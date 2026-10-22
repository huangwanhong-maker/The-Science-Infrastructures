"""A client for ECMWF's open data: the fields a manifest asks for, taken from the
per-step GRIB2 files of a forecast run by HTTP range requests.

The European Centre for Medium-Range Weather Forecasts publishes a subset of its
real-time forecasts, free and open, at https://data.ecmwf.int/forecasts/, with
copies on public clouds (Amazon S3, Azure, Google). Each run of a model (a
*cycle*, at 00, 06, 12 or 18 UTC) is a directory of GRIB2 files, one for each
step of the forecast::

    <root>/<yyyymmdd>/<HH>z/<model>/<resolution>/<stream>/<yyyymmddHH>0000-<step>h-<stream>-<type>.grib2

such as ``20261022/00z/ifs/0p25/oper/20261022000000-6h-oper-fc.grib2``. The
model is ``ifs`` (the physical model) or ``aifs-single`` (the AI model AIFS);
the stream ``oper`` (the medium-range forecast, formerly HRES) or ``enfo`` (the
ensemble); the type ``fc`` (a forecast), ``cf`` and ``pf`` (the control and
perturbed members) or ``ep`` (probabilities derived from the ensemble, with its
mean ``em`` and spread ``es``). The probabilities are not filed by step: one
file, named by step 240, holds every step up to 240 h, and another, named 360,
the rest.

Beside each file, an index (``.index``) lists its messages, one JSON object a
line, with the MARS keys and the message's place in the file::

    {"param": "2t", "levtype": "sfc", "step": "6", "type": "fc", "_offset": 23246183, "_length": 316751, ...}

with ``levelist`` (hPa) for pressure levels, and a range of hours (``"96-120"``)
as the step of a quantity over a period, such as a probability of 24-hour
precipitation. A file of the 0.25 degree forecast holds about 190 fields and
100 MB; the dozen fields an atlas draws are about 6 MB of it. The client reads
the index, selects the messages asked for, merges adjacent ones into runs of
bytes, and asks for all of them in one request (RFC 9110, section 14: several
ranges are answered as ``multipart/byteranges``). A server that answers only
single ranges, as Amazon S3 does, sends the whole file with status 200 instead;
the response is then closed unread and the runs are asked for one at a time.
Every message is checked by its framing (``GRIB`` ... ``7777``, edition 2, and
the length in section 0 equal to the index's) before it is written.

A manifest of kind ``ecmwf`` names the model, resolution, stream and type
(``source``), the hours of the runs to consider (``cycles``), and the fields
(``variables``: ``param``, ``levtype``, ``levelist``, and ``window`` for a
quantity over a period). A ``param`` may be a list of alternatives, the first in
the index being taken, since some streams give geopotential (``z``) where
others give geopotential height (``gh``). The run taken is the latest at one of
the ``cycles``, at or before the parameter ``cycle``, whose last file is
published: a run's files appear step by step over about an hour, so a run
without its last file is incomplete and the one before is taken. ECMWF keeps the
last 12 runs, about three days.

The fetched messages go to ``data/cache/<id>/<profile>/``, one file for each
file of the run (``fc-006h.grib2``, ``ep-240h.grib2``), with the selected index
entries, in file order, in the provenance record. The client needs only the
standard library; decoding is the product builder's (data/products/ecmwf_fields.py,
with ecCodes).

Licence: ECMWF's open data are published under CC BY 4.0 with ECMWF's Terms of
Use (https://apps.ecmwf.int/datasets/licences/general/): they may be used and
redistributed, commercially too, with the attribution "based on data and
products of the European Centre for Medium-Range Weather Forecasts (ECMWF)",
and an indication of any changes. ECMWF limits its portal to between 500 and
2000 simultaneous connections; this client makes one request at a time.

References: ECMWF, "ECMWF open data: real-time forecasts from IFS and AIFS"
(https://confluence.ecmwf.int/display/DAC/ECMWF+open+data%3A+real-time+forecasts+from+IFS+and+AIFS);
the ecmwf-opendata package (https://github.com/ecmwf/ecmwf-opendata), whose
index-and-range method this client follows.

Usage (the kind is registered in fetch.py's ACQUIRE table; this module also
runs alone, with fetch.py's arguments)::

    python data/acquisition/fetch.py ecmwf-ifs
    python data/acquisition/ecmwf_open.py ecmwf-aifs --set hours=48
    python data/acquisition/ecmwf_open.py ecmwf-ens --set cycle=2026-10-22T00:00
"""

from __future__ import annotations

import datetime as dt
import json
import re
import struct
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

UTC = dt.timezone.utc
USER_AGENT = "science-infrastructures-data/1 (+https://github.com/huangwanhong-maker/The-Science-Infrastructures)"
MAX_RANGES = 40        # byte ranges in one request; servers bound their number (Apache's MaxRanges: 200)
LOOKBACK = 6           # runs tried, latest first: at 00 and 12 UTC, three days, as long as ECMWF keeps them
CONTENT_RANGE = re.compile(rb"Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+|\*)", re.IGNORECASE)


class RangeError(OSError):
    """The server did not answer a range request with the bytes asked for."""


# ---- the layout of a run ------------------------------------------------------

def file_url(service: str, cycle: dt.datetime, source: dict, step: int) -> str:
    """The GRIB2 file of one step of a run (source: model, resolution, stream, type)."""
    s = source
    return (f"{service.rstrip('/')}/{cycle:%Y%m%d}/{cycle:%H}z/{s['model']}/{s['resolution']}/{s['stream']}/"
            f"{cycle:%Y%m%d%H}0000-{step}h-{s['stream']}-{s['type']}.grib2")


def index_url(url: str) -> str:
    return url[: -len(".grib2")] + ".index" if url.endswith(".grib2") else url + ".index"


def candidate_cycles(present: dt.datetime, hours=(0, 12), lookback: int = LOOKBACK) -> list[dt.datetime]:
    """The runs at the given hours (UTC) at or before present, latest first."""
    t = present.astimezone(UTC).replace(minute=0, second=0, microsecond=0)
    out = []
    while len(out) < lookback:
        if t.hour in hours:
            out.append(t)
        t -= dt.timedelta(hours=1)
    return out


def steps(first: int, last: int, every: int) -> list[int]:
    if every <= 0 or last < first:
        raise ValueError(f"no steps from {first} to {last} h every {every} h")
    return list(range(first, last + 1, every))


def file_step(step: int, file_steps: list[int] | None) -> int:
    """The step that names the file holding a step: itself, or for products filed
    together (the ensemble's probabilities) the first of file_steps at or after it."""
    if not file_steps:
        return step
    for s in sorted(file_steps):
        if s >= step:
            return s
    raise ValueError(f"no file holds step {step} h (files: {file_steps})")


def index_step(field: dict, step: int) -> str:
    """The step as the index writes it: '120', or for a quantity over a period of
    `window` hours ending at the step, '96-120'."""
    window = field.get("window")
    return f"{step - int(window)}-{step}" if window else str(step)


def signature(field: dict, kind: str) -> tuple:
    """What identifies a field in an index, whatever its alternatives."""
    params = field["param"] if isinstance(field["param"], list) else [field["param"]]
    return (tuple(params), field.get("levtype", "sfc"), str(field.get("levelist", "")), field.get("window"),
            field.get("type", kind))


# ---- indexes and byte ranges ----------------------------------------------------

def parse_index(text: str) -> list[dict]:
    """The entries of an index: one JSON object a line, each with _offset and _length."""
    out = []
    for n, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        entry = json.loads(line)
        if not isinstance(entry.get("_offset"), int) or not isinstance(entry.get("_length"), int):
            raise ValueError(f"index line {n} has no _offset or _length")
        out.append(entry)
    return out


def select(entries: list[dict], fields: list[dict], at: list[int], kind: str) -> list[dict]:
    """The index entries of each field at each step, in the order of the file.

    A field's param may be a list of alternatives (["gh", "z"]): the first that the
    index holds at a step is taken. Raises ValueError naming what is missing.
    """
    table = {}
    for e in entries:
        table[(e.get("param"), e.get("levtype"), str(e.get("levelist", "")), str(e.get("step")), e.get("type"))] = e
    chosen, missing = {}, []
    for field in fields:
        params, levtype, levelist, _, want_type = signature(field, kind)
        for step in at:
            label = index_step(field, step)
            hit = next((table[k] for p in params if (k := (p, levtype, levelist, label, want_type)) in table), None)
            if hit is None:
                where = f" at {levelist} hPa" if levelist else ""
                missing.append(f"{'/'.join(params)}{where} at step {label} h")
            else:
                chosen[hit["_offset"]] = hit
    if missing:
        raise ValueError("not in the index: " + ", ".join(missing))
    return [chosen[k] for k in sorted(chosen)]


def plan(entries: list[dict]) -> list[tuple[int, int]]:
    """Byte ranges (first, last), inclusive as HTTP counts them, that cover the
    entries, adjacent or overlapping messages merged into one range."""
    merged: list[list[int]] = []
    for a, b in sorted((e["_offset"], e["_offset"] + e["_length"]) for e in entries):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(a, b - 1) for a, b in merged]


def range_header(ranges: list[tuple[int, int]]) -> str:
    return "bytes=" + ",".join(f"{a}-{b}" for a, b in ranges)


def multipart_parts(body: bytes, boundary: str) -> dict[int, bytes]:
    """{first byte: bytes} of a multipart/byteranges body (RFC 9110, section 14.6).

    Each part is read by the length its Content-Range states, so that a boundary
    string inside the binary data cannot cut it short.
    """
    mark = b"--" + boundary.encode("ascii")
    parts, pos = {}, 0
    while True:
        at = body.find(mark, pos)
        if at < 0:
            raise RangeError("the multipart response ended without its closing boundary")
        pos = at + len(mark)
        if body[pos:pos + 2] == b"--":
            return parts
        head_end = body.find(b"\r\n\r\n", pos)
        if head_end < 0:
            raise RangeError("a part of the multipart response has no end of headers")
        m = CONTENT_RANGE.search(body, pos, head_end + 4)
        if not m:
            raise RangeError("a part of the multipart response has no Content-Range")
        first, last = int(m.group(1)), int(m.group(2))
        start = head_end + 4
        data = body[start:start + last - first + 1]
        if len(data) != last - first + 1:
            raise RangeError(f"the part for bytes {first}-{last} is cut short")
        parts[first] = data
        pos = start + len(data)


def grib_message(parts: dict[int, bytes], entry: dict) -> bytes:
    """One message, cut from the part that holds it and checked by its framing:
    'GRIB', edition 2, the total length of section 0 equal to the index's, '7777'."""
    offset, length = entry["_offset"], entry["_length"]
    for first, data in parts.items():
        if first <= offset and offset + length <= first + len(data):
            msg = data[offset - first: offset - first + length]
            break
    else:
        raise RangeError(f"no part of the response holds bytes {offset}-{offset + length - 1}")
    if msg[:4] != b"GRIB" or msg[-4:] != b"7777" or msg[7] != 2 or struct.unpack(">Q", msg[8:16])[0] != length:
        raise RangeError(f"bytes {offset}-{offset + length - 1} are not one GRIB2 message ({entry.get('param')})")
    return msg


# ---- requests --------------------------------------------------------------------

def _open(url: str, headers: dict | None = None, timeout: int = 120):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})}),
                                  timeout=timeout)


def retry(call, attempts: int = 3):
    """call(), again after a pause if the network fails; a missing file (HTTP 4xx)
    and a server that does not serve ranges are answers, not failures."""
    for attempt in range(1, attempts + 1):
        try:
            return call()
        except RangeError:
            raise
        except urllib.error.HTTPError as error:
            if error.code < 500 or attempt == attempts:
                raise
        except (urllib.error.URLError, OSError):
            if attempt == attempts:
                raise
        print(f"  attempt {attempt} failed; trying again")
        time.sleep(5 * attempt)


def read_text(url: str) -> str:
    with _open(url) as response:
        return response.read().decode("utf-8")


def read_ranges(url: str, ranges: list[tuple[int, int]], timeout: int = 300) -> dict[int, bytes]:
    """{first byte: bytes} for the ranges, in one request."""
    with _open(url, {"Range": range_header(ranges)}, timeout) as response:
        if response.status != 206:
            # The whole file is coming: close it unread.
            raise RangeError(f"the server answered {response.status} to a range request")
        content_type = response.headers.get("Content-Type", "")
        content_range = response.headers.get("Content-Range", "")
        body = response.read()
    m = re.match(r'multipart/byteranges;\s*boundary="?([^";]+)"?', content_type, re.IGNORECASE)
    if m:
        return multipart_parts(body, m.group(1))
    m = CONTENT_RANGE.match(f"Content-Range: {content_range}".encode())
    if not m:
        raise RangeError("a single-part answer without Content-Range")
    return {int(m.group(1)): body}


def fetch_messages(url: str, entries: list[dict]) -> tuple[bytes, int]:
    """The entries' messages, in file order, and the number of requests made."""
    ranges = plan(entries)
    parts, requests = {}, 0
    for i in range(0, len(ranges), MAX_RANGES):
        batch = ranges[i:i + MAX_RANGES]
        try:
            parts.update(retry(lambda: read_ranges(url, batch)))
            requests += 1
        except RangeError:
            if len(batch) == 1:
                raise
            for r in batch:                    # a server of single ranges only
                parts.update(retry(lambda: read_ranges(url, [r])))
                requests += 1
    return b"".join(grib_message(parts, e) for e in entries), requests


def latest_run(service: str, source: dict, cycles, present: dt.datetime, last_file: int) -> tuple[dt.datetime, str]:
    """The latest run whose last file is published, and that file's index."""
    tried = []
    for cycle in candidate_cycles(present, tuple(cycles)):
        url = index_url(file_url(service, cycle, source, last_file))
        try:
            return cycle, retry(lambda: read_text(url))
        except urllib.error.HTTPError as error:
            if error.code not in (403, 404):
                raise
            tried.append(f"{cycle:%Y-%m-%d %H}z")
    raise FileNotFoundError(f"no run with a file at {last_file} h among {', '.join(tried)}")


# ---- the kind 'ecmwf' of fetch.py ---------------------------------------------------

def acquire(manifest: dict, profile: str | None, overrides: dict, refresh: bool) -> bool:
    """The fields of a manifest of kind 'ecmwf' from the latest complete run,
    one file of messages for each file of the run."""
    import fetch                                # imported here: fetch.py imports this module
    values, tag = fetch.resolve(manifest, profile, overrides)
    folder = fetch.CACHE / manifest["id"] / tag
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    source, kind = manifest["source"], manifest["source"]["type"]
    fields = manifest["variables"]
    try:
        at = steps(int(manifest.get("first_step", 0)), int(values["hours"]), int(values["step_hours"]))
        by_file: dict[int, list[int]] = {}
        for s in at:
            by_file.setdefault(file_step(s, manifest.get("file_steps")), []).append(s)
    except ValueError as error:
        print(f"  {error}")
        return False
    wanted = json.loads(json.dumps([signature(f, kind) for f in fields]))     # as the provenance stores it
    present = fetch.resolve_time(values["cycle"], dt.datetime.now(UTC))

    # The run: the first server that answers says which is the latest complete one.
    cycle = service = None
    for service in [manifest["service"]] + [m["service"] for m in manifest.get("mirrors", [])]:
        try:
            cycle, last_index = latest_run(service, source, manifest.get("cycles", [0, 12]), present, max(by_file))
            break
        except (FileNotFoundError, urllib.error.URLError, OSError) as error:
            print(f"  {service} could not serve the run ({error})")
    if cycle is None:
        if (folder / "retrieval.json").exists():
            print("  keeping the cached retrieval")
            return True
        print("  REQUEST FAILED: no server has a complete run")
        return False
    print(f"  run {cycle:%Y-%m-%d %H:%M} UTC from {service}: {len(at)} step(s), {at[0]} to {at[-1]} h, "
          f"in {len(by_file)} file(s)")

    resolved = {"cycle": f"{cycle:%Y-%m-%dT%H:%M:%SZ}", "service": service, "steps": at, "fields": wanted}
    provenance = folder / "retrieval.json"
    if provenance.exists() and not refresh:
        old = fetch.load(provenance)
        same = {k: old.get("resolved", {}).get(k) for k in ("cycle", "steps", "fields")} == \
            {k: resolved[k] for k in ("cycle", "steps", "fields")}
        if same and all((folder / r["file"]).exists() for r in old.get("requests", [])):
            print(f"  cached in {folder.relative_to(fetch.DATA.parent)}")
            return True

    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    requests, total = [], 0
    for fstep, held in sorted(by_file.items()):
        url = file_url(service, cycle, source, fstep)
        name = f"{kind}-{fstep:03d}h.grib2"
        try:
            text = last_index if fstep == max(by_file) else retry(lambda: read_text(index_url(url)))
            entries = select(parse_index(text), fields, held, kind)
            data, count = fetch_messages(url, entries)
        except (ValueError, urllib.error.URLError, OSError) as error:
            print(f"  REQUEST FAILED: {url}: {error}")
            return False
        part = folder / (name + ".part")
        part.write_bytes(data)
        part.replace(folder / name)
        total += len(data)
        requests.append({"url": url, "file": name, "index": index_url(url), "steps": held, "requests": count,
                         "messages": [{k: e.get(k) for k in ("param", "levtype", "levelist", "step", "type")}
                                      | {"offset": e["_offset"], "length": e["_length"]} for e in entries]})
        print(f"  {name}: {len(entries)} field(s), {len(data) / 1e6:.1f} MB in {count} request(s)")
    fetch.write_provenance(folder, manifest, tag, values, resolved, requests, started)
    print(f"  {len(requests)} file(s), {total / 1e6:.1f} MB in {folder.relative_to(fetch.DATA.parent)}")
    return True


def main(argv: list[str]) -> int:
    """fetch.py's command line, with this kind registered."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import fetch
    fetch.ACQUIRE.setdefault("ecmwf", acquire)
    return fetch.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
