"""Acquire registered datasets into the local cache.

A dataset is described by a manifest in ``data/registry/`` (schema:
``data/registry/README.md``). There are six kinds:

``files``
    Fixed files, downloaded as they are, such as the DE440 kernel. Their
    checksums are recorded in the manifest and verified on every retrieval.
    A file may be a member of a zip archive on the web, read by range
    requests without downloading the archive.
``tap``
    A query to an astronomical archive through the Table Access Protocol, such
    as Gaia.
``griddap``
    A subset of a gridded dataset served by ERDDAP, such as a weather forecast
    or a sea-surface temperature analysis: which variables, which times, and
    how fine a grid.
``http``
    One or more web requests built from URL templates, such as the World
    Bank's indicators or the USGS earthquake catalogue.
``ncss``
    Subsets of a forecast model's collection on a THREDDS server, through its
    NetCDF Subset Service: groups of variables at chosen levels, times and
    grid spacing, such as the GFS at 500 hPa or NOAA's AI model AIGFS.
``opendap``
    Subsets of gridded datasets on an OPeNDAP (DAP2) server, by index ranges
    and strides, saved as netCDF classic files, such as TerraClimate's monthly
    climate normals.

Every kind except ``files`` declares parameters (which part of the data and how
much of it) and named profiles. A retrieval chooses a profile and may override
any parameter. Parameters of type ``time`` accept an ISO 8601 instant (UTC if
no offset is given) or an expression relative to the present or to the end of
the dataset: ``now``, ``today`` or ``last``, optionally followed by ``+`` or
``-`` and a number of hours (``h``) or days (``d``), as in ``now-1d``.

Files go to ``data/cache/<id>/``, or ``data/cache/<id>/<profile>/`` for the
parameterised kinds, with a provenance record, ``retrieval.json``, beside them.
The cache is not under version control.

Usage::

    python data/acquisition/fetch.py --list                        # datasets, parameters, profiles
    python data/acquisition/fetch.py                               # every 'files' dataset
    python data/acquisition/fetch.py jpl-de440s hyg-v41            # selected datasets
    python data/acquisition/fetch.py gaia-dr3                      # a query, default profile
    python data/acquisition/fetch.py gaia-dr3 --profile nearby
    python data/acquisition/fetch.py gaia-dr3 --set limit=50000 --set g_max=9.5
    python data/acquisition/fetch.py gfs --set end=now+3d --set stride=1
    python data/acquisition/fetch.py --record                      # record checksums of a first retrieval
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
REGISTRY = DATA / "registry"
CACHE = DATA / "cache"
CHUNK = 1 << 20

sys.path.insert(0, str(DATA / "catalogs"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import erddap  # noqa: E402
import ncss  # noqa: E402
import remotezip  # noqa: E402
import tap  # noqa: E402

UTC = dt.timezone.utc
RELATIVE_TIME = re.compile(r"^(now|today|last)(?:\s*([+-])\s*(\d+(?:\.\d+)?)\s*([hd]))?$")


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


# ---- fixed files ------------------------------------------------------------

def download(url: str, dest: Path, timeout: int = 60, attempts: int = 1) -> None:
    """Downloads url to dest, through a .part file; retries on failure, resuming
    where the last attempt stopped when the server allows it."""
    for attempt in range(1, attempts + 1):
        try:
            _download(url, dest, timeout)
            return
        except (urllib.error.URLError, OSError) as error:
            if attempt == attempts:
                raise
            print(f"  attempt {attempt} failed ({error}); trying again")
            time.sleep(10 * attempt)


def _download(url: str, dest: Path, timeout: int) -> None:
    """One attempt. A .part file left by an earlier attempt is continued by a
    range request (RFC 9110, section 14) if the server answers it (206), and
    begun again otherwise. A response shorter than its Content-Length is an
    error: the connection ended early, and the .part file is kept for the next
    attempt."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    headers = {"User-Agent": "laboratory-data-fetch/1 (+https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory)"}
    if have:
        headers["Range"] = f"bytes={have}-"
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=timeout) as response:
        resumed = have and response.status == 206
        if not resumed:
            have = 0
        total = have + int(response.headers.get("Content-Length") or 0)
        done = have
        with part.open("ab" if resumed else "wb") as out:
            for block in iter(lambda: response.read(CHUNK), b""):
                out.write(block)
                done += len(block)
                if total:
                    print(f"\r  {dest.name}: {done / total:6.1%}", end="", flush=True)
    print()
    if total and done != total:
        raise OSError(f"{dest.name}: the transfer ended at {done} of {total} bytes")
    part.replace(dest)


def extract_member(url: str, member: str, dest: Path, archives: dict) -> None:
    """Writes one member of a remote zip archive to dest, reading only its bytes."""
    if url not in archives:
        archives[url] = zipfile.ZipFile(remotezip.RemoteFile(url))
    data = archives[url].read(member)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    part.write_bytes(data)
    part.replace(dest)
    print(f"  {dest.name}: {len(data) / 1e6:.2f} MB from {member}")


def acquire_files(path: Path, manifest: dict, record: bool) -> bool:
    ident = manifest["id"]
    ok, changed, retrieved = True, False, []
    archives = {}
    for entry in manifest["files"]:
        dest = CACHE / ident / entry["file"]
        expected = entry.get("sha256")
        if dest.exists() and (expected is None or sha256_of(dest) == expected):
            print(f"  {entry['file']}: cached")
        elif "member" in entry:
            extract_member(entry["url"], entry["member"], dest, archives)
        else:
            download(entry["url"], dest, timeout=120, attempts=5)
        digest, size = sha256_of(dest), dest.stat().st_size
        if expected is not None and digest != expected:
            print(f"  {entry['file']}: CHECKSUM MISMATCH (expected {expected}, got {digest})")
            ok = False
        elif expected is None:
            if record:
                entry["sha256"], entry["size"] = digest, size
                changed = True
                print(f"  {entry['file']}: checksum recorded")
            else:
                print(f"  {entry['file']}: unverified (no recorded checksum); run with --record")
        retrieved.append({"file": entry["file"], "url": entry["url"], "size": size, "sha256": digest}
                         | ({"member": entry["member"]} if "member" in entry else {}))
    (CACHE / ident / "retrieval.json").write_text(json.dumps(
        {"dataset": ident, "retrieved_at": now(), "files": retrieved}, indent=2) + "\n", encoding="utf-8")
    if changed:
        path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return ok


# ---- queries ----------------------------------------------------------------

def resolve_time(expr: str, present: dt.datetime, last: dt.datetime | None = None) -> dt.datetime:
    """An instant (UTC) from an ISO 8601 text or a relative expression.

    'now' is the present rounded down to the hour, so that retrievals within the
    same hour ask for the same data; 'today' is 00:00 UTC of the present day;
    'last' is the last time of the dataset. A text without an offset is UTC.
    """
    m = RELATIVE_TIME.match(expr.strip().lower())
    if m:
        anchor, sign, amount, unit = m.groups()
        if anchor == "now":
            base = present.replace(minute=0, second=0, microsecond=0)
        elif anchor == "today":
            base = present.replace(hour=0, minute=0, second=0, microsecond=0)
        elif last is None:
            raise ValueError("'last' refers to the end of a dataset, which is not known here")
        else:
            base = last
        if sign:
            delta = dt.timedelta(hours=float(amount) * (24 if unit == "d" else 1))
            base = base + delta if sign == "+" else base - delta
        return base
    t = dt.datetime.fromisoformat(expr.strip().replace("Z", "+00:00"))
    return (t if t.tzinfo else t.replace(tzinfo=UTC)).astimezone(UTC)


def coerce(name: str, spec: dict, raw):
    kind = spec.get("type", "str")
    if kind == "time":
        try:
            epoch = dt.datetime(2000, 1, 1, tzinfo=UTC)
            resolve_time(str(raw), epoch, epoch)
        except ValueError:
            raise SystemExit(f"parameter {name}: expected an ISO 8601 time or an expression such as "
                             f"now-1d, today or last, got {raw!r}") from None
        return str(raw)
    try:
        value = int(raw) if kind == "int" else float(raw) if kind == "float" else str(raw)
    except ValueError:
        raise SystemExit(f"parameter {name}: expected {kind}, got {raw!r}") from None
    if "min" in spec and value < spec["min"] or "max" in spec and value > spec["max"]:
        raise SystemExit(f"parameter {name}: {value} is outside [{spec.get('min')}, {spec.get('max')}]")
    # A text parameter may be held to a pattern, as one that names a file on a server must be.
    if kind == "str" and spec.get("pattern") and not re.fullmatch(spec["pattern"], value):
        raise SystemExit(f"parameter {name}: {value!r} does not match {spec['pattern']}")
    return value


def resolve(manifest: dict, profile: str | None, overrides: dict) -> tuple[dict, str]:
    """Parameter values: defaults, then the profile, then overrides; and the cache tag."""
    spec = manifest["parameters"]
    profile = profile or manifest.get("default_profile")
    if profile and profile not in manifest.get("profiles", {}):
        raise SystemExit(f"{manifest['id']}: unknown profile {profile!r}; "
                         f"choose from {', '.join(manifest.get('profiles', {}))}")
    values = {k: coerce(k, s, s["default"]) for k, s in spec.items()}
    if profile:
        for k, v in manifest["profiles"][profile]["parameters"].items():
            values[k] = coerce(k, spec[k], v)
    for k, raw in overrides.items():
        if k not in spec:
            raise SystemExit(f"{manifest['id']}: unknown parameter {k!r}; choose from {', '.join(spec)}")
        values[k] = coerce(k, spec[k], raw)
    if overrides:
        digest = hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()[:8]
        tag = f"custom-{digest}"
    else:
        tag = profile or "default"
    return values, tag


def acquire_tap(manifest: dict, profile: str | None, overrides: dict, refresh: bool) -> bool:
    values, tag = resolve(manifest, profile, overrides)
    folder = CACHE / manifest["id"] / tag
    dest = folder / manifest.get("file", "result.csv")
    provenance = folder / "retrieval.json"
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    if dest.exists() and provenance.exists() and not refresh:
        if load(provenance).get("parameters") == values:
            print(f"  cached in {folder.relative_to(DATA.parent)}")
            return True
    adql = "\n".join(manifest["query"]).format(**values)
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    maxrec = values.get(manifest.get("maxrec_parameter", "limit"))
    try:
        size = tap.run_async(manifest["service"], adql, dest, fmt="csv", maxrec=maxrec)
    except tap.TapError as error:
        print(f"  QUERY FAILED: {error}")
        return False
    with dest.open("rb") as handle:
        rows = max(0, sum(1 for _ in handle) - 1)
    if maxrec and rows < maxrec and rows % 1000 == 0:
        print(f"  note: {rows} rows for a limit of {maxrec}; the service may cap results")
    provenance.write_text(json.dumps({
        "dataset": manifest["id"],
        "profile": tag,
        "parameters": values,
        "service": manifest["service"],
        "query": adql,
        "retrieved_at": now(),
        "elapsed_s": round(time.monotonic() - started, 1),
        "rows": rows,
        "file": dest.name,
        "size": size,
        "sha256": sha256_of(dest),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"  {rows} rows, {size / 1e6:.1f} MB in {folder.relative_to(DATA.parent)}")
    return True


# ---- gridded subsets (ERDDAP) and web requests ------------------------------

def fill(template, values: dict):
    """Fills '{name}' placeholders; a template that is one placeholder keeps the value's type."""
    if not isinstance(template, str):
        return template
    whole = re.fullmatch(r"\{(\w+)\}", template)
    return values[whole.group(1)] if whole else template.format(**values)


def cached(folder: Path, urls: list[str], files: list[str]) -> bool:
    provenance = folder / "retrieval.json"
    if not provenance.exists() or not all((folder / f).exists() for f in files):
        return False
    return [r["url"] for r in load(provenance).get("requests", [])] == urls


def write_provenance(folder: Path, manifest: dict, tag: str, values: dict, resolved: dict,
                     requests: list[dict], started: float) -> None:
    for r in requests:
        path = folder / r["file"]
        r["size"], r["sha256"] = path.stat().st_size, sha256_of(path)
    (folder / "retrieval.json").write_text(json.dumps({
        "dataset": manifest["id"],
        "profile": tag,
        "parameters": values,
        "resolved": resolved,
        "requests": requests,
        "retrieved_at": now(),
        "elapsed_s": round(time.monotonic() - started, 1),
    }, indent=2) + "\n", encoding="utf-8")


def griddap_request(manifest: dict, values: dict, present: dt.datetime) -> tuple[str, dict]:
    """The request URL for a griddap dataset, and the resolved axes for provenance."""
    service, dataset = manifest["service"], manifest["dataset"]
    meta = erddap.info(service, dataset)
    variables = [v.strip() for v in str(fill(manifest["variables"], values)).split(",") if v.strip()]
    unknown = [v for v in variables if v not in meta["variables"]]
    if unknown:
        raise erddap.ErddapError(f"unknown variable(s) {', '.join(unknown)}; "
                                 f"the dataset has {', '.join(meta['variables'])}")
    axes, resolved = [], {"variables": variables}
    for spec in manifest["axes"]:
        spec = {k: fill(v, values) for k, v in spec.items()}
        if spec["name"] != "time":
            axes.append(spec)
            resolved[spec["name"]] = {k: spec[k] for k in ("from", "to", "stride") if k in spec}
            continue
        first = erddap.parse_time(meta["global"]["time_coverage_start"])
        last = erddap.parse_time(meta["global"]["time_coverage_end"])
        start = resolve_time(str(spec["from"]), present, last)
        stop = resolve_time(str(spec.get("to", spec["from"])), present, last)
        start, stop = max(start, first), min(stop, last)
        if start > stop:
            raise erddap.ErddapError(f"the requested times lie outside the dataset "
                                     f"({erddap.iso(first)} to {erddap.iso(last)})")
        step = int(spec.get("step_hours", 1))
        times = erddap.time_values(service, dataset, start, stop)
        t0, t1, stride, count = erddap.select_times(times, step)
        axes.append({"name": "time", "from": erddap.iso(t0), "to": erddap.iso(t1), "stride": stride})
        resolved["time"] = {"first": erddap.iso(t0), "last": erddap.iso(t1), "count": count,
                            "step_hours": step, "dataset_end": erddap.iso(last)}
    return erddap.request_url(service, dataset, variables, axes), resolved


def acquire_griddap(manifest: dict, profile: str | None, overrides: dict, refresh: bool) -> bool:
    values, tag = resolve(manifest, profile, overrides)
    folder = CACHE / manifest["id"] / tag
    name = manifest.get("file", "data.nc")
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    started = time.monotonic()
    try:
        url, resolved = griddap_request(manifest, values, dt.datetime.now(UTC))
    except (erddap.ErddapError, urllib.error.URLError, OSError) as error:
        if (folder / name).exists():
            print(f"  the service could not be reached ({error}); keeping the cached retrieval")
            return True
        print(f"  REQUEST FAILED: {error}")
        return False
    if "time" in resolved:
        t = resolved["time"]
        print(f"  times {t['first']} to {t['last']}: {t['count']} step(s) of {t['step_hours']} h")
    if cached(folder, [url], [name]) and not refresh:
        print(f"  cached in {folder.relative_to(DATA.parent)}")
        return True
    folder.mkdir(parents=True, exist_ok=True)
    try:
        # ERDDAP may compute a large subset for minutes before sending a byte.
        download(url, folder / name, timeout=900, attempts=2)
    except (urllib.error.URLError, OSError) as error:
        print(f"  DOWNLOAD FAILED: {error}")
        return False
    write_provenance(folder, manifest, tag, values, resolved, [{"url": url, "file": name}], started)
    print(f"  {(folder / name).stat().st_size / 1e6:.1f} MB in {folder.relative_to(DATA.parent)}")
    return True


def followed(folder: Path, rule: dict) -> list[dict]:
    """Requests named by the items of a JSON file already fetched.

    rule: {"from": file, "items": key of the list in it, "where": {dotted path:
    value}, "url": template, "file": template}; templates index an item as
    str.format does, as in "{properties[url][geometry]}".
    """
    data = load(folder / rule["from"])
    out = []
    for item in data.get(rule["items"], []):
        def get(path, node=item):
            for key in path.split("."):
                node = node.get(key) if isinstance(node, dict) else None
            return node
        if all(get(k) == v for k, v in rule.get("where", {}).items()):
            url, file = rule["url"].format(**item), rule["file"].format(**item)
            # Names come from fetched content: keep them in the folder, and the requests secure.
            if Path(file).name != file or file.startswith(".") or not url.startswith("https://"):
                print(f"  linked request refused: {url} -> {file}")
                continue
            out.append({"url": url, "file": file})
    return out


def acquire_http(manifest: dict, profile: str | None, overrides: dict, refresh: bool) -> bool:
    values, tag = resolve(manifest, profile, overrides)
    folder = CACHE / manifest["id"] / tag
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    present = dt.datetime.now(UTC)
    filled = {k: resolve_time(v, present) if manifest["parameters"][k].get("type") == "time" else v
              for k, v in values.items()}
    # A time parameter with align_hours is moved down to a whole multiple of that
    # many hours since 00:00 UTC, as the cycles of a forecast model are.
    for k, spec in manifest["parameters"].items():
        if spec.get("type") == "time" and spec.get("align_hours"):
            filled[k] = ncss.align(filled[k], int(spec["align_hours"]))
    requests = [{"url": r["url"].format(**filled), "file": r["file"]} for r in manifest["requests"]]
    resolved = {k: erddap.iso(v) for k, v in filled.items() if isinstance(v, dt.datetime)}
    if cached(folder, [r["url"] for r in requests], [r["file"] for r in requests]) and not refresh:
        print(f"  cached in {folder.relative_to(DATA.parent)}")
        return True
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    for r in requests:
        try:
            download(r["url"], folder / r["file"], timeout=300, attempts=2)
        except (urllib.error.URLError, OSError) as error:
            print(f"  REQUEST FAILED: {r['url']}: {error}")
            return False
    for rule in manifest.get("follow", []):
        for r in followed(folder, rule):
            try:
                download(r["url"], folder / r["file"], timeout=300, attempts=2)
                requests.append(r)
            except (urllib.error.URLError, OSError) as error:
                print(f"  linked request failed, skipped: {r['url']}: {error}")
    write_provenance(folder, manifest, tag, values, resolved, requests, started)
    total = sum((folder / r["file"]).stat().st_size for r in requests)
    print(f"  {len(requests)} file(s), {total / 1e6:.2f} MB in {folder.relative_to(DATA.parent)}")
    return True


def acquire_ncss(manifest: dict, profile: str | None, overrides: dict, refresh: bool) -> bool:
    """Groups of variables from a THREDDS collection, one request (and file) per group.

    Times run from start to end every step_hours, aligned to 00:00 UTC; the
    collection's own step (time_step_hours in the manifest) sets the stride.
    """
    values, tag = resolve(manifest, profile, overrides)
    folder = CACHE / manifest["id"] / tag
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    present = dt.datetime.now(UTC)
    step = int(values.get("step_hours", manifest.get("time_step_hours", 3)))
    native = int(manifest.get("time_step_hours", 3))
    if step % native:
        print(f"  step_hours must be a multiple of {native}")
        return False
    start = ncss.align(resolve_time(values["start"], present), step, up=True)
    end = ncss.align(resolve_time(values["end"], present), step)
    stride = int(values.get("stride", 1))
    requests = []
    for group in manifest["groups"]:
        url = ncss.request_url(manifest["service"], manifest["dataset"], group["variables"], start, end,
                               step // native, stride, group.get("vertCoord"))
        requests.append({"url": url, "file": group["file"], "level": group.get("level")})
    resolved = {"start": erddap.iso(start), "end": erddap.iso(end), "step_hours": step}
    print(f"  times {resolved['start']} to {resolved['end']} every {step} h; {len(requests)} request(s)")
    if cached(folder, [r["url"] for r in requests], [r["file"] for r in requests]) and not refresh:
        print(f"  cached in {folder.relative_to(DATA.parent)}")
        return True
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    for r in requests:
        print(f"  {r['file']}: requesting (the server may take minutes to cut the subset)")
        try:
            download(r["url"], folder / r["file"], timeout=900, attempts=3)
        except (urllib.error.URLError, OSError) as error:
            print(f"  REQUEST FAILED: {r['file']}: {error}")
            return False
    write_provenance(folder, manifest, tag, values, resolved, requests, started)
    total = sum((folder / r["file"]).stat().st_size for r in requests)
    print(f"  {len(requests)} file(s), {total / 1e6:.1f} MB in {folder.relative_to(DATA.parent)}")
    return True


def column_indices(axis: dict, values: dict, present: dt.datetime) -> dict:
    """Indices of the grid point nearest (lat, lon), of the depths down to max_depth,
    and of the times from start to end on multiples of step_hours (UTC)."""
    import numpy as np
    lat, lon = float(values["lat"]), float(values["lon"])
    i_lat = int(np.argmin(np.abs(axis["lat"] - lat)))
    # Longitudes compared on the circle, whatever the axis's convention (0..360 or -180..180).
    gap = np.abs(((axis["lon"] - lon) + 180.0) % 360.0 - 180.0)
    i_lon = int(np.argmin(gap))
    depths = [i for i, z in enumerate(axis["depth"]) if z <= float(values.get("max_depth", 1e9)) + 1e-6]
    step = int(values.get("step_hours", 3))
    start = ncss.align(resolve_time(values["start"], present), step)
    end = ncss.align(resolve_time(values["end"], present), step)
    lo, hi = start.timestamp(), end.timestamp()
    times = [i for i, t in enumerate(axis["time"])
             if np.isfinite(t) and lo - 1 <= t <= hi + 1 and round(t) % (step * 3600) == 0]
    return {"lat": i_lat, "lon": i_lon, "depth": depths, "time": times,
            "start": start.isoformat().replace("+00:00", "Z"), "end": end.isoformat().replace("+00:00", "Z")}


def acquire_opendap_column(manifest: dict, values: dict, tag: str, refresh: bool) -> bool:
    """A column (time x depth) at the grid point nearest a place, for one variable of
    each of several OPeNDAP datasets, each saved as a netCDF classic file.

    The manifest's `column` names the coordinates (lat, lon, depth, time, and
    optionally run, the forecast reference time of each time). Times are asked
    for in short runs, several at once (`workers`), which suits servers that
    open one file per time step.
    """
    import numpy as np
    import opendap
    netcdf3 = opendap.netcdf3
    folder = CACHE / manifest["id"] / tag
    names = manifest["column"]
    base = manifest["service"].rstrip("/")
    present = dt.datetime.now(UTC)
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    try:
        axes, requests = {}, []
        for item in manifest["datasets"]:
            url = f"{base}/{item['path']}"
            if url not in axes:
                attrs = opendap.das(url)
                axis = {k: opendap.coordinate(url, names[k]) for k in ("lat", "lon", "depth")}
                axis["time"] = opendap.cf_epoch_seconds(opendap.coordinate(url, names["time"]), attrs[names["time"]]["units"])
                if names.get("run"):
                    axis["run"] = opendap.cf_epoch_seconds(opendap.coordinate(url, names["run"]), attrs[names["run"]]["units"])
                axes[url] = (axis, attrs, column_indices(axis, values, present))
            axis, attrs, idx = axes[url]
            if not idx["time"]:
                raise opendap.DapError(f"no times of {item['path']} between {idx['start']} and {idx['end']}")
            depth = f"[{idx['depth'][0]}:1:{idx['depth'][-1]}]"
            point = f"[{idx['lat']}][{idx['lon']}]"
            constraints = [f"[{a}:{s}:{b}]{depth}{point}" for a, s, b in opendap.runs(idx["time"], 4)]
            requests.append({"url": opendap.request_url(url, item["variable"], constraints[0]), "dataset": url,
                             "variable": item["variable"], "file": item["file"], "constraints": constraints})
    except (opendap.DapError, KeyError, urllib.error.URLError, OSError) as error:
        print(f"  REQUEST FAILED: {error}")
        return False
    if cached(folder, [r["url"] for r in requests], [r["file"] for r in requests]) and not refresh:
        print(f"  cached in {folder.relative_to(DATA.parent)}")
        return True
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    # Every variable's runs of times in one pool, so that all of them overlap.
    tasks = [(r["dataset"], r["variable"], c) for r in requests for c in r["constraints"]]
    print(f"  {len(tasks)} requests for {len(requests)} variable(s), {int(manifest.get('workers', 4))} at a time")
    try:
        results = iter(opendap.fetch_slabs(tasks, int(manifest.get("workers", 4))))
    except (opendap.DapError, urllib.error.URLError, OSError) as error:
        print(f"  REQUEST FAILED: {error}")
        return False
    for r in requests:
        axis, attrs, idx = axes[r["dataset"]]
        parts = [next(results) for _ in r.pop("constraints")]
        data = np.concatenate([p.reshape(p.shape[0], p.shape[1]) for p in parts]).astype(np.float32)
        va = {k: v for k, v in attrs.get(r["variable"], {}).items() if k not in ("coordinates", "NAVO_code") and not k.startswith("_")}
        for key in ("_FillValue", "missing_value"):
            fill = attrs.get(r["variable"], {}).get(key)
            if fill is not None and np.isfinite(fill):
                data[data == np.float32(fill)] = np.nan
        times = axis["time"][idx["time"]]
        variables = {
            "time": (["time"], times, {"units": "seconds since 1970-01-01 00:00:00 UTC", "standard_name": "time"}),
            "depth": (["depth"], axis["depth"][idx["depth"]], {"units": "m", "positive": "down", "standard_name": "depth"}),
            r["variable"]: (["time", "depth"], data, va),
        }
        if "run" in axis:
            variables["run"] = (["time"], axis["run"][idx["time"]], {"units": "seconds since 1970-01-01 00:00:00 UTC",
                                                                      "standard_name": "forecast_reference_time"})
        glob = {"source": f"{r['dataset']} (OPeNDAP), {r['variable']} at one grid point",
                "latitude": float(axis["lat"][idx["lat"]]), "longitude": float(axis["lon"][idx["lon"]]),
                "requested_latitude": float(values["lat"]), "requested_longitude": float(values["lon"])}
        netcdf3.write(folder / r["file"], {"time": len(times), "depth": len(idx["depth"])}, variables, glob)
        r["times"] = len(times)
    first = next(iter(axes.values()))
    resolved = {"start": first[2]["start"], "end": first[2]["end"],
                "grid_latitude": float(first[0]["lat"][first[2]["lat"]]),
                "grid_longitude": float(first[0]["lon"][first[2]["lon"]])}
    write_provenance(folder, manifest, tag, values, resolved, requests, started)
    print(f"  {len(requests)} file(s) in {time.monotonic() - started:.0f} s, in {folder.relative_to(DATA.parent)}")
    return True


def acquire_opendap(manifest: dict, profile: str | None, overrides: dict, refresh: bool) -> bool:
    """One variable from each of several OPeNDAP datasets, every stride-th point
    of the horizontal grid, each saved as a netCDF classic file (needs numpy);
    or, with a `column`, the column beneath one place (acquire_opendap_column)."""
    import opendap
    values, tag = resolve(manifest, profile, overrides)
    if manifest.get("column"):
        return acquire_opendap_column(manifest, values, tag, refresh)
    folder = CACHE / manifest["id"] / tag
    stride = int(values.get("stride", 1))
    print(f"  profile {tag}: " + ", ".join(f"{k}={v}" for k, v in values.items()))
    base = manifest["service"].rstrip("/")
    try:
        requests = []
        for item in manifest["datasets"]:
            url = f"{base}/{item['path']}"
            _, dims = opendap.dds(url)[item["variable"]]
            requests.append({"url": opendap.request_url(url, item["variable"], opendap.hyperslab(dims, stride)),
                             "dataset": url, "variable": item["variable"], "file": item["file"]})
    except (opendap.DapError, KeyError, urllib.error.URLError, OSError) as error:
        if (folder / "retrieval.json").exists():
            print(f"  the service could not be reached ({error}); keeping the cached retrieval")
            return True
        print(f"  REQUEST FAILED: {error}")
        return False
    if cached(folder, [r["url"] for r in requests], [r["file"] for r in requests]) and not refresh:
        print(f"  cached in {folder.relative_to(DATA.parent)}")
        return True
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    for r in requests:
        print(f"  {r['file']}: requesting {r['variable']} every {stride} point(s)")
        try:
            info = opendap.fetch_to_netcdf(r["dataset"], r["variable"], stride, folder / r["file"])
        except (opendap.DapError, urllib.error.URLError, OSError) as error:
            print(f"  REQUEST FAILED: {r['file']}: {error}")
            return False
        r["dimensions"] = info["dimensions"]
    write_provenance(folder, manifest, tag, values, {"stride": stride}, requests, started)
    total = sum((folder / r["file"]).stat().st_size for r in requests)
    print(f"  {len(requests)} file(s), {total / 1e6:.1f} MB in {folder.relative_to(DATA.parent)}")
    return True


ACQUIRE = {"tap": acquire_tap, "griddap": acquire_griddap, "http": acquire_http, "ncss": acquire_ncss,
           "opendap": acquire_opendap}


# ---- listing ----------------------------------------------------------------

def describe(manifest: dict) -> None:
    kind = manifest.get("kind", "files")
    print(f"{manifest['id']:30s} {kind:8s} {manifest['title']}")
    if kind == "files":
        total = sum(f.get("size") or 0 for f in manifest["files"])
        print(f"{'':39s}{len(manifest['files'])} file(s), {total / 1e6:.1f} MB")
        return
    print(f"{'':39s}parameters:")
    for name, spec in manifest["parameters"].items():
        bounds = f" {spec.get('min')}..{spec.get('max')}" if "min" in spec else ""
        print(f"{'':41s}{name} = {spec['default']} ({spec.get('type', 'str')}{bounds}) {spec['description']}")
    print(f"{'':39s}profiles:")
    for name, profile in manifest.get("profiles", {}).items():
        mark = " (default)" if name == manifest.get("default_profile") else ""
        print(f"{'':41s}{name}{mark}: {profile['description']}")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="*", help="dataset ids (default: every 'files' dataset)")
    parser.add_argument("--list", action="store_true", help="list datasets, their parameters and profiles")
    parser.add_argument("--profile", help="profile of a parameterised dataset")
    parser.add_argument("--set", action="append", default=[], metavar="NAME=VALUE",
                        help="override a parameter of a parameterised dataset (repeatable)")
    parser.add_argument("--refresh", action="store_true", help="query again even if a cached result exists")
    parser.add_argument("--record", action="store_true", help="record checksums of first retrieval")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace")     # consoles that cannot show ′, ° or CJK

    manifests = {p.stem: p for p in sorted(REGISTRY.glob("*.json"))}
    if args.list:
        for path in manifests.values():
            describe(load(path))
        return 0

    overrides = {}
    for item in args.set:
        name, sep, value = item.partition("=")
        if not sep:
            raise SystemExit(f"--set expects NAME=VALUE, got {item!r}")
        overrides[name.strip()] = value.strip()

    if args.datasets:
        unknown = [d for d in args.datasets if d not in manifests]
        if unknown:
            print(f"unknown dataset(s): {', '.join(unknown)}; see --list", file=sys.stderr)
            return 2
        chosen = args.datasets
    else:
        chosen = [k for k, p in manifests.items() if load(p).get("kind", "files") == "files"]

    ok = True
    for ident in chosen:
        manifest = load(manifests[ident])
        print(f"{ident}: {manifest['title']}")
        kind = manifest.get("kind", "files")
        if kind == "files":
            ok &= acquire_files(manifests[ident], manifest, args.record)
        else:
            ok &= ACQUIRE[kind](manifest, args.profile, overrides, args.refresh)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
