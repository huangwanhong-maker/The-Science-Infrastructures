"""Keep the forecasts: each run of a forecast product, as it was, in an archive.

A field product of a forecast (field_grid.py, ai_forecast.py) is rewritten each
time its source is fetched again, and the forecast it held is lost. A forecast
can only be scored once the weather it foretold has happened, days later; so
each run is kept, as the product was when it was built, in

    data/cache/archive/<product>/<run>/fields-<product>.json and .bin

with the archive's own record added to the JSON (``archive``), and an index of
the runs, data/cache/archive/index.json, which verification.py reads.

A run is named by its start, the instant of the analysis it was computed from,
written compactly (2026-10-22T06Z). The products do not all state it, and it is
found by the first of these that applies:

1. the product states it (``start``, the start of its own ``verification``, or
   ``source.run`` given as an instant, as ecmwf_fields.py writes it);
2. the product was computed here from an analysis (``source.run`` a record of
   the computation, as for Pangu-Weather): its first time is the analysis;
3. the retrieval it was built from is still in the cache, unchanged, and
   records the run each time comes from (THREDDS's ``reftime``): the latest;
4. otherwise, the latest time of the product at or before the retrieval time
   minus 4 hours, the time a run takes to be computed and published (a GFS run
   is complete on NOMADS about 3.5 to 4 hours after its analysis time). A
   product that begins after that instant (a series asked for from "now") has
   its time step extended backwards to it.

The rule of step 4 is an estimate: a run published late is named a cycle too
early, and the record says which step named each run (``basis``). Where the
source records the run of each time (step 3), the times after the start that
come from an earlier run (a best time series fetched while its latest run was
still being published) are listed in the record (``other_runs``), and are not
scored as that run's forecast.

The forecast products are, unless chosen: every field product of a model with
more than four times that is not a climatology. Archiving is idempotent: a run
already kept from the same retrieval is left as it is, and a later retrieval of
the same run replaces it (it is at least as complete). Runs older than the
retention (30 days unless set) are removed.

A run of every product every six hours takes about 130 MB (most of it the GFS
aloft), 16 GB in 30 days; ``--keep scored`` keeps only the variables that
verification.py scores (500 hPa height, 2 m temperature, sea-level pressure,
wave height), about a quarter of that.

Usage::

    python data/products/archive.py                                    # every forecast product
    python data/products/archive.py --products gfs-levels-forecast aigfs-forecast
    python data/products/archive.py --retention-days 14 --keep scored
    python data/products/archive.py --list                             # the runs kept
"""

from __future__ import annotations

import argparse
import datetime as dt
import functools
import json
import math
import shutil
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
from product_index import OUT  # noqa: E402

ARCHIVE = DATA / "cache" / "archive"
LATENCY = dt.timedelta(hours=4)       # from a run's analysis time to its publication
RETENTION_DAYS = 30
UTC = dt.timezone.utc


def parse_time(text: str) -> dt.datetime:
    """An ISO 8601 instant as an aware datetime in UTC (no offset means UTC)."""
    t = dt.datetime.fromisoformat(text.replace("Z", "+00:00"))
    return (t if t.tzinfo else t.replace(tzinfo=UTC)).astimezone(UTC)


def iso(t: dt.datetime) -> str:
    return t.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_name(t: dt.datetime) -> str:
    """A run's folder name: its start in UTC, to the hour (2026-10-22T06Z), or the minute if not on the hour."""
    t = t.astimezone(UTC)
    return t.strftime("%Y-%m-%dT%HZ") if t.minute == 0 else t.strftime("%Y-%m-%dT%H%MZ")


def parse_run_name(name: str) -> dt.datetime:
    clock = name[11:-1]
    return dt.datetime.strptime(name[:10], "%Y-%m-%d").replace(
        hour=int(clock[:2]), minute=int(clock[2:] or 0), tzinfo=UTC)


def forecast_products(index: dict) -> list[str]:
    """Field products of a model, with more than four times, that are not climatologies."""
    return [e["name"] for e in index.get("fields", [])
            if e.get("model") and not e.get("climatology") and (e.get("times") or {}).get("count", 0) > 4]


def start_by_rule(times: list[dt.datetime], retrieved: dt.datetime) -> dt.datetime:
    """Step 4: the latest time at or before the retrieval minus LATENCY, the time step extended backwards if need be."""
    limit = retrieved - LATENCY
    earlier = [t for t in times if t <= limit]
    if earlier:
        return max(earlier)
    step = times[1] - times[0] if len(times) > 1 else dt.timedelta(hours=6)
    return times[0] - step * math.ceil((times[0] - limit) / step)


def source_runs(meta: dict, data: Path = DATA) -> dict[str, str]:
    """For each time of the product, the run it comes from (THREDDS's reftime, ISO 8601), read from the
    files of the product's retrieval if that retrieval is still the product's own (the same retrieval
    time); empty if it is not, or records no runs."""
    source = meta.get("source") or {}
    if not source.get("dataset") or not source.get("profile"):
        return {}
    return dict(_reftimes(str(data / "cache" / source["dataset"] / source["profile"]), source.get("retrieved_at"),
                          tuple(meta["times"])))


@functools.lru_cache(maxsize=32)
def _reftimes(folder: str, retrieved_at: str | None, times: tuple[str, ...]) -> tuple[tuple[str, str], ...]:
    import field_grid
    import netcdf3
    folder = Path(folder)
    try:
        retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ()
    if retrieval.get("retrieved_at") != retrieved_at:
        return ()
    wanted = set(times)
    for request in sorted(retrieval.get("requests", []), key=lambda r: r.get("size") or 0):
        path = folder / request.get("file", "")
        if path.suffix != ".nc" or not path.exists():
            continue
        try:
            ds = netcdf3.Dataset.open(path)
            if "reftime" not in ds.variables or "time" not in ds.variables:
                continue
            pairs = tuple((t, r) for t, r in zip(field_grid.decode_times(ds, "time"), field_grid.decode_times(ds, "reftime"))
                          if t in wanted)
        except (OSError, ValueError, KeyError, SystemExit):
            continue
        if pairs:
            return pairs
    return ()


def start_from_retrieval(meta: dict, data: Path = DATA) -> dt.datetime | None:
    """Step 3: the latest run named by THREDDS's reftime in the files of the product's retrieval."""
    runs = source_runs(meta, data)
    return parse_time(max(runs.values())) if runs else None


def other_runs(meta: dict, start: dt.datetime, data: Path = DATA) -> list[str]:
    """Times after the start that come from another run: a best time series fetched while its latest run
    was still being published takes its later times from the run before. They are not that run's forecast."""
    return sorted(t for t, r in source_runs(meta, data).items() if parse_time(t) > start and parse_time(r) != start)


def start_of(meta: dict, data: Path = DATA) -> tuple[dt.datetime, str]:
    """A forecast product's start, and which rule found it (see the module's description)."""
    stated = meta.get("start") or (meta.get("verification") or {}).get("start")
    if stated:
        return parse_time(stated), "stated by the product"
    times = [parse_time(t) for t in meta["times"]]
    source = meta.get("source") or {}
    if isinstance(source.get("run"), str):
        return parse_time(source["run"]), "stated by its source (source.run)"
    if isinstance(source.get("run"), dict):
        return times[0], "computed here from the analysis at its first time"
    found = start_from_retrieval(meta, data)
    if found:
        return found, "the run of its source (reftime)"
    return start_by_rule(times, parse_time(source["retrieved_at"])), "the latest time at or before the retrieval minus 4 h"


def subset(meta: dict, codes: bytes, names: set[str]) -> tuple[dict, bytes]:
    """The product with only the named variables, and its binary renumbered to match."""
    n = meta["grid"]["nlon"] * meta["grid"]["nlat"] * 2
    out, blocks, offset = dict(meta), [], 0
    variables = []
    for v in meta["variables"]:
        if v["name"] not in names:
            continue
        offsets = []
        for o in v["offsets"]:
            blocks.append(codes[o:o + n])
            offsets.append(offset)
            offset += n
        variables.append({**v, "offsets": offsets})
    out["variables"] = variables
    out["vectors"] = [w for w in meta.get("vectors", []) if w["u"] in names and w["v"] in names]
    return out, b"".join(blocks)


def scored_variables(meta: dict) -> set[str]:
    import verification
    return {v["name"] for q in verification.QUANTITIES if (v := verification.find_variable(meta, q))}


def snapshot(name: str, products: Path = OUT, archive: Path = ARCHIVE, keep: str = "all",
             now: dt.datetime | None = None, data: Path = DATA) -> dict | None:
    """Keeps the current run of product `name`; returns its record, or None if it was kept already."""
    path = products / f"fields-{name}.json"
    meta = json.loads(path.read_text(encoding="utf-8"))
    retrieved = (meta.get("source") or {}).get("retrieved_at")
    if not retrieved or not meta.get("times"):
        print(f"  {name}: no retrieval time or no times; not archived")
        return None
    start, basis = start_of(meta, data)
    run = run_name(start)
    folder = archive / name / run
    kept = folder / path.name
    if kept.exists():
        old = json.loads(kept.read_text(encoding="utf-8"))
        if (old.get("source") or {}).get("retrieved_at") == retrieved:
            return None
    codes = (products / meta["file"]).read_bytes()
    if keep == "scored":
        meta, codes = subset(meta, codes, scored_variables(meta))
        if not meta["variables"]:
            print(f"  {name}: no variable that is scored; not archived")
            return None
    record = {"run": run, "start": iso(start), "basis": basis, "keep": keep,
              "archived_at": iso(now or dt.datetime.now(UTC))}
    foreign = other_runs(meta, start, data)
    if foreign:
        record["other_runs"] = foreign
    meta = {**meta, "archive": record}
    # Written beside the run, then moved into place, so a run is never half there.
    part = archive / name / f"{run}.part"
    if part.exists():
        shutil.rmtree(part)
    part.mkdir(parents=True)
    (part / meta["file"]).write_bytes(codes)
    (part / path.name).write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    if folder.exists():
        shutil.rmtree(folder)
    part.rename(folder)
    return record


def prune(archive: Path = ARCHIVE, now: dt.datetime | None = None, days: float = RETENTION_DAYS) -> list[Path]:
    """Removes the runs that started more than `days` before `now`; returns their folders."""
    limit = (now or dt.datetime.now(UTC)) - dt.timedelta(days=days)
    removed = []
    for folder in sorted(archive.glob("*/*Z")):
        try:
            start = parse_run_name(folder.name)
        except ValueError:
            continue
        if folder.is_dir() and start < limit:
            shutil.rmtree(folder)
            removed.append(folder)
    for product in archive.iterdir() if archive.exists() else []:
        if product.is_dir() and not any(product.iterdir()):
            product.rmdir()
    return removed


def rebuild_index(archive: Path = ARCHIVE, days: float = RETENTION_DAYS) -> dict:
    """The index of the archive, read from the runs that are there."""
    products = {}
    for path in sorted(archive.glob("*/*Z/fields-*.json")):
        meta = json.loads(path.read_text(encoding="utf-8"))
        rec = meta.get("archive") or {}
        entry = products.setdefault(path.parents[1].name, {"title": meta.get("title"), "model": meta.get("model"),
                                                           "runs": []})
        entry["title"], entry["model"] = meta.get("title"), meta.get("model")
        times = meta["times"]
        entry["runs"].append({
            "run": path.parent.name, "start": rec.get("start"), "basis": rec.get("basis"),
            "times": {"count": len(times), "first": times[0], "last": times[-1]},
            "retrieved_at": (meta.get("source") or {}).get("retrieved_at"),
            "archived_at": rec.get("archived_at"), "keep": rec.get("keep", "all"),
            "other_runs": rec.get("other_runs", []),
            "variables": [v["name"] for v in meta["variables"]],
            "path": path.relative_to(archive).as_posix(),
        })
    for entry in products.values():
        entry["runs"].sort(key=lambda r: r["start"] or "")
    index = {"about": "Runs of forecast products, kept as they were built, so that they can be scored "
                      "later (data/products/archive.py, verification.py). Paths are relative to this folder; "
                      "times are UTC.",
             "updated": iso(dt.datetime.now(UTC)), "retention_days": days, "products": products}
    archive.mkdir(parents=True, exist_ok=True)
    (archive / "index.json").write_text(json.dumps(index, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return index


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--products", nargs="*", help="names of field products (default: every forecast product)")
    parser.add_argument("--retention-days", type=float, default=RETENTION_DAYS, help="keep runs this many days")
    parser.add_argument("--keep", choices=["all", "scored"], default="all",
                        help="every variable of a run, or only those verification.py scores")
    parser.add_argument("--list", action="store_true", help="list the runs kept, and archive nothing")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace", line_buffering=True)

    if args.list:
        index = json.loads((ARCHIVE / "index.json").read_text(encoding="utf-8")) if (ARCHIVE / "index.json").exists() else {}
        for name, entry in (index.get("products") or {}).items():
            print(f"{name} ({entry['model']}): {len(entry['runs'])} run(s)")
            for r in entry["runs"]:
                print(f"  {r['run']:18s} {r['times']['count']:3d} times to {r['times']['last']}  "
                      f"retrieved {r['retrieved_at']}  ({r['basis']})")
        return 0

    index_path = OUT / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8")) if index_path.exists() else {}
    names = args.products if args.products else forecast_products(index)
    failed = []
    for name in names:
        if not (OUT / f"fields-{name}.json").exists():
            print(f"  {name}: no such product")
            failed.append(name)
            continue
        try:
            record = snapshot(name, keep=args.keep)
        except (OSError, ValueError, KeyError) as error:
            print(f"  {name}: FAILED ({type(error).__name__}: {error})")
            failed.append(name)
            continue
        if record:
            print(f"  {name}: kept the run of {record['start']} ({record['basis']})")
        else:
            print(f"  {name}: its run is kept already")
    removed = prune(days=args.retention_days)
    for folder in removed:
        print(f"  removed {folder.parent.name}/{folder.name}, older than {args.retention_days:g} days")
    index = rebuild_index(days=args.retention_days)
    runs = sum(len(e["runs"]) for e in index["products"].values())
    print(f"archive: {runs} run(s) of {len(index['products'])} product(s) in {ARCHIVE.relative_to(DATA.parent)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
