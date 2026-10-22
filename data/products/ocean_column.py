"""Build a column product: the ocean beneath a place, by depth and over time.

A column is what a current profiler moored at a place would record: at each
time, the eastward and northward current and the temperature at each depth,
from the surface to the sea floor. Here it comes from an ocean model's
analyses and forecast at the model's grid point nearest the place (a column
dataset such as data/registry/hycom-column.json, one netCDF classic file per
variable, dimensions time x depth).

The product says which times are forecasts: those after the latest run of the
model (its forecast reference time). It also gives the deepest level with
values, above which the sea floor lies.

Output: data/cache/products/columns/column-<dataset>-<profile>.json (schema
columns/1). Columns are made on demand (data/acquisition/on_demand.py), one
per place, so they are not listed in the index of products.

Usage::

    python data/acquisition/fetch.py hycom-column --set lat=32 --set lon=135
    python data/products/ocean_column.py --source hycom-column --profile custom-<hash>
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
import netcdf3  # noqa: E402
from product_index import OUT  # noqa: E402

COLUMNS = OUT / "columns"
DIGITS = {"m s-1": 3, "m/s": 3, "degC": 3, "psu": 3}


def _iso(seconds: float) -> str:
    import datetime as dt
    return dt.datetime.fromtimestamp(seconds, dt.timezone.utc).isoformat().replace("+00:00", "Z")


def build(folder: Path, manifest: dict) -> Path:
    """The column product from a retrieval folder; returns its path."""
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    columns, times, depths, runs, place, first_runs = [], None, None, [], {}, None
    for r in retrieval["requests"]:
        ds = netcdf3.Dataset.open(folder / r["file"])
        t = [float(x) for x in ds.values("time")]
        z = [float(x) for x in ds.values("depth")]
        if times is None:
            times, depths = t, z
            place = {"grid_latitude": ds.attributes.get("latitude"), "grid_longitude": ds.attributes.get("longitude"),
                     "latitude": ds.attributes.get("requested_latitude"), "longitude": ds.attributes.get("requested_longitude")}
        elif t != times or z != depths:
            raise SystemExit(f"{r['file']}: its times or depths differ from those of {retrieval['requests'][0]['file']}")
        if "run" in ds.variables:
            these = [float(x) for x in ds.values("run")]
            runs += these
            first_runs = first_runs or these
        values = list(ds.values(r["variable"]))
        a = ds.variables[r["variable"]].attributes
        units = a.get("units", "")
        digits = DIGITS.get(units, 4)
        nz = len(z)
        rows = [[None if math.isnan(v) else round(v, digits) for v in values[i * nz:(i + 1) * nz]] for i in range(len(t))]
        columns.append({"name": r["variable"], "standard_name": a.get("standard_name"), "long_name": a.get("long_name"),
                        "units": units, "values": rows})
    # The deepest level with a value at any time, in any variable.
    deepest = max((k for c in columns for row in c["values"] for k, v in enumerate(row) if v is not None), default=-1)
    # The newest run used: times after it are forecasts. A run time later than the
    # retrieval is a mislabelled file of the server's collection, and is ignored.
    import datetime as dt
    retrieved = dt.datetime.fromisoformat(retrieval["retrieved_at"]).timestamp()
    latest_run = max((x for x in runs if x <= retrieved), default=None)
    lon = place.get("grid_longitude")
    if lon is not None and lon > 180:
        place["grid_longitude"] = lon - 360
    name = f"column-{manifest['id']}-{retrieval['profile']}"
    product = {
        "schema": "columns/1", "name": name,
        "title": f"{manifest.get('model', manifest['title'])}: the column beneath "
                 f"{abs(place['grid_latitude']):.2f}° {'N' if place['grid_latitude'] >= 0 else 'S'}, "
                 f"{abs(place['grid_longitude']):.2f}° {'E' if place['grid_longitude'] >= 0 else 'W'}",
        "domain": manifest.get("domain", "oceanography"), "model": manifest.get("model"),
        "place": place,
        "times": [_iso(x) for x in times],
        "depths": depths,
        "deepest_level": deepest,
        "latest_run": _iso(latest_run) if latest_run is not None else None,
        # The run each time comes from (of the first variable), as the server's collection gives it.
        "runs": [_iso(x) for x in first_runs] if first_runs else None,
        "variables": columns,
        "units": {"times": "UTC", "depths": "m, positive down"},
        "source": {"dataset": manifest["id"], "profile": retrieval["profile"], "provider": manifest["provider"],
                   "credit": manifest.get("credit"), "licence": manifest["licence"], "citation": manifest["citation"],
                   "retrieved_at": retrieval["retrieved_at"], "parameters": retrieval["parameters"],
                   "resolved": retrieval.get("resolved", {})},
    }
    COLUMNS.mkdir(parents=True, exist_ok=True)
    path = COLUMNS / f"{name}.json"
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    return path


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default="hycom-column")
    parser.add_argument("--profile", required=True, help="the retrieval's cache folder, such as custom-<hash>")
    args = parser.parse_args(argv)
    folder = DATA / "cache" / args.source / args.profile
    if not (folder / "retrieval.json").exists():
        print(f"not fetched: {folder}", file=sys.stderr)
        return 1
    manifest = json.loads((DATA / "registry" / f"{args.source}.json").read_text(encoding="utf-8"))
    path = build(folder, manifest)
    product = json.loads(path.read_text(encoding="utf-8"))
    print(f"{path.relative_to(DATA.parent)}: {len(product['times'])} times x {len(product['depths'])} depths, "
          f"{len(product['variables'])} variables; deepest level {product['depths'][product['deepest_level']]:.0f} m; "
          f"forecast after {product['latest_run']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
