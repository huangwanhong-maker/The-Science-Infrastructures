"""Build field products from ECMWF's open forecasts: the GRIB2 messages fetched by
data/acquisition/ecmwf_open.py, decoded with ecCodes, sampled to a coarser grid,
and written through field_grid's writer.

Each variable of a product is described in the dataset's manifest: the field it
comes from (param, levtype, levelist) and its name, level, long name, CF standard
name and units in the product. Where the GFS gives the same quantity, these are
the GFS products' (air_temperature at "2 m", geopotential_height at "500 hPa",
precipitation_flux in kg m-2 s-1), so that the Earth Atlas offers ECMWF beside
the GFS under each quantity, and pairs the eastward and northward winds of a
level into a vector, as field_grid.build does.

Conversions, each to the product's SI units:

- Geopotential height is given in geopotential metres and kept as it is. Where a
  stream gives geopotential Phi (m2 s-2) instead, the height is Z = Phi / g0, with
  g0 = 9.80665 m s-2, the standard acceleration of gravity (3rd CGPM, 1901), by
  which WMO defines the geopotential metre.
- Precipitation is given as the total accumulated from the start of the run: by
  the IFS in metres of water, by AIFS in kg m-2. A depth of water d weighs
  rho_w d per unit area, with rho_w = 1000 kg m-3 (ECMWF's convention for its
  water equivalents, parameter 228), so 1 m is 1000 kg m-2. The mean rate over
  the interval (t - dt, t] is P(t) = (A(t) - A(t - dt)) / dt, in kg m-2 s-1, and
  is given at t, the end of the interval, as the GFS gives its rate. At the start
  of the run there is no interval before, and no value.
- Probabilities are in percent, as ECMWF gives them. A quantity over a period
  (24 h precipitation, the largest gust) is placed at the end of the period.

*Caveat.* The accumulations are packed with 12 bits, in steps of a power of two
that grow with the largest total: 2^-14 m (0.06 mm) at 6 h, 2^-12 m (0.24 mm)
once some total passes 0.5 m, from about 60 h. Each total is rounded to within
half a step, so a difference of two is uncertain by up to one step: 0.24 mm in
6 h, or 0.04 mm an hour, far below the forecast's own error. A difference may
then be slightly negative where no rain fell; it is set to zero.

The grid is sampled, not averaged: every stride-th point of the 0.25 degree grid,
from 90 N and 180 W, as the other 1 degree products are point values. Rows run
north to south and columns eastward from -180 degrees, as field_grid
normalises every product.

Decoding needs ecCodes (the Python package eccodes, Apache-2.0, whose wheels
carry the ecCodes library). Without it the builder stops and says how to
install it; fetching needs only the standard library.

Output: fields-<dataset>-<profile>.json and .bin in data/cache/products/, through
field_grid.build, from netCDF classic files, one per level, written beside the
fetched messages (data/cache/<dataset>/<profile>/netcdf/).

Usage::

    .venv/Scripts/python data/products/ecmwf_fields.py --source ecmwf-ifs          # fields-ecmwf-ifs-forecast
    .venv/Scripts/python data/products/ecmwf_fields.py --source ecmwf-aifs         # fields-ecmwf-aifs-forecast
    .venv/Scripts/python data/products/ecmwf_fields.py --source ecmwf-ens          # fields-ecmwf-ens-probabilities
    .venv/Scripts/python data/products/ecmwf_fields.py --source ecmwf-ifs --profile glance
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
import netcdf3  # noqa: E402
import field_grid  # noqa: E402

STANDARD_GRAVITY = 9.80665       # m s-2, g0 (3rd CGPM, 1901)
WATER_DENSITY = 1000.0           # kg m-3, rho_w of ECMWF's water equivalents
# Factors from the units ecCodes reports to the units of the product.
FACTORS = {
    ("m", "kg m-2"): WATER_DENSITY, ("kg m**-2", "kg m-2"): 1.0,
    ("m**2 s**-2", "m"): 1.0 / STANDARD_GRAVITY, ("gpm", "m"): 1.0, ("m", "m"): 1.0,
    ("K", "K"): 1.0, ("Pa", "Pa"): 1.0, ("m s**-1", "m s-1"): 1.0, ("%", "%"): 1.0,
}
GRID_KEYS = ("gridType", "Ni", "Nj", "latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees",
             "iDirectionIncrementInDegrees", "jDirectionIncrementInDegrees", "iScansNegatively", "jScansPositively")
INSTALL = ("decoding GRIB2 needs ecCodes, which this Python cannot load ({error}). Install it in the data "
           "library's environment:\n    .venv/Scripts/python -m pip install eccodes==2.50.0\n"
           "(.venv/bin/python on Unix; the wheel carries the ecCodes library, so nothing else is needed), "
           "or install data/products/requirements-models.txt")


def load_eccodes():
    try:
        import eccodes
        eccodes.codes_get_api_version()
    except (ImportError, RuntimeError, OSError) as error:
        raise SystemExit(INSTALL.format(error=error)) from None
    return eccodes


def factor(grib_units: str, units: str) -> float:
    try:
        return FACTORS[(grib_units, units)]
    except KeyError:
        raise SystemExit(f"no conversion from {grib_units!r} to {units!r}; add it to FACTORS") from None


def mean_rates(accumulations: list[np.ndarray], leads: list[int]) -> list[np.ndarray]:
    """Mean rates (per second) over the intervals between leads (hours), from totals
    accumulated since the start of the run; none at the start itself. Negative
    differences, which only the packing's rounding makes, are set to zero."""
    out = []
    for i, (a, lead) in enumerate(zip(accumulations, leads)):
        before, start = (accumulations[i - 1], leads[i - 1]) if i else (np.zeros_like(a), 0)
        if lead == start:
            out.append(np.full_like(a, np.nan))
        else:
            out.append(np.maximum(a - before, 0.0) / ((lead - start) * 3600.0))
    return out


def sampled_axes(grid: dict, stride: int) -> tuple[np.ndarray, np.ndarray]:
    """Latitudes and longitudes (degrees) of every stride-th point of a regular grid
    that runs once round the Earth, rows from the north."""
    if grid["gridType"] != "regular_ll" or grid["iScansNegatively"] or grid["jScansPositively"]:
        raise SystemExit(f"expected a regular latitude-longitude grid scanned from the north-west, got {grid}")
    ni, nj, di, dj = grid["Ni"], grid["Nj"], grid["iDirectionIncrementInDegrees"], grid["jDirectionIncrementInDegrees"]
    if abs(ni * di - 360.0) > 1e-6 or ni % stride or (nj - 1) % stride:
        raise SystemExit(f"a stride of {stride} does not divide the {ni} x {nj} grid evenly")
    lat = grid["latitudeOfFirstGridPointInDegrees"] - dj * np.arange(0, nj, stride)
    lon = (grid["longitudeOfFirstGridPointInDegrees"] + di * np.arange(0, ni, stride) + 180.0) % 360.0 - 180.0
    return lat, lon


def messages(path: Path, expected: list[dict], eccodes):
    """Each message of a fetched file, with its entry in the provenance (in file
    order), its decoded keys, and its values (rows north to south)."""
    out = []
    with path.open("rb") as handle:
        while (h := eccodes.codes_grib_new_from_file(handle)) is not None:
            try:
                if len(out) == len(expected):
                    raise SystemExit(f"{path.name}: more messages than the provenance lists")
                entry = expected[len(out)]
                info = {k: eccodes.codes_get(h, k) for k in ("shortName", "units", "startStep", "endStep", "dataDate", "dataTime")}
                if info["shortName"] != entry["param"]:
                    raise SystemExit(f"{path.name}: message {len(out) + 1} is {info['shortName']}, the index said {entry['param']}")
                if entry.get("levelist") and int(entry["levelist"]) != eccodes.codes_get(h, "level"):
                    raise SystemExit(f"{path.name}: {entry['param']} is at {eccodes.codes_get(h, 'level')} hPa, "
                                     f"the index said {entry['levelist']}")
                grid = {k: eccodes.codes_get(h, k) for k in GRID_KEYS}
                values = eccodes.codes_get_values(h).astype(np.float64)
                if eccodes.codes_get(h, "bitmapPresent"):
                    values[values == eccodes.codes_get(h, "missingValue")] = np.nan
                out.append((entry, info, grid, values.reshape(grid["Nj"], grid["Ni"])))
            finally:
                eccodes.codes_release(h)
    if len(out) != len(expected):
        raise SystemExit(f"{path.name}: {len(out)} messages, the provenance lists {len(expected)}")
    return out


def matches(variable: dict, entry: dict, lead: int, kind: str) -> bool:
    params = variable["param"] if isinstance(variable["param"], list) else [variable["param"]]
    window = variable.get("window")
    step = f"{lead - int(window)}-{lead}" if window else str(lead)
    return (entry["param"] in params and entry.get("levtype") == variable.get("levtype", "sfc")
            and str(entry.get("levelist") or "") == str(variable.get("levelist", ""))
            and entry.get("type") == variable.get("type", kind) and str(entry.get("step")) == step)


def write_groups(folder: Path, cycle: dt.datetime, leads: list[int], lat, lon, variables: list[tuple[dict, np.ndarray]],
                 title: str) -> list[tuple[Path, str | None]]:
    """One netCDF classic file per level, as field_grid.build reads them."""
    folder.mkdir(parents=True, exist_ok=True)
    groups: dict[str | None, list] = {}
    for spec, data in variables:
        groups.setdefault(spec.get("level"), []).append((spec, data))
    files = []
    for level, members in groups.items():
        out = {
            "time": (["time"], np.array(leads, np.float64), {"units": f"hours since {cycle:%Y-%m-%d %H:%M:%S}", "standard_name": "time"}),
            "lat": (["lat"], lat, {"units": "degrees_north", "standard_name": "latitude"}),
            "lon": (["lon"], lon, {"units": "degrees_east", "standard_name": "longitude"}),
        }
        for spec, data in members:
            attrs = {"long_name": spec["long_name"], "units": spec["units"]}
            for key in ("standard_name", "cell_methods"):
                if spec.get(key):
                    attrs[key] = spec[key]
            out[spec["name"]] = (["time", "lat", "lon"], data.astype(np.float32), attrs)
        path = folder / f"{(level or 'single').replace(' ', '')}.nc"
        netcdf3.write(path, {"time": len(leads), "lat": len(lat), "lon": len(lon)}, out,
                      {"title": title, "source": "data/products/ecmwf_fields.py"})
        files.append((path, level))
    return files


def build(source_id: str, profile: str | None = None, name: str | None = None) -> dict:
    eccodes = load_eccodes()
    manifest = json.loads((DATA / "registry" / f"{source_id}.json").read_text(encoding="utf-8"))
    tag = profile or manifest.get("default_profile", "default")
    folder = DATA / "cache" / source_id / tag
    if not (folder / "retrieval.json").exists():
        raise SystemExit(f"not fetched; run: python data/acquisition/fetch.py {source_id} --profile {tag}")
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    kind = manifest["source"]["type"]
    stride = int(retrieval["parameters"].get("stride", 4))
    cycle = dt.datetime.fromisoformat(retrieval["resolved"]["cycle"].replace("Z", "+00:00"))
    leads = [int(s) for s in retrieval["resolved"]["steps"]]
    specs = manifest["variables"]

    # Every message, placed by variable and lead.
    found: dict[tuple[int, int], tuple[dict, np.ndarray]] = {}
    lat = lon = grid0 = None
    for r in retrieval["requests"]:
        for entry, info, grid, values in messages(folder / r["file"], r["messages"], eccodes):
            run = dt.datetime.strptime(f"{info['dataDate']}{info['dataTime']:04d}", "%Y%m%d%H%M").replace(tzinfo=dt.timezone.utc)
            if run != cycle:
                raise SystemExit(f"{r['file']}: a message of the run of {run:%Y-%m-%d %H} UTC, expected {cycle:%Y-%m-%d %H} UTC")
            if grid0 is None:
                grid0, (lat, lon) = grid, sampled_axes(grid, stride)
            elif grid != grid0:
                raise SystemExit(f"{r['file']}: {entry['param']} is on another grid")
            lead = int(info["endStep"])
            hits = [i for i, spec in enumerate(specs) if matches(spec, entry, lead, kind)]
            for i in hits:
                found[(i, lead)] = (info, values[::stride, ::stride])
    missing = [f"{specs[i]['name']} at {lead} h" for i in range(len(specs)) for lead in leads if (i, lead) not in found]
    if missing:
        raise SystemExit("missing from the fetched files: " + ", ".join(missing))

    interval = leads[1] - leads[0] if len(leads) > 1 else leads[0]
    variables, changes = [], [f"sampled to {0.25 * stride:g}° (every {stride}th point of the 0.25° grid)"]
    for i, spec in enumerate(specs):
        infos = [found[(i, lead)][0] for lead in leads]
        if spec.get("derive") == "rate":
            if any(info["startStep"] != 0 for info in infos):
                raise SystemExit(f"{spec['name']}: expected totals accumulated from the start of the run")
            totals = [found[(i, lead)][1] * factor(info["units"], "kg m-2") for lead, info in zip(leads, infos)]
            data = np.stack(mean_rates(totals, leads))
            changes.append(f"{spec['param']} accumulated from the start turned into mean rates over {interval} h, "
                           f"placed at the end of each interval")
        else:
            data = np.stack([found[(i, lead)][1] * factor(info["units"], spec["units"]) for lead, info in zip(leads, infos)])
            if any(info["units"] != spec["units"] and factor(info["units"], spec["units"]) != 1.0 for info in infos):
                changes.append(f"{infos[0]['shortName']} converted from {infos[0]['units']} to {spec['units']}")
        variables.append(({**spec, "long_name": spec["long_name"].format(interval=interval)}, data))

    label = manifest["title"].split(":")[0]
    title = f"{label}, run of {cycle:%Y-%m-%d %H:%M} UTC"
    files = write_groups(folder / "netcdf", cycle, leads, lat, lon, variables, title)
    source = {"dataset": manifest["id"], "profile": tag, "provider": manifest["provider"], "credit": manifest.get("credit"),
              "licence": manifest["licence"], "citation": manifest["citation"],
              "retrieved_at": retrieval["retrieved_at"], "request": retrieval["requests"][0]["url"],
              "run": f"{cycle:%Y-%m-%dT%H:%M:%SZ}", "changes": "; ".join(changes),
              "decoded_with": f"ecCodes {eccodes.codes_get_api_version()}"}
    return field_grid.build(files, name or f"{source_id}-{tag}", title, manifest.get("domain", "meteorology"),
                            source, {}, manifest.get("model"))


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True, help="a registered dataset of kind 'ecmwf', already fetched")
    parser.add_argument("--profile", help="the profile (cache folder) of the retrieval; default: the dataset's")
    parser.add_argument("--name", help="product name (default: <source>-<profile>)")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace")
    product = build(args.source, args.profile, args.name)
    g = product["grid"]
    print(f"fields-{product['name']}: {len(product['variables'])} variable(s) "
          f"({', '.join(v['name'] for v in product['variables'])}), {len(product['times'])} time(s) "
          f"from {product['times'][0]}, {g['nlon']} x {g['nlat']} grid, vectors: {[v['name'] for v in product['vectors']] or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
