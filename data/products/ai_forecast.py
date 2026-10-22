"""Run an AI weather model from an analysis, and build its forecast as a field product.

The laboratory runs the model itself, on its own computer, from an analysis
the data library has fetched; the forecast is then a field product like any
other, which the Earth Atlas draws, compares with other models and verifies.

Pangu-Weather (Bi et al. 2023) is a 3D Earth-specific transformer trained on
ERA5. It maps the state of the atmosphere at one instant (five variables at
13 pressure levels and four at the surface, on a 0.25 degree grid) to the
state 1, 3, 6 or 24 hours later; the laboratory uses the 24 and 6 hour models.
A forecast of n hours is composed from them by hierarchical temporal
aggregation (Bi et al. 2023, Methods): the 24 hour model is applied as often
as it fits, from the analysis, and the 6 hour model fills the remainder from
the last 24 hour state, which keeps the number of steps, and so the
accumulated error, small (a lead of 66 hours is 24 + 24 + 6 + 6 + 6).

The analysis comes from another system than ERA5 (here the GFS); the model
takes it as it is, with the loss of accuracy its authors report for such
initial states. Geopotential is height times g = 9.80665 m s-2.

The model's parameters are licensed CC BY-NC-SA 4.0: the forecasts are for
research and learning, not for commercial use, and say so in their source.

Output: fields-<name>.json and .bin, every 6 hours, at 1 degree (every fourth
point of the model's grid, as the other 1 degree products are point values):
pressure at sea level, temperature at 2 m, wind at 10 m, geopotential height
at 500 hPa, temperature at 850 hPa, and wind at 850, 500 and 250 hPa. With
--verify, the root-mean-square error of 500 hPa height over 20-80 N against
the GFS's analyses (from the product fields-gfs-levels-forecast) is recorded,
with persistence for comparison.

Usage (in the laboratory's model environment, which has onnxruntime and ecCodes)::

    .venv/Scripts/python data/products/ai_forecast.py --source gfs-analysis-0p25 --hours 168
    .venv/Scripts/python data/products/ai_forecast.py --source gfs-analysis-0p25 --profile custom-272f93c4 --hours 72 --name pangu-hindcast --verify
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
import netcdf3  # noqa: E402
import field_grid  # noqa: E402
from product_index import OUT  # noqa: E402

GRAVITY = 9.80665
LEVELS = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50]     # hPa, the model's order
UPPER = ["gh", "q", "t", "u", "v"]           # ecCodes short names, in the model's order (Z, Q, T, U, V)
SURFACE = [("prmsl", "meanSea", 0), ("10u", "heightAboveGround", 10), ("10v", "heightAboveGround", 10),
           ("2t", "heightAboveGround", 2)]
STRIDE = 4                                    # 0.25 -> 1 degree
MODEL = "Pangu-Weather (laboratory)"


def read_analysis(path: Path) -> tuple[np.ndarray, np.ndarray, dt.datetime]:
    """The model's inputs from a GRIB2 analysis: (5, 13, 721, 1440) and (4, 721, 1440), float32, and its time."""
    import eccodes
    upper = np.full((5, 13, 721, 1440), np.nan, np.float32)
    surface = np.full((4, 721, 1440), np.nan, np.float32)
    when = None
    with open(path, "rb") as f:
        while (h := eccodes.codes_grib_new_from_file(f)) is not None:
            try:
                name, kind, level = (eccodes.codes_get(h, k) for k in ("shortName", "typeOfLevel", "level"))
                if eccodes.codes_get(h, "Ni") != 1440 or eccodes.codes_get(h, "Nj") != 721 \
                        or eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees") != 90.0:
                    raise SystemExit(f"{path}: expected the 0.25 degree grid from 90 N and 0 E")
                values = eccodes.codes_get_values(h).reshape(721, 1440).astype(np.float32)
                date, hhmm = eccodes.codes_get(h, "dataDate"), eccodes.codes_get(h, "dataTime")
                when = dt.datetime.strptime(f"{date}{hhmm:04d}", "%Y%m%d%H%M").replace(tzinfo=dt.timezone.utc)
                if kind == "isobaricInhPa" and name in UPPER and level in LEVELS:
                    if name == "gh":
                        values = values * GRAVITY                 # geopotential, m2 s-2
                    upper[UPPER.index(name), LEVELS.index(level)] = values
                for i, (sname, skind, slevel) in enumerate(SURFACE):
                    if name == sname and kind == skind and level == slevel:
                        surface[i] = values
            finally:
                eccodes.codes_release(h)
    missing = int(np.isnan(upper).any(axis=(2, 3)).sum() + np.isnan(surface).any(axis=(1, 2)).sum())
    if missing:
        raise SystemExit(f"{path}: {missing} of the model's 69 input fields are missing")
    return upper, surface, when


class Pangu:
    """The 24 and 6 hour models, and their composition into a forecast."""

    def __init__(self, folder: Path, threads: int):
        import onnxruntime as ort
        options = ort.SessionOptions()
        options.intra_op_num_threads = threads
        options.enable_cpu_mem_arena = False
        options.enable_mem_pattern = False
        options.enable_mem_reuse = False
        options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        self.models = {}
        for hours in (24, 6):
            path = folder / f"pangu_weather_{hours}.onnx"
            if not path.exists():
                raise SystemExit(f"{path} is missing; run: python data/acquisition/fetch.py pangu-weather")
            self.models[hours] = ort.InferenceSession(str(path), sess_options=options, providers=["CPUExecutionProvider"])

    def step(self, hours: int, upper: np.ndarray, surface: np.ndarray):
        out_upper, out_surface = self.models[hours].run(None, {"input": upper, "input_surface": surface})
        return out_upper, out_surface

    def forecast(self, upper, surface, total_hours: int, every: int = 6, on_step=None):
        """Yields (lead_hours, upper, surface) every `every` hours up to total_hours."""
        state24 = (upper, surface)
        state = (upper, surface)
        for lead in range(every, total_hours + 1, every):
            started = time.monotonic()
            if lead % 24 == 0:
                state24 = self.step(24, *state24)
                state = state24
            else:
                state = self.step(6, *state)
            if on_step:
                on_step(lead, time.monotonic() - started)
            yield lead, state[0], state[1]


# The fields kept for the atlas, by (file group, level): {name: (long name, standard name, units, values)}.
# Names and long names leave the level out; field_grid.build adds it.
def kept(upper: np.ndarray, surface: np.ndarray) -> dict:
    s = (slice(None, None, STRIDE), slice(None, None, STRIDE))
    lev = lambda p: LEVELS.index(p)
    return {
        ("msl", None): {"prmsl": ("Pressure at mean sea level", "air_pressure_at_sea_level", "Pa", surface[0][s])},
        ("2m", "2 m"): {"t": ("Air temperature", "air_temperature", "K", surface[3][s])},
        ("10m", "10 m"): {"u": ("Eastward wind", "eastward_wind", "m s-1", surface[1][s]),
                          "v": ("Northward wind", "northward_wind", "m s-1", surface[2][s])},
        ("850hPa", "850 hPa"): {"t": ("Temperature", "air_temperature", "K", upper[2, lev(850)][s]),
                                "u": ("Eastward wind", "eastward_wind", "m s-1", upper[3, lev(850)][s]),
                                "v": ("Northward wind", "northward_wind", "m s-1", upper[4, lev(850)][s])},
        ("500hPa", "500 hPa"): {"z": ("Geopotential height", "geopotential_height", "m", upper[0, lev(500)][s] / GRAVITY),
                                "u": ("Eastward wind", "eastward_wind", "m s-1", upper[3, lev(500)][s]),
                                "v": ("Northward wind", "northward_wind", "m s-1", upper[4, lev(500)][s])},
        ("250hPa", "250 hPa"): {"u": ("Eastward wind", "eastward_wind", "m s-1", upper[3, lev(250)][s]),
                                "v": ("Northward wind", "northward_wind", "m s-1", upper[4, lev(250)][s])},
    }


def write_groups(folder: Path, init: dt.datetime, leads: list[int], frames: list[dict]) -> list[tuple[Path, str | None]]:
    """One netCDF classic file per level group, as field_grid.build reads them."""
    lat = np.arange(90.0, -90.0 - 1e-9, -0.25 * STRIDE)
    lon = np.arange(0.0, 360.0 - 1e-9, 0.25 * STRIDE)
    files = []
    for key in frames[0]:
        group, level = key
        variables = {
            "time": (["time"], np.array(leads, np.float64), {"units": f"hours since {init:%Y-%m-%d %H:%M:%S}", "standard_name": "time"}),
            "lat": (["lat"], lat, {"units": "degrees_north", "standard_name": "latitude"}),
            "lon": (["lon"], lon, {"units": "degrees_east", "standard_name": "longitude"}),
        }
        for name, (long_name, standard, units, _) in frames[0][key].items():
            data = np.stack([f[key][name][3] for f in frames]).astype(np.float32)
            variables[name] = (["time", "lat", "lon"], data, {"long_name": long_name, "standard_name": standard, "units": units})
        path = folder / f"{group}.nc"
        netcdf3.write(path, {"time": len(leads), "lat": len(lat), "lon": len(lon)}, variables,
                      {"title": f"{MODEL} forecast from {init:%Y-%m-%d %H:%M} UTC", "source": "data/products/ai_forecast.py"})
        files.append((path, level))
    return files


def frame_of(upper, surface) -> dict:
    return {k: {n: (a, b, c, np.array(d, np.float32)) for n, (a, b, c, d) in v.items()} for k, v in kept(upper, surface).items()}


def rmse_band(a: np.ndarray, b: np.ndarray, lat: np.ndarray, south=20.0, north=80.0) -> float:
    band = (lat >= south) & (lat <= north)
    w = np.cos(np.radians(lat[band]))[:, None] * np.ones((1, a.shape[1]))
    d = (a[band] - b[band]) ** 2
    ok = np.isfinite(d)
    return float(np.sqrt((d[ok] * w[ok]).sum() / w[ok].sum()))


def verify(init: dt.datetime, leads: list[int], z500: list[np.ndarray]) -> list[dict]:
    """RMSE of 500 hPa height over 20-80 N against the GFS (1 degree), and of persistence."""
    meta_path = OUT / "fields-gfs-levels-forecast.json"
    if not meta_path.exists():
        print("  no fields-gfs-levels-forecast product to verify against")
        return []
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    var = next(v for v in meta["variables"] if v["name"] == "Geopotential_height_isobaric@500hPa")
    g = meta["grid"]
    codes = np.fromfile(OUT / meta["file"], "<i2")
    n = g["nlon"] * g["nlat"]
    def slice_at(i):
        q = codes[var["offsets"][i] // 2: var["offsets"][i] // 2 + n].reshape(g["nlat"], g["nlon"]).astype(np.float64)
        out = q * var["scale"] + var["offset"]
        out[q == -32768] = np.nan
        return out
    times = [dt.datetime.fromisoformat(t.replace("Z", "+00:00")) for t in meta["times"]]
    retrieved = dt.datetime.fromisoformat(meta["source"]["retrieved_at"])
    lat = g["lat0"] - np.arange(g["nlat"]) * g["dlat"]
    # Our grid starts at 0 E; the product's at its lon0 (-180): roll ours to match.
    shift = int(round((g["lon0"] % 360) / g["dlon"]))
    if init not in times:
        print(f"  the GFS product has no analysis at {init:%Y-%m-%d %H:%M}")
        return []
    start = slice_at(times.index(init))
    out = []
    for lead, z in zip(leads, z500):
        t = init + dt.timedelta(hours=lead)
        if t not in times:
            continue
        truth = slice_at(times.index(t))
        ours = np.roll(z, -shift, axis=1)
        out.append({"hours": lead, "model": round(rmse_band(ours, truth, lat), 1),
                    "persistence": round(rmse_band(start, truth, lat), 1),
                    "against": "analysis" if t <= retrieved - dt.timedelta(hours=6) else "forecast"})
    return out


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", default="gfs-analysis-0p25", help="a registered analysis dataset, already fetched")
    parser.add_argument("--profile", default="latest", help="its profile (cache folder)")
    parser.add_argument("--hours", type=int, default=168, help="length of the forecast, a multiple of 6")
    parser.add_argument("--name", default="pangu-forecast", help="product name")
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help="CPU threads for the model")
    parser.add_argument("--verify", action="store_true", help="verify 500 hPa height against the GFS's analyses")
    parser.add_argument("--label", default=MODEL, help="the model's name in the atlas, distinct for each product")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace")

    folder = DATA / "cache" / args.source / args.profile
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    upper, surface, init = read_analysis(folder / retrieval["requests"][0]["file"])
    print(f"analysis {init:%Y-%m-%d %H:%M} UTC from {args.source}/{args.profile}")
    weights = json.loads((DATA / "registry" / "pangu-weather.json").read_text(encoding="utf-8"))
    model = Pangu(DATA / "cache" / "pangu-weather", args.threads)

    leads, frames = [0], [frame_of(upper, surface)]
    started = time.monotonic()
    for lead, up, sf in model.forecast(upper, surface, args.hours,
                                       on_step=lambda h, s: print(f"  +{h:3d} h: {s:5.1f} s", flush=True)):
        leads.append(lead)
        frames.append(frame_of(up, sf))
    elapsed = time.monotonic() - started
    print(f"{len(leads) - 1} steps in {elapsed / 60:.1f} min")

    work = DATA / "cache" / "ai-forecasts" / args.name
    work.mkdir(parents=True, exist_ok=True)
    files = write_groups(work, init, leads, frames)
    source = {
        "dataset": args.source, "profile": args.profile,
        "provider": f"computed by the laboratory with {weights['title'].split(':')[0]} from the {retrieval['dataset']} analysis of {init:%Y-%m-%d %H:%M} UTC",
        "credit": f"Pangu-Weather (Bi et al. 2023), run by the laboratory from the GFS analysis of {init:%d %b %H} UTC",
        "licence": "Derived from Pangu-Weather, whose parameters are licensed CC BY-NC-SA 4.0: for research and learning, not for commercial use",
        "citation": weights["citation"],
        "retrieved_at": retrieval["retrieved_at"],
        "request": retrieval["requests"][0]["url"],
        "run": {"steps": len(leads) - 1, "minutes": round(elapsed / 60, 1), "threads": args.threads,
                "composition": "24 hour model as often as it fits, then the 6 hour model"},
    }
    product = field_grid.build(files, args.name, f"{args.label}, from {init:%Y-%m-%d %H:%M} UTC", "meteorology", source, {}, args.label)
    if args.verify:
        z500 = [f[("500hPa", "500 hPa")]["z"][3] for f in frames[1:]]
        table = verify(init, leads[1:], z500)
        for r in table:
            print(f"  +{r['hours']:3d} h  model {r['model']:6.1f} m   persistence {r['persistence']:6.1f} m   ({r['against']})")
        path = OUT / f"fields-{args.name}.json"
        meta = json.loads(path.read_text(encoding="utf-8"))
        meta["verification"] = {"quantity": "root-mean-square error of 500 hPa geopotential height, 20-80 N, m",
                                "against": "GFS (fields-gfs-levels-forecast)", "start": f"{init:%Y-%m-%dT%H:%M:%SZ}", "leads": table}
        path.write_text(json.dumps(meta, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    g = product["grid"]
    print(f"fields-{args.name}: {len(product['variables'])} variables, {len(product['times'])} times, {g['nlon']} x {g['nlat']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
