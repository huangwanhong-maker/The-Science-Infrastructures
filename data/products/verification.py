"""Score the forecasts against what happened: the scorecard of the Earth Atlas.

A forecast is worth what its errors are, measured against the weather it
foretold, once that weather has happened. The atlas keeps each run of each
forecast (archive.py) and scores every run whose valid times have since been
observed or analysed, by lead time: the hours from the run's start to the
instant foretold. Each quantity is scored against its own truth.

1. **500 hPa height**, over the bands 20-80 deg N and 20-80 deg S, against
   analyses. The truth at a valid time T is the analysis at T in the earliest
   run (archived, or a current product) whose time T is at or before its start,
   preferring the GFS's analyses: those of the GFS's own products, and the GFS
   analysis that began a run computed here (Pangu-Weather). Where no GFS
   analysis of T is held, another model's is taken (ECMWF's, as its forecast's
   first time); each score lists every truth it used. Two scores, both with
   weights cos(latitude) (an area weight on a regular grid):

   - the root-mean-square error, RMSE = sqrt(sum w (f - a)^2 / sum w);
   - the anomaly correlation in its centred form (Wilks 2011, section 8.6.4,
     "Anomaly correlation", the centred ACC):

         ACC = sum w (f' - <f'>)(a' - <a'>) / sqrt(sum w (f' - <f'>)^2 . sum w (a' - <a'>)^2)

     with anomalies f' = f - c and a' = a - c from the climatology c, and <.>
     the weighted mean over the band. The climatology is NCEP/NCAR
     Reanalysis 1's long-term monthly means of 1991-2020 (Kalnay et al. 1996;
     registry ncep-r1-climatology), each month's mean placed at the middle of
     the month and interpolated linearly in time between them, then
     bilinearly from its 2.5 deg grid to the analyses' grid.

   A field is one case; a lead's RMSE is the root of the mean of its cases'
   mean squared errors, and its ACC the mean of their ACCs. Persistence (the
   analysis at the run's start, held) is scored beside the models.

2. **Station quantities**: 2 m temperature (and, where the observations carry
   them, sea-level pressure and significant wave height) against the reports
   of an observations product (METAR for the air, NDBC buoys for the waves).
   The forecast is sampled at each station as the page's FieldProduct.sample
   does: bilinear in space (physis grid.js, renormalised round missing
   values), linear in time between the forecast's times, and only within the
   forecast's span. RMSE and bias (mean of forecast minus observation) over
   the reports at or after the run's start, with each report's lead rounded to
   the nearest 6 hours (halves upward). No correction is made for the
   difference between a station's height and the model's ground: the products
   carry no model orography (were one added, 6.5 K per km would be the
   correction of temperature). Reports outside a plausible range (2 m
   temperature 180-340 K, sea-level pressure 850-1100 hPa, wave height 0-25 m)
   are left out, and sea-level pressure is scored only at stations at or below
   500 m, since a METAR without the SLP group gives the altimeter setting
   (QNH), which departs from the sea-level pressure at high stations.
   Temperatures of buoys (measured 3-5 m above the sea) are not taken for 2 m.

3. **Tropical cyclones**: in each run's sea-level pressure, the meteorology
   domain's tracker (domains/meteorology/js/cyclones.js, through
   track_cyclones.mjs, as the page runs it, except that a track may begin up
   to 40 degrees from the equator, since a run may begin in the middle of a
   storm's life) finds the cyclones; each is matched to an observed storm
   when, at the first time at or after the run's start that both have, their
   centres are within 3 degrees of arc, nearest pairs first, one track to one
   storm. The observed storms are IBTrACS's best tracks; IBTrACS lags by days,
   so each is extended past its last analysis by GDACS's observed positions of
   the storm of the same name (the warning centres' positions as they are
   issued), and GDACS's storms that IBTrACS does not hold yet are added. The
   errors at each common time thereafter: the great-circle distance between
   the centres (km, haversine on the sphere of radius R1 = 6371.0088 km, the
   mean radius of GRS 80; Sinnott 1984; Moritz 2000), and the central pressure,
   as the mean absolute error and the bias (forecast minus observed, hPa).
   Central pressures of a 1 degree field are well above the true ones; the
   bias says by how much.

The forecasts are the runs kept in the archive and the current products; the
same run of a model held by two products (the GFS through ERDDAP and through
THREDDS) is scored once for each quantity. Before any archive exists, what can
be scored is what the current products already contain: a run that began
before the latest analysis held (a hindcast, or a product a day or more old
beside a newer one), scored at the leads whose valid times are analysed. The
scores then rest on one or two runs; each states its count of cases, runs and
the window of valid times. Scores with no case are left out.

The ACC's centred form removes each field's mean anomaly before correlating;
the uncentred form (Wilks's second) does not, and is larger when the forecast
and analysis share a mean anomaly over the band.

References:
  Wilks, D. S. (2011), Statistical Methods in the Atmospheric Sciences, 3rd ed.,
  Academic Press, section 8.6 (8.6.3 mean squared error; 8.6.4 anomaly correlation).
  Kalnay, E., et al. (1996), The NCEP/NCAR 40-Year Reanalysis Project, Bull.
  Amer. Meteor. Soc. 77, 437-471.
  Sinnott, R. W. (1984), Virtues of the haversine, Sky and Telescope 68(2), 159.
  Moritz, H. (2000), Geodetic Reference System 1980, J. Geodesy 74, 128-133.

Output: data/cache/products/verification-earth-atlas.json (kind `verification`).

Usage::

    python data/products/verification.py
    python data/products/verification.py --no-cyclones          # without Node
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
ROOT = DATA.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
from archive import ARCHIVE, forecast_products, iso, other_runs, parse_time, start_of  # noqa: E402
from product_index import OUT, register  # noqa: E402

MISSING = -32768
NAME = "earth-atlas"
EARTH_RADIUS_KM = 6371.0088         # R1 of GRS 80 (Moritz 2000)
MATCH_DEG = 3.0                      # a model's cyclone and a storm, at their first common time
GAP_DAYS = 5                         # the longest gap between a storm's IBTrACS and GDACS positions
LEAD_STEP_H = 6
BANDS = {"20–80° N": (20.0, 80.0), "20–80° S": (-80.0, -20.0)}
QUANTITIES = ("z500", "t2m", "mslp", "hs")
TITLES = {"z500": "500 hPa height", "t2m": "2 m temperature", "mslp": "Sea-level pressure",
          "hs": "Significant wave height", "cyclones": "Tropical cyclones"}
# A station quantity: the observed variable's standard name, its plausible range (SI), the pattern its level
# must match (screen height for 2 m temperature, not a buoy's 3-5 m) and the highest station (m) where it is used.
STATIONS = {"t2m": {"standard_name": "air_temperature", "level": r"\b2 m\b|screen", "range": (180.0, 340.0)},
            "mslp": {"standard_name": "air_pressure_at_sea_level", "range": (85000.0, 110000.0), "highest": 500.0},
            "hs": {"standard_name": "sea_surface_wave_significant_height", "range": (0.0, 25.0)}}
CLIMATOLOGY = DATA / "cache" / "ncep-r1-climatology" / "ltm"
RELIEF = OUT / "fields-etopo-2022-medium.json"


# ---- reading products ----------------------------------------------------------------

def find_variable(meta: dict, quantity: str) -> dict | None:
    """The variable of a field product that holds a quantity (not an ensemble's spread)."""
    for v in meta.get("variables", []):
        sn, level = v.get("standard_name"), v.get("level")
        if "standard_deviation" in (v.get("cell_methods") or ""):
            continue
        if quantity == "z500" and sn == "geopotential_height" and level == "500 hPa":
            return v
        if quantity == "t2m" and sn == "air_temperature" and (level == "2 m" or (level is None and "2m" in v["name"].lower())):
            return v
        if quantity == "mslp" and sn == "air_pressure_at_sea_level":
            return v
        if quantity == "hs" and sn == "sea_surface_wave_significant_height":
            return v
    return None


def model_name(meta: dict) -> str:
    """The model of a product; a hindcast is the same model run from an earlier analysis."""
    return re.sub(r"\s+hindcast\b", "", meta.get("model") or meta["name"])


class Run:
    """One run of a forecast: a field product (archived or current) and its start."""

    def __init__(self, path: Path, origin: str):
        self.path, self.origin = path, origin
        self.meta = json.loads(path.read_text(encoding="utf-8"))
        self.product = self.meta["name"]
        self.model = model_name(self.meta)
        rec = self.meta.get("archive")
        self.start = parse_time(rec["start"]) if rec else start_of(self.meta)[0]
        self.grid = self.meta["grid"]
        # Times after the start that come from an earlier run end the run's forecast (archive.other_runs).
        foreign = set(rec.get("other_runs") or []) if rec else set(other_runs(self.meta, self.start))
        times = self.meta["times"]
        cut = next((i for i, t in enumerate(times) if t in foreign), len(times))
        self.times = [parse_time(t) for t in times[:cut]]
        self.ms = np.array([t.timestamp() * 1000 for t in self.times])
        self.bin = path.parent / self.meta["file"]
        self.cache: dict = {}
        # Whose analyses its times at or before its start are: a run computed here (source.run) began from
        # the analysis of its source dataset (the GFS's, for Pangu-Weather); a model's product, from its own.
        source = self.meta.get("source") or {}
        self.analysis, self.analysis_note = self.model, f"{self.product}, at times at or before its run's start"
        if isinstance(source.get("run"), dict) and source.get("dataset"):
            manifest = DATA / "registry" / f"{source['dataset']}.json"
            if manifest.exists():
                self.analysis = json.loads(manifest.read_text(encoding="utf-8")).get("model") or self.model
                self.analysis_note = f"{source['dataset']}, the analysis that began {self.product}"

    def field(self, variable: dict, ti: int) -> np.ndarray:
        """Values of a variable at a time index, (nlat, nlon), NaN for no value."""
        key = (variable["name"], ti)
        if key not in self.cache:
            n = self.grid["nlon"] * self.grid["nlat"]
            q = np.fromfile(self.bin, "<i2", count=n, offset=variable["offsets"][ti]).reshape(self.grid["nlat"], self.grid["nlon"])
            out = q.astype(np.float64) * variable["scale"] + variable["offset"]
            out[q == MISSING] = np.nan
            if len(self.cache) > 64:
                self.cache.pop(next(iter(self.cache)))
            self.cache[key] = out
        return self.cache[key]

    def label(self) -> str:
        return f"{self.product} {self.start:%Y-%m-%dT%HZ} ({self.origin})"


def runs_available(archive: Path = ARCHIVE, products: Path = OUT) -> list[Run]:
    """The archived runs, then the current forecast products that are not archived already."""
    runs, seen = [], set()
    index_path = archive / "index.json"
    if index_path.exists():
        index = json.loads(index_path.read_text(encoding="utf-8"))
        for entry in index.get("products", {}).values():
            for r in entry["runs"]:
                path = archive / r["path"]
                if path.exists():
                    runs.append(Run(path, "archive"))
                    seen.add((runs[-1].product, r["start"], r["retrieved_at"]))
    pindex = json.loads((products / "index.json").read_text(encoding="utf-8")) if (products / "index.json").exists() else {}
    for name in forecast_products(pindex):
        path = products / f"fields-{name}.json"
        if not path.exists():
            continue
        run = Run(path, "current")
        if (run.product, iso(run.start), run.meta["source"].get("retrieved_at")) not in seen:
            runs.append(run)
    runs.sort(key=lambda r: (r.start, r.origin != "archive", r.product))
    return runs


# ---- interpolation, as physis's grid.js does it ------------------------------------------

def bilinear_points(values: np.ndarray, grid: dict, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Values at places (degrees): bilinear between the four surrounding grid points, the weights of
    those with a value renormalised if they carry at least half the weight; NaN outside the grid.
    The rule of foundations/physis/js/grid.js (gridPosition, bilinear), vectorised."""
    nlon, nlat = grid["nlon"], grid["nlat"]
    lon, lat = np.asarray(lon, np.float64), np.asarray(lat, np.float64)
    x = (lon - grid["lon0"]) / grid["dlon"]
    y = (grid["lat0"] - lat) / grid["dlat"]
    inside = (y >= -0.5) & (y <= nlat - 0.5)
    if grid.get("periodic"):
        x = np.mod(x, nlon)
    else:
        inside &= (x >= -0.5) & (x <= nlon - 0.5)
    x, y = np.where(inside, x, 0.0), np.clip(np.where(inside, y, 0.0), 0, nlat - 1)
    c0, r0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - c0, y - r0
    c1 = c0 + 1
    c1 = np.where(c1 >= nlon, c1 - nlon if grid.get("periodic") else nlon - 1, c1)
    c0 = np.maximum(c0, 0)
    r1 = np.minimum(r0 + 1, nlat - 1)
    total, weight = np.zeros(lon.shape), np.zeros(lon.shape)
    for r, c, w in ((r0, c0, (1 - fx) * (1 - fy)), (r0, c1, fx * (1 - fy)), (r1, c0, (1 - fx) * fy), (r1, c1, fx * fy)):
        v = values[r, c]
        ok = np.isfinite(v)
        total += np.where(ok, w * np.nan_to_num(v), 0.0)
        weight += np.where(ok, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(inside & (weight >= 0.5), total / weight, np.nan)


def lattice(grid: dict) -> tuple[np.ndarray, np.ndarray]:
    """Longitudes and latitudes of a grid's points, (nlat, nlon) each."""
    lon = grid["lon0"] + np.arange(grid["nlon"]) * grid["dlon"]
    lat = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    return np.meshgrid(lon, lat)


def same_grid(a: dict, b: dict) -> bool:
    return all(abs(a[k] - b[k]) < 1e-6 for k in ("lon0", "dlon", "lat0", "dlat")) and \
        a["nlon"] == b["nlon"] and a["nlat"] == b["nlat"]


def regrid(values: np.ndarray, source: dict, target: dict) -> np.ndarray:
    if same_grid(source, target):
        return values
    lon, lat = lattice(target)
    return bilinear_points(values, source, lon, lat)


def sample_run(run: Run, variable: dict, lon: np.ndarray, lat: np.ndarray, ms: np.ndarray) -> np.ndarray:
    """A forecast at places and instants (ms since 1970, UTC): bilinear in space, linear in time between
    the run's times; NaN outside its span (FieldProduct.sample, without its extension past the ends)."""
    out = np.full(len(ms), np.nan)
    t = run.ms
    inside = (ms >= t[0]) & (ms <= t[-1])
    if len(t) < 2 or not inside.any():
        return out
    i0 = np.clip(np.searchsorted(t, ms, side="right") - 1, 0, len(t) - 2)
    for k in np.unique(i0[inside]):
        sel = inside & (i0 == k)
        w = (ms[sel] - t[k]) / (t[k + 1] - t[k])
        a = bilinear_points(run.field(variable, k), run.grid, lon[sel], lat[sel])
        b = bilinear_points(run.field(variable, k + 1), run.grid, lon[sel], lat[sel])
        out[sel] = np.where(w == 0, a, a + (b - a) * w)
    return out


# ---- scores ---------------------------------------------------------------------------

def band_weights(grid: dict, south: float, north: float) -> np.ndarray:
    """cos(latitude) on the rows of a band (inclusive), zero elsewhere, (nlat, nlon)."""
    lat = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    w = np.where((lat >= south - 1e-9) & (lat <= north + 1e-9), np.cos(np.radians(lat)), 0.0)
    return np.repeat(w[:, None], grid["nlon"], axis=1)


def rmse(f: np.ndarray, a: np.ndarray, w: np.ndarray) -> float:
    ok = np.isfinite(f) & np.isfinite(a) & (w > 0)
    return float(np.sqrt((w[ok] * (f[ok] - a[ok]) ** 2).sum() / w[ok].sum()))


def bias(f: np.ndarray, a: np.ndarray, w: np.ndarray | None = None) -> float:
    w = np.ones_like(f) if w is None else w
    ok = np.isfinite(f) & np.isfinite(a) & (w > 0)
    return float((w[ok] * (f[ok] - a[ok])).sum() / w[ok].sum())


def acc(f: np.ndarray, a: np.ndarray, c: np.ndarray, w: np.ndarray) -> float:
    """The centred anomaly correlation (Wilks 2011, 8.6.4), weighted by w."""
    ok = np.isfinite(f) & np.isfinite(a) & np.isfinite(c) & (w > 0)
    w, fa, aa = w[ok], f[ok] - c[ok], a[ok] - c[ok]
    fa = fa - (w * fa).sum() / w.sum()
    aa = aa - (w * aa).sum() / w.sum()
    norm = math.sqrt((w * fa * fa).sum() * (w * aa * aa).sum())
    return float((w * fa * aa).sum() / norm) if norm > 0 else math.nan


def great_circle_km(lon1, lat1, lon2, lat2) -> np.ndarray:
    """Great-circle distance (km) on the sphere of radius R1, by the haversine (Sinnott 1984). Well
    conditioned for short distances, the haversine loses precision near antipodes (to about 0.2 m);
    the sphere departs from the ellipsoid by up to about 0.5 % of the distance."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(np.subtract(lon2, lon1)) / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.minimum(1.0, np.sqrt(h)))


def lead_hours(start: dt.datetime, t: dt.datetime) -> float:
    return (t - start).total_seconds() / 3600


class Tally:
    """Values of cases, gathered by score (quantity, metric, region, model) and lead."""

    def __init__(self):
        self.cases: dict = {}

    def add(self, key: tuple, lead: int, value: float, valid: dt.datetime, run: Run | None, against: str,
            weight: int = 1, until: dt.datetime | None = None) -> None:
        """A case's value (a sum over `weight` reports, for a station score), valid at `valid` (to `until`)."""
        if not math.isfinite(value):
            return
        until = until or valid
        s = self.cases.setdefault(key, {"leads": {}, "against": set(), "runs": set(), "first": valid, "last": until,
                                        "products": set()})
        s["leads"].setdefault(lead, []).append((value, weight))
        s["against"].add(against)
        if run:
            s["runs"].add((run.model, run.start))
            s["products"].add(run.product)
        s["first"], s["last"] = min(s["first"], valid), max(s["last"], until)


# ---- 1. 500 hPa height ------------------------------------------------------------------

class Climatology:
    """Monthly means of 500 hPa height (NCEP/NCAR R1, 1991-2020), linear in time between the months' middles."""

    def __init__(self, folder: Path = CLIMATOLOGY):
        import netcdf3
        retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
        ds = netcdf3.Dataset.open(folder / retrieval["requests"][0]["file"])
        levels = list(np.asarray(ds.array("level"), np.float64))
        z = np.asarray(ds.array("hgt"), np.float64)[:, levels.index(500.0)]       # (12, lat, lon)
        lat, lon = np.asarray(ds.array("lat"), np.float64), np.asarray(ds.array("lon"), np.float64)
        if lat[0] < lat[-1]:
            lat, z = lat[::-1], z[:, ::-1]
        self.months = z
        self.grid = {"lon0": float(lon[0]), "dlon": float(lon[1] - lon[0]), "nlon": len(lon),
                     "lat0": float(lat[0]), "dlat": float(lat[0] - lat[1]), "nlat": len(lat), "periodic": True}
        self.retrieved_at = retrieval["retrieved_at"]
        self.cache: dict = {}

    @staticmethod
    def weights(t: dt.datetime) -> tuple[int, int, float]:
        """The two months (0-11) whose middles bracket t, and the weight of the later one."""
        def middle(year: int, month: int) -> dt.datetime:       # month 1-12; months beyond wrap the year
            year, month = year + (month - 1) // 12, (month - 1) % 12 + 1
            a = dt.datetime(year, month, 1, tzinfo=dt.timezone.utc)
            b = dt.datetime(year + (month == 12), month % 12 + 1, 1, tzinfo=dt.timezone.utc)
            return a + (b - a) / 2
        for m in range(0, 13):                 # middles of December before, then January to December
            a, b = middle(t.year, m), middle(t.year, m + 1)
            if a <= t < b:
                return (m - 1) % 12, m % 12, (t - a) / (b - a)
        raise ValueError(t)

    def on(self, t: dt.datetime, grid: dict) -> np.ndarray:
        key = (t, tuple(grid[k] for k in ("lon0", "dlon", "nlon", "lat0", "dlat", "nlat")))
        if key not in self.cache:
            i, j, w = self.weights(t)
            self.cache[key] = regrid((1 - w) * self.months[i] + w * self.months[j], self.grid, grid)
        return self.cache[key]


def analyses(runs: list[Run], quantity: str) -> dict:
    """For each valid time (ms), the analysis that is its truth: (run, time index), from the earliest run
    whose time it is, at or before the run's start; the GFS's first."""
    truth = {}
    for run in sorted(runs, key=lambda r: (r.analysis != "GFS", r.start)):
        v = find_variable(run.meta, quantity)
        if not v:
            continue
        for i, t in enumerate(run.times):
            if t <= run.start and run.ms[i] not in truth:
                truth[run.ms[i]] = (run, i)
    return truth


def score_heights(runs: list[Run], tally: Tally, climatology: Climatology | None, log=print) -> None:
    truth = analyses(runs, "z500")
    if not truth:
        log("  500 hPa height: no analyses")
        return
    first_truth = next(iter(truth.values()))[0]
    target = first_truth.grid
    weights = {region: band_weights(target, *band) for region, band in BANDS.items()}

    def analysis(ms: float) -> tuple[np.ndarray, str]:
        run, i = truth[ms]
        return regrid(run.field(find_variable(run.meta, "z500"), i), run.grid, target), \
            f"{run.analysis} analyses ({run.analysis_note})"

    seen, starts = set(), {}
    for run in runs:
        v = find_variable(run.meta, "z500")
        if not v or (run.model, run.start) in seen:
            continue
        seen.add((run.model, run.start))
        n = 0
        for i, t in enumerate(run.times):
            lead = lead_hours(run.start, t)
            if lead <= 0 or lead % LEAD_STEP_H or run.ms[i] not in truth:
                continue
            a, against = analysis(run.ms[i])
            f = regrid(run.field(v, i), run.grid, target)
            c = climatology.on(t, target) if climatology else None
            for region, w in weights.items():
                tally.add(("z500", "rmse", region, run.model), int(lead), rmse(f, a, w) ** 2, t, run, against)
                if c is not None:
                    tally.add(("z500", "acc", region, run.model), int(lead), acc(f, a, c, w), t, run, against)
            starts.setdefault(run.start.timestamp() * 1000, set()).add(int(lead))
            n += 1
        if n:
            log(f"  500 hPa height: {run.label()}, {n} lead(s)")
    # Persistence: the analysis at each scored run's start, held.
    for start_ms, leads in sorted(starts.items()):
        if start_ms not in truth:
            continue
        p, _ = analysis(start_ms)
        start = dt.datetime.fromtimestamp(start_ms / 1000, dt.timezone.utc)
        for lead in sorted(leads):
            ms = start_ms + lead * 3600000
            if ms not in truth:
                continue
            a, against = analysis(ms)
            t = start + dt.timedelta(hours=lead)
            c = climatology.on(t, target) if climatology else None
            for region, w in weights.items():
                tally.add(("z500", "rmse", region, "Persistence"), lead, rmse(p, a, w) ** 2, t, None, against)
                if c is not None:
                    tally.add(("z500", "acc", region, "Persistence"), lead, acc(p, a, c, w), t, None, against)


# ---- 2. stations ----------------------------------------------------------------------

def observation_products(products: Path = OUT) -> list[dict]:
    index = json.loads((products / "index.json").read_text(encoding="utf-8")) if (products / "index.json").exists() else {}
    files = [products / e["file"] for e in index.get("observations", []) if e.get("file")]
    files += [p for p in sorted(products.glob("observations-*.json")) if p not in files]
    out = []
    for path in files:
        if path.exists():
            p = json.loads(path.read_text(encoding="utf-8"))
            if p.get("schema", "").startswith("observations/"):
                out.append(p)
    return out


def station_reports(product: dict, quantity: str) -> dict | None:
    """The reports of a quantity in an observations product: lon, lat, ms, value arrays, one per
    station and instant, within the quantity's plausible range."""
    spec = STATIONS[quantity]
    var = next((v for v in product["variables"] if v.get("standard_name") == spec["standard_name"]
                and (not spec.get("level") or re.search(spec["level"], v.get("level") or ""))), None)
    if not var:
        return None
    fields = product["fields"]
    try:
        k = {name: fields.index(name) for name in ("id", "lon", "lat", "time")}
        kv = fields.index(var["name"])
    except ValueError:
        return None
    ke = fields.index("elevation") if "elevation" in fields else None
    lo, hi = spec["range"]
    highest = spec.get("highest", math.inf)
    seen, lon, lat, ms, value = set(), [], [], [], []
    for row in product["rows"]:
        x = row[kv]
        if x is None or row[k["lon"]] is None or row[k["lat"]] is None or not lo <= x <= hi:
            continue
        if ke is not None and (row[ke] or 0) > highest:
            continue
        key = (row[k["id"]], row[k["time"]])
        if key in seen:
            continue
        seen.add(key)
        lon.append(row[k["lon"]]); lat.append(row[k["lat"]]); value.append(x)
        ms.append(parse_time(row[k["time"]]).timestamp() * 1000)
    if not value:
        return None
    return {"name": product["name"], "variable": var, "lon": np.array(lon, np.float64), "lat": np.array(lat, np.float64),
            "ms": np.array(ms, np.float64), "value": np.array(value, np.float64)}


def score_stations(runs: list[Run], tally: Tally, products: Path = OUT, log=print) -> list[str]:
    used = []
    for product in observation_products(products):
        for quantity in STATIONS:
            reports = station_reports(product, quantity)
            if not reports:
                continue
            against = (f"{product.get('title') or product['name']} ({'observations-' + product['name']}), "
                       f"{reports['variable']['name']}")
            region = product.get("title") or product["name"]
            seen = set()
            for run in runs:
                v = find_variable(run.meta, quantity)
                if not v or (run.model, run.start) in seen:
                    continue
                f = sample_run(run, v, reports["lon"], reports["lat"], reports["ms"])
                lead = (reports["ms"] - run.start.timestamp() * 1000) / 3600000
                ok = np.isfinite(f) & (lead >= 0)
                if not ok.any():
                    continue
                seen.add((run.model, run.start))
                step = np.floor(lead / LEAD_STEP_H + 0.5).astype(int) * LEAD_STEP_H      # halves go up
                for h in np.unique(step[ok]):
                    sel = ok & (step == h)
                    d = f[sel] - reports["value"][sel]
                    first, last = (dt.datetime.fromtimestamp(float(x) / 1000, dt.timezone.utc)
                                   for x in (reports["ms"][sel].min(), reports["ms"][sel].max()))
                    for metric, value in (("rmse", float((d * d).sum())), ("bias", float(d.sum()))):
                        tally.add((quantity, metric, region, run.model), int(h), value, first, run, against,
                                  weight=int(sel.sum()), until=last)
                log(f"  {TITLES[quantity]} at {product['name']}: {run.label()}, {int(ok.sum())} report(s)")
            used.append(product["name"])
    return sorted(set(used))


# ---- 3. tropical cyclones --------------------------------------------------------------------

def storms_of(path: Path) -> list[dict]:
    """The observed positions (and central pressures, hPa, where given) of each storm of a tracks product."""
    p = json.loads(path.read_text(encoding="utf-8"))
    fields = p["fields"]
    k_obs = fields.index("observed") if "observed" in fields else None
    k_p = next((fields.index(f) for f in fields if re.search(r"pres|mslp", f, re.I)), None)
    # Central pressures in Pa or hPa, told apart by their size (about 1e5 or 1e3).
    given = [q[k_p] for tr in p["tracks"] for q in tr["points"] if k_p is not None and q[k_p] is not None]
    scale = 0.01 if given and float(np.median(given)) > 2000 else 1.0
    out = []
    for tr in p["tracks"]:
        pts = sorted((q for q in tr["points"] if k_obs is None or q[k_obs]), key=lambda q: q[0])
        if not pts:
            continue
        out.append({"id": tr.get("id") or tr.get("name"), "name": tr.get("name"),
                    "t": np.array([q[0] for q in pts], np.float64), "lon": np.array([q[1] for q in pts], np.float64),
                    "lat": np.array([q[2] for q in pts], np.float64),
                    "hpa": np.array([q[k_p] * scale if k_p is not None and q[k_p] is not None else np.nan
                                     for q in pts], np.float64)})
    return out


def observed_storms(products: Path = OUT) -> tuple[list[dict], str]:
    """The storms observed: IBTrACS's best tracks, each extended past its last analysis by GDACS's observed
    positions of the storm of the same name (IBTrACS lags by days; GDACS relays the warning centres' positions
    as they are issued), with GDACS's storms that IBTrACS does not hold yet; GDACS's alone without IBTrACS."""
    index = json.loads((products / "index.json").read_text(encoding="utf-8")) if (products / "index.json").exists() else {}
    entries = [e for e in index.get("tracks", []) if (products / e["file"]).exists()]
    storms = {}
    for e in entries:
        if "ibtracs" in e["name"]:
            for s in storms_of(products / e["file"]):
                if s["id"] not in storms or len(s["t"]) > len(storms[s["id"]]["t"]):
                    storms[s["id"]] = s
    # Names recur from year to year: a GDACS storm joins the latest IBTrACS storm of its name, if that one
    # was observed within GAP_DAYS of the GDACS storm's first position.
    named = {}
    for s in sorted(storms.values(), key=lambda s: s["t"][-1]):
        named[str(s["name"]).upper()] = s
    extended = added = 0
    for e in entries:
        if "gdacs" not in e["name"]:
            continue
        for g in storms_of(products / e["file"]):
            s = named.get(re.sub(r"-\d{2}$", "", str(g["name"])).upper())
            if s is not None and g["t"][0] - s["t"][-1] > GAP_DAYS * 86400000:
                s = None
            if s is None:
                storms[g["id"]] = g
                added += 1
                continue
            later = g["t"] > s["t"][-1]
            if later.any():
                for k in ("t", "lon", "lat", "hpa"):
                    s[k] = np.concatenate([s[k], g[k][later]])
                extended += 1
    has_ibtracs = any("ibtracs" in e["name"] for e in entries)
    source = ("IBTrACS" + (f", extended by GDACS's observed positions ({extended} storm(s) extended, {added} added)"
                           if extended or added else "")) if has_ibtracs else "GDACS's observed positions"
    return [s for s in storms.values() if len(s["t"]) >= 2], source


def storm_at(storm: dict, t: float) -> tuple[float, float, float] | None:
    """A storm's position and pressure at an instant (ms), linear between its positions (across the date line)."""
    ts = storm["t"]
    if t < ts[0] or t > ts[-1]:
        return None
    k = min(max(int(np.searchsorted(ts, t, side="right")) - 1, 0), len(ts) - 2) if len(ts) > 1 else 0
    if len(ts) == 1 or ts[k + 1] == ts[k]:
        return float(storm["lon"][k]), float(storm["lat"][k]), float(storm["hpa"][k])
    w = (t - ts[k]) / (ts[k + 1] - ts[k])
    dl = storm["lon"][k + 1] - storm["lon"][k]
    dl = dl - 360 if dl > 180 else dl + 360 if dl < -180 else dl
    lon = (storm["lon"][k] + dl * w + 540) % 360 - 180
    return float(lon), float(storm["lat"][k] + (storm["lat"][k + 1] - storm["lat"][k]) * w), \
        float(storm["hpa"][k] + (storm["hpa"][k + 1] - storm["hpa"][k]) * w)


def model_tracks(runs: list[Run], log=print) -> dict:
    """The cyclones of each run's sea-level pressure, from track_cyclones.mjs; cached beside archived runs,
    and found again if the tracker's code has changed since."""
    node = shutil.which("node")
    script = Path(__file__).with_name("track_cyclones.mjs")
    tracker = hashlib.sha256(b"".join(p.read_bytes() for p in (
        script, ROOT / "domains" / "meteorology" / "js" / "cyclones.js"))).hexdigest()[:16]
    out, job = {}, []
    for run in runs:
        v = find_variable(run.meta, "mslp")
        if not v:
            continue
        cache = run.path.parent / f"tracks-{run.product}.json" if run.origin == "archive" else None
        if cache and cache.exists():
            kept = json.loads(cache.read_text(encoding="utf-8"))
            if kept.get("tracker") == tracker and kept.get("variable") == v["name"]:
                out[run.path] = kept["tracks"]
                continue
        job.append({"path": str(run.path), "variable": v["name"], "cache": cache})
    if job and not node:
        log("  cyclones: Node is not on the path; runs not yet tracked are left out")
        return out
    if job:
        request = {"relief": str(RELIEF) if RELIEF.exists() else None,
                   "products": [{"path": j["path"], "variable": j["variable"]} for j in job]}
        r = subprocess.run([node, str(script)], input=json.dumps(request),
                           cwd=ROOT, capture_output=True, text=True, encoding="utf-8", timeout=1800)
        if r.returncode != 0:
            log(f"  cyclones: the tracker failed: {r.stderr.strip()[-400:]}")
            return out
        found = json.loads(r.stdout)
        for j in job:
            res = found.get(j["path"], {})
            if "tracks" not in res:
                log(f"  cyclones: {j['path']}: {res.get('error')}")
                continue
            out[Path(j["path"])] = res["tracks"]
            if j["cache"]:
                j["cache"].write_text(json.dumps({"variable": j["variable"], "tracker": tracker,
                                                  "tracks": res["tracks"]}) + "\n", encoding="utf-8")
    return out


def match(run: Run, tracks: list[dict], storms: list[dict]) -> list[tuple[dict, dict]]:
    """Pairs (storm, model track): at their first common time at or after the run's start, within MATCH_DEG
    of arc; nearest pairs first, one track to one storm."""
    start = run.start.timestamp() * 1000
    pairs = []
    for i, tr in enumerate(tracks):
        for j, s in enumerate(storms):
            for p in tr["points"]:
                if p["t"] < start:
                    continue
                q = storm_at(s, p["t"])
                if q:
                    d = float(great_circle_km(p["lon"], p["lat"], q[0], q[1])) / EARTH_RADIUS_KM * 180 / math.pi
                    if d <= MATCH_DEG:
                        pairs.append((d, i, j))
                    break
    pairs.sort()
    used_t, used_s, out = set(), set(), []
    for _, i, j in pairs:
        if i in used_t or j in used_s:
            continue
        used_t.add(i); used_s.add(j)
        out.append((storms[j], tracks[i]))
    return out


def score_cyclones(runs: list[Run], tally: Tally, products: Path = OUT, log=print) -> str | None:
    storms, source = observed_storms(products)
    if not storms:
        log("  cyclones: no observed storms")
        return None
    tracks = model_tracks(runs, log)
    region = "worldwide"
    against = f"{source}: centres" + (", and IBTrACS's central pressures" if "IBTrACS" in source else "")
    seen = set()
    for run in runs:
        if run.path not in tracks or (run.model, run.start) in seen:
            continue
        seen.add((run.model, run.start))
        pairs = match(run, tracks[run.path], storms)
        for storm, tr in pairs:
            for p in tr["points"]:
                t = dt.datetime.fromtimestamp(p["t"] / 1000, dt.timezone.utc)
                lead = lead_hours(run.start, t)
                q = storm_at(storm, p["t"])
                if lead < 0 or lead % LEAD_STEP_H or not q or p["t"] > run.ms[-1]:
                    continue
                km = float(great_circle_km(p["lon"], p["lat"], q[0], q[1]))
                tally.add(("cyclones", "track_km", region, run.model), int(lead), km, t, run, against)
                if math.isfinite(q[2]):
                    e = p["pa"] / 100 - q[2]
                    tally.add(("cyclones", "pressure_hpa", region, run.model), int(lead), abs(e), t, run, against)
                    tally.add(("cyclones", "bias", region, run.model), int(lead), e, t, run, against)
        if pairs:
            log(f"  cyclones: {run.label()}: {', '.join(str(s['name']) for s, _ in pairs)}")
    return source


# ---- the product ----------------------------------------------------------------------------

UNITS = {("z500", "rmse"): "m", ("z500", "acc"): "1", ("t2m", "rmse"): "K", ("t2m", "bias"): "K",
         ("mslp", "rmse"): "hPa", ("mslp", "bias"): "hPa", ("hs", "rmse"): "m", ("hs", "bias"): "m",
         ("cyclones", "track_km"): "km", ("cyclones", "pressure_hpa"): "hPa", ("cyclones", "bias"): "hPa"}
DEFINITIONS = {
    ("z500", "rmse"): "root-mean-square error, cos(latitude) weighted; a lead's value is the root of its cases' mean squared error",
    ("z500", "acc"): "anomaly correlation, centred (Wilks 2011, 8.6.4), cos(latitude) weighted, anomalies from NCEP/NCAR R1 1991-2020 monthly means interpolated to the day; the mean of the cases",
    ("t2m", "rmse"): "root-mean-square of forecast minus report, over all reports at the lead (to the nearest 6 h); no height correction",
    ("t2m", "bias"): "mean of forecast minus report; no height correction (6.5 K per km of height between station and model ground is not applied)",
    ("mslp", "rmse"): "root-mean-square of forecast minus report", ("mslp", "bias"): "mean of forecast minus report",
    ("hs", "rmse"): "root-mean-square of forecast minus report", ("hs", "bias"): "mean of forecast minus report",
    ("cyclones", "track_km"): "mean great-circle distance between the model's and the observed centre",
    ("cyclones", "pressure_hpa"): "mean absolute error of the central pressure",
    ("cyclones", "bias"): "mean of the model's minus the observed central pressure",
}
DIGITS = {"m": 1, "1": 3, "K": 2, "hPa": 1, "km": 0}


def assemble(tally: Tally) -> list[dict]:
    scores = []
    for (quantity, metric, region, model), s in sorted(tally.cases.items(), key=lambda kv: (
            list(TITLES).index(kv[0][0]), kv[0][1], kv[0][2], kv[0][3] == "Persistence", kv[0][3])):
        units = UNITS[(quantity, metric)]
        leads = []
        for hours in sorted(s["leads"]):
            values = s["leads"][hours]
            count = sum(w for _, w in values)
            total = sum(v for v, _ in values)
            mean = total / count
            if metric == "rmse":
                mean = math.sqrt(mean)
            if quantity == "mslp":
                mean /= 100.0            # Pa to hPa
            leads.append({"hours": hours, "value": round(mean, DIGITS[units]), "count": count})
        scores.append({
            "quantity": quantity, "title": TITLES[quantity], "metric": metric, "region": region, "units": units,
            "model": model, "against": "; ".join(sorted(s["against"])), "definition": DEFINITIONS[(quantity, metric)],
            "runs": len(s["runs"]) if model != "Persistence" else None,
            "products": sorted(s["products"]) or None,
            "window": {"start": iso(s["first"]), "end": iso(s["last"])},
            "leads": leads,
        })
    return scores


def runs_summary(runs: list[Run]) -> dict:
    """For each product, its model, how many of its runs are archived or current, and their first and last starts."""
    out = {}
    for r in runs:
        e = out.setdefault(r.product, {"model": r.model, "archived": 0, "current": 0, "first": iso(r.start),
                                       "last": iso(r.start)})
        e["current" if r.origin == "current" else "archived"] += 1
        e["first"], e["last"] = min(e["first"], iso(r.start)), max(e["last"], iso(r.start))
    return out


def build(archive: Path = ARCHIVE, products: Path = OUT, cyclones: bool = True, log=print) -> dict:
    runs = runs_available(archive, products)
    log(f"{len(runs)} run(s): {sum(r.origin == 'archive' for r in runs)} archived, "
        f"{sum(r.origin == 'current' for r in runs)} current")
    tally = Tally()
    notes = []
    try:
        climatology = Climatology()
    except (OSError, ValueError, KeyError) as error:
        climatology = None
        notes.append(f"No anomaly correlation: the climatology could not be read ({type(error).__name__}); "
                     "fetch it with: python data/acquisition/fetch.py ncep-r1-climatology")
    score_heights(runs, tally, climatology, log)
    stations = score_stations(runs, tally, products, log)
    source = score_cyclones(runs, tally, products, log) if cyclones else None
    scores = assemble(tally)
    archived = sum(r.origin == "archive" for r in runs)
    notes += [
        "Lead: the hours from a run's start (the instant of the analysis it began from) to the instant foretold. "
        "A run's start is stated by the product, or read from its source, or estimated (data/products/archive.py).",
        "500 hPa height is scored against analyses: at each valid time, the analysis in the earliest run whose "
        "time it is, at or before that run's start, the GFS's first (including the GFS analysis that began a "
        "Pangu-Weather run); where no GFS analysis of that time is held, another model's (ECMWF's first time). "
        "Each score lists the truths it used. The GFS is thus scored mostly against its own analyses, which "
        "favours it at short leads.",
        "Station scores sample the forecast bilinearly at each station, linearly in time between its times, "
        "within its span; each report's lead is rounded to the nearest 6 h. No height correction is applied: "
        "the products carry no model orography.",
        "Persistence is the analysis at a run's start, held unchanged: the forecast a model must beat.",
    ]
    if not archived:
        notes.append("No archived run yet: these scores use only what the current products contain, runs that "
                     "began before the latest analysis held (a hindcast, or an older product beside a newer one), "
                     "at the leads already analysed. They rest on one or two runs; see each score's count and window.")
    if cyclones and source:
        notes.append(f"Tropical cyclones found by the meteorology domain's tracker in each run's sea-level pressure "
                     f"and matched to {source} within {MATCH_DEG:g} degrees at their first common time. Central "
                     "pressures of a 1 degree field are well above the true ones.")
    product = {
        "schema": "verification/1", "name": NAME, "title": "Scores of the Earth Atlas's forecasts",
        "updated": iso(dt.datetime.now(dt.timezone.utc)), "time_scale": "UTC",
        "scores": scores, "notes": notes,
        "sources": {
            "runs": runs_summary(runs),
            "climatology": {"dataset": "ncep-r1-climatology", "period": "1991-01-01/2020-12-31",
                            "credit": "NCEP/NCAR Reanalysis 1, NOAA PSL",
                            "retrieved_at": climatology.retrieved_at} if climatology else None,
            "observations": stations, "tracks": source,
        },
    }
    path = products / f"verification-{NAME}.json"
    path.write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    if products == OUT:
        register("verification", {"name": NAME, "file": path.name, "title": product["title"], "updated": product["updated"],
                                  "count": len(scores), "models": sorted({s["model"] for s in scores}),
                                  "quantities": sorted({s["quantity"] for s in scores})})
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no-cyclones", action="store_true", help="leave out the cyclones (they need Node)")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace", line_buffering=True)
    product = build(cyclones=not args.no_cyclones)
    print(f"verification-{NAME}: {len(product['scores'])} score(s)")
    for s in product["scores"]:
        row = "  ".join(f"{l['hours']}h {l['value']:g}" for l in s["leads"][:8])
        print(f"  {s['title']:24s} {s['metric']:12s} {s['region']:22s} {s['model']:26s} {row}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
