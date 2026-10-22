"""Checks of the archive, the verification and the scheduler: the scores against closed forms on
synthetic fields, station sampling against a linear field, the great-circle distance, the tracker
end to end on a synthetic cyclone, the archive's naming and retention in a temporary folder, and
the scheduler's due times; then, where the cache has them, the scores against the Pangu-Weather
hindcast's own record and the scorecard's shape.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

from checks import PRODUCTS, ROOT, assert_, check, skip

import archive  # noqa: E402
import verification as v  # noqa: E402
import schedule  # noqa: E402

UTC = dt.timezone.utc
T = lambda text: dt.datetime.fromisoformat(text.replace("Z", "+00:00"))  # noqa: E731
GRID = {"lon0": -180.0, "dlon": 1.0, "nlon": 360, "lat0": 90.0, "dlat": 1.0, "nlat": 181, "periodic": True}
lon, lat = v.lattice(GRID)
phi, lam = np.radians(lat), np.radians(lon)

# ---- 1. scores against closed forms ------------------------------------------------------
c = 5500 + 300 * np.cos(phi)                                   # a climatology
a = c + 80 * np.cos(phi) * np.sin(3 * lam) + 40 * np.sin(2 * phi)   # an analysis
for region, band in v.BANDS.items():
    w = v.band_weights(GRID, *band)
    check(f"z500 {region}: a field against itself has RMSE 0", v.rmse(a, a, w), 1e-12, "m")
    check(f"z500 {region}: a field against itself has ACC 1", abs(v.acc(a, a, c, w) - 1), 1e-12)
    check(f"z500 {region}: the negative anomaly has ACC -1", abs(v.acc(c - (a - c), a, c, w) + 1), 1e-12)
    check(f"z500 {region}: an offset of 10 m gives RMSE 10 m", abs(v.rmse(a + 10, a, w) - 10), 1e-9, "m")
    check(f"z500 {region}: an offset of 10 m gives bias 10 m", abs(v.bias(a + 10, a, w) - 10), 1e-9, "m")
    # The centred form: an anomaly scaled and shifted correlates perfectly (the uncentred form would not).
    check(f"z500 {region}: ACC is centred (2a' + 25 m has ACC 1)", abs(v.acc(c + 2 * (a - c) + 25, a, c, w) - 1), 1e-12)
# cos(latitude) is an area weight: an error of sin(latitude) over 20-80 N has the mean square
# (sin^3 80 - sin^3 20) / (3 (sin 80 - sin 20)) of the integral over the band's area; the sum over
# 1 degree rows, which counts the band's two end rows whole, departs from it by O(1 degree), 0.4 %.
s1, s2 = math.sin(math.radians(20)), math.sin(math.radians(80))
exact = math.sqrt((s2 ** 3 - s1 ** 3) / (3 * (s2 - s1)))
check("z500: the cos-weighted RMSE of an error sin(lat) over 20-80 N is the area integral's",
      abs(v.rmse(a + np.sin(phi), a, v.band_weights(GRID, 20, 80)) - exact), 5e-3)

# ---- 2. station sampling against a linear field -------------------------------------------
class Stub:
    """A run whose field is linear in longitude, latitude and time: bilinear and linear sampling are exact."""
    def __init__(self, grid, times):
        self.grid, self.times = grid, times
        self.ms = np.array([t.timestamp() * 1000 for t in times])
        self.glon, self.glat = v.lattice(grid)

    def field(self, variable, k):
        return 280 + 0.05 * self.glon - 0.3 * self.glat + 1.5 * k

regional = {"lon0": -20.0, "dlon": 0.5, "nlon": 81, "lat0": 60.0, "dlat": 0.5, "nlat": 61, "periodic": False}
times = [T("2026-10-22T00:00:00Z") + dt.timedelta(hours=6 * k) for k in range(5)]
run = Stub(regional, times)
rng = np.random.default_rng(1)
plon, plat = rng.uniform(-19.9, 19.9, 200), rng.uniform(30.1, 59.9, 200)
pms = rng.uniform(run.ms[0], run.ms[-1], 200)
expected = 280 + 0.05 * plon - 0.3 * plat + 1.5 * (pms - run.ms[0]) / 21600000
got = v.sample_run(run, None, plon, plat, pms)
check("stations: a linear field is sampled exactly (bilinear in space, linear in time)",
      float(np.max(np.abs(got - expected))), 1e-9, "K")
outside = v.sample_run(run, None, np.array([0.0, 0.0, 25.0]), np.array([45.0, 45.0, 45.0]),
                       np.array([run.ms[0] - 1, run.ms[-1] + 1, run.ms[0]]))
assert_("stations: no value before or after the forecast's span, nor outside its grid",
        bool(np.isnan(outside).all()), str(outside))
hole = run.field(None, 0).copy()
hole[0, 0] = np.nan                     # the north-west corner, a weight below one half for a place near the centre
x = v.bilinear_points(hole, regional, np.array([-19.6]), np.array([59.6]))[0]
check("stations: a missing neighbour's weight is shared among the others (physis's rule)",
      abs(x - (280 + 0.05 * -19.6 - 0.3 * 59.6)), 0.05, "K")
check("stations: a periodic grid wraps across 180 degrees",
      abs(v.bilinear_points(np.repeat(lat[:, :1], 360, axis=1), GRID, np.array([179.6]), np.array([10.25]))[0] - 10.25), 1e-9)

# ---- 3. the great-circle distance -------------------------------------------------------
R = v.EARTH_RADIUS_KM
check("cyclones: from the equator to the pole is a quarter of the great circle, pi R / 2",
      abs(float(v.great_circle_km(0, 0, 0, 90)) - math.pi * R / 2), 1e-9, "km")
check("cyclones: one degree of the equator is 2 pi R / 360 = 111.195 km",
      abs(float(v.great_circle_km(179.5, 0, -179.5, 0)) - 2 * math.pi * R / 360), 1e-9, "km")
# Near antipodes the haversine is ill-conditioned (arcsin at 1): a metre is its precision there.
check("cyclones: antipodes are pi R apart", abs(float(v.great_circle_km(30, 10, -150, -10)) - math.pi * R), 1e-3, "km")


def write_product(folder: Path, name: str, times: list[str], variables: dict, retrieved: str, model: str = "Test") -> Path:
    """A field product as field_grid.py writes it, from {name: (standard_name, level, values (nt, nlat, nlon))}."""
    folder.mkdir(parents=True, exist_ok=True)
    blocks, entries, offset = [], [], 0
    for vname, (standard, level, values) in variables.items():
        lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
        scale, off = (hi - lo) / 65000 or 1.0, (lo + hi) / 2
        q = np.round((values - off) / scale).astype("<i2")
        offsets = []
        for k in range(len(times)):
            offsets.append(offset)
            blocks.append(q[k].tobytes())
            offset += q[k].nbytes
        entries.append({"name": vname, "standard_name": standard, "level": level, "units": "1", "scale": scale,
                        "offset": off, "offsets": offsets})
    (folder / f"fields-{name}.bin").write_bytes(b"".join(blocks))
    meta = {"schema": "field-grid/1", "name": name, "title": name, "domain": "meteorology", "model": model,
            "source": {"retrieved_at": retrieved}, "grid": {**GRID, "crs": "", "order": ""}, "times": times,
            "time_scale": "UTC", "variables": entries, "vectors": [], "file": f"fields-{name}.bin"}
    path = folder / f"fields-{name}.json"
    path.write_text(json.dumps(meta), encoding="utf-8")
    return path


# ---- 4. the tracker, end to end, on a synthetic cyclone ------------------------------------------
if not shutil.which("node"):
    skip("cyclones: Node is not on the path")
else:
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        stimes = [f"2026-10-22T{h:02d}:00:00Z" for h in (0, 6, 12, 18)] + [f"2026-10-23T{h:02d}:00:00Z" for h in (0, 6, 12, 18)]
        centres = [(140.0 - 1.2 * k, 15.0 + 0.4 * k) for k in range(len(stimes))]
        p = np.stack([101300 - 4000 * np.exp(-(v.great_circle_km(clon, clat, lon, lat) / 250.0) ** 2)
                      for clon, clat in centres])
        path = write_product(tmp, "synthetic", stimes, {"prmsl": ("air_pressure_at_sea_level", None, p)},
                             "2026-10-22T00:00:00Z")
        r = v.Run(path, "current")
        r.start = T(stimes[0])
        tracks = v.model_tracks([r], log=lambda *_: None).get(path, [])
        storm = {"id": "S", "name": "SYNTHETIC", "t": r.ms.copy(), "lon": np.array([c[0] for c in centres]),
                 "lat": np.array([c[1] for c in centres]), "hpa": np.full(len(centres), 973.0)}
        pairs = v.match(r, tracks, [storm])
        assert_("cyclones: the tracker finds the synthetic cyclone and it is matched", len(pairs) == 1,
                f"{len(tracks)} track(s)")
        if pairs:
            _, tr = pairs[0]
            km = [float(v.great_circle_km(q["lon"], q["lat"], *v.storm_at(storm, q["t"])[:2])) for q in tr["points"]]
            check("cyclones: the found centres lie within 30 km of the true ones (on a 1 degree grid)", max(km), 30, "km")
            check("cyclones: the found central pressure is the true one", max(abs(q["pa"] / 100 - 973) for q in tr["points"]),
                  2.0, "hPa")
            assert_("cyclones: the track lasts the whole run", len(tr["points"]) == len(stimes), str(len(tr["points"])))

# ---- 5. the archive: names, rule, idempotence, retention -----------------------------------------
assert_("archive: a run is named compactly", archive.run_name(T("2026-10-22T06:00:00Z")) == "2026-10-22T06Z"
        and archive.run_name(T("2026-10-22T06:30:00Z")) == "2026-10-22T0630Z")
assert_("archive: a run's name reads back as its start",
        archive.parse_run_name("2026-10-22T0630Z") == T("2026-10-22T06:30:00Z")
        and archive.parse_run_name("2026-10-22T18Z") == T("2026-10-22T18:00:00Z"))
six = [T("2026-10-20T18:00:00Z") + dt.timedelta(hours=6 * k) for k in range(24)]
assert_("archive: the start is the latest time at or before the retrieval minus 4 h",
        archive.start_by_rule(six, T("2026-10-21T17:03:24Z")) == T("2026-10-21T12:00:00Z"),
        archive.start_by_rule(six, T("2026-10-21T17:03:24Z")))
later = [T("2026-10-22T00:00:00Z") + dt.timedelta(hours=6 * k) for k in range(28)]
assert_("archive: a product that begins after its run has its time step extended back to it",
        archive.start_by_rule(later, T("2026-10-21T19:34:30Z")) == T("2026-10-21T12:00:00Z"),
        archive.start_by_rule(later, T("2026-10-21T19:34:30Z")))
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    products, store = tmp / "products", tmp / "archive"
    times = [f"2026-10-22T{h:02d}:00:00Z" for h in (0, 6, 12, 18)] + ["2026-10-23T00:00:00Z", "2026-10-23T06:00:00Z"]
    z = np.stack([a + k for k in range(len(times))])
    t2 = np.stack([280 + 0 * a + k for k in range(len(times))])
    variables = {"z": ("geopotential_height", "500 hPa", z), "cape": ("atmosphere_convective_available_potential_energy_wrt_surface", None, t2)}
    write_product(products, "test-forecast", times, variables, "2026-10-22T10:30:00Z")
    first = archive.snapshot("test-forecast", products, store, data=tmp)
    again = archive.snapshot("test-forecast", products, store, data=tmp)
    assert_("archive: retrieved at 10:30, the run is 06Z (06:30 at the latest), kept in <product>/<run>/",
            first is not None and first["run"] == "2026-10-22T06Z"
            and (store / "test-forecast" / "2026-10-22T06Z" / "fields-test-forecast.bin").exists(), first)
    assert_("archive: archiving the same retrieval again changes nothing", again is None)
    write_product(products, "test-forecast", times, variables, "2026-10-22T10:00:00Z")
    replaced = archive.snapshot("test-forecast", products, store, data=tmp)
    kept = json.loads((store / "test-forecast" / "2026-10-22T06Z" / "fields-test-forecast.json").read_text(encoding="utf-8"))
    assert_("archive: a later retrieval of the same run replaces it",
            replaced is not None and kept["source"]["retrieved_at"] == "2026-10-22T10:00:00Z")
    write_product(products, "test-forecast", [t.replace("2026-10-2", "2026-10-1") for t in times], variables,
                  "2026-10-12T15:00:00Z")
    archive.snapshot("test-forecast", products, store, data=tmp)
    index = archive.rebuild_index(store)
    assert_("archive: the index lists both runs, in order",
            [r["run"] for r in index["products"]["test-forecast"]["runs"]] == ["2026-10-12T06Z", "2026-10-22T06Z"])
    removed = archive.prune(store, now=T("2026-10-23T00:00:00Z"), days=10)
    assert_("archive: retention removes the runs that started more than 10 days ago, and only those",
            [p.name for p in removed] == ["2026-10-12T06Z"]
            and (store / "test-forecast" / "2026-10-22T06Z").exists(), [p.name for p in removed])
    meta = json.loads((products / "fields-test-forecast.json").read_text(encoding="utf-8"))
    sub, codes = archive.subset(meta, (products / meta["file"]).read_bytes(), {"z"})
    full = np.frombuffer((products / meta["file"]).read_bytes(), "<i2")
    n = GRID["nlon"] * GRID["nlat"]
    same = all(np.array_equal(np.frombuffer(codes, "<i2")[sub["variables"][0]["offsets"][k] // 2:][:n],
                              full[meta["variables"][0]["offsets"][k] // 2:][:n]) for k in range(len(times)))
    assert_("archive: keeping only the scored variables keeps their values, renumbered",
            [w["name"] for w in sub["variables"]] == ["z"] and len(codes) == len(times) * n * 2 and same)

# ---- 6. the scheduler: due times, recipe steps, failures, the lock ---------------------------------
now = T("2026-10-22T15:05:00Z")
assert_("schedule: '1h' + 10 min was last due at 14:10", schedule.latest_slot("1h", 10, now) == T("2026-10-22T14:10:00Z"))
assert_("schedule: '1h' + 10 min is next due at 15:10", schedule.next_slot("1h", 10, now) == T("2026-10-22T15:10:00Z"))
assert_("schedule: a slot is due at its instant exactly",
        schedule.latest_slot("1h", 10, T("2026-10-22T15:10:00Z")) == T("2026-10-22T15:10:00Z"))
assert_("schedule: '6h' + 300 min falls at 05, 11, 17 and 23 UTC",
        [schedule.latest_slot("6h", 300, T(f"2026-10-22T{h}Z")).strftime("%d %H:%M") for h in ("04:59:00", "05:00:00", "16:59:59", "23:30:00")]
        == ["21 23:00", "22 05:00", "22 11:00", "22 23:00"])
assert_("schedule: '1d' + 420 min falls at 07:00 UTC each day",
        schedule.latest_slot("1d", 420, T("2026-10-22T06:00:00Z")) == T("2026-10-21T07:00:00Z")
        and schedule.next_slot("1d", 420, T("2026-10-22T06:00:00Z")) == T("2026-10-22T07:00:00Z"))
job = {"name": "j", "every": "6h", "offset_minutes": 300}
assert_("schedule: a job never run is due", schedule.is_due(job, None, now))
assert_("schedule: a job run since its latest instant is not due", not schedule.is_due(job, "2026-10-22T11:00:05Z", now))
assert_("schedule: a job last run before its latest instant is due", schedule.is_due(job, "2026-10-22T10:59:00Z", now))
assert_("schedule: missed instants give one run (due once, then not)",
        schedule.is_due(job, "2026-10-20T05:00:00Z", now) and not schedule.is_due(job, "2026-10-22T15:05:00Z", now))
assert_("schedule: a disabled job is never due", not schedule.is_due({**job, "enabled": False}, None, now))
try:
    schedule.period_seconds("6 hours")
    assert_("schedule: a malformed period is refused", False)
except ValueError:
    assert_("schedule: a malformed period is refused", True)
recipe = ["data/acquisition/recipe.py", "applications/earth-atlas/recipe.json"]
cmd, unknown = schedule.known_steps([*recipe, "gfs", "no-such-step", "fields-gfs"])
assert_("schedule: a recipe command keeps the steps its recipe has and drops the others",
        cmd == [*recipe, "gfs", "fields-gfs"] and unknown == ["no-such-step"], (cmd, unknown))
assert_("schedule: a recipe command with no known step is skipped",
        schedule.known_steps([*recipe, "no-such-step"]) == (None, ["no-such-step"]))
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    (tmp / "ok.py").write_text("import sys; open(sys.argv[1], 'a').write('ran\\n')\n", encoding="utf-8")
    (tmp / "fail.py").write_text("import sys; sys.exit(3)\n", encoding="utf-8")
    (tmp / "slow.py").write_text("import time; time.sleep(30)\n", encoding="utf-8")
    marks = tmp / "marks.txt"
    state = tmp / "state.json"
    jobs = {"name": "t", "every": "1h", "timeout_minutes": 0.05,
            "run": [["fail.py"], ["missing.py"], ["slow.py"], ["ok.py", str(marks)]]}
    lines = []
    status = schedule.run_job(tmp / "schedule.json", jobs, root=tmp, state_path=state, say=lines.append)
    entry = json.loads(state.read_text(encoding="utf-8"))["schedule.json"]["t"]
    assert_("schedule: a job goes on past failing, missing and overlong commands",
            status == "failed" and marks.exists() and len(entry["failed"]) == 3, (status, entry))
    assert_("schedule: the state records the job's start, end and status",
            entry["status"] == "failed" and entry["last_start"] <= entry["last_end"])
    again = schedule.run_job(tmp / "schedule.json", jobs, root=tmp, state_path=state, say=lines.append)
    assert_("schedule: a job just run is not run again by a second runner", again is None)
    one, two = schedule.Lock(tmp / "x.lock"), schedule.Lock(tmp / "x.lock")
    held = one.acquire(wait=False)
    refused = not two.acquire(wait=False)
    one.release()
    after = two.acquire(wait=False)
    two.release()
    assert_("schedule: the lock admits one holder at a time", held and refused and after)

# ---- 7. the climatology ------------------------------------------------------------------------
i, j, w = v.Climatology.weights(T("2027-01-01T00:00:00Z"))
assert_("climatology: on 1 January the means of December and January weigh one half each",
        (i, j) == (11, 0) and abs(w - 0.5) < 1e-12, (i, j, w))
i, j, w = v.Climatology.weights(T("2026-07-16T12:00:00Z"))
assert_("climatology: at the middle of July, July alone", (i, j) == (6, 7) and w == 0, (i, j, w))
try:
    clim = v.Climatology()
except (OSError, ValueError, KeyError):
    clim = None
    skip("climatology: not fetched (python data/acquisition/fetch.py ncep-r1-climatology)")
if clim:
    january = clim.on(T("2026-01-16T12:00:00Z"), GRID)
    band = v.band_weights(GRID, 20, 80)
    mean = float((january * band).sum() / band.sum())
    check("climatology: January's 500 hPa height over 20-80 N is near 5500 m (NCEP/NCAR R1, 1991-2020)",
          abs(mean - 5500), 30, "m")
    check("climatology: at the equator it is near 5860 m all year", abs(float(january[90].mean()) - 5860), 25, "m")

# ---- 8. the real products, where the cache has them ---------------------------------------------
hindcast = PRODUCTS / "fields-pangu-hindcast.json"
levels = PRODUCTS / "fields-gfs-levels-forecast.json"
if hindcast.exists() and levels.exists() and json.loads(hindcast.read_text(encoding="utf-8")).get("verification"):
    record = json.loads(hindcast.read_text(encoding="utf-8"))["verification"]
    runs = [v.Run(hindcast, "current"), v.Run(levels, "current")]
    tally = v.Tally()
    v.score_heights(runs, tally, None, log=lambda *_: None)
    mine = {l["hours"]: l["value"] for s in v.assemble(tally) if s["model"].startswith("Pangu") and s["region"] == "20–80° N"
            and s["metric"] == "rmse" for l in s["leads"]}
    theirs = {l["hours"]: l["model"] for l in record["leads"] if l["against"] == "analysis" and l["hours"] in mine}
    if theirs:
        check(f"z500: the Pangu hindcast's RMSE agrees with ai_forecast.py's own record at {len(theirs)} leads",
              max(abs(mine[h] - theirs[h]) for h in theirs), 0.15, "m")
    else:
        skip("z500: the hindcast and the GFS analyses have no lead in common")
else:
    skip("z500: no Pangu hindcast with its own verification to compare with")

card = PRODUCTS / f"verification-{v.NAME}.json"
if card.exists():
    p = json.loads(card.read_text(encoding="utf-8"))
    scores = p["scores"]
    shape = all({"quantity", "metric", "region", "units", "model", "against", "leads"} <= set(s) for s in scores)
    assert_(f"scorecard: {len(scores)} scores, each with its quantity, metric, region, units, model, truth and leads",
            p["schema"] == "verification/1" and shape and scores)
    assert_("scorecard: leads ascend in steps of 6 h, each with a count",
            all([l["hours"] for l in s["leads"]] == sorted(l["hours"] for l in s["leads"])
                and all(l["hours"] % 6 == 0 and l["count"] > 0 for l in s["leads"]) for s in scores))
    accs = [l["value"] for s in scores if s["metric"] == "acc" for l in s["leads"]]
    assert_("scorecard: every anomaly correlation lies in [-1, 1]", all(-1 <= x <= 1 for x in accs), str(accs[:5]))
    rm = [l["value"] for s in scores if s["metric"] in ("rmse", "track_km", "pressure_hpa") for l in s["leads"]]
    assert_("scorecard: errors are not negative", all(x >= 0 for x in rm))
else:
    skip("scorecard: not built (python data/products/verification.py)")
