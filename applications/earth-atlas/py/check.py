"""Checks the plug-ins the server runs, without a page or a server, each against
an answer known in closed form, and prints one JSON object (the Earth Atlas's
suite runs this and reads it); exits 1 when a check fails.

  python            the Python vorticity plug-in: a solid-body rotation u = U cos φ
                    has ζ = 2 U sin φ / a and δ = 0;
  round_trip        a product written as the page reads it, read back within half a quantum;
  observations      a Python DataAssimilation with a need of kind 'observations' and an
                    optional need, resolved from a library on disk: the set is passed, the
                    optional need left without a source, the background interpolated in
                    time to the instant, the observation drawn;

and, when a C++ compiler is at hand (c++, g++ or clang++, or g++ in WSL on
Windows), the C++ plug-ins through the file contract of api/cpp/plugin.hpp:

  cpp               vorticity.cpp: the same solid-body rotation;
  cpp_advection     advection.cpp, a ForecastModel: the cosine bell of Williamson et
                    al. (1992), test case 1, carried over the poles by a solid-body
                    rotation, returns to itself after one revolution of 12 days; a
                    uniform field stays uniform; the verification is the host's;
  cpp_parcels       parcels.cpp, Trajectories: a parcel in u = U cos φ keeps its latitude
                    and turns by U t / a;
  cpp_cressman      cressman.cpp, DataAssimilation, from a library on disk with an
                    observation set: one observation on a grid point is drawn exactly,
                    the background kept beyond the radius, and the reports outside the
                    window ignored.

    python applications/earth-atlas/py/check.py
"""

from __future__ import annotations

import datetime as dt
import json
import math
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "api" / "python"))
import host  # noqa: E402
from science import fields, hosting  # noqa: E402
from science import plugins as api  # noqa: E402

U, R, DEG = 20.0, 6371008.8, math.pi / 180
GRID = {"lon0": -180.0, "dlon": 2.5, "nlon": 144, "lat0": 90.0, "dlat": 2.5, "nlat": 73, "periodic": True}
T0 = dt.datetime(2026, 10, 20, tzinfo=dt.timezone.utc)
TIMES = [T0, T0 + dt.timedelta(hours=48)]
HOUR = dt.timedelta(hours=1)
LONS = GRID["lon0"] + np.arange(GRID["nlon"]) * GRID["dlon"]
LATS = GRID["lat0"] - np.arange(GRID["nlat"]) * GRID["dlat"]
LON, LAT = np.meshgrid(LONS, LATS)


def cell(lon: float, lat: float) -> tuple[int, int]:
    return round((GRID["lat0"] - lat) / GRID["dlat"]), round((lon - GRID["lon0"]) / GRID["dlon"]) % GRID["nlon"]


def rotation() -> fields.Field:
    u = (U * np.cos(LAT * DEG)).astype(np.float32)
    v = np.zeros_like(u)
    return api.make_product("test-zonal", "Solid-body rotation", GRID, TIMES, [
        {"name": "u", "long_name": "Eastward wind at 500 hPa", "standard_name": "eastward_wind", "units": "m s-1", "level": "500 hPa", "slices": [u, u]},
        {"name": "v", "long_name": "Northward wind at 500 hPa", "standard_name": "northward_wind", "units": "m s-1", "level": "500 hPa", "slices": [v, v]},
    ])


def vorticity_errors(product: fields.Field) -> dict:
    row = round((GRID["lat0"] - 45) / GRID["dlat"])
    zeta, delta = product.slice("zeta", 0), product.slice("delta", 0)
    return {"zeta_error": abs(float(zeta[row, 10]) - 2 * U * math.sin(45 * DEG) / R), "delta": abs(float(delta[row, 10]))}


def given(product: fields.Field, names: list[str], level: str | None = None) -> api.Input:
    return api.Input(product, names, [product.variables[n] for n in names], product.meta, product.meta["title"], level)


# ---- the library on disk: a background at two times, and an observation set --------------------


B_TIMES = [T0, T0 + 6 * HOUR]
AT = T0 + 3 * HOUR                     # the instant shown, between the background's times
STATION = (0.0, 40.0)                  # on a grid point


def background_at(when: dt.datetime) -> np.ndarray:
    """280 K + 5 K cos φ cos λ, warming by 1 K an hour."""
    return (280 + (when - T0) / HOUR + 5 * np.cos(LAT * DEG) * np.cos(LON * DEG)).astype(np.float32)


def library(base: Path) -> hosting.Library:
    products = base / "products"
    t2m = api.make_product("test-t2m", "Test background", GRID, B_TIMES, [
        {"name": "t2m", "long_name": "Temperature at 2 m", "standard_name": "air_temperature", "units": "K", "level": "2 m",
         "slices": [background_at(t) for t in B_TIMES]}], model="Test model")
    path = fields.write(t2m, products)
    b = float(background_at(AT)[cell(*STATION)])
    rows = [["S1", "Station on a grid point", STATION[0], STATION[1], 10.0, fields.iso(AT), b + 3.0],
            ["S2", "Station reporting too early", 60.0, -30.0, 10.0, fields.iso(AT - 5 * HOUR), 400.0]]
    observations = {"schema": "observations/1", "name": "test-stations", "title": "Test stations",
                    "fields": ["id", "name", "lon", "lat", "elevation", "time", "air_temperature"], "rows": rows}
    (products / "observations-test-stations.json").write_text(json.dumps(observations), encoding="utf-8")
    index = {"fields": [{"name": "test-t2m", "file": path.name}],
             "observations": [{"name": "test-stations", "file": "observations-test-stations.json", "title": "Test stations",
                               "variables": ["air_temperature"]}]}
    (products / "index.json").write_text(json.dumps(index), encoding="utf-8")
    return hosting.Library(products, base / "plugins", host.ROOT)


class Drawn(api.DataAssimilation):
    """A Python analysis that sets the cell of each observation to its value."""
    id = "check-observations"
    title = "A check of observation needs"

    def needs(self):
        return [{"name": "background", "label": "Background", "standard_name": "air_temperature", "level": "2 m"},
                {"name": "stations", "label": "Observed by", "kind": "observations", "variable": "air_temperature"},
                {"name": "truth", "label": "Or by a model", "standard_name": "air_temperature", "level": "2 m", "distinct": True, "optional": True}]

    def observations(self, inputs, when, params, context):
        set_, variable = inputs["stations"].observations, inputs["stations"].variable
        ix = {f: i for i, f in enumerate(set_["fields"])}
        self.seen = {"set": inputs["stations"].name, "truth": inputs["truth"], "when": when}
        return [{"lon": r[ix["lon"]], "lat": r[ix["lat"]], "value": r[ix[variable]]} for r in set_["rows"]
                if abs(fields.parse_time(r[ix["time"]]) - when) <= 1.5 * HOUR]

    def analyse(self, background, grid, observations, params, context):
        out = background.copy()
        for o in observations:
            out[cell(o["lon"], o["lat"])] = o["value"]
        return out


def observations_check(app_host: hosting.Host, lib: hosting.Library) -> dict:
    plugin = Drawn()
    needs = {n["name"]: n for n in hosting.resolve(plugin, lib)}
    served = [c["key"] for c in needs["stations"]["candidates"]]
    optional = [c["key"] for c in needs["truth"]["candidates"]]
    context = api.Context(time=AT)
    results = app_host.run_plugin(plugin, {"background": "test-t2m|2 m", "stations": "test-stations", "truth": ""}, {}, context, lib)
    product = fields.read(app_host.out / plugin.id / "fields-plugin-check-observations.json")
    r, c = cell(*STATION)
    far = cell(-120, -40)
    b = background_at(AT)
    analysed, increment = product.slice("t2m", 0), product.slice("t2m_increment", 0)
    q = product.variables["t2m"]["scale"] + product.variables["t2m_increment"]["scale"]
    out = {"served_by": served, "optional": optional, "time": fields.iso(product.times[0]), "truth_given": plugin.seen["truth"] is not None,
           "drawn_error": abs(float(analysed[r, c]) - float(b[r, c] + 3)),
           "background_error": abs(float(analysed[far]) - float(b[far])), "far_increment": abs(float(increment[far])),
           "quantum": q, "markers": len(next(x for x in results if x["kind"] == "sketch")["markers"])}
    out["ok"] = (served == ["test-stations"] and optional[-1] == "" and not out["truth_given"] and out["time"] == fields.iso(AT)
                 and out["drawn_error"] <= q and out["background_error"] <= q and out["far_increment"] <= q and out["markers"] == 1)
    return out


# ---- the C++ plug-ins ------------------------------------------------------------------------------


# Advection is checked on a 2° grid. On the 2.5° grid the bell moves exactly half a cell
# in each step of an hour (u0 Δt / a = 2π/288 = 1.25°), interpolating every time where
# cubic interpolation damps most, and l2 is 0.13; on 1.25° it moves a whole cell and the
# test says little of the interpolation.
GRID2 = {"lon0": -180.0, "dlon": 2.0, "nlon": 180, "lat0": 90.0, "dlat": 2.0, "nlat": 91, "periodic": True}
LON2, LAT2 = np.meshgrid(GRID2["lon0"] + np.arange(GRID2["nlon"]) * 2.0, GRID2["lat0"] - np.arange(GRID2["nlat"]) * 2.0)


def williamson(alpha: float):
    """Test case 1 of Williamson et al. (1992), on GRID2: u = u0 (cos φ cos α + sin φ cos λ sin α),
    v = −u0 sin λ sin α, u0 = 2πa / 12 days; and the cosine bell h = (h0/2)(1 + cos(π r / R))
    within r < R = a/3 of (3π/2, 0), here a temperature of 250 K + 50 K h/h0."""
    u0 = 2 * math.pi * R / (12 * 86400)
    lam, phi = LON2 * DEG, LAT2 * DEG
    u = (u0 * (np.cos(phi) * math.cos(alpha) + np.sin(phi) * np.cos(lam) * math.sin(alpha))).astype(np.float32)
    v = (-u0 * np.sin(lam) * math.sin(alpha)).astype(np.float32)
    lc, pc = 3 * math.pi / 2, 0.0
    r = R * np.arccos(np.clip(math.sin(pc) * np.sin(phi) + math.cos(pc) * np.cos(phi) * np.cos(lam - lc), -1, 1))
    bell = np.where(r < R / 3, 0.5 * (1 + np.cos(3 * math.pi * r / R)), 0.0)
    return u, v, (250 + 50 * bell).astype(np.float32)


def weighted(x: np.ndarray) -> float:
    """The root mean square over the sphere (GRID2), weighted by the cosine of latitude."""
    w = np.cos(LAT2 * DEG)
    return float(np.sqrt(np.sum(w * x.astype(np.float64) ** 2) / np.sum(w)))


def advection_check(binaries: Path, work: Path) -> dict:
    cpp = hosting.ProcessPlugin(host.plugin_path("plugins/advection.cpp"), work, binaries)
    out = {}
    days12 = [T0, T0 + 288 * HOUR]
    u, v, t = williamson(math.pi / 2 - 0.05)
    wind = api.make_product("test-rotation", "Solid-body rotation over the poles", GRID2, days12, [
        {"name": "u", "long_name": "Eastward wind at 10 m", "standard_name": "eastward_wind", "units": "m s-1", "level": "10 m", "slices": [u, u]},
        {"name": "v", "long_name": "Northward wind at 10 m", "standard_name": "northward_wind", "units": "m s-1", "level": "10 m", "slices": [v, v]}])
    for mode in ("cubic", "clipped"):
        bell = api.make_product("test-bell", "Cosine bell", GRID2, days12, [
            {"name": "t", "long_name": "Temperature at 2 m", "standard_name": "air_temperature", "units": "K", "level": "2 m", "slices": [t, t]}])
        start = time.perf_counter()
        results = cpp.run({"scalar": given(bell, ["t"], "2 m"), "wind": given(wind, ["u", "v"], "10 m")},
                          {"hours": 288, "interpolation": mode}, api.Context(time=T0))
        seconds = time.perf_counter() - start
        product = next(r["product"] for r in results if r["kind"] == "field")
        end = product.slice("t", len(product.times) - 1).astype(np.float64)
        error = end - t
        anomaly = t.astype(np.float64) - 250
        leads = product.meta.get("verification", {}).get("leads", [])
        own = api.rms_difference(end, t, GRID2, 20, 80)
        out[mode] = {"l2": weighted(error) / weighted(anomaly), "linf": float(np.max(np.abs(error)) / np.max(anomaly)),
                     "least": float(np.min(end) - 250), "times": len(product.times), "seconds": round(seconds, 2),
                     "verified_at": [r["hours"] for r in leads], "verification_error": abs(leads[-1]["model"] - own) if leads else None,
                     "persistence": leads[-1]["persistence"] if leads else None, "table": any(r["kind"] == "table" for r in results)}
    uniform = np.full_like(t, 280.0)
    flat = api.make_product("test-uniform", "Uniform", GRID2, days12, [
        {"name": "t", "long_name": "Temperature at 2 m", "standard_name": "air_temperature", "units": "K", "level": "2 m", "slices": [uniform, uniform]}])
    results = cpp.run({"scalar": given(flat, ["t"], "2 m"), "wind": given(wind, ["u", "v"], "10 m")}, {"hours": 48}, api.Context(time=T0))
    product = next(r["product"] for r in results if r["kind"] == "field")
    out["uniform_error"] = max(float(np.max(np.abs(product.slice("t", k) - 280.0))) for k in range(len(product.times)))
    c, m = out["cubic"], out["clipped"]
    # Tolerances: on 2°, cubic gives l2 0.065 and l∞ 0.047, clipped l2 0.082 and no undershoot (bilinear: l2 0.65).
    out["ok"] = (c["l2"] < 0.08 and c["linf"] < 0.06 and m["l2"] < 0.10 and m["least"] >= -1e-3 and c["times"] == 49
                 and c["verified_at"] == [288.0] and c["verification_error"] < 1e-6 and c["persistence"] == 0 and c["table"]
                 and out["uniform_error"] < 1e-4)
    return out


def parcels_check(binaries: Path, work: Path, zonal: fields.Field) -> dict:
    cpp = hosting.ProcessPlugin(host.plugin_path("plugins/parcels.cpp"), work, binaries)
    results = cpp.run({"wind": given(zonal, ["u", "v"], "500 hPa")}, {"hours": 24, "step": 1, "count": 1, "radius": 0},
                      api.Context(time=T0, pinned=[10, 40]))
    [sketch] = [r for r in results if r["kind"] == "sketch"]
    path = sketch["polylines"][0]["coords"]
    out = {"latitude_error": abs(path[-1] - 40), "longitude_error": abs(path[-2] - 10 - U * 86400 / R / DEG),
           "steps": len(path) // 2 - 1, "moving_starts": sketch["moving"][0]["times"][0] == int(T0.timestamp() * 1000),
           "label": [x["text"] for x in sketch["labels"]], "notes": sum(r["kind"] == "note" for r in results)}
    out["ok"] = (out["latitude_error"] < 1e-6 and out["longitude_error"] < 0.01 and out["steps"] == 24 and out["moving_starts"]
                 and out["label"] == ["start"] and out["notes"] == 1)
    return out


def cressman_check(app_host: hosting.Host, lib: hosting.Library) -> dict:
    plugin, _ = app_host.load_plugin("plugins/cressman.cpp")
    results = app_host.run_plugin(plugin, {"background": "test-t2m|2 m", "stations": "test-stations"},
                                  {"radius": 800, "passes": 3, "epsilon": 0, "window": 1.5}, api.Context(time=AT), lib)
    product = fields.read(app_host.out / plugin.id / "fields-plugin-cressman-cpp.json")
    b = background_at(AT)
    at, inside, beyond = cell(*STATION), cell(5, 40), cell(40, 40)      # 0, 426 and 3 380 km from the station
    analysed, increment = product.slice("t2m", 0), product.slice("t2m_increment", 0)
    q = product.variables["t2m"]["scale"] + product.variables["t2m_increment"]["scale"]
    [table] = [r for r in results if r["kind"] == "table"]
    [sketch] = [r for r in results if r["kind"] == "sketch"]
    out = {"drawn_error": abs(float(analysed[at]) - float(b[at] + 3)), "inside_increment_error": abs(float(increment[inside]) - 3),
           "beyond_increment": abs(float(increment[beyond])), "beyond_background_error": abs(float(analysed[beyond]) - float(b[beyond])),
           "quantum": q, "time": fields.iso(product.times[0]), "used": table["rows"][0][1], "markers": len(sketch["markers"]),
           "increment_is_difference": bool(product.variables["t2m_increment"].get("difference"))}
    out["ok"] = (out["drawn_error"] <= 2 * q and out["inside_increment_error"] <= 2 * q and out["beyond_increment"] <= q
                 and out["beyond_background_error"] <= 2 * q and out["time"] == fields.iso(AT) and out["used"] == "1"
                 and out["markers"] == 1 and out["increment_is_difference"])
    return out


def main() -> int:
    out, ok = {}, True
    zonal = rotation()
    wind = api.Input(zonal, ["u", "v"], zonal.meta["variables"], zonal.meta, "test", "500 hPa")
    context = api.Context(time=T0)

    # The Python plug-in, in this process.
    plugin, _ = host.load_plugin("plugins/vorticity.py")
    [result] = plugin.run({"wind": wind}, {}, context)
    e = vorticity_errors(result["product"])
    out["python"] = {**e, "ok": e["zeta_error"] < 2e-8 and e["delta"] < 1e-9}
    ok &= out["python"]["ok"]

    work = Path(tempfile.mkdtemp(prefix="check-"))
    try:
        # The round trip: written as the page reads it, read back within half a quantum.
        path = fields.write(result["product"], work / "round-trip")
        back = fields.read(path)
        q = back.variables["zeta"]["scale"] / 2
        err = float(np.nanmax(np.abs(back.slice("zeta", 0) - result["product"].slice("zeta", 0))))
        out["round_trip"] = {"max_error": err, "half_quantum": q, "ok": err <= q * 1.001}
        ok &= out["round_trip"]["ok"]

        # Needs of kind 'observations', and optional needs, for a Python plug-in. The results go
        # to the scratch directory; the C++ programs are built where the atlas's server builds them.
        lib = library(work)
        app_host = hosting.Host(host.APP, work / "products", work / "plugins", binaries=host.OUT / "bin")
        try:
            out["observations"] = observations_check(app_host, lib)
        except Exception as error:   # reported, and counted as a failure
            out["observations"] = {"ok": False, "error": f"{type(error).__name__}: {error}"[-600:]}
        ok &= out["observations"]["ok"]

        # The C++ plug-ins, compiled and run through the file contract, when a compiler is at hand.
        cxx = hosting.compiler()
        names = ("cpp", "cpp_advection", "cpp_parcels", "cpp_cressman")
        if not cxx:
            for name in names:
                out[name] = {"skipped": "no C++ compiler (c++, g++ or clang++, nor g++ in WSL)"}
        else:
            out["compiler"] = repr(cxx)
            checks = {
                "cpp": lambda: vorticity_cpp(wind, context),
                "cpp_advection": lambda: advection_check(host.OUT / "bin", work / "plugins"),
                "cpp_parcels": lambda: parcels_check(host.OUT / "bin", work / "plugins", zonal),
                "cpp_cressman": lambda: cressman_check(app_host, lib),
            }
            for name, run in checks.items():
                try:
                    out[name] = run()
                except Exception as error:   # reported, and counted as a failure
                    out[name] = {"ok": False, "error": f"{type(error).__name__}: {error}"[-600:]}
                ok &= out[name]["ok"]
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print(json.dumps(out))
    return 0 if ok else 1


def vorticity_cpp(wind: api.Input, context: api.Context) -> dict:
    cpp = host.ProcessPlugin(host.plugin_path("plugins/vorticity.cpp"))
    [result] = [r for r in cpp.run({"wind": wind}, {}, context) if r["kind"] == "field"]
    e = vorticity_errors(result["product"])
    return {**e, "ok": e["zeta_error"] < 2e-8 and e["delta"] < 1e-9}


if __name__ == "__main__":
    sys.exit(main())
