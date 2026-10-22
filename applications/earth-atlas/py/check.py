"""Checks the plug-ins the server runs, without a page or a server: the Python
vorticity plug-in, the round trip of a product through the files the page reads,
and, when a C++ compiler is at hand, the C++ vorticity plug-in through the file
contract. A solid-body rotation u = U cos φ has ζ = 2 U sin φ / a and δ = 0.
Prints one JSON object; exits 1 when a check fails. The Earth Atlas's suite runs it.

    python applications/earth-atlas/py/check.py
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fields  # noqa: E402
import host  # noqa: E402
import plugins as api  # noqa: E402

U, R, DEG = 20.0, 6371008.8, math.pi / 180
GRID = {"lon0": -180.0, "dlon": 2.5, "nlon": 144, "lat0": 90.0, "dlat": 2.5, "nlat": 73, "periodic": True}
T0 = dt.datetime(2026, 10, 20, tzinfo=dt.timezone.utc)
TIMES = [T0, T0 + dt.timedelta(hours=48)]


def rotation() -> fields.Field:
    lats = GRID["lat0"] - np.arange(GRID["nlat"]) * GRID["dlat"]
    u = np.repeat((U * np.cos(lats * DEG))[:, None], GRID["nlon"], axis=1).astype(np.float32)
    v = np.zeros_like(u)
    return api.make_product("test-zonal", "Solid-body rotation", GRID, TIMES, [
        {"name": "u", "long_name": "Eastward wind at 500 hPa", "standard_name": "eastward_wind", "units": "m s-1", "level": "500 hPa", "slices": [u, u]},
        {"name": "v", "long_name": "Northward wind at 500 hPa", "standard_name": "northward_wind", "units": "m s-1", "level": "500 hPa", "slices": [v, v]},
    ])


def vorticity_errors(product: fields.Field) -> dict:
    row = round((GRID["lat0"] - 45) / GRID["dlat"])
    zeta, delta = product.slice("zeta", 0), product.slice("delta", 0)
    return {"zeta_error": abs(float(zeta[row, 10]) - 2 * U * math.sin(45 * DEG) / R), "delta": abs(float(delta[row, 10]))}


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

    # The round trip: written as the page reads it, read back within half a quantum.
    work = Path(tempfile.mkdtemp(prefix="check-"))
    try:
        path = fields.write(result["product"], work)
        back = fields.read(path)
        q = back.variables["zeta"]["scale"] / 2
        err = float(np.nanmax(np.abs(back.slice("zeta", 0) - result["product"].slice("zeta", 0))))
        out["round_trip"] = {"max_error": err, "half_quantum": q, "ok": err <= q * 1.001}
        ok &= out["round_trip"]["ok"]
    finally:
        shutil.rmtree(work, ignore_errors=True)

    # The C++ plug-in, compiled and run through the file contract, when a compiler is at hand.
    if shutil.which("c++") or shutil.which("g++") or shutil.which("clang++"):
        try:
            cpp = host.ProcessPlugin(host.plugin_path("plugins/vorticity.cpp"))
            [result] = [r for r in cpp.run({"wind": wind}, {}, context) if r["kind"] == "field"]
            e = vorticity_errors(result["product"])
            out["cpp"] = {**e, "ok": e["zeta_error"] < 2e-8 and e["delta"] < 1e-9}
        except Exception as error:   # reported, and counted as a failure
            out["cpp"] = {"ok": False, "error": f"{type(error).__name__}: {error}"[-600:]}
        ok &= out["cpp"]["ok"]
    else:
        out["cpp"] = {"skipped": "no C++ compiler"}

    print(json.dumps(out))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
