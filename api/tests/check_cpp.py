"""The Python client's values for the checks of the C and C++ clients
(check_cpp.cpp, check_c.c), run by api/tests/run-tests.mjs: the same requests,
made by science.Client against the same service, written as one JSON object on
standard output, which the check programs read on their standard input and
hold their own answers to. The two clients read the same JSON from the same
service, so they must agree exactly: samples and series as doubles, grids as
float32.

    python api/tests/check_cpp.py http://127.0.0.1:<port>
"""

from __future__ import annotations

import datetime as dt
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "api" / "python"))
import science  # noqa: E402

UTC = dt.timezone.utc

# Places over the globe, the poles and both sides of the antimeridian, as the service's own checks take them.
LON = [139.69, -0.13, 179.9, -179.95, 0.0, 359.9, 12.3, -73.97, 151.2, 37.6]
LAT = [35.69, 51.51, -0.4, 10.2, 89.95, -89.99, -33.9, 40.78, -33.87, 55.75]


def plain(values) -> list:
    """Numbers as JSON can hold them: NaN as null."""
    return [None if (v is None or math.isnan(v)) else float(v) for v in values]


def iso(when: dt.datetime) -> str:
    return when.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def flat(grid) -> list:
    rows = grid.values.tolist() if hasattr(grid.values, "tolist") else grid.values
    return plain([v for row in rows for v in row])


def main(argv: list[str]) -> int:
    api = science.Client(argv[0])
    out: dict = {"version": api.info()["version"]}
    names = [f["name"] for f in api.fields()]
    out["fields"] = names

    if "gfs-forecast" in names:
        gfs = api.field("gfs-forecast")
        between = iso(gfs.times[2] + dt.timedelta(hours=2))
        out["sample"] = {"field": "gfs-forecast", "var": "tmp2m", "lon": LON, "lat": LAT, "time": between,
                         "values": plain(gfs.sample("tmp2m", LON, LAT, time=between))}
        # More places than a GET of the C++ client carries: it asks by POST.
        lon = [round(-179.5 + 0.599 * k, 3) for k in range(600)]
        lat = [round(-80.0 + 0.27 * k, 3) for k in range(600)]
        out["many"] = {"field": "gfs-forecast", "var": "rh2m", "lon": lon, "lat": lat, "time": between,
                       "values": plain(gfs.sample("rh2m", lon, lat, time=between))}
        s = gfs.series("tmp2m", 139.69, 35.69)
        out["series"] = {"field": "gfs-forecast", "var": "tmp2m", "lon": 139.69, "lat": 35.69, "units": s.units,
                         "times": [iso(t) for t in s.times], "times_unix": [t.timestamp() for t in s.times],
                         "values": plain(s.values)}
        g = gfs.grid("tmp2m", time=3, bbox=(170, -10, -170, 10), stride=2)
        out["grid"] = {"field": "gfs-forecast", "var": "tmp2m", "time": "3", "bbox": [170, -10, -170, 10], "stride": 2,
                       "grid": g.grid, "at": g.time, "values": flat(g)}
        b = gfs.grid("tmp2m", time=between, bbox=(122, 24, 150, 46))
        out["blend"] = {"field": "gfs-forecast", "var": "tmp2m", "time": between, "bbox": [122, 24, 150, 46], "stride": 1,
                        "grid": b.grid, "at": b.time, "values": flat(b)}

    # The domains' functions, through the bridge, the Sun's chain as the examples take it.
    ms = 1792735200000                     # 2026-10-23T06:00:00Z
    cases = [("astronomy", "sunPosition", [2460000.5]), ("astronomy", "julianDate", [ms]),
             ("astronomy", "subsolarPoint", [ms]), ("meteorology", "windFromDirection", [5, 0]),
             ("meteorology", "windFromDirection", [-3, 4]), ("geodesy", "distanceKm", [35.69, 139.69, 51.51, -0.13]),
             ("physis", "bilinear", [{"$typed": "Float32Array", "values": [1, 2, 3, None]},
                                     {"lon0": 0, "dlon": 1, "nlon": 2, "lat0": 1, "dlat": 1, "nlat": 2, "periodic": False}, 0.25, 0.75])]
    computed = []
    try:
        for unit, function, args in cases:
            computed.append({"unit": unit, "function": function, "args": args, "result": api.compute(unit, function, *args)})
        sub = computed[2]["result"]
        computed.append({"unit": "astronomy", "function": "sunAltAz", "args": [139.69, 35.69, sub],
                         "result": api.compute("astronomy", "sunAltAz", 139.69, 35.69, sub)})
        out["compute"] = computed
    except science.ApiError as error:
        out["compute_error"] = str(error).splitlines()[0][:300]       # the first line: a worker's error may carry its stack

    events = [e["name"] for e in api.events()]
    if events:
        out["events"] = {"name": events[0], "limit": 5, "events": api.events(events[0], limit=5)["events"]}
    out["verification"] = [v["name"] for v in api.verification()]
    sys.stdout.buffer.write(json.dumps(out, ensure_ascii=False, allow_nan=False).encode("utf-8") + b"\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
