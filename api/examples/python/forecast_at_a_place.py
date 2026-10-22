"""The forecast at a place over the coming days: the temperature at 2 m and the
wind at 10 m from GFS, every six hours, and the temperature tomorrow from every
model that forecasts it. Sampling is done by the service, as the Earth Atlas does
it (bilinear in space, linear in time).

    python api/examples/python/forecast_at_a_place.py [lon lat] [--url http://127.0.0.1:8800]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))    # the package, when not installed
import science  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("lon", type=float, nargs="?", default=139.69, help="degrees east (default: Tokyo)")
parser.add_argument("lat", type=float, nargs="?", default=35.69, help="degrees north")
parser.add_argument("--url", help="the service (default: $SCIENCE_API_URL or http://127.0.0.1:8800)")
args = parser.parse_args()

api = science.Client(args.url)
gfs = api.field("gfs-forecast")
temperature = gfs.series("tmp2m", args.lon, args.lat)
u = gfs.series("ugrd10m", args.lon, args.lat).values
v = gfs.series("vgrd10m", args.lon, args.lat).values

print(f"GFS at {args.lon:.2f}, {args.lat:.2f} (UTC)")
for when, t, ue, vn in zip(temperature.times, temperature.values, u, v):
    speed = math.hypot(ue, vn)
    # The direction the wind comes from, by the meteorology domain's own function.
    frm = api.compute("meteorology", "windFromDirection", float(ue), float(vn))
    print(f"  {when:%a %d %H:%M}   {t - 273.15:5.1f} °C   wind {speed:4.1f} m/s from {frm:3.0f}°")

print("\nTomorrow at this hour, 2 m temperature by model:")
for entry in api.fields():
    for variable in entry.get("variables", []):
        if variable.get("standard_name") == "air_temperature" and variable.get("level") == "2 m" and "stdDev" not in variable["name"]:
            value = api.field(entry["name"]).sample(variable["name"], args.lon, args.lat, time="now+1d")
            if not math.isnan(value):
                print(f"  {entry.get('model') or entry['name']:28s} {value - 273.15:5.1f} °C")
