"""The domains' algorithms from Python, through the compute bridge: where the Sun
stands over a place now (astronomy: the subsolar point, then the altitude and
azimuth), the distance between two places (geodesy), and the direction a wind
comes from (meteorology). The functions run as the pages run them; the service
lists every one at /api/v1/compute.

    python api/examples/python/sun_and_wind.py [lon lat] [--url ...]
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import science  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("lon", type=float, nargs="?", default=139.69)
parser.add_argument("lat", type=float, nargs="?", default=35.69)
parser.add_argument("--url")
args = parser.parse_args()

api = science.Client(args.url)
now = dt.datetime.now(dt.timezone.utc)
ms = now.timestamp() * 1000
sub = api.compute("astronomy", "subsolarPoint", ms)                 # {lon, lat}, degrees
sky = api.compute("astronomy", "sunAltAz", args.lon, args.lat, sub)  # {altitude, azimuth}, degrees
light = api.compute("astronomy", "daylight", sky["altitude"])
print(f"{now:%Y-%m-%d %H:%M} UTC: the Sun is overhead at {sub['lon']:.2f}, {sub['lat']:.2f};")
print(f"  at {args.lon:.2f}, {args.lat:.2f} it stands {sky['altitude']:.1f}° high, azimuth {sky['azimuth']:.1f}° ({light}).")

km = api.compute("geodesy", "distanceKm", args.lat, args.lon, 51.4779, -0.0015)
print(f"  Greenwich is {km:,.0f} km away (on the sphere of mean radius).")

for u, v in [(5, 0), (0, -5), (-3, 4)]:
    print(f"  a wind of u = {u}, v = {v} m/s comes from {api.compute('meteorology', 'windFromDirection', u, v):.0f}°")

sun = next(f for f in api.units("astronomy")["functions"] if f["name"] == "sunPosition")
print(f"\nastronomy.sunPosition({', '.join(sun['parameters'])}): {sun['doc']}")
