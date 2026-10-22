"""A sensor that posts its readings to the service, as a home weather station or a
microcontroller would: it defines its source once (the variables and the units it
measures in; the store converts them to SI), then posts a reading every few
seconds. Run sensor_stream.py in a second terminal to watch them arrive.

The readings here are made up (a daily cycle and some noise), so that the example
runs anywhere; a real device posts what it measures, the same way.

    python api/examples/python/sensor_post.py [--source example-garden] [--every 5] [--url ...] [--token ...]
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import science  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--source", default="example-garden")
parser.add_argument("--every", type=float, default=5.0, help="seconds between readings")
parser.add_argument("--count", type=int, default=0, help="stop after this many (0: never)")
parser.add_argument("--url")
parser.add_argument("--token", help="the service's token, when it has one (default: $SCIENCE_API_TOKEN)")
args = parser.parse_args()

api = science.Client(args.url, token=args.token)
define = {"title": "Example garden station", "variables": ["air_temperature:degree_C", "relative_humidity:%"],
          "stations": [{"id": "garden-1", "name": "The garden", "lon": 139.69, "lat": 35.69, "elevation": 40}]}
sent = 0
while not args.count or sent < args.count:
    now = dt.datetime.now(dt.timezone.utc)
    hour = now.hour + now.minute / 60 + 9                          # local solar-ish time in Tokyo
    t = 17 + 5 * math.sin((hour - 9) / 24 * 2 * math.pi) + random.gauss(0, 0.2)
    rh = 65 - 15 * math.sin((hour - 9) / 24 * 2 * math.pi) + random.gauss(0, 1)
    reading = {"station": "garden-1", "lon": 139.69, "lat": 35.69, "time": now,
               "values": {"air_temperature": round(t, 2), "relative_humidity": round(rh, 1)}}
    answer = api.post_observations(args.source, [reading], define=define if sent == 0 else None)
    print(f"{now:%H:%M:%S}  {t:5.2f} °C  {rh:4.1f} %   accepted {answer['accepted']}"
          + (f", refused: {answer['rejected']}" if answer["rejected"] else ""), flush=True)
    sent += 1
    time.sleep(args.every)
