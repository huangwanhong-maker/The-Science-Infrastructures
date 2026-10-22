"""A simulated weather sensor: the template a real device copies.

Every --every seconds it makes one reading at its place (--lon, --lat,
--elevation) and sends it: by HTTP POST to the application's API (--url, with
a Bearer token if --token is given), or, without a URL, straight into the
store (store.append), building the product after each. A reading the API does
not take is kept and sent again with the next, as a device without a network
should (at most a day of readings is kept).

The weather it reports is made up, but has the shape of the real thing:

- air temperature: a mean that falls with latitude and with height (the ICAO
  standard atmosphere's 6.5 K per km; ICAO Doc 7488, 1993) and follows the
  seasons, a diurnal cycle of 4 K amplitude with its maximum at 15:00 local
  solar time, a slow wandering of the weather (a first-order autoregressive
  process, 1.5 K, about a day's memory) and the sensor's own noise (0.1 K);
- relative humidity: from a dew point that wanders slowly below the mean
  temperature, so that it falls by day and rises by night, as it does;
- pressure at sea level: 1013.25 hPa, a wandering of the weather (8 hPa, two
  days' memory), and the atmosphere's solar semidiurnal tide, S2, of amplitude
  1.16 hPa cos^3(latitude), with maxima near 10:00 and 22:00 local solar time
  (Haurwitz 1956; Chapman and Lindzen 1970, chapter 1).

Local solar time is UTC plus the longitude over 15 degrees an hour plus the
equation of time, by Spencer's (1971) Fourier series as given by Duffie and
Beckman (2013, eq. 1.5.3), good to about half a minute.

The body of a POST is {"title", "variables", "readings": [reading, ...]}: the
description travels with the readings, so that a server can define a new
source on its first reading. A reading is {"station", "lon", "lat",
"elevation", "time" (ISO 8601, UTC), "values": {variable: value in SI}}.

Usage::

    python observation/sensors/simulate.py garden --lon 139.70 --lat 35.68
    python observation/sensors/simulate.py garden --lon 139.70 --lat 35.68 --every 10 --count 30
    python observation/sensors/simulate.py garden --lon 139.70 --lat 35.68 \\
        --url http://127.0.0.1:8800/api/v1/observations/garden --token $TOKEN
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import random
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import store  # noqa: E402
from observation_set import relative_humidity  # noqa: E402  (data/products, on store's path)

UTC = dt.timezone.utc
VARIABLES = [
    {"name": "air_temperature", "units": "K"},
    {"name": "relative_humidity", "units": "%"},
    {"name": "air_pressure_at_sea_level", "units": "Pa"},
]
BACKLOG = 1440                  # readings kept while the API cannot be reached


def equation_of_time(t: dt.datetime) -> float:
    """Apparent minus mean solar time, in minutes (Spencer 1971; Duffie and Beckman 2013, eq. 1.5.3)."""
    b = 2.0 * math.pi * (t.timetuple().tm_yday - 1) / 365.0
    return 229.2 * (0.000075 + 0.001868 * math.cos(b) - 0.032077 * math.sin(b)
                    - 0.014615 * math.cos(2 * b) - 0.04089 * math.sin(2 * b))


def solar_hour(t: dt.datetime, lon: float) -> float:
    """Local apparent solar time, in hours (0 to 24)."""
    hours = t.hour + t.minute / 60.0 + t.second / 3600.0
    return (hours + lon / 15.0 + equation_of_time(t) / 60.0) % 24.0


class Wander:
    """A first-order autoregressive process of standard deviation sigma and
    memory tau (seconds), stepped by the time elapsed."""

    def __init__(self, sigma: float, tau: float, rng: random.Random):
        self.sigma, self.tau, self.rng, self.x = sigma, tau, rng, rng.gauss(0.0, sigma)

    def step(self, seconds: float) -> float:
        a = math.exp(-seconds / self.tau)
        self.x = a * self.x + math.sqrt(1.0 - a * a) * self.sigma * self.rng.gauss(0.0, 1.0)
        return self.x


class Sensor:
    def __init__(self, lon: float, lat: float, elevation: float, seed: int | None = None):
        self.lon, self.lat, self.elevation = lon, lat, elevation
        self.rng = random.Random(seed)
        self.weather = Wander(1.5, 86400.0, self.rng)            # K
        self.dryness = Wander(1.5, 2 * 86400.0, self.rng)        # K, below the mean temperature
        self.synoptic = Wander(800.0, 2 * 86400.0, self.rng)     # Pa
        self.last = None

    def mean_temperature(self, t: dt.datetime) -> float:
        """K: warm at the equator, cold at the poles and on heights; summer in July
        in the north and in January in the south."""
        phi = math.radians(self.lat)
        season = math.cos(2.0 * math.pi * (t.timetuple().tm_yday - 197) / 365.0)
        return 300.0 - 30.0 * math.sin(phi) ** 2 + 12.0 * math.sin(phi) * season - 0.0065 * self.elevation

    def read(self, t: dt.datetime) -> dict:
        seconds = 60.0 if self.last is None else max(1.0, (t - self.last).total_seconds())
        self.last = t
        h = solar_hour(t, self.lon)
        mean = self.mean_temperature(t) + self.weather.step(seconds)
        temperature = mean + 4.0 * math.cos(2.0 * math.pi * (h - 15.0) / 24.0) + self.rng.gauss(0.0, 0.1)
        dew_point = min(temperature, mean - 6.0 - abs(self.dryness.step(seconds)))
        tide = 116.0 * math.cos(math.radians(self.lat)) ** 3 * math.cos(4.0 * math.pi * (h - 10.0) / 24.0)
        pressure = 101325.0 + self.synoptic.step(seconds) + tide + self.rng.gauss(0.0, 5.0)
        return {"air_temperature": round(temperature, 2),
                "relative_humidity": round(min(100.0, relative_humidity(temperature, dew_point)), 1),
                "air_pressure_at_sea_level": round(pressure, 1)}


def post(url: str, token: str | None, body: dict) -> None:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(url, data=json.dumps(body).encode("utf-8"), headers=headers, method="POST")
    with urllib.request.urlopen(request, timeout=15) as response:
        response.read()


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", help="the source's name: lowercase, digits and hyphens")
    parser.add_argument("--lon", type=float, required=True, help="longitude, degrees east (WGS 84)")
    parser.add_argument("--lat", type=float, required=True, help="latitude, degrees north (WGS 84)")
    parser.add_argument("--elevation", type=float, default=0.0, help="height above mean sea level, m")
    parser.add_argument("--station", help="the station's id (default: <source>-1)")
    parser.add_argument("--title", help="the source's title (default: 'Simulated sensor, <source>')")
    parser.add_argument("--every", type=float, default=60.0, help="seconds between readings (default 60)")
    parser.add_argument("--count", type=int, default=0, help="readings to make, then stop (default: no end)")
    parser.add_argument("--url", help="the API to post to, as http://127.0.0.1:8800/api/v1/observations/<source>")
    parser.add_argument("--token", help="a Bearer token for the API")
    parser.add_argument("--seed", type=int, help="seed of the random weather, for a repeatable run")
    args = parser.parse_args(argv)
    try:
        store._folder(args.source)
    except ValueError as error:
        raise SystemExit(str(error)) from None
    station = args.station or f"{args.source}-1"
    title = args.title or f"Simulated sensor, {args.source}"
    sensor = Sensor(args.lon, args.lat, args.elevation, args.seed)
    if not args.url and not any(s["name"] == args.source for s in store.sources()):
        store.define(args.source, title, VARIABLES, [{"id": station, "name": title, "lon": args.lon, "lat": args.lat,
                                                     "elevation": args.elevation}],
                     credit=f"{title} (simulated)", licence="Simulated data, free to use.")
    backlog, made = [], 0
    try:
        while True:
            t = dt.datetime.now(UTC).replace(microsecond=0)
            values = sensor.read(t)
            reading = {"station": station, "lon": args.lon, "lat": args.lat, "elevation": args.elevation,
                       "time": t.strftime("%Y-%m-%dT%H:%M:%SZ"), "values": values}
            shown = f"{values['air_temperature'] - 273.15:5.1f} C, {values['relative_humidity']:5.1f} %, " \
                    f"{values['air_pressure_at_sea_level'] / 100:7.1f} hPa"
            if args.url:
                backlog = (backlog + [reading])[-BACKLOG:]
                try:
                    post(args.url, args.token, {"title": title, "variables": VARIABLES, "readings": backlog})
                    print(f"{reading['time']}  {shown}  posted" + (f" with {len(backlog) - 1} kept" if len(backlog) > 1 else ""))
                    backlog = []
                except urllib.error.HTTPError as error:
                    detail = error.read().decode("utf-8", errors="replace")[:300]
                    if error.code in (400, 404, 409, 413, 422):
                        # Refused as it is: sending it again would not help.
                        print(f"{reading['time']}  {shown}  refused by the API ({error.code}): {detail}")
                        backlog = []
                    else:
                        print(f"{reading['time']}  {shown}  not posted ({error.code} {detail}); {len(backlog)} kept")
                except (urllib.error.URLError, OSError) as error:
                    print(f"{reading['time']}  {shown}  not posted ({error}); {len(backlog)} kept")
            else:
                rejected = []
                store.append(args.source, [reading], rejected)
                store.build_product(args.source)
                print(f"{reading['time']}  {shown}" + (f"  refused: {rejected[0][1]}" if rejected else ""))
            made += 1
            if args.count and made >= args.count:
                return 0
            time.sleep(args.every)
    except KeyboardInterrupt:
        return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
