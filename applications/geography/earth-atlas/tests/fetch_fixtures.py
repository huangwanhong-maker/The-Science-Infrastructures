"""Fetch reference subsolar points from the JPL Horizons API for the Earth Atlas's tests.

The subsolar point at a UTC instant is the place where the Sun stands at the
zenith: its latitude is the Sun's apparent geocentric declination, and its east
longitude is the Sun's apparent right ascension minus Greenwich apparent
sidereal time. Horizons gives the first two for a geocentric observer, and the
third as the local apparent sidereal time of an observer at longitude 0. (Its
"sub-solar point" quantity, asked of an observer at the Sun, is retarded by the
light-time and so differs by about 2 degrees of longitude.)

Writes tests/fixtures/subsolar.json.

Usage::

    python applications/geography/earth-atlas/tests/fetch_fixtures.py
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://ssd.jpl.nasa.gov/api/horizons.api"
OUT = Path(__file__).resolve().parent / "fixtures" / "subsolar.json"

# JD (UT) from 1955 to 2045: solstices, equinoxes and ordinary days, at several hours.
EPOCHS = [2435109.5, 2440000.25, 2444239.75, 2447892.5, 2451545.0, 2451623.8, 2451716.3, 2453000.1,
          2455000.6, 2457204.4, 2459000.9, 2460000.0, 2460390.25, 2460480.5, 2461000.3, 2461335.5,
          2462000.7, 2464000.2, 2467000.85, 2469800.45]


def query(params: dict) -> list[list[str]]:
    fixed = {"format": "json", "MAKE_EPHEM": "YES", "OBJ_DATA": "NO", "CSV_FORMAT": "YES",
             "EPHEM_TYPE": "OBSERVER", "TLIST_TYPE": "JD", "TIME_TYPE": "UT", "ANG_FORMAT": "DEG",
             "EXTRA_PREC": "YES", "COMMAND": "'10'", "TLIST": " ".join(f"'{t}'" for t in EPOCHS)}
    url = API + "?" + urllib.parse.urlencode({**fixed, **params})
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = json.load(response)
    if "error" in payload:
        raise RuntimeError(payload["error"])
    body = payload["result"].split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    return [[c.strip() for c in line.split(",")] for line in body.strip().splitlines()]


def main() -> None:
    radec = query({"CENTER": "'500@399'", "QUANTITIES": "'2'"})
    sidereal = query({"CENTER": "'coord@399'", "COORD_TYPE": "GEODETIC", "SITE_COORD": "'0,0,0'",
                      "QUANTITIES": "'7'"})
    points = []
    for jd, a, s in zip(EPOCHS, radec, sidereal):
        ra, dec = float(a[3]), float(a[4])
        gast = float(s[3]) * 15
        lon = (ra - gast + 540) % 360 - 180
        points.append({"jd_ut": jd, "utc": a[0], "ra_deg": ra, "dec_deg": dec, "gast_deg": gast,
                       "lon_deg": round(lon, 6), "lat_deg": round(dec, 6)})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "source": "JPL Horizons API (DE441): Sun, apparent RA and Dec from the geocentre; apparent "
                  "sidereal time at longitude 0. Subsolar east longitude = RA - GAST, latitude = Dec.",
        "points": points,
    }, indent=1) + "\n", encoding="utf-8")
    print(f"{len(points)} subsolar points written to {OUT}")


if __name__ == "__main__":
    main()
