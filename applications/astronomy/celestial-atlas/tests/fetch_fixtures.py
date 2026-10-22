"""Fetch reference values from the JPL Horizons API for the atlas's tests.

Writes tests/fixtures/horizons.json with:
  - barycentric state vectors (ICRF, km, km/s) of the Sun, the planets and the Moon
    at TDB epochs spanning the DE440s coverage;
  - geocentric astrometric right ascension and declination (ICRF, degrees) of
    several bodies at UTC epochs.

Horizons uses DE441 by default. DE441 and DE440 agree to well under a kilometre
for the planets within 1850-2150, so these values test the DE440 evaluator.

Usage::

    python tests/fetch_fixtures.py
"""

from __future__ import annotations

import datetime as dt
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://ssd.jpl.nasa.gov/api/horizons.api"
OUT = Path(__file__).resolve().parent / "fixtures" / "horizons.json"

TDB_EPOCHS = [2396800.5, 2415020.5, 2446895.5, 2451545.0, 2461335.75, 2488069.5, 2506200.5]
VECTOR_TARGETS = ["10", "199", "299", "399", "301", "4", "5", "6", "7", "8", "9"]
UTC_EPOCHS = [2446895.5, 2451545.0, 2461336.0]
OBSERVED_TARGETS = ["10", "301", "499", "599", "699"]


def query(params: dict) -> dict:
    fixed = {"format": "json", "MAKE_EPHEM": "YES", "OBJ_DATA": "NO", "CSV_FORMAT": "YES"}
    url = API + "?" + urllib.parse.urlencode({**fixed, **params})
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = json.load(response)
    if "error" in payload:
        raise RuntimeError(payload["error"])
    return payload


def rows(result: str) -> list[list[str]]:
    body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
    return [[c.strip() for c in line.split(",")] for line in body.strip().splitlines()]


def source(result: str) -> str:
    match = re.search(r"\{source:\s*([^}]+)\}", result)
    return match.group(1).strip() if match else "unknown"


def tlist(epochs: list[float]) -> str:
    return " ".join(f"'{jd}'" for jd in epochs)


def main() -> int:
    vectors = {}
    sources = set()
    for target in VECTOR_TARGETS:
        payload = query({
            "COMMAND": f"'{target}'", "EPHEM_TYPE": "VECTORS", "CENTER": "'500@0'",
            "TLIST": tlist(TDB_EPOCHS), "TLIST_TYPE": "JD", "TIME_TYPE": "TDB",
            "REF_PLANE": "FRAME", "REF_SYSTEM": "ICRF", "VEC_TABLE": "2",
            "OUT_UNITS": "KM-S", "VEC_LABELS": "NO",
        })
        sources.add(source(payload["result"]))
        vectors[target] = [
            {"jd_tdb": float(r[0]), "p": [float(x) for x in r[2:5]], "v": [float(x) for x in r[5:8]]}
            for r in rows(payload["result"])
        ]
        print(f"vectors {target}: {len(vectors[target])}")

    astrometric = {}
    for target in OBSERVED_TARGETS:
        payload = query({
            "COMMAND": f"'{target}'", "EPHEM_TYPE": "OBSERVER", "CENTER": "'500@399'",
            "TLIST": tlist(UTC_EPOCHS), "TLIST_TYPE": "JD", "TIME_TYPE": "UT",
            "QUANTITIES": "'1'", "ANG_FORMAT": "DEG", "EXTRA_PREC": "YES",
        })
        sources.add(source(payload["result"]))
        astrometric[target] = [
            {"jd_utc": jd, "ra_deg": float(r[3]), "dec_deg": float(r[4])}
            for jd, r in zip(UTC_EPOCHS, rows(payload["result"]))
        ]
        print(f"astrometric {target}: {len(astrometric[target])}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "source": "JPL Horizons API, " + API,
        "ephemerides": sorted(sources),
        "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "vectors_ssb_icrf_km": vectors,
        "astrometric_geocentric_icrf": astrometric,
    }, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
