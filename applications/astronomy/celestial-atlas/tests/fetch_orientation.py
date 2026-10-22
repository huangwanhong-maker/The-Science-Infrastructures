"""Fetch the sub-Earth points of the Sun, the Moon and the planets from JPL Horizons.

The sub-observer point of a body is where the line from its centre to the
observer pierces its surface; its longitude tests the body's rotation model
(the prime meridian W) and its latitude the direction of the pole. Horizons
computes it from the IAU rotation models (the Moon's from its mean-Earth frame)
with light-time, in planetodetic coordinates, with longitudes counted west on
prograde rotators and east on the others, as each table's header states.

Writes tests/fixtures/orientation.json.

Usage::

    python applications/astronomy/celestial-atlas/tests/fetch_orientation.py
"""

from __future__ import annotations

import datetime
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://ssd.jpl.nasa.gov/api/horizons.api"
OUT = Path(__file__).resolve().parent / "fixtures" / "orientation.json"
TARGETS = {"sun": "10", "mercury": "199", "venus": "299", "moon": "301", "mars": "499", "jupiter": "599",
           "saturn": "699", "uranus": "799", "neptune": "899", "pluto": "999"}
EPOCHS = [2451545.0, 2455000.25, 2457388.5, 2461335.5, 2465000.75]     # JD (TT)


def main() -> None:
    entries = []
    for key, target in TARGETS.items():
        params = {"format": "json", "COMMAND": f"'{target}'", "CENTER": "'500@399'", "MAKE_EPHEM": "YES",
                  "EPHEM_TYPE": "OBSERVER", "TLIST": " ".join(f"'{t}'" for t in EPOCHS), "TLIST_TYPE": "JD",
                  "TIME_TYPE": "TT", "QUANTITIES": "'14'", "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES",
                  "EXTRA_PREC": "YES", "OBJ_DATA": "NO"}
        with urllib.request.urlopen(API + "?" + urllib.parse.urlencode(params), timeout=60) as r:
            result = json.load(r)["result"]
        west = "West-longitude positive" in result
        frame = re.search(r"Target pole/equ : (\S+)", result).group(1)
        body = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0]
        for jd, line in zip(EPOCHS, body.strip().splitlines()):
            cells = [c.strip() for c in line.split(",")]
            entries.append({"body": key, "jd_tt": jd, "lon_deg": float(cells[3]), "lat_deg": float(cells[4]),
                            "west_positive": west, "frame": frame})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "source": "JPL Horizons API: observer sub-longitude and sub-latitude (quantity 14) of each body, "
                  "seen from the Earth's centre; planetodetic; Moon in MOON_ME.",
        "retrieved_at": datetime.date.today().isoformat(),
        "entries": entries,
    }, indent=1) + "\n", encoding="utf-8")
    print(f"{len(entries)} sub-Earth points written to {OUT}")


if __name__ == "__main__":
    main()
