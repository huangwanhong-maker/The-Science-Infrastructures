"""The astronomy domain's Horizons fixtures as plain text, for its C++ checks.

The C++ checks read no JSON; this script writes the values of fixtures/horizons.json (JPL
Horizons barycentric vectors) and fixtures/orientation.json (Horizons sub-Earth points) to
fixtures/horizons.txt and fixtures/orientation.txt, unchanged, one value per field. Run it when
the JSON fixtures change:

    python domains/astronomy/tests/cpp_fixtures.py
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent / "fixtures"


def main() -> None:
    h = json.loads((HERE / "horizons.json").read_text(encoding="utf-8"))
    rows = [f"# From fixtures/horizons.json: {h['source']}, {', '.join(h['ephemerides'])}, retrieved {h['retrieved_at']}.",
            "# ssb <NAIF id> <JD TDB> <x y z km> <vx vy vz km/s>, ICRF, relative to the solar-system barycentre."]
    for target, entries in h["vectors_ssb_icrf_km"].items():
        for e in entries:
            rows.append(" ".join(["ssb", target, repr(e["jd_tdb"])] + [repr(x) for x in e["p"]] + [repr(x) for x in e["v"]]))
    (HERE / "horizons.txt").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    o = json.loads((HERE / "orientation.json").read_text(encoding="utf-8"))
    rows = [f"# From fixtures/orientation.json: {o['source']} Retrieved {o['retrieved_at']}.",
            "# sub <body> <JD TT> <longitude deg> <latitude deg> <1 if west-positive> <frame>"]
    for e in o["entries"]:
        rows.append(" ".join(["sub", e["body"], repr(e["jd_tt"]), repr(e["lon_deg"]), repr(e["lat_deg"]), "1" if e["west_positive"] else "0", e["frame"]]))
    (HERE / "orientation.txt").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print("wrote fixtures/horizons.txt and fixtures/orientation.txt")


if __name__ == "__main__":
    main()
