"""Reference values of geodesics on WGS 84, for the checks of geodesy/geodesic.

Computed with GeographicLib's Python package (Karney), whose algorithms the C++ ports, at cases
chosen to reach every branch of the inverse problem: points anywhere, nearly antipodal pairs (the
astroid), poles, the equator (and equatorial lines near the antipode), meridians, short and very
short lines; and the direct problem from the same starting points. GeographicLib's own documented
test cases (its test_geodesic.py, values it gives with their sources) are written as they are.

Run with the reference environment (geographiclib 2.1):
    python -I domains/geodesy/tests/make_geodesic_fixtures.py
which writes domains/geodesy/tests/fixtures/geodesic-reference.txt.
"""

from __future__ import annotations

import math
import random
from pathlib import Path

import geographiclib
from geographiclib.geodesic import Geodesic
from geographiclib.test.test_geodesic import GeodesicTest

OUT = Path(__file__).resolve().parent / "fixtures" / "geodesic-reference.txt"
G = Geodesic.WGS84
ALL = Geodesic.ALL
UNROLL = Geodesic.ALL | Geodesic.LONG_UNROLL


def fmt(x: float) -> str:
    return repr(float(x))


def main() -> None:
    rnd = random.Random(20261028)
    pairs: list[tuple[str, float, float, float, float]] = []

    def uniform_lat() -> float:
        return math.degrees(math.asin(2 * rnd.random() - 1))

    for _ in range(300):
        pairs.append(("random", uniform_lat(), rnd.uniform(-180, 180), uniform_lat(), rnd.uniform(-180, 180)))
    for _ in range(120):  # nearly antipodal: the astroid and the strip near the cut
        lat1 = rnd.uniform(-89, 89)
        lon1 = rnd.uniform(-180, 180)
        dlat = 10 ** rnd.uniform(-9, 0) * rnd.choice([-1, 1])
        dlon = 10 ** rnd.uniform(-9, 0.5)
        pairs.append(("antipodal", lat1, lon1, -lat1 + dlat, lon1 + 180 - dlon))
    for _ in range(40):  # from or to a pole
        pole = rnd.choice([-90.0, 90.0])
        if rnd.random() < 0.5:
            pairs.append(("pole", pole, rnd.uniform(-180, 180), uniform_lat(), rnd.uniform(-180, 180)))
        else:
            pairs.append(("pole", uniform_lat(), rnd.uniform(-180, 180), pole, rnd.uniform(-180, 180)))
    for _ in range(40):  # along the equator, some near the antipode
        lon1 = rnd.uniform(-180, 180)
        dlon = rnd.choice([rnd.uniform(0, 179), 180 - 10 ** rnd.uniform(-6, 0.3)])
        pairs.append(("equator", 0.0, lon1, 0.0, lon1 + rnd.choice([-1, 1]) * dlon))
    for _ in range(40):  # meridians, and through the pole
        lon1 = rnd.uniform(-180, 180)
        pairs.append(("meridian", uniform_lat(), lon1, uniform_lat(), lon1 + rnd.choice([0.0, 180.0])))
    for _ in range(80):  # short and very short lines
        lat1 = uniform_lat() * 0.99
        lon1 = rnd.uniform(-180, 180)
        d = 10 ** rnd.uniform(-8, 0)
        a = rnd.uniform(0, 2 * math.pi)
        pairs.append(("short", lat1, lon1, max(-90.0, min(90.0, lat1 + d * math.cos(a))), lon1 + d * math.sin(a)))

    lines = [
        f"# GeographicLib {geographiclib.__version__} (Python, C. F. F. Karney), Geodesic.WGS84 (a = {G.a}, f = 1/298.257223563),",
        "# computed 2026-10-28 by domains/geodesy/tests/make_geodesic_fixtures.py (seed 20261028). Degrees and metres.",
        "# inv kind lat1 lon1 lat2 lon2 | s12 azi1 azi2 a12 m12 M12 M21 S12              (Inverse, ALL)",
        "# dir kind lat1 lon1 azi1 s12 | lat2 lon2 azi2 a12 m12 M12 M21 S12              (Direct, ALL | LONG_UNROLL)",
        "# doc lat1 lon1 azi1 lat2 lon2 azi2 s12 a12 m12 M12 M21 S12                      (GeographicLib's test cases)",
    ]
    for kind, lat1, lon1, lat2, lon2 in pairs:
        r = G.Inverse(lat1, lon1, lat2, lon2, ALL)
        lines.append("inv " + kind + " " + " ".join(fmt(x) for x in (lat1, lon1, lat2, lon2)) + " | " +
                     " ".join(fmt(r[k]) for k in ("s12", "azi1", "azi2", "a12", "m12", "M12", "M21", "S12")))
    for kind, lat1, lon1, lat2, lon2 in pairs[:400]:
        azi1 = rnd.uniform(-180, 180)
        s12 = 10 ** rnd.uniform(-3, 7.6)  # 1 mm to 40 000 km: past the antipode and round again
        r = G.Direct(lat1, lon1, azi1, s12, UNROLL)
        lines.append("dir " + kind + " " + " ".join(fmt(x) for x in (lat1, lon1, azi1, s12)) + " | " +
                     " ".join(fmt(r[k]) for k in ("lat2", "lon2", "azi2", "a12", "m12", "M12", "M21", "S12")))
    for case in GeodesicTest.testcases:
        lines.append("doc " + " ".join(fmt(x) for x in case))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="ascii")
    print(f"wrote {OUT}: {len(pairs)} inverse, {min(400, len(pairs))} direct, {len(GeodesicTest.testcases)} documented")


if __name__ == "__main__":
    main()
