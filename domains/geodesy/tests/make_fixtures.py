"""Reference values for the geodesy domain's C++ checks, computed with ERFA through pyerfa.

ERFA's gd2gc and gc2gd convert between geodetic and Earth-centred coordinates on WGS 84 and
GRS 80 (Fukushima 2006 for the inverse). This script evaluates them at points from the poles to
the equator, from 10 km below the ellipsoid to geostationary height, and writes
tests/fixtures/geodetic-reference.txt, which the C++ checks read. It needs pyerfa and numpy:

    python domains/geodesy/tests/make_fixtures.py
"""

from __future__ import annotations

import datetime
from pathlib import Path

import erfa
import numpy as np

OUT = Path(__file__).resolve().parent / "fixtures" / "geodetic-reference.txt"


def main() -> None:
    rows = [f"# ERFA {erfa.__version__} (SOFA {erfa.version.sofa_version}) via pyerfa, computed {datetime.date.today().isoformat()}",
            "# by domains/geodesy/tests/make_fixtures.py. gd2gc n lon lat h x y z; gc2gd n x y z lon lat h. Radians, metres."]
    rng = np.random.default_rng(20261028)
    points = [(0.0, 0.0, 0.0), (0.0, np.pi / 2, 0.0), (1.0, -np.pi / 2, 100.0), (2.0, 1e-9, 1.0), (-2.5, 0.7, -10000.0),
              (0.3, 0.9, 35786000.0), (3.1, -1.2, 8848.86), (-1.0, 1.5707963, 0.0), (0.5, -0.5, 400000.0)]
    for _ in range(40):
        points.append((rng.uniform(-np.pi, np.pi), rng.uniform(-np.pi / 2, np.pi / 2), rng.uniform(-5000, 2e7)))
    for n in (1, 2):  # 1 = WGS 84, 2 = GRS 80
        for lon, lat, h in points:
            x = erfa.gd2gc(n, lon, lat, h)
            rows.append(" ".join(["gd2gc", str(n)] + [repr(float(v)) for v in (lon, lat, h, *x)]))
            lo, la, hh = erfa.gc2gd(n, x)
            rows.append(" ".join(["gc2gd", str(n)] + [repr(float(v)) for v in (*x, lo, la, hh)]))
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT} ({len(rows)} lines)")


if __name__ == "__main__":
    main()
