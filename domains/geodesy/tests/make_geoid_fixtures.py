"""Reference values for the geodesy domain's geoid: EGM96 undulations from PROJ.

PROJ (through pyproj, with its network enabled) carries a point from WGS 84 ellipsoidal heights
(EPSG:4979) to EGM96 heights (EPSG:4326+5773) by NGA's EGM96 grid at 15′ (us_nga_egm96_15.tif from
cdn.proj.org), interpolated bilinearly; the undulation is N = h − H. This script evaluates it at
seeded random points and at a few named ones, and writes tests/fixtures/geoid-reference.txt, which
the C++ checks read beside the data library's own copy of the grid (GeographicLib's egm96-15.pgm).

    python domains/geodesy/tests/make_geoid_fixtures.py
"""

from __future__ import annotations

import datetime
from pathlib import Path

import numpy as np
import pyproj

OUT = Path(__file__).resolve().parent / "fixtures" / "geoid-reference.txt"


def main() -> None:
    pyproj.network.set_network_enabled(True)
    t = pyproj.Transformer.from_crs("EPSG:4979", "EPSG:4326+5773", always_xy=True)
    rng = np.random.default_rng(20261028)
    points = [(0.0, 0.0), (-3.009444, 16.775833), (138.7274, 35.3606), (78.5, 5.0), (-75.0, -20.0), (147.0, -6.0), (0.0, 89.9), (10.0, -89.9)]
    for _ in range(400):
        points.append((float(rng.uniform(-180, 180)), float(np.degrees(np.arcsin(rng.uniform(-0.999, 0.999))))))
    rows = [f"# pyproj {pyproj.__version__} (PROJ {pyproj.proj_version_str}), NGA's EGM96 15' grid (us_nga_egm96_15.tif), bilinear; computed {datetime.date.today().isoformat()}",
            "# by domains/geodesy/tests/make_geoid_fixtures.py. lon lat (degrees) N (m)."]
    for lon, lat in points:
        h = t.transform(lon, lat, 0.0)[2]
        rows.append(f"{lon!r} {lat!r} {-h!r}")
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(points)} points")


if __name__ == "__main__":
    main()
