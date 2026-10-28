"""Reference values for the geodesy domain's map projections and grids, from PROJ and GEOTRANS.

PROJ (through pyproj) gives the UTM coordinates of points in their zones (EPSG:326zz north, 327zz
south; PROJ's extended transverse Mercator of Poder and Engsager) with the meridian convergence and
point scale there (Proj.get_factors), and the UPS coordinates of the polar caps (EPSG:32661 north,
32761 south); NGA's GEOTRANS (through the mgrs package) gives MGRS references at 1 m. This script
evaluates them at seeded random points and at the grids' special places (zone edges, Norway,
Svalbard, the band edges, the poles), and writes tests/fixtures/projection-reference.txt, which the
C++ checks read. It needs pyproj and mgrs:

    python domains/geodesy/tests/make_projection_fixtures.py
"""

from __future__ import annotations

import datetime
from pathlib import Path

import mgrs
import numpy as np
import pyproj

OUT = Path(__file__).resolve().parent / "fixtures" / "projection-reference.txt"


def zone_of(lat: float, lon: float) -> int:
    if lat >= 56 and lat < 64 and 3 <= lon < 12:
        return 32
    if 72 <= lat < 84 and 0 <= lon < 42:
        return 31 if lon < 9 else 33 if lon < 21 else 35 if lon < 33 else 37
    return min(60, int((lon + 180) // 6) + 1)


def main() -> None:
    rng = np.random.default_rng(20261028)
    rows = [f"# pyproj {pyproj.__version__} (PROJ {pyproj.proj_version_str}) and mgrs (GEOTRANS), computed {datetime.date.today().isoformat()}",
            "# by domains/geodesy/tests/make_projection_fixtures.py. Degrees and metres; convergence in degrees as PROJ gives it.",
            "# utm zone north lat lon easting northing convergence scale | ups north lat lon easting northing | mgrs lat lon reference"]
    # UTM: random points in random zones, out to 4 degrees from the central meridian, and in the exceptions.
    points = []
    for _ in range(300):
        lat = float(rng.uniform(-79.9, 83.9))
        zone = int(rng.integers(1, 61))
        lon = float((6 * zone - 183) + rng.uniform(-4.0, 4.0))
        points.append((zone, lat, lon))
    for lat, lon in [(60, 5), (60, 10.5), (63.9, 3.1), (75, 8.9), (75, 9.1), (78, 20.9), (80, 32.9), (83.9, 41.9), (0, 0), (0, -0.0001),
                     (-79.9, 177), (83.99, -179.99), (45, 2.999999), (45, 3.000001), (-0.0001, 0.0001)]:
        points.append((zone_of(lat, lon), lat, lon))
    for zone, lat, lon in points:
        north = lat >= 0
        crs = pyproj.CRS.from_epsg((32600 if north else 32700) + zone)
        t = pyproj.Transformer.from_crs("EPSG:4326", crs, always_xy=True)
        e, n = t.transform(lon, lat)
        f = pyproj.Proj(crs).get_factors(lon, lat)
        rows.append(f"utm {zone} {int(north)} {lat!r} {lon!r} {e!r} {n!r} {f.meridian_convergence!r} {f.meridional_scale!r}")
    # UPS: the caps beyond UTM.
    for north, crs in ((1, 32661), (0, 32761)):
        t = pyproj.Transformer.from_crs("EPSG:4326", f"EPSG:{crs}", always_xy=True)
        for _ in range(60):
            lat = float(rng.uniform(84, 90)) if north else float(rng.uniform(-90, -80))
            lon = float(rng.uniform(-180, 180))
            e, n = t.transform(lon, lat)
            rows.append(f"ups {north} {lat!r} {lon!r} {e!r} {n!r}")
        lat = 90.0 if north else -90.0
        e, n = t.transform(0.0, lat)
        rows.append(f"ups {north} {lat!r} 0.0 {e!r} {n!r}")
    # MGRS: random points everywhere, and the special ones.
    m = mgrs.MGRS()
    specials = [(60, 5), (60, 10), (75, 10), (75, 20), (85, 0.5), (85, -0.5), (85, 90), (85, -90), (-81, 0.5), (-81, -0.5),
                (-81, 90), (-81, -90), (-85, -60), (0, 0.1), (0.0001, -0.1), (-50, -70), (48.8566, 2.3522), (35.6762, 139.6503)]
    for _ in range(400):
        lat = float(np.degrees(np.arcsin(rng.uniform(-1, 1))))
        lon = float(rng.uniform(-180, 180))
        specials.append((lat, lon))
    for lat, lon in specials:
        rows.append(f"mgrs {lat!r} {lon!r} {m.toMGRS(lat, lon, MGRSPrecision=5)}")
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(rows) - 3} rows")


if __name__ == "__main__":
    main()
