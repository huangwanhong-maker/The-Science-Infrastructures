"""Reference values for the navigation bridge's galactocentric frame, from astropy.

astropy.coordinates carries stars from the ICRS (barycentric, with proper motions and radial
velocities) into its Galactocentric frame with the default parameters ("v4.0": Sgr A* at RA
266.4051°, Dec −28.936175°, at 8.122 kpc, the Sun 20.8 pc above the plane, moving at (12.9, 245.6,
7.78) km/s). This script draws seeded random stars from 1 pc to 50 kpc, writes their ICRS
Cartesian positions and velocities as astropy computes them, and the Galactocentric ones, to
tests/fixtures/galactocentric-reference.txt, which the C++ checks of the frame read.

    python domains/navigation/tests/make_frames_fixtures.py
"""

from __future__ import annotations

import datetime
from pathlib import Path

import astropy
import astropy.coordinates as coord
import astropy.units as u
import numpy as np

OUT = Path(__file__).resolve().parent / "fixtures" / "galactocentric-reference.txt"


def main() -> None:
    rng = np.random.default_rng(20261028)
    n = 300
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    distance = 10 ** rng.uniform(0, np.log10(50000), n)       # pc
    pmra = rng.normal(0, 50, n)                                # mas/yr, μ_α cos δ
    pmdec = rng.normal(0, 50, n)
    rv = rng.normal(0, 60, n)                                  # km/s
    icrs = coord.ICRS(ra=ra * u.deg, dec=dec * u.deg, distance=distance * u.pc,
                      pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr, radial_velocity=rv * u.km / u.s)
    gc = icrs.transform_to(coord.Galactocentric())
    x = icrs.cartesian.xyz.to(u.m).value
    v = icrs.velocity.d_xyz.to(u.m / u.s).value
    X = gc.cartesian.xyz.to(u.m).value
    V = gc.velocity.d_xyz.to(u.m / u.s).value
    rows = [f"# astropy {astropy.__version__}, Galactocentric() with its default parameters (v4.0); computed {datetime.date.today().isoformat()}",
            "# by domains/navigation/tests/make_frames_fixtures.py. ICRS x y z (m) vx vy vz (m/s), then Galactocentric x y z (m) vx vy vz (m/s)."]
    for i in range(n):
        rows.append(" ".join(repr(float(q)) for q in (*x[:, i], *v[:, i], *X[:, i], *V[:, i])))
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {n} stars")


if __name__ == "__main__":
    main()
