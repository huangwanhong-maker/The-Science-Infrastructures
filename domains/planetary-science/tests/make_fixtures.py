"""Reference values for the Moon's crustal field (planetary-science/magnetism), from pyshtools.

The crustal field of T2015_449 (Tsunakawa et al. 2015, expanded by Wieczorek; the data library's
moon-magnetic-t2015) evaluated by pyshtools (Wieczorek and Meschede 2018) at mounts on and above the
Moon's sphere of 1737.4 km, in the Moon's mean-Earth frame: north = −B_θ, east = B_φ, down = −B_r.
magnetism/tests/check_magnetism.cpp evaluates the same points, and Cosmopoiesis's tests/check_runtime.py
mounts magnetometers at the same addresses and compares.

    .venv/Scripts/python -I domains/planetary-science/tests/make_fixtures.py      (needs pyshtools 4.13 or later)
"""

import datetime
from pathlib import Path

import pyshtools as pysh

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MODEL = ROOT / "data" / "cache" / "moon-magnetic-t2015" / "T2015_449.sh"
OUT = HERE / "fixtures" / "lunar-field.txt"

# (name, latitude °N, longitude °E, height m): the swirl of Reiner Gamma on the ground and 30 km over
# it, Tranquility Base, the Gerasimovich anomaly on the far side, and a point near the south pole.
MOUNTS = [
    ("reiner_gamma", 7.5, -59.0, 0.0),
    ("reiner_gamma_30km", 7.5, -59.0, 30e3),
    ("tranquility", 0.674, 23.473, 0.0),
    ("gerasimovich", -23.25, -123.25, 0.0),
    ("south", -80.0, 10.0, 10e3),
]


def main() -> None:
    c = pysh.SHMagCoeffs.from_file(str(MODEL), format="shtools", header=True, r0_index=0, file_units="T", units="T")
    lines = [f"# T2015_449 (moon-magnetic-t2015) by pyshtools {pysh.__version__}, {datetime.date.today().isoformat()}; r0 = {c.r0:.1f} m, degree {c.lmax}.",
             "# name latitude_deg longitude_deg height_m north_T east_T down_T (the Moon's mean-Earth frame)"]
    for name, lat, lon, h in MOUNTS:
        br, bt, bp = c.expand(lat=lat, lon=lon, r=c.r0 + h)
        lines.append(f"{name} {lat:.6f} {lon:.6f} {h:.1f} {-bt:.17e} {bp:.17e} {-br:.17e}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
