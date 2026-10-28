"""Reference values for the navigation bridge's SGP4 and TEME checks (orbits/tests/check_sgp4.cpp).

Writes, beside Vallado's verification files (SGP4-VER.TLE and tcppver.out, copied unaltered from
python-sgp4 2.27, which distributes them as Vallado does):

- sgp4-reference.txt: python-sgp4's compiled copy of Vallado's C++ (vallado_cpp) for every set of
  SGP4-VER.TLE at a few times, in the two cases tcppver.out does not cover: AFSPC's operation mode
  with WGS 72, and the improved mode with WGS 84.
- teme-reference.txt: the TEME → ITRF example of Vallado et al. 2006 (AIAA 2006-6753 rev. 2,
  Appendix C) computed with skyfield at the instant exactly (UT1 as two parts), with the
  rotation's rate taken as skyfield takes it (dθ_GMST82/dt).

Deterministic; run with the reference environment's Python (python-sgp4 2.27, skyfield 1.55):

    python -I domains/navigation/tests/make_sgp4_fixtures.py
"""

from __future__ import annotations

from pathlib import Path

import sgp4
from sgp4.api import WGS72, WGS84, Satrec, accelerated

HERE = Path(__file__).resolve().parent / "fixtures"
TIMES = (0.0, -720.0, 360.0, 1440.0, 4320.0)


def sets():
    lines = (HERE / "SGP4-VER.TLE").read_text(encoding="ascii").splitlines()
    for i, line in enumerate(lines):
        if line.startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            yield line, lines[i + 1]


def sgp4_reference() -> None:
    assert accelerated, "the compiled vallado_cpp is the reference"
    rows = [
        f"# python-sgp4 {sgp4.__version__} (vallado_cpp, Vallado's C++), computed 2026-10-28 by domains/navigation/tests/make_sgp4_fixtures.py",
        "# for each set of SGP4-VER.TLE (in its order): catalogue gravity opsmode tsince(min) error x y z (km) vx vy vz (km/s)",
        "# initialised by sgp4init with the TLE's elements and the epoch as Vallado's twoline2rv forms it, (jdsatepoch + jdsatepochF) - 2433281.5",
    ]
    for line1, line2 in sets():
        t = Satrec.twoline2rv(line1, line2, WGS72)
        epoch = (t.jdsatepoch + t.jdsatepochF) - 2433281.5
        for gravity, name, mode in ((WGS72, "wgs72", "a"), (WGS84, "wgs84", "i")):
            s = Satrec()
            s.sgp4init(gravity, mode, t.satnum, epoch, t.bstar, t.ndot, t.nddot, t.ecco, t.argpo, t.inclo, t.mo, t.no_kozai, t.nodeo)
            for tsince in TIMES:
                e, r, v = s.sgp4_tsince(tsince)
                values = " ".join(f"{x:.15e}" for x in (*r, *v))
                rows.append(f"{t.satnum_str} {name} {mode} {tsince:.1f} {e} {values}")
    (HERE / "sgp4-reference.txt").write_text("\n".join(rows) + "\n", encoding="ascii")


def teme_reference() -> None:
    import skyfield
    from numpy import array
    from skyfield.constants import DEG2RAD
    from skyfield.sgp4lib import TEME_to_ITRF

    arcsecond = DEG2RAD / 3600.0
    r_teme = array([5094.18016210, 6127.64465950, 6380.34453270])
    v_teme = array([-4.746131487, 0.785818041, 5.531931288])
    # 2004 April 6, 07:51:28.386009 UTC; UT1 − UTC = −0.4399619 s; x_p = −0.140682″, y_p = 0.333309″.
    d1 = 2453101.5
    d2 = (7 * 3600 + 51 * 60 + 28.386009 - 0.4399619) / 86400.0
    r, v = TEME_to_ITRF(d1, r_teme, v_teme * 86400.0, -0.140682 * arcsecond, 0.333309 * arcsecond, d2)
    v = v / 86400.0
    rows = [
        f"# skyfield {skyfield.__version__} TEME_to_ITRF (GMST 1982, AIAA 2006-6753 Appendix C), computed 2026-10-28",
        "# by domains/navigation/tests/make_sgp4_fixtures.py: the paper's example at its instant exactly (UT1 in two parts),",
        "# the rotation's rate dθ_GMST82/dt (skyfield's), no LOD. ut1 d1 d2; r_itrf (km); v_itrf (km/s)",
        f"ut1 {d1!r} {d2!r}",
        "r " + " ".join(f"{x:.15e}" for x in r),
        "v " + " ".join(f"{x:.15e}" for x in v),
    ]
    (HERE / "teme-reference.txt").write_text("\n".join(rows) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sgp4_reference()
    teme_reference()
    print("wrote", HERE / "sgp4-reference.txt", "and", HERE / "teme-reference.txt")
