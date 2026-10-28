"""Reference values for the metrology domain's C++ checks, computed with ERFA through pyerfa.

ERFA (BSD 3-clause) is the IAU's SOFA library under another name; its results are the
standard against which the domain's time scales and Earth orientation are checked. This
script evaluates ERFA's routines at dates spread over 1850-2200 (and around leap seconds) and
writes them, one per line, to tests/fixtures/erfa-reference.txt, which the C++ checks read.
It needs pyerfa (pip install pyerfa) and is run by hand when the fixtures are to be renewed:

    python domains/metrology/tests/make_fixtures.py

Lines: a name, then numbers in the order the C++ reads them (see check_time.cpp and
check_frames.cpp).
"""

from __future__ import annotations

import datetime
from pathlib import Path

import erfa
import numpy as np

OUT = Path(__file__).resolve().parent / "fixtures" / "erfa-reference.txt"

# Two-part TT dates (MJD zero point + MJD), from 1850 to 2200, with awkward fractions.
DATES = [
    (2400000.5, -3000.25), (2400000.5, 15020.0), (2400000.5, 33282.5), (2400000.5, 41317.0),
    (2400000.5, 50123.9999), (2400000.5, 51544.5), (2400000.5, 53736.0), (2400000.5, 57754.0),
    (2400000.5, 60000.123456789), (2400000.5, 61341.375), (2400000.5, 69807.0), (2400000.5, 88069.75),
    (2451545.0, 0.0), (2451545.0, -36525.0), (2451545.0, 73050.5),
]


def line(*values) -> str:
    out = []
    for v in values:
        if isinstance(v, str):
            out.append(v)
        elif isinstance(v, (int, np.integer)):
            out.append(str(int(v)))
        else:
            out.append(repr(float(v)))
    return " ".join(out)


def main() -> None:
    rows = [
        f"# ERFA {erfa.__version__} (SOFA {erfa.version.sofa_version}) via pyerfa, computed {datetime.date.today().isoformat()}",
        "# by domains/metrology/tests/make_fixtures.py. Units: radians, days, seconds.",
    ]
    for d1, d2 in DATES:
        dpsi, deps = erfa.nut06a(d1, d2)
        rows.append(line("nut06a", d1, d2, dpsi, deps))
        dpsi0, deps0 = erfa.nut00a(d1, d2)
        rows.append(line("nut00a", d1, d2, dpsi0, deps0))
        gamb, phib, psib, epsa = erfa.pfw06(d1, d2)
        rows.append(line("pfw06", d1, d2, gamb, phib, psib, epsa))
        rows.append(line("obl06", d1, d2, erfa.obl06(d1, d2)))
        r = erfa.pnm06a(d1, d2)
        rows.append(line("pnm06a", d1, d2, *r.ravel()))
        x, y = erfa.bpn2xy(r)
        rows.append(line("s06", d1, d2, x, y, erfa.s06(d1, d2, x, y)))
        rows.append(line("eect00", d1, d2, erfa.eect00(d1, d2)))
        rows.append(line("era00", d1, d2, erfa.era00(d1, d2)))
        rows.append(line("gmst06", d1, d2, d1, d2, erfa.gmst06(d1, d2, d1, d2)))
        rows.append(line("gmst82", d1, d2, erfa.gmst82(d1, d2)))
        rows.append(line("gst06a", d1, d2, d1, d2, erfa.gst06a(d1, d2, d1, d2)))
        # The classical apparent sidereal time, GMST + Δψ cos ε_A + complementary terms.
        ee = erfa.ee00(d1, d2, epsa, dpsi)
        rows.append(line("gast_classical", d1, d2, d1, d2, erfa.anp(erfa.gmst06(d1, d2, d1, d2) + ee)))
        rows.append(line("sp00", d1, d2, erfa.sp00(d1, d2)))
        rows.append(line("dtdb", d1, d2, erfa.dtdb(d1, d2, 0.0, 0.0, 0.0, 0.0)))
        g1, g2 = erfa.tttcg(d1, d2)
        rows.append(line("tttcg", d1, d2, g1, g2))
        b1, b2 = erfa.tdbtcb(d1, d2)
        rows.append(line("tdbtcb", d1, d2, b1, b2))
        xp, yp = 0.25 * erfa.DAS2R, 0.35 * erfa.DAS2R
        rc2t = erfa.c2t06a(d1, d2, d1, d2 + 0.3 / 86400, xp, yp)
        rows.append(line("c2t06a", d1, d2, d1, d2 + 0.3 / 86400, xp, yp, *rc2t.ravel()))
        y_, m_, dd, fd = erfa.jd2cal(d1, d2)
        rows.append(line("jd2cal", d1, d2, y_, m_, dd, fd))

    # Calendar dates to JD.
    for (y, m, d) in [(1850, 1, 1), (1900, 2, 28), (1960, 1, 1), (1972, 6, 30), (2000, 2, 29), (2016, 12, 31), (2026, 10, 28), (2199, 12, 31)]:
        djm0, djm = erfa.cal2jd(y, m, d)
        rows.append(line("cal2jd", y, m, d, djm0, djm))

    # TAI − UTC, including the drifting years and the edges of leap seconds.
    for (y, m, d, fd) in [(1960, 1, 1, 0.0), (1962, 6, 15, 0.5), (1965, 3, 1, 0.25), (1968, 2, 1, 0.0), (1971, 12, 31, 0.999),
                          (1972, 1, 1, 0.0), (1972, 6, 30, 0.5), (1972, 7, 1, 0.0), (1999, 1, 1, 0.0), (2016, 12, 31, 0.9),
                          (2017, 1, 1, 0.0), (2026, 10, 28, 0.5)]:
        rows.append(line("dat", y, m, d, fd, erfa.dat(y, m, d, fd)))

    # UTC to TAI and back, around the leap seconds of 1972, 2005 and 2016, and in the 1960s.
    for (u1, u2) in [(2441317.5, 0.0), (2441499.5, -0.000001), (2453735.5, 0.99999), (2453736.5, 0.00001), (2457753.5, 0.99999),
                     (2457753.5, 0.9999999), (2457754.5, 0.0), (2438000.5, 0.37), (2460000.5, 0.123456), (2453750.5, 0.892100694)]:
        t1, t2 = erfa.utctai(u1, u2)
        rows.append(line("utctai", u1, u2, t1, t2))
        v1, v2 = erfa.taiutc(t1, t2)
        rows.append(line("taiutc", t1, t2, v1, v2))

    # Date-time text to JD and back, UTC with leap seconds and TT.
    for (scale, y, m, d, h, mi, s) in [("UTC", 2016, 12, 31, 23, 59, 60.5), ("UTC", 1994, 6, 30, 23, 59, 60.13599),
                                         ("UTC", 2026, 10, 28, 12, 34, 56.789), ("TT", 2000, 1, 1, 12, 0, 0.0),
                                         ("UTC", 1972, 6, 30, 23, 59, 60.0), ("TT", 1850, 7, 4, 6, 30, 15.25)]:
        a, b = erfa.dtf2d(scale, y, m, d, h, mi, s)
        rows.append(line("dtf2d", scale, y, m, d, h, mi, s, a, b))
        yy, mm, dd, ihmsf = erfa.d2dtf(scale, 5, a, b)
        rows.append(line("d2dtf", scale, 5, a, b, yy, mm, dd, *ihmsf))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT} ({len(rows)} lines)")


if __name__ == "__main__":
    main()
