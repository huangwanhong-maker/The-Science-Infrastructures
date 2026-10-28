"""Transcribe the series of the IAU 2006/2000A precession-nutation model and of TDB - TT from ERFA.

The series are numerical tables published by the IERS and the IAU: the IAU 2000A nutation
(Mathews, Herring and Buffett 2002; IERS Conventions 2010, Tables 5.3a and 5.3b), the series
for the CIO locator s + XY/2 (Capitaine, Wallace and Chapront 2005; IERS Conventions 2010,
Table 5.2d), the complementary terms of the equation of the equinoxes (Capitaine and Gontier
1993; IERS Conventions 2010, Table 5.2e) and Fairhead and Bretagnon's series for TDB - TT
(Fairhead and Bretagnon 1990, A&A 229, 240). ERFA (BSD 3-clause; derived with permission
from the IAU's SOFA library) carries them as C arrays; this script reads those arrays from a
copy of ERFA's sources and writes them as C++ headers of the metrology domain, so that the
transcription is mechanical and can be repeated. Run it from the repository root:

    python domains/metrology/tools/transcribe_erfa.py <erfa/src directory>

It writes:
    domains/metrology/frames/include/metrology/frames/series_nutation.hpp
    domains/metrology/frames/include/metrology/frames/series_cio.hpp
    domains/metrology/time/include/metrology/time/series_tdb.hpp

The sources were taken from https://github.com/liberfa/erfa (branch master) on 2026-10-28.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

NOTICE = """// Transcribed mechanically by domains/metrology/tools/transcribe_erfa.py from ERFA's {file}
// (https://github.com/liberfa/erfa, retrieved 2026-10-28). ERFA is distributed under the following
// terms, which this transcription keeps:
//
//   Copyright (C) 2013-2021, NumFOCUS Foundation. All rights reserved. This library is derived,
//   with permission, from the International Astronomical Union's "Standards of Fundamental
//   Astronomy" library, available from http://www.iausofa.org. Redistribution and use in source
//   and binary forms, with or without modification, are permitted provided that the following
//   conditions are met: (1) redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer; (2) redistributions in binary
//   form must reproduce the above copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided with the distribution;
//   (3) neither the name of the Standards Of Fundamental Astronomy Board, the International
//   Astronomical Union nor the names of its contributors may be used to endorse or promote
//   products derived from this software without specific prior written permission. THIS
//   SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
//   IMPLIED WARRANTIES ARE DISCLAIMED.
"""


def block(text: str, start: str) -> str:
    """The text of the C initialiser that begins at `start`, up to its matching close."""
    i = text.index(start)
    i = text.index("{", i + len(start) - 1)
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1 : j]
    raise ValueError(f"unterminated initialiser {start}")


def strip_comments(text: str) -> str:
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)


def groups(body: str) -> list[list[str]]:
    """The innermost {...} groups of an initialiser, each as a list of its numbers."""
    return [[x.strip() for x in g.split(",") if x.strip()] for g in re.findall(r"\{([^{}]*)\}", body)]


def terms(body: str) -> list[tuple[list[str], str, str]]:
    """TERM groups {{n1, ..., n8}, s, c}."""
    out = []
    for m in re.finditer(r"\{\s*\{([^{}]*)\}\s*,\s*([^,{}]+?)\s*,\s*([^,{}]+?)\s*\}", body):
        out.append(([x.strip() for x in m.group(1).split(",")], m.group(2).strip(), m.group(3).strip()))
    return out


def number(x: str) -> str:
    """A C number as C++ writes it, unchanged except for spacing."""
    return re.sub(r"\s+", "", x)


def nutation(src: Path) -> str:
    text = strip_comments((src / "nut00a.c").read_text(encoding="utf-8"))
    xls = groups(block(text, "xls[] = {"))
    xpl = groups(block(text, "xpl[] = {"))
    assert len(xls) == 678 and all(len(g) == 11 for g in xls), (len(xls), set(len(g) for g in xls))
    assert len(xpl) == 687 and all(len(g) == 17 for g in xpl), (len(xpl), set(len(g) for g in xpl))
    lines = [
        "// The IAU 2000A nutation series: 678 luni-solar and 687 planetary terms (Mathews, Herring and",
        "// Buffett 2002, J. Geophys. Res. 107(B4), 2068; IERS Conventions 2010, Tables 5.3a and 5.3b).",
        "// Units: 0.1 microarcsecond, and the same per Julian century.",
        "//",
        NOTICE.format(file="nut00a.c").rstrip(),
        "",
        "#pragma once",
        "",
        "namespace metrology::frames::series {",
        "",
        "// Luni-solar terms: multipliers of l, l', F, D, Ω; then (A + A' t) sin + A'' cos in longitude,",
        "// (B + B' t) cos + B'' sin in obliquity.",
        "struct LuniSolarTerm {",
        "  int nl, nlp, nf, nd, nom;",
        "  double sp, spt, cp, ce, cet, se;",
        "};",
        "// Planetary terms: multipliers of l, F, D, Ω, the longitudes of Mercury to Neptune and the",
        "// general precession in longitude; then the sin and cos coefficients in longitude and obliquity.",
        "struct PlanetaryTerm {",
        "  int nl, nf, nd, nom, nme, nve, nea, nma, nju, nsa, nur, nne, npa;",
        "  int sp, cp, se, ce;",
        "};",
        "",
        "inline constexpr LuniSolarTerm luni_solar[678] = {",
    ]
    for g in xls:
        lines.append("    {" + ", ".join(g[:5]) + ", " + ", ".join(number(x) for x in g[5:]) + "},")
    lines += ["};", "", "inline constexpr PlanetaryTerm planetary[687] = {"]
    for g in xpl:
        lines.append("    {" + ", ".join(g) + "},")
    lines += ["};", "", "}  // namespace metrology::frames::series", ""]
    return "\n".join(lines)


def cio(src: Path) -> str:
    s06 = strip_comments((src / "s06.c").read_text(encoding="utf-8"))
    sp = [number(x) for x in block(s06, "sp[] = {").split(",") if x.strip()]
    assert len(sp) == 6
    tables = {}
    for name in ("s0", "s1", "s2", "s3", "s4"):
        tables[name] = terms(block(s06, f"{name}[] = {{"))
    eect = strip_comments((src / "eect00.c").read_text(encoding="utf-8"))
    e0 = terms(block(eect, "e0[] = {"))
    e1 = terms(block(eect, "e1[] = {"))
    counts = {k: len(v) for k, v in tables.items()}
    assert counts == {"s0": 33, "s1": 3, "s2": 25, "s3": 4, "s4": 1}, counts
    assert len(e0) == 33 and len(e1) == 1, (len(e0), len(e1))

    def emit(name: str, rows) -> list[str]:
        out = [f"inline constexpr ArgumentTerm {name}[{len(rows)}] = {{"]
        for nfa, s, c in rows:
            out.append("    {{" + ", ".join(nfa) + "}, " + number(s) + ", " + number(c) + "},")
        out.append("};")
        return out

    lines = [
        "// The series for the CIO locator, s + XY/2 (IAU 2006/2000A; Capitaine, Wallace and Chapront",
        "// 2005, A&A 432, 355; IERS Conventions 2010, Table 5.2d), and the complementary terms of the",
        "// equation of the equinoxes (Capitaine and Gontier 1993; IERS Conventions 2010, Table 5.2e).",
        "// Coefficients in arcseconds; arguments are integer combinations of l, l', F, D, Ω, L_Ve, L_E, p_A.",
        "//",
        NOTICE.format(file="s06.c and eect00.c").rstrip(),
        "",
        "#pragma once",
        "",
        "namespace metrology::frames::series {",
        "",
        "struct ArgumentTerm {",
        "  int nfa[8];",
        "  double s, c;  // sine and cosine coefficients, arcseconds",
        "};",
        "",
        "// The polynomial part of s + XY/2, arcseconds, in powers 0..5 of t (Julian centuries TT).",
        "inline constexpr double s06_polynomial[6] = {" + ", ".join(sp) + "};",
        "",
    ]
    for name in ("s0", "s1", "s2", "s3", "s4"):
        lines += emit(f"s06_{name}", tables[name]) + [""]
    lines += emit("eect00_e0", e0) + [""] + emit("eect00_e1", e1)
    lines += ["", "}  // namespace metrology::frames::series", ""]
    return "\n".join(lines)


def tdb(src: Path) -> str:
    text = strip_comments((src / "dtdb.c").read_text(encoding="utf-8"))
    rows = groups(block(text, "fairhd[787][3] = {"))
    assert len(rows) == 787 and all(len(r) == 3 for r in rows)
    lines = [
        "// Fairhead and Bretagnon's series for TDB - TT at the geocentre: 787 terms A sin(ω T + φ) in",
        "// powers 0..4 of T, Julian millennia TDB from J2000.0 (Fairhead and Bretagnon 1990, A&A 229, 240;",
        "// the series as adjusted in the IERS Conventions 2010, §10.1). Rows 0-473 multiply T^0,",
        "// 474-678 T^1, 679-763 T^2, 764-783 T^3, 784-786 T^4. A in seconds, ω in rad per millennium.",
        "//",
        NOTICE.format(file="dtdb.c").rstrip(),
        "",
        "#pragma once",
        "",
        "namespace metrology::time::series {",
        "",
        "struct FairheadTerm {",
        "  double amplitude, frequency, phase;",
        "};",
        "",
        "inline constexpr FairheadTerm fairhead[787] = {",
    ]
    for r in rows:
        lines.append("    {" + ", ".join(number(x) for x in r) + "},")
    lines += ["};", "", "}  // namespace metrology::time::series", ""]
    return "\n".join(lines)


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    out = {
        ROOT / "domains/metrology/frames/include/metrology/frames/series_nutation.hpp": nutation(src),
        ROOT / "domains/metrology/frames/include/metrology/frames/series_cio.hpp": cio(src),
        ROOT / "domains/metrology/time/include/metrology/time/series_tdb.hpp": tdb(src),
    }
    for path, text in out.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
        print(f"wrote {path.relative_to(ROOT)} ({text.count(chr(10))} lines)")


if __name__ == "__main__":
    main()
