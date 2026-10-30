"""Reference values for the geophysics domain's PREM (structure/tests/check_structure.cpp).

Fetches ObsPy's tabulation of PREM (obspy/taup/data/prem.nd, the isotropic model at named
discontinuities, as TauP reads it) at the commit that last changed it, checks its SHA-256, and writes
tests/fixtures/prem-obspy.txt: one line per row of the table, "depth_km vp_km_s vs_km_s rho_g_cm3"
(the Q columns left out), with the table's own section names (mantle, outer-core, inner-core) as
comment lines, unaltered below a header that states the source. Two rows of one depth are the two
sides of a discontinuity, the upper first. ObsPy's table carries PREM's crust to the surface (no
ocean): its first rows are the upper crust from 0 km.

ObsPy is LGPL-3.0; the numbers are PREM's (Dziewonski and Anderson 1981), as ObsPy's table gives
them, kept here as reference values with their source.

Deterministic; the standard library only:

    python -I domains/geophysics/tests/make_prem_fixtures.py
"""

from __future__ import annotations

import datetime
import hashlib
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "fixtures" / "prem-obspy.txt"
COMMIT = "92fd0283a2a4af909d072afc774d45663187ed0c"  # 2015-12-02, the last change to prem.nd
URL = f"https://raw.githubusercontent.com/obspy/obspy/{COMMIT}/obspy/taup/data/prem.nd"
SHA256 = "efd8960973133957f14fad0f98e993e0e9236b045d20dda6be61e4bb73917f45"


def main() -> int:
    with urllib.request.urlopen(URL, timeout=60) as r:
        data = r.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != SHA256:
        raise SystemExit(f"prem.nd has SHA-256 {digest}, not {SHA256}")
    lines = [
        "# PREM (Dziewonski and Anderson 1981), isotropic, as ObsPy's TauP tabulates it (obspy/taup/data/prem.nd).",
        f"# Source: {URL}",
        f"# SHA-256 {SHA256}; retrieved {datetime.date.today().isoformat()} by tests/make_prem_fixtures.py.",
        "# Columns: depth (km), vp (km/s), vs (km/s), density (g/cm3); the table's Q columns are left out.",
        "# Two rows of one depth are the two sides of a discontinuity, the upper first. The table's crust reaches the surface (no ocean).",
    ]
    for raw in data.decode("ascii").splitlines():
        parts = raw.split()
        if not parts:
            continue
        if len(parts) == 1:
            lines.append(f"# {parts[0]}")
            continue
        lines.append(" ".join(parts[:4]))
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"{OUT.relative_to(HERE.parents[2])}: {sum(1 for x in lines if not x.startswith('#'))} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
