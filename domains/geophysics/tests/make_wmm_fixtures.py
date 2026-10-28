"""Reference values for the geophysics domain's geomagnetism (geomagnetism/tests/check_wmm.cpp).

Reads NOAA's WMM2025 files from the data library's cache (python data/acquisition/fetch.py
noaa-wmm), checks each against the SHA-256 the registry records, and writes into tests/fixtures/:

- wmm2025-test-values.txt: the 100 test points distributed with the coefficients
  (WMM2025COF.zip, member WMM2025COF/WMM2025_TestValues.txt), unaltered below a header line;
- wmm2025-table6.txt: the 12 test points of the report's Table 6, with the grid variation
  (NCEI's WMM2025_TEST_VALUES.txt), unaltered below a header line;
- wmm2025-report.txt: values read out of the technical report's text (Chulliat et al. 2025,
  WMM_Report_2025.pdf) as "key value unit" lines: the high-precision example of Tables 3a and
  3b, the grid variations of Table 3d, the field above the poles (section 1.4), the
  geomagnetic poles (section 1.6) and the construction of the declination's error model
  (section 3.4). The report prints its digits in groups of five; the groups
  are joined, and the minus sign U+2212 read as '-'. Nothing is typed by hand: every number is
  matched in the report's text, and the script stops if one is not found.

Deterministic; run with an interpreter that has pypdf (the reference environment):

    python -I domains/geophysics/tests/make_wmm_fixtures.py
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / "data" / "cache" / "noaa-wmm"
REGISTRY = ROOT / "data" / "registry" / "noaa-wmm.json"
OUT = HERE / "fixtures"


def verified(name: str, manifest: dict) -> Path:
    path = CACHE / name
    if not path.exists():
        raise SystemExit(f"{path} is missing: python data/acquisition/fetch.py noaa-wmm")
    entry = next(e for e in manifest["files"] if e["file"] == name)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != entry["sha256"]:
        raise SystemExit(f"{name}: SHA-256 {digest} is not the registry's {entry['sha256']}")
    return path


def source(name: str, manifest: dict, retrieved: str) -> str:
    entry = next(e for e in manifest["files"] if e["file"] == name)
    where = entry["url"] + (f", member {entry['member']}" if "member" in entry else "")
    return f"NOAA NCEI, {where} (SHA-256 {entry['sha256'][:12]}...), retrieved {retrieved} through the data library (noaa-wmm)"


def number(text: str) -> str:
    """A number as the report prints it, its digit groups joined: '64 57402.34844 73705' -> '6457402.3484473705'."""
    return text.replace(" ", "")


def main() -> None:
    manifest = json.loads(REGISTRY.read_text(encoding="utf-8"))
    retrieved = json.loads((CACHE / "retrieval.json").read_text(encoding="utf-8"))["retrieved_at"][:10]
    OUT.mkdir(exist_ok=True)

    # ---- the two files of test values, unaltered ---------------------------------------------------
    for name, out, what in (("WMM2025_TestValues.txt", "wmm2025-test-values.txt", "WMM2025's 100 test points, distributed with its coefficients"),
                            ("WMM2025_TEST_VALUES.txt", "wmm2025-table6.txt", "WMM2025's test values of the technical report's Table 6")):
        body = verified(name, manifest).read_bytes().decode("ascii")
        header = (f"# {what}: {source(name, manifest, retrieved)}; copied unaltered below this line and the next by "
                  f"domains/geophysics/tests/make_wmm_fixtures.py.\n# Heights are above the WGS 84 ellipsoid; times are decimal years.\n")
        (OUT / out).write_text(header + body, encoding="ascii", newline="")
        print(f"wrote {out}")

    # ---- the report ---------------------------------------------------------------------------------
    reader = PdfReader(str(verified("WMM_Report_2025.pdf", manifest)))
    page = lambda printed: reader.pages[printed - 1].extract_text().replace("−", "-")  # noqa: E731
    flat = lambda printed: re.sub(r"\s+", " ", page(printed))  # noqa: E731
    rows: list[tuple[str, str, str]] = []

    def find(pattern: str, text: str, what: str) -> re.Match:
        m = re.search(pattern, text)
        if not m:
            raise SystemExit(f"the report's {what} was not found")
        return m

    # Table 3a, the example's place and time (page 17).
    p17 = page(17)
    for key, label, unit in (("time", "Time", "yr"), ("height", "Height-above-Ellipsoid", "km"), ("latitude", "Latitude", "deg"), ("longitude", "Longitude", "deg")):
        m = find(rf"(?m)^{label} (-?[\d ]+\.[\d ]+?) {unit}\s*$", p17, f"Table 3a row {label}")
        rows.append((f"table3a.{key}", number(m.group(1)), unit))
    # Table 3b, the computation for WMM2025 (page 18): 36 rows "n name value unit".
    found = 0
    for line in page(18).splitlines():
        m = re.match(r"^(\d{1,2}) (\S+) (-?[\d ]+\.[\d ]+?)\s+(nT/yr|nT|rad/yr|rad|m|yr)\s*$", line.strip())
        if m:
            rows.append((f"table3b.{m.group(2)}", number(m.group(3)), m.group(4)))
            found += 1
    if found != 36:
        raise SystemExit(f"Table 3b: {found} rows read, not 36")
    # Table 3d, the grid variations for WMM2025 (page 20), by grid: UPS north, UPS south, UTM 10, UTM 11.
    p20 = flat(20)
    t3d = find(r"Table 3d:(.*?)Table 3e:", p20, "Table 3d").group(1)
    for label in ("TrueN-to-GridN", "GridN-to-MagN", "TrueN-to-MagN"):
        values = re.findall(r"-?\d+\.\d{5} \d{5}", find(rf"{label} ((?:-?\d+\.\d{{5}} \d{{5}} ?){{4}})", t3d, f"Table 3d row {label}").group(1))
        for grid, v in zip(("UPS-north", "UPS-south", "UTM-10", "UTM-11"), values, strict=True):
            rows.append((f"table3d.{label}.{grid}", number(v), "deg"))
    # Section 1.4, the field above the poles on 1 January 2025 at r = 6 371 200 m, the pole given λ = 0.
    p17f = flat(17)
    m = find(r"On January 1, (\d{4}), directly above the North \(resp\. South\) Pole at ([\d,]+) meters", p17f, "section 1.4's instant and radius")
    rows.append(("poles.year", m.group(1), "yr"))
    rows.append(("poles.radius", m.group(2).replace(",", ""), "m"))
    m = find(r"predicted by the WMM lies in the half-plane of the ([\d.]+)°([EW]) \(resp\. ([\d.]+)°([EW])\) meridian", p17f, "section 1.4's meridians")
    rows.append(("poles.north.meridian", ("-" if m.group(2) == "W" else "") + m.group(1), "deg"))
    rows.append(("poles.south.meridian", ("-" if m.group(4) == "W" else "") + m.group(3), "deg"))
    m = find(r"predicted by the WMM are (-?[\d.]+) nT, (-?[\d.]+) nT, and (-?[\d.]+) nT respectively at the North Pole, (-?[\d.]+) nT, (-?[\d.]+) nT, and (-?[\d.]+) nT respectively at the South Pole",
             p17f, "section 1.4's components")
    for i, key in enumerate(("north.Xprime", "north.Yprime", "north.Zprime", "south.Xprime", "south.Yprime", "south.Zprime")):
        rows.append((f"poles.{key}", m.group(i + 1), "nT"))
    # Section 1.6, the geomagnetic poles of the centred dipole at 2025.0.
    m = find(r"Based on the WMM2025 coefficients for ([\d.]+) the geomagnetic north pole is at ([\d.]+)°([EW]) longitude and ([\d.]+)°([NS]) geocentric latitude \(([\d.]+)°([NS]) geodetic latitude\), "
             r"and the geomagnetic south pole is at ([\d.]+)°([EW]) longitude and ([\d.]+)°([NS]) geocentric latitude \(([\d.]+)°([NS]) geodetic latitude\)", p20, "section 1.6's geomagnetic poles")
    sign = lambda v, h: ("-" if h in "WS" else "") + v  # noqa: E731
    rows.append(("dipole.year", m.group(1), "yr"))
    rows.append(("dipole.north.longitude", sign(m.group(2), m.group(3)), "deg"))
    rows.append(("dipole.north.geocentric_latitude", sign(m.group(4), m.group(5)), "deg"))
    rows.append(("dipole.north.geodetic_latitude", sign(m.group(6), m.group(7)), "deg"))
    rows.append(("dipole.south.longitude", sign(m.group(8), m.group(9)), "deg"))
    rows.append(("dipole.south.geocentric_latitude", sign(m.group(10), m.group(11)), "deg"))
    rows.append(("dipole.south.geodetic_latitude", sign(m.group(12), m.group(13)), "deg"))
    m = find(r"The axis of the dipole is currently inclined at ([\d.]+)° to the Earth", p20, "section 1.6's tilt")
    rows.append(("dipole.tilt", m.group(1), "deg"))
    # Section 3.4, how the error model of D is built: its value where H is greatest on the surface.
    m = find(r"so that δD = ([\d.]+)° where H = (\d+) nT, the maximum value of H at the Earth.s surface in WMM2025", flat(63), "section 3.4's declination error")
    rows.append(("uncertainty.declination", m.group(1), "deg"))
    rows.append(("uncertainty.horizontal", m.group(2), "nT"))

    header = [f"# Values read out of the WMM2025 technical report (Chulliat et al. 2025, doi:10.25923/prbc-s316): {source('WMM_Report_2025.pdf', manifest, retrieved)};",
              "# written by domains/geophysics/tests/make_wmm_fixtures.py. key value unit; table3a/3b: the high-precision example (pp. 17-18), table3d: its grid variations (p. 20),",
              "# poles: the field above the geographic poles, the pole given longitude 0 (section 1.4, p. 17), dipole: the geomagnetic poles (section 1.6, p. 20),",
              "# uncertainty: the declination error where H is greatest on the surface (section 3.4, p. 63). W and S are negative."]
    (OUT / "wmm2025-report.txt").write_text("\n".join(header + [f"{k} {v} {u}" for k, v, u in rows]) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote wmm2025-report.txt: {len(rows)} values")


if __name__ == "__main__":
    main()
