"""The fixtures of the colour-vision checks (visualization/styles/cvd.js), made by a reference library.

colorspacious 1.1.2 (N. J. Smith and S. van der Walt; MIT licence) implements the simulation of Machado,
Oliveira and Fernandes (2009), "A physiologically-based model for simulation of color vision deficiency",
IEEE TVCG 15(6), 1291-1298: a matrix on linear sRGB for each deficiency and severity. This script writes
the three matrices at severity 1 (protanopia, deuteranopia, tritanopia) and the simulated sRGB of a set of
colours, unclamped, as colorspacious's cspace_convert gives them, to tests/fixtures/cvd.json.

Run from the repository root (colorspacious and numpy installed):
    python visualization/tests/make_cvd_fixtures.py
"""

import datetime
import json
from pathlib import Path

import colorspacious
from colorspacious import cspace_convert
from colorspacious.cvd import machado_et_al_2009_matrix

# The infrastructure's colours and the primaries: the applications' main colours, the palettes' stops.
COLOURS = [
    "#ca8dbe", "#ee8fbd", "#b892dd", "#c2477f", "#7d54ab", "#c9a85c", "#c4506e", "#a8577e", "#8e7cc3",
    "#2a1f3a", "#6b4a86", "#c08081", "#f0d9b5", "#5d5aa0", "#a99ad0", "#f2e6ee", "#d38e8a", "#9c4a45",
    "#f6e3b4", "#e9aa40", "#965208", "#ffffff", "#000000", "#ff0000", "#00ff00", "#0000ff", "#808080",
]
KINDS = {"protan": "protanomaly", "deutan": "deuteranomaly", "tritan": "tritanomaly"}


def rgb(hexcode: str) -> list[float]:
    return [int(hexcode[i:i + 2], 16) / 255 for i in (1, 3, 5)]


def main() -> None:
    out = {
        "source": f"colorspacious {colorspacious.__version__}: machado_et_al_2009_matrix(kind, 100), and cspace_convert from "
                  "{'name': 'sRGB1+CVD', 'cvd_type': kind, 'severity': 100} to 'sRGB1' (Machado, Oliveira and Fernandes 2009)",
        "retrieved": datetime.date.today().isoformat(),
        "matrices": {k: machado_et_al_2009_matrix(v, 100).tolist() for k, v in KINDS.items()},
        "colours": [],
    }
    for h in COLOURS:
        row = {"hex": h}
        for k, v in KINDS.items():
            row[k] = [float(x) for x in cspace_convert(rgb(h), {"name": "sRGB1+CVD", "cvd_type": v, "severity": 100}, "sRGB1")]
        out["colours"].append(row)
    path = Path(__file__).resolve().parent / "fixtures" / "cvd.json"
    path.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {path} ({len(COLOURS)} colours)")


if __name__ == "__main__":
    main()
