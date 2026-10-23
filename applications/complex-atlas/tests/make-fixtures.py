"""Writes the page's fixtures (tests/fixtures/) from the engine's own answers.

The fixtures let the page be seen without the engine (?fixtures in the address). Each fixture
complex is named by its key, which is how the page finds its analysis, gauge run or state sum.

Usage::

    python applications/complex-atlas/engine.py build
    python applications/complex-atlas/tests/make-fixtures.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))
import engine  # noqa: E402

FIXTURES = APP / "tests" / "fixtures"
# key: (family, parameters)
COMPLEXES = {
    "torus": ("words", {"words": "abAB"}),
    "surface": ("surface", {"genus": 2}),
    "projective-plane": ("nonorientable", {"crosscaps": 1}),
    "lattice-torus": ("lattice-torus", {"n": 4, "m": 4}),
    "boundary-simplex": ("boundary-simplex", {"n": 4}),
    "melon": ("melon", {"D": 3, "insertions": 2, "seed": 1}),
}


def write(name: str, answer: dict) -> None:
    if "error" in answer:
        raise SystemExit(f"{name}: {answer['error']}")
    (FIXTURES / f"{name}.json").write_text(json.dumps(answer, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


def main() -> int:
    st = engine.status()
    if not st["current"]:
        print("The engine is not built and current: python applications/complex-atlas/engine.py build")
        return 1
    kind = st["toolchain"]
    FIXTURES.mkdir(parents=True, exist_ok=True)
    write("catalogue", engine.ask(kind, {"op": "catalogue"}))
    for key, (family, params) in COMPLEXES.items():
        complex_ = engine.ask(kind, {"op": "build", "family": family, "params": params})["complex"]
        complex_["name"] = key
        write(f"build-{key}", {"complex": complex_})
        analysis = engine.ask(kind, {"op": "analyse", "complex": complex_})
        write(f"analyse-{key}", analysis)
        if key == "lattice-torus":
            loops = analysis["analysis"]["fundamental_group"]["simplified"]["loops"]
            write(f"gauge-{key}", engine.ask(kind, {"op": "gauge", "complex": complex_, "group": {"kind": "U1"}, "beta": 2.0,
                                                   "sweeps": 400, "thermalise": 100, "seed": 1, "update": "heatbath", "loops": loops}))
        if key in ("boundary-simplex", "melon"):
            write(f"turaev-viro-{key}", engine.ask(kind, {"op": "turaev_viro", "complex": complex_, "r": 4}))
    print(f"Fixtures written to {FIXTURES.relative_to(APP.parents[1]).as_posix()}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
