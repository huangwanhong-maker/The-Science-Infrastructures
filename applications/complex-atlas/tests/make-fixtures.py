"""Writes the page's fixtures (tests/fixtures/) from the engine's own answers.

The fixtures let the page be seen without the engine (?fixtures in the address). Each fixture
complex is named by its key, which is how the page finds its analysis, gauge run, state sum,
dimensions or spin foam; the generation of spacetime has its own keys (theories, expand-<theory>-
<order>, large-n-<theory>, condensate, random-geometry-<kind>, dimensions-<name>, spin-foam-<name>).

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
    generation(kind)
    print(f"Fixtures written to {FIXTURES.relative_to(APP.parents[1]).as_posix()}.")
    return 0


# The generation of spacetime (physis decision 0006): theories and their expansions, large-N
# phases, a condensate, random geometries with their dimensions, spin foams. Each operation the
# engine does not have yet is skipped, and its fixtures are left as they are.
EXPANSIONS = [("scalar", 1), ("scalar", 2), ("vector", 1), ("vector", 2), ("matrix", 1), ("matrix", 2), ("matrix", 3),
              ("tensor-melonic", 1), ("tensor-melonic", 2), ("tensor-wheel", 1), ("coloured-simplicial", 1),
              ("coloured-simplicial", 2), ("coloured-simplicial", 3), ("gft-boulatov", 1), ("gft-boulatov", 2), ("gft-ooguri", 1)]


def generation(kind: str) -> None:
    def ask(name, request):
        answer = engine.ask(kind, request)
        if "error" in answer:
            print(f"  skipped {name}: {answer['error']}")
            return None
        write(name, answer)
        return answer

    theories = ask("theories", {"op": "theories"})
    if theories:
        params = {t["id"]: {p["name"]: p["default"] for p in t.get("params", [])} for t in theories["theories"]}
        for theory, order in EXPANSIONS:
            ask(f"expand-{theory}-{order}", {"op": "expand", "theory": theory, "params": params.get(theory, {}), "order": order})
        for theory in ("matrix", "tensor-melonic", "tensor-wheel", "coloured-simplicial", "vector"):
            ask(f"large-n-{theory}", {"op": "large_n", "theory": theory, "params": params.get(theory, {})})
    ask("condensate", {"op": "condensate", "spins": [{"two_j": 1, "Q": 1, "E": 0}], "G": 1, "phi": [-1, 1], "points": 201})
    for geometry, size in (("quadrangulation", 400), ("melonic", 200)):
        answer = engine.ask(kind, {"op": "random_geometry", "kind": geometry, "size": size, "D": 3, "seed": 1})
        if "error" in answer:
            print(f"  skipped random-geometry-{geometry}: {answer['error']}")
            continue
        answer["complex"]["name"] = geometry
        write(f"random-geometry-{geometry}", answer)
        ask(f"dimensions-{geometry}", {"op": "dimensions", "complex": answer["complex"], "graph": "skeleton", "steps": 200, "starts": 64, "seed": 1})
    lattice = json.loads((FIXTURES / "build-lattice-torus.json").read_text(encoding="utf-8"))["complex"]
    ask("dimensions-lattice-torus", {"op": "dimensions", "complex": lattice, "graph": "skeleton", "steps": 100, "starts": 16, "seed": 1})
    simplex = json.loads((FIXTURES / "build-boundary-simplex.json").read_text(encoding="utf-8"))["complex"]
    melon = json.loads((FIXTURES / "build-melon.json").read_text(encoding="utf-8"))["complex"]
    for key, complex_ in (("boundary-simplex", simplex), ("simplex-4", {**simplex, "name": "simplex-4"}), ("melon", melon)):
        answer = ask(f"spin-foam-{key}", {"op": "spin_foam", "complex": complex_, "r": 5, "mode": "sample", "two_j": 2, "sweeps": 2000, "seed": 1})
        if answer:
            answer["foam"]["name"] = f"foam-{key}"
            write(f"spin-foam-{key}", answer)


if __name__ == "__main__":
    sys.exit(main())
