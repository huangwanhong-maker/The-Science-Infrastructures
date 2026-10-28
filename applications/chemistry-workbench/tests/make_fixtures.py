"""Reference values for the Chemistry Workbench's checks, made with RDKit (BSD), not typed.

The CIP descriptors of the two configurations of CHBrClF that the plug-in Symmetry and generation
builds: the engine's coordinates are read (engine.py ask), written as a molecule with 3D coordinates,
and RDKit assigns stereochemistry from them (AssignStereochemistryFrom3D, then its new CIP labeller,
rdCIPLabeler, Hanson et al. 2018, J. Chem. Inf. Model. 58, 1755). The test compares the engine's own
descriptor and its sign τ with these.

Usage (with RDKit installed, e.g. pip install rdkit into a scratch environment)::

    python applications/chemistry-workbench/tests/make_fixtures.py
"""

from __future__ import annotations

import datetime
import json
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler

APP = Path(__file__).resolve().parents[1]
ROOT = APP.parents[1]


def engine(request: dict) -> dict:
    r = subprocess.run([sys.executable, str(APP / "engine.py"), "ask", "symmetry-generation", json.dumps(request)],
                       capture_output=True, text=True, encoding="utf-8", cwd=ROOT)
    answer = json.loads(r.stdout)
    if not answer.get("ok"):
        raise RuntimeError(answer)
    return answer["result"]


def descriptor(molecule: dict) -> str:
    rw = Chem.RWMol()
    for atom in molecule["atoms"]:
        rw.AddAtom(Chem.Atom(atom["element"]))
    for a, b, order in molecule["bonds"]:
        rw.AddBond(int(a), int(b), Chem.BondType.SINGLE if order == 1 else Chem.BondType.DOUBLE)
    conf = Chem.Conformer(rw.GetNumAtoms())
    for i, atom in enumerate(molecule["atoms"]):
        conf.SetAtomPosition(i, atom["x"])
    mol = rw.GetMol()
    mol.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    rdCIPLabeler.AssignCIPLabels(mol)
    centre = mol.GetAtomWithIdx(0)
    return centre.GetProp("_CIPCode") if centre.HasProp("_CIPCode") else ""


def main() -> None:
    geometry = engine({"op": "grammar", "example": "stereocentre"})["geometry"]
    out = {
        "source": f"RDKit {rdBase.rdkitVersion}: AssignStereochemistryFrom3D and rdCIPLabeler on the engine's coordinates",
        "made": datetime.date.today().isoformat(),
        "plus": {"tau": geometry["tauPlus"], "cip": descriptor(geometry["plus"])},
        "minus": {"tau": geometry["tauMinus"], "cip": descriptor(geometry["minus"])},
    }
    target = APP / "tests" / "fixtures" / "cip.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
