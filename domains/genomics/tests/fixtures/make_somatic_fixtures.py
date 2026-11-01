"""Reference values for the checks of the somatic module and of sequence/contexts.hpp.

Usage (from any directory; Python 3.12 with scipy):

    python make_somatic_fixtures.py <path to sigprofilermatrixgenerator-1.4.0-py3-none-any.whl>

Writes, beside this script:

- somatic-classes.txt: the classes of mutational catalogues as SigProfilerMatrixGenerator 1.4.0
  (BSD-2-Clause, https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator) defines them. The
  wheel is read, never installed or imported: its SigProfilerMatrixGeneratorFunc.py is parsed and
  the lists `mutation_types` and `mutation_types_non_tsb` (the 78 doublet classes) and
  `indel_types` (the 83 indel classes, followed by its 'complex' and 'non_matching') are taken as
  written. The DBS78 and SBS96 matrices of MutationMatrixGenerator.py are pandas groupbys (sorted
  keys), so the doublets are listed sorted; the 96 classes are rebuilt as
  SigProfilerMatrixGeneratorFunc.py builds its 1536 pentanucleotide classes (perm(5, "ACGT"), a
  pyrimidine middle base, the three other bases) and reduced by name[1:8] as MutationMatrixGenerator
  reduces them, sorted.
- somatic-chi-square.txt: Pearson's chi-square of 2 x K tables of repeat lengths (normal, tumour),
  computed by scipy.stats.chi2_contingency(correction=False) on the columns where either count is
  positive: statistic, degrees of freedom and p-value. Tables drawn from a fixed numpy seed.
"""

import ast
import sys
import zipfile
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import chi2_contingency

HERE = Path(__file__).resolve().parent


def classes(wheel: Path) -> list[str]:
    with zipfile.ZipFile(wheel) as z:
        src = z.read("SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGeneratorFunc.py").decode("utf8")
    lists: dict[str, list[str]] = {}
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in ("mutation_types", "mutation_types_non_tsb", "indel_types") and isinstance(node.value, ast.List):
                values = [e.value for e in node.value.elts if isinstance(e, ast.Constant)]
                if name not in lists:
                    lists[name] = values
    dbs = sorted(lists["mutation_types"] + lists["mutation_types_non_tsb"])
    assert len(dbs) == 78, len(dbs)
    ids = [x for x in lists["indel_types"] if x not in ("complex", "non_matching")]
    # The 83 of the ID83 matrix (the 'Ins:M' classes of the list are SigProfiler's own, not ID83's).
    id83 = [x for x in ids if ":Ins:M:" not in x]
    assert len(id83) == 83, len(id83)
    bases = "ACGT"
    five = ["".join((a, b, c, d, e)) for a in bases for b in bases for c in bases for d in bases for e in bases]
    sbs1536 = []
    for mut in five:
        centre = mut[2]
        if centre in "CT":
            for base in bases:
                if base != centre:
                    sbs1536.append(mut[0:2] + "[" + centre + ">" + base + "]" + mut[3:])
    sbs96 = sorted({m[1:8] for m in sbs1536})
    assert len(sbs96) == 96
    return [f"SBS96\t{x}" for x in sbs96] + [f"DBS78\t{x}" for x in dbs] + [f"ID83\t{x}" for x in id83]


def chi_square_tables() -> list[str]:
    rng = np.random.default_rng(20261101)
    out = []
    for t in range(40):
        K = int(rng.integers(2, 9))
        centre = int(rng.integers(8, 25))
        lengths = list(range(centre - K // 2, centre - K // 2 + K))
        normal = rng.multinomial(int(rng.integers(15, 120)), rng.dirichlet(np.ones(K) * 2))
        shift = rng.dirichlet(np.ones(K) * (2 if t % 3 else 20))
        tumour = rng.multinomial(int(rng.integers(15, 120)), shift)
        cols = [i for i in range(K) if normal[i] + tumour[i] > 0]
        if len(cols) < 2:
            continue
        table = np.array([[normal[i] for i in cols], [tumour[i] for i in cols]])
        stat, p, dof, _ = chi2_contingency(table, correction=False)
        fmt = lambda v: ",".join(f"{lengths[i]}:{int(v[i])}" for i in range(K) if v[i] > 0)
        out.append(f"{fmt(normal)}\t{fmt(tumour)}\t{stat:.17g}\t{dof}\t{p:.17g}")
    return out


def main() -> None:
    wheel = Path(sys.argv[1])
    lines = classes(wheel)
    (HERE / "somatic-classes.txt").write_text(
        "# SigProfilerMatrixGenerator 1.4.0 (" + wheel.name + "), parsed by make_somatic_fixtures.py; retrieved 2026-11-01\n"
        + "\n".join(lines) + "\n", encoding="utf8", newline="\n")
    tables = chi_square_tables()
    (HERE / "somatic-chi-square.txt").write_text(
        f"# normal\ttumour\tstatistic\tdof\tp; scipy {scipy.__version__} chi2_contingency(correction=False), numpy {np.__version__}; 2026-11-01\n"
        + "\n".join(tables) + "\n", encoding="utf8", newline="\n")
    print(len(lines), "classes,", len(tables), "tables")


if __name__ == "__main__":
    main()
