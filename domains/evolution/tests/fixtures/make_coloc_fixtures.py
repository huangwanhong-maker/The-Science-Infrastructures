"""Reference values for the association module's coloc.hpp, from R's coloc package.

Writes coloc.txt: three seeded scenarios of two traits' estimates at 60 variants (beta and its variance):
a shared causal variant (`shared`), two distinct ones (`distinct`) and the second trait null (`one`);
and coloc.abf's posteriors PP.H0 … PP.H4 and SNP.PP.H4, the first trait quantitative (sdY = 1, so the
prior sd is 0.15) and the second case–control (prior sd 0.2), with coloc's default priors.

    RSCRIPT=<path to Rscript with coloc> python make_coloc_fixtures.py

Run 2026-11-01 with R 4.4.3 (conda-forge), coloc 5.2.3, numpy 2.5.4.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def run_r(script):
    rscript = os.environ.get("RSCRIPT")
    if not rscript:
        sys.exit("set RSCRIPT to an Rscript that has coloc")
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "script.R"
        path.write_text(script, encoding="utf-8")
        out = subprocess.run([rscript, str(path)], capture_output=True, text=True, encoding="utf-8", errors="replace")
        if out.returncode:
            sys.exit(out.stderr)
        return out.stdout


def scenario(rng, p, c1, c2):
    # Correlated z-scores: AR(1) LD, ρ = 0.8, signals at c1 and c2 (None for none).
    R = 0.8 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
    L = np.linalg.cholesky(R)
    def z(c, size):
        mean = np.zeros(p) if c is None else size * R[:, c]
        return mean + L @ rng.normal(size=p)
    v1 = rng.uniform(0.002, 0.004, size=p)
    v2 = rng.uniform(0.004, 0.008, size=p)
    return z(c1, 7.0) * np.sqrt(v1), v1, z(c2, 6.0) * np.sqrt(v2), v2


def main():
    rng = np.random.default_rng(2014)
    p = 60
    cases = {"shared": (20, 20), "distinct": (12, 45), "one": (30, None)}
    fmt = lambda a: ",".join(f"{v:.17g}" for v in a)
    lines = ["# evolution association: coloc.abf on seeded scenarios (see make_coloc_fixtures.py), retrieved 2026-11-01."]
    script = "suppressMessages(library(coloc))\nw <- function(name, v) cat(name, ':', sprintf('%.17g', as.numeric(v)), '\n')\n"
    for name, (c1, c2) in cases.items():
        b1, v1, b2, v2 = scenario(rng, p, c1, c2)
        for key, vec in (("b1", b1), ("v1", v1), ("b2", b2), ("v2", v2)):
            lines.append(f"{name}_{key} : " + " ".join(f"{v:.17g}" for v in vec))
        script += f"""
snp <- paste0("s", 1:{p})
d1 <- list(beta = c({fmt(b1)}), varbeta = c({fmt(v1)}), type = "quant", sdY = 1, snp = snp)
d2 <- list(beta = c({fmt(b2)}), varbeta = c({fmt(v2)}), type = "cc", snp = snp)
r <- suppressMessages(coloc.abf(d1, d2))
w("{name}_pp", r$summary[2:6]); w("{name}_snp_pp_h4", r$results$SNP.PP.H4[match(snp, r$results$snp)])
"""
    script += "cat('versions :', R.version.string, as.character(packageVersion('coloc')), '\n')\n"
    out = run_r(script)
    lines += [l.strip() for l in out.splitlines() if l.strip()]
    (HERE / "coloc.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
