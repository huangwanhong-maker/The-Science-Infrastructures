"""Reference values for the structural module's genotypes: SVTyper's own bayes_gt run on a grid of supports.

SVTyper (Chiang et al. 2015, Nature Methods 12, 966-968; MIT licence, copyright 2014 Colby Chiang) is
Python 2. Its functions log_choose and bayes_gt (svtyper/statistics.py at commit
3c5218a1ea148dc736126bf0019fb8e2a6828f9e, https://github.com/hall-lab/svtyper, retrieved 2026-11-01)
are copied below verbatim but for ``xrange`` written ``range``. For each (ref, alt, is_dup) the three
log10 likelihoods bayes_gt returns are written, with the genotype, GQ and SQ that singlesample.py's
bayesian_genotype derives from them (its lines quoted in ``qualities``).

Usage::

    python domains/genomics/tests/fixtures/make_structural_fixtures.py

writes structural-svtyper.tsv beside this script.
"""

from __future__ import annotations

import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMMIT = "3c5218a1ea148dc736126bf0019fb8e2a6828f9e"


# ---- svtyper/statistics.py at 3c5218a (xrange -> range) ----------------------------------------
# efficient combinatorial function to handle extremely large numbers
def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in range(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# return the genotype and log10 p-value
def bayes_gt(ref, alt, is_dup):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [1e-2, 0.2, 1/3.0]
    else:
        p_alt = [1e-3, 0.5, 0.9]

    total = ref + alt
    log_combo = log_choose(total, alt)

    lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return (lp_homref, lp_het, lp_homalt)
# ---- end of the copy ---------------------------------------------------------------------------


def qualities(gt_lplist):
    # singlesample.py, bayesian_genotype:
    #   best, second_best = sorted([(i, e) for i, e in enumerate(gt_lplist)], key=lambda(x): x[1], reverse=True)[0:2]
    #   sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log))
    #   phred_gq = min(-10 * (second_best[1] - best[1]), 200)
    ranked = sorted([(i, e) for i, e in enumerate(gt_lplist)], key=lambda x: x[1], reverse=True)[0:2]
    best, second = ranked
    gt_sum = sum(10 ** g for g in gt_lplist)
    sq = abs(-10 * (gt_lplist[0] - math.log(gt_sum, 10)))
    gq = int(min(-10 * (second[1] - best[1]), 200))
    return best[0], gq, sq


def main() -> None:
    rows = ["# SVTyper bayes_gt (svtyper/statistics.py at %s) and bayesian_genotype's GT, GQ, SQ" % COMMIT,
            "# ref\talt\tis_dup\tlp_homref\tlp_het\tlp_homalt\tgt\tgq\tsq"]
    for is_dup in (False, True):
        for ref, alt in [(0, 5), (5, 0), (10, 10), (20, 3), (3, 20), (15, 8), (8, 15), (30, 30), (1, 1), (40, 2), (2, 40), (12, 25), (60, 31)]:
            lp = bayes_gt(ref, alt, is_dup)
            gt, gq, sq = qualities(lp)
            rows.append("\t".join([str(ref), str(alt), "1" if is_dup else "0"] + ["%.12g" % x for x in lp] + [str(gt), str(gq), "%.12g" % sq]))
    (HERE / "structural-svtyper.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {len(rows) - 2} rows")


if __name__ == "__main__":
    main()
