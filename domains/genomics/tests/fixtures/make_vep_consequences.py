"""Ensembl VEP's consequence terms and their order of severity, for genomics/nomenclature's checks.

Reads modules/Bio/EnsEMBL/Variation/Utils/Constants.pm of ensembl-variation at release/116 (the
OVERLAP_CONSEQUENCES table: each term's SO accession, impact and rank) and writes
vep-consequences.tsv: term, SO accession, impact, rank, in the order of rank.

Usage: python make_vep_consequences.py   (network: raw.githubusercontent.com)
"""

from __future__ import annotations

import datetime as dt
import re
import sys
import urllib.request
from pathlib import Path

RELEASE = "release/116"
URL = f"https://raw.githubusercontent.com/Ensembl/ensembl-variation/{RELEASE}/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm"


def main() -> int:
    with urllib.request.urlopen(URL, timeout=60) as r:
        text = r.read().decode("utf-8")
    blocks = re.findall(r"'(\w+)' => Bio::EnsEMBL::Variation::OverlapConsequence->new_fast\(\{(.*?)\}\n\)", text, re.S)
    # The default consequence (intergenic_variant) is defined apart, as $DEFAULT_OVERLAP_CONSEQUENCE.
    default = re.search(r"\$DEFAULT_OVERLAP_CONSEQUENCE = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast\(\{(.*?)\}\n\)", text, re.S)
    if default:
        term = re.search(r"'SO_term' => '([^']+)'", default.group(1))
        blocks.append((term.group(1), default.group(1)))
    rows = []
    for name, body in blocks:
        acc = re.search(r"'SO_accession' => '([^']+)'", body)
        impact = re.search(r"'impact' => '([^']+)'", body)
        rank = re.search(r"'rank' => '(\d+)'", body)
        if acc and impact and rank:
            rows.append((name, acc.group(1), impact.group(1), int(rank.group(1))))
    rows.sort(key=lambda r: r[3])
    here = Path(__file__).resolve().parent
    out = here / "vep-consequences.tsv"
    with out.open("w", encoding="utf-8", newline="\n") as f:
        f.write(f"# Ensembl VEP's consequences, from ensembl-variation {RELEASE} Constants.pm ({URL}), retrieved {dt.date.today().isoformat()} by make_vep_consequences.py\n")
        f.write("# columns: SO term, SO accession, impact, rank (1 the most severe)\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print(f"{out.name}: {len(rows)} terms")
    return 0


if __name__ == "__main__":
    sys.exit(main())
