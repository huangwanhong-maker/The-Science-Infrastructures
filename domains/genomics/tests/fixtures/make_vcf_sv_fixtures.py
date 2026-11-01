"""The structural-variant examples of the VCF 4.3 specification, for genomics/nomenclature's checks.

Reads VCFv4.3.tex of samtools/hts-specs at a pinned commit and writes
  vcf43-structural.vcf   the "VCF STRUCTURAL VARIANT EXAMPLE" of §5.3, verbatim, its columns joined
                         by tabs (the specification aligns them with spaces);
  vcf43-breakends.tsv    every breakend table of §5.4 (its subsection, then CHROM POS ID REF ALT QUAL
                         FILTER INFO), the LaTeX ($, \\_, \\[, $<$) removed.

Usage: python make_vcf_sv_fixtures.py   (network: raw.githubusercontent.com)
"""

from __future__ import annotations

import datetime as dt
import re
import sys
import urllib.request
from pathlib import Path

COMMIT = "510c107ba89aa13113818c9e4fa3071283cb7f7b"  # hts-specs master, 2026-08-11
URL = f"https://raw.githubusercontent.com/samtools/hts-specs/{COMMIT}/VCFv4.3.tex"


def untex(cell: str) -> str:
    s = cell.strip()
    s = s.replace("$<$", "<").replace("$>$", ">").replace("\\_", "_").replace("$", "")
    s = s.replace("\\[", "[").replace("\\]", "]")
    s = re.sub(r"\\textbf\{([^}]*)\}", r"\1", s)
    return s.strip()


def main() -> int:
    with urllib.request.urlopen(URL, timeout=60) as r:
        tex = r.read().decode("utf-8")
    here = Path(__file__).resolve().parent
    today = dt.date.today().isoformat()
    # §5.3: the verbatim example.
    m = re.search(r"VCF STRUCTURAL VARIANT EXAMPLE\n\n(.*?)\\end\{verbatim\}", tex, re.S)
    lines = [l for l in m.group(1).splitlines() if l.strip()]
    vcf = []
    for l in lines:
        if l.startswith("##"):
            vcf.append(l.strip())
        else:
            cols = l.split()
            vcf.append("\t".join(cols))
    (here / "vcf43-structural.vcf").write_text("\n".join(vcf) + "\n", encoding="utf-8", newline="\n")
    # §5.4: the breakend tables.
    start = tex.index("\\subsection{Specifying complex rearrangements with breakends}")
    end = tex.index("\\subsection{Representing unspecified alleles")
    section = tex[start:end]
    rows = []
    current = "Breakends"
    for part in re.split(r"(\\subsubsection\{[^}]*\})", section):
        h = re.match(r"\\subsubsection\{([^}]*)\}", part)
        if h:
            current = h.group(1)
            continue
        for table in re.findall(r"\\begin\{tabular\}.*?\\end\{tabular\}", part, re.S):
            for row in table.split("\\\\"):
                cells = [untex(c) for c in row.split("&")]
                cells = [c for c in cells if c != ""] if len(cells) > 8 else cells
                if len(cells) != 8 or cells[0].startswith("\\#") or "CHROM" in cells[0] or "begin{tabular}" in cells[0]:
                    continue
                cells[0] = re.sub(r"^.*\}\s*", "", cells[0])
                if not re.match(r"^<?[A-Za-z0-9_]+>?$", cells[0]) or not cells[1].isdigit():
                    continue
                rows.append([current] + cells)
    with (here / "vcf43-breakends.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write(f"# The breakend tables of VCF 4.3 §5.4, from samtools/hts-specs {COMMIT} VCFv4.3.tex ({URL}), retrieved {today} by make_vcf_sv_fixtures.py\n")
        f.write("# columns: subsection, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"vcf43-structural.vcf: {len(vcf)} lines; vcf43-breakends.tsv: {len(rows)} records")
    return 0


if __name__ == "__main__":
    sys.exit(main())
