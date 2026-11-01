"""Reference values for the checks of the genomics domain's codon module (codon/tests/check_*.cpp).

Run with Biopython and codon-bias installed (the values recorded were made with Python 3.12,
Biopython 1.88 and codon-bias 0.5.0, 2026-11-01):

    python make_codon_fixtures.py <data/cache/ncbi-genomes>

It reads the RefSeq genomes the data library fetches (data/registry/ncbi-genomes.json) and writes:

- codon-ecoli-genes.tsv: for every CDS of Escherichia coli K-12 (NC_000913.3) that codon usage
  counts (the exclusions of codon/genes.hpp, applied here to Biopython's parse of the record), its
  key (locus tag, 0-based start and end, strand), its length in codons, its CAI by Biopython's
  Bio.SeqUtils.CodonAdaptationIndex against the ribosomal proteins (products "30S/50S ribosomal
  subunit protein …", the module's RIBOSOMAL_PROTEINS) under table 11, and its effective number of
  codons by codon-bias's EffectiveNumberOfCodons in Wright's form (robust=False, pseudocount=0,
  mean='unweighted') with a flag saying whether codon-bias's extra rule (it drops amino acids whose
  F is at most 1e-6, which Wright's definition keeps) left the value unchanged.
- codon-weights.tsv: the relative adaptiveness w of every codon by Biopython (Sharp and Li's 0.5 for
  codons absent from the reference set) for E. coli's ribosomal proteins under table 11 and for the
  thirteen genes of the human mitochondrion (NC_012920.1) under table 2; and the tRNA adaptiveness
  w of E. coli's codons by codon-bias's TrnaAdaptationIndex (s_values='dosReis') from the genome's
  elongator tRNA genes, with prokaryote=False: codon-bias's prokaryotic row (C, A, 0.89) lets every
  C34 anticodon read NNA as lysidine does, where only tRNA-Ile2's C34 is lysidine, so the check
  compares the other rules here and AUA by hand; and the pooled genome's N_c by codon-bias in Sun, Yang
  and Xia's form (robust=True, pseudocount=1, mean='weighted', families not split).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import AfterPosition, BeforePosition
from Bio.SeqUtils import CodonAdaptationIndex
from codonbias.scores import EffectiveNumberOfCodons, TrnaAdaptationIndex

HERE = Path(__file__).resolve().parent
RIBOSOMAL = re.compile(r"^(30S|50S) ribosomal (subunit )?protein [LS][0-9]+[A-Z]?$", re.IGNORECASE)


def coding_genes(record, default_table):
    """The CDSs codon usage counts, as codon/genes.hpp chooses them."""
    out = []
    for f in record.features:
        if f.type != "CDS":
            continue
        q = f.qualifiers
        if "pseudo" in q or "pseudogene" in q:
            continue
        if any(isinstance(p, (BeforePosition, AfterPosition)) for part in f.location.parts for p in (part.start, part.end)):
            continue
        if "ribosomal_slippage" in q or "exception" in q or "artificial_location" in q:
            continue
        excepts = q.get("transl_except", [])
        if any("aa:TERM" not in e for e in excepts):
            continue
        completed = bool(excepts)
        table = int(q.get("transl_table", [default_table])[0])
        x = str(f.extract(record.seq)).upper()
        frame = int(q.get("codon_start", ["1"])[0]) - 1
        x = x[frame:]
        if completed:
            x = x[: len(x) - len(x) % 3]
        if len(x) % 3:
            continue
        code = CodonTable.unambiguous_dna_by_id[table]
        if len(x) >= 3 and x[-3:] in code.stop_codons:
            x = x[:-3]
        codons = [x[i : i + 3] for i in range(0, len(x), 3)]
        if any(c in code.stop_codons for c in codons):
            continue
        tag = q.get("locus_tag", [""])[0]
        start = min(int(p.start) for p in f.location.parts)
        stop = max(int(p.end) for p in f.location.parts)
        strand = "-" if f.location.strand == -1 else "+"
        out.append({"key": f"{tag}|{start}|{stop}|{strand}", "cds": x, "product": q.get("product", [""])[0], "table": table})
    return out


def trna_counts(record):
    """Elongator tRNA genes by anticodon (initiators and selenocysteine's left out)."""
    counts = {}
    for f in record.features:
        if f.type != "tRNA" or "pseudo" in f.qualifiers:
            continue
        note = f.qualifiers.get("note", [""])[0]
        m = re.search(r"tRNA-(initiator )?([A-Za-z]{3,4})\(([ACGU]{3})\)", note)
        if not m or m.group(1) or m.group(2) in ("Sec", "fMet"):
            continue
        counts[m.group(3)] = counts.get(m.group(3), 0) + 1
    return pd.DataFrame({"anti_codon": list(counts), "GCN": list(counts.values())})


def main(cache: Path) -> None:
    ecoli = SeqIO.read(cache / "NC_000913.3.gb", "genbank")
    genes = coding_genes(ecoli, 11)
    reference = [g["cds"] for g in genes if RIBOSOMAL.search(g["product"])]
    table11 = CodonTable.unambiguous_dna_by_id[11]
    cai = CodonAdaptationIndex(reference, table11)
    wright = EffectiveNumberOfCodons(robust=False, pseudocount=0, mean="unweighted", genetic_code=11)
    rows = []
    for g in genes:
        F, P, N = wright._calc_F(g["cds"])
        deg = wright._aa_deg
        comparable = bool(np.all((N <= 1) | (F > 1e-6) | (deg < 2)))
        rows.append(f"{g['key']}\t{len(g['cds']) // 3}\t{cai.calculate(g['cds']):.15g}\t{wright.get_score(g['cds']):.15g}\t{int(comparable)}")
    with open(HERE / "codon-ecoli-genes.tsv", "w", newline="\n") as out:
        out.write("# Escherichia coli K-12 (NC_000913.3): CAI against its ribosomal proteins (Biopython 1.88) and Wright's N_c (codon-bias 0.5.0).\n")
        out.write(f"# {len(genes)} genes, {len(reference)} ribosomal proteins in the reference set. Made by make_codon_fixtures.py, 2026-11-01.\n")
        out.write("key\tcodons\tcai\tnc_wright\tcomparable\n")
        out.write("\n".join(rows) + "\n")

    lines = ["# set\tcodon\tvalue"]
    for codon in sorted(cai):
        lines.append(f"cai-ecoli-ribosomal\t{codon}\t{cai[codon]:.15g}")
    human = SeqIO.read(cache / "NC_012920.1.gb", "genbank")
    mito = coding_genes(human, 2)
    cai2 = CodonAdaptationIndex([g["cds"] for g in mito], CodonTable.unambiguous_dna_by_id[2])
    for codon in sorted(cai2):
        lines.append(f"cai-human-mitochondrion\t{codon}\t{cai2[codon]:.15g}")
    tai = TrnaAdaptationIndex(tGCN=trna_counts(ecoli), prokaryote=False, s_values="dosReis", genetic_code=11)
    for codon, w in sorted(tai.weights.items()):
        lines.append(f"tai-ecoli\t{codon}\t{w:.15g}")
    sun = EffectiveNumberOfCodons(robust=True, pseudocount=1, mean="weighted", genetic_code=11)
    pooled = "".join(g["cds"] for g in genes)
    lines.append(f"nc-sun-ecoli-pooled\t-\t{sun.get_score(pooled):.15g}")
    lines.append(f"genes-human-mitochondrion\t-\t{len(mito)}")
    with open(HERE / "codon-weights.tsv", "w", newline="\n") as out:
        out.write("# Codon weights by Biopython 1.88 (CAI's w, Sharp and Li 1987) and codon-bias 0.5.0 (tAI's w, dos Reis et al. 2004; N_c of Sun, Yang and Xia 2013).\n")
        out.write("# Made by make_codon_fixtures.py from data/cache/ncbi-genomes (NC_000913.3, NC_012920.1), 2026-11-01.\n")
        out.write("\n".join(lines) + "\n")
    print(len(genes), "E. coli genes,", len(reference), "ribosomal proteins,", len(mito), "mitochondrial genes")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
