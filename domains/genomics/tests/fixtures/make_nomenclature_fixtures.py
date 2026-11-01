"""Reference answers for genomics/nomenclature's checks, from three services that implement the
HGVS Nomenclature, SPDI and the Sequence Ontology's consequences:

  Mutalyzer 3 (https://mutalyzer.nl/api; Lefter et al. 2021, Bioinformatics 37, 2811-2817):
      description_extract (the description of an observed sequence against a reference) and
      normalize (a description normalised, with its protein and RNA predictions and genomic equivalent);
  NCBI Variation Services (https://api.ncbi.nlm.nih.gov/variation/v0; Holmes et al. 2020):
      the contextual SPDI (VOCA) of a genomic description, its VCF fields and its HGVS;
  Ensembl VEP's REST service (https://rest.ensembl.org; McLaren et al. 2016), release 116, RefSeq
      transcripts: consequence terms, HGVS c. and p.

It writes three tables, each with the query, the answer, the services' versions and the date:

  nomenclature-extractor.tsv  reference, observed, Mutalyzer's description (synthetic sequences:
                              hand-made cases of every kind, and 80 random ones of a seeded generator)
  nomenclature-variants.tsv   real variants of HBB, CFTR, BRCA1 and the mitochondrial genome:
                              Mutalyzer's c., p., r. and g.; NCBI's SPDI, VCF and HGVS; VEP's terms, c. and p.
  nomenclature-positions.tsv  genomic positions of the three genes' MANE transcripts (exon edges, the
                              middle of every intron, UTRs, flanks) and their c. numbering by VEP

Reads the data library's cache (data/cache/ncbi-human-genes, data/cache/ncbi-genomes) for the
reference bases at the positions it asks about. Usage: python make_nomenclature_fixtures.py
"""

from __future__ import annotations

import datetime as dt
import json
import random
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CACHE = ROOT / "data" / "cache"
MUTALYZER = "https://mutalyzer.nl/api"
NCBI = "https://api.ncbi.nlm.nih.gov/variation/v0"
VEP = "https://rest.ensembl.org"


def get(url: str, data: bytes | None = None, attempts: int = 3):
    for k in range(attempts):
        try:
            req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "the-science-infrastructures-fixtures"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            try:
                return json.loads(e.read().decode("utf-8"))
            except Exception:
                if k + 1 == attempts:
                    return {"error": f"HTTP {e.code}"}
        except Exception as e:  # network
            if k + 1 == attempts:
                return {"error": str(e)}
        time.sleep(2 * (k + 1))
    return {"error": "no answer"}


def q(s: str) -> str:
    return urllib.parse.quote(s, safe="")


# ---- the genomic records (for reference bases) ---------------------------------------------------------
def genbank_sequence(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    origin = text[text.index("\nORIGIN") :]
    return re.sub(r"[^acgtnACGTN]", "", origin.split("\n", 1)[1]).upper()


def mrna_exons(path: Path, transcript: str):
    """The exons (1-based inclusive, slice coordinates) of an mRNA feature, and its strand."""
    text = path.read_text(encoding="utf-8")
    for m in re.finditer(r"^     mRNA            (.*?)(?=^     \S)", text, re.S | re.M):
        block = m.group(0)
        if f'/transcript_id="{transcript}"' not in block:
            continue
        loc = "".join(l.strip() for l in block.split("/")[0].splitlines())
        loc = loc.replace("mRNA", "").strip()
        reverse = loc.startswith("complement(")
        spans = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc)]
        if reverse:
            spans = list(reversed(spans))
        return spans, reverse
    raise SystemExit(f"{transcript} not found in {path.name}")


GENES = [
    # (gene, transcript, chromosome accession, VEP's chromosome, file, offset)
    ("HBB", "NM_000518.5", "NC_000011.10", "11", "ncbi-human-genes/NC_000011.10_5222001-5232000.gb", 5222000),
    ("CFTR", "NM_000492.4", "NC_000007.14", "7", "ncbi-human-genes/NC_000007.14_117479001-117670000.gb", 117479000),
    ("BRCA1", "NM_007294.4", "NC_000017.11", "17", "ncbi-human-genes/NC_000017.11_43043001-43172000.gb", 43043000),
]

# Real variants: (gene, description given to Mutalyzer, what it is).
VARIANTS = [
    ("HBB", "NC_000011.10(NM_000518.5):c.20A>T", "HbS, sickle cell"),
    ("HBB", "NC_000011.10(NM_000518.5):c.19G>A", "HbC"),
    ("HBB", "NC_000011.10(NM_000518.5):c.79G>A", "HbE"),
    ("HBB", "NC_000011.10(NM_000518.5):c.9T>C", "synonymous His3"),
    ("HBB", "NC_000011.10(NM_000518.5):c.1A>G", "start lost"),
    ("HBB", "NC_000011.10(NM_000518.5):c.118C>T", "codon 39 nonsense, beta0-thalassaemia"),
    ("HBB", "NC_000011.10(NM_000518.5):c.27dup", "codons 8/9 +G frameshift"),
    ("HBB", "NC_000011.10(NM_000518.5):c.25_26del", "codon 8 -AA frameshift"),
    ("HBB", "NC_000011.10(NM_000518.5):c.126_129del", "codons 41/42 -TTCT frameshift"),
    ("HBB", "NC_000011.10(NM_000518.5):c.135del", "codon 44 -C frameshift"),
    ("HBB", "NC_000011.10(NM_000518.5):c.92+1G>A", "IVS-I-1 splice donor"),
    ("HBB", "NC_000011.10(NM_000518.5):c.92+5G>C", "IVS-I-5"),
    ("HBB", "NC_000011.10(NM_000518.5):c.93-21G>A", "IVS-I-110"),
    ("HBB", "NC_000011.10(NM_000518.5):c.316-197C>T", "IVS-II-654, deep intronic"),
    ("HBB", "NC_000011.10(NM_000518.5):c.-78A>G", "the -28 TATA box, upstream of the transcript"),
    ("HBB", "NC_000011.10(NM_000518.5):c.-50A>C", "the cap site, the first base of the 5' UTR"),
    ("HBB", "NC_000011.10(NM_000518.5):c.*110T>C", "the polyadenylation signal, 3' UTR"),
    ("HBB", "NC_000011.10(NM_000518.5):c.444A>T", "stop lost"),
    ("HBB", "NC_000011.10(NM_000518.5):c.442_444del", "the stop codon deleted"),
    ("HBB", "NC_000011.10(NM_000518.5):c.10_12del", "an in-frame deletion"),
    ("HBB", "NC_000011.10(NM_000518.5):c.10_11insGGC", "an in-frame insertion"),
    ("HBB", "NC_000011.10(NM_000518.5):c.20_22delinsTTT", "an in-frame delins"),
    ("HBB", "NC_000011.10(NM_000518.5):c.82_84inv", "an inversion"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.1521_1523del", "F508del, the 3' rule"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.1520_1522del", "F508del written 5' of its 3' place"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.1652G>A", "G551D"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.1624G>T", "G542X"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.3846G>A", "W1282X"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.350G>A", "R117H"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.3909C>G", "N1303K"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.3484C>T", "R1162X"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.489+1G>T", "621+1G>T splice donor"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.1585-1G>A", "1717-1G>A splice acceptor"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.3718-2477C>T", "3849+10kbC>T deep intronic"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.2051_2052delinsG", "2183AA>G frameshift"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.1210-12T[5]", "the 5T allele of the polythymidine tract"),
    ("CFTR", "NC_000007.14(NM_000492.4):c.3140-26A>G", "3272-26A>G"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.68_69del", "185delAG"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.5266dup", "5382insC"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.181T>G", "C61G"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.4327C>T", "R1443X"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.5096G>A", "R1699Q"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.4357+1G>A", "splice donor"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.213-11T>G", "splice region"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.5503C>T", "R1835X"),
    ("BRCA1", "NC_000017.11(NM_007294.4):c.1A>G", "start lost, minus strand"),
    ("MT", "NC_012920.1:m.3243A>G", "MELAS, MT-TL1"),
    ("MT", "NC_012920.1:m.8993T>G", "NARP, MT-ATP6 L156R"),
    ("MT", "NC_012920.1:m.11778G>A", "LHON, MT-ND4 R340H"),
    ("MT", "NC_012920.1:m.14484T>C", "LHON, MT-ND6 M64V, minus strand"),
    ("MT", "NC_012920.1:m.3460G>A", "LHON, MT-ND1 A52T"),
    ("MT", "NC_012920.1:m.1555A>G", "aminoglycoside deafness, MT-RNR1"),
    ("MT", "NC_012920.1:m.8344A>G", "MERRF, MT-TK"),
]

# Hand-made changes of every kind on short sequences: (reference, observed).
EXTRACTOR_CASES = [
    ("ATGCAGCAGCAGTTT", "ATGCAGCAGCAGCAGCAGTTT"),  # two copies of a tandem unit: a repeat
    ("ATGCAGCAGCAGTTT", "ATGCAGCAGCAGCAGTTT"),  # one copy: a duplication
    ("ATGCAGCAGCAGTTT", "ATGCAGTTT"),  # two units deleted
    ("ACGTTTTTACG", "ACGTTTTACG"),  # one of a T-run deleted, the 3' rule
    ("ACGTTTTTACG", "ACGTTTTTTACG"),  # a T-run lengthened by one: dup
    ("ACGTTTTTACG", "ACGTTTTTTTACG"),  # by two: a repeat
    ("ACGTACGTAAGG", "ACGTACGTCAGG"),  # a substitution
    ("AAAATTTCCCCCGGGG", "AAAATTGGGGGCCCGGGG"),  # a delins
    ("GCATCAGCCTAA", "GCACTGACCTAA"),  # an inversion (TCAG -> CTGA)
    ("ATGAGACT", "ATGACT"),  # SPDI's Table 2 deletion
    ("ATGACT", "ATGAGACT"),  # SPDI's Table 2 insertion
    ("TTGCATGCCAAGT", "TTGCATGCCATGCCAAGT"),  # a duplication of a non-repeated sequence
    ("ACGGTACCATTAG", "ACGGTAGGCCATTAG"),  # an insertion that is not a copy
    ("CCCATCTTTGGTGTT", "CCCATCTTGGTGTT"),  # CFTR's ATCTTTGGT, one T deleted
    ("CCCATCTTTGGTGTT", "CCCATCGGTGTT"),  # CTT deleted (F508del's context)
]


def random_case(rng: random.Random):
    """A short reference with a tandem repeat and one random change."""
    unit = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 4)))
    left = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 8)))
    right = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 8)))
    ref = left + unit * rng.randint(2, 5) + right
    kind = rng.choice(["del", "ins", "sub", "delins", "unit_ins", "unit_del", "inv"])
    p = rng.randint(1, len(ref) - 3)
    if kind == "del":
        L = rng.randint(1, 3)
        obs = ref[:p] + ref[p + L :]
    elif kind == "ins":
        obs = ref[:p] + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 4))) + ref[p:]
    elif kind == "sub":
        b = rng.choice([x for x in "ACGT" if x != ref[p]])
        obs = ref[:p] + b + ref[p + 1 :]
    elif kind == "delins":
        L = rng.randint(2, 4)
        obs = ref[:p] + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 4))) + ref[p + L :]
    elif kind == "unit_ins":
        at = len(left) + len(unit) * rng.randint(0, 2)
        obs = ref[:at] + unit * rng.randint(1, 3) + ref[at:]
    elif kind == "unit_del":
        at = len(left)
        obs = ref[:at] + ref[at + len(unit) :]
    else:
        L = rng.randint(3, 6)
        seg = ref[p : p + L]
        rc = seg[::-1].translate(str.maketrans("ACGT", "TGCA"))
        obs = ref[:p] + rc + ref[p + L :]
    return ref, obs


def versions():
    v = {}
    m = get(f"{MUTALYZER}/version")
    v["mutalyzer"] = f"Mutalyzer {m.get('mutalyzer')} (API {m.get('api')})"
    s = get(f"{VEP}/info/software?content-type=application/json")
    v["vep"] = f"Ensembl VEP REST, release {s.get('release')}"
    v["ncbi"] = "NCBI Variation Services v0 (var_service.yaml 0.0.0-dev)"
    return v


def main() -> int:
    today = dt.date.today().isoformat()
    ver = versions()
    header = f"# retrieved {today} by make_nomenclature_fixtures.py; {ver['mutalyzer']}; {ver['ncbi']}; {ver['vep']}\n"

    # ---- A. the extractor -------------------------------------------------------------------------------
    rng = random.Random(20261101)
    cases = list(EXTRACTOR_CASES) + [random_case(rng) for _ in range(80)]
    with (HERE / "nomenclature-extractor.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("# Mutalyzer 3's description_extract on short synthetic sequences (hand-made cases, then 80 random ones, seed 20261101)\n")
        f.write(header)
        f.write("# columns: reference, observed, Mutalyzer's description\n")
        for ref, obs in cases:
            if ref == obs:
                continue
            r = get(f"{MUTALYZER}/description_extract/?reference={ref}&observed={obs}")
            desc = r if isinstance(r, str) else json.dumps(r)
            f.write(f"{ref}\t{obs}\t{desc}\n")
            time.sleep(0.2)
    print("nomenclature-extractor.tsv written")

    # ---- B. real variants ---------------------------------------------------------------------------------
    rows = []
    for gene, desc, what in VARIANTS:
        m = get(f"{MUTALYZER}/normalize/{q(desc)}")
        if "custom" in m:
            m = m["custom"]
        norm = m.get("normalized_description", "")
        prot = (m.get("protein") or {}).get("description", "")
        rna = (m.get("rna") or {}).get("description", "")
        eq = (m.get("equivalent_descriptions") or {}).get("g") or []
        g = eq[0]["description"] if eq else (norm if ":m." in norm or ":g." in norm else "")
        errors = ";".join(e.get("code", "") for e in (m.get("errors") or []))
        time.sleep(0.3)
        spdi = vcf = ncbi_hgvs = ""
        if g:
            c = get(f"{NCBI}/hgvs/{q(g)}/contextuals")
            sp = ((c.get("data") or {}).get("spdis") or [])
            if sp:
                s = sp[0]
                spdi = f"{s['seq_id']}:{s['position']}:{s['deleted_sequence']}:{s['inserted_sequence']}"
                vf = get(f"{NCBI}/spdi/{q(spdi)}/vcf_fields").get("data") or {}
                if vf:
                    vcf = f"{vf.get('chrom')}:{vf.get('pos')}:{vf.get('ref')}:{vf.get('alt')}"
                h = get(f"{NCBI}/spdi/{q(spdi)}/hgvs").get("data") or {}
                ncbi_hgvs = h.get("hgvs", "")
            time.sleep(0.4)
        terms = hgvsc = hgvsp = ""
        if g:
            transcript = re.search(r"\((N[MR]_[0-9.]+)\)", desc)
            if gene == "MT":
                pos = re.search(r"m\.(\d+)([ACGT])>([ACGT])", g)
                v = get(f"{VEP}/vep/human/region/MT:{pos.group(1)}-{pos.group(1)}:1/{pos.group(3)}?hgvs=1&content-type=application/json") if pos else []
            else:
                v = get(f"{VEP}/vep/human/hgvs/{q(g)}?refseq=1&hgvs=1&content-type=application/json")
            if isinstance(v, list) and v:
                for t in v[0].get("transcript_consequences", []):
                    tid = t.get("transcript_id", "")
                    if (transcript and tid == transcript.group(1)) or (gene == "MT" and t.get("hgvsc") and "upstream" not in t.get("consequence_terms", [])[0] and "downstream" not in t.get("consequence_terms", [])[0]):
                        terms = ",".join(t.get("consequence_terms", []))
                        hgvsc = t.get("hgvsc", "") or ""
                        hgvsp = t.get("hgvsp", "") or ""
                        if gene == "MT":
                            terms = t.get("gene_symbol", "") + "|" + terms
                        break
            time.sleep(0.4)
        rows.append((gene, desc, what, norm, prot, rna, g, spdi, vcf, ncbi_hgvs, terms, hgvsc, hgvsp, errors))
        print(desc, "->", norm, prot, g, spdi, terms, hgvsp)
    with (HERE / "nomenclature-variants.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("# Real variants of HBB (NM_000518.5), CFTR (NM_000492.4), BRCA1 (NM_007294.4) and the mitochondrial genome (NC_012920.1)\n")
        f.write(header)
        f.write("# columns: gene, description asked, what it is, Mutalyzer normalized, Mutalyzer protein, Mutalyzer RNA, Mutalyzer genomic equivalent, "
                "NCBI contextual SPDI, NCBI VCF (chrom:pos:ref:alt), NCBI HGVS of the SPDI, VEP consequence terms (MT: gene|terms), VEP hgvsc, VEP hgvsp, Mutalyzer errors\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print("nomenclature-variants.tsv written")

    # ---- C. positions -------------------------------------------------------------------------------------
    out = []
    for gene, transcript, acc, chrom, file, offset in GENES:
        path = CACHE / file
        seq = genbank_sequence(path)
        exons, reverse = mrna_exons(path, transcript)
        positions = set()
        lo = min(a for a, b in exons)
        hi = max(b for a, b in exons)
        for a, b in exons:
            positions.update({a, b})
        # introns, in genomic order
        spans = sorted(exons)
        for (a1, b1), (a2, b2) in zip(spans, spans[1:]):
            i1, i2 = b1 + 1, a2 - 1
            mid = (i1 + i2) // 2
            positions.update({i1, i1 + 1, mid - 1, mid, mid + 1, i2 - 1, i2})
        positions.update({lo - 100, lo - 1, hi + 1, hi + 100})
        lines = []
        for p in sorted(positions):
            ref = seq[p - 1]
            alt = "A" if ref != "A" else "C"
            lines.append((p, ref, alt))
        for k in range(0, len(lines), 150):
            batch = lines[k : k + 150]
            body = json.dumps({"variants": [f"{chrom} {p + offset} . {ref} {alt} . . ." for p, ref, alt in batch]}).encode()
            res = get(f"{VEP}/vep/human/region?refseq=1&hgvs=1&content-type=application/json", data=body)
            if not isinstance(res, list):
                print("VEP:", res)
                continue
            for x in res:
                gpos = int(x["start"])
                for t in x.get("transcript_consequences", []):
                    if t.get("transcript_id") == transcript and t.get("hgvsc"):
                        cpos = re.match(r".*:c\.(.+?)[ACGT]>[ACGT]$", t["hgvsc"])
                        out.append((acc, transcript, str(gpos), cpos.group(1) if cpos else t["hgvsc"], ",".join(t.get("consequence_terms", []))))
            time.sleep(1)
    with (HERE / "nomenclature-positions.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("# c. numbering of genomic positions of three MANE transcripts (exon edges, every intron's ends and middle, UTRs, flanks), by Ensembl VEP's HGVS c.\n")
        f.write(header)
        f.write("# columns: chromosome, transcript, genomic position (1-based), c. position, VEP consequence terms of a substitution there\n")
        for r in out:
            f.write("\t".join(r) + "\n")
    print(f"nomenclature-positions.tsv: {len(out)} positions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
