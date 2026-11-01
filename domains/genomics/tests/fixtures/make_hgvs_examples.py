"""The examples of the HGVS Nomenclature, for the checks of genomics/nomenclature's grammar.

Reads the sources of https://hgvs-nomenclature.org at a pinned commit of
github.com/HGVSnomenclature/hgvs-nomenclature (version 21.1 of the Nomenclature), and writes
hgvs-examples.tsv: every description the recommendations give as an example, with the page and line
it is on and how the standard marks it:

  valid       written as an example (in syntax.yaml's `examples`, or in code on a recommendations page)
  position    a location alone (c.88+1, g.(123_456)_(789_1000)): the checks parse it from a location rule
  invalid     written in <code class="invalid">: the standard's own counter-examples
  in-context  a counter-example whose fault only the reference sequence or the transcript shows, with why
  excluded    an example the standard writes that is not a description in its own syntax, with the
              reason (a placeholder accession, ISCN's commas, a fragment of prose)

A code span is taken as a description when it begins with a coordinate type (g. m. o. c. n. r. p.),
optionally after a reference (an accession, its version and a selector, then ':'), or with '[' and
such a reference (an allele of descriptions). Placeholders of the prose (p.[variant1;variant2]) are
left out.

Usage: python make_hgvs_examples.py   (network: raw.githubusercontent.com)
"""

from __future__ import annotations

import datetime as dt
import html
import re
import sys
import urllib.request
from pathlib import Path

COMMIT = "78c03ac7e412331fb0dccc0d5ae5e408c5b7c26f"  # 2026-10-06, HGVS Nomenclature 21.1
BASE = f"https://raw.githubusercontent.com/HGVSnomenclature/hgvs-nomenclature/{COMMIT}/"
PAGES = [
    "docs/syntax.yaml",
    "docs/recommendations/general.md",
    "docs/recommendations/uncertain.md",
    "docs/background/numbering.md",
    *[f"docs/recommendations/DNA/{p}.md" for p in ("substitution", "deletion", "duplication", "insertion", "inversion", "delins", "repeated", "alleles", "complex", "other")],
    *[f"docs/recommendations/RNA/{p}.md" for p in ("substitution", "deletion", "duplication", "insertion", "inversion", "delins", "repeated", "alleles", "splicing", "adjoined_transcript")],
    *[f"docs/recommendations/protein/{p}.md" for p in ("substitution", "deletion", "duplication", "insertion", "delins", "frameshift", "extension", "repeated", "alleles")],
]

# Examples the standard writes that are not descriptions in its own syntax, with the reason.
EXCLUDED = {
    "NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]": "an insertion with no position (DNA/complex.md, the pericentric inversion of chromosome 2)",
    "NG_012232.1(NM_004006.2):c.123+51_123+52ins[[XXXXXX.y:g.393_1295];123+45_123+51]": "a placeholder accession (XXXXXX.y) in nested brackets",
    "c.IVS2+2T>G": "the IVS notation, retracted (numbering.md)",
    "c.IVS2-1G>T": "the IVS notation, retracted (numbering.md)",
    "p.2366Gln/Lys": "the retracted notation of polymorphisms (protein/substitution.md)",
    "p.2366Q/K": "the retracted notation of polymorphisms (protein/substitution.md)",
    "OPRM1:c.118AA": "a reader's suggestion for genotypes in prose, not adopted (DNA/alleles.md)",
    "OPRM1:c.118GA": "a reader's suggestion for genotypes in prose, not adopted (DNA/alleles.md)",
    "OPRM1:c.118GG": "a reader's suggestion for genotypes in prose, not adopted (DNA/alleles.md)",
    "OPRM1:c.118A": "a reader's suggestion for alleles in prose, not adopted (DNA/alleles.md)",
    "OPRM1:c.118G": "a reader's suggestion for alleles in prose, not adopted (DNA/alleles.md)",
    "GJB2:c.76AC": "a reader's suggestion for genotypes in prose, not adopted (DNA/alleles.md)",
    "chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter": "ISCN's own notation, with commas in numbers, which HGVS does not allow (DNA/complex.md)",
    "chr11:g.pter_15,825,272::chr2:g.8,247,757_cen_qter": "ISCN's own notation, with commas in numbers (DNA/complex.md)",
    "chr2:g.pter_cen_174500098::chr11:g.108111987_qter": "ISCN's own notation, a range through cen (DNA/complex.md)",
    "chr11:g.pter_cen_108111981::chr2:g.174500099_qter": "ISCN's own notation, a range through cen (DNA/complex.md)",
    "chr14:g.pter_cen_29,745,313::chr3:g.pter_36,969,141inv": "ISCN's own notation, with commas and cen (DNA/complex.md)",
    "chr14:g.29,745,314_qterinv::CATTTGTTCAAATTTAGTTCAAATGA::chr3:g.36,969,142_cen_qter": "ISCN's own notation, with commas and cen (DNA/complex.md)",
    "r.(100_150)delN[15]": "an RNA description with the DNA letter N (RNA is written in lower case, r.549_550insn; RNA/deletion.md)",
    "NM_002354.2:r.?_555::guaugauuuuuuaataa::NM_000251.2:r.212_?": "a 't' in an RNA linker, which is written with 'u' (RNA/adjoined_transcript.md)",
    "r.76[a>c];[=]": "an allele written inside a position, which DNA/alleles.md forbids (c.2376[G>C];[G>C]) and RNA/alleles.md recommends",
    "p.Gln[21]": "a residue without its position, in prose (RNA/repeated.md)",
    "p.Q[21]": "a residue without its position, in prose (RNA/repeated.md)",
    "p.Met1ext": "a fragment of p.Met1extMet-5, split by the page's markup (protein/extension.md)",
    "c.-": "a fragment of the numbering's prose",
    "c.*": "a fragment of the numbering's prose",
    "p.[": "a fragment of the alleles' prose",
}

# Counter-examples whose fault no text alone can show: it takes the reference sequence or the transcript.
IN_CONTEXT = {
    "NM_004006.1:c.-128354C>T": "the position lies outside the transcript NM_004006.1 (DNA/substitution.md)",
    "NM_000109.3:c.-401C>T": "the position lies outside the transcript NM_000109.3 (DNA/substitution.md)",
    "NM_004006.2:c.5690del": "the 3' rule: the deleted A is the last of an A-stretch ending at c.5697 (DNA/deletion.md)",
    "NC_000023.11(NM_004006.2):c.(?_-244)_(31+1_32-1)del": "c.-244 is the transcript's first nucleotide; a deletion beyond it is described in g. (DNA/deletion.md)",
    "NC_000023.11(NM_004006.2):c.(10086+1_10087-1)_(*2691_?)del": "c.*2691 is the transcript's last nucleotide (DNA/deletion.md)",
    "NC_000023.11(NM_004006.2):c.(?_-244)_(31+1_32-1)dup": "c.-244 is the transcript's first nucleotide (DNA/duplication.md)",
    "NC_000023.11(NM_004006.2):c.(10086+1_10087-1)_(*2691_?)dup": "c.*2691 is the transcript's last nucleotide (DNA/duplication.md)",
    "NC_000023.11(NM_004006.2):c.(-205839_-62966)_(*21568_*61692)dup": "a whole gene's duplication described beyond its transcript (DNA/duplication.md)",
    "NM_000492.3:c.1522_1524delTTT": "the 3' rule on CFTR's ATCTTTGGT: the change is c.1521_1523del (DNA/deletion.md)",
    "r.1033_1035del": "the 3' rule across an exon junction (RNA/deletion.md)",
    "r.6_7insu": "the inserted u repeats the one before it: r.7dup (RNA/duplication.md)",
    "NG_012232.1(NM_004006.2):r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]": "what the splicing of the transcript produces decides it (RNA/insertion.md)",
    "r.2949_2950ins[2950-30_2950-12;uuag]": "what the splicing of the transcript produces decides it (RNA/insertion.md)",
    "p.(Val89_Gln119del)": "the 3' rule on the protein: p.(Asp90_Val120del) (protein/deletion.md)",
}

# Written in plain code although the text around them says they are wrong.
INVALID_UNMARKED = {
    "NM_004006.2:c.[762_768del;767_774dup]": "general.md: descriptions removing part of a reference sequence and replacing it with part of the same sequence are not allowed",
    "c.[235A>T;237G>T]": "DNA/substitution.md: two changes in one codon, one nucleotide apart, are c.235_237delinsTAT",
    "p.[Lys79*;Lys79Asn]": "DNA/substitution.md: the wrong prediction that c.[235A>T;237G>T] leads to",
    "c.23ins24": "DNA/insertion.md: c.23ins24 is not correct",
    "p.Tyr4TerfsTer1": "protein/frameshift.md: not p.Tyr4TerfsTer1",
    "p.[Ser68_Ala74insSerGln];[Ser68_Ala74=]": "protein/alleles.md inserts between Ser68 and Ala74, which are not adjacent, against protein/insertion.md",
}

# Marked invalid by the page's markup although its text gives it as the right form.
MARKUP_LAPSES = {
    "r.124_500delins[AB053210.2:r.1289-365_1289-73inv]": "RNA/inversion.md: 'the insertion is in an inverted orientation, so inv should be used; r.124_500delins[AB053210.2:r.1289-365_1289-73inv]'",
}

REFERENCE = r"[A-Za-z][A-Za-z0-9_\-]*(?:\.[0-9]+)?(?:\([A-Za-z][A-Za-z0-9_\-]*(?:\.[0-9]+)?\))?"
DESCRIPTION = re.compile(rf"^(?:\[?{REFERENCE}:)?[gmocnrp]\..+|^\[{REFERENCE}:[gmocnrp]\..+")
# A location alone: c.88+1, g.(123_456)_(789_1000), p.His7.
POSITION = re.compile(r"^(?:[A-Za-z][A-Za-z0-9_.()\-]*:)?[gmocnrp]\.(?:[-*]?[0-9?(][0-9+*?()_\-]*|pter|qter|cen|(?:[A-Z][a-z]{2}|[A-Z*])[0-9]+)$")
# A statement about a whole sequence: c.=, p.0?, r.(?).
WHOLE = re.compile(r"^(?:[A-Za-z][A-Za-z0-9_.()\-]*:)?[gmocnrp]\.(?:=|\?|0|0\?|\(=\)|\(\?\))$")
# Words of the prose written in code as placeholders of a form (p.[variant1;variant2], g.(fragment-start_fragment-end)).
PLACEHOLDER = re.compile(r"variant|change|fragment|position|normal|first|last|probe|left|right|XXX|#|…|\.\.\.|\s")


def fetch(path: str) -> str:
    with urllib.request.urlopen(BASE + path, timeout=60) as r:
        return r.read().decode("utf-8")


def candidates(path: str, text: str):
    """(description, line, mark) for every code span or YAML example of a page."""
    for number, line in enumerate(text.splitlines(), 1):
        if path.endswith(".yaml"):
            m = re.match(r"^\s*-\s+(\S.*)$", line)
            if m and ":" in m.group(1) and re.search(r"[gmocnrp]\.", m.group(1)):
                yield m.group(1).strip(), number, "valid"
            continue
        for m in re.finditer(r'<code class="invalid">(.*?)</code>', line):
            yield html.unescape(m.group(1)).replace("\\_", "_").replace("\\*", "*"), number, "invalid"
        stripped = re.sub(r'<code class="invalid">.*?</code>', "", line)
        for m in re.finditer(r"`([^`]+)`", stripped):
            yield m.group(1).replace("\\_", "_").replace("\\*", "*"), number, "valid"


def main() -> int:
    here = Path(__file__).resolve().parent
    seen: dict[str, tuple[str, int, str]] = {}
    for path in PAGES:
        text = fetch(path)
        for desc, line, mark in candidates(path, text):
            desc = desc.strip()
            if desc not in EXCLUDED and (not DESCRIPTION.match(desc) or PLACEHOLDER.search(desc)):
                continue
            key = desc + "\t" + mark
            if key not in seen:
                seen[key] = (path.removeprefix("docs/"), line, mark)
    rows, done = [], set()
    for key, (page, line, mark) in seen.items():
        desc = key.split("\t")[0]
        if desc in EXCLUDED:
            status, reason = "excluded", EXCLUDED[desc]
        elif desc in IN_CONTEXT:
            status, reason = "in-context", IN_CONTEXT[desc]
        elif desc in INVALID_UNMARKED:
            status, reason = "invalid", "written as an example, but " + INVALID_UNMARKED[desc]
        elif mark == "valid" and re.match(r"^(?:[A-Za-z][A-Za-z0-9_.()\-]*:)?p\.[0-9]+$", desc):
            status, reason = "excluded", "the numbering's prose (p.1, p.2, ...): a residue's position is written with its residue"
        elif desc in MARKUP_LAPSES:
            status, reason = "valid", "marked invalid by the page's markup: " + MARKUP_LAPSES[desc]
        elif mark == "valid" and POSITION.match(desc) and not WHOLE.match(desc):
            status, reason = "position", ""
        else:
            status, reason = mark, ""
        if (desc, status) in done:
            continue
        done.add((desc, status))
        rows.append((desc, status, page, line, reason))
    today = dt.date.today().isoformat()
    out = here / "hgvs-examples.tsv"
    with out.open("w", encoding="utf-8", newline="\n") as f:
        f.write(f"# The examples of the HGVS Nomenclature 21.1, from github.com/HGVSnomenclature/hgvs-nomenclature at {COMMIT}\n")
        f.write(f"# retrieved {today} by make_hgvs_examples.py; columns: description, status (valid, position, invalid, in-context, excluded), page, line, reason\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    counts = {s: sum(1 for r in rows if r[1] == s) for s in ("valid", "position", "invalid", "in-context", "excluded")}
    print(f"{out.name}: {len(rows)} examples {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
