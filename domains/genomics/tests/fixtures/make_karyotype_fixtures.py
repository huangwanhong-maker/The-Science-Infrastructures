"""Makes the fixtures of the karyotype module's checks (domains/genomics/karyotype/tests).

Run from the repository root:

    python domains/genomics/tests/fixtures/make_karyotype_fixtures.py

It writes, beside itself:

1. karyotype-cytogps.tsv: the loss-gain-fusion (LGF) model of six karyotypes and the detailed forms
   of their derivatives, as the CytoGPS web service computed them on 2020-04-10 (Abrams et al. 2019,
   Bioinformatics 35, 5365-5366), shipped in the R package RCytoGPS 1.2.14 (CRAN, Apache-2.0) as
   inst/Examples/JSONfiles/CytoGPS_Result1.json and CytoGPS_Result2.json; the bands of each clone's
   loss, gain and fusion are given by name.
2. karyotype-iscn2009.tsv: the examples of ISCN 2009 with their page numbers and, where the book
   gives one, their detailed form, as quoted by the Karyotype Ontology (Warrender, Moorman and Lord
   2019, Bioinformatics 35, 5264-5270; github.com/jaydchan/tawny-karyotype at commit 3a875e5,
   src/ncl/karyotype/iscnexamples.clj, LGPL-3.0). The detailed forms the quotation garbled (listed in
   GARBLED below, each with what is wrong) are left out; the examples are written in ISCN 2009's
   order of breakpoints.
3. karyotype-ordulu2014.tsv: the sequencing descriptions of Ordulu et al. (2014), Am. J. Hum. Genet.
   94, 695-709, read from the article at PubMed Central (PMC4067557): each short system with its
   detailed system, on GRCh37/hg19.
4. karyotype-cll-corpus.txt: the karyotypes of chronic lymphocytic leukaemia as RCytoGPS ships them
   for its pre-cleaning (inst/PreClean/startKaryotypes.csv), one a line, for the parser's coverage of
   karyotypes as laboratories write them.

Only the standard library is used.
"""

from __future__ import annotations

import csv
import datetime as dt
import html.parser
import io
import json
import re
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CRAN = "https://cran.r-project.org/src/contrib/Archive/RCytoGPS/RCytoGPS_1.2.14.tar.gz"
CRAN_CURRENT = "https://cran.r-project.org/src/contrib/RCytoGPS_1.2.14.tar.gz"
ONTOLOGY = "https://raw.githubusercontent.com/jaydchan/tawny-karyotype/3a875e5ef7de9372aee70d9e6cc8c22fc5187544/src/ncl/karyotype/iscnexamples.clj"
ORDULU = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4067557/"
AGENT = {"User-Agent": "Mozilla/5.0 (Noetarium; fixtures of genomics/karyotype)"}

# Detailed forms the Karyotype Ontology's quotation garbled: the example is kept, its detailed form
# is not.
GARBLED = {
    "46,XX,t(3;7;7*)(q21;q22;p13)": "the detailed form names t(2;7;7*) and chromosome 2 for chromosome 3",
    "46,XX,t(2;7;5)(p21;q22;q23)": "7pter7q22 lacks its arrow",
    "46,X,t(X;22;1)(q24;q11.2;p33)": "1p33_>1pter for 1p33->1pter",
    "46,XX,t(3;9;9*;22)(p13;q22;q34;q11.2)": "9*ter for 9*qter",
    "46,XX,t(5;14;9)(q13q23;q24q21;p12p23)": "the detailed form names t(5;14;19) and is broken by a stray +",
    "46,XX,ins(5;2)(q31;p23p13)": "the two derivatives are joined by an arrow, not a semicolon",
    "46,XX,ins(2;2*)(p13;q21q31)": "the detailed form is that of ins(2)(p13q21q31)",
}
# Examples whose quotation is not a karyotype at all.
BROKEN = {
    "46,XY,der(1)(del(1)(p21p33)hsr(1)(p21)": "an unbalanced parenthesis",
}


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=120) as r:
        return r.read()


def rcytogps() -> tarfile.TarFile:
    for url in (CRAN_CURRENT, CRAN):
        try:
            return tarfile.open(fileobj=io.BytesIO(fetch(url)), mode="r:gz")
        except Exception:  # the current release moves to the archive when a newer one appears
            continue
    raise SystemExit("RCytoGPS 1.2.14 could not be fetched from CRAN")


def cytogps(tar: tarfile.TarFile, today: str) -> None:
    rows = []
    computed = ""
    for name in ("CytoGPS_Result1.json", "CytoGPS_Result2.json"):
        data = json.loads(tar.extractfile(f"RCytoGPS/inst/Examples/JSONfiles/{name}").read())
        bands = data["iscn2016_bands"]
        computed = data.get("date", computed)
        for item in data["output"]:
            clones = item.get("parsing_result", [])
            if not clones:
                rows.append((item["karyotype"], item["status"], "", "", "", "", ""))
            for k, clone in enumerate(clones):
                lgf = clone["loss_gain_fusion_computing"]
                pick = lambda key: ",".join(bands[i] for i, x in enumerate(lgf[key]) if x)
                rows.append((item["karyotype"], item["status"], str(k), pick("loss"), pick("gain"), pick("fusion"),
                             "|".join(clone.get("derivative_chromosome_detailed_systems", []))))
    with open(HERE / "karyotype-cytogps.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# CytoGPS's LGF model (computed {computed}) from RCytoGPS 1.2.14 (CRAN, Apache-2.0), retrieved {today}\n")
        f.write("# karyotype\tstatus\tclone\tloss\tgain\tfusion\tdetailed\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def corpus(tar: tarfile.TarFile) -> None:
    text = tar.extractfile("RCytoGPS/inst/PreClean/startKaryotypes.csv").read().decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text)))
    lines = [" ".join(r["Karyotype"].split()) for r in rows if r.get("Karyotype")]
    header = ("# Karyotypes of chronic lymphocytic leukaemia as RCytoGPS 1.2.14 ships them for its pre-cleaning\n"
              "# (inst/PreClean/startKaryotypes.csv, Apache-2.0), one a line; not all of them are ISCN.\n")
    with open(HERE / "karyotype-cll-corpus.txt", "w", encoding="utf-8", newline="\n") as f:
        f.write(header + "\n".join(lines) + "\n")


def iscn2009(today: str) -> None:
    text = fetch(ONTOLOGY).decode("utf-8")
    # Commented-out examples count as much as the others: strip the comment marks first.
    text = "\n".join(re.sub(r"^\s*;;", "", line) for line in text.splitlines())
    rows, seen = [], set()
    for m in re.finditer(r':label\s+"The\s+([^"]*?)\s+karyotype"(\s*:comment\s+"([^"]*)")?', text):
        short = " ".join(m.group(1).split())
        comment = " ".join((m.group(3) or "").split())
        if short in seen or short in BROKEN:
            continue
        seen.add(short)
        page = re.search(r"ISCN(\d{4}) pg (\d+)", comment)
        if not page:
            continue  # the ontology's own karyotypes, not quoted from ISCN
        tail = comment[comment.rfind("'") + 1:].strip() if comment.count("'") >= 2 else ""
        detailed = tail.split(" or ")[0].strip() if "->" in tail else ""
        if short in GARBLED:
            detailed = ""
        rows.append((page.group(1) if page else "", page.group(2) if page else "", short, detailed))
    with open(HERE / "karyotype-iscn2009.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# ISCN 2009 examples as quoted by the Karyotype Ontology (tawny-karyotype, commit 3a875e5), retrieved {today}\n")
        f.write("# edition\tpage\tshort\tdetailed\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


class Text(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []
        self.skip = 0

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.skip += 1
        if tag in ("p", "div", "br", "li", "tr", "h1", "h2", "h3", "h4"):
            self.parts.append("\n")

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.skip = max(0, self.skip - 1)

    def handle_data(self, data):
        if not self.skip:
            self.parts.append(data)


def ordulu(today: str) -> None:
    p = Text()
    p.feed(fetch(ORDULU).decode("utf-8", "replace"))
    lines = [" ".join(x.split()) for x in "".join(p.parts).split("\n")]
    lines = [x for x in lines if x]
    rows = []
    for i, line in enumerate(lines):
        if line.startswith("Short system: "):
            short = line[len("Short system: "):]
            detailed = ""
            for j in range(i + 1, min(i + 4, len(lines))):
                if lines[j].startswith("Detailed system: "):
                    detailed = lines[j][len("Detailed system: "):]
                    break
            rows.append((short, detailed))
    with open(HERE / "karyotype-ordulu2014.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# Ordulu et al. (2014), Am J Hum Genet 94:695-709, PMC4067557: short and detailed systems, retrieved {today}\n")
        f.write("# short\tdetailed\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def main() -> None:
    today = dt.date.today().isoformat()
    tar = rcytogps()
    cytogps(tar, today)
    corpus(tar)
    iscn2009(today)
    ordulu(today)
    print("wrote karyotype-cytogps.tsv, karyotype-iscn2009.tsv, karyotype-ordulu2014.tsv, karyotype-cll-corpus.txt")


if __name__ == "__main__":
    main()
