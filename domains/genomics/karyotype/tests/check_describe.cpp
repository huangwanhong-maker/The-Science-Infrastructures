// Checks of the descriptions written from results (describe.hpp): copy-number segments as arr, on
// GRCh38 (Conlin et al. 2025's examples, exactly) and GRCh37 (Ordulu et al. 2014's); normal results;
// whole chromosomes grouped; mosaic ranges, proportions and homozygosity; junctions as seq in the
// short and detailed systems, Ordulu et al.'s t(1;8), t(5;10) and del(12) exactly on GRCh37, and
// inversions, duplications, unbalanced derivatives and complex rearrangements recognised; and the
// circle closed: what is written, read back and interpreted, gives the segments and junctions it was
// written from.

#include <genomics/karyotype/describe.hpp>
#include <genomics/karyotype/interpret.hpp>

#include <calx/verification/check.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace genomics::karyotype;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-ideogram";
  if (!std::filesystem::exists(cache / "ideogram_9606_GCF_000001305.16_850_V1.tsv") || !std::filesystem::exists(cache / "ucsc-hg19-cytoBand.txt")) {
    s.skip("the ideograms are not in the data cache (python data/acquisition/fetch.py ncbi-ideogram): descriptions are not checked");
    return s.done();
  }
  const Ideogram g38 = read_ideogram(read_file(cache / "ideogram_9606_GCF_000001305.16_850_V1.tsv"), 850);
  const Ideogram g37 = read_cytoband(read_file(cache / "ucsc-hg19-cytoBand.txt"), 850, nullptr, "GRCh37");

  // ---- arr --------------------------------------------------------------------------------------------
  {
    // Conlin et al. (2025), Table 1 and text, GRCh38.
    s.require("1:147,105,904-147,917,509 at one copy is arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1, as Conlin et al. write it",
              describe_arr({{0, 147105903, 147917509, 1}}, g38) == "arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1",
              describe_arr({{0, 147105903, 147917509, 1}}, g38));
    s.require("7:73,265,388-74,726,130 at one copy is arr[GRCh38] 7q11.23(73,265,388_74,726,130)x1, as Conlin et al. write it",
              describe_arr({{6, 73265387, 74726130, 1}}, g38) == "arr[GRCh38] 7q11.23(73,265,388_74,726,130)x1");
    ArrOptions o37;
    o37.build = "hg19";
    s.require("Ordulu et al.'s deletion, 12:21,340,001-23,500,000 on GRCh37, is arr[hg19] 12p12.1(21,340,001_23,500,000)x1 (with ISCN 2016's underscore)",
              describe_arr({{11, 21340000, 23500000, 1}}, g37, o37) == "arr[hg19] 12p12.1(21,340,001_23,500,000)x1");
    ArrOptions female, male;
    female.male = false;
    male.male = true;
    s.require("a normal result is arr(X,1-22)x2 for a female and arr(X,Y)x1,(1-22)x2 for a male",
              describe_arr({}, g38, female) == "arr (X,1-22)x2" && describe_arr({}, g38, male) == "arr (X,Y)x1,(1-22)x2", describe_arr({}, g38, male));
    std::vector<CopySegment> whole;
    for (int c = 0; c < 22; c++)
      if (c != 20) whole.push_back({c, 0, g38.length[c], 3});
    whole.push_back({chr_x, 0, g38.length[chr_x], 2});
    s.require("whole chromosomes with the same copies grouped, the sex chromosomes first: arr[GRCh38] (X)x2,(1-20,22)x3 in a male",
              describe_arr(whole, g38, male) == "arr[GRCh38] (X)x2,(1-20,22)x3", describe_arr(whole, g38, male));
    std::vector<CopySegment> mixed{{1, 10000000, 20000000, 2, 3, 0.3}, {10, 230679, 10807544, 2}, {chr_x, 0, 2700000, 3}};
    mixed[1].homozygous = true;
    s.require("a mosaic range of copies with its proportion, a homozygous run, and the X first",
              describe_arr(mixed, g38, female) == "arr[GRCh38] Xp22.33(1_2,700,000)x3,2p25.1p24.1(10,000,001_20,000,000)x2~3[0.3],11p15.5p15.4(230,680_10,807,544)x2 hmz",
              describe_arr(mixed, g38, female));
    s.require("segments at their normal copies are left out", describe_arr({{4, 1000, 2000, 2}}, g38, female) == "arr (X,1-22)x2");
  }

  // ---- seq: Ordulu et al. (2014), GRCh37 ---------------------------------------------------------------------
  SeqOptions o;
  o.build = "GRCh37/hg19";
  {
    // t(1;8)(q23.3;q12.2): der(1) joins 1 up to base 164,789,896 to 8 from base 61,760,650; der(8) joins 8
    // up to 61,760,649 to 1 from 164,789,897 (Ordulu et al., Figure 3 and "General Principles").
    const SeqDescription t = describe_seq({{0, 164789896, true, 7, 61760649, true}, {7, 61760649, true, 0, 164789896, true}}, g37, o);
    s.require("Ordulu et al.'s t(1;8): seq[GRCh37/hg19] t(1;8)(q23.3;q12.2)", t.kind == "t" && t.short_form == "seq[GRCh37/hg19] t(1;8)(q23.3;q12.2)", t.short_form);
    s.require("and its detailed system, exactly: t(1;8)(1pter→1q23.3(164,789,896)::8q12.2(61,760,650)→8qter;8pter→8q12.2(61,760,649)::1q23.3(164,789,897)→1qter)",
              t.detailed == "seq[GRCh37/hg19] t(1;8)(1pter\xE2\x86\x92" "1q23.3(164,789,896)::8q12.2(61,760,650)\xE2\x86\x92" "8qter;8pter\xE2\x86\x92" "8q12.2(61,760,649)::1q23.3(164,789,897)\xE2\x86\x92" "1qter)",
              t.detailed);
    s.require("its junctions as Ordulu et al.'s rearrangement sequences: 1q23.3(+)(164,789,896)::8q12.2(+)(61,760,650) and 8q12.2(+)(61,760,649)::1q23.3(+)(164,789,897)",
              t.junctions.size() == 2 && t.junctions[0] == "1q23.3(+)(164,789,896)::8q12.2(+)(61,760,650)" && t.junctions[1] == "8q12.2(+)(61,760,649)::1q23.3(+)(164,789,897)",
              t.junctions.empty() ? "" : t.junctions[0]);
    // t(5;10)(p13.3;q21.3) with a 6 bp deletion on 10 and a 13 bp deletion on 5: der(5) = 10qter→10q21.3(67,539,997)::
    // 5p13.3(29,658,440)→5qter, der(10) = 10pter→10q21.3(67,539,990)::5p13.3(29,658,426)→5pter.
    const SeqDescription t5 = describe_seq({{9, 67539996, false, 4, 29658439, true}, {9, 67539990, true, 4, 29658426, false}}, g37, o);
    s.require("Ordulu et al.'s t(5;10) with its small deletions at both junctions: the short and detailed systems exactly",
              t5.short_form == "seq[GRCh37/hg19] t(5;10)(p13.3;q21.3)" &&
                  t5.detailed == "seq[GRCh37/hg19] t(5;10)(10qter\xE2\x86\x92" "10q21.3(67,539,997)::5p13.3(29,658,440)\xE2\x86\x92" "5qter;10pter\xE2\x86\x92"
                                 "10q21.3(67,539,990)::5p13.3(29,658,426)\xE2\x86\x92" "5pter)" &&
                  t5.acentric.size() == 2,
              t5.short_form + " | " + t5.detailed);
    const SeqDescription del = describe_seq({{11, 21340000, true, 11, 23500000, true}}, g37, o);
    s.require("Ordulu et al.'s deletion: seq[GRCh37/hg19] del(12)(p12.1p12.1), detailed del(12)(pter→p12.1(21,340,000)::p12.1(23,500,001)→qter)",
              del.kind == "del" && del.short_form == "seq[GRCh37/hg19] del(12)(p12.1p12.1)" &&
                  del.detailed == "seq[GRCh37/hg19] del(12)(pter\xE2\x86\x92p12.1(21,340,000)::p12.1(23,500,001)\xE2\x86\x92qter)",
              del.short_form + " | " + del.detailed);
  }

  // ---- seq: other patterns, GRCh38 -------------------------------------------------------------------------------
  {
    SeqOptions o38;
    const SeqDescription inv = describe_seq({{2, 50000000, true, 2, 120000000, false}, {2, 50000000, false, 2, 120000000, true}}, g38, o38);
    s.require("two junctions reversing 3:50,000,000-120,000,000 are a pericentric inversion, inv(3)(p21.31q13.33)",
              inv.kind == "inv" && inv.short_form == "seq[GRCh38] inv(3)(p21.31q13.33)", inv.short_form + " | " + inv.detailed);
    const SeqDescription dup = describe_seq({{6, 120000000, true, 6, 110000000, true}}, g38, o38);
    s.require("a junction from 7:120,000,000 back to 110,000,000 is a tandem duplication, dup(7)(q31.1q31.31)", dup.kind == "dup" && dup.short_form == "seq[GRCh38] dup(7)(q31.1q31.31)",
              dup.short_form);
    const SeqDescription der = describe_seq({{8, 130000000, true, 21, 23200000, true}}, g38, o38);
    s.require("one junction from 9 to 22, its partner absent, is an unbalanced derivative der(9)t(9;22)(q34.11;q11.23)",
              der.kind == "der" && der.short_form == "seq[GRCh38] der(9)t(9;22)(q34.11;q11.23)", der.short_form + " | " + der.detailed);
    // A translocation of three chromosomes, each break shared by two junctions: 2q, 3q and 8q exchanged in a cycle.
    const SeqDescription cx = describe_seq({{1, 100000000, true, 7, 100000000, true}, {2, 150000000, true, 1, 100000000, true}, {7, 100000000, true, 2, 150000000, true}}, g38, o38);
    s.require("the three junctions of a translocation of three chromosomes are one rearrangement, written (2,3,8)cx with its three derivatives",
              cx.kind == "cx" && cx.short_form == "seq[GRCh38] (2,3,8)cx" && cx.derivatives.size() == 3,
              cx.short_form + " | " + cx.detailed);
    const SeqDescription back = describe_seq({{6, 110000000, false, 6, 120000000, false}}, g38, o38);
    s.require("the same tandem duplication's junction read from its other end is the same dup", back.kind == "dup" && back.short_form == dup.short_form, back.short_form);
    const SeqDescription both = describe_seq({{0, 164789896, true, 7, 61760649, true}, {11, 21340000, true, 11, 23500000, true}, {7, 61760649, true, 0, 164789896, true}}, g37, o);
    s.require("a sample's junctions grouped into its rearrangements and described in ISCN's order: t(1;8)(q23.3;q12.2),del(12)(p12.1p12.1)",
              both.short_form == "seq[GRCh37/hg19] t(1;8)(q23.3;q12.2),del(12)(p12.1p12.1)" && both.kind == "t,del", both.short_form);
    const SeqDescription none = describe_seq({}, g38, o38);
    s.require("no junction, no description", none.kind.empty() && none.short_form.empty());
  }

  // ---- the circle closed ---------------------------------------------------------------------------------------------
  {
    SeqOptions o38;
    o38.unicode = false;
    const SeqDescription del = describe_seq({{11, 21340000, true, 11, 23500000, true}}, g38, o38);
    const Parsed p = parse(del.detailed);
    const Interpretation in = interpret(p.description, g38);
    const bool exact = p.ok() && in.techniques.size() == 1 && in.techniques[0].has_line && in.techniques[0].line.changes.size() == 1 &&
                       in.techniques[0].line.changes[0].start == 21340000 && in.techniques[0].line.changes[0].end == 23500000 &&
                       in.techniques[0].line.changes[0].copies == 1 && in.techniques[0].line.changes[0].inner_start == 21340000;
    s.require("a deletion written in the detailed system, read back and interpreted, loses exactly 12:21,340,001-23,500,000", exact, del.detailed + " " + p.errors());
    const SeqDescription t = describe_seq({{0, 164789896, true, 7, 61760649, true}, {7, 61760649, true, 0, 164789896, true}}, g38, o38);
    const Parsed q = parse(t.detailed);
    const Interpretation iq = interpret(q.description, g38);
    s.require("a translocation written in the detailed system, read back and interpreted, is balanced, with its two derivatives",
              q.ok() && iq.techniques.size() == 1 && iq.techniques[0].line.balanced && iq.techniques[0].line.fusions.size() == 2, t.detailed + " " + q.errors());
    const std::string arr = describe_arr({{0, 147105903, 147917509, 1}}, g38);
    const Interpretation ia = interpret(parse(arr).description, g38);
    s.require("an arr description read back gives its segment and holds its bands", ia.techniques.size() == 1 && ia.techniques[0].changes.size() == 1 &&
                                                                                        ia.techniques[0].changes[0].start == 147105903 && ia.techniques[0].diagnostics.empty());
    const Interpretation bad = interpret(parse("arr[GRCh38] 1p36.33(147,105,904_147,917,509)x1").description, g38);
    s.require("an arr description whose bands do not hold its coordinates is warned of", bad.techniques.size() == 1 && bad.techniques[0].diagnostics.size() == 1 &&
                                                                                             bad.techniques[0].diagnostics[0].message.find("1q21.1q21.2") != std::string::npos);
  }
  return s.done();
}
