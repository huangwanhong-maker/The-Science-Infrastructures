// Checks of ISCN's grammar (grammar.hpp), its typed reading (iscn.hpp) and its canonical form
// (format.hpp): the examples of ISCN 2009 (as the Karyotype Ontology quotes them), of Ordulu et al.
// 2014, of ISCN 2020 (as Miller and Madan 2021 list its changes), of the GenQA assessment and of
// Conlin et al. 2025 accepted and read; the canonical form a fixed point; ISCN's order of
// abnormalities, of the chromosomes of translocations and of breakpoints imposed; ISCN 2016's order
// turned into 2020's; faults of syntax located with what could follow, and faults of sense named;
// a corpus of laboratory karyotypes read as far as it is ISCN.

#include <genomics/karyotype/format.hpp>
#include <genomics/karyotype/iscn.hpp>

#include <calx/verification/check.hpp>

#include <filesystem>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics::karyotype;

namespace {

std::vector<std::vector<std::string>> read_tsv(const std::filesystem::path& p) {
  std::ifstream f(p);
  std::vector<std::vector<std::string>> rows;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string c;
    while (std::getline(ss, c, '\t')) cols.push_back(c);
    rows.push_back(cols);
  }
  return rows;
}

std::string canonical(std::string_view text, Edition e = Edition::iscn2020) {
  const Parsed p = parse(text, e);
  return p.accepted() ? format(p.description) : "(rejected: " + p.errors() + ")";
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const std::filesystem::path fixtures(GENOMICS_FIXTURES);

  // ---- the grammar compiles, and is read without ambiguity of structure ---------------------------------
  s.require("ISCN's grammar compiles on morphe's BNF reader", iscn_bnf().rule_names().size() > 50, std::to_string(iscn_bnf().rule_names().size()) + " rules");

  // ---- examples of ISCN 2009, Ordulu 2014, ISCN 2020, GenQA, Conlin 2025 ---------------------------------------
  {
    const auto rows = read_tsv(fixtures / "karyotype-iscn2009.tsv");
    std::size_t accepted = 0, sound = 0, detailed = 0, stable = 0;
    std::string failed, unstable;
    for (const auto& r : rows) {
      if (r.size() < 3) continue;
      const Parsed p = parse(r[2], Edition::iscn2016);
      accepted += p.accepted();
      sound += p.ok();
      if (!p.ok() && failed.size() < 400) failed += r[2] + " (" + p.errors() + "); ";
      if (r.size() > 3 && !r[3].empty()) detailed += parse(r[3], Edition::iscn2016).ok();
      const std::string f = canonical(r[2], Edition::iscn2016);
      stable += canonical(f) == f;
      if (canonical(f) != f && unstable.size() < 400) unstable += f + " → " + canonical(f) + "; ";
    }
    std::size_t with_detail = 0;
    for (const auto& r : rows) with_detail += r.size() > 3 && !r[3].empty();
    s.require("all " + std::to_string(rows.size()) + " examples of ISCN 2009 (Karyotype Ontology) are accepted and mean something", accepted == rows.size() && sound == rows.size(), failed);
    s.require("all " + std::to_string(with_detail) + " of their detailed forms are accepted and mean something", detailed == with_detail);
    s.require("the canonical form of each is a fixed point of parsing and writing", stable == rows.size(), unstable);
  }
  {
    const auto rows = read_tsv(fixtures / "karyotype-ordulu2014.tsv");
    std::size_t ok = 0;
    std::string failed;
    for (const auto& r : rows)
      for (const auto& t : r) {
        if (t.find("{") != std::string::npos) continue;  // homologies in braces are not read (open question 3)
        const Parsed p = parse(t);
        ok += p.ok();
        if (!p.ok()) failed += t + " (" + p.errors() + "); ";
      }
    std::size_t total = 0;
    for (const auto& r : rows)
      for (const auto& t : r) total += t.find("{") == std::string::npos;
    s.require("the " + std::to_string(total) + " sequencing descriptions of Ordulu et al. (2014) without homologies are accepted, with their arr and ish parts",
              ok == total, failed);
  }
  {
    // ISCN 2020's new examples (Miller and Madan 2021), GenQA's assessment (Hastings, Tabiner and Deans,
    // ESHG 2022 poster P15.001.B), Conlin et al. (2025) Table 1, ISCN 2024's summary (Hastings et al. 2025).
    const char* examples[] = {
        "46,XX,inv(2)(p23p13)", "46,Y,del(X)(p21p11.4)", "46,XY,der(9)inv(9)(p23p13)del(9)(q22q23)", "47,XX,+7,der(7)t(1;7)(q12;p22)x2",
        "47,XY,+der(4)t(4;11)(q21;q23),t(4;11)(q21;q23)", "46,XX,der(1)del(1)(p34p22)ins(1;17)(p22;q11q25)", "47,XX,t(9;22;6)(q34;q11.2;p21),+der(22)t(9;22;6)",
        "46,XX,dup(1)(p34p31)", "46,XX,dup(1)(p31p34)", "46,XY,hsr(21)(q22)", "46,XX,ins(2)(p13q31q21)", "46,XX,ins(5;2)(p14;q32q22)", "46,XX,r(7)(p15q31)",
        "46,XX,r(20)(p13q13.3)", "46,XY,del(20)(q11.2-13.1q13.3)", "46,XX,del(4)(p15.2)", "81<3n>,XXX,+8", "mos 45,X[25]/47,XXX[12]/46,XX[13]",
        "47,XXYc[5]",
        "49~55,XY,-3,+6,-7,+10,+17,+18,+21c,+21,+r,+mar1,+mar2x2,3~9dmin[cp10]/47,XY,+21c[10]", "47,XX,+der(22)t(11;22)(q23;q11)dmat",
        "47,XY,t(4;11)(q21;q23),+der(4)t(4;11)[1]/48,idem,+19[9]/46,XY[10]", "arr[GRCh38] (X)x2,(Y)x1,(1-20,22)x3",
        "46,XY,der(6)t(3;6)(p13;q21)", "46,XY,der(6)(6pter→6q21::3p13→3pter)", "arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1",
        "arr[GRCh38] 1q21.1q21.2(147029796x2,147105904_147917509x1,148354661x2)", "46,XY,rec(2)dup(2q)inv(2)(p23q32)dpat",
        "47,XX,t(9;22)(q34;q11.2),+21c[15]/47,XX,+21c[5]", "arr (X,1-22)x2", "seq[GRCh38] t(9;10)(q34.3;q24.2)", "46,U", "seq (13)x3"};
    std::size_t ok = 0;
    std::string failed;
    for (const char* e : examples) {
      const Parsed p = parse(e);
      ok += p.ok();
      if (!p.ok()) failed += std::string(e) + " (" + p.errors() + "); ";
    }
    s.require("the examples of ISCN 2020 and 2024, GenQA's assessment and Conlin et al. are accepted and mean something", ok == std::size(examples), failed);
  }

  // ---- the canonical form --------------------------------------------------------------------------------
  {
    const char* fixed[] = {"46,XX", "47,XY,+21", "46,XY,t(9;22)(q34;q11.2)", "47,XX,+21c[15]/47,XX,t(9;22)(q34;q11.2),+21c[5]",
                           "45,X,-Y", "46,X,i(X)(q10)", "mos 47,XXY[10]/46,XY[20]", "chi 46,XX[20]//46,XY[10]", "45~48,XY,+8,+mar[cp10]",
                           "arr[GRCh38] 7q11.23(73,265,388_74,726,130)x1", "46,XX,del(5)(q13q33)", "46,XY,der(22)t(9;22)(q34;q11.2)",
                           "46,XX,ins(5;2)(p14;q32q22)", "47,XY,+mar1 c", "46,XX,inv dup(1)(q25q22)", "46,XY,del(20)(q11.2-13.1q13.3)"};
    std::size_t same = 0;
    std::string differ;
    for (const char* f : fixed) {
      const std::string c = canonical(f);
      same += c == f;
      if (c != f) differ += std::string(f) + " → " + c + "; ";
    }
    s.require("descriptions already in ISCN's form are written back unchanged", same == std::size(fixed), differ);
    s.require("abnormalities in ISCN's order: sex chromosomes, then autosomes by number; numerical before structural on one chromosome",
              canonical("48,XY,+21,del(5)(q13q33),+8,t(22;9)(q11.2;q34)") == "48,XY,del(5)(q13q33),+8,t(9;22)(q34;q11.2),+21",
              canonical("48,XY,+21,del(5)(q13q33),+8,t(22;9)(q11.2;q34)"));
    s.require("on one chromosome: +3 before inv(3), −10 before der(10), a constitutional +21c before an acquired +21",
              canonical("47,XX,inv(3)(q21q26.2),+3") == "47,XX,+3,inv(3)(q21q26.2)" &&
                  canonical("45,XY,der(10)t(10;17)(q22;p12),-10") == "45,XY,-10,der(10)t(10;17)(q22;p12)" &&
                  canonical("48,XY,+21,+21c") == "48,XY,+21c,+21");
    s.require("structural abnormalities of one chromosome alphabetically, the stemline's before an additional one: der before t, i(8) before +i(8)",
              canonical("47,XY,t(4;11)(q21;q23),+der(4)t(4;11)(q21;q23)") == "47,XY,+der(4)t(4;11)(q21;q23),t(4;11)(q21;q23)" &&
                  canonical("47,XX,+i(8)(q10),i(8)(q10)") == "47,XX,i(8)(q10),+i(8)(q10)");
    s.require("the sex chromosomes' abnormalities first, X before Y, unknown chromosomes and markers last",
              canonical("47,XY,+mar,t(Y;1)(q11.2;p13),inv(X)(p21q26),+der(?)t(?;3)(?;p21)") ==
                  "47,XY,inv(X)(p21q26),t(Y;1)(q11.2;p13),+der(?)t(?;3)(?;p21),+mar",
              canonical("47,XY,+mar,t(Y;1)(q11.2;p13),inv(X)(p21q26),+der(?)t(?;3)(?;p21)"));
    s.require("a translocation of three chromosomes rotated to begin with the lowest, its order kept: t(9;22;6) is t(6;9;22)",
              canonical("46,XX,t(9;22;6)(q34;q11.2;p21)") == "46,XX,t(6;9;22)(p21;q34;q11.2)", canonical("46,XX,t(9;22;6)(q34;q11.2;p21)"));
    s.require("the breakpoints of an inversion and of a deletion from pter to qter (ISCN 2020): inv(2)(p13p23) is inv(2)(p23p13)",
              canonical("46,XX,inv(2)(p13p23)") == "46,XX,inv(2)(p23p13)" && canonical("46,Y,del(X)(p11.4p21)") == "46,Y,del(X)(p21p11.4)" &&
                  canonical("46,XX,inv(3)(q26.2q21)") == "46,XX,inv(3)(q21q26.2)");
    s.require("the signs written in Unicode (×, →, −) read as their ASCII stand-ins",
              canonical("46,XX.arr[hg19] 12p12.1(21,340,001–23,500,000)×1") == "46,XX.arr[hg19] 12p12.1(21,340,001_23,500,000)x1" &&
                  canonical("45,XY,\xE2\x88\x92" "7") == "45,XY,-7");
  }

  // ---- ISCN 2016's order turned into 2020's ------------------------------------------------------------------
  {
    // Miller and Madan (2021): ISCN 2020 §9.2.5 and §9.2.9 list the bands from pter, ISCN 2016 from the
    // centromere; the examples dup(1)(p31p34) → dup(1)(p34p31) direct, ins(2)(p13q21q31) → ins(2)(p13q31q21),
    // ins(5;2)(p14;q22q32) → ins(5;2)(p14;q32q22), the q arm unchanged.
    s.require("ISCN 2016's direct dup(1)(p31p34) is ISCN 2020's dup(1)(p34p31); dup(1)(q22q25) is unchanged",
              canonical("46,XX,dup(1)(p31p34)", Edition::iscn2016) == "46,XX,dup(1)(p34p31)" && canonical("46,XX,dup(1)(q22q25)", Edition::iscn2016) == "46,XX,dup(1)(q22q25)");
    s.require("ISCN 2016's direct insertions into p, ins(2)(p13q21q31) and ins(5;2)(p14;q22q32), are ISCN 2020's ins(2)(p13q31q21) and ins(5;2)(p14;q32q22); into q unchanged",
              canonical("46,XX,ins(2)(p13q21q31)", Edition::iscn2016) == "46,XX,ins(2)(p13q31q21)" &&
                  canonical("46,XY,ins(5;2)(p14;q22q32)", Edition::iscn2016) == "46,XY,ins(5;2)(p14;q32q22)" &&
                  canonical("46,XX,ins(5;2)(q31;p13p23)", Edition::iscn2016) == "46,XX,ins(5;2)(q31;p13p23)");
  }

  // ---- the typed reading -----------------------------------------------------------------------------------------
  {
    const Parsed p = parse("49~55<2n>,XYc,+8,der(1)t(1;3)(p32;q21)t(1;11)(q25;q13)x2,+mar2x3,3~9dmin,inc[cp10]/50,idem,+9[5]");
    const Description& d = p.description;
    const CellLine& a = d.lines.at(0);
    s.require("a cell line's count range, ploidy, constitutional sex chromosomes, composite count of cells and inc are read",
              p.ok() && a.count == 49 && a.count_hi == 55 && a.ploidy == 2 && a.sex == "XY" && a.sex_constitutional && a.composite && a.cells == 10 && a.incomplete);
    const Abnormality& der = a.abnormalities.at(1);
    s.require("a derivative with its two components, its multiplier, and its chromosome 1",
              der.type == Type::der && der.components.size() == 2 && der.components[1].type == Type::t && der.multiplier == 2 && der.chromosome() == 0 &&
                  der.components[0].groups.at(1).at(0).label() == "q21");
    const Abnormality& mar = a.abnormalities.at(2);
    const Abnormality& dm = a.abnormalities.at(3);
    s.require("markers and double minutes: +mar2x3 is marker 2 three times, 3~9dmin three to nine double minutes",
              mar.type == Type::marker && mar.marker_number == 2 && mar.multiplier == 3 && dm.type == Type::dmin && dm.count == 3 && dm.count_hi == 9);
    s.require("idem refers to the stemline", d.lines.at(1).reference == CellLine::Reference::idem && d.lines.at(1).separator == "/");
    const Parsed q = parse("seq[GRCh37/hg19] t(1;8)(1pter→1q23.3(164,789,896)::8q12.2(61,760,650)→8qter;8pter→8q12.2(61,760,649)::1q23.3(164,789,897)→1qter)dn");
    const Abnormality& t = q.description.techniques.at(0).structural.at(0);
    const LongPart& seg = t.long_form.at(0).parts.at(0);
    s.require("a sequencing result's detailed form: two chromosomes, segments with nucleotide positions, its build and inheritance",
              q.ok() && q.description.techniques[0].build == "GRCh37/hg19" && t.long_form.size() == 2 && seg.to.has_position && seg.to.position == 164789896 &&
                  seg.from.kind == Breakpoint::Kind::pter && t.inheritance == "dn");
    const Parsed r = parse("arr[GRCh38] 2p25.3p11.2(204,193_86,768,824)x2~3[0.3],11p15.5p15.4(230,680_10,807,544)x2 hmz,(13)x3");
    const auto& e = r.description.techniques.at(0).arr;
    s.require("arr entries: a mosaic range of copies with its proportion, a homozygous run, a whole chromosome",
              r.ok() && e.size() == 3 && e[0].copies == 2 && e[0].copies_hi == 3 && e[0].proportion > 0.29 && e[0].start == 204193 && e[1].flags.at(0) == "hmz" &&
                  e[2].kind == ArrEntry::Kind::whole && e[2].wholes.at(0).first == 12 && e[2].copies == 3);
  }

  // ---- faults ---------------------------------------------------------------------------------------------------------
  {
    const Parsed a = parse("46,XX,dul(1)(q21q31)");
    s.require("a misspelt abbreviation: stopped after \"46,XX,du\" at column 9, 'p' (of 'dup') expected",
              !a.accepted() && a.syntax.furthest == 8 && a.errors().find("column 9") != std::string::npos && a.errors().find("'p' (of 'dup')") != std::string::npos,
              a.errors());
    const Parsed b = parse("46,XY,t(4;11)(q21;q23");
    s.require("an unclosed parenthesis: the text ends, and ')' could have followed", !b.accepted() && b.errors().find("ends") != std::string::npos && b.errors().find("')'") != std::string::npos, b.errors());
    const Parsed c = parse("46,XY,t(9;22)(q34)");
    s.require("a translocation of two chromosomes with one group of breakpoints is accepted by the grammar and refused by sense",
              c.accepted() && !c.ok() && c.errors().find("2 chromosomes but 1 group") != std::string::npos, c.errors());
    const Parsed d = parse("47,XX,+23");
    s.require("chromosome 23 is named as not a human chromosome", d.accepted() && !d.ok() && d.errors().find("no chromosome 23") != std::string::npos, d.errors());
    const Parsed e = parse("46,XX,del(5)(q13q21q31)");
    s.require("a deletion with three breakpoints is refused", !e.ok() && e.errors().find("one breakpoint (terminal) or two") != std::string::npos, e.errors());
    const Parsed f = parse("46,XX,del(5)(q13)t(5;9)(q31;p13)");
    s.require("abnormalities after a deletion are refused: only a derivative is made of others", !f.ok() && f.errors().find("only a derivative") != std::string::npos, f.errors());
    const Parsed g = parse("47,XY,+12[5]/48,sdl2,+8[3]");
    s.require("sdl2 with no second sideline before it is refused", !g.ok() && g.errors().find("sdl2") != std::string::npos, g.errors());
    const Parsed h = parse("46,XY,ins(5;2)(p14;q22)");
    s.require("an insertion without the two ends of its segment is refused", !h.ok() && h.errors().find("insertion") != std::string::npos, h.errors());
  }

  // ---- a corpus of laboratory karyotypes --------------------------------------------------------------------------------
  {
    std::ifstream f(fixtures / "karyotype-cll-corpus.txt");
    std::string line;
    std::size_t n = 0, raw = 0, ok = 0, stable = 0;
    std::map<std::string, std::size_t> why;
    while (std::getline(f, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      n++;
      raw += parse(line, Edition::iscn2016).ok();
      // RCytoGPS's own pre-cleaning: the cell lines a laboratory calls nonclonal are not ISCN.
      std::string clean;
      std::stringstream parts(line);
      std::string part;
      while (std::getline(parts, part, '/'))
        if (part.rfind("nonclonal", 0) != 0) clean += (clean.empty() ? "" : "/") + part;
      const Parsed p = parse(clean, Edition::iscn2016);
      if (!p.ok()) {
        const std::string e = p.errors();
        why[e.substr(0, e.find(';')).substr(0, 60)]++;
        continue;
      }
      ok++;
      const std::string c = format(p.description);
      stable += canonical(c) == c;
    }
    const double share = n ? static_cast<double>(ok) / static_cast<double>(n) : 0;
    std::string reasons;
    for (const auto& [k, v] : why) reasons += std::to_string(v) + "× " + k + "; ";
    s.require("of " + std::to_string(n) + " CLL karyotypes as laboratories wrote them (RCytoGPS's corpus), " + std::to_string(raw) + " are ISCN as written and " +
                  std::to_string(ok) + " (" + std::to_string(static_cast<int>(share * 100 + 0.5)) + "%) once their nonclonal cell lines are set aside, as RCytoGPS's pre-cleaning does",
              n > 500 && share > 0.5, reasons.substr(0, 1500));
    s.require("every one of them read is written back to a fixed point", stable == ok);
  }
  return s.done();
}
