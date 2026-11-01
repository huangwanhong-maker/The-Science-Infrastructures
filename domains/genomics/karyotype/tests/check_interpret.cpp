// Checks of the interpretation (interpret.hpp) on GRCh38's 850-band ideogram: for every kind of
// abnormality the derivative chromosomes as words of oriented segments, the copy number at every
// band, the losses and gains from where to where with their extents, the fusions, the count of
// chromosomes; the net imbalances ISCN 2009 states for its examples (pp. 80-82, as quoted by the
// Karyotype Ontology) and ISCN 2020 for i(X)(q10); the detailed system written from the short
// equal to ISCN 2009's own, for its 27 examples that give one; invariants (a balanced
// rearrangement keeps the multiset of segments, the count agrees with the chromosomes made).

#include <genomics/karyotype/interpret.hpp>

#include <calx/verification/check.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics::karyotype;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

struct Result {
  Parsed parsed;
  Interpretation in;
  const LineInterpretation& line(std::size_t i = 0) const { return in.lines.at(i); }
};

Result run(const std::string& text, const Ideogram& g, Edition e = Edition::iscn2020) {
  Result r;
  r.parsed = parse(text, e);
  r.in = interpret(r.parsed.description, g);
  return r;
}

// The changes of a line as "5q13→5q33 x1", in their order along the genome.
std::string changes(const LineInterpretation& L) {
  std::string s;
  for (const auto& c : L.changes)
    s += (s.empty() ? "" : ", ") + chromosome_name(c.chromosome) + c.from_band + "→" + chromosome_name(c.chromosome) + c.to_band + " x" + std::to_string(c.copies);
  return s;
}

std::vector<std::string> abnormal_long_forms(const LineInterpretation& L) {
  std::vector<std::string> out;
  LongFormOptions o;
  o.unicode = false;
  for (const auto& ch : L.chromosomes)
    if (ch.abnormal && !ch.marker) out.push_back(ch.label + ": " + long_form(ch, o));
  return out;
}

std::string joined(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& x : v) s += (s.empty() ? "" : " | ") + x;
  return s;
}

bool diagnostics_clean(const Result& r) {
  if (!r.parsed.ok()) return false;
  for (const auto& L : r.in.lines)
    for (const auto& d : L.diagnostics)
      if (d.level != Diagnostic::Level::note) return false;
  return true;
}

std::string diagnostics(const Result& r) {
  std::string s = r.parsed.errors();
  for (const auto& L : r.in.lines)
    for (const auto& d : L.diagnostics) s += (s.empty() ? "" : "; ") + std::string(level_name(d.level)) + ": " + d.message;
  return s;
}

// The copies at a band, as "min" or "min-max".
std::string copies_at(const LineInterpretation& L, const Ideogram& g, const std::string& band) {
  const int c = chromosome_index(band.substr(0, band.find_first_of("pq")));
  const auto r = g.find(c, band.substr(band.find_first_of("pq")));
  if (!r) return "?";
  int lo = 1 << 30, hi = -1;
  for (std::size_t b = r->first; b <= r->last; b++) { lo = std::min(lo, L.bands[b].min); hi = std::max(hi, L.bands[b].max); }
  return lo == hi ? std::to_string(lo) : std::to_string(lo) + "-" + std::to_string(hi);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-ideogram";
  const std::filesystem::path path = cache / "ideogram_9606_GCF_000001305.16_850_V1.tsv";
  if (!std::filesystem::exists(path)) {
    s.skip("the ideogram is not in the data cache (python data/acquisition/fetch.py ncbi-ideogram): the interpretation is not checked");
    return s.done();
  }
  const Ideogram g = read_ideogram(read_file(path), 850);

  // ---- normal and numerical ------------------------------------------------------------------------
  {
    const Result r = run("46,XX", g);
    const auto& L = r.line();
    bool flat = true;
    for (const auto& b : L.bands) flat = flat && b.min == b.max && b.min == b.normal;
    s.require("46,XX: 46 chromosomes, every band at its normal copy number (2, and 0 on Y), balanced, no fusion",
              L.counted == 46 && L.changes.empty() && L.balanced && flat && L.fusions.empty() && copies_at(L, g, "Yq11.23") == "0");
  }
  {
    const Result r = run("47,XX,+21", g);
    const auto& L = r.line();
    s.require("47,XX,+21: chromosome 21 gained whole, every one of its bands at 3 copies, from 21p13 to 21q22.3",
              L.counted == 47 && changes(L) == "21p13→21q22.3 x3" && copies_at(L, g, "21p") == "3" && copies_at(L, g, "21q") == "3" && L.lgf.gain.size() == g.end_of(20) - g.begin_of(20),
              changes(L));
  }
  {
    s.require("45,X: one X against a female's two, the whole X lost", changes(run("45,X", g).line()) == "Xp22.33→Xq28 x1");
    const auto xxy = run("47,XXY", g);
    s.require("47,XXY: two X against a male's one, a gain of X", changes(xxy.line()) == "Xp22.33→Xq28 x2" && xxy.line().male && xxy.line().counted == 47);
    const auto xy = run("45,X,-Y", g);
    s.require("45,X,-Y: the loss of Y in a male, the X written already present; 45 chromosomes",
              changes(xy.line()) == "Yp11.32→Yq12 x0" && xy.line().counted == 45 && diagnostics_clean(xy), diagnostics(xy));
    const auto kl = run("46,XXYc,-X", g);
    s.require("46,XXYc,-X (ISCN 2009 p. 57): Klinefelter cells that lost an X hold XY, normal for a male; 46 chromosomes",
              kl.line().changes.empty() && kl.line().counted == 46, changes(kl.line()));
    const auto tr = run("71,XXX,+8,+10", g);
    s.require("71,XXX,+8,+10 (ISCN 2009 p. 55): near-triploid, four copies of 8 and 10 against three",
              tr.line().ploidy == 3 && changes(tr.line()) == "8p23.3→8q24.3 x4, 10p15.3→10q26.3 x4" && tr.line().counted == 71, changes(tr.line()));
    const auto hap = run("26,X,+4,+6,+21", g);
    s.require("26,X,+4,+6,+21 (ISCN 2009 p. 55): near-haploid, two copies of 4, 6 and 21 against one",
              hap.line().ploidy == 1 && hap.line().changes.size() == 3 && hap.line().counted == 26, changes(hap.line()));
  }

  // ---- deletions, duplications, inversions -------------------------------------------------------------------
  {
    const Result r = run("46,XX,del(5)(q13q33)", g);
    const auto& L = r.line();
    const auto& c = L.changes.at(0);
    const auto q13 = g.find(4, "q13"), q33 = g.find(4, "q33");
    s.require("del(5)(q13q33): lost from 5q13 to 5q33, certainly between them, possibly from the start of 5q13 to the end of 5q33",
              changes(L) == "5q13→5q33 x1" && c.inner_start == q13->end && c.inner_end == q33->start && c.outer_start == q13->start && c.outer_end == q33->end,
              changes(L));
    s.require("del(5)(q13q33): the bands strictly between at 1 copy, 5q13 and 5q33 at 1-2, the rest at 2; one fusion, 5q13::5q33",
              copies_at(L, g, "5q21") == "1" && copies_at(L, g, "5q13") == "1-2" && copies_at(L, g, "5q33") == "1-2" && copies_at(L, g, "5q35") == "2" &&
                  L.fusions.size() == 1 && L.lgf.fusion_labels == std::vector<std::string>{"5q13", "5q33"},
              copies_at(L, g, "5q13") + " " + copies_at(L, g, "5q21"));
    s.require("del(5)(q13): a terminal deletion loses from 5q13 to the telomere, 5q35.3", changes(run("46,XX,del(5)(q13)", g).line()) == "5q13→5q35.3 x1");
    s.require("del(5)(p13): on the short arm, from the telomere 5p15.33 to 5p13", changes(run("46,XX,del(5)(p13)", g).line()) == "5p15.33→5p13 x1");
    const auto d = run("46,XX,dup(1)(q22q25)", g);
    s.require("dup(1)(q22q25): gained from 1q22 to 1q25, written pter→q25::q22→qter", changes(d.line()) == "1q22→1q25 x3" &&
                                                                                          joined(abnormal_long_forms(d.line())) == "dup(1)(q22q25): 1pter->1q25::1q22->1qter");
    const auto inv = run("46,XX,inv(3)(p13q21)", g);
    s.require("inv(3)(p13q21): balanced, pter→p13::q21→p13::q21→qter, two fusions", inv.line().balanced && inv.line().fusions.size() == 2 &&
                                                                                        joined(abnormal_long_forms(inv.line())) == "inv(3)(p13q21): 3pter->3p13::3q21->3p13::3q21->3qter");
    const auto same = run("46,XX,del(5)(q13q13)", g);
    s.require("del(5)(q13q13): two breaks in one band, a loss within 5q13 only, nothing certain", same.line().changes.size() == 1 &&
                                                                                                   same.line().changes[0].inner_start >= same.line().changes[0].inner_end &&
                                                                                                   copies_at(same.line(), g, "5q12") == "2" && copies_at(same.line(), g, "5q14") == "2");
  }

  // ---- translocations and derivatives ----------------------------------------------------------------------------
  {
    const Result r = run("46,XY,t(9;22)(q34;q11.2)", g);
    const auto& L = r.line();
    s.require("t(9;22)(q34;q11.2): balanced; der(9) = 9pter→9q34::22q11.2→22qter and der(22) = 22pter→22q11.2::9q34→9qter (ISCN; the Philadelphia chromosome)",
              L.balanced && L.counted == 46 && joined(abnormal_long_forms(L)) ==
                                                     "der(9)t(9;22)(q34;q11.2): 9pter->9q34::22q11.2->22qter | der(22)t(9;22)(q34;q11.2): 22pter->22q11.2::9q34->9qter",
              joined(abnormal_long_forms(L)));
    s.require("t(9;22): the two derivatives share their breaks, so each band holds two copies exactly; fusions at 9q34 and 22q11.2",
              copies_at(L, g, "9q34") == "2" && copies_at(L, g, "22q11.2") == "2" && L.fusions.size() == 2 &&
                  L.lgf.fusion_labels == std::vector<std::string>{"22q11.2", "9q34"} && L.lgf.loss.empty() && L.lgf.gain.empty());
    const auto ph = run("46,XY,der(22)t(9;22)(q34;q11.2)", g);
    s.require("der(22)t(9;22) alone: 9q34→qter gained, 22q11.2→qter lost", changes(ph.line()) == "9q34→9q34.3 x3, 22q11.2→22q13.33 x1", changes(ph.line()));
    const auto d6 = run("46,XY,der(6)t(3;6)(p13;q21)", g);
    s.require("der(6)t(3;6)(p13;q21) is der(6)(6pter→6q21::3p13→3pter), as Conlin et al. (2025) write it; 3p gained, 6q lost",
              joined(abnormal_long_forms(d6.line())) == "der(6)t(3;6)(p13;q21): 6pter->6q21::3p13->3pter" &&
                  changes(d6.line()) == "3p26.3→3p13 x3, 6q21→6q27 x1",
              joined(abnormal_long_forms(d6.line())) + " / " + changes(d6.line()));
    const auto d1 = run("46,XX,der(1)t(1;3)(p32;q21)", g);
    s.require("der(1)t(1;3)(p32;q21) is read from the end of 3q, its centric segment 1p32→1qter forward: 3qter→3q21::1p32→1qter",
              joined(abnormal_long_forms(d1.line())) == "der(1)t(1;3)(p32;q21): 3qter->3q21::1p32->1qter", joined(abnormal_long_forms(d1.line())));
    const auto three = run("46,XX,der(1)t(1;3)(p32;q21)t(1;11)(q25;q13)", g);
    s.require("der(1)t(1;3)(p32;q21)t(1;11)(q25;q13) (ISCN 2009 p. 63): both ends of 1 replaced, 3qter→3q21::1p32→1q25::11q13→11qter",
              joined(abnormal_long_forms(three.line())) == "der(1)t(1;3)(p32;q21)t(1;11)(q25;q13): 3qter->3q21::1p32->1q25::11q13->11qter",
              joined(abnormal_long_forms(three.line())));
    const auto chain = run("46,XY,der(1)t(1;3)(p32;q21)t(3;7)(q28;q11.2)", g);
    s.require("der(1)t(1;3)(p32;q21)t(3;7)(q28;q11.2) (ISCN 2009 p. 63): the second translocation acts on the 3 that the first attached",
              joined(abnormal_long_forms(chain.line())) == "der(1)t(1;3)(p32;q21)t(3;7)(q28;q11.2): 7qter->7q11.2::3q28->3q21::1p32->1qter",
              joined(abnormal_long_forms(chain.line())));
    const auto ider = run("46,XX,ider(22)(q10)t(9;22)(q34;q11.2)", g);
    s.require("ider(22)(q10)t(9;22) (ISCN 2009 p. 64): the long arm of the Philadelphia chromosome doubled",
              joined(abnormal_long_forms(ider.line())) == "ider(22)(q10)t(9;22)(q34;q11.2): 9qter->9q34::22q11.2->22q10::22q10->22q11.2::9q34->9qter",
              joined(abnormal_long_forms(ider.line())));
    const auto plus = run("47,XX,t(9;22)(q34;q11.2),+der(22)t(9;22)", g);
    s.require("47,XX,t(9;22)(q34;q11.2),+der(22)t(9;22) (ISCN 2009 p. 66): the extra Ph's breakpoints taken from the translocation; 22pter→22q11.2 and 9q34→9qter gained",
              plus.line().counted == 47 && changes(plus.line()) == "9q34→9q34.3 x3, 22p13→22q11.2 x3", changes(plus.line()) + " " + diagnostics(plus));
    const auto m = run("45,XY,-10,der(10)t(10;17)(q22;p12)", g);
    s.require("45,XY,-10,der(10)t(10;17)(q22;p12) (ISCN 2009 p. 60): no normal 10 left, 45 chromosomes; 10q22→qter lost twice over, 17p gained",
              m.line().counted == 45 && changes(m.line()) == "10p15.3→10q22 x1, 10q22→10q26.3 x0, 17p13.3→17p12 x3" && diagnostics_clean(m),
              changes(m.line()) + " " + diagnostics(m));
    const auto way = run("46,XX,t(2;7;5)(p21;q22;q23)", g);
    s.require("t(2;7;5)(p21;q22;q23) (ISCN 2009 p. 78): balanced, each chromosome receiving the segment of the one before it",
              way.line().balanced && joined(abnormal_long_forms(way.line())) ==
                                         "der(2)t(2;7;5)(p21;q22;q23): 5qter->5q23::2p21->2qter | der(7)t(2;7;5)(p21;q22;q23): 7pter->7q22::2p21->2pter | "
                                         "der(5)t(2;7;5)(p21;q22;q23): 5pter->5q23::7q22->7qter",
              joined(abnormal_long_forms(way.line())));
  }

  // ---- whole arms: the net imbalances ISCN 2009 states (pp. 80-82) ---------------------------------------------------
  {
    struct Case { const char* text; const char* changes; int count; const char* source; };
    const Case cases[] = {
        {"45,XX,der(1;3)(p10;q10)", "1q11→1q44 x1, 3p26.3→3p11.1 x1", 45, "monosomy for 1q and for 3p"},
        {"46,XX,+1,der(1;3)(p10;q10)", "1p36.33→1p11.1 x3, 3p26.3→3p11.1 x1", 46, "trisomy 1p, monosomy 3p"},
        {"46,XX,der(1;3)(p10;q10),+3", "1q11→1q44 x1, 3q11.1→3q29 x3", 46, "monosomy 1q, trisomy 3q"},
        {"47,XX,+der(1;3)(p10;q10)", "1p36.33→1p11.1 x3, 3q11.1→3q29 x3", 47, "trisomy 1p and 3q"},
        {"44,XX,-1,der(1;3)(p10;q10)", "1p36.33→1p11.1 x1, 1q11→1q44 x0, 3p26.3→3p11.1 x1", 44, "nullisomy 1q, monosomy 1p and 3p"},
        {"45,XX,der(13;21)(q10;q10)", "13p13→13p11.1 x1, 21p13→21p11.1 x1", 45, "loss of the short arms of 13 and 21"},
        {"46,XX,der(13;21)(q10;q10),+21", "13p13→13p11.1 x1, 21q11.1→21q22.3 x3", 46, "loss of 13p, trisomy 21q"},
        {"46,XX,+21,der(21;21)(q10;q10)", "21p13→21p11.1 x1, 21q11.1→21q22.3 x3", 46, "two 21q in the derivative and one normal 21"},
        {"46,X,i(X)(q10)", "Xp22.33→Xp11.1 x1, Xq11.1→Xq28 x3", 46, "a single copy of Xp and three of Xq (ISCN 2020 §9.2.11)"},
    };
    for (const auto& c : cases) {
      const Result r = run(c.text, g);
      // A change from the centromere is written from the band beside it on its side (1q11, 3p11.1).
      const std::string got = changes(r.line());
      s.require(std::string(c.text) + ": " + c.source + " (" + std::to_string(c.count) + " chromosomes)",
                got == c.changes && r.line().counted == c.count && diagnostics_clean(r), got + " / " + std::to_string(r.line().counted) + " " + diagnostics(r));
    }
  }

  // ---- rings, isochromosomes, dicentrics, insertions, markers, recombinants ---------------------------------------------
  {
    const auto ring = run("46,XX,r(7)(p22q36)", g);
    s.require("r(7)(p22q36): pter→p22 and q36→qter lost, the ring closed at 7p22::7q36",
              changes(ring.line()) == "7p22.3→7p22 x1, 7q36→7q36.3 x1" && ring.line().fusions.size() == 1 && ring.line().fusions[0].ring &&
                  joined(abnormal_long_forms(ring.line())) == "r(7)(p22q36): ::7p22->7q36::",
              changes(ring.line()) + " " + joined(abnormal_long_forms(ring.line())));
    const auto i17 = run("46,XX,i(17)(q10)", g);
    s.require("i(17)(q10): 17p lost, 17q gained, qter→q10::q10→qter", changes(i17.line()) == "17p13.3→17p11.1 x1, 17q11.1→17q25.3 x3" &&
                                                                         joined(abnormal_long_forms(i17.line())) == "i(17)(q10): 17qter->17q10::17q10->17qter",
              changes(i17.line()));
    const auto dic = run("45,XX,dic(13;15)(q22;q24)", g);
    s.require("dic(13;15)(q22;q24): 13pter→13q22::15q24→15pter; 13q22→qter and 15q24→qter lost; 45 chromosomes",
              joined(abnormal_long_forms(dic.line())) == "dic(13;15)(q22;q24): 13pter->13q22::15q24->15pter" && dic.line().counted == 45 &&
                  changes(dic.line()) == "13q22→13q34 x1, 15q24→15q26.3 x1",
              joined(abnormal_long_forms(dic.line())) + " " + changes(dic.line()));
    const auto idic = run("46,X,idic(Y)(q12)", g);
    s.require("idic(Y)(q12): Ypter→Yq12::Yq12→Ypter, two copies of Y from pter to q12", joined(abnormal_long_forms(idic.line())) == "idic(Y)(q12): Ypter->Yq12::Yq12->Ypter",
              joined(abnormal_long_forms(idic.line())));
    const auto ins = run("46,XY,ins(5;2)(p14;q32q22)", g);
    s.require("ins(5;2)(p14;q32q22) (ISCN 2020; 2009's ins(5;2)(p14;q22q32)): balanced, 5pter→5p14::2q32→2q22::5p14→5qter and 2pter→2q22::2q32→2qter",
              ins.line().balanced && joined(abnormal_long_forms(ins.line())) ==
                                          "der(5)ins(5;2)(p14;q32q22): 5pter->5p14::2q32->2q22::5p14->5qter | der(2)ins(5;2)(p14;q32q22): 2pter->2q22::2q32->2qter",
              joined(abnormal_long_forms(ins.line())));
    const auto mar = run("47,XY,+mar", g);
    s.require("47,XY,+mar: a marker counted, of unknown content, no change of copy number known", mar.line().counted == 47 && mar.line().unknown == 1 && mar.line().balanced);
    const auto dmin = run("46,XY,9~50dmin", g);
    s.require("9~50dmin: double minutes are not counted among the chromosomes", dmin.line().counted == 46 && diagnostics_clean(dmin), diagnostics(dmin));
    const auto rec = run("46,XY,rec(2)dup(2q)inv(2)(p23q32)dpat", g);
    s.require("rec(2)dup(2q)inv(2)(p23q32): the recombinant of a pericentric inversion, 2qter→2q32::2p23→2qter; 2q32→qter gained, 2pter→p23 lost",
              joined(abnormal_long_forms(rec.line())) == "rec(2)dup(2q)inv(2)(p23q32): 2qter->2q32::2p23->2qter" && changes(rec.line()) == "2p25.3→2p23 x1, 2q32→2q37.3 x3",
              joined(abnormal_long_forms(rec.line())) + " " + changes(rec.line()));
    const auto fis = run("47,XY,-10,+fis(10)(p10),+fis(10)(q10)", g);
    s.require("47,XY,-10,+fis(10)(p10),+fis(10)(q10): a centric fission keeps the copy number, 47 chromosomes", fis.line().balanced && fis.line().counted == 47);
  }

  // ---- clones -----------------------------------------------------------------------------------------------------------
  {
    const auto r = run("47,XY,+12[10]/49,idem,+18,+19[5]/46,XY[20]", g);
    s.require("idem: the second clone holds the stemline's +12 and its own +18 and +19, 49 chromosomes; the normal clone none",
              r.in.lines.size() == 3 && r.line(1).counted == 49 && r.line(1).changes.size() == 3 && r.line(2).changes.empty() && diagnostics_clean(r),
              changes(r.line(1)) + " " + diagnostics(r));
    const auto sl = run("45,XX,dic(1;20)(q32;q13.3)[5]/45,sl,+1,-dic(1;20)[3]", g);
    s.require("sl,−dic(1;20): the sideline lacks the stemline's dicentric but not the chromosomes it replaced: one 1 and 20, +1 restoring 1",
              sl.line(1).counted == 45 && changes(sl.line(1)) == "20p13→20q13.33 x1", changes(sl.line(1)) + " " + diagnostics(sl));
    const auto mos = run("mos 47,XXY[10]/46,XY[20]", g);
    s.require("mos 47,XXY[10]/46,XY[20]: two cell lines, a gain of X in the first", mos.in.lines.size() == 2 && changes(mos.line(0)) == "Xp22.33→Xq28 x2" && mos.line(1).changes.empty());
  }

  // ---- what does not add up --------------------------------------------------------------------------------------------
  {
    const auto r = run("46,XX,+21", g);
    bool warned = false;
    for (const auto& d : r.line().diagnostics) warned = warned || d.message.find("give 47 chromosomes, and the count is 46") != std::string::npos;
    s.require("46,XX,+21: the count written disagrees with the 47 chromosomes described, and is warned of", warned, diagnostics(r));
    const auto b = run("46,XX,del(5)(q37)", g);
    bool band = false;
    for (const auto& d : b.line().diagnostics) band = band || d.message.find("no band 5q37") != std::string::npos;
    s.require("del(5)(q37): no such band in the ideogram, said so", band, diagnostics(b));
  }

  // ---- the detailed system written from the short: ISCN 2009's own --------------------------------------------------------
  {
    std::ifstream f(std::filesystem::path(GENOMICS_FIXTURES) / "karyotype-iscn2009.tsv");
    std::string line;
    std::size_t n = 0, same = 0;
    std::string differ;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> cols;
      std::stringstream ss(line);
      std::string c;
      while (std::getline(ss, c, '\t')) cols.push_back(c);
      if (cols.size() < 4 || cols[3].empty()) continue;
      n++;
      std::string want;
      for (std::size_t i = 0; i < cols[3].size(); i++)
        if (!(cols[3][i] == ' ' && i > 0 && cols[3][i - 1] == ';')) want += cols[3][i];
      const Result r = run(cols[2], g, Edition::iscn2016);
      LongFormOptions o;
      o.unicode = false;
      const std::string got = r.in.lines.empty() ? "" : detailed(r.line(), r.parsed.description.lines.at(0), o);
      same += got == want;
      if (got != want) differ += "\n    " + cols[2] + ": " + got + " vs " + want;
    }
    s.require("the detailed system written from the short is ISCN 2009's own for its " + std::to_string(n) + " examples that give one", n >= 25 && same == n, differ);
  }

  // ---- invariants over every example of ISCN 2009 -------------------------------------------------------------------------
  {
    std::ifstream f(std::filesystem::path(GENOMICS_FIXTURES) / "karyotype-iscn2009.tsv");
    std::string line;
    // Examples whose quotation does not add up to its own count (a + or a sex chromosome lost in
    // the quotation, or in the book): reported, not counted against the interpretation.
    const std::set<std::string> known = {"47,XY,dic(15;15)(q12;q12)", "50,XX,+1,+3,+8,+r,+mar", "46,XY,t(X;15;18)(p11.1;p11.1;q11.1)",
                                         "50,XX,+1,del(1)(p13),+dup(1)(q21q32),+inv(1)(p31q41),+8,r(10)(p12q25),-21"};
    std::size_t n = 0, counted = 0, balanced_ok = 0, balanced_n = 0;
    std::string off;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> cols;
      std::stringstream ss(line);
      std::string c;
      while (std::getline(ss, c, '\t')) cols.push_back(c);
      if (cols.size() < 3) continue;
      const Result r = run(cols[2], g, Edition::iscn2016);
      if (!r.parsed.ok() || r.in.lines.empty()) continue;
      n++;
      bool agree = true;
      for (const auto& L : r.in.lines)
        for (const auto& d : L.diagnostics)
          if (d.message.find("chromosomes, and the count is") != std::string::npos) agree = false;
      if (!agree && known.count(cols[2])) { n--; continue; }
      counted += agree;
      if (!agree) off += cols[2] + "; ";
      // Rearrangements that only move material (t, inv, ins between two chromosomes, whole-arm t,
      // fis in pairs) keep the copy number.
      const auto& L = r.line();
      bool only_moves = !L.abnormalities.empty();
      for (const auto& a : L.abnormalities)
        only_moves = only_moves && a.sign == 0 && (a.type == Type::t || a.type == Type::inv || (a.type == Type::ins && a.chromosomes.size() == 2 && !a.chromosomes[1].unknown));
      if (only_moves) { balanced_n++; balanced_ok += L.balanced; }
    }
    s.require("the chromosomes described agree with the count written for all " + std::to_string(n) + " examples of ISCN 2009 whose quotation adds up (" +
                  std::to_string(known.size()) + " do not)",
              counted == n, off);
    s.require("every one of its " + std::to_string(balanced_n) + " examples made only of translocations, inversions and insertions is balanced", balanced_n > 10 && balanced_ok == balanced_n);
  }
  return s.done();
}
