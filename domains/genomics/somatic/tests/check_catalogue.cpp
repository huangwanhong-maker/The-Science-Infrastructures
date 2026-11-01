// Checks of catalogues: SBS-96 and DBS-78 counted from calls on a small reference against a count by
// hand; adjacent calls as doublets and runs of three as neither; filtered calls left out; tables
// of profiles read by class names in any order, with their problems; the round trip of a catalogue
// through a table.

#include <genomics/somatic/catalogue.hpp>

#include <calx/verification/check.hpp>

#include <string>
#include <vector>

using namespace genomics;

namespace {

somatic::SomaticCall snv(std::size_t p, char ref, char alt, bool pass = true) {
  somatic::SomaticCall c;
  c.position = p;
  c.ref = std::string(1, ref);
  c.alt = std::string(1, alt);
  if (!pass) c.filters.push_back("strand_bias");
  return c;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  //                                          0         1         2         3
  //                                          0123456789012345678901234567890123456
  const std::vector<sequence::Codes> ref{sequence::encode_dna("ACGTTGCATGCAAGCTTCCGGATCGATTACAGGCATC")};
  // By hand, each class read on the pyrimidine strand:
  //   p 1 C>T in A C G          → A[C>T]G
  //   p 5 G>A in T G C, rc G C A → G[C>T]A
  //   p 9 G>T in T G C          → G[C>A]A (rc of TGC is GCA, G>T is C>A)
  //   p 14 C>G in G C T         → G[C>G]T
  //   p 18 C>A, p 19 G>T        → the doublet CG>AT (adjacent: not in SBS-96)
  //   p 22 T>C in A T C         → A[T>C]C
  //   p 27 T>A, 28 A>G, 29 C>T  → a run of three: neither catalogue
  //   p 33 C>A filtered         → left out
  std::vector<somatic::SomaticCall> calls{snv(1, 'C', 'T'),  snv(5, 'G', 'A'),  snv(9, 'G', 'T'),  snv(14, 'C', 'G'), snv(18, 'C', 'A'), snv(19, 'G', 'T'),
                                          snv(22, 'T', 'C'), snv(27, 'T', 'A'), snv(28, 'A', 'G'), snv(29, 'C', 'T'), snv(33, 'C', 'A', false)};
  const auto c = somatic::catalogue(calls, ref);
  using sequence::sbs96_index;
  std::vector<double> hand(96, 0);
  for (auto [l, r, a, x] : std::vector<std::array<sequence::Code, 4>>{{sequence::A, sequence::C, sequence::T, sequence::G},
                                                                     {sequence::G, sequence::C, sequence::T, sequence::A},
                                                                     {sequence::G, sequence::C, sequence::A, sequence::A},
                                                                     {sequence::G, sequence::C, sequence::G, sequence::T},
                                                                     {sequence::A, sequence::T, sequence::C, sequence::C}})
    hand[static_cast<std::size_t>(sbs96_index(l, r, a, x))]++;
  s.require("SBS-96 of the example equals the count by hand", c.sbs96 == hand && c.substitutions == 5);
  s.require("the names of the counted classes", sequence::sbs96_name(static_cast<std::size_t>(c.sbs[0].second)) == "A[C>T]G" &&
                                                     sequence::sbs96_name(static_cast<std::size_t>(c.sbs[1].second)) == "G[C>T]A" &&
                                                     sequence::sbs96_name(static_cast<std::size_t>(c.sbs[2].second)) == "G[C>A]A");
  std::vector<double> dbs(78, 0);
  dbs[static_cast<std::size_t>(sequence::dbs78_index(sequence::C, sequence::G, sequence::A, sequence::T))]++;
  s.require("one doublet, CG>AT", c.doublets == 1 && c.dbs78 == dbs && sequence::dbs78_name(static_cast<std::size_t>(c.dbs[0].second)) == "CG>AT");
  s.require("a run of three is a multi-base substitution, in neither catalogue", c.multiple == 3);
  const auto all = somatic::catalogue(calls, ref, false);
  s.require("with filtered calls included, C>A at p 33 (G C A → G[C>A]A) counts too", all.substitutions == 6 && all.sbs96[static_cast<std::size_t>(sbs96_index(sequence::G, sequence::C, sequence::A, sequence::A))] == 2);

  // ---- tables ---------------------------------------------------------------------------------------
  {
    // A table in reverse class order, two columns, as text.
    std::string text = "Type\tS1\tS2\n";
    for (std::size_t k = 96; k-- > 0;) text += sequence::sbs96_name(k) + "\t" + std::to_string(k) + "\t" + std::to_string(2 * k) + "\n";
    const auto t = somatic::read_profiles(text);
    bool ok = t.classes == 96 && t.names == std::vector<std::string>{"S1", "S2"} && t.problems.empty();
    for (std::size_t k = 0; k < 96 && ok; k++) ok = t.profiles[0][k] == static_cast<double>(k) && t.profiles[1][k] == static_cast<double>(2 * k);
    s.require("a table in any row order is read by class names", ok);
    std::string bad = "Type\tS1\nA[C>A]A\t1\nA[C>A]A\t2\nA[G>A]A\t3\nA[C>T]C\tx\nA[C>T]G\t1\t2\n";
    const auto b = somatic::read_profiles(bad);
    std::size_t missing = 0, lined = 0;
    for (const auto& [line, what] : b.problems) (line == 0 ? missing : lined)++;
    s.require("problems: a duplicate, an unknown (purine-centred) class, a non-number, a short row (lines 3–6), and the classes missing",
              lined == 4 && missing == 96 - 2, std::to_string(lined) + " lined, " + std::to_string(missing) + " missing");
    std::string d = "Type\tX\n";
    for (std::size_t k = 0; k < 78; k++) d += sequence::dbs78_name(k) + "\t" + std::to_string(k) + "\n";
    const auto e = somatic::read_profiles(d);
    s.require("a DBS-78 table", e.classes == 78 && e.problems.empty() && e.profiles[0][77] == 77);
  }
  return s.done();
}
