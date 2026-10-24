// Checks of IUPAC patterns found by one Aho–Corasick automaton: the words a pattern stands for
// (their number the product of its letters' sets, in order, and the faults refused); and every
// occurrence of nine patterns on both strands of 300 kb of random sequence with runs of N, against
// their definition by brute force, position by position with iupac_matches: forward where the
// window reads the pattern, reverse where its reverse complement does and it is not forward. A
// palindromic site is found once; one strand alone gives the forward occurrences only.

#include <genomics/index/patterns.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

using namespace genomics;
using index::PatternHit;
using index::PatternSet;
using sequence::Code;
using sequence::Codes;

namespace {

constexpr std::uint32_t CHECK = 0x50415454;  // "PATT"

// Whether the window of x at i reads p on the forward strand, or on the reverse one.
bool reads_forward(const Codes& x, std::size_t i, const std::string& p) {
  for (std::size_t j = 0; j < p.size(); j++)
    if (!sequence::iupac_matches(p[j], x[i + j])) return false;
  return true;
}
bool reads_reverse(const Codes& x, std::size_t i, const std::string& p) {
  const std::size_t L = p.size();
  for (std::size_t j = 0; j < L; j++) {
    const Code b = x[i + L - 1 - j];
    if (b > 3 || !sequence::iupac_matches(p[j], sequence::complement(b))) return false;
  }
  return true;
}

std::vector<PatternHit> brute(const Codes& x, const std::vector<std::string>& patterns, bool both) {
  std::vector<PatternHit> out;
  for (std::size_t i = 0; i < x.size(); i++)
    for (std::uint32_t p = 0; p < patterns.size(); p++) {
      const std::size_t L = patterns[p].size();
      if (i + L > x.size()) continue;
      if (reads_forward(x, i, patterns[p])) out.push_back({i, i + L, p, false});
      else if (both && reads_reverse(x, i, patterns[p])) out.push_back({i, i + L, p, true});
    }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the words of a pattern ----------------------------------------------------------------------------------
  {
    const auto one = index::iupac_words("GAATTC"), two = index::iupac_words("ggwcc"), four = index::iupac_words("CCNGG");
    const auto rybn = index::iupac_words("RYBN");
    s.require("GAATTC is one word, GGWCC two (GGACC, GGTCC), CCNGG four, RYBN 2·2·3·4 = 48", one.size() == 1 && two.size() == 2 && four.size() == 4 && rybn.size() == 48 &&
              sequence::decode_dna(two[0]) == "GGACC" && sequence::decode_dna(two[1]) == "GGTCC");
    bool ordered = true, matching = true;
    for (std::size_t k = 1; k < rybn.size(); k++) ordered = ordered && rybn[k - 1] < rybn[k];
    for (const auto& w : rybn)
      for (std::size_t j = 0; j < 4; j++) matching = matching && sequence::iupac_matches("RYBN"[j], w[j]);
    s.require("the words are distinct, in order, and each reads its pattern letter by letter", ordered && matching);
    int refused = 0;
    for (const char* bad : {"", "GAXTC", "GA-TC"}) {
      try { index::iupac_words(bad); } catch (const std::invalid_argument&) { refused++; }
    }
    try { index::iupac_words("NNNNNNNNNNN", 1u << 20); } catch (const std::invalid_argument&) { refused++; }
    s.require("an empty pattern, a letter outside IUPAC's, and eleven Ns (4^11 words) beyond the bound are refused", refused == 4);
  }

  // ---- occurrences against their definition -------------------------------------------------------------------------
  {
    Codes x = simulate::random_sequence(300000, 0.45, 20261024, CHECK);
    calx::random::Stream r(20261024, CHECK, 1, 0);
    for (int run = 0; run < 40; run++) {  // runs of N, 1 to 20 bases
      const std::size_t at = r.below(x.size() - 20), n = 1 + r.below(20);
      for (std::size_t i = 0; i < n; i++) x[at + i] = sequence::N;
    }
    const std::vector<std::string> patterns = {"GAATTC", "GGATCC", "GGWCC", "CCNGG", "TATAWAW", "RGATCY", "AAGCTT", "ACGNNNT", "GCC"};
    const PatternSet both(patterns, true), one(patterns, false);
    const auto got = both.find(x), want = brute(x, patterns, true);
    s.require("nine patterns on both strands of 300 kb with runs of N: every occurrence, and only those, as brute force finds them", got == want,
              std::to_string(got.size()) + " and " + std::to_string(want.size()));
    std::size_t reverse = 0;
    for (const auto& h : got) reverse += h.reverse;
    s.require("both strands give reverse occurrences of the patterns that are not palindromes", reverse > 0);
    bool palindromes_once = true;
    for (const auto& h : got) palindromes_once = palindromes_once && !(h.reverse && (h.pattern == 0 || h.pattern == 1 || h.pattern == 2 || h.pattern == 3 || h.pattern == 6));
    s.require("GAATTC, GGATCC, GGWCC, CCNGG and AAGCTT, closed under reverse complementation, are found once a site, forward", palindromes_once);
    s.require("one strand alone gives exactly the forward occurrences", one.find(x) == brute(x, patterns, false));
    s.require("the automaton holds each pattern's words and the reverse words not among them",
              both.forward_words(2) == 2 && both.forward_words(3) == 4 && both.forward_words(4) == 4 && both.words() > one.words() && one.words() == 1 + 1 + 2 + 4 + 4 + 4 + 1 + 64 + 1);
  }

  return s.done();
}
