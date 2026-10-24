// Patterns of DNA written in IUPAC letters (R, Y, W, N, …), every occurrence on both strands of
// a sequence found in one pass by one Aho–Corasick automaton (morphe::formal::AhoCorasick; Aho
// and Corasick 1975, Communications of the ACM 18, 333–340): restriction sites, promoter boxes,
// any small set of motifs.
//
// - A pattern stands for the set of words its letters allow: each letter the set of bases it
//   names (Cornish-Bowden 1985, Nucleic Acids Research 13, 3021–3030; sequence/alphabet.hpp), the
//   words the product of the sets, |W(p)| = Π_i |set(p_i)|. GAATTC is one word, GGWCC two, CCNGG
//   four. The automaton holds every word, so patterns standing for more than `max_words` words in
//   all are refused (ten Ns are about a million).
// - The reverse strand. A window x of the forward strand is an occurrence of p on the reverse
//   strand when its reverse complement is a word of p, that is when x ∈ rc(W(p)). The words of
//   W(p) are entered as forward words of p and those of rc(W(p)) \ W(p) as its reverse words, in
//   one automaton, so that one pass over the forward strand finds both and a window is reported
//   at most once for each pattern: on the forward strand when it is a word of p, even when its
//   reverse complement is one too. A palindromic site (GAATTC, or GGWCC, whose set is closed
//   under reverse complementation) is therefore counted once, as a restriction map counts it.
// - Coordinates are the forward strand's, [begin, end). A window holding N matches no pattern,
//   whatever its letters: N in a sequence is an unknown base, not every base (N is outside the
//   automaton's alphabet and returns it to its root).
//
// Construction O(4 Σ |w|) over the words entered; a pass O(n + occurrences).

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/formal/aho_corasick.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::index {

using sequence::Code;
using sequence::Codes;

// The words a pattern of IUPAC letters stands for, in lexicographic order of their codes.
// std::invalid_argument for an empty pattern, a letter that is not IUPAC's, or more than
// max_words words.
inline std::vector<Codes> iupac_words(std::string_view pattern, std::size_t max_words = std::size_t{1} << 20) {
  if (pattern.empty()) throw std::invalid_argument("genomics::index::iupac_words: an empty pattern");
  std::size_t count = 1;
  for (char c : pattern) {
    const std::uint8_t m = sequence::iupac_mask(c);
    if (m == 0) throw std::invalid_argument(std::string("genomics::index::iupac_words: '") + c + "' is not an IUPAC letter of DNA");
    count *= static_cast<std::size_t>(std::popcount(m));
    if (count > max_words) throw std::invalid_argument("genomics::index::iupac_words: " + std::string(pattern) + " stands for more than " + std::to_string(max_words) + " words");
  }
  std::vector<Codes> out{Codes{}};
  for (char c : pattern) {
    const std::uint8_t m = sequence::iupac_mask(c);
    std::vector<Codes> next;
    next.reserve(out.size() * static_cast<std::size_t>(std::popcount(m)));
    for (const auto& w : out)
      for (Code b = 0; b < 4; b++)
        if ((m >> b) & 1) {
          next.push_back(w);
          next.back().push_back(b);
        }
    out = std::move(next);
  }
  return out;  // lexicographic: each step extends the words in order by the bases in order
}

struct PatternHit {
  std::size_t begin = 0, end = 0;  // forward-strand coordinates
  std::uint32_t pattern = 0;       // its index among the patterns
  bool reverse = false;            // the window's reverse complement is the word
  friend bool operator==(const PatternHit&, const PatternHit&) = default;
};

class PatternSet {
 public:
  // The patterns, each in IUPAC letters (case is ignored); both_strands adds the reverse words.
  explicit PatternSet(const std::vector<std::string>& patterns, bool both_strands = true, std::size_t max_words = std::size_t{1} << 20)
      : patterns_(patterns), both_(both_strands) {
    std::vector<Codes> words;
    std::size_t total = 0;
    forward_.assign(patterns.size(), 0);
    for (std::uint32_t p = 0; p < patterns.size(); p++) {
      const auto w = iupac_words(patterns[p], max_words);
      total += w.size() * (both_strands ? 2 : 1);
      if (total > max_words) throw std::invalid_argument("genomics::index::PatternSet: the patterns stand for more than " + std::to_string(max_words) + " words");
      forward_[p] = w.size();
      for (const auto& x : w) {
        words.push_back(x);
        tags_.push_back({p, false});
      }
      if (!both_strands) continue;
      const std::set<Codes> own(w.begin(), w.end());
      for (const auto& x : w) {
        Codes r = sequence::reverse_complement(x);
        if (own.count(r)) continue;  // a word of p itself: its windows are forward occurrences
        words.push_back(std::move(r));
        tags_.push_back({p, true});
      }
    }
    automaton_ = morphe::formal::AhoCorasick(words, 4);
  }

  std::size_t patterns() const { return patterns_.size(); }
  const std::string& pattern(std::size_t p) const { return patterns_.at(p); }
  // The words a pattern stands for, and the words in the automaton in all (both strands).
  std::size_t forward_words(std::size_t p) const { return forward_.at(p); }
  std::size_t words() const { return tags_.size(); }
  bool both_strands() const { return both_; }

  // Every occurrence in a sequence of DNA codes, in order of position, then of pattern.
  std::vector<PatternHit> find(const Codes& text) const {
    std::vector<PatternHit> out;
    for (const auto& m : automaton_.matches(text)) {
      const Tag& t = tags_[m.pattern];
      out.push_back({m.begin, m.end, t.pattern, t.reverse});
    }
    std::sort(out.begin(), out.end(), [](const PatternHit& a, const PatternHit& b) { return a.begin != b.begin ? a.begin < b.begin : a.pattern < b.pattern; });
    return out;
  }

 private:
  struct Tag {
    std::uint32_t pattern = 0;
    bool reverse = false;
  };
  std::vector<std::string> patterns_;
  bool both_ = true;
  std::vector<std::size_t> forward_;
  std::vector<Tag> tags_;
  morphe::formal::AhoCorasick automaton_;
};

}  // namespace genomics::index
