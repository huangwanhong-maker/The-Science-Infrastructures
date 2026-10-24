// Exact repeats of a genome: maximal and supermaximal repeats and the longest repeated
// substring, on one strand or on both, by morphe's suffix array and LCP array.
//
// For a word x, a maximal pair is two occurrences of a factor α that cannot be extended to the
// left or to the right together (the letters before them differ, or one is at the start; the
// letters after them differ, or one is at the end); α is a maximal repeat when it has a maximal
// pair, and a supermaximal repeat when it is a maximal repeat contained in no other (Gusfield
// 1997, Algorithms on Strings, Trees, and Sequences, §7.12).
//
// - On the suffix array SA with its LCP array, the internal nodes of the suffix tree are the
//   lcp-intervals [i, j] with value ℓ > 0 (LCP[k] ≥ ℓ inside, = ℓ somewhere, < ℓ at both
//   borders), enumerated bottom-up with a stack in O(n) (Kasai et al. 2001, §5; Abouelhoda,
//   Kurtz and Ohlebusch 2004, J. Discrete Algorithms 2, 53–86, Algorithm 4.4). Such an interval
//   is right-maximal by construction, and α = x[SA[i] .. SA[i] + ℓ) is a maximal repeat exactly
//   when its node is left-diverse, when the letters x[SA[k] − 1] (the Burrows–Wheeler letters;
//   a suffix at the start counts as a letter of its own) are not all equal on [i, j] (Gusfield,
//   Theorem 7.12.1); with the runs of the BWT precomputed, that is a test in O(1).
// - It is supermaximal exactly when the node's children are all leaves (no lcp-interval nested
//   inside: LCP[k] = ℓ on the whole of (i, j]) with pairwise distinct left letters (Gusfield,
//   §7.12).
// - Both strands: the word x #₀ rc(x) #₁ (morphe::formal::concatenate, its separators unique
//   letters), whose repeats include the inverted ones, a copy on each strand. Each double-strand
//   repeat appears there as α and as rc(α) with the same places; it is reported once, as the
//   smaller of the two, its occurrences given on the forward strand with their strands.
// - A repeat containing N is not reported (an assembly's runs of N would be its longest repeats).
//
// Complexity: O(n) to build the suffix array (SA-IS) and the LCP array, and O(n) for the
// traversal, with O(n) words of memory (32-bit indices: genomes below 2³¹ bases on both
// strands); reporting costs the occurrences reported.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/formal/suffix_array.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace genomics::annotation {

using sequence::Code;
using sequence::Codes;

struct Occurrence {
  std::size_t position = 0;  // the first base, forward strand
  bool reverse = false;      // the repeat reads on the reverse strand there
  bool operator<(const Occurrence& o) const { return position != o.position ? position < o.position : reverse < o.reverse; }
  bool operator==(const Occurrence& o) const { return position == o.position && reverse == o.reverse; }
};

struct Repeat {
  std::size_t length = 0;
  std::vector<Occurrence> occurrences;  // every occurrence of the word, sorted
  bool supermaximal = false;
};

struct RepeatOptions {
  std::size_t min_length = 1;
  bool both_strands = false;
  bool supermaximal_only = false;
};

namespace detail {

// The text the suffix array is built on, and how to map its positions back.
struct RepeatText {
  morphe::formal::Word text;
  std::size_t n = 0;  // the genome's length
  bool both = false;
  std::vector<std::uint32_t> n_before;  // the number of N codes before each position (prefix sums)
  std::uint32_t n_code = 0;             // the code of N in text
  Occurrence map(std::size_t p, std::size_t length) const {
    if (!both || p < n) return {p, false};
    const std::size_t o = p - (n + 1);
    return {n - o - length, true};
  }
  bool has_n(std::size_t p, std::size_t length) const { return n_before[p + length] != n_before[p]; }
};

inline RepeatText repeat_text(const Codes& x, bool both) {
  RepeatText t;
  t.n = x.size();
  t.both = both;
  for (Code c : x)
    if (c > 4) throw std::invalid_argument("genomics::annotation::repeats: DNA codes 0 to 4");
  if (both) {
    const auto c = morphe::formal::concatenate(std::vector<Codes>{x, sequence::reverse_complement(x)}, 5);
    t.text = c.text;
    t.n_code = 4 + 2;
  } else {
    t.text.assign(x.begin(), x.end());
    t.n_code = 4;
  }
  t.n_before.assign(t.text.size() + 1, 0);
  for (std::size_t i = 0; i < t.text.size(); i++) t.n_before[i + 1] = t.n_before[i] + (t.text[i] == t.n_code);
  return t;
}

// Whether the word at p is not larger than its reverse complement (both strands).
inline bool canonical(const RepeatText& t, std::size_t p, std::size_t length) {
  const std::uint32_t shift = t.both ? 2 : 0;
  for (std::size_t i = 0; i < length; i++) {
    const std::uint32_t a = t.text[p + i] - shift, b = 3 - (t.text[p + length - 1 - i] - shift);
    if (a != b) return a < b;
  }
  return true;
}

}  // namespace detail

// The maximal repeats of x of at least min_length bases (with their supermaximality), longest
// first, then by first occurrence.
inline std::vector<Repeat> maximal_repeats(const Codes& x, const RepeatOptions& o = {}) {
  std::vector<Repeat> out;
  if (x.size() < 2) return out;
  const detail::RepeatText t = detail::repeat_text(x, o.both_strands);
  const std::size_t N = t.text.size();
  const std::vector<std::uint32_t> sa = morphe::formal::suffix_array<std::uint32_t>(t.text, o.both_strands ? 4 + 1 + 2 : 5);
  const std::vector<std::uint32_t> lcp = morphe::formal::lcp_array(t.text, sa);
  // The BWT letters, a suffix at the start getting a letter of its own, and the end of each run.
  constexpr std::uint32_t START = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> bwt(N), run_end(N);
  for (std::size_t k = 0; k < N; k++) bwt[k] = sa[k] ? t.text[sa[k] - 1] : START;
  run_end[N - 1] = static_cast<std::uint32_t>(N - 1);
  for (std::size_t k = N - 1; k-- > 0;) run_end[k] = bwt[k] == bwt[k + 1] && bwt[k] != START ? run_end[k + 1] : static_cast<std::uint32_t>(k);

  struct Interval {
    std::uint32_t lcp, lb;
    bool child;
  };
  std::vector<Interval> stack{{0, 0, false}};
  const auto report = [&](const Interval& v, std::size_t rb) {
    if (v.lcp < o.min_length || v.lcp == 0) return;
    if (run_end[v.lb] >= rb) return;  // not left-diverse
    const std::size_t p = sa[v.lb], len = v.lcp;
    if (t.has_n(p, len)) return;
    if (t.both && !detail::canonical(t, p, len)) return;
    bool super = !v.child;
    if (super)
      for (std::size_t a = v.lb; a <= rb && super; a++)
        for (std::size_t b = a + 1; b <= rb; b++)
          if (bwt[a] == bwt[b] && bwt[a] != START) { super = false; break; }
    if (o.supermaximal_only && !super) return;
    Repeat r;
    r.length = len;
    r.supermaximal = super;
    for (std::size_t k = v.lb; k <= rb; k++) r.occurrences.push_back(t.map(sa[k], len));
    std::sort(r.occurrences.begin(), r.occurrences.end());
    r.occurrences.erase(std::unique(r.occurrences.begin(), r.occurrences.end()), r.occurrences.end());
    out.push_back(std::move(r));
  };
  for (std::size_t i = 1; i <= N; i++) {
    const std::uint32_t h = i < N ? lcp[i] : 0;
    std::uint32_t lb = static_cast<std::uint32_t>(i - 1);
    bool pending = false;
    while (h < stack.back().lcp) {
      const Interval last = stack.back();
      stack.pop_back();
      report(last, i - 1);
      lb = last.lb;
      if (h <= stack.back().lcp) stack.back().child = true;
      else pending = true;
    }
    if (h > stack.back().lcp) stack.push_back({h, lb, pending});
  }
  std::sort(out.begin(), out.end(), [](const Repeat& a, const Repeat& b) {
    return a.length != b.length ? a.length > b.length : a.occurrences.front() < b.occurrences.front();
  });
  return out;
}

inline std::vector<Repeat> supermaximal_repeats(const Codes& x, RepeatOptions o = {}) {
  o.supermaximal_only = true;
  return maximal_repeats(x, o);
}

// The longest repeated substring (on both strands, an inverted repeat counts): its length and
// two of its occurrences; length 0 when no base repeats.
inline Repeat longest_repeat(const Codes& x, bool both_strands = false) {
  Repeat best;
  if (x.size() < 2) return best;
  const detail::RepeatText t = detail::repeat_text(x, both_strands);
  const std::vector<std::uint32_t> sa = morphe::formal::suffix_array<std::uint32_t>(t.text, both_strands ? 4 + 1 + 2 : 5);
  const std::vector<std::uint32_t> lcp = morphe::formal::lcp_array(t.text, sa);
  std::size_t at = 0;
  for (std::size_t k = 1; k < sa.size(); k++)
    if (lcp[k] > best.length && !t.has_n(sa[k], lcp[k])) {
      best.length = lcp[k];
      at = k;
    }
  if (best.length) {
    best.occurrences = {t.map(sa[at - 1], best.length), t.map(sa[at], best.length)};
    std::sort(best.occurrences.begin(), best.occurrences.end());
  }
  return best;
}

}  // namespace genomics::annotation
