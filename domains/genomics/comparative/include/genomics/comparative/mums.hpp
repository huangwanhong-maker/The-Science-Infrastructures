// Maximal unique matches of two genomes, the anchors of whole-genome alignment (Delcher, Kasif,
// Fleischmann, Peterson, White and Salzberg 1999, "Alignment of whole genomes", Nucleic Acids
// Research 27, 2369–2376: MUMmer).
//
// - A maximal unique match (MUM) of genomes A and B is a word that occurs exactly once in A and
//   exactly once in B and cannot be extended at either end, the bases before (and after) its two
//   occurrences differing or one of them missing (ibid.; that a unique match contained in a longer
//   one is always an extension of it is the theory's Proposition 4.2).
//   On the reverse strand, a MUM of A and the reverse complement of B. Bases N match nothing:
//   each genome is cut into its runs free of N, and no match crosses one.
// - The suffix array of the runs of both genomes, each run ending in a separator of its own (morphe's
//   `concatenate` and SA-IS, with Kasai's LCP array: morphe/formal/suffix_array.hpp). A word occurs
//   exactly once in A and once in B when exactly two suffixes, one from each, begin with it; those
//   are neighbours in the suffix array, at an index i with LCP[i] ≥ m (the minimum length) and
//   LCP[i − 1] < LCP[i] > LCP[i + 1], their common prefix of length LCP[i] being right-maximal by
//   definition; it is a MUM when it is also left-maximal (theory, Proposition 4.3). One pass over
//   the array finds every MUM: O(n) time for n = |A| + |B| bases with SA-IS, O(n) words of memory
//   (the text, the suffix array, the LCP array and the ranks: about 16 bytes a base).
// - MUMs never nest: no MUM's interval in A (or in B) contains another's (theory, Proposition 4.4),
//   so their starts in A are distinct, and so are their starts in B.
//
// Coordinates are 0-based. A forward MUM has A[a, a + length) = B[b, b + length); a reverse one
// A[a, a + length) = the reverse complement of B[b, b + length), b on B's forward strand. Genomes
// shorter than 2^32 − 1 bases together (32-bit suffix arrays).

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/formal/suffix_array.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::comparative {

using sequence::Code;
using sequence::Codes;

struct Match {
  std::size_t a = 0, b = 0;  // starts in A and in B (B's forward strand)
  std::size_t length = 0;
  bool reverse = false;
  std::size_t a_end() const { return a + length; }
  std::size_t b_end() const { return b + length; }
  friend bool operator==(const Match&, const Match&) = default;
  friend bool operator<(const Match& x, const Match& y) {
    if (x.reverse != y.reverse) return x.reverse < y.reverse;
    if (x.a != y.a) return x.a < y.a;
    if (x.b != y.b) return x.b < y.b;
    return x.length < y.length;
  }
};

namespace detail {

// The runs of x free of N, as words, with where each begins in x.
inline void runs(const Codes& x, std::vector<Codes>& words, std::vector<std::size_t>& offsets) {
  std::size_t i = 0;
  while (i < x.size()) {
    if (x[i] > 3) {
      if (x[i] != sequence::N) throw std::invalid_argument("genomics::comparative: DNA codes 0 to 4");
      i++;
      continue;
    }
    std::size_t j = i;
    while (j < x.size() && x[j] <= 3) j++;
    words.emplace_back(x.begin() + static_cast<std::ptrdiff_t>(i), x.begin() + static_cast<std::ptrdiff_t>(j));
    offsets.push_back(i);
    i = j;
  }
}

}  // namespace detail

// The MUMs of A and B (of A and the reverse complement of B when `reverse`) of at least
// min_length bases, in increasing a.
inline std::vector<Match> maximal_unique_matches(const Codes& A, const Codes& B, std::size_t min_length, bool reverse = false) {
  if (min_length == 0) throw std::invalid_argument("genomics::comparative::maximal_unique_matches: min_length >= 1");
  const Codes rc = reverse ? sequence::reverse_complement(B) : Codes{};
  const Codes& Bs = reverse ? rc : B;
  std::vector<Codes> words;
  std::vector<std::size_t> offsets;
  detail::runs(A, words, offsets);
  const std::size_t runs_a = words.size();
  detail::runs(Bs, words, offsets);
  std::vector<Match> out;
  if (runs_a == 0 || words.size() == runs_a) return out;
  const morphe::formal::Concatenation c = morphe::formal::concatenate(words, 4);
  const std::vector<std::uint32_t> sa = morphe::formal::generalised_suffix_array<std::uint32_t>(c);
  const std::vector<std::uint32_t> lcp = morphe::formal::lcp_array(c.text, sa);
  const std::size_t n = sa.size();
  for (std::size_t i = 1; i < n; i++) {
    const std::size_t l = lcp[i];
    if (l < min_length) continue;
    if (lcp[i - 1] >= l || (i + 1 < n && lcp[i + 1] >= l)) continue;  // a third suffix shares it
    std::size_t p = sa[i - 1], q = sa[i];
    auto [wp, op] = c.locate(p);
    auto [wq, oq] = c.locate(q);
    if ((wp < runs_a) == (wq < runs_a)) continue;  // both in one genome: a repeat, not a match
    if (op > 0 && oq > 0 && c.text[p - 1] == c.text[q - 1]) continue;  // extends to the left
    if (wp >= runs_a) {
      std::swap(p, q);
      std::swap(wp, wq);
      std::swap(op, oq);
    }
    Match m;
    m.a = offsets[wp] + op;
    m.length = l;
    const std::size_t b = offsets[wq] + oq;  // on the strand of Bs
    m.b = reverse ? B.size() - b - l : b;
    m.reverse = reverse;
    out.push_back(m);
  }
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace genomics::comparative
