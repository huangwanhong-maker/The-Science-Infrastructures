// Edit distance by bit vectors: Myers' algorithm, with blocks for patterns of any length, a
// band that follows the cells that can still matter, and a traceback to a CIGAR.
//
// The edit distance with unit costs is the sum over the paths of the edit graph in the (min, +)
// semiring with weights 0 for a match and 1 for a mismatch or a gap. With the query (the
// pattern, length m) as rows i = 0 .. m and the target (the text, length n) as columns
// j = 0 .. n, its matrix C has the boundary C(i, 0) = i and, in the two modes kept here,
//   global      C(0, j) = j: both sequences end to end, the answer C(m, n);
//   semiglobal  C(0, j) = 0: the query end to end, the target's ends free (approximate string
//               matching, Sellers 1980), the answers C(m, j) for every end j.
// A match is two equal codes below `sigma` (4 for DNA): N, or any code ≥ sigma, matches
// nothing, not even itself, as align::edit_distance counts it.
//
// The schedule. Adjacent cells of C differ by −1, 0 or +1 (Ukkonen 1985), so a column is its
// vertical differences, two bits per row: Pv (the difference is +1) and Mv (−1). Myers (1999,
// J. ACM 46, 395–415; in the formulation of Hyyrö 2001, "Explaining and extending
// the bit-parallel approximate string matching algorithm of Myers", Tech. Rep. A-2001-10,
// Univ. Tampere) computes the next column from the previous one and the match vector Eq of
// the text's symbol with a handful of word operations; the one addition carries the effect of
// a match down a run of rows, which is how a whole column of 64 cells is one step. A pattern
// longer than 64 is cut into blocks of 64 rows, each passing to the next the horizontal
// difference at its last row (Myers 1999, the blocks; Hyyrö 2003, "A bit-vector algorithm for computing
// Levenshtein and Damerau edit distances", Nordic J. Computing 10, 29–39, for the carry of
// a horizontal −1 into a block: Eq |= 1).
//
// The band. For a threshold k, only the cells that can lie on a path of cost ≤ k are needed:
// in semiglobal mode those with C ≤ k, in global mode those with C + h ≤ k, where
// h(i, j) = |(m − i) − (n − j)| is the least cost from (i, j) to (m, n). Two facts make the
// band of blocks that holds them move by at most one row down per column, and never up at its
// top (theory note fast-alignment.md §3): C is non-decreasing along diagonals
// (C(i, j) ≥ C(i − 1, j − 1), Ukkonen 1985), and h is consistent. So a column computes
// the blocks [first, last], adds the block below when the last row of the band was within k
// in the previous column, and drops blocks whose every cell is beyond k (by the bound
// C ≥ score − (rows below)) at the bottom, and in global mode at the top. Cells outside the
// band are taken as their upper bounds (vertical differences of +1), which can only raise the
// values of cells beyond k; every cell within k is exact. This is the band of Ukkonen (1985,
// Information and Control 64, 100–118) in blocks, as Edlib computes it (Šošić and Šikić 2017,
// Bioinformatics 33, 1394–1395), with conditions proved here rather than tuned. The distance
// itself is found by doubling k from 64 until it is reached (Ukkonen 1985).
//
// The traceback keeps the band's (Pv, Mv) and the score of each block's last row for every
// column; a cell's value is that score less the differences below it (two popcounts), so the
// path is read back from (m, n) (global) or from the best end (semiglobal) in O(n + m) steps,
// preferring a pair, then a query base against a gap (I), then a target base (D).
//
// Complexity: O(⌈m/64⌉ n) word operations without a band; with the band, O(⌈k/64⌉ n) for the
// threshold k (≤ 2d for the distance d, by doubling). Memory: (σ + 1)⌈m/64⌉ words for the
// match vectors; the traceback 20 bytes per block of the band per column (5000 × 5000 at 10%
// divergence: about 2 MB; whole genomes need the wavefront algorithm's BiWFA instead).
// Distances fit in 32 bits (m + n < 2³¹).

#pragma once

#include <genomics/align/types.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::align {

// A match of the whole query ending at target position `end` (exclusive), with its distance.
struct EditHit {
  std::size_t end = 0;
  std::int64_t distance = 0;
  friend bool operator==(const EditHit&, const EditHit&) = default;
};

// The work of a run, against the (m + 1)(n + 1) cells of the full matrix: the passes made (one
// per threshold tried) and the blocks of 64 cells computed.
struct MyersStats {
  std::uint64_t passes = 0, blocks = 0;
  std::uint64_t cells() const { return 64 * blocks; }
};

// The query prepared for Myers' algorithm: for each symbol c < σ and each block b, the word
// whose bit t is set when the query's base 64b + t is c; symbol σ (every code ≥ σ) matches
// nothing.
class MyersPattern {
 public:
  using Word = std::uint64_t;
  static constexpr std::size_t W = 64;

  explicit MyersPattern(const Codes& pattern, int sigma = 4) : m_(pattern.size()), sigma_(sigma) {
    if (sigma < 1 || sigma > 255) throw std::invalid_argument("genomics::align::MyersPattern: 1 <= sigma <= 255");
    blocks_ = (m_ + W - 1) / W;
    peq_.assign(static_cast<std::size_t>(sigma + 1) * blocks_, 0);
    for (std::size_t i = 0; i < m_; i++)
      if (pattern[i] < sigma) peq_[pattern[i] * blocks_ + i / W] |= Word{1} << (i % W);
  }
  std::size_t size() const { return m_; }
  std::size_t blocks() const { return blocks_; }
  int sigma() const { return sigma_; }
  // The match vectors of a text symbol, one word per block.
  const Word* eq(Code c) const { return peq_.data() + static_cast<std::size_t>(c < sigma_ ? c : sigma_) * blocks_; }
  // The last row of block b, as a row index of C (1-based rows), and the bit that holds it.
  std::size_t bottom(std::size_t b) const { return std::min((b + 1) * W, m_); }
  Word last_bit(std::size_t b) const { return Word{1} << ((bottom(b) - 1) % W); }

 private:
  std::size_t m_ = 0, blocks_ = 0;
  int sigma_ = 4;
  std::vector<Word> peq_;  // symbol-major: peq_[c · blocks + b]
};

namespace detail {

using MyersWord = std::uint64_t;

// One block of one column (Myers 1999; Hyyrö 2003): the vertical
// differences (Pv, Mv) of the previous column become this column's, given the match vector Eq
// and the horizontal difference hin ∈ {−1, 0, +1} entering at the block's top; returns the
// horizontal difference at the block's last row (the bit `last`).
inline int myers_block(MyersWord& Pv, MyersWord& Mv, MyersWord Eq, int hin, MyersWord last) {
  const MyersWord hneg = hin < 0 ? 1 : 0, hpos = hin > 0 ? 1 : 0;
  const MyersWord Xv = Eq | Mv;
  Eq |= hneg;
  const MyersWord Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
  MyersWord Ph = Mv | ~(Xh | Pv);
  MyersWord Mh = Pv & Xh;
  const int hout = (Ph & last) ? 1 : (Mh & last) ? -1 : 0;
  Ph = (Ph << 1) | hpos;
  Mh = (Mh << 1) | hneg;
  Pv = Mh | ~(Xv | Ph);
  Mv = Ph & Xv;
  return hout;
}

// The band of every column, kept for the traceback: column j holds blocks [first[j], last[j]]
// at offset[j] (none when first > last).
struct MyersColumns {
  std::vector<MyersWord> P, M;
  std::vector<std::int32_t> S;  // C at each block's last row
  std::vector<std::size_t> offset;
  std::vector<std::int64_t> first, last;
};

inline constexpr std::int64_t MYERS_FAR = std::numeric_limits<std::int64_t>::max() / 4;

// One pass at threshold k over the target t[0, n). For each column j whose band reaches the
// last row with C(m, j) ≤ k, calls end(j, C(m, j)). Keeps the band when `keep` is given.
// Returns false when the band empties (no path of cost ≤ k remains).
template <class End>
bool myers_pass(const MyersPattern& p, const Code* t, std::size_t n, bool global, std::int64_t k, MyersColumns* keep, End&& end, MyersStats* stats = nullptr) {
  using Word = MyersWord;
  if (stats) stats->passes++;
  const std::size_t m = p.size(), B = p.blocks();
  const std::int64_t mm = static_cast<std::int64_t>(m), nn = static_cast<std::int64_t>(n);
  std::vector<Word> P(B, ~Word{0}), M(B, 0);
  std::vector<std::int64_t> S(B);
  for (std::size_t b = 0; b < B; b++) S[b] = static_cast<std::int64_t>(p.bottom(b));
  // The least cost from (i, j) to the end; zero in semiglobal mode.
  const auto h = [&](std::int64_t i, std::int64_t j) -> std::int64_t {
    if (!global) return 0;
    const std::int64_t d = (mm - i) - (nn - j);
    return d < 0 ? -d : d;
  };
  // A lower bound of C + h over the cells of block b in column j, from C ≥ S − (rows below).
  const auto bound = [&](std::size_t b, std::int64_t j) -> std::int64_t {
    const std::int64_t lo = static_cast<std::int64_t>(b * MyersPattern::W) + 1, hi = static_cast<std::int64_t>(p.bottom(b));
    const std::int64_t base = S[b] - hi;  // C(i) ≥ base + i
    if (!global) return base + lo;
    // min over i in [lo, hi] of i + |c0 − i|, with c0 = m − n + j.
    const std::int64_t c0 = mm - nn + j;
    const std::int64_t g = lo <= c0 ? c0 : 2 * lo - c0;
    return base + g;
  };
  // The exact least C + h over a block, cell by cell (used every few columns).
  const auto exact = [&](std::size_t b, std::int64_t j) -> std::int64_t {
    const std::int64_t lo = static_cast<std::int64_t>(b * MyersPattern::W) + 1, hi = static_cast<std::int64_t>(p.bottom(b));
    std::int64_t c = S[b], best = c + h(hi, j);
    for (std::int64_t i = hi; i > lo; i--) {
      const std::size_t bit = static_cast<std::size_t>(i - lo);
      c -= static_cast<std::int64_t>((P[b] >> bit) & 1) - static_cast<std::int64_t>((M[b] >> bit) & 1);
      best = std::min(best, c + h(i - 1, j));
    }
    return best;
  };
  std::int64_t first = 0, last = std::min<std::int64_t>(static_cast<std::int64_t>(B) - 1, k / static_cast<std::int64_t>(MyersPattern::W));
  const auto store = [&](std::size_t j) {
    if (!keep) return;
    keep->offset[j] = keep->P.size();
    keep->first[j] = first;
    keep->last[j] = last;
    for (std::int64_t b = first; b <= last; b++) {
      keep->P.push_back(P[static_cast<std::size_t>(b)]);
      keep->M.push_back(M[static_cast<std::size_t>(b)]);
      keep->S.push_back(static_cast<std::int32_t>(S[static_cast<std::size_t>(b)]));
    }
  };
  if (keep) {
    keep->P.clear();
    keep->M.clear();
    keep->S.clear();
    keep->offset.assign(n + 1, 0);
    keep->first.assign(n + 1, 1);
    keep->last.assign(n + 1, 0);
  }
  store(0);
  if (last == static_cast<std::int64_t>(B) - 1 && mm <= k) end(std::size_t{0}, mm);
  // Row 0 belongs to the band too: in semiglobal mode always (C(0, j) = 0), in global mode while
  // j + |m − n + j| ≤ k (that sum never decreases). With no block left, the band is row 0 alone
  // (first = 0, last = −1), and a block grows below it again when row 0 allows.
  const auto row0_dead = [&](std::int64_t j) {
    if (!global) return false;
    const std::int64_t d = mm - nn + j;
    return j + (d < 0 ? -d : d) > k;
  };
  for (std::size_t j = 1; j <= n; j++) {
    const std::int64_t jj = static_cast<std::int64_t>(j);
    const Word* eq = p.eq(t[j - 1]);
    // C(bottom of the band, j − 1), and its row.
    const bool rows = last >= first;
    const std::int64_t prev_bottom = rows ? S[static_cast<std::size_t>(last)] : (global ? jj - 1 : 0);
    const std::int64_t prev_row = rows ? static_cast<std::int64_t>(p.bottom(static_cast<std::size_t>(last))) : 0;
    // Above the band: row 0 (C(0, j) − C(0, j − 1) = 1 in global mode, 0 in semiglobal), or a
    // region beyond k taken at its upper bound, +1.
    int hin = (global || first > 0) ? 1 : 0;
    for (std::int64_t b = first; b <= last; b++) {
      const std::size_t u = static_cast<std::size_t>(b);
      hin = myers_block(P[u], M[u], eq[u], hin, p.last_bit(u));
      S[u] += hin;
    }
    // The block below, when a cell of its first row can be within k (diagonal monotonicity:
    // C(r + 1, j) ≥ C(r, j − 1), with the same h).
    if (last + 1 < static_cast<std::int64_t>(B) && prev_bottom + h(prev_row, jj - 1) <= k) {
      last++;
      const std::size_t u = static_cast<std::size_t>(last);
      P[u] = ~Word{0};
      M[u] = 0;
      S[u] = prev_bottom + static_cast<std::int64_t>(p.bottom(u) - u * MyersPattern::W);
      S[u] += myers_block(P[u], M[u], eq[u], hin, p.last_bit(u));
    }
    if (stats && last >= first) stats->blocks += static_cast<std::uint64_t>(last - first + 1);
    const bool strong = (j % 8) == 0;
    while (last >= first && (bound(static_cast<std::size_t>(last), jj) > k || (strong && exact(static_cast<std::size_t>(last), jj) > k))) last--;
    const bool dead0 = row0_dead(jj);
    if (global)
      while (first <= last && (first > 0 || dead0) && (bound(static_cast<std::size_t>(first), jj) > k || (strong && exact(static_cast<std::size_t>(first), jj) > k)))
        first++;
    if (last < first && (first > 0 || dead0)) {
      if (keep) {
        keep->first.resize(j);
        keep->last.resize(j);
        keep->offset.resize(j);
      }
      return false;
    }
    store(j);
    if (last == static_cast<std::int64_t>(B) - 1 && S[B - 1] <= k) end(j, S[B - 1]);
  }
  return true;
}

// C(i, j) from the kept band; MYERS_FAR outside it.
inline std::int64_t myers_cell(const MyersPattern& p, const MyersColumns& c, bool global, std::size_t i, std::size_t j) {
  if (i == 0) return global ? static_cast<std::int64_t>(j) : 0;
  if (j >= c.first.size()) return MYERS_FAR;
  const std::int64_t b = static_cast<std::int64_t>((i - 1) / MyersPattern::W);
  if (b < c.first[j] || b > c.last[j]) return MYERS_FAR;
  const std::size_t at = c.offset[j] + static_cast<std::size_t>(b - c.first[j]);
  const std::size_t base = static_cast<std::size_t>(b) * MyersPattern::W, top = p.bottom(static_cast<std::size_t>(b)) - base;
  const std::size_t from = i - base;  // the bit of row i + 1
  MyersWord mask = top == 64 ? ~MyersWord{0} : ((MyersWord{1} << top) - 1);
  mask &= from >= 64 ? 0 : (~MyersWord{0} << from);
  return c.S[at] - std::popcount(c.P[at] & mask) + std::popcount(c.M[at] & mask);
}

// The path into (i, j) read back through the kept band; returns the CIGAR and sets the start.
inline Cigar myers_traceback(const MyersPattern& p, const MyersColumns& c, const Codes& query, const Codes& target, bool global, std::size_t i,
                             std::size_t j, std::size_t& target_begin) {
  Cigar rev;
  const int sigma = p.sigma();
  while (i > 0 || (global && j > 0)) {
    if (i == 0) { push(rev, 'D', static_cast<std::uint32_t>(j)); j = 0; break; }
    if (j == 0) { push(rev, 'I', static_cast<std::uint32_t>(i)); i = 0; break; }
    const std::int64_t v = myers_cell(p, c, global, i, j);
    const bool same = query[i - 1] == target[j - 1] && query[i - 1] < sigma;
    if (myers_cell(p, c, global, i - 1, j - 1) + (same ? 0 : 1) == v) { push(rev, 'M'); i--; j--; continue; }
    if (myers_cell(p, c, global, i - 1, j) + 1 == v) { push(rev, 'I'); i--; continue; }
    if (myers_cell(p, c, global, i, j - 1) + 1 == v) { push(rev, 'D'); j--; continue; }
    throw std::logic_error("genomics::align::myers: the traceback lost its path");
  }
  target_begin = j;
  return reversed(std::move(rev));
}

inline void myers_mode(Mode mode, const char* who) {
  if (mode != Mode::global && mode != Mode::semiglobal) throw std::invalid_argument(std::string(who) + ": Myers' algorithm is kept for the global and semiglobal modes");
}

}  // namespace detail

// Every end position in the target of a match of the whole query with at most k edits
// (semiglobal), in increasing order of the end, with its distance. O(⌈k/64⌉ n).
inline std::vector<EditHit> myers_search(const MyersPattern& p, const Codes& target, std::int64_t k, MyersStats* stats = nullptr) {
  std::vector<EditHit> hits;
  if (k < 0) return hits;
  if (p.size() == 0) {
    for (std::size_t j = 0; j <= target.size(); j++) hits.push_back({j, 0});
    return hits;
  }
  detail::myers_pass(p, target.data(), target.size(), false, k, nullptr, [&](std::size_t j, std::int64_t d) { hits.push_back({j, d}); }, stats);
  return hits;
}
inline std::vector<EditHit> myers_search(const Codes& query, const Codes& target, std::int64_t k, int sigma = 4, MyersStats* stats = nullptr) {
  return myers_search(MyersPattern(query, sigma), target, k, stats);
}

// The edit distance of the query and the target in the global or semiglobal mode (for
// semiglobal, the least over the target's ends), or −1 when it exceeds max_k (max_k < 0: no
// limit). By doubling the threshold from 64.
inline std::int64_t myers_distance(const MyersPattern& p, const Codes& target, Mode mode = Mode::global, std::int64_t max_k = -1, MyersStats* stats = nullptr) {
  detail::myers_mode(mode, "genomics::align::myers_distance");
  const bool global = mode == Mode::global;
  const std::int64_t m = static_cast<std::int64_t>(p.size()), n = static_cast<std::int64_t>(target.size());
  const std::int64_t ceiling = global ? std::max(m, n) : m;
  if (m == 0 || n == 0) {
    const std::int64_t d = global ? std::max(m, n) : m;
    return max_k >= 0 && d > max_k ? -1 : d;
  }
  std::int64_t k = std::min<std::int64_t>(64, ceiling);
  if (max_k >= 0) k = std::min(k, max_k);
  for (;;) {
    std::int64_t best = -1;
    if (!global || std::abs(m - n) <= k)
      detail::myers_pass(
          p, target.data(), target.size(), global, k, nullptr,
          [&](std::size_t j, std::int64_t d) {
            if (!global || static_cast<std::int64_t>(j) == n) best = best < 0 ? d : std::min(best, d);
          },
          stats);
    if (best >= 0) return best;
    if (k >= ceiling || (max_k >= 0 && k >= max_k)) return -1;
    k = std::min(2 * k, ceiling);
    if (max_k >= 0) k = std::min(k, max_k);
  }
}
inline std::int64_t myers_distance(const Codes& query, const Codes& target, Mode mode = Mode::global, std::int64_t max_k = -1, int sigma = 4,
                                   MyersStats* stats = nullptr) {
  return myers_distance(MyersPattern(query, sigma), target, mode, max_k, stats);
}

// An optimal alignment under unit costs, with its CIGAR (M, I, D); the score is minus the
// distance. Semiglobal: the leftmost end of least distance. Empty when the distance exceeds
// max_k (max_k < 0: no limit).
inline std::optional<Alignment> myers_align(const Codes& query, const Codes& target, Mode mode = Mode::global, std::int64_t max_k = -1, int sigma = 4,
                                            MyersStats* stats = nullptr) {
  detail::myers_mode(mode, "genomics::align::myers_align");
  const bool global = mode == Mode::global;
  const MyersPattern p(query, sigma);
  const std::int64_t m = static_cast<std::int64_t>(query.size()), n = static_cast<std::int64_t>(target.size());
  Alignment a;
  a.query_begin = 0;
  a.query_end = query.size();
  if (m == 0 || n == 0) {
    const std::int64_t d = global ? std::max(m, n) : m;
    if (max_k >= 0 && d > max_k) return std::nullopt;
    if (global) { push(a.cigar, 'D', static_cast<std::uint32_t>(n)); push(a.cigar, 'I', static_cast<std::uint32_t>(m)); a.target_end = target.size(); }
    else push(a.cigar, 'I', static_cast<std::uint32_t>(m));
    a.score = -d;
    return a;
  }
  const std::int64_t d = myers_distance(p, target, mode, max_k, stats);
  if (d < 0) return std::nullopt;
  detail::MyersColumns c;
  std::size_t end_j = target.size();
  bool found = false;
  detail::myers_pass(
      p, target.data(), target.size(), global, std::max<std::int64_t>(d, 0), &c,
      [&](std::size_t j, std::int64_t v) {
        if (global || found || v != d) return;
        end_j = j;
        found = true;
      },
      stats);
  std::size_t tb = 0;
  a.cigar = detail::myers_traceback(p, c, query, target, global, query.size(), end_j, tb);
  a.target_begin = tb;
  a.target_end = end_j;
  a.score = -d;
  return a;
}

}  // namespace genomics::align
