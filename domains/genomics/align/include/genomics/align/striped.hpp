// Local alignment on vector lanes: Farrar's striped Smith–Waterman with affine gaps, its score
// in 8-bit lanes first and in 16- or 32-bit lanes when 8 bits overflow, the end of the best
// alignment, its start from the reversed prefixes, and its path.
//
// The sum. Smith and Waterman's local alignment is the (max, +) sum over the paths of the edit
// graph that may begin and end anywhere, with Gotoh's three states for affine gaps (gotoh.hpp):
//   E(i, j) = max(H(i, j − 1) − (o + e), E(i, j − 1) − e)       a target base against a gap
//   F(i, j) = max(H(i − 1, j) − (o + e), F(i − 1, j) − e)       a query base against a gap
//   H(i, j) = max(0, H(i − 1, j − 1) + s(q_i, t_j), E(i, j), F(i, j))
// with the query as rows. Its optimum is gotoh(…, Mode::local)'s.
//
// The schedule (Farrar 2007, Bioinformatics 23, 156–161). A column of the
// matrix (one target base) is computed a vector at a time, the query cut into W stripes: lane l
// of segment s holds the query position s + l·S, S = ⌈m/W⌉, so that the diagonal dependency
// H(i − 1, j − 1) is the previous column's segment s − 1 in the same lane (and for s = 0, its
// last segment shifted up one lane). The vertical dependency F, which runs down the column,
// crosses from lane to lane only at the wrap; it is first taken as 0 and then corrected by the
// "lazy F" loop, which carries F from each lane's last segment into the next lane's first and
// stops as soon as no lane can raise an H (F ≤ H − (o + e) in every lane): rarely more than one
// pass, because a gap that long seldom pays. The corrected H also opens gaps E into the next
// column (Farrar's update of E in the loop), so every transition of Gotoh's recurrence is kept
// and the scores are Gotoh's exactly.
// - The query profile (Rognes and Seeberg 2000, Bioinformatics 16, 699–706): for each target
//   symbol, the segments of scores s(q_i, c) in the striped order, so the inner loop loads a
//   vector and never gathers.
// - Lanes of 8 bits first, as SSW (Zhao, Lee, Garrison and Marth 2013, PLoS ONE 8, e82138):
//   unsigned bytes with the profile biased by β = max(1, −min s) so that it is never negative
//   (H + s + β saturating, then − β saturating at 0, which is the local floor itself); when the
//   best score reaches 255 − β a byte may have saturated, and the pass is redone in 16-bit lanes,
//   and in 32-bit lanes past 32767 − max s. Saturating arithmetic is calx::cpu::Lanes' adds and
//   subs; the width is the build's widest (64 bytes of AVX-512: 64 lanes of 8 bits).
// - The end: the first cell of the best score in row-major order (the smallest query end, then
//   the smallest target end), gotoh.hpp's rule, found from the columns whose maximum reaches the
//   best. The start: the same kernel on the reversed prefixes q[0, i) and t[0, j), whose best
//   alignment ends, reversed, exactly at (i, j) (any better-placed end would be an earlier cell
//   of the best score in the forward matrix). The path: the global alignment of the region,
//   whose optimum is the local one's, by banded_global (banded.hpp).
//
// Complexity: O(⌈m/W⌉ n) vector operations per pass, plus the lazy-F passes; the profile
// σ ⌈m/W⌉ W bytes per width. Gap extension must be at least 1 (the lazy-F loop's progress), and
// codes must be below the matrix's size. A score of 0 is the empty alignment at (0, 0).

#pragma once

#include <genomics/align/banded.hpp>
#include <genomics/align/types.hpp>

#include <calx/cpu/lanes.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpsabi"

namespace genomics::align {

// The work of a run: the width of the last pass (8, 16 or 32 bits), the passes made, the cells
// they computed (columns × striped query length), and the lazy-F corrections.
struct StripedStats {
  int bits = 0;
  std::uint64_t passes = 0, cells = 0, lazy_f = 0;
};

namespace detail {

template <class T>
using StripedV = calx::cpu::Lanes<T>;

// The striped query profile of one lane type: for each symbol c, S vectors whose lane l of
// segment s holds s(q[s + l·S], c) (+ β in 8-bit lanes), and −β past the query's end.
template <class T>
struct StripedTable {
  std::size_t segments = 0;
  std::vector<StripedV<T>> profile;  // [c · S + s]

  void build(const Codes& query, const Matrix& matrix, int bias, bool biased) {
    using V = StripedV<T>;
    constexpr std::size_t W = V::width;
    const std::size_t m = query.size();
    segments = std::max<std::size_t>(1, (m + W - 1) / W);
    profile.assign(static_cast<std::size_t>(matrix.size) * segments, V::zero());
    for (int c = 0; c < matrix.size; c++)
      for (std::size_t s = 0; s < segments; s++) {
        V v = V::zero();
        for (std::size_t l = 0; l < W; l++) {
          const std::size_t i = s + l * segments;
          const int score = i < m ? matrix(query[i], static_cast<Code>(c)) : -bias;
          v.set(l, static_cast<T>(biased ? score + bias : score));
        }
        profile[static_cast<std::size_t>(c) * segments + s] = v;
      }
  }
};

struct StripedBest {
  std::int64_t score = 0;
  std::size_t i = 0, j = 0;  // the cell (query, target), 0-based
  bool overflow = false;
};

// One pass of Farrar's kernel over the target. Biased: unsigned lanes with the profile biased
// by `bias`; otherwise signed lanes floored at 0. `limit`: a best score at which the lanes may
// have saturated.
template <class T, bool Biased>
StripedBest striped_pass(const StripedTable<T>& table, std::size_t m, const Code* target, std::size_t n, int bias, int gap_open, int gap_extend,
                         std::int64_t limit, StripedStats* stats) {
  using V = StripedV<T>;
  constexpr std::size_t W = V::width;
  const std::size_t S = table.segments;
  std::vector<V> H1(S, V::zero()), H2(S, V::zero()), E(S, V::zero());
  V* Hload = H1.data();
  V* Hstore = H2.data();
  const V vbias = V::splat(static_cast<T>(Biased ? bias : 0)), vgo = V::splat(static_cast<T>(gap_open + gap_extend)), vge = V::splat(static_cast<T>(gap_extend));
  const V zero = V::zero();
  // H + s, floored at 0: by the bias in unsigned lanes, by a maximum in signed ones.
  const auto add = [&](V h, V p) {
    if constexpr (Biased) return subs(adds(h, p), vbias);
    else return max(adds(h, p), zero);
  };
  const auto dec = [&](V x, V g) {
    if constexpr (Biased) return subs(x, g);
    else return max(subs(x, g), zero);
  };
  // The first query position of the column holding `value`, or m when none does.
  const auto first_row = [&](T value) {
    std::size_t best = m;
    const V target_value = V::splat(value);
    for (std::size_t s = 0; s < S; s++) {
      std::uint64_t bits = (Hstore[s] == target_value).bits();
      while (bits) {
        const std::size_t l = static_cast<std::size_t>(std::countr_zero(bits));
        const std::size_t i = s + l * S;
        if (i < m) { best = std::min(best, i); break; }
        bits &= bits - 1;
      }
    }
    return best;
  };
  StripedBest r;
  T best = 0;
  std::uint64_t lazy = 0;
  for (std::size_t j = 0; j < n; j++) {
    const V* P = table.profile.data() + static_cast<std::size_t>(target[j]) * S;
    V vF = zero, vMax = zero;
    V vH = Hstore[S - 1].shift_up(0);  // H(i − 1, j − 1) for segment 0: the last segment, one lane up
    std::swap(Hload, Hstore);
    for (std::size_t s = 0; s < S; s++) {
      vH = add(vH, P[s]);
      const V e = E[s];
      vH = max(vH, e);
      vH = max(vH, vF);
      vMax = max(vMax, vH);
      Hstore[s] = vH;
      const V hg = dec(vH, vgo);
      E[s] = max(dec(e, vge), hg);
      vF = max(dec(vF, vge), hg);
      vH = Hload[s];
    }
    // Lazy F: carry F across the lanes until it can raise no H.
    vF = vF.shift_up(0);
    std::size_t s = 0;
    while ((vF > dec(Hstore[s], vgo)).any()) {
      const V h = max(Hstore[s], vF);
      Hstore[s] = h;
      vMax = max(vMax, h);
      E[s] = max(E[s], dec(h, vgo));
      vF = dec(vF, vge);
      lazy++;
      if (++s == S) {
        s = 0;
        vF = vF.shift_up(0);
      }
    }
    const V vb = V::splat(best);
    if ((vMax > vb).any()) {
      best = vMax.hmax();
      r.i = first_row(best);
      r.j = j;
      if (static_cast<std::int64_t>(best) >= limit) {
        r.overflow = true;
        break;
      }
    } else if (best > 0 && (vMax == vb).any()) {
      const std::size_t i = first_row(best);
      if (i < r.i) { r.i = i; r.j = j; }
    }
  }
  r.score = best;
  if (stats) {
    stats->passes++;
    stats->cells += static_cast<std::uint64_t>(n) * S * W;
    stats->lazy_f += lazy;
    stats->bits = static_cast<int>(8 * sizeof(T));
  }
  return r;
}

}  // namespace detail

// The query prepared for the striped kernel: its profiles in 8-bit and 16-bit lanes (the
// 32-bit one is made when a pass needs it).
class StripedProfile {
 public:
  StripedProfile(const Codes& query, const Scoring& s) : query_(query), scoring_(s) {
    const Matrix& M = s.matrix;
    if (M.size <= 0 || M.score.size() != static_cast<std::size_t>(M.size) * static_cast<std::size_t>(M.size))
      throw std::invalid_argument("genomics::align::StripedProfile: the substitution matrix is not square");
    if (s.gap_open < 0 || s.gap_extend < 1) throw std::invalid_argument("genomics::align::StripedProfile: gaps need open >= 0 and extend >= 1");
    for (Code c : query)
      if (c >= M.size) throw std::invalid_argument("genomics::align::StripedProfile: a query code outside the matrix");
    bias_ = std::max(1, -M.min_score());
    max_ = std::max(0, M.max_score());
    // 8 bits when the biased scores and the gap costs fit in a byte.
    byte_ = max_ + bias_ <= 255 && s.gap_open + s.gap_extend <= 255;
    if (byte_) p8_.build(query, M, bias_, true);
    word_ = max_ + bias_ <= 32767 && s.gap_open + s.gap_extend <= 32767;
    if (word_) p16_.build(query, M, bias_, false);
  }
  const Codes& query() const { return query_; }
  const Scoring& scoring() const { return scoring_; }

  // The best local score and its end, by the narrowest lanes that do not overflow.
  Endpoint end(const Codes& target, StripedStats* stats = nullptr) const {
    for (Code c : target)
      if (c >= scoring_.matrix.size) throw std::invalid_argument("genomics::align::striped: a target code outside the matrix");
    const std::size_t m = query_.size(), n = target.size();
    if (m == 0 || n == 0) return {};
    const int o = scoring_.gap_open, e = scoring_.gap_extend;
    detail::StripedBest b;
    b.overflow = true;
    if (byte_) b = detail::striped_pass<std::uint8_t, true>(p8_, m, target.data(), n, bias_, o, e, 255 - bias_, stats);
    if (b.overflow && word_) b = detail::striped_pass<std::int16_t, false>(p16_, m, target.data(), n, bias_, o, e, 32767 - max_, stats);
    if (b.overflow) {
      if (p32_.profile.empty()) p32_.build(query_, scoring_.matrix, bias_, false);
      b = detail::striped_pass<std::int32_t, false>(p32_, m, target.data(), n, bias_, o, e, std::numeric_limits<std::int32_t>::max() - max_, stats);
      if (b.overflow) throw std::overflow_error("genomics::align::striped: a score beyond 32 bits");
    }
    if (b.score <= 0) return {};
    return {b.score, b.i + 1, b.j + 1};
  }

 private:
  Codes query_;
  Scoring scoring_;
  int bias_ = 1, max_ = 0;
  bool byte_ = false, word_ = false;
  detail::StripedTable<std::uint8_t> p8_;
  detail::StripedTable<std::int16_t> p16_;
  mutable detail::StripedTable<std::int32_t> p32_;
};

// The best local score and the end of the first cell that holds it (row-major), as gotoh_end
// in Mode::local gives them.
inline Endpoint striped_end(const StripedProfile& profile, const Codes& target, StripedStats* stats = nullptr) { return profile.end(target, stats); }

// The best local alignment with its path: the end by the striped kernel, the start by the
// kernel on the reversed prefixes, the path by the banded global alignment of the region.
inline Alignment striped_local(const StripedProfile& profile, const Codes& target, StripedStats* stats = nullptr) {
  const Endpoint end = profile.end(target, stats);
  Alignment a;
  if (end.score <= 0) return a;
  const Codes& query = profile.query();
  const Codes rq(query.rbegin() + static_cast<std::ptrdiff_t>(query.size() - end.query_end), query.rend());
  const Codes rt(target.rbegin() + static_cast<std::ptrdiff_t>(target.size() - end.target_end), target.rend());
  const Endpoint back = StripedProfile(rq, profile.scoring()).end(rt, stats);
  if (back.score != end.score) throw std::logic_error("genomics::align::striped_local: the reversed prefixes lost the score");
  a.query_begin = end.query_end - back.query_end;
  a.target_begin = end.target_end - back.target_end;
  a.query_end = end.query_end;
  a.target_end = end.target_end;
  const Codes qr(query.begin() + static_cast<std::ptrdiff_t>(a.query_begin), query.begin() + static_cast<std::ptrdiff_t>(a.query_end));
  const Codes tr(target.begin() + static_cast<std::ptrdiff_t>(a.target_begin), target.begin() + static_cast<std::ptrdiff_t>(a.target_end));
  const Alignment g = banded_global(qr, tr, profile.scoring());
  if (g.score != end.score) throw std::logic_error("genomics::align::striped_local: the region's global optimum is not the local one");
  a.cigar = g.cigar;
  a.score = end.score;
  return a;
}
inline Alignment striped_local(const Codes& query, const Codes& target, const Scoring& s, StripedStats* stats = nullptr) {
  return striped_local(StripedProfile(query, s), target, stats);
}

}  // namespace genomics::align

#pragma GCC diagnostic pop
