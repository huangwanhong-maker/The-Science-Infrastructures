// One query against many targets, one target per lane: inter-sequence vectorisation of the
// Smith–Waterman score with affine gaps, for scoring a batch (a database search, the candidates
// of a mapper).
//
// The schedule (Rognes 2011, BMC Bioinformatics 12, 221, SWIPE; after Wozniak 1997, CABIOS
// 13, 145–150, and Rognes and Seeberg 2000). Where Farrar's striped kernel (striped.hpp) puts the
// cells of one column in the lanes, here lane l holds the same cell (i, j) of W different
// alignments: the query against target l. Every dependency of Gotoh's recurrence then stays in
// its lane, so a column is computed row by row with no correction loop:
//   H(i, j) = max(0, H(i − 1, j − 1) + s(q_i, t_l[j]), E(i, j), F(i, j))
//   E(i, j + 1) = max(H(i, j) − (o + e), E(i, j) − e),   F(i + 1, j) = max(H(i, j) − (o + e), F(i, j) − e)
// The price is the score vector: lane l needs s(q_i, t_l[j]), a different target base in each
// lane. For each column the kernel builds, once, the vectors S_a = (s(a, t_0[j]), …, s(a,
// t_{W−1}[j])) for every query symbol a (Rognes 2011's "score profile"), σW scalar writes against
// the m vector rows that use them. Targets are taken in order of length, so a batch wastes
// little on the shorter ones; past its end a target scores −16384 a base, which only lowers H.
//
// Lanes of 16 bits (calx::cpu::Lanes<int16_t>, 32 on AVX-512), saturating; a lane whose best
// reaches 32767 − max s may have saturated and is rescored alone by the striped kernel, whose
// widths go to 32 bits. Scores only: the ends and the paths are striped_local's.
//
// Complexity: O(m · Σ_batches max length) vector operations, plus σW per column; memory O(m)
// vectors. Gap extension at least 1, codes below the matrix's size.

#pragma once

#include <genomics/align/striped.hpp>
#include <genomics/align/types.hpp>

#include <calx/cpu/lanes.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpsabi"

namespace genomics::align {

// The work of a batch run: the batches, the cells computed (lanes × rows × columns, padding
// included), the targets rescored for saturation.
struct InterStats {
  std::uint64_t batches = 0, cells = 0, rescored = 0;
};

// The best local score of the query against each target, in the order of the targets.
inline std::vector<std::int64_t> inter_scores(const Codes& query, const std::vector<Codes>& targets, const Scoring& s, InterStats* stats = nullptr) {
  using V = calx::cpu::Lanes<std::int16_t>;
  constexpr std::size_t W = V::width;
  const Matrix& M = s.matrix;
  const int sigma = M.size;
  if (sigma <= 0 || M.score.size() != static_cast<std::size_t>(sigma) * static_cast<std::size_t>(sigma))
    throw std::invalid_argument("genomics::align::inter_scores: the substitution matrix is not square");
  if (s.gap_open < 0 || s.gap_extend < 1 || s.gap_open + s.gap_extend > 16000) throw std::invalid_argument("genomics::align::inter_scores: gaps need open >= 0, extend >= 1");
  if (M.max_score() > 1000 || M.min_score() < -1000) throw std::invalid_argument("genomics::align::inter_scores: scores beyond ±1000");
  for (Code c : query)
    if (c >= sigma) throw std::invalid_argument("genomics::align::inter_scores: a query code outside the matrix");
  for (const auto& t : targets)
    for (Code c : t)
      if (c >= sigma) throw std::invalid_argument("genomics::align::inter_scores: a target code outside the matrix");
  const std::size_t m = query.size();
  std::vector<std::int64_t> out(targets.size(), 0);
  if (m == 0 || targets.empty()) return out;
  const std::int16_t PAD = -16384;
  const std::int64_t limit = 32767 - std::max(0, M.max_score());
  // The targets by length, so that a batch is of like lengths.
  std::vector<std::size_t> order(targets.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return targets[a].size() < targets[b].size(); });
  const V zero = V::zero(), voe = V::splat(static_cast<std::int16_t>(s.gap_open + s.gap_extend)), ve = V::splat(static_cast<std::int16_t>(s.gap_extend));
  std::vector<V> Hcol(m + 1), E(m + 1), S(static_cast<std::size_t>(sigma));
  std::vector<std::int16_t> lanes(static_cast<std::size_t>(sigma) * W);
  std::vector<std::size_t> rescore;
  for (std::size_t first = 0; first < order.size(); first += W) {
    const std::size_t count = std::min(W, order.size() - first);
    std::size_t longest = 0;
    for (std::size_t l = 0; l < count; l++) longest = std::max(longest, targets[order[first + l]].size());
    std::fill(Hcol.begin(), Hcol.end(), zero);
    std::fill(E.begin(), E.end(), zero);
    V best = zero;
    for (std::size_t j = 0; j < longest; j++) {
      // The score vectors of this column, one per query symbol.
      for (std::size_t l = 0; l < W; l++) {
        const Codes* t = l < count ? &targets[order[first + l]] : nullptr;
        const bool in = t && j < t->size();
        const Code c = in ? (*t)[j] : 0;
        for (int a = 0; a < sigma; a++) lanes[static_cast<std::size_t>(a) * W + l] = in ? static_cast<std::int16_t>(M(static_cast<Code>(a), c)) : PAD;
      }
      for (int a = 0; a < sigma; a++) S[static_cast<std::size_t>(a)] = V::load(lanes.data() + static_cast<std::size_t>(a) * W);
      V diag = Hcol[0], F = zero;
      for (std::size_t i = 1; i <= m; i++) {
        const V up = Hcol[i];
        V h = max(adds(diag, S[query[i - 1]]), zero);
        h = max(h, E[i]);
        h = max(h, F);
        best = max(best, h);
        Hcol[i] = h;
        const V hg = max(subs(h, voe), zero);
        E[i] = max(max(subs(E[i], ve), zero), hg);
        F = max(max(subs(F, ve), zero), hg);
        diag = up;
      }
    }
    for (std::size_t l = 0; l < count; l++) {
      const std::int64_t v = best[l];
      if (v >= limit) rescore.push_back(order[first + l]);
      else out[order[first + l]] = v;
    }
    if (stats) {
      stats->batches++;
      stats->cells += static_cast<std::uint64_t>(longest) * m * W;
    }
  }
  if (!rescore.empty()) {
    const StripedProfile p(query, s);
    for (std::size_t k : rescore) out[k] = p.end(targets[k]).score;
    if (stats) stats->rescored += rescore.size();
  }
  return out;
}

}  // namespace genomics::align

#pragma GCC diagnostic pop
