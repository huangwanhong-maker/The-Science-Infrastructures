// Banded alignment with affine gaps and traceback: Gotoh's recurrences restricted to the
// vertices (i, j) of the edit graph whose diagonal j − i lies in [low, high] (Chao, Pearson and
// Miller 1992, CABIOS 8, 481–487; the band of Fickett 1984, Nucleic Acids Res. 12, 175–179, and
// Ukkonen 1985, Information and Control 64, 100–118). It is gotoh.hpp's schedule on a sub-graph:
// the same kernel, the same modes and the same tie-breaking, over the paths that stay in the band.
//
// When it equals the full alignment. The banded optimum is the maximum over the paths inside
// the band, so it never exceeds gotoh's, and equals it exactly when some optimal path of the
// whole edit graph stays inside the band (then the two may still return different co-optimal
// paths). For global alignment a band can be certified after the fact (banded_is_optimal): a
// path that leaves [low, high] has at least
//   g_out = min(2 (high + 1) − (n − m), (n − m) − 2 (low − 1))
// gap columns (it must go out to the diagonal and come back to n − m), and a path with g gap
// columns scores at most U(g) = s_max (m + n − g) / 2 − o − g e. If U(g_out) ≤ the banded score,
// no path outside the band beats it, and the banded alignment is optimal. banded_global doubles
// the band until the certificate holds (the iteration of Ukkonen 1985 and Fickett 1984), so it is
// always optimal and costs O(m w) for the final half-width w.
//
// Complexity: O(m w) time and O(m w) bytes of directions for a band of w = high − low + 1
// diagonals, O(n) for the scores.

#pragma once

#include <genomics/align/gotoh.hpp>

#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace genomics::align {

// An optimal alignment among the paths whose diagonals j − i stay in [diagonal_low,
// diagonal_high]. Throws std::invalid_argument when the band holds no path the mode allows
// (for global, unless low ≤ 0 ≤ high and low ≤ n − m ≤ high).
inline Alignment banded(const Codes& query, const Codes& target, const Scoring& s, Mode m, std::int64_t diagonal_low, std::int64_t diagonal_high) {
  detail::validate(query, target, s, "genomics::align::banded");
  if (diagonal_low > diagonal_high) throw std::invalid_argument("genomics::align::banded: an empty band");
  const std::int64_t mm = static_cast<std::int64_t>(query.size()), nn = static_cast<std::int64_t>(target.size());
  // Diagonals beyond the edit graph change nothing; clipping keeps the arithmetic small.
  const std::int64_t low = std::max(diagonal_low, -mm), high = std::min(diagonal_high, nn);
  if (low > high) throw std::invalid_argument("genomics::align::banded: the band misses the edit graph");
  detail::Trace trace;
  const auto best = detail::fill<true>(query, target, s, m, m, low, high, &trace);
  if (!best.found) throw std::invalid_argument("genomics::align::banded: no path of the mode lies in the band");
  return detail::traceback(trace, best);
}

// Whether a global alignment of score `score` found in the band [low, high] is optimal over the
// whole edit graph, by the bound above (sufficient, not necessary).
inline bool banded_is_optimal(std::size_t m, std::size_t n, const Scoring& s, std::int64_t low, std::int64_t high, std::int64_t score) {
  const std::int64_t d = static_cast<std::int64_t>(n) - static_cast<std::int64_t>(m);
  const std::int64_t g_out = std::min(2 * (high + 1) - d, d - 2 * (low - 1));
  const std::int64_t total = static_cast<std::int64_t>(m + n);
  if (g_out > total) return true;  // no path leaves the band
  const std::int64_t smax = std::max(0, s.matrix.max_score());
  // 2 U(g) = smax (m + n − g) − 2 (o + g e), linear in g, so its maximum over g ≥ g_out is at
  // g_out or at m + n; compared doubled, without rounding.
  const auto twice_bound = [&](std::int64_t g) { return smax * (total - g) - 2 * (s.gap_open + g * s.gap_extend); };
  return std::max(twice_bound(g_out), twice_bound(total)) <= 2 * score;
}

// An optimal global alignment by bands of doubling half-width around the diagonals 0 .. n − m,
// starting from `half_width`, until banded_is_optimal certifies the result.
inline Alignment banded_global(const Codes& query, const Codes& target, const Scoring& s, std::int64_t half_width = 16) {
  const std::int64_t d = static_cast<std::int64_t>(target.size()) - static_cast<std::int64_t>(query.size());
  std::int64_t w = std::max<std::int64_t>(1, half_width);
  const std::int64_t widest = static_cast<std::int64_t>(query.size() + target.size());
  for (;;) {
    const std::int64_t low = std::min<std::int64_t>(0, d) - w, high = std::max<std::int64_t>(0, d) + w;
    Alignment a = banded(query, target, s, Mode::global, low, high);
    if (w >= widest || banded_is_optimal(query.size(), target.size(), s, low, high, a.score)) return a;
    w *= 2;
  }
}

}  // namespace genomics::align
