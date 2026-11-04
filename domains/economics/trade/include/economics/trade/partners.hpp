// The partners of trade: how each economy's exports spread over its partners, how intensely it
// trades with each, and how much of a pair's trade goes both ways, on nomos's measures of
// concentration and specialisation (nomos::measurement).
//
// x_ij ≥ 0 are the exports of i to j, a square matrix with x_ii = 0; X_i = Σ_j x_ij the exports of
// i, M_j = Σ_i x_ij the imports of j, W = Σ x the world's trade.
//
// - Partner concentration. With s_ij = x_ij / X_i, the Herfindahl–Hirschman index H_i = Σ_j s_ij²
//   (Hirschman 1945, on the concentration of a country's trade by partner; Herfindahl 1950), in
//   [1/(n−1), 1], and its numbers equivalent 1/H_i, the number of equal partners that would give
//   the same index (Adelman 1969).
// - Trade intensity. I_ij = (x_ij / X_i) / (M_j / W): j's share of i's exports over j's share of
//   the world's imports (Brown 1949; Kojima 1964). I > 1: i sends more to j than j's size as a
//   buyer predicts, from nearness, history, language or policy. It is the location quotient of the
//   bilateral matrix with partners as categories, so its mean over exporters, weighted by their
//   exports, is 1 for each partner. (Kojima's variant divides by M_j / (W − M_i), leaving out i's
//   own imports; not here.)
// - Two-way trade. Grubel and Lloyd's GL_ij = 1 − |x_ij − x_ji| / (x_ij + x_ji) (Grubel and Lloyd
//   1975) for the pair: the share of the pair's gross trade that is reciprocal exchange in the
//   circulation module's split (economics::net_transfer), since x_ij + x_ji − |x_ij − x_ji| =
//   2 min(x_ij, x_ji). The world's reciprocal share is the aggregate GL, weighted by gross trade.
//
// Undefined values (an economy without exports, a pair without trade) are NaN.

#pragma once

#include <calx/linalg/rectangular.hpp>
#include <nomos/measurement/concentration.hpp>
#include <nomos/measurement/specialisation.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace economics::trade {

using calx::linalg::Rectangular;

namespace detail {
inline void require_square(const Rectangular& x) {
  if (x.rows == 0 || x.rows != x.cols) throw std::invalid_argument("economics::trade: a square matrix of exports, exporters by importers");
}
}  // namespace detail

inline std::vector<double> exports(const Rectangular& x) {
  std::vector<double> X(x.rows, 0.0);
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t j = 0; j < x.cols; j++) X[i] += x(i, j);
  return X;
}

inline std::vector<double> imports(const Rectangular& x) {
  std::vector<double> M(x.cols, 0.0);
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t j = 0; j < x.cols; j++) M[j] += x(i, j);
  return M;
}

struct PartnerConcentration {
  std::vector<double> hhi;          // H_i = Σ_j s_ij²
  std::vector<double> equivalent;   // 1 / H_i
  std::vector<std::size_t> first;   // i's largest partner
  std::vector<double> first_share;  // its share of i's exports
};

inline PartnerConcentration partner_concentration(const Rectangular& x) {
  detail::require_square(x);
  constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
  const std::size_t n = x.rows;
  PartnerConcentration P;
  P.hhi.assign(n, NaN);
  P.equivalent.assign(n, NaN);
  P.first.assign(n, n);
  P.first_share.assign(n, NaN);
  for (std::size_t i = 0; i < n; i++) {
    std::vector<double> row(x.a.begin() + static_cast<std::ptrdiff_t>(i * n), x.a.begin() + static_cast<std::ptrdiff_t>((i + 1) * n));
    double total = 0;
    for (double v : row) total += v;
    if (!(total > 0)) continue;
    P.hhi[i] = nomos::measurement::herfindahl_hirschman(row);
    P.equivalent[i] = 1 / P.hhi[i];
    const auto top = static_cast<std::size_t>(std::max_element(row.begin(), row.end()) - row.begin());
    P.first[i] = top;
    P.first_share[i] = row[top] / total;
  }
  return P;
}

// I_ij, NaN where i exports nothing or j imports nothing.
inline Rectangular intensity(const Rectangular& x) {
  detail::require_square(x);
  return nomos::measurement::location_quotients(x);
}

// GL_ij for every pair (symmetric), NaN where the pair does not trade.
inline Rectangular two_way(const Rectangular& x) {
  detail::require_square(x);
  const std::size_t n = x.rows;
  Rectangular g(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      g(i, j) = i != j && x(i, j) + x(j, i) > 0 ? nomos::measurement::grubel_lloyd(x(i, j), x(j, i)) : std::numeric_limits<double>::quiet_NaN();
  return g;
}

// The world's share of trade that goes both ways: Σ (x_ij + x_ji) GL_ij over pairs, over W.
inline double reciprocal_share(const Rectangular& x) {
  detail::require_square(x);
  double both = 0, total = 0;
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t j = 0; j < x.cols; j++) {
      total += x(i, j);
      if (i < j) both += 2 * std::min(x(i, j), x(j, i));
    }
  return total > 0 ? both / total : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace economics::trade
