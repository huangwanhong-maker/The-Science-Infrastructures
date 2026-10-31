// Terrain (theory/street-oceanography.md §2.1–2.2): a field of components in unlike units, read by each
// arrival through a weighting; the order that survives every weighting; and the length over which the
// uses of a street change, with the compression ratio built on it.
//
// Contraction. τ_a = ⟨w_a, T⟩ (Huang, Street Oceanography, eq. 3). A Contraction keeps its weighting: no
// scalar terrain exists here without the w that made it (Theorem 2.2: no nonzero w is invariant under
// changes of the components' units).
//
// Dominance. T ≼ T′ for every non-negative weighting exactly when T_i ≤ T′_i for every i (Theorem 2.2b);
// covering() returns the Hasse diagram of that order among items (holdings, segments, streets): i ⋖ j when
// i is dominated by j with nothing strictly between. Where the signs of the components' effects are not
// fixed (Theorem 2.2c) nothing is comparable, and the caller says so by not asking.
//
// Use correlation. For points with uses, ρ(r) = (P_same(r) − P_∞)/(1 − P_∞), P_same(r) the share of pairs at
// separation in the bin of r that share a use, P_∞ = Σ_u p_u²; L_grad = ∫₀^{r₀} ρ dr up to its first zero r₀
// (Corrsin's integral scale; Tennekes and Lumley 1972, §8.1), by the trapezoidal rule on the bins' centres,
// ρ(0) = 1 by definition, each bin placed at the mean separation of its pairs. Positions in one dimension
// (s along the street) or two (x, y in a local plane).
// Accuracy: the bins' width bounds the resolution of L_grad; pairs are counted exactly, O(n²).
//
// Compression. C = L_walk / L_grad (Street Oceanography, eq. 4); L_walk is declared (Remark 2.5).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace urban::street {

struct Weighting {
  std::string name;              // whose reading: "a visitor on a budget", "a resident", …
  std::vector<double> w;         // one weight per component, in the inverse of its unit
};

struct Contraction {
  Weighting weighting;           // the w that made these numbers
  std::vector<double> value;     // τ per item
  std::vector<std::size_t> rank; // items by τ, highest first
};

// τ_i = ⟨w, T_i⟩ for items i (rows of terrain), and their ranking.
inline Contraction contract(const std::vector<std::vector<double>>& terrain, const Weighting& w) {
  Contraction c{w, {}, {}};
  for (const auto& row : terrain) {
    if (row.size() != w.w.size()) throw std::invalid_argument("urban::street::contract: one weight per component");
    double t = 0;
    for (std::size_t i = 0; i < row.size(); i++) t += w.w[i] * row[i];
    c.value.push_back(t);
  }
  c.rank.resize(c.value.size());
  std::iota(c.rank.begin(), c.rank.end(), std::size_t{0});
  std::stable_sort(c.rank.begin(), c.rank.end(), [&](std::size_t a, std::size_t b) { return c.value[a] > c.value[b]; });
  return c;
}

// Whether a is dominated by b: a_i ≤ b_i for every component, and a ≠ b.
inline bool dominated(const std::vector<double>& a, const std::vector<double>& b) {
  bool strict = false;
  for (std::size_t i = 0; i < a.size(); i++) {
    if (a[i] > b[i]) return false;
    if (a[i] < b[i]) strict = true;
  }
  return strict;
}

struct Cover { std::size_t below = 0, above = 0; };

// The covering pairs of the dominance order (its Hasse diagram), O(n³) in the items.
inline std::vector<Cover> covering(const std::vector<std::vector<double>>& terrain) {
  const std::size_t n = terrain.size();
  std::vector<std::vector<char>> less(n, std::vector<char>(n, 0));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) less[i][j] = i != j && dominated(terrain[i], terrain[j]);
  std::vector<Cover> out;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      if (!less[i][j]) continue;
      bool between = false;
      for (std::size_t k = 0; k < n && !between; k++) between = less[i][k] && less[k][j];
      if (!between) out.push_back({i, j});
    }
  return out;
}

// The items no other item dominates (the maximal elements; the Pareto front).
inline std::vector<std::size_t> maximal(const std::vector<std::vector<double>>& terrain) {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < terrain.size(); i++) {
    bool dom = false;
    for (std::size_t j = 0; j < terrain.size() && !dom; j++) dom = dominated(terrain[i], terrain[j]);
    if (!dom) out.push_back(i);
  }
  return out;
}

// The share of ordered pairs of items that two weightings rank the same way (pairs tied under either are
// left out): 1 when they agree on every pair, as they must on dominated pairs.
inline double agreement(const Contraction& a, const Contraction& b) {
  std::size_t same = 0, total = 0;
  for (std::size_t i = 0; i < a.value.size(); i++)
    for (std::size_t j = i + 1; j < a.value.size(); j++) {
      const double da = a.value[i] - a.value[j], db = b.value[i] - b.value[j];
      if (da == 0 || db == 0) continue;
      total++;
      if ((da > 0) == (db > 0)) same++;
    }
  return total ? static_cast<double>(same) / static_cast<double>(total) : 1.0;
}

// ---- the use correlation and the length of the terrain ------------------------------------------------

struct UseCorrelation {
  std::vector<double> r;          // bin centres, m (r[0] = 0 is the definition ρ(0) = 1)
  std::vector<double> rho;        // ρ(r); NaN in bins with no pairs
  std::vector<std::size_t> pairs; // pairs per bin
  double p_infinity = 0;          // Σ p_u²
  double length = 0;              // L_grad, m
  double first_zero = 0;          // r₀, m (the last bin's edge when ρ stays positive: then L_grad is a lower bound)
  bool single_regime = false;     // one use only: P_∞ = 1, ρ undefined, L_grad = 0
  bool reached_zero = false;
};

// Points (x, y) with uses (y empty for positions along a street). Bins of width `bin` up to r_max.
inline UseCorrelation use_correlation(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::size_t>& use,
                                      double bin, double r_max) {
  const std::size_t n = x.size();
  if (use.size() != n || (!y.empty() && y.size() != n)) throw std::invalid_argument("urban::street::use_correlation: one position and one use per point");
  if (!(bin > 0 && r_max > bin)) throw std::invalid_argument("urban::street::use_correlation: a positive bin narrower than the range");
  UseCorrelation out;
  if (n < 2) { out.single_regime = true; return out; }
  std::size_t kinds = 0;
  for (auto u : use) kinds = std::max(kinds, u + 1);
  std::vector<double> share(kinds, 0);
  for (auto u : use) share[u] += 1.0 / static_cast<double>(n);
  for (double p : share) out.p_infinity += p * p;
  if (out.p_infinity >= 1 - 1e-12) { out.single_regime = true; return out; }
  const std::size_t bins = static_cast<std::size_t>(std::ceil(r_max / bin));
  std::vector<double> same(bins, 0), all(bins, 0), sum_r(bins, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const double dx = x[i] - x[j], dy = y.empty() ? 0.0 : y[i] - y[j];
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d >= r_max) continue;
      const std::size_t b = std::min(bins - 1, static_cast<std::size_t>(d / bin));
      all[b] += 1;
      sum_r[b] += d;
      if (use[i] == use[j]) same[b] += 1;
    }
  out.r.push_back(0);
  out.rho.push_back(1);
  out.pairs.push_back(0);
  for (std::size_t b = 0; b < bins; b++) {
    // Each bin at the mean separation of its pairs (its centre when empty): exact for positions on a lattice.
    out.r.push_back(all[b] > 0 ? sum_r[b] / all[b] : (static_cast<double>(b) + 0.5) * bin);
    out.pairs.push_back(static_cast<std::size_t>(all[b]));
    out.rho.push_back(all[b] > 0 ? (same[b] / all[b] - out.p_infinity) / (1 - out.p_infinity) : std::numeric_limits<double>::quiet_NaN());
  }
  // ∫ ρ dr by trapezoids over the bins with pairs, up to the first sign change (linearly interpolated).
  double area = 0, r_prev = 0, rho_prev = 1;
  out.first_zero = out.r.back() + 0.5 * bin;
  for (std::size_t k = 1; k < out.r.size(); k++) {
    if (!std::isfinite(out.rho[k])) continue;
    const double r = out.r[k], rho = out.rho[k];
    if (rho <= 0) {
      const double r0 = r_prev + (r - r_prev) * rho_prev / (rho_prev - rho);
      area += 0.5 * rho_prev * (r0 - r_prev);
      out.first_zero = r0;
      out.reached_zero = true;
      break;
    }
    area += 0.5 * (rho_prev + rho) * (r - r_prev);
    r_prev = r, rho_prev = rho;
  }
  out.length = area;
  return out;
}

// C = L_walk / L_grad; NaN for a single regime (L_grad = 0).
inline double compression(double walk, const UseCorrelation& c) {
  if (c.single_regime || !(c.length > 0)) return std::numeric_limits<double>::quiet_NaN();
  return walk / c.length;
}

}  // namespace urban::street
