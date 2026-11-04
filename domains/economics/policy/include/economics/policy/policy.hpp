// Decisions on the perceived regime of an economy: a policymaker's POMDP whose state is the regime, whose observations
// are the year's growth, and whose belief is the regimes' filter; a fund that supports economies under a budget; a
// budget on support itself.
//
// The data give the regimes (economics/regimes: the transition matrix and each regime's law of growth). What policy does
// and what it is worth are assumptions, declared in `PolicyAssumptions` and carried into every result:
// - actions: wait (0), support (1: fiscal or monetary easing), restrain (2: tightening);
// - support moves a share `support_effect` of a recession's persistence to normal growth, and in a boom a share
//   `overheating` of the boom's return to normal to the boom itself; restraint moves a share `restraint_effect` of a
//   boom's persistence to normal growth, and in a recession a share `deepening` of the recovery back to recession;
// - a year's welfare by regime, less the action's cost.
// The observation is the next year's growth, in bins with edges `edges` (in the regimes' units: each economy's own
// deviations from its normal when the regimes are standardised, else per cent), each regime's probability of a bin
// that of its fitted normal law. The POMDP is then tyche's; its belief, updated by Bayes's rule on the binned growth, is
// the regimes' filter coarsened to the bins.
//
// The fund's model keeps two actions (wait, support), each economy one component of tyche's weakly coupled problem,
// support costing one unit of a year's budget (Adelman and Mersereau 2008); the budget on support is a constrained POMDP
// whose cost is one each year support is given (Altman 1999).

#pragma once

#include <tyche/decision/pomdp.hpp>
#include <tyche/processes/hidden_markov.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace economics {

struct PolicyAssumptions {
  std::vector<double> welfare{-4.0, 1.0, -0.5};  // a year in contraction, normal growth, boom
  std::array<double, 3> cost{0.0, 0.15, 0.2};    // wait, support, restrain
  double support_effect = 0.5, overheating = 0.2;
  double restraint_effect = 0.35, deepening = 0.2;
  std::vector<double> edges{-1.5, -0.75, 0, 0.75, 1.5};  // bins of growth, in the regimes' units (each economy's deviations)
  double discount = 0.9;
};

enum PolicyAction : std::size_t { wait = 0, support = 1, restrain = 2 };

namespace detail {
inline double normal_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Moves a share `f` of row i's probability of `from` to `to`.
inline void shift(std::vector<double>& row, std::size_t from, std::size_t to, double f) {
  const double m = f * row[from];
  row[from] -= m, row[to] += m;
}
}  // namespace detail

// The probability of each growth bin in each regime.
inline std::vector<std::vector<double>> growth_bins(const tyche::processes::GaussianHmm& R, const std::vector<double>& edges) {
  std::vector<std::vector<double>> Z(R.states, std::vector<double>(edges.size() + 1));
  for (std::size_t k = 0; k < R.states; k++)
    for (std::size_t o = 0; o <= edges.size(); o++) {
      const double lo = o == 0 ? 0.0 : detail::normal_cdf((edges[o - 1] - R.mean[k]) / R.sd[k]);
      const double hi = o == edges.size() ? 1.0 : detail::normal_cdf((edges[o] - R.mean[k]) / R.sd[k]);
      Z[k][o] = hi - lo;
    }
  return Z;
}

inline std::size_t growth_bin(double g, const std::vector<double>& edges) {
  std::size_t o = 0;
  while (o < edges.size() && g >= edges[o]) o++;
  return o;
}

// The policymaker's POMDP on three regimes (recession 0, normal 1, boom 2); `actions` 2 keeps wait and support (the
// fund's model), 3 all three.
inline tyche::decision::Pomdp policy_pomdp(const tyche::processes::GaussianHmm& R, const PolicyAssumptions& A, std::size_t actions = 3) {
  if (R.states != 3) throw std::invalid_argument("economics::policy_pomdp: three regimes (recession, normal, boom)");
  if (A.welfare.size() != 3) throw std::invalid_argument("economics::policy_pomdp: a welfare for each regime");
  const std::size_t S = 3, O = A.edges.size() + 1;
  tyche::decision::Pomdp P(S, actions, O, A.discount);
  const auto Z = growth_bins(R, A.edges);
  for (std::size_t a = 0; a < actions; a++)
    for (std::size_t s = 0; s < S; s++) {
      std::vector<double> row(S);
      for (std::size_t s2 = 0; s2 < S; s2++) row[s2] = R.a(s, s2);
      if (a == support && s == 0) detail::shift(row, 0, 1, A.support_effect);
      if (a == support && s == 2) detail::shift(row, 1, 2, A.overheating);
      if (a == restrain && s == 2) detail::shift(row, 2, 1, A.restraint_effect);
      if (a == restrain && s == 0) detail::shift(row, 1, 0, A.deepening);
      for (std::size_t s2 = 0; s2 < S; s2++) P.transition(a, s, s2) = row[s2];
      for (std::size_t o = 0; o < O; o++) P.observe(a, s, o) = Z[s][o];
      P.reward(s, a) = A.welfare[s] - A.cost[a];
    }
  P.initial = {1.0 / 3, 1.0 / 3, 1.0 / 3};
  return P;
}

// The cost of support, one each year it is given (for the budget on support).
inline std::vector<double> support_cost(const tyche::decision::Pomdp& P) {
  std::vector<double> C(P.actions * P.states, 0.0);
  for (std::size_t s = 0; s < P.states; s++) C[support * P.states + s] = 1;
  return C;
}

// An economy's belief about its regime after its last observed year, carried `ahead` years by the chain (waiting).
inline std::vector<double> belief_ahead(const tyche::processes::GaussianHmm& R, const std::vector<double>& filtered, int ahead = 1) {
  std::vector<double> b = filtered, n(R.states);
  for (int t = 0; t < ahead; t++) {
    for (std::size_t j = 0; j < R.states; j++) {
      n[j] = 0;
      for (std::size_t i = 0; i < R.states; i++) n[j] += b[i] * R.a(i, j);
    }
    b = n;
  }
  return b;
}

}  // namespace economics
