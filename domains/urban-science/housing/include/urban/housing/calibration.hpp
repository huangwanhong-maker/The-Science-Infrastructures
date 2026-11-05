// What the aggregates reveal of the households' preferences: the inverse demand of the centre's segments, read from
// the rents of new lettings and the stock (ADR 0075).
//
// With the market's logit (market.hpp) and one class, a segment's occupied stock satisfies
//   ln O_jq = ln M_q − ln(1 + Σ_k e^{V_kq}) + δ_j − α ln r_jq + ξ_jq
// in period q. The period's terms are common to its segments, so that
//   ln r_jq = a_q + b_j − (1/α) ln S_jq + e_jq,
// a regression of the log rent on the log stock with period and segment effects, whose slope is −1/α: by how much
// a segment's rent falls, against the others' in the same period, as its stock grows. Its period effects a_q carry the
// demand common to the centre (the households who consider it, and the outside option's draw): a_q α is the log of
// the demand index, up to a constant. tyche's least squares gives the slope with heteroskedasticity-consistent (HC1)
// and segment-clustered (CR1) errors.
//
// The assumption that identifies α: a segment's stock is predetermined, decided years before it is let (towers are
// planned, approved and built over three to six years), so that it is not correlated with the period's shock to that
// segment's demand ξ_jq beyond what the segment's effect holds. Developers build where demand has grown, so the
// assumption is a model, declared. Overlapping years of the moving annual medians make e serially correlated; the
// clustered error allows it.

#pragma once

#include <urban/housing/market.hpp>

#include <tyche/regression/linear.hpp>

#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace urban::housing {

struct RentPanel {
  std::vector<std::size_t> segment, period;
  std::vector<double> log_rent, log_stock;
  std::size_t segments = 0, periods = 0;
};

struct InverseDemand {
  double slope = 0;             // ∂ ln r / ∂ ln S
  double se_hc1 = 0, se_cluster = 0;
  double alpha = 0;             // −1 / slope
  double alpha_se = 0;          // by the delta method, from the clustered error
  std::vector<double> segment_effect;  // b_j
  std::vector<double> period_effect;   // a_q, the first period's 0 (periods without observations NaN)
  double r_squared = 0;
  std::size_t observations = 0, clusters = 0;
};

inline InverseDemand inverse_demand(const RentPanel& P) {
  const std::size_t n = P.log_rent.size(), J = P.segments, Q = P.periods;
  detail::require(n == P.log_stock.size() && n == P.segment.size() && n == P.period.size() && n > 0, "a panel of equal columns");
  // The periods that hold observations; the first of them is the reference.
  std::vector<int> used(Q, 0);
  for (const std::size_t q : P.period) used[q] = 1;
  std::vector<long> column(Q, -1);
  std::size_t k = 1 + J;
  bool first = true;
  for (std::size_t q = 0; q < Q; q++)
    if (used[q]) {
      if (first) { first = false; continue; }
      column[q] = static_cast<long>(k++);
    }
  calx::linalg::Rectangular X(n, k);
  for (std::size_t i = 0; i < n; i++) {
    X(i, 0) = P.log_stock[i];
    X(i, 1 + P.segment[i]) = 1;
    if (column[P.period[i]] >= 0) X(i, static_cast<std::size_t>(column[P.period[i]])) = 1;
  }
  tyche::regression::OlsOptions o;
  o.covariance = tyche::regression::Covariance::hc1;
  const auto hc = tyche::regression::ols(X, P.log_rent, o);
  o.covariance = tyche::regression::Covariance::cluster;
  for (const std::size_t j : P.segment) o.clusters.push_back(static_cast<std::int64_t>(j));
  const auto cl = tyche::regression::ols(X, P.log_rent, o);
  InverseDemand out;
  out.slope = hc.beta[0], out.se_hc1 = hc.se[0], out.se_cluster = cl.se[0];
  out.alpha = -1 / out.slope, out.alpha_se = out.se_cluster / (out.slope * out.slope);
  out.segment_effect.assign(hc.beta.begin() + 1, hc.beta.begin() + 1 + static_cast<long>(J));
  out.period_effect.assign(Q, std::nan(""));
  for (std::size_t q = 0; q < Q; q++)
    if (used[q]) out.period_effect[q] = column[q] >= 0 ? hc.beta[static_cast<std::size_t>(column[q])] : 0.0;
  out.r_squared = hc.r_squared, out.observations = n, out.clusters = cl.clusters;
  return out;
}

// The quarterly hazard h whose four quarters give the year's share p of moves: 1 − (1 − h)^4 = p.
inline double quarterly_hazard(double yearly_share) { return 1 - std::pow(1 - yearly_share, 0.25); }

}  // namespace urban::housing
