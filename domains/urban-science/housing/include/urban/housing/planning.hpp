// The planner's decision: how much new supply to permit, when demand is seen only through rents (ADR 0071 §8, the
// design's mixed observability).
//
// A MOMDP (tyche decision; Ong, Png, Hsu and Lee 2010) whose state is (x, y):
// - x, seen: the stock, in units of supply above a reference (0 … X − 1);
// - y, hidden: the demand, in units of households who consider the centre (0 … Y − 1), a Markov chain T_Y estimated
//   from the history of demand (the study reads it from the rents);
// - the action a: the units permitted this period (0 … A − 1), completed by the next: x' = min(x + a, X − 1), the lag
//   of construction;
// - the observation: the band of the rent at which the market clears at (x', y'), read with error (a neighbouring band
//   with probability `confusion` each side), so that the rent reveals demand only in part and only through the stock;
// - the reward: the period's welfare at (x, y) (the market's surplus, its rents and the cost of empty dwellings, as
//   the caller computes them from urban::housing::market) less the cost of the units permitted.
// Beliefs live on {x} × Δ(Y). The solvers are tyche's, on `to_pomdp` of this model: HSVI's interval for the optimum,
// and the simple policies below evaluated by simulation.
//
// Caveat. The planner here permits; it does not build. Its welfare weighs a household's surplus as one dollar, and
// the costs of construction are the caller's assumptions.

#pragma once

#include <urban/housing/market.hpp>

#include <tyche/decision/pomdp.hpp>
#include <tyche/decision/structured.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace urban::housing {

struct PlannerSpec {
  std::size_t supply = 0, demand = 0, actions = 0, bands = 3;
  std::vector<double> demand_transition;  // T_Y(y' | y) at [y Y + y']
  std::vector<double> welfare;            // at [x Y + y], a period
  std::vector<double> rent;               // the clearing rent at [x Y + y]
  std::vector<double> cuts;               // bands − 1 increasing thresholds of rent
  double confusion = 0.1;
  double cost = 0;                        // a unit permitted
  double discount = 0.9;
  std::size_t x0 = 0;
  std::vector<double> y0;                 // the belief about demand at the start
};

inline std::size_t band_of(const PlannerSpec& P, double rent) {
  std::size_t b = 0;
  while (b < P.cuts.size() && rent >= P.cuts[b]) b++;
  return b;
}

inline tyche::decision::Momdp planner(const PlannerSpec& P) {
  const std::size_t X = P.supply, Y = P.demand, A = P.actions, O = P.bands;
  detail::require(X > 0 && Y > 0 && A > 0 && O > 0, "the planner needs supply, demand, actions and bands");
  detail::require(P.demand_transition.size() == Y * Y && P.welfare.size() == X * Y && P.rent.size() == X * Y, "T_Y is Y·Y, welfare and rent X·Y");
  detail::require(P.cuts.size() + 1 == O, "bands − 1 cuts");
  tyche::decision::Momdp M;
  M.visible = X, M.hidden = Y, M.actions = A, M.observations = O, M.discount = P.discount, M.x0 = P.x0, M.y0 = P.y0;
  const std::size_t S = X * Y;
  M.T.assign(A * S * S, 0.0), M.Z.assign(A * S * O, 0.0), M.R.assign(A * S, 0.0);
  for (std::size_t a = 0; a < A; a++)
    for (std::size_t x = 0; x < X; x++)
      for (std::size_t y = 0; y < Y; y++) {
        const std::size_t s = x * Y + y, x2 = std::min(x + a, X - 1);
        for (std::size_t y2 = 0; y2 < Y; y2++) M.T[(a * S + s) * S + x2 * Y + y2] = P.demand_transition[y * Y + y2];
        M.R[a * S + s] = P.welfare[s] - P.cost * static_cast<double>(a);
        // What is read of the state entered (the same (x, y) entered by any action).
        const std::size_t band = band_of(P, P.rent[s]);
        for (std::size_t o = 0; o < O; o++) {
          const double up = band > 0 ? P.confusion : 0, down = band + 1 < O ? P.confusion : 0;
          M.Z[(a * S + s) * O + o] = o == band ? 1 - up - down : (o + 1 == band ? up : (o == band + 1 ? down : 0.0));
        }
      }
  return M;
}

// The seen supply and the belief about demand of a belief of to_pomdp(M) (its mass lies on one x).
inline tyche::decision::MomdpBelief split(const tyche::decision::Momdp& M, std::span<const double> b) {
  tyche::decision::MomdpBelief out{0, std::vector<double>(M.hidden, 0.0)};
  double best = -1;
  for (std::size_t x = 0; x < M.visible; x++) {
    double mass = 0;
    for (std::size_t y = 0; y < M.hidden; y++) mass += b[x * M.hidden + y];
    if (mass > best) best = mass, out.x = x;
  }
  for (std::size_t y = 0; y < M.hidden; y++) out.y[y] = b[out.x * M.hidden + y] / best;
  return out;
}

// Simple policies to set against the plan.
// - never: permit nothing; steady: the same units every period;
// - follow the rents: permit the most when the expected band of rent (under the belief, at the supply seen) is the
//   highest, nothing when it is the lowest, one unit between (a planner who reacts to rents as they are read);
// - certainty equivalent: act as the MDP that sees the state would at the most likely demand.
inline std::size_t policy_follow_rents(const PlannerSpec& P, const tyche::decision::MomdpBelief& b) {
  double band = 0;
  for (std::size_t y = 0; y < P.demand; y++) band += b.y[y] * static_cast<double>(band_of(P, P.rent[b.x * P.demand + y]));
  const double top = static_cast<double>(P.bands - 1);
  if (band >= top - 0.5) return P.actions - 1;
  if (band <= 0.5) return 0;
  return std::min<std::size_t>(1, P.actions - 1);
}

inline std::size_t policy_certainty_equivalent(const tyche::decision::MdpSolution& mdp, std::size_t Y, const tyche::decision::MomdpBelief& b) {
  const std::size_t y = static_cast<std::size_t>(std::max_element(b.y.begin(), b.y.end()) - b.y.begin());
  return mdp.policy[b.x * Y + y];
}

}  // namespace urban::housing
