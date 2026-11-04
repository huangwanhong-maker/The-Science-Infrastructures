// Planning the maintenance of a city's surfaces under a budget: each surface a POMDP on its fitted deterioration, the
// city their weakly coupled sum.
//
// A surface (per square metre) is a POMDP (tyche decision) on the five classes:
// - actions: do nothing (0), inspect (1), maintain (2: one class better), renew (3: class 5), each followed by a year's
//   wear exp(Q), Q the type's fitted chain;
// - observations: an inspection reads the class entered, or a neighbouring class with probability `confusion` on each
//   side; the other actions read nothing (each of the five symbols equally);
// - reward: the service of the class held (good 1 … failed −1, a year) less the action's cost;
// - the belief now: from the last rating, read as an inspection, carried by the wear of the years since.
// The city is the sum of its surfaces, each weighted by its area, joined only by the year's budget: tyche's Lagrangian
// bound (Adelman and Mersereau 2008) is an upper bound on the value any plan within the budget can reach, and its
// ranked plan, which actions to take on which surfaces this year, keeps the budget (Madanat and Ben-Akiva 1994;
// Papakonstantinou and Shinozuka 2014 plan inspection and repair of one structure so).
//
// Caveat. The service values, costs and the noise of an inspection are assumptions, declared in `Assumptions` and
// carried into every result: the data give the deterioration, not the City's values or prices.

#pragma once

#include <urban/maintenance/deterioration.hpp>

#include <tyche/decision/pomdp.hpp>
#include <tyche/decision/structured.hpp>
#include <tyche/processes/panel_markov.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace urban::maintenance {

struct Assumptions {
  std::array<double, 5> service{1.0, 0.8, 0.5, 0.2, -1.0};  // a year a square metre, in classes 5 … 1
  std::array<double, 4> cost{0.0, 0.02, 0.3, 1.2};          // a square metre: nothing, inspect, maintain, renew
  double confusion = 0.1;                                   // an inspection's chance of each neighbouring class
  double discount = 0.95;                                   // a year
};

enum Action : std::size_t { nothing = 0, inspect = 1, maintain = 2, renew = 3 };

inline tyche::decision::Pomdp surface_model(const Deterioration& D, const Assumptions& A) {
  tyche::decision::Pomdp P(5, 4, 5, A.discount);
  std::vector<std::vector<double>> wear(5);
  for (std::size_t s = 0; s < 5; s++) wear[s] = tyche::processes::transition_row(D.generator(), s, 1.0);
  for (std::size_t a = 0; a < 4; a++)
    for (std::size_t s = 0; s < 5; s++) {
      const std::size_t after = a == maintain ? (s == 0 ? 0 : s - 1) : a == renew ? 0 : s;
      for (std::size_t s2 = 0; s2 < 5; s2++) P.transition(a, s, s2) = wear[after][s2];
      for (std::size_t o = 0; o < 5; o++) {
        double z = 1.0 / 5;
        if (a == inspect) {
          const double up = s > 0 ? A.confusion : 0, down = s + 1 < 5 ? A.confusion : 0;  // read one class better, worse
          z = o == s ? 1 - up - down : (o + 1 == s ? up : (o == s + 1 ? down : 0.0));
        }
        P.observe(a, s, o) = z;
      }
      P.reward(s, a) = A.service[s] - A.cost[a];
    }
  return P;
}

// The belief now of a surface last rated in class `cls`, `since` years ago: the rating read as an inspection (from a
// uniform prior), then the years' wear.
inline std::vector<double> belief_now(const Deterioration& D, int cls, double since, const Assumptions& A) {
  std::vector<double> b(5, 0.0);
  const std::size_t o = state_of(cls);
  double total = 0;
  for (std::size_t s = 0; s < 5; s++) {
    const double up = s > 0 ? A.confusion : 0, down = s + 1 < 5 ? A.confusion : 0;
    const double z = o == s ? 1 - up - down : (o + 1 == s ? up : (o == s + 1 ? down : 0.0));
    total += (b[s] = z);
  }
  for (double& x : b) x /= total;
  if (since <= 0) return b;
  std::vector<double> out(5, 0.0);
  for (std::size_t s = 0; s < 5; s++)
    if (b[s] > 0) {
      const auto row = tyche::processes::transition_row(D.generator(), s, since);
      for (std::size_t k = 0; k < 5; k++) out[k] += b[s] * row[k];
    }
  return out;
}

struct CityPlan {
  double lambda = 0, upper = 0;       // the Lagrangian bound on the discounted value of any plan within the budget
  double budget = 0, spent = 0;       // this year's
  std::vector<std::size_t> action;    // one per surface
  std::map<std::string, std::array<int, 4>> counts;  // by type: nothing, inspect, maintain, renew
  double nothing_value = 0;           // the value of doing nothing, exactly (a feasible plan)
  int evaluations = 0;
};

inline CityPlan plan_city(const std::vector<Rated>& surfaces, const std::map<std::string, Deterioration>& fitted, const Assumptions& A, double budget,
                          double epsilon = 0.05, int max_trials = 150) {
  std::map<std::string, tyche::decision::Pomdp> models;
  for (const auto& [type, D] : fitted) models.emplace(type, surface_model(D, A));
  std::vector<tyche::decision::Component> C;
  std::vector<std::size_t> kept;
  CityPlan out;
  out.budget = budget;
  const std::vector<double> cost(A.cost.begin(), A.cost.end());
  for (std::size_t i = 0; i < surfaces.size(); i++) {
    const auto& r = surfaces[i];
    const auto it = models.find(r.type);
    if (it == models.end() || !(r.area > 0)) continue;
    C.push_back({&it->second, cost, belief_now(fitted.at(r.type), r.cls, r.since, A), r.area});
    kept.push_back(i);
  }
  const auto sol = tyche::decision::lagrangian_bound(C, budget, {.epsilon = epsilon, .max_trials = max_trials});
  const auto act = tyche::decision::coupled_actions(sol, C, budget);
  out.lambda = sol.lambda, out.upper = sol.upper, out.evaluations = sol.evaluations;
  out.action.assign(surfaces.size(), nothing);
  for (std::size_t k = 0; k < C.size(); k++) {
    out.action[kept[k]] = act[k];
    out.spent += C[k].weight * cost[act[k]];
    out.counts[surfaces[kept[k]].type][act[k]]++;
  }
  // Doing nothing for ever: Σ_i w_i b_i·α, α = R(·, nothing) + γ T α, exactly.
  std::map<std::string, std::vector<double>> alpha;
  for (const auto& [type, P] : models) {
    std::vector<double> v(5, 0.0), next(5);
    for (int it = 0; it < 20000; it++) {
      double change = 0;
      for (std::size_t s = 0; s < 5; s++) {
        double e = 0;
        for (std::size_t s2 = 0; s2 < 5; s2++) e += P.transition(nothing, s, s2) * v[s2];
        next[s] = P.reward(s, nothing) + P.discount * e;
        change = std::max(change, std::fabs(next[s] - v[s]));
      }
      v.swap(next);
      if (change < 1e-12) break;
    }
    alpha[type] = v;
  }
  for (std::size_t k = 0; k < C.size(); k++) {
    const auto& a = alpha[surfaces[kept[k]].type];
    for (std::size_t s = 0; s < 5; s++) out.nothing_value += C[k].weight * C[k].belief[s] * a[s];
  }
  return out;
}

}  // namespace urban::maintenance
