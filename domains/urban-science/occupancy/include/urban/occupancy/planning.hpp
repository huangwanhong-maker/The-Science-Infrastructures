// The planner's interventions on vacant space as a weakly coupled POMDP over the buildings (the design's §1.6 and §10;
// theory/occupancy.md, §5), on tyche's `decision`.
//
// A building (the planner's model). Its hidden state is its lump, occupied (O), vacant and lettable (V) or vacant under
// works (W), whose law is the model's lumped chain with its neighbours and the demand held at the last census
// (Proposition 2 of model.hpp; the held neighbourhood is the decoupling, declared: the agent-based model itself keeps
// the neighbours' feedback, and the policies are evaluated in it too). A year is a step:
// - actions: nothing (0), inspect (1: the next census surveys it for certain), incentive (2: if it is V, it is offered
//   to occupiers at the extra rate κ λ_0 for the year);
// - transitions: exp(Q_a), Q_a the lumped generator with the building enrolled (a = 2) or not;
// - observations: the census surveys a building in a year with chance ρ (CLUE surveys every commercial property at least
//   once every two years: ρ = 1/2 is the declared default), and then reads its lump; else it reads nothing (symbol 3);
//   an inspection reads it surely;
// - reward: minus the expected years in V within the year from the state at its start, ∫_0^1 exp(t Q_a)(s, V) dt
//   (Simpson's rule on 64 intervals, its error below 10⁻⁹ for rates under 10 a year);
// - costs (the budget's units): nothing 0, inspect 0.1, incentive 1 (an incentive-year), declared.
// Works are not the planner's target: an incentive does nothing to a building under renovation or construction.
//
// The city. The buildings' sum, joined only by the year's budget B: tyche's Lagrangian bound (Adelman and Mersereau
// 2008) is an upper bound on the discounted value of any plan within the budget, and its ranked policy keeps the
// budget. Buildings whose lumped models agree are one group, solved once. Against it, four simple rules, each enrolling
// buildings last seen vacant and lettable while the budget allows: every vacancy (in a random order), the longest
// vacant first, the most vacant neighbourhood first; and doing nothing, valued exactly.
//
// Two evaluations of a policy:
// - in the planner's own decoupled model (episodes of its lumped chains, the bound's problem): the policy's value lies
//   below the bound, and the two make an interval for the best plan's;
// - in the agent-based model (Gillespie over the uses, the neighbours' feedback kept): the years of lettable vacancy
//   the policy averts, against doing nothing on common random numbers.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <string>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>
#include <nomos/agents/engine.hpp>
#include <tyche/decision/pomdp.hpp>
// tyche's structured.hpp (its line 164, CoupledGroup built from two of its five members) draws -Wmissing-field-initializers
// under -Wextra; benign, and reported to tyche. Silenced here only, around the include.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#include <tyche/decision/structured.hpp>
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#include <tyche/processes/panel_markov.hpp>
#include <urban/occupancy/model.hpp>

namespace urban::occupancy {

namespace td = tyche::decision;

struct PlanAssumptions {
  double survey = 0.5;                       // ρ, a census's chance of surveying a building in a year
  std::array<double, 3> cost{0.0, 0.1, 1.0}; // nothing, inspect, incentive (incentive-years)
  double discount = 0.95;                    // a year
};

enum Action : std::size_t { nothing = 0, inspect = 1, incentive = 2 };
constexpr std::size_t not_surveyed = 3;

using Generator3 = std::array<std::array<double, lumps>, lumps>;

inline std::vector<double> lump_row(const Generator3& Q, std::size_t s, double t) {
  calx::linalg::Matrix q(lumps);
  for (std::size_t a = 0; a < lumps; a++)
    for (std::size_t b = 0; b < lumps; b++) q(a, b) = Q[a][b];
  return tyche::processes::transition_row(q, s, t);
}

// The POMDP of a building from its two lumped generators (not enrolled, enrolled).
inline td::Pomdp building_model(const Generator3& Q0, const Generator3& Q1, const PlanAssumptions& A) {
  td::Pomdp P(lumps, 3, lumps + 1, A.discount);
  for (std::size_t a = 0; a < 3; a++) {
    const auto& Q = a == incentive ? Q1 : Q0;
    for (std::size_t s = 0; s < lumps; s++) {
      const auto row = lump_row(Q, s, 1.0);
      for (std::size_t s2 = 0; s2 < lumps; s2++) P.transition(a, s, s2) = row[s2];
      const int m = 64;  // Simpson's rule for the years in V
      double v = 0;
      for (int k = 0; k <= m; k++) {
        const double w = (k == 0 || k == m) ? 1 : (k % 2 ? 4 : 2);
        v += w * lump_row(Q, s, static_cast<double>(k) / m)[vacant];
      }
      P.reward(s, a) = -v / (3.0 * m);
    }
    for (std::size_t s2 = 0; s2 < lumps; s2++) {
      if (a == inspect) P.observe(a, s2, s2) = 1;
      else P.observe(a, s2, s2) = A.survey, P.observe(a, s2, not_surveyed) = 1 - A.survey;
    }
  }
  P.initial = {1, 0, 0};
  return P;
}

inline td::Pomdp building_model(const Model& M, const State& s, std::size_t i, const PlanAssumptions& A) {
  return building_model(lumped_generator(M, s, i, false), lumped_generator(M, s, i, true), A);
}

// The belief a year (or `years`) after the lump was seen, doing nothing.
inline td::Belief belief_after(const td::Pomdp& P, std::size_t lump, int years = 1) {
  td::Belief b(lumps, 0.0);
  b[lump] = 1;
  for (int k = 0; k < years; k++) b = td::predict(P, b, nothing);
  return b;
}

// What the planner knows of each building.
struct Knowledge {
  std::vector<int> lump;               // last lump seen (−1 none)
  std::vector<int> year;               // when
  std::vector<int> since;              // the first year of the run of V sightings it is in (when lump is V)
  std::vector<td::Belief> belief;      // the POMDP's, by building
};

struct CityPlan {
  std::vector<td::Pomdp> models;                // distinct
  std::vector<std::size_t> model_of;            // by building (models.size(): not planned)
  std::vector<std::size_t> planned;             // the buildings planned, in order of their components
  std::vector<td::Component> components;        // one a planned building, its belief at the start
  td::CoupledSolution solution;
  double budget = 0;
  double nothing_value = 0;                     // Σ_i b_i·α_nothing, exactly (a feasible plan)
  std::vector<std::vector<double>> alpha_nothing;  // by model: the value of doing nothing for ever from each lump
  PlanAssumptions assumptions;
  // The components and the solution point into `models`: a plan is moved, never copied.
  CityPlan() = default;
  CityPlan(CityPlan&&) = default;
  CityPlan& operator=(CityPlan&&) = default;
  CityPlan(const CityPlan&) = delete;
  CityPlan& operator=(const CityPlan&) = delete;
};

// The city's plan from the configuration s (the last census, the demand held), for the buildings marked, with their
// beliefs now. The models of buildings with equal lumped generators (to 10⁻¹²) are shared.
inline CityPlan plan_city(const Model& M, const State& s, const std::vector<char>& planned, const std::vector<td::Belief>& beliefs, double budget,
                          const PlanAssumptions& A = {}, const td::CoupledOptions& opt = {.epsilon = 0.05, .max_trials = 200}) {
  CityPlan C;
  C.budget = budget, C.assumptions = A;
  const std::size_t n = s.use.size();
  std::map<std::vector<long long>, std::size_t> seen;
  std::vector<std::size_t> of(n, std::numeric_limits<std::size_t>::max());
  for (std::size_t i = 0; i < n; i++) {
    if (!planned[i]) continue;
    const auto Q0 = lumped_generator(M, s, i, false), Q1 = lumped_generator(M, s, i, true);
    std::vector<long long> key;
    for (const auto* Q : {&Q0, &Q1})
      for (const auto& r : *Q)
        for (double v : r) key.push_back(std::llround(v * 1e12));
    const auto it = seen.find(key);
    if (it != seen.end()) of[i] = it->second;
    else {
      of[i] = C.models.size();
      seen[key] = of[i];
      C.models.push_back(building_model(Q0, Q1, A));
    }
  }
  C.model_of.assign(n, C.models.size());
  const std::vector<double> cost(A.cost.begin(), A.cost.end());
  for (std::size_t i = 0; i < n; i++) {
    if (!planned[i]) continue;
    C.model_of[i] = of[i];
    C.planned.push_back(i);
    C.components.push_back({&C.models[of[i]], cost, beliefs[i], 1.0});
  }
  C.solution = td::lagrangian_bound(C.components, budget, opt);
  // Doing nothing for ever, exactly: α = R(·, nothing) + γ T α, by iteration to 10⁻¹².
  std::vector<std::vector<double>> alpha(C.models.size());
  for (std::size_t m = 0; m < C.models.size(); m++) {
    const auto& P = C.models[m];
    std::vector<double> v(lumps, 0.0), next(lumps);
    for (int it = 0; it < 100000; it++) {
      double change = 0;
      for (std::size_t a = 0; a < lumps; a++) {
        double e = 0;
        for (std::size_t b = 0; b < lumps; b++) e += P.transition(nothing, a, b) * v[b];
        next[a] = P.reward(a, nothing) + P.discount * e;
        change = std::max(change, std::fabs(next[a] - v[a]));
      }
      v.swap(next);
      if (change < 1e-12) break;
    }
    alpha[m] = v;
  }
  for (std::size_t k = 0; k < C.components.size(); k++)
    for (std::size_t a = 0; a < lumps; a++) C.nothing_value += C.components[k].belief[a] * alpha[of[C.planned[k]]][a];
  C.alpha_nothing = std::move(alpha);
  return C;
}

// ---- policies ----------------------------------------------------------------------------------------------------

enum class Rule { nothing, every_vacancy, longest_vacant, by_neighbours, pomdp };

inline const char* rule_name(Rule r) {
  switch (r) {
    case Rule::nothing: return "nothing";
    case Rule::every_vacancy: return "every vacancy";
    case Rule::longest_vacant: return "the longest vacant";
    case Rule::by_neighbours: return "by neighbours";
    case Rule::pomdp: return "the POMDP's ranked policy";
  }
  return "";
}

constexpr std::uint32_t order_stream = 0x6f636372u;   // "occr": a rule's random order
constexpr std::uint32_t census_stream = 0x63656e73u;  // "cens": who a census surveys
constexpr std::uint32_t start_stream = 0x73747274u;   // "strt": a decoupled episode's first states
constexpr std::uint32_t step_stream = 0x73746570u;    // "step": a decoupled episode's draws

// One action for each planned building (in the order of C.planned) for the year.
inline std::vector<std::size_t> decide(Rule rule, const CityPlan& C, const Knowledge& K, const Space& X, int year, std::uint64_t seed) {
  const std::size_t m = C.planned.size();
  std::vector<std::size_t> act(m, nothing);
  if (rule == Rule::nothing) return act;
  if (rule == Rule::pomdp) {
    auto comps = C.components;
    for (std::size_t k = 0; k < m; k++) comps[k].belief = K.belief[C.planned[k]];
    return td::coupled_actions(C.solution, comps, C.budget);
  }
  struct Candidate {
    double priority, key;
    std::size_t k;
  };
  std::vector<Candidate> cand;
  for (std::size_t k = 0; k < m; k++) {
    const std::size_t i = C.planned[k];
    if (K.lump[i] != static_cast<int>(vacant)) continue;
    calx::random::Stream rng(seed, order_stream, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(year));
    double priority = 0;
    if (rule == Rule::longest_vacant) priority = year - K.since[i];
    else if (rule == Rule::by_neighbours) {
      const auto& nb = X.neighbours[i];
      double v = 0;
      for (const std::size_t j : nb) v += K.lump[j] == static_cast<int>(vacant) || K.lump[j] == static_cast<int>(works);
      priority = nb.empty() ? 0 : v / static_cast<double>(nb.size());
    }
    cand.push_back({priority, rng.uniform(), k});
  }
  std::sort(cand.begin(), cand.end(), [](const Candidate& a, const Candidate& b) { return a.priority != b.priority ? a.priority > b.priority : a.key < b.key; });
  double spent = 0;
  const double c = C.assumptions.cost[incentive];
  for (const auto& x : cand) {
    if (spent + c > C.budget + 1e-12) break;
    act[x.k] = incentive, spent += c;
  }
  return act;
}

// The planner learns a year's census: the lump read where the census or an inspection surveyed (o < 3).
inline void learn(Knowledge& K, const CityPlan& C, std::size_t k, std::size_t a, std::size_t o, int year) {
  const std::size_t i = C.planned[k];
  K.belief[i] = td::update(C.models[C.model_of[i]], K.belief[i], a, o).belief;
  if (o == not_surveyed) return;
  if (o == vacant && K.lump[i] != static_cast<int>(vacant)) K.since[i] = year;
  K.lump[i] = static_cast<int>(o), K.year[i] = year;
}

struct Evaluation {
  double mean = 0, error = 0;     // the discounted value (minus the discounted years in V), and its standard error;
                                  // in the decoupled model, of the policy for `horizon` years and then nothing, exactly
  double years = 0, years_error = 0;  // undiscounted years in V over the horizon
  double spent = 0;               // the mean budget spent a year
  double truncation = 0;          // the agent-based model's: the discounted value beyond the horizon is at most this
  int episodes = 0, horizon = 0;
};

// In the planner's decoupled model: each episode draws every building's lump from its belief now, then runs the lumped
// chains under the policy, the reward the expected years in V (the POMDP's own), observations by the census; at the
// horizon each building's value of doing nothing from its lump then (α_nothing) is added, so that the value is that of
// a feasible plan (the policy for `horizon` years, then nothing) with no truncation: a lower bound on V*, in
// expectation, beside the Lagrangian bound above it.
inline Evaluation evaluate_decoupled(const CityPlan& C, Rule rule, const Space& X, const Knowledge& K0, int start, int horizon, int episodes, std::uint64_t seed) {
  Evaluation E;
  E.episodes = episodes, E.horizon = horizon;
  const std::size_t m = C.planned.size();
  const double gamma = C.assumptions.discount;
  E.truncation = 0;
  double s1 = 0, s2 = 0, y1 = 0, y2 = 0, spent = 0;
  for (int e = 0; e < episodes; e++) {
    Knowledge K = K0;
    std::vector<std::size_t> state(m);
    for (std::size_t k = 0; k < m; k++) {
      calx::random::Stream r(seed + static_cast<std::uint64_t>(e), start_stream, static_cast<std::uint32_t>(C.planned[k]), 0);
      state[k] = calx::random::categorical(r, C.components[k].belief);
    }
    double g = 0, years = 0, disc = 1;
    for (int t = 0; t < horizon; t++) {
      const int year = start + t;
      const auto act = decide(rule, C, K, X, year, seed + static_cast<std::uint64_t>(e));
      for (std::size_t k = 0; k < m; k++) {
        const auto& P = C.models[C.model_of[C.planned[k]]];
        const double r = P.reward(state[k], act[k]);
        g += disc * r, years -= r, spent += C.assumptions.cost[act[k]];
        calx::random::Stream rng(seed + static_cast<std::uint64_t>(e), step_stream, static_cast<std::uint32_t>(C.planned[k]), static_cast<std::uint32_t>(t));
        const std::size_t s2n = calx::random::categorical(rng, std::span<const double>(P.T).subspan((act[k] * lumps + state[k]) * lumps, lumps));
        const std::size_t o = calx::random::categorical(rng, std::span<const double>(P.Z).subspan((act[k] * lumps + s2n) * (lumps + 1), lumps + 1));
        state[k] = s2n;
        learn(K, C, k, act[k], o, year + 1);
      }
      disc *= gamma;
    }
    for (std::size_t k = 0; k < m; k++) g += disc * C.alpha_nothing[C.model_of[C.planned[k]]][state[k]];
    s1 += g, s2 += g * g, y1 += years, y2 += years * years;
  }
  const double n = episodes;
  E.mean = s1 / n, E.error = n > 1 ? std::sqrt(std::max(0.0, (s2 - s1 * s1 / n) / (n - 1) / n)) : 0;
  E.years = y1 / n, E.years_error = n > 1 ? std::sqrt(std::max(0.0, (y2 - y1 * y1 / n) / (n - 1) / n)) : 0;
  E.spent = spent / (n * horizon);
  return E;
}

// In the agent-based model: from the configuration s0 at the census of `start − 1` (the planner's knowledge K0), a
// year without intervention, then `horizon` years under the policy; the years in V measured exactly from each run's
// changes; the census at each year's end surveys by ρ (or by an inspection). Episode e runs on seed + e; the same
// seeds for every policy, so that their differences are paired.
struct AbmRun {
  std::vector<double> vacant_by_year;   // the mean share of the planned buildings in V at each year's start
};

inline Evaluation evaluate_abm(const Model& M, const CityPlan& C, Rule rule, const State& s0, const Knowledge& K0, int start, int horizon, int episodes,
                               std::uint64_t seed, AbmRun* run = nullptr) {
  Evaluation E;
  E.episodes = episodes, E.horizon = horizon;
  const std::size_t m = C.planned.size(), n = s0.use.size();
  const double gamma = C.assumptions.discount;
  E.truncation = static_cast<double>(m) * std::pow(gamma, horizon) / (1 - gamma);
  double s1 = 0, s2 = 0, y1 = 0, y2 = 0, spent = 0;
  if (run) run->vacant_by_year.assign(static_cast<std::size_t>(horizon + 1), 0.0);
  const auto in_v = [&](std::size_t u) { return M.uses->lump[u] == vacant; };
  for (int e = 0; e < episodes; e++) {
    const std::uint64_t se = seed + static_cast<std::uint64_t>(e);
    State s = s0;
    std::fill(s.enrolled.begin(), s.enrolled.end(), 0);
    s.changes.clear();
    nomos::agents::gillespie(M, s, start - 1.0, static_cast<double>(start), year_seed(se, start - 1));
    s.changes.clear();
    Knowledge K = K0;
    double g = 0, years = 0, disc = 1;
    for (int t = 0; t < horizon; t++) {
      const int year = start + t;
      if (run) {
        double v = 0;
        for (const std::size_t i : C.planned) v += in_v(s.use[i]);
        run->vacant_by_year[static_cast<std::size_t>(t)] += v / static_cast<double>(m) / episodes;
      }
      const auto act = decide(rule, C, K, *M.space, year, se);
      std::fill(s.enrolled.begin(), s.enrolled.end(), 0);
      for (std::size_t k = 0; k < m; k++) {
        if (act[k] == incentive) s.enrolled[C.planned[k]] = 1;
        spent += C.assumptions.cost[act[k]];
      }
      // The year, and the years in V measured from its changes.
      std::vector<double> since(n, static_cast<double>(year));
      std::vector<char> v(n, 0);
      for (const std::size_t i : C.planned) v[i] = in_v(s.use[i]);
      s.changes.clear();
      nomos::agents::Log log;
      nomos::agents::gillespie(M, s, static_cast<double>(year), year + 1.0, year_seed(se, year), &log);
      double held = 0;
      std::vector<char> planned_now(n, 0);
      for (const std::size_t i : C.planned) planned_now[i] = 1;
      for (std::size_t q = 0; q < s.changes.size(); q++) {
        const auto& c = s.changes[q];
        if (!planned_now[c.building]) continue;
        const double tq = log.events[q].time;
        if (v[c.building]) held += tq - since[c.building];
        since[c.building] = tq, v[c.building] = in_v(c.to);
      }
      for (const std::size_t i : C.planned)
        if (v[i]) held += year + 1.0 - since[i];
      g -= disc * held, years += held;
      // The census at the year's end.
      for (std::size_t k = 0; k < m; k++) {
        const std::size_t i = C.planned[k];
        calx::random::Stream rng(se, census_stream, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(year));
        const bool surveyed = act[k] == inspect || rng.uniform() < C.assumptions.survey;
        learn(K, C, k, act[k], surveyed ? static_cast<std::size_t>(M.uses->lump[s.use[i]]) : not_surveyed, year + 1);
      }
      disc *= gamma;
    }
    if (run) {
      double v = 0;
      for (const std::size_t i : C.planned) v += in_v(s.use[i]);
      run->vacant_by_year[static_cast<std::size_t>(horizon)] += v / static_cast<double>(m) / episodes;
    }
    s1 += g, s2 += g * g, y1 += years, y2 += years * years;
  }
  const double N = episodes;
  E.mean = s1 / N, E.error = N > 1 ? std::sqrt(std::max(0.0, (s2 - s1 * s1 / N) / (N - 1) / N)) : 0;
  E.years = y1 / N, E.years_error = N > 1 ? std::sqrt(std::max(0.0, (y2 - y1 * y1 / N) / (N - 1) / N)) : 0;
  E.spent = spent / (N * horizon);
  return E;
}

}  // namespace urban::occupancy
