// Perception and decision, drawn (ADR 0074): the engine's answers for the Decisions view. Each kind computes what one
// figure draws, from the World Bank's growth of real GDP (the request's "growth" and "economies" paths) and the
// economics domain's regimes and policy, on tyche's decision and processes and nomos's agents.
//
//   {"op": "decide", "kind": "regimes"}                     the pooled regimes and every economy's reading
//   {"op": "decide", "kind": "policy", "params": {...}}     the policymaker's POMDP solved: bounds, the policy over the
//                                                            simplex of beliefs, along its edges, and an economy's years
//   {"op": "decide", "kind": "online", "params": {...}}     POMCP from an economy's belief, simulations against bounds
//   {"op": "decide", "kind": "fund", "params": {...}}       a fund's support among economies under a year's budget
//   {"op": "decide", "kind": "budget", "params": {...}}     a budget on support: the frontier, the dual, the mixture
//   {"op": "decide", "kind": "coordination", "params": {}}  two economies' policymakers, each seeing its own growth
//   {"op": "decide", "kind": "herding", "params": {...}}    Kirman's herding: a run, and its law against the exact one
//   {"op": "decide", "kind": "sorting", "params": {...}}    Schelling's sorting: frames and segregation by tick

#pragma once

#include <economics/policy/policy.hpp>
#include <economics/regimes/regimes.hpp>
#include <nomos/agents/choice.hpp>
#include <nomos/agents/engine.hpp>
#include <science/json.hpp>
#include <tyche/decision/alpha.hpp>
#include <tyche/decision/constrained.hpp>
#include <tyche/decision/multiagent.hpp>
#include <tyche/decision/online.hpp>
#include <tyche/decision/point_based.hpp>
#include <tyche/decision/structured.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace decide {

using science::Json;
using Arr = Json::Array;
using Obj = Json::Object;
namespace td = tyche::decision;

inline Json num(double x, double unit = 1e-4) { return std::isfinite(x) ? Json(std::round(x / unit) * unit) : Json(); }
inline Json nums(const std::vector<double>& v, double unit = 1e-4) {
  Arr a;
  for (const double x : v) a.push_back(num(x, unit));
  return Json(std::move(a));
}

inline Json load(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error(path + ": cannot be read");
  std::ostringstream s;
  s << in.rdbuf();
  return Json::parse(s.str());
}

inline double param(const Json& request, const char* key, double fallback) {
  const Json* p = request.find("params");
  if (!p || !p->is_object()) return fallback;
  const Json* v = p->find(key);
  return v && v->is_number() ? v->as_number() : fallback;
}
inline std::string param_text(const Json& request, const char* key, const std::string& fallback) {
  const Json* p = request.find("params");
  if (!p || !p->is_object()) return fallback;
  const Json* v = p->find(key);
  return v && v->is_string() ? v->as_string() : fallback;
}

// The World Bank's growth of real GDP by economy, aggregates left out, each series 1961 … its last year.
inline std::vector<economics::GrowthSeries> growth_series(const Json& request) {
  const Json economies = load(request["economies"].as_string());
  std::set<std::string> real;
  for (const auto& e : economies[1].as_array())
    if (e["region"]["id"].as_string() != "NA" && !e["id"].as_string().empty()) real.insert(e["id"].as_string());
  const Json values = load(request["growth"].as_string());
  std::map<std::string, std::map<int, double>> by;
  for (const auto& r : values[1].as_array()) {
    const std::string code = r["countryiso3code"].is_string() ? r["countryiso3code"].as_string() : "";
    if (!real.count(code) || !r["value"].is_number()) continue;
    by[code][std::stoi(r["date"].as_string())] = r["value"].as_number();
  }
  std::vector<economics::GrowthSeries> out;
  for (const auto& [code, years] : by) {
    economics::GrowthSeries g{code, years.begin()->first, {}};
    for (int y = years.begin()->first; y <= years.rbegin()->first; y++) {
      const auto it = years.find(y);
      g.growth.push_back(it == years.end() ? std::nan("") : it->second);
    }
    out.push_back(std::move(g));
  }
  return out;
}

inline Json model_json(const economics::Regimes& R) {
  Arr A;
  for (std::size_t i = 0; i < 3; i++) {
    Arr row;
    for (std::size_t j = 0; j < 3; j++) row.push_back(num(R.model.a(i, j)));
    A.push_back(Json(std::move(row)));
  }
  return Json(Obj{{"standardised", Json(R.standardised)}, {"mean", nums(R.model.mean)}, {"sd", nums(R.model.sd)}, {"transition", Json(std::move(A))},
                  {"duration", nums({R.duration(0), R.duration(1), R.duration(2)})}, {"stationary", nums(R.stationary())},
                  {"log_likelihood", num(R.log_likelihood, 0.01)}, {"iterations", Json(R.iterations)}, {"converged", Json(R.converged)},
                  {"observations", Json(static_cast<double>(R.observations))}, {"economies", Json(static_cast<double>(R.readings.size()))}});
}

inline const economics::RegimeReading* reading_of(const economics::Regimes& R, const std::string& code) {
  for (const auto& r : R.readings)
    if (r.economy == code) return &r;
  return nullptr;
}

// The World Bank's names of its economies (the page renames them by the data library's naming policy).
inline Json economy_names(const Json& request) {
  const Json economies = load(request["economies"].as_string());
  Obj names;
  for (const auto& e : economies[1].as_array())
    if (e["region"]["id"].as_string() != "NA" && !e["id"].as_string().empty()) names.push_back({e["id"].as_string(), Json(e["name"].as_string())});
  return Json(std::move(names));
}

inline Json regimes(const Json& request) {
  const auto R = economics::perceive_regimes(growth_series(request), 3, 10, true, true);
  Arr readings;
  for (const auto& r : R.readings) {
    Arr f, m;
    for (std::size_t t = 0; t < r.growth.size(); t++) f.push_back(nums(r.filtered[t], 1e-3)), m.push_back(nums(r.smoothed[t], 1e-3));
    Arr path;
    for (const auto k : r.path) path.push_back(Json(static_cast<double>(k)));
    readings.push_back(Json(Obj{{"economy", Json(r.economy)}, {"first", Json(r.first_year)}, {"growth", nums(r.growth, 0.01)}, {"centre", num(r.centre, 0.01)}, {"scale", num(r.scale, 0.01)}, {"filtered", Json(std::move(f))},
                                {"smoothed", Json(std::move(m))}, {"path", Json(std::move(path))}}));
  }
  return Json(Obj{{"model", model_json(R)}, {"readings", Json(std::move(readings))}, {"names", economy_names(request)}});
}

inline economics::PolicyAssumptions assumptions(const Json& request) {
  economics::PolicyAssumptions A;
  A.welfare = {param(request, "welfare_recession", A.welfare[0]), param(request, "welfare_normal", A.welfare[1]), param(request, "welfare_boom", A.welfare[2])};
  A.cost = {0.0, param(request, "cost_support", A.cost[1]), param(request, "cost_restrain", A.cost[2])};
  A.support_effect = param(request, "support_effect", A.support_effect);
  A.restraint_effect = param(request, "restraint_effect", A.restraint_effect);
  A.discount = param(request, "discount", A.discount);
  return A;
}

inline Json assumptions_json(const economics::PolicyAssumptions& A) {
  return Json(Obj{{"welfare", nums(A.welfare)}, {"cost", nums({A.cost[0], A.cost[1], A.cost[2]})}, {"support_effect", num(A.support_effect)},
                  {"overheating", num(A.overheating)}, {"restraint_effect", num(A.restraint_effect)}, {"deepening", num(A.deepening)},
                  {"edges", nums(A.edges)}, {"discount", num(A.discount)}});
}

inline Json vectors_json(const td::ValueFunction& V) {
  Arr a;
  for (const auto& v : V.vectors) a.push_back(Json(Obj{{"values", nums(v.values)}, {"action", Json(static_cast<double>(v.action))}}));
  return Json(std::move(a));
}

// The belief an economy holds, after a year's growth, for the next: its filtered law that year carried one year by the
// chain (the last year of its record when `year` is 0 or beyond it).
inline std::vector<double> belief_for(const economics::Regimes& R, const std::string& code, int year = 0) {
  const auto* r = reading_of(R, code);
  if (!r) throw std::invalid_argument("no growth record for " + code);
  const int last = r->first_year + static_cast<int>(r->growth.size()) - 1;
  const int y = year <= 0 || year > last ? last : std::max(year, r->first_year);
  return economics::belief_ahead(R.model, r->filtered[static_cast<std::size_t>(y - r->first_year)]);
}

inline Json policy(const Json& request) {
  const auto R = economics::perceive_regimes(growth_series(request), 3, 10, true, true);
  const auto A = assumptions(request);
  const auto P = economics::policy_pomdp(R.model, A);
  const auto H = td::hsvi(P, {.epsilon = param(request, "epsilon", 0.02), .max_trials = 600});
  const auto Q = td::qmdp(P), F = td::fast_informed_bound(P), L = td::blind_lower_bound(P);
  // The policy over the simplex: a triangular grid of beliefs (n + 1)(n + 2)/2.
  const int n = static_cast<int>(param(request, "grid", 30));
  Arr grid;
  for (int i = 0; i <= n; i++)
    for (int j = 0; i + j <= n; j++) {
      const std::vector<double> b{static_cast<double>(i) / n, static_cast<double>(j) / n, static_cast<double>(n - i - j) / n};
      grid.push_back(Json(Arr{num(b[0]), num(b[1]), num(b[2]), Json(static_cast<double>(H.lower.action(b))), num(H.lower.value(b)), num(H.upper.value(b))}));
    }
  // Along each edge of the simplex: (recession, normal), (normal, boom), (recession, boom).
  Arr edges;
  const std::size_t pairs[3][2] = {{0, 1}, {1, 2}, {0, 2}};
  for (const auto& e : pairs) {
    Arr pts;
    for (int k = 0; k <= 100; k++) {
      std::vector<double> b(3, 0.0);
      const double p = k / 100.0;
      b[e[0]] = 1 - p, b[e[1]] = p;
      pts.push_back(Json(Arr{num(p), num(H.lower.value(b)), num(H.upper.value(b)), num(F.value(b)), num(Q.value(b)), num(L.value(b)), Json(static_cast<double>(H.lower.action(b)))}));
    }
    edges.push_back(Json(Obj{{"from", Json(static_cast<double>(e[0]))}, {"to", Json(static_cast<double>(e[1]))}, {"points", Json(std::move(pts))}}));
  }
  Arr trace;
  const std::size_t step = std::max<std::size_t>(1, H.trace.size() / 300);
  for (std::size_t t = 0; t < H.trace.size(); t += step) trace.push_back(Json(Arr{Json(static_cast<double>(t)), num(H.trace[t].first), num(H.trace[t].second)}));
  trace.push_back(Json(Arr{Json(static_cast<double>(H.trials)), num(H.lower_value), num(H.upper_value)}));
  // An economy's years: its filtered belief and the action the policy takes on it.
  const std::string code = param_text(request, "economy", "USA");
  Json years = Json();
  if (const auto* r = reading_of(R, code)) {
    Arr ys;
    for (std::size_t t = 0; t < r->growth.size(); t++) {
      const auto b = r->filtered[t];
      ys.push_back(Json(Arr{Json(r->first_year + static_cast<int>(t)), num(r->growth[t], 0.01), nums(b, 1e-3), Json(static_cast<double>(H.lower.action(b)))}));
    }
    years = Json(Obj{{"economy", Json(code)}, {"years", Json(std::move(ys))}});
  }
  Arr T;
  for (std::size_t a = 0; a < P.actions; a++) {
    Arr rows;
    for (std::size_t s = 0; s < 3; s++) {
      Arr row;
      for (std::size_t s2 = 0; s2 < 3; s2++) row.push_back(num(P.transition(a, s, s2)));
      rows.push_back(Json(std::move(row)));
    }
    T.push_back(Json(std::move(rows)));
  }
  Arr Zs;
  for (const auto& z : economics::growth_bins(R.model, A.edges)) Zs.push_back(nums(z));
  return Json(Obj{{"model", model_json(R)}, {"assumptions", assumptions_json(A)}, {"transition", Json(std::move(T))}, {"bins", Json(std::move(Zs))},
                  {"hsvi", Json(Obj{{"lower", num(H.lower_value)}, {"upper", num(H.upper_value)}, {"trials", Json(H.trials)}, {"converged", Json(H.converged)},
                                    {"vectors", vectors_json(H.lower)}, {"points", Json(static_cast<double>(H.upper.points.size()))}, {"trace", Json(std::move(trace))}})},
                  {"bounds", Json(Obj{{"qmdp", vectors_json(Q)}, {"fib", vectors_json(F)}, {"blind", vectors_json(L)}})},
                  {"grid", Json(Obj{{"n", Json(n)}, {"points", Json(std::move(grid))}})}, {"edges", Json(std::move(edges))}, {"economy", std::move(years)}});
}

inline Json online(const Json& request) {
  const auto R = economics::perceive_regimes(growth_series(request), 3, 10, true, true);
  const auto A = assumptions(request);
  auto P = economics::policy_pomdp(R.model, A);
  const std::string code = param_text(request, "economy", "USA");
  const int year = static_cast<int>(param(request, "year", 2009));
  P.initial = belief_for(R, code, year);
  const auto H = td::hsvi(P, {.epsilon = 0.02, .max_trials = 300});
  const auto L = td::blind_lower_bound(P);
  // Q at the belief from HSVI's bounds, by one step of lookahead.
  std::vector<double> qlo(P.actions), qhi(P.actions);
  for (std::size_t a = 0; a < P.actions; a++) {
    double lo = td::expected_reward(P, P.initial, a), hi = lo;
    const auto po = td::observation_law(P, P.initial, a);
    for (std::size_t o = 0; o < P.observations; o++)
      if (po[o] > 0) {
        const auto nb = td::update(P, P.initial, a, o).belief;
        lo += P.discount * po[o] * H.lower.value(nb), hi += P.discount * po[o] * H.upper.value(nb);
      }
    qlo[a] = lo, qhi[a] = hi;
  }
  const td::TabularModel<td::Pomdp> M{P};
  Arr ladder;
  for (const int sims : {30, 100, 300, 1000, 3000, 10000, 30000}) {
    td::Pomcp<td::TabularModel<td::Pomdp>> planner(M, {.simulations = sims, .reward_span = P.reward_span(), .particles = 2000, .seed = 1000u + static_cast<std::uint64_t>(sims)});
    planner.set_leaf([&](const std::size_t& s) { double v = -1e300; for (const auto& b : L.vectors) v = std::max(v, b.values[s]); return v; });
    const std::size_t a = planner.plan();
    ladder.push_back(Json(Obj{{"simulations", Json(sims)}, {"q", nums(planner.values())}, {"action", Json(static_cast<double>(a))}}));
  }
  return Json(Obj{{"economy", Json(code)}, {"year", Json(year)}, {"belief", nums(P.initial)}, {"q_lower", nums(qlo)}, {"q_upper", nums(qhi)},
                  {"hsvi_action", Json(static_cast<double>(H.lower.action(P.initial)))}, {"ladder", Json(std::move(ladder))}});
}

inline Json fund(const Json& request) {
  const auto R = economics::perceive_regimes(growth_series(request), 3, 10, true, true);
  const auto A = assumptions(request);
  const auto P = economics::policy_pomdp(R.model, A, 2);
  const double budget = param(request, "budget", 20);
  const int year = static_cast<int>(param(request, "year", 2020));
  const std::vector<double> cost{0, 1};
  struct Candidate { std::string code; double now, growth, z; std::vector<double> ahead; };
  std::vector<Candidate> pool;
  for (const auto& r : R.readings) {
    const int last = r.first_year + static_cast<int>(r.growth.size()) - 1;
    if (year < r.first_year || year > last) continue;  // economies with no record that year
    const auto t = static_cast<std::size_t>(year - r.first_year);
    // A belief is an estimate: kept to six places, so that economies equally sure of their regime tie exactly.
    auto ahead = economics::belief_ahead(R.model, r.filtered[t]);
    double sum = 0;
    for (auto& v : ahead) sum += (v = std::round(v * 1e6) / 1e6);
    for (auto& v : ahead) v /= sum;
    // The fall in the economy's own deviations, before the model's clip at ±6.
    pool.push_back({r.economy, r.filtered[t][0], r.growth[t], (r.growth[t] - r.centre) / r.scale, std::move(ahead)});
  }
  // The model cannot tell apart economies of equal belief, and the allocation keeps the order of equal gains:
  // among them the budget goes to the deepest falls in each economy's own deviations, a declared rule.
  const auto depth = [](double z) { return std::isfinite(z) ? z : 1e300; };
  std::stable_sort(pool.begin(), pool.end(), [&](const Candidate& a, const Candidate& b) { return depth(a.z) < depth(b.z); });
  std::vector<td::Component> C;
  int sure = 0;
  for (const auto& c : pool) C.push_back({&P, cost, c.ahead, 1.0}), sure += c.now >= 0.995;
  const auto sol = td::lagrangian_bound(C, budget, {.epsilon = 0.02, .max_trials = 150});
  const auto act = td::coupled_actions(sol, C, budget);
  // The bound against λ, by the fast informed bound alone (each value a valid bound).
  Arr curve;
  const double reach = std::max(0.5, 4 * sol.lambda);
  for (int k = 0; k <= 40; k++) {
    const double lambda = reach * k / 40;
    const auto F = td::fast_informed_bound(td::detail::relaxed(P, cost, lambda));
    double total = lambda * budget / (1 - P.discount);
    for (const auto& c : C) total += F.value(c.belief);
    curve.push_back(Json(Arr{num(lambda), num(total, 0.01)}));
  }
  // Waiting everywhere, exactly: a plan within any budget.
  const auto wait = td::blind_lower_bound(P);
  double waiting = 0;
  for (const auto& c : C) waiting += td::dot(wait.vectors[0].values, c.belief);
  Arr economies;
  int supported = 0;
  for (std::size_t i = 0; i < C.size(); i++) {
    supported += act[i] == 1;
    economies.push_back(Json(Arr{Json(pool[i].code), num(pool[i].now, 1e-3), num(C[i].belief[0], 1e-3), Json(static_cast<double>(act[i])),
                                 num(pool[i].growth, 0.01), num(pool[i].z, 0.01)}));
  }
  return Json(Obj{{"budget", num(budget)}, {"year", Json(year)}, {"lambda", num(sol.lambda)}, {"upper", num(sol.upper, 0.01)}, {"waiting", num(waiting, 0.01)}, {"supported", Json(supported)}, {"sure", Json(sure)},
                  {"curve", Json(std::move(curve))}, {"economies", Json(std::move(economies))}, {"evaluations", Json(sol.evaluations)}});
}

inline Json budget(const Json& request) {
  const auto R = economics::perceive_regimes(growth_series(request), 3, 10, true, true);
  const auto A = assumptions(request);
  auto P = economics::policy_pomdp(R.model, A, 2);
  const std::string code = param_text(request, "economy", "USA");
  const int year = static_cast<int>(param(request, "year", 2009));
  P.initial = belief_for(R, code, year);
  const auto Cst = economics::support_cost(P);
  const double c = param(request, "budget", 1.0);
  calx::random::Stream rng(20261104, 0x62756467u, 0, 0);  // "budg"
  auto B = td::random_beliefs(P, 150, 20, rng);
  B.push_back(P.initial);
  const auto solve = [&](double lambda) { return td::scalarised_infinite(P, Cst, lambda, B); };
  Arr frontier;
  for (int k = -1; k <= 34; k++) {
    const double lambda = k < 0 ? 0.0 : 0.001 * std::pow(1.3, k);
    const auto s = solve(lambda);
    frontier.push_back(Json(Arr{num(lambda), num(s.reward), num(s.cost)}));
  }
  const auto sol = td::lagrangian_dual(solve, c, 12.0, 30);
  auto scal = P;
  for (std::size_t a = 0; a < P.actions; a++)
    for (std::size_t s = 0; s < P.states; s++) scal.reward(s, a) -= sol.lambda * Cst[a * P.states + s];
  scal.initial = P.initial;
  const double dual = td::hsvi(scal, {.epsilon = 0.05, .max_trials = 300}).upper_value + sol.lambda * c;
  const auto free = td::hsvi(P, {.epsilon = 0.05, .max_trials = 300});
  return Json(Obj{{"economy", Json(code)}, {"year", Json(year)}, {"belief", nums(P.initial)}, {"budget", num(c)}, {"frontier", Json(std::move(frontier))},
                  {"lambda", num(sol.lambda)}, {"weight", num(sol.weight)}, {"reward", num(sol.reward)}, {"cost", num(sol.cost)},
                  {"low", Json(Arr{num(sol.low.reward), num(sol.low.cost)})}, {"high", Json(Arr{num(sol.high.reward), num(sol.high.cost)})},
                  {"dual", num(dual)}, {"free", Json(Arr{num(free.lower_value), num(free.upper_value)})}, {"feasible", Json(sol.feasible)}});
}

// Two neighbouring economies, each a policymaker seeing only its own growth (below or above its normal), sharing one welfare,
// support by both shortening both recessions more than either alone (a declared spillover). States: the pair of
// regimes, recession (0) or not (1); actions: wait or support each; horizon 3.
inline Json coordination(const Json& request) {
  const auto R = economics::perceive_regimes(growth_series(request), 3, 10, true, true);
  const double spill = param(request, "spillover", 0.25);
  // Two regimes from the three: recession, and the rest weighted by the stationary law.
  const auto pi = R.stationary();
  const double stay = R.model.a(0, 0), rec_from_rest = (pi[1] * R.model.a(1, 0) + pi[2] * R.model.a(2, 0)) / (pi[1] + pi[2]);
  const auto Zfull = economics::growth_bins(R.model, {0.0});  // below or above the economy's own normal
  const double low_in_rec = Zfull[0][0], low_in_rest = (pi[1] * Zfull[1][0] + pi[2] * Zfull[2][0]) / (pi[1] + pi[2]);
  td::DecPomdp D;
  D.states = 4, D.actions = {2, 2}, D.observations = {2, 2}, D.discount = 1;
  const std::size_t S = 4, JA = 4, JO = 4;
  D.T.assign(JA * S * S, 0.0), D.Z.assign(JA * S * JO, 0.0), D.R.assign(JA * S, 0.0);
  const double welfare[2] = {-3.0, 1.0}, cost = 0.4;
  for (std::size_t a0 = 0; a0 < 2; a0++)
    for (std::size_t a1 = 0; a1 < 2; a1++) {
      const std::size_t ja = a0 + 2 * a1;
      for (std::size_t s = 0; s < S; s++) {
        const std::size_t r0 = s % 2, r1 = s / 2;  // 0 recession, 1 not
        D.R[ja * S + s] = welfare[r0] + welfare[r1] - cost * (a0 + a1);
        // Each economy's chance of being in recession next year.
        const auto next_rec = [&](std::size_t r, std::size_t mine, std::size_t other) {
          if (r == 1) return rec_from_rest;
          double p = stay * (1 - 0.35 * mine);
          if (mine && other) p *= 1 - spill;
          return p;
        };
        const double p0 = next_rec(r0, a0, a1), p1 = next_rec(r1, a1, a0);
        for (std::size_t s2 = 0; s2 < S; s2++) {
          const std::size_t q0 = s2 % 2, q1 = s2 / 2;
          D.T[(ja * S + s) * S + s2] = (q0 == 0 ? p0 : 1 - p0) * (q1 == 0 ? p1 : 1 - p1);
          for (std::size_t o0 = 0; o0 < 2; o0++)
            for (std::size_t o1 = 0; o1 < 2; o1++) {
              const double l0 = q0 == 0 ? low_in_rec : low_in_rest, l1 = q1 == 0 ? low_in_rec : low_in_rest;
              D.Z[(ja * S + s2) * JO + o0 + 2 * o1] = (o0 == 0 ? l0 : 1 - l0) * (o1 == 0 ? l1 : 1 - l1);
            }
        }
      }
    }
  const double p = param(request, "shock", 0.6);  // each economy's chance of contraction at the start (a shared shock)
  D.initial = {p * p, (1 - p) * p, p * (1 - p), (1 - p) * (1 - p)};
  const int h = 3;
  td::JointPolicy best;
  for (const auto& t : td::all_trees(2, 2, h)) {
    const auto br = td::best_response(D, {t, t}, 1);
    if (br.value > best.value) best = br;
  }
  const auto J = td::jesp(D, h);
  // Each acting alone, the other always waiting: the best tree against waiting.
  td::PolicyTree waiting{2, h, std::vector<std::size_t>(td::PolicyTree::nodes(2, h), 0)};
  const auto alone = td::best_response(D, {waiting, waiting}, 0);
  const auto tree_json = [](const td::PolicyTree& q) {
    Arr a;
    for (const auto x : q.action) a.push_back(Json(static_cast<double>(x)));
    return Json(std::move(a));
  };
  return Json(Obj{{"horizon", Json(h)}, {"shock", num(p)}, {"spillover", num(spill)}, {"recession_stay", num(stay)}, {"recession_from_rest", num(rec_from_rest)},
                  {"low_in_recession", num(low_in_rec)}, {"low_otherwise", num(low_in_rest)},
                  {"optimum", Json(Obj{{"value", num(best.value)}, {"trees", Json(Arr{tree_json(best.trees[0]), tree_json(best.trees[1])})}})},
                  {"jesp", Json(Obj{{"value", num(J.value)}, {"trees", Json(Arr{tree_json(J.trees[0]), tree_json(J.trees[1])})}})},
                  {"alone", num(alone.value)}, {"never", num(td::value(D, {waiting, waiting}))}});
}

inline Json herding(const Json& request) {
  const std::size_t N = static_cast<std::size_t>(param(request, "agents", 40));
  const nomos::agents::Kirman M{param(request, "epsilon", 0.004), param(request, "delta", 0.3)};
  const int T = static_cast<int>(param(request, "time", 4000));
  // The law is read from a run fifty times as long as the one drawn: the market switches seldom.
  const int long_run = 50 * T;
  std::vector<int> x(N, 0);
  for (std::size_t i = 0; i < N / 2; i++) x[i] = 1;
  std::vector<double> series, seen(N + 1, 0.0);
  for (int t = 0; t < long_run; t++) {
    nomos::agents::gillespie(M, x, t, t + 1, 20261104 + static_cast<std::uint64_t>(t));
    const int k = nomos::agents::Kirman::ones(x);
    if (t < T) series.push_back(static_cast<double>(k) / static_cast<double>(N));
    seen[static_cast<std::size_t>(k)] += 1.0 / long_run;
  }
  return Json(Obj{{"agents", Json(static_cast<double>(N))}, {"epsilon", num(M.epsilon)}, {"delta", num(M.delta)}, {"threshold", num(M.delta / (static_cast<double>(N) - 1))},
                  {"long_run", Json(static_cast<double>(long_run))}, {"series", nums(series, 1e-3)}, {"seen", nums(seen, 1e-5)}, {"exact", nums(M.stationary(N), 1e-5)}});
}

inline Json sorting(const Json& request) {
  const std::size_t side = static_cast<std::size_t>(param(request, "side", 40));
  const nomos::agents::Schelling M{{side, side, true, nomos::agents::Neighbourhood::moore}, param(request, "threshold", 1.0 / 3)};
  const std::size_t n = static_cast<std::size_t>(param(request, "agents", 0.8 * side * side));
  auto x = M.place(n, 7);
  const std::set<int> keep{0, 1, 2, 4, 8, 30};
  Arr frames, segregation, content;
  for (std::uint32_t t = 0; t <= 30; t++) {
    if (keep.count(static_cast<int>(t))) {
      Arr cells;
      for (const int c : x.cell) cells.push_back(Json(c));
      frames.push_back(Json(Obj{{"tick", Json(static_cast<double>(t))}, {"cells", Json(std::move(cells))}}));
    }
    segregation.push_back(num(M.segregation(x))), content.push_back(num(M.share_content(x)));
    if (t < 30) nomos::agents::tick(M, x, nomos::agents::Schedule::random_sequential, 7, t);
  }
  return Json(Obj{{"side", Json(static_cast<double>(side))}, {"agents", Json(static_cast<double>(n))}, {"threshold", num(M.threshold)},
                  {"frames", Json(std::move(frames))}, {"segregation", Json(std::move(segregation))}, {"content", Json(std::move(content))}});
}

inline Json run(const Json& request) {
  const std::string kind = request["kind"].as_string();
  if (kind == "regimes") return regimes(request);
  if (kind == "policy") return policy(request);
  if (kind == "online") return online(request);
  if (kind == "fund") return fund(request);
  if (kind == "budget") return budget(request);
  if (kind == "coordination") return coordination(request);
  if (kind == "herding") return herding(request);
  if (kind == "sorting") return sorting(request);
  throw std::invalid_argument("no decision \"" + kind + "\"");
}

}  // namespace decide
