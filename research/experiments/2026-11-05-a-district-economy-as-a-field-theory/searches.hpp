// The searches compared at equal budgets (README §5.3): study (1)'s black-box searches of heuresis, with its settings
// (../2026-11-05-a-district-economy-without-a-model/without_a_model.cpp: GA, GA from the plan as built, separable CMA-ES
// with restarts, random search), and two more (simulated annealing from the plan as built, the categorical UMDA), all
// on the shared problem's mixed encoding and kernels, seeing the value alone; and the heat bath on the plan field, which uses the action's
// structure (its locality in space and time, and the people's linear response) to value every state of a frontage
// exactly at once.
//
// The common currency is the evaluation: one call of the black box's objective, the field's passers for a plan and the
// shared problem's evaluate(). The heat bath's local work is converted at the measured cost of one such call; its full
// rebuilds (after each sweep, and for each link tried) each count as one. The flows' cache, which every search shares,
// is warmed before any run is timed.

#pragma once

#include <calx/scheduler/thread_pool.hpp>
#include <heuresis/annealing/annealing.hpp>
#include <heuresis/distributions/distributions.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>

#include "plan.hpp"

namespace field {

namespace hs = heuresis::search;

struct Run {
  std::string algorithm;
  std::uint64_t seed = 0;
  double value = 0, evaluations = 0, seconds = 0;
  std::vector<std::array<double, 2>> curve;  // (evaluations or their equivalents, the best value so far)
  Plan plan;
};

template <class G>
Run finish(const std::string& name, std::uint64_t seed, const hs::Outcome<G>& out, const std::function<Plan(const G&)>& plan_of, double seconds) {
  Run r;
  r.algorithm = name, r.seed = seed, r.seconds = seconds, r.evaluations = out.evaluations;
  if (out.best) r.plan = plan_of(out.best->genotype), r.value = -out.best->evaluation.value();
  for (const auto& p : out.progress) r.curve.push_back({p.evaluations, -p.value});
  return r;
}

// Study (1)'s GA: μ = λ = 100, two elites, crossover 0.9, tournaments of two; from random plans, or from the plan as
// built mutated five times.
inline Run run_ga(const Model& M, const economy::Kernels& K, std::uint64_t seed, double budget, calx::scheduler::ThreadPool& pool, bool from_built) {
  const auto T0 = std::chrono::steady_clock::now();
  const Plan built = M.problem->as_built();
  heuresis::genetic::Operators<Plan> ops{[&](calx::random::Stream& s) {
                                           if (!from_built) return K.sample(s);
                                           Plan p = built;
                                           for (int k = 0; k < 5; k++) p = K.mutate(p, s);
                                           return p;
                                         },
                                         [&](const Plan& x, calx::random::Stream& s) { return K.mutate(x, s); },
                                         [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return K.crossover(a, b, s); }};
  heuresis::genetic::GaOptions o;
  o.mu = 100, o.lambda = 100, o.elites = 2, o.crossover_rate = 0.9, o.tournament = 2;
  heuresis::genetic::GeneticAlgorithm<Plan> ga(ops, o, hs::Draws(seed));
  const auto out = hs::run(ga, [&](const Plan& p) { return -M.evaluate(p).value; }, hs::Budget{budget}, hs::Pool(pool));
  return finish<Plan>(from_built ? "GA from as built" : "GA", seed, out, [](const Plan& p) { return p; }, std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count());
}

// Study (1)'s CMA-ES: separable, on the canonical relaxation [0, 1]^621 (bounds by transformation), σ₀ = 0.3, each
// restart from a new random mean 0.1 + 0.8 u until the budget is spent.
// With pool null the evaluations run on the calling thread (keyed draws: the run is the same either way).
inline Run run_cma(const Model& M, std::uint64_t seed, double budget, calx::scheduler::ThreadPool* pool) {
  const auto T0 = std::chrono::steady_clock::now();
  const Problem& P = *M.problem;
  const std::size_t D = P.dimension();
  Run r;
  r.algorithm = "CMA-ES", r.seed = seed;
  double spent = 0, best = -1e300;
  std::uint64_t serial = 0, restart = 0;
  while (spent < budget) {
    heuresis::strategies::CmaOptions o;
    o.separable = true, o.bounds = heuresis::variation::Box::cube(D, 0, 1);
    auto ms = hs::Draws(seed)(hs::tag("study.cma.mean"), restart, 0);
    heuresis::variation::Real mean(D);
    for (auto& x : mean) x = 0.1 + 0.8 * ms.uniform();
    heuresis::strategies::Cma cma(mean, 0.3, o, hs::Draws(seed * 1000 + restart), serial);
    const auto f = [&](const heuresis::variation::Real& u) { return -M.evaluate(P.decode(u)).value; };
    const auto out = pool ? hs::run(cma, f, hs::Budget{budget - spent}, hs::Pool(*pool)) : hs::run(cma, f, hs::Budget{budget - spent});
    for (const auto& p : out.progress)
      if (-p.value > best) best = -p.value, r.curve.push_back({spent + p.evaluations, best});
    if (out.best && -out.best->evaluation.value() >= best - 1e-9) r.plan = P.decode(out.best->genotype), best = std::max(best, -out.best->evaluation.value());
    spent += out.evaluations, serial = cma.next_serial(), restart++;
    if (out.evaluations == 0) break;
  }
  r.value = best, r.evaluations = spent;
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
  return r;
}

inline Run run_sa(const Model& M, const economy::Kernels& K, std::uint64_t seed, double budget, double T0) {
  const auto t0 = std::chrono::steady_clock::now();
  const Plan start = M.problem->as_built();
  heuresis::annealing::SimulatedAnnealing<Plan> sa([start](calx::random::Stream&) { return start; }, [&](const Plan& x, calx::random::Stream& s) { return K.mutate(x, s); },
                                                    heuresis::annealing::Schedule::modified_lam(T0, static_cast<std::uint64_t>(budget)), hs::Draws(seed));
  const auto out = hs::run(sa, [&](const Plan& p) { return -M.evaluate(p).value; }, hs::Budget{budget});
  return finish<Plan>("SA", seed, out, [](const Plan& p) { return p; }, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}

// The encoding as a product of categorical variables: kinds (6), hour genes (24 or 25 values), links (2).
inline heuresis::distributions::ProductCategorical categorical_family(const Problem& P) {
  heuresis::distributions::ProductCategorical F;
  for (std::size_t f = 0; f < P.size(); f++) F.sizes.push_back(economy::kinds);
  for (std::size_t i = 0; i < economy::hour_genes * P.size(); i++) F.sizes.push_back(i % 2 == 0 ? 24 : 25);
  for (std::size_t j = 0; j < P.links.size(); j++) F.sizes.push_back(2);
  return F;
}
inline Plan plan_of_categories(const Problem& P, const std::vector<std::uint32_t>& x) {
  Plan p;
  const std::size_t n = P.size(), H = economy::hour_genes * n;
  for (std::size_t f = 0; f < n; f++) p.kind.push_back(static_cast<std::int64_t>(x[f]));
  for (std::size_t i = 0; i < H; i++) p.hours.push_back(static_cast<std::int64_t>(x[n + i]));
  for (std::size_t j = 0; j < P.links.size(); j++) p.links.push_back(static_cast<std::uint8_t>(x[n + H + j]));
  return p;
}

inline Run run_umda(const Model& M, std::uint64_t seed, double budget, calx::scheduler::ThreadPool& pool) {
  namespace hd = heuresis::distributions;
  const auto T0 = std::chrono::steady_clock::now();
  const auto F = categorical_family(*M.problem);
  hd::DistributionSearch<hd::ProductCategorical, hd::CategoricalUmda> umda(F, F.uniform(), 100, hd::CategoricalUmda{30, true}, hs::Draws(seed));
  const Problem& P = *M.problem;
  const auto out = hs::run(umda, [&](const std::vector<std::uint32_t>& x) { return -M.evaluate(plan_of_categories(P, x)).value; }, hs::Budget{budget}, hs::Pool(pool));
  return finish<std::vector<std::uint32_t>>("UMDA", seed, out, [&](const std::vector<std::uint32_t>& x) { return plan_of_categories(P, x); },
                                            std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count());
}

inline Run run_random(const Model& M, const economy::Kernels& K, std::uint64_t seed, double budget, calx::scheduler::ThreadPool& pool) {
  const auto T0 = std::chrono::steady_clock::now();
  hs::RandomSearch<Plan> rs([&](calx::random::Stream& s) { return K.sample(s); }, hs::Draws(seed), 32);
  const auto out = hs::run(rs, [&](const Plan& p) { return -M.evaluate(p).value; }, hs::Budget{budget}, hs::Pool(pool));
  return finish<Plan>("random", seed, out, [](const Plan& p) { return p; }, std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count());
}

// The heat bath on the plan field: sweeps over the frontages in a random order, each drawing its state from its exact
// conditional law at T (README §5.2), then each link from its own; T falls geometrically from T0 to T1 over `sweeps`,
// then T = 0 (iterated conditional modes) until a sweep changes nothing.
struct HeatBathOptions {
  std::size_t sweeps = 40, zero_sweeps = 20;
  double T0 = 3000, T1 = 3;
  bool links = true;
};

inline Run run_heat_bath(const Model& M, Plan start, std::uint64_t seed, const HeatBathOptions& o, double evaluation_seconds, const std::string& name) {
  const auto t0 = std::chrono::steady_clock::now();
  Run r;
  r.algorithm = name, r.seed = seed;
  PlanState S(M, std::move(start));
  LocalWork W;
  double full = 1;
  calx::random::Stream z(seed, 0x6865, 0x6174, 0);
  double best = S.value;
  Plan best_plan = S.plan;
  r.curve.push_back({full, best});
  const std::size_t n = M.problem->size(), L = M.problem->links.size();
  for (std::size_t sweep = 0; sweep < o.sweeps + o.zero_sweeps; sweep++) {
    const bool hot = sweep < o.sweeps;
    const double T = hot ? o.T0 * std::pow(o.T1 / o.T0, o.sweeps > 1 ? static_cast<double>(sweep) / static_cast<double>(o.sweeps - 1) : 0.0) : 0.0;
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    for (std::size_t i = n; i > 1; i--) std::swap(order[i - 1], order[z.below(i)]);
    double moved = 0;
    for (auto f : order) moved += std::fabs(S.update(f, T, z, W));
    S.rebuild(), full++;
    if (o.links) {
      std::vector<std::size_t> lo(L);
      std::iota(lo.begin(), lo.end(), std::size_t{0});
      for (std::size_t i = L; i > 1; i--) std::swap(lo[i - 1], lo[z.below(i)]);
      for (auto j : lo) {
        PlanState Q = S;
        Q.plan.links[j] ^= 1;
        Q.rebuild(), full++;
        const double d = Q.value - S.value;
        const bool take = T > 0 ? z.uniform() < 1 / (1 + std::exp(-d / T)) : d > 1e-6;
        if (take) S = std::move(Q), moved += std::fabs(d);
      }
    }
    if (S.value > best) best = S.value, best_plan = S.plan;
    r.curve.push_back({full + W.seconds / evaluation_seconds, best});
    if (!hot && moved < 1e-6) break;
  }
  r.value = best, r.plan = best_plan, r.evaluations = full + W.seconds / evaluation_seconds;
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return r;
}

}  // namespace field
