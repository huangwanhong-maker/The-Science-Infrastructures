// A district's economy without a model of its flow: study (1) of the shared problem (README.md; the problem in
// ../2026-11-05-a-district-economy-the-problem). The people come from the data alone (flows.hpp); the searches treat
// the value as a black box.
//
//   without_a_model <problem-inputs.json> <problem.json> <sensors.json> <results directory> [--threads 4] [--scale 1] [--seeds 5] [--profile]
//
// --scale multiplies the large budget (20,000 evaluations at 1); --profile times 200 random plans and stops.
//
// Runs, every one on the same value, at equal budgets over seeds 1…5:
//   large budget (20,000): the genetic algorithm on the mixed encoding (GA), separable CMA-ES on the canonical
//     relaxation, MAP-Elites over the spending's evening and weekend shares, NSGA-II on value, spread and mix, and
//     random search;
//   small budget (400): the same five, and Bayesian optimisation on a random ten-dimensional embedding of the
//     relaxation (REMBO), which cannot afford the large one.
// Every evaluation of every run is kept (its value, spread, mix, access, evening and weekend shares, by ticket), so
// that each algorithm is read by the same measures: its best value, the front of all it evaluated, and the
// repertoire of all it evaluated.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <calx/scheduler/thread_pool.hpp>
#include <heuresis/diversity/diversity.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/models/models.hpp>
#include <heuresis/multiobjective/multiobjective.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>

#include "economy.hpp"
#include "flows.hpp"

namespace {

using economy::Json;
using economy::numbers;
using economy::Plan;
namespace hs = heuresis::search;
namespace hv = heuresis::variation;
namespace mo = heuresis::multiobjective;
namespace dv = heuresis::diversity;
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}
void write_file(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  f << text;
  if (!f) throw std::runtime_error("cannot write " + path);
}
double r(double x, double d) { return std::isfinite(x) ? std::round(x * d) / d : std::nan(""); }

// What one evaluation measured.
struct Row {
  double value = 0, spread = 0, mix = 0, access = 0, evening = 0, weekend = 0;
};

// Every evaluation of a run, by ticket, and its best plan.
struct Tape {
  std::mutex m;
  std::map<hs::Ticket, Row> rows;
  double best = -std::numeric_limits<double>::infinity();
  Plan best_plan;
  hs::Ticket best_ticket = 0;
  void add(hs::Ticket t, const Row& row, const Plan& p) {
    std::lock_guard<std::mutex> lock(m);
    rows[t] = row;
    if (row.value > best || (row.value == best && t < best_ticket)) best = row.value, best_plan = p, best_ticket = t;
  }
};

struct Run {
  Run(std::string a, std::size_t b, std::uint64_t s) : algorithm(std::move(a)), budget(b), seed(s) {}
  std::string algorithm;
  std::size_t budget = 0;
  std::uint64_t seed = 0;
  double seconds = 0, evaluations = 0;
  std::vector<Row> rows;  // in ticket order
  double best = 0;
  Plan best_plan;
};

// The grid of the repertoire: evening share × weekend share of the spending.
constexpr double grid_hi = 0.6;
constexpr std::size_t grid_n = 12;
std::size_t cell_of(double evening, double weekend) {
  const auto c = [](double v) { return std::min<std::size_t>(grid_n - 1, static_cast<std::size_t>(std::max(0.0, v) / grid_hi * grid_n)); };
  return c(weekend) * grid_n + c(evening);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: without_a_model <problem-inputs.json> <problem.json> <sensors.json> <results directory> [--threads N] [--scale s] [--seeds k] [--profile]\n";
    return 2;
  }
  std::size_t threads = 4, seeds_n = 5;
  double scale = 1;
  bool profile = false;
  for (int i = 5; i < argc; i++) {
    if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) threads = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--scale") == 0 && i + 1 < argc) scale = std::stod(argv[++i]);
    else if (std::strcmp(argv[i], "--seeds") == 0 && i + 1 < argc) seeds_n = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--profile") == 0) profile = true;
  }
  const std::string out = argv[4];
  const auto T0 = Clock::now();
  const Json in = Json::parse(read_file(argv[1]));
  const Json pj = Json::parse(read_file(argv[2]));
  const Json sj = Json::parse(read_file(argv[3]));
  const economy::Problem P(in, &pj);
  const std::size_t n = P.size(), L = P.links.size(), D = P.dimension();

  // ---- the people, from the data alone ---------------------------------------------------------------------------
  const auto counters = without::counters_of(sj);
  const auto KF = without::fit_kriging(counters);
  const auto& K = KF.best;
  std::printf("kriging: %zu counters; l = %.0f m, nugget %.2f, s2 = %.3f; leave-one-out RMSE %.3f in ln(1 + count), R2 %.3f, median weekly ratio %.3f\n", counters.size(), K.length, K.nugget, K.scale,
              KF.loo_rmse, KF.loo_r2, KF.loo_weekly_median_ratio);
  const auto E = without::establishments_of(in);
  const auto hours = without::typical_hours(P);
  std::vector<without::Elasticity> elasticities;
  for (double radius : {25.0, 50.0, 100.0}) elasticities.push_back(without::fit_elasticity(counters, E, hours, radius));
  const auto& rule = elasticities[1];
  for (const auto& e : elasticities)
    std::printf("elasticity within %.0f m: %.3f (se %.3f, t %.2f; %zu rows, %zu counters)\n", e.radius, e.slope, e.se, e.t, e.rows, e.clusters);
  const without::Flows F(P, K, E, rule);
  const auto measure = [&](const Plan& p) { return P.evaluate(p, F(p)); };
  const Plan built = P.as_built();
  const auto O0 = measure(built);
  std::printf("as built: V = %.0f A$/week, spending %.0f, footfall %.0f people-passings/week, spread %.3f (%.1f h), mix %.3f, access %.1f m, evening %.3f, weekend %.3f\n", O0.value, O0.spending, O0.footfall,
              O0.spread, O0.effective_hours, O0.mix, O0.access, O0.evening_share, O0.weekend_share);
  std::fflush(stdout);

  calx::scheduler::ThreadPool pool(threads);
  const economy::Kernels ker{&P};
  {
    const auto t = Clock::now();
    const hs::Draws draws(99);
    std::vector<Plan> ps;
    for (std::size_t i = 0; i < 200; i++) {
      auto s = draws(hs::tag("study.timing"), 0, i);
      ps.push_back(ker.sample(s));
    }
    std::vector<double> v(ps.size());
    calx::scheduler::parallel_for(pool, ps.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t i = b; i < e; i++) v[i] = measure(ps[i]).value;
    });
    std::sort(v.begin(), v.end());
    std::printf("200 random plans: %.2f ms each on %zu threads; value median %.0f, best %.0f A$/week\n", 1e3 * seconds_since(t) / 200, threads, v[100], v.back());
    // Again, with the distances of their links cached: the cost of an evaluation once a search has met its links.
    const auto t2 = Clock::now();
    calx::scheduler::parallel_for(pool, ps.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t i = b; i < e; i++) v[i] = measure(ps[i]).value;
    });
    std::printf("the same plans again, their distances cached: %.2f ms each\n", 1e3 * seconds_since(t2) / 200);
  }
  if (profile) return 0;

  // ---- the searches ---------------------------------------------------------------------------------------------------
  const std::size_t large = static_cast<std::size_t>(std::max(1000.0, std::round(20000 * scale))), small = 400;
  std::vector<std::uint64_t> seeds;
  for (std::size_t s = 1; s <= seeds_n; s++) seeds.push_back(s);
  std::vector<Run> runs;
  // The objective, with the tape: f(genotype, ticket) for each representation.
  const auto row_of = [](const economy::Outcome& o) { return Row{o.value, o.spread, o.mix, o.access, o.evening_share, o.weekend_share}; };
  const auto finish = [&](Run& run, Tape& tape) {
    for (const auto& [t, row] : tape.rows) run.rows.push_back(row);
    run.evaluations = static_cast<double>(run.rows.size());
    run.best = tape.best, run.best_plan = tape.best_plan;
    std::printf("  %-10s budget %5zu seed %llu: best %.0f A$/week, %zu evaluations, %.0f s\n", run.algorithm.c_str(), run.budget, static_cast<unsigned long long>(run.seed), run.best, run.rows.size(), run.seconds);
    std::fflush(stdout);
  };
  const Plan* most_spread_plan = nullptr;
  std::vector<Plan> nsga_front_plans;  // seed 1, large budget: its final population's non-dominated plans
  std::vector<std::pair<std::size_t, hs::Individual<Plan>>> me_cells;  // MAP-Elites seed 1, large budget
  (void)most_spread_plan;

  for (const std::size_t budget : {small, large}) {
    const bool is_small = budget == small;
    for (auto seed : seeds) {
      // GA on the mixed encoding.
      {
        Tape tape;
        Run run{"GA", budget, seed};
        heuresis::genetic::GaOptions o;
        o.mu = is_small ? 20 : 100, o.lambda = o.mu, o.elites = 2, o.crossover_rate = 0.9, o.tournament = 2;
        heuresis::genetic::GeneticAlgorithm<Plan> ga({[&](calx::random::Stream& s) { return ker.sample(s); }, [&](const Plan& x, calx::random::Stream& s) { return ker.mutate(x, s); },
                                                      [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s); }},
                                                     o, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(ga, [&](const Plan& p, hs::Ticket t) { const auto oc = measure(p); tape.add(t, row_of(oc), p); return hs::scalar(-oc.value); }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool));
        run.seconds = seconds_since(t0);
        finish(run, tape);
        runs.push_back(std::move(run));
      }
      // The same GA started from the plan as built, as a planner would start: its first generation the plan as built
      // mutated five times over (outside the comparison of equal knowledge, §2.3).
      {
        Tape tape;
        Run run{"GA from as built", budget, seed};
        heuresis::genetic::GaOptions o;
        o.mu = is_small ? 20 : 100, o.lambda = o.mu, o.elites = 2, o.crossover_rate = 0.9, o.tournament = 2;
        heuresis::genetic::GeneticAlgorithm<Plan> ga({[&](calx::random::Stream& s) {
                                                        Plan p = built;
                                                        for (int k = 0; k < 5; k++) p = ker.mutate(p, s);
                                                        return p;
                                                      },
                                                      [&](const Plan& x, calx::random::Stream& s) { return ker.mutate(x, s); },
                                                      [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s); }},
                                                     o, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(ga, [&](const Plan& p, hs::Ticket t) { const auto oc = measure(p); tape.add(t, row_of(oc), p); return hs::scalar(-oc.value); }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool));
        run.seconds = seconds_since(t0);
        finish(run, tape);
        runs.push_back(std::move(run));
      }
      // Separable CMA-ES on the canonical relaxation, restarted from a new random mean when it stops.
      {
        Tape tape;
        Run run{"CMA-ES", budget, seed};
        const auto t0 = Clock::now();
        std::uint64_t serial = 0, restart = 0;
        double spent = 0;
        while (spent < static_cast<double>(budget)) {
          heuresis::strategies::CmaOptions o;
          o.separable = true, o.bounds = hv::Box::cube(D, 0, 1);
          auto ms = hs::Draws(seed)(hs::tag("study.cma.mean"), restart, 0);
          hv::Real mean(D);
          for (auto& x : mean) x = 0.1 + 0.8 * ms.uniform();
          heuresis::strategies::Cma cma(mean, 0.3, o, hs::Draws(seed * 1000 + restart), serial);
          const auto oc = hs::run(cma, [&](const hv::Real& u, hs::Ticket t) { const Plan p = P.decode(u); const auto e = measure(p); tape.add(t, row_of(e), p); return hs::scalar(-e.value); },
                                  hs::Budget{static_cast<double>(budget) - spent}, hs::Pool(pool));
          spent += oc.evaluations, serial = cma.next_serial(), restart++;
          if (oc.evaluations == 0) break;
        }
        run.seconds = seconds_since(t0);
        finish(run, tape);
        runs.push_back(std::move(run));
      }
      // MAP-Elites over the evening and weekend shares of the spending.
      {
        Tape tape;
        Run run{"MAP-Elites", budget, seed};
        dv::MapElites<Plan, dv::Grid>::Options o;
        o.initial = is_small ? 100 : 500, o.variation_probability = 0.5, o.breadth = is_small ? 8 : 32;
        dv::MapElites<Plan, dv::Grid> me({[&](calx::random::Stream& s) { return ker.sample(s); }, [&](const Plan& x, calx::random::Stream& s) { return ker.mutate(x, s); },
                                          [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s).first; }},
                                         dv::Grid({0.0, 0.0}, {grid_hi, grid_hi}, {grid_n, grid_n}), o, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(me, [&](const Plan& p, hs::Ticket t) {
                  const auto oc = measure(p);
                  tape.add(t, row_of(oc), p);
                  return hs::evaluation({-oc.value}, {}, {oc.weekend_share, oc.evening_share});
                }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool), o.breadth);
        run.seconds = seconds_since(t0);
        if (!is_small && seed == 1)
          for (auto c : me.archive().filled()) me_cells.push_back({c, *me.archive().elite(c)});
        finish(run, tape);
        runs.push_back(std::move(run));
      }
      // NSGA-II on value, spread and mix.
      {
        Tape tape;
        Run run{"NSGA-II", budget, seed};
        const std::size_t mu = is_small ? 20 : 100;
        mo::Nsga2<Plan> alg({[&](calx::random::Stream& s) { return ker.sample(s); }, [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s); },
                             [&](const Plan& x, calx::random::Stream& s) { return ker.mutate(x, s); }},
                            mu, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(alg, [&](const Plan& p, hs::Ticket t) { const auto oc = measure(p); tape.add(t, row_of(oc), p); return hs::evaluation({-oc.value, -oc.spread, -oc.mix}); },
                hs::Budget{static_cast<double>(budget)}, hs::Pool(pool));
        run.seconds = seconds_since(t0);
        if (!is_small && seed == 1) {
          std::vector<mo::Point> pts;
          for (const auto& x : alg.population()) pts.push_back(x.evaluation.objectives);
          for (auto i : mo::non_dominated(pts)) nsga_front_plans.push_back(alg.population()[i].genotype);
        }
        finish(run, tape);
        runs.push_back(std::move(run));
      }
      // Random search.
      {
        Tape tape;
        Run run{"random", budget, seed};
        hs::RandomSearch<Plan> rs([&](calx::random::Stream& s) { return ker.sample(s); }, hs::Draws(seed), 32);
        const auto t0 = Clock::now();
        hs::run(rs, [&](const Plan& p, hs::Ticket t) { const auto oc = measure(p); tape.add(t, row_of(oc), p); return hs::scalar(-oc.value); }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool), 32);
        run.seconds = seconds_since(t0);
        finish(run, tape);
        runs.push_back(std::move(run));
      }
      // Bayesian optimisation on a random embedding (REMBO, Wang et al. 2016): z ∈ [−√d, √d]^d, u = (clip(A z, −1, 1) + 1)/2.
      if (is_small) {
        Tape tape;
        Run run{"BO (REMBO)", budget, seed};
        const std::size_t d = 10;
        std::vector<double> A(D * d);
        auto as = hs::Draws(seed)(hs::tag("study.rembo.embedding"), 0, 0);
        for (auto& x : A) x = as.normal();
        const auto embed = [&](const hv::Real& z) {
          std::vector<double> u(D);
          for (std::size_t i = 0; i < D; i++) {
            double v = 0;
            for (std::size_t j = 0; j < d; j++) v += A[i * d + j] * z[j];
            u[i] = (std::clamp(v, -1.0, 1.0) + 1) / 2;
          }
          return u;
        };
        heuresis::models::BayesianOptions o;
        o.initial = 2 * d + 2, o.breadth = 4;
        heuresis::models::BayesianOptimisation bo(hv::Box::cube(d, -std::sqrt(static_cast<double>(d)), std::sqrt(static_cast<double>(d))), o, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(bo, [&](const hv::Real& z, hs::Ticket t) { const Plan p = P.decode(embed(z)); const auto e = measure(p); tape.add(t, row_of(e), p); return hs::scalar(-e.value); },
                hs::Budget{static_cast<double>(budget)}, hs::Pool(pool), 4);
        run.seconds = seconds_since(t0);
        finish(run, tape);
        runs.push_back(std::move(run));
      }
    }
  }

  // ---- reading the runs -------------------------------------------------------------------------------------------
  // The best value so far, along each run.
  std::ostringstream any;
  any << "algorithm,budget,seed,evaluations,best_value\n";
  std::ostringstream rcsv;
  rcsv << "algorithm,budget,seed,evaluations,seconds,best_value,best_spread,best_mix,best_access\n";
  for (const auto& run : runs) {
    double best = -std::numeric_limits<double>::infinity();
    const std::size_t step = run.budget == small ? 10 : 100;
    for (std::size_t i = 0; i < run.rows.size(); i++) {
      best = std::max(best, run.rows[i].value);
      if ((i + 1) % step == 0 || i + 1 == run.rows.size()) any << run.algorithm << "," << run.budget << "," << run.seed << "," << i + 1 << "," << r(best, 10) << "\n";
    }
    const auto ob = measure(run.best_plan);
    rcsv << run.algorithm << "," << run.budget << "," << run.seed << "," << run.evaluations << "," << r(run.seconds, 100) << "," << r(run.best, 10) << "," << r(ob.spread, 1e5) << "," << r(ob.mix, 1e5) << "," << r(ob.access, 100) << "\n";
  }
  write_file(out + "/anytime.csv", any.str());
  write_file(out + "/runs.csv", rcsv.str());

  // Fronts of all that each run evaluated, in value, spread and mix; hypervolumes on the union's ideal and nadir.
  std::map<const Run*, std::vector<mo::Point>> fronts;
  std::vector<double> lo(3, 1e300), hi(3, -1e300);
  for (const auto& run : runs) {
    if (run.budget != large) continue;
    std::vector<mo::Point> pts;
    for (const auto& row : run.rows) pts.push_back({-row.value, -row.spread, -row.mix});
    std::vector<mo::Point> f;
    for (auto i : mo::non_dominated(pts)) f.push_back(pts[i]);
    for (const auto& p : f)
      for (std::size_t k = 0; k < 3; k++) lo[k] = std::min(lo[k], p[k]), hi[k] = std::max(hi[k], p[k]);
    fronts[&run] = std::move(f);
  }
  // The value's axis: from the union's best to the plan as built, so that plans worse than as built fall outside.
  hi[0] = std::max(lo[0] + 1, -O0.value);
  std::ostringstream hvc, fcsv;
  hvc << "algorithm,seed,front,hypervolume\n";
  fcsv << "algorithm,seed,value,spread,mix\n";
  const std::vector<double> ref(3, 1.1);
  for (const auto& [run, f] : fronts) {
    std::vector<mo::Point> q;
    for (const auto& p : f) {
      mo::Point x(3);
      for (std::size_t k = 0; k < 3; k++) x[k] = (p[k] - lo[k]) / (hi[k] - lo[k]);
      if (x[0] <= 1.1) q.push_back(x);
      fcsv << run->algorithm << "," << run->seed << "," << r(-p[0], 10) << "," << r(-p[1], 1e6) << "," << r(-p[2], 1e6) << "\n";
    }
    hvc << run->algorithm << "," << run->seed << "," << f.size() << "," << (q.empty() ? 0.0 : mo::hypervolume(q, ref)) << "\n";
  }
  write_file(out + "/hypervolume.csv", hvc.str());
  write_file(out + "/fronts.csv", fcsv.str());

  // Repertoires of all that each run evaluated, on MAP-Elites' grid: coverage, and the QD-score Σ max(0, best value) (A$k/week).
  std::ostringstream qd;
  qd << "algorithm,budget,seed,coverage,positive_cells,qd_score_kaud\n";
  for (const auto& run : runs) {
    std::vector<double> best(grid_n * grid_n, -std::numeric_limits<double>::infinity());
    for (const auto& row : run.rows) {
      const auto c = cell_of(row.evening, row.weekend);
      best[c] = std::max(best[c], row.value);
    }
    double filled = 0, positive = 0, score = 0;
    for (double b : best)
      if (std::isfinite(b)) filled += 1, positive += b > 0, score += std::max(0.0, b) / 1000;
    qd << run.algorithm << "," << run.budget << "," << run.seed << "," << r(filled / (grid_n * grid_n), 1e4) << "," << positive << "," << r(score, 10) << "\n";
  }
  write_file(out + "/repertoire-scores.csv", qd.str());
  {
    Json cells = Json::array();
    for (const auto& [c, x] : me_cells) {
      const auto o = measure(x.genotype);
      cells.push(Json::object({{"cell", static_cast<double>(c)}, {"weekend_bin", static_cast<double>(c / grid_n)}, {"evening_bin", static_cast<double>(c % grid_n)},
                               {"value", r(o.value, 10)}, {"evening_share", r(o.evening_share, 1e5)}, {"weekend_share", r(o.weekend_share, 1e5)}, {"spread", r(o.spread, 1e5)}, {"mix", r(o.mix, 1e5)},
                               {"frontages_by_kind", numbers(o.frontages_by_kind, 0)}, {"spending_by_hour", numbers(o.spending_by_hour, 1)}}));
    }
    write_file(out + "/repertoire.json", Json::object({{"grid", Json::object({{"descriptors", Json::array({"evening share of spending (18-24 h)", "weekend share of spending"})}, {"low", 0.0}, {"high", grid_hi}, {"resolution", static_cast<double>(grid_n)}})},
                                                       {"cells", cells}})
                                              .dump());
  }

  // ---- the plans --------------------------------------------------------------------------------------------------
  const Run* best_run = nullptr;
  for (const auto& run : runs)
    if (!best_run || run.best > best_run->best) best_run = &run;
  const Plan best = best_run->best_plan;
  // From NSGA-II's final front (seed 1, the large budget): the most spread and the most mixed of the plans worth at
  // least the plan as built; and the front itself.
  Plan spread_plan = best, mixed_plan = best;
  {
    double bs = -1, bm = -1;
    std::ostringstream nf;
    nf << "value,spread,mix,access,evening,weekend\n";
    for (const auto& p : nsga_front_plans) {
      const auto o = measure(p);
      nf << r(o.value, 10) << "," << r(o.spread, 1e6) << "," << r(o.mix, 1e6) << "," << r(o.access, 100) << "," << r(o.evening_share, 1e6) << "," << r(o.weekend_share, 1e6) << "\n";
      if (o.value >= O0.value && o.spread > bs) bs = o.spread, spread_plan = p;
      if (o.value >= O0.value && o.mix > bm) bm = o.mix, mixed_plan = p;
    }
    write_file(out + "/nsga-front.csv", nf.str());
  }
  const auto plan_entry = [&](const std::string& key, const std::string& title, const std::string& from, const Plan& p) {
    const auto flows = F(p);
    const auto o = P.evaluate(p, flows);
    Json fr = Json::array();
    for (std::size_t f = 0; f < n; f++) {
      double passers = 0;
      for (std::size_t h = 0; h < economy::week; h++) passers += flows[f * economy::week + h];
      fr.push(Json::object({{"kind", economy::kind_name(p.kind[f])}, {"hours", Json::array({static_cast<double>(p.hours[4 * f]), static_cast<double>(p.hours[4 * f + 1]), static_cast<double>(p.hours[4 * f + 2]), static_cast<double>(p.hours[4 * f + 3])})},
                            {"open_hours_week", p.weekly_hours(f)}, {"passers_week", r(passers, 1)}, {"spending_week", r(o.spending_by_frontage[f], 1)}, {"value_week", r(o.value_by_frontage[f], 1)}}));
    }
    return Json::object({{"key", key}, {"title", title}, {"from", from}, {"outcome", economy::outcome_json(o, true)}, {"plan", economy::plan_json(p)}, {"frontages", fr}});
  };
  Json plans = Json::array();
  plans.push(plan_entry("as-built", "the district as built", "CLUE 2024 kinds, OpenStreetMap's median hours, no link", built));
  plans.push(plan_entry("best", "the best value found", best_run->algorithm + ", budget " + std::to_string(best_run->budget) + ", seed " + std::to_string(best_run->seed), best));
  plans.push(plan_entry("spread", "the most spread of NSGA-II's final front worth at least the plan as built", "NSGA-II, budget 20,000, seed 1", spread_plan));
  plans.push(plan_entry("mixed", "the most mixed of NSGA-II's final front worth at least the plan as built", "NSGA-II, budget 20,000, seed 1", mixed_plan));
  // The best plan of each algorithm at the large budget, seed 1.
  for (const auto& run : runs)
    if (run.budget == large && run.seed == 1) plans.push(plan_entry("best-" + run.algorithm, "the best of " + run.algorithm, run.algorithm + ", budget 20,000, seed 1", run.best_plan));
  write_file(out + "/plans.json", Json::object({{"plans", plans}}).dump());

  // ---- where the model-free plan is blind ---------------------------------------------------------------------------
  Json blind = Json::object();
  {
    // 1. Links: what the value says of each (only its cost: the flows do not move), and what the network says: the
    //    share of the people's shortest walks between frontages (pairs weighted by their passers as built) that the
    //    link would carry, and how many frontages' own segment would lose or gain such walks by a tenth or more.
    std::vector<double> mass(n, 0.0);
    for (std::size_t f = 0; f < n; f++)
      for (std::size_t h = 0; h < economy::week; h++) mass[f] += F.base[f * economy::week + h];
    const auto routes = [&](const hv::Bits& made, std::vector<double>* own, double* through_links) {
      auto Dd = P.with_links(made);
      study::WalkGraph G;
      std::map<std::size_t, std::size_t> node;
      std::map<std::uint64_t, std::size_t> edge;
      std::vector<char> is_link;
      for (std::size_t i = 0; i < Dd.complex.size(); i++)
        if (Dd.complex.cells[i].dim == 0) node[i] = G.nodes++;
      for (std::size_t i = 0; i < Dd.complex.size(); i++) {
        const auto& c = Dd.complex.cells[i];
        if (c.dim != 1) continue;
        const auto ends = Dd.complex.ends(i);
        edge[c.id] = G.edges.size();
        G.edges.push_back({node.at(ends[0]), node.at(ends[1]), Dd.length(c.id)});
        const auto nm = Dd.name.count(c.id) ? Dd.name.at(c.id) : std::string();
        is_link.push_back(nm == "a passage" || (nm.size() > 16 && nm.substr(nm.size() - 16) == " carried through"));
      }
      G.finish();
      std::vector<study::OnSegment> at;
      for (const auto& f : P.frontages) {
        double bd = 1e300, s = 0;
        std::uint64_t seg = 0;
        for (const auto& [id, l] : Dd.line) {
          if (is_link[edge.at(id)]) continue;
          const auto [along, dist] = physis::computational_geometry::nearest_on_polyline(l, f.anchor);
          if (dist < bd) bd = dist, seg = id, s = along;
        }
        at.push_back({edge.at(seg), s});
      }
      std::vector<std::vector<double>> dist(n);
      for (std::size_t f = 0; f < n; f++) dist[f] = study::from_sources(G, {at[f]});
      std::vector<double> load(G.edges.size(), 0.0);
      double total = 0, linked = 0;
      for (std::size_t f = 0; f < n; f++)
        for (std::size_t g = f + 1; g < n; g++) {
          const double w = mass[f] * mass[g];
          const double dfg = study::to_target(G, dist[f], {at[f]}, at[g]);
          total += w;
          bool via = false;
          for (std::size_t e = 0; e < G.edges.size(); e++) {
            const auto& E2 = G.edges[e];
            const double a = std::min(dist[f][E2.u] + E2.length + dist[g][E2.v], dist[f][E2.v] + E2.length + dist[g][E2.u]);
            if (a <= dfg + 1e-6) {
              load[e] += w;
              via = via || is_link[e];
            }
          }
          if (via) linked += w;
        }
      own->assign(n, 0.0);
      for (std::size_t f = 0; f < n; f++) (*own)[f] = load[at[f].edge] / total;
      *through_links = linked / total;
    };
    std::vector<double> own0;
    double through0 = 0;
    routes(hv::Bits(L, 0), &own0, &through0);
    Json links = Json::array();
    for (std::size_t k = 0; k < L; k++) {
      hv::Bits made(L, 0);
      made[k] = 1;
      std::vector<double> own;
      double through = 0;
      routes(made, &own, &through);
      std::size_t changed = 0;
      Json change = Json::array();
      for (std::size_t f = 0; f < n; f++) {
        const double c = own0[f] > 0 ? own[f] / own0[f] - 1 : 0;
        change.push(r(c, 1e4));
        changed += std::fabs(c) >= 0.1;
      }
      Plan with = best;
      with.links.assign(L, 0);
      with.links[k] = 1;
      Plan without = best;
      without.links.assign(L, 0);
      const auto ow = measure(with), on = measure(without);
      links.push(Json::object({{"link", static_cast<double>(k)}, {"kind", P.links[k].kind}, {"joins", Json::array({P.links[k].from, P.links[k].to})}, {"length_m", r(P.links[k].length, 10)},
                               {"share_of_walks_through_it", r(through, 1e5)}, {"frontages_changed_by_a_tenth", static_cast<double>(changed)}, {"own_segment_change", change},
                               {"value_change_aud_week", r(ow.value - on.value, 10)}, {"access_change_m", r(ow.access - on.access, 100)}}));
    }
    Json own0j = Json::array();
    for (double v : own0) own0j.push(r(v, 1e5));
    blind.set("links", links);
    blind.set("own_segment_share_without_links", own0j);
    // 2. Conservation: the people the rule adds or takes away.
    const auto ob = measure(best);
    blind.set("footfall", Json::object({{"as_built", r(O0.footfall, 1)}, {"best", r(ob.footfall, 1)}, {"change", r(ob.footfall / O0.footfall - 1, 1e5)}}));
    // 3. Lanes: the frontages on a laneway, by kriging, against the laneway counters.
    Json lanes = Json::array();
    for (std::size_t f = 0; f < n; f++)
      if (P.frontages[f].way != "street") {
        double w = 0;
        for (std::size_t h = 0; h < economy::week; h++) w += F.base[f * economy::week + h];
        lanes.push(Json::object({{"frontage", static_cast<double>(f)}, {"street", P.frontages[f].street}, {"way", P.frontages[f].way}, {"kriged_week", r(w, 1)}}));
      }
    Json lane_counters = Json::array();
    for (const auto& c : counters)
      if (c.way == "laneway") {
        double w = 0;
        for (double v : c.week) w += v;
        lane_counters.push(Json::object({{"id", static_cast<double>(c.id)}, {"name", c.name}, {"week", r(w, 1)}}));
      }
    blind.set("lanes", Json::object({{"frontages", lanes}, {"laneway_counters", lane_counters}}));
    // 4. Weeks that are not typical: a public holiday on the Monday, and a wet week (1 mm in each daylight hour),
    //    by the median responses of the counters in the district's box.
    std::vector<double> hol, rain;
    const auto& box = in["district"]["district_box"];
    for (const auto& s : sj["sensors"].as_array()) {
      const double x = s["x"].as_number(), y = s["y"].as_number();
      if (x < box["x0"].as_number() || x > box["x1"].as_number() || y < box["y0"].as_number() || y > box["y1"].as_number()) continue;
      if (s["response"].has("holiday")) hol.push_back(s["response"]["holiday"]["total"].as_number());
      if (s["response"].has("rain")) rain.push_back(s["response"]["rain"]["total"].as_number());
    }
    std::sort(hol.begin(), hol.end()), std::sort(rain.begin(), rain.end());
    const double bh = hol[hol.size() / 2], br = rain[rain.size() / 2];
    Json weeks = Json::array();
    for (const auto* pl : {&built, &best}) {
      auto f0 = F(*pl), fh = f0, fr = f0;
      for (std::size_t f = 0; f < n; f++)
        for (std::size_t h = 0; h < economy::week; h++) {
          if (h < 24) fh[f * economy::week + h] *= std::exp(bh);
          if (economy::hour_of_day(h) >= 7 && economy::hour_of_day(h) < 19) fr[f * economy::week + h] *= std::exp(br * 1.0);
        }
      weeks.push(Json::object({{"plan", pl == &built ? "as built" : "best"}, {"typical", r(P.evaluate(*pl, f0).value, 1)}, {"holiday_monday", r(P.evaluate(*pl, fh).value, 1)}, {"wet_week", r(P.evaluate(*pl, fr).value, 1)}}));
    }
    blind.set("weeks", Json::object({{"holiday_log_response", r(bh, 1e5)}, {"rain_log_response_per_mm", r(br, 1e5)}, {"counters", static_cast<double>(hol.size())}, {"values", weeks}}));
  }
  write_file(out + "/blindness.json", blind.dump());

  // ---- the people, for drawing --------------------------------------------------------------------------------------
  {
    Json grid = Json::array();
    for (const auto& g : KF.grid) grid.push(Json::array({g[0], g[1], r(g[2], 100)}));
    Json loo = Json::array();
    for (std::size_t i = 0; i < counters.size(); i++) {
      double po = 0, ob = 0;
      for (std::size_t h = 0; h < economy::week; h++) po += std::expm1(KF.loo[i][h]), ob += counters[i].week[h];
      loo.push(Json::object({{"id", static_cast<double>(counters[i].id)}, {"name", counters[i].name}, {"observed_week", r(ob, 1)}, {"predicted_week", r(po, 1)}}));
    }
    Json base = Json::array();
    for (std::size_t f = 0; f < n; f++) {
      std::vector<double> w(F.base.begin() + static_cast<std::ptrdiff_t>(f * economy::week), F.base.begin() + static_cast<std::ptrdiff_t>((f + 1) * economy::week));
      base.push(numbers(w, 1));
    }
    // The kriged passers at the middle of each segment of the district, at four hours, keyed by the segment's id.
    Json segs = Json::object();
    const std::size_t hours4[] = {8, 13, 19, 5 * 24 + 14};  // Monday 08:00, 13:00, 19:00; Saturday 14:00
    for (const auto id : P.base.of(1)) {
      const auto& line = P.base.line.at(id);
      const auto mid = physis::computational_geometry::point_at(line, physis::computational_geometry::polyline_length(line) / 2);
      const auto m = K.predict({mid[0], mid[1]});
      Json v = Json::array();
      for (auto h : hours4) v.push(r(std::expm1(m[h]), 1));
      segs.set(std::to_string(id), v);
    }
    Json el = Json::array();
    for (const auto& e : elasticities) el.push(Json::object({{"radius_m", e.radius}, {"slope", r(e.slope, 1e5)}, {"se", r(e.se, 1e5)}, {"t", r(e.t, 100)}, {"r2_within", r(e.r2_within, 1e4)}, {"rows", static_cast<double>(e.rows)}, {"counters", static_cast<double>(e.clusters)}}));
    write_file(out + "/people.json", Json::object({{"kriging", Json::object({{"length_m", K.length}, {"nugget_share", K.nugget}, {"scale", r(K.scale, 1e5)}, {"log_likelihood", r(K.log_likelihood, 100)},
                                                                          {"grid", grid}, {"loo_rmse", r(KF.loo_rmse, 1e4)}, {"loo_r2", r(KF.loo_r2, 1e4)}, {"loo_weekly_median_ratio", r(KF.loo_weekly_median_ratio, 1e4)}, {"loo", loo}})},
                                                   {"elasticity", el}, {"rule_radius_m", rule.radius}, {"passers_as_built", base}, {"segments_hours", Json::array({"Monday 08:00", "Monday 13:00", "Monday 19:00", "Saturday 14:00"})},
                                                   {"segments", segs}})
                                         .dump());
  }

  // ---- summary ----------------------------------------------------------------------------------------------------
  {
    Json rj = Json::array();
    for (const auto& run : runs) rj.push(Json::object({{"algorithm", run.algorithm}, {"budget", static_cast<double>(run.budget)}, {"seed", static_cast<double>(run.seed)}, {"best_value", r(run.best, 10)}, {"seconds", r(run.seconds, 10)}}));
    write_file(out + "/summary.json", Json::object({{"as_built", economy::outcome_json(O0)}, {"best", economy::outcome_json(measure(best))}, {"best_from", best_run->algorithm + " budget " + std::to_string(best_run->budget) + " seed " + std::to_string(best_run->seed)},
                                                    {"budgets", Json::array({static_cast<double>(small), static_cast<double>(large)})}, {"seeds", static_cast<double>(seeds_n)}, {"threads", static_cast<double>(threads)},
                                                    {"normalisation", Json::object({{"low", numbers(lo, 4)}, {"high", numbers(hi, 4)}})}, {"runs", rj}, {"seconds", r(seconds_since(T0), 1)}})
                                         .dump());
  }
  std::printf("done in %.0f s\n", seconds_since(T0));
  return 0;
}
