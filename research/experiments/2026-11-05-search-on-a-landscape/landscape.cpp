// Search on a landscape: how six algorithms of heuresis move over Rosenbrock's, Rastrigin's and Branin's functions
// in two dimensions.
//
// Every algorithm is driven by hand through the ask–tell protocol (search/protocol.hpp): ask a generation, evaluate
// it, tell it, and read the algorithm's state before and after. The run's best and its improvements are kept by
// search::Elite, exactly as search::run keeps them, so that analysis/ reads the runs as it reads any outcome.
//
// Written to the directory given as the first argument:
//   trajectories.json  seed 0 of every algorithm on every function: each generation's state and candidates, 4 s.f.
//   landscapes.json    f − f* of each function on a 141 × 141 grid of its box, 4 s.f.
//   anytime.csv        the best f − f* after b evaluations, quartiles over 31 seeds (analysis::anytime_quantile)
//   spread.csv         the distance of the best point to the nearest minimiser and the search's spread, quartiles
//   runs.csv           one line per run: what it spent, how it ended, its runtimes to three targets

#include <heuresis/analysis/analysis.hpp>
#include <heuresis/benchmarks/benchmarks.hpp>
#include <heuresis/distributions/distributions.hpp>
#include <heuresis/models/models.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>
#include <heuresis/swarms/swarms.hpp>

#include <tyche/descriptive/summary.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <numbers>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace hs = heuresis::search;
namespace hv = heuresis::variation;
namespace hb = heuresis::benchmarks;
namespace ha = heuresis::analysis;
namespace st = heuresis::strategies;
namespace sw = heuresis::swarms;
namespace hd = heuresis::distributions;
namespace hm = heuresis::models;
using Real = hv::Real;

namespace {

constexpr std::size_t seeds = 31;
constexpr double budget = 2000;      // evaluations, for every algorithm but Bayesian optimisation
constexpr double bo_budget = 100;    // Bayesian optimisation's, whose evaluations are meant to be expensive
constexpr double precision = 1e-12;  // a run ends once f − f* ≤ 10⁻¹²
const std::vector<std::size_t> gp_counts = {6, 15, 30, 60, 100};
constexpr std::size_t gp_grid = 50, landscape_grid = 141;

struct Landscape {
  std::string name;
  std::function<double(std::span<const double>)> f;
  hv::Box box;
  double fstar;
  std::vector<Real> minimisers;
  Real start;     // the initial point of the (1+1)-ES, CMA-ES and the cross-entropy method
  double sigma0;  // their initial step size
};

std::vector<Landscape> landscapes() {
  const double pi = std::numbers::pi;
  return {
      {"rosenbrock", [](std::span<const double> x) { return hb::rosenbrock(x); }, {{-2, -1}, {2, 3}}, 0.0, {{1, 1}}, {-1.5, -0.5}, 0.5},
      {"rastrigin", [](std::span<const double> x) { return hb::rastrigin(x); }, hv::Box::cube(2, -5.12, 5.12), 0.0, {{0, 0}}, {3.5, -3.0}, 2.0},
      {"branin", [](std::span<const double> x) { return hb::branin(x); }, {{-5, 0}, {10, 15}}, 0.397887357729738, {{-pi, 12.275}, {pi, 2.275}, {3 * pi, 2.475}}, {2.5, 12.0}, 3.0},
  };
}

const std::vector<std::string> algorithms = {"(1+1)-ES", "CMA-ES", "DE", "PSO", "CE", "BO", "random"};

double nearest(const Landscape& L, std::span<const double> x) {
  double d = INFINITY;
  for (const auto& m : L.minimisers) {
    double s = 0;
    for (std::size_t i = 0; i < x.size(); i++) s += (x[i] - m[i]) * (x[i] - m[i]);
    d = std::min(d, std::sqrt(s));
  }
  return d;
}

// The root mean square distance of points from their centroid.
double rms_spread(const std::vector<Real>& pts) {
  if (pts.empty()) return NAN;
  const std::size_t n = pts[0].size();
  Real c(n, 0.0);
  for (const auto& p : pts)
    for (std::size_t i = 0; i < n; i++) c[i] += p[i] / static_cast<double>(pts.size());
  double s = 0;
  for (const auto& p : pts)
    for (std::size_t i = 0; i < n; i++) s += (p[i] - c[i]) * (p[i] - c[i]);
  return std::sqrt(s / static_cast<double>(pts.size()));
}

// One generation: its state before it was drawn, what it drew, and the search after it was told.
struct Generation {
  double evaluations = 0, best = INFINITY, distance = INFINITY, spread = NAN;
  Real centre;                   // the parent (ES), the mean (CMA-ES, CE) the generation was drawn about
  double sigma = NAN;            // the step size it was drawn with (ES, CMA-ES)
  std::vector<double> cov;       // entries (1,1), (1,2), (2,2): CMA-ES's C (its law is N(m, σ²C)), CE's covariance
  std::vector<Real> points;      // the candidates
  std::vector<double> values;    // their values
  std::vector<Real> velocities;  // PSO: the velocities that brought the particles to the candidates
  std::vector<Real> population;  // DE: the population after selection
  std::vector<double> population_values;
};

struct GpSnapshot {
  std::size_t count = 0;
  std::vector<double> mean;  // the surrogate's mean on a gp_grid × gp_grid grid of the box, row by row in x₂
  std::vector<double> sd;    // its standard deviation
};

struct Trace {
  std::string ending;
  ha::Run run;
  std::vector<Generation> generations;
  std::vector<GpSnapshot> gp;
};

// Asks, evaluates, tells, until the budget or f* + precision; `pre` reads the state a generation is drawn from,
// `post` the state after it is told. The best and its improvements are kept by search::Elite as search::run does.
template <class A, class Pre, class Post>
Trace drive(A& a, const Landscape& L, double limit, Pre&& pre, Post&& post) {
  Trace t;
  hs::Elite<Real> elite;
  const double target = L.fstar + precision;
  double evals = 0;
  try {
    for (;;) {
      if (elite.value() <= target) { t.ending = "target"; break; }
      if (evals >= limit) { t.ending = "budget"; break; }
      if (a.stopped()) { t.ending = "stopped"; break; }
      Generation g;
      pre(a, g);
      auto c = a.ask(static_cast<std::size_t>(std::ceil(limit - evals)));
      if (c.empty()) { t.ending = a.stopped() ? "stopped" : "stalled"; break; }
      std::vector<hs::Told> told(c.size());
      for (std::size_t i = 0; i < c.size(); i++) {
        const double v = L.f(c[i].genotype);
        evals += 1;
        told[i] = {c[i].ticket, hs::scalar(v)};
        elite.offer(hs::Individual<Real>{c[i].ticket, c[i].genotype, told[i].evaluation, c[i].origin}, evals, evals);
        g.points.push_back(c[i].genotype);
        g.values.push_back(v);
      }
      a.tell(std::span<const hs::Told>(told));
      g.evaluations = evals;
      g.best = elite.value();
      g.distance = nearest(L, elite.best()->genotype);
      post(a, g);
      t.generations.push_back(std::move(g));
    }
  } catch (const std::exception& e) {
    t.ending = std::string("error: ") + e.what();
  }
  t.run = {elite.progress(), evals, evals};
  return t;
}

const char* stop_name(st::CmaStop r) {
  switch (r) {
    case st::CmaStop::none: return "none";
    case st::CmaStop::generations: return "generations";
    case st::CmaStop::tol_x: return "TolX";
    case st::CmaStop::tol_fun: return "TolFun";
    case st::CmaStop::condition: return "ConditionCov";
    case st::CmaStop::no_effect_axis: return "NoEffectAxis";
    case st::CmaStop::no_effect_coordinate: return "NoEffectCoord";
    case st::CmaStop::stagnation: return "Stagnation";
    case st::CmaStop::tol_up_sigma: return "TolUpSigma";
  }
  return "?";
}

Trace run_one(const Landscape& L, const std::string& name, std::uint64_t seed) {
  const hs::Draws draws(seed);
  const double rn = std::sqrt(static_cast<double>(L.box.size()));
  auto none = [](auto&, Generation&) {};
  if (name == "(1+1)-ES") {
    st::OnePlusOneEs a(L.start, L.sigma0, st::OnePlusOneOptions{}, draws);
    return drive(
        a, L, budget,
        [&](st::OnePlusOneEs& e, Generation& g) {
          g.centre = e.parent() ? e.parent()->genotype : L.start;
          g.sigma = e.sigma();
        },
        [&](st::OnePlusOneEs& e, Generation& g) { g.spread = e.sigma() * rn; });
  }
  if (name == "CMA-ES") {
    st::Cma a(L.start, L.sigma0, st::CmaOptions{}, draws);
    Trace t = drive(
        a, L, budget,
        [](st::Cma& c, Generation& g) {
          const auto& s = c.state();
          g.centre = s.mean;
          g.sigma = s.sigma;
          g.cov = {s.C(0, 0), s.C(0, 1), s.C(1, 1)};
        },
        [](st::Cma& c, Generation& g) {
          const auto& s = c.state();
          g.spread = s.sigma * std::sqrt(s.C(0, 0) + s.C(1, 1));
        });
    if (t.ending == "stopped") t.ending += std::string(" (") + stop_name(a.stop_reason()) + ")";
    return t;
  }
  if (name == "DE") {
    st::DifferentialEvolution a(L.box, st::DeOptions{}, draws);
    return drive(a, L, budget, none, [](st::DifferentialEvolution& d, Generation& g) {
      for (const auto& x : d.population()) {
        g.population.push_back(x.genotype);
        g.population_values.push_back(x.evaluation.value());
      }
      g.spread = rms_spread(g.population);
    });
  }
  if (name == "PSO") {
    sw::ParticleSwarm a(L.box, sw::PsoOptions{}, draws);
    return drive(a, L, budget, none, [](sw::ParticleSwarm& p, Generation& g) {
      g.velocities = p.velocities();
      g.spread = rms_spread(p.positions());
    });
  }
  if (name == "CE") {
    hd::DistributionSearch a(hd::Gaussian{2}, hd::Gaussian::isotropic(L.start, L.sigma0), 50, hd::CrossEntropy{10, 0.7, 0.8, 5}, draws);
    using A = decltype(a);
    return drive(
        a, L, budget,
        [](A& d, Generation& g) {
          const auto& p = d.parameters();
          g.centre = p.m;
          g.cov = {p.C(0, 0), p.C(0, 1), p.C(1, 1)};
        },
        [](A& d, Generation& g) {
          const auto& p = d.parameters();
          g.spread = std::sqrt(p.C(0, 0) + p.C(1, 1));
        });
  }
  if (name == "BO") {
    hm::BayesianOptimisation a(L.box, hm::BayesianOptions{}, draws);
    std::vector<double> told;  // every value told so far, to undo the surrogate's standardisation
    std::vector<GpSnapshot> snapshots;
    Trace t = drive(a, L, bo_budget, none, [&](hm::BayesianOptimisation& b, Generation& g) {
      told.insert(told.end(), g.values.begin(), g.values.end());
      if (seed != 0 || std::find(gp_counts.begin(), gp_counts.end(), told.size()) == gp_counts.end()) return;
      // The surrogate's standardisation (models/bayesian.hpp): mean 0 and the sample standard deviation 1.
      double mu = 0, s2 = 0;
      for (double v : told) mu += v / static_cast<double>(told.size());
      for (double v : told) s2 += (v - mu) * (v - mu) / static_cast<double>(told.size() - 1);
      const double sd = s2 > 0 ? std::sqrt(s2) : 1;
      const hm::GaussianProcess gp = b.surrogate();
      GpSnapshot snap{told.size(), {}, {}};
      for (std::size_t j = 0; j < gp_grid; j++)
        for (std::size_t i = 0; i < gp_grid; i++) {
          const std::vector<double> u{static_cast<double>(i) / (gp_grid - 1), static_cast<double>(j) / (gp_grid - 1)};
          const auto p = gp.predict(u);
          snap.mean.push_back(mu + sd * p.mean - L.fstar);
          snap.sd.push_back(sd * std::sqrt(std::max(0.0, p.variance)));
        }
      snapshots.push_back(std::move(snap));
    });
    t.gp = std::move(snapshots);
    return t;
  }
  // random search, the baseline
  hs::RandomSearch<Real> a([box = L.box](calx::random::Stream& s) { return hv::uniform_point(box, s); }, draws);
  return drive(a, L, budget, none, none);
}

std::string num(double v) {
  if (!std::isfinite(v)) return "null";
  char b[32];
  std::snprintf(b, sizeof b, "%.4g", v);
  return b;
}
std::string full(double v) {
  if (!std::isfinite(v)) return "";
  char b[40];
  std::snprintf(b, sizeof b, "%.10g", v);
  return b;
}
std::string point(const Real& x) {
  std::string s = "[";
  for (std::size_t i = 0; i < x.size(); i++) s += (i ? "," : "") + num(x[i]);
  return s + "]";
}
std::string points(const std::vector<Real>& xs) {
  std::string s = "[";
  for (std::size_t i = 0; i < xs.size(); i++) s += (i ? "," : "") + point(xs[i]);
  return s + "]";
}
std::string numbers(const std::vector<double>& v) {
  std::string s = "[";
  for (std::size_t i = 0; i < v.size(); i++) s += (i ? "," : "") + num(v[i]);
  return s + "]";
}

std::string heuresis_version() {
#ifdef HEURESIS_COMMIT
  return HEURESIS_COMMIT;
#else
  return "6b3cc88";
#endif
}

// The grid of budgets 1, …, B, ten to a decade, and B itself.
std::vector<double> budgets(double B) {
  std::vector<double> b;
  for (int k = 0;; k++) {
    const double v = std::round(std::pow(10.0, k / 10.0));
    if (v >= B) break;
    if (b.empty() || v > b.back()) b.push_back(v);
  }
  b.push_back(B);
  return b;
}

// Values as gaps f − f*, which keep their meaning at four significant digits when f* is not 0.
std::vector<double> gaps(std::vector<double> v, double fstar) {
  for (auto& x : v) x -= fstar;
  return v;
}

double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  return tyche::descriptive::quantile_sorted(v, q, 7);
}

}  // namespace

int main(int argc, char** argv) {
  const std::filesystem::path out = argc > 1 ? argv[1] : "results";
  std::filesystem::create_directories(out);
  const auto L = landscapes();
  const auto t0 = std::chrono::steady_clock::now();

  // Every (function, algorithm, seed), on threads; each run is a pure function of its seed.
  struct Job {
    std::size_t f, a, seed;
  };
  std::vector<Job> jobs;
  for (std::size_t a = 0; a < algorithms.size(); a++)  // Bayesian optimisation first: its runs are the longest
    for (std::size_t f = 0; f < L.size(); f++)
      for (std::size_t s = 0; s < seeds; s++) jobs.push_back({f, a, s});
  std::stable_sort(jobs.begin(), jobs.end(), [](const Job& x, const Job& y) { return (algorithms[x.a] == "BO") > (algorithms[y.a] == "BO"); });
  std::vector<Trace> traces(L.size() * algorithms.size() * seeds);
  auto slot = [&](std::size_t f, std::size_t a, std::size_t s) -> Trace& { return traces[(f * algorithms.size() + a) * seeds + s]; };
  std::atomic<std::size_t> next{0};
  std::mutex io;
  const unsigned threads = std::max(1u, std::thread::hardware_concurrency());
  std::vector<std::thread> pool;
  for (unsigned w = 0; w < threads; w++)
    pool.emplace_back([&] {
      for (std::size_t k; (k = next++) < jobs.size();) {
        const Job& j = jobs[k];
        slot(j.f, j.a, j.seed) = run_one(L[j.f], algorithms[j.a], j.seed);
        if (j.seed == seeds - 1) {
          std::lock_guard lock(io);
          std::cerr << L[j.f].name << " " << algorithms[j.a] << " done\n";
        }
      }
    });
  for (auto& t : pool) t.join();

  // runs.csv
  {
    std::ofstream o(out / "runs.csv");
    o << "function,algorithm,seed,evaluations,final,ending,rt_1e-2,rt_1e-4,rt_1e-8\n";
    for (std::size_t f = 0; f < L.size(); f++)
      for (std::size_t a = 0; a < algorithms.size(); a++)
        for (std::size_t s = 0; s < seeds; s++) {
          const Trace& t = slot(f, a, s);
          const double final_value = t.run.progress.empty() ? INFINITY : t.run.progress.back().value - L[f].fstar;
          o << L[f].name << ',' << algorithms[a] << ',' << s << ',' << t.run.evaluations << ',' << full(final_value) << ',' << t.ending;
          for (double e : {1e-2, 1e-4, 1e-8}) {
            const auto r = ha::runtime(t.run, L[f].fstar + e);
            o << ',' << (r ? full(*r) : "");
          }
          o << '\n';
        }
  }

  // anytime.csv: quartiles of the best f − f* after b evaluations, over the seeds.
  {
    std::ofstream o(out / "anytime.csv");
    o << "function,algorithm,budget,q25,q50,q75\n";
    for (std::size_t f = 0; f < L.size(); f++)
      for (std::size_t a = 0; a < algorithms.size(); a++) {
        std::vector<ha::Run> runs;
        for (std::size_t s = 0; s < seeds; s++) runs.push_back(slot(f, a, s).run);
        const auto b = budgets(algorithms[a] == "BO" ? bo_budget : budget);
        const auto q25 = ha::anytime_quantile(runs, b, 0.25), q50 = ha::anytime_quantile(runs, b, 0.5), q75 = ha::anytime_quantile(runs, b, 0.75);
        for (std::size_t k = 0; k < b.size(); k++)
          o << L[f].name << ',' << algorithms[a] << ',' << b[k] << ',' << full(q25[k] - L[f].fstar) << ',' << full(q50[k] - L[f].fstar) << ',' << full(q75[k] - L[f].fstar) << '\n';
      }
  }

  // spread.csv: the state after the last generation completed within b evaluations, quartiles over the seeds.
  {
    std::ofstream o(out / "spread.csv");
    o << "function,algorithm,budget,distance_q25,distance_q50,distance_q75,spread_q25,spread_q50,spread_q75\n";
    for (std::size_t f = 0; f < L.size(); f++)
      for (std::size_t a = 0; a < algorithms.size(); a++) {
        if (algorithms[a] == "random") continue;
        for (double b : budgets(algorithms[a] == "BO" ? bo_budget : budget)) {
          std::vector<double> d, sp;
          for (std::size_t s = 0; s < seeds; s++) {
            const auto& gs = slot(f, a, s).generations;
            const Generation* last = nullptr;
            for (const auto& g : gs)
              if (g.evaluations <= b) last = &g;
            if (!last) continue;
            d.push_back(last->distance);
            if (std::isfinite(last->spread)) sp.push_back(last->spread);
          }
          if (d.size() * 2 < seeds) continue;  // fewer than half the runs have completed a generation
          o << L[f].name << ',' << algorithms[a] << ',' << b << ',' << full(quantile(d, 0.25)) << ',' << full(quantile(d, 0.5)) << ',' << full(quantile(d, 0.75));
          if (sp.size() * 2 >= seeds) o << ',' << full(quantile(sp, 0.25)) << ',' << full(quantile(sp, 0.5)) << ',' << full(quantile(sp, 0.75));
          else o << ",,,";
          o << '\n';
        }
      }
  }

  // landscapes.json: f − f* on a grid of each box, row by row in x₂.
  {
    std::ofstream o(out / "landscapes.json");
    o << "{\"grid\":" << landscape_grid << ",\"note\":\"values: f - f* on a grid of the box, row by row in x2, 4 significant digits\",\"functions\":{";
    for (std::size_t f = 0; f < L.size(); f++) {
      const auto& b = L[f].box;
      std::vector<double> v;
      for (std::size_t j = 0; j < landscape_grid; j++)
        for (std::size_t i = 0; i < landscape_grid; i++) {
          const Real x{b.lower[0] + b.width(0) * static_cast<double>(i) / (landscape_grid - 1), b.lower[1] + b.width(1) * static_cast<double>(j) / (landscape_grid - 1)};
          v.push_back(L[f].f(x) - L[f].fstar);
        }
      o << (f ? "," : "") << '"' << L[f].name << "\":{\"lower\":" << point(b.lower) << ",\"upper\":" << point(b.upper) << ",\"fstar\":" << full(L[f].fstar)
        << ",\"minimisers\":" << points(L[f].minimisers) << ",\"values\":" << numbers(v) << '}';
    }
    o << "}}\n";
  }

  // trajectories.json: seed 0 of every algorithm on every function.
  {
    std::ofstream o(out / "trajectories.json");
    o << "{\"meta\":{\"heuresis\":\"" << heuresis_version() << "\",\"seed\":0,\"budget\":" << budget << ",\"bo_budget\":" << bo_budget << ",\"precision\":" << precision
      << ",\"digits\":4,\"note\":\"Every value is a gap f - f*. Each generation: e = evaluations after it, gap = the best gap so far; pts and f its candidates and their gaps; "
         "m, s, C the state it was drawn from (ES: the parent and sigma; CMA-ES: the mean, sigma and C, whose distribution is N(m, s^2 C); CE: the mean and covariance); v (PSO) "
         "the velocities after the generation (from the second generation on, the steps that brought the particles to pts; in the first, their initial velocities); pop and pf "
         "(DE) the population after selection and its gaps. C lists the entries (1,1), (1,2), (2,2). gp: Bayesian optimisation's surrogate mean (as a gap) and standard deviation "
         "on a grid of the box after the given number of evaluations, row by row in x2.\"},\"functions\":{";
    for (std::size_t f = 0; f < L.size(); f++) {
      o << (f ? "," : "") << '"' << L[f].name << "\":{\"lower\":" << point(L[f].box.lower) << ",\"upper\":" << point(L[f].box.upper) << ",\"fstar\":" << full(L[f].fstar)
        << ",\"minimisers\":" << points(L[f].minimisers) << ",\"start\":" << point(L[f].start) << ",\"sigma0\":" << num(L[f].sigma0) << ",\"algorithms\":{";
      bool first = true;
      for (std::size_t a = 0; a < algorithms.size(); a++) {
        if (algorithms[a] == "random") continue;
        const Trace& t = slot(f, a, 0);
        o << (first ? "" : ",") << '"' << algorithms[a] << "\":{\"ending\":\"" << t.ending << "\",\"generations\":[";
        first = false;
        for (std::size_t k = 0; k < t.generations.size(); k++) {
          const Generation& g = t.generations[k];
          o << (k ? "," : "") << "{\"e\":" << g.evaluations << ",\"gap\":" << num(g.best - L[f].fstar) << ",\"pts\":" << points(g.points) << ",\"f\":" << numbers(gaps(g.values, L[f].fstar));
          if (!g.centre.empty()) o << ",\"m\":" << point(g.centre);
          if (std::isfinite(g.sigma)) o << ",\"s\":" << num(g.sigma);
          if (!g.cov.empty()) o << ",\"C\":" << numbers(g.cov);
          if (!g.velocities.empty()) o << ",\"v\":" << points(g.velocities);
          if (!g.population.empty()) o << ",\"pop\":" << points(g.population) << ",\"pf\":" << numbers(gaps(g.population_values, L[f].fstar));
          o << '}';
        }
        o << ']';
        if (!t.gp.empty()) {
          o << ",\"gp\":{\"grid\":" << gp_grid << ",\"snapshots\":[";
          for (std::size_t k = 0; k < t.gp.size(); k++)
            o << (k ? "," : "") << "{\"count\":" << t.gp[k].count << ",\"mean\":" << numbers(t.gp[k].mean) << ",\"sd\":" << numbers(t.gp[k].sd) << '}';
          o << "]}";
        }
        o << '}';
      }
      o << "}}";
    }
    o << "}}\n";
  }

  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::cerr << "done in " << secs << " s on " << threads << " threads\n";
  return 0;
}
