// Benchmarking continuous search: ten algorithms of heuresis on eight shifted and rotated functions, in the manner of
// COCO (Hansen, Auger, Ros, Mersmann, Tušar and Brockhoff 2021, "COCO: a platform for comparing continuous optimizers
// in a black-box setting", Optim. Methods Softw. 36, 114–144).
//
// Each function has 15 instances in each of the dimensions 2, 5, 10 and 20: its optimum moved by a shift uniform in
// [−4, 4]ⁿ and its coordinates turned by a rotation drawn from the Haar measure (benchmarks::Rotation, keyed by the
// instance), so that an algorithm cannot profit from an optimum at the centre or from the coordinate axes. The
// separable ellipsoid is only shifted; Schwefel's function, whose structure is its separability, only has its
// coordinates' signs flipped. Every algorithm gets 10⁴·n evaluations and stops at f − f* ≤ 10⁻⁸; its runtimes to the
// eleven targets f − f* ∈ {10², 10¹, …, 10⁻⁸} are read from its progress by heuresis::analysis.
//
// Written to the directory given as the first argument:
//   runs.csv    one line per run: what it spent, the gap it reached, how it ended, its runtime to each target
//   ecdf.csv    the empirical distribution of runtimes over (run, target) pairs (analysis::runtime_distribution),
//               per dimension pooled over the functions and per function, at budgets 10^(k/20)·n
//   ert.csv     the expected running time (analysis::expected_running_time) to each target, with the success rate and
//               Clopper and Pearson's 95% interval (analysis::clopper_pearson)
//   timing.csv  the processor time spent per algorithm and dimension

#include <heuresis/analysis/analysis.hpp>
#include <heuresis/benchmarks/benchmarks.hpp>
#include <heuresis/distributions/distributions.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>
#include <heuresis/swarms/swarms.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
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
using Real = hv::Real;
using Objective = std::function<double(std::span<const double>)>;

namespace {

constexpr std::uint64_t suite_seed = 20261105;
constexpr std::size_t instances = 15;
const std::vector<std::size_t> dimensions = {2, 5, 10, 20};
constexpr double budget_per_n = 1e4;
const std::vector<double> targets = {1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
constexpr std::size_t batch_limit = 1000;  // the most candidates asked at once (a generation is asked in parts)

const std::vector<std::string> functions = {"sphere", "ellipsoid", "rotated ellipsoid", "rosenbrock", "rastrigin", "ackley", "griewank", "schwefel"};
const std::vector<std::string> algorithms = {"random", "(1+1)-ES", "CMA-ES", "IPOP-CMA-ES", "BIPOP-CMA-ES", "DE", "SHADE", "PSO", "SPSO 2011", "CE"};

// Schwefel's function is a sum of g(x) = 418.98… − x sin √|x|; its minimum per coordinate, by golden sections on
// [400, 440], where g is unimodal (its minimiser 420.968746…).
double schwefel_coordinate_minimum() {
  auto g = [](double x) { const double v[1] = {x}; return hb::schwefel(std::span<const double>(v, 1)); };
  double a = 400, b = 440;
  const double r = (std::sqrt(5.0) - 1) / 2;
  double c = b - r * (b - a), d = a + r * (b - a);
  for (int k = 0; k < 200; k++) {
    if (g(c) < g(d)) b = d;
    else a = c;
    c = b - r * (b - a);
    d = a + r * (b - a);
  }
  return g(0.5 * (a + b));
}
const double schwefel_min = schwefel_coordinate_minimum();

struct Instance {
  std::size_t function = 0, n = 0, index = 0;
  Objective f;
  double fstar = 0;
  Real minimiser;
};

Instance instance(std::size_t fi, std::size_t n, std::size_t k) {
  Instance I{fi, n, k, {}, 0, {}};
  const std::uint64_t key = 100 * fi + k;
  calx::random::Stream s = hs::Draws(suite_seed)(hs::tag("experiment.benchmark.shift"), key, n);
  Real shift(n);
  for (auto& v : shift) v = s.uniform(-4, 4);
  const hb::Rotation R(n, suite_seed, key);
  const std::string& name = functions[fi];
  if (name == "sphere") I.f = hb::transformed([](std::span<const double> y) { return hb::sphere(y); }, shift, R);
  else if (name == "ellipsoid") I.f = hb::transformed([](std::span<const double> y) { return hb::ellipsoid(y); }, shift);
  else if (name == "rotated ellipsoid") I.f = hb::transformed([](std::span<const double> y) { return hb::ellipsoid(y); }, shift, R);
  else if (name == "rosenbrock") {
    // Rosenbrock's minimiser (1, …, 1) moved to the origin, so that the instance's minimiser is the shift.
    I.f = hb::transformed([](std::span<const double> y) { Real z(y.begin(), y.end()); for (auto& v : z) v += 1; return hb::rosenbrock(z); }, shift, R);
  } else if (name == "rastrigin") I.f = hb::transformed([](std::span<const double> y) { return hb::rastrigin(y); }, shift, R);
  else if (name == "ackley") I.f = hb::transformed([](std::span<const double> y) { return hb::ackley(y); }, shift, R);
  else if (name == "griewank") {
    // At its classical scale: the box [−5, 5]ⁿ is [−500, 500]ⁿ of Griewank's coordinates.
    I.f = hb::transformed([](std::span<const double> y) { Real z(y.begin(), y.end()); for (auto& v : z) v *= 100; return hb::griewank(z); }, shift, R);
  } else {
    // Schwefel's function on [−500, 500]ⁿ, scaled into [−5, 5]ⁿ, each coordinate's sign flipped with probability ½.
    // Outside its domain −z sin √|z| falls without bound, below f*, so a search on ℝⁿ that left the box would "solve"
    // it there; as in COCO's f20, a penalty 10⁴ Σ max(0, |xᵢ| − 5)² outside the box removes those points. It is 0 in
    // the box and exceeds the most the function can fall below f* outside it (Schwefel's g(z) is below 0 only from
    // |z| ≈ 665 on, where the penalty is above 2.7 × 10⁴ against g ≥ 419 − |z|).
    Real sign(n);
    calx::random::Stream t = hs::Draws(suite_seed)(hs::tag("experiment.benchmark.signs"), key, n);
    for (auto& v : sign) v = t.uniform() < 0.5 ? -1.0 : 1.0;
    I.f = [sign](std::span<const double> x) {
      Real z(x.size());
      double penalty = 0;
      for (std::size_t i = 0; i < x.size(); i++) {
        z[i] = 100 * sign[i] * x[i];
        const double out = std::max(0.0, std::fabs(x[i]) - 5);
        penalty += 1e4 * out * out;
      }
      return hb::schwefel(z) + penalty;
    };
    I.fstar = static_cast<double>(n) * schwefel_min;
    I.minimiser.resize(n);
    for (std::size_t i = 0; i < n; i++) I.minimiser[i] = sign[i] * 4.20968746;
    return I;
  }
  I.minimiser = shift;
  return I;
}

struct Result {
  ha::Run run;  // its progress as gaps f − f*, kept at the first improvement below each target and the last
  std::string ending;
  double seconds = 0;
};

// The progress thinned to what runtimes to the targets need: the first improvement at or below each, and the last.
std::vector<hs::Progress> thin(const std::vector<hs::Progress>& p, double fstar) {
  std::vector<hs::Progress> out;
  std::size_t next = 0;
  for (std::size_t i = 0; i < p.size(); i++) {
    hs::Progress q = p[i];
    q.value -= fstar;
    bool keep = i + 1 == p.size();
    while (next < targets.size() && q.value <= targets[next]) { keep = true; next++; }
    if (keep) out.push_back(q);
  }
  return out;
}

// Asks, evaluates, tells, until 10⁴·n evaluations or f − f* ≤ 10⁻⁸, keeping the best and its improvements with
// search::Elite as search::run does; an algorithm that throws ends its run, and the reason is kept.
template <class A>
Result drive(A& a, const Instance& I) {
  Result r;
  hs::Elite<Real> elite;
  const double limit = budget_per_n * static_cast<double>(I.n), target = I.fstar + targets.back();
  double evals = 0;
  try {
    for (;;) {
      if (elite.value() <= target) { r.ending = "target"; break; }
      if (evals >= limit) { r.ending = "budget"; break; }
      if (a.stopped()) { r.ending = "stopped"; break; }
      auto c = a.ask(std::min<std::size_t>(batch_limit, static_cast<std::size_t>(std::ceil(limit - evals))));
      if (c.empty()) { r.ending = a.stopped() ? "stopped" : "stalled"; break; }
      std::vector<hs::Told> told(c.size());
      for (std::size_t i = 0; i < c.size(); i++) {
        evals += 1;
        told[i] = {c[i].ticket, hs::scalar(I.f(c[i].genotype))};
        elite.offer(hs::Individual<Real>{c[i].ticket, std::move(c[i].genotype), told[i].evaluation, std::move(c[i].origin)}, evals, evals);
      }
      a.tell(std::span<const hs::Told>(told));
    }
  } catch (const std::exception& e) {
    r.ending = std::string("error: ") + e.what();
  }
  r.run = {thin(elite.progress(), I.fstar), evals, evals};
  return r;
}

// The processor time of the calling thread, so that a run's cost does not depend on what else the machine runs.
double thread_seconds() {
  timespec t{};
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return static_cast<double>(t.tv_sec) + 1e-9 * static_cast<double>(t.tv_nsec);
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

Result run_one(const Instance& I, const std::string& name) {
  const std::size_t n = I.n;
  const std::uint64_t seed = hs::mix(suite_seed ^ hs::mix(1000003ull * I.function + 1009ull * n + I.index));
  const hs::Draws draws(seed);
  const hv::Box box = hv::Box::cube(n, -5, 5), start_box = hv::Box::cube(n, -4, 4);
  calx::random::Stream s = draws(hs::tag("experiment.benchmark.start"), 0, 0);
  const Real start = hv::uniform_point(start_box, s);
  const double sigma0 = 2;
  const double t0 = thread_seconds();
  Result r;
  if (name == "random") {
    hs::RandomSearch<Real> a([box](calx::random::Stream& t) { return hv::uniform_point(box, t); }, draws);
    r = drive(a, I);
  } else if (name == "(1+1)-ES") {
    st::OnePlusOneOptions o;
    o.min_sigma = 1e-12;  // stopped once its step has collapsed
    st::OnePlusOneEs a(start, sigma0, o, draws);
    r = drive(a, I);
  } else if (name == "CMA-ES") {
    st::Cma a(start, sigma0, st::CmaOptions{}, draws);
    r = drive(a, I);
    if (r.ending == "stopped") r.ending += std::string(" (") + stop_name(a.stop_reason()) + ")";
  } else if (name == "IPOP-CMA-ES" || name == "BIPOP-CMA-ES") {
    st::RestartOptions o;
    o.regime = name == "IPOP-CMA-ES" ? st::RestartRegime::ipop : st::RestartRegime::bipop;
    o.start = start_box;
    o.sigma0 = sigma0;
    o.max_large_runs = 30;  // the budget, not the number of restarts, ends a run
    st::CmaRestarts a(o, draws);
    r = drive(a, I);
  } else if (name == "DE") {
    st::DifferentialEvolution a(box, st::DeOptions{}, draws);
    r = drive(a, I);
  } else if (name == "SHADE") {
    st::Shade a(box, st::ShadeOptions{}, draws);
    r = drive(a, I);
  } else if (name == "PSO") {
    sw::PsoOptions o;
    o.rule = sw::PsoRule::inertia;
    o.topology = sw::PsoTopology::ring;
    sw::ParticleSwarm a(box, o, draws);
    r = drive(a, I);
  } else if (name == "SPSO 2011") {
    sw::ParticleSwarm a(box, sw::PsoOptions{}, draws);
    r = drive(a, I);
  } else {  // the cross-entropy method with Kroese, Porotsky and Rubinstein's dynamic smoothing
    const std::size_t lambda = std::max<std::size_t>(50, 10 * n), elite = std::max<std::size_t>(10, lambda / 5);
    hd::DistributionSearch a(hd::Gaussian{n}, hd::Gaussian::isotropic(start, sigma0), lambda, hd::CrossEntropy{elite, 0.7, 0.8, 5}, draws);
    r = drive(a, I);
  }
  r.seconds = thread_seconds() - t0;
  return r;
}

std::string full(double v) {
  if (!std::isfinite(v)) return "";
  char b[40];
  std::snprintf(b, sizeof b, "%.6g", v);
  return b;
}

}  // namespace

int main(int argc, char** argv) {
  const std::filesystem::path out = argc > 1 ? argv[1] : "results";
  std::filesystem::create_directories(out);
  const auto t0 = std::chrono::steady_clock::now();
  std::cerr << "Schwefel's minimum per coordinate: " << full(schwefel_min) << "\n";

  // Every (function, dimension, instance, algorithm), the largest dimensions first, on threads.
  struct Job {
    std::size_t f, d, k, a;
  };
  std::vector<Job> jobs;
  for (std::size_t d = dimensions.size(); d-- > 0;)
    for (std::size_t a = 0; a < algorithms.size(); a++)
      for (std::size_t f = 0; f < functions.size(); f++)
        for (std::size_t k = 0; k < instances; k++) jobs.push_back({f, d, k, a});
  auto index = [&](std::size_t f, std::size_t d, std::size_t k, std::size_t a) { return ((f * dimensions.size() + d) * instances + k) * algorithms.size() + a; };
  std::vector<Result> results(jobs.size());
  std::vector<Instance> problems(functions.size() * dimensions.size() * instances);
  for (std::size_t f = 0; f < functions.size(); f++)
    for (std::size_t d = 0; d < dimensions.size(); d++)
      for (std::size_t k = 0; k < instances; k++) problems[(f * dimensions.size() + d) * instances + k] = instance(f, dimensions[d], k);
  // Each instance's minimiser is checked: f(x*) − f* must be far below the last target.
  for (const auto& I : problems) {
    const double gap = I.f(I.minimiser) - I.fstar;
    if (!(std::fabs(gap) < 1e-9)) {
      std::cerr << "instance " << functions[I.function] << " n=" << I.n << " #" << I.index << ": f(x*) − f* = " << gap << "\n";
      return 1;
    }
  }
  std::atomic<std::size_t> next{0}, finished{0};
  std::mutex io;
  const unsigned threads = std::max(1u, std::thread::hardware_concurrency());
  std::vector<std::thread> pool;
  for (unsigned w = 0; w < threads; w++)
    pool.emplace_back([&] {
      for (std::size_t j; (j = next++) < jobs.size();) {
        const Job& b = jobs[j];
        results[index(b.f, b.d, b.k, b.a)] = run_one(problems[(b.f * dimensions.size() + b.d) * instances + b.k], algorithms[b.a]);
        const std::size_t done = ++finished;
        if (done % 400 == 0) {
          std::lock_guard lock(io);
          std::cerr << done << "/" << jobs.size() << " runs, " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
        }
      }
    });
  for (auto& t : pool) t.join();

  // runs.csv
  {
    std::ofstream o(out / "runs.csv");
    o << "function,n,instance,algorithm,evaluations,gap,ending";
    for (double t : targets) o << ",rt_" << full(t);
    o << '\n';
    for (std::size_t f = 0; f < functions.size(); f++)
      for (std::size_t d = 0; d < dimensions.size(); d++)
        for (std::size_t k = 0; k < instances; k++)
          for (std::size_t a = 0; a < algorithms.size(); a++) {
            const Result& r = results[index(f, d, k, a)];
            o << functions[f] << ',' << dimensions[d] << ',' << k << ',' << algorithms[a] << ',' << r.run.evaluations << ','
              << (r.run.progress.empty() ? "" : full(r.run.progress.back().value)) << ',' << r.ending;
            for (double t : targets) {
              const auto rt = ha::runtime(r.run, t);
              o << ',' << (rt ? full(*rt) : "");
            }
            o << '\n';
          }
  }

  // ecdf.csv: per dimension over all functions ("all"), and per function.
  {
    std::ofstream o(out / "ecdf.csv");
    o << "function,n,algorithm,budget_per_n,fraction\n";
    std::vector<double> per_n;
    for (int k = 0; k <= 80; k++) per_n.push_back(std::pow(10.0, k / 20.0));
    for (std::size_t d = 0; d < dimensions.size(); d++)
      for (std::size_t a = 0; a < algorithms.size(); a++) {
        std::vector<double> b;
        for (double v : per_n) b.push_back(v * static_cast<double>(dimensions[d]));
        std::vector<ha::Run> all;
        for (std::size_t f = 0; f <= functions.size(); f++) {
          std::vector<ha::Run> runs;
          if (f < functions.size())
            for (std::size_t k = 0; k < instances; k++) runs.push_back(results[index(f, d, k, a)].run);
          else
            runs = all;
          const auto F = ha::runtime_distribution(runs, targets, b);
          for (std::size_t i = 0; i < b.size(); i++) o << (f < functions.size() ? functions[f] : "all") << ',' << dimensions[d] << ',' << algorithms[a] << ',' << full(per_n[i]) << ',' << full(F[i]) << '\n';
          if (f < functions.size()) all.insert(all.end(), runs.begin(), runs.end());
        }
      }
  }

  // ert.csv
  {
    std::ofstream o(out / "ert.csv");
    o << "function,n,algorithm,target,ert,successes,runs,rate,lower,upper\n";
    for (std::size_t f = 0; f < functions.size(); f++)
      for (std::size_t d = 0; d < dimensions.size(); d++)
        for (std::size_t a = 0; a < algorithms.size(); a++) {
          std::vector<ha::Run> runs;
          for (std::size_t k = 0; k < instances; k++) runs.push_back(results[index(f, d, k, a)].run);
          for (double t : targets) {
            const auto e = ha::expected_running_time(runs, t);
            const auto p = ha::clopper_pearson(e.successes, e.runs);
            o << functions[f] << ',' << dimensions[d] << ',' << algorithms[a] << ',' << full(t) << ',' << full(e.value) << ',' << e.successes << ',' << e.runs << ','
              << full(p.estimate) << ',' << full(p.lower) << ',' << full(p.upper) << '\n';
          }
        }
  }

  // timing.csv
  {
    std::ofstream o(out / "timing.csv");
    o << "algorithm,n,cpu_seconds,evaluations\n";
    for (std::size_t a = 0; a < algorithms.size(); a++)
      for (std::size_t d = 0; d < dimensions.size(); d++) {
        double secs = 0, evals = 0;
        for (std::size_t f = 0; f < functions.size(); f++)
          for (std::size_t k = 0; k < instances; k++) {
            secs += results[index(f, d, k, a)].seconds;
            evals += results[index(f, d, k, a)].run.evaluations;
          }
        o << algorithms[a] << ',' << dimensions[d] << ',' << full(secs) << ',' << full(evals) << '\n';
      }
  }

  std::cerr << "done in " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s on " << threads << " threads\n";
  return 0;
}
