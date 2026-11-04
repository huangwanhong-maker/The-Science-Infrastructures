// Fronts in three objectives: how NSGA-II, NSGA-III, MOEA/D (Tchebycheff and PBI) and SMS-EMOA approach the
// Pareto fronts of DTLZ1, DTLZ2 and DTLZ7 (three objectives) and of ZDT1 and ZDT3 (two), measured by the exact
// hypervolume and by IGD against evaluations, over 15 seeds.
//
// Every algorithm is heuresis's, driven by hand through the ask–tell protocol so that the state can be measured at
// fixed numbers of evaluations: at each checkpoint the algorithm holds what it had after its last completed update
// (for the generational ones, the last generation told whole). The fronts of ZDT3 and DTLZ7, and their
// hypervolumes, are not in heuresis's benchmarks; they are computed here from the definitions (see `records`):
// on g = 1 a point of either front is a record of a function of one variable, so each front is a set of records
// (ZDT3) or a product of two such sets (DTLZ7), and each hypervolume is a one-dimensional integral of a running
// minimum or maximum.
//
// Usage: fronts <results directory>

#include <heuresis/benchmarks/multiobjective.hpp>
#include <heuresis/multiobjective/multiobjective.hpp>
#include <heuresis/search/draws.hpp>
#include <heuresis/variation/real.hpp>

#include <tyche/descriptive/summary.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numbers>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace mo = heuresis::multiobjective;
namespace hb = heuresis::benchmarks;
using heuresis::search::Draws;
using heuresis::search::Told;
using heuresis::variation::Box;
using heuresis::variation::Real;
using mo::Point;

namespace {

constexpr std::size_t seeds = 15;
constexpr double near = 0.01;  // a point within this Euclidean distance of the true front is "on" it

// ---------------------------------------------------------------------------------------------------------------
// The fronts of ZDT3 and DTLZ7.
//
// ZDT3 at g = 1: f = (t, z(t)), z(t) = 1 − √t − t sin(10πt). The point at t is non-dominated among the curve's
// points exactly when z(t) < z(s) for every s < t: the front is the set of strict records of z from below.
// DTLZ7 (m = 3) at g = 1: f = (t₁, t₂, 6 − φ(t₁) − φ(t₂)), φ(t) = t(1 + sin 3πt). Since f₃ is separable, (t₁, t₂) is
// non-dominated exactly when each tᵢ is a strict record of φ from above (if a′ ≤ a with φ(a′) ≥ φ(a), not both
// equal, then (a′, b) dominates (a, b); conversely a dominating (a′, b′) makes a′ or b′ such a witness). The front
// is T* × T*, four patches (T* is two intervals). With r = (1, 1) and r = (1, 1, 7):
//   HV(ZDT3) = ∫₀¹ (1 − Ψ(y)) dy,       Ψ(y) = min_{t ≤ y} z(t);
//   HV(DTLZ7) = 1 + 2 ∫₀¹ Φ(y) dy,      Φ(y) = max_{t ≤ y} φ(t),
// the measure above the front's lower envelope at each (y₁) or (y₁, y₂).
double zdt3_curve(double t) { return 1 - std::sqrt(t) - t * std::sin(10 * std::numbers::pi * t); }
double dtlz7_phi(double t) { return t * (1 + std::sin(3 * std::numbers::pi * t)); }

// The strict records on the grid t_i = i/(N − 1): of z from below (lower = true) or of φ from above.
std::vector<double> records(std::size_t N, const std::function<double(double)>& h, bool lower) {
  std::vector<double> out;
  double best = lower ? INFINITY : -INFINITY;
  for (std::size_t i = 0; i < N; i++) {
    const double t = static_cast<double>(i) / static_cast<double>(N - 1), v = h(t);
    if (lower ? v < best : v > best) { out.push_back(t); best = v; }
  }
  return out;
}

// ∫₀¹ of the running extreme of h by the midpoint rule on N cells (the running extreme is continuous and monotone).
double running_integral(std::size_t N, const std::function<double(double)>& h, bool lower) {
  double best = lower ? INFINITY : -INFINITY, s = 0;
  for (std::size_t i = 0; i < N; i++) {
    const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(N), v = h(t);
    best = lower ? std::min(best, v) : std::max(best, v);
    s += best;
  }
  // the running extreme is attained up to the cell's midpoint; the error is O(1/N)
  return s / static_cast<double>(N);
}

// A dense, sorted record set and the nearest front point to a, within `near`: the distance and the component
// (segment or patch) of the nearest point, or (∞, −1).
struct Records {
  std::vector<double> t;
  std::vector<int> segment;  // consecutive grid points share a segment
};
Records dense_records(std::size_t N, const std::function<double(double)>& h, bool lower) {
  Records r;
  double best = lower ? INFINITY : -INFINITY;
  std::size_t last = 0;
  int seg = -1;
  for (std::size_t i = 0; i < N; i++) {
    const double t = static_cast<double>(i) / static_cast<double>(N - 1), v = h(t);
    if (lower ? v < best : v > best) {
      if (seg < 0 || i != last + 1) seg++;
      r.t.push_back(t);
      r.segment.push_back(seg);
      last = i;
      best = v;
    }
  }
  return r;
}

std::pair<double, int> nearest_zdt3(const Records& R, const Point& a) {
  double best = INFINITY;
  int seg = -1;
  const auto lo = std::lower_bound(R.t.begin(), R.t.end(), a[0] - near), hi = std::upper_bound(R.t.begin(), R.t.end(), a[0] + near);
  for (auto it = lo; it != hi; ++it) {
    const double t = *it, d = std::hypot(a[0] - t, a[1] - zdt3_curve(t));
    if (d < best) { best = d; seg = R.segment[static_cast<std::size_t>(it - R.t.begin())]; }
  }
  return {best, seg};
}

std::pair<double, int> nearest_dtlz7(const Records& R, const Point& a) {
  double best = INFINITY;
  int patch = -1;
  const auto lo1 = std::lower_bound(R.t.begin(), R.t.end(), a[0] - near), hi1 = std::upper_bound(R.t.begin(), R.t.end(), a[0] + near);
  const auto lo2 = std::lower_bound(R.t.begin(), R.t.end(), a[1] - near), hi2 = std::upper_bound(R.t.begin(), R.t.end(), a[1] + near);
  for (auto i = lo1; i != hi1; ++i)
    for (auto j = lo2; j != hi2; ++j) {
      const double f3 = 6 - dtlz7_phi(*i) - dtlz7_phi(*j);
      const double d = std::sqrt((a[0] - *i) * (a[0] - *i) + (a[1] - *j) * (a[1] - *j) + (a[2] - f3) * (a[2] - f3));
      if (d < best) { best = d; patch = R.segment[static_cast<std::size_t>(i - R.t.begin())] + 2 * R.segment[static_cast<std::size_t>(j - R.t.begin())]; }
    }
  return {best, patch};
}

// ---------------------------------------------------------------------------------------------------------------
// The problems.

struct Problem {
  std::string name;
  std::size_t m = 0, n = 0;
  std::function<std::vector<double>(const Real&)> f;
  Point ref;                      // the hypervolume's reference point
  double hv_true = 0;             // the front's hypervolume
  std::string hv_how;             // how it is known
  double hv_set = NAN;            // the hypervolume of a natural set of about μ points on the front
  double igd_set = NAN;           // and its IGD
  std::string hv_set_how;
  std::vector<Point> igd_ref;     // the reference set of IGD
  double budget = 0;              // evaluations
  std::size_t mu = 0;             // the population of NSGA-II, NSGA-III and SMS-EMOA
  std::vector<Point> directions;  // Das–Dennis: NSGA-III's reference directions, MOEA/D's weights
  std::size_t divisions = 0;
  // The Euclidean distance from a point to the true front, and the component of the front nearest it.
  std::function<std::pair<double, int>(const Point&)> distance;
  double scale = 0;  // the front's intercept, for associating near points with directions (0: not measured)
};

std::vector<Point> normalised(std::vector<Point> w) {
  for (auto& p : w) {
    double s = 0;
    for (double v : p) s += v * v;
    for (double& v : p) v /= std::sqrt(s);
  }
  return w;
}

std::vector<Problem> problems(Records& zdt3_records, Records& dtlz7_records, std::ostringstream& notes) {
  std::vector<Problem> out;
  zdt3_records = dense_records(20001, zdt3_curve, true);
  dtlz7_records = dense_records(20001, dtlz7_phi, false);
  {
    Problem p;
    p.name = "DTLZ1";
    p.m = 3;
    p.n = 3 + 5 - 1;
    p.f = [](const Real& x) { return hb::dtlz1(x, 3); };
    p.ref = {1, 1, 1};
    p.hv_true = hb::dtlz1_hypervolume(3);
    p.hv_how = "closed form 1 - 2^-m/m!";
    for (auto w : mo::das_dennis(3, 99)) {
      for (double& v : w) v *= 0.5;
      p.igd_ref.push_back(w);
    }
    p.budget = 40000;
    p.mu = 92;
    p.divisions = 12;
    p.directions = mo::das_dennis(3, 12);
    p.distance = [](const Point& a) { return std::make_pair(std::fabs(a[0] + a[1] + a[2] - 0.5) / std::sqrt(3.0), 0); };
    p.scale = 0.5;
    {
      std::vector<Point> s;
      for (auto w : p.directions) {
        for (double& v : w) v *= 0.5;
        s.push_back(w);
      }
      p.hv_set = mo::hypervolume(s, p.ref);
      p.igd_set = mo::inverted_generational_distance(s, p.igd_ref);
      p.hv_set_how = "the 91 Das-Dennis directions (p = 12) where they meet the front";
    }
    out.push_back(p);
  }
  {
    Problem p;
    p.name = "DTLZ2";
    p.m = 3;
    p.n = 3 + 10 - 1;
    p.f = [](const Real& x) { return hb::dtlz2(x, 3); };
    p.ref = {1, 1, 1};
    p.hv_true = hb::dtlz2_hypervolume(3);
    p.hv_how = "closed form 1 - pi^(m/2)/(2^m Gamma(m/2+1))";
    p.igd_ref = normalised(mo::das_dennis(3, 99));
    p.budget = 40000;
    p.mu = 92;
    p.divisions = 12;
    p.directions = mo::das_dennis(3, 12);
    p.distance = [](const Point& a) { return std::make_pair(std::fabs(std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) - 1), 0); };
    p.scale = 1;
    p.hv_set = mo::hypervolume(normalised(p.directions), p.ref);
    p.igd_set = mo::inverted_generational_distance(normalised(p.directions), p.igd_ref);
    p.hv_set_how = "the 91 Das-Dennis directions (p = 12) where they meet the front";
    out.push_back(p);
  }
  {
    Problem p;
    p.name = "DTLZ7";
    p.m = 3;
    p.n = 3 + 20 - 1;
    p.f = [](const Real& x) { return hb::dtlz7(x, 3); };
    p.ref = {1, 1, 7};
    const double I = running_integral(4000000, dtlz7_phi, false);
    p.hv_true = 1 + 2 * I;
    p.hv_how = "1 + 2 * integral of the running maximum of t(1 + sin 3 pi t), midpoint rule on 4e6 cells";
    const auto T = records(151, dtlz7_phi, false);
    for (double a : T)
      for (double b : T) p.igd_ref.push_back({a, b, 6 - dtlz7_phi(a) - dtlz7_phi(b)});
    p.budget = 40000;
    p.mu = 92;
    p.divisions = 12;
    p.directions = mo::das_dennis(3, 12);
    const Records* R = &dtlz7_records;
    p.distance = [R](const Point& a) { return nearest_dtlz7(*R, a); };
    out.push_back(p);
    // The integral's convergence, and the exact hypervolume of sampled fronts, which converge to it from below.
    notes << "  \"dtlz7_integral\": {\"cells_1e5\": " << std::setprecision(12) << 1 + 2 * running_integral(100000, dtlz7_phi, false)
          << ", \"cells_1e6\": " << 1 + 2 * running_integral(1000000, dtlz7_phi, false) << ", \"cells_4e6\": " << p.hv_true << "},\n";
    for (std::size_t N : {201, 801, 2001}) {
      const auto TT = records(N, dtlz7_phi, false);
      std::vector<Point> s;
      for (double a : TT)
        for (double b : TT) s.push_back({a, b, 6 - dtlz7_phi(a) - dtlz7_phi(b)});
      notes << "  \"dtlz7_sampled_front_grid_" << N << "\": {\"points\": " << s.size() << ", \"hypervolume\": " << mo::hypervolume(s, p.ref) << "},\n";
    }
  }
  {
    Problem p;
    p.name = "ZDT1";
    p.m = 2;
    p.n = 30;
    p.f = [](const Real& x) { return hb::zdt1(x); };
    p.ref = {1, 1};
    p.hv_true = hb::zdt1_hypervolume();
    p.hv_how = "closed form 2/3";
    for (int i = 0; i < 1000; i++) p.igd_ref.push_back(hb::zdt1_front(i / 999.0));
    p.budget = 25000;
    p.mu = 100;
    p.divisions = 99;
    p.directions = mo::das_dennis(2, 99);
    {
      std::vector<Point> s;
      for (int i = 0; i < 100; i++) s.push_back(hb::zdt1_front(i / 99.0));
      p.hv_set = mo::hypervolume(s, p.ref);
      p.igd_set = mo::inverted_generational_distance(s, p.igd_ref);
      p.hv_set_how = "100 points of the front evenly spaced in f1";
    }
    p.distance = [](const Point& a) {
      // nearest point of f₂ = 1 − √f₁, by a fine scan near a
      double best = INFINITY;
      for (int i = -200; i <= 200; i++) {
        const double t = std::clamp(a[0] + i * near / 200, 0.0, 1.0);
        best = std::min(best, std::hypot(a[0] - t, a[1] - (1 - std::sqrt(t))));
      }
      return std::make_pair(best, 0);
    };
    out.push_back(p);
  }
  {
    Problem p;
    p.name = "ZDT3";
    p.m = 2;
    p.n = 30;
    p.f = [](const Real& x) { return hb::zdt3(x); };
    p.ref = {1, 1};
    const double I = running_integral(4000000, zdt3_curve, true);
    p.hv_true = 1 - I;
    p.hv_how = "integral of 1 minus the running minimum of 1 - sqrt t - t sin 10 pi t, midpoint rule on 4e6 cells";
    for (double t : records(4001, zdt3_curve, true)) p.igd_ref.push_back({t, zdt3_curve(t)});
    p.budget = 25000;
    p.mu = 100;
    p.divisions = 99;
    p.directions = mo::das_dennis(2, 99);
    const Records* R = &zdt3_records;
    p.distance = [R](const Point& a) { return nearest_zdt3(*R, a); };
    out.push_back(p);
    notes << "  \"zdt3_integral\": {\"cells_1e5\": " << 1 - running_integral(100000, zdt3_curve, true) << ", \"cells_1e6\": " << 1 - running_integral(1000000, zdt3_curve, true)
          << ", \"cells_4e6\": " << p.hv_true << "},\n";
    for (std::size_t N : {1001, 10001, 100001}) {
      std::vector<Point> s;
      for (double t : records(N, zdt3_curve, true)) s.push_back({t, zdt3_curve(t)});
      notes << "  \"zdt3_sampled_front_grid_" << N << "\": {\"points\": " << s.size() << ", \"hypervolume\": " << mo::hypervolume(s, p.ref) << "},\n";
    }
  }
  // The separable characterisation of DTLZ7's front, checked against non-dominated sorting of a grid.
  {
    std::vector<Point> grid;
    for (int i = 0; i <= 60; i++)
      for (int j = 0; j <= 60; j++) {
        const double a = i / 60.0, b = j / 60.0;
        grid.push_back({a, b, 6 - dtlz7_phi(a) - dtlz7_phi(b)});
      }
    const auto nd = mo::non_dominated(grid);
    const auto T = records(61, dtlz7_phi, false);
    notes << "  \"dtlz7_check\": {\"grid\": \"61 x 61\", \"non_dominated\": " << nd.size() << ", \"records_squared\": " << T.size() * T.size() << "},\n";
    std::vector<Point> curve;
    for (int i = 0; i <= 2000; i++) curve.push_back({i / 2000.0, zdt3_curve(i / 2000.0)});
    notes << "  \"zdt3_check\": {\"grid\": 2001, \"non_dominated\": " << mo::non_dominated(curve).size() << ", \"records\": " << records(2001, zdt3_curve, true).size() << "},\n";
    notes << "  \"zdt3_segments\": " << zdt3_records.segment.back() + 1 << ", \"dtlz7_intervals\": " << dtlz7_records.segment.back() + 1 << ",\n";
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// Runs.

const std::vector<std::string> algorithms = {"NSGA-II", "NSGA-III", "MOEA/D-Tch", "MOEA/D-PBI", "SMS-EMOA"};

struct Run {
  std::vector<double> hv, igd;   // at the checkpoints
  std::vector<Point> population;  // the final population's objectives
  std::vector<Point> front;       // its non-dominated points
  double seconds = 0;
};

std::vector<double> checkpoints(double budget) {
  std::vector<double> c;
  for (double e = 500; e <= budget + 0.5; e += 500) c.push_back(e);
  return c;
}

template <class A, class Population>
Run drive(A& alg, const Problem& p, Population&& population) {
  Run r;
  const auto cps = checkpoints(p.budget);
  double evals = 0;
  std::size_t next = 0;
  const auto start = std::chrono::steady_clock::now();
  while (next < cps.size()) {
    const auto limit = static_cast<std::size_t>(cps[next] - evals);
    auto c = alg.ask(limit);
    if (c.empty()) throw std::runtime_error("an algorithm offered nothing");
    std::vector<Told> told;
    told.reserve(c.size());
    for (const auto& x : c) told.push_back({x.ticket, heuresis::search::evaluation(p.f(x.genotype))});
    evals += static_cast<double>(c.size());
    alg.tell(told);
    while (next < cps.size() && evals >= cps[next]) {
      const auto front = alg.front();
      r.hv.push_back(mo::hypervolume(front, p.ref));
      r.igd.push_back(mo::inverted_generational_distance(front, p.igd_ref));
      next++;
    }
  }
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  for (const auto& x : population(alg))
    if (!x.evaluation.failed()) r.population.push_back(x.evaluation.objectives);
  r.front = alg.front();
  return r;
}

Run run_one(const Problem& p, const std::string& algorithm, std::uint64_t seed) {
  const Box box = Box::cube(p.n, 0, 1);
  const Draws draws(seed);
  if (algorithm == "NSGA-II") {
    mo::Nsga2<Real> a(mo::real_kernels(box, 20, 20, 0.9), p.mu, draws);
    return drive(a, p, [](const auto& x) { return x.population(); });
  }
  if (algorithm == "NSGA-III") {
    // Deb and Jain 2014, §VI: SBX with η_c = 30 and probability 1, polynomial mutation η_m = 20 at rate 1/n.
    mo::Nsga3<Real> a(mo::real_kernels(box, 30, 20, 1.0), p.directions, p.mu, draws);
    return drive(a, p, [](const auto& x) { return x.population(); });
  }
  if (algorithm == "SMS-EMOA") {
    mo::SmsEmoa<Real> a(mo::real_kernels(box, 20, 20, 0.9), p.mu, draws);
    return drive(a, p, [](const auto& x) { return x.population(); });
  }
  mo::Moead<Real>::Options o;
  o.neighbours = 20;
  o.scalarisation = algorithm == "MOEA/D-PBI" ? mo::Scalarisation::pbi : mo::Scalarisation::tchebycheff;
  o.theta = 5;
  mo::Moead<Real> a(mo::real_kernels(box, 20, 20, 0.9), p.directions, o, draws);
  return drive(a, p, [](const auto& x) { return x.solutions(); });
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: fronts <results directory>\n";
    return 2;
  }
  const std::filesystem::path dir = argv[1];
  std::filesystem::create_directories(dir);
  Records zdt3_records, dtlz7_records;
  std::ostringstream notes;
  notes << std::setprecision(12);
  const auto P = problems(zdt3_records, dtlz7_records, notes);

  // Tasks: every problem, algorithm and seed; SMS-EMOA's first, being the longest.
  struct Task {
    std::size_t p, a, s;
  };
  std::vector<Task> tasks;
  for (std::size_t a = algorithms.size(); a-- > 0;)
    for (std::size_t p = 0; p < P.size(); p++)
      for (std::size_t s = 0; s < seeds; s++) tasks.push_back({p, a, s});
  std::vector<Run> runs(P.size() * algorithms.size() * seeds);
  auto slot = [&](std::size_t p, std::size_t a, std::size_t s) -> Run& { return runs[(p * algorithms.size() + a) * seeds + s]; };
  std::atomic<std::size_t> next{0};
  std::mutex log;
  const unsigned threads = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
  const auto start = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < threads; t++)
    pool.emplace_back([&] {
      for (std::size_t i; (i = next++) < tasks.size();) {
        const auto [p, a, s] = tasks[i];
        slot(p, a, s) = run_one(P[p], algorithms[a], 1 + s);
        if (s == 0) {
          std::lock_guard<std::mutex> g(log);
          std::cerr << P[p].name << " " << algorithms[a] << " seed 1: " << slot(p, a, s).seconds << " s, HV " << slot(p, a, s).hv.back() << "\n";
        }
      }
    });
  for (auto& t : pool) t.join();
  const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  // The curves: quartiles over the seeds of HV / HV* and of IGD at each checkpoint.
  {
    std::ofstream out(dir / "curves.csv");
    out << "problem,algorithm,evaluations,hv_ratio_q1,hv_ratio_median,hv_ratio_q3,igd_q1,igd_median,igd_q3\n";
    out << std::setprecision(8);
    for (std::size_t p = 0; p < P.size(); p++) {
      const auto cps = checkpoints(P[p].budget);
      for (std::size_t a = 0; a < algorithms.size(); a++)
        for (std::size_t k = 0; k < cps.size(); k++) {
          std::vector<double> hv, igd;
          for (std::size_t s = 0; s < seeds; s++) {
            hv.push_back(slot(p, a, s).hv[k] / P[p].hv_true);
            igd.push_back(slot(p, a, s).igd[k]);
          }
          namespace td = tyche::descriptive;
          out << P[p].name << "," << algorithms[a] << "," << cps[k] << "," << td::quantile(hv, 0.25) << "," << td::median(hv) << "," << td::quantile(hv, 0.75) << ","
              << td::quantile(igd, 0.25) << "," << td::median(igd) << "," << td::quantile(igd, 0.75) << "\n";
        }
    }
  }

  // Each run's end: HV, IGD, GD, how many of its non-dominated points lie on the true front, how they share out
  // among the front's components, and (DTLZ1, DTLZ2) how many of the 91 reference directions they occupy.
  std::vector<std::size_t> median_seed(P.size() * algorithms.size());
  {
    std::ofstream out(dir / "final.csv");
    out << "problem,algorithm,seed,hv,hv_ratio,igd,gd,front_size,on_front,components,occupied,seconds\n";
    out << std::setprecision(8);
    for (std::size_t p = 0; p < P.size(); p++)
      for (std::size_t a = 0; a < algorithms.size(); a++) {
        std::vector<std::pair<double, std::size_t>> order;
        for (std::size_t s = 0; s < seeds; s++) {
          const Run& r = slot(p, a, s);
          std::size_t on = 0;
          std::vector<int> comps(P[p].name == "DTLZ7" ? 4 : P[p].name == "ZDT3" ? zdt3_records.segment.back() + 1 : 1, 0);
          std::vector<Point> near_points;
          for (const auto& f : r.front) {
            const auto [d, c] = P[p].distance(f);
            if (d <= near) {
              on++;
              if (c >= 0 && static_cast<std::size_t>(c) < comps.size()) comps[static_cast<std::size_t>(c)]++;
              near_points.push_back(f);
            }
          }
          long occupied = -1;
          if (P[p].scale > 0 && !near_points.empty()) {
            for (auto& q : near_points)
              for (double& v : q) v /= P[p].scale;
            std::vector<std::size_t> line;
            std::vector<double> dist;
            mo::associate(near_points, P[p].directions, line, dist);
            std::sort(line.begin(), line.end());
            occupied = static_cast<long>(std::unique(line.begin(), line.end()) - line.begin());
          } else if (P[p].scale > 0) {
            occupied = 0;
          }
          std::string components;
          for (std::size_t c = 0; c < comps.size(); c++) components += (c ? "|" : "") + std::to_string(comps[c]);
          out << P[p].name << "," << algorithms[a] << "," << 1 + s << "," << r.hv.back() << "," << r.hv.back() / P[p].hv_true << "," << r.igd.back() << ","
              << mo::generational_distance(r.front, P[p].igd_ref) << "," << r.front.size() << "," << on << "," << components << "," << occupied << "," << r.seconds << "\n";
          order.push_back({r.hv.back(), s});
        }
        std::sort(order.begin(), order.end());
        median_seed[p * algorithms.size() + a] = order[seeds / 2].second;
      }
  }

  // The final population of the median seed (by final HV) of every problem and algorithm.
  {
    std::ofstream out(dir / "populations.csv");
    out << "problem,algorithm,seed,f1,f2,f3\n";
    out << std::setprecision(7);
    for (std::size_t p = 0; p < P.size(); p++)
      for (std::size_t a = 0; a < algorithms.size(); a++) {
        const std::size_t s = median_seed[p * algorithms.size() + a];
        for (const auto& f : slot(p, a, s).population) out << P[p].name << "," << algorithms[a] << "," << 1 + s << "," << f[0] << "," << f[1] << "," << (f.size() > 2 ? f[2] : 0.0) << "\n";
      }
  }

  // The true fronts, sampled for drawing, and the reference directions.
  {
    std::ofstream out(dir / "fronts.csv");
    out << "problem,component,f1,f2,f3\n";
    out << std::setprecision(7);
    for (auto w : mo::das_dennis(3, 40)) out << "DTLZ1,0," << 0.5 * w[0] << "," << 0.5 * w[1] << "," << 0.5 * w[2] << "\n";
    for (auto w : normalised(mo::das_dennis(3, 40))) out << "DTLZ2,0," << w[0] << "," << w[1] << "," << w[2] << "\n";
    const Records T = dense_records(81, dtlz7_phi, false);
    for (std::size_t i = 0; i < T.t.size(); i++)
      for (std::size_t j = 0; j < T.t.size(); j++)
        out << "DTLZ7," << T.segment[i] + 2 * T.segment[j] << "," << T.t[i] << "," << T.t[j] << "," << 6 - dtlz7_phi(T.t[i]) - dtlz7_phi(T.t[j]) << "\n";
    for (int i = 0; i <= 400; i++) out << "ZDT1,0," << i / 400.0 << "," << 1 - std::sqrt(i / 400.0) << ",0\n";
    const Records Z = dense_records(4001, zdt3_curve, true);
    for (std::size_t i = 0; i < Z.t.size(); i += 2) out << "ZDT3," << Z.segment[i] << "," << Z.t[i] << "," << zdt3_curve(Z.t[i]) << ",0\n";
  }
  {
    std::ofstream out(dir / "directions.csv");
    out << "problem,w1,w2,w3\n";
    for (const auto& w : mo::das_dennis(3, 12)) out << "three objectives," << w[0] << "," << w[1] << "," << w[2] << "\n";
  }
  {
    std::ofstream out(dir / "reference.json");
    out << std::setprecision(12) << "{\n";
    out << notes.str();
    out << "  \"seeds\": " << seeds << ", \"near\": " << near << ", \"threads\": " << threads << ", \"wall_seconds\": " << total << ",\n";
    out << "  \"problems\": [\n";
    for (std::size_t p = 0; p < P.size(); p++) {
      out << "    {\"name\": \"" << P[p].name << "\", \"m\": " << P[p].m << ", \"n\": " << P[p].n << ", \"budget\": " << P[p].budget << ", \"mu\": " << P[p].mu
          << ", \"weights\": " << P[p].directions.size() << ", \"divisions\": " << P[p].divisions << ", \"reference\": [";
      for (std::size_t i = 0; i < P[p].ref.size(); i++) out << (i ? ", " : "") << P[p].ref[i];
      out << "], \"hv_true\": " << P[p].hv_true << ", \"hv_how\": \"" << P[p].hv_how << "\", ";
      if (!std::isnan(P[p].hv_set)) out << "\"hv_set\": " << P[p].hv_set << ", \"hv_set_ratio\": " << P[p].hv_set / P[p].hv_true << ", \"igd_set\": " << P[p].igd_set << ", \"hv_set_how\": \"" << P[p].hv_set_how << "\", ";
      out << "\"igd_reference_points\": " << P[p].igd_ref.size() << "}" << (p + 1 < P.size() ? "," : "") << "\n";
    }
    out << "  ]\n}\n";
  }
  std::cerr << "done in " << total << " s\n";
  return 0;
}
