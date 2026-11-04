// Repertoires of an arm: quality-diversity on heuresis's planar arm of twelve joints.
//
// The arm (heuresis::benchmarks::planar_arm) has d = 12 links of length 1/d from a base at the origin; its genotype
// is the twelve joint angles in [−π, π], its descriptor the tip mapped from [−1, 1]² to [0, 1]², and its objective
// the variance of the angles (an arm that bends evenly is the better one). The question is how much of the
// reachable disc each search illuminates, and how good the arm it keeps for each place is: MAP-Elites on a grid with
// Gaussian mutation at three step sizes and with the iso+line variation, CVT-MAP-Elites on Voronoi cells of sites
// drawn in the reachable disc, novelty search, random search, and CMA-ES with IPOP restarts minimising the variance
// alone. Every evaluation of every search is also offered to a passive grid archive of 50 × 50 cells, so that all
// are measured on the same cells: coverage (the share of the cells meeting the reachable disc that hold an elite)
// and the QD-score with offset π² (an upper bound of the variance on the box), divided by π² times those cells.
//
// Usage: repertoires <results directory>

#include <heuresis/benchmarks/combinatorial.hpp>
#include <heuresis/diversity/diversity.hpp>
#include <heuresis/search/random_search.hpp>
#include <heuresis/search/run.hpp>
#include <heuresis/strategies/restarts.hpp>
#include <heuresis/variation/real.hpp>

#include <tyche/descriptive/summary.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numbers>
#include <string>
#include <thread>
#include <vector>

namespace dv = heuresis::diversity;
using heuresis::search::Draws;
using heuresis::search::Told;
using heuresis::search::tag;
using heuresis::variation::Box;
using heuresis::variation::Real;

namespace {

constexpr std::size_t joints = 12;
constexpr std::size_t seeds = 15;
constexpr double budget = 100000;
constexpr std::size_t resolution = 50;
constexpr std::size_t sites_count = 1000;
const double offset = std::numbers::pi * std::numbers::pi;  // the variance of angles in [−π, π] is at most π²

const Box box = Box::cube(joints, -std::numbers::pi, std::numbers::pi);
const dv::Grid grid({0, 0}, {1, 1}, {resolution, resolution});

heuresis::search::Evaluation arm(const Real& x) { return heuresis::benchmarks::planar_arm(x); }

// Whether a grid cell meets the open disc of radius ½ about (½, ½), the arm's reachable set in descriptor space.
bool reachable(std::size_t c) {
  const auto x = grid.centre(c);
  const double h = 0.5 / static_cast<double>(resolution);
  const double dx = std::max({x[0] - h - 0.5, 0.0, 0.5 - (x[0] + h)}), dy = std::max({x[1] - h - 0.5, 0.0, 0.5 - (x[1] + h)});
  return std::hypot(dx, dy) < 0.5;
}

// The CVT sites: uniform in the disc of radius ½ (rejection from the square, one keyed stream). heuresis asks the
// caller for the sites until tyche has k-means; uniform sites give cells of equal volume in expectation only.
std::vector<std::vector<double>> disc_sites() {
  std::vector<std::vector<double>> s;
  calx::random::Stream r = Draws(2026)(tag("experiment.arm.sites"), 0, 0);
  while (s.size() < sites_count) {
    const double x = r.uniform(), y = r.uniform();
    if (std::hypot(x - 0.5, y - 0.5) < 0.5) s.push_back({x, y});
  }
  return s;
}

dv::Kernels<Real> kernels(const std::string& variation, double sigma) {
  dv::Kernels<Real> k;
  k.initialise = [](calx::random::Stream& s) { return heuresis::variation::uniform_point(box, s); };
  if (variation == "gaussian") {
    k.mutate = [sigma](const Real& x, calx::random::Stream& s) {
      Real y = heuresis::variation::gaussian_mutation(x, sigma, s);
      heuresis::variation::repair(y, box, heuresis::variation::Repair::reflect, s);
      return y;
    };
  } else {
    // iso+line (Vassiliades and Mouret 2018): x + σ_iso N(0, I) + σ_line N(0, 1)(y − x), repaired by reflection; the
    // mutation after it is the identity.
    k.vary = [](const Real& a, const Real& b, calx::random::Stream& s) {
      Real y = heuresis::variation::iso_line(a, b, 0.05, 0.3, s);
      heuresis::variation::repair(y, box, heuresis::variation::Repair::reflect, s);
      return y;
    };
    k.mutate = [](const Real& x, calx::random::Stream&) { return x; };
  }
  return k;
}

struct Algorithm {
  std::string name;
  std::string kind;       // map-elites, cvt, novelty, random, cma
  std::string variation;  // gaussian, iso+line
  double sigma = 0;
};

const std::vector<Algorithm> algorithms = {
    {"MAP-Elites, iso+line", "map-elites", "iso+line", 0},
    {"MAP-Elites, Gaussian 0.03", "map-elites", "gaussian", 0.03},
    {"MAP-Elites, Gaussian 0.1", "map-elites", "gaussian", 0.1},
    {"MAP-Elites, Gaussian 0.3", "map-elites", "gaussian", 0.3},
    {"CVT-MAP-Elites, iso+line", "cvt", "iso+line", 0},
    {"novelty search, Gaussian 0.1", "novelty", "gaussian", 0.1},
    {"random search", "random", "", 0},
    {"CMA-ES, IPOP", "cma", "", 0},
};

std::vector<double> checkpoints() {
  std::vector<double> c = {100, 200, 500};
  for (double e = 1000; e <= budget + 0.5; e += 1000) c.push_back(e);
  return c;
}

struct Elite {
  std::size_t cell = 0;
  double variance = 0;
  Real angles;
};

struct Run {
  std::vector<double> coverage, qd, mean_variance, own_coverage;  // at the checkpoints
  std::vector<Elite> elites;      // the passive grid archive at the end
  std::vector<Elite> cvt_elites;  // CVT-MAP-Elites' own archive at the end (cell = site index)
  double best_variance = INFINITY;
  std::size_t novelty_archive = 0, restarts = 0;
  double seconds = 0;
};

std::size_t reachable_cells() {
  std::size_t r = 0;
  for (std::size_t c = 0; c < grid.cells(); c++) r += reachable(c);
  return r;
}

template <class A, class Own>
Run drive(A& alg, Own&& own_coverage) {
  Run r;
  dv::Archive<Real, dv::Grid> passive(grid);
  const auto cps = checkpoints();
  const double R = static_cast<double>(reachable_cells());
  double evals = 0;
  std::size_t next = 0;
  const auto start = std::chrono::steady_clock::now();
  auto measure = [&] {
    double total = 0;
    for (auto c : passive.filled()) total += passive.elite(c)->evaluation.value();
    r.coverage.push_back(static_cast<double>(passive.size()) / R);
    r.qd.push_back(passive.qd_score(offset) / (offset * R));
    r.mean_variance.push_back(passive.size() ? total / static_cast<double>(passive.size()) : NAN);
    r.own_coverage.push_back(own_coverage(alg));
  };
  while (next < cps.size()) {
    const auto limit = static_cast<std::size_t>(cps[next] - evals);
    auto c = alg.ask(limit);
    if (c.empty()) {
      // the search stopped (CMA-ES out of restarts): the remaining checkpoints keep its last state
      while (next < cps.size()) { measure(); next++; }
      break;
    }
    std::vector<Told> told;
    told.reserve(c.size());
    for (auto& x : c) {
      auto e = arm(x.genotype);
      r.best_variance = std::min(r.best_variance, e.value());
      told.push_back({x.ticket, e});
      passive.offer({x.ticket, x.genotype, std::move(e), x.origin});
    }
    evals += static_cast<double>(c.size());
    alg.tell(told);
    while (next < cps.size() && evals >= cps[next]) { measure(); next++; }
  }
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  for (auto c : passive.filled()) r.elites.push_back({c, passive.elite(c)->evaluation.value(), passive.elite(c)->genotype});
  return r;
}

Run run_one(const Algorithm& a, std::uint64_t seed, const std::vector<std::vector<double>>& sites) {
  const Draws draws(seed);
  auto none = [](const auto&) { return NAN; };
  if (a.kind == "map-elites") {
    dv::MapElites<Real, dv::Grid>::Options o;
    o.initial = 1000;
    o.variation_probability = a.variation == "iso+line" ? 1.0 : 0.0;
    dv::MapElites<Real, dv::Grid> me(kernels(a.variation, a.sigma), grid, o, draws);
    return drive(me, [](const auto& m) { return m.archive().coverage(); });
  }
  if (a.kind == "cvt") {
    dv::MapElites<Real, dv::Voronoi>::Options o;
    o.initial = 1000;
    o.variation_probability = 1.0;
    dv::MapElites<Real, dv::Voronoi> me(kernels(a.variation, a.sigma), dv::Voronoi(sites), o, draws);
    Run r = drive(me, [](const auto& m) { return m.archive().coverage(); });
    for (auto c : me.archive().filled()) r.cvt_elites.push_back({c, me.archive().elite(c)->evaluation.value(), me.archive().elite(c)->genotype});
    return r;
  }
  if (a.kind == "novelty") {
    dv::NoveltySearch<Real>::Options o;
    o.population = 100;
    o.neighbours = 15;
    o.threshold = 0.02;
    dv::NoveltySearch<Real> ns(kernels(a.variation, a.sigma), o, draws);
    Run r = drive(ns, none);
    r.novelty_archive = ns.archive().size();
    return r;
  }
  if (a.kind == "random") {
    heuresis::search::RandomSearch<Real> rs([](calx::random::Stream& s) { return heuresis::variation::uniform_point(box, s); }, draws);
    return drive(rs, none);
  }
  heuresis::strategies::RestartOptions o;
  o.regime = heuresis::strategies::RestartRegime::ipop;
  o.start = box;
  o.sigma0 = 1.0;  // a sixth of the range 2π
  o.max_large_runs = 64;
  o.cma.bounds = box;
  heuresis::strategies::CmaRestarts cma(o, draws);
  Run r = drive(cma, none);
  r.restarts = cma.runs();
  return r;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: repertoires <results directory>\n";
    return 2;
  }
  const std::filesystem::path dir = argv[1];
  std::filesystem::create_directories(dir);
  const auto sites = disc_sites();
  const std::size_t R = reachable_cells();

  std::vector<Run> runs(algorithms.size() * seeds);
  auto slot = [&](std::size_t a, std::size_t s) -> Run& { return runs[a * seeds + s]; };
  std::atomic<std::size_t> next{0};
  std::mutex log;
  const unsigned threads = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
  const auto start = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < threads; t++)
    pool.emplace_back([&] {
      for (std::size_t i; (i = next++) < runs.size();) {
        const std::size_t a = i / seeds, s = i % seeds;
        slot(a, s) = run_one(algorithms[a], 1 + s, sites);
        if (s == 0) {
          std::lock_guard<std::mutex> g(log);
          std::cerr << algorithms[a].name << " seed 1: " << slot(a, s).seconds << " s, coverage " << slot(a, s).coverage.back() << "\n";
        }
      }
    });
  for (auto& t : pool) t.join();
  const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  namespace td = tyche::descriptive;
  const auto cps = checkpoints();

  {
    std::ofstream out(dir / "curves.csv");
    out << "algorithm,evaluations,coverage_q1,coverage_median,coverage_q3,qd_q1,qd_median,qd_q3,mean_variance_median,own_coverage_median\n";
    out << std::setprecision(7);
    for (std::size_t a = 0; a < algorithms.size(); a++)
      for (std::size_t k = 0; k < cps.size(); k++) {
        std::vector<double> c, q, v, o;
        for (std::size_t s = 0; s < seeds; s++) {
          c.push_back(slot(a, s).coverage[k]);
          q.push_back(slot(a, s).qd[k]);
          v.push_back(slot(a, s).mean_variance[k]);
          o.push_back(slot(a, s).own_coverage[k]);
        }
        out << "\"" << algorithms[a].name << "\"," << cps[k] << "," << td::quantile(c, 0.25) << "," << td::median(c) << "," << td::quantile(c, 0.75) << "," << td::quantile(q, 0.25) << ","
            << td::median(q) << "," << td::quantile(q, 0.75) << "," << td::median(v) << ",";
        if (std::isnan(o[0])) out << "";
        else out << td::median(o);
        out << "\n";
      }
  }

  std::vector<std::size_t> median_seed(algorithms.size());
  {
    std::ofstream out(dir / "final.csv");
    out << "algorithm,seed,coverage_10k,coverage,qd_10k,qd,mean_variance,own_coverage,best_variance,cells_filled,novelty_archive,restarts,seconds\n";
    out << std::setprecision(7);
    const std::size_t k10 = static_cast<std::size_t>(std::find(cps.begin(), cps.end(), 10000.0) - cps.begin());
    for (std::size_t a = 0; a < algorithms.size(); a++) {
      std::vector<std::pair<double, std::size_t>> order;
      for (std::size_t s = 0; s < seeds; s++) {
        const Run& r = slot(a, s);
        out << "\"" << algorithms[a].name << "\"," << 1 + s << "," << r.coverage[k10] << "," << r.coverage.back() << "," << r.qd[k10] << "," << r.qd.back() << "," << r.mean_variance.back() << ",";
        if (!std::isnan(r.own_coverage.back())) out << r.own_coverage.back();
        out << "," << r.best_variance << "," << r.elites.size() << "," << r.novelty_archive << "," << r.restarts << "," << r.seconds << "\n";
        order.push_back({r.qd.back(), s});
      }
      std::sort(order.begin(), order.end());
      median_seed[a] = order[seeds / 2].second;
    }
  }

  // The passive archives of the median seed (by final QD-score): each elite's cell and variance.
  {
    std::ofstream out(dir / "elites.csv");
    out << "algorithm,seed,i,j,variance\n";
    out << std::setprecision(6);
    for (std::size_t a = 0; a < algorithms.size(); a++) {
      if (algorithms[a].name == "MAP-Elites, Gaussian 0.03" || algorithms[a].name == "MAP-Elites, Gaussian 0.3") continue;
      const std::size_t s = median_seed[a];
      for (const auto& e : slot(a, s).elites) out << "\"" << algorithms[a].name << "\"," << 1 + s << "," << e.cell % resolution << "," << e.cell / resolution << "," << e.variance << "\n";
    }
  }
  // CVT-MAP-Elites' own archive of the median seed, on its sites.
  {
    std::ofstream out(dir / "cvt.csv");
    out << "site,x,y,filled,variance\n";
    out << std::setprecision(6);
    const std::size_t a = 4, s = median_seed[a];
    std::vector<double> v(sites.size(), NAN);
    for (const auto& e : slot(a, s).cvt_elites) v[e.cell] = e.variance;
    for (std::size_t j = 0; j < sites.size(); j++) out << j << "," << sites[j][0] << "," << sites[j][1] << "," << !std::isnan(v[j]) << "," << (std::isnan(v[j]) ? 0.0 : v[j]) << "\n";
  }
  // Arms for chosen places: the elites of MAP-Elites (iso+line, median seed) whose cells hold points on rings about
  // the base, and the zero-variance arms (every angle equal), which trace one curve.
  {
    std::ofstream out(dir / "arms.csv");
    out << "ring,direction,i,j,variance";
    for (std::size_t q = 0; q < joints; q++) out << ",theta" << q + 1;
    out << "\n";
    out << std::setprecision(7);
    const std::size_t s = median_seed[0];
    std::vector<const Elite*> by_cell(grid.cells(), nullptr);
    for (const auto& e : slot(0, s).elites) by_cell[e.cell] = &e;
    const std::vector<double> rings = {0.06, 0.16, 0.26, 0.36, 0.46};
    for (std::size_t ri = 0; ri < rings.size(); ri++)
      for (std::size_t k = 0; k < 16; k++) {
        const double phi = 2 * std::numbers::pi * (static_cast<double>(k) + 0.5 * static_cast<double>(ri % 2)) / 16;
        const std::vector<double> d = {0.5 + rings[ri] * std::cos(phi), 0.5 + rings[ri] * std::sin(phi)};
        const std::size_t c = grid.cell(d);
        if (!by_cell[c]) continue;
        out << ri << "," << k << "," << c % resolution << "," << c / resolution << "," << by_cell[c]->variance;
        for (double t : by_cell[c]->angles) out << "," << t;
        out << "\n";
      }
  }
  {
    std::ofstream out(dir / "reference.json");
    out << std::setprecision(10) << "{\n";
    out << "  \"joints\": " << joints << ", \"angle_box\": [" << -std::numbers::pi << ", " << std::numbers::pi << "], \"seeds\": " << seeds << ", \"budget\": " << budget << ",\n";
    out << "  \"grid\": [" << resolution << ", " << resolution << "], \"reachable_cells\": " << R << ", \"cvt_sites\": " << sites.size() << ", \"qd_offset\": " << offset << ",\n";
    out << "  \"map_elites_initial\": 1000, \"iso_line\": [0.05, 0.3], \"novelty\": {\"population\": 100, \"neighbours\": 15, \"threshold\": 0.02},\n";
    out << "  \"cma\": {\"regime\": \"ipop\", \"sigma0\": 1.0, \"bounds\": \"transform\", \"max_large_runs\": 64},\n";
    out << "  \"median_seeds\": {";
    for (std::size_t a = 0; a < algorithms.size(); a++) out << (a ? ", " : "") << "\"" << algorithms[a].name << "\": " << 1 + median_seed[a];
    out << "},\n  \"threads\": " << threads << ", \"wall_seconds\": " << total << "\n}\n";
  }
  std::cerr << "done in " << total << " s\n";
  return 0;
}
