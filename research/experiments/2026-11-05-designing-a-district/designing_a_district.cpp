// Designing a district: the redesign of a few blocks of Melbourne's Hoddle Grid as a heterogeneous, many-objective
// search, on the search foundation heuresis and the measures of the architecture and urban-science domains
// (README.md; the problem in district.hpp).
//
//   designing_a_district <results/district-inputs.json> <results directory> [--threads 4] [--scale 1] [--profile]
//
// --scale multiplies the budgets (at 1: 210 generations of 210 for the five-objective searches, of 100 for NSGA-II,
// 84,000 evaluations for MAP-Elites); --profile times 400 random designs and stops.
//
// Runs, each on the same problem and its constraints:
//   1. NSGA-III on the five objectives, three seeds;
//   2. MOEA/D (Tchebycheff, the same 210 directions) on the five objectives, three seeds;
//   3. random search with the same budget, three seeds (the baseline);
//   4. NSGA-II on two of them (floor area and winter sun), three seeds;
//   5. MAP-Elites over plot ratio and open share, winter sun its quality.
// Every random draw is keyed by its purpose (heuresis::search::Draws), so each run is the same on any number of threads.
// Writes the results/ tables from which README.md and figures.py draw everything.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <calx/scheduler/thread_pool.hpp>
#include <heuresis/diversity/diversity.hpp>
#include <heuresis/multiobjective/multiobjective.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/variation/variation.hpp>

#include "district.hpp"

namespace {

using namespace study;
namespace mo = heuresis::multiobjective;
namespace dv = heuresis::diversity;
using Individual = hs::Individual<Design>;
using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// ---- the kernels of the mixed space (composed by hand: heuresis variation/README.md §6) -----------------------------

struct Operators {
  const Problem* P = nullptr;
  hv::IntegerBox storeys, uses;
  std::size_t n = 0, k = 0;
  double crossover_probability = 0.9;
  std::int64_t creep = 6;

  explicit Operators(const Problem& p) : P(&p), storeys(p.storey_box()), uses(p.use_box()), n(p.parcels.size()), k(p.links.size()) {}

  Design sample(calx::random::Stream& s) const {
    Design d;
    d.storeys = hv::uniform_integers(storeys, s);
    d.uses = hv::uniform_integers(uses, s);
    d.links = hv::random_bits(k, s);
    return d;
  }
  // Storeys crept by at most `creep` at two parcels in n on average; a use redrawn at one in n; a link flipped at
  // one in k.
  Design mutate(const Design& x, calx::random::Stream& s) const {
    Design y;
    y.storeys = hv::creep_mutation(x.storeys, storeys, creep, 2.0 / static_cast<double>(n), s);
    y.uses = hv::reset_mutation(x.uses, uses, 1.0 / static_cast<double>(n), s);
    y.links = x.links;
    if (k) hv::standard_bit_mutation(y.links, 1.0 / static_cast<double>(k), s);
    return y;
  }
  // A parcel's storeys and use travel together; parcels are in order along the area, so the two-point crossover
  // exchanges a stretch of the area between two cuts across it. Links cross uniformly.
  std::pair<Design, Design> crossover(const Design& a, const Design& b, calx::random::Stream& s) const {
    if (s.uniform() >= crossover_probability) return {a, b};
    std::vector<std::array<std::int64_t, 2>> za(n), zb(n);
    for (std::size_t i = 0; i < n; i++) za[i] = {a.storeys[i], a.uses[i]}, zb[i] = {b.storeys[i], b.uses[i]};
    auto [ca, cb] = hv::two_point_crossover(za, zb, s);
    auto [pa, pb] = hv::uniform_crossover(a.links, b.links, s);
    Design x, y;
    for (std::size_t i = 0; i < n; i++) x.storeys.push_back(ca[i][0]), x.uses.push_back(ca[i][1]), y.storeys.push_back(cb[i][0]), y.uses.push_back(cb[i][1]);
    x.links = std::move(pa), y.links = std::move(pb);
    return {x, y};
  }
  mo::Kernels<Design> multiobjective() const {
    return {[this](calx::random::Stream& s) { return sample(s); },
            [this](const Design& a, const Design& b, calx::random::Stream& s) { return crossover(a, b, s); },
            [this](const Design& x, calx::random::Stream& s) { return mutate(x, s); }};
  }
  dv::Kernels<Design> repertoire() const {
    return {[this](calx::random::Stream& s) { return sample(s); },
            [this](const Design& x, calx::random::Stream& s) { return mutate(x, s); },
            [this](const Design& a, const Design& b, calx::random::Stream& s) { return crossover(a, b, s).first; }};
  }
};

// ---- fronts ------------------------------------------------------------------------------------------------------

// The feasible individuals no other feasible one dominates.
std::vector<Individual> feasible_front(const std::vector<Individual>& pop) {
  std::vector<std::size_t> ok;
  for (std::size_t i = 0; i < pop.size(); i++)
    if (pop[i].evaluation.feasible()) ok.push_back(i);
  std::vector<Individual> out;
  for (auto i : ok) {
    bool dominated = false;
    for (auto j : ok)
      if (j != i && mo::dominates(pop[j].evaluation.objectives, pop[i].evaluation.objectives)) { dominated = true; break; }
    if (!dominated) {
      bool twin = false;  // one of equal objective vectors kept
      for (const auto& x : out) twin = twin || x.evaluation.objectives == pop[i].evaluation.objectives;
      if (!twin) out.push_back(pop[i]);
    }
  }
  return out;
}

std::vector<mo::Point> points_of(const std::vector<Individual>& xs, std::size_t m = 0) {
  std::vector<mo::Point> out;
  for (const auto& x : xs) {
    auto f = x.evaluation.objectives;
    if (m) f.resize(m);
    out.push_back(f);
  }
  return out;
}

struct Checkpoint {
  std::string algorithm;
  std::uint64_t seed = 0;
  double evaluations = 0, seconds = 0;
  std::vector<mo::Point> front;
};

struct RunOutcome {
  std::string algorithm;
  std::uint64_t seed = 0;
  std::vector<Individual> front;
  std::vector<Checkpoint> checkpoints;
  double evaluations = 0, seconds = 0;
};

template <class Algorithm, class Objective>
RunOutcome run_generational(const std::string& name, std::uint64_t seed, Algorithm& alg, Objective&& f, std::size_t mu, std::size_t generations, std::size_t every,
                            calx::scheduler::ThreadPool& pool, const std::function<std::vector<Individual>(const Algorithm&)>& population) {
  RunOutcome out;
  out.algorithm = name, out.seed = seed;
  const auto t0 = Clock::now();
  for (std::size_t g = 0; g < generations; g += every) {
    const std::size_t chunk = std::min(every, generations - g);
    const auto o = hs::run(alg, f, hs::Budget{static_cast<double>(mu * chunk)}, hs::Pool(pool));
    out.evaluations += o.evaluations;
    out.checkpoints.push_back({name, seed, out.evaluations, seconds_since(t0), points_of(feasible_front(population(alg)))});
  }
  out.front = feasible_front(population(alg));
  out.seconds = seconds_since(t0);
  return out;
}

// Merges individuals into a feasible non-dominated archive.
void merge(std::vector<Individual>& archive, const std::vector<Individual>& more) {
  std::vector<Individual> all = archive;
  for (const auto& x : more)
    if (x.evaluation.feasible()) all.push_back(x);
  archive = feasible_front(all);
}

// ---- output helpers ---------------------------------------------------------------------------------------------

Json numbers(const std::vector<double>& v, int digits = 6) {
  Json a = Json::array();
  const double p = std::pow(10.0, digits);
  for (double x : v) a.push(std::isfinite(x) ? Json(std::round(x * p) / p) : Json());
  return a;
}
Json ints(const hv::Integers& v) {
  Json a = Json::array();
  for (auto x : v) a.push(Json(static_cast<double>(x)));
  return a;
}
Json bits(const hv::Bits& v) {
  Json a = Json::array();
  for (auto x : v) a.push(Json(static_cast<double>(x)));
  return a;
}
Json design_json(const Design& d) { return Json::object({{"storeys", ints(d.storeys)}, {"uses", ints(d.uses)}, {"links", bits(d.links)}}); }
Json measures_json(const Measures& m) {
  return Json::object({{"floor_area_m2", std::round(m.floor_area)}, {"winter_sun_m2h", std::round(m.sun)}, {"heating_kwh_per_m2_day", std::round(m.heat * 1e5) / 1e5},
                       {"walk_m", std::round(m.access * 10) / 10}, {"walk_to_shop_m", std::round(m.access_shop * 10) / 10}, {"walk_to_open_space_m", std::round(m.access_open * 10) / 10},
                       {"mix_nats", std::round(m.mix * 1e5) / 1e5}, {"use_entropy_nats", std::round(m.mix_entropy * 1e5) / 1e5}, {"departure_from_mixing_nats", std::round(m.mix_departure * 1e5) / 1e5},
                       {"open_share", std::round(m.open_share * 1e5) / 1e5}, {"dwelling_share", std::round(m.dwelling_share * 1e5) / 1e5}, {"shop_share", std::round(m.shop_share * 1e5) / 1e5},
                       {"plot_ratio", std::round(m.plot_ratio * 1e4) / 1e4}, {"floor_by_use_m2", numbers({m.floor_by_use[0], m.floor_by_use[1], m.floor_by_use[2]}, 0)},
                       {"heating_kwh_day", std::round(m.heat_total)}, {"buildings", static_cast<double>(m.built)}, {"open_parcels", static_cast<double>(m.open)},
                       {"links", static_cast<double>(m.links)}});
}
Json ring_json(const Ring& r) {
  Json a = Json::array();
  for (const auto& p : r) a.push(Json::array({std::round(p[0] * 100) / 100, std::round(p[1] * 100) / 100}));
  return a;
}
void write(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  f << text;
  if (!f) throw std::runtime_error("cannot write " + path);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: designing_a_district <district-inputs.json> <results directory> [--threads N] [--scale s]\n";
    return 2;
  }
  const std::string inputs = argv[1], out = argv[2];
  std::size_t threads = 4;
  double scale = 1;
  bool profile = false;
  for (int i = 3; i < argc; i++) {
    if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) threads = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--scale") == 0 && i + 1 < argc) scale = std::stod(argv[++i]);
    else if (std::strcmp(argv[i], "--profile") == 0) profile = true;
  }
  const auto T0 = Clock::now();
  std::ifstream in(inputs, std::ios::binary);
  std::stringstream buffer;
  buffer << in.rdbuf();
  const Json J = Json::parse(buffer.str());
  const Problem P(J);
  const Operators ops(P);
  const std::size_t n = P.parcels.size(), K = P.links.size();
  std::printf("district: %zu parcels, %zu candidate links, %zu ground samples, %zu places; %zu segments dropped by planarisation\n", n, K, P.shade.points.size(),
              P.place_count, P.dropped.size());
  const auto& R = P.reference;
  std::printf("as built: floor %.0f m2, sun %.0f m2h, heating %.4f kWh/m2/day, walk %.1f m (shop %.1f, open %.1f), mix %.4f nats; open %.3f, homes %.3f, shops %.3f, plot ratio %.2f\n",
              R.floor_area, R.sun, R.heat, R.access, R.access_shop, R.access_open, R.mix, R.open_share, R.dwelling_share, R.shop_share, R.plot_ratio);
  std::printf("prepared in %.1f s\n", seconds_since(T0));
  std::fflush(stdout);

  calx::scheduler::ThreadPool pool(threads);
  // The cost of one evaluation, and the bounds of a random design.
  {
    const auto t = Clock::now();
    const hs::Draws draws(99);
    std::vector<Design> ds;
    for (std::size_t i = 0; i < 400; i++) {
      auto s = draws(hs::tag("study.timing"), 0, i);
      ds.push_back(ops.sample(s));
    }
    std::vector<hs::Evaluation> e(ds.size());
    calx::scheduler::parallel_for(pool, ds.size(), 1, [&](std::size_t b, std::size_t en, std::size_t) {
      for (std::size_t i = b; i < en; i++) e[i] = P.evaluate(ds[i]);
    });
    std::size_t feasible = 0, failed = 0;
    for (const auto& x : e) feasible += x.feasible(), failed += x.failed();
    std::printf("400 random designs: %.2f ms each on %zu threads; %zu feasible, %zu failed; per evaluation:", 1e3 * seconds_since(t) / 400, threads, feasible, failed);
    const char* stage[] = {"massing", "heating", "sun", "walking", "mix"};
    for (std::size_t k = 0; k < 5; k++) std::printf(" %s %.2f ms", stage[k], static_cast<double>(P.stage_ns[k].load()) / 1e6 / 401);
    std::printf("\n");
    std::fflush(stdout);
  }
  if (profile) return 0;

  const std::size_t mu = 210, generations = static_cast<std::size_t>(std::max(5.0, std::round(210 * scale))), every = 5;
  const auto references = mo::das_dennis(5, 6);
  std::vector<RunOutcome> runs;
  const std::vector<std::uint64_t> seeds{1, 2, 3};
  const auto objective = [&P](const Design& d) { return P.evaluate(d); };

  for (auto seed : seeds) {
    mo::Nsga3<Design> alg(ops.multiobjective(), references, mu, hs::Draws(seed));
    runs.push_back(run_generational<mo::Nsga3<Design>>("NSGA-III", seed, alg, objective, mu, generations, every, pool, [](const mo::Nsga3<Design>& a) { return a.population(); }));
    std::printf("NSGA-III seed %llu: %.0f evaluations in %.0f s, front %zu\n", static_cast<unsigned long long>(seed), runs.back().evaluations, runs.back().seconds, runs.back().front.size());
    std::fflush(stdout);
  }
  for (auto seed : seeds) {
    mo::Moead<Design>::Options o;
    o.neighbours = 20, o.replacements = 2, o.scalarisation = mo::Scalarisation::tchebycheff;
    mo::Moead<Design> alg(ops.multiobjective(), references, o, hs::Draws(seed));
    runs.push_back(run_generational<mo::Moead<Design>>("MOEA/D", seed, alg, objective, mu, generations, every, pool, [](const mo::Moead<Design>& a) { return a.solutions(); }));
    std::printf("MOEA/D seed %llu: %.0f evaluations in %.0f s, front %zu\n", static_cast<unsigned long long>(seed), runs.back().evaluations, runs.back().seconds, runs.back().front.size());
    std::fflush(stdout);
  }
  for (auto seed : seeds) {
    hs::RandomSearch<Design> alg([&ops](calx::random::Stream& s) { return ops.sample(s); }, hs::Draws(seed));
    RunOutcome r;
    r.algorithm = "random", r.seed = seed;
    const auto t0 = Clock::now();
    std::vector<Individual> archive;
    for (std::size_t g = 0; g < generations; g += every) {
      hs::Record<Design> rec;
      const std::size_t chunk = std::min(every, generations - g);
      const auto o = hs::run(alg, objective, hs::Budget{static_cast<double>(mu * chunk)}, hs::Pool(pool), mu, &rec);
      r.evaluations += o.evaluations;
      merge(archive, rec.told);
      r.checkpoints.push_back({"random", seed, r.evaluations, seconds_since(t0), points_of(archive)});
    }
    r.front = archive, r.seconds = seconds_since(t0);
    std::printf("random seed %llu: %.0f evaluations in %.0f s, archive %zu\n", static_cast<unsigned long long>(seed), r.evaluations, r.seconds, r.front.size());
    std::fflush(stdout);
    runs.push_back(std::move(r));
  }

  // Normalisation for the hypervolume: the ideal and the nadir of the union of the final fronts, the reference 1.1.
  std::vector<double> lo(5, std::numeric_limits<double>::infinity()), hi(5, -std::numeric_limits<double>::infinity());
  {
    std::vector<Individual> all;
    for (const auto& r : runs) all.insert(all.end(), r.front.begin(), r.front.end());
    for (const auto& x : feasible_front(all))
      for (std::size_t i = 0; i < 5; i++) lo[i] = std::min(lo[i], x.evaluation.objectives[i]), hi[i] = std::max(hi[i], x.evaluation.objectives[i]);
  }
  const auto normalise = [&](const std::vector<mo::Point>& pts) {
    std::vector<mo::Point> o;
    for (const auto& p : pts) {
      mo::Point q(5);
      for (std::size_t i = 0; i < 5; i++) q[i] = (p[i] - lo[i]) / (hi[i] - lo[i]);
      o.push_back(q);
    }
    return o;
  };
  const std::vector<double> ref(5, 1.1);
  std::ostringstream hvcsv;
  hvcsv << "algorithm,seed,evaluations,seconds,front,hypervolume\n";
  std::map<std::string, std::vector<double>> final_hv;
  {
    const auto t = Clock::now();
    for (const auto& r : runs)
      for (std::size_t c = 0; c < r.checkpoints.size(); c++) {
        const auto& cp = r.checkpoints[c];
        const bool last = c + 1 == r.checkpoints.size();
        const double v = mo::hypervolume(normalise(cp.front), ref);
        hvcsv << cp.algorithm << "," << cp.seed << "," << cp.evaluations << "," << cp.seconds << "," << cp.front.size() << "," << v << "\n";
        if (last) final_hv[r.algorithm].push_back(v);
      }
    std::printf("hypervolumes in %.0f s\n", seconds_since(t));
  }
  write(out + "/hypervolume.csv", hvcsv.str());

  // NSGA-II on two objectives: floor area and winter sun.
  std::vector<RunOutcome> two;
  const auto objective_two = [&P](const Design& d) { return P.evaluate_two(d); };
  const std::size_t mu2 = 100;
  for (auto seed : seeds) {
    mo::Nsga2<Design> alg(ops.multiobjective(), mu2, hs::Draws(100 + seed));
    two.push_back(run_generational<mo::Nsga2<Design>>("NSGA-II", seed, alg, objective_two, mu2, generations, every, pool, [](const mo::Nsga2<Design>& a) { return a.population(); }));
    std::printf("NSGA-II seed %llu: %.0f evaluations in %.0f s, front %zu\n", static_cast<unsigned long long>(seed), two.back().evaluations, two.back().seconds, two.back().front.size());
    std::fflush(stdout);
  }

  // MAP-Elites over plot ratio and open share; winter sun the quality.
  const double pr_max = 32, open_min = 0.05, open_max = 0.65;
  const std::size_t res_pr = 16, res_open = 12;
  dv::MapElites<Design, dv::Grid>::Options meo;
  meo.initial = 1000, meo.variation_probability = 0.5, meo.breadth = 64;
  dv::MapElites<Design, dv::Grid> me(ops.repertoire(), dv::Grid({0.0, open_min}, {pr_max, open_max}, {res_pr, res_open}), meo, hs::Draws(7));
  const auto objective_me = [&P](const Design& d) { return P.evaluate_repertoire(d); };
  const double me_budget = std::max(2000.0, std::round(84000 * scale));
  std::ostringstream mecsv;
  mecsv << "evaluations,seconds,coverage,qd_score,best\n";
  {
    const auto t0 = Clock::now();
    double spent = 0;
    while (spent < me_budget) {
      const double chunk = std::min(2000.0, me_budget - spent);
      const auto o = hs::run(me, objective_me, hs::Budget{chunk}, hs::Pool(pool), 64);
      spent += o.evaluations;
      double best = std::numeric_limits<double>::infinity();
      for (auto c : me.archive().filled())
        if (me.archive().elite(c)->evaluation.feasible()) best = std::min(best, me.archive().elite(c)->evaluation.objectives[0]);
      mecsv << spent << "," << seconds_since(t0) << "," << me.archive().coverage() << "," << me.archive().qd_score(0.0) << "," << best << "\n";
    }
    std::printf("MAP-Elites: %.0f evaluations in %.0f s, coverage %.3f, QD-score %.2f\n", spent, seconds_since(t0), me.archive().coverage(), me.archive().qd_score(0.0));
  }
  write(out + "/map-elites-progress.csv", mecsv.str());

  // ---- the results ------------------------------------------------------------------------------------------------
  // Fronts: every algorithm's final front, with its measures.
  const auto first = [&](const std::string& name) -> const RunOutcome& {
    for (const auto& r : runs)
      if (r.algorithm == name && r.seed == 1) return r;
    throw std::logic_error("no run " + name);
  };
  std::ostringstream fcsv;
  fcsv << "algorithm,seed,index,floor_area_m2,winter_sun_m2h,heating_kwh_per_m2_day,walk_m,walk_to_shop_m,walk_to_open_space_m,mix_nats,open_share,dwelling_share,shop_share,plot_ratio,links,f1,f2,f3,f4,f5,n1,n2,n3,n4,n5\n";
  const auto front_rows = [&](const std::string& alg, std::uint64_t seed, const std::vector<Individual>& front) {
    for (std::size_t i = 0; i < front.size(); i++) {
      const auto m = P.measure(front[i].genotype);
      const auto& f = front[i].evaluation.objectives;
      char line[1024];
      std::snprintf(line, sizeof line, "%s,%llu,%zu,%.1f,%.1f,%.5f,%.2f,%.2f,%.2f,%.5f,%.5f,%.5f,%.5f,%.4f,%zu", alg.c_str(), static_cast<unsigned long long>(seed), i, m.floor_area, m.sun, m.heat,
                    m.access, m.access_shop, m.access_open, m.mix, m.open_share, m.dwelling_share, m.shop_share, m.plot_ratio, m.links);
      fcsv << line;
      for (std::size_t k = 0; k < f.size(); k++) fcsv << "," << f[k];
      for (std::size_t k = f.size(); k < 5; k++) fcsv << ",";
      if (f.size() == 5)
        for (std::size_t k = 0; k < 5; k++) fcsv << "," << (f[k] - lo[k]) / (hi[k] - lo[k]);
      else
        fcsv << ",,,,,";
      fcsv << "\n";
    }
  };
  for (const auto& r : runs) front_rows(r.algorithm, r.seed, r.front);
  for (const auto& r : two) front_rows(r.algorithm, r.seed, r.front);
  write(out + "/fronts.csv", fcsv.str());

  // The chosen designs, from NSGA-III's first front: the extreme of each measure, and the balanced design (the least
  // largest normalised objective, the Tchebycheff point nearest the ideal).
  const auto& N3 = first("NSGA-III").front;
  struct Chosen { std::string key, title; std::size_t index; };
  std::vector<Chosen> chosen;
  const char* keys[] = {"most-floor", "most-sun", "least-heating", "shortest-walks", "most-mixed"};
  const char* titles[] = {"the most floor", "the most winter sun", "the least heating", "the shortest walks", "the most mixed"};
  for (std::size_t k = 0; k < 5; k++) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < N3.size(); i++)
      if (N3[i].evaluation.objectives[k] < N3[best].evaluation.objectives[k]) best = i;
    chosen.push_back({keys[k], titles[k], best});
  }
  {
    std::size_t best = 0;
    double bv = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < N3.size(); i++) {
      double v = -std::numeric_limits<double>::infinity();
      for (std::size_t k = 0; k < 5; k++) v = std::max(v, (N3[i].evaluation.objectives[k] - lo[k]) / (hi[k] - lo[k]));
      if (v < bv) bv = v, best = i;
    }
    chosen.push_back({"balanced", "a balanced design", best});
  }

  // designs.json: keyed by the city model's building ids and CLUE's base properties, for the workbench's GIS.
  {
    Json parcels = Json::array();
    for (const auto& p : P.parcels)
      parcels.push(Json::object({{"building", p.id}, {"clue_property", p.clue_property.empty() ? Json() : Json(p.clue_property)}, {"name", p.name}, {"footprint_m2", std::round(p.area * 10) / 10},
                                 {"storeys_bound", static_cast<double>(p.bound)}, {"as_built", Json::object({{"storeys", static_cast<double>(p.floors0)}, {"use", use_name(p.use0)}, {"clue_use", p.clue_use}, {"height_source", p.height_source}})}}));
    Json links = Json::array();
    for (const auto& q : P.links)
      links.push(Json::object({{"kind", q.kind}, {"from", numbers({q.a[0], q.a[1]}, 2)}, {"to", numbers({q.b[0], q.b[1]}, 2)}, {"length_m", std::round(q.length * 10) / 10},
                               {"joins", Json::array({q.from, q.to})}}));
    Json designs = Json::array();
    const auto one = [&](const std::string& key, const std::string& title, const std::string& from, const Design& d) {
      const auto m = P.measure(d);
      Json by = Json::object();
      for (std::size_t i = 0; i < n; i++) {
        const bool open = d.uses[i] == open_space;
        by.set(P.parcels[i].id, Json::object({{"storeys", open ? 0.0 : static_cast<double>(d.storeys[i])}, {"height_m", open ? 0.0 : std::round(P.A.storey * static_cast<double>(d.storeys[i]) * 10) / 10},
                                               {"use", use_name(d.uses[i])}}));
      }
      Json opened = Json::array();
      for (std::size_t k = 0; k < K; k++)
        if (d.links[k]) opened.push(Json(static_cast<double>(k)));
      designs.push(Json::object({{"key", key}, {"title", title}, {"from", from}, {"measures", measures_json(m)}, {"objectives_relative_to_as_built", numbers(P.objectives(m), 5)},
                                 {"parcels", by}, {"links_made", opened}}));
    };
    one("as-built", "the district as built", "CLUE's latest census (floors, predominant use); no link made", P.existing());
    for (const auto& c : chosen) one(c.key, c.title, "NSGA-III, seed 1, final front, member " + std::to_string(c.index), N3[c.index].genotype);
    Json doc = Json::object({{"study", "research/experiments/2026-11-05-designing-a-district"},
                             {"about", "Designs of the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street (Melbourne's Hoddle Grid), keyed by the city model's building ids (OpenStreetMap ways) and CLUE's base properties, for drawing on the Urban Science Workbench's GIS. A study of methods on a real place, not a planning proposal."},
                             {"frame", J["frame"]}, {"origin", J["origin"]}, {"storey_m", P.A.storey},
                             {"conventions", "storeys 0 and use 'open space': the parcel's building removed and its ground given to the public; links: a blind lane carried through to the street beyond it, or a passage through a block, each an arcade at ground level between the points 'from' and 'to' (local metres), taking 4 m of ground floor along its length; links_made lists the indices of the links made; heights = storeys x storey_m"},
                             {"parcels", parcels}, {"links", links}, {"designs", designs}});
    write(out + "/designs.json", doc.dump());
  }

  // The fronts' designs (NSGA-III seed 1; NSGA-II seed 1), for the figures.
  {
    Json a = Json::array();
    for (const auto& x : N3) a.push(Json::object({{"design", design_json(x.genotype)}, {"objectives", numbers(x.evaluation.objectives, 6)}}));
    Json b = Json::array();
    for (const auto& r : two)
      if (r.seed == 1)
        for (const auto& x : r.front) b.push(Json::object({{"design", design_json(x.genotype)}, {"objectives", numbers(x.evaluation.objectives, 6)}}));
    Json ch = Json::array();
    for (const auto& c : chosen) ch.push(Json::object({{"key", c.key}, {"title", c.title}, {"index", static_cast<double>(c.index)}}));
    write(out + "/front-designs.json", Json::object({{"nsga3", a}, {"nsga2", b}, {"chosen", ch}}).dump());
  }

  // The repertoire.
  {
    Json cells = Json::array();
    const auto& ar = me.archive();
    for (auto c : ar.filled()) {
      const auto& x = *ar.elite(c);
      const auto m = P.measure(x.genotype);
      cells.push(Json::object({{"cell", static_cast<double>(c)}, {"i", static_cast<double>(c % res_pr)}, {"j", static_cast<double>(c / res_pr)}, {"feasible", x.evaluation.feasible()},
                               {"quality", x.evaluation.objectives[0]}, {"descriptors", numbers(x.evaluation.descriptors, 5)}, {"measures", measures_json(m)}, {"design", design_json(x.genotype)}}));
    }
    write(out + "/map-elites.json", Json::object({{"grid", Json::object({{"plot_ratio", numbers({0, pr_max}, 3)}, {"open_share", numbers({open_min, open_max}, 3)}, {"resolution", Json::array({static_cast<double>(res_pr), static_cast<double>(res_open)})}})},
                                                  {"quality", "winter sun on the public ground relative to the district as built, negated (minimised)"},
                                                  {"constraints", "homes at least 20% of the floor, shops between 5% and 20% of it; open share is a descriptor here, not a constraint"},
                                                  {"budget", me_budget}, {"coverage", ar.coverage()}, {"qd_score", ar.qd_score(0.0)}, {"cells", cells}})
                                       .dump());
  }

  // Sun maps of the as-built district and the chosen designs; the district's cells and links, for drawing.
  {
    Json pts = Json::array(), owner = Json::array();
    for (std::size_t p = 0; p < P.shade.points.size(); p++) pts.push(numbers({P.shade.points[p][0], P.shade.points[p][1]}, 2)), owner.push(Json(static_cast<double>(P.shade.owner[p])));
    Json maps = Json::object();
    maps.set("as-built", numbers(P.sun_map(P.existing()), 2));
    for (const auto& c : chosen) maps.set(c.key, numbers(P.sun_map(N3[c.index].genotype), 2));
    Json suns = Json::array();
    for (const auto& s : P.shade.suns) suns.push(numbers({s.altitude, s.azimuth, s.weight}, 6));
    write(out + "/sun-maps.json", Json::object({{"spacing_m", P.A.sample}, {"points", pts}, {"owner", owner}, {"suns", suns}, {"hours", maps}}).dump());

    Json blocks = Json::array(), holes = Json::array(), segments = Json::array();
    for (const auto id : P.base.of(2)) blocks.push(Json::object({{"id", std::to_string(id)}, {"ring", ring_json(P.base.ring(id))}, {"area_m2", std::round(P.base.area(id))}}));
    for (const auto& r : P.base.hole_rings) holes.push(ring_json(r));
    for (const auto id : P.base.of(1)) segments.push(Json::object({{"id", std::to_string(id)}, {"type", P.base.cell(id).type}, {"name", P.base.name.at(id)}, {"line", ring_json(P.base.line.at(id))}}));
    Json order = Json::array(), bounds = Json::array(), place = Json::array();
    for (const auto& p : P.parcels) order.push(p.id), bounds.push(Json(static_cast<double>(p.bound))), place.push(Json(static_cast<double>(p.place)));
    Json links = Json::array();
    for (const auto& q : P.links) links.push(Json::object({{"kind", q.kind}, {"a", numbers({q.a[0], q.a[1]}, 2)}, {"b", numbers({q.b[0], q.b[1]}, 2)}, {"length_m", q.length}, {"joins", Json::array({q.from, q.to})}}));
    Json dropped = Json::array();
    for (const auto& d : P.dropped) dropped.push(d);
    const auto topo = morphe::topology::admissible(P.base.complex);
    write(out + "/district-cells.json", Json::object({{"parcel_order", order}, {"storeys_bound", bounds}, {"place", place}, {"places", static_cast<double>(P.place_count)},
                                                      {"blocks", blocks}, {"holes", holes}, {"segments", segments}, {"links", links}, {"dropped", dropped},
                                                      {"regular", topo.regular}, {"cells", Json::array({static_cast<double>(P.base.of(0).size()), static_cast<double>(P.base.of(1).size()), static_cast<double>(P.base.of(2).size())})},
                                                      {"land_m2", P.land}, {"open_space_m2", P.open_space_area}})
                                           .dump());
  }

  // The summary.
  {
    Json hvj = Json::object();
    for (const auto& [alg, v] : final_hv) hvj.set(alg, numbers(v, 6));
    Json runs_j = Json::array();
    for (const auto& r : runs) runs_j.push(Json::object({{"algorithm", r.algorithm}, {"seed", static_cast<double>(r.seed)}, {"evaluations", r.evaluations}, {"seconds", std::round(r.seconds)}, {"front", static_cast<double>(r.front.size())}}));
    for (const auto& r : two) runs_j.push(Json::object({{"algorithm", r.algorithm}, {"seed", static_cast<double>(r.seed)}, {"evaluations", r.evaluations}, {"seconds", std::round(r.seconds)}, {"front", static_cast<double>(r.front.size())}}));
    // Two-objective hypervolumes: NSGA-II's fronts, and NSGA-III's fronts projected on the same two objectives.
    std::vector<mo::Point> all2;
    for (const auto& r : two) for (const auto& p : points_of(r.front)) all2.push_back(p);
    for (const auto& r : runs) if (r.algorithm == "NSGA-III") for (const auto& p : points_of(r.front, 2)) all2.push_back(p);
    std::vector<double> lo2(2, 1e300), hi2(2, -1e300);
    for (auto i : mo::non_dominated(all2))
      for (std::size_t k = 0; k < 2; k++) lo2[k] = std::min(lo2[k], all2[i][k]), hi2[k] = std::max(hi2[k], all2[i][k]);
    const auto hv2 = [&](const std::vector<mo::Point>& pts) {
      std::vector<mo::Point> q;
      for (const auto& p : pts) q.push_back({(p[0] - lo2[0]) / (hi2[0] - lo2[0]), (p[1] - lo2[1]) / (hi2[1] - lo2[1])});
      return mo::hypervolume_2d(q, std::vector<double>{1.1, 1.1});
    };
    Json a2 = Json::array(), b2 = Json::array();
    for (const auto& r : two) a2.push(hv2(points_of(r.front)));
    for (const auto& r : runs) if (r.algorithm == "NSGA-III") b2.push(hv2(points_of(r.front, 2)));
    Json chosen_j = Json::array();
    for (const auto& c : chosen) {
      const auto m = P.measure(N3[c.index].genotype);
      chosen_j.push(Json::object({{"key", c.key}, {"title", c.title}, {"measures", measures_json(m)}, {"objectives", numbers(N3[c.index].evaluation.objectives, 5)}}));
    }
    Json parcels_j = Json::object({{"count", static_cast<double>(n)}, {"joined_to_clue", 0.0}});
    double joined = 0;
    for (const auto& p : P.parcels) joined += !p.clue_property.empty();
    parcels_j.set("joined_to_clue", joined);
    std::array<double, 4> u0{};
    for (const auto& p : P.parcels) u0[static_cast<std::size_t>(p.use0)] += 1;
    Json summary = Json::object({
        {"parcels", parcels_j}, {"uses_as_built", numbers({u0[0], u0[1], u0[2], u0[3]}, 0)}, {"candidate_links", static_cast<double>(K)},
        {"ground_samples", static_cast<double>(P.shade.points.size())}, {"land_m2", std::round(P.land)}, {"existing_open_space_m2", std::round(P.open_space_area)},
        {"as_built", measures_json(R)}, {"normalisation", Json::object({{"low", numbers(lo, 6)}, {"high", numbers(hi, 6)}, {"reference", numbers(ref, 3)}})},
        {"hypervolume_final", hvj}, {"hypervolume_two_objectives", Json::object({{"NSGA-II", a2}, {"NSGA-III projected", b2}, {"low", numbers(lo2, 6)}, {"high", numbers(hi2, 6)}})},
        {"runs", runs_j}, {"chosen", chosen_j}, {"threads", static_cast<double>(threads)}, {"generations", static_cast<double>(generations)}, {"population", static_cast<double>(mu)},
        {"population_two_objectives", static_cast<double>(mu2)}, {"map_elites_budget", me_budget}, {"seconds", std::round(seconds_since(T0))}});
    write(out + "/summary.json", summary.dump());
  }
  std::printf("done in %.0f s\n", seconds_since(T0));
  return 0;
}
