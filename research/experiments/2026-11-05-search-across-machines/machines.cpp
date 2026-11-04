// Search across machines: how information crosses a structured population, what synchronisation costs on
// heterogeneous workers, and how an archipelago's topology and migration interval decide a deceptive search.
//
//   machines <results/> [part: all | takeover | synchronisation | islands] [threads] [pilot]
//
// (a) Takeover. A cellular population on a graph under selection alone (heuresis::populations::CellularEA with no
//     variation): one best individual at cell 0, the rest of the cells holding distinct worse values. With each
//     cell taking the best of its closed neighbourhood, cell i holds after g generations the best of its ball of
//     radius g, so the takeover time is the eccentricity of cell 0, exactly; with binary tournaments in the
//     neighbourhood the wave is slower and random. Eccentricities by breadth-first search (morphe::combinatorics).
// (b) Synchronisation. 32 simulated workers (heuresis::search::run_asynchronous, the discrete-event master–worker
//     model) evaluate a real-coded GA on the 10-dimensional Rastrigin function. An evaluation lasts a lognormal
//     time of mean 1 and log-spread σ, drawn from the stream keyed by its ticket; in the straggler case worker 0
//     takes ten times as long. Compared: the generational GA (it waits for each generation's slowest
//     evaluation), the steady-state GA (it keeps every worker busy), and 32 asynchronous islands on machines of
//     their own (heuresis::populations::run_islands_asynchronous). The protocol's identity is shown alongside:
//     the generational GA's states, generation by generation, on the serial driver, a thread pool, generations
//     cut into parts, and the workers at every σ.
// (c) Islands. 16 islands of a GA on concatenated deceptive traps of order 5 (synchronous migration,
//     heuresis::populations::run_islands), over five topologies and six migration intervals, and one panmictic
//     GA of the same total population; the success rate (the optimum found within the budget) with Clopper and
//     Pearson's interval (heuresis::analysis).

#include <heuresis/analysis/runtimes.hpp>
#include <heuresis/benchmarks/benchmarks.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/populations/populations.hpp>
#include <heuresis/search/search.hpp>

#include <morphe/combinatorics/graph.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using heuresis::populations::CellSelection;
using heuresis::populations::Topology;
using heuresis::search::Budget;
using heuresis::search::Draws;
using heuresis::search::Evaluation;
using heuresis::search::tag;
using heuresis::variation::Bits;
using heuresis::variation::Real;

constexpr double inf = std::numeric_limits<double>::infinity();

// ---- (a) takeover ------------------------------------------------------------------------------------------

std::vector<std::size_t> distances_from(const Topology& t, std::size_t s) {
  morphe::combinatorics::Multigraph G;
  G.vertices = t.size();
  for (std::size_t i = 0; i < t.size(); i++)
    for (auto j : t.out[i])
      if (i < j) G.edges.push_back({i, j});
  return morphe::combinatorics::distances(G, s);
}

struct Wave {
  std::size_t takeover = 0;               // generations until every cell holds the best
  std::vector<std::vector<int>> values;   // the cells' values after each generation, cells by distance from 0
};

constexpr std::uint32_t initial_tag = tag("experiment.machines.initial");

Wave takeover(const Topology& t, CellSelection selection, std::uint64_t seed, bool keep, std::size_t cap = 100000) {
  const std::size_t N = t.size();
  const Draws draws(seed);
  calx::random::Stream s = draws(initial_tag, 0, 0);
  auto p = heuresis::variation::random_permutation(N, s);
  std::swap(p[0], p[static_cast<std::size_t>(std::find(p.begin(), p.end(), 0u) - p.begin())]);
  std::vector<double> init(p.begin(), p.end());
  heuresis::populations::CellularEA<double> ea(init, t, selection, {}, {}, 1, draws);
  const auto d = distances_from(t, 0);
  std::vector<std::size_t> order(N);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return d[a] < d[b]; });
  auto f = [](const double& g) { return g; };
  Wave w;
  for (std::size_t g = 0; g <= cap; g++) {
    auto c = ea.ask(heuresis::search::unlimited);
    auto e = heuresis::search::Serial{}(std::span<const heuresis::search::Candidate<double>>(c), f);
    std::vector<heuresis::search::Told> told(c.size());
    for (std::size_t i = 0; i < c.size(); i++) told[i] = {c[i].ticket, e[i]};
    ea.tell(told);
    const auto& cells = ea.cells();
    if (keep) {
      std::vector<int> row(N);
      for (std::size_t r = 0; r < N; r++) row[r] = static_cast<int>(cells[order[r]].genotype);
      w.values.push_back(std::move(row));
    }
    if (std::all_of(cells.begin(), cells.end(), [](const auto& x) { return x.genotype == 0; })) {
      w.takeover = g;  // generation 0 is the initial population
      return w;
    }
  }
  throw std::runtime_error("takeover not reached");
}

void part_takeover(const std::string& out) {
  struct Named {
    std::string name, kind;
    Topology t;
  };
  // Heat maps: four graphs of 256 cells.
  std::vector<Named> maps{{"ring of 256", "ring", Topology::ring(256)}, {"torus 16 × 16", "torus", Topology::torus(16, 16)}, {"hypercube of dimension 8", "hypercube", Topology::hypercube(8)}, {"complete graph of 256", "complete", Topology::complete(256)}};
  {
    std::ofstream f(out + "/takeover_waves.json");
    f << "{\n";
    for (std::size_t m = 0; m < maps.size(); m++) {
      const auto d = distances_from(maps[m].t, 0);
      std::vector<std::size_t> sorted(d.begin(), d.end());
      std::sort(sorted.begin(), sorted.end());
      f << "  \"" << maps[m].kind << "\": {\"name\": \"" << maps[m].name << "\", \"eccentricity\": " << sorted.back() << ", \"distance\": [";
      for (std::size_t i = 0; i < sorted.size(); i++) f << (i ? "," : "") << sorted[i];
      f << "]";
      for (auto sel : {CellSelection::best, CellSelection::tournament}) {
        const Wave w = takeover(maps[m].t, sel, 1, true);
        f << ", \"" << (sel == CellSelection::best ? "best" : "tournament") << "\": {\"takeover\": " << w.takeover << ", \"values\": [";
        for (std::size_t g = 0; g < w.values.size(); g++) {
          f << (g ? "," : "") << "[";
          for (std::size_t r = 0; r < w.values[g].size(); r++) f << (r ? "," : "") << w.values[g][r];
          f << "]";
        }
        f << "]}";
        std::cerr << maps[m].name << (sel == CellSelection::best ? ", best: " : ", tournament: ") << w.takeover << " generations (eccentricity " << sorted.back() << ")\n";
      }
      f << "}" << (m + 1 < maps.size() ? "," : "") << "\n";
    }
    f << "}\n";
  }
  // Takeover times against the eccentricity.
  std::vector<Named> graphs;
  for (std::size_t k : {16, 32, 64, 128, 256}) graphs.push_back({"ring " + std::to_string(k), "ring", Topology::ring(k)});
  for (auto [r, c] : std::vector<std::pair<std::size_t, std::size_t>>{{4, 4}, {6, 9}, {8, 8}, {8, 32}, {16, 16}, {32, 32}}) graphs.push_back({"torus " + std::to_string(r) + "x" + std::to_string(c), "torus", Topology::torus(r, c)});
  for (auto [r, c] : std::vector<std::pair<std::size_t, std::size_t>>{{8, 8}, {16, 16}, {32, 32}}) graphs.push_back({"Moore torus " + std::to_string(r) + "x" + std::to_string(c), "torus (Moore)", Topology::torus(r, c, true)});
  for (std::size_t dim : {4, 5, 6, 7, 8, 9, 10}) graphs.push_back({"hypercube " + std::to_string(dim), "hypercube", Topology::hypercube(dim)});
  for (std::size_t k : {16, 64, 256, 1024}) graphs.push_back({"complete " + std::to_string(k), "complete", Topology::complete(k)});
  std::ofstream f(out + "/takeover_times.csv");
  f << "graph,kind,cells,eccentricity,selection,seed,takeover\n";
  std::size_t exact = 0;
  for (const auto& g : graphs) {
    const auto d = distances_from(g.t, 0);
    const std::size_t ecc = *std::max_element(d.begin(), d.end());
    const Wave b = takeover(g.t, CellSelection::best, 0, false);
    exact += b.takeover == ecc;
    f << g.name << ',' << g.kind << ',' << g.t.size() << ',' << ecc << ",best,0," << b.takeover << '\n';
    double mean = 0;
    for (std::uint64_t s = 0; s < 20; s++) {
      const Wave w = takeover(g.t, CellSelection::tournament, s, false);
      f << g.name << ',' << g.kind << ',' << g.t.size() << ',' << ecc << ",tournament," << s << ',' << w.takeover << '\n';
      mean += static_cast<double>(w.takeover) / 20;
    }
    std::fprintf(stderr, "  %-20s ecc %3zu  best %3zu  tournament %.1f\n", g.name.c_str(), ecc, b.takeover, mean);
  }
  std::cerr << "best-of-neighbourhood takeover equal to the eccentricity on " << exact << " of " << graphs.size() << " graphs\n";
}

// ---- (b) synchronisation -----------------------------------------------------------------------------------

constexpr std::size_t W = 32;           // workers (and islands)
constexpr std::size_t D = 20;           // Rastrigin's dimension
constexpr double B = 32000;             // evaluations

heuresis::genetic::Operators<Real> real_operators() {
  const auto box = heuresis::variation::Box::cube(D, -5.12, 5.12);
  return {heuresis::genetic::uniform_point(box), heuresis::genetic::polynomial(box, 20, 1.0 / D), heuresis::genetic::simulated_binary(box, 15)};
}

double rastrigin(const Real& x) { return heuresis::benchmarks::rastrigin(x); }

struct Spread {
  std::string name;
  double sigma = 0;
  bool straggler = false;
};

// A lognormal duration of mean 1 (log-spread σ); worker 0 ten times slower in the straggler case.
double duration_of(const Spread& sp, std::size_t worker, calx::random::Stream& s) {
  double d = sp.sigma > 0 ? std::exp(sp.sigma * s.normal() - sp.sigma * sp.sigma / 2) : 1.0;
  if (sp.straggler && worker == 0) d *= 10;
  return d;
}

heuresis::genetic::GaOptions ga_options(std::size_t mu) {
  heuresis::genetic::GaOptions o;
  o.mu = mu;
  o.lambda = mu;
  o.crossover_rate = 0.9;
  o.tournament = 2;
  o.elites = 1;
  return o;
}

struct Timed {
  double time = 0, best = inf, evaluations = 0;
};

Timed run_scheme(const std::string& scheme, const Spread& sp, std::uint64_t seed, const Budget& budget) {
  using namespace heuresis;
  Timed r;
  if (scheme == "islands") {
    std::vector<genetic::GeneticAlgorithm<Real>> islands;
    for (std::size_t i = 0; i < W; i++) islands.emplace_back(real_operators(), ga_options(8), Draws(seed, static_cast<std::uint32_t>(i)));
    populations::IslandClock clock;
    clock.duration = [sp](const Evaluation&, std::size_t island, calx::random::Stream& s) { return duration_of(sp, island, s); };
    clock.seed = seed;
    const auto o = populations::run_islands_asynchronous(islands, Topology::ring(W), populations::Migration{5, 1}, rastrigin, budget, clock);
    r.time = o.time;
    r.best = o.best->evaluation.value();
    r.evaluations = o.evaluations;
    return r;
  }
  search::Workers workers;
  workers.count = W;
  workers.duration = [sp](const Evaluation&, std::size_t worker, calx::random::Stream& s) { return duration_of(sp, worker, s); };
  workers.seed = seed;
  if (scheme == "generational") {
    genetic::GeneticAlgorithm<Real> a(real_operators(), ga_options(W), Draws(seed));
    const auto o = search::run_asynchronous(a, rastrigin, budget, workers);
    r.time = o.time;
    r.best = o.best->evaluation.value();
    r.evaluations = o.evaluations;
    return r;
  }
  if (scheme == "steady-state") {
    genetic::SteadyStateGA<Real> a(real_operators(), ga_options(W), Draws(seed));
    const auto o = search::run_asynchronous(a, rastrigin, budget, workers);
    r.time = o.time;
    r.best = o.best->evaluation.value();
    r.evaluations = o.evaluations;
    return r;
  }
  throw std::invalid_argument("unknown scheme " + scheme);
}

// An optimiser that watches another: after every tell, the best of its population, and its generation.
template <class A>
class Watched {
 public:
  using Genotype = typename A::Genotype;
  explicit Watched(A& a) : a_(&a) {}
  std::vector<heuresis::search::Candidate<Genotype>> ask(std::size_t limit) { return a_->ask(limit); }
  void tell(std::span<const heuresis::search::Told> told) {
    a_->tell(told);
    tells_ += told.size();
    const auto& pop = a_->population();
    if (pop.empty()) return;
    double b = inf;
    for (const auto& x : pop) b = std::min(b, x.evaluation.value());
    if constexpr (requires { a_->generation(); }) {
      if (a_->generation() != generation_) {  // one value per generation
        generation_ = a_->generation();
        trace_.push_back(b);
      }
    } else {
      trace_.resize(tells_, b);  // one value per tell: trace_[i] after i + 1 tells
    }
  }
  bool stopped() const { return a_->stopped(); }
  const std::vector<double>& trace() const { return trace_; }

 private:
  A* a_;
  std::vector<double> trace_;
  std::uint64_t generation_ = 0;
  std::size_t tells_ = 0;
};

void part_synchronisation(const std::string& out, calx::scheduler::ThreadPool& pool, bool pilot) {
  std::vector<Spread> spreads;
  for (double s : {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0}) {
    char name[32];
    std::snprintf(name, sizeof name, "%.2f", s);
    spreads.push_back({name, s, false});
  }
  spreads.push_back({"straggler", 0.0, true});
  const std::vector<std::string> schemes{"generational", "steady-state", "islands"};
  const std::size_t seeds = pilot ? 3 : 10;
  // Two values of Rastrigin's function: a loose one every scheme reaches, a tight one.
  const std::vector<double> targets{10.0, 1.0};
  struct Job {
    std::size_t scheme, spread;
    std::uint64_t seed;
  };
  std::vector<Job> jobs;
  for (std::size_t a = 0; a < schemes.size(); a++)
    for (std::size_t s = 0; s < spreads.size(); s++)
      for (std::uint64_t k = 0; k < seeds; k++) jobs.push_back({a, s, k});
  std::vector<Timed> full(jobs.size());
  std::vector<std::vector<Timed>> to_target(targets.size(), std::vector<Timed>(jobs.size()));
  calx::scheduler::parallel_for(pool, jobs.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t j = b; j < e; j++) {
      const auto& J = jobs[j];
      full[j] = run_scheme(schemes[J.scheme], spreads[J.spread], J.seed, Budget{B});
      for (std::size_t k = 0; k < targets.size(); k++) {
        Budget bt{B};
        bt.target = targets[k];
        to_target[k][j] = run_scheme(schemes[J.scheme], spreads[J.spread], J.seed, bt);
      }
    }
  });
  {
    std::ofstream f(out + "/synchronisation.csv");
    f << "scheme,spread,sigma,straggler,seed,makespan,best,evaluations,time_to_10,evaluations_to_10,time_to_1,evaluations_to_1\n";
    f.precision(10);
    for (std::size_t j = 0; j < jobs.size(); j++) {
      const auto& J = jobs[j];
      f << schemes[J.scheme] << ',' << spreads[J.spread].name << ',' << spreads[J.spread].sigma << ',' << (spreads[J.spread].straggler ? 1 : 0) << ',' << J.seed << ',' << full[j].time << ',' << full[j].best << ',' << full[j].evaluations;
      for (std::size_t k = 0; k < targets.size(); k++) {
        const bool reached = to_target[k][j].best <= targets[k];
        f << ',' << (reached ? to_target[k][j].time : inf) << ',' << (reached ? to_target[k][j].evaluations : inf);
      }
      f << '\n';
    }
  }
  for (std::size_t a = 0; a < schemes.size(); a++)
    for (std::size_t s = 0; s < spreads.size(); s++) {
      std::vector<double> t, q;
      std::size_t hit = 0;
      for (std::size_t j = 0; j < jobs.size(); j++)
        if (jobs[j].scheme == a && jobs[j].spread == s) {
          t.push_back(full[j].time);
          q.push_back(full[j].best);
          hit += to_target[1][j].best <= targets[1];
        }
      std::sort(t.begin(), t.end());
      std::sort(q.begin(), q.end());
      std::fprintf(stderr, "  %-13s %-9s makespan %8.1f  best (median) %7.3f  reached %zu/%zu\n", schemes[a].c_str(), spreads[s].name.c_str(), t[t.size() / 2], q[q.size() / 2], hit, t.size());
    }

  // The protocol's identity: the generational GA's states on every placement; the steady state's beside it.
  struct Placement {
    std::string name;
    std::function<std::vector<double>(bool)> trace;  // argument: steady state
  };
  auto watched = [&](bool steady, auto&& drive) {
    using namespace heuresis;
    if (steady) {
      genetic::SteadyStateGA<Real> a(real_operators(), ga_options(W), Draws(0));
      Watched<genetic::SteadyStateGA<Real>> w(a);
      drive(w);
      // Thin the steady state's trace to one value per W tells, as a generation's worth.
      std::vector<double> t;
      for (std::size_t i = W - 1; i < w.trace().size(); i += W) t.push_back(w.trace()[i]);
      return t;
    }
    genetic::GeneticAlgorithm<Real> a(real_operators(), ga_options(W), Draws(0));
    Watched<genetic::GeneticAlgorithm<Real>> w(a);
    drive(w);
    return w.trace();
  };
  calx::scheduler::ThreadPool four(4);
  std::vector<Placement> placements;
  placements.push_back({"serial", [&](bool st) { return watched(st, [&](auto& w) { heuresis::search::run(w, rastrigin, Budget{B}); }); }});
  placements.push_back({"pool of 4 threads", [&](bool st) { return watched(st, [&](auto& w) { heuresis::search::run(w, rastrigin, Budget{B}, heuresis::search::Pool(four)); }); }});
  placements.push_back({"parts of 7", [&](bool st) { return watched(st, [&](auto& w) { heuresis::search::run(w, rastrigin, Budget{B}, heuresis::search::Serial{}, 7); }); }});
  for (const Spread& sp : {Spread{"workers at sigma 0.5", 0.5, false}, Spread{"workers at sigma 2", 2.0, false}, Spread{"workers with a straggler", 0.0, true}}) {
    placements.push_back({sp.name, [&, sp](bool st) {
                            return watched(st, [&](auto& w) {
                              heuresis::search::Workers workers;
                              workers.count = W;
                              workers.duration = [sp](const Evaluation&, std::size_t worker, calx::random::Stream& s) { return duration_of(sp, worker, s); };
                              workers.seed = 0;
                              heuresis::search::run_asynchronous(w, rastrigin, Budget{B}, workers);
                            });
                          }});
  }
  // The traces at generations 0–19 and every tenth after (bit for bit, as 17 digits), and the largest
  // difference from the serial run over every generation.
  std::ofstream f(out + "/identity.csv");
  f << "scheme,placement,generation,best\n";
  f.precision(17);
  std::ofstream summary(out + "/identity_summary.csv");
  summary << "scheme,placement,generations,largest_difference_from_serial,identical\n";
  summary.precision(17);
  for (bool steady : {false, true}) {
    const char* scheme = steady ? "steady-state" : "generational";
    std::vector<std::vector<double>> traces;
    for (const auto& p : placements) traces.push_back(p.trace(steady));
    for (std::size_t k = 0; k < traces.size(); k++) {
      const std::size_t len = std::min(traces[k].size(), traces[0].size());
      double worst = 0;
      bool same = traces[k].size() == traces[0].size();
      for (std::size_t g = 0; g < len; g++) {
        worst = std::max(worst, std::fabs(traces[k][g] - traces[0][g]));
        same = same && traces[k][g] == traces[0][g];
      }
      summary << scheme << ',' << placements[k].name << ',' << traces[k].size() << ',' << worst << ',' << (same ? 1 : 0) << '\n';
      std::cerr << "  " << scheme << ", " << placements[k].name << ": " << traces[k].size() << " generations, largest difference from the serial run " << worst << (same ? " (identical)" : "") << '\n';
      for (std::size_t g = 0; g < traces[k].size(); g++)
        if (g < 20 || g % 10 == 0 || g + 1 == traces[k].size()) f << scheme << ',' << placements[k].name << ',' << g << ',' << traces[k][g] << '\n';
    }
  }
}

// ---- (c) islands on deceptive traps ------------------------------------------------------------------------

constexpr std::size_t TRAP_N = 40, TRAP_K = 5, ISLANDS = 16, ISLAND_MU = 20;
constexpr double TRAP_BUDGET = 16 * 20 * 400;  // 400 generations of every island

double trap(const Bits& x) { return heuresis::benchmarks::trap(x, TRAP_K); }

heuresis::genetic::Operators<Bits> bit_operators() {
  return {heuresis::genetic::random_bits(TRAP_N), heuresis::genetic::bit_flip(), heuresis::genetic::k_point_crossover<std::uint8_t>(2)};
}

struct IslandRun {
  bool success = false;
  double evaluations = 0, best = inf;
};

IslandRun run_archipelago(const Topology& t, std::size_t interval, std::uint64_t seed, calx::scheduler::ThreadPool* pool = nullptr) {
  using namespace heuresis;
  std::vector<genetic::GeneticAlgorithm<Bits>> islands;
  for (std::size_t i = 0; i < t.size(); i++) islands.emplace_back(bit_operators(), ga_options(ISLAND_MU), Draws(seed, static_cast<std::uint32_t>(i)));
  Budget b{TRAP_BUDGET};
  b.target = 0;
  const auto o = populations::run_islands(islands, t, populations::Migration{interval, 1}, trap, b, pool);
  return {o.best->evaluation.value() <= 0, o.evaluations, o.best->evaluation.value()};
}

IslandRun run_panmictic(std::uint64_t seed) {
  using namespace heuresis;
  genetic::GeneticAlgorithm<Bits> a(bit_operators(), ga_options(ISLANDS * ISLAND_MU), Draws(seed));
  Budget b{TRAP_BUDGET};
  b.target = 0;
  const auto o = search::run(a, trap, b);
  return {o.best->evaluation.value() <= 0, o.evaluations, o.best->evaluation.value()};
}

void part_islands(const std::string& out, calx::scheduler::ThreadPool& pool, bool pilot) {
  struct Named {
    std::string name;
    Topology t;
  };
  const std::vector<Named> topologies{{"isolated", Topology::isolated(ISLANDS)}, {"ring", Topology::ring(ISLANDS)}, {"torus 4 × 4", Topology::torus(4, 4)}, {"hypercube", Topology::hypercube(4)}, {"complete", Topology::complete(ISLANDS)}};
  const std::vector<std::size_t> intervals{1, 2, 5, 10, 25, 50};
  const std::size_t seeds = pilot ? 8 : 40;
  struct Job {
    int topology;  // −1: the panmictic GA
    std::size_t interval;
    std::uint64_t seed;
  };
  std::vector<Job> jobs;
  for (std::uint64_t s = 0; s < seeds; s++) jobs.push_back({-1, 0, s});
  // Isolated islands do not depend on the interval: they are run once.
  for (std::size_t t = 0; t < topologies.size(); t++)
    for (auto iv : intervals)
      if (t > 0 || iv == intervals[0])
        for (std::uint64_t s = 0; s < seeds; s++) jobs.push_back({static_cast<int>(t), iv, s});
  std::vector<IslandRun> runs(jobs.size());
  calx::scheduler::parallel_for(pool, jobs.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t j = b; j < e; j++) runs[j] = jobs[j].topology < 0 ? run_panmictic(jobs[j].seed) : run_archipelago(topologies[jobs[j].topology].t, jobs[j].interval, jobs[j].seed);
  });
  {
    std::ofstream f(out + "/islands.csv");
    f << "topology,interval,seed,success,evaluations,best\n";
    for (std::size_t j = 0; j < jobs.size(); j++)
      f << (jobs[j].topology < 0 ? std::string("panmictic") : topologies[jobs[j].topology].name) << ',' << jobs[j].interval << ',' << jobs[j].seed << ',' << (runs[j].success ? 1 : 0) << ',' << runs[j].evaluations << ',' << runs[j].best << '\n';
  }
  std::ofstream f(out + "/islands_summary.csv");
  f << "topology,interval,runs,successes,rate,lower,upper,ert\n";
  auto summarise = [&](int t, std::size_t iv, const std::string& name) {
    std::size_t k = 0, n = 0;
    double spent = 0;
    for (std::size_t j = 0; j < jobs.size(); j++)
      if (jobs[j].topology == t && jobs[j].interval == iv) {
        n++;
        k += runs[j].success;
        spent += runs[j].evaluations;
      }
    const auto p = heuresis::analysis::clopper_pearson(k, n);
    const double ert = k ? spent / static_cast<double>(k) : inf;  // expected running time (Auger and Hansen 2005)
    f << name << ',' << iv << ',' << n << ',' << k << ',' << p.estimate << ',' << p.lower << ',' << p.upper << ',' << ert << '\n';
    std::fprintf(stderr, "  %-12s interval %2zu: %2zu/%zu  ERT %.0f\n", name.c_str(), iv, k, n, ert);
  };
  summarise(-1, 0, "panmictic");
  for (std::size_t t = 0; t < topologies.size(); t++)
    for (auto iv : intervals)
      if (t > 0 || iv == intervals[0]) summarise(static_cast<int>(t), iv, topologies[t].name);
  // The synchronous archipelago on four threads is the same run.
  const IslandRun one = run_archipelago(topologies[1].t, 5, 3), many = run_archipelago(topologies[1].t, 5, 3, &pool);
  std::cerr << "ring, interval 5, seed 3: one thread " << one.evaluations << " evaluations, best " << one.best << "; the pool " << many.evaluations << ", best " << many.best << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: machines <results/> [all | takeover | synchronisation | islands] [threads] [pilot]\n";
    return 2;
  }
  const std::string out = argv[1], part = argc > 2 ? argv[2] : "all";
  const std::size_t threads = argc > 3 ? std::stoul(argv[3]) : 4;
  const bool pilot = argc > 4 && std::string(argv[4]) == "pilot";
  calx::scheduler::ThreadPool pool(threads);
  const auto t0 = std::chrono::steady_clock::now();
  auto lap = [&](const char* what) { std::cerr << what << ": " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n"; };
  if (part == "all" || part == "takeover") {
    part_takeover(out);
    lap("takeover");
  }
  if (part == "all" || part == "synchronisation") {
    part_synchronisation(out, pool, pilot);
    lap("synchronisation");
  }
  if (part == "all" || part == "islands") {
    part_islands(out, pool, pilot);
    lap("islands");
  }
  return 0;
}
