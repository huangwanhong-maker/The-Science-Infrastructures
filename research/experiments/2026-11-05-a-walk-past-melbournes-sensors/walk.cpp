// A walk past Melbourne's sensors: the shortest closed walk through the City of Melbourne's pedestrian
// counters, by every permutation search heuresis holds, on geodesic distances from the geodesy domain.
//
//   walk <sensors.csv> <results/> [budget] [seeds] [threads] [restarts] [kicks] [algorithms] [ant seeds]
//
// The problem. The n sensors of the Pedestrian Counting System (sensors.csv of the data library's
// melbourne-pedestrian-counts: only their latitudes and longitudes on WGS 84 are read) are the cities of a
// symmetric travelling-salesman problem; the distance between two sensors is the length of the geodesic between
// them on the WGS 84 ellipsoid (Karney 2013, geodesy::geodesic). A tour is a permutation read as a cycle; its
// length, in metres, is minimised.
//
// The best known tour. Iterated local search with white-box 2-opt (heuresis::local::two_opt_descent, which
// evaluates a move from four distances) and the double-bridge kick (heuresis::local::double_bridge), accepting a
// kicked optimum that is not worse, from `restarts` random starts with `kicks` kicks each. It is not proved optimal.
//
// The comparison. Each algorithm is a black box on the protocol: it sees a tour's length and nothing else
// (the ant colony also receives the inverse distances as its heuristic, as the method requires). Every run is
// driven by heuresis::search::run for `budget` evaluations; its anytime curve (each improvement of the best) is
// read at log-spaced checkpoints by heuresis::analysis::best_at. Runs are independent and run on a calx thread
// pool, one run per task; their results do not depend on the pool.
//
// Outputs (results/): sensors.csv (positions on WGS 84 and in the local ENU frame of the data library's city
// model), best_known.json, convergence.csv, final.csv, tours.json, timing.json.

#include <heuresis/analysis/runtimes.hpp>
#include <heuresis/annealing/annealing.hpp>
#include <heuresis/distributions/distributions.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/local/local.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/swarms/swarms.hpp>

#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>
#include <geodesy/geodesic/geodesic.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using heuresis::search::Budget;
using heuresis::search::Draws;
using heuresis::search::tag;
using heuresis::variation::Permutation;

struct Sensor {
  int id = 0;
  double latitude = 0, longitude = 0;  // degrees, WGS 84
  double east = 0, north = 0;          // metres, local ENU at the city model's origin
};

std::vector<std::string> split(const std::string& line, char sep) {
  std::vector<std::string> out;
  std::string field;
  std::istringstream in(line);
  while (std::getline(in, field, sep)) out.push_back(field);
  if (!line.empty() && line.back() == sep) out.emplace_back();
  return out;
}

// The sensors' positions only: location_id, latitude, longitude (semicolons, a byte-order mark).
std::vector<Sensor> read_sensors(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::string line;
  std::getline(in, line);
  if (line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF) line = line.substr(3);
  if (!line.empty() && line.back() == '\r') line.pop_back();
  const auto header = split(line, ';');
  auto column = [&](const std::string& name) {
    const auto it = std::find(header.begin(), header.end(), name);
    if (it == header.end()) throw std::runtime_error("no column " + name);
    return static_cast<std::size_t>(it - header.begin());
  };
  const std::size_t c_id = column("location_id"), c_lat = column("latitude"), c_lon = column("longitude");
  std::vector<Sensor> out;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    const auto f = split(line, ';');
    if (f.size() != header.size()) throw std::runtime_error("a row of " + std::to_string(f.size()) + " fields: " + line);
    out.push_back({std::stoi(f[c_id]), std::stod(f[c_lat]), std::stod(f[c_lon])});
  }
  return out;
}

struct Problem {
  std::size_t n = 0;
  std::vector<double> d;  // geodesic distances (m), row by row
  double operator()(std::size_t a, std::size_t b) const { return d[a * n + b]; }
  double length(const Permutation& p) const {
    double L = 0;
    for (std::size_t i = 0; i < n; i++) L += d[p[i] * n + p[(i + 1) % n]];
    return L;
  }
};

// One run's record: its anytime curve read at the checkpoints, its best tour.
struct Result {
  std::string algorithm;
  std::uint64_t seed = 0;
  std::vector<double> at;  // best length after each checkpoint
  double evaluations = 0, best = 0;
  Permutation tour;
  double seconds = 0;
};

template <class A>
Result drive(const std::string& name, std::uint64_t seed, A& a, const Problem& P, double budget, const std::vector<double>& checkpoints, std::size_t batch = heuresis::search::unlimited) {
  const auto t0 = std::chrono::steady_clock::now();
  auto f = [&P](const Permutation& p) { return P.length(p); };
  auto o = heuresis::search::run(a, f, Budget{budget}, heuresis::search::Serial{}, batch);
  Result r;
  r.algorithm = name;
  r.seed = seed;
  const heuresis::analysis::Run run = heuresis::analysis::run_of(o);
  for (double c : checkpoints) r.at.push_back(heuresis::analysis::best_at(run, c));
  r.evaluations = o.evaluations;
  r.best = o.best->evaluation.value();
  r.tour = o.best->genotype;
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return r;
}

constexpr std::uint32_t start_tag = tag("experiment.melbourne.start");
constexpr std::uint32_t kick_tag = tag("experiment.melbourne.kick");
constexpr std::uint32_t pilot_tag = tag("experiment.melbourne.pilot");

// The GA's mutation: an inversion with probability 0.2, else an insertion with probability 0.2, else none.
Permutation mutate(const Permutation& p, calx::random::Stream& s) {
  const double u = s.uniform();
  if (u < 0.2) return heuresis::variation::inversion_mutation(p, s);
  if (u < 0.4) return heuresis::variation::insertion_mutation(p, s);
  return p;
}

struct Gumbel {
  Permutation operator()(const heuresis::distributions::PlackettLuce& f, const std::vector<double>& w, calx::random::Stream& r) const { return f.sample_by_gumbel(w, r); }
};

Result run_algorithm(const std::string& name, std::uint64_t seed, const Problem& P, const std::vector<double>& eta, double budget, const std::vector<double>& checkpoints, double t0_sa) {
  using namespace heuresis;
  const Draws draws(seed);
  const std::size_t n = P.n;
  auto sample = [n](calx::random::Stream& s) { return variation::random_permutation(n, s); };
  if (name == "random") {
    search::RandomSearch<Permutation> a(sample, draws);
    return drive(name, seed, a, P, budget, checkpoints, 1000);
  }
  if (name == "GA-OX" || name == "GA-PMX" || name == "GA-CX") {
    genetic::Crossover<Permutation> x = name == "GA-OX" ? genetic::Crossover<Permutation>(variation::order_crossover)
                                      : name == "GA-PMX" ? genetic::Crossover<Permutation>(variation::partially_mapped_crossover)
                                                         : genetic::Crossover<Permutation>(variation::cycle_crossover);
    genetic::GaOptions o;
    o.mu = 100;
    o.lambda = 100;
    o.crossover_rate = 0.9;
    o.tournament = 3;
    o.elites = 2;
    genetic::GeneticAlgorithm<Permutation> a({sample, mutate, x}, o, draws);
    return drive(name, seed, a, P, budget, checkpoints);
  }
  if (name == "MMAS") {
    swarms::AcoOptions o;
    o.variant = swarms::AcoVariant::max_min;
    o.ants = 25;
    o.alpha = 1;
    o.beta = 3;
    o.rho = 0.02;
    o.p_best = 0.05;
    swarms::AntColony a(n, eta, o, draws);
    return drive(name, seed, a, P, budget, checkpoints);
  }
  if (name == "ILS") {
    calx::random::Stream s = draws(start_tag, 0, 0);
    local::IteratedLocalSearch<Permutation> a(variation::random_permutation(n, s), local::two_opt(), local::double_bridge, local::Pivot::first, local::Acceptance::not_worse, draws);
    return drive(name, seed, a, P, budget, checkpoints);
  }
  if (name == "SA") {
    // Geometric cooling from T₀ to T₀/10⁴ over the budget, 100 steps at each level.
    const std::uint64_t level = 100;
    const double alpha = std::pow(1e-4, static_cast<double>(level) / budget);
    annealing::SimulatedAnnealing<Permutation> a(sample, variation::inversion_mutation, annealing::Schedule::geometric(t0_sa, alpha, level), draws);
    return drive(name, seed, a, P, budget, checkpoints);
  }
  if (name == "Plackett-Luce") {
    distributions::PlackettLuceRule rule{20};
    rule.iterations = 10;
    rule.smoothing = 0.3;
    rule.floor = 0.01;
    distributions::DistributionSearch<distributions::PlackettLuce, distributions::PlackettLuceRule, Gumbel> a(distributions::PlackettLuce{}, distributions::PlackettLuce::uniform(n), 200, rule, draws, Gumbel{});
    return drive(name, seed, a, P, budget, checkpoints);
  }
  throw std::invalid_argument("unknown algorithm " + name);
}

std::string tour_json(const Permutation& p) {
  std::string s = "[";
  for (std::size_t i = 0; i < p.size(); i++) s += (i ? "," : "") + std::to_string(p[i]);
  return s + "]";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: walk <sensors.csv> <results/> [budget] [seeds] [threads] [restarts] [kicks] [algorithms]\n";
    return 2;
  }
  const std::string sensors_path = argv[1], out = argv[2];
  const double budget = argc > 3 ? std::stod(argv[3]) : 1e6;
  const std::size_t seeds = argc > 4 ? std::stoul(argv[4]) : 10;
  const std::size_t threads = argc > 5 ? std::stoul(argv[5]) : 4;
  const std::size_t restarts = argc > 6 ? std::stoul(argv[6]) : 32;
  const std::size_t kicks = argc > 7 ? std::stoul(argv[7]) : 20000;
  std::vector<std::string> algorithms{"random", "Plackett-Luce", "GA-CX", "GA-PMX", "GA-OX", "MMAS", "SA", "ILS"};
  if (argc > 8 && std::string(argv[8]) != "all") algorithms = split(argv[8], ',');
  // The ant colony's construction costs O(n²) powers per ant, two orders more than the others' per evaluation;
  // it may be given fewer seeds at the same budget.
  const std::size_t ant_seeds = argc > 9 ? std::stoul(argv[9]) : seeds;

  // The sensors, their geodesics and their places in the city model's frame.
  std::vector<Sensor> sensors = read_sensors(sensors_path);
  const std::size_t n = sensors.size();
  const auto wgs84 = geodesy::figure::Ellipsoid::wgs84();
  const auto geodesic = geodesy::geodesic::Geodesic::wgs84();
  constexpr double deg = std::numbers::pi / 180;
  const double lat0 = -37.8085, lon0 = 144.9525;  // the origin of the data library's city model of Melbourne
  const auto place = geodesy::figure::Place::at(wgs84, {lat0 * deg, lon0 * deg, 0});
  for (auto& s : sensors) {
    const auto enu = place.enu_of(wgs84.cartesian({s.latitude * deg, s.longitude * deg, 0}));
    s.east = enu[0];
    s.north = enu[1];
  }
  Problem P;
  P.n = n;
  P.d.assign(n * n, 0);
  double worst_relative = 0, worst_absolute = 0, asymmetry = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const double s = geodesic.inverse_degrees(sensors[i].latitude, sensors[i].longitude, sensors[j].latitude, sensors[j].longitude).s12;
      const double back = geodesic.inverse_degrees(sensors[j].latitude, sensors[j].longitude, sensors[i].latitude, sensors[i].longitude).s12;
      asymmetry = std::max(asymmetry, std::fabs(s - back));
      P.d[i * n + j] = P.d[j * n + i] = s;
      const double chord = std::hypot(sensors[i].east - sensors[j].east, sensors[i].north - sensors[j].north);
      worst_absolute = std::max(worst_absolute, std::fabs(chord - s));
      worst_relative = std::max(worst_relative, std::fabs(chord - s) / s);
    }
  std::vector<double> eta(n * n, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (i != j) eta[i * n + j] = 1 / P(i, j);
  {
    std::ofstream f(out + "/sensors.csv");
    f << "index,location_id,latitude,longitude,east_m,north_m\n";
    f.precision(10);
    for (std::size_t i = 0; i < n; i++) f << i << ',' << sensors[i].id << ',' << sensors[i].latitude << ',' << sensors[i].longitude << ',' << sensors[i].east << ',' << sensors[i].north << '\n';
  }
  std::cerr << n << " sensors; geodesic against the local frame's chord: at most " << worst_absolute << " m (" << worst_relative << " relative); asymmetry " << asymmetry << " m\n";

  calx::scheduler::ThreadPool pool(threads);

  // The best known tour: white-box ILS from `restarts` random starts.
  const auto tb = std::chrono::steady_clock::now();
  std::vector<double> restart_best(restarts);
  std::vector<std::size_t> restart_kick(restarts);
  std::vector<Permutation> restart_tour(restarts);
  calx::scheduler::parallel_for(pool, restarts, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t r = b; r < e; r++) {
      const Draws draws(100000 + r);
      calx::random::Stream s = draws(start_tag, 0, 0);
      Permutation x = heuresis::variation::random_permutation(n, s);
      heuresis::local::two_opt_descent(x, P);
      double lx = P.length(x);
      Permutation best = x;
      double lb = lx;
      std::size_t at = 0;
      for (std::size_t k = 0; k < kicks; k++) {
        calx::random::Stream ks = draws(kick_tag, k, 0);
        Permutation y = heuresis::local::double_bridge(x, ks);
        heuresis::local::two_opt_descent(y, P);
        const double ly = P.length(y);
        if (ly <= lx) {
          x = std::move(y);
          lx = ly;
          if (lx < lb - 1e-9) {
            lb = lx;
            best = x;
            at = k + 1;
          }
        }
      }
      restart_best[r] = lb;
      restart_kick[r] = at;
      restart_tour[r] = best;
    }
  });
  const std::size_t rb = static_cast<std::size_t>(std::min_element(restart_best.begin(), restart_best.end()) - restart_best.begin());
  const double best_known = restart_best[rb];
  const Permutation& best_tour = restart_tour[rb];
  std::size_t reached = 0;
  for (double v : restart_best)
    if (v <= best_known + 1e-6) reached++;
  double chord_length = 0;
  for (std::size_t i = 0; i < n; i++) {
    const auto& a = sensors[best_tour[i]];
    const auto& c = sensors[best_tour[(i + 1) % n]];
    chord_length += std::hypot(a.east - c.east, a.north - c.north);
  }
  const double best_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tb).count();
  std::cerr << "best known " << best_known << " m, reached by " << reached << " of " << restarts << " restarts (" << best_seconds << " s)\n";
  {
    std::ofstream f(out + "/best_known.json");
    f.precision(12);
    f << "{\n  \"sensors\": " << n << ",\n  \"length_m\": " << best_known << ",\n  \"length_in_local_frame_m\": " << chord_length << ",\n  \"restarts\": " << restarts << ",\n  \"kicks\": " << kicks
      << ",\n  \"restarts_reaching_best\": " << reached << ",\n  \"restart_lengths_m\": [";
    for (std::size_t r = 0; r < restarts; r++) f << (r ? ", " : "") << restart_best[r];
    f << "],\n  \"restart_kick_of_best\": [";
    for (std::size_t r = 0; r < restarts; r++) f << (r ? ", " : "") << restart_kick[r];
    f << "],\n  \"geodesic_vs_chord_max_abs_m\": " << worst_absolute << ",\n  \"geodesic_vs_chord_max_rel\": " << worst_relative << ",\n  \"geodesic_asymmetry_max_m\": " << asymmetry
      << ",\n  \"seconds\": " << best_seconds << ",\n  \"tour\": " << tour_json(best_tour) << "\n}\n";
  }

  // The temperature of annealing: the mean |Δ| of an inversion on random tours (a pilot of 2000 moves).
  double t0_sa = 0;
  {
    const Draws draws(424242);
    for (std::size_t k = 0; k < 2000; k++) {
      calx::random::Stream s = draws(pilot_tag, k, 0);
      Permutation p = heuresis::variation::random_permutation(n, s);
      t0_sa += std::fabs(P.length(heuresis::variation::inversion_mutation(p, s)) - P.length(p));
    }
    t0_sa /= 2000;
  }
  std::cerr << "annealing from T0 = " << t0_sa << " m\n";

  // The comparison.
  std::vector<double> checkpoints;
  for (double c = 100; c < budget * 1.0001; c *= std::pow(10.0, 1.0 / 16)) checkpoints.push_back(std::round(c));
  if (checkpoints.back() != budget) checkpoints.push_back(budget);
  struct Job {
    std::string algorithm;
    std::uint64_t seed;
  };
  std::vector<Job> jobs;
  // The slow ones first, so that the pool's dynamic schedule ends together.
  for (const char* a : {"MMAS", "Plackett-Luce", "ILS", "GA-OX", "GA-PMX", "GA-CX", "SA", "random"})
    if (std::find(algorithms.begin(), algorithms.end(), a) != algorithms.end())
      for (std::uint64_t s = 0; s < (std::string(a) == "MMAS" ? ant_seeds : seeds); s++) jobs.push_back({a, s});
  std::vector<Result> results(jobs.size());
  const auto tc = std::chrono::steady_clock::now();
  calx::scheduler::parallel_for(pool, jobs.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t j = b; j < e; j++) {
      results[j] = run_algorithm(jobs[j].algorithm, jobs[j].seed, P, eta, budget, checkpoints, t0_sa);
      std::fprintf(stderr, "  %-14s seed %2llu: %10.1f m (%+.2f%%) in %.1f s\n", results[j].algorithm.c_str(), static_cast<unsigned long long>(results[j].seed), results[j].best, 100 * (results[j].best / best_known - 1), results[j].seconds);
    }
  });
  const double compare_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tc).count();

  {
    std::ofstream f(out + "/convergence.csv");
    f << "algorithm,seed,evaluations,best_m\n";
    f.precision(10);
    for (const auto& r : results)
      for (std::size_t c = 0; c < checkpoints.size(); c++)
        if (std::isfinite(r.at[c])) f << r.algorithm << ',' << r.seed << ',' << checkpoints[c] << ',' << r.at[c] << '\n';
  }
  {
    std::ofstream f(out + "/final.csv");
    f << "algorithm,seed,evaluations,best_m,excess,seconds\n";
    f.precision(10);
    for (const auto& r : results) f << r.algorithm << ',' << r.seed << ',' << r.evaluations << ',' << r.best << ',' << r.best / best_known - 1 << ',' << r.seconds << '\n';
  }
  {
    // Each algorithm's best tour over its seeds at the budget.
    std::map<std::string, const Result*> best;
    for (const auto& r : results)
      if (!best.count(r.algorithm) || r.best < best[r.algorithm]->best) best[r.algorithm] = &r;
    std::ofstream f(out + "/tours.json");
    f.precision(10);
    f << "{\n  \"budget\": " << budget << ",\n  \"best_known_m\": " << best_known << ",\n  \"best_known\": " << tour_json(best_tour) << ",\n  \"algorithms\": {";
    bool first = true;
    for (const auto& [name, r] : best) {
      f << (first ? "" : ",") << "\n    \"" << name << "\": {\"seed\": " << r->seed << ", \"length_m\": " << r->best << ", \"tour\": " << tour_json(r->tour) << "}";
      first = false;
    }
    f << "\n  }\n}\n";
  }
  {
    std::map<std::string, double> seconds;
    for (const auto& r : results) seconds[r.algorithm] += r.seconds;
    std::ofstream f(out + "/timing.json");
    f << "{\n  \"threads\": " << threads << ",\n  \"best_known_seconds\": " << best_seconds << ",\n  \"comparison_seconds\": " << compare_seconds << ",\n  \"t0_annealing_m\": " << t0_sa
      << ",\n  \"seeds\": " << seeds << ",\n  \"ant_seeds\": " << ant_seeds << ",\n  \"budget\": " << budget << ",\n  \"cpu_seconds_per_algorithm\": {";
    bool first = true;
    for (const auto& [name, s] : seconds) {
      f << (first ? "" : ",") << "\n    \"" << name << "\": " << s;
      first = false;
    }
    f << "\n  }\n}\n";
  }
  std::cerr << "comparison: " << compare_seconds << " s\n";
  return 0;
}
