// Rugged landscapes: Kauffman's NK landscapes of n = 20, enumerated whole (2²⁰ points), from smooth (k = 0) to
// uncorrelated (k = 19), five instances of each k.
//
// For each landscape, exactly: the global optimum; the local optima under one-bit flips and their basins under
// steepest descent; the one-step correlation ρ(1) of the random walk and its correlation length; the
// fitness–distance correlation over every point. Then seven searches of heuresis, each with the same budget of
// evaluations on every landscape (10 seeds per instance): random search, the (1+1) EA, a generational GA, simulated
// annealing on the modified Lam schedule, tabu search, the compact GA and UMDA; for each, the gap between the best
// value it found and the enumerated optimum, the rank of that value among all 2²⁰, and whether it found the optimum.
//
// The objective is heuresis's NK as written, 1 − F(x) (minimised), F the mean of the n tables' entries.
//
// Usage: rugged <results directory>

#include <heuresis/annealing/simulated_annealing.hpp>
#include <heuresis/benchmarks/pseudo_boolean.hpp>
#include <heuresis/distributions/bernoulli.hpp>
#include <heuresis/distributions/distribution_search.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/landscapes/landscapes.hpp>
#include <heuresis/local/tabu.hpp>
#include <heuresis/search/random_search.hpp>
#include <heuresis/search/run.hpp>

#include <tyche/descriptive/summary.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace ls = heuresis::landscapes;
using heuresis::search::Budget;
using heuresis::search::Draws;
using heuresis::search::tag;
using heuresis::variation::Bits;

namespace {

constexpr std::size_t n = 20;
constexpr std::size_t instances = 5;
constexpr std::size_t seeds = 10;
constexpr std::uint64_t nk_seed = 2026;
const std::vector<std::size_t> ks = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19};
const std::vector<double> budgets = {5000, 50000};
const std::vector<std::string> algorithms = {"random", "(1+1) EA", "GA", "SA", "tabu", "cGA", "UMDA"};
// The local optima of these landscapes (instance 0) are written out whole, for the picture of their basins.
const std::vector<std::size_t> pictured = {2, 6, 10};
constexpr std::size_t walk_length = 200000, max_lag = 20;

struct Landscape {
  std::size_t k = 0, instance = 0;
  double optimum_fitness = 0;  // F* = 1 − min f
  std::size_t local_optima = 0;
  double rho1_exact = 0, rho1_walk = 0, length_exact = 0, length_walk = 0;
  std::vector<double> rho;  // the walk's autocorrelation, lags 0..max_lag
  double fdc = 0;
  double global_basin = 0, largest_basin = 0;  // shares of the space
  double mean_optimum_distance = 0;            // mean Hamming distance of the local optima to the global optimum
  double optima_fitness_basin_correlation = 0;  // Spearman-free: Pearson between F and log basin size
  // The local optima, if pictured: F, basin size, Hamming distance to the global optimum.
  std::vector<std::array<double, 3>> optima;
  double seconds = 0;
};

struct Result {
  std::size_t k = 0, instance = 0, seed = 0, algorithm = 0, budget = 0;
  double gap = 0, rank = 0;  // F* − F(best), and the number of points strictly better than the best found
  bool success = false;
};

double run_algorithm(std::size_t a, double B, std::uint64_t seed, const ls::Hypercube& h) {
  auto f = [&h](const Bits& x) { return h[ls::Hypercube::index(x)]; };
  const Draws draws(seed);
  const Budget budget{B};
  auto best = [](const auto& outcome) { return outcome.best->evaluation.value(); };
  switch (a) {
    case 0: {
      heuresis::search::RandomSearch<Bits> rs([](calx::random::Stream& s) { return heuresis::variation::random_bits(n, s); }, draws);
      return best(heuresis::search::run(rs, f, budget, heuresis::search::Serial{}, 1000));
    }
    case 1: {
      auto ea = heuresis::genetic::one_plus_one_ea(n, draws);
      return best(heuresis::search::run(ea, f, budget));
    }
    case 2: {
      namespace hg = heuresis::genetic;
      hg::Operators<Bits> op{hg::random_bits(n), hg::bit_flip(), hg::uniform_crossover<std::uint8_t>()};
      hg::GaOptions o;
      o.mu = 50;
      o.crossover_rate = 0.9;
      o.tournament = 2;
      o.survival = hg::Survival::generational;
      o.elites = 1;
      hg::GeneticAlgorithm<Bits> ga(op, o, draws);
      return best(heuresis::search::run(ga, f, budget));
    }
    case 3: {
      namespace ha = heuresis::annealing;
      auto flip = [](const Bits& x, calx::random::Stream& s) {
        Bits y = x;
        heuresis::variation::flip_exactly(y, 1, s);
        return y;
      };
      ha::SimulatedAnnealing<Bits> sa([](calx::random::Stream& s) { return heuresis::variation::random_bits(n, s); }, flip, ha::Schedule::modified_lam(0.05, static_cast<std::uint64_t>(B)), draws);
      return best(heuresis::search::run(sa, f, budget));
    }
    case 4: {
      calx::random::Stream s = draws(tag("experiment.nk.tabu.start"), 0, 0);
      heuresis::local::TabuSearch<Bits> tabu(heuresis::variation::random_bits(n, s), heuresis::local::one_flip(), 5, draws);
      return best(heuresis::search::run(tabu, f, budget));
    }
    case 5: {
      namespace hd = heuresis::distributions;
      hd::DistributionSearch cga(hd::ProductBernoulli{}, hd::ProductBernoulli::uniform(n), 2, hd::Compact{50}, draws);
      return best(heuresis::search::run(cga, f, budget));
    }
    default: {
      namespace hd = heuresis::distributions;
      hd::DistributionSearch umda(hd::ProductBernoulli{}, hd::ProductBernoulli::uniform(n), 100, hd::Umda{50}, draws);
      return best(heuresis::search::run(umda, f, budget));
    }
  }
}

void analyse(Landscape& L, std::vector<Result>& results) {
  const auto start = std::chrono::steady_clock::now();
  ls::Hypercube h;
  {
    const heuresis::benchmarks::NK nk(n, L.k, nk_seed, L.instance);
    h = ls::Hypercube::tabulate(n, nk);
  }
  const auto values = h.values();
  const std::size_t N = h.size();
  const std::size_t opt = static_cast<std::size_t>(std::min_element(values.begin(), values.end()) - values.begin());
  const double fstar = values[opt];
  L.optimum_fitness = 1 - fstar;
  L.local_optima = ls::local_minima(h).size();
  L.rho1_exact = ls::one_step_correlation(h);
  L.length_exact = ls::correlation_length(L.rho1_exact);
  {
    calx::random::Stream s = Draws(nk_seed + L.instance)(tag("experiment.nk.walk.start"), L.k, 0);
    auto flip = [](const Bits& x, calx::random::Stream& r) {
      Bits y = x;
      heuresis::variation::flip_exactly(y, 1, r);
      return y;
    };
    const auto walk = ls::random_walk(heuresis::variation::random_bits(n, s), [&h](const Bits& x) { return h[ls::Hypercube::index(x)]; }, flip, walk_length, Draws(nk_seed + 100 * L.k + L.instance));
    L.rho = ls::autocorrelation(walk, max_lag);
    L.rho1_walk = L.rho[1];
    L.length_walk = ls::correlation_length(L.rho1_walk);
  }
  {
    std::vector<double> d(N);
    for (std::size_t m = 0; m < N; m++) d[m] = static_cast<double>(std::popcount(m ^ opt));
    L.fdc = ls::fitness_distance_correlation(values, d);
  }
  {
    const auto b = ls::basins(h);
    std::size_t largest = 0, global = 0;
    double dist = 0;
    std::vector<double> F, logs;
    for (std::size_t q = 0; q < b.optima.size(); q++) {
      largest = std::max(largest, b.sizes[q]);
      if (b.optima[q] == opt) global = b.sizes[q];
      dist += static_cast<double>(std::popcount(b.optima[q] ^ opt));
      F.push_back(1 - values[b.optima[q]]);
      logs.push_back(std::log(static_cast<double>(b.sizes[q])));
      if (L.instance == 0 && std::find(pictured.begin(), pictured.end(), L.k) != pictured.end())
        L.optima.push_back({1 - values[b.optima[q]], static_cast<double>(b.sizes[q]), static_cast<double>(std::popcount(b.optima[q] ^ opt))});
    }
    L.global_basin = static_cast<double>(global) / static_cast<double>(N);
    L.largest_basin = static_cast<double>(largest) / static_cast<double>(N);
    L.mean_optimum_distance = dist / static_cast<double>(b.optima.size());
    L.optima_fitness_basin_correlation = b.optima.size() > 2 ? tyche::descriptive::correlation(F, logs) : NAN;
  }
  std::vector<double> sorted(values.begin(), values.end());
  std::sort(sorted.begin(), sorted.end());
  for (std::size_t bi = 0; bi < budgets.size(); bi++)
    for (std::size_t a = 0; a < algorithms.size(); a++)
      for (std::size_t s = 0; s < seeds; s++) {
        const std::uint64_t seed = 1000000ull * L.k + 1000ull * L.instance + 10ull * bi + s + 1;
        const double v = run_algorithm(a, budgets[bi], seed, h);
        Result r;
        r.k = L.k;
        r.instance = L.instance;
        r.seed = s + 1;
        r.algorithm = a;
        r.budget = bi;
        r.gap = v - fstar;  // = F* − F(best)
        r.rank = static_cast<double>(std::lower_bound(sorted.begin(), sorted.end(), v) - sorted.begin());
        r.success = v == fstar;
        results.push_back(r);
      }
  L.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: rugged <results directory>\n";
    return 2;
  }
  const std::filesystem::path dir = argv[1];
  std::filesystem::create_directories(dir);
  std::vector<Landscape> L;
  for (auto k : ks)
    for (std::size_t i = 0; i < instances; i++) {
      Landscape l;
      l.k = k;
      l.instance = i;
      L.push_back(l);
    }
  // The k = 19 landscapes (160 MB of tables each) first, so that no more than a few are held at once at the end.
  std::vector<std::size_t> order(L.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return L[a].k > L[b].k; });
  std::vector<std::vector<Result>> results(L.size());
  std::atomic<std::size_t> next{0};
  std::mutex log;
  const unsigned threads = std::max(1u, std::min(6u, std::thread::hardware_concurrency()));
  const auto start = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < threads; t++)
    pool.emplace_back([&] {
      for (std::size_t i; (i = next++) < order.size();) {
        analyse(L[order[i]], results[order[i]]);
        std::lock_guard<std::mutex> g(log);
        std::cerr << "k = " << L[order[i]].k << " instance " << L[order[i]].instance << ": " << L[order[i]].local_optima << " local optima, " << L[order[i]].seconds << " s\n";
      }
    });
  for (auto& t : pool) t.join();
  const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  namespace td = tyche::descriptive;

  {
    std::ofstream out(dir / "landscapes.csv");
    out << "k,instance,optimum_fitness,local_optima,rho1_exact,rho1_walk,length_exact,length_walk,fdc,global_basin,largest_basin,mean_optimum_distance,fitness_basin_correlation\n";
    out << std::setprecision(8);
    for (const auto& l : L)
      out << l.k << "," << l.instance << "," << l.optimum_fitness << "," << l.local_optima << "," << l.rho1_exact << "," << l.rho1_walk << "," << l.length_exact << "," << l.length_walk << ","
          << l.fdc << "," << l.global_basin << "," << l.largest_basin << "," << l.mean_optimum_distance << "," << l.optima_fitness_basin_correlation << "\n";
  }
  {
    std::ofstream out(dir / "autocorrelation.csv");
    out << "k,lag,rho_mean\n";
    out << std::setprecision(6);
    for (auto k : ks)
      for (std::size_t s = 0; s <= max_lag; s++) {
        double m = 0;
        for (const auto& l : L)
          if (l.k == k) m += l.rho[s];
        out << k << "," << s << "," << m / static_cast<double>(instances) << "\n";
      }
  }
  {
    std::ofstream out(dir / "optima.csv");
    out << "k,instance,fitness,basin,distance\n";
    out << std::setprecision(7);
    for (const auto& l : L)
      for (const auto& o : l.optima) out << l.k << "," << l.instance << "," << o[0] << "," << o[1] << "," << o[2] << "\n";
  }
  {
    // Over the 5 instances × 10 seeds of each k: the quartiles of the gap, the median rank, the share of runs that
    // found the optimum.
    std::ofstream out(dir / "performance.csv");
    out << "k,algorithm,budget,gap_q1,gap_median,gap_q3,gap_mean,rank_median,rank_q3,success\n";
    out << std::setprecision(7);
    for (auto k : ks)
      for (std::size_t bi = 0; bi < budgets.size(); bi++)
        for (std::size_t a = 0; a < algorithms.size(); a++) {
          std::vector<double> gap, rank;
          double success = 0;
          for (std::size_t li = 0; li < L.size(); li++)
            if (L[li].k == k)
              for (const auto& r : results[li])
                if (r.algorithm == a && r.budget == bi) {
                  gap.push_back(r.gap);
                  rank.push_back(r.rank);
                  success += r.success;
                }
          out << k << "," << algorithms[a] << "," << budgets[bi] << "," << td::quantile(gap, 0.25) << "," << td::median(gap) << "," << td::quantile(gap, 0.75) << "," << td::mean(gap) << ","
              << td::median(rank) << "," << td::quantile(rank, 0.75) << "," << success / static_cast<double>(gap.size()) << "\n";
        }
  }
  {
    // Per instance: how many of the ten seeds found the optimum, and the median gap.
    std::ofstream out(dir / "instances.csv");
    out << "k,instance,algorithm,budget,successes,gap_median\n";
    out << std::setprecision(7);
    for (std::size_t li = 0; li < L.size(); li++)
      for (std::size_t bi = 0; bi < budgets.size(); bi++)
        for (std::size_t a = 0; a < algorithms.size(); a++) {
          std::vector<double> gap;
          int successes = 0;
          for (const auto& r : results[li])
            if (r.algorithm == a && r.budget == bi) {
              gap.push_back(r.gap);
              successes += r.success;
            }
          out << L[li].k << "," << L[li].instance << "," << algorithms[a] << "," << budgets[bi] << "," << successes << "," << td::median(gap) << "\n";
        }
  }
  {
    std::ofstream out(dir / "reference.json");
    out << std::setprecision(10) << "{\n";
    out << "  \"n\": " << n << ", \"k\": [";
    for (std::size_t i = 0; i < ks.size(); i++) out << (i ? ", " : "") << ks[i];
    out << "], \"instances\": " << instances << ", \"nk_seed\": " << nk_seed << ", \"seeds_per_instance\": " << seeds << ",\n";
    out << "  \"budgets\": [5000, 50000], \"walk_length\": " << walk_length << ",\n";
    out << "  \"algorithms\": {\"random\": \"uniform bit strings\", \"(1+1) EA\": \"standard bit mutation 1/n\", \"GA\": \"mu 50, uniform crossover 0.9, bit flip 1/n, binary tournament, generational, 1 elite\","
           " \"SA\": \"one-bit flips, modified Lam schedule, T0 0.05, steps = budget\", \"tabu\": \"one-flip neighbourhood, tenure 5, random start\", \"cGA\": \"K 50\", \"UMDA\": \"lambda 100, mu 50, borders 1/n\"},\n";
    out << "  \"threads\": " << threads << ", \"wall_seconds\": " << total << "\n}\n";
  }
  std::cerr << "done in " << total << " s\n";
  return 0;
}
