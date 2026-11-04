// Quantum search, as the search foundation sees it: the law of a device, what the law costs in oracle calls,
// and a quantum-inspired algorithm that is a product law in disguise.
//
//   quantum <results/> [threads]
//
// Nothing here simulates a quantum computer. heuresis::quantum::EmulatedGrover draws measurements from the exact
// law that r Grover iterations leave on a register (a marked point with probability sin²((2r + 1)θ), sin²θ = t/N,
// uniformly within the marked and within the unmarked points); to know which points are marked it evaluates every
// point once, classically, a cost reported apart. The speed-up measured is in oracle calls (Grover iterations), the
// unit in which Dürr and Høyer state it, not in seconds.
//
// (1) The law. On N = 4096 points with t marked, the closed form against the emulator's frequencies.
// (2) Grover adaptive search (Dürr and Høyer 1996; Baritompa, Bulger and Wood 2005) on N = 2⁶ … 2²⁰ points with a
//     unique minimum: the oracle calls until the threshold holds the minimum, over 200 runs at each size, against
//     Dürr and Høyer's bound and against classical random search (heuresis::search::RandomSearch), whose
//     expectation is N evaluations. On a QUBO of 20 bits, the threshold's descent, measurement by measurement.
// (3) The quantum-inspired evolutionary algorithm (Han and Kim 2002) on the same QUBO, beside PBIL: its
//     observations against the product-Bernoulli family's reference draws on the same streams, and its angles.

#include <heuresis/benchmarks/benchmarks.hpp>
#include <heuresis/distributions/distributions.hpp>
#include <heuresis/quantum/quantum.hpp>
#include <heuresis/search/search.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <string>
#include <vector>

namespace {

using heuresis::search::Budget;
using heuresis::search::Draws;
using heuresis::search::tag;
using heuresis::variation::Bits;

constexpr double inf = std::numeric_limits<double>::infinity();

// ---- (1) the law ---------------------------------------------------------------------------------------------

void part_law(const std::string& out) {
  const std::uint64_t N = 4096;
  heuresis::quantum::EmulatedGrover device(N, [](std::uint64_t i) { return static_cast<double>(i); });
  const std::vector<std::uint64_t> marked{1024, 256, 64, 16, 4, 1};
  std::ofstream law(out + "/grover_law.csv");
  law << "N,t,r,probability\n";
  law.precision(17);
  for (auto t : marked)
    for (std::uint64_t r = 0; r <= 60; r++) law << N << ',' << t << ',' << r << ',' << heuresis::quantum::grover_success(t, N, r) << '\n';
  std::ofstream em(out + "/grover_emulated.csv");
  em << "N,t,r,measurements,hits,frequency,probability,z\n";
  em.precision(10);
  const std::uint64_t K = 4000;
  const Draws draws(2026);
  constexpr std::uint32_t purpose = tag("experiment.quantum.law");
  double worst = 0;
  std::size_t cells = 0;
  for (std::size_t m = 0; m < marked.size(); m++)
    for (std::uint64_t r = 0; r <= 60; r += 2) {
      const double p = heuresis::quantum::grover_success(marked[m], N, r);
      std::uint64_t hits = 0;
      for (std::uint64_t k = 0; k < K; k++) {
        calx::random::Stream s = draws(purpose, m * 1000 + r, k);
        const auto got = device.measure(static_cast<double>(marked[m]), r, s);
        hits += got.index < marked[m];
      }
      const double freq = static_cast<double>(hits) / static_cast<double>(K);
      const double sd = std::sqrt(std::max(p * (1 - p), 1e-12) / static_cast<double>(K));
      const double z = (freq - p) / sd;
      if (p > 1e-9 && p < 1 - 1e-9) {
        worst = std::max(worst, std::fabs(z));
        cells++;
      }
      em << N << ',' << marked[m] << ',' << r << ',' << K << ',' << hits << ',' << freq << ',' << p << ',' << z << '\n';
    }
  std::cerr << "the emulator against the law: the largest |z| over " << cells << " cells of " << K << " measurements is " << worst << '\n';
}

// ---- (2) Grover adaptive search ------------------------------------------------------------------------------

Bits index_bits(std::uint64_t i, std::size_t n) {
  Bits x(n);
  for (std::size_t b = 0; b < n; b++) x[b] = static_cast<std::uint8_t>((i >> b) & 1u);
  return x;
}

void part_adaptive(const std::string& out, calx::scheduler::ThreadPool& pool, const heuresis::benchmarks::Qubo& qubo) {
  using namespace heuresis;
  constexpr std::uint32_t values_tag = tag("experiment.quantum.values");
  struct Size {
    std::uint64_t N;
    std::vector<double> calls, measurements, random;
  };
  std::vector<Size> sizes;
  for (unsigned k = 6; k <= 20; k++) sizes.push_back({std::uint64_t{1} << k, {}, {}, {}});
  const std::size_t runs = 200, random_runs = 100;
  const auto t0 = std::chrono::steady_clock::now();
  for (auto& S : sizes) {
    // A landscape with a unique minimum: the values a random permutation of 0..N−1.
    calx::random::Stream vs = Draws(11)(values_tag, S.N, 0);
    std::vector<std::uint32_t> values(S.N);
    std::iota(values.begin(), values.end(), 0u);
    for (std::size_t i = S.N; i > 1; i--) std::swap(values[i - 1], values[vs.below(i)]);
    // One emulator for every run: measure() reads it and changes nothing, so the runs may share it.
    quantum::EmulatedGrover device(S.N, [&](std::uint64_t i) { return static_cast<double>(values[i]); });
    S.calls.resize(runs);
    S.measurements.resize(runs);
    calx::scheduler::parallel_for(pool, runs, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t r = b; r < e; r++) {
        quantum::GroverAdaptiveSearch<quantum::EmulatedGrover> gas(device, {}, Draws(r));
        Budget budget;
        budget.target = 0;  // the minimum
        const auto o = search::run(gas, [&](const std::uint64_t& i) { return static_cast<double>(values[i]); }, budget);
        if (o.ending != search::Ending::target) throw std::runtime_error("Grover adaptive search did not reach the minimum");
        S.calls[r] = static_cast<double>(gas.oracle_calls());
        S.measurements[r] = o.evaluations;
      }
    });
    if (S.N <= (std::uint64_t{1} << 14)) {
      S.random.resize(random_runs);
      calx::scheduler::parallel_for(pool, random_runs, 1, [&](std::size_t b, std::size_t e, std::size_t) {
        for (std::size_t r = b; r < e; r++) {
          const std::uint64_t n = S.N;
          search::RandomSearch<std::uint64_t> rs([n](calx::random::Stream& s) { return s.below(n); }, Draws(r));
          Budget budget;
          budget.target = 0;
          const auto o = search::run(rs, [&](const std::uint64_t& i) { return static_cast<double>(values[i]); }, budget, search::Serial{}, 256);
          S.random[r] = o.progress.back().evaluations;  // the evaluation at which the minimum was drawn
        }
      });
    }
    double mc = 0;
    for (double c : S.calls) mc += c / static_cast<double>(runs);
    std::fprintf(stderr, "  N = %7llu: mean oracle calls %8.1f (bound %8.1f), random search %.0f\n", static_cast<unsigned long long>(S.N), mc, quantum::dh_expected_bound(S.N),
                 S.random.empty() ? 0.0 : std::accumulate(S.random.begin(), S.random.end(), 0.0) / static_cast<double>(S.random.size()));
  }
  std::cerr << "  (" << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s)\n";
  {
    std::ofstream f(out + "/gas_scaling.csv");
    f << "N,method,run,oracle_calls,measurements\n";
    for (const auto& S : sizes) {
      for (std::size_t r = 0; r < S.calls.size(); r++) f << S.N << ",grover," << r << ',' << S.calls[r] << ',' << S.measurements[r] << '\n';
      for (std::size_t r = 0; r < S.random.size(); r++) f << S.N << ",random," << r << ",0," << S.random[r] << '\n';
    }
    std::ofstream g(out + "/gas_bound.csv");
    g << "N,dh_expected_bound,dh_cutoff\n";
    for (const auto& S : sizes) g << S.N << ',' << quantum::dh_expected_bound(S.N) << ',' << quantum::dh_cutoff(S.N) << '\n';
  }

  // The threshold's descent on a QUBO of 20 bits, measurement by measurement.
  const std::size_t n = qubo.n();
  const std::uint64_t N = std::uint64_t{1} << n;
  auto f = [&](const std::uint64_t& i) { return qubo.value(index_bits(i, n)); };
  quantum::EmulatedGrover device(N, [&](std::uint64_t i) { return f(i); });
  std::ofstream trace(out + "/gas_descent.csv");
  trace << "run,measurement,oracle_calls,threshold,rank\n";
  trace.precision(12);
  std::vector<double> to_min;
  for (std::uint64_t run = 0; run < 100; run++) {
    quantum::GroverAdaptiveSearch<quantum::EmulatedGrover> gas(device, {}, Draws(1000 + run));
    double best = inf;
    for (std::uint64_t m = 0;; m++) {
      auto c = gas.ask(1);
      const double v = f(c[0].genotype);
      const search::Told told{c[0].ticket, search::scalar(v)};
      gas.tell(std::span<const search::Told>(&told, 1));
      best = std::min(best, v);
      if (run < 12) trace << run << ',' << m + 1 << ',' << gas.oracle_calls() << ',' << gas.threshold() << ',' << device.marked(gas.threshold()) + 1 << '\n';
      if (best <= device.minimum()) {
        to_min.push_back(static_cast<double>(gas.oracle_calls()));
        break;
      }
    }
  }
  std::sort(to_min.begin(), to_min.end());
  std::fprintf(stderr, "QUBO of %zu bits: minimum %.6f; oracle calls to it over 100 runs: median %.0f, mean %.0f (bound %.0f); the classical emulation evaluated %llu points\n", n, device.minimum(), to_min[to_min.size() / 2],
               std::accumulate(to_min.begin(), to_min.end(), 0.0) / static_cast<double>(to_min.size()), quantum::dh_expected_bound(N), static_cast<unsigned long long>(device.classical_evaluations()));
  std::ofstream g(out + "/gas_qubo.csv");
  g << "run,oracle_calls\n";
  for (std::size_t r = 0; r < to_min.size(); r++) g << r << ',' << to_min[r] << '\n';
}

// ---- (3) the quantum-inspired EA beside PBIL ----------------------------------------------------------------

void part_inspired(const std::string& out, calx::scheduler::ThreadPool& pool, const heuresis::benchmarks::Qubo& qubo, const Bits& optimum, double minimum) {
  using namespace heuresis;
  const std::size_t n = qubo.n(), generations = 400, individuals = 10, seeds = 20;
  const distributions::ProductBernoulli family;
  struct Track {
    std::vector<double> best;                       // best so far after each generation
    std::vector<std::vector<double>> agreement;     // P(xᵢ = x*ᵢ) after each generation (seed 0 only)
    std::uint64_t bits = 0, differing = 0;
  };
  std::vector<Track> qiea(seeds), pbil(seeds);
  calx::scheduler::parallel_for(pool, seeds, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t seed = b; seed < e; seed++) {
      // QIEA by hand: each generation's observations against the product family's reference draws.
      quantum::QuantumInspiredOptions qo;
      qo.individuals = individuals;
      qo.global_period = 20;
      const Draws draws(seed);
      quantum::QuantumInspired q(n, qo, draws);
      Track& T = qiea[seed];
      double best = inf;
      for (std::size_t g = 0; g < generations; g++) {
        const auto angles = q.angles();  // before the generation is observed
        auto c = q.ask(search::unlimited);
        std::vector<search::Told> told(c.size());
        for (std::size_t k = 0; k < c.size(); k++) {
          calx::random::Stream s = draws(quantum::QuantumInspired::observe_tag, g, k);
          const Bits reference = family.sample(quantum::observation_probabilities(angles[k]), s);
          for (std::size_t i = 0; i < n; i++) T.differing += reference[i] != c[k].genotype[i];
          T.bits += n;
          const double v = qubo.value(c[k].genotype);
          best = std::min(best, v);
          told[k] = {c[k].ticket, search::scalar(v)};
        }
        q.tell(told);
        T.best.push_back(best);
        if (seed == 0) {
          const auto p = quantum::observation_probabilities(q.angles()[0]);
          std::vector<double> a(n);
          for (std::size_t i = 0; i < n; i++) a[i] = optimum[i] ? p[i] : 1 - p[i];
          T.agreement.push_back(std::move(a));
        }
      }
      // PBIL with the same budget: λ = 10 a generation, the best one learnt at rate 0.1 (Baluja 1994).
      distributions::DistributionSearch<distributions::ProductBernoulli, distributions::Pbil> pb(family, distributions::ProductBernoulli::uniform(n), individuals, distributions::Pbil{1, 0.1}, draws);
      Track& P = pbil[seed];
      best = inf;
      for (std::size_t g = 0; g < generations; g++) {
        auto c = pb.ask(search::unlimited);
        std::vector<search::Told> told(c.size());
        for (std::size_t k = 0; k < c.size(); k++) {
          const double v = qubo.value(c[k].genotype);
          best = std::min(best, v);
          told[k] = {c[k].ticket, search::scalar(v)};
        }
        pb.tell(told);
        P.best.push_back(best);
        if (seed == 0) {
          const auto& p = pb.parameters();
          std::vector<double> a(n);
          for (std::size_t i = 0; i < n; i++) a[i] = optimum[i] ? p[i] : 1 - p[i];
          P.agreement.push_back(std::move(a));
        }
      }
    }
  });
  std::uint64_t bits = 0, differing = 0;
  std::size_t q_found = 0, p_found = 0;
  for (std::size_t s = 0; s < seeds; s++) {
    bits += qiea[s].bits;
    differing += qiea[s].differing;
    q_found += qiea[s].best.back() <= minimum + 1e-9;
    p_found += pbil[s].best.back() <= minimum + 1e-9;
  }
  std::fprintf(stderr, "QIEA's observations against the product-Bernoulli reference draws: %llu of %llu bits differ\n", static_cast<unsigned long long>(differing), static_cast<unsigned long long>(bits));
  std::fprintf(stderr, "the QUBO's minimum found in %zu of %zu QIEA runs and %zu of %zu PBIL runs (%zu evaluations)\n", q_found, seeds, p_found, seeds, generations * individuals);
  // The anytime curves, every generation to the 50th and every fifth after; the evaluations to the minimum.
  std::ofstream f(out + "/inspired_convergence.csv");
  f << "algorithm,seed,generation,evaluations,excess\n";
  f.precision(8);
  for (std::size_t s = 0; s < seeds; s++)
    for (std::size_t g = 0; g < generations; g++) {
      if (g + 1 > 50 && (g + 1) % 5) continue;
      f << "QIEA," << s << ',' << g + 1 << ',' << (g + 1) * individuals << ',' << qiea[s].best[g] - minimum << '\n';
      f << "PBIL," << s << ',' << g + 1 << ',' << (g + 1) * individuals << ',' << pbil[s].best[g] - minimum << '\n';
    }
  std::ofstream h(out + "/inspired_hitting.csv");
  h << "algorithm,seed,evaluations_to_minimum\n";
  for (std::size_t s = 0; s < seeds; s++)
    for (const Track* T : {&qiea[s], &pbil[s]}) {
      std::size_t g = 0;
      while (g < T->best.size() && T->best[g] > minimum + 1e-9) g++;
      h << (T == &qiea[s] ? "QIEA," : "PBIL,") << s << ',' << (g < T->best.size() ? static_cast<double>((g + 1) * individuals) : inf) << '\n';
    }
  std::ofstream a(out + "/inspired_agreement.json");
  a.precision(6);
  auto matrix = [&](const std::vector<std::vector<double>>& m) {
    a << "[";
    for (std::size_t g = 0; g < m.size(); g++) {
      a << (g ? "," : "") << "[";
      for (std::size_t i = 0; i < m[g].size(); i++) a << (i ? "," : "") << m[g][i];
      a << "]";
    }
    a << "]";
  };
  a << "{\n  \"bits\": " << n << ",\n  \"generations\": " << generations << ",\n  \"individuals\": " << individuals << ",\n  \"observed_bits\": " << bits << ",\n  \"differing_bits\": " << differing
    << ",\n  \"qiea_found\": " << q_found << ",\n  \"pbil_found\": " << p_found << ",\n  \"seeds\": " << seeds << ",\n  \"minimum\": " << minimum << ",\n  \"qiea\": ";
  matrix(qiea[0].agreement);
  a << ",\n  \"pbil\": ";
  matrix(pbil[0].agreement);
  a << "\n}\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: quantum <results/> [threads]\n";
    return 2;
  }
  const std::string out = argv[1];
  const std::size_t threads = argc > 2 ? std::stoul(argv[2]) : 4;
  calx::scheduler::ThreadPool pool(threads);
  const auto t0 = std::chrono::steady_clock::now();
  part_law(out);
  // A QUBO of 20 bits, its minimum by enumerating its 2²⁰ points.
  const auto qubo = heuresis::benchmarks::Qubo::random(20, 7);
  double minimum = inf;
  std::uint64_t argmin = 0;
  for (std::uint64_t i = 0; i < (std::uint64_t{1} << 20); i++) {
    const double v = qubo.value(index_bits(i, 20));
    if (v < minimum) minimum = v, argmin = i;
  }
  const Bits optimum = index_bits(argmin, 20);
  std::fprintf(stderr, "QUBO(20, seed 7): minimum %.9f at %llu\n", minimum, static_cast<unsigned long long>(argmin));
  part_adaptive(out, pool, qubo);
  part_inspired(out, pool, qubo, optimum, minimum);
  std::ofstream f(out + "/qubo.json");
  f.precision(12);
  f << "{\n  \"bits\": 20,\n  \"seed\": 7,\n  \"minimum\": " << minimum << ",\n  \"argmin\": " << argmin << ",\n  \"optimum\": [";
  for (std::size_t i = 0; i < optimum.size(); i++) f << (i ? "," : "") << static_cast<int>(optimum[i]);
  f << "]\n}\n";
  std::cerr << "total " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
  return 0;
}
