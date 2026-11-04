// The genetic algorithm as a population (README.md).
//
// A generational genetic algorithm without crossover, whose parents are drawn by fitness-proportional selection
// (heuresis::selection::roulette) and whose children copy them, is the haploid Wright–Fisher model with genic
// selection: given k carriers of a variant of fitness 1 + s among M individuals of fitness 1, the number of
// carriers among the children is Binomial(M, ψ), ψ = k(1 + s)/(k(1 + s) + M − k) (Fisher 1930; Wright 1931).
// This program shows it by simulation with heuresis's operators, against the evolution domain's Wright–Fisher
// replicates and Kimura's diffusion (evolution::population), and against the exact Markov chain:
//
//   law:       one generation's law of the carriers' count under roulette, universal sampling and tournaments;
//   fixation:  the probability that one advantageous (or deleterious) variant fixes, across M and 2Ms;
//   drift:     the decay of the expected heterozygosity 2p(1 − p) of a neutral label, by scheme, against the
//              exact rate 1 − Σ_j q_j² (q_j the probability that a child's parent is j);
//   takeover:  the time for the best to fill the population under tournaments, truncation and roulette;
//   balance:   OneMax under mutation and selection: the equilibrium mean against the infinite population.
//
// Every number of the record is written to <results>/ by this program.
//
//   ga_population <results directory> [law|fixation|drift|takeover|balance|all]

#include <heuresis/benchmarks/pseudo_boolean.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/selection/selection.hpp>
#include <heuresis/variation/variation.hpp>

#include <evolution/population/wright_fisher.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

namespace hs = heuresis::search;
namespace sel = heuresis::selection;
namespace hv = heuresis::variation;
namespace hg = heuresis::genetic;
namespace ep = evolution::population;

namespace {

std::filesystem::path out_dir;
using Pop = std::vector<std::uint8_t>;  // 1: carries the variant (or the label, or is the best)

double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

constexpr std::uint32_t LAW = hs::tag("study.law"), FIX = hs::tag("study.fixation"), DRIFT = hs::tag("study.drift"), ARRANGE = hs::tag("study.arrange"),
                        TAKE = hs::tag("study.takeover"), BALANCE = hs::tag("study.balance");

enum class Scheme { roulette, sus, tournament2, tournament4, tournament2_fresh_ties, truncation_half, truncation_quarter };
const char* name_of(Scheme s) {
  switch (s) {
    case Scheme::roulette: return "roulette";
    case Scheme::sus: return "universal sampling";
    case Scheme::tournament2: return "tournament 2";
    case Scheme::tournament4: return "tournament 4";
    case Scheme::tournament2_fresh_ties: return "tournament 2 fresh ties";
    case Scheme::truncation_half: return "truncation 1/2";
    case Scheme::truncation_quarter: return "truncation 1/4";
  }
  return "";
}

// The M parents of the next generation from fitnesses w (larger better), by the scheme; every draw from r.
// Tournaments and truncation read the population through the order w_i > w_j, ties to the earlier index (the
// library's rule), or, for the fresh-ties tournament, through keys drawn anew for each tournament.
std::vector<std::size_t> parents(Scheme scheme, const std::vector<double>& w, calx::random::Stream& r) {
  const std::size_t M = w.size();
  std::vector<std::size_t> p(M);
  const sel::Less better = [&w](std::size_t i, std::size_t j) { return w[i] > w[j]; };
  switch (scheme) {
    case Scheme::roulette:
      for (auto& x : p) x = sel::roulette(w, r);
      break;
    case Scheme::sus: p = sel::universal_sampling(w, M, r); break;
    case Scheme::tournament2:
    case Scheme::tournament4:
      for (auto& x : p) x = sel::tournament(M, scheme == Scheme::tournament2 ? 2 : 4, better, r);
      break;
    case Scheme::tournament2_fresh_ties: {
      std::vector<double> key(M);
      const sel::Less keyed = [&](std::size_t i, std::size_t j) { return w[i] > w[j] || (w[i] == w[j] && key[i] < key[j]); };
      for (auto& x : p) {
        for (auto& k : key) k = r.uniform();
        x = sel::tournament(M, 2, keyed, r);
      }
      break;
    }
    case Scheme::truncation_half:
    case Scheme::truncation_quarter: {
      const std::size_t mu = scheme == Scheme::truncation_half ? M / 2 : M / 4;
      const auto best = sel::truncation(M, mu, better);
      for (auto& x : p) x = best[r.below(mu)];
      break;
    }
  }
  return p;
}

// Fitnesses 1 + s·g of a population of two kinds.
std::vector<double> two_kinds(const Pop& g, double s) {
  std::vector<double> w(g.size());
  for (std::size_t i = 0; i < g.size(); i++) w[i] = 1 + s * g[i];
  return w;
}

// One generation: children copy their parents.
std::size_t next_generation(Scheme scheme, Pop& g, double s, calx::random::Stream& r) {
  const auto p = parents(scheme, two_kinds(g, s), r);
  Pop next(g.size());
  std::size_t k = 0;
  for (std::size_t c = 0; c < g.size(); c++) k += (next[c] = g[p[c]]);
  g = std::move(next);
  return k;
}

double log_choose(double n, double k) {
  using calx::special_functions::log_gamma;
  return log_gamma(n + 1) - log_gamma(k + 1) - log_gamma(n - k + 1);
}
double binomial_pmf(std::size_t n, double p, std::size_t j) {
  if (p <= 0) return j == 0;
  if (p >= 1) return j == n;
  return std::exp(log_choose(double(n), double(j)) + double(j) * std::log(p) + double(n - j) * std::log1p(-p));
}

// ---- one generation's law --------------------------------------------------------------------------------------

void law(calx::scheduler::ThreadPool& pool) {
  const std::size_t M = 40, k = 10, G = 200000;
  const double s = 0.5;
  const double psi = double(k) * (1 + s) / (double(k) * (1 + s) + double(M - k));
  const double p0 = double(k) / double(M);
  FILE* f = std::fopen((out_dir / "one_generation.csv").string().c_str(), "w");
  std::fprintf(f, "scheme,count,frequency,theory\n");
  FILE* fs = std::fopen((out_dir / "one_generation_summary.csv").string().c_str(), "w");
  std::fprintf(fs, "scheme,M,k,s,generations,mean,variance,theory_mean,theory_variance,total_variation\n");
  for (Scheme scheme : {Scheme::roulette, Scheme::sus, Scheme::tournament2}) {
    const std::size_t chunks = 64;
    std::vector<std::vector<double>> parts(chunks, std::vector<double>(M + 1, 0));
    calx::scheduler::parallel_for(pool, chunks, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t c = b; c < e; c++)
        for (std::size_t g = c * G / chunks; g < (c + 1) * G / chunks; g++) {
          Pop pop(M, 0);
          std::fill(pop.begin(), pop.begin() + k, 1);
          calx::random::Stream r = hs::Draws(11)(LAW, std::uint64_t(scheme), g);
          parts[c][next_generation(scheme, pop, s, r)] += 1;
        }
    });
    std::vector<double> h(M + 1, 0);
    for (const auto& part : parts)
      for (std::size_t j = 0; j <= M; j++) h[j] += part[j];
    // The theory: roulette Binomial(M, ψ); a tournament of two Binomial(M, 1 − (1 − p)²); universal sampling ⌊Mψ⌋
    // or ⌈Mψ⌉ copies, the latter with probability Mψ − ⌊Mψ⌋.
    std::vector<double> t(M + 1, 0);
    if (scheme == Scheme::roulette)
      for (std::size_t j = 0; j <= M; j++) t[j] = binomial_pmf(M, psi, j);
    else if (scheme == Scheme::tournament2)
      for (std::size_t j = 0; j <= M; j++) t[j] = binomial_pmf(M, 1 - (1 - p0) * (1 - p0), j);
    else {
      const double e = double(M) * psi, fl = std::floor(e);
      t[std::size_t(fl)] = 1 - (e - fl);
      if (e > fl) t[std::size_t(fl) + 1] = e - fl;
    }
    double mean = 0, var = 0, tm = 0, tv = 0, tvd = 0;
    for (std::size_t j = 0; j <= M; j++) {
      h[j] /= double(G);
      mean += double(j) * h[j];
      tm += double(j) * t[j];
      tvd += 0.5 * std::fabs(h[j] - t[j]);
    }
    for (std::size_t j = 0; j <= M; j++) var += (double(j) - mean) * (double(j) - mean) * h[j], tv += (double(j) - tm) * (double(j) - tm) * t[j];
    for (std::size_t j = 0; j <= M; j++) std::fprintf(f, "%s,%zu,%.6f,%.6f\n", name_of(scheme), j, h[j], t[j]);
    std::fprintf(fs, "%s,%zu,%zu,%.2f,%zu,%.5f,%.5f,%.5f,%.5f,%.5f\n", name_of(scheme), M, k, s, G, mean, var, tm, tv, tvd);
    std::printf("law %s: mean %.4f (theory %.4f), variance %.4f (theory %.4f), TV %.4f\n", name_of(scheme), mean, tm, var, tv, tvd);
  }
  std::fclose(f);
  std::fclose(fs);
}

// ---- fixation ------------------------------------------------------------------------------------------------

// The exact probability that one copy fixes in the haploid Wright–Fisher chain with genic selection: u solves
// u_i = Σ_j P_ij u_j on the transient states 1 … M − 1, u_0 = 0, u_M = 1, P_ij = Binomial(M, ψ_i)(j).
double exact_fixation(std::size_t M, double s) {
  const std::size_t n = M - 1;
  calx::linalg::Matrix A = calx::linalg::Matrix::identity(n), B(n);
  for (std::size_t i = 1; i < M; i++) {
    const double psi = double(i) * (1 + s) / (double(i) * (1 + s) + double(M - i));
    for (std::size_t j = 1; j < M; j++) A(i - 1, j - 1) -= binomial_pmf(M, psi, j);
    B(i - 1, 0) = binomial_pmf(M, psi, M);
  }
  const auto X = calx::linalg::solve(A, B);
  return X(0, 0);
}

void fixation(calx::scheduler::ThreadPool& pool) {
  const std::vector<std::size_t> sizes{20, 50, 100, 200};
  const std::vector<double> scaled{-2, -1, 0, 1, 2, 4, 8, 16};  // S = 2Ms
  FILE* f = std::fopen((out_dir / "fixation.csv").string().c_str(), "w");
  std::fprintf(f, "M,S,s,replicates,fixed_ga,u_ga,se_ga,mean_fixation_generations,fixed_domain,u_domain,u_kimura,u_diffusion,u_exact,wall_s\n");
  for (std::size_t mi = 0; mi < sizes.size(); mi++)
    for (std::size_t si = 0; si < scaled.size(); si++) {
      const std::size_t M = sizes[mi];
      const double S = scaled[si], s = S / (2.0 * double(M));
      const std::size_t R = S < 0 ? 40000 : 20000;
      const std::size_t chunks = 128;
      std::vector<std::size_t> fixed(chunks, 0);
      std::vector<double> generations(chunks, 0);
      const hs::Draws draws(1000 + 100 * mi + si);
      const auto t0 = std::chrono::steady_clock::now();
      calx::scheduler::parallel_for(pool, chunks, 1, [&](std::size_t b, std::size_t e, std::size_t) {
        for (std::size_t c = b; c < e; c++)
          for (std::size_t rep = c * R / chunks; rep < (c + 1) * R / chunks; rep++) {
            Pop pop(M, 0);
            pop[0] = 1;
            std::size_t k = 1, g = 0;
            while (k > 0 && k < M) {
              calx::random::Stream r = draws(FIX, g, rep);
              k = next_generation(Scheme::roulette, pop, s, r);
              g++;
            }
            if (k == M) fixed[c]++, generations[c] += double(g);
          }
      });
      const double wall = seconds_since(t0);
      const std::size_t nf = std::accumulate(fixed.begin(), fixed.end(), std::size_t{0});
      const double gf = std::accumulate(generations.begin(), generations.end(), 0.0);
      const double u = double(nf) / double(R);
      // The domain's Wright–Fisher replicates: N = M/2 diploids (2N = M gametes), genic selection h = ½, the
      // allele's advantage s/2·2 = s (fitnesses 1 + 2s, 1 + s, 1), from one copy.
      ep::WrightFisherSpec w;
      w.individuals = M / 2;
      w.s = 2 * s;
      w.h = 0.5;
      w.p0 = 1.0 / double(M);
      w.generations = 200000;
      const auto dom = ep::wright_fisher_replicates(w, 5000 + 100 * mi + si, R, pool, w.generations);
      const double kim = ep::kimura_fixation(double(M) / 2, s, 1.0 / double(M));
      const double diff = ep::fixation_probability(double(M) / 2, 2 * s, 0.5, 1.0 / double(M));
      const double ex = exact_fixation(M, s);
      std::fprintf(f, "%zu,%.0f,%.6f,%zu,%zu,%.6e,%.6e,%.2f,%zu,%.6e,%.6e,%.6e,%.6e,%.2f\n", M, S, s, R, nf, u, std::sqrt(u * (1 - u) / double(R)), nf ? gf / double(nf) : 0.0, dom.fixed,
                   dom.fixation(), kim, diff, ex, wall);
      std::fflush(f);
      std::printf("fixation M=%zu S=%+.0f: GA %.5f, domain %.5f, Kimura %.5f, exact %.5f (%.1f s)\n", M, S, u, dom.fixation(), kim, ex, wall);
    }
  std::fclose(f);
}

// ---- drift -------------------------------------------------------------------------------------------------

// The coalescence probability per generation c = Σ_j q_j², q_j the probability that a child's parent is j, for
// neutral fitness: roulette 1/M; universal sampling 0 (each one copy); truncation to μ 1/μ; a tournament of two
// with ties to the earlier index (q_j = (2(M − j) − 1)/M², j = 0 … M − 1) (4M² − 1)/(3M³); with fresh ties 1/M.
// Then E[2p'(1 − p')] = (1 − c)·2p(1 − p) exactly, the arrangement of the labels being exchangeable.
double coalescence(Scheme scheme, std::size_t M) {
  const double m = double(M);
  switch (scheme) {
    case Scheme::roulette: return 1 / m;
    case Scheme::sus: return 0;
    case Scheme::tournament2: return (4 * m * m - 1) / (3 * m * m * m);
    case Scheme::tournament2_fresh_ties: return 1 / m;
    case Scheme::truncation_half: return 2 / m;
    default: return 0;
  }
}

void drift(calx::scheduler::ThreadPool& pool) {
  const std::size_t M = 50, R = 4000, T = 300;
  FILE* f = std::fopen((out_dir / "drift.csv").string().c_str(), "w");
  std::fprintf(f, "scheme,generation,mean_H,theory_H,coalescence\n");
  for (Scheme scheme : {Scheme::roulette, Scheme::sus, Scheme::tournament2, Scheme::tournament2_fresh_ties, Scheme::truncation_half}) {
    const std::size_t chunks = 64;
    std::vector<std::vector<double>> H(chunks, std::vector<double>(T + 1, 0));
    calx::scheduler::parallel_for(pool, chunks, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t c = b; c < e; c++)
        for (std::size_t rep = c * R / chunks; rep < (c + 1) * R / chunks; rep++) {
          // Half the population carries the label, at positions drawn uniformly (Fisher–Yates).
          Pop pop(M, 0);
          std::fill(pop.begin(), pop.begin() + M / 2, 1);
          calx::random::Stream a = hs::Draws(21)(ARRANGE, std::uint64_t(scheme), rep);
          for (std::size_t i = M - 1; i > 0; i--) std::swap(pop[i], pop[a.below(i + 1)]);
          std::size_t k = M / 2;
          for (std::size_t g = 0; g <= T; g++) {
            const double p = double(k) / double(M);
            H[c][g] += 2 * p * (1 - p);
            if (g == T) break;
            calx::random::Stream r = hs::Draws(22)(DRIFT, g, rep * 8 + std::uint64_t(scheme));
            k = next_generation(scheme, pop, 0.0, r);
          }
        }
    });
    const double cc = coalescence(scheme, M);
    for (std::size_t g = 0; g <= T; g++) {
      double h = 0;
      for (const auto& part : H) h += part[g];
      std::fprintf(f, "%s,%zu,%.6f,%.6f,%.6f\n", name_of(scheme), g, h / double(R), 0.5 * std::pow(1 - cc, double(g)), cc);
    }
    std::printf("drift %s done\n", name_of(scheme));
  }
  // The domain's neutral Wright–Fisher replicates: N = 25 diploids, 2N = 50 gametes, from p = ½.
  ep::WrightFisherSpec w;
  w.individuals = M / 2;
  w.p0 = 0.5;
  w.generations = T;
  const auto dom = ep::wright_fisher_replicates(w, 23, R, pool, 1);
  for (std::size_t g = 0; g <= T; g++) std::fprintf(f, "domain Wright-Fisher,%zu,%.6f,%.6f,%.6f\n", g, dom.h_s[g], 0.5 * ep::heterozygosity_decay(double(M) / 2, double(g)), 1 / double(M));
  std::fclose(f);
}

// ---- takeover ----------------------------------------------------------------------------------------------

void takeover(calx::scheduler::ThreadPool& pool) {
  const std::size_t M = 500, R = 1000, T = 200;
  struct Case {
    Scheme scheme;
    double s;
    std::string label;
  };
  const std::vector<Case> cases{{Scheme::tournament2, 1.0, "tournament 2"},       {Scheme::tournament4, 1.0, "tournament 4"},
                                {Scheme::truncation_half, 1.0, "truncation 1/2"}, {Scheme::truncation_quarter, 1.0, "truncation 1/4"},
                                {Scheme::roulette, 1.0, "roulette s = 1"},        {Scheme::roulette, 4.0, "roulette s = 4"},
                                {Scheme::tournament2, 4.0, "tournament 2 s = 4"}, {Scheme::truncation_half, 4.0, "truncation 1/2 s = 4"}};
  FILE* f = std::fopen((out_dir / "takeover.csv").string().c_str(), "w");
  std::fprintf(f, "case,generation,mean_fraction_fixed,deterministic\n");
  FILE* ft = std::fopen((out_dir / "takeover_times.csv").string().c_str(), "w");
  std::fprintf(ft, "case,replicates,fixed,mean_time,sd_time,deterministic_time,closed_form\n");
  std::vector<std::vector<std::size_t>> trajectories_t2(2);  // tournament 2 at s = 1 and s = 4, to compare
  std::vector<std::vector<std::size_t>> trajectories_tr(2);
  for (std::size_t ci = 0; ci < cases.size(); ci++) {
    const Case& cs = cases[ci];
    std::vector<std::vector<std::size_t>> counts(R);  // each replicate's count by generation
    calx::scheduler::parallel_for(pool, R, 8, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t rep = b; rep < e; rep++) {
        Pop pop(M, 0);
        pop[0] = 1;
        std::size_t k = 1;
        counts[rep].push_back(k);
        for (std::size_t g = 0; g < T && k > 0 && k < M; g++) {
          // The same streams for every case of one scheme, so that s's effect (or none) is seen draw for draw.
          calx::random::Stream r = hs::Draws(31)(TAKE, g, rep * 8 + std::uint64_t(cs.scheme));
          k = next_generation(cs.scheme, pop, cs.s, r);
          counts[rep].push_back(k);
        }
      }
    });
    // Order-based schemes see s only through the order: the runs at s = 1 and s = 4 must agree count for count.
    if (cs.scheme == Scheme::tournament2) {
      auto& slot = trajectories_t2[cs.s == 1.0 ? 0 : 1];
      for (const auto& c : counts) slot.insert(slot.end(), c.begin(), c.end());
    }
    if (cs.scheme == Scheme::truncation_half) {
      auto& slot = trajectories_tr[cs.s == 1.0 ? 0 : 1];
      for (const auto& c : counts) slot.insert(slot.end(), c.begin(), c.end());
    }
    std::vector<double> mean(T + 1, 0);
    std::size_t fixed = 0;
    double sum = 0, sum2 = 0;
    for (const auto& c : counts) {
      if (c.back() != M) continue;
      fixed++;
      const double t = double(c.size() - 1);
      sum += t, sum2 += t * t;
      for (std::size_t g = 0; g <= T; g++) mean[g] += double(g < c.size() ? c[g] : M) / double(M);
    }
    // The deterministic recursion of the best's share from 1/M: tournament 1 − (1 − p)^k, truncation min(1, p/τ),
    // roulette p(1 + s)/(1 + sp); its takeover time the first generation with fewer than ½ others expected.
    auto step = [&](double p) {
      switch (cs.scheme) {
        case Scheme::tournament2: return 1 - (1 - p) * (1 - p);
        case Scheme::tournament4: return 1 - std::pow(1 - p, 4);
        case Scheme::truncation_half: return std::min(1.0, 2 * p);
        case Scheme::truncation_quarter: return std::min(1.0, 4 * p);
        default: return p * (1 + cs.s) / (1 + cs.s * p);
      }
    };
    std::vector<double> det(T + 1);
    det[0] = 1.0 / double(M);
    double det_time = -1;
    for (std::size_t g = 0; g <= T; g++) {
      if (g) det[g] = step(det[g - 1]);
      if (det_time < 0 && double(M) * (1 - det[g]) < 0.5) det_time = double(g);
    }
    // Goldberg and Deb's (1991) closed forms: tournaments (ln M + ln ln M)/ln k, truncation ln M/ln(1/τ);
    // roulette between two values, from odds 1/(M − 1) to M − 1: 2 ln(M − 1)/ln(1 + s).
    double closed = 0;
    const double m = double(M);
    switch (cs.scheme) {
      case Scheme::tournament2: closed = (std::log(m) + std::log(std::log(m))) / std::log(2.0); break;
      case Scheme::tournament4: closed = (std::log(m) + std::log(std::log(m))) / std::log(4.0); break;
      case Scheme::truncation_half: closed = std::log(m) / std::log(2.0); break;
      case Scheme::truncation_quarter: closed = std::log(m) / std::log(4.0); break;
      default: closed = 2 * std::log(m - 1) / std::log(1 + cs.s);
    }
    for (std::size_t g = 0; g <= T; g++) std::fprintf(f, "%s,%zu,%.6f,%.6f\n", cs.label.c_str(), g, fixed ? mean[g] / double(fixed) : 0.0, det[g]);
    const double mt = fixed ? sum / double(fixed) : 0, sd = fixed > 1 ? std::sqrt((sum2 - sum * sum / double(fixed)) / double(fixed - 1)) : 0;
    std::fprintf(ft, "%s,%zu,%zu,%.3f,%.3f,%.0f,%.3f\n", cs.label.c_str(), R, fixed, mt, sd, det_time, closed);
    std::printf("takeover %s: %zu of %zu fixed, mean %.2f generations (deterministic %.0f, closed form %.2f)\n", cs.label.c_str(), fixed, R, mt, det_time, closed);
  }
  std::fclose(f);
  std::fclose(ft);
  FILE* fj = std::fopen((out_dir / "takeover_invariance.json").string().c_str(), "w");
  std::fprintf(fj, "{\"tournament 2, s = 1 and s = 4 identical\": %s, \"truncation 1/2, s = 1 and s = 4 identical\": %s}\n", trajectories_t2[0] == trajectories_t2[1] ? "true" : "false",
               trajectories_tr[0] == trajectories_tr[1] ? "true" : "false");
  std::fclose(fj);
}

// ---- mutation–selection balance on OneMax ------------------------------------------------------------------

// The infinite population on the classes 0 … n of the number of ones: selection (roulette: x_k ← k x_k / Σ j x_j;
// a tournament of two: the winner's class law, P(≥ c) = 1 − (1 − G(c))², G(c) = Σ_{j ≥ c} x_j), then standard bit
// mutation at rate μ (from class k: a ~ Bin(k, μ) ones lost, b ~ Bin(n − k, μ) gained), iterated to its fixed point
// from the binomial class law of random strings. Exact for an infinite asexual population (Eigen 1971's
// quasispecies, on the classes of OneMax).
double infinite_population(std::size_t n, double mu, bool tournament) {
  std::vector<std::vector<double>> P(n + 1, std::vector<double>(n + 1, 0));
  for (std::size_t k = 0; k <= n; k++)
    for (std::size_t a = 0; a <= k; a++) {
      const double pa = binomial_pmf(k, mu, a);
      if (pa < 1e-300) continue;
      for (std::size_t b = 0; b <= n - k; b++) P[k][k - a + b] += pa * binomial_pmf(n - k, mu, b);
    }
  std::vector<double> x(n + 1);
  for (std::size_t k = 0; k <= n; k++) x[k] = binomial_pmf(n, 0.5, k);
  double last = -1;
  for (int it = 0; it < 200000; it++) {
    std::vector<double> y(n + 1, 0);
    if (!tournament) {
      double z = 0;
      for (std::size_t k = 0; k <= n; k++) z += double(k) * x[k];
      for (std::size_t k = 0; k <= n; k++) y[k] = double(k) * x[k] / z;
    } else {
      std::vector<double> G(n + 2, 0);
      for (std::size_t k = n + 1; k-- > 0;) G[k] = G[k + 1] + x[k];
      for (std::size_t k = 0; k <= n; k++) y[k] = (1 - (1 - G[k]) * (1 - G[k])) - (1 - (1 - G[k + 1]) * (1 - G[k + 1]));
    }
    std::vector<double> z(n + 1, 0);
    for (std::size_t k = 0; k <= n; k++)
      if (y[k] > 0)
        for (std::size_t j = 0; j <= n; j++) z[j] += y[k] * P[k][j];
    x = z;
    double m = 0;
    for (std::size_t k = 0; k <= n; k++) m += double(k) * x[k];
    if (std::fabs(m - last) < 1e-13) break;
    last = m;
  }
  double m = 0;
  for (std::size_t k = 0; k <= n; k++) m += double(k) * x[k];
  return m / double(n);
}

void balance(calx::scheduler::ThreadPool& pool) {
  const std::size_t n = 100, T = 2000, from = 1000, seeds = 4;
  const std::vector<double> chis{0.25, 0.5, 1, 2, 4, 8};
  const std::vector<std::size_t> sizes{20, 100, 500};
  struct Job {
    bool tournament;
    std::size_t M, ci;
    std::uint64_t seed;
  };
  std::vector<Job> jobs;
  for (bool t : {false, true})
    for (std::size_t M : sizes)
      for (std::size_t ci = 0; ci < chis.size(); ci++)
        for (std::uint64_t sd = 1; sd <= seeds; sd++) jobs.push_back({t, M, ci, sd});
  std::vector<std::vector<double>> traj(jobs.size());
  const auto t0 = std::chrono::steady_clock::now();
  calx::scheduler::parallel_for(pool, jobs.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t j = b; j < e; j++) {
      const Job job = jobs[j];
      const double mu = chis[job.ci] / double(n);
      auto& tr = traj[j];
      if (!job.tournament) {
        // The GA by roulette, built from heuresis's parts: weights the number of ones, children copied and
        // mutated by standard bit mutation, generational, no elites.
        std::vector<hv::Bits> pop(job.M);
        for (std::size_t i = 0; i < job.M; i++) {
          calx::random::Stream s = hs::Draws(job.seed)(BALANCE, 0, i);
          pop[i] = hv::random_bits(n, s);
        }
        for (std::size_t g = 0; g < T; g++) {
          std::vector<double> w(job.M);
          double mean = 0;
          for (std::size_t i = 0; i < job.M; i++) mean += (w[i] = double(n) - heuresis::benchmarks::onemax(pop[i]));
          tr.push_back(mean / double(job.M * n));
          calx::random::Stream r = hs::Draws(job.seed)(BALANCE, g + 1, 1u << 20);
          std::vector<hv::Bits> next(job.M);
          for (std::size_t c = 0; c < job.M; c++) {
            next[c] = pop[sel::roulette(w, r)];
            hv::standard_bit_mutation(next[c], mu, r);
          }
          pop = std::move(next);
        }
      } else {
        // The library's genetic algorithm: tournaments of two, no crossover, standard bit mutation, generational
        // replacement without elites, minimising OneMax's zeros.
        hg::GaOptions o;
        o.mu = job.M;
        o.crossover_rate = 0;
        o.tournament = 2;
        o.elites = 0;
        hg::GeneticAlgorithm<hv::Bits> ga({hg::random_bits(n), hg::bit_flip(mu), {}}, o, hs::Draws(job.seed + 100));
        for (std::size_t g = 0; g <= T; g++) {
          auto cs = ga.ask(hs::unlimited);
          std::vector<hs::Told> told;
          for (const auto& c : cs) told.push_back({c.ticket, hs::scalar(heuresis::benchmarks::onemax(c.genotype))});
          ga.tell(told);
          double mean = 0;
          for (const auto& x : ga.population()) mean += double(n) - x.evaluation.value();
          tr.push_back(mean / double(job.M * n));
        }
        tr.resize(T);
      }
    }
  });
  const double wall = seconds_since(t0);
  FILE* f = std::fopen((out_dir / "balance.csv").string().c_str(), "w");
  std::fprintf(f, "scheme,M,chi,seed,equilibrium\n");
  FILE* ftr = std::fopen((out_dir / "balance_trajectories.csv").string().c_str(), "w");
  std::fprintf(ftr, "scheme,M,chi,generation,mean_fraction\n");
  for (std::size_t j = 0; j < jobs.size(); j++) {
    const Job job = jobs[j];
    const double eq = std::accumulate(traj[j].begin() + from, traj[j].end(), 0.0) / double(T - from);
    std::fprintf(f, "%s,%zu,%.2f,%llu,%.6f\n", job.tournament ? "tournament 2" : "roulette", job.M, chis[job.ci], static_cast<unsigned long long>(job.seed), eq);
    if (job.seed == 1 && job.M == 100)
      for (std::size_t g = 0; g < T; g += 4) std::fprintf(ftr, "%s,%zu,%.2f,%zu,%.5f\n", job.tournament ? "tournament 2" : "roulette", job.M, chis[job.ci], g, traj[j][g]);
  }
  std::fclose(f);
  std::fclose(ftr);
  FILE* fi = std::fopen((out_dir / "balance_infinite.csv").string().c_str(), "w");
  std::fprintf(fi, "scheme,chi,infinite,linkage_equilibrium\n");
  for (bool t : {false, true})
    for (double chi : {0.1, 0.15, 0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.6, 8.0}) {
      const double inf = infinite_population(n, chi / double(n), t);
      // For roulette, with the classes' variance binomial (independent loci): Var/k̄ = μ(2k̄ − n) gives
      // f = (1 + χ)/(1 + 2χ).
      std::fprintf(fi, "%s,%.2f,%.6f,%.6f\n", t ? "tournament 2" : "roulette", chi, inf, t ? NAN : (1 + chi) / (1 + 2 * chi));
    }
  std::fclose(fi);
  std::printf("balance: %zu runs in %.1f s\n", jobs.size(), wall);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: ga_population <results directory> [law|fixation|drift|takeover|balance|all]\n");
    return 2;
  }
  out_dir = argv[1];
  std::filesystem::create_directories(out_dir);
  const std::string mode = argc > 2 ? argv[2] : "all";
  calx::scheduler::ThreadPool pool(std::max(1u, std::thread::hardware_concurrency()));
  const auto t0 = std::chrono::steady_clock::now();
  if (mode == "law" || mode == "all") law(pool);
  if (mode == "drift" || mode == "all") drift(pool);
  if (mode == "takeover" || mode == "all") takeover(pool);
  if (mode == "balance" || mode == "all") balance(pool);
  if (mode == "fixation" || mode == "all") fixation(pool);
  std::printf("total %.1f s\n", seconds_since(t0));
  return 0;
}
