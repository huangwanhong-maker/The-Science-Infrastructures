// Wright–Fisher populations forward in time: the frequency of an allele A₁ (against A₂) in D
// demes of N diploid individuals, under drift, selection, mutation and migration; trajectories
// for drawing, and many replicates in parallel; with the diffusion results they are checked
// against.
//
// One generation in each deme, from the frequency p of A₁ among its 2N gametes:
//
//   1. selection on the zygotes formed at random, with fitnesses 1 + s, 1 + hs, 1 for
//      A₁A₁, A₁A₂, A₂A₂ (h the dominance of A₁): p ← (p²(1 + s) + pq(1 + hs)) / w̄;
//   2. mutation A₁ → A₂ at rate μ and A₂ → A₁ at rate ν: p ← p(1 − μ) + (1 − p)ν;
//   3. migration, the island model (Wright 1931, Genetics 16, 97–159): a fraction m of each
//      deme's gametes comes from the pool of all demes, p ← (1 − m)p + m p̄;
//   4. drift: 2N gametes drawn, p ← Binomial(2N, p)/(2N) (Fisher 1930; Wright 1931).
//
// The diffusion approximation of this chain (Kimura 1962, Genetics 47, 713–719) gives the
// probability that A₁, at frequency p, is fixed:
//
//   u(p) = ∫₀ᵖ G(x) dx / ∫₀¹ G(x) dx,   G(x) = exp(−2Ns(2hx + (1 − 2h)x²)),
//
// which for genic selection (h = ½, the allele's own advantage s_K = s/2) is Kimura's
// u(p) = (1 − e^{−4N s_K p}) / (1 − e^{−4N s_K}), and p without selection. Without selection
// and mutation the expected heterozygosity 2pq decays exactly as (1 − 1/(2N))^t, and among
// many islands the equilibrium F_ST = 1 − H_S/H_T is about 1/(1 + 4Nm), exactly
// (1 − m)²/(2N − (2N − 1)(1 − m)²) for identity by descent with an infinite pool (Wright 1931).
//
// Draws: replicate r uses the stream (seed, WRIGHT_FISHER, r), the demes in order each
// generation, so a replicate is the same whichever thread runs it; replicates are run as a
// calx task graph of shards of replicates, each summing its own, merged in shard order: the
// summary is the same on any number of threads.
//
// Complexity: O(D) per generation (one binomial variate per deme, O(1) expected by Hörmann's
// BTRS), O(R · T · D) for R replicates of T generations. A replicate without mutation stops when
// every deme has fixed or lost A₁ and migration cannot change it.

#pragma once

#include <calx/quadrature/gauss_kronrod.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/task_graph.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace evolution::population {

inline constexpr std::uint32_t WRIGHT_FISHER = 0x57524648;

struct WrightFisherSpec {
  std::size_t individuals = 100;   // N, diploid, in each deme: 2N gametes
  std::size_t demes = 1;
  double s = 0, h = 0.5;           // fitnesses 1 + s, 1 + hs, 1 of A₁A₁, A₁A₂, A₂A₂
  double mu = 0, nu = 0;           // mutation A₁ → A₂ and A₂ → A₁, per gamete and generation
  double m = 0;                    // the fraction of each deme's gametes from the pool of all demes
  double p0 = 0.5;                 // the starting frequency of A₁ in every deme
  std::size_t generations = 1000;  // the most generations run
};

namespace detail {

inline void check(const WrightFisherSpec& w) {
  if (w.individuals == 0 || w.demes == 0) throw std::invalid_argument("evolution::population::wright_fisher: N and the demes are positive");
  if (!(w.s > -1) || !(1 + w.h * w.s > 0)) throw std::invalid_argument("evolution::population::wright_fisher: fitnesses are positive");
  if (w.mu < 0 || w.mu > 1 || w.nu < 0 || w.nu > 1 || w.m < 0 || w.m > 1 || w.p0 < 0 || w.p0 > 1)
    throw std::invalid_argument("evolution::population::wright_fisher: rates and frequencies in [0, 1]");
}

// One generation of every deme, in place; returns whether nothing can change any more.
inline bool generation(const WrightFisherSpec& w, std::vector<double>& p, std::vector<double>& next, calx::random::Stream& r) {
  const std::size_t D = p.size();
  const double two_n = 2.0 * static_cast<double>(w.individuals);
  double pool = 0;
  for (std::size_t d = 0; d < D; d++) {
    double x = p[d], y = 1 - x;
    if (w.s != 0) {
      const double a = x * x * (1 + w.s), b = x * y * (1 + w.h * w.s), c = y * y;
      x = (a + b) / (a + 2 * b + c);
    }
    x = x * (1 - w.mu) + (1 - x) * w.nu;
    next[d] = x;
    pool += x;
  }
  pool /= static_cast<double>(D);
  bool settled = w.mu == 0 && w.nu == 0;
  for (std::size_t d = 0; d < D; d++) {
    const double x = std::clamp((1 - w.m) * next[d] + w.m * pool, 0.0, 1.0);
    p[d] = static_cast<double>(r.binomial(static_cast<std::int64_t>(two_n), x)) / two_n;
    settled = settled && (p[d] == 0 || p[d] == 1) && (w.m == 0 || p[d] == p[0]);
  }
  return settled;
}

}  // namespace detail

// The frequencies of every deme, generation by generation (generation-major, from generation
// 0), until the spec's last generation or until nothing can change.
struct Trajectory {
  std::size_t demes = 0;
  std::vector<double> p;  // (generations + 1) × demes
  std::size_t generations() const { return demes ? p.size() / demes - 1 : 0; }
  double at(std::size_t generation, std::size_t deme) const { return p[generation * demes + deme]; }
};

inline Trajectory wright_fisher_trajectory(const WrightFisherSpec& w, std::uint64_t seed, std::uint32_t replicate = 0) {
  detail::check(w);
  calx::random::Stream r(seed, WRIGHT_FISHER, replicate, 0);
  Trajectory t;
  t.demes = w.demes;
  std::vector<double> p(w.demes, w.p0), next(w.demes);
  t.p.reserve((w.generations + 1) * w.demes);
  t.p.insert(t.p.end(), p.begin(), p.end());
  for (std::size_t g = 0; g < w.generations; g++) {
    const bool settled = detail::generation(w, p, next, r);
    t.p.insert(t.p.end(), p.begin(), p.end());
    if (settled) break;
  }
  return t;
}

// Means over replicates, at the generations 0, every, 2·every, … ≤ the spec's last.
struct WrightFisherSummary {
  std::size_t replicates = 0;
  std::size_t fixed = 0, lost = 0;   // replicates settled with A₁ in every deme, or in none
  double mean_settling = 0;          // the mean generation at which they settled
  std::size_t every = 1;
  std::vector<double> mean_p;        // of the mean over demes
  std::vector<double> h_s, h_s_sq;   // of H_S, the mean over demes of 2pq, and of its square
  std::vector<double> h_t;           // of H_T = 2 p̄ (1 − p̄), p̄ the mean over demes
  calx::scheduler::Trace trace;      // where and when each shard ran

  double fixation() const { return replicates ? static_cast<double>(fixed) / static_cast<double>(replicates) : 0; }
  // 1 − Σ H_S / Σ H_T over the records [first, last): F_ST as a ratio of averages.
  double fst(std::size_t first, std::size_t last) const {
    double hs = 0, ht = 0;
    for (std::size_t k = first; k < std::min(last, h_t.size()); k++) {
      hs += h_s[k];
      ht += h_t[k];
    }
    return ht > 0 ? 1 - hs / ht : 0;
  }
};

inline WrightFisherSummary wright_fisher_replicates(const WrightFisherSpec& w, std::uint64_t seed, std::size_t replicates,
                                                    calx::scheduler::ThreadPool& pool, std::size_t every = 1, std::size_t shard = 64) {
  detail::check(w);
  if (every == 0 || shard == 0) throw std::invalid_argument("evolution::population::wright_fisher_replicates: positive steps");
  const std::size_t records = w.generations / every + 1, shards = (replicates + shard - 1) / shard;
  struct Part {
    std::size_t fixed = 0, lost = 0;
    double settling = 0;
    std::vector<double> mean_p, h_s, h_s_sq, h_t;
  };
  std::vector<Part> parts(shards);
  calx::scheduler::TaskGraph graph;
  std::vector<std::size_t> ids;
  for (std::size_t k = 0; k < shards; k++) {
    ids.push_back(graph.add(
        "replicates " + std::to_string(k * shard),
        [&, k](std::size_t) {
          Part& part = parts[k];
          part.mean_p.assign(records, 0);
          part.h_s.assign(records, 0);
          part.h_s_sq.assign(records, 0);
          part.h_t.assign(records, 0);
          std::vector<double> p(w.demes), next(w.demes);
          const double D = static_cast<double>(w.demes);
          for (std::size_t rep = k * shard; rep < std::min(replicates, (k + 1) * shard); rep++) {
            calx::random::Stream r(seed, WRIGHT_FISHER, static_cast<std::uint32_t>(rep), 0);
            std::fill(p.begin(), p.end(), w.p0);
            const auto record = [&](std::size_t at) {
              double mean = 0, hs = 0;
              for (const double x : p) {
                mean += x;
                hs += 2 * x * (1 - x);
              }
              mean /= D;
              hs /= D;
              part.mean_p[at] += mean;
              part.h_s[at] += hs;
              part.h_s_sq[at] += hs * hs;
              part.h_t[at] += 2 * mean * (1 - mean);
            };
            record(0);
            bool settled = false;
            std::size_t g = 0;
            while (g < w.generations && !settled) {
              settled = detail::generation(w, p, next, r);
              g++;
              if (g % every == 0) record(g / every);
            }
            if (settled) {
              // Nothing changes any more: the remaining records repeat the last state.
              for (std::size_t at = g / every + 1; at < records; at++) record(at);
              part.settling += static_cast<double>(g);
              bool all = true, none = true;
              for (const double x : p) {
                all = all && x == 1;
                none = none && x == 0;
              }
              part.fixed += all;
              part.lost += none;
            }
          }
        },
        {}, static_cast<double>(std::min(replicates, (k + 1) * shard) - k * shard)));
  }
  WrightFisherSummary out;
  graph.add(
      "merge",
      [&](std::size_t) {
        out.replicates = replicates;
        out.every = every;
        out.mean_p.assign(records, 0);
        out.h_s.assign(records, 0);
        out.h_s_sq.assign(records, 0);
        out.h_t.assign(records, 0);
        double settling = 0;
        for (const Part& part : parts) {
          out.fixed += part.fixed;
          out.lost += part.lost;
          settling += part.settling;
          for (std::size_t i = 0; i < records; i++) {
            out.mean_p[i] += part.mean_p[i];
            out.h_s[i] += part.h_s[i];
            out.h_s_sq[i] += part.h_s_sq[i];
            out.h_t[i] += part.h_t[i];
          }
        }
        const double R = static_cast<double>(std::max<std::size_t>(1, replicates));
        for (std::size_t i = 0; i < records; i++) {
          out.mean_p[i] /= R;
          out.h_s[i] /= R;
          out.h_s_sq[i] /= R;
          out.h_t[i] /= R;
        }
        const std::size_t settled = out.fixed + out.lost;
        out.mean_settling = settled ? settling / static_cast<double>(settled) : 0;
      },
      ids);
  out.trace = graph.run(pool);
  return out;
}

// ---- the diffusion -------------------------------------------------------------------------------------

// Kimura's fixation probability for genic selection, s_K the advantage of the allele (the
// spec's s = 2s_K with h = ½), effective size N.
inline double kimura_fixation(double N, double s_k, double p) {
  const double a = 4 * N * s_k;
  if (std::fabs(a) < 1e-12) return p;
  return -std::expm1(-a * p) / -std::expm1(-a);
}

// The diffusion's fixation probability for fitnesses 1 + s, 1 + hs, 1, by quadrature.
inline double fixation_probability(double N, double s, double h, double p) {
  if (p <= 0) return 0;
  if (p >= 1) return 1;
  const auto f = [&](double x) { return -2 * N * s * (2 * h * x + (1 - 2 * h) * x * x); };
  // The largest exponent on [0, 1], subtracted so that G never overflows.
  double top = std::max(f(0), f(1));
  if (1 - 2 * h != 0) {
    const double vertex = -h / (1 - 2 * h);
    if (vertex > 0 && vertex < 1) top = std::max(top, f(vertex));
  }
  const auto G = [&](double x) { return std::exp(f(x) - top); };
  const double below = calx::quadrature::integrate(G, 0, p, 0, 1e-12).value;
  const double above = calx::quadrature::integrate(G, p, 1, 0, 1e-12).value;
  return below / (below + above);
}

// E[2p_t q_t] / (2p₀q₀) = (1 − 1/(2N))^t, without selection or mutation.
inline double heterozygosity_decay(double N, double t) { return std::pow(1 - 1 / (2 * N), t); }

// The equilibrium F_ST of the island model: 1/(1 + 4Nm), and the exact fixed point of
// F′ = (1 − m)² (1/(2N) + (1 − 1/(2N)) F).
inline double island_fst(double N, double m) { return 1 / (1 + 4 * N * m); }
inline double island_fst_exact(double N, double m) {
  const double k = (1 - m) * (1 - m);
  return k / (2 * N - (2 * N - 1) * k);
}

}  // namespace evolution::population
