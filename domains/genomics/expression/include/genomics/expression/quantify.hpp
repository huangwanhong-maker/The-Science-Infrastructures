// Transcript abundances from equivalence-class counts by the EM algorithm (Bray, Pimentel,
// Melsted and Pachter 2016, Nature Biotechnology 34, 525–527, Online Methods; Patro, Duggal,
// Love, Irizarry and Kingsford 2017, Nature Methods 14, 417–419, for the EM over classes), with
// kallisto's convergence rule, TPM, and the bootstrap of the counts.
//
// - The likelihood (Bray et al., Online Methods; Pachter 2011, arXiv:1104.3889): with α_t the
//   probability that a fragment comes from transcript t and l̃_t its effective length
//   (transcriptome.hpp), a fragment of class c (n_c of them) has probability Σ_{t∈c} α_t / l̃_t,
//   so that
//     L(α) = Π_c (Σ_{t∈c} α_t / l̃_t)^{n_c},   α in the simplex.
//   log L is concave in α, a sum of logarithms of linear forms.
// - The EM (Dempster, Laird and Rubin 1977, J. R. Statist. Soc. B 39, 1–38): the latent
//   variable is the transcript of each fragment; given α its posterior within class c is
//   (α_t / l̃_t) / Σ_{u∈c} α_u / l̃_u, and the new α is the expected fraction of the fragments,
//     α'_t = (1/N) Σ_{c∋t} n_c (α_t / l̃_t) / Σ_{u∈c} α_u / l̃_u.
//   It is a fixed point iteration on the simplex whose fixed points are the stationary points
//   of L, and it never decreases L (Dempster et al., Theorem 1; theory note, Proposition 5.2).
//   Held here as estimated counts N α, as kallisto holds them, starting from N/T each.
// - Convergence (kallisto's EMAlgorithm.h, release 0.46): at least 50 rounds and at most
//   10 000; it stops when every transcript of estimated count above 0.01 changes by less than
//   one per cent in a round. At the end counts below 10⁻⁸ are set to 0.
// - TPM: ρ_t ∝ α_t / l̃_t, the fraction of the molecules, times 10⁶ (Li and Dewey 2011, BMC
//   Bioinformatics 12, 323).
// - The bootstrap (Bray et al., Online Methods): the class counts are resampled from the
//   multinomial law of N fragments with probabilities n_c / N, by conditional binomials
//   (calx::random), and the EM is run again; each replicate from its own stream, in parallel.
//
// A round costs O(Σ_c |c|) for the classes with fragments: the denominators of the classes in
// parallel, then each transcript's sum over its classes in parallel, each slot written by one
// task, so the estimates are the same on any number of threads. Below 2^18 incidences a round
// runs on the calling thread, where waking the pool would cost more than the round.
//
// *Remark.* The EM, the bootstrap and the likelihood are general statistics; they move to the
// inference foundation when it exists.

#pragma once

#include <genomics/expression/pseudoalign.hpp>
#include <genomics/expression/transcriptome.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::expression {

struct EmOptions {
  int min_rounds = 50, max_rounds = 10000;
  double change_limit = 1e-2;  // transcripts of estimated count above this ...
  double change = 1e-2;        // ... must change by less than this fraction for the EM to stop
  double zero_below = 1e-8;    // estimated counts below this are set to 0 at the end
  bool trace = false;          // keep log L after every round
};

struct Abundances {
  std::vector<double> counts;            // estimated fragments per transcript, N α
  std::vector<double> tpm;
  std::vector<double> effective_length;
  std::size_t fragments = 0;             // N, the fragments with a class
  int rounds = 0;
  bool converged = false;
  std::vector<double> log_likelihood;    // log L after each round, when traced
};

// TPM from estimated counts and effective lengths.
inline std::vector<double> tpm(const std::vector<double>& counts, const std::vector<double>& effective) {
  std::vector<double> out(counts.size(), 0.0);
  double total = 0;
  for (std::size_t t = 0; t < counts.size(); t++) total += out[t] = effective[t] > 0 ? counts[t] / effective[t] : 0.0;
  if (total > 0)
    for (double& x : out) x *= 1e6 / total;
  return out;
}

// log L(α) at estimated counts N α (normalised here to the simplex).
inline double log_likelihood(const ClassTable& classes, const std::vector<std::uint32_t>& n, const std::vector<double>& effective, const std::vector<double>& counts) {
  double N = 0;
  for (double x : counts) N += x;
  if (!(N > 0)) return -INFINITY;
  double ll = 0;
  for (std::uint32_t c = 0; c < n.size(); c++) {
    if (n[c] == 0) continue;
    double s = 0;
    for (std::uint32_t t : classes.members(c)) s += counts[t] / (N * effective[t]);
    ll += static_cast<double>(n[c]) * std::log(s);
  }
  return ll;
}

// Each class's fragments divided equally among its transcripts: the baseline the EM improves on.
inline std::vector<double> equal_split(const ClassTable& classes, const std::vector<std::uint32_t>& n, std::size_t transcripts) {
  std::vector<double> out(transcripts, 0.0);
  for (std::uint32_t c = 0; c < n.size(); c++) {
    if (n[c] == 0) continue;
    const auto m = classes.members(c);
    for (std::uint32_t t : m) out[t] += static_cast<double>(n[c]) / static_cast<double>(m.size());
  }
  return out;
}

inline Abundances quantify(const ClassTable& classes, const std::vector<std::uint32_t>& n, const std::vector<double>& effective, const EmOptions& o,
                           calx::scheduler::ThreadPool& pool) {
  if (n.size() > classes.size()) throw std::invalid_argument("genomics::expression::quantify: more counts than classes");
  const std::size_t T = effective.size();
  for (double l : effective)
    if (!(l > 0)) throw std::invalid_argument("genomics::expression::quantify: effective lengths must be positive");
  Abundances a;
  a.effective_length = effective;
  // The classes with fragments, and for each transcript the classes it is in.
  std::vector<std::uint32_t> used;
  std::vector<std::uint32_t> degree(T + 1, 0);
  for (std::uint32_t c = 0; c < n.size(); c++) {
    if (n[c] == 0) continue;
    used.push_back(c);
    a.fragments += n[c];
    for (std::uint32_t t : classes.members(c)) {
      if (t >= T) throw std::invalid_argument("genomics::expression::quantify: a class names a transcript without an effective length");
      degree[t + 1]++;
    }
  }
  for (std::size_t t = 0; t < T; t++) degree[t + 1] += degree[t];
  std::vector<std::uint32_t> incident(degree[T]), fill(degree.begin(), degree.end() - 1);
  for (std::uint32_t k = 0; k < used.size(); k++)
    for (std::uint32_t t : classes.members(used[k])) incident[fill[t]++] = k;
  const double N = static_cast<double>(a.fragments);
  std::vector<double> alpha(T, T ? N / static_cast<double>(T) : 0.0), next(T, 0.0), weight(used.size(), 0.0);
  if (a.fragments == 0) {
    a.counts.assign(T, 0.0);
    a.tpm.assign(T, 0.0);
    a.converged = true;
    return a;
  }
  // A round is two parallel loops; below about 2^18 incidences the threads' wake-ups cost more
  // than the round, so it runs on the calling thread (the same arithmetic, slot for slot).
  const std::size_t grain = 256;
  calx::scheduler::ThreadPool serial(1);
  calx::scheduler::ThreadPool& runner = degree[T] >= (1u << 18) ? pool : serial;
  for (int round = 1; round <= o.max_rounds; round++) {
    // n_c / Σ_{u∈c} α_u / l̃_u for each class with fragments.
    calx::scheduler::parallel_for(runner, used.size(), grain, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t k = b; k < e; k++) {
        double d = 0;
        for (std::uint32_t t : classes.members(used[k])) d += alpha[t] / effective[t];
        weight[k] = d > 0 ? static_cast<double>(n[used[k]]) / d : 0.0;
      }
    });
    calx::scheduler::parallel_for(runner, T, grain, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t t = b; t < e; t++) {
        double s = 0;
        for (std::uint32_t q = degree[t]; q < degree[t + 1]; q++) s += weight[incident[q]];
        next[t] = alpha[t] / effective[t] * s;
      }
    });
    bool stop = true;
    for (std::size_t t = 0; t < T && stop; t++)
      if (next[t] > o.change_limit && std::fabs(next[t] - alpha[t]) / next[t] > o.change) stop = false;
    std::swap(alpha, next);
    a.rounds = round;
    if (o.trace) a.log_likelihood.push_back(log_likelihood(classes, n, effective, alpha));
    if (stop && round >= o.min_rounds) {
      a.converged = true;
      break;
    }
  }
  for (double& x : alpha)
    if (x < o.zero_below) x = 0;
  a.counts = alpha;
  a.tpm = tpm(alpha, effective);
  return a;
}

inline Abundances quantify(const Pseudoalignment& p, const std::vector<double>& effective, const EmOptions& o, calx::scheduler::ThreadPool& pool) {
  return quantify(p.classes, p.counts, effective, o, pool);
}

// ---- the bootstrap ------------------------------------------------------------------------------------
struct Bootstrap {
  std::size_t replicates = 0;
  std::vector<double> mean, sd;  // of the estimated counts, per transcript
};

// The class counts resampled from Multinomial(N, n_c / N).
inline std::vector<std::uint32_t> resample(const std::vector<std::uint32_t>& n, calx::random::Stream& r) {
  std::int64_t left = 0;
  for (std::uint32_t x : n) left += x;
  double mass = static_cast<double>(left);
  std::vector<std::uint32_t> out(n.size(), 0);
  for (std::size_t c = 0; c < n.size() && left > 0; c++) {
    if (n[c] == 0) continue;
    const double p = std::min(1.0, static_cast<double>(n[c]) / mass);
    const std::int64_t x = r.binomial(left, p);
    out[c] = static_cast<std::uint32_t>(x);
    left -= x;
    mass -= static_cast<double>(n[c]);
  }
  return out;
}

inline Bootstrap bootstrap(const ClassTable& classes, const std::vector<std::uint32_t>& n, const std::vector<double>& effective, std::size_t replicates,
                           std::uint64_t seed, const EmOptions& o, calx::scheduler::ThreadPool& pool) {
  const std::size_t T = effective.size();
  std::vector<std::vector<double>> runs(replicates);
  calx::scheduler::parallel_for(pool, replicates, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    calx::scheduler::ThreadPool one(1);
    for (std::size_t k = b; k < e; k++) {
      calx::random::Stream r(seed, BOOTSTRAP, static_cast<std::uint32_t>(k), 0);
      runs[k] = quantify(classes, resample(n, r), effective, o, one).counts;
    }
  });
  Bootstrap out;
  out.replicates = replicates;
  out.mean.assign(T, 0.0);
  out.sd.assign(T, 0.0);
  if (replicates == 0) return out;
  for (const auto& x : runs)
    for (std::size_t t = 0; t < T; t++) out.mean[t] += x[t] / static_cast<double>(replicates);
  for (const auto& x : runs)
    for (std::size_t t = 0; t < T; t++) out.sd[t] += (x[t] - out.mean[t]) * (x[t] - out.mean[t]);
  for (double& v : out.sd) v = replicates > 1 ? std::sqrt(v / static_cast<double>(replicates - 1)) : 0.0;
  return out;
}

}  // namespace genomics::expression
