// The copying model of Li and Stephens: a haplotype as an imperfect mosaic of the haplotypes of
// a reference panel, a hidden Markov model whose states are the panel's haplotypes; its
// likelihood, the posterior probabilities of copying each haplotype, and the most probable
// copying path, each in O(N K) time for N sites and K haplotypes.
//
// The model (Li and Stephens 2003, Genetics 165, 2213–2233: the conditional distribution of
// one haplotype given K others). A panel of K haplotypes h_1 … h_K at N biallelic sites with genetic positions
// g_1 ≤ … ≤ g_N; the target x_1 … x_N copies, at each site j, one haplotype Z_j of the panel:
//
// - Z_1 is uniform on {1, …, K};
// - between sites j − 1 and j, at distance d_j = g_j − g_{j−1} on the map, the copying process
//   jumps with probability r_j = 1 − exp(−ρ d_j / K) to a haplotype drawn uniformly (possibly
//   the one it copied) and otherwise stays:
//     P(Z_j = l | Z_{j−1} = k) = (1 − r_j)[k = l] + r_j / K,
//   the chain observed at the sites of a Poisson process of jumps of rate ρ/K per unit of the
//   map. With ρ = 4N_e c and the map in Morgans, ρ is the population-scaled recombination rate;
//   on the coalescent's chromosome [0, 1), ρ is the coalescent's own (population/coalescent.hpp);
// - x_j is h_{Z_j}(j) with probability 1 − ε and the other allele with probability ε. In Li and
//   Stephens's form a mutation, with probability θ/(K + θ), makes the allele uniform on the two,
//   so ε = θ/(2(K + θ)); the value of θ they suggest, θ̃ = (Σ_{m=1}^{K−1} 1/m)⁻¹, is the
//   default here;
// - a missing allele (−1) of the target is emitted with probability 1 by every state, which
//   marginalises over it: imputation (imputation.hpp) reads the posterior at such sites.
//
// The panel is complete (alleles 0 and 1); it is a population::Haplotypes, sites × sequences, so
// that the K alleles of a site are contiguous and every step below is one pass over them. The
// map is the panel's positions unless another is given.
//
// The algorithms, each the HMM recursion of Rabiner (1989, Proc. IEEE 77, 257–286, §III and
// §V.A) with the sum over predecessors made constant-time by the uniform jump:
//
// - forward: α_j(l) = e_j(l) Σ_k α_{j−1}(k) P(k → l) = e_j(l) ((1 − r_j) α_{j−1}(l) + (r_j/K) Σ_k α_{j−1}(k)),
//   and Σ_k α̂_{j−1}(k) = 1 for the scaled α̂, so each site costs O(K); log P(x) = Σ_j log c_j
//   for the scales c_j;
// - backward: β_{j−1}(k) = (1 − r_j) e_j(k) β_j(k) + (r_j/K) Σ_l e_j(l) β_j(l), scaled by its
//   own constants, so that it computes log P(x) independently of the forward pass;
// - posterior: γ_j(k) = P(Z_j = k | x) ∝ α̂_j(k) β̂_j(k);
// - Viterbi: δ_j(l) = log e_j(l) + max(log(1 − r_j + r_j/K) + δ_{j−1}(l), log(r_j/K) + max_k δ_{j−1}(k)).
//   The second term includes k = l with the smaller weight r_j/K ≤ 1 − r_j + r_j/K, which never
//   wins over the first, so the maximum is exact; ties stay, and jump to the smallest k;
// - forward_dense: the same likelihood by the general O(K²) recursion over the explicit
//   transition matrix of each site, the reference the checks hold the fast forward to.
//
// *Remark (genomics::hmm).* The genomics domain's Model is a homogeneous chain whose
// emissions depend on the state alone, run in O(T S²). The copying chain is inhomogeneous: its
// jumps change with the map, its emissions with the panel's alleles at each site; and its
// uniform jumps make a step O(K), not O(K²). So the recursions are written here, reusing the
// genomics domain's Path and log helpers. On an evenly spaced map the chain is homogeneous, and
// with the site folded into the symbol (symbol 2j + x_j, emitted by state k with the copying
// probability of site j) genomics::hmm computes the same likelihood, posterior and Viterbi path:
// the check `copying` holds the two to each other.
//
// Complexity: O(N K) time for each of forward, backward, posterior, Viterbi and log_likelihood;
// O(N K) memory for the stored passes and the posterior, O(K) for log_likelihood, N K bytes
// for Viterbi's back-pointers. forward_dense costs O(N K²) and is the reference, not the path.
//
// Accuracy: every step of the scaled passes is normalised, so no range is lost on long
// chromosomes; the checks find the fast and dense forwards, the brute-force sum over all K^N
// paths and genomics::hmm agreeing within 1e-12 relative. A target that no path can emit (ε = 0
// and an allele no haplotype carries) has log P = −∞, zero posteriors and an empty path.
//
// Limits: the model's ρ and θ are parameters, not estimated here. The conditional
// distribution approximates the coalescent's; Li and Stephens's likelihood of a whole sample,
// the product of such conditionals as haplotypes are added one by one, depends on the order of
// adding. Only the conditional is written here.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <genomics/hmm/model.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace evolution::haplotypes {

using genomics::hmm::NEG_INF;
using genomics::hmm::Path;
using genomics::hmm::safe_log;
using population::Haplotypes;

// The purpose of this module's streams (the first identity word of calx::random::Stream).
inline constexpr std::uint32_t MOSAIC = 0x4d4f5341;  // "MOSA"

// Li and Stephens's θ̃ = (Σ_{m=1}^{K−1} 1/m)⁻¹ for a panel of K haplotypes (infinite for K < 2).
inline double li_stephens_theta(std::size_t K) {
  double a = 0;
  for (std::size_t m = K; m-- > 1;) a += 1.0 / static_cast<double>(m);
  return a > 0 ? 1 / a : std::numeric_limits<double>::infinity();
}

// ε = θ/(2(K + θ)), the probability of miscopying an allele (1/2 in the limit θ → ∞).
inline double miscopying(double theta, std::size_t K) {
  if (!std::isfinite(theta)) return 0.5;
  return theta / (2 * (static_cast<double>(K) + theta));
}

struct CopyingModel {
  double rho = 1;       // ρ: the copying process jumps at rate ρ/K per unit of the map
  double theta = -1;    // θ: ε = θ/(2(K + θ)); negative for Li and Stephens's θ̃(K)
  double epsilon = -1;  // ε itself when it is ≥ 0, overriding θ
};

// The copying chain on a panel: its K states, N sites, the jump probability r_j before each
// site (r_0 = 1: the first state is uniform) and ε.
struct Chain {
  std::size_t states = 0, sites = 0;
  std::vector<double> jump;
  double epsilon = 0;

  // P(Z_j = k | Z_{j−1} = k) and P(Z_j = l | Z_{j−1} = k) for l ≠ k.
  double stay(std::size_t j) const { return 1 - jump[j] + jump[j] / static_cast<double>(states); }
  double move(std::size_t j) const { return jump[j] / static_cast<double>(states); }
  double transition(std::size_t j, std::size_t k, std::size_t l) const { return k == l ? stay(j) : move(j); }
};

namespace detail {

[[noreturn]] inline void fail(const char* where, const std::string& why) {
  throw std::invalid_argument(std::string("evolution::haplotypes::") + where + ": " + why);
}

inline void check_panel(const Haplotypes& panel, const char* where) {
  if (panel.sequences == 0 || panel.sites == 0) fail(where, "a panel of at least one haplotype and one site");
  if (panel.alleles.size() != panel.sites * panel.sequences) fail(where, "the panel's alleles are sites × sequences");
  for (const std::int8_t a : panel.alleles)
    if (a != 0 && a != 1) fail(where, "the panel is complete: alleles 0 and 1");
}

inline void check_target(const Chain& c, std::span<const std::int8_t> x, const char* where) {
  if (x.size() != c.sites) fail(where, "one allele of the target for each site of the panel");
  for (const std::int8_t a : x)
    if (a < -1 || a > 1) fail(where, "an allele of the target is 0, 1 or −1 (missing)");
}

// e_j(k) for every k: 1 − ε where h_k(j) = x, ε elsewhere; 1 when x is missing.
inline void emissions(const Haplotypes& panel, std::size_t j, std::int8_t x, double eps, double* out) {
  const std::size_t K = panel.sequences;
  if (x < 0) {
    std::fill(out, out + K, 1.0);
    return;
  }
  const std::int8_t* h = panel.alleles.data() + j * K;
  const double same = 1 - eps;
  for (std::size_t k = 0; k < K; k++) out[k] = h[k] == x ? same : eps;
}

}  // namespace detail

// The chain of the model on a panel, the map being the panel's positions unless one is given.
inline Chain copying_chain(const Haplotypes& panel, const CopyingModel& m = {}, std::span<const double> map = {}) {
  detail::check_panel(panel, "copying_chain");
  if (map.empty()) map = panel.positions;
  if (map.size() != panel.sites) detail::fail("copying_chain", "one genetic position for each site");
  if (!(m.rho >= 0) || !std::isfinite(m.rho)) detail::fail("copying_chain", "ρ >= 0");
  Chain c;
  c.states = panel.sequences;
  c.sites = panel.sites;
  c.jump.resize(c.sites);
  c.jump[0] = 1;
  const double K = static_cast<double>(c.states);
  for (std::size_t j = 1; j < c.sites; j++) {
    const double d = map[j] - map[j - 1];
    if (!(d >= 0)) detail::fail("copying_chain", "the genetic map is non-decreasing");
    c.jump[j] = -std::expm1(-m.rho * d / K);
  }
  c.epsilon = m.epsilon >= 0 ? m.epsilon : miscopying(m.theta >= 0 ? m.theta : li_stephens_theta(c.states), c.states);
  if (!(c.epsilon >= 0 && c.epsilon <= 0.5)) detail::fail("copying_chain", "ε in [0, ½]: miscopying no likelier than copying");
  return c;
}

// ---- forward and backward ---------------------------------------------------------------------------

struct Forward {
  std::size_t states = 0, sites = 0;
  std::vector<double> alpha;  // α̂_j(k), sites × states, each row summing to one
  std::vector<double> scale;  // c_j
  double log_likelihood = NEG_INF;
  double at(std::size_t j, std::size_t k) const { return alpha[j * states + k]; }
};
struct Backward {
  std::size_t states = 0, sites = 0;
  std::vector<double> beta;   // β̂_j(k), sites × states, each row summing to one
  std::vector<double> scale;  // d_j
  double log_likelihood = NEG_INF;
  double at(std::size_t j, std::size_t k) const { return beta[j * states + k]; }
};

namespace detail {

// One forward pass, calling row(j, α̂_j) on each normalised row; returns log P(x).
template <class Row>
double forward_pass(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x, Row&& row) {
  const std::size_t K = c.states, N = c.sites;
  std::vector<double> a(K), e(K);
  double logp = 0;
  for (std::size_t j = 0; j < N; j++) {
    emissions(panel, j, x[j], c.epsilon, e.data());
    double total = 0;
    if (j == 0) {
      const double u = 1 / static_cast<double>(K);
      for (std::size_t k = 0; k < K; k++) total += a[k] = u * e[k];
    } else {
      const double stay = 1 - c.jump[j], move = c.move(j);
      for (std::size_t k = 0; k < K; k++) total += a[k] = e[k] * (stay * a[k] + move);
    }
    if (!(total > 0)) {
      row(j, static_cast<const std::vector<double>&>(a), 0.0);
      return NEG_INF;
    }
    const double inv = 1 / total;
    for (double& v : a) v *= inv;
    logp += std::log(total);
    row(j, static_cast<const std::vector<double>&>(a), total);
  }
  return logp;
}

}  // namespace detail

inline Forward forward(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x) {
  detail::check_target(c, x, "forward");
  Forward f;
  f.states = c.states;
  f.sites = c.sites;
  f.alpha.assign(c.sites * c.states, 0);
  f.scale.assign(c.sites, 0);
  f.log_likelihood = detail::forward_pass(panel, c, x, [&](std::size_t j, const std::vector<double>& a, double total) {
    f.scale[j] = total;
    if (total > 0) std::copy(a.begin(), a.end(), f.alpha.begin() + static_cast<std::ptrdiff_t>(j * c.states));
  });
  return f;
}

// log P(x) by the forward recursion in O(K) memory.
inline double log_likelihood(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x) {
  detail::check_target(c, x, "log_likelihood");
  return detail::forward_pass(panel, c, x, [](std::size_t, const std::vector<double>&, double) {});
}

inline Backward backward(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x) {
  detail::check_target(c, x, "backward");
  const std::size_t K = c.states, N = c.sites;
  Backward b;
  b.states = K;
  b.sites = N;
  b.beta.assign(N * K, 0);
  b.scale.assign(N, 0);
  // β_N(k) = 1, normalised.
  const double uniform = 1 / static_cast<double>(K);
  std::fill(b.beta.begin() + static_cast<std::ptrdiff_t>((N - 1) * K), b.beta.end(), uniform);
  b.scale[N - 1] = static_cast<double>(K);
  double logp = std::log(static_cast<double>(K));
  std::vector<double> e(K), w(K);
  for (std::size_t j = N - 1; j > 0; j--) {
    detail::emissions(panel, j, x[j], c.epsilon, e.data());
    const double* next = &b.beta[j * K];
    double W = 0;
    for (std::size_t k = 0; k < K; k++) W += w[k] = e[k] * next[k];
    const double stay = 1 - c.jump[j], jumped = c.move(j) * W;
    double* row = &b.beta[(j - 1) * K];
    double d = 0;
    for (std::size_t k = 0; k < K; k++) d += row[k] = stay * w[k] + jumped;
    b.scale[j - 1] = d;
    if (!(d > 0)) return b;
    const double inv = 1 / d;
    for (std::size_t k = 0; k < K; k++) row[k] *= inv;
    logp += std::log(d);
  }
  detail::emissions(panel, 0, x[0], c.epsilon, e.data());
  double z = 0;
  for (std::size_t k = 0; k < K; k++) z += uniform * e[k] * b.beta[k];
  b.log_likelihood = logp + safe_log(z);
  return b;
}

// ---- posteriors ------------------------------------------------------------------------------------

struct Posterior {
  std::size_t states = 0, sites = 0;
  std::vector<double> gamma;  // γ_j(k) = P(Z_j = k | x), sites × states
  double log_likelihood = NEG_INF;
  double at(std::size_t j, std::size_t k) const { return gamma[j * states + k]; }
  std::span<const double> site(std::size_t j) const { return {gamma.data() + j * states, states}; }
};

inline Posterior posterior(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x) {
  const Forward f = forward(panel, c, x);
  Posterior p;
  p.states = c.states;
  p.sites = c.sites;
  p.log_likelihood = f.log_likelihood;
  p.gamma.assign(c.sites * c.states, 0);
  if (f.log_likelihood == NEG_INF) return p;
  const Backward b = backward(panel, c, x);
  const std::size_t K = c.states;
  for (std::size_t j = 0; j < c.sites; j++) {
    double s = 0;
    for (std::size_t k = 0; k < K; k++) s += p.gamma[j * K + k] = f.alpha[j * K + k] * b.beta[j * K + k];
    for (std::size_t k = 0; k < K; k++) p.gamma[j * K + k] /= s;
  }
  return p;
}

// ---- Viterbi -------------------------------------------------------------------------------------------

// The most probable copying path and log P(x, Z*).
inline Path viterbi(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x) {
  detail::check_target(c, x, "viterbi");
  const std::size_t K = c.states, N = c.sites;
  Path out;
  std::vector<std::uint8_t> stayed(N * K, 1);  // the back-pointer at (j, l): l itself, or best[j − 1]
  std::vector<std::uint32_t> best(N, 0);       // argmax_k δ_j(k), the smallest
  std::vector<double> v(K), e(K);
  const double le_same = safe_log(1 - c.epsilon), le_other = safe_log(c.epsilon);
  const auto log_emit = [&](std::size_t j, std::size_t k) {
    if (x[j] < 0) return 0.0;
    return panel.alleles[j * K + k] == x[j] ? le_same : le_other;
  };
  const double l0 = -std::log(static_cast<double>(K));
  double top = NEG_INF;
  for (std::size_t k = 0; k < K; k++) {
    v[k] = l0 + log_emit(0, k);
    if (v[k] > top) {
      top = v[k];
      best[0] = static_cast<std::uint32_t>(k);
    }
  }
  for (std::size_t j = 1; j < N; j++) {
    const double ls = safe_log(c.stay(j)), lm = safe_log(c.move(j)), jumped = lm + top;
    std::uint8_t* bp = &stayed[j * K];
    double next_top = NEG_INF;
    for (std::size_t k = 0; k < K; k++) {
      const double kept = ls + v[k];
      double s;
      if (kept >= jumped) s = kept;
      else {
        s = jumped;
        bp[k] = 0;
      }
      v[k] = s + log_emit(j, k);
      if (v[k] > next_top) {
        next_top = v[k];
        best[j] = static_cast<std::uint32_t>(k);
      }
    }
    top = next_top;
  }
  out.log_probability = top;
  if (top == NEG_INF) return out;
  out.states.resize(N);
  out.states[N - 1] = best[N - 1];
  for (std::size_t j = N - 1; j > 0; j--) {
    const std::uint32_t l = out.states[j];
    out.states[j - 1] = stayed[j * K + l] ? l : best[j - 1];
  }
  return out;
}

// log P(x, Z) of a given copying path.
inline double path_log_probability(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x, std::span<const std::uint32_t> path) {
  detail::check_target(c, x, "path_log_probability");
  if (path.size() != c.sites) detail::fail("path_log_probability", "one state for each site");
  const std::size_t K = c.states;
  const auto emit = [&](std::size_t j, std::uint32_t k) {
    if (x[j] < 0) return 0.0;
    return safe_log(panel.alleles[j * K + k] == x[j] ? 1 - c.epsilon : c.epsilon);
  };
  double s = -std::log(static_cast<double>(K)) + emit(0, path[0]);
  for (std::size_t j = 1; j < c.sites; j++) s += safe_log(c.transition(j, path[j - 1], path[j])) + emit(j, path[j]);
  return s;
}

// ---- the dense reference ---------------------------------------------------------------------------------

// log P(x) by the general forward recursion over each site's K × K transition matrix: O(N K²).
inline double forward_dense(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> x) {
  detail::check_target(c, x, "forward_dense");
  const std::size_t K = c.states, N = c.sites;
  std::vector<double> a(K), next(K), e(K), T(K * K);
  detail::emissions(panel, 0, x[0], c.epsilon, e.data());
  double total = 0, logp = 0;
  for (std::size_t k = 0; k < K; k++) total += a[k] = e[k] / static_cast<double>(K);
  if (!(total > 0)) return NEG_INF;
  for (double& v : a) v /= total;
  logp += std::log(total);
  for (std::size_t j = 1; j < N; j++) {
    for (std::size_t k = 0; k < K; k++)
      for (std::size_t l = 0; l < K; l++) T[k * K + l] = c.transition(j, k, l);
    detail::emissions(panel, j, x[j], c.epsilon, e.data());
    total = 0;
    for (std::size_t l = 0; l < K; l++) {
      double s = 0;
      for (std::size_t k = 0; k < K; k++) s += a[k] * T[k * K + l];
      total += next[l] = s * e[l];
    }
    if (!(total > 0)) return NEG_INF;
    for (std::size_t l = 0; l < K; l++) a[l] = next[l] / total;
    logp += std::log(total);
  }
  return logp;
}

// ---- simulation from the model ------------------------------------------------------------------------

struct Mosaic {
  std::vector<std::uint32_t> path;  // Z_1 … Z_N
  std::vector<std::int8_t> alleles;  // x_1 … x_N
};

// A target drawn from the model: the copying path, then each allele miscopied with probability ε.
inline Mosaic sample_mosaic(const Haplotypes& panel, const Chain& c, calx::random::Stream& r) {
  detail::check_panel(panel, "sample_mosaic");
  Mosaic m;
  m.path.resize(c.sites);
  m.alleles.resize(c.sites);
  std::uint32_t z = 0;
  for (std::size_t j = 0; j < c.sites; j++) {
    if (j == 0 || r.bernoulli(c.jump[j])) z = static_cast<std::uint32_t>(r.below(c.states));
    m.path[j] = z;
    const std::int8_t h = panel.alleles[j * c.states + z];
    m.alleles[j] = r.bernoulli(c.epsilon) ? static_cast<std::int8_t>(1 - h) : h;
  }
  return m;
}

}  // namespace evolution::haplotypes
