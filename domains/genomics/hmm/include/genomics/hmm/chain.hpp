// Hidden Markov chains whose emissions are given step by step and whose transitions may change from
// step to step: the HMMs of observations that are not symbols of a finite alphabet (genotype
// likelihoods, allele counts, depths) taken at uneven distances (sites along a genome).
//
// A chain has S states and T steps: log initial probabilities log π_k, log emissions
// log e_t(k) = log P(x_t | π_t = k) given for every step (the observation enters only through
// them), and log transitions log a_t(k, l) = log P(π_{t+1} = l | π_t = k) for each step t → t + 1,
// or one matrix for all of them (`stationary`). The probability of the observations is
// P(x) = Σ_π π_{π₁} Π_t e_t(π_t) Π_t a_t(π_t, π_{t+1}), the same sum over the paths of the trellis
// as model.hpp's, with the weights read from the chain instead of from a table of symbols
// (Rabiner 1989, Proc. IEEE 77, 257–286, §III; Durbin et al. 1998, §3.2). As there:
//
// - viterbi: the tropical recurrence δ_{t+1}(l) = log e_{t+1}(l) + max_k [δ_t(k) + log a_t(k, l)],
//   with back-pointers, ties to the smallest k (Rabiner eqs. 32–35);
// - forward and backward in the log semiring (log-sum-exp), which keeps the range whatever the
//   emissions (they need not sum to one over observations: likelihoods up to a factor per step,
//   as genotype likelihoods are, leave the posterior unchanged and shift log P(x) by Σ_t log c_t);
// - posterior: γ_t(k) = exp(α_t(k) + β_t(k) − log P(x)), exact (Durbin et al. §3.2), and its
//   decoding;
// - path_log_probability of a given path, for checks.
//
// bcftools roh (Narasimhan et al. 2016, Bioinformatics 32, 1749–1751) is such a chain over the
// variant sites of a sample, two states, emissions from genotype likelihoods and transitions
// from the distance between sites; genomics::structural builds it (allelic.hpp).
//
// Complexity O(T S²) for each algorithm, memory O(T S) (one byte per back-pointer when S ≤ 256).
// Accuracy: forward and backward agree with the enumeration of every path to rounding (checks).
//
// *Remark.* Like model.hpp this is statistics, and moves with it to the inference foundation.

#pragma once

#include <genomics/hmm/model.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::hmm {

struct Chain {
  std::size_t states = 0, steps = 0;
  std::vector<double> log_initial;     // S
  std::vector<double> log_emission;    // T × S
  std::vector<double> log_transition;  // (T − 1) × S × S (row k, column l), or S × S when stationary
  bool stationary = false;

  Chain() = default;
  Chain(std::size_t S, std::size_t T, bool one_transition = false)
      : states(S), steps(T), log_initial(S, NEG_INF), log_emission(T * S, 0), log_transition((one_transition ? 1 : (T > 0 ? T - 1 : 0)) * S * S, NEG_INF), stationary(one_transition) {
    if (S == 0) throw std::invalid_argument("genomics::hmm::Chain: at least one state");
  }
  double e(std::size_t t, std::size_t k) const { return log_emission[t * states + k]; }
  double& e(std::size_t t, std::size_t k) { return log_emission[t * states + k]; }
  // log a_t(k, l): the transition from step t to step t + 1.
  double a(std::size_t t, std::size_t k, std::size_t l) const { return log_transition[(stationary ? 0 : t) * states * states + k * states + l]; }
  double& a(std::size_t t, std::size_t k, std::size_t l) { return log_transition[(stationary ? 0 : t) * states * states + k * states + l]; }

  void validate() const {
    const std::size_t S = states, T = steps;
    if (log_initial.size() != S || log_emission.size() != T * S || log_transition.size() != (stationary ? 1 : (T > 0 ? T - 1 : 0)) * S * S)
      throw std::invalid_argument("genomics::hmm::Chain: the sizes of the parameters do not agree with the states and steps");
  }
};

namespace detail {

template <class Pointer>
Path chain_viterbi(const Chain& c) {
  const std::size_t S = c.states, T = c.steps;
  Path out;
  std::vector<Pointer> back(T * S);
  std::vector<double> v(S), w(S);
  for (std::size_t k = 0; k < S; k++) v[k] = c.log_initial[k] + c.e(0, k);
  for (std::size_t t = 1; t < T; t++) {
    for (std::size_t l = 0; l < S; l++) {
      double best = NEG_INF;
      std::size_t arg = 0;
      for (std::size_t k = 0; k < S; k++) {
        const double x = v[k] + c.a(t - 1, k, l);
        if (x > best) { best = x; arg = k; }
      }
      w[l] = best + c.e(t, l);
      back[t * S + l] = static_cast<Pointer>(arg);
    }
    std::swap(v, w);
  }
  const std::size_t arg = static_cast<std::size_t>(std::max_element(v.begin(), v.end()) - v.begin());
  out.log_probability = v[arg];
  if (out.log_probability == NEG_INF) return out;
  out.states.resize(T);
  out.states[T - 1] = static_cast<std::uint32_t>(arg);
  for (std::size_t t = T - 1; t > 0; t--) out.states[t - 1] = static_cast<std::uint32_t>(back[t * S + out.states[t]]);
  return out;
}

}  // namespace detail

// The most probable path and its log probability (log P(x, π*)).
inline Path viterbi(const Chain& c) {
  c.validate();
  if (c.steps == 0) { Path p; p.log_probability = 0; return p; }
  if (c.states <= 256) return detail::chain_viterbi<std::uint8_t>(c);
  return detail::chain_viterbi<std::uint32_t>(c);
}

// log α_t(k) = log P(x_1..t, π_t = k), T × S.
inline std::vector<double> log_forward(const Chain& c, double* log_likelihood = nullptr) {
  c.validate();
  const std::size_t S = c.states, T = c.steps;
  std::vector<double> alpha(T * S, NEG_INF);
  if (T == 0) { if (log_likelihood) *log_likelihood = 0; return alpha; }
  for (std::size_t k = 0; k < S; k++) alpha[k] = c.log_initial[k] + c.e(0, k);
  for (std::size_t t = 1; t < T; t++)
    for (std::size_t l = 0; l < S; l++) {
      double s = NEG_INF;
      for (std::size_t k = 0; k < S; k++) s = log_add(s, alpha[(t - 1) * S + k] + c.a(t - 1, k, l));
      alpha[t * S + l] = s + c.e(t, l);
    }
  if (log_likelihood) {
    double s = NEG_INF;
    for (std::size_t k = 0; k < S; k++) s = log_add(s, alpha[(T - 1) * S + k]);
    *log_likelihood = s;
  }
  return alpha;
}

// log β_t(k) = log P(x_{t+1..T} | π_t = k), T × S.
inline std::vector<double> log_backward(const Chain& c, double* log_likelihood = nullptr) {
  c.validate();
  const std::size_t S = c.states, T = c.steps;
  std::vector<double> beta(T * S, NEG_INF);
  if (T == 0) { if (log_likelihood) *log_likelihood = 0; return beta; }
  for (std::size_t k = 0; k < S; k++) beta[(T - 1) * S + k] = 0;
  for (std::size_t t = T - 1; t > 0; t--)
    for (std::size_t k = 0; k < S; k++) {
      double s = NEG_INF;
      for (std::size_t l = 0; l < S; l++) s = log_add(s, c.a(t - 1, k, l) + c.e(t, l) + beta[t * S + l]);
      beta[(t - 1) * S + k] = s;
    }
  if (log_likelihood) {
    double s = NEG_INF;
    for (std::size_t k = 0; k < S; k++) s = log_add(s, c.log_initial[k] + c.e(0, k) + beta[k]);
    *log_likelihood = s;
  }
  return beta;
}

inline double log_likelihood(const Chain& c) {
  double l = NEG_INF;
  (void)log_forward(c, &l);
  return l;
}

// The posterior of the states at every step.
inline Posterior posterior(const Chain& c) {
  Posterior p;
  p.states = c.states;
  const std::vector<double> alpha = log_forward(c, &p.log_likelihood);
  p.gamma.assign(c.steps * c.states, 0);
  if (p.log_likelihood == NEG_INF) return p;
  const std::vector<double> beta = log_backward(c);
  for (std::size_t i = 0; i < p.gamma.size(); i++) p.gamma[i] = std::exp(alpha[i] + beta[i] - p.log_likelihood);
  return p;
}

inline double path_log_probability(const Chain& c, const std::vector<std::uint32_t>& path) {
  if (path.size() != c.steps) throw std::invalid_argument("genomics::hmm::path_log_probability: one state for each step");
  if (path.empty()) return 0;
  double s = c.log_initial[path[0]] + c.e(0, path[0]);
  for (std::size_t t = 1; t < path.size(); t++) s += c.a(t - 1, path[t - 1], path[t]) + c.e(t, path[t]);
  return s;
}

// The chain of a symbol model over a sequence: the same probabilities, read step by step.
inline Chain chain_of(const Model& m, const Codes& x) {
  if (m.ends()) throw std::invalid_argument("genomics::hmm::chain_of: a model without an end state");
  Chain c(m.states, x.size(), true);
  for (std::size_t k = 0; k < m.states; k++) c.log_initial[k] = safe_log(m.initial[k]);
  for (std::size_t t = 0; t < x.size(); t++)
    for (std::size_t k = 0; k < m.states; k++) c.e(t, k) = safe_log(m.emit(k, x[t]));
  for (std::size_t k = 0; k < m.states; k++)
    for (std::size_t l = 0; l < m.states; l++) c.a(0, k, l) = safe_log(m.a(k, l));
  return c;
}

}  // namespace genomics::hmm
