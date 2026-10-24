// Hidden Markov models on a finite alphabet, and their three problems: the probability of a
// sequence, its most probable path, and the estimation of the parameters from sequences.
//
// A model has S emitting states, K symbols, initial probabilities π_k (the transitions a_{0k}
// out of a silent begin state), transitions a_kl, emissions e_k(b), and optionally the
// transitions a_{k0} into a silent end state (Durbin, Eddy, Krogh and Mitchison 1998,
// Biological Sequence Analysis, §3.2): with an end, Σ_l a_kl + a_{k0} = 1 and a path's
// probability includes its ending; without one (Rabiner's convention), Σ_l a_kl = 1 and the
// length is given. A code ≥ K in a sequence is missing data (N in DNA): every state emits it
// with probability 1, which marginalises over the symbol.
//
// The probability of x = x_1 … x_T is the sum over the state paths π of the product of their
// weights, P(x) = Σ_π a_{0π_1} Π_t e_{π_t}(x_t) a_{π_t π_{t+1}} (· a_{π_T 0}); the most
// probable path is the same expression with the sum replaced by a maximum. Both are one
// dynamic programme over the trellis in a semiring of log weights (path_sum<R>): LogSum
// (⊕ = log-sum-exp, ⊗ = +) gives log P(x), Tropical (⊕ = max, ⊗ = +) the Viterbi score
// (Durbin et al. §3.2; Pachter and Sturmfels 2005, Algebraic Statistics for Computational
// Biology, §2.2, for HMMs as tropical and ordinary polynomials).
//
// - viterbi: the tropical recurrence δ_t(l) = log e_l(x_t) + max_k (δ_{t−1}(k) + log a_kl) with
//   back-pointers, ties to the smallest k (Viterbi 1967; Rabiner 1989, Proc. IEEE 77, 257–286,
//   eqs. 32–35).
// - forward and backward in probability space with scaling: α̂_t = α_t / Π_{s≤t} c_s normalised
//   at every t (Rabiner 1989, §V.A), log P(x) = Σ_t log c_t (+ log Σ_k α̂_T(k) a_k0);
//   the backward variables are normalised at every t by their own constants, so that
//   backward() computes log P(x) independently of forward().
// - posterior: γ_t(k) = P(π_t = k | x) ∝ α̂_t(k) β̂_t(k), normalised over k at each t, which is
//   exact whatever the scales (Durbin et al. §3.2); posterior decoding takes argmax_k γ_t(k).
// - Baum–Welch: EM for the parameters from many sequences (Baum, Petrie, Soules and Weiss
//   1970; Rabiner 1989 eqs. 40a–c and 109–110 for several sequences; Durbin et al. §3.3), with
//   ξ_t(k, l) ∝ α̂_t(k) a_kl e_l(x_{t+1}) β̂_{t+1}(l) normalised over (k, l). Pseudocounts r are
//   added to the expected counts of the parameters the model allows (a zero stays a zero, so
//   that a topology is kept): this is EM for the maximum a posteriori estimate under Dirichlet
//   priors with parameters r + 1 (Dempster, Laird and Rubin 1977), so the quantity that never
//   decreases is the log posterior log P(x | θ) + Σ r log θ, reported beside the likelihood;
//   with r = 0 it is the likelihood itself. Sequences are taken in chunks of fixed size whose
//   counts are merged in order, so the result is the same on any number of threads.
// - sampling of a path and its symbols from a calx::random::Stream, of a given length or to
//   the end state.
//
// Complexity: O(T S²) time for each of viterbi, forward, backward and an iteration of
// Baum–Welch per sequence of length T, on flat arrays with the transition matrix transposed so
// that the inner sums run over contiguous memory; memory O(T S) for the posterior and for
// Viterbi's back-pointers (one byte each when S ≤ 256), O(S) for log_likelihood. path_sum<LogSum>
// costs an exponential and a logarithm per term and is the reference, not the fast path.
//
// Accuracy: the scaled recursions lose no range (each step is normalised), and agree with
// log-space sums to rounding (the checks: 1e-10 relative on enumerations). A sequence the
// model cannot emit has log P = −∞, an empty Viterbi path and zero posteriors.
//
// *Remark.* The machinery names no object in nature: it is statistics, and moves to the
// planned inference foundation when it exists (CLAUDE.md §8, open decision 3). It is written
// for any number of states and symbols so that the move is mechanical.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::hmm {

using sequence::Code;
using sequence::Codes;

inline constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

// The purpose of the streams of this module (the first identity word of calx::random::Stream).
inline constexpr std::uint32_t SAMPLE = 0x484d4d53;  // "HMMS"

inline double safe_log(double p) { return p > 0 ? std::log(p) : NEG_INF; }

// log(e^a + e^b) without overflow, −∞ the identity.
inline double log_add(double a, double b) {
  if (a < b) std::swap(a, b);
  if (b == NEG_INF) return a;
  return a + std::log1p(std::exp(b - a));
}

// ---- semirings of log weights ----------------------------------------------------------------
// The tropical semiring (max, +) and the log semiring (log-sum-exp, +), on ℝ ∪ {−∞}: the two
// readings of one recurrence. Both are commutative semirings with zero −∞ and one 0, and the
// map log: (ℝ≥0, +, ×) → (ℝ ∪ {−∞}, ⊕, +) is an isomorphism onto the log semiring; the tropical
// semiring is its limit under log_h with base h → 0 (Maslov dequantisation).
struct Tropical {
  static constexpr double zero() { return NEG_INF; }
  static constexpr double one() { return 0; }
  static double plus(double a, double b) { return std::max(a, b); }
  static double times(double a, double b) { return a + b; }
};
struct LogSum {
  static constexpr double zero() { return NEG_INF; }
  static constexpr double one() { return 0; }
  static double plus(double a, double b) { return log_add(a, b); }
  static double times(double a, double b) { return a + b; }
};

// ---- the model -------------------------------------------------------------------------------------
struct Model {
  std::size_t states = 0, symbols = 0;
  std::vector<double> initial;     // π_k = a_{0k}, S
  std::vector<double> transition;  // a_kl, S × S, row-major
  std::vector<double> emission;    // e_k(b), S × K, row-major
  std::vector<double> end;         // a_{k0}, S; empty when the end is not modelled

  Model() = default;
  Model(std::size_t S, std::size_t K, bool with_end = false)
      : states(S), symbols(K), initial(S, 0), transition(S * S, 0), emission(S * K, 0), end(with_end ? S : 0, 0) {
    if (S == 0 || K == 0) throw std::invalid_argument("genomics::hmm::Model: at least one state and one symbol");
  }

  bool ends() const { return !end.empty(); }
  double a(std::size_t k, std::size_t l) const { return transition[k * states + l]; }
  double& a(std::size_t k, std::size_t l) { return transition[k * states + l]; }
  double e(std::size_t k, std::size_t b) const { return emission[k * symbols + b]; }
  double& e(std::size_t k, std::size_t b) { return emission[k * symbols + b]; }
  // The emission of a code, 1 for a missing one (a code ≥ K).
  double emit(std::size_t k, Code b) const { return b < symbols ? emission[k * symbols + b] : 1.0; }

  // Throws std::invalid_argument unless the sizes agree, every probability is in [0, 1] and
  // every distribution sums to one within `tolerance`.
  void validate(double tolerance = 1e-9) const {
    const auto fail = [](const std::string& m) { throw std::invalid_argument("genomics::hmm::Model: " + m); };
    if (initial.size() != states || transition.size() != states * states || emission.size() != states * symbols || (ends() && end.size() != states))
      fail("the sizes of the parameters do not agree with the states and symbols");
    const auto in_range = [](const std::vector<double>& v) {
      for (double p : v)
        if (!(p >= 0 && p <= 1)) return false;
      return true;
    };
    if (!in_range(initial) || !in_range(transition) || !in_range(emission) || !in_range(end)) fail("a probability outside [0, 1]");
    double s = 0;
    for (double p : initial) s += p;
    if (std::fabs(s - 1) > tolerance) fail("the initial probabilities sum to " + std::to_string(s));
    for (std::size_t k = 0; k < states; k++) {
      double r = ends() ? end[k] : 0, q = 0;
      for (std::size_t l = 0; l < states; l++) r += a(k, l);
      for (std::size_t b = 0; b < symbols; b++) q += e(k, b);
      if (std::fabs(r - 1) > tolerance) fail("the transitions out of state " + std::to_string(k) + " sum to " + std::to_string(r));
      if (std::fabs(q - 1) > tolerance) fail("the emissions of state " + std::to_string(k) + " sum to " + std::to_string(q));
    }
  }

  // Scales every distribution to sum to one (a row of zeros is left as it is).
  void normalise() {
    const auto scale = [](double* p, std::size_t n, double extra, double* q) {
      double s = extra;
      for (std::size_t i = 0; i < n; i++) s += p[i];
      if (!(s > 0)) return;
      for (std::size_t i = 0; i < n; i++) p[i] /= s;
      if (q) *q /= s;
    };
    scale(initial.data(), states, 0, nullptr);
    for (std::size_t k = 0; k < states; k++) {
      scale(&transition[k * states], states, ends() ? end[k] : 0, ends() ? &end[k] : nullptr);
      scale(&emission[k * symbols], symbols, 0, nullptr);
    }
  }
};

// ---- the generic recurrence ---------------------------------------------------------------------
// ⊕ over every state path of ⊗ of its log weights, in the semiring R: log P(x) for LogSum, the
// log probability of the best path for Tropical. O(T S²) operations of R.
template <class R>
double path_sum(const Model& m, const Codes& x) {
  const std::size_t S = m.states, T = x.size();
  if (T == 0) return m.ends() ? R::zero() : R::one();
  std::vector<double> v(S), w(S);
  for (std::size_t k = 0; k < S; k++) v[k] = R::times(safe_log(m.initial[k]), safe_log(m.emit(k, x[0])));
  for (std::size_t t = 1; t < T; t++) {
    for (std::size_t l = 0; l < S; l++) {
      double s = R::zero();
      for (std::size_t k = 0; k < S; k++) s = R::plus(s, R::times(v[k], safe_log(m.a(k, l))));
      w[l] = R::times(s, safe_log(m.emit(l, x[t])));
    }
    std::swap(v, w);
  }
  double s = R::zero();
  for (std::size_t k = 0; k < S; k++) s = R::plus(s, m.ends() ? R::times(v[k], safe_log(m.end[k])) : v[k]);
  return s;
}

namespace detail {

inline std::vector<double> transposed(const Model& m) {
  const std::size_t S = m.states;
  std::vector<double> t(S * S);
  for (std::size_t k = 0; k < S; k++)
    for (std::size_t l = 0; l < S; l++) t[l * S + k] = m.a(k, l);
  return t;
}

}  // namespace detail

// ---- Viterbi ------------------------------------------------------------------------------------------
struct Path {
  std::vector<std::uint32_t> states;  // π_1 … π_T; empty when x cannot be emitted
  double log_probability = NEG_INF;   // log P(x, π*), natural logarithm
};

namespace detail {

template <class Pointer>
Path viterbi_with(const Model& m, const Codes& x) {
  const std::size_t S = m.states, T = x.size();
  Path out;
  std::vector<double> logT(S * S), logE(S * (m.symbols + 1));
  {
    const auto t = transposed(m);
    for (std::size_t i = 0; i < S * S; i++) logT[i] = safe_log(t[i]);
    for (std::size_t k = 0; k < S; k++) {
      for (std::size_t b = 0; b < m.symbols; b++) logE[k * (m.symbols + 1) + b] = safe_log(m.e(k, b));
      logE[k * (m.symbols + 1) + m.symbols] = 0;  // missing
    }
  }
  const auto le = [&](std::size_t k, Code b) { return logE[k * (m.symbols + 1) + std::min<std::size_t>(b, m.symbols)]; };
  std::vector<Pointer> back(T * S);
  std::vector<double> v(S), w(S);
  for (std::size_t k = 0; k < S; k++) v[k] = safe_log(m.initial[k]) + le(k, x[0]);
  for (std::size_t t = 1; t < T; t++) {
    Pointer* bp = &back[t * S];
    for (std::size_t l = 0; l < S; l++) {
      const double* row = &logT[l * S];
      double best = NEG_INF;
      std::size_t arg = 0;
      for (std::size_t k = 0; k < S; k++) {
        const double c = v[k] + row[k];
        if (c > best) { best = c; arg = k; }
      }
      w[l] = best + le(l, x[t]);
      bp[l] = static_cast<Pointer>(arg);
    }
    std::swap(v, w);
  }
  double best = NEG_INF;
  std::size_t arg = 0;
  for (std::size_t k = 0; k < S; k++) {
    const double c = v[k] + (m.ends() ? safe_log(m.end[k]) : 0.0);
    if (c > best) { best = c; arg = k; }
  }
  out.log_probability = best;
  if (best == NEG_INF) return out;
  out.states.resize(T);
  out.states[T - 1] = static_cast<std::uint32_t>(arg);
  for (std::size_t t = T - 1; t > 0; t--) out.states[t - 1] = static_cast<std::uint32_t>(back[t * S + out.states[t]]);
  return out;
}

}  // namespace detail

// The most probable path of x and its log probability.
inline Path viterbi(const Model& m, const Codes& x) {
  if (x.empty()) {
    Path p;
    p.log_probability = m.ends() ? NEG_INF : 0;
    return p;
  }
  if (m.states <= 256) return detail::viterbi_with<std::uint8_t>(m, x);
  if (m.states <= 65536) return detail::viterbi_with<std::uint16_t>(m, x);
  return detail::viterbi_with<std::uint32_t>(m, x);
}

// The log probability of a given path (and x along it): log P(x, π).
inline double path_log_probability(const Model& m, const Codes& x, const std::vector<std::uint32_t>& path) {
  if (path.size() != x.size()) throw std::invalid_argument("genomics::hmm::path_log_probability: the path and the sequence differ in length");
  if (x.empty()) return m.ends() ? NEG_INF : 0;
  double s = safe_log(m.initial[path[0]]) + safe_log(m.emit(path[0], x[0]));
  for (std::size_t t = 1; t < x.size(); t++) s += safe_log(m.a(path[t - 1], path[t])) + safe_log(m.emit(path[t], x[t]));
  if (m.ends()) s += safe_log(m.end[path.back()]);
  return s;
}

// ---- forward and backward ------------------------------------------------------------------------
struct Forward {
  std::vector<double> alpha;  // α̂_t(k), T × S, each row summing to one
  std::vector<double> scale;  // c_t
  double log_likelihood = NEG_INF;
};
struct Backward {
  std::vector<double> beta;   // β̂_t(k), T × S, each row summing to one
  std::vector<double> scale;  // d_t
  double log_likelihood = NEG_INF;
};

namespace detail {

// One forward pass; `store` keeps α̂ (T × S), otherwise only O(S) memory is used.
inline double forward_pass(const Model& m, const Codes& x, const std::vector<double>& aT, std::vector<double>* alpha, std::vector<double>* scales) {
  const std::size_t S = m.states, T = x.size();
  if (T == 0) return m.ends() ? NEG_INF : 0;
  std::vector<double> v(S), w(S);
  double logp = 0, c = 0;
  for (std::size_t k = 0; k < S; k++) c += v[k] = m.initial[k] * m.emit(k, x[0]);
  const auto finish_row = [&](std::vector<double>& row, std::size_t t) {
    if (scales) (*scales)[t] = c;
    if (!(c > 0)) return false;
    const double inv = 1 / c;
    for (auto& p : row) p *= inv;
    logp += std::log(c);
    if (alpha) std::copy(row.begin(), row.end(), alpha->begin() + static_cast<std::ptrdiff_t>(t * S));
    return true;
  };
  if (alpha) alpha->assign(T * S, 0);
  if (scales) scales->assign(T, 0);
  if (!finish_row(v, 0)) return NEG_INF;
  for (std::size_t t = 1; t < T; t++) {
    c = 0;
    for (std::size_t l = 0; l < S; l++) {
      const double* col = &aT[l * S];
      double s = 0;
      for (std::size_t k = 0; k < S; k++) s += v[k] * col[k];
      c += w[l] = s * m.emit(l, x[t]);
    }
    std::swap(v, w);
    if (!finish_row(v, t)) return NEG_INF;
  }
  if (m.ends()) {
    double z = 0;
    for (std::size_t k = 0; k < S; k++) z += v[k] * m.end[k];
    logp += safe_log(z);
  }
  return logp;
}

}  // namespace detail

inline Forward forward(const Model& m, const Codes& x) {
  Forward f;
  f.log_likelihood = detail::forward_pass(m, x, detail::transposed(m), &f.alpha, &f.scale);
  return f;
}

// log P(x) by the forward recursion in O(S) memory.
inline double log_likelihood(const Model& m, const Codes& x) { return detail::forward_pass(m, x, detail::transposed(m), nullptr, nullptr); }

inline Backward backward(const Model& m, const Codes& x) {
  const std::size_t S = m.states, T = x.size();
  Backward b;
  if (T == 0) {
    b.log_likelihood = m.ends() ? NEG_INF : 0;
    return b;
  }
  b.beta.assign(T * S, 0);
  b.scale.assign(T, 0);
  double logp = 0;
  // β_T(k) = a_k0 (or 1), normalised.
  {
    double d = 0;
    for (std::size_t k = 0; k < S; k++) d += b.beta[(T - 1) * S + k] = m.ends() ? m.end[k] : 1.0;
    b.scale[T - 1] = d;
    if (!(d > 0)) return b;
    for (std::size_t k = 0; k < S; k++) b.beta[(T - 1) * S + k] /= d;
    logp += std::log(d);
  }
  std::vector<double> u(S);
  for (std::size_t t = T - 1; t > 0; t--) {
    const double* next = &b.beta[t * S];
    for (std::size_t l = 0; l < S; l++) u[l] = m.emit(l, x[t]) * next[l];
    double d = 0;
    double* row = &b.beta[(t - 1) * S];
    for (std::size_t k = 0; k < S; k++) {
      const double* ak = &m.transition[k * S];
      double s = 0;
      for (std::size_t l = 0; l < S; l++) s += ak[l] * u[l];
      d += row[k] = s;
    }
    b.scale[t - 1] = d;
    if (!(d > 0)) return b;
    for (std::size_t k = 0; k < S; k++) row[k] /= d;
    logp += std::log(d);
  }
  double z = 0;
  for (std::size_t k = 0; k < S; k++) z += m.initial[k] * m.emit(k, x[0]) * b.beta[k];
  b.log_likelihood = logp + safe_log(z);
  return b;
}

// ---- posteriors ---------------------------------------------------------------------------------------
struct Posterior {
  std::size_t states = 0;
  std::vector<double> gamma;  // γ_t(k) = P(π_t = k | x), T × S
  double log_likelihood = NEG_INF;
  double at(std::size_t t, std::size_t k) const { return gamma[t * states + k]; }
};

inline Posterior posterior(const Model& m, const Codes& x) {
  const std::size_t S = m.states, T = x.size();
  Posterior p;
  p.states = S;
  const Forward f = forward(m, x);
  p.log_likelihood = f.log_likelihood;
  p.gamma.assign(T * S, 0);
  if (f.log_likelihood == NEG_INF) return p;
  const Backward b = backward(m, x);
  for (std::size_t t = 0; t < T; t++) {
    double s = 0;
    for (std::size_t k = 0; k < S; k++) s += p.gamma[t * S + k] = f.alpha[t * S + k] * b.beta[t * S + k];
    for (std::size_t k = 0; k < S; k++) p.gamma[t * S + k] /= s;
  }
  return p;
}

// The state of largest posterior at each position (ties to the smallest).
inline std::vector<std::uint32_t> posterior_decode(const Posterior& p) {
  const std::size_t S = p.states, T = S ? p.gamma.size() / S : 0;
  std::vector<std::uint32_t> out(T, 0);
  for (std::size_t t = 0; t < T; t++) {
    const double* row = &p.gamma[t * S];
    out[t] = static_cast<std::uint32_t>(std::max_element(row, row + S) - row);
  }
  return out;
}

// ---- Baum–Welch -----------------------------------------------------------------------------------------
struct Pseudocounts {
  double initial = 0, transition = 0, emission = 0, end = 0;
};
struct TrainOptions {
  int max_iterations = 200;
  double tolerance = 1e-8;  // stop when the log posterior rises by less than this (nats) per symbol
  Pseudocounts pseudocounts;
  bool initial = true, transition = true, emission = true, end = true;  // which parameters are re-estimated
  std::size_t chunk = 8;  // sequences per shard of the E-step (fixed, so results do not depend on threads)
};
struct Training {
  std::vector<double> log_likelihood;  // Σ log P(x | θ_i) before each update, and after the last
  std::vector<double> log_posterior;   // the same plus Σ r log θ: the quantity EM never decreases
  int iterations = 0;
  bool converged = false;
};

namespace detail {

struct Counts {
  std::vector<double> initial, transition, emission, end;
  double log_likelihood = 0;
  void reset(const Model& m) {
    initial.assign(m.states, 0);
    transition.assign(m.states * m.states, 0);
    emission.assign(m.states * m.symbols, 0);
    end.assign(m.states, 0);
    log_likelihood = 0;
  }
  void add(const Counts& o) {
    for (std::size_t i = 0; i < initial.size(); i++) initial[i] += o.initial[i];
    for (std::size_t i = 0; i < transition.size(); i++) transition[i] += o.transition[i];
    for (std::size_t i = 0; i < emission.size(); i++) emission[i] += o.emission[i];
    for (std::size_t i = 0; i < end.size(); i++) end[i] += o.end[i];
    log_likelihood += o.log_likelihood;
  }
};

// Adds the expected counts of one sequence.
inline void expected_counts(const Model& m, const Codes& x, Counts& c) {
  const std::size_t S = m.states, T = x.size();
  if (T == 0) return;
  const Forward f = forward(m, x);
  c.log_likelihood += f.log_likelihood;
  if (f.log_likelihood == NEG_INF) return;
  const Backward b = backward(m, x);
  std::vector<double> g(S), xi(S * S), u(S);
  for (std::size_t t = 0; t < T; t++) {
    double s = 0;
    for (std::size_t k = 0; k < S; k++) s += g[k] = f.alpha[t * S + k] * b.beta[t * S + k];
    for (std::size_t k = 0; k < S; k++) g[k] /= s;
    if (t == 0)
      for (std::size_t k = 0; k < S; k++) c.initial[k] += g[k];
    if (t == T - 1)
      for (std::size_t k = 0; k < S; k++) c.end[k] += g[k];
    if (x[t] < m.symbols)
      for (std::size_t k = 0; k < S; k++) c.emission[k * m.symbols + x[t]] += g[k];
    if (t + 1 < T) {
      for (std::size_t l = 0; l < S; l++) u[l] = m.emit(l, x[t + 1]) * b.beta[(t + 1) * S + l];
      double z = 0;
      for (std::size_t k = 0; k < S; k++) {
        const double ak = f.alpha[t * S + k];
        const double* row = &m.transition[k * S];
        for (std::size_t l = 0; l < S; l++) z += xi[k * S + l] = ak * row[l] * u[l];
      }
      const double inv = 1 / z;
      for (std::size_t i = 0; i < S * S; i++) c.transition[i] += xi[i] * inv;
    }
  }
}

// Σ r log θ over the parameters that receive pseudocounts.
inline double log_prior(const Model& m, const TrainOptions& o) {
  double s = 0;
  const auto add = [&](const std::vector<double>& v, double r, bool on) {
    if (!on || r == 0) return;
    for (double p : v)
      if (p > 0) s += r * std::log(p);
  };
  add(m.initial, o.pseudocounts.initial, o.initial);
  add(m.transition, o.pseudocounts.transition, o.transition);
  add(m.emission, o.pseudocounts.emission, o.emission);
  add(m.end, o.pseudocounts.end, o.end && m.ends());
  return s;
}

}  // namespace detail

// Re-estimates m from the sequences until the log posterior rises by less than the tolerance,
// or for max_iterations. The E-step runs on `pool` when one is given.
inline Training baum_welch(Model& m, const std::vector<Codes>& sequences, const TrainOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  m.validate(1e-6);
  const std::size_t S = m.states, K = m.symbols;
  std::size_t symbols = 0;
  for (const auto& x : sequences) symbols += x.size();
  const std::size_t chunk = std::max<std::size_t>(1, o.chunk), shards = (sequences.size() + chunk - 1) / chunk;
  std::vector<detail::Counts> parts(shards);
  detail::Counts total;
  const auto e_step = [&] {
    const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
      detail::Counts& c = parts[begin / chunk];
      c.reset(m);
      for (std::size_t i = begin; i < end; i++) detail::expected_counts(m, sequences[i], c);
    };
    if (pool) calx::scheduler::parallel_for(*pool, sequences.size(), chunk, body);
    else
      for (std::size_t s = 0; s < shards; s++) body(s * chunk, std::min(sequences.size(), (s + 1) * chunk), 0);
    total.reset(m);
    for (const auto& c : parts) total.add(c);
  };
  // Pseudocounts only where the model allows a transition or an emission.
  const auto estimate = [](std::vector<double>& p, const std::vector<double>& n, std::size_t row, double r, double* extra_p, double extra_n, double extra_r) {
    const std::size_t rows = row ? p.size() / row : 0;
    for (std::size_t k = 0; k < rows; k++) {
      double s = 0;
      for (std::size_t j = 0; j < row; j++)
        if (p[k * row + j] > 0) s += n[k * row + j] + r;
      const bool with_extra = extra_p && extra_p[k] > 0;
      if (with_extra) s += extra_n + extra_r;
      if (!(s > 0)) continue;
      for (std::size_t j = 0; j < row; j++)
        if (p[k * row + j] > 0) p[k * row + j] = (n[k * row + j] + r) / s;
      if (with_extra) extra_p[k] = (extra_n + extra_r) / s;
    }
  };
  Training out;
  const double tol = o.tolerance * static_cast<double>(std::max<std::size_t>(1, symbols));
  for (int it = 0; it < o.max_iterations; it++) {
    e_step();
    out.log_likelihood.push_back(total.log_likelihood);
    out.log_posterior.push_back(total.log_likelihood + detail::log_prior(m, o));
    if (total.log_likelihood == NEG_INF) break;
    if (it > 0 && out.log_posterior[it] - out.log_posterior[it - 1] < tol) {
      out.converged = true;
      out.iterations = it;
      return out;
    }
    if (o.initial) estimate(m.initial, total.initial, S, o.pseudocounts.initial, nullptr, 0, 0);
    if (o.transition || (o.end && m.ends())) {
      // A row of transitions and its end are one distribution; a part not re-estimated keeps its
      // mass, and the other part is re-estimated within the rest.
      for (std::size_t k = 0; k < S; k++) {
        double* row = &m.transition[k * S];
        const double* n = &total.transition[k * S];
        double* e = m.ends() ? &m.end[k] : nullptr;
        if (o.transition && (!e || o.end)) {
          std::vector<double> p(row, row + S), c(n, n + S);
          double ep = e ? *e : 0;
          estimate(p, c, S, o.pseudocounts.transition, e ? &ep : nullptr, total.end[k], o.pseudocounts.end);
          std::copy(p.begin(), p.end(), row);
          if (e) *e = ep;
        } else if (o.transition) {
          double s = 0, kept = e ? *e : 0;
          for (std::size_t l = 0; l < S; l++)
            if (row[l] > 0) s += n[l] + o.pseudocounts.transition;
          if (s > 0)
            for (std::size_t l = 0; l < S; l++)
              if (row[l] > 0) row[l] = (1 - kept) * (n[l] + o.pseudocounts.transition) / s;
        } else if (e && o.end) {
          // Only the end: its expected frequency among the row's transitions.
          double s = total.end[k] + o.pseudocounts.end;
          for (std::size_t l = 0; l < S; l++) s += n[l];
          if (s > 0 && *e > 0) {
            const double ne = (total.end[k] + o.pseudocounts.end) / s, scale = (1 - ne) / (1 - *e);
            for (std::size_t l = 0; l < S; l++) row[l] *= scale;
            *e = ne;
          }
        }
      }
    }
    if (o.emission) estimate(m.emission, total.emission, K, o.pseudocounts.emission, nullptr, 0, 0);
    out.iterations = it + 1;
  }
  e_step();
  out.log_likelihood.push_back(total.log_likelihood);
  out.log_posterior.push_back(total.log_likelihood + detail::log_prior(m, o));
  const std::size_t n = out.log_posterior.size();
  out.converged = n >= 2 && out.log_posterior[n - 1] - out.log_posterior[n - 2] < tol;
  return out;
}

// ---- sampling -------------------------------------------------------------------------------------------
struct Sample {
  std::vector<std::uint32_t> states;
  Codes symbols;
};

namespace detail {

inline std::size_t draw(calx::random::Stream& r, const double* p, std::size_t n, double total) {
  double u = r.uniform() * total;
  for (std::size_t i = 0; i < n; i++) {
    if (u < p[i]) return i;
    u -= p[i];
  }
  for (std::size_t i = n; i-- > 0;)
    if (p[i] > 0) return i;
  return n - 1;
}

}  // namespace detail

// A path of exactly `length` states and its symbols; the end state, if any, is not consulted
// (the transitions out of a state are drawn conditioned on not ending).
inline Sample sample(const Model& m, std::size_t length, calx::random::Stream& r) {
  Sample s;
  s.states.reserve(length);
  s.symbols.reserve(length);
  for (std::size_t t = 0; t < length; t++) {
    std::size_t k;
    if (t == 0) k = detail::draw(r, m.initial.data(), m.states, 1.0);
    else {
      const double* row = &m.transition[s.states.back() * m.states];
      double total = 0;
      for (std::size_t l = 0; l < m.states; l++) total += row[l];
      k = detail::draw(r, row, m.states, total);
    }
    s.states.push_back(static_cast<std::uint32_t>(k));
    s.symbols.push_back(static_cast<Code>(detail::draw(r, &m.emission[k * m.symbols], m.symbols, 1.0)));
  }
  return s;
}

// A path to the end state (the model must have one), of at most max_length states.
inline Sample sample_to_end(const Model& m, calx::random::Stream& r, std::size_t max_length) {
  if (!m.ends()) throw std::invalid_argument("genomics::hmm::sample_to_end: the model has no end state");
  Sample s;
  std::vector<double> row(m.states + 1);
  std::size_t k = detail::draw(r, m.initial.data(), m.states, 1.0);
  while (s.states.size() < max_length) {
    s.states.push_back(static_cast<std::uint32_t>(k));
    s.symbols.push_back(static_cast<Code>(detail::draw(r, &m.emission[k * m.symbols], m.symbols, 1.0)));
    std::copy(&m.transition[k * m.states], &m.transition[k * m.states] + m.states, row.begin());
    row[m.states] = m.end[k];
    const std::size_t next = detail::draw(r, row.data(), m.states + 1, 1.0);
    if (next == m.states) break;
    k = next;
  }
  return s;
}

}  // namespace genomics::hmm
