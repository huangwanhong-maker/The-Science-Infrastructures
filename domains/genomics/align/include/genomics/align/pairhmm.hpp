// The pair hidden Markov model of alignment (Durbin, Eddy, Krogh and Mitchison 1998, Biological
// Sequence Analysis, ch. 4): its most probable alignment (Viterbi), the probability of the two
// sequences summed over all alignments (forward and backward, in log space), the posterior
// probability that two letters are aligned, and the alignment of maximum expected accuracy.
//
// The model (Durbin et al. 1998, §4.1, Fig. 4.2). States M (a pair x_i ◇ y_j emitted with
// probability p(x_i, y_j)), X (x_i against a gap, emitted with q(x_i)) and Y (y_j against a gap,
// q(y_j)); a begin that behaves as M and an end. Transitions: M → M 1 − 2δ − τ, M → X and M → Y δ,
// X → X and Y → Y ε, X → M and Y → M 1 − ε − τ, every state → end τ; X and Y never follow each
// other. In alignment.hpp's terms it is the three-state scheme with these transitions and
// emissions, so that
//
//   Viterbi  is the sum over paths in the tropical semiring (max, +) of the log probabilities,
//   forward  is the sum in the probability semiring (+, ×) (here as log-sum-exp in log space),
//
// and the checks hold each to morphe's generic programme with those weights.
//
// Recurrences (Durbin et al. 1998, §4.2, the Viterbi algorithm for pair HMMs; §4.3, the forward
// algorithm; §4.4, the backward algorithm and the posterior of an aligned pair, P(x_i ◇ y_j | x, y)
// = f^M(i, j) b^M(i, j) / P(x, y)), with f^M(0, 0) = 1 standing for the begin. Maximum expected
// accuracy (Holmes and Durbin 1998, J. Comput. Biol. 5, 493–504; Do, Mahabhashyam, Brudno and
// Batzoglou 2005, ProbCons, Genome Res. 15, 330–340): the alignment maximising the sum of the
// posteriors of its aligned pairs, by A(i, j) = max(A(i − 1, j − 1) + P(i, j), A(i − 1, j),
// A(i, j − 1)).
//
// Numbers. Every quantity is a natural logarithm; sums of probabilities are log-sum-exp,
// log(e^a + e^b) = max + log1p(e^{−|a − b|}), exact to rounding, so that sequences of any length
// neither underflow nor overflow. Complexity: O(mn) time each; Viterbi and the tables of forward
// and backward keep O(mn) doubles (three per vertex), the posterior one.

#pragma once

#include <genomics/align/types.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::align {

struct PairHmm {
  int size = 0;           // the alphabet: codes 0 .. size − 1
  std::vector<double> p;  // size × size, p[a · size + b] = P(a ◇ b) in M (x's letter a, y's b)
  std::vector<double> q;  // size, the emission of X and Y
  double delta = 0.2, epsilon = 0.1, tau = 0.1;

  // DNA with P(a ◇ a) = identity / 4 and the other pairs (1 − identity) / 12, q uniform.
  static PairHmm dna(double identity, double delta, double epsilon, double tau) {
    PairHmm h;
    h.size = 4;
    h.p.assign(16, (1 - identity) / 12);
    for (int a = 0; a < 4; a++) h.p[static_cast<std::size_t>(a * 5)] = identity / 4;
    h.q.assign(4, 0.25);
    h.delta = delta;
    h.epsilon = epsilon;
    h.tau = tau;
    return h;
  }
  // The target frequencies of a substitution matrix, p_ab ∝ q_a q_b e^{λ s_ab} (Altschul 1991,
  // J. Mol. Biol. 219, 555–565: a matrix scaled by λ is the log-odds of these frequencies),
  // normalised; with λ the Karlin–Altschul root (statistics.hpp) they sum to 1 already.
  static PairHmm from_matrix(const Matrix& m, std::span<const double> background, double lambda, double delta, double epsilon, double tau) {
    PairHmm h;
    h.size = m.size;
    h.p.assign(static_cast<std::size_t>(m.size * m.size), 0.0);
    h.q.assign(background.begin(), background.begin() + m.size);
    double qt = 0, pt = 0;
    for (double x : h.q) qt += x;
    for (auto& x : h.q) x /= qt;
    for (int a = 0; a < m.size; a++)
      for (int b = 0; b < m.size; b++) pt += h.p[static_cast<std::size_t>(a * m.size + b)] = h.q[static_cast<std::size_t>(a)] * h.q[static_cast<std::size_t>(b)] * std::exp(lambda * m(static_cast<Code>(a), static_cast<Code>(b)));
    for (auto& x : h.p) x /= pt;
    h.delta = delta;
    h.epsilon = epsilon;
    h.tau = tau;
    return h;
  }

  void validate(const Codes& x, const Codes& y) const {
    if (size <= 0 || p.size() != static_cast<std::size_t>(size * size) || q.size() != static_cast<std::size_t>(size)) throw std::invalid_argument("genomics::align::PairHmm: emissions of the wrong size");
    if (!(delta > 0 && epsilon > 0 && tau > 0 && 2 * delta + tau < 1 && epsilon + tau < 1)) throw std::invalid_argument("genomics::align::PairHmm: need δ, ε, τ > 0, 2δ + τ < 1, ε + τ < 1");
    for (Code c : x)
      if (c >= size) throw std::invalid_argument("genomics::align::PairHmm: a code of x outside the alphabet");
    for (Code c : y)
      if (c >= size) throw std::invalid_argument("genomics::align::PairHmm: a code of y outside the alphabet");
  }
  double log_p(Code a, Code b) const { return std::log(p[static_cast<std::size_t>(a) * static_cast<std::size_t>(size) + b]); }
  double log_q(Code a) const { return std::log(q[a]); }
};

inline constexpr double LOG_ZERO = -std::numeric_limits<double>::infinity();

// log(e^a + e^b), exactly to rounding, −∞ the identity.
inline double log_sum_exp(double a, double b) {
  if (a == LOG_ZERO) return b;
  if (b == LOG_ZERO) return a;
  return a > b ? a + std::log1p(std::exp(b - a)) : b + std::log1p(std::exp(a - b));
}
inline double log_sum_exp(double a, double b, double c) { return log_sum_exp(log_sum_exp(a, b), c); }

// Log-space tables of the three states over the vertices (i, j), row-major, (m + 1)(n + 1) each.
struct PairTables {
  std::size_t m = 0, n = 0;
  std::vector<double> M, X, Y;
  double log_probability = LOG_ZERO;  // log P(x, y)
  std::size_t at(std::size_t i, std::size_t j) const { return i * (n + 1) + j; }
};

namespace detail {

struct LogTransitions {
  double mm, mx, xm, xx, end;
  explicit LogTransitions(const PairHmm& h)
      : mm(std::log(1 - 2 * h.delta - h.tau)), mx(std::log(h.delta)), xm(std::log(1 - h.epsilon - h.tau)), xx(std::log(h.epsilon)), end(std::log(h.tau)) {}
};

}  // namespace detail

// ---- Viterbi ----------------------------------------------------------------------------------------
struct PairViterbi {
  Alignment alignment;  // global; its score field unused (0)
  double log_probability = LOG_ZERO;
};

inline PairViterbi viterbi(const PairHmm& h, const Codes& x, const Codes& y) {
  h.validate(x, y);
  const detail::LogTransitions t(h);
  PairTables v;
  v.m = x.size();
  v.n = y.size();
  const std::size_t cells = (v.m + 1) * (v.n + 1);
  v.M.assign(cells, LOG_ZERO);
  v.X.assign(cells, LOG_ZERO);
  v.Y.assign(cells, LOG_ZERO);
  v.M[0] = 0;
  for (std::size_t i = 0; i <= v.m; i++)
    for (std::size_t j = 0; j <= v.n; j++) {
      if (i == 0 && j == 0) continue;
      const std::size_t c = v.at(i, j);
      if (i && j) {
        const std::size_t d = v.at(i - 1, j - 1);
        v.M[c] = h.log_p(x[i - 1], y[j - 1]) + std::max({t.mm + v.M[d], t.xm + v.X[d], t.xm + v.Y[d]});
      }
      if (i) {
        const std::size_t u = v.at(i - 1, j);
        v.X[c] = h.log_q(x[i - 1]) + std::max(t.mx + v.M[u], t.xx + v.X[u]);
      }
      if (j) {
        const std::size_t l = v.at(i, j - 1);
        v.Y[c] = h.log_q(y[j - 1]) + std::max(t.mx + v.M[l], t.xx + v.Y[l]);
      }
    }
  PairViterbi out;
  const std::size_t e = v.at(v.m, v.n);
  out.log_probability = t.end + std::max({v.M[e], v.X[e], v.Y[e]});
  // Traceback: M before X before Y on ties. The candidates are recomputed by the same expressions
  // as in the fill, so the one the maximum chose compares equal exactly.
  const auto first_max = [](double a, double b, double c) { return a >= b && a >= c ? 0 : (b >= c ? 1 : 2); };
  int state = first_max(v.M[e], v.X[e], v.Y[e]);
  std::size_t i = v.m, j = v.n;
  Cigar rev;
  while (i > 0 || j > 0) {
    if (state == 0) {
      push(rev, 'M');
      const std::size_t d = v.at(i - 1, j - 1);
      state = first_max(t.mm + v.M[d], t.xm + v.X[d], t.xm + v.Y[d]);
      i--;
      j--;
    } else if (state == 1) {
      push(rev, 'I');
      const std::size_t u = v.at(i - 1, j);
      state = t.mx + v.M[u] >= t.xx + v.X[u] ? 0 : 1;
      i--;
    } else {
      push(rev, 'D');
      const std::size_t l = v.at(i, j - 1);
      state = t.mx + v.M[l] >= t.xx + v.Y[l] ? 0 : 2;
      j--;
    }
  }
  out.alignment.query_end = v.m;
  out.alignment.target_end = v.n;
  out.alignment.cigar = reversed(std::move(rev));
  return out;
}

// ---- forward and backward ------------------------------------------------------------------------------
inline PairTables forward(const PairHmm& h, const Codes& x, const Codes& y) {
  h.validate(x, y);
  const detail::LogTransitions t(h);
  PairTables f;
  f.m = x.size();
  f.n = y.size();
  const std::size_t cells = (f.m + 1) * (f.n + 1);
  f.M.assign(cells, LOG_ZERO);
  f.X.assign(cells, LOG_ZERO);
  f.Y.assign(cells, LOG_ZERO);
  f.M[0] = 0;
  for (std::size_t i = 0; i <= f.m; i++)
    for (std::size_t j = 0; j <= f.n; j++) {
      if (i == 0 && j == 0) continue;
      const std::size_t c = f.at(i, j);
      if (i && j) {
        const std::size_t d = f.at(i - 1, j - 1);
        f.M[c] = h.log_p(x[i - 1], y[j - 1]) + log_sum_exp(t.mm + f.M[d], t.xm + f.X[d], t.xm + f.Y[d]);
      }
      if (i) {
        const std::size_t u = f.at(i - 1, j);
        f.X[c] = h.log_q(x[i - 1]) + log_sum_exp(t.mx + f.M[u], t.xx + f.X[u]);
      }
      if (j) {
        const std::size_t l = f.at(i, j - 1);
        f.Y[c] = h.log_q(y[j - 1]) + log_sum_exp(t.mx + f.M[l], t.xx + f.Y[l]);
      }
    }
  const std::size_t e = f.at(f.m, f.n);
  f.log_probability = t.end + log_sum_exp(f.M[e], f.X[e], f.Y[e]);
  return f;
}

// b^s(i, j): the log probability of what follows (i, j) in state s; log P(x, y) = b^M(0, 0).
inline PairTables backward(const PairHmm& h, const Codes& x, const Codes& y) {
  h.validate(x, y);
  const detail::LogTransitions t(h);
  PairTables b;
  b.m = x.size();
  b.n = y.size();
  const std::size_t cells = (b.m + 1) * (b.n + 1);
  b.M.assign(cells, LOG_ZERO);
  b.X.assign(cells, LOG_ZERO);
  b.Y.assign(cells, LOG_ZERO);
  for (std::size_t i = b.m + 1; i-- > 0;)
    for (std::size_t j = b.n + 1; j-- > 0;) {
      const std::size_t c = b.at(i, j);
      if (i == b.m && j == b.n) {
        b.M[c] = b.X[c] = b.Y[c] = t.end;
        continue;
      }
      double pair = LOG_ZERO, xs = LOG_ZERO, ys = LOG_ZERO;  // what follows: a pair, an X step, a Y step
      if (i < b.m && j < b.n) pair = h.log_p(x[i], y[j]) + b.M[b.at(i + 1, j + 1)];
      if (i < b.m) xs = h.log_q(x[i]) + b.X[b.at(i + 1, j)];
      if (j < b.n) ys = h.log_q(y[j]) + b.Y[b.at(i, j + 1)];
      b.M[c] = log_sum_exp(t.mm + pair, t.mx + xs, t.mx + ys);
      b.X[c] = log_sum_exp(t.xm + pair, t.xx + xs);
      b.Y[c] = log_sum_exp(t.xm + pair, t.xx + ys);
    }
  b.log_probability = b.M[0];
  return b;
}

// P(x_i ◇ y_j | x, y) for 1 ≤ i ≤ m, 1 ≤ j ≤ n, as an m × n row-major matrix (entry (i − 1, j − 1)).
inline std::vector<double> posterior(const PairHmm& h, const Codes& x, const Codes& y) {
  const PairTables f = forward(h, x, y), b = backward(h, x, y);
  std::vector<double> post(x.size() * y.size(), 0.0);
  for (std::size_t i = 1; i <= x.size(); i++)
    for (std::size_t j = 1; j <= y.size(); j++) post[(i - 1) * y.size() + (j - 1)] = std::exp(f.M[f.at(i, j)] + b.M[b.at(i, j)] - f.log_probability);
  return post;
}

// ---- maximum expected accuracy --------------------------------------------------------------------------
struct ExpectedAccuracy {
  Alignment alignment;   // global; its score field unused (0)
  double expected_pairs = 0;  // the sum of the posteriors of its aligned pairs
  double accuracy = 0;        // expected_pairs / min(m, n), ProbCons's expected accuracy
};

inline ExpectedAccuracy mea_alignment(const std::vector<double>& post, std::size_t m, std::size_t n) {
  if (post.size() != m * n) throw std::invalid_argument("genomics::align::mea_alignment: the posterior is not m × n");
  std::vector<double> A((m + 1) * (n + 1), 0.0);
  const auto at = [n](std::size_t i, std::size_t j) { return i * (n + 1) + j; };
  for (std::size_t i = 1; i <= m; i++)
    for (std::size_t j = 1; j <= n; j++) A[at(i, j)] = std::max({A[at(i - 1, j - 1)] + post[(i - 1) * n + (j - 1)], A[at(i - 1, j)], A[at(i, j - 1)]});
  ExpectedAccuracy out;
  out.expected_pairs = A[at(m, n)];
  out.accuracy = std::min(m, n) ? out.expected_pairs / static_cast<double>(std::min(m, n)) : 0;
  std::size_t i = m, j = n;
  Cigar rev;
  while (i > 0 || j > 0) {
    if (i && j && A[at(i, j)] == A[at(i - 1, j - 1)] + post[(i - 1) * n + (j - 1)]) { push(rev, 'M'); i--; j--; }
    else if (i && A[at(i, j)] == A[at(i - 1, j)]) { push(rev, 'I'); i--; }
    else { push(rev, 'D'); j--; }
  }
  out.alignment.query_end = m;
  out.alignment.target_end = n;
  out.alignment.cigar = reversed(std::move(rev));
  return out;
}

}  // namespace genomics::align
