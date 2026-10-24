// The spin foam of a history computed exactly (theory note §6.6): physis's contraction of the ℤ₂
// lattice gauge theory (physis/gauge/contraction.hpp) along the order of time, and the couplings by
// maximum a posteriori from the exact likelihood and its derivatives (calx jets, through physis).
//
// The domain contributes what is the record's, and nothing of the method:
//
//   the theory, in temporal gauge: its links the rungs, its 2-cells the lenses and the triangles
//     with their couplings, its fixed links the attested signs when the record is conditioned on;
//   the order of time: at each occasion, each tie's lens from its last meeting as its new rung
//     arrives, then the occasion's triangles; the old rung of a tie leaves as its new one enters,
//     so a cut holds one rung for each open tie, and its state lives on the ℤ/2 cycles of the
//     network of that moment (physis's Proposition: the prior's states are spin networks);
//   cuts at dates, and the counts of frustrated triangles read there;
//   the statistics of the couplings: the Jeffreys prior on each coupling's rate of change, the
//     maximum a posteriori and its Laplace errors. They belong to the inference foundation when it
//     exists (README §4); until then they are here.

#pragma once

#include <history/construction/complex.hpp>
#include <history/foam/model.hpp>
#include <calx/linalg/dense.hpp>
#include <physis/gauge/contraction.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace history {

using physis::gauge::Contraction;
using physis::gauge::ContractionLimits;
using physis::gauge::Intractable;

// The theory of a model in temporal gauge: links the rungs, 2-cells the cells.
inline physis::gauge::Z2Theory theory_of(const HistoryModel& M, bool conditioned) {
  physis::gauge::Z2Theory T;
  T.links = M.rung.size();
  for (const auto& c : M.cells) {
    T.plaquettes.push_back(c.rungs);
    T.beta.push_back(c.beta);
  }
  if (conditioned) T.fixed = M.attested;
  return T;
}

// The order of time: cells by the occasion of their latest rung; there, the lenses in the order of
// their new rungs, then the triangles.
inline std::vector<std::size_t> time_order(const HistoryComplex& K, const HistoryModel& M) {
  auto when = [&](std::size_t r) { return K.vertex[K.rungs[M.rung[r][0]][M.rung[r][1]].occasion]; };
  std::vector<std::tuple<std::size_t, int, std::size_t, std::size_t>> key;
  for (std::size_t c = 0; c < M.cells.size(); c++) {
    std::size_t latest = M.cells[c].rungs.front();
    for (std::size_t r : M.cells[c].rungs)
      if (when(r) > when(latest) || (when(r) == when(latest) && r > latest)) latest = r;
    key.push_back({when(latest), M.cells[c].triangle ? 1 : 0, latest, c});
  }
  std::sort(key.begin(), key.end());
  std::vector<std::size_t> order;
  for (const auto& k : key) order.push_back(std::get<3>(k));
  return order;
}

// A cut after the vertex `after` of K (the occasion K.order[after]): the rungs it reads stay live.
struct Cut {
  std::size_t after = 0;
  std::vector<std::size_t> keep;
};

// The exact ensemble of a model: the prior (conditioned = false) or the posterior. Two-points are
// pairs of rungs; counters[q], the triangles counted at cuts[q]; `networks`, how many of the
// heaviest spin networks to give at each cut.
inline Contraction contract(const HistoryComplex& K, const HistoryModel& M, bool conditioned,
                            const std::vector<std::array<std::size_t, 2>>& two_points = {},
                            const std::vector<Cut>& cuts = {}, const std::vector<Counter>& counters = {},
                            std::size_t networks = 0, const ContractionLimits& limits = {}) {
  const auto order = time_order(K, M);
  auto when = [&](std::size_t c) {
    std::size_t v = 0;
    for (std::size_t r : M.cells[c].rungs) v = std::max(v, K.vertex[K.rungs[M.rung[r][0]][M.rung[r][1]].occasion]);
    return v;
  };
  std::vector<physis::gauge::ContractionCut> at;
  for (std::size_t q = 0; q < cuts.size(); q++) {
    physis::gauge::ContractionCut c;
    while (c.after < order.size() && when(order[c.after]) <= cuts[q].after) c.after++;
    c.keep = cuts[q].keep;
    if (q < counters.size())
      for (const auto& t : counters[q].triangles) c.products.push_back({t[0], t[1], t[2]});
    at.push_back(std::move(c));
  }
  return physis::gauge::contract(theory_of(M, conditioned), order, two_points, at, networks, limits);
}

// log Z of the prior or the posterior alone, by one forward pass: the test of reach, and the
// amplitude; `largest` receives the largest state.
inline double log_z(const HistoryComplex& K, const HistoryModel& M, bool conditioned, std::size_t& largest, const ContractionLimits& limits = {}) {
  return physis::gauge::log_z(theory_of(M, conditioned), time_order(K, M), limits, largest);
}

// ---- the couplings by maximum a posteriori (theory note §6.1) ------------------------------------

// The parameters of the weight: one persistence β for each colour with lenses, and the balance β_Δ
// if there are triangles; each cell's parameter.
struct Parameters {
  std::vector<std::string> names;
  std::vector<std::size_t> of_cell;
};

inline Parameters parameters_of(const HistoryComplex& K, const HistoryModel& M) {
  Parameters P;
  std::map<std::string, std::size_t> index;
  auto id = [&](const std::string& name) {
    const auto it = index.find(name);
    if (it != index.end()) return it->second;
    index[name] = P.names.size();
    P.names.push_back(name);
    return P.names.size() - 1;
  };
  for (std::size_t c = 0; c < M.cells.size(); c++)
    P.of_cell.push_back(M.cells[c].triangle ? id("balance") : id(K.record.ties[K.lenses[c].tie].colour));
  return P;
}

struct Fit {
  std::vector<std::string> names;
  std::vector<double> beta, error;      // the MAP and its standard error (Laplace)
  std::vector<double> start;            // the KT estimates it started from
  std::vector<std::size_t> cells;       // cells of each parameter
  std::vector<double> gradient;         // at the MAP: zero
  double log_likelihood = 0;            // log P(attested signs | β̂) = log Z_post − log Z_prior
  double log_posterior = 0;             // with the prior, −Σ log cosh β
  std::size_t iterations = 0;
  bool converged = false;
};

namespace transfer_detail {

inline HistoryModel with_betas(HistoryModel M, const Parameters& P, const std::vector<double>& beta) {
  for (std::size_t c = 0; c < M.cells.size(); c++) M.cells[c].beta = beta[P.of_cell[c]];
  return M;
}

struct Evaluation {
  double objective = 0, log_likelihood = 0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

// Φ(β) = log Z_post − log Z_prior − Σ log cosh β_p, with its gradient (order ≥ 1) and Hessian
// (order 2), by physis's jets with N = the number of parameters.
template <std::size_t N>
inline Evaluation evaluate_n(const HistoryComplex& K, const HistoryModel& M, const Parameters& P, const std::vector<double>& beta,
                             const ContractionLimits& limits, int order) {
  const auto prior = theory_of(M, false), post = theory_of(M, true);
  const auto sequence = time_order(K, M);
  Evaluation e;
  if (order == 0) {
    e.log_likelihood = physis::gauge::log_z(post, sequence, limits) - physis::gauge::log_z(prior, sequence, limits);
  } else {
    std::array<double, N> par{};
    for (std::size_t p = 0; p < N; p++) par[p] = beta[p];
    auto run = [&](const physis::gauge::Z2Theory& T) {
      if (order == 2) return physis::gauge::log_z_derivatives<N, 2>(T, sequence, P.of_cell, par, limits);
      return physis::gauge::log_z_derivatives<N, 1>(T, sequence, P.of_cell, par, limits);
    };
    const auto a = run(prior), b = run(post);
    e.log_likelihood = b.log_z - a.log_z;
    e.gradient.assign(N, 0);
    for (std::size_t p = 0; p < N; p++) e.gradient[p] = b.gradient[p] - a.gradient[p] - std::tanh(beta[p]);
    if (order == 2) {
      e.hessian.assign(N, std::vector<double>(N, 0));
      for (std::size_t p = 0; p < N; p++)
        for (std::size_t q = 0; q < N; q++) e.hessian[p][q] = b.hessian[p][q] - a.hessian[p][q];
      for (std::size_t p = 0; p < N; p++) e.hessian[p][p] -= 1 - std::tanh(beta[p]) * std::tanh(beta[p]);
    }
  }
  e.objective = e.log_likelihood;
  for (double b : beta) e.objective -= physis::gauge::log_cosh(b);
  return e;
}

inline Evaluation evaluate(const HistoryComplex& K, const HistoryModel& M0, const Parameters& P, const std::vector<double>& beta,
                           const ContractionLimits& limits = {}, int order = 1) {
  const HistoryModel M = with_betas(M0, P, beta);
  switch (beta.size()) {
    case 1: return evaluate_n<1>(K, M, P, beta, limits, order);
    case 2: return evaluate_n<2>(K, M, P, beta, limits, order);
    case 3: return evaluate_n<3>(K, M, P, beta, limits, order);
    case 4: return evaluate_n<4>(K, M, P, beta, limits, order);
    case 5: return evaluate_n<5>(K, M, P, beta, limits, order);
    case 6: return evaluate_n<6>(K, M, P, beta, limits, order);
    case 7: return evaluate_n<7>(K, M, P, beta, limits, order);
    case 8: return evaluate_n<8>(K, M, P, beta, limits, order);
    default: throw std::invalid_argument("history: between one and eight couplings can be fitted");
  }
}

// The Newton step from the gradient and Hessian (calx linalg), damped toward the gradient until it
// rises: a step s with (−H + μ) s = g and g · s > 0.
inline std::vector<double> newton_step(const std::vector<std::vector<double>>& H, const std::vector<double>& g) {
  const std::size_t m = g.size();
  for (double mu = 0;; mu = mu == 0 ? 1e-6 : 10 * mu) {
    calx::linalg::Matrix A(m), B(m);
    for (std::size_t p = 0; p < m; p++) {
      for (std::size_t q = 0; q < m; q++) A(p, q) = -H[p][q];
      A(p, p) += mu;
      B(p, 0) = g[p];
    }
    try {
      const auto X = calx::linalg::solve(A, B);
      std::vector<double> s(m);
      double rise = 0;
      for (std::size_t p = 0; p < m; p++) { s[p] = X(p, 0); rise += s[p] * g[p]; }
      if (rise > 0 || mu > 1e8) return s;
    } catch (const std::runtime_error&) {
      if (mu > 1e8) return g;
    }
  }
}

}  // namespace transfer_detail

// The maximum a posteriori couplings: the exact log-likelihood of the attested signs,
// log Z_post(β) − log Z_prior(β), plus the Jeffreys prior of each coupling's rate of change
// p = (1 − tanh β)/2, which is −log cosh β up to a constant. Its gradient is
// Σ_{f of p} (⟨h_f⟩_post − ⟨h_f⟩_prior) − tanh β_p and its Hessian Cov_post − Cov_prior − sech² β,
// both exact. Newton from the KT estimates, a full step tried first and halved until Φ rises; the
// standard errors from the inverse of minus the Hessian at the maximum (Laplace). When the cells
// are independent (no rung shared), the maximum is the KT estimate (k + ½)/(n + 1) of each
// coupling (Proposition 6.4).
inline Fit fit_couplings(const HistoryComplex& K, const HistoryModel& M, const Couplings& kt, const ContractionLimits& limits = {}) {
  using namespace transfer_detail;
  const Parameters P = parameters_of(K, M);
  const std::size_t m = P.names.size();
  Fit fit;
  fit.names = P.names;
  fit.cells.assign(m, 0);
  for (std::size_t p : P.of_cell) fit.cells[p]++;
  std::vector<double> beta(m);
  for (std::size_t p = 0; p < m; p++) beta[p] = P.names[p] == "balance" ? kt.balance : kt.of(P.names[p]);
  fit.start = beta;
  Evaluation e = evaluate(K, M, P, beta, limits, 2);
  for (fit.iterations = 0; fit.iterations < 40; fit.iterations++) {
    double gmax = 0;
    for (double g : e.gradient) gmax = std::max(gmax, std::fabs(g));
    if (gmax < 1e-8) { fit.converged = true; break; }
    const auto step = newton_step(e.hessian, e.gradient);
    auto at = [&](double length) {
      std::vector<double> next(m);
      for (std::size_t p = 0; p < m; p++) next[p] = beta[p] + length * step[p];
      return next;
    };
    Evaluation full = evaluate(K, M, P, at(1), limits, 2);
    if (full.objective >= e.objective - 1e-12) {
      beta = at(1);
      e = std::move(full);
      continue;
    }
    bool improved = false;
    for (double length = 0.5; length > 1e-10; length /= 2)
      if (evaluate(K, M, P, at(length), limits, 0).objective >= e.objective - 1e-12) {
        beta = at(length);
        e = evaluate(K, M, P, beta, limits, 2);
        improved = true;
        break;
      }
    if (!improved) { fit.converged = true; break; }
  }
  fit.beta = beta;
  fit.gradient = e.gradient;
  fit.log_likelihood = e.log_likelihood;
  fit.log_posterior = e.objective;
  fit.error.assign(m, std::nan(""));
  calx::linalg::Matrix A(m);
  for (std::size_t p = 0; p < m; p++)
    for (std::size_t q = 0; q < m; q++) A(p, q) = -e.hessian[p][q];
  try {
    const auto inverse = calx::linalg::solve(A, calx::linalg::Matrix::identity(m));
    for (std::size_t p = 0; p < m; p++)
      if (inverse(p, p) > 0) fit.error[p] = std::sqrt(inverse(p, p));
  } catch (const std::runtime_error&) {
  }
  return fit;
}

}  // namespace history
