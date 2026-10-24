// Stocks redistributed by a column-stochastic matrix (theory note §3.3).
//
// When each position passes the shares M_ij of its stock to the others in each period (columns
// summing to 1, so nothing is created or destroyed), the stocks follow k(t + 1) = M k(t). The
// total is conserved; for an irreducible, aperiodic M the stocks tend to the stationary vector π
// with M π = π and Σ π = Σ k(0) (Perron and Frobenius), at the rate of M's second eigenvalue. A
// drift of stocks under a fixed rule says nothing about who exploits whom: the user's paper on
// loop deficits makes this point with M = [[0.9, 0.6], [0.1, 0.4]] (its Counterexample 8.2),
// whose stocks drift from (1, 1) to (12/7, 2/7) with no loop and no deficit at all.

#pragma once

#include <calx/linalg/dense.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace economics {

inline void require_column_stochastic(const calx::linalg::Matrix& M, double tolerance = 1e-12) {
  for (std::size_t j = 0; j < M.n; j++) {
    double s = 0;
    for (std::size_t i = 0; i < M.n; i++) {
      if (M(i, j) < 0) throw std::invalid_argument("economics: a stochastic matrix has no negative share");
      s += M(i, j);
    }
    if (std::fabs(s - 1) > tolerance) throw std::invalid_argument("economics: each column of a stochastic matrix sums to 1");
  }
}

// k(0), k(1), …, k(steps).
inline std::vector<std::vector<double>> stock_path(const calx::linalg::Matrix& M, std::vector<double> k, std::size_t steps) {
  require_column_stochastic(M);
  if (k.size() != M.n) throw std::invalid_argument("economics::stock_path: one stock per position");
  std::vector<std::vector<double>> path{k};
  for (std::size_t t = 0; t < steps; t++) {
    std::vector<double> next(M.n, 0.0);
    for (std::size_t i = 0; i < M.n; i++)
      for (std::size_t j = 0; j < M.n; j++) next[i] += M(i, j) * k[j];
    k = next;
    path.push_back(k);
  }
  return path;
}

// The stationary stocks with the given total: (M − I) π = 0 with one equation replaced by Σ π.
// Unique when M is irreducible; calx's solve throws on a singular system otherwise.
inline std::vector<double> stationary(const calx::linalg::Matrix& M, double total) {
  require_column_stochastic(M);
  const std::size_t n = M.n;
  calx::linalg::Matrix A(n), B(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) A(i, j) = M(i, j) - (i == j ? 1.0 : 0.0);
  for (std::size_t j = 0; j < n; j++) A(n - 1, j) = 1;
  B(n - 1, 0) = total;
  const auto X = calx::linalg::solve(A, B);
  std::vector<double> pi(n);
  for (std::size_t i = 0; i < n; i++) pi[i] = X(i, 0);
  return pi;
}

}  // namespace economics
