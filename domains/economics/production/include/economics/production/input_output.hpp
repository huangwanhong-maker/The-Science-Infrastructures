// Input–output analysis: production as a circular flow (theory note §4; Leontief 1936, 1986;
// Miller and Blair 2009, Input–Output Analysis, 2nd ed., ch. 2).
//
// An input–output table of n sectors records z_ij, the deliveries of sector i to sector j in a
// period, with each sector's final demand y_i, value added v_j and gross output x:
//
//   x_i = Σ_j z_ij + y_i   (where each sector's output goes),
//   x_j = Σ_i z_ij + v_j   (what each sector's output costs).
//
// The technical coefficients a_ij = z_ij / x_j make the quantity model x = A x + y, solved by the
// Leontief inverse L = (I − A)⁻¹: x = L y, and the total requirement of sector i per unit of final
// demand for j is L_ij (Miller and Blair, eq. 2.11). Its column sums are the output multipliers.
// The price model is the dual, p = p A + c with c_j = v_j / x_j the value added per unit of
// output, so p = c L: in a table in money every price is 1 (eq. 2.39), and a change in value
// added propagates through L as cost-push. The economy is productive when L ≥ 0 (for a
// non-negative A, the Hawkins–Simon condition: every leading principal minor of I − A positive).
//
// The deliveries z_ij are a flow on the network of sectors, and its Hodge decomposition (physis,
// through economics::circulation) separates what is driven from sector to sector from what
// circulates round the cycles of production, the circular flow Quesnay drew and Leontief
// measured.

#pragma once

#include <economics/circulation/flows.hpp>
#include <economics/transfer/network.hpp>
#include <calx/linalg/dense.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace economics {

struct InputOutputTable {
  std::vector<std::string> sectors;
  calx::linalg::Matrix z;                // deliveries from row to column
  std::vector<double> final_demand;      // y
  std::vector<double> value_added;       // v

  std::size_t n() const { return sectors.size(); }
  // Gross output by rows: x_i = Σ_j z_ij + y_i.
  std::vector<double> output() const {
    std::vector<double> x(n(), 0.0);
    for (std::size_t i = 0; i < n(); i++) {
      for (std::size_t j = 0; j < n(); j++) x[i] += z(i, j);
      x[i] += final_demand.at(i);
    }
    return x;
  }
  // The largest |row total − column total| over the sectors: zero for a balanced table.
  double imbalance() const {
    const auto x = output();
    double worst = 0;
    for (std::size_t j = 0; j < n(); j++) {
      double col = value_added.at(j);
      for (std::size_t i = 0; i < n(); i++) col += z(i, j);
      worst = std::max(worst, std::fabs(col - x[j]));
    }
    return worst;
  }
};

struct LeontiefModel {
  calx::linalg::Matrix A, L;
  std::vector<double> output, multipliers, value_added_coefficients, prices;
  std::vector<double> backward, forward;  // Rasmussen's normalised linkages: column and row sums of L over their mean
  bool productive = false;                // L ≥ 0
};

inline LeontiefModel leontief(const InputOutputTable& T) {
  using calx::linalg::Matrix;
  const std::size_t n = T.n();
  if (T.z.n != n || T.final_demand.size() != n || T.value_added.size() != n) throw std::invalid_argument("economics::leontief: one row, column and entry per sector");
  LeontiefModel M;
  M.output = T.output();
  M.A = Matrix(n);
  for (std::size_t j = 0; j < n; j++) {
    if (!(M.output[j] > 0)) throw std::invalid_argument("economics::leontief: every sector has positive output");
    for (std::size_t i = 0; i < n; i++) M.A(i, j) = T.z(i, j) / M.output[j];
  }
  M.L = calx::linalg::solve(Matrix::identity(n) - M.A, Matrix::identity(n));
  M.productive = true;
  for (double v : M.L.a) M.productive = M.productive && v >= -1e-12;
  M.multipliers.assign(n, 0.0);
  std::vector<double> rows(n, 0.0);
  double mean = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      M.multipliers[j] += M.L(i, j);
      rows[i] += M.L(i, j);
      mean += M.L(i, j);
    }
  mean /= static_cast<double>(n * n);
  for (std::size_t k = 0; k < n; k++) {
    M.backward.push_back(M.multipliers[k] / static_cast<double>(n) / mean);
    M.forward.push_back(rows[k] / static_cast<double>(n) / mean);
  }
  M.value_added_coefficients.resize(n);
  for (std::size_t j = 0; j < n; j++) M.value_added_coefficients[j] = T.value_added[j] / M.output[j];
  M.prices.assign(n, 0.0);
  for (std::size_t j = 0; j < n; j++)
    for (std::size_t i = 0; i < n; i++) M.prices[j] += M.value_added_coefficients[i] * M.L(i, j);
  return M;
}

// The output that meets a final demand: x = L y.
inline std::vector<double> output_for(const LeontiefModel& M, const std::vector<double>& y) {
  std::vector<double> x(M.L.n, 0.0);
  for (std::size_t i = 0; i < M.L.n; i++)
    for (std::size_t j = 0; j < M.L.n; j++) x[i] += M.L(i, j) * y.at(j);
  return x;
}

// The prices under value-added coefficients c: p = c L (cost-push).
inline std::vector<double> prices_for(const LeontiefModel& M, const std::vector<double>& c) {
  std::vector<double> p(M.L.n, 0.0);
  for (std::size_t j = 0; j < M.L.n; j++)
    for (std::size_t i = 0; i < M.L.n; i++) p[j] += c.at(i) * M.L(i, j);
  return p;
}

// The deliveries as a transfer network of sectors, one leg for each positive z_ij (i ≠ j), its
// flow z_ij and its transport 1 (a table in money registers a delivery at its value).
inline TransferNetwork delivery_network(const InputOutputTable& T) {
  TransferNetwork N;
  N.positions = T.sectors;
  for (std::size_t i = 0; i < T.n(); i++)
    for (std::size_t j = 0; j < T.n(); j++)
      if (i != j && T.z(i, j) > 0) N.legs.push_back({i, j, 1.0, T.z(i, j), "deliveries"});
  return N;
}

}  // namespace economics
