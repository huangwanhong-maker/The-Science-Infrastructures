// Checks of input–output analysis (input_output.hpp).
//
// Miller and Blair's two-sector example (2009, §2.3): deliveries [[150, 500], [200, 100]], final
// demand (350, 1700), outputs (1000, 2000); A = [[0.15, 0.25], [0.20, 0.05]] and
// L = (I − A)⁻¹ = (1/0.7575) [[0.95, 0.25], [0.20, 0.85]]; L y = x; every price of a table in
// money is 1; the output multipliers are L's column sums; a rise in one sector's value added
// raises every price by its row of L (cost-push). A table that cannot meet its demand is not
// productive. The deliveries of a three-sector table with a cycle split into driven flow and
// circulation as least squares by hand says, their net inflows the sectors' inputs less their
// deliveries.

#include <economics/production/input_output.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;
  InputOutputTable T;
  T.sectors = {"agriculture", "manufacturing"};
  T.z = calx::linalg::Matrix(2);
  T.z(0, 0) = 150, T.z(0, 1) = 500, T.z(1, 0) = 200, T.z(1, 1) = 100;
  T.final_demand = {350, 1700};
  T.value_added = {650, 1400};
  s.check("Miller and Blair's table balances, rows and columns", T.imbalance(), 1e-9);
  const auto M = leontief(T);
  const double A[2][2] = {{0.15, 0.25}, {0.20, 0.05}};
  const double d = 0.7575, L[2][2] = {{0.95 / d, 0.25 / d}, {0.20 / d, 0.85 / d}};
  double ea = 0, el = 0;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      ea = std::max(ea, std::fabs(M.A(i, j) - A[i][j]));
      el = std::max(el, std::fabs(M.L(i, j) - L[i][j]));
    }
  s.check("... its technical coefficients A", ea, 1e-12);
  s.check("... its Leontief inverse (1/0.7575) [[0.95, 0.25], [0.20, 0.85]]", el, 1e-12);
  const auto x = output_for(M, T.final_demand);
  s.check("... L y is the output (1000, 2000)", std::fabs(x[0] - 1000) + std::fabs(x[1] - 2000), 1e-9);
  s.check("... every price of a table in money is 1", std::fabs(M.prices[0] - 1) + std::fabs(M.prices[1] - 1), 1e-12);
  s.check("... the output multipliers are L's column sums", std::fabs(M.multipliers[0] - (L[0][0] + L[1][0])) + std::fabs(M.multipliers[1] - (L[0][1] + L[1][1])), 1e-12);
  s.require("... and it is productive", M.productive);
  // Cost-push: raise agriculture's value added per unit by 0.1.
  auto c = M.value_added_coefficients;
  c[0] += 0.1;
  const auto p = prices_for(M, c);
  s.check("a rise of 0.1 in agriculture's value added raises each price by 0.1 L₀ⱼ", std::fabs(p[0] - (1 + 0.1 * L[0][0])) + std::fabs(p[1] - (1 + 0.1 * L[0][1])), 1e-12);

  // A table that cannot meet its demand.
  InputOutputTable U;
  U.sectors = {"a", "b"};
  U.z = calx::linalg::Matrix(2);
  U.z(0, 0) = 50, U.z(0, 1) = 60, U.z(1, 0) = 70, U.z(1, 1) = 40;
  U.final_demand = {-80, -80};
  U.value_added = {-90, -70};
  s.require("a table whose sectors use more than they make is not productive", !leontief(U).productive);

  // The circulation of deliveries in three sectors with a cycle.
  InputOutputTable C;
  C.sectors = {"s", "t", "u"};
  C.z = calx::linalg::Matrix(3);
  C.z(0, 1) = 40, C.z(1, 2) = 40, C.z(2, 0) = 40, C.z(0, 2) = 25;
  C.final_demand = {100, 80, 60};
  C.value_added = {0, 0, 0};
  const auto N = delivery_network(C);
  const auto flow = circulation(N);
  s.check("the deliveries' three shares add to 1", flow.driven + flow.local + flow.global - 1, 1e-12);
  // Least squares by hand: α_t − α_s = α_u − α_t = 5, so the driven part is (5, 5, −10, 10) on
  // s → t, t → u, u → s, s → u, of ‖·‖² = 250 out of 3 · 40² + 25² = 5425; the rest circulates.
  s.check("... the driven share is 250 / 5425 (least squares by hand)", flow.driven - 250.0 / 5425.0, 1e-12);
  s.check("... the attraction rises by 5 from s to t and from t to u",
          std::fabs(flow.attraction[1] - flow.attraction[0] - 5) + std::fabs(flow.attraction[2] - flow.attraction[1] - 5), 1e-12);
  double inflow = 0;
  for (std::size_t j = 0; j < 3; j++) {
    double in = 0, out = 0;
    for (std::size_t i = 0; i < 3; i++)
      if (i != j) { in += C.z(i, j); out += C.z(j, i); }
    inflow = std::max(inflow, std::fabs(flow.net_inflow[j] - (in - out)));
  }
  s.check("... each sector's net inflow is its inputs less its deliveries", inflow, 1e-12);
  return s.done();
}
