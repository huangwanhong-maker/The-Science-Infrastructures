// Checks of circulation (flows.hpp) and stocks (stocks.hpp).
//
// Flows: equal flows round an unfilled directed triangle are wholly global circulation, with every
// position balanced, and wholly local once the triangle is a face; a flow d₀ α is wholly driven,
// with attraction α; a round trip A → B → A is circulation round the hole it makes; the shares add
// to 1; balances by kind are inflow − outflow; reciprocal exchange and net transfer add up to the
// flow, and the net network keeps every net inflow. Pairing: its three terms add up to Σ F log U, and
// a driven flow meets no deficit. Conductances: with the gross volumes, the driven flow is Ohm's law
// for the attraction and carries every net inflow; with any conductances a driven flow meets no
// deficit. Stocks: the paper's M = [[0.9, 0.6], [0.1, 0.4]] from (1, 1)
// follows k(t) = π − (5/7)(0.3)ᵗ (1, −1) with π = (12/7, 2/7), its total conserved.

#include <economics/circulation/flows.hpp>
#include <economics/circulation/stocks.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <set>
#include <utility>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x63697263u, 0, 0);  // "circ"

  TransferNetwork T;
  T.positions = {"a", "b", "c"};
  T.legs = {{0, 1, 1.2, 5, "money"}, {1, 2, 0.5, 5, "money"}, {2, 0, 1.0, 5, "money"}};
  const auto open = circulation(T);
  s.check("equal flows round an unfilled triangle: wholly global circulation", open.global - 1, 1e-12);
  double unbalanced = 0;
  for (double x : open.net_inflow) unbalanced = std::max(unbalanced, std::fabs(x));
  s.check("... every position balanced", unbalanced, 1e-12);
  T.faces = {{{0, 1}, {1, 1}, {2, 1}}};
  s.check("... and wholly local once the triangle is a face", circulation(T).local - 1, 1e-12);

  // A driven flow.
  TransferNetwork D;
  const std::vector<double> alpha{0.0, 2.0, 5.0, 1.0};
  D.positions = {"p", "q", "r", "s"};
  for (std::size_t u = 0; u < 4; u++)
    for (std::size_t v = u + 1; v < 4; v++) D.legs.push_back({u, v, 1.0, alpha[v] - alpha[u], "goods"});
  const auto driven = circulation(D);
  s.check("a flow d₀ α is wholly driven", driven.driven - 1, 1e-12);
  double attr = 0;
  for (std::size_t v = 0; v < 4; v++) attr = std::max(attr, std::fabs((driven.attraction[v] - driven.attraction[0]) - (alpha[v] - alpha[0])));
  s.check("... its attraction is α up to a constant", attr, 1e-12);

  // A round trip.
  TransferNetwork R;
  R.positions = {"A", "B"};
  R.legs = {{0, 1, 1.0, 3, "money"}, {1, 0, 1.0, 3, "money"}};
  s.check("a round trip A → B → A is circulation round the hole it makes", circulation(R).global - 1, 1e-12);

  // A random flow: shares add to 1; balances by kind.
  TransferNetwork X;
  for (std::size_t v = 0; v < 9; v++) X.positions.push_back("x" + std::to_string(v));
  for (std::size_t k = 0; k < 25; k++) {
    const std::size_t u = z.below(9), v = z.below(9);
    if (u != v) X.legs.push_back({u, v, std::exp(z.uniform() - 0.5), z.uniform() * 10, k % 2 ? "money" : "goods"});
  }
  const auto rc = circulation(X);
  s.check("the three shares of a random flow add to 1", rc.driven + rc.local + rc.global - 1, 1e-10);
  const auto B = balances(X);
  double sum = 0;
  for (const auto& [kind, b] : B)
    for (std::size_t v = 0; v < 9; v++) sum += b[v];
  s.check("balances by kind sum to zero over the positions (nothing is created by moving it)", sum, 1e-9);
  double agree = 0;
  for (std::size_t v = 0; v < 9; v++) agree = std::max(agree, std::fabs(B.at("money")[v] + B.at("goods")[v] - rc.net_inflow[v]));
  s.check("... and add up to the net inflow", agree, 1e-9);

  // Reciprocal exchange and net transfer: reciprocal + net = total, the net network has one leg
  // per unbalanced pair, and its net inflows are those of the flows.
  const auto NT = net_transfer(X);
  s.check("reciprocal exchange and net transfer add up to the whole flow", NT.reciprocal + NT.net - NT.total, 1e-9);
  const auto nc = circulation(NT.network);
  double same = 0;
  for (std::size_t v = 0; v < 9; v++) same = std::max(same, std::fabs(nc.net_inflow[v] - rc.net_inflow[v]));
  s.check("... and the net network has the same net inflow at every position", same, 1e-9);
  s.require("... with one leg per pair, along the net, its gross volume no less", [&] {
    std::set<std::pair<std::size_t, std::size_t>> seen;
    for (std::size_t e = 0; e < NT.network.legs.size(); e++) {
      const auto& l = NT.network.legs[e];
      if (!seen.insert(std::minmax(l.from, l.to)).second || !(l.flow >= 0) || !(NT.gross[e] >= l.flow)) return false;
    }
    return NT.gross.size() == NT.network.legs.size();
  }());

  // Pairing.
  const auto P = pairing(X);
  s.check("the pairing's three terms add up to ⟨F, log U⟩", P.prices + P.local + P.global - P.total, 1e-10);
  TransferNetwork Y = D;
  Y.legs[0].ratio = 1.7;  // a deficit round the triangles through p and q
  Y.faces = {{{0, 1}, {3, 1}, {1, -1}}};  // p → q → r → p
  const auto Pd = pairing(Y);
  s.check("a driven flow meets no deficit: its local and global pairings vanish", std::fabs(Pd.local) + std::fabs(Pd.global), 1e-12);

  // Conductances. The net network's gross volumes as conductances: the driven net flow on each
  // link is its gross volume times the difference in attraction, and it carries every net inflow.
  const auto gc = circulation(NT.network, NT.gross);
  double ohm = 0;
  for (std::size_t e = 0; e < NT.network.legs.size(); e++) {
    const auto& l = NT.network.legs[e];
    ohm = std::max(ohm, std::fabs(gc.parts.driven[e] - NT.gross[e] * (gc.attraction[l.to] - gc.attraction[l.from])));
  }
  const Cochain gd = physis::discretisation::d0_transpose(NT.network.complex(), gc.parts.driven);
  for (std::size_t v = 0; v < 9; v++) ohm = std::max(ohm, std::fabs(gd[v] - rc.net_inflow[v]));
  s.check("with the gross volumes as conductances, the driven flow is gross × the difference in attraction, carrying every net inflow", ohm, 1e-9);
  s.check("... the energy shares add to 1", gc.driven + gc.local + gc.global - 1, 1e-10);
  // A flow c d₀ α with conductances c meets no deficit either.
  TransferNetwork Yc = Y;
  Cochain c;
  for (auto& l : Yc.legs) {
    c.push_back(std::exp(z.uniform() - 0.5));
    l.flow = c.back() * (alpha[l.to] - alpha[l.from]);
  }
  const auto Pc = pairing(Yc, c);
  s.check("... and with any conductances, a driven flow c d₀ α meets no deficit, its pairings adding up",
          std::fabs(Pc.local) + std::fabs(Pc.global) + std::fabs(Pc.prices + Pc.local + Pc.global - Pc.total), 1e-12);

  // Stocks.
  calx::linalg::Matrix M(2);
  M(0, 0) = 0.9, M(0, 1) = 0.6, M(1, 0) = 0.1, M(1, 1) = 0.4;
  const auto path = stock_path(M, {1.0, 1.0}, 12);
  double worst = 0, total = 0;
  for (std::size_t t = 0; t < path.size(); t++) {
    const double c = (5.0 / 7.0) * std::pow(0.3, static_cast<double>(t));
    worst = std::max({worst, std::fabs(path[t][0] - (12.0 / 7.0 - c)), std::fabs(path[t][1] - (2.0 / 7.0 + c))});
    total = std::max(total, std::fabs(path[t][0] + path[t][1] - 2));
  }
  s.check("the paper's stock drift follows π − (5/7)(0.3)ᵗ (1, −1)", worst, 1e-12);
  s.check("... with its total conserved", total, 1e-12);
  const auto pi = stationary(M, 2.0);
  s.check("... towards the stationary stocks (12/7, 2/7)", std::fabs(pi[0] - 12.0 / 7.0) + std::fabs(pi[1] - 2.0 / 7.0), 1e-12);
  return s.done();
}
