// Checks of transfer networks (network.hpp).
//
// The worked numbers of the user's paper on loop deficits (Huang 2026, The Geometry of Value):
// the triangle of transports 1.2, 0.5, 1.0 has holonomy 0.60 and deficit 0.511, unchanged when
// re-described by g = (1, 2.5, 0.4), which sends its legs to 3.0, 0.08, 2.5; the network of five
// positions and six legs has chord holonomies 0.975 and 1.0725 in the tree gauge. Then: deficits
// change sign with the direction of a loop and add under concatenation; every chord holonomy is
// unchanged by a random re-description, while the log prices shift by log g; a network built
// from prices has no deficit and gives its prices back; on a directed cycle with steady flows the
// retentions multiply to the holonomy.

#include <economics/transfer/network.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x7472616eu, 0, 0);  // "tran"

  // The triangle.
  TransferNetwork T;
  T.positions = {"a", "b", "c"};
  T.legs = {{0, 1, 1.2, 10, "money"}, {1, 2, 0.5, 12, "money"}, {2, 0, 1.0, 6, "money"}};
  const Path round{{0, 1}, {1, 1}, {2, 1}};
  s.check("the paper's triangle: holonomy 0.60", holonomy(T, round) - 0.60, 1e-12);
  s.check("... deficit −ln 0.6 = 0.511", deficit(T, round) - (-std::log(0.6)), 1e-12);
  const TransferNetwork G = regauged(T, {1.0, 2.5, 0.4});
  s.require("... re-described by g = (1, 2.5, 0.4): legs 3.0, 0.08, 2.5",
            std::fabs(G.legs[0].ratio - 3.0) < 1e-12 && std::fabs(G.legs[1].ratio - 0.08) < 1e-12 && std::fabs(G.legs[2].ratio - 2.5) < 1e-12);
  s.check("... and its holonomy unchanged", holonomy(G, round) - 0.60, 1e-12);
  const Path reversed{{2, -1}, {1, -1}, {0, -1}};
  s.check("the deficit of the reversed loop is the negative", deficit(T, reversed) + deficit(T, round), 1e-12);
  Path twice = round;
  twice.insert(twice.end(), round.begin(), round.end());
  s.check("deficits add under concatenation", deficit(T, twice) - 2 * deficit(T, round), 1e-12);
  // Retention on the directed triangle with steady flows 10, 12, 6.
  double product = 1;
  for (double r : retention(T, round)) product *= r;
  s.check("retentions on a directed cycle multiply to the holonomy", product - 0.60, 1e-12);

  // The five positions and six legs: tree r→1, r→2, 1→3, 2→4; chords 1→2, 3→4.
  TransferNetwork F;
  F.positions = {"r", "1", "2", "3", "4"};
  F.legs = {{0, 1, 1.3}, {0, 2, 0.8}, {1, 3, 0.9}, {2, 4, 1.5}, {1, 2, 0.6}, {3, 4, 1.1}};
  const Circuits C = circuits(F);
  s.require("the paper's five positions: two chords, 1→2 and 3→4", C.chords == std::vector<std::size_t>{4, 5});
  s.check("... chord 1→2 carries its cycle's holonomy 0.975", C.holonomy[0] - 0.975, 1e-12);
  s.check("... chord 3→4 carries 1.0725", C.holonomy[1] - 1.0725, 1e-12);
  // Invariance under a random re-description.
  std::vector<double> g(F.positions.size());
  for (double& x : g) x = std::exp(z.uniform() * 4 - 2);
  const Circuits Cg = circuits(regauged(F, g));
  double worst = 0;
  for (std::size_t i = 0; i < C.holonomy.size(); i++) worst = std::max(worst, std::fabs(Cg.holonomy[i] - C.holonomy[i]));
  s.check("chord holonomies are unchanged by a random re-description", worst, 1e-12);
  double shift = 0;
  for (std::size_t v = 0; v < g.size(); v++)
    shift = std::max(shift, std::fabs((Cg.log_price[v] - C.log_price[v]) - (std::log(g[v]) - std::log(g[0]))));
  s.check("... while the log prices shift by log g", shift, 1e-12);

  // A network built from prices has no deficit, and gives its prices back.
  TransferNetwork P;
  const std::size_t n = 12;
  std::vector<double> phi(n);
  for (std::size_t v = 0; v < n; v++) {
    P.positions.push_back("p" + std::to_string(v));
    phi[v] = z.uniform() * 3 - 1.5;
  }
  for (std::size_t k = 0; k < 30; k++) {
    const std::size_t u = z.below(n), v = z.below(n);
    if (u != v) P.legs.push_back({u, v, std::exp(phi[v] - phi[u])});
  }
  for (std::size_t v = 1; v < n; v++) P.legs.push_back({v - 1, v, std::exp(phi[v] - phi[v - 1])});  // connected
  const Circuits Cp = circuits(P);
  s.check("a network built from prices has no deficit round any loop", Cp.defect, 1e-12);
  double back = 0;
  for (std::size_t v = 0; v < n; v++) back = std::max(back, std::fabs((Cp.log_price[v] - Cp.log_price[0]) - (phi[v] - phi[0])));
  s.check("... and its prices are recovered up to a constant", back, 1e-12);
  return s.done();
}
