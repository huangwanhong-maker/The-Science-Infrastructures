// Checks of trade (partners.hpp) on a world of four economies whose numbers are worked by hand, and
// on a random one: the partner concentration of a row is its Herfindahl index, the intensities
// average to one for each partner when weighted by exports, the two-way share of a pair is
// Grubel–Lloyd's, and the world's reciprocal share is circulation's reciprocal exchange.

#include <economics/circulation/flows.hpp>
#include <economics/trade/partners.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;

  // x_ij: A sends 3 to B and 1 to C; B sends 1 to A and 1 to D; C sends 2 to A; D sends 4 to B.
  const auto x = trade::Rectangular::from_rows({{0, 3, 1, 0}, {1, 0, 0, 1}, {2, 0, 0, 0}, {0, 4, 0, 0}});
  const auto P = trade::partner_concentration(x);
  s.check("partners: A's exports, 3/4 and 1/4, have H = 0.625", P.hhi[0] - 0.625, 1e-15);
  s.check("... and 1.6 equal partners", P.equivalent[0] - 1.6, 1e-14);
  s.require("... its largest partner is B, with three quarters", P.first[0] == 1 && std::fabs(P.first_share[0] - 0.75) < 1e-15);
  s.check("partners: C sells to A alone, H = 1", P.hhi[2] - 1, 1e-15);

  const auto I = trade::intensity(x);
  // W = 12; A's share to B is 3/4, B's share of world imports 7/12: I = (3/4)/(7/12) = 9/7.
  s.check("intensity: I_AB = (3/4) / (7/12) = 9/7", I(0, 1) - 9.0 / 7, 1e-14);
  // B sends half its exports to D, which takes 1/12 of the world's imports: I_BD = (1/2)/(1/12) = 6.
  s.check("intensity: I_BD = (1/2) / (1/12) = 6", I(1, 3) - 6, 1e-14);

  const auto G = trade::two_way(x);
  s.check("two-way: A and B trade 3 and 1, GL = 1 − 2/4 = 0.5", G(0, 1) - 0.5, 1e-15);
  s.require("two-way: symmetric, and NaN for a pair that does not trade", G(1, 0) == G(0, 1) && std::isnan(G(2, 3)));

  // A random world: the identities hold for any matrix.
  calx::random::Stream rng(20261025, 0x74726164u, 0, 0);  // "trad"
  const std::size_t n = 12;
  trade::Rectangular y(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (i != j && rng.uniform() < 0.7) y(i, j) = std::exp(3 * rng.uniform());
  const auto X = trade::exports(y), M = trade::imports(y);
  double W = 0;
  for (double v : X) W += v;
  const auto Iy = trade::intensity(y);
  double worst = 0;
  for (std::size_t j = 0; j < n; j++) {
    if (!(M[j] > 0)) continue;
    double mean = 0;
    for (std::size_t i = 0; i < n; i++)
      if (X[i] > 0) mean += X[i] / W * Iy(i, j);
    worst = std::max(worst, std::fabs(mean - 1));
  }
  s.check("intensity: for every partner, the exporters' intensities weighted by their exports average to 1", worst, 1e-12);

  TransferNetwork N;
  for (std::size_t i = 0; i < n; i++) N.positions.push_back(std::to_string(i));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++)
      if (y(i, j) > 0) N.legs.push_back({i, j, 1.0, y(i, j), "goods"});
  const auto T = net_transfer(N);
  s.check("two-way: the world's reciprocal share is circulation's reciprocal exchange over the whole", trade::reciprocal_share(y) - T.reciprocal / T.total, 1e-12);

  return s.done();
}
