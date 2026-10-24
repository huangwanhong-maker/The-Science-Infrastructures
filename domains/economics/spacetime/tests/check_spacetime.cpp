// Checks of the spacetime complex of a record of transactions (complex.hpp).
//
// A market built from a price system (log prices c_x + r_x t, every exchange their ratio) is flat:
// every lens and triangle has zero curvature, no loop has a deficit, and the Hodge decomposition
// of its transports is wholly exact. A lens of two agents exchanging a → b at p and q has the
// closed form (r_a − r_b)(t_q − t_p) + log U_q − log U_p; a triangle of exchanges 2, 3 and 5 at
// one event has log(2 · 3 / 5). Changing one exchange's rate by e^ε changes the curvature by ±ε on
// exactly the 2-cells through it. Reach: viable along a chain of exchanges, unknown where one link
// is unknown, lost where one is not viable; and making one exchange less viable never makes any
// presence more reachable (the user's Proposition 5.1).

#include <economics/circulation/flows.hpp>
#include <economics/spacetime/complex.hpp>
#include <physis/discretisation/hodge.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <set>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x73706163u, 0, 0);  // "spac"

  // A flat market: three agents, four dates, every pair exchanging at most dates.
  Transactions T;
  T.agents = {"a", "b", "c"};
  T.holding_rate = {0.02, 0.05, -0.01};
  const std::vector<double> c0{0.0, 0.3, -0.4};
  for (double t : {0.0, 0.5, 1.2, 2.0}) T.events.push_back({t, {0, 1, 2}, ""});
  auto log_price = [&](std::size_t x, std::size_t v) { return c0[x] + T.holding_rate[x] * T.events[v].time; };
  for (std::size_t v = 0; v < 4; v++)
    for (std::size_t x = 0; x < 3; x++)
      for (std::size_t y = x + 1; y < 3; y++)
        if ((v + x + y) % 4 != 3) T.exchanges.push_back({v, (v % 2) ? y : x, (v % 2) ? x : y, 0.0, 1.0});
  for (auto& e : T.exchanges) e.ratio = std::exp(log_price(e.to, e.event) - log_price(e.from, e.event));
  const Spacetime S = spacetime(T);
  double kmax = 0;
  for (double k : S.curvature) kmax = std::max(kmax, std::fabs(k));
  s.require("a flat market has lenses and triangles", S.lenses > 0 && S.curvature.size() > S.lenses, std::to_string(S.lenses) + " lenses, " + std::to_string(S.curvature.size() - S.lenses) + " triangles");
  s.check("... every one of zero curvature", kmax, 1e-12);
  s.check("... no loop with a deficit", circuits(S.network).defect, 1e-12);
  const auto H = physis::discretisation::hodge_decomposition(S.network.complex(), S.network.transport());
  s.check("... its transports wholly exact (a price system)", H.coexact_norm2 + H.harmonic_norm2, 1e-20);

  // The lens closed form.
  {
    Transactions L;
    L.agents = {"a", "b"};
    L.holding_rate = {0.04, 0.01};
    L.events = {{0.0, {0, 1}, "p"}, {0.75, {0, 1}, "q"}};
    L.exchanges = {{0, 0, 1, 1.30, 0}, {1, 0, 1, 1.45, 0}};
    const Spacetime W = spacetime(L);
    const double expected = (0.04 - 0.01) * 0.75 + std::log(1.45) - std::log(1.30);
    s.check("a lens a → b at p and q: (r_a − r_b) Δt + log U_q − log U_p", W.curvature[0] - expected, 1e-12);
  }
  // Rates for each event: a lens over two months at different rates sums them.
  {
    Transactions L;
    L.agents = {"a", "b"};
    L.events = {{0.0, {0, 1}, "p"}, {0.5, {0, 1}, "m"}, {1.0, {0, 1}, "q"}};
    L.holding_rates = {{0.04, 0.06, 0.0}, {0.01, 0.02, 0.0}};
    L.exchanges = {{0, 0, 1, 1.30, 0}, {2, 0, 1, 1.45, 0}};
    const Spacetime W = spacetime(L);
    const double expected = (0.04 - 0.01) * 0.5 + (0.06 - 0.02) * 0.5 + std::log(1.45) - std::log(1.30);
    s.check("rates for each event: a lens over two periods sums their differences", W.curvature[0] - expected, 1e-12);
  }
  // A triangle at one event.
  {
    Transactions R;
    R.agents = {"a", "b", "c"};
    R.events = {{0.0, {0, 1, 2}, ""}};
    R.exchanges = {{0, 0, 1, 2.0, 0}, {0, 1, 2, 3.0, 0}, {0, 0, 2, 5.0, 0}};
    const Spacetime W = spacetime(R);
    s.require("three exchanges at one event make one triangle", W.lenses == 0 && W.curvature.size() == 1);
    s.check("... whose curvature is the triangular arbitrage log(2 · 3 / 5)", W.curvature[0] - std::log(1.2), 1e-12);
  }
  // Stokes: a change of e^ε in one exchange moves the curvature of exactly the 2-cells through it.
  {
    Transactions P = T;
    const std::size_t k = 4;
    const double eps = 0.03;
    P.exchanges[k].ratio *= std::exp(eps);
    const Spacetime W = spacetime(P);
    // The rung of exchange k, and the 2-cells whose words contain it.
    std::size_t rung = 0;
    for (std::size_t e = W.refined.temporal; e < W.refined.cw.edges.size(); e++) {
      const auto [t, m] = W.refined.rung_of_edge[e - W.refined.temporal];
      if (W.exchange_of[t][m] == k) rung = e;
    }
    std::set<std::size_t> through;
    for (std::size_t f = 0; f < W.refined.cw.faces.size(); f++)
      for (const auto& step : W.refined.cw.faces[f].boundary)
        if (step.edge == rung) through.insert(f);
    bool exact = true;
    for (std::size_t f = 0; f < W.curvature.size(); f++) {
      const double d = std::fabs(W.curvature[f] - S.curvature[f]);
      exact = exact && (through.count(f) ? std::fabs(d - eps) < 1e-12 : d < 1e-12);
    }
    s.require("changing one exchange by e^ε moves by ±ε the curvature of exactly the 2-cells through it", exact && !through.empty(),
              std::to_string(through.size()) + " 2-cells");
  }
  // Reach.
  {
    Transactions R;
    R.agents = {"a", "b", "c"};
    R.events = {{0.0, {0}, "e0"}, {1.0, {0, 1}, "e1"}, {2.0, {1, 2}, "e2"}, {3.0, {2}, "e3"}};
    R.exchanges = {{1, 0, 1, 1.0, 1}, {2, 1, 2, 1.0, 1}};
    const Spacetime W = spacetime(R);
    const std::size_t c3 = W.refined.presence_of.at({2, 3});
    using V = Viability;
    s.require("value from a at e0 reaches c at e3 along a viable chain", reach(W, R, 0, 0)[c3] == V::yes);
    s.require("... is unknown there if one link is unknown", reach(W, R, 0, 0, {V::yes, V::unknown})[c3] == V::unknown);
    s.require("... and lost if one link is not viable", reach(W, R, 0, 0, {V::no, V::unknown})[c3] == V::no);
    s.require("... and never reaches back in time", reach(W, R, 2, 2)[W.refined.presence_of.at({0, 1})] == V::no);
    // Monotonicity on the flat market, with random viabilities.
    bool monotone = true;
    for (int trial = 0; trial < 40; trial++) {
      std::vector<V> v(T.exchanges.size());
      for (auto& x : v) x = static_cast<V>(z.below(3));
      const auto before = reach(S, T, 0, 0, v);
      const std::size_t k = z.below(v.size());
      if (v[k] == V::no) continue;
      v[k] = static_cast<V>(static_cast<int>(v[k]) - 1);
      const auto after = reach(S, T, 0, 0, v);
      for (std::size_t p = 0; p < before.size(); p++) monotone = monotone && static_cast<int>(after[p]) <= static_cast<int>(before[p]);
    }
    s.require("making one exchange less viable never makes any presence more reachable", monotone);
  }
  return s.done();
}
