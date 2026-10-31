// Checks of the flow module (theory/street-oceanography.md §6):
//
// - Theorem 6.1: for random stays the person-time integral equals the sum of the stays' parts in the window
//   exactly (to rounding), and Little's gap vanishes when every stay lies inside the window.
// - The cordon budget accumulates and its residual is the miscount.
// - The street day (Representation 6.5): the crowd conserved to rounding on both configurations; opening the
//   carriageway lowers the peak density and does not lengthen the time on the street; red crossings under a
//   heavy demand make queues on an ordinary day.
// - Theorem 6.3 on a district of 4 × 3 blocks: an attraction alone satisfies detailed balance (σ = 0, every
//   affinity 0) and the Hodge decomposition returns it as the exact part with the attraction as potential;
//   circulation round the blocks makes σ > 0 with a coexact part, and adding the attraction leaves every cycle
//   affinity unchanged; circulation round a hole is harmonic.
// - The window of the street day: the crowd still conserved, a red crossing passes nobody, the walkers' speeds
//   within [0, v_f], the signal green half the time; a street of its own (blocks of given lengths, a crossing
//   without a signal) as long as its parts, its crowd conserved, its unsignalled crossing always passing.
// - Currents: the edges' terms of the entropy production sum to Schnakenberg's σ from physis on a driven
//   ring; a balanced edge has none; a current meter of constructed counts finds its flood and ebb hours and a
//   residual of zero; the sign test of handedness against its binomial law.
// - Transit: a district of one intersection with an exit is exponential with mean 1/exit rate; on a grid the
//   survival starts at 1 and falls.

#include <urban/flow/budget.hpp>
#include <urban/flow/currents.hpp>
#include <urban/flow/kinematic.hpp>
#include <urban/flow/movement.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace urban::flow;

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(8, 8, 8, 8);

  // ---- person-time ----
  {
    std::vector<double> a, d;
    for (int i = 0; i < 500; i++) {
      const double t = rng.uniform(-100, 1000);
      a.push_back(t), d.push_back(t + rng.exponential() * 60);
    }
    const auto p = person_time(a, d, 900);
    s.check("person-time: ∫N dt equals the sum of the stays in the window (relative)", (p.integral - p.sum_of_stays) / p.sum_of_stays, 1e-12);
    std::vector<double> a2, d2;
    for (int i = 0; i < 300; i++) {
      const double t = rng.uniform(10, 800);
      a2.push_back(t), d2.push_back(t + rng.uniform(1, 50));
    }
    const auto q = person_time(a2, d2, 900);
    s.check("person-time: Little's L = λW exactly when every stay lies in the window", q.little_gap, 1e-12);
    const auto b = cordon_budget({5, 3, 0, 2}, {1, 2, 4, 1}, 10, 60, 12);
    s.check("cordon: N accumulates and the residual is the miscount (12 counted, 12 accumulated)", b.residual, 1e-12);
    s.check("cordon: person-time Σ N_k Δt", b.person_time - 60 * (10 + 14 + 15 + 11), 1e-9);
  }

  // ---- the street day ----
  {
    StreetDay day;
    day.blocks = 3, day.block_length = 125, day.from = 11, day.to = 14, day.snapshot = 120;
    day.rule_from = 11, day.rule_to = 14;   // the rule over the whole run, so that the two days differ only by it
    day.demand.assign(24, 0);
    for (int h = 10; h < 20; h++) day.demand[static_cast<std::size_t>(h)] = 12;   // people per second: heavy
    const auto ordinary = street_day(day);
    day.rule = true;
    const auto opened = street_day(day);
    s.check("street day: the crowd is conserved on an ordinary day (relative to those entered)", ordinary.budget_residual / std::max(1.0, ordinary.entered.back()), 1e-10);
    s.check("street day: and with the carriageway opened", opened.budget_residual / std::max(1.0, opened.entered.back()), 1e-10);
    s.require("street day: opening the carriageway lowers the peak density", opened.max_density < ordinary.max_density,
              std::to_string(ordinary.max_density) + " → " + std::to_string(opened.max_density));
    s.require("street day: and does not lengthen the time on the street", opened.mean_time <= ordinary.mean_time * 1.001,
              std::to_string(ordinary.mean_time) + " → " + std::to_string(opened.mean_time));
    s.require("street day: red crossings make queues under heavy demand", ordinary.queue_spill > 0, std::to_string(ordinary.queue_spill));
    // The window: four minutes at two seconds, on the ordinary day.
    StreetDay w = day;
    w.rule = false, w.window_from = 12 * 3600, w.window_to = 12 * 3600 + 240, w.window_step = 2;
    const auto o = street_day(w);
    s.check("street day: the crowd is conserved with a window kept", o.budget_residual / std::max(1.0, o.entered.back()), 1e-10);
    s.require("street day: the window has its 121 snapshots", o.window_times.size() == 121 && std::fabs(o.window_times.back() - w.window_to) < 1e-6,
              std::to_string(o.window_times.size()));
    double red_flow = 0, lo = 1e9, hi = -1e9;
    std::size_t greens = 0, reds = 0;
    for (std::size_t k = 0; k < o.window_times.size(); k++) {
      for (std::size_t c = 0; c < o.crossings.size(); c++) {
        std::size_t at = 0;
        while (at < o.interfaces.size() && std::fabs(o.interfaces[at] - o.crossings[c]) > 1e-9) at++;
        if (o.window_green[k][c]) greens++;
        else reds++, red_flow = std::max(red_flow, std::fabs(o.window_flow[k][at]));
      }
      for (double v : o.window_speed[k]) lo = std::min(lo, v), hi = std::max(hi, v);
    }
    s.check("street day: a red crossing passes nobody (people per second at its stop line)", red_flow, 1e-6);
    s.require("street day: the walkers' speeds lie within [0, v_f]", lo >= 0 && hi <= w.free_speed + 1e-12, std::to_string(lo) + ", " + std::to_string(hi));
    s.require("street day: the signal is green about half the time", std::fabs(static_cast<double>(greens) / static_cast<double>(greens + reds) - 0.5) < 0.02,
              std::to_string(greens) + " green, " + std::to_string(reds) + " red");
    // A street of its own: blocks of 90, 160 and 70 m, the middle crossing unsignalled.
    StreetDay own = w;
    own.block_lengths = {90, 160, 70}, own.signalled = {1, 0, 1};
    const auto q = street_day(own);
    s.check("street day: a street of its own is as long as its blocks and crossings", q.interfaces.back() - (90 + 160 + 70 + 3 * own.crossing), 1e-9);
    s.check("street day: and conserves its crowd", q.budget_residual / std::max(1.0, q.entered.back()), 1e-10);
    bool free_middle = true;
    for (const auto& g : q.window_green) free_middle = free_middle && g[1] == 1;
    s.require("street day: an unsignalled crossing always passes walkers", free_middle);
  }

  // ---- currents past a counting line ----
  {
    calx::linalg::Matrix k(4);
    for (std::size_t i = 0; i < 4; i++) k(i, (i + 1) % 4) = 2.0 + 0.3 * static_cast<double>(i), k((i + 1) % 4, i) = 0.7 + 0.1 * static_cast<double>(i);
    const auto irr = physis::statistical_mechanics::irreversibility(k);
    double sum = 0;
    for (std::size_t i = 0; i < 4; i++) sum += edge_entropy_production(irr.pi[i] * k(i, (i + 1) % 4), irr.pi[(i + 1) % 4] * k((i + 1) % 4, i));
    s.check("currents: the edges' terms sum to Schnakenberg's σ (physis) on a driven ring", sum - irr.entropy_production, 1e-13);
    s.check("currents: a balanced edge has no term", edge_entropy_production(5, 5), 0);
    s.check("currents: (3 − 1) ln 3", edge_entropy_production(3, 1) - 2 * std::log(3.0), 1e-15);
    // Four weeks of counts toward the north peaking at 08:00 and toward the south at 17:00 (civil): the net
    // current is −300 sin(2π(τ − 12.5)/24) sin(2π·4.5/24), greatest at 06:30 and least at 18:30, its mean zero.
    urban::rhythm::Record d1, d2;
    d1.t0 = d2.t0 = urban::rhythm::parse_utc("2025-06-01T14:00:00Z");
    d1.offsets = d2.offsets = {{d1.t0 - 86400, 10.0}};
    const double two_pi = 2 * std::numbers::pi;
    for (std::size_t i = 0; i < 28 * 24; i++) {
      const double tau = d1.civil_hour(d1.time(i)) + 0.5;
      d1.y.push_back(200 + 150 * std::cos(two_pi * (tau - 8) / 24));
      d2.y.push_back(200 + 150 * std::cos(two_pi * (tau - 17) / 24));
    }
    const auto m = current_meter(d1, d2);
    s.require("currents: the flood at 06:30 and the ebb at 18:30 on the mean weekday", m.flood_at == 6.5 && m.ebb_at == 18.5,
              std::to_string(m.flood_at) + ", " + std::to_string(m.ebb_at));
    s.check("currents: no residual current", m.residual, 1e-9);
    s.check("currents: the flood's rate", m.flood - 300 * std::sin(two_pi * 4.5 / 24), 1e-9);
    s.require("currents: the exchange is irreversible", m.sigma > 0, std::to_string(m.sigma));
    // Handedness: four footpaths, each as left-hand traffic runs, one against; a footpath counting across its
    // side does not enter. P(X ≥ 3), X ~ binomial(4, ½) = 5/16.
    using C = Compass;
    const auto hand = handedness({C::north, C::south, C::west, C::east, C::north}, {{{C::east, C::west}}, {{C::east, C::west}}, {{C::north, C::south}}, {{C::north, C::south}}, {{C::north, C::south}}},
                                 {+3.0, -2.0, +1.0, +4.0, +9.0});
    s.require("currents: handedness counts the footpaths along their side, 3 of 4 left-handed", hand.footpaths == 4 && hand.left == 3,
              std::to_string(hand.left) + " of " + std::to_string(hand.footpaths));
    s.check("currents: the sign test's P(X ≥ 3) = 5/16", hand.p_left - 5.0 / 16, 1e-12);
    s.check("currents: and two-sided, P(X ≥ 3) + P(X ≤ 1) = 10/16", hand.p_two_sided - 10.0 / 16, 1e-12);
  }

  // ---- circulation (Theorem 6.3) ----
  {
    District d;
    d.nx = 4, d.ny = 3;
    const std::size_t V = (d.nx + 1) * (d.ny + 1);
    d.attraction.resize(V);
    for (auto& v : d.attraction) v = rng.uniform(-1, 1);
    const auto c0 = circulation(d);
    s.check("circulation: an attraction alone produces no entropy (detailed balance)", c0.irreversibility.entropy_production, 1e-10);
    s.require("circulation: and Kolmogorov's criterion holds", c0.irreversibility.detailed_balance);
    s.check("circulation: its log ratios are exact (no coexact part)", std::sqrt(c0.hodge.coexact_norm2 + c0.hodge.harmonic_norm2), 1e-9);
    double mean = 0;
    for (double v : d.attraction) mean += v / static_cast<double>(V);
    double err = 0;
    for (std::size_t i = 0; i < V; i++) err = std::max(err, std::fabs(c0.hodge.potential[i] - (d.attraction[i] - mean)));
    s.check("circulation: the Hodge potential is the attraction, up to a constant", err, 1e-8);
    d.block_circulation.assign(d.nx * d.ny, 0);
    for (auto& g : d.block_circulation) g = rng.uniform(0.2, 1.0);
    const auto c1 = circulation(d);
    s.require("circulation: circulation round the blocks produces entropy", c1.irreversibility.entropy_production > 1e-6 && !c1.irreversibility.detailed_balance,
              std::to_string(c1.irreversibility.entropy_production));
    s.require("circulation: with a coexact part", c1.hodge.coexact_norm2 > 1e-6);
    s.check("circulation: the attraction leaves every cycle affinity unchanged", c1.attraction_free_affinity_gap, 1e-9);
    double sigma = 0;
    for (const auto& cy : c1.irreversibility.cycles) sigma += cy.current * cy.affinity;
    s.check("circulation: σ = Σ J(C) A(C) (Schnakenberg)", (sigma - c1.irreversibility.entropy_production) / c1.irreversibility.entropy_production, 1e-9);
    District hole;
    hole.nx = 3, hole.ny = 3, hole.holes = {4}, hole.hole_circulation = {1.0};
    const auto c2 = circulation(hole);
    s.require("circulation: circulation round a hole is harmonic", c2.hodge.harmonic_norm2 > 1e-6, std::to_string(c2.hodge.harmonic_norm2));
  }

  // ---- transit ----
  {
    District one;
    one.nx = 1, one.ny = 1, one.exits = {0, 1, 2, 3}, one.exit_rate = 1.0 / 120, one.rate = 1.0 / 60;
    const auto t = transit(one, 1200);
    s.check("transit: every intersection an exit at rate 1/120 s gives an exponential of mean 120 s", t.mean - 120, 1e-6, "s");
    District g;
    g.nx = 4, g.ny = 3, g.exits = {0, 19}, g.entry = 7;
    const auto u = transit(g, 7200);
    bool falls = std::fabs(u.survival[0] - 1) < 1e-12;
    for (std::size_t i = 1; i < u.survival.size(); i++) falls = falls && u.survival[i] <= u.survival[i - 1] + 1e-12;
    s.require("transit: the survival starts at 1 and falls", falls);
  }
  return s.done();
}
