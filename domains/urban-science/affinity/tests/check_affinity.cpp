// Checks of the affinity module (theory/street-oceanography.md §5):
//
// - Catalysis: the yield from physis's absorption probabilities against a Monte Carlo of the same process
//   (jump chain, 40 000 visitors): agreement within four standard errors; the yield falls with separation;
//   for a pure walk the coupling length grows as √D (doubling D multiplies ξ_c by √2 within 15%); a search
//   speed lengthens it; the dwell is phase-type with finite moments.
// - The coupling window: both stores open all day equals the time-homogeneous yield (the day long against
//   the dwell); the ring store closed through the visitor's day gives no yield.
// - The plume: its masses sum to one at every time; P(bought) rises to the yield.
// - Co-location: uses attracted within cells give CLQ > 1 against both nulls; uses that only favour the same
//   cells give CLQ > 1 against the free null and not against the stratified one.
// - Spacing: a periodic use of period 50 m is found at 50 m; a Turing example's wavelength is 2π/√k²_fastest.

#include <urban/affinity/catalysis.hpp>
#include <urban/affinity/colocation.hpp>
#include <urban/affinity/spacing.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace urban::affinity;

int main() {
  calx::verification::Suite s;

  // ---- catalysis ----
  Pathway p;
  p.length = 300, p.spacing = 10, p.diffusivity = 1.0, p.bridal = 100, p.ring = 140;
  {
    const double y = pathway_yield(p);
    // Monte Carlo of the jump chain: from each state jump to q with probability k(s, q)/Σ k(s, ·).
    const auto k = pathway_rates(p);
    const std::size_t S = k.n, bought = S - 2, left = S - 1;
    calx::random::Stream rng(1, 2, 3, 4);
    const std::size_t visitors = 40000;
    std::size_t hits = 0;
    for (std::size_t v = 0; v < visitors; v++) {
      std::size_t st = start_state(p);
      while (st != bought && st != left) {
        double total = 0;
        for (std::size_t q = 0; q < S; q++) total += k(st, q);
        double u = rng.uniform() * total;
        std::size_t q = 0;
        for (; q < S; q++) {
          u -= k(st, q);
          if (u <= 0 && k(st, q) > 0) break;
        }
        st = std::min(q, S - 1);
      }
      hits += st == bought;
    }
    const double mc = static_cast<double>(hits) / visitors, se = std::sqrt(mc * (1 - mc) / visitors);
    s.check("catalysis: the absorption probability against 40 000 simulated visitors (4 standard errors)", (y - mc) / (4 * se), 1.0);
    s.require("catalysis: the yield is a probability", y > 0 && y < 1, std::to_string(y));
    // The plume: its masses sum to one at every time, P(bought) never falls, and after six hours (twelve mean
    // stays on the street) it is the yield.
    const auto pl = pathway_plume(p, 6 * 3600, 360);
    double worst = 0;
    bool rising = true;
    for (std::size_t k2 = 0; k2 < pl.times.size(); k2++) {
      double m = pl.bought[k2] + pl.left[k2];
      for (std::size_t i = 0; i < pl.position.size(); i++) m += pl.errand[k2][i] + pl.need[k2][i];
      worst = std::max(worst, std::fabs(m - 1));
      if (k2 > 0 && pl.bought[k2] < pl.bought[k2 - 1] - 1e-15) rising = false;
    }
    s.check("plume: the visitor's masses sum to one at every time", worst, 1e-10);
    s.require("plume: the chance of having bought never falls", rising);
    s.check("plume: after six hours P(bought) is the yield", pl.bought.back() - y, 1e-4);
    s.require("plume: it starts at the bridal store's door with the errand", pl.errand[0][start_state(p) / 2] == 1 && pl.bought[0] == 0);
  }
  {
    // A long street, so that the ends are far from the stores: the yield beyond the neighbourhood decays as a
    // killed walk's hitting probability, ξ_c = −Δs / ln λ, λ + 1/λ = 2 + ℓΔs²/D.
    Pathway lp = p;
    lp.length = 1500, lp.bridal = 500;
    std::vector<double> d;
    for (double x = 10; x <= 300; x += 10) d.push_back(x);
    const auto c1 = coupling(lp, d);
    bool falls = true;
    for (std::size_t i = 1; i < c1.yield.size(); i++) falls = falls && c1.yield[i] <= c1.yield[i - 1] + 1e-12;
    s.require("catalysis: the yield falls with the separation of the stores", falls);
    const double xi1 = walk_coupling_length(lp.diffusivity, lp.leave_rate, lp.spacing);
    s.check("catalysis: the coupling length is a killed walk's e-folding length (D = 1 m²/s, relative)", (c1.coupling_length - xi1) / xi1, 2e-3);
    Pathway q = lp;
    q.diffusivity = 2.0;
    const auto c2 = coupling(q, d);
    const double xi2 = walk_coupling_length(q.diffusivity, q.leave_rate, q.spacing);
    s.check("catalysis: and with D = 2 m²/s, near √2 times longer (relative)", (c2.coupling_length - xi2) / xi2, 2e-3);
    s.check("catalysis: the closed form tends to √(D/ℓ) on a fine street (Δs = 0.1 m, relative)",
            walk_coupling_length(1, 1.0 / 1800, 0.1) / std::sqrt(1800.0) - 1, 1e-5);
    Pathway w = lp;
    w.search_speed = 0.2;
    const auto c3 = coupling(w, d);
    s.require("catalysis: a search speed toward the ring store lengthens the coupling", c3.coupling_length > c1.coupling_length || std::isnan(c3.coupling_length),
              std::to_string(c1.coupling_length) + " → " + std::to_string(c3.coupling_length));
    s.require("catalysis: the dwell is phase-type with a finite mean and median ≤ p90", c1.dwell_mean > 0 && c1.dwell_median <= c1.dwell_p90, std::to_string(c1.dwell_mean));
  }
  {
    Pathway f = p;
    f.leave_rate = 1.0 / 900;
    const double homogeneous = pathway_yield(f);
    // Arrivals every half hour through the day; the last ones are cut by the day's end, a loss of order
    // e^{−ℓ · 30 min} of one arrival in 48.
    const auto all_day = coupling_window(f, 0.0, 24.0, 0.0, 24.0, {0.0}, 120, 0.5);
    s.check("window: both stores open through the day equal the time-homogeneous yield (the day's end aside)", all_day.yield[0] - homogeneous, 5e-3);
    const auto shut = coupling_window(f, 10.0, 11.0, 0.0, 1.0, {0.0}, 120, 0.5);
    s.check("window: the ring store closed through the visitor's day gives no yield", shut.yield[0], 1e-12);
    const auto c = capture(homogeneous, 0.3, 2000, 3000);
    s.check("capture: B's contribution Y v_J beside its own takings p_B v_B", c.ratio - homogeneous * 3000 / 600, 1e-12);
  }

  // ---- co-location ----
  {
    calx::random::Stream rng(9, 9, 9, 9);
    std::vector<double> x, y;
    std::vector<std::size_t> u;
    // Twenty cells of 100 m. In each, bridal shops and jewellers (attraction within the cell): each bridal
    // shop has a jeweller within 5 m; other uses scattered.
    for (int c = 0; c < 20; c++) {
      const double cx = 100.0 * (c % 5), cy = 100.0 * (c / 5);
      for (int k = 0; k < 3; k++) {
        const double bx = cx + rng.uniform(10, 90), by = cy + rng.uniform(10, 90);
        x.push_back(bx), y.push_back(by), u.push_back(0);
        x.push_back(bx + rng.uniform(-3, 3)), y.push_back(by + rng.uniform(-3, 3)), u.push_back(1);
      }
      for (int k = 0; k < 10; k++) x.push_back(cx + rng.uniform(0, 100)), y.push_back(cy + rng.uniform(0, 100)), u.push_back(2);
    }
    const auto a = affinity(x, y, u, {"bridal", "jeweller", "other"}, 100, 0, 499, 3);
    s.require("co-location: bridal → jeweller attracted against the free null", a.free.clq[0 * 3 + 1] > 2 && a.free.p_upper[1] < 0.01, std::to_string(a.free.clq[1]));
    s.require("co-location: and against the stratified null", a.stratified.p_upper[1] < 0.01, std::to_string(a.stratified.p_upper[1]));
    // Uses that only favour the same cells: in half of the cells a use A, in the other half B and C mixed.
    std::vector<double> x2, y2;
    std::vector<std::size_t> u2;
    for (int c = 0; c < 20; c++) {
      const double cx = 100.0 * (c % 5), cy = 100.0 * (c / 5);
      for (int k = 0; k < 12; k++) {
        x2.push_back(cx + rng.uniform(0, 100)), y2.push_back(cy + rng.uniform(0, 100));
        u2.push_back(c % 2 == 0 ? 0 : (rng.uniform() < 0.5 ? 1 : 2));
      }
    }
    const auto b = affinity(x2, y2, u2, {"a", "b", "c"}, 100, 0, 499, 5);
    s.require("co-location: B → C attracted against the free null (they share cells)", b.free.p_upper[1 * 3 + 2] < 0.01, std::to_string(b.free.clq[5]));
    s.require("co-location: but not against the stratified null (independent response to terrain)", b.stratified.p_upper[1 * 3 + 2] > 0.01,
              std::to_string(b.stratified.p_upper[5]));
  }

  // ---- spacing ----
  {
    std::vector<double> s_pos;
    std::vector<std::size_t> u;
    calx::random::Stream rng(4, 4, 4, 4);
    for (int i = 0; i < 120; i++) {
      const double pos = 5.0 * i + rng.uniform(-1, 1);
      s_pos.push_back(pos);
      u.push_back(std::fmod(pos, 50.0) < 12 ? 1 : 0);
    }
    const auto o = observed_spacing(s_pos, u, 1, 600);
    s.check("spacing: a use of period 50 m is found at 50 m", o.best - 50, 2.0, "m");
    s.require("spacing: with a small false-alarm probability", o.false_alarm < 1e-3, std::to_string(o.false_alarm));
    // Schnakenberg-like kinetics at a steady state: J = [[1, 1], [−2, −1.5]], d_u = 1, d_v = 20.
    const auto t = predicted_spacing(1, 1, -2, -1.5, 1, 20);
    s.require("spacing: the example is Turing-unstable", t.turing.unstable);
    s.check("spacing: its wavelength is 2π/√k²_fastest", t.wavelength - 2 * std::numbers::pi / std::sqrt(t.turing.k2_fastest), 1e-12);
  }
  return s.done();
}
