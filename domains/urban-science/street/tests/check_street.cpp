// Checks of the street module (theory/street-oceanography.md §1–2):
//
// - Theorem 2.2: on random terrains, every non-negative weighting ranks every dominated pair the same way;
//   the covering pairs are the transitive reduction of the dominance order (their closure is the order, by
//   brute force); a change of a component's unit leaves the order and changes a contraction's ranking.
// - The use correlation of periodic blocks of length b (two uses alternating) is the triangle 1 − 2r/b, so
//   L_grad = b/4; a street of one use is a single regime with no compression ratio.
// - Theorems 2.7 and 2.8 on a street's holdings: M = between + within to rounding, and fifty steps of
//   exchange blind to kind never raise M.
// - The opening forcing's lines against a fine numerical transform of the sampled forcing, and the closed
//   form for one holding open Monday to Friday, 9 to 17.
// - Kendall's τ_b: 1 for one order, −1 for its reverse; the ascent and descent of holdings that open and
//   close in one order (τ = 1, small p) and of holdings that do not.
// - The visits: every arrival in an hour with something open goes somewhere; the arrivals in hours with
//   nothing open are counted unserved; a holding open the same hours every day, visited at a constant rate,
//   has no line power off the day's harmonics.

#include <urban/street/hours.hpp>
#include <urban/street/mixing.hpp>
#include <urban/street/street.hpp>
#include <urban/street/terrain.hpp>
#include <urban/street/visits.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <vector>

using namespace urban::street;

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(2026, 10, 30, 1);

  // ---- terrain: dominance and contractions ----
  {
    std::vector<std::vector<double>> T(24, std::vector<double>(3));
    for (auto& row : T)
      for (auto& v : row) v = std::floor(rng.uniform(0, 4));
    bool agree = true;
    for (int trial = 0; trial < 200; trial++) {
      Weighting w{"random", {rng.uniform(), rng.uniform(), rng.uniform()}};
      const auto c = contract(T, w);
      for (std::size_t i = 0; i < T.size(); i++)
        for (std::size_t j = 0; j < T.size(); j++)
          if (dominated(T[i], T[j]) && !(c.value[i] <= c.value[j] + 1e-12)) agree = false;
    }
    s.require("terrain: 200 non-negative weightings agree on every dominated pair (Theorem 2.2b)", agree);
    const auto cover = covering(T);
    const std::size_t n = T.size();
    std::vector<std::vector<char>> closure(n, std::vector<char>(n, 0));
    for (const auto& c : cover) closure[c.below][c.above] = 1;
    for (std::size_t k = 0; k < n; k++)
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++)
          if (closure[i][k] && closure[k][j]) closure[i][j] = 1;
    bool same = true;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) same = same && (closure[i][j] != 0) == dominated(T[i], T[j]);
    s.require("terrain: the closure of the covering pairs is the dominance order", same, std::to_string(cover.size()) + " covers");
    for (auto m : maximal(T))
      for (std::size_t j = 0; j < n; j++) same = same && !dominated(T[m], T[j]);
    s.require("terrain: the maximal elements are dominated by nothing", same);
    // A change of unit: two holdings, one wide and cheap, one narrow and dear.
    const std::vector<std::vector<double>> pair{{12, 1}, {6, 3}};
    Weighting w{"width in metres, price band", {1, 1}};
    Weighting w2{"width in tens of metres, price band", {0.1, 1}};
    const auto a = contract(pair, w), b = contract(pair, w2);
    s.require("terrain: a change of the width's unit reverses a contraction's ranking of incomparable holdings", a.rank[0] != b.rank[0] && agreement(a, b) == 0.0);
    s.require("terrain: and they are incomparable in the dominance order", !dominated(pair[0], pair[1]) && !dominated(pair[1], pair[0]));
  }

  // ---- the use correlation and the compression ratio ----
  {
    const double b = 8;
    std::vector<double> x;
    std::vector<std::size_t> u;
    for (int i = 0; i < 400; i++) x.push_back(i), u.push_back((i / static_cast<int>(b)) % 2);
    const auto c = use_correlation(x, {}, u, 1.0, 40);
    s.check("terrain length: alternating blocks of 8 have L_grad = b/4 = 2 m (binned at 1 m)", c.length - b / 4, 0.3, "m");
    s.require("terrain length: their correlation reaches zero near b/2", c.reached_zero && std::fabs(c.first_zero - b / 2) < 1.0, std::to_string(c.first_zero));
    s.check("compression: 400 m of walking over L_grad", compression(400, c) - 400 / c.length, 1e-12);
    std::vector<std::size_t> one(400, 0);
    const auto c1 = use_correlation(x, {}, one, 1.0, 40);
    s.require("compression: a street of one use is a single regime, with no ratio", c1.single_regime && std::isnan(compression(400, c1)));
  }

  // ---- mixing: the chain rule and exchange ----
  {
    Street st;
    st.id = "test", st.length = 600, st.uses = {"food", "fashion", "jewellery", "office"};
    for (int i = 0; i < 160; i++) {
      Holding h;
      h.id = std::to_string(i), h.s = rng.uniform(0, 600);
      h.use = h.s < 200 ? (rng.uniform() < 0.7 ? 0 : 3) : h.s < 400 ? (rng.uniform() < 0.6 ? 1 : 2) : static_cast<std::size_t>(rng.below(4));
      st.holdings.push_back(h);
    }
    const auto m = street_mixing(st, 25, 100);
    s.check("mixing: M = between + within across scales (Theorem 2.7)", m.scales.between + m.scales.within - m.mixing.departure, 1e-12);
    s.require("mixing: a street of three regimes departs from mixing", m.mixing.departure > 0.05, std::to_string(m.mixing.departure));
    auto counts = m.counts;
    double prev = m.mixing.departure;
    bool never_up = true;
    for (int k = 0; k < 50; k++) {
      counts = exchange_step(counts, 0.2);
      const double now = physis::statistical_mechanics::mixing(counts).departure;
      never_up = never_up && now <= prev + 1e-12;
      prev = now;
    }
    s.require("mixing: fifty steps of exchange blind to kind never raise M (Theorem 2.8)", never_up, std::to_string(prev));
  }

  // ---- hours ----
  {
    Street st;
    st.length = 100, st.uses = {"office"};
    Holding h;
    h.id = "office", h.s = 50, h.declared = true;
    for (int d = 0; d < 5; d++) h.hours.push_back({24.0 * d + 9, 24.0 * d + 17});
    st.holdings.push_back(h);
    const auto lines = week_lines(st, 21);
    s.check("hours: the mean line is the hours open over 168", lines[0].amplitude - 40.0 / 168, 1e-14);
    // F̂(7) in closed form: 5 days in one phase, (e^{−iν9} − e^{−iν17})/(iν)/168 each, ν = 2π/24.
    const double nu = 2 * std::numbers::pi / 24;
    const std::complex<double> one_day = (std::exp(std::complex<double>(0, -nu * 9)) - std::exp(std::complex<double>(0, -nu * 17))) / std::complex<double>(0, nu);
    s.check("hours: the daily line in closed form", std::abs(lines[7].value - 5.0 * one_day / 168.0), 1e-14);
    const auto f = opening_forcing(st, 168 * 60);
    double err = 0;
    for (std::size_t m = 1; m <= 21; m++) {
      std::complex<double> c{0, 0};
      for (std::size_t k = 0; k < f.size(); k++)
        c += f[k] * std::exp(std::complex<double>(0, -2 * std::numbers::pi * static_cast<double>(m) * (static_cast<double>(k) + 0.5) / static_cast<double>(f.size())));
      err = std::max(err, std::abs(c / static_cast<double>(f.size()) - lines[m].value));
    }
    s.check("hours: the exact lines against the sampled forcing's transform (minutes)", err, 2e-4);
    s.require("hours: the daily harmonics are flagged daily", lines[7].daily && lines[14].daily && !lines[1].daily);
  }
  {
    s.check("Kendall: one order has τ = 1", kendall_tau({1, 2, 3, 4}, {10, 20, 30, 40}) - 1, 1e-15);
    s.check("Kendall: its reverse has τ = −1", kendall_tau({1, 2, 3, 4}, {40, 30, 20, 10}) + 1, 1e-15);
    Street st;
    st.length = 100, st.uses = {"shop"};
    for (int i = 0; i < 10; i++) {
      Holding h;
      h.id = std::to_string(i), h.s = 10.0 * i, h.declared = true;
      h.hours.push_back({8.0 + 0.25 * i, 17.0 + 0.3 * i});   // Monday: first to open, first to close
      st.holdings.push_back(h);
    }
    const auto ad = ascent_descent(st, 0, 1999);
    s.require("ascent and descent: one order of opening and closing has τ = 1 and p < 0.01", std::fabs(ad.tau - 1) < 1e-12 && ad.p < 0.01, std::to_string(ad.p));
    s.check("ascent and descent: the asymmetry of dispersions 0.25 h and 0.30 h steps", ad.asymmetry - (0.3 - 0.25) / (0.3 + 0.25), 1e-12);
    Street split = st;
    split.holdings[0].hours = {{11.5, 14.5}, {17.5, 22.0}};
    const auto sd = ascent_descent(split, 0, 9);
    s.require("ascent and descent: a holding open twice a day opens at its first opening and closes at its last closing",
              sd.open[0] == 11.5 && sd.close[0] == 22.0, std::to_string(sd.open[0]) + "–" + std::to_string(sd.close[0]));
  }

  // ---- visits ----
  {
    Street st;
    st.length = 300, st.uses = {"cafe", "bridal", "jeweller"};
    const double day_hours[3][2] = {{7, 18}, {11, 19}, {10, 20}};
    for (int i = 0; i < 9; i++) {
      Holding h;
      h.id = std::to_string(i), h.s = 30.0 * i + 15, h.use = static_cast<std::size_t>(i % 3), h.declared = true;
      for (int d = 0; d < 7; d++) h.hours.push_back({24.0 * d + day_hours[i % 3][0], 24.0 * d + day_hours[i % 3][1]});
      st.holdings.push_back(h);
    }
    ArrivalType worker{"office worker", std::vector<double>(168, 20.0), {2.0, -1.0, 0.0}, 0.01, 0};
    ArrivalType couple{"couple", std::vector<double>(168, 5.0), {0.0, 2.0, 1.5}, 0.003, 150};
    const auto v = visits(st, {worker, couple});
    double served = 0, presence_open = 0;
    for (std::size_t k = 0; k < 168; k++) {
      bool any = false;
      for (const auto& h : st.holdings) any = any || open_at(h, k + 0.5);
      if (any) presence_open += 20.0;
      for (std::size_t h = 0; h < st.holdings.size(); h++) served += v.visits[0][h][k];
    }
    s.check("visits: every arrival in an open hour goes to some holding", served - presence_open, 1e-9);
    s.check("visits: the arrivals in closed hours are unserved", v.unserved[0] - (168 * 20.0 - presence_open), 1e-9);
    s.check("visits: a holding open the same hours each day, at constant presence, has no power off the daily lines", v.off_daily[1], 1e-12);
  }

  return s.done();
}
