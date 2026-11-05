// Checks of the households of the centre (housing/market.hpp, households.hpp, calibration.hpp, planning.hpp), on small
// cases whose answers are known exactly.
//
// - The logit: the market's choice probabilities against their closed form; a household's choice among (building,
//   type) options, with ln n for a building's n vacant dwellings, equal to the choice among the dwellings one by one;
//   its draws' frequencies within four standard errors of those probabilities.
// - Clearing: with one class the rents that clear are explicit, and the solver finds them; with several classes and
//   types closed to some, demand equals the occupied stock; Berry's inversion recovers the segments' values; Roy's
//   identity holds for the surplus; in the agents' market, rents formed by excess demand settle where the vacancy is
//   natural.
// - Schelling's limit: with no rent, no type and no distance, a step utility at the threshold and the best option
//   taken, the households on a lattice move exactly as nomos::agents::Schelling's agents, tick by tick, draw for draw,
//   at thresholds 1/3 and 1/2; and sorting rises with the threshold.
// - Reproducibility: a run of a small city depends on its seed alone, and its event log, replayed on the state it
//   began from, gives the state it ended in.
// - The inverse demand: a panel made from a known slope and effects, with noise, gives the slope back within three
//   standard errors.
// - The planner's MOMDP: its beliefs equal those of the POMDP it is; with demand known, HSVI's interval holds the
//   MDP's value, which a two-state case gives in closed form; with demand hidden, HSVI's interval holds the value of
//   exact value iteration.

#include <urban/housing/calibration.hpp>
#include <urban/housing/households.hpp>
#include <urban/housing/market.hpp>
#include <urban/housing/planning.hpp>

#include <nomos/agents/choice.hpp>
#include <nomos/agents/engine.hpp>
#include <nomos/agents/space.hpp>

#include <tyche/decision/alpha.hpp>
#include <tyche/decision/point_based.hpp>
#include <tyche/decision/pomdp.hpp>
#include <tyche/decision/structured.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace uh = urban::housing;
namespace td = tyche::decision;
using calx::verification::relative;

namespace {

std::string fmt(const char* f, double a, double b = 0, double c = 0) {
  char buf[256];
  std::snprintf(buf, sizeof buf, f, a, b, c);
  return buf;
}

// A small market: three segments (apartments and houses in two areas, a student segment), two or three classes.
uh::Market small_market(bool several) {
  uh::Market m;
  m.segments = {{uh::apartment, 0, 400, 1.0, 0.5}, {uh::house, 0, 120, 0.4, 0.5}, {uh::apartment, 1, 300, 0.7, 2.0}, {uh::student, 0, 150, 0.2, 0.8}};
  const double ninf = -std::numeric_limits<double>::infinity();
  if (several) {
    m.classes = {{0.5, 600, {0.0, 0.0, 0.5}}, {0.3, 1500, {0.6, 0.2, ninf}}, {0.2, 3000, {1.0, 0.3, ninf}}};
  } else {
    m.classes = {{1.0, 1200, {0.0, 0.0, 0.0}}};
  }
  m.households = 5000, m.alpha = 2.5, m.mean_income = 1200, m.distance_cost = 0.3, m.vacancy = 0.05;
  return m;
}

// A lattice city: one apartment a cell, its graph the lattice's neighbourhood.
uh::City lattice_city(const nomos::agents::Lattice& L) {
  uh::City c;
  c.areas = 1, c.zones = 1;
  std::vector<std::pair<std::size_t, std::size_t>> edges;
  for (std::size_t i = 0; i < L.cells(); i++) {
    c.buildings.push_back({0, 0, static_cast<double>(i % L.width), static_cast<double>(i / L.width), 0});
    for (const std::size_t m : L.neighbours(i))
      if (m > i) edges.push_back({i, m});
  }
  c.graph = nomos::agents::Graph(L.cells(), edges);
  return c;
}

// A small city: buildings in zones and two areas, mixed capacities.
struct SmallCity {
  uh::City city;
  std::vector<std::array<int, uh::types>> capacity;
};
SmallCity small_city() {
  SmallCity s;
  s.city.areas = 2, s.city.zones = 5;
  calx::random::Stream rng(20261105, 0x63697479u, 0, 0);  // "city"
  for (std::size_t b = 0; b < 40; b++) {
    const std::size_t zone = b % 5, area = zone < 3 ? 0 : 1;
    s.city.buildings.push_back({area, zone, 100.0 * static_cast<double>(b % 8), 100.0 * static_cast<double>(b / 8), 0.2 + 0.1 * static_cast<double>(zone)});
    s.capacity.push_back({b % 4 == 0 ? 3 : 0, 4 + static_cast<int>(rng.below(20)), b % 10 == 3 ? 12 : 0});
  }
  return s;
}

uh::Preferences small_preferences(std::size_t segments) {
  uh::Preferences p;
  p.base.assign(segments, 0.0);
  for (std::size_t j = 0; j < segments; j++) p.base[j] = 40.0 + 0.2 * static_cast<double>(j % 2);  // α ln 400 ≈ 42
  p.income = {500, 1400, 3000};
  p.alpha = 7.0, p.mean_income = 1400, p.income_power = 0.25;
  p.taste_owner = {0.8, 0.0, -std::numeric_limits<double>::infinity()};
  p.taste_renter = {0.2, 0.0, 0.4};
  p.family_house = 0.5, p.distance_cost = 0.4, p.neighbourhood = 1.0, p.stay = 1.5, p.outside = 0;
  return p;
}

// A household of an income class, tenure and size, not yet resident.
uh::Household household(std::uint8_t income, std::uint8_t tenure, std::uint8_t size) {
  uh::Household h;
  h.income = income, h.tenure = tenure, h.size = size, h.active = true;
  return h;
}

std::vector<uh::Household> newcomers(std::size_t n, std::uint64_t seed, std::uint32_t t) {
  calx::random::Stream rng(seed, 0x6e657772u, t, 0);  // "newr"
  std::vector<uh::Household> out;
  for (std::size_t k = 0; k < n; k++) {
    const auto income = static_cast<std::uint8_t>(rng.below(3));
    const std::uint8_t tenure = rng.uniform() < 0.3 ? uh::owner : uh::renter, size = rng.uniform() < 0.5 ? 1 : 2;
    out.push_back(household(income, tenure, size));
  }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the logit ----------------------------------------------------------------------------------------------
  {
    const auto m = small_market(false);
    const std::vector<double> r = {500, 650, 420, 300};
    const auto p = uh::choice(m, 0, r);
    double denom = 1;
    std::vector<double> e(4);
    for (std::size_t j = 0; j < 4; j++) denom += (e[j] = std::exp(m.segments[j].base - m.alpha * std::log(r[j]) - m.distance_cost * m.segments[j].distance));
    double worst = std::fabs(p[4] - 1 / denom);
    for (std::size_t j = 0; j < 4; j++) worst = std::max(worst, std::fabs(p[j] - e[j] / denom));
    s.check("the market's choice probabilities against e^V / (1 + Σ e^V)", worst, 1e-15);
  }
  {
    // A newcomer among two buildings, three and one vacant dwellings, and the outside option.
    uh::City city;
    city.areas = 1, city.zones = 2;
    city.buildings = {{0, 0, 0, 0, 0.0}, {0, 1, 0, 0, 1.0}};
    uh::Preferences pref;
    pref.base = {0.0, 0.7, 0.0};
    pref.income = {1000};
    pref.alpha = 0, pref.distance_cost = 0.5, pref.outside = 0.2;
    pref.taste_renter = {0, 0, 0};
    const uh::Simulation sim(city, pref, {});
    auto st = sim.empty({{0, 3, 0}, {0, 1, 0}}, {400, 400, 400});
    st.households.push_back(household(0, uh::renter, 1));
    sim.begin(st);
    const auto O = sim.choices(st, 0);
    const auto P = nomos::agents::logit(O.utility, 1.0);
    // The dwellings one by one: three of utility 0.7, one of 0.7 − 0.5, and the outside option 0.2.
    const double e0 = std::exp(0.7), e1 = std::exp(0.2), eo = std::exp(0.2), D = 3 * e0 + e1 + eo;
    std::map<std::uint32_t, double> want = {{0 * 3 + 1, 3 * e0 / D}, {1 * 3 + 1, e1 / D}};
    double worst = std::fabs(P.back() - eo / D);
    for (std::size_t k = 0; k < O.option.size(); k++) worst = std::max(worst, std::fabs(P[k] - want[O.option[k]]));
    s.check("a choice among (building, type) with ln n equals the choice among its dwellings one by one", worst, 1e-15);
    // Its draws, over independent streams.
    const int n = 200000;
    std::map<std::uint32_t, int> got;
    int outside = 0;
    for (int k = 0; k < n; k++) {
      auto x = st;
      calx::random::Stream rng(99, 0x61676e74u, 0, static_cast<std::uint32_t>(k));
      sim.update(0, x, x, rng, nullptr);
      if (x.households[0].resident) got[x.households[0].building * 3 + x.households[0].type]++;
      else outside++;
    }
    double z = std::fabs(outside - n * eo / D) / std::sqrt(n * eo / D * (1 - eo / D));
    for (const auto& [o, q] : want) z = std::max(z, std::fabs(got[o] - n * q) / std::sqrt(n * q * (1 - q)));
    s.check("its draws' frequencies, in standard errors, over 200,000 streams", z, 4.0);
  }

  {
    // A zone of two buildings at different distances: the option (zone, apartment) carries ln Σ n_b e^{−θ d_b} and its
    // building is drawn in proportion, which is the choice among its five dwellings one by one.
    uh::City city;
    city.areas = 1, city.zones = 1;
    city.buildings = {{0, 0, 0, 0, 0.2}, {0, 0, 0, 0, 1.4}};
    uh::Preferences pref;
    pref.base = {0.0, 0.3, 0.0};
    pref.income = {1000};
    pref.alpha = 0, pref.distance_cost = 0.8, pref.outside = 0.0;
    pref.taste_renter = {0, 0, 0};
    const uh::Simulation sim(city, pref, {});
    auto st = sim.empty({{0, 2, 0}, {0, 3, 0}}, {400, 400, 400});
    st.households.push_back(household(0, uh::renter, 1));
    sim.begin(st);
    const double e0 = std::exp(0.3 - 0.8 * 0.2), e1 = std::exp(0.3 - 0.8 * 1.4), D = 2 * e0 + 3 * e1 + 1;
    const auto C = sim.choices(st, 0);
    const auto P = nomos::agents::logit(C.utility, 1.0);
    s.check("a zone's option: its probability equals the sum over its dwellings one by one", std::fabs(P[0] - (2 * e0 + 3 * e1) / D) + std::fabs(P[1] - 1 / D), 1e-15);
    // The two steps (a segment, then an option within it) give the same law.
    const auto Ps = sim.segment_choice(st, st.households[0], false);
    s.check("the choice in two steps, a segment then an option: the segment's probability", std::fabs(Ps[uh::apartment] - P[0]) + std::fabs(Ps.back() - P[1]), 1e-15);
    const int n = 200000;
    int in0 = 0, in1 = 0, none = 0;
    for (int k = 0; k < n; k++) {
      auto x = st;
      calx::random::Stream rng(98, 0x61676e74u, 0, static_cast<std::uint32_t>(k));
      sim.update(0, x, x, rng, nullptr);
      if (!x.households[0].resident) none++;
      else (x.households[0].building == 0 ? in0 : in1)++;
    }
    const auto z = [&](int got, double p) { return std::fabs(got - n * p) / std::sqrt(n * p * (1 - p)); };
    s.check("and its building drawn within the zone: frequencies in standard errors", std::max({z(in0, 2 * e0 / D), z(in1, 3 * e1 / D), z(none, 1 / D)}), 4.0);
  }

  // ---- clearing -------------------------------------------------------------------------------------------------
  {
    auto m = small_market(false);
    m.classes[0].taste = {0, 0, 0};
    const auto c = uh::clear(m);
    double O = 0;
    for (std::size_t j = 0; j < m.segments.size(); j++) O += m.occupied(j);
    double worst = 0;
    for (std::size_t j = 0; j < m.segments.size(); j++) {
      const auto& g = m.segments[j];
      const double lr = (g.base - m.distance_cost * g.distance - std::log(m.occupied(j) / (m.households - O))) / m.alpha;
      worst = std::max(worst, std::fabs(std::log(c.rents[j]) - lr));
    }
    s.require("one class: the solver converges", c.converged, std::to_string(c.iterations) + " iterations");
    s.check("one class: the clearing log rents against ln(O_j / (M − Σ O)) = δ_j − α ln r_j − θ d_j", worst, 1e-10);
  }
  {
    const auto m = small_market(true);
    const auto c = uh::clear(m);
    const auto d = uh::demand(m, c.rents);
    double worst = 0;
    for (std::size_t j = 0; j < m.segments.size(); j++) worst = std::max(worst, relative(d[j], m.occupied(j)));
    s.require("three classes, student apartments open to one: the solver converges", c.converged,
              std::to_string(c.iterations) + " iterations; rents " + fmt("%.1f, %.1f, %.1f", c.rents[0], c.rents[1], c.rents[2]) + fmt(", %.1f", c.rents[3]));
    s.check("three classes: demand equals the occupied stock, relatively", worst, 1e-10);
    // Berry's inversion recovers the bases from the rents and the stock.
    auto blank = m;
    for (auto& g : blank.segments) g.base = 0;
    std::vector<double> occ;
    for (std::size_t j = 0; j < m.segments.size(); j++) occ.push_back(m.occupied(j));
    const auto inv = uh::invert(blank, c.rents, occ);
    double dev = 0;
    for (std::size_t j = 0; j < m.segments.size(); j++) dev = std::max(dev, std::fabs(inv.base[j] - m.segments[j].base));
    s.check("Berry's inversion recovers the segments' values from the clearing rents", dev, 1e-9);
    // Roy's identity: ∂ surplus / ∂ ln r_j = −D_j.
    double roy = 0;
    for (std::size_t j = 0; j < m.segments.size(); j++) {
      auto up = c.rents, dn = c.rents;
      const double h = 1e-5;
      up[j] *= std::exp(h), dn[j] *= std::exp(-h);
      const double deriv = (uh::surplus_log_rent(m, up) - uh::surplus_log_rent(m, dn)) / (2 * h);
      roy = std::max(roy, relative(-deriv, d[j]));
    }
    s.check("Roy's identity: the surplus's derivative in ln r_j is −D_j, relatively", roy, 1e-7);
  }
  {
    // Rents formed by excess demand settle where the tick's demand meets the dwellings offered; the vacancy that
    // remains is turnover times the quarter a letting takes.
    auto sc = small_city();
    auto pref = small_preferences(sc.city.segments());
    uh::Dynamics dyn;
    dyn.exit = {0.01, 0.03}, dyn.hazard = {0.02, 0.08}, dyn.adjustment = 0.3, dyn.letting_lag = true;
    const uh::Simulation sim(sc.city, pref, dyn);
    auto st = sim.empty(sc.capacity, std::vector<double>(sc.city.segments(), 500.0));
    double gap = 0, vac = 0, drift = 0;
    int n = 0;
    std::vector<double> at200;
    for (std::uint32_t t = 0; t < 400; t++) {
      sim.begin(st);
      sim.arrive(st, newcomers(40, 11, t));
      sim.run(st, 11);
      if (t == 199) at200 = st.rent;
      if (t >= 200) {
        const auto cnt = uh::segment_counts(sc.city, st);
        double cap = 0, occ = 0, g = 0;
        int held = 0;
        for (std::size_t j = 0; j < cnt.capacity.size(); j++) {
          cap += cnt.capacity[j], occ += cnt.occupied[j];
          if (cnt.capacity[j] > 0) g += std::log((st.demand[j] + 1) / (st.offered[j] + 1)), held++;
        }
        gap += g / held;
        vac += 1 - occ / cap, n++;
      }
    }
    for (std::size_t j = 0; j < st.rent.size(); j++) drift = std::max(drift, std::fabs(std::log(st.rent[j] / at200[j])) / 200);
    s.check("rents formed by excess demand: the mean of ln((D + 1)/(V + 1)) over 200 quarters", gap / n, 0.05);
    s.check("and the rents settle: their largest drift a quarter over those 200 quarters", drift, 0.002);
    s.require("the vacancy that remains is turnover's", vac / n > 0 && vac / n < 0.2, fmt("mean vacancy %.4f", vac / n));
  }

  // ---- Schelling's limit --------------------------------------------------------------------------------------------
  for (const double f : {1.0 / 3, 0.5}) {
    const nomos::agents::Lattice L{20, 20, true, nomos::agents::Neighbourhood::moore};
    nomos::agents::Schelling S{L, f};
    auto x = S.place(300, 5);
    const auto city = lattice_city(L);
    uh::Preferences pref;
    pref.base.assign(city.segments(), 0.0);
    pref.income = {1, 1};
    pref.alpha = 0, pref.neighbourhood = 1, pref.step = true, pref.threshold = f;
    pref.scale = std::numeric_limits<double>::infinity(), pref.outside = -std::numeric_limits<double>::infinity();
    pref.taste_renter = {0, 0, 0};
    uh::Dynamics dyn;
    dyn.forced_by_neighbourhood = true;
    const uh::Simulation sim(city, pref, dyn);
    std::vector<std::array<int, uh::types>> cap(L.cells(), {0, 1, 0});
    auto st = sim.empty(cap, std::vector<double>(city.segments(), 1.0));
    for (std::size_t i = 0; i < x.where.size(); i++) {
      st.households.push_back(household(static_cast<std::uint8_t>(x.cell[x.where[i]]), uh::renter, 1));
      sim.place(st, i, x.where[i], uh::apartment);
    }
    bool same = true;
    int ticks = 0;
    for (std::uint32_t t = 0; t < 30 && same; t++, ticks++) {
      nomos::agents::tick(S, x, nomos::agents::Schedule::random_sequential, 2026, t);
      sim.step(st, 2026);
      std::vector<int> cells(L.cells(), -1);
      for (const auto& h : st.households) cells[h.building] = h.income;
      same = cells == x.cell;
    }
    s.require("Schelling's limit at f = " + fmt("%.3f", f) + ": the households' cells equal nomos's Schelling's, tick by tick", same,
              std::to_string(ticks) + " ticks; segregation " + fmt("%.3f", S.segregation(x)) + ", content " + fmt("%.3f", S.share_content(x)));
  }
  {
    // Sorting rises with the threshold (the mean share of like neighbours after 40 ticks).
    std::vector<double> seg;
    for (const double f : {0.0, 1.0 / 3, 0.5}) {
      const nomos::agents::Lattice L{30, 30, true, nomos::agents::Neighbourhood::moore};
      nomos::agents::Schelling S{L, f};
      const auto city = lattice_city(L);
      uh::Preferences pref;
      pref.base.assign(city.segments(), 0.0);
      pref.income = {1, 1};
      pref.alpha = 0, pref.neighbourhood = 1, pref.step = true, pref.threshold = f;
      pref.scale = std::numeric_limits<double>::infinity(), pref.outside = -std::numeric_limits<double>::infinity();
      uh::Dynamics dyn;
      dyn.forced_by_neighbourhood = true;
      const uh::Simulation sim(city, pref, dyn);
      auto x = S.place(720, 3);
      auto st = sim.empty(std::vector<std::array<int, uh::types>>(L.cells(), {0, 1, 0}), std::vector<double>(1 * uh::types, 1.0));
      for (std::size_t i = 0; i < x.where.size(); i++) {
        st.households.push_back(household(static_cast<std::uint8_t>(x.cell[x.where[i]]), uh::renter, 1));
        sim.place(st, i, x.where[i], uh::apartment);
      }
      for (std::uint32_t t = 0; t < 40; t++) sim.step(st, 77);
      std::fill(x.cell.begin(), x.cell.end(), -1);
      for (std::size_t i = 0; i < st.households.size(); i++) x.where[i] = st.households[i].building, x.cell[x.where[i]] = st.households[i].income;
      seg.push_back(S.segregation(x));
    }
    s.require("sorting rises with the threshold: like neighbours at f = 0, 1/3, 1/2", seg[0] < seg[1] && seg[1] < seg[2],
              fmt("%.3f, %.3f, %.3f", seg[0], seg[1], seg[2]));
  }

  // ---- reproducibility and the record -----------------------------------------------------------------------------
  {
    auto sc = small_city();
    auto pref = small_preferences(sc.city.segments());
    uh::Dynamics dyn;
    dyn.hazard = {0.03, 0.1}, dyn.burden = 0.6, dyn.adjustment = 0.2;
    const uh::Simulation sim(sc.city, pref, dyn);
    const auto run = [&](std::uint64_t seed, nomos::agents::Log* log, uh::State* start) {
      auto st = sim.empty(sc.capacity, std::vector<double>(sc.city.segments(), 450.0));
      for (std::uint32_t t = 0; t < 60; t++) {
        if (t == 20 && start) *start = st;
        sim.begin(st);
        sim.arrive(st, newcomers(30, seed, t));
        if (t == 30) sim.set_capacity(st, 7, uh::apartment, 0, t >= 20 ? log : nullptr);
        sim.run(st, seed, t >= 20 ? log : nullptr);
      }
      return st;
    };
    const auto a = run(7, nullptr, nullptr), b = run(7, nullptr, nullptr), c = run(8, nullptr, nullptr);
    const auto key = [](const uh::State& st) {
      std::vector<long> k;
      for (const auto& h : st.households) k.push_back(h.resident ? static_cast<long>(h.building * 3 + h.type) : -1);
      for (const double r : st.rent) k.push_back(std::lround(r * 1e6));
      return k;
    };
    s.require("a run depends on its seed alone: seed 7 twice the same, seed 8 different", key(a) == key(b) && key(a) != key(c));
    // The log replayed on the state at tick 20 gives the state at the end.
    nomos::agents::Log log;
    uh::State start;
    const auto end = run(7, &log, &start);
    std::vector<long> where(end.households.size(), -1);
    for (std::size_t h = 0; h < start.households.size(); h++)
      if (start.households[h].resident) where[h] = static_cast<long>(start.households[h].building * 3 + start.households[h].type);
    for (const auto& e : log.events) {
      if (e.rule == uh::rule::entered) where[e.agent] = static_cast<long>(e.touched[0]);
      else if (e.rule == uh::rule::moved) where[e.agent] = static_cast<long>(e.touched[1]);
      else if (e.rule == uh::rule::left || e.rule == uh::rule::evicted) where[e.agent] = -1;
    }
    std::vector<long> truth(end.households.size(), -1);
    for (std::size_t h = 0; h < end.households.size(); h++)
      if (end.households[h].resident) truth[h] = static_cast<long>(end.households[h].building * 3 + end.households[h].type);
    s.require("the event log, replayed on the state at tick 20, gives the state at tick 60", where == truth, std::to_string(log.events.size()) + " events");
  }

  // ---- the inverse demand -------------------------------------------------------------------------------------------
  {
    uh::RentPanel P;
    P.segments = 6, P.periods = 40;
    calx::random::Stream rng(31, 0x70616e6cu, 0, 0);  // "panl"
    const double slope = -0.15;
    for (std::size_t j = 0; j < 6; j++)
      for (std::size_t q = 0; q < 40; q++) {
        const double ls = 7 + 0.3 * static_cast<double>(j) + 0.02 * static_cast<double>(q) * static_cast<double>(j + 1) + 0.05 * rng.normal();
        P.segment.push_back(j), P.period.push_back(q), P.log_stock.push_back(ls);
        P.log_rent.push_back(6 + 0.01 * static_cast<double>(q) + 0.1 * static_cast<double>(j) + slope * ls + 0.02 * rng.normal());
      }
    const auto F = uh::inverse_demand(P);
    s.check("the inverse demand recovers its slope from a known panel, in standard errors (HC1)", (F.slope - slope) / F.se_hc1, 3.0);
  }

  // ---- the planner's MOMDP --------------------------------------------------------------------------------------------
  {
    // Demand known (one level): an MDP on the supply, and two levels solvable by hand.
    uh::PlannerSpec P;
    P.supply = 2, P.demand = 1, P.actions = 2, P.bands = 2;
    P.demand_transition = {1.0};
    P.welfare = {1.0, 3.0};
    P.rent = {600, 400};
    P.cuts = {500};
    P.confusion = 0.1, P.cost = 4.0, P.discount = 0.8, P.x0 = 0, P.y0 = {1.0};
    const auto M = uh::planner(P);
    const auto Q = td::to_pomdp(M);
    td::validate(Q);
    const auto mdp = td::solve_mdp(Q, 1e-12);
    // V(1) = 3 / (1 − γ) whatever is done (supply cannot fall); V(0) = max(1 + γ V(0), 1 − 4 + γ V(1)).
    const double V1 = 3 / (1 - 0.8), V0 = std::max(1 / (1 - 0.8), 1 - 4 + 0.8 * V1);
    s.check("demand known: the MDP's value against the closed form, V(0) = max(1/(1 − γ), 1 − c + γ 3/(1 − γ))", std::fabs(mdp.value[0] - V0) + std::fabs(mdp.value[1] - V1), 1e-9);
    const auto H = td::hsvi(Q, {.epsilon = 1e-6});
    s.require("demand known: HSVI's interval holds the MDP's value", H.lower_value <= V0 + 1e-9 && V0 <= H.upper_value + 1e-9,
              fmt("[%.6f, %.6f] against %.6f", H.lower_value, H.upper_value, V0));
  }
  {
    // Demand hidden: two levels of supply, two of demand, two bands of rent.
    uh::PlannerSpec P;
    P.supply = 2, P.demand = 2, P.actions = 2, P.bands = 2;
    P.demand_transition = {0.8, 0.2, 0.3, 0.7};
    P.welfare = {2.0, 0.5, 1.0, 3.0};  // (x, y): low supply suits low demand, high supply high demand
    P.rent = {450, 700, 350, 520};
    P.cuts = {500};
    P.confusion = 0.15, P.cost = 1.0, P.discount = 0.6, P.x0 = 0, P.y0 = {0.5, 0.5};
    const auto M = uh::planner(P);
    const auto Q = td::to_pomdp(M);
    td::validate(Q);
    // The MOMDP's belief equals the POMDP's.
    td::MomdpBelief b{0, {0.5, 0.5}};
    td::Belief bp = Q.initial;
    double dev = 0;
    for (int step = 0; step < 6; step++) {
      const std::size_t a = static_cast<std::size_t>(step % 2), x2 = std::min<std::size_t>(b.x + a, 1), o = static_cast<std::size_t>(step % 3 == 0);
      b = td::update(M, b, a, x2, o);
      bp = td::update(Q, bp, a, x2 * M.observations + o).belief;
      for (std::size_t y = 0; y < 2; y++) dev = std::max(dev, std::fabs(b.y[y] - bp[x2 * 2 + y]));
    }
    s.check("the MOMDP's belief over demand equals its POMDP's, over six steps", dev, 1e-14);
    const int horizon = 40;
    const auto V = td::exact_value_iteration(Q, horizon);
    double rmax = 0;
    for (const double r : Q.R) rmax = std::max(rmax, std::fabs(r));
    const double tail = std::pow(Q.discount, horizon) * rmax / (1 - Q.discount);
    const double exact = V.back().value(Q.initial);
    const auto H = td::hsvi(Q, {.epsilon = 1e-6});
    s.require("demand hidden: HSVI's interval holds exact value iteration's V*(b₀) (horizon 40, its tail added)",
              H.lower_value <= exact + tail + 1e-9 && exact - tail - 1e-9 <= H.upper_value && H.converged,
              fmt("[%.7f, %.7f] against %.7f", H.lower_value, H.upper_value, exact) + fmt(" ± %.1e", tail));
  }

  return s.done();
}
