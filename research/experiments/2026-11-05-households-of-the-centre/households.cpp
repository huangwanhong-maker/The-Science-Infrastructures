// The households of the centre: a case study of urban-science's housing on the City of Melbourne (ADR 0075).
//
//   households <housing-melbourne.json> <results directory> [calibrate|history|shock|schelling|planner|all, or several joined] [reuse]
//   (reuse: κ and π as an earlier run chose them, from results/trials.json)
//
// (1) Calibration and history: the inverse demand of the centre's segments from Homes Victoria's rents and CLUE's
//     stock; the households' classes, tenure and moves from the 2021 Census; the agents placed in 2002, their values
//     calibrated at 2002's rents, then run quarter by quarter to 2025Q3 with CLUE's stock and the demand that the
//     rents reveal; rents of new lettings against the report, and the households of 2021 against the Census.
// (2) A supply shock: the Development Activity Monitor's pipeline of towers or of student apartments, 2026-2030.
// (3) Schelling's limit: exactly on a lattice; in the city, sorting by block with prices off, against the Census.
// (4) The planner's permits under uncertain demand: a MOMDP, HSVI's interval and simple policies, and the cycles of
//     rents each makes.
// Results go to results/ as JSON; figures.py draws them.

#include "centre.hpp"

#include <urban/housing/planning.hpp>

#include <nomos/agents/choice.hpp>
#include <nomos/agents/space.hpp>
#include <tyche/decision/alpha.hpp>
#include <tyche/decision/point_based.hpp>
#include <tyche/decision/pomdp.hpp>
#include <tyche/decision/simulate.hpp>
#include <tyche/decision/structured.hpp>

#include <calx/signal/spectrum.hpp>

#include <chrono>
#include <functional>
#include <cstdio>
#include <iostream>
#include <thread>

using namespace centre;
namespace td = tyche::decision;

namespace {

std::string results;

void save(const std::string& name, const Json& j) {
  write_text(results + "/" + name, j.dump() + "\n");
  std::printf("  wrote %s\n", name.c_str());
}

double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

// ---- (1a) the inverse demand ------------------------------------------------------------------------------------

// The groups whose suburbs lie within the municipality (Flemington, most of South Yarra and of Port Melbourne lie
// outside it, so that CLUE's stock is not theirs).
const std::vector<std::string> inner_groups = {"CBD-St Kilda Rd", "Southbank", "Docklands", "Carlton-Parkville", "North Melbourne-West Melbourne", "East Melbourne"};

struct Calibration {
  uh::InverseDemand fit;
  std::vector<std::string> segment_name;     // the panel's segments
  std::vector<double> metro;                 // the outside option's rent: the metropolis's other regions, by report quarter
  std::vector<double> demand;                // the demand index by report quarter, exp(α (a_q − ln r̄_q)), 1 at the first
  double hazard_renter = 0, hazard_owner = 0, moved_share = 0;
  double exit_renter = 0, exit_owner = 0, from_outside = 0, within = 0;
  CensusTotals census;
  std::array<double, classes> class_share{};
  double owner_share = 0, owner_share_house = 0, owner_share_flat = 0, one_person = 0;
  double owner_share_newcomers = 0;  // the owners' share of the newcomers that keeps the Census's share among residents
  double harmonic_income = 0;
  Json panel;
};

// The outside option's rent: the lettings-weighted geometric mean of the medians of all properties in the metropolis's
// regions other than Inner Melbourne (which holds the centre), by report quarter.
const std::vector<std::string> metro_regions = {"Inner Eastern Melbourne", "Southern Melbourne", "Outer Western Melbourne", "North Western Melbourne",
                                                "North Eastern Melbourne", "Outer Eastern Melbourne", "South Eastern Melbourne"};
std::vector<double> metro_index(const Data& D) {
  const std::size_t Q = D.quarters.size();
  std::vector<double> w;
  for (const auto& r : metro_regions) {
    const auto& s = D.rents.at(r).at("All properties");
    double total = 0;
    for (std::size_t q = 0; q < Q; q++)
      if (std::isfinite(s.count[q])) total += s.count[q];
    w.push_back(total);
  }
  std::vector<double> out(Q, std::nan(""));
  for (std::size_t q = 0; q < Q; q++) {
    double num = 0, den = 0;
    for (std::size_t k = 0; k < metro_regions.size(); k++) {
      const double m = D.rents.at(metro_regions[k]).at("All properties").median[q];
      if (std::isfinite(m) && m > 0) num += w[k] * std::log(m), den += w[k];
    }
    if (den > 0) out[q] = std::exp(num / den);
  }
  return out;
}

Calibration calibrate(const Data& D) {
  Calibration C;
  const std::size_t Q = D.quarters.size();
  // Stock by group and type at a time.
  const auto stock = [&](const std::string& group, std::size_t type, double time) {
    double s = 0;
    for (const auto& b : D.buildings)
      if (!b.development && D.group_of_area[b.area] == group) s += stock_at(b, type, time);
    return s;
  };
  uh::RentPanel P;
  Json rows = Json::array();
  std::size_t seg = 0;
  for (const auto& g : inner_groups)
    for (const std::size_t type : {uh::apartment, uh::house}) {
      double lettings = 0;
      const auto idx = rent_index(D, g, type == uh::apartment ? flats : houses, &lettings);
      if (lettings < 80) continue;  // a year's lettings: too few houses let in most groups
      std::size_t n = 0;
      for (std::size_t q = 12; q < Q; q++) {  // 2003Q1 on: the moving year lies within CLUE's censuses
        if (!std::isfinite(idx[q])) continue;
        const double s = stock(g, type, quarter_time(static_cast<int>(q)) - 0.5);  // the middle of the year ending in q
        if (!(s > 0)) continue;
        P.segment.push_back(seg), P.period.push_back(q), P.log_rent.push_back(std::log(idx[q])), P.log_stock.push_back(std::log(s));
        rows.push(Json::array({static_cast<double>(seg), static_cast<double>(q), std::log(idx[q]), std::log(s)}));
        n++;
      }
      if (n) C.segment_name.push_back(g + (type == uh::apartment ? " flats" : " houses")), seg++;
    }
  P.segments = seg, P.periods = Q;
  C.fit = uh::inverse_demand(P);
  C.metro = metro_index(D);
  C.demand.assign(Q, std::nan(""));
  double first = std::nan("");
  for (std::size_t q = 0; q < Q; q++)
    if (std::isfinite(C.fit.period_effect[q])) {
      const double a = C.fit.period_effect[q] - std::log(C.metro[q]);
      if (!std::isfinite(first)) first = a;
      C.demand[q] = std::exp(C.fit.alpha * (a - first));
    }
  C.panel = Json::object({{"columns", Json::array({"segment", "quarter", "log_rent", "log_stock"})}, {"rows", rows}});
  // The Census.
  C.census = census_totals(D);
  const auto& T = C.census;
  double h = 0;
  for (std::size_t c = 0; c < classes; c++) h += T.households[c];
  for (std::size_t c = 0; c < classes; c++) C.class_share[c] = T.households[c] / h;
  C.owner_share = T.owners / (T.owners + T.renters);
  C.owner_share_house = T.owners_house / (T.owners_house + T.renters_house);
  C.owner_share_flat = T.owners_flat / (T.owners_flat + T.renters_flat);
  C.one_person = T.one_person / T.all_sizes;
  // Moves (G44): of the persons whose address a year before is stated or known to differ, those who came from another
  // SA2, from overseas or from an address not stated are taken as arrivals from outside the centre, whose exits match
  // them in a steady state; those who moved within their SA2 as moves within it (both assumptions: an SA2 is a part of
  // the centre). Owners are taken to leave and to move a quarter as often as renters (an assumption).
  {
    double same = 0, in_sa2 = 0, other_sa2 = 0, overseas = 0, unstated = 0;
    const auto& c44 = D.census_columns_g44;
    const auto g = [](double v) { return std::isfinite(v) ? v : 0.0; };
    for (const auto& a : D.sa1) {
      same += g(a.g44[column(c44, "Sme_Usl_ad_1_yr_ago_as_2021_P")]), in_sa2 += g(a.g44[column(c44, "Dif_Usl_ad_1_yr_ago_Sme_SA2_P")]);
      other_sa2 += g(a.g44[column(c44, "Dif_Us_ad_1_ago_Dif_SA2_Tot_P")]), overseas += g(a.g44[column(c44, "Difnt_Usl_add_1_yr_ago_OS_P")]);
      unstated += g(a.g44[column(c44, "Difnt_Usl_add_1_yr_ago_Ns_P")]);
    }
    const double all = same + in_sa2 + other_sa2 + overseas + unstated;
    C.moved_share = (all - same) / all, C.from_outside = (other_sa2 + overseas + unstated) / all, C.within = in_sa2 / all;
  }
  const double scale = 1 / ((1 - C.owner_share) + C.owner_share / 4);
  C.exit_renter = uh::quarterly_hazard(C.from_outside) * scale, C.exit_owner = C.exit_renter / 4;
  C.hazard_renter = uh::quarterly_hazard(C.within) * scale, C.hazard_owner = C.hazard_renter / 4;
  // Owners leave a quarter as often as renters, so that newcomers hold fewer owners than residents do: in a steady
  // state, inflow by tenure is proportional to stock times exit.
  C.owner_share_newcomers = C.owner_share * C.exit_owner / (C.owner_share * C.exit_owner + (1 - C.owner_share) * C.exit_renter);
  double inv = 0;
  for (std::size_t c = 0; c < classes; c++) inv += C.class_share[c] / T.median_income[c];
  C.harmonic_income = 1 / inv;
  return C;
}

// ---- (1b) the agents --------------------------------------------------------------------------------------------

struct Params {
  uh::Preferences pref;
  uh::Dynamics dyn;
  double arrivals = 1000;                 // newcomers a quarter at tick 0
  std::vector<double> demand_by_tick;     // the demand index at each tick (1 at tick 0)
  std::vector<double> reference_by_tick;  // the outside option's rent at each tick
  std::vector<double> rent0;              // the rents of 2002Q3 by segment
};

struct Run {
  uh::State state, census;                // at the end, and at 2021Q3
  std::vector<std::vector<double>> rents; // [tick][segment], the rents let at
  std::vector<std::vector<int>> lettings; // [tick][segment]
  std::vector<uh::TickRecord> record;
  Params params;
  double burn_acceptance = 0, burn_gap = 0;
};

uh::Household draw_household(calx::random::Stream& rng, const Calibration& C, double owner_share, bool student) {
  uh::Household h;
  const double u = rng.uniform();
  h.income = static_cast<std::uint8_t>(u < C.class_share[0] ? 0 : (u < C.class_share[0] + C.class_share[1] ? 1 : 2));
  h.tenure = !student && rng.uniform() < owner_share ? uh::owner : uh::renter;
  h.size = student || rng.uniform() < C.one_person ? 1 : 2;
  h.active = true;
  return h;
}

std::vector<uh::Household> newcomers(const Calibration& C, std::size_t n, std::uint64_t seed, std::uint32_t t) {
  calx::random::Stream rng(seed, 0x61727276u, t, 0);  // "arrv"
  std::vector<uh::Household> out;
  out.reserve(n);
  for (std::size_t k = 0; k < n; k++) out.push_back(draw_household(rng, C, C.owner_share_newcomers, false));
  return out;
}

// Rents of 2002Q3 by segment: apartments the group's flats index; houses the group's houses where the report
// publishes them, else the group's flats scaled by Inner Melbourne's ratio of houses to flats; student apartments the
// group's one-bedroom flats (an assumption: the report holds no student apartments).
std::vector<double> initial_rents(const Data& D, const uh::City& city, int quarter) {
  std::vector<double> r(city.segments(), 0.0);
  const auto inner_flats = rent_index(D, "Inner Melbourne", flats), inner_houses = rent_index(D, "Inner Melbourne", houses);
  for (std::size_t a = 0; a < city.areas; a++) {
    const auto& g = D.group_of_area[a];
    const auto f = rent_index(D, g, flats), hs = rent_index(D, g, houses), one = rent_index(D, g, {"1 bedroom flat"});
    const auto at = [&](const std::vector<double>& v) {  // the quarter's value, or the nearest published
      for (int d = 0; d < 12; d++) {
        if (quarter + d < static_cast<int>(v.size()) && std::isfinite(v[static_cast<std::size_t>(quarter + d)])) return v[static_cast<std::size_t>(quarter + d)];
        if (quarter - d >= 0 && std::isfinite(v[static_cast<std::size_t>(quarter - d)])) return v[static_cast<std::size_t>(quarter - d)];
      }
      return std::nan("");
    };
    double fl = at(f);
    if (!std::isfinite(fl)) fl = at(inner_flats);
    double ho = at(hs);
    if (!std::isfinite(ho)) ho = fl * at(inner_houses) / at(inner_flats);
    double st = at(one);
    if (!std::isfinite(st)) st = fl;
    r[uh::apartment * city.areas + a] = fl, r[uh::house * city.areas + a] = ho, r[uh::student * city.areas + a] = st;
  }
  return r;
}

Params base_params(const Data& D, const uh::City& city, const Calibration& C, double power, double kappa) {
  Params P;
  auto& p = P.pref;
  p.income.assign(C.census.median_income.begin(), C.census.median_income.end());
  p.alpha = C.fit.alpha, p.income_power = power;
  // The mean income at which α_c = α: the one that makes the classes' mean coefficient α (Σ_c w_c (ȳ/y_c)^π = 1).
  {
    double m = 0;
    for (std::size_t c = 0; c < classes; c++) m += C.class_share[c] * std::pow(1 / C.census.median_income[c], power);
    p.mean_income = std::pow(1 / m, 1 / std::max(power, 1e-12));
    if (power == 0) p.mean_income = 1;
  }
  p.taste_owner = {0.0, 0.0, -std::numeric_limits<double>::infinity()};
  p.taste_renter = {0.0, 0.0, 0.0};
  p.family_house = 0, p.student_single_only = true;
  p.distance_cost = 0.3;  // per km (an assumption: not identified by the aggregates)
  p.neighbourhood = 0, p.stay = 0, p.outside = 0;
  P.rent0 = initial_rents(D, city, quarter_of_tick(0));
  P.reference_by_tick.resize(ticks);
  for (int t = 0; t < ticks; t++) P.reference_by_tick[static_cast<std::size_t>(t)] = C.metro[static_cast<std::size_t>(quarter_of_tick(t))];
  p.base.assign(city.segments(), 0.0);
  for (std::size_t j = 0; j < city.segments(); j++) p.base[j] = C.fit.alpha * std::log(P.rent0[j] / P.reference_by_tick[0]) - 4;
  P.dyn.hazard = {C.hazard_owner, C.hazard_renter};
  P.dyn.exit = {C.exit_owner, C.exit_renter};
  P.dyn.adjustment = kappa, P.dyn.cap = 0.15, P.dyn.letting_lag = true;
  P.demand_by_tick.resize(ticks);
  double last = 1;
  for (int t = 0; t < ticks; t++) {
    const double d = C.demand[static_cast<std::size_t>(quarter_of_tick(t))];
    if (std::isfinite(d)) last = d;
    P.demand_by_tick[static_cast<std::size_t>(t)] = last;
  }
  // The index is 1 at its first quarter (2003Q1); at tick 0 (2002Q3), the first value is held.
  const double d0 = P.demand_by_tick[0];
  for (auto& d : P.demand_by_tick) d /= d0;
  P.arrivals = 0;  // set by the burn-in
  return P;
}

// The households of 2002Q3, placed in its stock at 96 % occupancy (an assumption), their classes, tenure (by the
// Census's share of owners among houses and among flats) and sizes drawn as the 2021 Census holds them (an
// assumption for 2002).
uh::State initial_state(const uh::Simulation& sim, const Data& D, const Calibration& C, const Params& P, std::uint64_t seed) {
  auto st = sim.empty(capacities(D, 0), P.rent0);
  calx::random::Stream rng(seed, 0x696e6974u, 0, 0);  // "init"
  for (std::size_t b = 0; b < D.buildings.size(); b++)
    for (std::size_t t = 0; t < uh::types; t++) {
      const int n = static_cast<int>(std::floor(0.96 * st.capacity[b][t]));
      for (int k = 0; k < n; k++) {
        st.households.push_back(draw_household(rng, C, t == uh::house ? C.owner_share_house : C.owner_share_flat, t == uh::student));
        sim.place(st, st.households.size() - 1, b, t);
      }
    }
  return st;
}

double owners_in_houses(const uh::State& st) {
  double o = 0, n = 0;
  for (const auto& h : st.households)
    if (h.resident && h.type == uh::house) n++, o += h.tenure == uh::owner;
  return n > 0 ? o / n : 0;
}

// The burn-in: rents held at 2002Q3's, the segments' values moved until each segment's demand meets its offer
// (Berry's contraction read on the flow of lettings), and the owners' taste for houses until owners hold the Census's
// share of houses. Newcomers a quarter are twenty times the residents who leave (an assumption: one in twenty of the
// households who look at the centre in a quarter settle there, as the planner's market holds it, so that the centre is a
// small part of the households' choices, as the inverse demand assumes, and its population is steady).
Run burn_in(const Data& D, const uh::City& city, const Calibration& C, Params P, std::uint64_t seed) {
  {
    double n = 0;
    for (const auto& b : D.buildings) n += std::floor(0.96 * b.stock[0][0]) + std::floor(0.96 * b.stock[1][0]) + std::floor(0.96 * b.stock[2][0]);
    P.arrivals = 20 * n * (C.owner_share * C.exit_owner + (1 - C.owner_share) * C.exit_renter);
  }
  uh::Dynamics frozen = P.dyn;
  frozen.adjustment = 0;
  Run R;
  for (int round = 0; round < 14; round++) {
    const uh::Simulation sim(city, P.pref, frozen);
    auto st = round == 0 ? initial_state(sim, D, C, P, seed) : R.state;
    st.rent = P.rent0, st.reference_rent = P.reference_by_tick[0];
    std::vector<double> dem(city.segments(), 0.0), off(city.segments(), 0.0);
    double arrived = 0, entered = 0;
    const int T = 6;
    for (int k = 0; k < T; k++) {
      sim.compact(st);
      sim.begin(st);
      sim.arrive(st, newcomers(C, static_cast<std::size_t>(std::lround(P.arrivals)), seed, static_cast<std::uint32_t>(1000 + round * T + k)));
      sim.run(st, seed + 7919 * static_cast<std::uint64_t>(round));
      for (std::size_t j = 0; j < city.segments(); j++) dem[j] += st.demands.back()[j] / T, off[j] += st.offers.back()[j] / T;
      arrived += st.ticks.back().arrivals, entered += st.ticks.back().entered;
    }
    // Drop the households who left, and the record of the burn-in.
    double gap = 0;
    int n = 0;
    for (std::size_t j = 0; j < city.segments(); j++) {
      if (off[j] <= 0 && dem[j] <= 0) continue;
      const double g = std::log((off[j] + 1) / (dem[j] + 1));
      P.pref.base[j] += std::clamp(g, -1.0, 1.0);
      gap += std::fabs(g), n++;
    }
    const double target = C.owner_share_house, got = owners_in_houses(st);
    if (got > 0 && got < 1) P.pref.taste_owner[uh::house] += std::log(target / (1 - target)) - std::log(got / (1 - got));
    const double acc = arrived > 0 ? entered / arrived : 0.5;
    R.burn_acceptance = acc, R.burn_gap = gap / std::max(n, 1);
    R.state = std::move(st);
  }
  // A fresh record from tick 0.
  auto& st = R.state;
  std::vector<uh::Household> kept;
  for (const auto& h : st.households)
    if (h.resident) kept.push_back(h);
  const uh::Simulation sim(city, P.pref, P.dyn);
  auto fresh = sim.empty(capacities(D, 0), P.rent0);
  for (auto h : kept) {
    h.resident = false;
    fresh.households.push_back(h);
    sim.place(fresh, fresh.households.size() - 1, h.building, h.type);
  }
  for (auto& h : fresh.households) h.since = -1;
  fresh.reference_rent = P.reference_by_tick[0];
  R.state = std::move(fresh);
  R.params = P;
  return R;
}

// The history, 2002Q3 to 2025Q3, from the burn-in's state.
Run history(const Data& D, const uh::City& city, const Calibration& C, const Run& burn, std::uint64_t seed) {
  Run R;
  R.params = burn.params;
  const auto& P = R.params;
  const uh::Simulation sim(city, P.pref, P.dyn);
  R.state = burn.state;
  auto& st = R.state;
  for (int t = 0; t < ticks; t++) {
    st.reference_rent = P.reference_by_tick[static_cast<std::size_t>(t)];
    sim.compact(st);
    sim.begin(st);
    if (t > 0) {
      const auto cap = capacities(D, t);
      for (std::size_t b = 0; b < cap.size(); b++)
        for (std::size_t k = 0; k < uh::types; k++)
          if (cap[b][k] != st.capacity[b][k]) sim.set_capacity(st, b, k, cap[b][k]);
    }
    sim.arrive(st, newcomers(C, static_cast<std::size_t>(std::lround(P.arrivals * P.demand_by_tick[static_cast<std::size_t>(t)])), seed, static_cast<std::uint32_t>(t)));
    sim.run(st, seed);
    if (t == census_tick) R.census = st;
  }
  R.rents = st.rents, R.lettings = st.lettings, R.record = st.ticks;
  return R;
}

// ---- measures ---------------------------------------------------------------------------------------------------

// A group's moving annual median rent of new lettings in a type, by tick (from tick 3), over its areas' segments.
std::vector<double> lettings_median(const Data& D, const uh::City& city, const Run& R, const std::string& group, std::size_t type) {
  std::vector<double> out(R.rents.size(), std::nan(""));
  for (std::size_t t = 3; t < R.rents.size(); t++) {
    std::vector<std::pair<double, double>> v;
    double total = 0;
    for (std::size_t s = t - 3; s <= t; s++)
      for (std::size_t a = 0; a < city.areas; a++)
        if (D.group_of_area[a] == group) {
          const std::size_t j = type * city.areas + a;
          if (R.lettings[s][j] > 0) v.push_back({R.rents[s][j], static_cast<double>(R.lettings[s][j])}), total += R.lettings[s][j];
        }
    if (v.empty()) continue;
    std::sort(v.begin(), v.end());
    double run = 0;
    for (const auto& [r, w] : v)
      if ((run += w) >= total / 2) { out[t] = r; break; }
  }
  return out;
}

// The centre's index: the groups' medians weighted by their lettings, geometrically.
std::vector<double> centre_index(const std::vector<std::vector<double>>& by_group, const std::vector<double>& weights) {
  std::vector<double> out(by_group[0].size(), std::nan(""));
  for (std::size_t t = 0; t < out.size(); t++) {
    double num = 0, den = 0;
    for (std::size_t g = 0; g < by_group.size(); g++)
      if (std::isfinite(by_group[g][t])) num += weights[g] * std::log(by_group[g][t]), den += weights[g];
    if (den > 0) out[t] = std::exp(num / den);
  }
  return out;
}

// Each building's SA1 (the SA1 whose ring holds its point, else the nearest SA1's point; −1 beyond 2 km).
std::vector<int> sa1_of_buildings(const Data& D) {
  std::vector<int> out(D.buildings.size(), -1);
  for (std::size_t b = 0; b < D.buildings.size(); b++) {
    const double x = D.buildings[b].x, y = D.buildings[b].y;
    for (std::size_t s = 0; s < D.sa1.size() && out[b] < 0; s++)
      for (const auto& r : D.sa1[s].rings)
        if (inside(r, x, y)) { out[b] = static_cast<int>(s); break; }
    if (out[b] < 0) {
      double best = 2000;
      for (std::size_t s = 0; s < D.sa1.size(); s++) {
        const double d = std::hypot(D.sa1[s].x - x, D.sa1[s].y - y);
        if (d < best) best = d, out[b] = static_cast<int>(s);
      }
    }
  }
  return out;
}

// Duncan and Duncan's dissimilarity of two groups over units: ½ Σ |a_i/A − b_i/B| (Duncan and Duncan 1955, "A
// methodological analysis of segregation indexes", Am. Sociol. Rev. 20, 210–217).
double dissimilarity(const std::vector<double>& a, const std::vector<double>& b) {
  double A = 0, B = 0, d = 0;
  for (std::size_t i = 0; i < a.size(); i++) A += a[i], B += b[i];
  for (std::size_t i = 0; i < a.size(); i++) d += std::fabs(a[i] / A - b[i] / B);
  return d / 2;
}

// Households by SA1 and income class, and renters and owners, of a state.
struct BySa1 {
  std::vector<std::array<double, classes>> income;
  std::vector<double> renters, owners, households, moved;
  std::vector<std::vector<double>> rent;  // renters' current rents
};
BySa1 by_sa1(const Data& D, const uh::City& city, const uh::State& st, const std::vector<int>& sa1_of, int tick) {
  BySa1 out;
  const std::size_t S = D.sa1.size();
  out.income.assign(S, {0, 0, 0}), out.renters.assign(S, 0), out.owners.assign(S, 0), out.households.assign(S, 0), out.moved.assign(S, 0), out.rent.assign(S, {});
  for (const auto& h : st.households) {
    if (!h.resident || sa1_of[h.building] < 0) continue;
    const auto s = static_cast<std::size_t>(sa1_of[h.building]);
    out.income[s][h.income]++, out.households[s]++;
    (h.tenure == uh::renter ? out.renters : out.owners)[s]++;
    if (h.since > tick - 4) out.moved[s]++;
    if (h.tenure == uh::renter) out.rent[s].push_back(st.rent[city.segment(h.building, h.type)]);
  }
  return out;
}

double median(std::vector<double> v) {
  if (v.empty()) return std::nan("");
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double census_d(const Data& D) {
  std::vector<double> lo, hi;
  for (const auto& s : D.sa1) {
    const auto c = income_classes(D, s);
    lo.push_back(c[0]), hi.push_back(c[2]);
  }
  return dissimilarity(lo, hi);
}

double model_d(const BySa1& m) {
  std::vector<double> lo, hi;
  for (const auto& c : m.income) lo.push_back(c[0]), hi.push_back(c[2]);
  return dissimilarity(lo, hi);
}

Json vec(const std::vector<double>& v) {
  Json a = Json::array();
  for (const double x : v) a.push(std::isfinite(x) ? Json(x) : Json());
  return a;
}

// The fit of a history to the report: the log error of the centre's flats index, 2003Q1 to 2025Q3.
struct RentFit {
  std::vector<std::vector<double>> model, observed;  // by inner group, by tick
  std::vector<double> centre_model, centre_observed, weights;
  double rmse = 0, rmse_groups = 0;
};
RentFit rent_fit(const Data& D, const uh::City& city, const Run& R) {
  RentFit F;
  for (const auto& g : inner_groups) {
    double lettings = 0;
    const auto idx = rent_index(D, g, flats, &lettings);
    std::vector<double> obs(ticks);
    for (int t = 0; t < ticks; t++) obs[static_cast<std::size_t>(t)] = idx[static_cast<std::size_t>(quarter_of_tick(t))];
    F.observed.push_back(obs), F.model.push_back(lettings_median(D, city, R, g, uh::apartment)), F.weights.push_back(lettings);
  }
  F.centre_model = centre_index(F.model, F.weights), F.centre_observed = centre_index(F.observed, F.weights);
  double e = 0, eg = 0;
  int n = 0, ng = 0;
  for (int t = 2; t < ticks; t++) {
    const auto u = static_cast<std::size_t>(t);
    if (std::isfinite(F.centre_model[u]) && std::isfinite(F.centre_observed[u])) e += std::pow(std::log(F.centre_model[u] / F.centre_observed[u]), 2), n++;
    for (std::size_t g = 0; g < F.model.size(); g++)
      if (std::isfinite(F.model[g][u]) && std::isfinite(F.observed[g][u])) eg += std::pow(std::log(F.model[g][u] / F.observed[g][u]), 2), ng++;
  }
  F.rmse = std::sqrt(e / std::max(n, 1)), F.rmse_groups = std::sqrt(eg / std::max(ng, 1));
  return F;
}

// Runs f(seed) for each seed on up to three threads.
template <class F>
void parallel(const std::vector<std::uint64_t>& seeds, F&& f) {
  for (std::size_t i = 0; i < seeds.size(); i += 3) {
    std::vector<std::thread> pool;
    for (std::size_t k = i; k < std::min(seeds.size(), i + 3); k++) pool.emplace_back([&, k] { f(k, seeds[k]); });
    for (auto& t : pool) t.join();
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "households <housing-melbourne.json> <results directory> [calibrate|history|shock|schelling|planner|all]\n");
    return 2;
  }
  results = argv[2];
  const std::string task = argc > 3 ? argv[3] : "all";
  const auto t0 = std::chrono::steady_clock::now();
  const Data D = load(argv[1]);
  const uh::City city = make_city(D);
  std::printf("the centre: %zu buildings (%zu of the pipeline), %zu blocks, %zu areas, %zu SA1s, %zu developments\n", D.buildings.size(),
              static_cast<std::size_t>(std::count_if(D.buildings.begin(), D.buildings.end(), [](const BuildingRecord& b) { return b.development; })), D.blocks.size(),
              D.areas.size(), D.sa1.size(), D.developments.size());

  // ---- (1a) ----
  const Calibration C = calibrate(D);
  std::printf("inverse demand: slope %.4f (HC1 %.4f, clustered %.4f), α = %.2f ± %.2f, %zu observations in %zu segments, R² %.3f\n", C.fit.slope, C.fit.se_hc1,
              C.fit.se_cluster, C.fit.alpha, C.fit.alpha_se, C.fit.observations, C.fit.clusters, C.fit.r_squared);
  std::printf("census: classes %.3f %.3f %.3f (median incomes %.0f %.0f %.0f), owners %.3f (houses %.3f, flats %.3f), one person %.3f, moved %.3f (from outside %.3f, "
              "within %.3f) → exits %.4f %.4f, moves %.4f %.4f a quarter\n",
              C.class_share[0], C.class_share[1], C.class_share[2], C.census.median_income[0], C.census.median_income[1], C.census.median_income[2], C.owner_share,
              C.owner_share_house, C.owner_share_flat, C.one_person, C.moved_share, C.from_outside, C.within, C.exit_owner, C.exit_renter, C.hazard_owner, C.hazard_renter);
  {
    Json seg = Json::array();
    for (const auto& s : C.segment_name) seg.push(s);
    Json q = Json::array();
    for (const auto& s : D.quarters) q.push(s);
    save("calibration.json",
         Json::object({{"metro_index", vec(C.metro)},
                       {"inverse_demand", Json::object({{"slope", C.fit.slope}, {"se_hc1", C.fit.se_hc1}, {"se_cluster", C.fit.se_cluster}, {"alpha", C.fit.alpha},
                                                         {"alpha_se", C.fit.alpha_se}, {"r_squared", C.fit.r_squared}, {"observations", static_cast<double>(C.fit.observations)},
                                                         {"segments", seg}, {"segment_effect", vec(C.fit.segment_effect)}, {"period_effect", vec(C.fit.period_effect)}})},
                       {"panel", C.panel},
                       {"quarters", q},
                       {"demand_index", vec(C.demand)},
                       {"census", Json::object({{"class_share", vec({C.class_share[0], C.class_share[1], C.class_share[2]})},
                                                {"median_income", vec({C.census.median_income[0], C.census.median_income[1], C.census.median_income[2]})},
                                                {"harmonic_income", C.harmonic_income},
                                                {"owner_share", C.owner_share}, {"owner_share_house", C.owner_share_house}, {"owner_share_flat", C.owner_share_flat},
                                                {"one_person", C.one_person}, {"moved_share", C.moved_share}, {"hazard_owner", C.hazard_owner}, {"hazard_renter", C.hazard_renter},
                                                {"from_outside", C.from_outside}, {"within", C.within}, {"exit_owner", C.exit_owner}, {"exit_renter", C.exit_renter},
                                                {"owner_share_newcomers", C.owner_share_newcomers},
                                                {"occupied", C.census.occupied}, {"unoccupied", C.census.unoccupied}, {"d_low_high", census_d(D)}})}}));
  }
  if (task == "calibrate") return 0;

  // ---- (1b) the agents: κ by the rents, η by the Census's sorting ----
  const auto sa1_of = sa1_of_buildings(D);
  const double d_census = census_d(D);
  struct Trial {
    double power, kappa, rmse, rmse_groups, d, acceptance, gap, arrivals, seconds;
  };
  std::vector<Trial> trials;
  std::map<double, Run> burns;  // by π
  const auto trial = [&](double power, double kappa, std::uint64_t seed) {
    const auto t1 = std::chrono::steady_clock::now();
    if (!burns.count(power)) burns[power] = burn_in(D, city, C, base_params(D, city, C, power, kappa), seed);
    Run b = burns[power];
    b.params.dyn.adjustment = kappa;
    const auto R = history(D, city, C, b, seed);
    const auto F = rent_fit(D, city, R);
    const auto m = by_sa1(D, city, R.census, sa1_of, census_tick);
    Trial x{power, kappa, F.rmse, F.rmse_groups, model_d(m), b.burn_acceptance, b.burn_gap, b.params.arrivals, seconds_since(t1)};
    std::printf("  π %.2f κ %.2f: rents' log RMSE %.4f (groups %.4f), D %.3f (Census %.3f); burn-in acceptance %.3f, gap %.4f, arrivals %.0f a quarter; %.0f s\n", power,
                kappa, x.rmse, x.rmse_groups, x.d, d_census, x.acceptance, x.gap, x.arrivals, x.seconds);
    std::fflush(stdout);
    trials.push_back(x);
    return R;
  };
  double kappa = 0.2, power = 0.5;
  const bool reuse = argc > 4 && std::string(argv[4]) == "reuse";
  if (reuse) {
    // The values chosen by an earlier run, and its trials kept.
    const Json tr = read_json(results + "/trials.json");
    kappa = tr["kappa"].as_number(), power = tr["power"].as_number();
    for (const auto& x : tr["trials"].as_array())
      trials.push_back({x["power"].as_number(), x["kappa"].as_number(), x["rmse"].as_number(), x["rmse_groups"].as_number(), x["d"].as_number(), x["acceptance"].as_number(),
                        x["gap"].as_number(), x["arrivals"].as_number(), x["seconds"].as_number()});
  } else {
    double best = 1e9;
    for (const double k : {0.1, 0.2, 0.4, 0.8, 1.6}) {
      trial(0.5, k, 1);
      if (trials.back().rmse < best) best = trials.back().rmse, kappa = k;
    }
    best = 1e9;
    for (const double e : {0.0, 0.25, 0.35, 0.5, 0.75}) {
      if (e != 0.5) trial(e, kappa, 1);
      const auto& x = e == 0.5 ? *std::find_if(trials.begin(), trials.end(), [&](const Trial& y) { return y.power == 0.5 && y.kappa == kappa; }) : trials.back();
      if (std::fabs(x.d - d_census) < best) best = std::fabs(x.d - d_census), power = e;
    }
  }
  std::printf("chosen: κ = %.2f, π = %.2f\n", kappa, power);
  {
    Json tr = Json::array();
    for (const auto& x : trials)
      tr.push(Json::object({{"power", x.power}, {"kappa", x.kappa}, {"rmse", x.rmse}, {"rmse_groups", x.rmse_groups}, {"d", x.d}, {"acceptance", x.acceptance},
                            {"gap", x.gap}, {"arrivals", x.arrivals}, {"seconds", x.seconds}}));
    save("trials.json", Json::object({{"trials", tr}, {"kappa", kappa}, {"power", power}, {"d_census", d_census}}));
  }

  // The history at the chosen values, four seeds.
  const std::vector<std::uint64_t> seeds = {11, 12, 13, 14};
  std::vector<Run> runs(seeds.size());
  {
    const auto burn = burn_in(D, city, C, base_params(D, city, C, power, kappa), 1);
    parallel(seeds, [&](std::size_t k, std::uint64_t seed) {
      Run b = burn;
      b.params.dyn.adjustment = kappa;
      runs[k] = history(D, city, C, b, seed);
    });
    std::printf("history: %zu seeds, %.0f s so far\n", seeds.size(), seconds_since(t0));
  }
  {
    // Rents.
    Json groups = Json::array();
    std::vector<RentFit> fits;
    for (const auto& R : runs) fits.push_back(rent_fit(D, city, R));
    for (std::size_t g = 0; g < inner_groups.size(); g++) {
      Json models = Json::array();
      for (const auto& F : fits) models.push(vec(F.model[g]));
      groups.push(Json::object({{"group", inner_groups[g]}, {"observed", vec(fits[0].observed[g])}, {"model", models}, {"weight", fits[0].weights[g]}}));
    }
    Json centre_m = Json::array();
    Json rmse = Json::array();
    for (const auto& F : fits) centre_m.push(vec(F.centre_model)), rmse.push(Json::array({F.rmse, F.rmse_groups}));
    // Student apartments and houses of the model, and the segments' rents by tick (the first seed).
    Json seg = Json::array();
    for (std::size_t j = 0; j < city.segments(); j++) {
      std::vector<double> r;
      for (const auto& x : runs[0].rents) r.push_back(x[j]);
      seg.push(Json::object({{"type", static_cast<double>(j / city.areas)}, {"area", D.areas[j % city.areas]}, {"rent", vec(r)}}));
    }
    Json areas = Json::array();
    for (const auto& a : D.areas) areas.push(a);
    Json ticksq = Json::array();
    for (int t = 0; t < ticks; t++) ticksq.push(quarter_label(quarter_of_tick(t)));
    // The record of moves.
    Json rec = Json::array();
    for (const auto& r : runs[0].record)
      rec.push(Json::array({static_cast<double>(r.arrivals), static_cast<double>(r.entered), static_cast<double>(r.refused), static_cast<double>(r.moved),
                            static_cast<double>(r.left), static_cast<double>(r.evicted), static_cast<double>(r.searched)}));
    save("history.json", Json::object({{"ticks", ticksq}, {"groups", groups}, {"centre_model", centre_m}, {"centre_observed", vec(fits[0].centre_observed)}, {"rmse", rmse},
                                        {"segments", seg}, {"areas", areas}, {"record_columns", Json::array({"arrivals", "entered", "refused", "moved", "left", "evicted", "searched"})},
                                        {"record", rec}, {"kappa", kappa}, {"power", power}, {"arrivals0", runs[0].params.arrivals}, {"metro", vec(C.metro)},
                                        {"reference_by_tick", vec(runs[0].params.reference_by_tick)},
                                        {"demand_by_tick", vec(runs[0].params.demand_by_tick)}, {"base", vec(runs[0].params.pref.base)},
                                        {"taste_owner_house", runs[0].params.pref.taste_owner[uh::house]}}));
  }
  {
    // The Census of 2021 against the model's 2021Q3, by SA1 (each seed) and by small area.
    Json rows = Json::array();
    std::vector<BySa1> ms;
    for (const auto& R : runs) ms.push_back(by_sa1(D, city, R.census, sa1_of, census_tick));
    Json dm = Json::array();
    for (const auto& m : ms) dm.push(model_d(m));
    for (std::size_t s = 0; s < D.sa1.size(); s++) {
      const auto& a = D.sa1[s];
      const auto c = income_classes(D, a);
      const auto& c37 = D.census_columns_g37;
      const auto g = [](double v) { return std::isfinite(v) ? v : 0.0; };
      Json model = Json::array();
      for (const auto& m : ms)
        model.push(Json::array({m.income[s][0], m.income[s][1], m.income[s][2], m.renters[s], m.owners[s], m.moved[s], median(m.rent[s])}));
      Json ring = Json::array();
      for (const auto& p : a.rings[0]) ring.push(Json::array({p[0], p[1]}));
      rows.push(Json::object({{"sa1", a.code}, {"area", D.areas[a.area]}, {"x", a.x}, {"y", a.y}, {"ring", ring},
                              {"census", Json::array({c[0], c[1], c[2], g(a.g37[column(c37, "R_Tot_Total")]), g(a.g37[column(c37, "O_OR_Total")]) + g(a.g37[column(c37, "O_MTG_Total")]),
                                                      g(a.g44[column(D.census_columns_g44, "Difnt_Usl_add_1_yr_ago_Tot_P")]), a.g02[column(D.census_columns_g02, "Median_rent_weekly")],
                                                      g(a.g36[column(D.census_columns_g36, "Unoccupied_PDs_Dwgs")]), g(a.g36[column(D.census_columns_g36, "OPDs_Tot_OPDs_Dwellings")])})},
                              {"model", model}}));
    }
    // Vacancy at 2021Q3 by seed.
    Json vac = Json::array();
    for (const auto& R : runs) {
      const auto cnt = uh::segment_counts(city, R.census);
      double cap = 0, occ = 0;
      for (std::size_t j = 0; j < cnt.capacity.size(); j++) cap += cnt.capacity[j], occ += cnt.occupied[j];
      vac.push(Json::array({cap, occ}));
    }
    save("census-2021.json", Json::object({{"columns_census", Json::array({"low", "middle", "high", "renters", "owners", "persons_moved", "median_rent", "unoccupied", "occupied"})},
                                            {"columns_model", Json::array({"low", "middle", "high", "renters", "owners", "moved_in_year", "median_rent"})},
                                            {"rows", rows}, {"d_census", d_census}, {"d_model", dm}, {"capacity_occupied", vac}}));
  }
  {
    // The stock and its places, for the maps.
    Json b = Json::array();
    for (const auto& r : D.buildings) {
      if (r.development) continue;
      b.push(vec({r.x, r.y, static_cast<double>(r.area), static_cast<double>(r.zone), static_cast<double>(r.stock[0][0]), static_cast<double>(r.stock[1][0]),
                  static_cast<double>(r.stock[2][0]), static_cast<double>(r.stock[0][years - 1]), static_cast<double>(r.stock[1][years - 1]),
                  static_cast<double>(r.stock[2][years - 1])}));
    }
    Json blocks = Json::array();
    for (const auto& k : D.blocks) {
      Json ring = Json::array();
      for (const auto& q : k.rings[0]) ring.push(Json::array({q[0], q[1]}));
      blocks.push(Json::object({{"block", static_cast<double>(k.id)}, {"area", static_cast<double>(k.area)}, {"x", k.x}, {"y", k.y}, {"ring", ring}}));
    }
    Json totals = Json::array();
    for (int y = 0; y < years; y++) {
      std::vector<double> by(D.areas.size() * uh::types, 0.0);
      for (const auto& r : D.buildings)
        if (!r.development)
          for (std::size_t t = 0; t < uh::types; t++) by[t * D.areas.size() + r.area] += r.stock[t][static_cast<std::size_t>(y)];
      totals.push(vec(by));
    }
    Json areas = Json::array();
    for (const auto& a : D.areas) areas.push(a);
    const auto an = anchor();
    save("stock.json", Json::object({{"columns", Json::array({"x", "y", "area", "zone", "house_2002", "apartment_2002", "student_2002", "house_2024", "apartment_2024", "student_2024"})},
                                      {"buildings", b}, {"blocks", blocks}, {"areas", areas}, {"totals_by_year_type_area", totals}, {"anchor", vec({an[0], an[1]})}}));
  }

  // ---- (1c) the frequencies of rents: the report's and the model's quarterly growth ----
  {
    const auto obs = rent_fit(D, city, runs[0]);
    Json spectra = Json::array();
    const auto spectrum_of = [&](const std::vector<double>& index) {
      std::vector<double> g;
      for (int t = 4; t < ticks; t++) {
        const auto u = static_cast<std::size_t>(t);
        if (std::isfinite(index[u]) && std::isfinite(index[u - 1])) g.push_back(std::log(index[u] / index[u - 1]));
      }
      const auto M = calx::signal::multitaper(g, 4.0, 3.0, 5);  // cycles a year
      return Json::object({{"frequency", vec(M.frequency)}, {"density", vec(M.density)}, {"n", static_cast<double>(g.size())}});
    };
    spectra.push(Json::object({{"name", "the report, the centre's flats"}, {"spectrum", spectrum_of(obs.centre_observed)}}));
    for (std::size_t k = 0; k < runs.size(); k++) spectra.push(Json::object({{"name", "the model, seed " + std::to_string(seeds[k])}, {"spectrum", spectrum_of(rent_fit(D, city, runs[k]).centre_model)}}));
    save("rent-spectra.json", Json::object({{"spectra", spectra}, {"note", "multitaper (NW 3, 5 tapers) of the quarterly growth of the log index, 2003Q3-2025Q3; frequency in cycles a year"}}));
  }

  const auto wants = [&](const char* name) { return task == "all" || task.find(name) != std::string::npos; };
  if (wants("shock")) {
  // ---- (2) supply shocks, 2025Q4 to 2030Q4 ----
  struct Scenario {
    std::string name;
    bool towers_uc, towers_approved, students;
  };
  const std::vector<Scenario> scenarios = {{"baseline", false, false, false},
                                           {"towers under construction", true, false, false},
                                           {"every approved tower", true, true, false},
                                           {"student apartments", false, false, true}};
  const int horizon = 114;  // ticks 93 … 113: 2025Q4 … 2030Q4
  // Demand and the metropolis's rent beyond the report: their mean growth of 2010-2019 (an assumption).
  double g_demand = 0, g_metro = 0;
  {
    const auto q = [](int y) { return static_cast<std::size_t>((y - 2000) * 4 + 1); };  // each year's June quarter
    g_demand = std::log(C.demand[q(2019)] / C.demand[q(2010)]) / 9, g_metro = std::log(C.metro[q(2019)] / C.metro[q(2010)]) / 9;
  }
  std::printf("beyond 2025Q3: demand grows %.4f and the metropolis's rents %.4f a year (their means of 2010-2019)\n", g_demand, g_metro);
  const auto extra_at = [&](const Scenario& S, int t) {
    std::vector<std::array<double, uh::types>> x(D.buildings.size(), {0.0, 0.0, 0.0});
    const double time = quarter_time(quarter_of_tick(t));
    const auto ramp = [&](double a, double b) { return std::clamp((time - a) / (b - a), 0.0, 1.0); };
    for (const auto& d : D.developments) {
      double share = 0, apartments = 0, students = 0;
      if (d.status == "COMPLETED") share = ramp(d.year, d.year + 1), apartments = d.dwellings, students = d.student;
      else if (d.status == "UNDER CONSTRUCTION") share = ramp(2026, 2028), apartments = S.towers_uc ? d.dwellings : 0, students = S.students ? d.student : 0;
      else if (d.status == "APPROVED") share = ramp(2027, 2030), apartments = S.towers_approved ? d.dwellings : 0, students = S.students ? d.student : 0;
      x[d.building][uh::apartment] += share * apartments, x[d.building][uh::student] += share * students;
    }
    return x;
  };
  std::vector<int> is_new(D.buildings.size(), 0);
  for (const auto& d : D.developments)
    if (d.status != "COMPLETED") is_new[d.building] = 1;
  struct Future {
    Run run;
    std::vector<std::array<double, classes>> residents;  // by tick from 92
    std::vector<std::array<double, classes>> blocks;     // at the end, by block
    std::array<double, classes> in_new{};                // residents of the pipeline's buildings at the end
    std::vector<uh::TickRecord> record;
  };
  std::vector<std::vector<Future>> futures(scenarios.size(), std::vector<Future>(runs.size()));
  const auto count = [&](const uh::State& st) {
    std::array<double, classes> c{};
    for (const auto& h : st.households)
      if (h.resident) c[h.income]++;
    return c;
  };
  for (std::size_t s = 0; s < scenarios.size(); s++) {
    parallel(seeds, [&](std::size_t k, std::uint64_t seed) {
      Future F;
      F.run = runs[k];
      const auto& P = F.run.params;
      const uh::Simulation sim(city, P.pref, P.dyn);
      auto& st = F.run.state;
      F.residents.push_back(count(st));
      const double d92 = P.demand_by_tick[ticks - 1], m92 = P.reference_by_tick[ticks - 1];
      for (int t = ticks; t < horizon; t++) {
        const double years_on = (t - (ticks - 1)) / 4.0;
        st.reference_rent = m92 * std::exp(g_metro * years_on);
        sim.compact(st);
        sim.begin(st);
        const auto extra = extra_at(scenarios[s], t);
        const auto cap = capacities(D, t, &extra);
        for (std::size_t b = 0; b < cap.size(); b++)
          for (std::size_t q = 0; q < uh::types; q++)
            if (cap[b][q] != st.capacity[b][q]) sim.set_capacity(st, b, q, cap[b][q]);
        sim.arrive(st, newcomers(C, static_cast<std::size_t>(std::lround(P.arrivals * d92 * std::exp(g_demand * years_on))), seed, static_cast<std::uint32_t>(t)));
        sim.run(st, seed);
        F.residents.push_back(count(st));
      }
      F.blocks.assign(D.blocks.size(), {0, 0, 0});
      for (const auto& h : st.households)
        if (h.resident) {
          F.blocks[D.buildings[h.building].zone][h.income]++;
          if (is_new[h.building]) F.in_new[h.income]++;
        }
      F.run.rents = st.rents, F.run.lettings = st.lettings;
      F.record.assign(st.ticks.end() - (horizon - ticks), st.ticks.end());
      futures[s][k] = std::move(F);
    });
    std::printf("scenario %s: %.0f s so far\n", scenarios[s].name.c_str(), seconds_since(t0));
    std::fflush(stdout);
  }
  {
    Json out = Json::array();
    for (std::size_t s = 0; s < scenarios.size(); s++) {
      Json seeds_out = Json::array();
      for (const auto& F : futures[s]) {
        Json groups = Json::array();
        for (const auto& g : inner_groups) {
          std::vector<double> m = lettings_median(D, city, F.run, g, uh::apartment);
          groups.push(vec(std::vector<double>(m.begin() + (ticks - 1), m.end())));
        }
        std::vector<std::vector<double>> by_group;
        std::vector<double> weights;
        for (const auto& g : inner_groups) {
          double lettings = 0;
          (void)rent_index(D, g, flats, &lettings);
          by_group.push_back(lettings_median(D, city, F.run, g, uh::apartment)), weights.push_back(lettings);
        }
        const auto centre_m = centre_index(by_group, weights);
        // The student apartments' rent, geometrically over areas weighted by their 2024 stock (the same weights in every
        // scenario), by tick.
        std::vector<double> student_rent;
        std::vector<double> student_weight(city.areas, 0.0);
        for (const auto& r : D.buildings)
          if (!r.development) student_weight[r.area] += r.stock[uh::student][years - 1];
        for (int t = ticks - 1; t < horizon; t++) {
          double num = 0, den = 0;
          for (std::size_t a = 0; a < city.areas; a++) {
            const std::size_t j = uh::student * city.areas + a;
            const double cap = student_weight[a];
            if (cap > 0) num += cap * std::log(F.run.rents[static_cast<std::size_t>(t)][j]), den += cap;
          }
          student_rent.push_back(den > 0 ? std::exp(num / den) : std::nan(""));
        }
        Json res = Json::array();
        for (const auto& c : F.residents) res.push(vec({c[0], c[1], c[2]}));
        Json blocks = Json::array();
        for (const auto& c : F.blocks) blocks.push(vec({c[0], c[1], c[2]}));
        Json rec = Json::array();
        for (const auto& r : F.record) rec.push(vec({static_cast<double>(r.arrivals), static_cast<double>(r.entered), static_cast<double>(r.refused), static_cast<double>(r.moved),
                                                     static_cast<double>(r.left), static_cast<double>(r.evicted)}));
        seeds_out.push(Json::object({{"centre", vec(std::vector<double>(centre_m.begin() + (ticks - 1), centre_m.end()))}, {"groups", groups}, {"student", vec(student_rent)},
                                     {"residents", res}, {"blocks", blocks}, {"in_new", vec({F.in_new[0], F.in_new[1], F.in_new[2]})}, {"record", rec}}));
      }
      double added_apartments = 0, added_students = 0;
      {
        const auto x = extra_at(scenarios[s], horizon - 1), x0 = extra_at(scenarios[0], horizon - 1);
        for (std::size_t b = 0; b < x.size(); b++) added_apartments += x[b][uh::apartment] - x0[b][uh::apartment], added_students += x[b][uh::student] - x0[b][uh::student];
      }
      out.push(Json::object({{"name", scenarios[s].name}, {"added_apartments", added_apartments}, {"added_students", added_students}, {"seeds", seeds_out}}));
    }
    Json tq = Json::array();
    for (int t = ticks - 1; t < horizon; t++) tq.push(quarter_label(quarter_of_tick(t)));
    Json blocks = Json::array();
    for (const auto& b : D.blocks) {
      Json ring = Json::array();
      for (const auto& p : b.rings[0]) ring.push(Json::array({p[0], p[1]}));
      blocks.push(Json::object({{"block", static_cast<double>(b.id)}, {"area", D.areas[b.area]}, {"x", b.x}, {"y", b.y}, {"ring", ring}}));
    }
    Json devs = Json::array();
    for (const auto& d : D.developments)
      devs.push(Json::object({{"status", d.status}, {"year", static_cast<double>(d.year)}, {"x", d.x}, {"y", d.y}, {"dwellings", d.dwellings}, {"student", d.student},
                              {"block", static_cast<double>(D.blocks[D.buildings[d.building].zone].id)}}));
    save("shocks.json", Json::object({{"ticks", tq}, {"scenarios", out}, {"blocks", blocks}, {"developments", devs}, {"g_demand", g_demand}, {"g_metro", g_metro}}));
  }
  }

  // ---- (3) Schelling's limit ----
  if (wants("schelling")) {
    // (a) On a lattice: the households' limit and nomos's Schelling, tick by tick.
    Json lattices = Json::array();
    for (const double f : {1.0 / 3, 0.5}) {
      const nomos::agents::Lattice L{40, 40, true, nomos::agents::Neighbourhood::moore};
      nomos::agents::Schelling S{L, f};
      auto x = S.place(1200, 5);
      uh::City lc;
      lc.areas = 1, lc.zones = 1;
      std::vector<std::pair<std::size_t, std::size_t>> edges;
      for (std::size_t i = 0; i < L.cells(); i++) {
        lc.buildings.push_back({0, 0, static_cast<double>(i % L.width), static_cast<double>(i / L.width), 0});
        for (const std::size_t m : L.neighbours(i))
          if (m > i) edges.push_back({i, m});
      }
      lc.graph = nomos::agents::Graph(L.cells(), edges);
      uh::Preferences pref;
      pref.base.assign(lc.segments(), 0.0);
      pref.income = {1, 1};
      pref.neighbourhood = 1, pref.step = true, pref.threshold = f;
      pref.scale = std::numeric_limits<double>::infinity(), pref.outside = -std::numeric_limits<double>::infinity();
      uh::Dynamics dyn;
      dyn.forced_by_neighbourhood = true;
      const uh::Simulation sim(lc, pref, dyn);
      auto st = sim.empty(std::vector<std::array<int, uh::types>>(L.cells(), {0, 1, 0}), std::vector<double>(lc.segments(), 1.0));
      for (std::size_t i = 0; i < x.where.size(); i++) {
        uh::Household h;
        h.income = static_cast<std::uint8_t>(x.cell[x.where[i]]), h.tenure = uh::renter, h.size = 1, h.active = true;
        st.households.push_back(h);
        sim.place(st, i, x.where[i], uh::apartment);
      }
      std::vector<double> seg_n, seg_h, moves;
      seg_n.push_back(S.segregation(x)), seg_h.push_back(S.segregation(x));
      Json start = Json(std::vector<double>(x.cell.begin(), x.cell.end()));
      int differ = 0;
      for (std::uint32_t t = 0; t < 40; t++) {
        nomos::agents::tick(S, x, nomos::agents::Schedule::random_sequential, 2026, t);
        sim.step(st, 2026);
        std::vector<int> cells(L.cells(), -1);
        for (const auto& h : st.households) cells[h.building] = h.income;
        differ += cells != x.cell;
        auto y = x;
        std::fill(y.cell.begin(), y.cell.end(), -1);
        for (std::size_t i = 0; i < st.households.size(); i++) y.where[i] = st.households[i].building, y.cell[y.where[i]] = st.households[i].income;
        seg_n.push_back(S.segregation(x)), seg_h.push_back(S.segregation(y)), moves.push_back(st.ticks.back().moved);
      }
      lattices.push(Json::object({{"threshold", f}, {"segregation_nomos", vec(seg_n)}, {"segregation_housing", vec(seg_h)}, {"moves", vec(moves)},
                                  {"ticks_differing", static_cast<double>(differ)}, {"start", start}, {"end", Json(std::vector<double>(x.cell.begin(), x.cell.end()))},
                                  {"width", 40.0}}));
      std::printf("lattice f = %.3f: %d ticks of 40 differ; segregation %.3f → %.3f\n", f, differ, seg_n.front(), seg_n.back());
    }
    // (b) In the city of 2021Q3: sorting by block with prices, types and distance off (the step at f, the best option),
    // and the calibrated model with a taste η for like neighbours.
    const auto& base = runs[0].census;
    Json limits = Json::array();
    for (const double f : {1.0 / 3, 0.5}) {
      uh::Preferences pref;
      pref.base.assign(city.segments(), 0.0);
      pref.income.assign(classes, 1.0);
      pref.taste_owner = {0, 0, -std::numeric_limits<double>::infinity()};
      pref.neighbourhood = 1, pref.step = true, pref.threshold = f;
      pref.scale = std::numeric_limits<double>::infinity(), pref.outside = -std::numeric_limits<double>::infinity();
      uh::Dynamics dyn;
      dyn.forced_by_neighbourhood = true;
      const uh::Simulation sim(city, pref, dyn);
      auto st = sim.empty(base.capacity, base.rent);
      for (const auto& h0 : base.households)
        if (h0.resident) {
          auto h = h0;
          h.resident = false;
          st.households.push_back(h);
          sim.place(st, st.households.size() - 1, h0.building, h0.type);
        }
      std::vector<double> d, content, moved;
      const auto measure = [&] {
        d.push_back(model_d(by_sa1(D, city, st, sa1_of, 0)));
        double c = 0, n = 0;
        for (const auto& h : st.households)
          if (h.resident) n++, c += sim.like(st, D.buildings[h.building].zone, h.income, true, h.building) >= f;
        content.push_back(c / n);
      };
      measure();
      for (std::uint32_t t = 0; t < 30; t++) {
        sim.step(st, 99);
        measure();
        moved.push_back(st.ticks.back().moved);
        if (st.ticks.back().moved == 0) break;
      }
      std::vector<std::array<double, classes>> blocks(D.blocks.size(), {0, 0, 0});
      for (const auto& h : st.households)
        if (h.resident) blocks[D.buildings[h.building].zone][h.income]++;
      Json bl = Json::array();
      for (const auto& c : blocks) bl.push(vec({c[0], c[1], c[2]}));
      limits.push(Json::object({{"threshold", f}, {"d", vec(d)}, {"content", vec(content)}, {"moved", vec(moved)}, {"blocks", bl}}));
      std::printf("the city's Schelling limit at f = %.3f: D %.3f → %.3f in %zu ticks, content %.3f → %.3f\n", f, d.front(), d.back(), moved.size(), content.front(), content.back());
    }
    // The calibrated model's blocks at 2021Q3 (seed 1) and the Census's SA1s, for the maps.
    std::vector<std::array<double, classes>> blocks(D.blocks.size(), {0, 0, 0});
    for (const auto& h : base.households)
      if (h.resident) blocks[D.buildings[h.building].zone][h.income]++;
    Json bl = Json::array();
    for (const auto& c : blocks) bl.push(vec({c[0], c[1], c[2]}));
    // A taste for like neighbours added to the calibrated model: D at 2021Q3 as η grows.
    Json etas = Json::array();
    for (const double e : {0.0, 2.0, 4.0, 8.0}) {
      auto P = base_params(D, city, C, power, kappa);
      P.pref.neighbourhood = e;
      const auto b = burn_in(D, city, C, P, 1);
      Run bb = b;
      bb.params.dyn.adjustment = kappa;
      const auto R = history(D, city, C, bb, 11);
      const double dm = model_d(by_sa1(D, city, R.census, sa1_of, census_tick));
      const auto F = rent_fit(D, city, R);
      etas.push(Json::object({{"eta", e}, {"d", dm}, {"rmse", F.rmse}}));
      std::printf("  η %.1f: D %.3f, rents' log RMSE %.4f\n", e, dm, F.rmse);
      std::fflush(stdout);
    }
    save("schelling.json", Json::object({{"lattices", lattices}, {"city_limits", limits}, {"model_blocks", bl}, {"etas", etas}, {"d_census", d_census},
                                         {"d_model", model_d(by_sa1(D, city, base, sa1_of, census_tick))}}));
  }
  // ---- (4) the planner's permits under uncertain demand ----
  if (wants("planner")) {
    // The market of the planner: the centre's apartments as one segment, the Census's three classes, α and π as
    // calibrated; its stock in units of 10 % of 2024's apartments (about two years of the centre's completions), its
    // demand in the same units at the rent of 2025Q3.
    double S0 = 0;
    for (const auto& b : D.buildings)
      if (!b.development) S0 += b.stock[uh::apartment][years - 1];
    const double unit = 0.10 * S0;
    const std::size_t X = 7, Y = 11, A = 3, y_now = 3;
    const auto obs = rent_fit(D, city, runs[0]);
    const double r_ref = obs.centre_observed[ticks - 1], metro_ref = C.metro[static_cast<std::size_t>(quarter_of_tick(ticks - 1))];
    uh::Market m;
    m.segments = {{uh::apartment, 0, S0, 0.0, 0.0}};
    const auto& pref = runs[0].params.pref;
    for (std::size_t c = 0; c < classes; c++) m.classes.push_back({C.class_share[c], C.census.median_income[c], {0, 0, 0}});
    m.alpha = pref.alpha, m.mean_income = pref.mean_income, m.income_power = pref.income_power, m.vacancy = 0.05, m.reference_rent = metro_ref;
    const double O_ref = m.occupied(0), M_ref = 20 * O_ref;  // one in twenty of those who consider the centre live there (as the agents)
    m.households = M_ref;
    m.segments[0].base = uh::invert(m, std::vector<double>{r_ref}, std::vector<double>{O_ref}).base[0];
    uh::PlannerSpec P;
    P.supply = X, P.demand = Y, P.actions = A, P.bands = 3;
    P.welfare.assign(X * Y, 0.0), P.rent.assign(X * Y, 0.0);
    Json table = Json::array();
    for (std::size_t x = 0; x < X; x++)
      for (std::size_t y = 0; y < Y; y++) {
        auto mm = m;
        mm.segments[0].supply = S0 + unit * static_cast<double>(x);
        mm.households = M_ref * (S0 + unit * (static_cast<double>(y) - static_cast<double>(y_now))) / S0;
        const auto cl = uh::clear(mm, {r_ref});
        const double r = cl.rents[0];
        const double welfare = 52 * (uh::surplus_dollars(mm, cl.rents) + r * mm.occupied(0)) / 1e6;  // A$ million a year
        P.rent[x * Y + y] = r, P.welfare[x * Y + y] = welfare;
        table.push(vec({static_cast<double>(x), static_cast<double>(y), r, welfare}));
      }
    // The demand's transitions: the yearly change of the demand index (June to June, 2003-2025) in units, rounded and
    // held within ±2.
    std::vector<double> steps(5, 0.0);
    std::vector<double> changes;
    for (int yr = 2004; yr <= 2025; yr++) {
      const auto q = [](int y2) { return static_cast<std::size_t>((y2 - 2000) * 4 + 1); };
      if (!std::isfinite(C.demand[q(yr)]) || !std::isfinite(C.demand[q(yr - 1)])) continue;
      // A change of demand by k units is a change of ln M by ln((S0 + k u)/S0) ≈ k u / S0.
      const double k = std::log(C.demand[q(yr)] / C.demand[q(yr - 1)]) / std::log(1 + unit / S0);
      changes.push_back(k);
      steps[static_cast<std::size_t>(std::clamp<long>(std::lround(k), -2, 2) + 2)] += 1;
    }
    const double nyears = static_cast<double>(changes.size());
    for (auto& v : steps) v /= nyears;
    P.demand_transition.assign(Y * Y, 0.0);
    for (std::size_t y = 0; y < Y; y++)
      for (int k = -2; k <= 2; k++) {
        const long y2 = std::clamp<long>(static_cast<long>(y) + k, 0, static_cast<long>(Y) - 1);
        P.demand_transition[y * Y + static_cast<std::size_t>(y2)] += steps[static_cast<std::size_t>(k + 2)];
      }
    P.cuts = {0.99 * r_ref, 1.01 * r_ref};  // rents below, near and above 2025Q3's, within 1 % (an assumption)
    P.confusion = 0.15;                    // a reading's error each side (an assumption)
    // A$ million a unit permitted: A$410,000 a dwelling, and a discount of 0.93 a year (assumptions, chosen so that a unit
    // pays its way at some demands and not at others, a problem whose answer depends on what is hidden: at A$380,000 the
    // MDP permits at every demand, at A$430,000 almost never).
    P.cost = unit * 410000 / 1e6;
    P.discount = 0.93;
    P.x0 = 0, P.y0.assign(Y, 0.0), P.y0[y_now - 1] = 0.25, P.y0[y_now] = 0.5, P.y0[y_now + 1] = 0.25;
    const auto M = uh::planner(P);
    auto Q = td::to_pomdp(M);
    // Rewards less their least (values shift by a constant, policies do not).
    const double rmin = *std::min_element(Q.R.begin(), Q.R.end());
    for (auto& r : Q.R) r -= rmin;
    td::validate(Q);
    const auto t4 = std::chrono::steady_clock::now();
    const auto H = td::hsvi(Q, {.epsilon = 5.0, .max_trials = 4000});
    std::printf("planner: %zu states, HSVI [%.1f, %.1f] (A$ million, shifted by %.1f a year) after %d trials, %.0f s\n", Q.states, H.lower_value, H.upper_value, rmin, H.trials,
                seconds_since(t4));
    const auto mdp = td::solve_mdp(Q);
    const auto qm = td::qmdp(Q);
    using Policy = std::function<std::size_t(std::span<const double>)>;
    std::vector<std::pair<std::string, Policy>> policies = {
        {"never permit", [](std::span<const double>) { return std::size_t{0}; }},
        {"one unit a year", [](std::span<const double>) { return std::size_t{1}; }},
        {"follow the rents", [&](std::span<const double> b) { return uh::policy_follow_rents(P, uh::split(M, b)); }},
        {"certainty equivalent", [&](std::span<const double> b) { return uh::policy_certainty_equivalent(mdp, Y, uh::split(M, b)); }},
        {"QMDP", [&](std::span<const double> b) { return qm.action(b); }},
        {"HSVI", [&](std::span<const double> b) { return H.lower.action(b); }},
    };
    Json pol = Json::array();
    for (const auto& [name, pi] : policies) {
      const auto E = td::evaluate(Q, pi, 3000, 150, 4242);
      // The rents each policy makes: the true state's clearing rent, a year at a time, over 400 episodes of 120
      // years, and its spectrum (the periodogram of ln r after 20 years, averaged).
      std::vector<double> mean_density;
      std::vector<double> freq;
      std::vector<double> path;   // one episode's rents, for the figure
      std::vector<double> actions(A, 0.0);
      for (int e = 0; e < 400; e++) {
        calx::random::Stream rng(777, 0x63796c65u, static_cast<std::uint32_t>(e), 0);  // "cyle"
        std::size_t s = calx::random::categorical(rng, Q.initial);
        td::Belief b = Q.initial;
        std::vector<double> lr;
        for (int t = 0; t < 120; t++) {
          const std::size_t a = pi(b);
          actions[a]++;
          const std::size_t s2 = calx::random::categorical(rng, std::span<const double>(Q.T).subspan((a * Q.states + s) * Q.states, Q.states));
          const std::size_t o = calx::random::categorical(rng, std::span<const double>(Q.Z).subspan((a * Q.states + s2) * Q.observations, Q.observations));
          b = td::update(Q, b, a, o).belief;
          s = s2;
          lr.push_back(std::log(P.rent[s]));
        }
        if (e == 0) path = lr;
        const std::vector<double> tail(lr.begin() + 20, lr.end());
        const auto sp = calx::signal::periodogram(tail, 1.0);
        if (mean_density.empty()) mean_density.assign(sp.density.size(), 0.0), freq = sp.frequency;
        for (std::size_t j = 0; j < sp.density.size(); j++) mean_density[j] += sp.density[j] / 400;
      }
      for (auto& v : actions) v /= 400.0 * 120;
      pol.push(Json::object({{"name", name}, {"value", E.mean + rmin / (1 - Q.discount)}, {"error", E.error}, {"truncation", E.truncation}, {"frequency", vec(freq)},
                             {"density", vec(mean_density)}, {"path", vec(path)}, {"actions", vec(actions)}}));
      std::printf("  %-22s %.1f ± %.1f (A$ million); permits 0/1/2: %.2f %.2f %.2f\n", name.c_str(), E.mean + rmin / (1 - Q.discount), E.error, actions[0], actions[1], actions[2]);
      std::fflush(stdout);
    }
    // The plan's action at each supply for beliefs over demand concentrated at each level (and spread), from HSVI.
    Json plan = Json::array();
    for (std::size_t x = 0; x < X; x++) {
      Json row = Json::array();
      for (std::size_t y = 0; y < Y; y++) {
        td::Belief b(Q.states, 0.0);
        b[x * Y + y] = 1;
        row.push(static_cast<double>(H.lower.action(b)));
      }
      plan.push(row);
    }
    save("planner.json", Json::object({{"unit", unit}, {"S0", S0}, {"r_ref", r_ref}, {"M_ref", M_ref}, {"y_now", static_cast<double>(y_now)}, {"table", table},
                                       {"steps", vec(steps)}, {"changes", vec(changes)}, {"cuts", vec(P.cuts)}, {"cost", P.cost}, {"discount", P.discount},
                                       {"hsvi", vec({H.lower_value + rmin / (1 - Q.discount), H.upper_value + rmin / (1 - Q.discount), static_cast<double>(H.trials),
                                                     H.converged ? 1.0 : 0.0})},
                                       {"policies", pol}, {"plan", plan}, {"X", static_cast<double>(X)}, {"Y", static_cast<double>(Y)}}));
  }
  std::printf("done in %.0f s\n", seconds_since(t0));
  return 0;
}
