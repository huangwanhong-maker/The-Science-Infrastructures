// The engine's operations on representations (Street Oceanography §24): the bridal store and the ring store
// as a catalytic pathway, the crowd on a street as a kinematic wave, with a rule that opens its carriageway,
// and movement in a district of blocks with its circulation and transit. Their parameters are declared by
// the request and returned with the result.

#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <urban/affinity/catalysis.hpp>
#include <urban/flow/budget.hpp>
#include <urban/flow/kinematic.hpp>
#include <urban/flow/movement.hpp>

#include "data.hpp"
#include "ops_street.hpp"

namespace so {

inline Json op_pathway(const Json& r) {
  urban::affinity::Pathway p;
  p.length = num(r, "length", 600), p.spacing = num(r, "spacing", 10), p.diffusivity = num(r, "diffusivity", 1.0);
  p.search_speed = num(r, "search", 0), p.leave_rate = 1.0 / (60 * num(r, "leaveMinutes", 30));
  p.store_residence = 60 * num(r, "storeMinutes", 15), p.catalysis = 1.0 / (60 * num(r, "catalysisMinutes", 20));
  p.purchase = 1.0 / (60 * num(r, "purchaseMinutes", 30)), p.bridal = num(r, "bridal", 300), p.ring = p.bridal + p.spacing;
  if (p.length / p.spacing > 120) throw std::invalid_argument("at most 120 frontages");
  std::vector<double> sep;
  for (double d = p.spacing; d <= std::min(p.bridal, p.length - p.bridal) + 1e-9; d += p.spacing) sep.push_back(d);
  const auto c = urban::affinity::coupling(p, sep);
  urban::affinity::Pathway walk = p;
  walk.search_speed = 0;
  const auto c0 = urban::affinity::coupling(walk, sep);
  std::vector<double> offsets;
  for (double o = -6; o <= 6 + 1e-9; o += 1) offsets.push_back(o);
  urban::affinity::Pathway near = p;
  near.ring = p.bridal + num(r, "windowSeparation", 30);
  const auto w = urban::affinity::coupling_window(near, num(r, "bridalOpen", 11), num(r, "bridalClose", 19), num(r, "ringOpen", 11), num(r, "ringClose", 19), offsets, 300, 0.5);
  const auto cap = urban::affinity::capture(c.adjacent, num(r, "pBridal", 0.3), num(r, "vBridal", 2000), num(r, "vRing", 3000));
  return Json::object({
      {"parameters", Json::object({{"length", p.length}, {"spacing", p.spacing}, {"diffusivity", p.diffusivity}, {"search", p.search_speed},
                                   {"leaveMinutes", 1 / p.leave_rate / 60}, {"storeMinutes", p.store_residence / 60}, {"catalysisMinutes", 1 / p.catalysis / 60},
                                   {"purchaseMinutes", 1 / p.purchase / 60}, {"bridal", p.bridal}})},
      {"coupling", Json::object({{"separation", numbers(c.separation, 1)}, {"yield", numbers(c.yield, 6)}, {"adjacent", c.adjacent},
                                 {"length", c.coupling_length}, {"half", c.half_length}, {"dwellMean", c.dwell_mean}, {"dwellMedian", c.dwell_median},
                                 {"dwellP90", c.dwell_p90}})},
      {"walk", Json::object({{"yield", numbers(c0.yield, 6)}, {"adjacent", c0.adjacent}, {"length", c0.coupling_length}, {"half", c0.half_length}})},
      {"scales", Json::object({{"walkClosedForm", urban::affinity::walk_coupling_length(p.diffusivity, p.leave_rate, p.spacing)},
                               {"diffusive", std::sqrt(p.diffusivity / p.leave_rate)}, {"ballistic", p.search_speed / p.leave_rate}})},
      {"window", Json::object({{"offset", numbers(w.offset, 2)}, {"yield", numbers(w.yield, 6)}, {"bridalOpen", w.b_open}, {"bridalClose", w.b_close},
                               {"ringOpen", w.j_open}, {"ringClose", w.j_close}, {"separation", near.ring - near.bridal}})},
      {"capture", Json::object({{"own", cap.own}, {"created", cap.created}, {"ratio", cap.ratio}})}});
}

// A street of its own from a request: blockLengths (m, each ending at a cross street), signalled (per crossing,
// 1 or 0), demand (24 values: walkers entering at the street's start each second, hour by hour of the day).
// Without them, the declared street of `blocks` blocks of blockLength and a declared day's shape times peak.
inline void street_of(const Json& r, urban::flow::StreetDay& d) {
  if (r["blockLengths"].is_array()) {
    d.block_lengths.clear();
    for (const auto& v : r["blockLengths"].as_array()) d.block_lengths.push_back(v.as_number());
    if (d.block_lengths.empty() || d.block_lengths.size() > 24) throw std::invalid_argument("1 to 24 blocks");
  }
  if (r["signalled"].is_array()) {
    d.signalled.clear();
    for (const auto& v : r["signalled"].as_array()) d.signalled.push_back(v.as_number() != 0 ? 1 : 0);
  }
  if (r["demand"].is_array()) {
    d.demand.clear();
    for (const auto& v : r["demand"].as_array()) d.demand.push_back(std::max(0.0, v.is_number() ? v.as_number() : 0.0));
    if (d.demand.size() != 24) throw std::invalid_argument("a demand for each hour of the day");
  }
}

inline Json op_streetday(const Json& r) {
  urban::flow::StreetDay d;
  d.blocks = static_cast<std::size_t>(num(r, "blocks", 8));
  d.block_length = num(r, "blockLength", 125), d.footway = num(r, "footway", 12), d.carriageway = num(r, "carriageway", 27);
  d.cycle = num(r, "cycle", 120), d.green = num(r, "green", 0.5);
  d.rule_from = 12, d.rule_to = num(r, "ruleTo", 18), d.from = num(r, "from", 10), d.to = num(r, "to", 20), d.snapshot = num(r, "snapshot", 300);
  d.cell = num(r, "cell", 12.5);
  // The demand (people per second entering at Kyōbashi), declared: a weekday or a weekend afternoon.
  const double peak = num(r, "peak", 6);
  d.demand.assign(24, 0);
  const double shape[24] = {0, 0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.2, 0.35, 0.55, 0.75, 0.9, 1.0, 1.0, 0.95, 0.8, 0.6, 0.4, 0.25, 0.1, 0.05, 0};
  for (int h = 0; h < 24; h++) d.demand[static_cast<std::size_t>(h)] = peak * shape[h];
  if (d.blocks < 1 || d.blocks > 12) throw std::invalid_argument("1 to 12 blocks");
  street_of(r, d);
  auto to_json = [&](const urban::flow::StreetDayResult& o) {
    Json snaps = Json::array();
    for (const auto& row : o.density) snaps.push(numbers(row, 4));
    return Json::object({{"times", numbers(o.times, 0)}, {"density", snaps}, {"entered", numbers(o.entered, 1)}, {"left", numbers(o.left, 1)},
                         {"budgetResidual", o.budget_residual}, {"maxDensity", o.max_density}, {"maxDensityRuleHours", o.max_density_in_rule_hours},
                         {"meanTime", o.mean_time}, {"queue", o.queue_spill},
                         {"throughput", o.throughput}});
  };
  d.rule = false;
  const auto ordinary = urban::flow::street_day(d);
  d.rule = true;
  const auto opened = urban::flow::street_day(d);
  return Json::object({{"position", numbers(ordinary.position, 1)}, {"ordinary", to_json(ordinary)}, {"rule", to_json(opened)}, {"demand", numbers(d.demand, 3)},
                       {"parameters", Json::object({{"blocks", static_cast<double>(d.blocks)}, {"blockLength", d.block_length}, {"footway", d.footway},
                                                    {"carriageway", d.carriageway}, {"cycle", d.cycle}, {"green", d.green}, {"peak", peak},
                                                    {"ruleFrom", d.rule_from}, {"ruleTo", d.rule_to}, {"jam", d.jam}, {"freeSpeed", d.free_speed}})}});
}

inline Json op_district(const Json& r) {
  urban::flow::District d;
  d.nx = static_cast<std::size_t>(num(r, "nx", 5)), d.ny = static_cast<std::size_t>(num(r, "ny", 4));
  if (d.nx < 1 || d.ny < 1 || d.nx * d.ny > 64) throw std::invalid_argument("a district of 1 to 64 blocks");
  if (r["holes"].is_array())
    for (const auto& h : r["holes"].as_array()) d.holes.push_back(static_cast<std::size_t>(h.as_number()));
  const std::size_t V = (d.nx + 1) * (d.ny + 1);
  // An attraction: a station at an intersection, falling off with walking distance (grid steps).
  const double station = num(r, "attraction", 2.0);
  const auto at = static_cast<std::size_t>(num(r, "station", 0));
  d.attraction.assign(V, 0);
  for (std::size_t v = 0; v < V; v++) {
    const double di = std::fabs(static_cast<double>(v % (d.nx + 1)) - static_cast<double>(at % (d.nx + 1)));
    const double dj = std::fabs(static_cast<double>(v / (d.nx + 1)) - static_cast<double>(at / (d.nx + 1)));
    d.attraction[v] = station * std::exp(-(di + dj) / 2);
  }
  d.block_circulation.assign(d.nx * d.ny, num(r, "circulation", 0.0));
  for (auto h : d.holes) (void)h, d.hole_circulation.push_back(num(r, "holeCirculation", 0.0));
  d.exits = {0, V - 1};
  d.entry = at;
  const auto c = urban::flow::circulation(d);
  const auto t = urban::flow::transit(d, 3600);
  Json edges = Json::array();
  for (std::size_t e = 0; e < c.complex.edges.size(); e++)
    edges.push(Json::object({{"tail", static_cast<double>(c.complex.edges[e][0])}, {"head", static_cast<double>(c.complex.edges[e][1])}, {"phi", c.phi[e]},
                             {"exact", c.hodge.exact[e]}, {"coexact", c.hodge.coexact[e]}, {"harmonic", c.hodge.harmonic[e]}, {"current", c.current[e]}}));
  Json cycles = Json::array();
  for (const auto& cy : c.irreversibility.cycles) cycles.push(Json::object({{"affinity", cy.affinity}, {"current", cy.current}}));
  return Json::object({{"nx", static_cast<double>(d.nx)}, {"ny", static_cast<double>(d.ny)}, {"holes", numbers(std::vector<double>(d.holes.begin(), d.holes.end()), 0)},
                       {"edges", edges}, {"potential", numbers(c.hodge.potential)}, {"pi", numbers(c.irreversibility.pi)},
                       {"entropyProduction", c.irreversibility.entropy_production}, {"detailedBalance", c.irreversibility.detailed_balance},
                       {"norms", Json::object({{"exact", c.hodge.exact_norm2}, {"coexact", c.hodge.coexact_norm2}, {"harmonic", c.hodge.harmonic_norm2}})},
                       {"affinityGapWithoutAttraction", c.attraction_free_affinity_gap}, {"cycles", cycles},
                       {"transit", Json::object({{"time", numbers(t.time, 0)}, {"density", numbers(t.density, 9)}, {"survival", numbers(t.survival, 6)},
                                                 {"mean", t.mean}, {"median", t.median}, {"p90", t.p90}})}});
}

inline Json op_person_time(const Json& r) {
  // Little's identity on stays drawn from the transit law above is a theorem; here, on the visits the
  // request lists (arrive, depart in seconds) over a window.
  std::vector<double> a, d;
  for (const auto& v : r["stays"].as_array()) a.push_back(v[0].as_number()), d.push_back(v[1].as_number());
  const auto p = urban::flow::person_time(a, d, num(r, "window", 3600));
  return Json::object({{"integral", p.integral}, {"sumOfStays", p.sum_of_stays}, {"arrivalRate", p.arrival_rate}, {"meanPresent", p.mean_present},
                       {"meanStay", p.mean_stay}, {"littleGap", p.little_gap}, {"units", static_cast<double>(p.units)}});
}

}  // namespace so
