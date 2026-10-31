// The engine's operations of the street read as a sea, which the map draws:
//
//   tides                every outdoor sensor's tide: its readings (high and low water, range, form, the year's
//                        swing) and its predicted mean week, from the calendars' constituents fitted to the counts
//                        since October 2024 (the domain's rhythm/tide.hpp)
//   sea {date}           a civil date: each sensor's predicted level through the day (quarter-hours), its observed
//                        counts and its counts each way if the date is in the record, and the day's holiday and weather
//   currents             every sensor that counts each way: its mean week each way, the net current, its
//                        irreversibility, flood and ebb (the domain's flow/currents.hpp)
//   streetwave {…}       the crowd on a street of its own blocks and signals, entering at its measured rate, in a
//                        fine window about an hour: densities, flows, the walkers' speeds, the signals (the
//                        domain's flow/kinematic.hpp; a representation)
//   plume {…}            a visitor's law after leaving a store that raises a need, along its street (the domain's
//                        affinity/catalysis.hpp; a representation)

#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <urban/affinity/catalysis.hpp>
#include <urban/flow/currents.hpp>
#include <urban/flow/kinematic.hpp>
#include <urban/rhythm/tide.hpp>

#include "data.hpp"
#include "ops_models.hpp"
#include "ops_street.hpp"

namespace so {

namespace ur = urban::rhythm;

// The tides are fitted from the live record's start; a sensor needs eight weeks of it.
inline double tide_from() { return ur::parse_utc("2024-10-01"); }

struct SensorTide {
  long long id = 0;
  ur::Record recent;   // the record from tide_from() on, kept: the records' cache holds only a few
  ur::Tide tide;
  ur::TideReadings readings;
};

// The part of a record from tide_from() on, copied (the records' cache may let go of the record itself).
inline ur::Record recent_part(const ur::Record& r) {
  ur::Record out;
  out.id = r.id, out.step = r.step, out.offsets = r.offsets;
  std::size_t first = r.size();
  for (std::size_t i = 0; i < r.size(); i++)
    if (r.time(i) >= tide_from()) { first = i; break; }
  out.t0 = first < r.size() ? r.time(first) : r.t0;
  if (first < r.size()) out.y.assign(r.y.begin() + static_cast<std::ptrdiff_t>(first), r.y.end());
  return out;
}

// Every outdoor sensor's tide, fitted on as many threads as the machine has (each fit is independent).
inline const std::vector<SensorTide>& sensor_tides() {
  static std::vector<SensorTide> list;
  static bool done = false;
  static std::mutex lock;
  std::lock_guard<std::mutex> g(lock);
  if (done) return list;
  std::vector<std::pair<long long, ur::Record>> work;
  for (const auto& s : sensors()) {
    if (s.location_type != "Outdoor") continue;
    try {
      auto rec = recent_part(record(s.id));
      std::size_t n = 0;
      for (double v : rec.y) n += std::isfinite(v);
      if (n >= 56 * 24) work.emplace_back(s.id, std::move(rec));
    } catch (...) {
    }
  }
  std::vector<SensorTide> out(work.size());
  std::vector<char> ok(work.size(), 0);
  std::atomic<std::size_t> next{0};
  const std::size_t threads = std::max<std::size_t>(1, std::min<std::size_t>(work.size(), std::thread::hardware_concurrency()));
  std::vector<std::thread> pool;
  for (std::size_t w = 0; w < threads; w++)
    pool.emplace_back([&] {
      for (std::size_t i = next++; i < work.size(); i = next++) {
        try {
          out[i].id = work[i].first;
          out[i].tide = ur::tide(work[i].second, tide_from(), 1e18);
          out[i].readings = ur::readings(out[i].tide);
          out[i].recent = std::move(work[i].second);
          ok[i] = 1;
        } catch (...) {
        }
      }
    });
  for (auto& t : pool) t.join();
  for (std::size_t i = 0; i < out.size(); i++)
    if (ok[i]) list.push_back(std::move(out[i]));
  done = true;
  return list;
}

inline Json day_json(const ur::TideDay& d) {
  return Json::object({{"high", std::round(d.high)}, {"highAt", std::round(d.high_at * 100) / 100}, {"low", std::round(d.low)},
                       {"lowAt", std::round(d.low_at * 100) / 100}, {"mean", std::round(d.mean)}, {"range", std::round(d.range())}});
}

inline Json op_tides() {
  Json list = Json::array();
  for (const auto& t : sensor_tides()) {
    const auto& s = sensor(t.id);
    const auto& R = t.readings;
    std::vector<double> week;
    for (int k = 0; k < 168; k++) week.push_back(ur::level(t.tide, k + 0.5, ur::of_the_week));
    list.push(Json::object({{"id", static_cast<double>(t.id)}, {"name", s.name}, {"lat", s.lat}, {"lon", s.lon}, {"directions", strings(s.directions)},
                            {"from", ur::format_date(t.tide.first)}, {"to", ur::format_date(t.tide.last)}, {"hours", static_cast<double>(t.tide.observations)},
                            {"constituents", static_cast<double>(t.tide.constituents.size())}, {"varianceExplained", R.variance_explained},
                            {"mean", R.mean}, {"weekday", day_json(R.weekday)}, {"weekend", day_json(R.weekend)},
                            {"diurnal", R.diurnal}, {"semidiurnal", R.semidiurnal}, {"form", std::isfinite(R.form) ? Json(R.form) : Json()},
                            {"formClass", R.form_class}, {"annual", R.annual}, {"week", numbers(week, 1)}}));
  }
  return list;
}

inline long long civil_day_of(const std::string& date) {
  int y = 0, m = 0, d = 0;
  if (std::sscanf(date.c_str(), "%d-%d-%d", &y, &m, &d) != 3) throw std::invalid_argument("a date YYYY-MM-DD");
  return ur::days_from_civil(y, static_cast<unsigned>(m), static_cast<unsigned>(d));
}

inline Json op_sea(const Json& r) {
  const std::string date = r["date"].text();
  const long long day = civil_day_of(date);
  Json list = Json::array();
  std::size_t observed_sensors = 0;
  for (const auto& t : sensor_tides()) {
    const auto& rec = t.recent;
    const auto curve = ur::day_curve(t.tide, day, 4);
    const auto obs = ur::observed_day(rec, day);
    std::vector<double> o(obs.begin(), obs.end()), w1(24, std::nan("")), w2(24, std::nan(""));
    bool any = false;
    for (double v : o) any = any || std::isfinite(v);
    observed_sensors += any;
    try {
      const auto& d = directions(t.id);
      const auto a = ur::observed_day(d.way[0], day), b = ur::observed_day(d.way[1], day);
      w1.assign(a.begin(), a.end()), w2.assign(b.begin(), b.end());
    } catch (...) {
    }
    list.push(Json::object({{"id", static_cast<double>(t.id)}, {"level", numbers(curve, 1)}, {"observed", numbers(o, 0)},
                            {"way1", numbers(w1, 0)}, {"way2", numbers(w2, 0)}}));
  }
  // The day's holiday and weather (Melbourne's civil day, from the first sensor's offsets).
  std::string holiday;
  try {
    for (const auto& h : holidays())
      if (h.date == date) holiday = h.name;
  } catch (...) {
  }
  Json weather_day = Json();
  try {
    const auto& w = weather();
    const auto& rec = record(sensors().front().id);
    double rain = 0, tmax = -1e9, tmin = 1e9;
    std::size_t n = 0;
    const double start = static_cast<double>(day) * 86400 - 12 * 3600;
    for (long long k = static_cast<long long>(std::floor((start - w.t0) / w.step)); k < static_cast<long long>(std::floor((start - w.t0) / w.step)) + 40; k++) {
      if (k < 0 || k >= static_cast<long long>(w.temperature.size())) continue;
      const double t = w.t0 + w.step * static_cast<double>(k);
      if (rec.civil_day(t) != day) continue;
      const auto i = static_cast<std::size_t>(k);
      if (std::isfinite(w.precipitation[i])) rain += w.precipitation[i];
      if (std::isfinite(w.temperature[i])) tmax = std::max(tmax, w.temperature[i]), tmin = std::min(tmin, w.temperature[i]);
      n++;
    }
    if (n >= 20) weather_day = Json::object({{"rain", std::round(rain * 10) / 10}, {"tmax", std::round(tmax * 10) / 10}, {"tmin", std::round(tmin * 10) / 10}});
  } catch (...) {
  }
  return Json::object({{"date", date}, {"sensors", list}, {"observed", static_cast<double>(observed_sensors)}, {"holiday", holiday}, {"weather", weather_day}});
}

// A sensor's side of its street, from its name ("Town Hall (West)"), and a compass direction by its name.
inline bool compass_of(const std::string& word, urban::flow::Compass& c) {
  using C = urban::flow::Compass;
  if (word == "North") return c = C::north, true;
  if (word == "South") return c = C::south, true;
  if (word == "East") return c = C::east, true;
  if (word == "West") return c = C::west, true;
  return false;
}
inline bool side_of(const std::string& name, urban::flow::Compass& c) {
  for (const char* w : {"North", "South", "East", "West"})
    if (name.find(std::string("(") + w + ")") != std::string::npos) return compass_of(w, c);
  return false;
}

inline Json op_currents() {
  static Json cached;
  static bool done = false;
  if (done) return cached;
  Json list = Json::array();
  std::vector<urban::flow::Compass> sides;
  std::vector<std::array<urban::flow::Compass, 2>> towards;
  std::vector<double> residuals;
  Json footpaths = Json::array();
  for (const auto& s : sensors()) {
    if (s.location_type != "Outdoor" || s.directions.size() != 2) continue;
    try {
      const auto& d = directions(s.id);
      const auto m = urban::flow::current_meter(d.way[0], d.way[1]);
      std::size_t weeks = 1000000;
      for (double w : m.weeks) weeks = std::min(weeks, static_cast<std::size_t>(w));
      if (weeks < 8) continue;
      list.push(Json::object({{"id", static_cast<double>(s.id)}, {"name", s.name}, {"lat", s.lat}, {"lon", s.lon}, {"toward", strings({d.toward[0], d.toward[1]})},
                              {"weeks", static_cast<double>(weeks)}, {"way1", numbers(m.way[0], 1)}, {"way2", numbers(m.way[1], 1)},
                              {"net", numbers(m.net, 1)}, {"irreversibility", numbers(m.irreversibility, 2)}, {"residual", m.residual}, {"transport", m.transport},
                              {"sigma", m.sigma}, {"floodAt", m.flood_at}, {"flood", m.flood}, {"ebbAt", m.ebb_at}, {"ebb", m.ebb}}));
      urban::flow::Compass side, a, b;
      if (side_of(s.name, side) && compass_of(d.toward[0], a) && compass_of(d.toward[1], b)) {
        sides.push_back(side), towards.push_back({a, b}), residuals.push_back(m.residual);
        footpaths.push(static_cast<double>(s.id));
      }
    } catch (...) {
    }
  }
  const auto h = urban::flow::handedness(sides, towards, residuals);
  cached = Json::object({{"sensors", list}, {"handedness", Json::object({{"footpaths", static_cast<double>(h.footpaths)}, {"left", static_cast<double>(h.left)},
                                                                         {"pLeft", h.p_left}, {"pTwoSided", h.p_two_sided}, {"named", footpaths}})}});
  done = true;
  return cached;
}

// ---- the crowd on a street, second by second -------------------------------------------------------------------

inline Json op_streetwave(const Json& r) {
  urban::flow::StreetDay d;
  d.blocks = static_cast<std::size_t>(num(r, "blocks", 8));
  d.block_length = num(r, "blockLength", 125), d.footway = num(r, "footway", 12), d.carriageway = num(r, "carriageway", 27);
  d.cycle = num(r, "cycle", 120), d.green = num(r, "green", 0.5);
  d.rule = r["rule"].is_bool() && r["rule"].as_bool();
  d.rule_from = num(r, "ruleFrom", 12), d.rule_to = num(r, "ruleTo", 18), d.snapshot = 600, d.cell = num(r, "cell", 12.5);
  const double peak = num(r, "peak", 8);
  d.demand.assign(24, 0);
  const double shape[24] = {0, 0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.2, 0.35, 0.55, 0.75, 0.9, 1.0, 1.0, 0.95, 0.8, 0.6, 0.4, 0.25, 0.1, 0.05, 0};
  for (int h = 0; h < 24; h++) d.demand[static_cast<std::size_t>(h)] = peak * shape[h];
  const double span = std::clamp(num(r, "window", 240), 30.0, 600.0);
  d.window_step = std::clamp(num(r, "step", 2), 0.5, 10.0);
  if (d.blocks < 1 || d.blocks > 12) throw std::invalid_argument("1 to 12 blocks");
  street_of(r, d);
  // From the morning, when the street is empty, to the window's end.
  d.from = std::clamp(num(r, "from", 6), 0.0, 23.0);
  const double at = std::clamp(num(r, "hour", 15), d.from, 24.0 - span / 3600);
  d.to = at + span / 3600 + 1e-6;
  d.window_from = at * 3600, d.window_to = at * 3600 + span;
  const auto o = urban::flow::street_day(d);
  Json density = Json::array(), speed = Json::array(), flow = Json::array(), green = Json::array();
  for (std::size_t k = 0; k < o.window_times.size(); k++) {
    density.push(numbers(o.window_density[k], 4));
    speed.push(numbers(o.window_speed[k], 3));
    flow.push(numbers(o.window_flow[k], 4));
    Json g = Json::array();
    for (char c : o.window_green[k]) g.push(static_cast<double>(c));
    green.push(g);
  }
  return Json::object({{"position", numbers(o.position, 2)}, {"interfaces", numbers(o.interfaces, 2)}, {"crossings", numbers(o.crossings, 2)},
                       {"crossingLength", o.crossing_length}, {"length", o.interfaces.empty() ? 0.0 : o.interfaces.back()},
                       {"times", numbers(o.window_times, 1)}, {"density", density}, {"speed", speed}, {"flow", flow}, {"width", numbers(o.window_width, 1)},
                       {"green", green}, {"rule", d.rule}, {"freeSpeed", d.free_speed}, {"jam", d.jam}, {"footway", d.footway}, {"carriageway", d.carriageway},
                       {"budgetResidual", o.budget_residual}, {"entered", o.entered.empty() ? 0.0 : o.entered.back()}, {"cycle", d.cycle},
                       {"signalled", numbers([&] { std::vector<double> v; for (std::size_t c = 0; c < o.crossings.size(); c++) v.push_back(d.signalled.empty() || d.signalled[c] ? 1 : 0); return v; }(), 0)}});
}

// ---- a visitor's plume ------------------------------------------------------------------------------------------

inline Json op_plume(const Json& r) {
  urban::affinity::Pathway p;
  p.length = num(r, "length", 600);
  p.spacing = std::max(num(r, "spacing", 10), p.length / 120);
  p.diffusivity = num(r, "diffusivity", 1.0), p.search_speed = num(r, "search", 0.15);
  p.leave_rate = 1.0 / (60 * num(r, "leaveMinutes", 30)), p.store_residence = 60 * num(r, "storeMinutes", 15);
  p.catalysis = 1.0 / (60 * num(r, "catalysisMinutes", 20)), p.purchase = 1.0 / (60 * num(r, "purchaseMinutes", 30));
  p.bridal = std::clamp(num(r, "bridal", 300), 0.0, p.length), p.ring = std::clamp(num(r, "ring", 330), 0.0, p.length);
  if (std::lround(p.bridal / p.spacing) == std::lround(p.ring / p.spacing)) throw std::invalid_argument("the two stores at one frontage");
  const double horizon = std::clamp(num(r, "horizon", 3600), 60.0, 6 * 3600.0);
  const auto steps = static_cast<std::size_t>(std::clamp(num(r, "steps", 60), 1.0, 240.0));
  const auto pl = urban::affinity::pathway_plume(p, horizon, steps);
  Json errand = Json::array(), need = Json::array();
  for (std::size_t k = 0; k < pl.times.size(); k++) errand.push(numbers(pl.errand[k], 6)), need.push(numbers(pl.need[k], 6));
  return Json::object({{"times", numbers(pl.times, 1)}, {"position", numbers(pl.position, 2)}, {"errand", errand}, {"need", need},
                       {"bought", numbers(pl.bought, 6)}, {"left", numbers(pl.left, 6)}, {"yield", urban::affinity::pathway_yield(p)},
                       {"spacing", p.spacing}, {"bridal", p.bridal}, {"ring", p.ring}, {"length", p.length}});
}

}  // namespace so
