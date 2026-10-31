// The engine's operations on terrain, hours and affinity: Melbourne's sensors and their surroundings
// (CLUE), and OpenStreetMap's named streets with their declared hours.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <urban/affinity/colocation.hpp>
#include <urban/affinity/spacing.hpp>
#include <urban/street/hours.hpp>
#include <urban/street/mixing.hpp>
#include <urban/street/terrain.hpp>
#include <urban/street/visits.hpp>

#include "data.hpp"

namespace so {

inline Json numbers(const std::vector<double>& v, int digits = 6) {
  Json a = Json::array();
  const double scale = std::pow(10.0, digits);
  for (double x : v) a.push(std::isfinite(x) ? Json(std::round(x * scale) / scale) : Json());
  return a;
}

inline Json strings(const std::vector<std::string>& v) {
  Json a = Json::array();
  for (const auto& s : v) a.push(s);
  return a;
}

inline Json mixing_json(const urban::street::StreetMixing& m) {
  return Json::object({{"place", m.place}, {"stretch", m.stretch}, {"places", static_cast<double>(m.counts.rows)},
                       {"entropy", m.mixing.entropy}, {"conditional", m.mixing.conditional}, {"departure", m.mixing.departure},
                       {"normalised", m.mixing.normalised}, {"between", m.scales.between}, {"within", m.scales.within},
                       {"local", numbers(m.mixing.local)}});
}

inline Json correlation_json(const urban::street::UseCorrelation& c, double walk) {
  Json pairs = Json::array();
  for (auto p : c.pairs) pairs.push(static_cast<double>(p));
  return Json::object({{"r", numbers(c.r, 2)}, {"rho", numbers(c.rho)}, {"pairs", pairs}, {"pInfinity", c.p_infinity},
                       {"length", c.length}, {"firstZero", c.first_zero}, {"reachedZero", c.reached_zero}, {"singleRegime", c.single_regime},
                       {"walk", walk}, {"compression", urban::street::compression(walk, c)}});
}

// ---- the catalogue -----------------------------------------------------------------------------------------

inline Json op_status() {
  Json out = Json::object();
  auto have = [&](const char* name) {
    try { read_json(name); return true; } catch (...) { return false; }
  };
  for (const char* p : {"street-sensors.json", "street-weather-melbourne.json", "street-holidays-victoria.json", "street-terrain-melbourne.json", "street-holdings.json"})
    out.set(p, have(p));
  return out;
}

inline Json op_sensors() {
  Json list = Json::array();
  std::set<long long> with_terrain;
  try {
    for (const auto& [id, v] : terrain().by_sensor) if (!v.empty()) with_terrain.insert(id);
  } catch (...) {}
  for (const auto& s : sensors())
    list.push(Json::object({{"id", static_cast<double>(s.id)}, {"code", s.code}, {"name", s.name}, {"lat", s.lat}, {"lon", s.lon},
                            {"locationType", s.location_type}, {"first", s.first}, {"last", s.last}, {"hours", static_cast<double>(s.hours)},
                            {"missing", static_cast<double>(s.missing)}, {"terrain", with_terrain.count(s.id) > 0}}));
  return list;
}

inline Json op_streets() {
  Json list = Json::array();
  for (const auto& n : streets()) {
    std::size_t declared = 0;
    for (const auto& h : n.street.holdings) declared += h.declared;
    list.push(Json::object({{"id", n.street.id}, {"name", n.street.name}, {"city", n.city}, {"area", n.area}, {"length", n.street.length},
                            {"holdings", static_cast<double>(n.street.holdings.size())}, {"declared", static_cast<double>(declared)}, {"note", n.note}}));
  }
  return list;
}

// ---- a named street: holdings, mixing, terrain length, hours ---------------------------------------------------

inline Json op_street(const Json& r) {
  const auto& n = named_street(r["street"].text("melbourne-swanston"));
  const auto& st = n.street;
  const double walk = num(r, "walk", 400);
  Json holdings = Json::array();
  for (std::size_t i = 0; i < st.holdings.size(); i++) {
    const auto& h = st.holdings[i];
    Json hours = Json::array();
    for (const auto& [a, b] : h.hours) hours.push(Json::array({a, b}));
    holdings.push(Json::object({{"id", h.id}, {"name", h.name}, {"s", h.s}, {"offset", h.offset}, {"use", static_cast<double>(h.use)},
                                {"kind", n.kinds[i]}, {"hours", hours}, {"declared", h.declared}, {"raw", n.raw_hours[i]}, {"status", n.status[i]},
                                {"lat", n.places[i][0]}, {"lon", n.places[i][1]}}));
  }
  Json axis = Json::array();
  for (const auto& p : n.axis) axis.push(Json::array({p[0], p[1]}));
  Json out = Json::object({{"id", st.id}, {"name", st.name}, {"city", n.city}, {"timezone", n.timezone}, {"area", n.area}, {"length", st.length},
                           {"uses", strings(st.uses)}, {"holdings", holdings}, {"note", n.note}, {"axis", axis}});
  // Mixing along the street and its scales.
  if (st.holdings.size() >= 4) {
    out.set("mixing", mixing_json(urban::street::street_mixing(st, num(r, "place", 25), num(r, "stretch", 100))));
    std::vector<double> s;
    std::vector<std::size_t> u;
    for (const auto& h : st.holdings) s.push_back(h.s), u.push_back(h.use);
    out.set("correlation", correlation_json(urban::street::use_correlation(s, {}, u, num(r, "bin", 5), std::min(300.0, std::max(50.0, st.length / 2))), walk));
  }
  // Hours: the forcing, its lines, the open fraction, ascent and descent.
  std::size_t declared = 0;
  for (const auto& h : st.holdings) declared += h.declared;
  out.set("declared", static_cast<double>(declared));
  if (declared > 0) {
    out.set("forcing", numbers(urban::street::opening_forcing(st, 168), 3));
    Json by_use = Json::object();
    for (std::size_t u = 0; u < st.uses.size(); u++) {
      bool any = false;
      for (const auto& h : st.holdings) any = any || (h.use == u && h.declared);
      if (any) by_use.set(st.uses[u], numbers(urban::street::opening_forcing(st, 168, {u}), 3));
    }
    out.set("forcingByUse", by_use);
    Json lines = Json::array();
    for (const auto& l : urban::street::week_lines(st, 42))
      lines.push(Json::object({{"m", static_cast<double>(l.m)}, {"frequency", l.frequency}, {"amplitude", l.amplitude},
                               {"phase", std::arg(l.value)}, {"daily", l.daily}}));
    out.set("lines", lines);
    std::vector<double> sgrid, theta;
    const double ds = std::max(5.0, st.length / 120);
    for (double x = 0; x <= st.length; x += ds) sgrid.push_back(x);
    for (int k = 0; k < 168 * 2; k++) theta.push_back(0.25 + 0.5 * k);
    Json raster = Json::array();
    for (const auto& row : urban::street::open_fraction(st, sgrid, theta, num(r, "kernel", 25))) raster.push(numbers(row, 3));
    out.set("openFraction", Json::object({{"s", numbers(sgrid, 1)}, {"theta", numbers(theta, 2)}, {"values", raster}}));
    Json days = Json::array();
    for (int d = 0; d < 7; d++) {
      const auto a = urban::street::ascent_descent(st, d, 999, 1);
      days.push(Json::object({{"day", d}, {"holdings", static_cast<double>(a.holdings)}, {"meanOpen", a.mean_open}, {"meanClose", a.mean_close},
                              {"sdOpen", a.sd_open}, {"sdClose", a.sd_close}, {"asymmetry", a.asymmetry}, {"tau", a.tau}, {"p", a.p},
                              {"open", numbers(a.open, 3)}, {"close", numbers(a.close, 3)}}));
    }
    out.set("ascentDescent", days);
  }
  return out;
}

// The visits of kinds of people to a named street's holdings (Representation 2.14), with declared types.
inline Json op_visits(const Json& r) {
  const auto& n = named_street(r["street"].text("melbourne-swanston"));
  std::vector<urban::street::ArrivalType> types;
  for (const auto& t : r["types"].as_array()) {
    urban::street::ArrivalType a;
    a.name = t["name"].text();
    for (const auto& v : t["presence"].as_array()) a.presence.push_back(v.as_number());
    a.affinity.assign(n.street.uses.size(), 0.0);
    if (t["affinity"].is_object())
      for (const auto& [use, w] : t["affinity"].as_object())
        for (std::size_t u = 0; u < n.street.uses.size(); u++)
          if (n.street.uses[u] == use) a.affinity[u] = w.as_number();
    a.distance_cost = num(t, "distanceCost", 0.004);
    a.entry = num(t, "entry", 0);
    types.push_back(std::move(a));
  }
  const auto v = urban::street::visits(n.street, types);
  Json byHolding = Json::array(), weekly = Json::array();
  for (std::size_t h = 0; h < n.street.holdings.size(); h++) {
    byHolding.push(Json::object({{"totals", numbers(v.by_holding[h], 2)}, {"share", numbers(v.share[h], 4)}, {"offDaily", v.off_daily[h]}}));
  }
  // The street's weekly visits by type (summed over holdings), for drawing.
  for (std::size_t a = 0; a < types.size(); a++) {
    std::vector<double> w(168, 0);
    for (std::size_t h = 0; h < n.street.holdings.size(); h++)
      for (std::size_t k = 0; k < 168; k++) w[k] += v.visits[a][h][k];
    weekly.push(numbers(w, 2));
  }
  return Json::object({{"types", strings(v.types)}, {"holdings", byHolding}, {"weekly", weekly}, {"unserved", numbers(v.unserved, 2)}});
}

// ---- a sensor's surroundings (CLUE) --------------------------------------------------------------------------------

inline std::vector<Establishment> about_sensor(long long id, int year) {
  const auto& t = terrain();
  const auto it = t.by_sensor.find(id);
  if (it == t.by_sensor.end()) throw std::invalid_argument("no surroundings of sensor " + std::to_string(id) + " in CLUE");
  std::vector<Establishment> out;
  for (const auto& e : it->second)
    if (e.year == year) out.push_back(e);
  return out;
}

inline Json op_terrain(const Json& r) {
  const auto& t = terrain();
  const long long id = static_cast<long long>(num(r, "sensor", 1));
  const int year = static_cast<int>(num(r, "year", latest_year(t)));
  const double walk = num(r, "walk", 400);
  const auto es = about_sensor(id, year);
  std::vector<double> x, y;
  std::vector<std::size_t> u;
  std::vector<double> composition(t.uses.size(), 0);
  Json points = Json::array();
  for (const auto& e : es) {
    const std::size_t k = t.use_index(e.use);
    x.push_back(e.x), y.push_back(e.y), u.push_back(k);
    composition[k] += 1;
    points.push(Json::array({std::round(e.x * 10) / 10, std::round(e.y * 10) / 10, static_cast<double>(k)}));
  }
  Json out = Json::object({{"sensor", static_cast<double>(id)}, {"name", sensor(id).name}, {"year", year}, {"uses", strings(t.uses)},
                           {"points", points}, {"composition", numbers(composition, 0)}});
  Json years = Json::array();
  for (int yy : t.years) years.push(yy);
  out.set("years", years);
  if (es.size() >= 8) {
    out.set("mixing", mixing_json(urban::street::plane_mixing(x, y, u, t.uses.size(), num(r, "cell", 25), num(r, "block", 100))));
    out.set("correlation", correlation_json(urban::street::use_correlation(x, y, u, num(r, "bin", 10), 200), walk));
  }
  // The record of the surroundings: M and the composition, year by year (Definition 7.1).
  Json series = Json::array();
  for (int yy : t.years) {
    const auto ey = about_sensor(id, yy);
    if (ey.size() < 8) continue;
    std::vector<double> xx, yv;
    std::vector<std::size_t> uu;
    std::vector<double> comp(t.uses.size(), 0);
    for (const auto& e : ey) xx.push_back(e.x), yv.push_back(e.y), uu.push_back(t.use_index(e.use)), comp[t.use_index(e.use)] += 1;
    const auto m = urban::street::plane_mixing(xx, yv, uu, t.uses.size(), 25, 100);
    series.push(Json::object({{"year", yy}, {"establishments", static_cast<double>(ey.size())}, {"departure", m.mixing.departure},
                              {"normalised", m.mixing.normalised}, {"composition", numbers(comp, 0)}}));
  }
  out.set("record", series);
  return out;
}

// The terrain of every sensor's surroundings as a vector of components in unlike units, its dominance order and
// the rankings of named weightings (Theorem 2.2).
inline Json op_dominance(const Json& r) {
  const auto& t = terrain();
  const int year = static_cast<int>(num(r, "year", latest_year(t)));
  const std::vector<std::string> components{"food and drink (count)", "retail (count)", "offices and professional (count)", "departure from mixing (nats)",
                                            "terrain length L_grad (m)"};
  std::vector<std::vector<double>> T;
  Json items = Json::array();
  for (const auto& s : sensors()) {
    if (s.location_type != "Outdoor" || !t.by_sensor.count(s.id)) continue;
    const auto es = about_sensor(s.id, year);
    if (es.size() < 20) continue;
    double food = 0, retail = 0, office = 0;
    std::vector<double> x, y;
    std::vector<std::size_t> u;
    for (const auto& e : es) {
      food += e.use == "food-and-drink";
      retail += e.use == "clothing-and-fashion" || e.use == "jewellery" || e.use == "bridal-and-wedding" || e.use == "other-retail";
      office += e.use == "offices-and-professional" || e.use == "finance";
      x.push_back(e.x), y.push_back(e.y), u.push_back(t.use_index(e.use));
    }
    const auto m = urban::street::plane_mixing(x, y, u, t.uses.size(), 25, 100);
    const auto c = urban::street::use_correlation(x, y, u, 10, 200);
    T.push_back({food, retail, office, m.mixing.departure, c.length});
    items.push(Json::object({{"id", static_cast<double>(s.id)}, {"name", s.name}, {"terrain", numbers(T.back(), 4)}}));
  }
  Json covers = Json::array();
  for (const auto& c : urban::street::covering(T)) covers.push(Json::array({static_cast<double>(c.below), static_cast<double>(c.above)}));
  Json maximal = Json::array();
  for (auto m : urban::street::maximal(T)) maximal.push(static_cast<double>(m));
  const std::vector<urban::street::Weighting> ws{
      {"a visitor with an afternoon to spend", {1.0, 1.0, 0.0, 10.0, 0.0}},
      {"an office worker at lunch", {1.0, 0.0, 0.5, 0.0, 0.0}},
      {"a buyer of rings", {0.0, 2.0, 0.0, 0.0, 0.05}},
      {"the same buyer, the length in kilometres", {0.0, 2.0, 0.0, 0.0, 50.0}}};
  Json rankings = Json::array();
  std::vector<urban::street::Contraction> cs;
  for (const auto& w : ws) {
    cs.push_back(urban::street::contract(T, w));
    Json rank = Json::array();
    for (auto k : cs.back().rank) rank.push(static_cast<double>(k));
    rankings.push(Json::object({{"name", w.name}, {"w", numbers(w.w)}, {"value", numbers(cs.back().value, 4)}, {"rank", rank}}));
  }
  Json agree = Json::array();
  for (const auto& a : cs) {
    Json row = Json::array();
    for (const auto& b : cs) row.push(urban::street::agreement(a, b));
    agree.push(row);
  }
  return Json::object({{"year", year}, {"components", strings(components)}, {"items", items}, {"covers", covers}, {"maximal", maximal},
                       {"rankings", rankings}, {"agreement", agree}});
}

// Every establishment of a census year about the sensors, once each, as offsets from its sensor (the page
// places them in its own frame from the sensor's latitude and longitude).
inline Json op_establishments(const Json& r) {
  const auto& t = terrain();
  const int year = static_cast<int>(num(r, "year", latest_year(t)));
  std::set<std::tuple<long long, long long, std::string, std::string>> seen;
  const double lat0 = -37.8136, lon0 = 144.9631, R = 6371008.8, deg = std::numbers::pi / 180;
  Json points = Json::array(), places = Json::object();
  for (const auto& s : sensors()) {
    if (!t.by_sensor.count(s.id)) continue;
    places.set(std::to_string(s.id), Json::array({s.lat, s.lon}));
    const double ox = (s.lon - lon0) * deg * R * std::cos(lat0 * deg), oy = (s.lat - lat0) * deg * R;
    for (const auto& e : t.by_sensor.at(s.id)) {
      if (e.year != year) continue;
      if (!seen.insert({std::llround(ox + e.x), std::llround(oy + e.y), e.name, e.anzsic}).second) continue;
      points.push(Json::array({static_cast<double>(s.id), std::round(e.x * 10) / 10, std::round(e.y * 10) / 10, static_cast<double>(t.use_index(e.use)), e.name}));
    }
  }
  return Json::object({{"year", year}, {"uses", strings(t.uses)}, {"sensors", places}, {"points", points}});
}

// ---- affinity: co-location and spacing ---------------------------------------------------------------------------

inline Json colocation_json(const urban::affinity::Affinity& a) {
  auto mat = [&](const std::vector<double>& v) { return numbers(v, 4); };
  Json counts = Json::array();
  for (auto c : a.free.count) counts.push(static_cast<double>(c));
  return Json::object({{"uses", strings(a.uses)}, {"count", counts}, {"clq", mat(a.free.clq)}, {"pUpper", mat(a.free.p_upper)}, {"pLower", mat(a.free.p_lower)},
                       {"stratifiedPUpper", mat(a.stratified.p_upper)}, {"stratifiedPLower", mat(a.stratified.p_lower)},
                       {"stratifiedNullMean", mat(a.stratified.null_mean)}, {"cell", a.cell}, {"radius", a.radius},
                       {"permutations", static_cast<double>(a.free.permutations)}});
}

inline Json op_affinity(const Json& r) {
  const std::string where = r["where"].text("melbourne");
  const double cell = num(r, "cell", 100), radius = num(r, "radius", 0);
  const auto perms = static_cast<std::size_t>(num(r, "permutations", 499));
  if (where == "melbourne") {
    const auto& t = terrain();
    const int year = static_cast<int>(num(r, "year", latest_year(t)));
    // Every sensor's surroundings in one plane about the CBD (−37.8136, 144.9631), duplicates dropped.
    const double lat0 = -37.8136, lon0 = 144.9631, R = 6371008.8, deg = std::numbers::pi / 180;
    std::set<std::tuple<long long, long long, std::string, std::string>> seen;
    std::vector<double> x, y;
    std::vector<std::size_t> u;
    for (const auto& s : sensors()) {
      if (!t.by_sensor.count(s.id)) continue;
      const double ox = (s.lon - lon0) * deg * R * std::cos(lat0 * deg), oy = (s.lat - lat0) * deg * R;
      for (const auto& e : t.by_sensor.at(s.id)) {
        if (e.year != year) continue;
        const double X = ox + e.x, Y = oy + e.y;
        if (!seen.insert({std::llround(X), std::llround(Y), e.name, e.anzsic}).second) continue;
        x.push_back(X), y.push_back(Y), u.push_back(t.use_index(e.use));
      }
    }
    auto a = urban::affinity::affinity(x, y, u, t.uses, cell, radius, perms, 1);
    Json out = colocation_json(a);
    out.set("where", "Melbourne: the surroundings of every pedestrian sensor (CLUE, " + std::to_string(year) + ")");
    out.set("points", static_cast<double>(x.size()));
    return out;
  }
  const auto& n = named_street(where);
  std::vector<double> x, y;
  std::vector<std::size_t> u;
  for (const auto& h : n.street.holdings) x.push_back(h.s), y.push_back(h.offset), u.push_back(h.use);
  auto a = urban::affinity::affinity(x, y, u, n.street.uses, cell, radius, perms, 1);
  Json out = colocation_json(a);
  out.set("where", n.street.name + " (OpenStreetMap)");
  out.set("points", static_cast<double>(x.size()));
  return out;
}

inline Json op_spacing(const Json& r) {
  Json out = Json::object();
  const auto p = urban::affinity::predicted_spacing(num(r, "j11", 1), num(r, "j12", 1), num(r, "j21", -2), num(r, "j22", -1.5), num(r, "du", 1), num(r, "dv", 20));
  out.set("predicted", Json::object({{"stable", p.turing.stable_homogeneous}, {"unstable", p.turing.unstable}, {"criticalRatio", p.turing.critical_ratio},
                                     {"k2Low", p.turing.k2_low}, {"k2High", p.turing.k2_high}, {"k2Fastest", p.turing.k2_fastest}, {"wavelength", p.wavelength}}));
  if (r["street"].is_string()) {
    const auto& n = named_street(r["street"].text());
    std::vector<double> s;
    std::vector<std::size_t> u;
    for (const auto& h : n.street.holdings) s.push_back(h.s), u.push_back(h.use);
    std::size_t which = 0;
    const std::string use = r["use"].text("food-and-drink");
    for (std::size_t k = 0; k < n.street.uses.size(); k++)
      if (n.street.uses[k] == use) which = k;
    const auto o = urban::affinity::observed_spacing(s, u, which, n.street.length);
    out.set("observed", Json::object({{"street", n.street.name}, {"use", n.street.uses.empty() ? "" : n.street.uses[which]}, {"wavelength", numbers(o.wavelength, 2)},
                                      {"power", numbers(o.power, 5)}, {"best", o.best}, {"bestPower", o.best_power}, {"falseAlarm", o.false_alarm},
                                      {"holdings", static_cast<double>(o.holdings)}, {"ofUse", static_cast<double>(o.of_use)}}));
  }
  return out;
}

}  // namespace so
