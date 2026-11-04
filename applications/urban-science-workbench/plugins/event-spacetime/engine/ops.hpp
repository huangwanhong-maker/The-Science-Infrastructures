// The operations of the plug-in "Event Spacetime" (the infrastructure's ADR 0067): the urban-science domain's events
// module on a district the page sends (its walking network and its buildings, from the city model) and the record of
// the district's events the data library compiles (city-events-<area>.json, under --data). Every number is computed
// here; the page draws.
//
// The record is assembled from the volumes and interfaces of the district (urban::events::city_network) and the
// dated facts of the product, each attached to the cell it names by identifier (A3 of the design):
//   births        a building's construction in the window begins its coupling, and each of its interfaces begins at
//                 the later of its two sides' births; a building constructed before the window is alive at ∂₀;
//   alterations   a refurbishment, a development completed in a building that stood, touch the building;
//   institutions  a registration on the heritage register touches the building it names;
//   regimes       a street event closes the cross-sections of the street spaces it names, and opens them at its end;
//   observations  each census of a building touches it, a local tick; a census of a building before its
//                 construction (the site's earlier building) is left out and counted.
// An event that the record cannot hold (it touches a cell outside its life) is left out and counted with its reason.

#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <architecture/space/massing.hpp>
#include <science/json.hpp>
#include <urban/events/complex.hpp>
#include <urban/events/memory.hpp>
#include <urban/events/network.hpp>
#include <urban/events/record.hpp>
#include <urban/landscape/district.hpp>

namespace es {

using science::Json;
namespace ue = urban::events;
namespace sm = physis::statistical_mechanics;
namespace ul = urban::landscape;
namespace as = architecture::space;

inline std::string& data_dir() {
  static std::string dir = "data/cache/products";
  return dir;
}

inline double num(const Json& j, const char* key, double fallback) { return j[key].is_number() ? j[key].as_number() : fallback; }
inline Json round(double x, int digits = 3) {
  if (!std::isfinite(x)) return Json();
  const double p = std::pow(10.0, digits);
  return Json(std::round(x * p) / p);
}
inline Json xy(const std::array<double, 2>& p) { return Json::array({round(p[0], 2), round(p[1], 2)}); }
inline std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

// "YYYY-MM-DD" (or "YYYY") as an MJD; NaN when it is not a date.
inline double mjd_of(const Json& v) {
  if (v.is_number()) return ue::mjd(static_cast<int>(v.as_number()));
  if (!v.is_string()) return std::nan("");
  const std::string s = v.as_string();
  int y = 0, m = 1, d = 1;
  if (std::sscanf(s.c_str(), "%d-%d-%d", &y, &m, &d) < 1 || y < 1000) return std::nan("");
  return ue::mjd(y, m, d);
}
inline std::string date_of(double mjd) {
  // The civil date of an MJD (Fliegel and Van Flandern's inverse).
  long l = static_cast<long>(std::floor(mjd)) + 2400001 + 68569;
  const long n = 4 * l / 146097;
  l -= (146097 * n + 3) / 4;
  const long i = 4000 * (l + 1) / 1461001;
  l = l - 1461 * i / 4 + 31;
  const long j = 80 * l / 2447, day = l - 2447 * j / 80;
  l = j / 11;
  const long month = j + 2 - 12 * l, year = 100 * (n - 49) + i + l;
  char buf[16];
  std::snprintf(buf, sizeof buf, "%04ld-%02ld-%02ld", year, month, day);
  return buf;
}

inline Json read_json(const std::string& name) {
  const std::string path = data_dir() + "/" + name;
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("no product " + name + " in the data library's cache: run the workbench's recipe (Build the engines and fetch the data)");
  std::stringstream ss;
  ss << f.rdbuf();
  return Json::parse(ss.str());
}

// A columnar table of the product: {"columns": [...], "rows": [[...]]}.
struct Table {
  std::map<std::string, std::size_t> col;
  const Json* rows = nullptr;
  explicit Table(const Json& t) {
    const auto& c = t["columns"];
    if (c.is_array())
      for (std::size_t k = 0; k < c.as_array().size(); k++) col[c[k].text()] = k;
    rows = &t["rows"];
  }
  std::size_t size() const { return rows && rows->is_array() ? rows->as_array().size() : 0; }
  const Json& at(std::size_t r, const std::string& name) const {
    static const Json none;
    const auto it = col.find(name);
    return it == col.end() ? none : (*rows)[r][it->second];
  }
};

// ---- the district and its buildings -------------------------------------------------------------------------

inline ul::District district_of(const Json& r) {
  ul::Network N;
  for (const auto& j : r["junctions"].as_array()) N.junctions.push_back({j[0].as_number(), j[1].as_number()});
  for (const auto& s : r["segments"].as_array()) {
    ul::Segment g;
    g.u = static_cast<std::size_t>(s["u"].as_number()), g.v = static_cast<std::size_t>(s["v"].as_number());
    for (const auto& p : s["line"].as_array()) g.line.push_back({p[0].as_number(), p[1].as_number()});
    g.kind = s["kind"].text("street"), g.name = s["name"].text();
    N.segments.push_back(std::move(g));
  }
  return ul::district(ul::planarised(N));
}

inline std::vector<as::Footprint> footprints_of(const Json& r) {
  std::vector<as::Footprint> F;
  for (const auto& b : r["buildings"].as_array()) {
    as::Footprint f;
    f.key = b["key"].text(), f.name = b["name"].text();
    for (const auto& q : b["ring"].as_array()) f.ring.push_back({q[0].as_number(), q[1].as_number()});
    // Counter-clockwise, as the massing expects.
    double s = 0;
    for (std::size_t i = 0; i < f.ring.size(); i++) {
      const auto& p = f.ring[i];
      const auto& q = f.ring[(i + 1) % f.ring.size()];
      s += p[0] * q[1] - q[0] * p[1];
    }
    if (s < 0) std::reverse(f.ring.begin(), f.ring.end());
    f.base = num(b, "base", 0), f.top = num(b, "top", 10);
    F.push_back(std::move(f));
  }
  if (F.size() > 4000) throw std::invalid_argument("at most 4000 buildings");
  return F;
}

// ---- the record ------------------------------------------------------------------------------------------------

struct Built {
  std::string key;
  ue::CityNetwork network;
  ue::Record record;
  ue::EventComplex complex;
  std::map<std::string, std::size_t> dropped;  // reason ↦ events left out
  std::map<std::string, std::size_t> kinds;    // type/kind ↦ events kept
  std::size_t massing_skipped = 0;
  Json sources = Json::array();
  // The panels of uses (stage E4): each base property's use at each census; each coupling's use by its property seen
  // at most censuses; each premises' class where one establishment holds it; the ANZSIC codes.
  std::vector<ue::Observation> building_panel, premises_panel;
  std::map<std::size_t, std::map<int, std::string>> use_at;
  std::vector<ue::AnzsicCode> anzsic;
};

inline std::shared_ptr<Built>& last() {
  static std::shared_ptr<Built> b;
  return b;
}

// The request without its operation and time: what determines the record built from it.
inline std::string key_of(const Json& r) {
  Json k = Json::object();
  for (const auto& [name, v] : r.as_object())
    if (name != "op" && name != "t" && name != "id") k.set(name, v);
  return k.dump();
}

inline std::shared_ptr<Built> build(const Json& r) {
  const std::string area = r["area"].text("melbourne");
  const auto D = district_of(r);
  const auto F = footprints_of(r);
  const auto M = as::massing(F);
  auto B = std::make_shared<Built>();
  B->massing_skipped = M.skipped.size();
  B->network = ue::city_network(M, F, D);
  const auto& N = B->network;
  // The record given with the request (the checks give a small one), else the data library's for the area.
  const Json P = r["record"].is_object() ? r["record"] : read_json("city-events-" + area + ".json");
  B->sources = P["sources"];
  auto& R = B->record;
  R.window = {mjd_of(P["window"]["from"]), mjd_of(P["window"]["to"])};
  if (!std::isfinite(R.window.from) || !std::isfinite(R.window.to)) throw std::runtime_error("the record's window is not a pair of dates");
  R.couplings = N.couplings, R.primitives = N.primitives;
  // The properties' buildings.
  std::map<std::string, std::size_t> of_property;
  {
    const Table T(P["properties"]);
    for (std::size_t i = 0; i < T.size(); i++) {
      const auto b = T.at(i, "building").text();
      const auto it = N.building.find(b);
      if (!b.empty() && it != N.building.end()) of_property[T.at(i, "id").text()] = it->second;
    }
  }
  const auto coupling_of = [&](const Json& match) -> long {
    if (match["property"].is_string()) {
      const auto it = of_property.find(match["property"].as_string());
      return it == of_property.end() ? -1 : static_cast<long>(it->second);
    }
    if (match["building"].is_string()) {
      const auto it = N.building.find(match["building"].as_string());
      return it == N.building.end() ? -1 : static_cast<long>(it->second);
    }
    return -1;
  };
  std::map<std::string, std::vector<std::size_t>> streets;  // a street's name, lower case ↦ its street spaces
  for (std::size_t c = 0; c < N.couplings.size(); c++)
    if (N.couplings[c].id.rfind("street:", 0) == 0 && !N.couplings[c].name.empty()) streets[lower(N.couplings[c].name)].push_back(c);
  std::map<std::size_t, std::vector<std::size_t>> cross;  // a street space ↦ its cross-sections
  for (std::size_t p = 0; p < N.primitives.size(); p++)
    if (N.primitives[p].type == "cross-section") cross[N.primitives[p].a].push_back(p);
  // Births first: each building's construction in the window, the earliest one recorded.
  const Table E(P["events"]);
  std::map<std::size_t, std::size_t> birth;  // coupling ↦ its birth event in R.events
  const auto drop = [&](const std::string& why) { B->dropped[why]++; };
  const auto in_window = [&](double t) { return t >= R.window.from && t <= R.window.to; };
  std::vector<std::size_t> later;  // rows of the events table to read after the births
  for (std::size_t i = 0; i < E.size(); i++) {
    const std::string type = E.at(i, "type").text(), kind = E.at(i, "kind").text();
    if (type != "birth") { later.push_back(i); continue; }
    const double from = mjd_of(E.at(i, "from")), to = mjd_of(E.at(i, "to"));
    const long c = coupling_of(E.at(i, "match"));
    if (c < 0) { drop("names no building of the district"); continue; }
    if (!std::isfinite(from)) { drop("no date"); continue; }
    if (!in_window(from)) continue;  // before the window: the building is alive at ∂₀
    const auto cu = static_cast<std::size_t>(c);
    if (birth.count(cu)) {
      // A second account of the same building's making (a development completed): kept as an alteration.
      later.push_back(i);
      continue;
    }
    ue::Event e;
    e.id = E.at(i, "id").text("birth " + std::to_string(i)), e.type = "birth", e.kind = kind;
    e.when = {from, std::isfinite(to) ? to : from};
    e.begins_couplings = {cu};
    e.source = E.at(i, "source").text(), e.provenance = E.at(i, "provenance").text("recorded");
    birth[cu] = R.events.size();
    R.events.push_back(e);
  }
  // Each interface begins with the later of its two sides.
  for (std::size_t p = 0; p < N.primitives.size(); p++) {
    long at = -1;
    for (auto c : {N.primitives[p].a, N.primitives[p].b}) {
      const auto it = birth.find(c);
      if (it == birth.end()) continue;
      if (at < 0 || R.events[it->second].when.from > R.events[static_cast<std::size_t>(at)].when.from ||
          (R.events[it->second].when.from == R.events[static_cast<std::size_t>(at)].when.from && static_cast<long>(it->second) > at))
        at = static_cast<long>(it->second);
    }
    if (at >= 0) R.events[static_cast<std::size_t>(at)].begins_primitives.push_back(p);
  }
  const auto born_by = [&](std::size_t c, double t) {
    const auto it = birth.find(c);
    return it == birth.end() || R.events[it->second].when.from <= t;
  };
  for (const std::size_t i : later) {
    const std::string type = E.at(i, "type").text(), kind = E.at(i, "kind").text();
    const double from = mjd_of(E.at(i, "from")), to = mjd_of(E.at(i, "to"));
    if (!std::isfinite(from)) { drop("no date"); continue; }
    if (!in_window(from)) { drop("outside the window"); continue; }
    ue::Event e;
    e.id = E.at(i, "id").text(type + " " + std::to_string(i)), e.kind = kind;
    e.when = {from, std::isfinite(to) && to >= from ? to : from};
    e.source = E.at(i, "source").text(), e.provenance = E.at(i, "provenance").text("recorded");
    const auto& match = E.at(i, "match");
    if (type == "regime") {
      // A street event under a permit: the street spaces it names, by name, taken into the event's use at its start
      // and given back at its end. Where the source says the street was closed to vehicles ("closure": true in the
      // event's detail), the event meets the street spaces' cross-sections, whose permeability it changes; Melbourne's
      // permits do not say, so their events are states of the street spaces alone.
      std::set<std::size_t> named;
      std::vector<std::string> names;
      if (match["street"].is_string()) names.push_back(match["street"].as_string());
      if (match["streets"].is_array())
        for (const auto& s : match["streets"].as_array()) names.push_back(s.text());
      for (const auto& nm : names) {
        const auto it = streets.find(lower(nm));
        if (it != streets.end()) named.insert(it->second.begin(), it->second.end());
      }
      if (named.empty()) { drop("names no street of the district"); continue; }
      e.type = "regime";
      const bool closure = E.at(i, "detail")["closure"].is_bool() && E.at(i, "detail")["closure"].as_bool();
      for (auto c : named) {
        if (closure) for (auto p : cross[c]) e.primitives.push_back(p);
        else e.couplings.push_back(c);
      }
      ue::Event end = e;
      e.when = {from, from};
      R.events.push_back(e);
      if (std::isfinite(to) && to > from && in_window(to)) {
        end.id = e.id + " (end)", end.kind = closure ? "reopened" : "ended", end.when = {to, to};
        R.events.push_back(end);
      }
      continue;
    }
    const long c = coupling_of(match);
    if (c < 0) { drop("names no building of the district"); continue; }
    if (!born_by(static_cast<std::size_t>(c), from)) { drop("before the building's construction"); continue; }
    e.type = type == "birth" ? "alteration" : type;
    e.couplings = {static_cast<std::size_t>(c)};
    R.events.push_back(e);
  }
  // The censuses: one observation of each building at each census year, a local tick.
  {
    const Table C(P["censuses"]);
    std::set<std::pair<std::size_t, int>> seen;
    for (std::size_t i = 0; i < C.size(); i++) {
      const auto it = of_property.find(C.at(i, "property").text());
      if (it == of_property.end()) continue;
      const auto& yv = C.at(i, "year");
      const int year = yv.is_number() ? static_cast<int>(yv.as_number()) : 0;
      if (year <= 0) continue;
      const auto span = ue::year_span(year);
      if (!in_window(span.from)) continue;
      if (!born_by(it->second, span.from)) { drop("a census before the building's construction"); continue; }
      if (!seen.insert({it->second, year}).second) continue;
      ue::Event e;
      e.id = "census " + std::to_string(year) + " " + N.couplings[it->second].id, e.type = "observation", e.kind = "census";
      e.when = span, e.couplings = {it->second}, e.source = "melbourne-clue-buildings";
      R.events.push_back(e);
    }
  }
  // The panels of uses.
  {
    const Table C(P["censuses"]);
    std::map<std::string, std::map<int, std::string>> by_property;
    for (std::size_t i = 0; i < C.size(); i++) {
      const std::string prop = C.at(i, "property").text();
      if (!of_property.count(prop)) continue;
      const auto& yv = C.at(i, "year");
      const std::string use = C.at(i, "predominant_use").text();
      if (!yv.is_number() || use.empty()) continue;
      by_property[prop].emplace(static_cast<int>(yv.as_number()), use);  // a property's first row of a census
    }
    std::map<std::size_t, std::pair<std::size_t, std::string>> best;
    for (const auto& [prop, ys] : by_property) {
      for (const auto& [y, u] : ys) B->building_panel.push_back({prop, u, y});
      auto& b = best[of_property.at(prop)];
      if (ys.size() > b.first) b = {ys.size(), prop};
    }
    for (const auto& [c, b] : best) B->use_at[c] = by_property.at(b.second);
    const Table Es(P["establishments"]);
    std::map<std::pair<std::string, int>, std::vector<std::string>> held;
    for (std::size_t i = 0; i < Es.size(); i++) {
      if (!of_property.count(Es.at(i, "property").text())) continue;
      const auto& pid = Es.at(i, "property_id");
      const auto& yv = Es.at(i, "year");
      if (!yv.is_number()) continue;
      const std::string premises = pid.is_number() ? std::to_string(static_cast<long long>(pid.as_number())) : pid.text();
      held[{premises, static_cast<int>(yv.as_number())}].push_back(Es.at(i, "anzsic4").text());
    }
    for (const auto& [k, codes] : held)
      if (codes.size() == 1 && !codes[0].empty()) B->premises_panel.push_back({"premises " + k.first, codes[0], k.second});
    if (P["anzsic"].is_object())
      for (const auto& [code, info] : P["anzsic"].as_object())
        B->anzsic.push_back({code, info["division"].text(), info["subdivision"].text(), info["group"].text()});
  }
  // What the record cannot hold is left out, with its reason.
  for (int pass = 0; pass < 4; pass++) {
    const auto bad = ue::problems(R);
    if (bad.empty()) break;
    std::set<std::string> ids;
    for (const auto& b : bad)
      if (b.rfind("event ", 0) == 0) ids.insert(b.substr(6, b.find(':') - 6));
    if (ids.empty()) throw std::runtime_error("the record is not sound: " + bad.front());
    const auto before = R.events.size();
    R.events.erase(std::remove_if(R.events.begin(), R.events.end(), [&](const ue::Event& e) { return ids.count(e.id) > 0; }), R.events.end());
    B->dropped["not held by the record (outside a cell's life)"] += before - R.events.size();
  }
  for (const auto& e : R.events) B->kinds[e.type + (e.kind.empty() ? "" : "/" + e.kind)]++;
  B->complex = ue::construct(R);
  B->key = key_of(r);
  return B;
}

// The record of a request: the one last built when the request names the same district and record, or none is
// named (a request of {t} alone); else built afresh.
inline std::shared_ptr<Built> record_of(const Json& r) {
  auto& L = last();
  if (r["junctions"].is_array()) {
    if (!L || L->key != key_of(r)) L = build(r);
  } else if (!L) {
    throw std::runtime_error("no record built yet: ask for the record first, or send the district with the time");
  }
  return L;
}

inline Json topology_json(const ue::Topology& T) {
  Json o = Json::object({{"vertices", static_cast<double>(T.cells[0])}, {"edges", static_cast<double>(T.cells[1])}, {"faces", static_cast<double>(T.cells[2])},
                         {"euler", static_cast<double>(T.euler)}, {"b0", static_cast<double>(T.b0)}});
  if (T.b1) o.set("b1", static_cast<double>(*T.b1));
  if (T.b2) o.set("b2", static_cast<double>(*T.b2));
  return o;
}

// The district's record and its event complex, for drawing: the couplings with their lives and the events they
// attend, the interfaces with their lives and meetings, the events (the censuses counted, not listed), the
// topology of K and K̂, and what was left out.
inline Json op_record(const Json& r) {
  const auto B = record_of(r);
  const auto& R = B->record;
  const auto& X = B->complex;
  const auto& N = B->network;
  // Each coupling's vertices: how many events it attends, and, for a building, the years of the censuses that
  // observed it (the other events are listed with whom they attend).
  std::vector<std::vector<double>> census_years(R.couplings.size());
  for (std::size_t v = 1; v + 1 < X.time.size(); v++) {
    const auto& e = R.events[static_cast<std::size_t>(X.event_of[v])];
    if (e.type == "observation")
      for (auto c : X.encounters.present[v]) census_years[c].push_back(std::stod(date_of(e.when.from).substr(0, 4)));
  }
  Json couplings = Json::array();
  for (std::size_t c = 0; c < R.couplings.size(); c++) {
    const auto& k = R.couplings[c];
    Json years = Json::array();
    for (double y : census_years[c]) years.push(y);
    couplings.push(Json::object({{"id", k.id}, {"type", k.type}, {"name", k.name}, {"at", xy(k.at)}, {"base", round(k.base, 1)}, {"top", round(k.top, 1)},
                                 {"life", Json::array({round(X.time[X.life_c[c][0]], 2), round(X.time[X.life_c[c][1]], 2)})},
                                 {"vertices", static_cast<double>(X.K.worldline[c].size())}, {"censuses", years}}));
  }
  Json primitives = Json::array();
  for (std::size_t p = 0; p < R.primitives.size(); p++) {
    const auto& q = R.primitives[p];
    Json meets = Json::array();
    for (auto v : X.encounters.pairs[p].meetings) meets.push(round(X.time[v], 2));
    primitives.push(Json::object({{"id", q.id}, {"type", q.type}, {"a", static_cast<double>(q.a)}, {"b", static_cast<double>(q.b)}, {"area", round(q.area, 1)},
                                  {"life", Json::array({round(X.time[X.life_p[p][0]], 2), round(X.time[X.life_p[p][1]], 2)})}, {"meets", meets}}));
  }
  Json events = Json::array();
  std::map<int, std::size_t> censuses;
  for (std::size_t v = 1; v + 1 < X.time.size(); v++) {
    const auto& e = R.events[static_cast<std::size_t>(X.event_of[v])];
    if (e.type == "observation") { censuses[std::stoi(date_of(e.when.from).substr(0, 4))]++; continue; }
    Json who = Json::array(), met = Json::array();
    std::array<double, 2> mid{};
    for (auto c : X.encounters.present[v]) who.push(static_cast<double>(c)), mid[0] += R.couplings[c].at[0], mid[1] += R.couplings[c].at[1];
    const double n = static_cast<double>(std::max<std::size_t>(1, X.encounters.present[v].size()));
    mid[0] /= n, mid[1] /= n;
    std::set<std::size_t> m(e.primitives.begin(), e.primitives.end());
    m.insert(e.begins_primitives.begin(), e.begins_primitives.end()), m.insert(e.ends_primitives.begin(), e.ends_primitives.end());
    for (auto p : m) met.push(static_cast<double>(p));
    events.push(Json::object({{"vertex", static_cast<double>(v)}, {"id", e.id}, {"type", e.type}, {"kind", e.kind}, {"from", round(e.when.from, 2)}, {"to", round(e.when.to, 2)},
                              {"date", date_of(e.when.from)}, {"at", xy(mid)}, {"attend", who}, {"meets", met}, {"source", e.source}, {"provenance", e.provenance}}));
  }
  Json census = Json::array();
  for (const auto& [y, n] : censuses) census.push(Json::object({{"year", static_cast<double>(y)}, {"buildings", static_cast<double>(n)}}));
  Json dropped = Json::object(), kinds = Json::object();
  for (const auto& [k, n] : B->dropped) dropped.set(k, static_cast<double>(n));
  for (const auto& [k, n] : B->kinds) kinds.set(k, static_cast<double>(n));
  const auto TK = ue::topology(X.K.cw), TR = ue::topology(X.refined.cw);
  return Json::object({
      {"window", Json::object({{"from", round(R.window.from, 2)}, {"to", round(R.window.to, 2)}, {"fromDate", date_of(R.window.from)}, {"toDate", date_of(R.window.to)}})},
      {"couplings", couplings}, {"primitives", primitives}, {"events", events}, {"censuses", census},
      {"network", Json::object({{"buildings", static_cast<double>(N.buildings)}, {"streets", static_cast<double>(N.streets)}, {"junctions", static_cast<double>(N.junctions)},
                                {"partyWalls", static_cast<double>(N.party_walls)}, {"frontages", static_cast<double>(N.frontages)}, {"crossSections", static_cast<double>(N.cross_sections)},
                                {"facades", static_cast<double>(N.facades)}, {"facadesFronting", static_cast<double>(N.facades_fronting)},
                                {"massingSkipped", static_cast<double>(B->massing_skipped)}})},
      {"topology", Json::object({{"K", topology_json(TK)}, {"refined", topology_json(TR)}, {"triangles", static_cast<double>(X.refined.triangles.size())}})},
      {"kinds", kinds}, {"dropped", dropped}, {"sources", B->sources}});
}

// ---- what the city remembers (stage E4) -------------------------------------------------------------------------

inline Json numbers(const std::vector<double>& v, int digits = 6) {
  Json a = Json::array();
  for (double x : v) a.push(round(x, digits));
  return a;
}

inline Json memory_json(const ue::UseClassification& U, const ue::UseMemory& M) {
  Json uses = Json::array();
  std::vector<std::size_t> order(U.uses.size());
  for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return U.C.law[a] > U.C.law[b]; });
  for (std::size_t k = 0; k < std::min<std::size_t>(order.size(), 12); k++)
    uses.push(Json::object({{"use", U.uses[order[k]]}, {"share", round(U.C.law[order[k]], 4)}, {"family", U.block_names[1][order[k]]}}));
  Json levels = Json::array();
  for (const auto& l : M.levels) levels.push(l);
  Json bound = Json::array();
  for (bool b : M.fit.at_bound) bound.push(b);
  return Json::object({{"levels", levels}, {"lambda", numbers(M.fit.lambda)}, {"lambdaError", numbers(M.fit.lambda_error)}, {"mu", numbers(M.fit.mu)},
                       {"muError", numbers(M.fit.mu_error)}, {"halfLife", numbers(M.half_life, 3)}, {"halfLifeError", numbers(M.half_life_error, 3)},
                       {"agree", numbers(M.agree, 0)}, {"transitions", M.transitions}, {"units", static_cast<double>(M.units)}, {"uses", static_cast<double>(M.uses)},
                       {"logL", round(M.fit.log_likelihood, 3)}, {"flatLogL", round(M.flat_log_likelihood, 3)}, {"ratio", round(M.ratio, 3)},
                       {"degrees", static_cast<double>(M.degrees)}, {"converged", M.fit.converged}, {"atBound", bound}, {"gradient", numbers(M.fit.gradient, 9)}, {"iterations", static_cast<double>(M.fit.iterations)}, {"top", uses}});
}

// What the district remembers of its uses: the buildings' uses (CLUE's space uses, families by CLUE's own prefixes) and
// the premises' classes (ANZSIC, where one establishment holds a premises), fitted level by level; each building's
// spans between censuses with their expected unseen changes and their modes' weights; the buildings' uses as chains on
// the event complex, their probability by the contraction against the product of their transitions; the rates of the
// record's kinds of events.
inline Json op_memory(const Json& r) {
  const auto B = record_of(r);
  const auto& R = B->record;
  const auto& X = B->complex;
  const auto units = [](const std::vector<ue::Observation>& O) { std::set<std::string> u; for (const auto& o : O) u.insert(o.unit); return u.size(); };
  const auto Ub = ue::prefix_classification(ue::use_counts(B->building_panel));
  const auto Mb = ue::remember(Ub, ue::transitions_of(Ub, B->building_panel), units(B->building_panel));
  const sm::NestedChain Xb(Ub.C, Mb.fit.lambda);
  Json premises;
  try {
    const auto Up = ue::anzsic_classification(B->anzsic, ue::use_counts(B->premises_panel));
    const auto Mp = ue::remember(Up, ue::transitions_of(Up, B->premises_panel), units(B->premises_panel));
    premises = memory_json(Up, Mp);
  } catch (const std::exception& e) {
    premises = Json::object({{"refused", e.what()}});
  }
  // Each building's spans between its censuses.
  Json spans = Json::array();
  for (const auto& [c, ys] : B->use_at) {
    std::vector<ue::Observation> O;
    for (const auto& [y, u] : ys) O.push_back({"", u, y});
    Json list = Json::array();
    for (const auto& sp : ue::spans(Ub, Xb, O))
      list.push(Json::object({{"from", static_cast<double>(sp.from)}, {"to", static_cast<double>(sp.to)}, {"first", sp.first}, {"last", sp.last},
                              {"expected", round(sp.expected_changes, 4)}, {"weights", numbers(sp.weights, 4)}}));
    spans.push(Json::object({{"coupling", static_cast<double>(c)}, {"spans", list}}));
  }
  // The buildings' uses as chains on the event complex.
  std::map<std::size_t, std::map<std::size_t, std::size_t>> seen;
  for (std::size_t v = 1; v + 1 < X.time.size(); v++) {
    const auto& e = R.events[static_cast<std::size_t>(X.event_of[v])];
    if (e.type != "observation" || e.couplings.empty()) continue;
    const std::size_t c = e.couplings[0];
    const auto it = B->use_at.find(c);
    if (it == B->use_at.end()) continue;
    const auto u = it->second.find(std::stoi(date_of(e.when.from).substr(0, 4)));
    if (u == it->second.end() || !Ub.index.count(u->second)) continue;
    seen[c][v] = Ub.of(u->second);
  }
  const auto chains = ue::chains_on_complex(X, Xb, seen);
  const auto Z = sm::contract(chains);
  double closed = 0;
  for (const auto& [c, m] : seen) {
    auto prev = m.begin();
    closed += std::log(Ub.C.law[prev->second]);
    for (auto it = std::next(m.begin()); it != m.end(); prev = it, ++it)
      closed += std::log(Xb.kernel(prev->second, it->second, (X.time[it->first] - X.time[prev->first]) / 365.2425));
  }
  Json rates = Json::array();
  for (const auto& e : ue::event_rates(R, X))
    rates.push(Json::object({{"type", e.type}, {"kind", e.kind}, {"of", e.of}, {"events", e.events}, {"exposure", round(e.exposure, 1)},
                             {"rate", round(e.rate, 6)}, {"error", round(e.error, 6)}}));
  return Json::object({{"buildings", memory_json(Ub, Mb)}, {"premises", premises}, {"spans", spans},
                       {"onComplex", Json::object({{"logZ", round(Z.log_z, 6)}, {"product", round(closed, 6)}, {"cells", static_cast<double>(chains.cells.size())},
                                                   {"pieces", static_cast<double>(chains.pieces.size())}, {"largest", static_cast<double>(Z.largest)},
                                                   {"buildings", static_cast<double>(seen.size())}})},
                       {"rates", rates}});
}

// The spin network of the moment at t (an MJD), of the record the request names, or of the one last built.
inline Json op_cut(const Json& r) {
  const auto B = record_of(r);
  const auto S = ue::cut(B->complex, num(r, "t", B->record.window.from));
  Json c = Json::array(), p = Json::array();
  for (auto x : S.couplings) c.push(static_cast<double>(x));
  for (auto x : S.primitives) p.push(static_cast<double>(x));
  return Json::object({{"t", round(S.t, 2)}, {"date", date_of(S.t)}, {"couplings", c}, {"primitives", p}});
}

}  // namespace es
