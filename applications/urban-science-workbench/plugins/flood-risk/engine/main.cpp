// The engine of the plug-in "Flood risk": what the Freetown Flood Risk Assessment's maps reach (the urban-science domain's exposure, its 153,435
// modelled buildings against each depth map) and what its loss tables say of the year (tyche's exceedance), with the assessment's own
// figures beside what is computed here from its event loss table.
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each ({"ok": true, "result": …} or
// {"ok": false, "error": …}): with --serve until its input ends, as the workbench's server keeps it; otherwise the one request given as an
// argument. --data names the directory of the data library's products, of which it reads freetown-flood.
//
//   summary                       the assessment's annual loss, the same from its events and from its curves; its curves, and those read from its events
//   exposed {peril, rp, threshold?, res?}
//                                 what a map reaches: buildings and value, by band of depth, by ward, on H3 cells of a resolution; each ward's share of
//                                 its own value
//   status {}                     the engine answers

#include <tyche/descriptive/exceedance.hpp>
#include <urban/exposure/exposure.hpp>

#include <physis/computational_geometry/h3.hpp>
#include <science/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace {

using science::Json;
namespace H = physis::computational_geometry::h3;
namespace td = tyche::descriptive;
namespace ue = urban::exposure;

std::string& data_dir() { static std::string d = "."; return d; }
constexpr double DEG = std::numbers::pi / 180;
double number_or(const Json& j, const char* key, double fallback) { return j.has(key) && j[key].is_number() ? j[key].as_number() : fallback; }

std::string slurp(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("the data library has no " + path + ": run the recipe");
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

const Json& product() {
  static const Json p = Json::parse(slurp(data_dir() + "/freetown-flood.json"));
  return p;
}

std::vector<double> numbers(const Json& a) {
  std::vector<double> out;
  out.reserve(a.size());
  for (std::size_t k = 0; k < a.size(); k++) out.push_back(a[k].as_number());
  return out;
}
Json array_of(const std::vector<double>& v) {
  Json out = Json::array();
  for (const double x : v) out.push(x);
  return out;
}

Json op_summary() {
  const Json& P = product();
  const auto events = numbers(P["events"]["losses"]);
  const double years = P["events"]["catalogue_years"].as_number();
  const double from_events = td::average_annual_loss(events, years);
  const auto rp = numbers(P["curves"]["oep"]["return_period"]), oep = numbers(P["curves"]["oep"]["loss"]), aep = numbers(P["curves"]["aep"]["loss"]);
  std::vector<double> rp_up, oep_up, aep_up;
  for (std::size_t k = 0; k < rp.size(); k++) if (rp[k] > 1) { rp_up.push_back(rp[k]); oep_up.push_back(oep[k]); aep_up.push_back(aep[k]); }
  const auto area_oep = td::curve_area(rp_up, oep_up, 0, true), area_aep = td::curve_area(rp_up, aep_up, 0, true);
  // The events in order, taken as the years' largest: each year's second event is not in the published curve either; where they differ, it shows.
  const auto read = td::occurrence_curve(events, years, rp_up);
  Json out = Json::object({
      {"reported", P["reported_aal"]["all"]},
      {"from_events", from_events},
      {"catalogue_years", years},
      {"events", static_cast<double>(events.size())},
      {"area_oep", Json::object({{"within", area_oep.within}, {"below", area_oep.below}, {"beyond", area_oep.beyond}, {"total", area_oep.total()}})},
      {"area_aep", Json::object({{"within", area_aep.within}, {"below", area_aep.below}, {"beyond", area_aep.beyond}, {"total", area_aep.total()}})},
      {"return_period", array_of(rp_up)},
      {"oep", array_of(oep_up)},
      {"aep", array_of(aep_up)},
      {"read_from_events", array_of(read)},
      {"by_occupancy", P["reported_aal"]["by_occupancy"]},
      {"credit", P["credit"]},
  });
  return out;
}

// The wards, and which each building falls in, found once.
struct Buildings {
  std::vector<ue::Asset> assets;
  std::vector<double> ward_value;
};

bool inside(double x, double y, const Json& ring) {
  bool hit = false;
  for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
    const double xi = ring[i][0].as_number(), yi = ring[i][1].as_number(), xj = ring[j][0].as_number(), yj = ring[j][1].as_number();
    if ((yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) hit = !hit;
  }
  return hit;
}

const Buildings& buildings() {
  static Buildings B;
  if (!B.assets.empty()) return B;
  const Json& P = product();
  const std::string raw = slurp(data_dir() + "/" + P["buildings"]["file"].as_string());
  const auto* v = reinterpret_cast<const double*>(raw.data());
  const std::size_t n = raw.size() / (4 * sizeof(double));
  struct Box { double w, s, e, n; std::size_t ward; const Json* outer; };
  std::vector<Box> boxes;
  const Json& wards = P["wards"];
  B.ward_value.assign(wards.size(), 0.0);
  for (std::size_t k = 0; k < wards.size(); k++) {
    const Json& polys = wards[k]["polygons"];
    for (std::size_t p = 0; p < polys.size(); p++) {
      const Json& r = polys[p]["outer"];
      Box b{1e9, 1e9, -1e9, -1e9, k, &r};
      for (std::size_t i = 0; i < r.size(); i++) { b.w = std::min(b.w, r[i][0].as_number()); b.e = std::max(b.e, r[i][0].as_number()); b.s = std::min(b.s, r[i][1].as_number()); b.n = std::max(b.n, r[i][1].as_number()); }
      boxes.push_back(b);
    }
  }
  B.assets.reserve(n);
  for (std::size_t k = 0; k < n; k++) {
    ue::Asset a{v[4 * k], v[4 * k + 1], v[4 * k + 2], -1};
    for (const Box& b : boxes) {
      if (a.lon < b.w || a.lon > b.e || a.lat < b.s || a.lat > b.n) continue;
      if (inside(a.lon, a.lat, *b.outer)) { a.zone = static_cast<int>(b.ward); B.ward_value[b.ward] += a.value; break; }
    }
    B.assets.push_back(a);
  }
  return B;
}

ue::Grid grid_of(const std::string& peril, int rp) {
  const Json& maps = product()["maps"];
  for (std::size_t k = 0; k < maps.size(); k++) {
    const Json& m = maps[k];
    if (m["peril"].text() != peril || static_cast<int>(m["return_period"].as_number()) != rp) continue;
    ue::Grid g;
    g.south = m["box"][0].as_number(), g.west = m["box"][1].as_number(), g.north = m["box"][2].as_number(), g.east = m["box"][3].as_number();
    g.width = static_cast<std::size_t>(m["width"].as_number());
    g.height = static_cast<std::size_t>(m["height"].as_number());
    g.unit = product()["depth_code"]["metres_per_unit"].as_number();
    const std::string raw = slurp(data_dir() + "/" + m["file"].as_string());
    g.values.assign(raw.begin(), raw.end());
    return g;
  }
  throw std::invalid_argument("no map of that peril and return period");
}

Json tally(const ue::Tally& t) { return Json::object({{"count", static_cast<double>(t.count)}, {"value", std::round(t.value)}}); }

Json op_exposed(const Json& r) {
  const std::string peril = r["peril"].text("FLRF");
  const int rp = static_cast<int>(number_or(r, "rp", 100));
  const double threshold = number_or(r, "threshold", 0.15);
  const int res = static_cast<int>(number_or(r, "res", 10));
  const ue::Grid g = grid_of(peril, rp);
  const Buildings& B = buildings();
  const Json& wards = product()["wards"];
  const std::vector<double> bands = {0.5, 1.0, 2.0, 4.0};
  const auto E = ue::exposed(g, B.assets, threshold, bands, {}, {}, wards.size());
  // On cells.
  std::map<std::uint64_t, std::pair<double, double>> cells;
  for (const auto& a : B.assets) {
    const double v = g.at(a.lat, a.lon);
    if (v < threshold || v <= 0) continue;
    auto& c = cells[H::cell_at({a.lat * DEG, a.lon * DEG}, res)];
    c.first += 1, c.second += a.value;
  }
  Json ids = Json::array(), rings = Json::array(), count = Json::array(), value = Json::array();
  for (const auto& [cell, c] : cells) {
    ids.push(H::to_string(cell));
    Json ring = Json::array();
    for (const auto& p : H::boundary(cell)) ring.push(Json::array({std::round(p.lng / DEG * 1e6) / 1e6, std::round(p.lat / DEG * 1e6) / 1e6}));
    rings.push(std::move(ring));
    count.push(c.first);
    value.push(std::round(c.second));
  }
  Json by_ward = Json::array();
  for (std::size_t k = 0; k < wards.size(); k++) by_ward.push(Json::object({{"count", static_cast<double>(E.zones[k].count)}, {"value", std::round(E.zones[k].value)}, {"of", std::round(B.ward_value[k])}}));
  Json by_band = Json::array();
  for (const auto& t : E.by_band) by_band.push(tally(t));
  return Json::object({{"peril", peril}, {"rp", rp}, {"threshold", threshold}, {"reached", tally(E.reached)}, {"all", tally(E.all)}, {"outside", static_cast<double>(E.outside)}, {"bands", array_of(bands)}, {"by_band", std::move(by_band)}, {"by_ward", std::move(by_ward)}, {"res", res}, {"cells", std::move(ids)}, {"rings", std::move(rings)}, {"count", std::move(count)}, {"value", std::move(value)}});
}

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "summary") return op_summary();
  if (op == "exposed") return op_exposed(r);
  if (op == "status") return Json::object({{"ok", true}});
  throw std::invalid_argument("no operation \"" + op + "\"");
}

std::string respond(const std::string& line) {
  Json reply = Json::object();
  try {
    const Json r = Json::parse(line);
    if (r.has("id")) reply.set("id", r["id"]);
    reply.set("ok", true);
    reply.set("result", answer(r));
  } catch (const std::exception& e) {
    reply = Json::object({{"ok", false}, {"error", e.what()}});
  }
  return reply.dump();
}

}  // namespace

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  bool serve = false;
  std::string request;
  for (int i = 1; i < argc; i++) {
    if (std::strcmp(argv[i], "--serve") == 0) serve = true;
    else if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) data_dir() = argv[++i];
    else request = argv[i];
  }
  if (serve) {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::cout << respond(line) << "\n" << std::flush;
    }
    return 0;
  }
  if (request.empty()) std::getline(std::cin, request);
  std::cout << respond(request) << "\n";
  return 0;
}
