// The engine of the plug-in "The territory": what is known of the world's places, binned on a global grid of hexagons
// (physis's H3) and combined (nomos's composite measures). Any quantity at a place can be moved onto the same cells, the idea
// of giga-spatial and of Space2Stats, so that a population and a building and a flood can be compared where they stand.
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each ({"ok": true, "result": …} or
// {"ok": false, "error": …}): with --serve until its input ends, as the workbench's server keeps it; otherwise the one
// request given as an argument. --data names the directory of the data library's products, of which it reads
// fields-ghsl-population, regions-countries and areas-in-dispute.
//
//   earth_people {res}                     the world's people on the cells of a resolution: each cell's ring, people, area and
//                                          the economy most of its people live in (the regions product, after the naming policy); a
//                                          cell holding any part of an area in dispute has no economy, so that no claimant's colour falls on it
//   bin_points   {res, points: [[lat, lon, weight], ...]}
//                                          any points' weights summed on the cells of a resolution, with each cell's ring and area
//   composite    {factors: [{values, weight, invert?}], low?, high?}
//                                          a composite of factors on the same cells: each scaled to 0..5 between its low and high
//                                          percentiles, a missing value left out of that cell's mean
//   status       {}                        the engine answers

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <nomos/measurement/composite.hpp>
#include <physis/computational_geometry/h3.hpp>
#include <science/json.hpp>

namespace {

using science::Json;
namespace H = physis::computational_geometry::h3;

std::string& data_dir() { static std::string d = "."; return d; }

double number_or(const Json& j, const char* key, double fallback) { return j.has(key) && j[key].is_number() ? j[key].as_number() : fallback; }
constexpr double DEG = std::numbers::pi / 180;

std::string slurp(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("the data library has no " + path + ": run the recipe");
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// A cell's ring, in degrees [lon, lat], unwrapped so that a cell across the antimeridian stays one polygon.
Json ring_of(std::uint64_t cell) {
  Json out = Json::array();
  double previous = 1e9;
  for (const auto& p : H::boundary(cell)) {
    double lon = p.lng / DEG;
    if (previous < 1e8) {
      while (lon - previous > 180) lon -= 360;
      while (lon - previous < -180) lon += 360;
    }
    previous = lon;
    out.push(Json::array({std::round(lon * 1e4) / 1e4, std::round(p.lat / DEG * 1e4) / 1e4}));
  }
  return out;
}

struct Binned { double people = 0; std::map<int, double> by_region; };

const std::map<std::uint64_t, Binned>& people_on(int res) {
  static std::map<int, std::map<std::uint64_t, Binned>> kept;
  auto it = kept.find(res);
  if (it != kept.end()) return it->second;
  const Json meta = Json::parse(slurp(data_dir() + "/fields-ghsl-population.json"));
  const Json& g = meta["grid"];
  const int nlon = static_cast<int>(g["nlon"].as_number()), nlat = static_cast<int>(g["nlat"].as_number());
  const double lon0 = g["lon0"].as_number(), dlon = g["dlon"].as_number(), lat0 = g["lat0"].as_number(), dlat = g["dlat"].as_number();
  const Json& v = meta["variables"][0];
  const double scale = v["scale"].as_number(), offset = v["offset"].as_number();
  const std::string raw = slurp(data_dir() + "/" + meta["file"].as_string());
  const auto* q = reinterpret_cast<const std::int16_t*>(raw.data());
  const Json regions = Json::parse(slurp(data_dir() + "/regions-countries.json"));
  const Json& rg = regions["grid"];
  const int rlon = static_cast<int>(rg["nlon"].as_number()), rlat = static_cast<int>(rg["nlat"].as_number());
  const double r0 = rg["lon0"].as_number(), rd = rg["dlon"].as_number(), s0 = rg["lat0"].as_number(), sd = rg["dlat"].as_number();
  const std::string rraw = slurp(data_dir() + "/" + regions["file"].as_string());
  const auto* rid = reinterpret_cast<const std::uint16_t*>(rraw.data());
  auto& out = kept[res];
  for (int j = 0; j < nlat; j++) {
    const double lat = lat0 - j * dlat;
    for (int i = 0; i < nlon; i++) {
      const std::int16_t value = q[j * nlon + i];
      if (value == -32768) continue;
      const double people = value * scale + offset;
      if (people <= 0.5) continue;
      const double lon = lon0 + i * dlon;
      Binned& b = out[H::cell_at({lat * DEG, lon * DEG}, res)];
      b.people += people;
      const int ri = std::clamp(static_cast<int>(std::lround((lon - r0) / rd)), 0, rlon - 1), rj = std::clamp(static_cast<int>(std::lround((s0 - lat) / sd)), 0, rlat - 1);
      if (const int id = rid[rj * rlon + ri]) b.by_region[id] += people;
    }
  }
  return out;
}

bool in_rings(double lon, double lat, const Json& rings) {
  bool inside = false;
  for (std::size_t k = 0; k < rings.size(); k++) {
    const Json& ring = rings[k];
    for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
      const double xi = ring[i][0].as_number(), yi = ring[i][1].as_number(), xj = ring[j][0].as_number(), yj = ring[j][1].as_number();
      if ((yi > lat) != (yj > lat) && lon < (xj - xi) * (lat - yi) / (yj - yi) + xi) inside = !inside;
    }
  }
  return inside;
}

// The cells of a resolution that hold any part of an area in dispute: those of its points, and of a grid of points within it finer than any cell asked for.
std::set<std::uint64_t> disputed_cells(int res) {
  std::set<std::uint64_t> out;
  const Json areas = Json::parse(slurp(data_dir() + "/areas-in-dispute.json"));
  const Json& polygons = areas["polygons"];
  constexpr double step = 0.05;
  for (std::size_t p = 0; p < polygons.size(); p++) {
    const Json& rings = polygons[p];
    double w = 1e9, s = 1e9, e = -1e9, n = -1e9;
    for (std::size_t k = 0; k < rings.size(); k++) for (std::size_t i = 0; i < rings[k].size(); i++) {
      const double lon = rings[k][i][0].as_number(), lat = rings[k][i][1].as_number();
      w = std::min(w, lon), e = std::max(e, lon), s = std::min(s, lat), n = std::max(n, lat);
      out.insert(H::cell_at({lat * DEG, lon * DEG}, res));
    }
    for (double lat = s; lat <= n; lat += step) for (double lon = w; lon <= e; lon += step)
      if (in_rings(lon, lat, rings)) out.insert(H::cell_at({lat * DEG, lon * DEG}, res));
  }
  return out;
}

Json op_earth_people(const Json& r) {
  const int res = static_cast<int>(number_or(r, "res", 3));
  if (res < 0 || res > 4) throw std::invalid_argument("the world's people are binned at resolutions 0 to 4");
  const Json regions = Json::parse(slurp(data_dir() + "/regions-countries.json"));
  const auto& bins = people_on(res);
  const std::set<std::uint64_t> disputed = disputed_cells(res);
  Json cells = Json::array(), rings = Json::array(), people = Json::array(), area = Json::array(), code = Json::array(), name = Json::array();
  double total = 0;
  for (const auto& [cell, b] : bins) {
    cells.push(H::to_string(cell));
    rings.push(ring_of(cell));
    people.push(std::round(b.people));
    area.push(std::round(H::cell_area_m2(cell) / 1e6));
    total += b.people;
    int best = 0;
    double most = 0;
    for (const auto& [id, n] : b.by_region) if (n > most) { most = n; best = id; }
    if (best && !disputed.count(cell)) {
      const Json& reg = regions["regions"][best - 1];
      code.push(reg["code"].text());
      name.push(reg["name"].is_string() ? reg["name"] : Json());
    } else {
      code.push(Json());
      name.push(Json());
    }
  }
  return Json::object({{"res", res}, {"total", total}, {"cells", std::move(cells)}, {"rings", std::move(rings)}, {"people", std::move(people)}, {"area_km2", std::move(area)}, {"code", std::move(code)}, {"name", std::move(name)}});
}

Json op_bin_points(const Json& r) {
  const int res = static_cast<int>(number_or(r, "res", 9));
  if (res < 0 || res > 15) throw std::invalid_argument("a resolution is 0 to 15");
  std::map<std::uint64_t, std::pair<double, int>> bins;
  const Json& points = r["points"];
  for (std::size_t k = 0; k < points.size(); k++) {
    const Json& p = points[k];
    auto& b = bins[H::cell_at({p[0].as_number() * DEG, p[1].as_number() * DEG}, res)];
    b.first += p.size() > 2 ? p[2].as_number() : 1.0;
    b.second++;
  }
  Json cells = Json::array(), rings = Json::array(), weight = Json::array(), count = Json::array(), area = Json::array();
  for (const auto& [cell, b] : bins) {
    cells.push(H::to_string(cell));
    rings.push(ring_of(cell));
    weight.push(b.first);
    count.push(b.second);
    area.push(H::cell_area_m2(cell));
  }
  return Json::object({{"res", res}, {"cells", std::move(cells)}, {"rings", std::move(rings)}, {"weight", std::move(weight)}, {"count", std::move(count)}, {"area_m2", std::move(area)}});
}

Json op_composite(const Json& r) {
  const Json& fs = r["factors"];
  const double low = number_or(r, "low", 0.02), high = number_or(r, "high", 0.98);
  std::vector<nomos::measurement::Factor> factors;
  for (std::size_t k = 0; k < fs.size(); k++) {
    nomos::measurement::Factor f;
    const Json& vs = fs[k]["values"];
    f.values.reserve(vs.size());
    for (std::size_t i = 0; i < vs.size(); i++) f.values.push_back(vs[i].is_number() ? vs[i].as_number() : std::nan(""));
    f.weight = number_or(fs[k], "weight", 1);
    f.invert = fs[k].has("invert") && fs[k]["invert"].is_bool() && fs[k]["invert"].as_bool();
    factors.push_back(std::move(f));
  }
  const auto c = nomos::measurement::composite(factors, low, high, 5.0);
  Json out = Json::array();
  for (const double v : c) out.push(std::isfinite(v) ? Json(std::round(v * 1000) / 1000) : Json());
  return Json::object({{"index", std::move(out)}});
}

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "earth_people") return op_earth_people(r);
  if (op == "bin_points") return op_bin_points(r);
  if (op == "composite") return op_composite(r);
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
