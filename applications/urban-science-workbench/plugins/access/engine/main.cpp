// The engine of the plug-in "Reach and use": how far one can go from a place over the ways drawn there (the urban-science domain's
// access, on morphe's shortest paths), and how the city's walkers share themselves among its places from day to day (its usage, the Urban
// Space Usage Index, on physis's H3 and Melbourne's counters).
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each ({"ok": true, "result": …} or
// {"ok": false, "error": …}): with --serve until its input ends, as the workbench's server keeps it; otherwise the one request given as an
// argument. --data names the directory of the data library's products, of which it reads street-sensors and street-counts.
//
//   reach  {ways: [[speed, x0, y0, x1, y1, ...], ...], from: [[x, y], ...], walk?, stitch?, limits?: [s, ...]}
//                         the time to each point of each way (m and s in a local frame); the metres within each limit
//   usage  {date, weeks?, res?}
//                         on a date (YYYY-MM-DD, Melbourne's), the day's share at each H3 cell holding counters, and its departure from the same
//                         weekday over the `weeks` before it, for the panel reporting on all of them
//   days   {cell, date, weeks?, res?}
//                         a cell's index and departure over the year before a date
//   status {}             the engine answers

#include <urban/access/access.hpp>
#include <urban/usage/usage.hpp>

#include <physis/computational_geometry/h3.hpp>
#include <science/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace {

using science::Json;
namespace H = physis::computational_geometry::h3;
namespace ua = urban::access;
namespace uu = urban::usage;

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

Json op_reach(const Json& r) {
  std::vector<ua::Way> ways;
  const Json& given = r["ways"];
  ways.reserve(given.size());
  for (std::size_t k = 0; k < given.size(); k++) {
    const Json& w = given[k];
    ua::Way way;
    way.speed = w[0].as_number();
    way.xy.reserve(w.size() - 1);
    for (std::size_t i = 1; i < w.size(); i++) way.xy.push_back(w[i].as_number());
    ways.push_back(std::move(way));
  }
  ua::Network n = ua::network_of(ways, number_or(r, "snap", 1.5));
  const double walk = number_or(r, "walk", 1.34), stitch = number_or(r, "stitch", 8);
  const std::size_t joined = ua::stitch_ends(n, stitch, walk);
  std::vector<std::pair<double, double>> from;
  for (std::size_t k = 0; k < r["from"].size(); k++) from.push_back({r["from"][k][0].as_number(), r["from"][k][1].as_number()});
  const ua::Reach R = ua::reach(n, from, walk);
  std::vector<double> limits = {300, 600, 900, 1800, 2700, 3600};
  if (r.has("limits")) {
    limits.clear();
    for (std::size_t k = 0; k < r["limits"].size(); k++) limits.push_back(r["limits"][k].as_number());
  }
  const auto within = ua::length_within(n, R.time, limits);
  Json times = Json::array();
  for (const auto& points : n.points) {
    Json line = Json::array();
    for (const auto p : points) line.push(std::isfinite(R.time[p]) ? Json(std::round(R.time[p] * 10) / 10) : Json());
    times.push(std::move(line));
  }
  Json lim = Json::array(), len = Json::array();
  for (std::size_t k = 0; k < limits.size(); k++) { lim.push(limits[k]); len.push(std::round(within[k])); }
  return Json::object({{"times", std::move(times)}, {"limits", std::move(lim)}, {"metres", std::move(len)}, {"points", static_cast<double>(n.x.size())}, {"stitched", static_cast<double>(joined)}, {"gap", R.start_gap}});
}

// Civil days from Howard Hinnant's algorithm, to count days without a calendar library.
long days_of(int y, unsigned m, unsigned d) {
  y -= m <= 2;
  const long era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + static_cast<long>(doe) - 719468;
}
long days_of(const std::string& date) {
  if (date.size() < 10) throw std::invalid_argument("a date is YYYY-MM-DD");
  return days_of(std::stoi(date.substr(0, 4)), static_cast<unsigned>(std::stoi(date.substr(5, 2))), static_cast<unsigned>(std::stoi(date.substr(8, 2))));
}
std::string date_of(long z) {
  z += 719468;
  const long era = (z >= 0 ? z : z - 146096) / 146097;
  const unsigned doe = static_cast<unsigned>(z - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const long y = static_cast<long>(yoe) + era * 400;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  const unsigned d = doy - (153 * mp + 2) / 5 + 1;
  const unsigned m = mp + (mp < 10 ? 3 : -9);
  char buf[16];
  std::snprintf(buf, sizeof buf, "%04ld-%02u-%02u", y + (m <= 2), m, d);
  return buf;
}

// Melbourne's counters by civil day: each sensor's walkers on each whole day of its record (a day missing an hour is no record).
struct Record {
  struct Sensor { int id; std::string name; double lat, lon; std::uint64_t cell = 0; };
  std::vector<Sensor> sensors;
  long first = 0, last = 0;
  uu::Matrix daily;     // [sensor][day]
  int res = -1;
};

Record& record() {
  static Record R;
  if (!R.sensors.empty()) return R;
  const Json meta = Json::parse(slurp(data_dir() + "/street-sensors.json"));
  const Json& ss = meta["sensors"];
  // The zone's offset is taken as the standard time's, +10 h, so that a day of counts runs 00:00 to 24:00 standard time: a daylight-saving
  // hour moves a day's edge by an hour of night, which the shares scarcely see.
  constexpr long offset = 10 * 3600;
  std::vector<std::vector<std::pair<long, double>>> kept;
  long lo = std::numeric_limits<long>::max(), hi = std::numeric_limits<long>::min();
  for (std::size_t k = 0; k < ss.size(); k++) {
    const Json& s = ss[k];
    if (s["location_type"].text() != "Outdoor" || !s.has("series") || !s["series"].is_string()) continue;
    const Json series = Json::parse(slurp(data_dir() + "/" + s["series"].as_string()));
    const std::string t0 = series["t0_utc"].as_string();
    const long start = days_of(t0) * 86400 + std::stoi(t0.substr(11, 2)) * 3600 + offset;
    const Json& counts = series["counts"];
    std::map<long, std::pair<double, int>> days;
    for (std::size_t i = 0; i < counts.size(); i++) {
      const long at = start + static_cast<long>(i) * 3600;
      auto& d = days[at >= 0 ? at / 86400 : -((-at + 86399) / 86400)];
      if (counts[i].is_number()) { d.first += counts[i].as_number(); d.second++; }
    }
    std::vector<std::pair<long, double>> whole;
    for (const auto& [day, v] : days) if (v.second == 24) { whole.push_back({day, v.first}); lo = std::min(lo, day); hi = std::max(hi, day); }
    if (whole.empty()) continue;
    R.sensors.push_back({static_cast<int>(s["id"].as_number()), s["name"].text(), s["lat"].as_number(), s["lon"].as_number()});
    kept.push_back(std::move(whole));
  }
  if (R.sensors.empty()) throw std::runtime_error("no counters in the data library");
  R.first = lo, R.last = hi;
  const auto days = static_cast<std::size_t>(hi - lo + 1);
  R.daily.assign(R.sensors.size(), std::vector<double>(days, std::numeric_limits<double>::quiet_NaN()));
  for (std::size_t k = 0; k < kept.size(); k++) for (const auto& [d, v] : kept[k]) R.daily[k][static_cast<std::size_t>(d - lo)] = v;
  return R;
}

// The counts summed on the cells of a resolution: only sensors that report on every day asked for, so that a changing panel moves no share.
struct Cells {
  std::vector<std::uint64_t> ids;
  std::vector<std::vector<std::size_t>> members;
};

Cells cells_of(Record& R, int res, const std::vector<std::size_t>& panel) {
  std::map<std::uint64_t, std::vector<std::size_t>> by;
  for (const auto k : panel) by[H::cell_at({R.sensors[k].lat * DEG, R.sensors[k].lon * DEG}, res)].push_back(k);
  Cells c;
  for (auto& [id, m] : by) { c.ids.push_back(id); c.members.push_back(std::move(m)); }
  return c;
}

Json ring_of(std::uint64_t cell) {
  Json out = Json::array();
  for (const auto& p : H::boundary(cell)) out.push(Json::array({std::round(p.lng / DEG * 1e6) / 1e6, std::round(p.lat / DEG * 1e6) / 1e6}));
  return out;
}

Json op_usage(const Json& r) {
  Record& R = record();
  const int res = static_cast<int>(number_or(r, "res", 9));
  const int weeks = static_cast<int>(number_or(r, "weeks", 12));
  const long asked = std::clamp(days_of(r["date"].as_string()), R.first + 7L * weeks, R.last);
  const auto at = [&](long d) { return static_cast<std::size_t>(d - R.first); };
  // The panel: reporting on that day and on the same weekday in each of the weeks before. The record has gaps (the archive ends in 2022, the
  // live feed keeps two years): a date inside one takes the nearest before it that a good part of the counters report on.
  std::vector<std::size_t> panel;
  long day = asked;
  for (; day >= R.first + 7L * weeks; day--) {
    panel.clear();
    for (std::size_t k = 0; k < R.sensors.size(); k++) {
      bool all = std::isfinite(R.daily[k][at(day)]);
      for (int w = 1; w <= weeks && all; w++) all = std::isfinite(R.daily[k][at(day - 7L * w)]);
      if (all) panel.push_back(k);
    }
    if (panel.size() >= 12) break;
  }
  if (panel.size() < 12) throw std::runtime_error("no good part of the counters report before " + date_of(asked));
  const bool moved = day != days_of(r["date"].as_string());
  const Cells C = cells_of(R, res, panel);
  const auto T = static_cast<std::size_t>(weeks) + 1;
  uu::Matrix counts(C.ids.size(), std::vector<double>(T, 0.0));
  for (std::size_t c = 0; c < C.ids.size(); c++) for (const auto k : C.members[c]) for (int w = 0; w <= weeks; w++) counts[c][static_cast<std::size_t>(w)] += R.daily[k][at(day - 7L * w)];
  const auto I = uu::shares(counts);
  std::vector<bool> base(T, true);
  base[0] = false;
  const auto B = uu::baseline_of(I, base);
  const auto Z = uu::departures(I, B);
  Json cells = Json::array(), rings = Json::array(), share = Json::array(), z = Json::array(), walkers = Json::array(), names = Json::array(), mean = Json::array();
  for (std::size_t c = 0; c < C.ids.size(); c++) {
    cells.push(H::to_string(C.ids[c]));
    rings.push(ring_of(C.ids[c]));
    share.push(I[c][0]);
    mean.push(B.mean[c]);
    z.push(std::isfinite(Z[c][0]) ? Json(std::round(Z[c][0] * 100) / 100) : Json());
    walkers.push(counts[c][0]);
    Json n = Json::array();
    for (const auto k : C.members[c]) n.push(R.sensors[k].name);
    names.push(std::move(n));
  }
  return Json::object({{"date", date_of(day)}, {"moved", moved}, {"res", res}, {"weeks", weeks}, {"panel", static_cast<double>(panel.size())}, {"cells", std::move(cells)}, {"rings", std::move(rings)}, {"share", std::move(share)}, {"mean", std::move(mean)}, {"z", std::move(z)}, {"walkers", std::move(walkers)}, {"names", std::move(names)}, {"first", date_of(R.first + 7L * weeks)}, {"last", date_of(R.last)}});
}

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "reach") return op_reach(r);
  if (op == "usage") return op_usage(r);
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
