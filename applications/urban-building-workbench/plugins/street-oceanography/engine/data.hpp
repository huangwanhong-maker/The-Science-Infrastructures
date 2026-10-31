// The data the engine reads: the data library's products under the directory it is given (--data), and
// nothing else. Each is read once and kept.
//
//   street-sensors.json             Melbourne's pedestrian sensors and their records
//   street-counts/<id>.json         a sensor's hourly counts in UTC, with Melbourne's offsets, and each direction's
//   street-weather-melbourne.json   hourly temperature, precipitation and wind (Open-Meteo)
//   street-holidays-victoria.json   Victoria's public holidays
//   street-terrain-melbourne.json   CLUE's establishments about each sensor, by census year, with their uses
//   street-holdings.json            OpenStreetMap's holdings along named streets, with their declared hours
//
// A product that is missing is refused with the recipe step that makes it.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <science/json.hpp>
#include <urban/rhythm/record.hpp>
#include <urban/street/street.hpp>

namespace so {

using science::Json;

inline std::string& data_dir() {
  static std::string dir = "data/cache/products";
  return dir;
}

inline Json read_json(const std::string& name) {
  const std::string path = data_dir() + "/" + name;
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("no product " + name + " in the data library's cache: run the workbench's recipe (Build the engines and fetch the data)");
  std::stringstream ss;
  ss << f.rdbuf();
  return Json::parse(ss.str());
}

inline double num(const Json& j, const char* key, double fallback = std::nan("")) { return j[key].is_number() ? j[key].as_number() : fallback; }

// ---- sensors and their records ------------------------------------------------------------------------

struct Sensor {
  long long id = 0;
  std::string code, name, location_type, series, first, last;
  std::vector<std::string> directions;   // toward which d1 and d2 count (the locations dataset's names), or none
  double lat = 0, lon = 0;
  long long hours = 0, missing = 0;
};

inline const std::vector<Sensor>& sensors() {
  static std::vector<Sensor> list;
  static bool loaded = false;
  if (!loaded) {
    const Json j = read_json("street-sensors.json");
    for (const auto& s : j["sensors"].as_array()) {
      Sensor x;
      x.id = static_cast<long long>(num(s, "id", 0));
      x.code = s["code"].text(), x.name = s["name"].text(), x.location_type = s["location_type"].text();
      x.series = s["series"].text(), x.first = s["first_utc"].text(), x.last = s["last_utc"].text();
      x.lat = num(s, "lat", 0), x.lon = num(s, "lon", 0);
      x.hours = static_cast<long long>(num(s, "hours", 0)), x.missing = static_cast<long long>(num(s, "missing", 0));
      if (s["directions"].is_array())
        for (const auto& d : s["directions"].as_array()) x.directions.push_back(d.text());
      list.push_back(x);
    }
    loaded = true;
  }
  return list;
}

inline const Sensor& sensor(long long id) {
  for (const auto& s : sensors())
    if (s.id == id) return s;
  throw std::invalid_argument("no sensor " + std::to_string(id));
}

inline const urban::rhythm::Record& record(long long id) {
  static std::map<long long, urban::rhythm::Record> cache;
  static std::mutex lock;
  std::lock_guard<std::mutex> g(lock);
  auto it = cache.find(id);
  if (it != cache.end()) return it->second;
  const auto& s = sensor(id);
  const Json j = read_json(s.series.empty() ? "street-counts/" + std::to_string(id) + ".json" : s.series);
  urban::rhythm::Record r;
  r.id = std::to_string(id);
  r.t0 = urban::rhythm::parse_utc(j["t0_utc"].text());
  r.step = num(j, "step_s", 3600);
  for (const auto& v : j["counts"].as_array()) r.y.push_back(v.is_number() ? v.as_number() : std::nan(""));
  for (const auto& o : j["offsets"].as_array()) r.offsets.emplace_back(urban::rhythm::parse_utc(o[0].text()), o[1].as_number());
  if (r.offsets.empty()) r.offsets.emplace_back(r.t0 - 86400, 10.0);
  if (r.offsets.front().first > r.t0) r.offsets.insert(r.offsets.begin(), {r.t0 - 86400, r.offsets.front().second});
  if (cache.size() > 24) cache.clear();
  return cache.emplace(id, std::move(r)).first->second;
}

// A sensor's counts each way (the live record's directions): two records on the same hours, d1 toward
// directions[0] and d2 toward directions[1]. Sensors without them are refused.
struct Directions {
  std::string toward[2];
  urban::rhythm::Record way[2];
};

inline const Directions& directions(long long id) {
  static std::map<long long, Directions> cache;
  static std::mutex lock;
  std::lock_guard<std::mutex> g(lock);
  auto it = cache.find(id);
  if (it != cache.end()) return it->second;
  const auto& s = sensor(id);
  if (s.directions.size() != 2) throw std::invalid_argument("sensor " + std::to_string(id) + " has no named directions");
  const Json j = read_json(s.series.empty() ? "street-counts/" + std::to_string(id) + ".json" : s.series);
  const Json& d = j["directions"];
  if (!d.is_object()) throw std::invalid_argument("sensor " + std::to_string(id) + " has no counts each way");
  Directions out;
  for (int w = 0; w < 2; w++) {
    auto& r = out.way[w];
    out.toward[w] = s.directions[static_cast<std::size_t>(w)];
    r.id = std::to_string(id) + (w ? ":2" : ":1");
    r.t0 = urban::rhythm::parse_utc(d["t0_utc"].text());
    r.step = num(d, "step_s", 3600);
    for (const auto& v : d[w ? "d2" : "d1"].as_array()) r.y.push_back(v.is_number() ? v.as_number() : std::nan(""));
    for (const auto& o : j["offsets"].as_array()) r.offsets.emplace_back(urban::rhythm::parse_utc(o[0].text()), o[1].as_number());
    if (r.offsets.empty() || r.offsets.front().first > r.t0) r.offsets.insert(r.offsets.begin(), {r.t0 - 86400, r.offsets.empty() ? 10.0 : r.offsets.front().second});
  }
  return cache.emplace(id, std::move(out)).first->second;
}

// ---- weather and holidays, aligned to a record -------------------------------------------------------------

struct Weather { double t0 = 0, step = 3600; std::vector<double> temperature, precipitation, wind; };

inline const Weather& weather() {
  static Weather w;
  static bool loaded = false;
  if (!loaded) {
    const Json j = read_json("street-weather-melbourne.json");
    w.t0 = urban::rhythm::parse_utc(j["t0_utc"].text());
    w.step = num(j, "step_s", 3600);
    auto read = [&](const char* key, std::vector<double>& out) {
      for (const auto& v : j[key].as_array()) out.push_back(v.is_number() ? v.as_number() : std::nan(""));
    };
    read("temperature_c", w.temperature);
    read("precipitation_mm", w.precipitation);
    read("wind_m_s", w.wind);
    loaded = true;
  }
  return w;
}

// A weather series on the record's hours (NaN outside the weather's span).
inline std::vector<double> aligned(const std::vector<double>& series, const Weather& w, const urban::rhythm::Record& r) {
  std::vector<double> out(r.size(), std::nan(""));
  for (std::size_t i = 0; i < r.size(); i++) {
    const double k = (r.time(i) - w.t0) / w.step;
    const long long j = std::llround(k);
    if (j >= 0 && j < static_cast<long long>(series.size()) && std::fabs(k - static_cast<double>(j)) < 0.01) out[i] = series[static_cast<std::size_t>(j)];
  }
  return out;
}

struct Holiday { long long day = 0; std::string date, name; };

inline const std::vector<Holiday>& holidays() {
  static std::vector<Holiday> list;
  static bool loaded = false;
  if (!loaded) {
    const Json j = read_json("street-holidays-victoria.json");
    for (const auto& h : j["holidays"].as_array()) {
      Holiday x;
      x.date = h["date"].text(), x.name = h["name"].text();
      x.day = static_cast<long long>(std::floor(urban::rhythm::parse_utc(x.date) / 86400.0));
      list.push_back(x);
    }
    loaded = true;
  }
  return list;
}

inline std::vector<long long> holiday_days() {
  std::vector<long long> d;
  for (const auto& h : holidays()) d.push_back(h.day);
  std::sort(d.begin(), d.end());
  d.erase(std::unique(d.begin(), d.end()), d.end());
  return d;
}

// ---- terrain about the sensors (CLUE) -------------------------------------------------------------------------

struct Establishment { int year = 0; double x = 0, y = 0; std::string anzsic, use, name; };

struct Terrain {
  std::vector<std::string> uses;                                   // the use classes, by index
  std::vector<int> years;
  std::map<long long, std::vector<Establishment>> by_sensor;
  std::size_t use_index(const std::string& u) const {
    for (std::size_t i = 0; i < uses.size(); i++)
      if (uses[i] == u) return i;
    return uses.empty() ? 0 : uses.size() - 1;
  }
};

inline const Terrain& terrain() {
  static Terrain t;
  static bool loaded = false;
  if (!loaded) {
    const Json j = read_json("street-terrain-melbourne.json");
    for (const auto& u : j["uses"].as_array()) t.uses.push_back(u["id"].text());
    for (const auto& y : j["years"].as_array()) t.years.push_back(static_cast<int>(y.as_number()));
    for (const auto& [key, v] : j["sensors"].as_object()) {
      auto& list = t.by_sensor[std::stoll(key)];
      for (const auto& e : v["establishments"].as_array()) {
        Establishment x;
        x.year = static_cast<int>(num(e, "year", 0));
        x.x = num(e, "x_m", 0), x.y = num(e, "y_m", 0);
        x.anzsic = e["anzsic4"].is_number() ? std::to_string(static_cast<long long>(e["anzsic4"].as_number())) : e["anzsic4"].text();
        x.use = e["use"].text(), x.name = e["name"].text();
        list.push_back(std::move(x));
      }
    }
    loaded = true;
  }
  return t;
}

inline int latest_year(const Terrain& t) { return t.years.empty() ? 0 : *std::max_element(t.years.begin(), t.years.end()); }

// ---- holdings along named streets (OpenStreetMap) --------------------------------------------------------------

struct NamedStreet {
  urban::street::Street street;
  std::string city, timezone, note;
  std::vector<std::array<double, 2>> axis;      // lat, lon
  std::vector<std::string> raw_hours, status, kinds;
  std::vector<std::array<double, 2>> places;    // each holding's lat, lon (NaN when the product has none)
  std::string area;                              // the city model's area the street lies in
};

inline const std::vector<NamedStreet>& streets() {
  static std::vector<NamedStreet> list;
  static bool loaded = false;
  if (!loaded) {
    const Json j = read_json("street-holdings.json");
    // One set of use classes across streets: those named by the holdings, in order of appearance.
    std::vector<std::string> uses;
    auto use_of = [&](const std::string& u) {
      for (std::size_t i = 0; i < uses.size(); i++)
        if (uses[i] == u) return i;
      uses.push_back(u);
      return uses.size() - 1;
    };
    for (const auto& s : j["streets"].as_array()) {
      NamedStreet n;
      n.street.id = s["id"].text(), n.street.name = s["name"].text(), n.street.length = num(s, "length_m", 0);
      n.city = s["city"].text(), n.timezone = s["timezone"].text(), n.note = s["note"].text(), n.area = s["area"].text();
      for (const auto& p : s["axis"].as_array()) n.axis.push_back({p[0].as_number(), p[1].as_number()});
      for (const auto& h : s["holdings"].as_array()) {
        urban::street::Holding x;
        x.id = h["osm"].text(), x.name = h["name"].text();
        x.s = num(h, "s_m", 0), x.offset = num(h, "offset_m", 0);
        x.use = use_of(h["use"].text("other"));
        std::vector<std::array<double, 3>> dm;
        for (const auto& iv : h["hours"].as_array()) dm.push_back({iv[0].as_number(), iv[1].as_number(), iv[2].as_number()});
        x.hours = urban::street::week_intervals(dm);
        const std::string st = h["hours_status"].text("none");
        x.declared = (st == "ok" || st == "partial") && !x.hours.empty();
        std::string kind;
        if (h["tags"].is_object())
          for (const auto& [k, v] : h["tags"].as_object()) kind += (kind.empty() ? "" : ", ") + k + "=" + v.text();
        n.raw_hours.push_back(h["hours_raw"].text());
        n.status.push_back(st);
        n.kinds.push_back(kind);
        n.places.push_back({num(h, "lat"), num(h, "lon")});
        n.street.holdings.push_back(std::move(x));
      }
      list.push_back(std::move(n));
    }
    for (auto& n : list) {
      n.street.uses = uses;
      double longest = n.street.length;
      for (const auto& h : n.street.holdings) longest = std::max(longest, h.s);
      n.street.length = std::max(longest, 1.0);
    }
    loaded = true;
  }
  return list;
}

inline const NamedStreet& named_street(const std::string& id) {
  for (const auto& s : streets())
    if (s.street.id == id) return s;
  throw std::invalid_argument("no street " + id);
}

}  // namespace so
