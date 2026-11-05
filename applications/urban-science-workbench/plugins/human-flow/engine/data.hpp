// The records the engine of "Human flow" reads: the data library's products under the directory it is given
// (--data), and nothing else. Each is read once and kept.
//
//   street-sensors.json             Melbourne's pedestrian sensors: their names, places and directions
//   street-counts/<id>.json         a sensor's hourly counts (UTC), the offsets of civil time, and each way's counts
//   street-holidays-victoria.json   Victoria's public holidays (left out of the counters' weeks)
//   street-terrain-melbourne.json   CLUE's establishments within 200 m of each sensor, by census year (the attractors)
//   city-events-melbourne.json      the Hoddle Grid's record of events: the City's event permits among them
//
// The counters' weeks are those of the district-economy problem (research/experiments/2026-11-05-a-district-economy-
// the-problem, calibrate.cpp), so that the walkers stand on the calibration of its flow study: over the reference year,
// the 52 weeks from Monday 27 October 2025 (civil), with Victoria's public holidays left out, a counter's mean week is
// its tide's week (urban::rhythm::tide, the week's constituents alone at the middle of each hour) and its mean week each
// way its current meter (urban::flow::current_meter). A counter covering less than 80% of the year is not used.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <science/json.hpp>
#include <urban/flow/currents.hpp>
#include <urban/rhythm/record.hpp>
#include <urban/rhythm/tide.hpp>

namespace hf {

using science::Json;
namespace ur = urban::rhythm;

inline std::string& data_dir() {
  static std::string dir = "data/cache/products";
  return dir;
}

inline bool product_exists(const std::string& name) {
  std::ifstream f(data_dir() + "/" + name, std::ios::binary);
  return static_cast<bool>(f);
}

inline Json read_json(const std::string& name) {
  const std::string path = data_dir() + "/" + name;
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("no product " + name + " in the data library's cache: run the workbench's recipe (the data lamp)");
  std::stringstream ss;
  ss << f.rdbuf();
  return Json::parse(ss.str());
}

inline double num(const Json& j, const char* key, double fallback = std::nan("")) { return j[key].is_number() ? j[key].as_number() : fallback; }

// The reference year (the problem's README §5; calibrate.cpp): Monday 2025-10-27 00:00 AEDT for 52 weeks.
inline double reference_start() { static const double t = ur::parse_utc("2025-10-26T13:00:00Z"); return t; }
inline constexpr std::size_t reference_weeks = 52;
inline constexpr double least_coverage = 0.8;

// ---- sensors ---------------------------------------------------------------------------------------------------

struct Sensor {
  long long id = 0;
  std::string name, series, location_type;
  std::vector<std::string> directions;   // the compass names d1 and d2 count toward, or none
  double lat = 0, lon = 0;
};

inline const std::map<long long, Sensor>& sensors() {
  static std::map<long long, Sensor> list;
  static std::once_flag once;
  std::call_once(once, [] {
    const Json j = read_json("street-sensors.json");
    for (const auto& s : j["sensors"].as_array()) {
      Sensor x;
      x.id = static_cast<long long>(num(s, "id", 0));
      x.name = s["name"].text(), x.series = s["series"].text(), x.location_type = s["location_type"].text();
      x.lat = num(s, "lat", 0), x.lon = num(s, "lon", 0);
      if (s["directions"].is_array())
        for (const auto& d : s["directions"].as_array()) x.directions.push_back(d.text());
      list[x.id] = std::move(x);
    }
  });
  return list;
}

inline const std::set<long long>& holiday_days() {
  static std::set<long long> days;
  static std::once_flag once;
  std::call_once(once, [] {
    const Json j = read_json("street-holidays-victoria.json");
    for (const auto& h : j["holidays"].as_array()) days.insert(static_cast<long long>(std::floor(ur::parse_utc(h["date"].text()) / ur::day_s)));
  });
  return days;
}

// A counter's weeks over the reference year: the mean count at each hour of the civil week (the tide's week), and,
// where it counts each way, the mean each way (its current meter). `covered` is the share of the year's hours counted.
struct Weeks {
  bool usable = false, directed = false;
  double covered = 0;
  std::array<double, 168> gross{};
  std::array<std::array<double, 168>, 2> way{};
  std::array<std::string, 2> toward;
  // The counts of each civil hour of a date (24), from the whole record (NaN where missing): for a date's surge.
  ur::Record record;
};

inline ur::Record window(const ur::Record& rec, double from, std::size_t hours) {
  ur::Record w;
  w.id = rec.id, w.t0 = from, w.step = rec.step, w.offsets = rec.offsets;
  w.y.assign(hours, std::nan(""));
  const long long i0 = static_cast<long long>(std::llround((from - rec.t0) / rec.step));
  for (std::size_t i = 0; i < hours; i++) {
    const long long j = i0 + static_cast<long long>(i);
    if (j >= 0 && j < static_cast<long long>(rec.size())) w.y[i] = rec.y[static_cast<std::size_t>(j)];
  }
  return w;
}

inline const Weeks& weeks(long long id) {
  static std::map<long long, Weeks> cache;
  static std::mutex lock;
  {
    std::lock_guard<std::mutex> g(lock);
    if (auto it = cache.find(id); it != cache.end()) return it->second;
  }
  Weeks out;
  const auto& S = sensors();
  const auto it = S.find(id);
  const std::string file = it == S.end() ? "" : (it->second.series.empty() ? "street-counts/" + std::to_string(id) + ".json" : it->second.series);
  if (it != S.end() && product_exists(file)) {
    const Json j = read_json(file);
    ur::Record full;
    full.id = std::to_string(id), full.t0 = ur::parse_utc(j["t0_utc"].text()), full.step = num(j, "step_s", 3600);
    for (const auto& v : j["counts"].as_array()) full.y.push_back(v.is_number() ? v.as_number() : std::nan(""));
    for (const auto& o : j["offsets"].as_array()) full.offsets.push_back({ur::parse_utc(o[0].text()), o[1].as_number()});
    if (full.offsets.empty() || full.offsets.front().first > full.t0) full.offsets.insert(full.offsets.begin(), {full.t0 - 86400, full.offsets.empty() ? 10.0 : full.offsets.front().second});
    out.record = full;
    const std::size_t hours = reference_weeks * 168;
    ur::Record typical = window(full, reference_start(), hours);
    std::size_t have = 0;
    for (std::size_t i = 0; i < typical.size(); i++) {
      if (std::isfinite(typical.y[i])) have++;
      if (holiday_days().count(typical.civil_day(typical.time(i)))) typical.y[i] = std::nan("");
    }
    out.covered = static_cast<double>(have) / static_cast<double>(hours);
    if (out.covered >= least_coverage) {
      const auto tide = ur::tide(typical, reference_start(), reference_start() + static_cast<double>(hours) * 3600.0);
      for (std::size_t h = 0; h < 168; h++) out.gross[h] = std::max(0.0, ur::level(tide, static_cast<double>(h) + 0.5, ur::of_the_week));
      out.usable = true;
      const Json& d = j["directions"];
      if (it->second.directions.size() == 2 && d.is_object()) {
        ur::Record d1, d2;
        d1.id = d2.id = full.id;
        d1.t0 = d2.t0 = ur::parse_utc(d["t0_utc"].text());
        d1.step = d2.step = num(d, "step_s", 3600);
        d1.offsets = d2.offsets = full.offsets;
        for (const auto& v : d["d1"].as_array()) d1.y.push_back(v.is_number() ? v.as_number() : std::nan(""));
        for (const auto& v : d["d2"].as_array()) d2.y.push_back(v.is_number() ? v.as_number() : std::nan(""));
        for (std::size_t i = 0; i < d1.size(); i++)
          if (holiday_days().count(d1.civil_day(d1.time(i)))) d1.y[i] = d2.y[i] = std::nan("");
        try {
          const auto cm = urban::flow::current_meter(d1, d2, reference_start(), reference_start() + static_cast<double>(hours) * 3600.0);
          for (std::size_t h = 0; h < 168; h++) out.way[0][h] = cm.way[0][h], out.way[1][h] = cm.way[1][h];
          out.toward = {it->second.directions[0], it->second.directions[1]};
          out.directed = true;
        } catch (const std::exception&) {
          out.directed = false;
        }
      }
    }
  }
  std::lock_guard<std::mutex> g(lock);
  return cache.emplace(id, std::move(out)).first->second;
}

// ---- CLUE's establishments about the sensors ---------------------------------------------------------------------

struct Establishment {
  double x = 0, y = 0;   // m east and north of its sensor
  std::string use, name, anzsic;
  int year = 0;
};

struct Terrain {
  int latest = 0;
  std::map<long long, std::vector<Establishment>> by_sensor;
};

inline const Terrain& terrain() {
  static Terrain t;
  static std::once_flag once;
  std::call_once(once, [] {
    const Json j = read_json("street-terrain-melbourne.json");
    for (const auto& y : j["years"].as_array()) t.latest = std::max(t.latest, static_cast<int>(y.as_number()));
    for (const auto& [key, v] : j["sensors"].as_object()) {
      auto& list = t.by_sensor[std::stoll(key)];
      for (const auto& e : v["establishments"].as_array()) {
        const int year = static_cast<int>(num(e, "year", 0));
        if (year != t.latest) continue;
        Establishment x;
        x.year = year, x.x = num(e, "x_m", 0), x.y = num(e, "y_m", 0), x.use = e["use"].text(), x.name = e["name"].text();
        x.anzsic = e["anzsic4"].is_number() ? std::to_string(static_cast<long long>(e["anzsic4"].as_number())) : e["anzsic4"].text();
        list.push_back(std::move(x));
      }
    }
  });
  return t;
}

// ---- the City's event permits (in the record of events) --------------------------------------------------------------

struct Permit {
  std::string id, kind, from, to, location;
  std::vector<std::string> streets;
};

inline const std::vector<Permit>& permits() {
  static std::vector<Permit> list;
  static std::once_flag once;
  std::call_once(once, [] {
    const Json j = read_json("city-events-melbourne.json");
    const auto& cols = j["events"]["columns"].as_array();
    std::map<std::string, std::size_t> at;
    for (std::size_t k = 0; k < cols.size(); k++) at[cols[k].text()] = k;
    for (const auto& r : j["events"]["rows"].as_array()) {
      if (r[at["source"]].text() != "melbourne-event-permits" || r[at["type"]].text() != "regime") continue;
      Permit p;
      p.id = r[at["id"]].text(), p.kind = r[at["kind"]].text(), p.from = r[at["from"]].text(), p.to = r[at["to"]].text();
      p.location = r[at["detail"]]["location"].text();
      for (const auto& s : r[at["match"]]["streets"].as_array()) p.streets.push_back(s.text());
      list.push_back(std::move(p));
    }
  });
  return list;
}

}  // namespace hf
