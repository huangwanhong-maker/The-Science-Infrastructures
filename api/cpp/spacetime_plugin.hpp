// Plug-ins over data in space and time in C++ (ADR 0079): the class of the Python plug-ins of api/python/science/spacetime.py, with the same
// records, parameters and results, for a program that a server compiles and runs. The header belongs to no application: spacetime.py's Host
// compiles a plug-in with -I api/cpp and the include directories its application gives (the Urban Science Workbench gives its domains and the
// foundations, so that a plug-in may use them). A plug-in derives from Analysis, describes itself, runs, and ends with SCIENCE_SPACETIME_PLUGIN(Type):
//
//     #include "spacetime_plugin.hpp"
//     namespace st = science::spacetime;
//     struct Late : st::Analysis {
//       science::Json describe() const override {
//         return science::Json::object({{"id", "late"}, {"title", "The last ones"}, {"description", "..."},
//                                       {"params", science::Json::array({st::param("n", "How many", 100, "integer")})}});
//       }
//       std::vector<science::Json> run(const std::vector<st::Record>& records, const st::Params& p, const st::Context&, st::Progress&) override {
//         st::Events out = records.at(0).events;
//         // ... keep what is wanted ...
//         return {out.to_json("The last ones")};
//       }
//     };
//     SCIENCE_SPACETIME_PLUGIN(Late)
//
// The program is run twice. With --describe it prints its description as JSON. With --run <dir> it reads dir/input.json, {"records", "params", "context"}, and
// writes dir/output.json, {"results": [...]}, reporting progress on standard output as lines "progress\t<fraction>\t<text>". Records and results have the forms
// spacetime.py sets out: degrees on WGS 84, time in ms UTC, a missing value NaN here and null in the files. A run that throws prints its reason on standard
// error and exits with 1; the host reports it.

#pragma once

#include "science/json.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace science::spacetime {

inline constexpr double missing = std::numeric_limits<double>::quiet_NaN();
inline constexpr double day_ms = 86400000.0;

inline std::vector<double> column(const Json& a) {
  std::vector<double> out;
  out.reserve(a.size());
  for (std::size_t k = 0; k < a.size(); k++) out.push_back(a[k].is_number() ? a[k].as_number() : missing);
  return out;
}
inline Json json_of(const std::vector<double>& v) {
  Json out = Json::array();
  for (const double x : v) out.push(std::isfinite(x) ? Json(x) : Json());
  return out;
}
using Values = std::map<std::string, std::vector<double>>;
inline Json json_of(const Values& v) {
  Json out = Json::object();
  for (const auto& [name, col] : v) out.set(name, json_of(col));
  return out;
}

struct Events {
  std::string title;
  std::vector<double> time, lat, lon;
  Values values;
  std::vector<std::string> labels;
  std::size_t size() const { return time.size(); }
  void push(double t, double la, double lo) { time.push_back(t); lat.push_back(la); lon.push_back(lo); }
  Json to_json(const std::string& name, const Json& show = {}, const std::string& note = "") const {
    Json out = Json::object({{"kind", "events"}, {"title", name}, {"time", json_of(time)}, {"lat", json_of(lat)}, {"lon", json_of(lon)}, {"values", json_of(values)}});
    if (!labels.empty()) out.set("labels", Json(labels));
    if (!show.is_null()) out.set("show", show);
    if (!note.empty()) out.set("note", note);
    return out;
  }
};

struct Track {
  std::string id, name;
  std::vector<double> time, lat, lon;
  Values values;
};

// A place with a series: its id, its ring of [lon, lat] (degrees), a name.
struct Cell {
  std::string id, name;
  std::vector<std::pair<double, double>> ring;
};

// The cells' values at each time: values[name][time][cell].
inline Json cells_result(const std::string& title, const std::vector<double>& time, const std::vector<Cell>& cells,
                         const std::map<std::string, std::vector<std::vector<double>>>& values, const Json& show = {}, const std::string& note = "") {
  Json items = Json::array();
  for (const Cell& c : cells) {
    Json ring = Json::array();
    for (const auto& [lon, lat] : c.ring) ring.push(Json::array({lon, lat}));
    Json item = Json::object({{"id", c.id}, {"ring", std::move(ring)}});
    if (!c.name.empty()) item.set("name", c.name);
    items.push(std::move(item));
  }
  Json vals = Json::object();
  for (const auto& [name, frames] : values) {
    Json f = Json::array();
    for (const auto& row : frames) f.push(json_of(row));
    vals.set(name, std::move(f));
  }
  Json out = Json::object({{"kind", "cells"}, {"title", title}, {"time", json_of(time)}, {"cells", std::move(items)}, {"values", std::move(vals)}});
  if (!show.is_null()) out.set("show", show);
  if (!note.empty()) out.set("note", note);
  return out;
}

inline Json tracks_result(const std::string& title, const std::vector<Track>& tracks, const Json& show = {}, const std::string& note = "") {
  Json items = Json::array();
  for (const Track& t : tracks)
    items.push(Json::object({{"id", t.id}, {"name", t.name}, {"time", json_of(t.time)}, {"lat", json_of(t.lat)}, {"lon", json_of(t.lon)}, {"values", json_of(t.values)}}));
  Json out = Json::object({{"kind", "tracks"}, {"title", title}, {"tracks", std::move(items)}});
  if (!show.is_null()) out.set("show", show);
  if (!note.empty()) out.set("note", note);
  return out;
}

inline Json table(const std::string& title, const std::vector<std::string>& columns, const std::vector<std::vector<Json>>& rows, const std::string& note = "") {
  Json r = Json::array();
  for (const auto& row : rows) r.push(Json(Json::Array(row.begin(), row.end())));
  return Json::object({{"kind", "table"}, {"title", title}, {"columns", Json(columns)}, {"rows", r}, {"note", note}});
}
inline Json note(const std::string& text) { return Json::object({{"kind", "note"}, {"text", text}}); }

// A record as a program reads it. `raw` holds the whole; `events` and `tracks` are read out for those kinds.
struct Record {
  std::string kind, title;
  Events events;
  std::vector<Track> tracks;
  Json raw;
  explicit Record(const Json& r) : kind(r["kind"].text()), title(r["title"].text()), raw(r) {
    if (kind == "events") {
      events.title = title;
      events.time = column(r["time"]), events.lat = column(r["lat"]), events.lon = column(r["lon"]);
      if (r["values"].is_object()) for (const auto& [name, col] : r["values"].as_object()) events.values[name] = column(col);
    } else if (kind == "tracks") {
      const Json& ts = r["tracks"];
      for (std::size_t k = 0; k < ts.size(); k++) {
        Track t;
        t.id = ts[k]["id"].is_string() ? ts[k]["id"].as_string() : ts[k]["id"].dump();
        t.name = ts[k]["name"].text();
        t.time = column(ts[k]["time"]), t.lat = column(ts[k]["lat"]), t.lon = column(ts[k]["lon"]);
        if (ts[k]["values"].is_object()) for (const auto& [name, col] : ts[k]["values"].as_object()) t.values[name] = column(col);
        tracks.push_back(std::move(t));
      }
    }
  }
};

struct Params {
  Json given;
  double number(const std::string& name, double fallback = 0) const { const Json& v = given[name]; return v.is_number() ? v.as_number() : fallback; }
  int integer(const std::string& name, int fallback = 0) const { return static_cast<int>(number(name, fallback)); }
  std::string text(const std::string& name, const std::string& fallback = "") const { return given[name].text(fallback); }
  bool flag(const std::string& name, bool fallback = false) const { const Json& v = given[name]; return v.is_bool() ? v.as_bool() : fallback; }
};

struct Context {
  double lat = missing, lon = missing, range_m = missing, now_ms = missing, window_ms = missing;
  explicit Context(const Json& c) {
    const auto read = [&](const char* k, double& x) { if (c[k].is_number()) x = c[k].as_number(); };
    read("lat", lat), read("lon", lon), read("range_m", range_m), read("now_ms", now_ms), read("window_ms", window_ms);
  }
};

struct Progress {
  void operator()(double fraction, const std::string& text = "") const {
    std::printf("progress\t%.4f\t%s\n", fraction, text.c_str());
    std::fflush(stdout);
  }
};

// A parameter: kind "number", "integer", "choice" (with options) or "boolean".
inline Json param(const std::string& name, const std::string& label, Json default_value, const std::string& kind = "number", Json minimum = {}, Json maximum = {}, Json options = {}, const std::string& help = "") {
  Json p = Json::object({{"name", name}, {"label", label}, {"kind", kind}, {"default", std::move(default_value)}});
  if (!minimum.is_null()) p.set("min", std::move(minimum));
  if (!maximum.is_null()) p.set("max", std::move(maximum));
  if (!options.is_null()) p.set("options", std::move(options));
  if (!help.empty()) p.set("help", help);
  return p;
}

class Analysis {
 public:
  virtual ~Analysis() = default;
  virtual Json describe() const = 0;
  virtual std::vector<Json> run(const std::vector<Record>& records, const Params& params, const Context& context, Progress& progress) = 0;
};

inline int main_of(Analysis& plugin, int argc, char** argv) {
  const std::string what = argc > 1 ? argv[1] : "";
  try {
    if (what == "--describe") {
      Json d = plugin.describe();
      d.set("kind", "spacetime-analysis");
      d.set("language", "cpp");
      std::printf("%s\n", d.dump().c_str());
      return 0;
    }
    if (what == "--run" && argc > 2) {
      const std::string dir = argv[2];
      std::ifstream in(dir + "/input.json", std::ios::binary);
      if (!in) throw std::runtime_error(dir + "/input.json cannot be read");
      std::ostringstream text;
      text << in.rdbuf();
      const Json request = Json::parse(text.str());
      std::vector<Record> records;
      for (std::size_t k = 0; k < request["records"].size(); k++) records.emplace_back(request["records"][k]);
      Progress progress;
      const auto results = plugin.run(records, Params{request["params"]}, Context(request["context"]), progress);
      std::ofstream out(dir + "/output.json", std::ios::binary);
      out << Json::object({{"results", Json(Json::Array(results.begin(), results.end()))}}).dump();
      return out.good() ? 0 : 1;
    }
  } catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }
  std::fprintf(stderr, "usage: %s --describe | --run <dir>\n", argc > 0 ? argv[0] : "plugin");
  return 2;
}

}  // namespace science::spacetime

#define SCIENCE_SPACETIME_PLUGIN(Type)                           \
  int main(int argc, char** argv) {                              \
    Type plugin;                                                 \
    return science::spacetime::main_of(plugin, argc, argv);      \
  }
