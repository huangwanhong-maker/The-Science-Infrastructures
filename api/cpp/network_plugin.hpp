// Plug-ins over networks of places in C++ (ADR 0046): the class of the Python
// plug-ins of api/python/science/networks.py, with the same network, parameters
// and results, for a program that a server compiles and runs. Like the Python
// class, the header belongs to no application: networks.py's Host compiles a
// plug-in with -I api/cpp and the include directories its application gives
// (the Economy Atlas gives the economics domain and the foundations, so that a
// plug-in may use the domain itself). A plug-in derives from Analysis, describes
// itself, runs, and ends with SCIENCE_NETWORK_PLUGIN(Type):
//
//     #include "network_plugin.hpp"
//     namespace net = science::network;
//     struct Mine : net::Analysis {
//       science::Json describe() const override {
//         return science::Json::object({{"id", "mine"}, {"title", "Mine"}, {"description", "..."},
//                                       {"params", science::Json::array({net::param("k", "k", 3, "integer")})}});
//       }
//       std::vector<science::Json> run(const net::Network& n, const net::Params& p, net::Progress& progress) override {
//         std::map<std::string, double> b;
//         const auto balance = n.balance();
//         for (std::size_t i = 0; i < n.size(); i++) b[n.places[i].id] = balance[i];
//         return {net::values("balance", b)};
//       }
//     };
//     SCIENCE_NETWORK_PLUGIN(Mine)
//
// The program is run twice. With --describe it prints its description as JSON
// (its kind, "network-analysis", and language, "cpp", added). With --run <dir>
// it reads dir/input.json and writes dir/output.json, reporting progress on
// standard output as lines "progress\t<fraction>\t<text>":
//
//   input.json, written by the host
//     {"network": {"title", "year", "unit",
//                  "places": [{"id", "name", "latitude", "longitude", "group"}],   degrees on WGS 84
//                  "flows": [{"from", "to", "value", "kind"?}],                     place ids
//                  "decomposition": {...} | null},
//      "params": {name: value}}                       defaulted from the description and typed by the host
//   output.json, written by the program
//     {"results": [result, ...]}, each one of
//       {"kind": "values", "title", "unit", "scale": "diverging"|"sequential", "note", "values": {id: number}}
//       {"kind": "arcs", "title", "unit", "note", "arcs": [{"from", "to", "value"}]}
//       {"kind": "sketch", "title", "polylines": [], "markers": [{lon, lat, colour, size}], "labels": [{lon, lat, text}], "moving": []}
//       {"kind": "table", "title", "columns": [...], "rows": [[...]], "note"}
//       {"kind": "note", "text"}
//
// A run that throws prints its reason on standard error and exits with 1; the
// host reports it. Distances are great-circle on the sphere of mean radius
// 6371.0088 km, as networks.py takes them.

#pragma once

#include "science/json.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <fstream>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace science::network {

inline constexpr double earth_radius_km = 6371.0088;

struct Place {
  std::string id, name, group;
  double latitude = 0, longitude = 0;   // degrees
};

struct Flow {
  std::size_t from = 0, to = 0;   // indices of places
  double value = 0;
  std::string kind;
};

// The places and the flows between them, as a plug-in reads them.
struct Network {
  std::string title, unit;
  std::optional<int> year;
  std::vector<Place> places;
  std::vector<Flow> flows;
  std::map<std::string, std::size_t> index;   // id -> place
  Json decomposition;                          // the application's own reading, or null

  static Network from(const Json& j) {
    Network n;
    n.title = j["title"].text();
    n.unit = j["unit"].text();
    if (j["year"].is_number()) n.year = static_cast<int>(j["year"].as_number());
    for (const auto& p : j["places"].as_array()) {
      Place q;
      q.id = id_of(p["id"]);
      q.name = p["name"].text(q.id);
      q.group = p["group"].text();
      q.latitude = p["latitude"].as_number();
      q.longitude = p["longitude"].as_number();
      if (!n.index.emplace(q.id, n.places.size()).second) throw std::invalid_argument("two places share the id " + q.id);
      n.places.push_back(std::move(q));
    }
    for (const auto& f : j["flows"].as_array()) {
      Flow g;
      g.from = n.at(id_of(f["from"]), "from");
      g.to = n.at(id_of(f["to"]), "to");
      g.value = f["value"].as_number();
      g.kind = f["kind"].text();
      n.flows.push_back(std::move(g));
    }
    n.decomposition = j["decomposition"];
    return n;
  }

  std::size_t size() const { return places.size(); }

  // F[i][j], the flow from place i to place j, flows between the same two places added.
  std::vector<std::vector<double>> matrix() const {
    std::vector<std::vector<double>> F(size(), std::vector<double>(size(), 0.0));
    for (const auto& f : flows) F[f.from][f.to] += f.value;
    return F;
  }
  std::vector<double> exports() const {
    std::vector<double> x(size(), 0.0);
    for (const auto& f : flows) x[f.from] += f.value;
    return x;
  }
  std::vector<double> imports() const {
    std::vector<double> m(size(), 0.0);
    for (const auto& f : flows) m[f.to] += f.value;
    return m;
  }
  // The net inflow of each place, imports − exports.
  std::vector<double> balance() const {
    auto b = imports();
    const auto x = exports();
    for (std::size_t i = 0; i < size(); i++) b[i] -= x[i];
    return b;
  }
  // The great-circle distance, km: atan2(|u × v|, u · v) on the unit vectors of the places.
  double distance_km(std::size_t i, std::size_t j) const {
    constexpr double deg = std::numbers::pi / 180;
    auto unit = [&](const Place& p) {
      const double f = p.latitude * deg, l = p.longitude * deg;
      return std::array<double, 3>{std::cos(f) * std::cos(l), std::cos(f) * std::sin(l), std::sin(f)};
    };
    const auto u = unit(places.at(i)), v = unit(places.at(j));
    const double c0 = u[1] * v[2] - u[2] * v[1], c1 = u[2] * v[0] - u[0] * v[2], c2 = u[0] * v[1] - u[1] * v[0];
    return earth_radius_km * std::atan2(std::sqrt(c0 * c0 + c1 * c1 + c2 * c2), u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
  }

 private:
  // An id as text: a string, or a number written as JSON writes it.
  static std::string id_of(const Json& j) { return j.is_string() ? j.as_string() : Json::number_text(j.as_number()); }
  std::size_t at(const std::string& id, const char* end) const {
    const auto it = index.find(id);
    if (it == index.end()) throw std::invalid_argument(std::string("a flow ") + end + " '" + id + "': no place has that id");
    return it->second;
  }
};

// The parameters of a run: the values the host gave, else the defaults of the description.
class Params {
 public:
  Params() = default;
  Params(Json given, Json described) : given_(std::move(given)), described_(std::move(described)) {}

  double number(const std::string& name) const { return value(name).as_number(); }
  long long integer(const std::string& name) const { return std::llround(value(name).as_number()); }
  std::string text(const std::string& name) const {
    const Json& v = value(name);
    return v.is_string() ? v.as_string() : Json::number_text(v.as_number());
  }
  bool flag(const std::string& name) const {
    const Json& v = value(name);
    return v.is_bool() ? v.as_bool() : v.is_number() ? v.as_number() != 0 : (v.text() == "true" || v.text() == "1");
  }

 private:
  Json given_, described_;
  const Json& value(const std::string& name) const {
    if (const Json* v = given_.find(name); v && !v->is_null()) return *v;
    if (described_.is_array())
      for (const auto& p : described_.as_array())
        if (p["name"].text() == name) return p["default"];
    throw std::invalid_argument("no parameter '" + name + "'");
  }
};

// Progress, as the line the host reads: "progress\t<fraction>\t<text>".
class Progress {
 public:
  void operator()(double fraction, const std::string& text = "") const {
    std::string clean;
    for (char c : text) clean += (c == '\t' || c == '\n' || c == '\r') ? ' ' : c;
    std::printf("progress\t%.4f\t%s\n", std::fmin(1.0, std::fmax(0.0, fraction)), clean.c_str());
    std::fflush(stdout);
  }
};

// ---- results -------------------------------------------------------------------------------------------

struct Arc {
  std::string from, to;
  double value = 0;
};
struct Marker {
  double lon = 0, lat = 0;
  std::string colour = "#c4506e";
  double size = 4;
};
struct Label {
  double lon = 0, lat = 0;
  std::string text;
};

// A number for each place, by id; NaN is left out. scale: "diverging" about 0, or "sequential".
inline Json values(const std::string& title, const std::map<std::string, double>& values, const std::string& unit = "",
                   const std::string& scale = "diverging", const std::string& note = "") {
  Json v = Json::object();
  for (const auto& [id, x] : values)
    if (std::isfinite(x)) v.set(id, x);
  return Json::object({{"kind", "values"}, {"title", title}, {"unit", unit}, {"scale", scale}, {"note", note}, {"values", v}});
}

inline Json arcs(const std::string& title, const std::vector<Arc>& arcs, const std::string& unit = "", const std::string& note = "") {
  Json a = Json::array();
  for (const auto& x : arcs)
    if (std::isfinite(x.value)) a.push(Json::object({{"from", x.from}, {"to", x.to}, {"value", x.value}}));
  return Json::object({{"kind", "arcs"}, {"title", title}, {"unit", unit}, {"note", note}, {"arcs", a}});
}

inline Json table(const std::string& title, const std::vector<std::string>& columns, const std::vector<std::vector<Json>>& rows,
                  const std::string& note = "") {
  Json r = Json::array();
  for (const auto& row : rows) r.push(Json(Json::Array(row.begin(), row.end())));
  return Json::object({{"kind", "table"}, {"title", title}, {"columns", Json(columns)}, {"rows", r}, {"note", note}});
}

inline Json note(const std::string& text) { return Json::object({{"kind", "note"}, {"text", text}}); }

// A sketch over the map: markers and labels in longitude and latitude, degrees.
inline Json sketch(const std::string& title, const std::vector<Marker>& markers = {}, const std::vector<Label>& labels = {}) {
  Json m = Json::array(), l = Json::array();
  for (const auto& x : markers) m.push(Json::object({{"lon", x.lon}, {"lat", x.lat}, {"colour", x.colour}, {"size", x.size}}));
  for (const auto& x : labels) l.push(Json::object({{"lon", x.lon}, {"lat", x.lat}, {"text", x.text}}));
  return Json::object({{"kind", "sketch"}, {"title", title}, {"polylines", Json::array()}, {"markers", m}, {"labels", l},
                       {"moving", Json::array()}});
}

// A parameter: kind "number", "integer", "choice" (with options) or "boolean"; null bounds are left out.
inline Json param(const std::string& name, const std::string& label, Json default_value, const std::string& kind = "number",
                  Json minimum = {}, Json maximum = {}, Json step = {}, Json options = {}, const std::string& help = "") {
  Json p = Json::object({{"name", name}, {"label", label}, {"kind", kind}, {"default", std::move(default_value)}});
  if (!minimum.is_null()) p.set("min", std::move(minimum));
  if (!maximum.is_null()) p.set("max", std::move(maximum));
  if (!step.is_null()) p.set("step", std::move(step));
  if (!options.is_null()) p.set("options", std::move(options));
  if (!help.empty()) p.set("help", help);
  return p;
}

// ---- the class and the program ---------------------------------------------------------------------------

// An analysis of a network of places: describe() gives {id, title, description, params}; run() answers results.
class Analysis {
 public:
  virtual ~Analysis() = default;
  virtual Json describe() const = 0;
  virtual std::vector<Json> run(const Network& network, const Params& params, Progress& progress) = 0;
};

inline std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error(path + ": cannot be read");
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

// The program: --describe, or --run <dir> (the contract above).
inline int main_of(Analysis& plugin, int argc, char** argv) {
  try {
    const std::string what = argc > 1 ? argv[1] : "";
    if (what == "--describe") {
      Json d = plugin.describe();
      d.set("kind", "network-analysis");
      d.set("language", "cpp");
      if (!d.has("params")) d.set("params", Json::array());
      std::printf("%s\n", d.dump().c_str());
      return 0;
    }
    if (what == "--run" && argc > 2) {
      const std::string dir = argv[2];
      const Json request = Json::parse(read_file(dir + "/input.json"));
      const Network network = Network::from(request["network"]);
      const Params params(request["params"], plugin.describe()["params"]);
      Progress progress;
      Json results(plugin.run(network, params, progress));
      std::ofstream out(dir + "/output.json", std::ios::binary);
      out << Json::object({{"results", std::move(results)}}).dump();
      if (!out) throw std::runtime_error(dir + "/output.json: cannot be written");
      return 0;
    }
    std::fprintf(stderr, "usage: %s --describe | --run <dir>\n", argc > 0 ? argv[0] : "plugin");
    return 2;
  } catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 1;
  }
}

}  // namespace science::network

#define SCIENCE_NETWORK_PLUGIN(Type)                         \
  int main(int argc, char** argv) {                          \
    Type plugin;                                             \
    return ::science::network::main_of(plugin, argc, argv);  \
  }
