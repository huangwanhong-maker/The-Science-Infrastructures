// Plug-ins in C++ (ADR 0029): the classes of the JavaScript and Python plug-ins
// (applications/earth-atlas/src/plugins/plugin.js, api/python/science/plugins.py),
// with the same names, hooks and results, for a program that a server compiles
// and runs. Like the Python classes, the header belongs to no application: the
// hosting (api/python/science/hosting.py) compiles a plug-in with -I api/cpp.
// A plug-in derives from one of the classes, sets its identity, needs,
// parameters and outputs in its constructor, writes the hooks of its class, and
// ends with PLUGIN_MAIN(Type):
//
//   Plugin            identity, parameters, needs and run(): the general case;
//   Diagnostic        a field computed from fields, at every time of the first need: compute();
//   ForecastModel<S>  a state S advanced in time from an analysis: grid(), initial(), advance(), record(), verify();
//   DataAssimilation  a background and observations to an analysis: observations(), analyse();
//   Trajectories      parcels carried through a flow: seeds(), velocity().
//
//     #include "plugin.hpp"
//     struct Mine : plugin::Diagnostic { Mine() { id = "mine"; ... } std::map<std::string, std::vector<float>> compute(...) override; };
//     PLUGIN_MAIN(Mine)
//
// The program is run twice. With --describe it prints its description as JSON;
// with --run <dir> it reads dir/input.txt and the files it names, runs, and
// writes dir/output.txt and its files, reporting progress on standard output as
// lines "progress <fraction> <text>". Every line of both files is tab-separated,
// its first cell naming its type:
//
//   input.txt, written by the host
//     time <ms>                                 the instant shown, ms since 1970-01-01T00:00:00Z (UTC)
//     pinned <lon> <lat>                        the chosen place, degrees; absent when none is chosen
//     param <name> <value>                      a parameter, as text (a choice's value need not be a number)
//     need <name> <product> <model> <level> <ntimes> <lon0> <dlon> <nlon> <lat0> <dlat> <nlat> <periodic 0|1>
//     times <need> <ISO 8601> ...               the need's times, UTC
//     data <need> <variable> <file> [<standard_name> <units> <long_name>]
//                                               one variable: float32, ntimes × nlat × nlon, row-major,
//                                               rows north to south, NaN where there is no value
//     observations <need> <set> <title> <variable> <file>
//                                               an observation set for a need of kind "observations": the
//                                               file holds its rows, tab-separated, the first line the
//                                               column names (id, lon, lat, time, ...), an empty cell for none
//   output.txt, written by the program
//     class <Plugin|Diagnostic|ForecastModel|DataAssimilation|Trajectories>
//     start <ms>                                a forecast's start, or an analysis's time
//     grid <lon0> <dlon> <nlon> <lat0> <dlat> <nlat> <periodic>   the grid of the fields; else the first need's
//     times <ISO 8601> ...                      the times of the fields; else the first need's
//     field <name> <long_name> <standard_name> <units> <level> <difference 0|1> <file>
//                                               one output variable, float32 as above
//     vector <name> <u> <v> <units> <standard_name> <level>    a pair of the fields that forms a vector
//     show <surface|flow|lines> <name>          what the page selects when the product is shown
//     verify <output> <need> <label> <unit>     a forecast's output, and the need holding its truth
//     observation <lon> <lat> <value> <error> <id>             an observation an analysis used
//     sketch <title>                            a sketch; the polyline, marker, label and moving lines after it are its own
//     polyline <colour> <width> <alpha> <dashed 0|1> <lon lat lon lat ...>
//     marker <lon> <lat> <colour> <size>
//     label <lon> <lat> <text>
//     moving <colour> <size> <ms ms ...> <lon lat lon lat ...> a point drawn where it is at the time shown
//     table <title> <note> <column> ...         a table; the row lines after it are its own
//     row <cell> ...
//     note <text>
//     error <text>                              the run failed, for this reason
//
// The generic parts of a run that need none of the program's hooks are the
// host's, written once in Python for every language: a forecast's verification
// against the input's later times beside persistence, and an analysis's
// increment and the sketch of its observations. What calls the hooks, the loops
// of the four classes and the parcel step, is here. Units are SI; longitudes and
// latitudes are in degrees, as in the products. Only the standard library is
// used, so any C++20 compiler builds a plug-in.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace plugin {

constexpr double PI = 3.14159265358979323846, DEG = PI / 180;
// The IUGG mean radius R1 = (2a + b) / 3 of WGS 84 (Moritz 2000), as domains/geodesy.
constexpr double EARTH_RADIUS_M = 6371008.8;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

// ---- grids, needs, outputs and parameters ------------------------------------------------

struct Grid {
  double lon0 = 0, dlon = 1; int nlon = 0; double lat0 = 90, dlat = 1; int nlat = 0; bool periodic = true;
  int size() const { return nlon * nlat; }
  double lon(int column) const { return lon0 + column * dlon; }
  double lat(int row) const { return lat0 - row * dlat; }
};

// A field needed, by standard names and level; or, with kind "observations", an
// observation set carrying `variable`. `optional` lets the user choose none, and
// `distinct` asks for a source other than the needs before it.
struct Need {
  std::string name, label; std::vector<std::string> standard_names; std::string level; int times = 0;
  bool distinct = false, optional = false; std::string kind, variable;
};
inline Need observed(std::string name, std::string label, std::string variable, bool optional = false) {
  Need n; n.name = std::move(name); n.label = std::move(label); n.kind = "observations"; n.variable = std::move(variable); n.optional = optional;
  return n;
}

struct Output { std::string name, long_name, standard_name, units, level; bool difference = false; };
struct Vector { std::string name, u, v, units, standard_name, level; };

// A parameter the panel offers: a number within bounds, or one of `choices`
// ({value, label} pairs, the value text), whose default is `text`.
struct Parameter {
  std::string name, label; double value = 0, min = 0, max = 0, step = 0; std::string unit;
  std::vector<std::pair<std::string, std::string>> choices; std::string text;
};
inline Parameter choice(std::string name, std::string label, std::string value, std::vector<std::pair<std::string, std::string>> choices) {
  Parameter p; p.name = std::move(name); p.label = std::move(label); p.text = std::move(value); p.choices = std::move(choices);
  return p;
}

// A number, when the whole text is one.
inline std::optional<double> to_number(const std::string& text) {
  if (text.empty()) return std::nullopt;
  char* end = nullptr;
  const double x = std::strtod(text.c_str(), &end);
  if (end != text.c_str() + text.size()) return std::nullopt;
  return x;
}

// The values of a run's parameters: as numbers (the map), and every one as text.
struct Params : std::map<std::string, double> {
  std::map<std::string, std::string> text;
  void set(const std::string& name, const std::string& value) {
    text[name] = value;
    if (auto x = to_number(value)) (*this)[name] = *x; else erase(name);
  }
  double number(const std::string& name, double fallback = NaN) const { auto it = find(name); return it == end() ? fallback : it->second; }
  std::string choice(const std::string& name, const std::string& fallback = "") const { auto it = text.find(name); return it == text.end() ? fallback : it->second; }
};

// ---- time: ms since 1970 UTC and ISO 8601 ------------------------------------------------

// Days since 1970-01-01 of a date of the proleptic Gregorian calendar, and back
// (Hinnant, "chrono-Compatible Low-Level Date Algorithms": days_from_civil, civil_from_days).
inline long long days_from_civil(long long y, unsigned m, unsigned d) {
  y -= m <= 2;
  const long long era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = unsigned(y - era * 400);
  const unsigned doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + (long long)doe - 719468;
}
inline void civil_from_days(long long z, long long& y, unsigned& m, unsigned& d) {
  z += 719468;
  const long long era = (z >= 0 ? z : z - 146096) / 146097;
  const unsigned doe = unsigned(z - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  d = doy - (153 * mp + 2) / 5 + 1;
  m = mp < 10 ? mp + 3 : mp - 9;
  y = (long long)yoe + era * 400 + (m <= 2);
}

// An instant written in ISO 8601 (2026-10-22T15:00:00Z, with seconds, fractions
// and an offset optional), as ms since 1970 UTC; NaN when it is not one.
inline double parse_time(std::string s) {
  std::replace(s.begin(), s.end(), ' ', 'T');
  int y = 0, mo = 0, d = 0, h = 0, mi = 0; double sec = 0; char zone[16] = "";
  const int n = std::sscanf(s.c_str(), "%d-%d-%dT%d:%d:%lf%15s", &y, &mo, &d, &h, &mi, &sec, zone);
  if (n < 3 || mo < 1 || mo > 12 || d < 1 || d > 31) return NaN;
  double ms = (double(days_from_civil(y, unsigned(mo), unsigned(d))) * 86400 + h * 3600 + mi * 60 + sec) * 1000;
  if (zone[0] == '+' || zone[0] == '-') {
    int zh = 0, zm = 0;
    std::sscanf(zone + 1, "%d:%d", &zh, &zm);
    ms -= (zone[0] == '+' ? 1 : -1) * (zh * 3600 + zm * 60) * 1000.0;
  }
  return ms;
}

// ms since 1970 UTC as ISO 8601 to the second, 2026-10-22T15:00:00Z.
inline std::string iso(double ms) {
  const long long s = (long long)std::floor(ms / 1000 + 1e-6), days = (s >= 0 ? s : s - 86399) / 86400, rest = s - days * 86400;
  long long y; unsigned m, d;
  civil_from_days(days, y, m, d);
  char buf[64];
  std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02lld:%02lld:%02lldZ", y, m, d, rest / 3600, rest / 60 % 60, rest % 60);
  return buf;
}

// ---- sampling: the page's rules (foundations/physis/js/grid.js) --------------------------

// Bilinear interpolation on the grid, as the page and the Python host do it:
// columns from lon0 by dlon (wrapping round when periodic), rows from lat0
// southward by dlat; a place within half a cell of the first or last row is
// clamped onto it, and one beyond, or beyond half a cell of the edge columns of
// a grid that does not wrap, has no value; a place whose valid neighbours carry
// at least half the weight takes their weighted mean. NaN where there is none.
inline double bilinear(const float* values, const Grid& g, double lon, double lat) {
  double x = (lon - g.lon0) / g.dlon, y = (g.lat0 - lat) / g.dlat;
  if (!std::isfinite(x) || !std::isfinite(y) || y < -0.5 || y > g.nlat - 0.5) return NaN;
  if (g.periodic) x = std::fmod(std::fmod(x, g.nlon) + g.nlon, g.nlon);
  else if (x < -0.5 || x > g.nlon - 0.5) return NaN;
  y = std::clamp(y, 0.0, double(g.nlat - 1));
  int c0 = int(std::floor(x)), r0 = int(std::floor(y));
  const double fx = x - c0, fy = y - r0;
  int c1 = c0 + 1;
  const int r1 = std::min(r0 + 1, g.nlat - 1);
  if (c1 >= g.nlon) c1 = g.periodic ? c1 - g.nlon : g.nlon - 1;
  c0 = std::max(c0, 0);
  double total = 0, weight = 0;
  const int rows[4] = {r0, r0, r1, r1}, cols[4] = {c0, c1, c0, c1};
  const double ws[4] = {(1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy};
  for (int k = 0; k < 4; k++) {
    const double v = values[size_t(rows[k]) * g.nlon + cols[k]];
    if (!std::isnan(v)) { total += ws[k] * v; weight += ws[k]; }
  }
  return weight >= 0.5 ? total / weight : NaN;
}
inline double bilinear(const std::vector<float>& values, const Grid& g, double lon, double lat) { return bilinear(values.data(), g, lon, lat); }

// The two time indices around an instant and the weight of the later one,
// clamped to the span; {-1, -1, 0} for no times.
struct Bracket { int i0 = -1, i1 = -1; double w = 0; };
inline Bracket time_bracket(const std::vector<double>& times, double ms) {
  const int n = int(times.size());
  if (n == 0) return {};
  if (ms <= times.front()) return {0, 0, 0};
  if (ms >= times.back()) return {n - 1, n - 1, 0};
  const int hi = int(std::upper_bound(times.begin(), times.end(), ms) - times.begin()), lo = hi - 1;
  return {lo, hi, (ms - times[lo]) / (times[hi] - times[lo])};
}

// The index of the last time at or before `ms` (within a second), or 0.
inline int index_at_or_before(const std::vector<double>& times, double ms) {
  int found = 0;
  for (int i = 0; i < int(times.size()); i++) if (times[i] <= ms + 1000) found = i;
  return found;
}

// ---- inputs and the context ---------------------------------------------------------------

struct Variable { std::string name, standard_name, units, long_name; };

// An observation set as the library holds it (schema observations/1): its column
// names and rows, cells as text, an empty cell where there is no value.
struct Observations {
  std::vector<std::string> fields;
  std::vector<std::vector<std::string>> rows;
  int column(const std::string& name) const {
    auto it = std::find(fields.begin(), fields.end(), name);
    return it == fields.end() ? -1 : int(it - fields.begin());
  }
  const std::string& text(size_t row, int column) const {
    static const std::string none;
    return column < 0 || size_t(column) >= rows[row].size() ? none : rows[row][column];
  }
  double number(size_t row, int column) const { auto x = to_number(text(row, column)); return x ? *x : NaN; }
  double time(size_t row, int column) const { const std::string& t = text(row, column); return t.empty() ? NaN : parse_time(t); }
};

// One resolved need: a field product's grid, times and variables as float32
// blocks; or, for a need of kind "observations", the set chosen (`product` its
// name, `model` its title).
struct Input {
  std::string product, model, level, kind = "field";
  Grid grid;
  std::vector<std::string> times;            // ISO 8601, UTC
  std::vector<double> times_ms;              // the same, ms since 1970
  std::vector<std::string> names;            // the variables, in the order of the need's standard names
  std::vector<Variable> variables;
  std::vector<std::vector<float>> data;      // per variable: ntimes × nlat × nlon
  std::string variable;                      // observations: the variable asked for
  Observations observations;

  const float* slice(int variable, int index) const {
    const int last = std::max(0, int(times.size()) - 1);
    index = std::max(0, std::min(last, index));
    return data[variable].data() + size_t(index) * grid.size();
  }
  // A variable's value at a place and instant: bilinear in space, linear in time.
  double sample(int variable, double ms, double lon, double lat) const {
    const Bracket b = time_bracket(times_ms, ms);
    if (b.i0 < 0) return bilinear(slice(variable, 0), grid, lon, lat);
    const double a = bilinear(slice(variable, b.i0), grid, lon, lat);
    if (b.w == 0) return a;
    const double c = bilinear(slice(variable, b.i1), grid, lon, lat);
    return a + (c - a) * b.w;
  }
};
// The inputs by need name; an optional need left without a source is absent.
using Inputs = std::map<std::string, Input>;

struct LonLat { double lon = 0, lat = 0; };

// What a plug-in sees while it runs: the instant shown, the chosen place, and a progress hook.
struct Context {
  double time = 0;                           // ms since 1970 UTC
  std::optional<LonLat> pinned;              // the chosen place, or none
  void progress(double fraction, const std::string& text = "") const {
    std::string t = text;
    std::replace(t.begin(), t.end(), '\t', ' ');
    std::replace(t.begin(), t.end(), '\n', ' ');
    std::cout << "progress\t" << fraction << '\t' << t << std::endl;
  }
};

// ---- results ------------------------------------------------------------------------------

struct Polyline { std::vector<double> coords; std::string colour = "#e8c3dc"; double width = 1.3, alpha = 0.85; bool dashed = false; };
struct Marker { double lon = 0, lat = 0; std::string colour = "#f6c47e"; double size = 3; };
struct Label { double lon = 0, lat = 0; std::string text; };
struct Moving { std::vector<double> times, coords; std::string colour = "#f6c47e"; double size = 3; };
// A sketch over the globe: coords are lon, lat, lon, lat, ... in degrees; a moving point's times in ms.
struct Sketch { std::string title; std::vector<Polyline> polylines; std::vector<Marker> markers; std::vector<Label> labels; std::vector<Moving> moving; };
// A table for the panel: headings, then rows of cells, the first cell of each a heading.
struct Table { std::string title; std::vector<std::string> columns; std::vector<std::vector<std::string>> rows; std::string note; };
struct Observation { double lon = 0, lat = 0, value = NaN, error = NaN; std::string id; };
struct Verify { std::string output, need, label, unit; };

// The field product of a run (named plugin-<id>): its grid, times and variables.
struct Fields {
  Grid grid;
  std::vector<double> times;                                  // ms
  std::vector<Output> outputs;
  std::map<std::string, std::vector<float>> values;           // per output: ntimes × nlat × nlon
  std::vector<Vector> vectors;
  std::map<std::string, std::string> show;
};

struct Results {
  Fields fields;                             // no outputs: no field product
  std::vector<Sketch> sketches;
  std::vector<Table> tables;
  std::vector<std::string> notes;
  double start = NaN;                        // a forecast's start, or an analysis's time, ms
  Verify verify;                             // a forecast's verification, when its output is named
  std::vector<Observation> observations;     // the observations an analysis used
};

// ---- numerical helpers, as plugin.js's ----------------------------------------------------

inline double wrap_longitude(double lon) { return std::fmod(std::fmod(lon + 180, 360) + 360, 360) - 180; }

// Points spread evenly over the sphere (a Fibonacci lattice).
inline std::vector<LonLat> lattice(int n) {
  const double golden = PI * (3 - std::sqrt(5.0));
  std::vector<LonLat> out;
  for (int i = 0; i < n; i++) {
    const double z = 1 - (2.0 * i + 1) / n;
    out.push_back({wrap_longitude(i * golden / DEG), std::asin(z) / DEG});
  }
  return out;
}

// A variable's values at one time, sampled bilinearly onto another grid.
inline std::vector<float> sample_on_grid(const Input& input, int variable, int index, const Grid& grid) {
  std::vector<float> out(grid.size());
  const float* s = input.slice(variable, index);
  for (int r = 0; r < grid.nlat; r++)
    for (int c = 0; c < grid.nlon; c++) out[size_t(r) * grid.nlon + c] = float(bilinear(s, input.grid, grid.lon(c), grid.lat(r)));
  return out;
}

// A flow: (lon, lat, ms) -> {u, v} in m/s, or none outside the data.
using Flow = std::function<std::optional<std::array<double, 2>>(double lon, double lat, double ms)>;

// One fourth-order Runge–Kutta step of a parcel on the sphere, from (lon, lat) in
// degrees at time ms, over dt seconds, through `flow`: dλ/dt = u / (a cos φ),
// dφ/dt = v / a. Returns the new place, or none when the flow has no value on the
// way; a parcel carried over a pole comes down the other side (plugin.js's parcelStep).
inline std::optional<LonLat> parcel_step(const Flow& flow, double lon, double lat, double ms, double dt) {
  const auto rate = [&](double lo, double la, double t) -> std::optional<std::array<double, 2>> {
    const auto uv = flow(wrap_longitude(lo), std::max(-89.99, std::min(89.99, la)), t);
    if (!uv) return std::nullopt;
    return std::array<double, 2>{(*uv)[0] / (EARTH_RADIUS_M * std::cos(la * DEG)) / DEG, (*uv)[1] / EARTH_RADIUS_M / DEG};
  };
  const auto k1 = rate(lon, lat, ms); if (!k1) return std::nullopt;
  const auto k2 = rate(lon + 0.5 * dt * (*k1)[0], lat + 0.5 * dt * (*k1)[1], ms + 500 * dt); if (!k2) return std::nullopt;
  const auto k3 = rate(lon + 0.5 * dt * (*k2)[0], lat + 0.5 * dt * (*k2)[1], ms + 500 * dt); if (!k3) return std::nullopt;
  const auto k4 = rate(lon + dt * (*k3)[0], lat + dt * (*k3)[1], ms + 1000 * dt); if (!k4) return std::nullopt;
  double nlon = lon + dt / 6 * ((*k1)[0] + 2 * (*k2)[0] + 2 * (*k3)[0] + (*k4)[0]);
  double nlat = lat + dt / 6 * ((*k1)[1] + 2 * (*k2)[1] + 2 * (*k3)[1] + (*k4)[1]);
  if (nlat > 90) { nlat = 180 - nlat; nlon += 180; } else if (nlat < -90) { nlat = -180 - nlat; nlon += 180; }
  return LonLat{wrap_longitude(nlon), nlat};
}

// ---- the classes ----------------------------------------------------------------------------

inline std::string quoted(const std::string& s);

struct Plugin {
  std::string id = "plugin", title = "A plug-in", description, domain = "meteorology";
  std::vector<Need> needs;
  std::vector<Parameter> parameters;
  // The context of the run, for hooks whose signature has none (Diagnostic::compute).
  Context context;
  virtual ~Plugin() = default;
  virtual std::string kind() const { return "Plugin"; }
  // The algorithm, in the general case: any results.
  virtual Results run(const Inputs&, const Params&, const Context&) { return {}; }
  // The members of the description that a class adds (JSON, each preceded by a comma).
  virtual std::string described() const { return ""; }
  std::string product_name() const { return "plugin-" + id; }
  const Input& first(const Inputs& inputs) const {
    if (needs.empty()) throw std::runtime_error(id + ": no needs, and the class needs a first one");
    auto it = inputs.find(needs[0].name);
    if (it == inputs.end() || it->second.kind != "field") throw std::runtime_error(id + ": the first need, " + needs[0].name + ", is not a field given");
    return it->second;
  }
};

inline std::string described_outputs(const std::vector<Output>& outputs, const std::vector<Vector>& vectors, const std::map<std::string, std::string>& show) {
  std::ostringstream o;
  o << ",\"outputs\":[";
  for (size_t i = 0; i < outputs.size(); i++) {
    const auto& v = outputs[i];
    o << (i ? "," : "") << "{\"name\":" << quoted(v.name) << ",\"long_name\":" << quoted(v.long_name)
      << ",\"standard_name\":" << quoted(v.standard_name) << ",\"units\":" << quoted(v.units);
    if (!v.level.empty()) o << ",\"level\":" << quoted(v.level);
    if (v.difference) o << ",\"difference\":true";
    o << "}";
  }
  o << "],\"vectors\":[";
  for (size_t i = 0; i < vectors.size(); i++) {
    const auto& v = vectors[i];
    o << (i ? "," : "") << "{\"name\":" << quoted(v.name) << ",\"u\":" << quoted(v.u) << ",\"v\":" << quoted(v.v)
      << ",\"units\":" << quoted(v.units) << ",\"standard_name\":" << quoted(v.standard_name) << ",\"level\":" << quoted(v.level) << "}";
  }
  o << "],\"show\":{";
  size_t k = 0;
  for (const auto& [key, value] : show) o << (k++ ? "," : "") << quoted(key) << ":" << quoted(value);
  o << "}";
  return o.str();
}

// A field computed from fields, at every time of the first need, on its grid.
struct Diagnostic : Plugin {
  std::vector<Output> outputs;
  std::vector<Vector> vectors;
  std::map<std::string, std::string> show;   // what to select: {"surface", name}, {"flow", vector}, {"lines", name}
  std::string kind() const override { return "Diagnostic"; }
  // {output name: values on the first need's grid} for time `index` of the first need.
  virtual std::map<std::string, std::vector<float>> compute(const Inputs& inputs, const Params& params, int index) = 0;
  std::string described() const override { return described_outputs(outputs, vectors, show); }

  Results run(const Inputs& inputs, const Params& params, const Context& context) override {
    const Input& first = this->first(inputs);
    const int ntimes = std::max<int>(1, int(first.times.size()));
    Results r;
    Fields& f = r.fields;
    f.grid = first.grid; f.times = first.times_ms; f.outputs = outputs; f.vectors = vectors; f.show = show;
    for (auto& o : f.outputs) if (o.level.empty()) o.level = first.level;
    if (f.show.empty() && !outputs.empty()) f.show = {{"surface", outputs[0].name}};
    for (const auto& o : outputs) f.values[o.name].reserve(size_t(ntimes) * first.grid.size());
    for (int k = 0; k < ntimes; k++) {
      auto values = compute(inputs, params, k);
      for (const auto& o : outputs) {
        auto it = values.find(o.name);
        if (it == values.end() || int(it->second.size()) != first.grid.size()) throw std::runtime_error("compute() did not return " + o.name + " on the grid");
        f.values[o.name].insert(f.values[o.name].end(), it->second.begin(), it->second.end());
      }
      context.progress(double(k + 1) / ntimes, "time " + std::to_string(k + 1) + " of " + std::to_string(ntimes));
    }
    return r;
  }
};

// What every forecast model shows of itself, whatever its state.
struct Forecast : Plugin {
  double dt = 1800;      // the step, in seconds
  double every = 6;      // hours between recorded times
  std::vector<Output> outputs;
  std::vector<Vector> vectors;
  std::map<std::string, std::string> show;
  Forecast() { parameters = {{"hours", "Forecast length", 72, 6, 240, 6, "h"}}; }
  std::string kind() const override { return "ForecastModel"; }
  // The grid of record(): called after initial(), so that it may be the input's.
  virtual Grid grid() const = 0;
  // The output to verify, the need holding the truth, a label and a unit; an empty output for none.
  virtual Verify verify() const { return {}; }
  std::string described() const override {
    std::ostringstream o;
    o << described_outputs(outputs, vectors, show) << ",\"dt\":" << dt << ",\"every\":" << every;
    const Verify v = verify();
    if (!v.output.empty())
      o << ",\"verify\":{\"output\":" << quoted(v.output) << ",\"need\":" << quoted(v.need) << ",\"label\":" << quoted(v.label) << ",\"unit\":" << quoted(v.unit) << "}";
    return o.str();
  }
};

// A state advanced in time from an analysis: the first need gives the start, at
// its time at or before the instant shown; the state is recorded at the start and
// every `every` hours. The host verifies the recorded output against the input's
// later times beside persistence (RMSE over 20–80° N, weighted by the cosine of
// latitude), as plugins.py's ForecastModel.run does, and tables the result.
// Steps and records are counted with rounding half to even, as Python's round().
template <class State>
struct ForecastModel : Forecast {
  // The state at time index `index` of the inputs.
  virtual State initial(const Inputs& inputs, int index, const Params& params, const Context& context) = 0;
  // One step of `seconds`.
  virtual void advance(State& state, double seconds, const Params& params) = 0;
  // {output name: values on grid()}: the state as fields.
  virtual std::map<std::string, std::vector<float>> record(const State& state, const Params& params) = 0;

  Results run(const Inputs& inputs, const Params& params, const Context& context) override {
    const Input& first = this->first(inputs);
    if (first.times_ms.empty()) throw std::runtime_error(id + ": the first need has no times to start from");
    const int i0 = index_at_or_before(first.times_ms, context.time);
    const double t0 = first.times_ms[i0];
    const double hours = params.number("hours", 72);
    const long steps = long(std::nearbyint(hours * 3600 / dt));
    const long record_every = std::max(1L, long(std::nearbyint(every * 3600 / dt)));
    Results r;
    Fields& f = r.fields;
    State state = initial(inputs, i0, params, context);
    f.grid = grid();
    const auto take = [&]() {
      auto values = record(state, params);
      for (const auto& o : outputs) {
        auto it = values.find(o.name);
        if (it == values.end() || int(it->second.size()) != f.grid.size()) throw std::runtime_error("record() did not return " + o.name + " on grid()");
        auto& all = f.values[o.name];
        all.insert(all.end(), it->second.begin(), it->second.end());
      }
    };
    take();
    f.times.push_back(t0);
    for (long s = 1; s <= steps; s++) {
      advance(state, dt, params);
      if (s % record_every == 0) { take(); f.times.push_back(t0 + s * dt * 1000); }
      if (s % 4 == 0) context.progress(double(s) / steps, "+" + std::to_string(long(std::lround(s * dt / 3600))) + " h");
    }
    f.outputs = outputs; f.vectors = vectors; f.show = show;
    if (f.show.empty() && !outputs.empty()) f.show = {{"surface", outputs[0].name}};
    r.start = t0;
    r.verify = verify();
    return r;
  }
};

// A background field and observations to an analysis, at the instant shown
// (clamped to the background's span), the background interpolated linearly in
// time to it. The program returns the analysis and the observations used; the
// host adds the increment (marked a difference) and the sketch of the
// observations, warm where above the background, as plugins.py does.
struct DataAssimilation : Plugin {
  std::string kind() const override { return "DataAssimilation"; }
  // What is observed at time `t` (ms), the analysis's time.
  virtual std::vector<Observation> observations(const Inputs&, double /* t */, const Params&, const Context&) { return {}; }
  // The analysis on the background's grid (NaN where the background is missing).
  virtual std::vector<float> analyse(const std::vector<float>& background, const Grid&, const std::vector<Observation>&, const Params&, const Context&) { return background; }

  Results run(const Inputs& inputs, const Params& params, const Context& context) override {
    const Input& first = this->first(inputs);
    const auto& times = first.times_ms;
    const double t = times.empty() ? context.time : std::max(times.front(), std::min(times.back(), context.time));
    const Bracket b = time_bracket(times, t);
    const float* a0 = first.slice(0, std::max(0, b.i0));
    const float* a1 = first.slice(0, std::max(0, b.i1));
    std::vector<float> background(first.grid.size());
    for (size_t i = 0; i < background.size(); i++) background[i] = b.w > 0 ? float(double(a0[i]) + (double(a1[i]) - a0[i]) * b.w) : a0[i];
    Results r;
    r.observations = observations(inputs, t, params, context);
    std::vector<float> analysis = analyse(background, first.grid, r.observations, params, context);
    if (int(analysis.size()) != first.grid.size()) throw std::runtime_error("analyse() did not return a field on the background's grid");
    const Variable v = first.variables.empty() ? Variable{first.names.empty() ? "x" : first.names[0], "", "", ""} : first.variables[0];
    r.fields.grid = first.grid;
    r.fields.times = {t};
    r.fields.outputs = {{v.name, (v.long_name.empty() ? v.name : v.long_name) + ", analysed", v.standard_name, v.units, first.level, false}};
    r.fields.values[v.name] = std::move(analysis);
    r.start = t;
    return r;
  }
};

// Parcels carried through a flow: the first need is a vector field. Override
// seeds(); velocity() samples the need by default. The result is a sketch of the
// paths, with the parcels moving along them as the time shown changes.
struct Seed { double lon = 0, lat = 0; std::string colour; };

struct Trajectories : Plugin {
  std::string colour = "#e8c3dc";
  Trajectories() {
    parameters = {{"hours", "Hours ahead", 48, 6, 240, 6, "h"}, {"step", "Time step", 1, 0.25, 6, 0.25, "h"}};
  }
  std::string kind() const override { return "Trajectories"; }
  // Where the parcels start.
  virtual std::vector<Seed> seeds(const Inputs&, const Params&, const Context&) { return {}; }
  // The flow: the first need's two variables, bilinear in space and linear in time.
  virtual Flow velocity(const Inputs& inputs, const Params&) {
    const Input* w = &this->first(inputs);
    return [w](double lon, double lat, double ms) -> std::optional<std::array<double, 2>> {
      const double a = w->sample(0, ms, lon, lat), b = w->sample(1, ms, lon, lat);
      if (std::isnan(a) || std::isnan(b)) return std::nullopt;
      return std::array<double, 2>{a, b};
    };
  }

  Results run(const Inputs& inputs, const Params& params, const Context& context) override {
    const Input& w = this->first(inputs);
    const auto& times = w.times_ms;
    if (times.empty()) throw std::runtime_error(id + ": the flow has no times");
    const Flow flow = velocity(inputs, params);
    const double t_end = times.back(), t0 = std::max(times.front(), std::min(context.time, t_end));
    const double dt = params.number("step", 1) * 3600;
    const long steps = std::max(0L, std::min(long(std::nearbyint(params.number("hours", 48) * 3600 / dt)), long(std::floor((t_end - t0) / 1000 / dt))));
    const std::vector<Seed> starts = seeds(inputs, params, context);
    Sketch sketch;
    for (size_t s = 0; s < starts.size(); s++) {
      const Seed& seed = starts[s];
      Polyline path; Moving point;
      path.coords = {seed.lon, seed.lat}; point.times = {t0};
      double lon = seed.lon, lat = seed.lat;
      for (long k = 0; k < steps; k++) {
        const auto next = parcel_step(flow, lon, lat, t0 + k * dt * 1000, dt);
        if (!next) break;
        lon = next->lon; lat = next->lat;
        path.coords.push_back(lon); path.coords.push_back(lat);
        point.times.push_back(t0 + (k + 1) * dt * 1000);
      }
      path.colour = seed.colour.empty() ? colour : seed.colour;
      point.coords = path.coords;
      sketch.polylines.push_back(std::move(path));
      sketch.moving.push_back(std::move(point));
      if (s % 8 == 7) context.progress(double(s + 1) / starts.size(), "parcel " + std::to_string(s + 1) + " of " + std::to_string(starts.size()));
    }
    sketch.title = std::to_string(starts.size()) + " parcels through the wind of " + w.model + ", from " + iso(t0) + " for " +
                   std::to_string(long(std::lround(steps * dt / 3600))) + " h";
    Results r;
    r.sketches.push_back(std::move(sketch));
    return r;
  }
};

// ---- the description, as JSON, for --describe ------------------------------------------------

inline std::string quoted(const std::string& s) {
  std::string out = "\"";
  for (unsigned char c : s) {
    if (c == '"') out += "\\\""; else if (c == '\\') out += "\\\\"; else if (c == '\n') out += "\\n";
    else if (c < 0x20) { char buf[8]; std::snprintf(buf, sizeof buf, "\\u%04x", c); out += buf; } else out += char(c);
  }
  return out + "\"";
}

// A number as JSON and in output.txt: %.15g (times in ms exactly), and null for NaN or an infinity.
inline std::string number(double x) {
  if (!std::isfinite(x)) return "null";
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.15g", x);
  return buf;
}

// A choice's value: a number when the text is one, else a string.
inline std::string value_json(const std::string& text) { auto x = to_number(text); return x ? number(*x) : quoted(text); }

inline std::string describe(const Plugin& p) {
  std::ostringstream o;
  o << "{\"id\":" << quoted(p.id) << ",\"title\":" << quoted(p.title) << ",\"description\":" << quoted(p.description)
    << ",\"domain\":" << quoted(p.domain) << ",\"class\":" << quoted(p.kind()) << ",\"parameters\":{";
  for (size_t i = 0; i < p.parameters.size(); i++) {
    const auto& q = p.parameters[i];
    o << (i ? "," : "") << quoted(q.name) << ":{\"label\":" << quoted(q.label);
    if (!q.choices.empty()) {
      o << ",\"default\":" << value_json(q.text) << ",\"choices\":[";
      for (size_t k = 0; k < q.choices.size(); k++) o << (k ? "," : "") << "[" << value_json(q.choices[k].first) << "," << quoted(q.choices[k].second) << "]";
      o << "]}";
    } else {
      o << ",\"default\":" << number(q.value) << ",\"min\":" << number(q.min) << ",\"max\":" << number(q.max)
        << ",\"step\":" << number(q.step) << ",\"unit\":" << quoted(q.unit) << "}";
    }
  }
  o << "},\"needs\":[";
  for (size_t i = 0; i < p.needs.size(); i++) {
    const auto& n = p.needs[i];
    o << (i ? "," : "") << "{\"name\":" << quoted(n.name) << ",\"label\":" << quoted(n.label);
    if (n.kind == "observations") o << ",\"kind\":\"observations\",\"variable\":" << quoted(n.variable);
    else {
      o << ",\"standard_names\":[";
      for (size_t k = 0; k < n.standard_names.size(); k++) o << (k ? "," : "") << quoted(n.standard_names[k]);
      o << "]";
      if (!n.level.empty()) o << ",\"level\":" << quoted(n.level);
      if (n.times) o << ",\"times\":" << n.times;
      if (n.distinct) o << ",\"distinct\":true";
    }
    if (n.optional) o << ",\"optional\":true";
    o << "}";
  }
  o << "]" << p.described() << "}\n";
  return o.str();
}

// ---- the files of a run ---------------------------------------------------------------------

inline std::vector<std::string> split(const std::string& line, char by = '\t') {
  std::vector<std::string> out;
  std::string part;
  std::istringstream in(line);
  while (std::getline(in, part, by)) out.push_back(part);
  if (!line.empty() && line.back() == by) out.push_back("");
  return out;
}

inline std::vector<float> read_f32(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  in.seekg(0, std::ios::end);
  std::vector<float> out(size_t(in.tellg()) / sizeof(float));
  in.seekg(0);
  in.read(reinterpret_cast<char*>(out.data()), std::streamsize(out.size() * sizeof(float)));
  return out;
}

inline void write_f32(const std::string& path, const std::vector<float>& v) {
  std::ofstream out(path, std::ios::binary);
  out.write(reinterpret_cast<const char*>(v.data()), std::streamsize(v.size() * sizeof(float)));
  if (!out) throw std::runtime_error("cannot write " + path);
}

inline Observations read_observations(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  Observations set;
  std::string line;
  if (std::getline(in, line)) set.fields = split(line);
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto cells = split(line);
    cells.resize(set.fields.size());
    set.rows.push_back(std::move(cells));
  }
  return set;
}

struct Run { Inputs inputs; Params params; Context context; std::string first; };

// The run's inputs from dir/input.txt; parameters not given take the plug-in's defaults.
inline Run load(const std::string& dir, const Plugin& p) {
  Run run;
  std::ifstream in(dir + "/input.txt");
  if (!in) throw std::runtime_error("no input.txt in " + dir);
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto f = split(line);
    if (f.empty()) continue;
    if (f[0] == "time" && f.size() > 1) run.context.time = std::stod(f[1]);
    else if (f[0] == "pinned" && f.size() > 2) run.context.pinned = LonLat{std::stod(f[1]), std::stod(f[2])};
    else if (f[0] == "param" && f.size() > 2) run.params.set(f[1], f[2]);
    else if (f[0] == "need" && f.size() > 12) {
      Input& inp = run.inputs[f[1]];
      if (run.first.empty()) run.first = f[1];
      inp.product = f[2]; inp.model = f[3]; inp.level = f[4];
      inp.grid = {std::stod(f[6]), std::stod(f[7]), std::stoi(f[8]), std::stod(f[9]), std::stod(f[10]), std::stoi(f[11]), f[12] == "1"};
    } else if (f[0] == "times" && f.size() > 1) {
      Input& inp = run.inputs[f[1]];
      inp.times.assign(f.begin() + 2, f.end());
      inp.times_ms.clear();
      for (const auto& t : inp.times) inp.times_ms.push_back(parse_time(t));
    } else if (f[0] == "data" && f.size() > 3) {
      Input& inp = run.inputs[f[1]];
      inp.names.push_back(f[2]);
      inp.variables.push_back({f[2], f.size() > 4 ? f[4] : "", f.size() > 5 ? f[5] : "", f.size() > 6 ? f[6] : ""});
      inp.data.push_back(read_f32(dir + "/" + f[3]));
    } else if (f[0] == "observations" && f.size() > 5) {
      Input& inp = run.inputs[f[1]];
      inp.kind = "observations"; inp.product = f[2]; inp.model = f[3]; inp.variable = f[4];
      inp.observations = read_observations(dir + "/" + f[5]);
    }
  }
  for (const auto& q : p.parameters) {
    if (run.params.text.count(q.name)) continue;
    if (!q.choices.empty()) run.params.set(q.name, q.text);
    else { run.params[q.name] = q.value; run.params.text[q.name] = number(q.value); }
  }
  return run;
}

// A cell of output.txt: tabs and line breaks become spaces.
inline std::string cell(std::string s) {
  for (char& c : s) if (c == '\t' || c == '\n' || c == '\r') c = ' ';
  return s;
}

inline std::string joined(const std::vector<double>& xs) {
  std::string out;
  for (size_t i = 0; i < xs.size(); i++) { if (i) out += ' '; out += std::isfinite(xs[i]) ? number(xs[i]) : "NaN"; }
  return out;
}

// The results as dir/output.txt and its float32 files.
inline void write(const std::string& dir, const Plugin& p, const Results& r) {
  std::ofstream out(dir + "/output.txt");
  out << "class\t" << p.kind() << '\n';
  if (std::isfinite(r.start)) out << "start\t" << number(r.start) << '\n';
  const Fields& f = r.fields;
  if (!f.outputs.empty()) {
    const Grid& g = f.grid;
    out << "grid\t" << number(g.lon0) << '\t' << number(g.dlon) << '\t' << g.nlon << '\t' << number(g.lat0) << '\t'
        << number(g.dlat) << '\t' << g.nlat << '\t' << (g.periodic ? 1 : 0) << '\n';
    out << "times";
    for (double t : f.times) out << '\t' << iso(t);
    out << '\n';
    const size_t expected = size_t(std::max<size_t>(1, f.times.size())) * g.size();
    for (size_t k = 0; k < f.outputs.size(); k++) {
      const Output& o = f.outputs[k];
      auto it = f.values.find(o.name);
      if (it == f.values.end() || it->second.size() != expected) throw std::runtime_error("no values, or not times × grid of them, for " + o.name);
      const std::string file = "out." + std::to_string(k) + ".f32";
      write_f32(dir + "/" + file, it->second);
      out << "field\t" << cell(o.name) << '\t' << cell(o.long_name) << '\t' << cell(o.standard_name) << '\t' << cell(o.units) << '\t'
          << cell(o.level) << '\t' << (o.difference ? 1 : 0) << '\t' << file << '\n';
    }
    for (const auto& v : f.vectors)
      out << "vector\t" << cell(v.name) << '\t' << cell(v.u) << '\t' << cell(v.v) << '\t' << cell(v.units) << '\t' << cell(v.standard_name) << '\t' << cell(v.level) << '\n';
    for (const auto& [key, value] : f.show) out << "show\t" << cell(key) << '\t' << cell(value) << '\n';
  }
  if (!r.verify.output.empty())
    out << "verify\t" << cell(r.verify.output) << '\t' << cell(r.verify.need) << '\t' << cell(r.verify.label) << '\t' << cell(r.verify.unit) << '\n';
  for (const auto& o : r.observations)
    out << "observation\t" << number(o.lon) << '\t' << number(o.lat) << '\t' << number(o.value) << '\t' << number(o.error) << '\t' << cell(o.id) << '\n';
  for (const auto& s : r.sketches) {
    out << "sketch\t" << cell(s.title) << '\n';
    for (const auto& l : s.polylines)
      out << "polyline\t" << cell(l.colour) << '\t' << number(l.width) << '\t' << number(l.alpha) << '\t' << (l.dashed ? 1 : 0) << '\t' << joined(l.coords) << '\n';
    for (const auto& m : s.markers) out << "marker\t" << number(m.lon) << '\t' << number(m.lat) << '\t' << cell(m.colour) << '\t' << number(m.size) << '\n';
    for (const auto& l : s.labels) out << "label\t" << number(l.lon) << '\t' << number(l.lat) << '\t' << cell(l.text) << '\n';
    for (const auto& m : s.moving) out << "moving\t" << cell(m.colour) << '\t' << number(m.size) << '\t' << joined(m.times) << '\t' << joined(m.coords) << '\n';
  }
  for (const auto& t : r.tables) {
    out << "table\t" << cell(t.title) << '\t' << cell(t.note);
    for (const auto& c : t.columns) out << '\t' << cell(c);
    out << '\n';
    for (const auto& row : t.rows) {
      out << "row";
      for (const auto& c : row) out << '\t' << cell(c);
      out << '\n';
    }
  }
  for (const auto& n : r.notes) out << "note\t" << cell(n) << '\n';
  if (!out) throw std::runtime_error("cannot write " + dir + "/output.txt");
}

template <class T>
int run_main(int argc, char** argv) {
  T p;
  try {
    if (argc >= 2 && std::string(argv[1]) == "--describe") { std::cout << describe(p); return 0; }
    if (argc >= 3 && std::string(argv[1]) == "--run") {
      const std::string dir = argv[2];
      Run run = load(dir, p);
      p.context = run.context;
      const Results r = p.run(run.inputs, run.params, run.context);
      write(dir, p, r);
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cerr << "usage: " << argv[0] << " --describe | --run <dir>\n";
  return 2;
}

}  // namespace plugin

#define PLUGIN_MAIN(Type) int main(int argc, char** argv) { return plugin::run_main<Type>(argc, argv); }
