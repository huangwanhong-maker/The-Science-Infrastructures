// The toolbox's requests and answers (ADR 0048): a tool reads the columns its roles name and its
// parameters, and answers blocks the page draws (the kinds are listed in ../tools.hpp).

#pragma once

#include <science/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace science::toolbox {

using Arr = Json::Array;
using Obj = Json::Object;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

// ---- reading the request -------------------------------------------------------------------------------

struct Request {
  const Json& json;
  explicit Request(const Json& j) : json(j) {}

  const Json& params() const { static const Json empty = Json(Obj{}); const Json* p = json.find("params"); return p ? *p : empty; }
  double number(const char* name, double fallback) const {
    const Json* v = params().find(name);
    return v && v->is_number() ? v->as_number() : fallback;
  }
  std::string text(const char* name, const std::string& fallback) const {
    const Json* v = params().find(name);
    return v && v->is_string() ? v->as_string() : fallback;
  }
  bool flag(const char* name, bool fallback) const {
    const Json* v = params().find(name);
    return v && v->is_bool() ? v->as_bool() : fallback;
  }
  // The column a role names (one), or every column it names.
  std::vector<std::string> role(const char* name) const {
    std::vector<std::string> out;
    const Json* roles = json.find("roles");
    const Json* v = roles ? roles->find(name) : nullptr;
    if (!v) return out;
    if (v->is_string()) out.push_back(v->as_string());
    else if (v->is_array())
      for (const auto& x : v->as_array())
        if (x.is_string()) out.push_back(x.as_string());
    return out;
  }
  std::string one(const char* name) const {
    const auto r = role(name);
    if (r.empty()) throw std::invalid_argument(std::string("choose a column for ") + name);
    return r[0];
  }
  const Json& column(const std::string& name) const {
    const Json* cols = json.find("columns");
    const Json* c = cols ? cols->find(name) : nullptr;
    if (!c || !c->is_array()) throw std::invalid_argument("the table has no column " + name);
    return *c;
  }
  // A column of numbers, NaN where a cell is empty or not a number.
  std::vector<double> numbers(const std::string& name) const {
    std::vector<double> out;
    for (const auto& x : column(name).as_array()) out.push_back(x.is_number() ? x.as_number() : NaN);
    return out;
  }
  std::vector<std::string> texts(const std::string& name) const {
    std::vector<std::string> out;
    for (const auto& x : column(name).as_array()) out.push_back(x.is_string() ? x.as_string() : x.is_number() ? Json(x.as_number()).dump() : "");
    return out;
  }
};

// Rows where every one of the columns is a number (listwise deletion), and the columns kept.
struct Complete {
  std::vector<std::size_t> rows;
  std::vector<std::vector<double>> columns;
};
inline Complete complete(const std::vector<std::vector<double>>& columns) {
  Complete C;
  if (columns.empty()) return C;
  C.columns.resize(columns.size());
  for (std::size_t i = 0; i < columns[0].size(); i++) {
    bool ok = true;
    for (const auto& c : columns) ok = ok && i < c.size() && std::isfinite(c[i]);
    if (!ok) continue;
    C.rows.push_back(i);
    for (std::size_t j = 0; j < columns.size(); j++) C.columns[j].push_back(columns[j][i]);
  }
  return C;
}

// ---- writing the answer ------------------------------------------------------------------------------

inline Json num(double x) { return std::isfinite(x) ? Json(x) : Json(nullptr); }
inline Json nums(const std::vector<double>& v) {
  Arr a;
  for (double x : v) a.push_back(num(x));
  return Json(std::move(a));
}
inline Json strs(const std::vector<std::string>& v) {
  Arr a;
  for (const auto& x : v) a.push_back(Json(x));
  return Json(std::move(a));
}
inline Json figures(const std::string& title, const std::vector<std::pair<std::string, double>>& items) {
  Arr a;
  for (const auto& [k, v] : items) a.push_back(Json(Arr{Json(k), num(v)}));
  return Json(Obj{{"kind", Json("figures")}, {"title", Json(title)}, {"items", Json(std::move(a))}});
}
inline Json table(const std::string& title, const std::vector<std::string>& columns, Arr rows, const std::string& note = "") {
  Obj o{{"kind", Json("table")}, {"title", Json(title)}, {"columns", strs(columns)}, {"rows", Json(std::move(rows))}};
  if (!note.empty()) o.push_back({"note", Json(note)});
  return Json(std::move(o));
}
struct Series {
  std::string name;
  std::vector<double> y, lower, upper;
  bool dashed = false;
};
inline Json lines(const std::string& title, const std::vector<double>& x, const std::string& xlabel, const std::string& ylabel, const std::vector<Series>& series) {
  Arr s;
  for (const auto& q : series) {
    Obj o{{"name", Json(q.name)}, {"y", nums(q.y)}};
    if (!q.lower.empty()) o.push_back({"lower", nums(q.lower)});
    if (!q.upper.empty()) o.push_back({"upper", nums(q.upper)});
    if (q.dashed) o.push_back({"dashed", Json(true)});
    s.push_back(Json(std::move(o)));
  }
  return Json(Obj{{"kind", Json("lines")}, {"title", Json(title)}, {"x", nums(x)}, {"xlabel", Json(xlabel)}, {"ylabel", Json(ylabel)}, {"series", Json(std::move(s))}});
}
inline Json bars(const std::string& title, const std::vector<std::string>& labels, const std::vector<std::pair<std::string, std::vector<double>>>& series,
                 bool stacked = false, double band = NaN) {
  Arr s;
  for (const auto& [name, values] : series) s.push_back(Json(Obj{{"name", Json(name)}, {"values", nums(values)}}));
  Obj o{{"kind", Json("bars")}, {"title", Json(title)}, {"labels", strs(labels)}, {"series", Json(std::move(s))}, {"stacked", Json(stacked)}};
  if (std::isfinite(band)) o.push_back({"band", Json(Arr{Json(-band), Json(band)})});
  return Json(std::move(o));
}
inline Json note(const std::string& text) { return Json(Obj{{"kind", Json("note")}, {"text", Json(text)}}); }
inline Json answer(const std::string& title, Arr blocks) { return Json(Obj{{"title", Json(title)}, {"blocks", Json(std::move(blocks))}}); }
// Figures that are words: a winner's name, a verdict.
inline Json words(const std::string& title, const std::vector<std::pair<std::string, std::string>>& items) {
  Arr a;
  for (const auto& [k, v] : items) a.push_back(Json(Arr{Json(k), Json(v)}));
  return Json(Obj{{"kind", Json("figures")}, {"title", Json(title)}, {"items", Json(std::move(a))}});
}
inline Json matrix(const std::string& title, const std::vector<std::string>& rows, const std::vector<std::string>& columns, const std::vector<std::vector<double>>& values) {
  Arr v;
  for (const auto& r : values) v.push_back(nums(r));
  return Json(Obj{{"kind", Json("matrix")}, {"title", Json(title)}, {"rows", strs(rows)}, {"columns", strs(columns)}, {"values", Json(std::move(v))}});
}
// Points, with a line through them if given (a fit, or the diagonal).
inline Json scatter(const std::string& title, const std::vector<double>& x, const std::vector<double>& y, const std::string& xlabel, const std::string& ylabel,
                    const std::vector<double>& line_x = {}, const std::vector<double>& line_y = {}) {
  Obj o{{"kind", Json("scatter")}, {"title", Json(title)}, {"x", nums(x)}, {"y", nums(y)}, {"xlabel", Json(xlabel)}, {"ylabel", Json(ylabel)}};
  if (!line_x.empty()) o.push_back({"line", Json(Obj{{"x", nums(line_x)}, {"y", nums(line_y)}})});
  return Json(std::move(o));
}
// Values for the map's places, by id (ISO alpha-3 for the library's economies).
inline Json map_values(const std::string& title, const std::map<std::string, double>& values, const std::string& scale, const std::string& unit = "") {
  Obj o;
  for (const auto& [k, x] : values)
    if (std::isfinite(x)) o.push_back({k, Json(x)});
  return Json(Obj{{"kind", Json("map")}, {"title", Json(title)}, {"values", Json(std::move(o))}, {"scale", Json(scale)}, {"unit", Json(unit)}});
}

// A value for each row of the table given, which the page can draw where its rows are places. A row the tool left out
// (a missing value) has no entry.
inline Json rows_of(const std::string& title, const std::vector<std::size_t>& rows, const std::vector<double>& values, const std::string& scale, const std::string& unit = "") {
  Arr at;
  for (std::size_t i : rows) at.push_back(Json(static_cast<double>(i)));
  return Json(Obj{{"kind", Json("rows")}, {"title", Json(title)}, {"rows", Json(std::move(at))}, {"values", nums(values)}, {"scale", Json(scale)}, {"unit", Json(unit)}});
}

// ---- reading what the page typed ---------------------------------------------------------------------------

// Labels as groups 0 … G − 1, in the order they first appear, and the names of the groups.
struct Groups {
  std::vector<std::size_t> of;
  std::vector<std::string> names;
};
inline Groups groups(const std::vector<std::string>& labels) {
  Groups g;
  std::map<std::string, std::size_t> index;
  for (const auto& l : labels) {
    const auto [it, fresh] = index.try_emplace(l, g.names.size());
    if (fresh) g.names.push_back(l);
    g.of.push_back(it->second);
  }
  return g;
}

inline std::vector<std::string> lines_of(const std::string& text) {
  std::vector<std::string> out;
  std::string line;
  for (char c : text + "\n") {
    if (c == '\n' || c == '\r') {
      const auto a = line.find_first_not_of(" \t"), b = line.find_last_not_of(" \t");
      if (a != std::string::npos) out.push_back(line.substr(a, b - a + 1));
      line.clear();
    } else line += c;
  }
  return out;
}

// Every number in a line, in order: "3,3  0,5" is 3 3 0 5; "-1.5e2" one number.
inline std::vector<double> numbers_in(const std::string& line) {
  std::vector<double> out;
  std::size_t i = 0;
  while (i < line.size()) {
    const char c = line[i];
    const bool starts = std::isdigit(static_cast<unsigned char>(c)) || ((c == '-' || c == '+' || c == '.') && i + 1 < line.size() &&
                        (std::isdigit(static_cast<unsigned char>(line[i + 1])) || line[i + 1] == '.'));
    if (!starts) { i++; continue; }
    std::size_t used = 0;
    try {
      out.push_back(std::stod(line.substr(i), &used));
    } catch (const std::exception&) {
      used = 1;
    }
    i += std::max<std::size_t>(used, 1);
  }
  return out;
}

// The words of a text, split at blanks, commas and the given separators.
inline std::vector<std::string> words_in(const std::string& text, const std::string& separators = "") {
  std::vector<std::string> out;
  std::string w;
  for (char c : text + " ") {
    if (std::isspace(static_cast<unsigned char>(c)) || c == ',' || separators.find(c) != std::string::npos) {
      if (!w.empty()) out.push_back(w);
      w.clear();
    } else w += c;
  }
  return out;
}

inline std::string fixed(double x, int digits = 3) {
  if (!std::isfinite(x)) return "—";
  char b[64];
  std::snprintf(b, sizeof b, "%.*f", digits, x);
  return b;
}

}  // namespace science::toolbox
