// JSON values (RFC 8259): read from text and written back, objects keeping the order of
// their keys. A record of the history domain is a file of JSON lines (codebook §2).
//
// A copy of the neuroscience domain's snn/json.hpp, in this domain's namespace. That header
// says it moves to a shared place when a second library needs JSON; this is the second, and
// the place is an open question of ADR 0041 (§4): a domain may not use the API's JSON, an
// edge, and no foundation yet holds formats. Until it is decided, the two copies are kept
// identical but for the namespace.

#pragma once

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace history {

class Json {
 public:
  using Array = std::vector<Json>;
  using Object = std::vector<std::pair<std::string, Json>>;

  Json() : v_(nullptr) {}
  Json(std::nullptr_t) : v_(nullptr) {}
  Json(bool b) : v_(b) {}
  Json(double d) : v_(d) {}
  Json(int i) : v_(static_cast<double>(i)) {}
  Json(long i) : v_(static_cast<double>(i)) {}
  Json(long long i) : v_(static_cast<double>(i)) {}
  Json(unsigned i) : v_(static_cast<double>(i)) {}
  Json(unsigned long i) : v_(static_cast<double>(i)) {}
  Json(unsigned long long i) : v_(static_cast<double>(i)) {}
  Json(const char* s) : v_(std::string(s)) {}
  Json(std::string s) : v_(std::move(s)) {}
  Json(std::string_view s) : v_(std::string(s)) {}
  Json(Array a) : v_(std::make_shared<Array>(std::move(a))) {}
  Json(Object o) : v_(std::make_shared<Object>(std::move(o))) {}
  // Values: a copy owns its own arrays and objects.
  Json(const Json& o) : v_(o.v_) { deepen(); }
  Json(Json&&) noexcept = default;
  Json& operator=(const Json& o) {
    if (this != &o) { v_ = o.v_; deepen(); }
    return *this;
  }
  Json& operator=(Json&&) noexcept = default;

  static Json array() { return Json(Array{}); }
  static Json object() { return Json(Object{}); }

  bool is_null() const { return std::holds_alternative<std::nullptr_t>(v_); }
  bool is_bool() const { return std::holds_alternative<bool>(v_); }
  bool is_number() const { return std::holds_alternative<double>(v_); }
  bool is_string() const { return std::holds_alternative<std::string>(v_); }
  bool is_array() const { return std::holds_alternative<std::shared_ptr<Array>>(v_); }
  bool is_object() const { return std::holds_alternative<std::shared_ptr<Object>>(v_); }

  bool as_bool() const { return get<bool>("a boolean"); }
  double as_number() const { return get<double>("a number"); }
  const std::string& as_string() const { return get<std::string>("a string"); }
  const Array& as_array() const { return *get<std::shared_ptr<Array>>("an array"); }
  Array& as_array() { return *get<std::shared_ptr<Array>>("an array"); }
  const Object& as_object() const { return *get<std::shared_ptr<Object>>("an object"); }
  Object& as_object() { return *get<std::shared_ptr<Object>>("an object"); }

  // Object access: the value of a key, or nullptr.
  const Json* find(std::string_view key) const {
    if (!is_object()) return nullptr;
    for (const auto& [k, v] : as_object()) if (k == key) return &v;
    return nullptr;
  }
  const Json& at(std::string_view key) const {
    const Json* v = find(key);
    if (!v) throw std::invalid_argument("missing key \"" + std::string(key) + "\"");
    return *v;
  }
  bool has(std::string_view key) const { return find(key) != nullptr; }
  // Sets a key, replacing it if present.
  Json& set(std::string key, Json value) {
    if (!is_object()) v_ = std::make_shared<Object>();
    for (auto& [k, v] : as_object()) if (k == key) { v = std::move(value); return v; }
    as_object().emplace_back(std::move(key), std::move(value));
    return as_object().back().second;
  }
  void push(Json value) {
    if (!is_array()) v_ = std::make_shared<Array>();
    as_array().push_back(std::move(value));
  }
  std::size_t size() const { return is_array() ? as_array().size() : is_object() ? as_object().size() : 0; }

  static Json parse(std::string_view text) {
    Parser p{text, 0};
    Json v = p.value();
    p.space();
    if (p.i != text.size()) p.fail("trailing characters");
    return v;
  }

  std::string dump(int indent = -1) const {
    std::string out;
    write(out, indent, 0);
    return out;
  }

 private:
  void deepen() {
    if (auto* a = std::get_if<std::shared_ptr<Array>>(&v_)) *a = std::make_shared<Array>(**a);
    else if (auto* o = std::get_if<std::shared_ptr<Object>>(&v_)) *o = std::make_shared<Object>(**o);
  }
  template <class T>
  const T& get(const char* what) const {
    if (const T* p = std::get_if<T>(&v_)) return *p;
    throw std::invalid_argument(std::string("JSON: expected ") + what);
  }
  template <class T>
  T& get(const char* what) {
    if (T* p = std::get_if<T>(&v_)) return *p;
    throw std::invalid_argument(std::string("JSON: expected ") + what);
  }

  static void write_string(std::string& out, const std::string& s) {
    out += '"';
    for (unsigned char c : s) {
      switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
          if (c < 0x20) { char buf[8]; std::snprintf(buf, sizeof buf, "\\u%04x", c); out += buf; }
          else out += static_cast<char>(c);
      }
    }
    out += '"';
  }

  void write(std::string& out, int indent, int depth) const {
    auto newline = [&](int d) {
      if (indent < 0) return;
      out += '\n';
      out.append(static_cast<std::size_t>(indent * d), ' ');
    };
    if (is_null()) out += "null";
    else if (is_bool()) out += as_bool() ? "true" : "false";
    else if (is_number()) {
      const double d = as_number();
      if (!std::isfinite(d)) { out += "null"; return; }
      char buf[32];
      if (d == std::floor(d) && std::fabs(d) < 9.007199254740992e15) std::snprintf(buf, sizeof buf, "%.0f", d);
      else std::snprintf(buf, sizeof buf, "%.17g", d);
      out += buf;
    } else if (is_string()) write_string(out, as_string());
    else if (is_array()) {
      const auto& a = as_array();
      out += '[';
      for (std::size_t i = 0; i < a.size(); i++) {
        if (i) out += ',';
        if (indent >= 0 && !a[i].is_number()) newline(depth + 1);
        a[i].write(out, indent, depth + 1);
      }
      if (indent >= 0 && !a.empty() && !a.back().is_number()) newline(depth);
      out += ']';
    } else {
      const auto& o = as_object();
      out += '{';
      for (std::size_t i = 0; i < o.size(); i++) {
        if (i) out += ',';
        newline(depth + 1);
        write_string(out, o[i].first);
        out += indent >= 0 ? ": " : ":";
        o[i].second.write(out, indent, depth + 1);
      }
      if (!o.empty()) newline(depth);
      out += '}';
    }
  }

  struct Parser {
    std::string_view s;
    std::size_t i;
    [[noreturn]] void fail(const char* why) const {
      throw std::invalid_argument("JSON at " + std::to_string(i) + ": " + why);
    }
    void space() { while (i < s.size() && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r' || s[i] == '\t')) i++; }
    bool take(char c) { space(); if (i < s.size() && s[i] == c) { i++; return true; } return false; }
    void expect(char c) { if (!take(c)) fail((std::string("expected '") + c + "'").c_str()); }
    Json value() {
      space();
      if (i >= s.size()) fail("unexpected end");
      const char c = s[i];
      if (c == '{') {
        i++;
        Object o;
        if (take('}')) return Json(std::move(o));
        do {
          space();
          std::string k = string();
          expect(':');
          o.emplace_back(std::move(k), value());
        } while (take(','));
        expect('}');
        return Json(std::move(o));
      }
      if (c == '[') {
        i++;
        Array a;
        if (take(']')) return Json(std::move(a));
        do a.push_back(value()); while (take(','));
        expect(']');
        return Json(std::move(a));
      }
      if (c == '"') return Json(string());
      if (s.substr(i, 4) == "true") { i += 4; return Json(true); }
      if (s.substr(i, 5) == "false") { i += 5; return Json(false); }
      if (s.substr(i, 4) == "null") { i += 4; return Json(nullptr); }
      return Json(number());
    }
    double number() {
      const std::size_t start = i;
      if (i < s.size() && (s[i] == '-' || s[i] == '+')) i++;
      while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.' || s[i] == 'e' || s[i] == 'E' || s[i] == '-' || s[i] == '+')) i++;
      if (start == i) fail("unexpected character");
      try { return std::stod(std::string(s.substr(start, i - start))); } catch (...) { fail("bad number"); }
    }
    static void utf8(std::string& out, std::uint32_t cp) {
      if (cp < 0x80) out += static_cast<char>(cp);
      else if (cp < 0x800) { out += static_cast<char>(0xC0 | (cp >> 6)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
      else if (cp < 0x10000) { out += static_cast<char>(0xE0 | (cp >> 12)); out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
      else { out += static_cast<char>(0xF0 | (cp >> 18)); out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F)); out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
    }
    std::uint32_t hex4() {
      if (i + 4 > s.size()) fail("short \\u escape");
      std::uint32_t v = 0;
      for (int k = 0; k < 4; k++) {
        const char c = s[i++];
        v <<= 4;
        if (c >= '0' && c <= '9') v |= c - '0';
        else if (c >= 'a' && c <= 'f') v |= c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v |= c - 'A' + 10;
        else fail("bad \\u escape");
      }
      return v;
    }
    std::string string() {
      if (i >= s.size() || s[i] != '"') fail("expected a string");
      i++;
      std::string out;
      while (i < s.size() && s[i] != '"') {
        char c = s[i++];
        if (c != '\\') { out += c; continue; }
        if (i >= s.size()) fail("unterminated escape");
        c = s[i++];
        switch (c) {
          case '"': out += '"'; break;
          case '\\': out += '\\'; break;
          case '/': out += '/'; break;
          case 'b': out += '\b'; break;
          case 'f': out += '\f'; break;
          case 'n': out += '\n'; break;
          case 'r': out += '\r'; break;
          case 't': out += '\t'; break;
          case 'u': {
            std::uint32_t cp = hex4();
            if (cp >= 0xD800 && cp < 0xDC00 && i + 6 <= s.size() && s[i] == '\\' && s[i + 1] == 'u') {
              i += 2;
              const std::uint32_t lo = hex4();
              cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
            }
            utf8(out, cp);
            break;
          }
          default: fail("bad escape");
        }
      }
      if (i >= s.size()) fail("unterminated string");
      i++;
      return out;
    }
  };

  std::variant<std::nullptr_t, bool, double, std::string, std::shared_ptr<Array>, std::shared_ptr<Object>> v_;
};

}  // namespace history
