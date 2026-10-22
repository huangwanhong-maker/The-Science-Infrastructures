// The JSON of the API's clients in C and C++ (api/README.md §3.5): a value
// (null, a boolean, a number, a string, an array or an object), its parser and
// its writer, on the standard library alone. It reads every answer of the
// service and writes every body the service takes; it is not meant as a general
// JSON library, and makes the choices the API needs:
//
// - Numbers are doubles. Integers are exact to 2^53, which holds every count,
//   index and event id of the API. They are read by strtod, which rounds
//   correctly, and written in the fewest significant digits (15 to 17) that read
//   back to the same double, as Python's repr writes them; so a number passes
//   through the service and the clients unchanged. The C locale's decimal point
//   is accounted for, so that a program that sets a locale still reads "1.5".
// - NaN and the infinities have no JSON: they are written as null, as the
//   service writes them, and number() reads null as NaN, so that a list of
//   values with gaps reads as doubles with NaN where there is no value.
// - Strings are UTF-8. Escapes are decoded, \uXXXX to UTF-8 with surrogate
//   pairs joined into one code point; a lone surrogate becomes U+FFFD, the
//   replacement character, rather than an error. Other bytes pass unchanged.
// - Objects keep the order of their keys, as the service writes them; a key
//   looked up that is absent reads as null. Lookup is linear, which suits
//   objects of tens of keys, as the API's are.
// - The parser follows RFC 8259 strictly (no comments, no trailing commas, no
//   leading zeros, no control characters in strings), and nests at most 512
//   levels deep, so that no answer can exhaust the stack.
//
//     science::Json j = science::Json::parse(R"({"values": [1.5, null], "units": "K"})");
//     j["values"][1].number();                          // NaN
//     science::Json body = science::Json::object({{"var", "tmp2m"}, {"lon", science::Json::array({139.69, 135.5})}});
//     body.dump();                                      // {"var":"tmp2m","lon":[139.69,135.5]}

#pragma once

#include <clocale>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace science {

// A text that is not JSON, or a value read as what it is not.
class JsonError : public std::runtime_error {
 public:
  explicit JsonError(const std::string& what) : std::runtime_error(what) {}
};

class Json {
 public:
  using Array = std::vector<Json>;
  using Member = std::pair<std::string, Json>;
  using Object = std::vector<Member>;
  enum class Type { null, boolean, number, string, array, object };

  Json() = default;
  Json(std::nullptr_t) {}
  Json(bool b) : value_(std::in_place_type<bool>, b) {}
  template <class T, std::enable_if_t<std::is_arithmetic_v<T> && !std::is_same_v<T, bool>, int> = 0>
  Json(T x) : value_(std::in_place_type<double>, static_cast<double>(x)) {}
  Json(const char* s) : value_(std::in_place_type<std::string>, s ? s : "") {}
  Json(std::string s) : value_(std::in_place_type<std::string>, std::move(s)) {}
  Json(std::string_view s) : value_(std::in_place_type<std::string>, s) {}
  Json(Array a) : value_(std::in_place_type<Array>, std::move(a)) {}
  Json(Object o) : value_(std::in_place_type<Object>, std::move(o)) {}
  // A list of numbers or strings, as an array.
  template <class T, std::enable_if_t<std::is_arithmetic_v<T> || std::is_convertible_v<const T&, std::string>, int> = 0>
  Json(const std::vector<T>& xs) : value_(std::in_place_type<Array>) {
    auto& a = std::get<Array>(value_);
    a.reserve(xs.size());
    for (const auto& x : xs) a.emplace_back(x);
  }

  static Json array(std::initializer_list<Json> xs) { return Json(Array(xs)); }
  static Json object(std::initializer_list<Member> xs) { return Json(Object(xs)); }
  static Json array() { return Json(Array{}); }
  static Json object() { return Json(Object{}); }

  // The value a text holds; JsonError, with the place, when it is not JSON.
  static Json parse(std::string_view text);
  // The value as compact JSON text.
  std::string dump() const {
    std::string out;
    write(out);
    return out;
  }

  Type type() const { return static_cast<Type>(value_.index()); }
  bool is_null() const { return type() == Type::null; }
  bool is_bool() const { return type() == Type::boolean; }
  bool is_number() const { return type() == Type::number; }
  bool is_string() const { return type() == Type::string; }
  bool is_array() const { return type() == Type::array; }
  bool is_object() const { return type() == Type::object; }

  bool as_bool() const { return std::get<bool>(checked(Type::boolean)); }
  // A number; null reads as NaN, the API's "no value".
  double as_number() const {
    if (is_null()) return std::numeric_limits<double>::quiet_NaN();
    return std::get<double>(checked(Type::number));
  }
  // A string, an array or an object. Of a temporary (the answer of a call) they
  // are moved out, so that `for (auto& f : api.fields().as_array())` holds no
  // reference into a value already gone.
  const std::string& as_string() const& { return std::get<std::string>(checked(Type::string)); }
  std::string as_string() && { return std::move(std::get<std::string>(checked(Type::string))); }
  const Array& as_array() const& { return std::get<Array>(checked(Type::array)); }
  Array& as_array() & { return std::get<Array>(checked(Type::array)); }
  Array as_array() && { return std::move(std::get<Array>(checked(Type::array))); }
  const Object& as_object() const& { return std::get<Object>(checked(Type::object)); }
  Object& as_object() & { return std::get<Object>(checked(Type::object)); }
  Object as_object() && { return std::move(std::get<Object>(checked(Type::object))); }

  // A string, or `fallback` when the value is not one (absent keys read as null).
  std::string text(const std::string& fallback = {}) const { return is_string() ? as_string() : fallback; }
  // The elements of an array, or the members of an object; 0 for any other value.
  std::size_t size() const {
    if (is_array()) return std::get<Array>(value_).size();
    if (is_object()) return std::get<Object>(value_).size();
    return 0;
  }

  // The member of that key, or nullptr.
  const Json* find(std::string_view key) const {
    if (!is_object()) return nullptr;
    for (const auto& m : std::get<Object>(value_))
      if (m.first == key) return &m.second;
    return nullptr;
  }
  bool has(std::string_view key) const { return find(key) != nullptr; }
  // The member of that key; null when it is absent or this is not an object.
  // Of a temporary, the member is moved out, as above.
  const Json& operator[](std::string_view key) const& {
    const Json* found = find(key);
    return found ? *found : none();
  }
  Json operator[](std::string_view key) && { return std::move(const_cast<Json&>((*this)[key])); }
  const Json& operator[](const char* key) const& { return (*this)[std::string_view(key)]; }
  Json operator[](const char* key) && { return std::move(*this)[std::string_view(key)]; }
  // The element at i; null beyond the end or when this is not an array.
  const Json& operator[](std::size_t i) const& {
    if (!is_array() || i >= std::get<Array>(value_).size()) return none();
    return std::get<Array>(value_)[i];
  }
  Json operator[](std::size_t i) && { return std::move(const_cast<Json&>((*this)[i])); }
  const Json& operator[](int i) const& { return i < 0 ? none() : (*this)[static_cast<std::size_t>(i)]; }
  Json operator[](int i) && { return i < 0 ? Json() : std::move(*this)[static_cast<std::size_t>(i)]; }

  // Sets a member, replacing one of the same key; a null value becomes an object first.
  Json& set(const std::string& key, Json value) {
    if (is_null()) value_ = Object{};
    for (auto& m : as_object())
      if (m.first == key) {
        m.second = std::move(value);
        return *this;
      }
    as_object().emplace_back(key, std::move(value));
    return *this;
  }
  // Appends an element; a null value becomes an array first.
  Json& push(Json value) {
    if (is_null()) value_ = Array{};
    as_array().push_back(std::move(value));
    return *this;
  }

  // Equality of values: numbers by ==, so NaN never equals; objects whatever the order of their keys.
  bool operator==(const Json& o) const {
    if (type() != o.type()) return false;
    if (!is_object()) return value_ == o.value_;
    const auto& a = std::get<Object>(value_);
    if (a.size() != o.size()) return false;
    for (const auto& m : a) {
      const Json* other = o.find(m.first);
      if (!other || !(m.second == *other)) return false;
    }
    return true;
  }
  bool operator!=(const Json& o) const { return !(*this == o); }

  // JSON for a number: the fewest digits that read back to the same double; null for NaN and the infinities.
  static std::string number_text(double x) {
    if (!std::isfinite(x)) return "null";
    char buf[40];
    for (int digits = 15; digits <= 17; digits++) {
      std::snprintf(buf, sizeof buf, "%.*g", digits, x);
      if (std::strtod(buf, nullptr) == x) break;
    }
    // The C locale may write a comma; JSON's point is a full stop.
    const char point = *std::localeconv()->decimal_point;
    if (point != '.')
      for (char* c = buf; *c; c++)
        if (*c == point) *c = '.';
    return buf;
  }

  // A string as JSON, with its quotes.
  static void quote(std::string& out, std::string_view s) {
    static const char* hex = "0123456789abcdef";
    out += '"';
    for (unsigned char c : s) {
      switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\b': out += "\\b"; break;
        case '\f': out += "\\f"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
          if (c < 0x20) {
            out += "\\u00";
            out += hex[c >> 4];
            out += hex[c & 15];
          } else {
            out += static_cast<char>(c);
          }
      }
    }
    out += '"';
  }

  void write(std::string& out) const {
    switch (type()) {
      case Type::null: out += "null"; break;
      case Type::boolean: out += std::get<bool>(value_) ? "true" : "false"; break;
      case Type::number: out += number_text(std::get<double>(value_)); break;
      case Type::string: quote(out, std::get<std::string>(value_)); break;
      case Type::array: {
        out += '[';
        bool first = true;
        for (const auto& x : std::get<Array>(value_)) {
          if (!first) out += ',';
          first = false;
          x.write(out);
        }
        out += ']';
        break;
      }
      case Type::object: {
        out += '{';
        bool first = true;
        for (const auto& m : std::get<Object>(value_)) {
          if (!first) out += ',';
          first = false;
          quote(out, m.first);
          out += ':';
          m.second.write(out);
        }
        out += '}';
        break;
      }
    }
  }

 private:
  using Value = std::variant<std::nullptr_t, bool, double, std::string, Array, Object>;
  Value value_;

  static const Json& none() {
    static const Json null;
    return null;
  }
  static const char* name(Type t) {
    static const char* names[] = {"null", "a boolean", "a number", "a string", "an array", "an object"};
    return names[static_cast<int>(t)];
  }
  const Value& checked(Type want) const {
    if (type() != want) throw JsonError(std::string("the value is ") + name(type()) + ", not " + name(want));
    return value_;
  }
  Value& checked(Type want) {
    static_cast<const Json&>(*this).checked(want);
    return value_;
  }

  class Parser;
};

// ---- the parser ------------------------------------------------------------------------------------

class Json::Parser {
 public:
  explicit Parser(std::string_view text) : s_(text) {}

  Json document() {
    Json value = parse(0);
    space();
    if (at_ < s_.size()) fail("text after the value");
    return value;
  }

 private:
  static constexpr int MAX_DEPTH = 512;
  std::string_view s_;
  std::size_t at_ = 0;

  [[noreturn]] void fail(const std::string& why) const {
    throw JsonError("not JSON: " + why + " at byte " + std::to_string(at_));
  }
  void space() {
    while (at_ < s_.size() && (s_[at_] == ' ' || s_[at_] == '\t' || s_[at_] == '\n' || s_[at_] == '\r')) at_++;
  }
  bool word(std::string_view w) {
    if (s_.substr(at_, w.size()) != w) return false;
    at_ += w.size();
    return true;
  }

  // A value inside `depth` arrays and objects.
  Json parse(int depth) {
    space();
    if (at_ >= s_.size()) fail("the text ends where a value should be");
    char c = s_[at_];
    if (c == '{') return object(depth);
    if (c == '[') return array(depth);
    if (c == '"') return Json(string());
    if (c == '-' || (c >= '0' && c <= '9')) return Json(number());
    if (word("true")) return Json(true);
    if (word("false")) return Json(false);
    if (word("null")) return Json();
    fail(std::string("unexpected '") + c + "'");
  }

  Json object(int depth) {
    if (depth >= MAX_DEPTH) fail("nested deeper than 512 levels");
    at_++;  // {
    Object members;
    space();
    if (at_ < s_.size() && s_[at_] == '}') {
      at_++;
      return Json(std::move(members));
    }
    for (;;) {
      space();
      if (at_ >= s_.size() || s_[at_] != '"') fail("a key of an object must be a string");
      std::string key = string();
      space();
      if (at_ >= s_.size() || s_[at_] != ':') fail("a key must be followed by ':'");
      at_++;
      members.emplace_back(std::move(key), parse(depth + 1));
      space();
      if (at_ < s_.size() && s_[at_] == ',') {
        at_++;
        continue;
      }
      if (at_ < s_.size() && s_[at_] == '}') {
        at_++;
        return Json(std::move(members));
      }
      fail("an object's members are separated by ',' and end with '}'");
    }
  }

  Json array(int depth) {
    if (depth >= MAX_DEPTH) fail("nested deeper than 512 levels");
    at_++;  // [
    Array items;
    space();
    if (at_ < s_.size() && s_[at_] == ']') {
      at_++;
      return Json(std::move(items));
    }
    for (;;) {
      items.push_back(parse(depth + 1));
      space();
      if (at_ < s_.size() && s_[at_] == ',') {
        at_++;
        continue;
      }
      if (at_ < s_.size() && s_[at_] == ']') {
        at_++;
        return Json(std::move(items));
      }
      fail("an array's elements are separated by ',' and end with ']'");
    }
  }

  // RFC 8259 §6: -? (0 | [1-9][0-9]*) (. [0-9]+)? ([eE] [+-]? [0-9]+)?, then strtod.
  double number() {
    std::size_t start = at_;
    auto digits = [&] {
      std::size_t from = at_;
      while (at_ < s_.size() && s_[at_] >= '0' && s_[at_] <= '9') at_++;
      return at_ - from;
    };
    if (s_[at_] == '-') at_++;
    if (at_ < s_.size() && s_[at_] == '0') {
      at_++;
      if (at_ < s_.size() && s_[at_] >= '0' && s_[at_] <= '9') fail("a number with a leading zero");
    } else if (!digits()) {
      fail("a number without digits");
    }
    if (at_ < s_.size() && s_[at_] == '.') {
      at_++;
      if (!digits()) fail("a number with no digits after its point");
    }
    if (at_ < s_.size() && (s_[at_] == 'e' || s_[at_] == 'E')) {
      at_++;
      if (at_ < s_.size() && (s_[at_] == '+' || s_[at_] == '-')) at_++;
      if (!digits()) fail("a number with no digits in its exponent");
    }
    std::string text(s_.substr(start, at_ - start));
    const char point = *std::localeconv()->decimal_point;
    if (point != '.')
      for (char& c : text)
        if (c == '.') c = point;
    return std::strtod(text.c_str(), nullptr);  // beyond the range of doubles: ±HUGE_VAL, as Python reads it
  }

  unsigned hex4() {
    if (at_ + 4 > s_.size()) fail("a \\u escape needs four hexadecimal digits");
    unsigned v = 0;
    for (int k = 0; k < 4; k++) {
      char c = s_[at_++];
      v <<= 4;
      if (c >= '0' && c <= '9') v |= unsigned(c - '0');
      else if (c >= 'a' && c <= 'f') v |= unsigned(c - 'a' + 10);
      else if (c >= 'A' && c <= 'F') v |= unsigned(c - 'A' + 10);
      else fail("a \\u escape needs four hexadecimal digits");
    }
    return v;
  }

  static void utf8(std::string& out, unsigned cp) {
    if (cp < 0x80) {
      out += static_cast<char>(cp);
    } else if (cp < 0x800) {
      out += static_cast<char>(0xC0 | (cp >> 6));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
      out += static_cast<char>(0xE0 | (cp >> 12));
      out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    } else {
      out += static_cast<char>(0xF0 | (cp >> 18));
      out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F));
      out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
      out += static_cast<char>(0x80 | (cp & 0x3F));
    }
  }

  std::string string() {
    at_++;  // "
    std::string out;
    for (;;) {
      if (at_ >= s_.size()) fail("a string without its closing quote");
      unsigned char c = static_cast<unsigned char>(s_[at_++]);
      if (c == '"') return out;
      if (c < 0x20) fail("a control character inside a string");
      if (c != '\\') {
        out += static_cast<char>(c);
        continue;
      }
      if (at_ >= s_.size()) fail("a string ends inside an escape");
      char e = s_[at_++];
      switch (e) {
        case '"': out += '"'; break;
        case '\\': out += '\\'; break;
        case '/': out += '/'; break;
        case 'b': out += '\b'; break;
        case 'f': out += '\f'; break;
        case 'n': out += '\n'; break;
        case 'r': out += '\r'; break;
        case 't': out += '\t'; break;
        case 'u': {
          unsigned cp = hex4();
          if (cp >= 0xD800 && cp <= 0xDBFF) {
            // A high surrogate joins the low one that follows (RFC 8259 §7); alone, it is U+FFFD.
            if (s_.substr(at_, 2) == "\\u") {
              std::size_t back = at_;
              at_ += 2;
              unsigned low = hex4();
              if (low >= 0xDC00 && low <= 0xDFFF) {
                cp = 0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
              } else {
                at_ = back;  // not a pair: the next escape is read on its own
                cp = 0xFFFD;
              }
            } else {
              cp = 0xFFFD;
            }
          } else if (cp >= 0xDC00 && cp <= 0xDFFF) {
            cp = 0xFFFD;
          }
          utf8(out, cp);
          break;
        }
        default: fail(std::string("an unknown escape \\") + e);
      }
    }
  }
};

inline Json Json::parse(std::string_view text) { return Parser(text).document(); }

}  // namespace science
