// Addresses: the language in which a point of the Solar System, or of the Galaxy, is written at an
// instant, read back, and answered in any frame and chart of the catalogue.
//
// An address names a point by the frame it is given in, a chart of that frame, its coordinates, a
// chain of local steps from it, and optionally an instant:
//
//   address  ::= point [ "@" instant ]
//   point    ::= base { "/" step }
//   base     ::= BODY [ "(" coords ")" ]                    a body's centre, or a point on its surface
//              | FRAME [ ":" CHART ] [ "(" coords ")" ]     the frame's origin, or coordinates in a chart
//   step     ::= ( "enu" | "ned" | "aer" ) "(" coords ")"   from the local horizon at the point so far
//   coords   ::= coord { "," coord }
//   coord    ::= STRING                                     a word, as an MGRS reference
//              | NUMBER [ UNIT | DMS | HMS ] [ HEMISPHERE ]
//   DMS      ::= "°" [ NUMBER "′" [ NUMBER "″" ] ]
//   HMS      ::= "h" [ NUMBER "m" [ NUMBER "s" ] ]
//   instant  ::= the rest of the text, in metrology's language of times ("2026-10-28T12:00:00Z",
//                "2026-10-28T12:00:00 TDB", "JD 2461000.5 TT", "MJD 61000 TAI"), or "now"
//
// The terminals are lexical classes: BODY a body's key (sun, mercury, …, earth, moon, …, pluto, and
// ssb, the barycentre, and galaxy, the Galaxy's centre), FRAME any other identifier (a frame of the
// catalogue, case not significant: ITRS, MOON_ME, IAU_MARS, GALACTOCENTRIC), CHART an identifier, UNIT
// one of m km au ly pc kpc Mpc deg rad arcmin arcsec, HEMISPHERE one of N S E W. Every alternative
// begins with a distinct terminal, so the grammar is LL(1) and every address has one parse.
//
//   earth(48.8583701 N, 2.2944813 E, 330 m)                        the Eiffel Tower's top, WGS 84
//   ITRS:geodetic(48.8583701 N, 2.2944813 E, 330 m) @ 2026-10-28T12:00:00Z
//   ITRS:mgrs("31UDQ4825211954", 330 m)
//   moon(43.31 S, 11.36 W)                                         Tycho, on the Moon's sphere
//   earth(35.36 N, 138.73 E)/aer(45°, 30°, 10 km)                  a point seen from Mount Fuji
//   ICRS:equatorial(266.4051°, 28.936175° S, 8.122 kpc)            the Galaxy's centre
//   GALACTOCENTRIC:cylindrical(8.2 kpc, 180°, 20.8 pc)
//   mars                                                           the centre of Mars
//   ICRS:equatorial(6 h 45 m 8.9 s, 16° 42′ 58″ S)                 the direction of Sirius (no distance)
//
// The canonical text of an address (print) writes each number in its shortest form that reads back
// to the same double (std::to_chars), each unit and hemisphere once, single spaces, frames in
// capitals; print(parse(w)) = w for every canonical w, and parse(print(a)) = a for every address.
// It is the language of STC-S (IVOA, Rots 2007; the Space-Time Coordinates of the Virtual
// Observatory) in spirit: a frame, a reference position and coordinates, in one line.
//
// Charts (coordinates on a frame): cartesian (x, y, z); spherical, and its names equatorial (right
// ascension and declination: frames with equatorial axes), ecliptic (ECLIPTIC_J2000) and galactic
// (GALACTIC), as (longitude, latitude, distance), the distance optional for a direction; cylindrical
// (ρ, φ, z); on a body's fixed frame, geodetic (WGS 84 in the ITRS) or planetocentric (the body's
// mean sphere, Archinal et al. 2018; selenographic its name on the Moon), as (latitude, longitude,
// height), a hemisphere letter fixing which is which, so that the order is free when they are
// given; and in the ITRS utm (zone, easting, northing, height) and mgrs (reference, height). Steps:
// enu (east, north, up), ned (north, east, down), aer (azimuth clockwise from north, elevation,
// range) from the horizon of the body's figure at the point so far.
//
// Resolution is passive (Cosmopoiesis formalism §6.1): a chart's inverse gives Cartesian
// coordinates in the frame; the frame graph carries them to another frame at the instant; the
// target chart gives coordinates again. Positions are geometric (no light time, no aberration).

#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>
#include <geodesy/projection/mgrs.hpp>
#include <geodesy/projection/utm.hpp>
#include <navigation/frames/catalogue.hpp>
#include <navigation/frames/galactocentric.hpp>

namespace navigation::frames {

// ---- the syntax -----------------------------------------------------------------------------------

struct Coordinate {
  double value = 0;        // in `unit`
  std::string unit;        // canonical: "", "m", "km", "au", "ly", "pc", "kpc", "Mpc", "deg", "rad", "arcmin", "arcsec", "h"
  char hemisphere = 0;     // 'N', 'S', 'E', 'W', or 0
  std::string word;        // a quoted word (an MGRS reference); then value and unit are unused
  bool is_word() const { return !word.empty(); }
  bool operator==(const Coordinate&) const = default;
};

struct Step {
  std::string chart;  // enu, ned, aer
  std::vector<Coordinate> coords;
  bool operator==(const Step&) const = default;
};

struct Address {
  bool body = true;          // a body's key, or a frame
  std::string name;          // the body's key (lower case) or the frame's name (upper case)
  std::string chart;         // a frame's chart ("" with coordinates means cartesian); for a body ""
  bool has_coords = false;
  std::vector<Coordinate> coords;
  std::vector<Step> steps;
  std::string instant;       // the text after "@", trimmed; "" when none
  bool operator==(const Address&) const = default;
};

class SyntaxError : public std::invalid_argument {
 public:
  SyntaxError(const std::string& what, std::size_t column, std::string expected)
      : std::invalid_argument(what + " at column " + std::to_string(column + 1) + " (expected " + expected + ")"),
        column_(column), expected_(std::move(expected)) {}
  std::size_t column() const { return column_; }          // 0-based, in bytes
  const std::string& expected() const { return expected_; }
 private:
  std::size_t column_;
  std::string expected_;
};

namespace address_detail {

inline const std::vector<std::string_view>& body_keys() {
  static const std::vector<std::string_view> k = {"sun", "mercury", "venus", "earth", "moon", "mars", "jupiter",
                                                  "saturn", "uranus", "neptune", "pluto", "ssb", "galaxy"};
  return k;
}
inline bool is_body(std::string_view s) {
  for (auto k : body_keys())
    if (k == s) return true;
  return false;
}
struct UnitInfo {
  std::string_view spelled, canonical;
  double factor;   // to metres or radians
  bool length;
};
inline const std::vector<UnitInfo>& units() {
  static const std::vector<UnitInfo> u = {
      {"m", "m", 1.0, true}, {"km", "km", 1e3, true}, {"au", "au", AU_M, true}, {"ly", "ly", LIGHT_YEAR_M, true},
      {"pc", "pc", PARSEC_M, true}, {"kpc", "kpc", 1e3 * PARSEC_M, true}, {"Mpc", "Mpc", 1e6 * PARSEC_M, true},
      {"deg", "deg", std::numbers::pi / 180, false}, {"rad", "rad", 1.0, false},
      {"arcmin", "arcmin", std::numbers::pi / 10800, false}, {"arcsec", "arcsec", std::numbers::pi / 648000, false},
      {"h", "h", std::numbers::pi / 12, false}};
  return u;
}
inline const UnitInfo* unit(std::string_view s) {
  for (const auto& u : units())
    if (u.spelled == s) return &u;
  return nullptr;
}

enum class Kind { name, number, word, symbol, end };
struct Token {
  Kind kind = Kind::end;
  std::string text;   // a name, a word, a symbol ("(", ")", ",", ":", "/", "@", "°", "′", "″")
  double number = 0;
  std::size_t column = 0;
};

class Lexer {
 public:
  explicit Lexer(std::string_view s) : s_(s) {}
  std::string_view source() const { return s_; }
  std::size_t position() const { return i_; }
  Token next() {
    while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t')) i_++;
    Token t;
    t.column = i_;
    if (i_ >= s_.size()) return t;
    const unsigned char c = static_cast<unsigned char>(s_[i_]);
    // The symbols of degrees, minutes and seconds of arc, in UTF-8.
    if (s_.substr(i_, 2) == "\xC2\xB0") { i_ += 2; return sym(t, "°"); }
    if (s_.substr(i_, 3) == "\xE2\x80\xB2") { i_ += 3; return sym(t, "′"); }
    if (s_.substr(i_, 3) == "\xE2\x80\xB3") { i_ += 3; return sym(t, "″"); }
    if (std::string_view("(),:/@").find(char(c)) != std::string_view::npos) { i_++; return sym(t, std::string(1, char(c))); }
    if (c == '"') {
      const std::size_t end = s_.find('"', i_ + 1);
      if (end == std::string_view::npos) throw SyntaxError("an unclosed quotation", i_, "a closing \"");
      t.kind = Kind::word;
      t.text = std::string(s_.substr(i_ + 1, end - i_ - 1));
      i_ = end + 1;
      if (t.text.empty()) throw SyntaxError("an empty word", t.column, "a word between the quotation marks");
      return t;
    }
    if (std::isdigit(c) || c == '.' || ((c == '-' || c == '+') && i_ + 1 < s_.size() &&
                                        (std::isdigit(static_cast<unsigned char>(s_[i_ + 1])) || s_[i_ + 1] == '.'))) {
      std::size_t j = i_ + ((c == '-' || c == '+') ? 1 : 0);
      bool digits = false;
      while (j < s_.size() && std::isdigit(static_cast<unsigned char>(s_[j]))) { j++; digits = true; }
      if (j < s_.size() && s_[j] == '.') { j++; while (j < s_.size() && std::isdigit(static_cast<unsigned char>(s_[j]))) { j++; digits = true; } }
      if (!digits) throw SyntaxError("not a number", i_, "digits");
      // An exponent only when digits follow it, so that "2.35E" is 2.35 and the hemisphere E.
      if (j < s_.size() && (s_[j] == 'e' || s_[j] == 'E')) {
        std::size_t k = j + 1;
        if (k < s_.size() && (s_[k] == '+' || s_[k] == '-')) k++;
        if (k < s_.size() && std::isdigit(static_cast<unsigned char>(s_[k]))) {
          while (k < s_.size() && std::isdigit(static_cast<unsigned char>(s_[k]))) k++;
          j = k;
        }
      }
      std::string text(s_.substr(i_, j - i_));
      if (text[0] == '+') text.erase(0, 1);
      double v = 0;
      const auto r = std::from_chars(text.data(), text.data() + text.size(), v);
      if (r.ec != std::errc() || r.ptr != text.data() + text.size()) throw SyntaxError("not a number", i_, "a number");
      t.kind = Kind::number;
      t.number = v;
      t.text = std::string(s_.substr(i_, j - i_));
      i_ = j;
      return t;
    }
    if (std::isalpha(c) || c == '_') {
      std::size_t j = i_;
      while (j < s_.size() && (std::isalnum(static_cast<unsigned char>(s_[j])) || s_[j] == '_')) j++;
      t.kind = Kind::name;
      t.text = std::string(s_.substr(i_, j - i_));
      i_ = j;
      return t;
    }
    throw SyntaxError("an unexpected character", i_, "a name, a number, a word in quotation marks, or one of ( ) , : / @");
  }

 private:
  static Token sym(Token t, std::string s) { t.kind = Kind::symbol; t.text = std::move(s); return t; }
  std::string_view s_;
  std::size_t i_ = 0;
};

class Parser {
 public:
  explicit Parser(std::string_view s) : lex_(s) { advance(); }

  Address address() {
    Address a = point();
    if (is("@")) {
      const std::size_t at = tok_.column + 1;
      std::string rest(lex_.source().substr(at));
      const auto b = rest.find_first_not_of(" \t"), e = rest.find_last_not_of(" \t");
      if (b == std::string::npos) throw SyntaxError("an empty instant", at, "a time after @");
      a.instant = rest.substr(b, e - b + 1);
      return a;
    }
    if (tok_.kind != Kind::end) fail("text after the point", "@ and an instant, or the end");
    return a;
  }

 private:
  Address point() {
    Address a;
    if (tok_.kind != Kind::name) fail("no point", "a body (earth, moon, mars, …) or a frame (ITRS, ICRS, …)");
    std::string name = tok_.text;
    std::string lower = name;
    for (char& c : lower) c = char(std::tolower(static_cast<unsigned char>(c)));
    advance();
    if (is_body(lower)) {
      a.body = true;
      a.name = lower;
      if (is(":")) fail("a body has no charts of its own", "( and coordinates, / and a step, @, or the end; or a frame, as ITRS:geodetic(…)");
    } else {
      a.body = false;
      for (char& c : name) c = char(std::toupper(static_cast<unsigned char>(c)));
      a.name = name;
      if (is(":")) {
        advance();
        if (tok_.kind != Kind::name) fail("no chart", "a chart: cartesian, spherical, equatorial, geodetic, planetocentric, utm, mgrs, …");
        a.chart = tok_.text;
        for (char& c : a.chart) c = char(std::tolower(static_cast<unsigned char>(c)));
        advance();
        if (!is("(")) fail("a chart without coordinates", "( and coordinates");
      }
    }
    if (is("(")) {
      advance();
      a.has_coords = true;
      a.coords = coords();
      expect(")", "a comma and a coordinate, or )");
    }
    while (is("/")) {
      advance();
      if (tok_.kind != Kind::name) fail("no step", "enu, ned or aer");
      Step st;
      st.chart = tok_.text;
      for (char& c : st.chart) c = char(std::tolower(static_cast<unsigned char>(c)));
      if (st.chart != "enu" && st.chart != "ned" && st.chart != "aer") fail("an unknown step", "enu, ned or aer");
      advance();
      expect("(", "( and coordinates");
      st.coords = coords();
      expect(")", "a comma and a coordinate, or )");
      a.steps.push_back(std::move(st));
    }
    return a;
  }

  std::vector<Coordinate> coords() {
    std::vector<Coordinate> out;
    out.push_back(coord());
    while (is(",")) { advance(); out.push_back(coord()); }
    return out;
  }

  Coordinate coord() {
    Coordinate c;
    if (tok_.kind == Kind::word) { c.word = tok_.text; advance(); return c; }
    if (tok_.kind != Kind::number) fail("no coordinate", "a number, or a word in quotation marks");
    const bool negative = !tok_.text.empty() && tok_.text[0] == '-';
    c.value = tok_.number;
    advance();
    if (is("°")) {
      // Degrees, minutes and seconds of arc: the sign of the degrees applies to the whole.
      advance();
      double v = std::fabs(c.value);
      if (tok_.kind == Kind::number) {
        const double m = tok_.number; advance();
        expect("′", "′ after the minutes");
        v += m / 60;
        if (tok_.kind == Kind::number) {
          const double sec = tok_.number; advance();
          expect("″", "″ after the seconds");
          v += sec / 3600;
        }
      }
      c.value = negative ? -v : v;
      c.unit = "deg";
    } else if (tok_.kind == Kind::name && tok_.text == "h") {
      // Hours, minutes and seconds of time, of right ascension.
      advance();
      double v = std::fabs(c.value);
      if (tok_.kind == Kind::number) {
        const double m = tok_.number; advance();
        if (!(tok_.kind == Kind::name && tok_.text == "m")) fail("minutes of time without m", "m after the minutes");
        advance();
        v += m / 60;
        if (tok_.kind == Kind::number) {
          const double sec = tok_.number; advance();
          if (!(tok_.kind == Kind::name && tok_.text == "s")) fail("seconds of time without s", "s after the seconds");
          advance();
          v += sec / 3600;
        }
      }
      c.value = negative ? -v : v;
      c.unit = "h";
    } else if (tok_.kind == Kind::name && unit(tok_.text)) {
      c.unit = std::string(unit(tok_.text)->canonical);
      advance();
    }
    if (tok_.kind == Kind::name && tok_.text.size() == 1 && std::string_view("NSEW").find(tok_.text[0]) != std::string_view::npos) {
      c.hemisphere = tok_.text[0];
      advance();
    }
    return c;
  }

  bool is(std::string_view s) const { return tok_.kind == Kind::symbol && tok_.text == s; }
  void expect(std::string_view s, const std::string& expected) {
    if (!is(s)) fail("a missing " + std::string(s), expected);
    advance();
  }
  void advance() { tok_ = lex_.next(); }
  [[noreturn]] void fail(const std::string& what, const std::string& expected) const {
    const std::string seen = tok_.kind == Kind::end ? "the end" : "\"" + tok_.text + "\"";
    throw SyntaxError(what + ", seeing " + seen, tok_.column, expected);
  }

  Lexer lex_;
  Token tok_;
};

// The shortest decimal that reads back to the same double.
inline std::string shortest(double v) {
  char buf[64];
  const auto r = std::to_chars(buf, buf + sizeof buf, v);
  return std::string(buf, r.ptr);
}
inline std::string fixed(double v, int decimals) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*f", decimals, v);
  return buf;
}
inline std::string significant(double v, int digits) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*g", digits, v);
  return buf;
}

}  // namespace address_detail

inline Address parse_address(std::string_view text) { return address_detail::Parser(text).address(); }

// How numbers are written: shortest (canonical), or rounded for display.
struct Precision {
  int angle_decimals = -1;      // decimals of degrees; −1: shortest
  int length_significant = -1;  // significant digits of lengths; −1: shortest
};

inline std::string print_coordinate(const Coordinate& c, Precision p = {}) {
  using namespace address_detail;
  if (c.is_word()) return "\"" + c.word + "\"";
  const bool angle = c.unit == "deg" || c.unit == "rad" || c.unit == "arcmin" || c.unit == "arcsec" || c.unit == "h";
  std::string s;
  if (angle && p.angle_decimals >= 0 && c.unit == "deg") s = fixed(c.value, p.angle_decimals);
  else if (!angle && !c.unit.empty() && p.length_significant > 0) s = significant(c.value, p.length_significant);
  else s = shortest(c.value);
  if (c.unit == "deg") s += "°";
  else if (!c.unit.empty()) s += " " + c.unit;
  if (c.hemisphere) s += std::string(" ") + c.hemisphere;
  return s;
}

inline std::string print_address(const Address& a, Precision p = {}) {
  std::string s = a.name;
  if (!a.body && !a.chart.empty()) s += ":" + a.chart;
  const auto list = [&](const std::vector<Coordinate>& cs) {
    std::string out = "(";
    for (std::size_t i = 0; i < cs.size(); i++) out += (i ? ", " : "") + print_coordinate(cs[i], p);
    return out + ")";
  };
  if (a.has_coords) s += list(a.coords);
  for (const auto& st : a.steps) s += "/" + st.chart + list(st.coords);
  if (!a.instant.empty()) s += " @ " + a.instant;
  return s;
}

// ---- the meaning ----------------------------------------------------------------------------------

// A body an address may name: its fixed and inertial frames in the catalogue and its figure.
struct BodyFrames {
  std::string key, name;
  std::string fixed;      // body-fixed frame ("" for the barycentre and the Galaxy)
  std::string inertial;   // centred, with ICRS axes (or the galactocentric frame for the Galaxy)
  std::optional<geodesy::figure::Ellipsoid> figure;  // the figure to stand on, if any
  bool surface = false;
  std::string chart;      // its surface chart: geodetic (the Earth) or planetocentric
};

inline std::optional<BodyFrames> body_frames(const Catalogue& c, std::string_view key) {
  const auto& g = c.graph();
  BodyFrames b;
  b.key = std::string(key);
  if (key == "ssb") { b.name = "the Solar System's barycentre"; b.inertial = "ICRS"; return b; }
  if (key == "galaxy") { b.name = "the Galaxy's centre"; b.inertial = "GALACTOCENTRIC"; return b; }
  const auto known = astronomy::ephemerides::body(key);
  if (!known) return std::nullopt;
  b.name = std::string(known->name);
  b.surface = true;
  if (key == "earth") {
    b.fixed = "ITRS"; b.inertial = "GCRS"; b.figure.emplace(geodesy::figure::Ellipsoid::wgs84()); b.chart = "geodetic";
  } else if (key == "moon") {
    b.fixed = c.has_moon_pa() ? "MOON_ME" : "IAU_MOON"; b.inertial = "MOON_ICRS"; b.figure.emplace(Catalogue::moon_sphere()); b.chart = "planetocentric";
  } else {
    std::string up(key);
    for (char& ch : up) ch = char(std::toupper(static_cast<unsigned char>(ch)));
    b.fixed = "IAU_" + up;
    b.inertial = up + "_ICRS";
    b.figure.emplace(geodesy::figure::Ellipsoid::sphere(b.name + ", IAU mean sphere", known->radius));
    b.chart = "planetocentric";
  }
  if (g.find(b.fixed) < 0 || g.find(b.inertial) < 0) return std::nullopt;
  return b;
}

// The body whose fixed frame a frame is, if any.
inline std::optional<BodyFrames> body_of_fixed_frame(const Catalogue& c, const std::string& frame) {
  for (auto key : address_detail::body_keys()) {
    const auto b = body_frames(c, key);
    if (b && b->surface && (b->fixed == frame || (key == "moon" && (frame == "MOON_ME" || frame == "IAU_MOON" || frame == "MOON_PA"))))
      return b;
  }
  return std::nullopt;
}

// A point (or a direction) of space, in a frame of the catalogue.
struct Event {
  std::string frame;
  Vec3 position{0, 0, 0};   // m; a unit vector when `direction`
  bool direction = false;
  std::string body;         // the body whose surface or fixed frame it was given on, or ""
};

class ResolveError : public std::invalid_argument {
 public:
  using std::invalid_argument::invalid_argument;
};

namespace address_detail {

inline double length_of(const Coordinate& c, const char* what) {
  if (c.is_word()) throw ResolveError(std::string(what) + " is a number, not a word");
  if (c.hemisphere) throw ResolveError(std::string(what) + " takes no hemisphere");
  if (c.unit.empty() || c.unit == "m") return c.value;
  const UnitInfo* u = unit(c.unit);
  if (!u || !u->length) throw ResolveError(std::string(what) + " is a length, not in " + c.unit);
  return c.value * u->factor;
}
inline double angle_of(const Coordinate& c, const char* what) {
  if (c.is_word()) throw ResolveError(std::string(what) + " is a number, not a word");
  double a = c.value * std::numbers::pi / 180;  // degrees by default
  if (!c.unit.empty() && c.unit != "deg") {
    const UnitInfo* u = unit(c.unit);
    if (!u || u->length) throw ResolveError(std::string(what) + " is an angle, not in " + c.unit);
    a = c.value * u->factor;
  }
  if (c.hemisphere == 'S' || c.hemisphere == 'W') a = -a;
  return a;
}
// Latitude and longitude, by their hemisphere letters when given (in either order), else in
// the order latitude, longitude (ISO 6709); then an optional third coordinate.
inline std::array<double, 2> latitude_longitude(const std::vector<Coordinate>& cs) {
  if (cs.size() < 2) throw ResolveError("a place needs a latitude and a longitude");
  const auto ns = [](const Coordinate& c) { return c.hemisphere == 'N' || c.hemisphere == 'S'; };
  const auto ew = [](const Coordinate& c) { return c.hemisphere == 'E' || c.hemisphere == 'W'; };
  int ilat = 0, ilon = 1;
  if (ns(cs[1]) || ew(cs[0])) { ilat = 1; ilon = 0; }
  if (ew(cs[ilat]) || ns(cs[ilon])) throw ResolveError("a latitude's hemisphere is N or S, a longitude's E or W");
  const double lat = angle_of(cs[ilat], "a latitude"), lon = angle_of(cs[ilon], "a longitude");
  if (std::fabs(lat) > std::numbers::pi / 2 + 1e-12) throw ResolveError("a latitude lies within 90° of the equator");
  return {lat, lon};
}
inline Vec3 from_spherical(double lon, double lat, double r) {
  return {r * std::cos(lat) * std::cos(lon), r * std::cos(lat) * std::sin(lon), r * std::sin(lat)};
}

}  // namespace address_detail

// The charts an address may give in a frame, the first of them its default for answers.
inline std::vector<std::string> charts_of(const Catalogue& c, const std::string& frame) {
  if (c.graph().find(frame) < 0) return {};
  std::vector<std::string> out;
  if (const auto b = body_of_fixed_frame(c, frame)) {
    out.push_back(b->chart);
    if (b->key == "earth") { out.push_back("utm"); out.push_back("mgrs"); }
    if (b->key == "moon") out.push_back("selenographic");
    out.push_back("spherical");
  } else if (frame == "ECLIPTIC_J2000") {
    out.push_back("ecliptic");
    out.push_back("spherical");
  } else if (frame == "GALACTIC") {
    out.push_back("galactic");
    out.push_back("spherical");
  } else if (frame == "GALACTOCENTRIC") {
    out.push_back("cartesian");
    out.push_back("cylindrical");
    out.push_back("spherical");
    return out;
  } else if (frame == "EARTH_MOON_ROTATING" || frame.find("_ENU") != std::string::npos || frame.find("_NED") != std::string::npos ||
             frame == "TIRS" || frame == "PEF" || frame == "MOON_PA") {
    out.push_back("cartesian");
    out.push_back("spherical");
    return out;
  } else {
    out.push_back("equatorial");
    out.push_back("spherical");
  }
  out.push_back("cartesian");
  out.push_back("cylindrical");
  return out;
}

// The event an address names, in the frame it was given in (steps apply in the body's fixed frame).
inline Event resolve(const Address& a, const Catalogue& c) {
  using namespace address_detail;
  namespace gf = geodesy::figure;
  const auto& g = c.graph();
  Event e;
  if (a.body) {
    const auto b = body_frames(c, a.name);
    if (!b) throw ResolveError("the catalogue does not hold " + a.name);
    if (!a.has_coords) {
      e.frame = b->inertial;  // its centre
      e.body = b->surface ? b->key : "";
    } else {
      if (!b->surface) throw ResolveError(b->name + " has no surface to give coordinates on; use a frame, as " + b->inertial + ":cartesian(…)");
      const auto ll = latitude_longitude(a.coords);
      const double h = a.coords.size() > 2 ? length_of(a.coords[2], "a height") : 0.0;
      if (a.coords.size() > 3) throw ResolveError("a place on a body has at most three coordinates");
      e.frame = b->fixed;
      e.position = b->figure->cartesian(gf::Geodetic{ll[0], ll[1], h});
      e.body = b->key;
    }
  } else {
    if (g.find(a.name) < 0) throw ResolveError("the catalogue has no frame " + a.name);
    e.frame = a.name;
    const auto fixed_body = body_of_fixed_frame(c, a.name);
    if (fixed_body) e.body = fixed_body->key;
    std::string chart = a.chart.empty() ? (a.has_coords ? "cartesian" : "") : a.chart;
    const auto allowed = charts_of(c, a.name);
    if (!chart.empty() && std::find(allowed.begin(), allowed.end(), chart) == allowed.end()) {
      std::string list;
      for (const auto& x : allowed) list += (list.empty() ? "" : ", ") + x;
      throw ResolveError("the chart " + chart + " is not one of " + a.name + "'s: " + list);
    }
    const auto& cs = a.coords;
    if (chart.empty()) {
      e.position = {0, 0, 0};
    } else if (chart == "cartesian") {
      if (cs.size() != 3) throw ResolveError("cartesian coordinates are three lengths");
      e.position = {length_of(cs[0], "x"), length_of(cs[1], "y"), length_of(cs[2], "z")};
    } else if (chart == "spherical" || chart == "equatorial" || chart == "ecliptic" || chart == "galactic") {
      if (cs.size() < 2 || cs.size() > 3) throw ResolveError("spherical coordinates are a longitude, a latitude and, for a point, a distance");
      if (cs[0].hemisphere == 'N' || cs[0].hemisphere == 'S' || cs[1].hemisphere == 'E' || cs[1].hemisphere == 'W')
        throw ResolveError("spherical coordinates are (longitude, latitude[, distance]) in that order");
      const double lon = angle_of(cs[0], "a longitude"), lat = angle_of(cs[1], "a latitude");
      if (std::fabs(lat) > std::numbers::pi / 2 + 1e-12) throw ResolveError("a latitude lies within 90° of the equator");
      if (cs.size() == 3) e.position = from_spherical(lon, lat, length_of(cs[2], "a distance"));
      else { e.position = from_spherical(lon, lat, 1.0); e.direction = true; }
    } else if (chart == "cylindrical") {
      if (cs.size() != 3) throw ResolveError("cylindrical coordinates are ρ, φ and z");
      const double rho = length_of(cs[0], "ρ"), phi = angle_of(cs[1], "φ"), z = length_of(cs[2], "z");
      e.position = {rho * std::cos(phi), rho * std::sin(phi), z};
    } else if (chart == "geodetic" || chart == "planetocentric" || chart == "selenographic") {
      if (!fixed_body) throw ResolveError(chart + " coordinates belong to a body's fixed frame");
      const auto ll = latitude_longitude(cs);
      const double h = cs.size() > 2 ? length_of(cs[2], "a height") : 0.0;
      e.position = fixed_body->figure->cartesian(gf::Geodetic{ll[0], ll[1], h});
    } else if (chart == "utm") {
      if (cs.size() < 3) throw ResolveError("UTM coordinates are a zone with its hemisphere, an easting and a northing");
      if (cs[0].is_word() || (cs[0].hemisphere != 'N' && cs[0].hemisphere != 'S') || !cs[0].unit.empty() || cs[0].value != std::floor(cs[0].value))
        throw ResolveError("a UTM zone is written as its number and hemisphere, as 31 N; UPS is zone 0");
      const int zone = int(cs[0].value);
      const geodesy::projection::UniversalGrids grids(gf::Ellipsoid::wgs84());
      const gf::Geodetic where = grids.reverse(zone, cs[0].hemisphere == 'N', length_of(cs[1], "an easting"), length_of(cs[2], "a northing"));
      const double h = cs.size() > 3 ? length_of(cs[3], "a height") : 0.0;
      e.position = gf::Ellipsoid::wgs84().cartesian({where.latitude, where.longitude, h});
    } else if (chart == "mgrs") {
      if (cs.empty() || !cs[0].is_word()) throw ResolveError("an MGRS reference is written in quotation marks, as \"31UDQ4825211954\"");
      const geodesy::projection::Mgrs mgrs{geodesy::projection::UniversalGrids(gf::Ellipsoid::wgs84())};
      const auto parsed = mgrs.reverse(cs[0].word);
      const double h = cs.size() > 1 ? length_of(cs[1], "a height") : 0.0;
      e.position = gf::Ellipsoid::wgs84().cartesian({parsed.geodetic.latitude, parsed.geodetic.longitude, h});
    }
  }
  // Steps from the local horizon at the point so far, on its body's figure.
  for (const auto& st : a.steps) {
    if (e.direction) throw ResolveError("a step starts from a point, not a direction");
    if (e.body.empty()) throw ResolveError("a step needs a point on a body (its horizon); this point is in " + e.frame);
    const auto b = body_frames(c, e.body);
    if (!b || !b->surface) throw ResolveError("a step needs a body with a figure");
    if (st.coords.size() != 3) throw ResolveError("a step has three coordinates");
    // The point in the body's fixed frame, its foot on the figure, and the horizon there.
    if (e.frame != b->fixed) throw ResolveError("a step starts from a point given on the body's surface or in its fixed frame " + b->fixed);
    const gf::Geodetic foot = b->figure->geodetic(e.position);
    const auto place = gf::Place::at(*b->figure, foot);
    Vec3 enu;
    if (st.chart == "enu") enu = {length_of(st.coords[0], "east"), length_of(st.coords[1], "north"), length_of(st.coords[2], "up")};
    else if (st.chart == "ned") enu = {length_of(st.coords[1], "east"), length_of(st.coords[0], "north"), -length_of(st.coords[2], "down")};
    else {
      const double az = angle_of(st.coords[0], "an azimuth"), el = angle_of(st.coords[1], "an elevation"), range = length_of(st.coords[2], "a range");
      enu = gf::enu_from_look({az, el, range});
    }
    e.position = place.ecef_of(enu);
  }
  return e;
}

// An event carried to another frame at an instant.
inline Event carry(const Event& e, const std::string& frame, const Catalogue& c, const Instant& t) {
  const auto& g = c.graph();
  Event out = e;
  out.frame = frame;
  const auto x = g.transform(e.frame, frame, t);
  out.position = e.direction ? x.direction(e.position) : x.position(e.position);
  const auto fixed_body = body_of_fixed_frame(c, frame);
  out.body = fixed_body ? fixed_body->key : (e.body);
  if (!fixed_body && frame != e.frame) out.body = "";
  return out;
}

// An event in a frame and chart: the named values and the address that gives it.
struct Value {
  std::string name;
  double value = 0;
  std::string unit;  // "m", "deg", …
};
struct Expression {
  std::string frame, chart;
  std::vector<Value> values;
  Address address;
};

namespace address_detail {
// A length in a unit fitting its size, for answers: m below 1000 km, km below 0.01 au, au below
// 0.1 pc, pc below 1 kpc, kpc below 1 Mpc, else Mpc.
inline const char* length_unit(double metres) {
  const double m = std::fabs(metres);
  return m < 1e6 ? "m" : m < 0.01 * AU_M ? "km" : m < 0.1 * PARSEC_M ? "au" : m < 1e3 * PARSEC_M ? "pc" : m < 1e6 * PARSEC_M ? "kpc" : "Mpc";
}
// A length in a unit: its own size's (length_unit), or a unit shared by the coordinates of one point.
inline Coordinate length_coordinate(double metres, const char* u = nullptr) {
  Coordinate c;
  if (!u) u = length_unit(metres);
  c.unit = u;
  c.value = metres / unit(u)->factor;
  return c;
}
inline Coordinate angle_coordinate(double radians, char positive = 0, char negative = 0) {
  Coordinate c;
  c.unit = "deg";
  double d = radians * 180 / std::numbers::pi;
  if (positive) { c.hemisphere = d < 0 ? negative : positive; d = std::fabs(d); }
  c.value = d;
  return c;
}
}  // namespace address_detail

inline Expression express(const Event& event, const std::string& frame, std::string chart, const Catalogue& c, const Instant& t) {
  using namespace address_detail;
  namespace gf = geodesy::figure;
  const auto allowed = charts_of(c, frame);
  if (allowed.empty()) throw ResolveError("the catalogue has no frame " + frame);
  if (chart.empty()) chart = allowed.front();
  if (std::find(allowed.begin(), allowed.end(), chart) == allowed.end()) throw ResolveError("the chart " + chart + " is not one of " + frame + "'s");
  const Event e = carry(event, frame, c, t);
  Expression x;
  x.frame = frame;
  x.chart = chart;
  x.address.body = false;
  x.address.name = frame;
  x.address.chart = chart;
  x.address.has_coords = true;
  const Vec3& p = e.position;
  const auto push = [&](const char* name, double value, const char* unit_name, const Coordinate& coord) {
    x.values.push_back({name, value, unit_name});
    x.address.coords.push_back(coord);
  };
  if (chart == "cartesian") {
    if (e.direction) throw ResolveError("a direction has no Cartesian coordinates; ask for a spherical chart");
    const char* u = length_unit(physis::geometry::norm(p));  // one unit for the point, by its distance
    push("x_m", p[0], "m", length_coordinate(p[0], u));
    push("y_m", p[1], "m", length_coordinate(p[1], u));
    push("z_m", p[2], "m", length_coordinate(p[2], u));
  } else if (chart == "spherical" || chart == "equatorial" || chart == "ecliptic" || chart == "galactic") {
    const double r = physis::geometry::norm(p);
    double lon = std::atan2(p[1], p[0]);
    if (lon < 0) lon += 2 * std::numbers::pi;
    const double lat = r > 0 ? std::asin(std::clamp(p[2] / r, -1.0, 1.0)) : 0.0;
    const char* lon_name = chart == "equatorial" ? "right_ascension_deg" : "longitude_deg";
    const char* lat_name = chart == "equatorial" ? "declination_deg" : "latitude_deg";
    push(lon_name, lon * 180 / std::numbers::pi, "deg", angle_coordinate(lon));
    Coordinate la = angle_coordinate(lat);
    push(lat_name, lat * 180 / std::numbers::pi, "deg", la);
    if (!e.direction) push("distance_m", r, "m", length_coordinate(r));
  } else if (chart == "cylindrical") {
    if (e.direction) throw ResolveError("a direction has no cylindrical coordinates");
    const double rho = std::hypot(p[0], p[1]), phi = std::atan2(p[1], p[0]);
    const char* u = length_unit(physis::geometry::norm(p));
    push("rho_m", rho, "m", length_coordinate(rho, u));
    push("phi_deg", phi * 180 / std::numbers::pi, "deg", angle_coordinate(phi));
    push("z_m", p[2], "m", length_coordinate(p[2], u));
  } else if (chart == "geodetic" || chart == "planetocentric" || chart == "selenographic") {
    if (e.direction) throw ResolveError("a direction is not a place");
    const auto b = body_of_fixed_frame(c, frame);
    const gf::Geodetic gd = b->figure->geodetic(p);
    push("latitude_deg", gd.latitude * 180 / std::numbers::pi, "deg", angle_coordinate(gd.latitude, 'N', 'S'));
    push("longitude_deg", gd.longitude * 180 / std::numbers::pi, "deg", angle_coordinate(gd.longitude, 'E', 'W'));
    Coordinate h;
    h.value = gd.height;
    h.unit = "m";
    push("height_m", gd.height, "m", h);
  } else if (chart == "utm" || chart == "mgrs") {
    if (e.direction) throw ResolveError("a direction is not a place");
    const gf::Geodetic gd = gf::Ellipsoid::wgs84().geodetic(p);
    const geodesy::projection::UniversalGrids grids(gf::Ellipsoid::wgs84());
    const auto u = grids.forward(gd.latitude, gd.longitude);
    Coordinate h;
    h.value = gd.height;
    h.unit = "m";
    if (chart == "utm") {
      Coordinate zone;
      zone.value = u.zone;
      zone.hemisphere = u.north ? 'N' : 'S';
      x.values.push_back({"zone", double(u.zone), ""});
      x.address.coords.push_back(zone);
      Coordinate east, north;
      east.unit = north.unit = "m";
      east.value = u.easting;
      north.value = u.northing;
      push("easting_m", u.easting, "m", east);
      push("northing_m", u.northing, "m", north);
      push("height_m", gd.height, "m", h);
    } else {
      const geodesy::projection::Mgrs mgrs{geodesy::projection::UniversalGrids(gf::Ellipsoid::wgs84())};
      Coordinate ref;
      ref.word = mgrs.format(u, 5);
      x.address.coords.push_back(ref);
      x.values.push_back({"mgrs", 0, ref.word});
      push("height_m", gd.height, "m", h);
    }
  }
  return x;
}

// A point seen from an observer standing on a body: ENU, NED or azimuth, elevation and range from
// the observer's horizon on the body's figure.
inline Expression express_from(const Event& point, const Event& observer, const std::string& chart, const Catalogue& c, const Instant& t) {
  using namespace address_detail;
  namespace gf = geodesy::figure;
  if (observer.body.empty()) throw ResolveError("an observer stands on a body");
  const auto b = body_frames(c, observer.body);
  if (!b || !b->surface) throw ResolveError("an observer stands on a body with a figure");
  const Event o = carry(observer, b->fixed, c, t);
  const Event p = carry(point, b->fixed, c, t);
  const auto place = gf::Place::at(*b->figure, b->figure->geodetic(o.position));
  Expression x;
  x.frame = b->fixed;
  x.chart = chart;
  // The address of the point from the observer: the observer's place and a step.
  x.address = express(observer, b->fixed, b->chart, c, t).address;
  Step st;
  st.chart = chart;
  if (p.direction) {
    const Vec3 enu = place.to_enu * p.position;
    const auto look = gf::look_angles(enu);
    if (chart != "aer") throw ResolveError("a direction is seen as an azimuth and an elevation (aer)");
    x.values = {{"azimuth_deg", look.azimuth * 180 / std::numbers::pi, "deg"}, {"elevation_deg", look.elevation * 180 / std::numbers::pi, "deg"}};
    st.coords = {angle_coordinate(look.azimuth), angle_coordinate(look.elevation)};
  } else {
    const Vec3 enu = place.enu_of(p.position);
    if (chart == "enu") {
      x.values = {{"east_m", enu[0], "m"}, {"north_m", enu[1], "m"}, {"up_m", enu[2], "m"}};
      const char* u = length_unit(physis::geometry::norm(enu));
      st.coords = {length_coordinate(enu[0], u), length_coordinate(enu[1], u), length_coordinate(enu[2], u)};
    } else if (chart == "ned") {
      x.values = {{"north_m", enu[1], "m"}, {"east_m", enu[0], "m"}, {"down_m", -enu[2], "m"}};
      const char* u = length_unit(physis::geometry::norm(enu));
      st.coords = {length_coordinate(enu[1], u), length_coordinate(enu[0], u), length_coordinate(-enu[2], u)};
    } else if (chart == "aer") {
      const auto look = gf::look_angles(enu);
      x.values = {{"azimuth_deg", look.azimuth * 180 / std::numbers::pi, "deg"}, {"elevation_deg", look.elevation * 180 / std::numbers::pi, "deg"},
                  {"range_m", look.range, "m"}};
      st.coords = {angle_coordinate(look.azimuth), angle_coordinate(look.elevation), length_coordinate(look.range)};
    } else {
      throw ResolveError("a point is seen from an observer in enu, ned or aer");
    }
  }
  x.address.steps.push_back(st);
  return x;
}

// The grammar, as the page shows it.
inline std::string address_grammar() {
  return "address  ::= point [ \"@\" instant ]\n"
         "point    ::= base { \"/\" step }\n"
         "base     ::= BODY [ \"(\" coords \")\" ] | FRAME [ \":\" CHART ] [ \"(\" coords \")\" ]\n"
         "step     ::= ( \"enu\" | \"ned\" | \"aer\" ) \"(\" coords \")\"\n"
         "coords   ::= coord { \",\" coord }\n"
         "coord    ::= STRING | NUMBER [ UNIT | DMS | HMS ] [ HEMISPHERE ]\n"
         "DMS      ::= \"°\" [ NUMBER \"′\" [ NUMBER \"″\" ] ]\n"
         "HMS      ::= \"h\" [ NUMBER \"m\" [ NUMBER \"s\" ] ]\n"
         "instant  ::= ISO 8601 with Z or a scale | JD or MJD with a scale | now\n";
}

}  // namespace navigation::frames
