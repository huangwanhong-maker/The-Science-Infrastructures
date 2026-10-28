// The two-line element set (TLE): the mean elements of an Earth satellite as the United States
// Space Force's catalogue publishes them, for the SGP4/SDP4 model (sgp4.hpp).
//
// A set is two lines of 69 columns, fixed by the punch cards they were designed for (Hoots and
// Roehrich 1980, Spacetrack Report #3, §12; Vallado, Crawford, Hujsak and Kelso 2006, "Revisiting
// Spacetrack Report #3", AIAA 2006-6753, Appendix A), often preceded by a line naming the object
// (a "three-line" set):
//
//   1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
//   2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
//
// line 1: the catalogue number (Alpha-5 beyond 99999: a letter for the leading two digits, I and O
// skipped), its classification, the international designator, the epoch (two-digit year, and the
// day of the year with its fraction, 1.0 being January 1 at 0h UTC), ṅ/2 (rev/day²), n̈/6 (rev/day³,
// mantissa and exponent with an implied decimal point), B* (per Earth radius, likewise), the
// ephemeris type and the element set's number; line 2: the inclination, the right ascension of
// the ascending node, the eccentricity (with an implied leading decimal point), the argument of
// perigee and the mean anomaly (degrees), the mean motion (rev/day) and the revolution number. The
// last column of each line is a checksum: the sum of the line's digits, each minus sign counting
// one, modulo 10. Years 57–99 are 1957–1999 and 00–56 are 2000–2056 (Vallado's convention, "temp
// fix for years from 1957-2056" in sgp4io.cpp).
//
// The elements are mean elements of SGP4's theory, in its TEME frame; they are not osculating
// elements and are meaningful only through SGP4. They are given here in SI units and radians, the
// epoch as a two-part Julian date in UTC (the TLE's own scale): d1 the day's 0h, d2 its fraction
// as the TLE writes it.
//
// The reading follows sgp4io.cpp (Vallado 2007, twoline2rv) as python-sgp4 2.27 transcribes it
// (sgp4/io.py; © Brandon Rhodes, MIT licence): the same columns, checks and conversions. A line
// may carry text past column 69 (Vallado's verification file puts the times of its runs there);
// it is ignored.

#pragma once

#include <cmath>
#include <cstdlib>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <metrology/time/julian.hpp>

namespace navigation::orbits {

namespace detail {
inline std::string trim(std::string_view s) {
  std::size_t a = 0, b = s.size();
  while (a < b && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
  while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
  return std::string(s.substr(a, b - a));
}
inline double number(std::string_view s, const char* what) {
  const std::string t = trim(s);
  if (t.empty()) return 0;
  char* end = nullptr;
  const double v = std::strtod(t.c_str(), &end);
  if (end == t.c_str() || *end != '\0') throw std::invalid_argument(std::string("TLE: ") + what + " is not a number: \"" + t + "\"");
  return v;
}
inline long integer(std::string_view s, const char* what) {
  const std::string t = trim(s);
  if (t.empty()) return 0;
  char* end = nullptr;
  const long v = std::strtol(t.c_str(), &end, 10);
  if (end == t.c_str() || *end != '\0') throw std::invalid_argument(std::string("TLE: ") + what + " is not an integer: \"" + t + "\"");
  return v;
}
// A field with an implied decimal point and an exponent, columns "SMMMMMEE": ±0.MMMMM × 10^EE.
inline double exponential(std::string_view line, std::size_t sign, const char* what) {
  std::string m = std::string(1, line[sign] == '-' ? '-' : '+') + "0." + std::string(line.substr(sign + 1, 5));
  for (char& c : m)
    if (c == ' ') c = '0';
  return number(m, what) * std::pow(10.0, static_cast<double>(integer(line.substr(sign + 6, 2), what)));
}
inline std::string_view without_end(std::string_view s) {
  while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ')) s.remove_suffix(1);
  return s;
}
}  // namespace detail

// The checksum of a line: its first 68 columns' digits, each '-' counting one, modulo 10.
inline int tle_checksum(std::string_view line) {
  int sum = 0;
  for (std::size_t i = 0; i < std::min<std::size_t>(68, line.size()); i++) {
    const char c = line[i];
    if (c >= '0' && c <= '9') sum += c - '0';
    else if (c == '-') sum += 1;
  }
  return sum % 10;
}

// The Alpha-5 catalogue number: five digits, or a letter (A = 10 … Z = 33, I and O skipped) and
// four digits, up to Z9999 = 339999.
inline long alpha5_number(std::string_view s) {
  const std::string t = detail::trim(s);
  if (t.empty()) return 0;
  const char c = t[0];
  if (c >= '0' && c <= '9') return detail::integer(t, "the catalogue number");
  if (c < 'A' || c > 'Z' || c == 'I' || c == 'O') throw std::invalid_argument("TLE: not an Alpha-5 catalogue number: \"" + t + "\"");
  const long lead = (c - 'A' + 10) - (c > 'I' ? 1 : 0) - (c > 'O' ? 1 : 0);
  return lead * 10000 + detail::integer(t.substr(1), "the catalogue number");
}

struct Tle {
  std::string name;              // the line naming the object in a three-line set, else empty
  std::string line1, line2;      // as given, without their line endings
  std::string catalogue;         // columns 3–7 as written (Alpha-5 kept)
  long number = 0;               // the catalogue number, Alpha-5 decoded
  char classification = 'U';     // U unclassified, C classified, S secret
  std::string designator;        // international designator: launch year, number and piece
  int epoch_year = 0;            // four digits
  double epoch_day = 0;          // day of the year and its fraction, 1.0 = January 1, 0h UTC
  metrology::time::JulianDate epoch;  // UTC: d1 the day's 0h, d2 its fraction
  // The mean elements (SGP4's, in TEME), SI and radians.
  double mean_motion_dot_half = 0;    // ṅ/2, rad/s² (the TLE's rev/day²)
  double mean_motion_ddot_sixth = 0;  // n̈/6, rad/s³ (the TLE's rev/day³)
  double bstar = 0;              // B*, per Earth radius (SGP4's drag term, in its own unit)
  int ephemeris_type = 0;
  long element_number = 0;
  double inclination = 0;        // rad
  double raan = 0;               // right ascension of the ascending node, rad
  double eccentricity = 0;
  double argument_of_perigee = 0;  // rad
  double mean_anomaly = 0;       // rad
  double mean_motion = 0;        // Kozai mean motion, rad/s
  long revolution = 0;           // revolutions at epoch
  // Each line's checksum as written (−1 when the column is not a digit) and as computed.
  int checksum1 = -1, checksum2 = -1, computed1 = 0, computed2 = 0;
  bool checksums_ok() const { return checksum1 == computed1 && checksum2 == computed2; }

  // The epoch in SGP4's own count: days from 1949 December 31, 0h UTC (JD 2433281.5).
  double days_since_1950() const { return (epoch.d1 - 2433281.5) + epoch.d2; }
};

// A set from its two lines (and its name, if a three-line set gave one). Malformed lines are
// refused with the column that is wrong; a wrong checksum is reported in the set, not refused.
inline Tle parse_tle(std::string_view first, std::string_view second, std::string_view name = {}) {
  using detail::exponential;
  using detail::integer;
  using detail::number;
  constexpr double deg = std::numbers::pi / 180.0;
  constexpr double rev_per_day = 2.0 * std::numbers::pi / 86400.0;  // rad/s for 1 rev/day
  const std::string_view a = detail::without_end(first), b = detail::without_end(second);
  auto refuse = [](int n, std::string_view line, const char* what) {
    throw std::invalid_argument("TLE line " + std::to_string(n) + " " + what + ": \"" + std::string(line) + "\"");
  };
  if (a.size() < 64 || a.substr(0, 2) != "1 ") refuse(1, a, "does not begin \"1 \" or is shorter than 64 columns");
  for (std::size_t i : {8u, 32u, 43u, 52u, 61u, 63u})
    if (a[i] != ' ') refuse(1, a, ("needs a space in column " + std::to_string(i + 1)).c_str());
  if (a[23] != '.' || a[34] != '.') refuse(1, a, "needs its decimal points in columns 24 and 35");
  if (b.size() < 68 || b.substr(0, 2) != "2 ") refuse(2, b, "does not begin \"2 \" or is shorter than 68 columns");
  for (std::size_t i : {7u, 16u, 25u, 33u, 42u, 51u})
    if (b[i] != ' ') refuse(2, b, ("needs a space in column " + std::to_string(i + 1)).c_str());
  for (std::size_t i : {11u, 20u, 37u, 46u})
    if (b[i] != '.') refuse(2, b, ("needs a decimal point in column " + std::to_string(i + 1)).c_str());
  if (a.substr(2, 5) != b.substr(2, 5)) throw std::invalid_argument("TLE: the catalogue numbers of lines 1 and 2 differ");

  Tle t;
  t.name = detail::trim(name);
  if (!t.name.empty() && t.name.rfind("0 ", 0) == 0) t.name = detail::trim(t.name.substr(2));  // the "0 NAME" form
  t.line1 = std::string(a);
  t.line2 = std::string(b);
  t.catalogue = std::string(a.substr(2, 5));
  t.number = alpha5_number(t.catalogue);
  t.classification = a[7] == ' ' ? 'U' : a[7];
  t.designator = detail::trim(a.substr(9, 8));
  const long yy = integer(a.substr(18, 2), "the epoch's year");
  t.epoch_year = static_cast<int>(yy < 57 ? 2000 + yy : 1900 + yy);
  // The day of the year: its whole days and its fraction read apart, so that the fraction keeps
  // the digits the TLE gives.
  const std::string day_text = detail::trim(a.substr(20, 12));
  const auto point = day_text.find('.');
  const long whole_day = integer(day_text.substr(0, point), "the epoch's day");
  const double fraction = point == std::string::npos ? 0.0 : number("0" + day_text.substr(point), "the epoch's fraction of a day");
  t.epoch_day = static_cast<double>(whole_day) + fraction;
  const metrology::time::JulianDate jan1 = metrology::time::calendar_to_jd(t.epoch_year, 1, 1);
  t.epoch = {jan1.d1 + jan1.d2 + static_cast<double>(whole_day - 1), fraction};
  t.mean_motion_dot_half = number(a.substr(33, 10), "ṅ/2") * rev_per_day / 86400.0;
  t.mean_motion_ddot_sixth = exponential(a, 44, "n̈/6") * rev_per_day / (86400.0 * 86400.0);
  t.bstar = exponential(a, 53, "B*");
  t.ephemeris_type = a[62] == ' ' ? 0 : static_cast<int>(integer(a.substr(62, 1), "the ephemeris type"));
  t.element_number = integer(a.substr(64, 4), "the element set's number");

  t.inclination = number(b.substr(8, 8), "the inclination") * deg;
  t.raan = number(b.substr(17, 8), "the node") * deg;
  std::string e = "0." + std::string(b.substr(26, 7));
  for (char& c : e)
    if (c == ' ') c = '0';
  t.eccentricity = number(e, "the eccentricity");
  t.argument_of_perigee = number(b.substr(34, 8), "the argument of perigee") * deg;
  t.mean_anomaly = number(b.substr(43, 8), "the mean anomaly") * deg;
  t.mean_motion = number(b.substr(52, 11), "the mean motion") * rev_per_day;
  t.revolution = integer(b.substr(63, 5), "the revolution number");

  t.computed1 = tle_checksum(a);
  t.computed2 = tle_checksum(b);
  t.checksum1 = a.size() > 68 && a[68] >= '0' && a[68] <= '9' ? a[68] - '0' : -1;
  t.checksum2 = b.size() > 68 && b[68] >= '0' && b[68] <= '9' ? b[68] - '0' : -1;
  return t;
}

// Every set in a text of two- or three-line sets (CelesTrak's files, Space-Track's): a line
// beginning "1 " followed by one beginning "2 " is a set, and a line just before it that is
// neither is its name. Lines beginning '#' are comments.
inline std::vector<Tle> read_tles(std::string_view text) {
  std::vector<std::string_view> lines;
  std::size_t start = 0;
  while (start <= text.size()) {
    std::size_t end = text.find('\n', start);
    if (end == std::string_view::npos) end = text.size();
    lines.push_back(detail::without_end(text.substr(start, end - start)));
    start = end + 1;
  }
  std::vector<Tle> out;
  for (std::size_t i = 0; i + 1 < lines.size(); i++) {
    if (lines[i].substr(0, 2) == "1 " && lines[i + 1].substr(0, 2) == "2 ") {
      std::string_view name;
      if (i > 0 && !lines[i - 1].empty() && lines[i - 1][0] != '#' && lines[i - 1].substr(0, 2) != "1 " && lines[i - 1].substr(0, 2) != "2 ")
        name = lines[i - 1];
      out.push_back(parse_tle(lines[i], lines[i + 1], name));
      i++;
    }
  }
  return out;
}

}  // namespace navigation::orbits
