// Chronology: the dates of a record and the calendars they are written in (codebook §4).
//
// A date is EDTF (ISO 8601-2:2019, level 1) in a subset: a year, a month or a day, with an
// optional qualifier on the whole (? uncertain, ~ approximate, % both), unspecified digits X
// in the year, or an interval of two such. It is written as the source writes it, in one of
// three calendars: the Gregorian (proleptic before 1582), the Julian, and Old Style as England
// wrote it until 1752 (Julian, the legal year beginning on 25 March). Every date becomes the
// span of Julian day numbers it may denote, and occasions are ordered by them.
//
// The day numbers are the chronological Julian day numbers (JDN): the integer count of days
// whose day 0 is 1 January 4713 BC of the Julian calendar (year −4712 in astronomical
// numbering), so that JDN 2451545 is 1 January 2000 of the Gregorian calendar. The
// conversions are those of Fliegel and Van Flandern (1968) and Richards, in the Explanatory
// Supplement to the Astronomical Almanac (3rd ed., 2013), §15.11.3; years are astronomical
// (year 0 is 1 BC). Valid for JDN ≥ 0.

#pragma once

#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace history {

enum class Calendar { gregorian, julian, julian_march };

inline std::optional<Calendar> calendar_of(std::string_view name) {
  if (name.empty() || name == "gregorian") return Calendar::gregorian;
  if (name == "julian") return Calendar::julian;
  if (name == "julian-march") return Calendar::julian_march;
  return std::nullopt;
}

inline std::string calendar_name(Calendar c) {
  return c == Calendar::gregorian ? "gregorian" : c == Calendar::julian ? "julian" : "julian-march";
}

// Floor division, so that the formulae hold for every year from −4712 on.
inline long floor_div(long a, long b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }

// The JDN of a day of the Gregorian calendar (Explanatory Supplement 2013, §15.11.3).
inline long jdn_gregorian(long year, int month, int day) {
  const long a = floor_div(14 - month, 12);
  const long y = year + 4800 - a, m = month + 12 * a - 3;
  return day + floor_div(153 * m + 2, 5) + 365 * y + floor_div(y, 4) - floor_div(y, 100) + floor_div(y, 400) - 32045;
}

// The JDN of a day of the Julian calendar, the year beginning on 1 January.
inline long jdn_julian(long year, int month, int day) {
  const long a = floor_div(14 - month, 12);
  const long y = year + 4800 - a, m = month + 12 * a - 3;
  return day + floor_div(153 * m + 2, 5) + 365 * y + floor_div(y, 4) - 32083;
}

struct Ymd {
  long year = 0;
  int month = 1, day = 1;
  friend bool operator==(const Ymd&, const Ymd&) = default;
};

// The inverses (Richards's algorithm; the Julian one is the Gregorian without the centuries).
inline Ymd gregorian_of(long jdn) {
  const long a = jdn + 32044, b = floor_div(4 * a + 3, 146097), c = a - floor_div(146097 * b, 4);
  const long d = floor_div(4 * c + 3, 1461), e = c - floor_div(1461 * d, 4), m = floor_div(5 * e + 2, 153);
  return {100 * b + d - 4800 + floor_div(m, 10), static_cast<int>(m + 3 - 12 * floor_div(m, 10)), static_cast<int>(e - floor_div(153 * m + 2, 5) + 1)};
}

inline Ymd julian_of(long jdn) {
  const long c = jdn + 32082, d = floor_div(4 * c + 3, 1461), e = c - floor_div(1461 * d, 4), m = floor_div(5 * e + 2, 153);
  return {d - 4800 + floor_div(m, 10), static_cast<int>(m + 3 - 12 * floor_div(m, 10)), static_cast<int>(e - floor_div(153 * m + 2, 5) + 1)};
}

// The JDN of a day as written in a calendar. In Old Style a day from 1 January to 24 March
// belongs to the legal year before: "1726-03-20" written Old Style is 20 March 1727, Julian.
inline long jdn_of(Calendar c, long year, int month, int day) {
  if (c == Calendar::gregorian) return jdn_gregorian(year, month, day);
  if (c == Calendar::julian_march && (month < 3 || (month == 3 && day < 25))) year += 1;
  return jdn_julian(year, month, day);
}

inline int days_in_month(Calendar c, long year, int month) {
  const bool gregorian = c == Calendar::gregorian;
  const long next_y = month == 12 ? year + 1 : year;
  const int next_m = month == 12 ? 1 : month + 1;
  return static_cast<int>((gregorian ? jdn_gregorian(next_y, next_m, 1) : jdn_julian(next_y, next_m, 1)) -
                          (gregorian ? jdn_gregorian(year, month, 1) : jdn_julian(year, month, 1)));
}

// A day as ISO 8601 text in the Gregorian calendar: "1808-10-02".
inline std::string iso_gregorian(long jdn) {
  const Ymd g = gregorian_of(jdn);
  char buf[32];
  std::snprintf(buf, sizeof buf, "%s%04ld-%02d-%02d", g.year < 0 ? "-" : "", g.year < 0 ? -g.year : g.year, g.month, g.day);
  return buf;
}

// The decimal Gregorian year of the beginning of a day, for drawing: 1808.75 for 2 October 1808.
inline double decimal_year(long jdn) {
  const Ymd g = gregorian_of(jdn);
  const long start = jdn_gregorian(g.year, 1, 1), end = jdn_gregorian(g.year + 1, 1, 1);
  return static_cast<double>(g.year) + static_cast<double>(jdn - start) / static_cast<double>(end - start);
}

// ---- The historians' convention ------------------------------------------------------------------
//
// A date that names no calendar, as a reader gives one (the cut of a moment, the end of a year, a
// year on the axis of time), is read as historians and astronomers read it (Meeus 1998, ch. 7;
// Explanatory Supplement 2013, §15.11.1): in the Julian calendar before the Gregorian reform,
// which began on 15 October 1582 (JDN 2299161), and in the Gregorian from then on. A record dated
// in the Julian calendar before the reform, as the canons of antiquity and of the Renaissance are,
// then has its years end where the reader's do. Read in the proleptic Gregorian, the Julian year
// 337 BC begins five days before the reader's 338 BC ends, and every occasion dated by its year
// alone falls in the year before.

inline constexpr long GREGORIAN_REFORM = 2299161;  // 15 October 1582, Gregorian

// The JDN of a day in the historians' convention. The ten days the reform removed, 5 to 14
// October 1582, are read in the Julian calendar, where they are 15 to 24 October.
inline long jdn_historical(long year, int month, int day) {
  const long g = jdn_gregorian(year, month, day);
  return g >= GREGORIAN_REFORM ? g : jdn_julian(year, month, day);
}

inline Ymd historical_of(long jdn) { return jdn >= GREGORIAN_REFORM ? gregorian_of(jdn) : julian_of(jdn); }

// The decimal year of the beginning of a day in the historians' convention, for drawing: years
// before the reform are Julian years, so that −337.5 is the middle of 338 BC as Athens dated it.
inline double historical_year(long jdn) {
  const Ymd h = historical_of(jdn);
  const long start = jdn_historical(h.year, 1, 1), end = jdn_historical(h.year + 1, 1, 1);
  return static_cast<double>(h.year) + static_cast<double>(jdn - start) / static_cast<double>(end - start);
}

struct Date {
  std::string text;                 // as written
  Calendar calendar = Calendar::gregorian;
  long first = 0, last = 0;         // the first and last days it may denote (JDN), inclusive
  int precision = 0;                // of its first part: 1 a year, 2 a month, 3 a day
  bool interval = false, uncertain = false, approximate = false;
};

namespace detail {

struct Part {
  long first = 0, last = 0;
  int precision = 0;
  bool uncertain = false, approximate = false;
};

inline long digits(std::string_view s, std::string_view what) {
  if (s.empty()) throw std::invalid_argument(std::string(what) + " is empty");
  long v = 0;
  for (char ch : s) {
    if (ch < '0' || ch > '9') throw std::invalid_argument(std::string(what) + " \"" + std::string(s) + "\" is not a number");
    v = 10 * v + (ch - '0');
  }
  return v;
}

// One part of a date: [-]YYYY[-MM[-DD]] with a qualifier at its end.
inline Part part(std::string_view s, Calendar c) {
  Part p;
  while (!s.empty() && (s.back() == '?' || s.back() == '~' || s.back() == '%')) {
    if (s.back() != '~') p.uncertain = true;
    if (s.back() != '?') p.approximate = true;
    s.remove_suffix(1);
  }
  bool negative = false;
  if (!s.empty() && s.front() == '-') { negative = true; s.remove_prefix(1); }
  const auto dash1 = s.find('-');
  const std::string_view ys = s.substr(0, dash1);
  if (ys.size() < 4) throw std::invalid_argument("a year has at least four digits: \"" + std::string(ys) + "\"");
  // Unspecified digits in the year (EDTF level 1): 190X is 1900 to 1909.
  std::string lo(ys), hi(ys);
  for (std::size_t i = 0; i < lo.size(); i++)
    if (lo[i] == 'X') { lo[i] = '0'; hi[i] = '9'; }
  const bool unspecified = lo != hi;
  long y0 = digits(lo, "the year"), y1 = digits(hi, "the year");
  if (negative) { y0 = -y0; y1 = -y1; std::swap(y0, y1); }
  if (dash1 == std::string_view::npos) {
    p.precision = 1;
    if (c == Calendar::julian_march) {
      // The legal year Y runs from 25 March of Y to 24 March of Y + 1.
      p.first = jdn_julian(y0, 3, 25);
      p.last = jdn_julian(y1 + 1, 3, 24);
    } else {
      p.first = jdn_of(c, y0, 1, 1);
      p.last = jdn_of(c, y1, 12, 31);
    }
    return p;
  }
  if (unspecified) throw std::invalid_argument("unspecified digits are allowed only in a year alone");
  const std::string_view rest = s.substr(dash1 + 1);
  const auto dash2 = rest.find('-');
  const long month = digits(rest.substr(0, dash2), "the month");
  if (month < 1 || month > 12) throw std::invalid_argument("no month " + std::to_string(month));
  const int m = static_cast<int>(month);
  if (dash2 == std::string_view::npos) {
    p.precision = 2;
    // In Old Style January and February of the legal year Y fall in the year Y + 1.
    const long historical = c == Calendar::julian_march && m < 3 ? y0 + 1 : y0;
    const int dim = days_in_month(c == Calendar::gregorian ? c : Calendar::julian, historical, m);
    if (c == Calendar::julian_march && m == 3) {
      // March of a legal year is split by its new year: 25-31 March of Y, 1-24 March of Y + 1.
      p.first = jdn_julian(y0, 3, 25);
      p.last = jdn_julian(y0 + 1, 3, 24);
    } else {
      p.first = jdn_of(c, y0, m, 1);
      p.last = jdn_of(c, y0, m, dim);
    }
    return p;
  }
  const long day = digits(rest.substr(dash2 + 1), "the day");
  const int dim = days_in_month(c == Calendar::gregorian ? c : Calendar::julian, c == Calendar::julian_march && (m < 3 || (m == 3 && day < 25)) ? y0 + 1 : y0, m);
  if (day < 1 || day > dim) throw std::invalid_argument("no day " + std::to_string(day) + " in that month");
  p.precision = 3;
  p.first = p.last = jdn_of(c, y0, m, static_cast<int>(day));
  return p;
}

}  // namespace detail

// Parses a date as written in a calendar; throws std::invalid_argument with the reason.
inline Date parse_date(std::string_view text, Calendar c = Calendar::gregorian) {
  Date d;
  d.text = std::string(text);
  d.calendar = c;
  if (text.empty()) throw std::invalid_argument("the date is empty");
  // A slash after the first character separates an interval (a leading '-' is a sign).
  const auto slash = text.find('/', 1);
  const detail::Part a = detail::part(text.substr(0, slash), c);
  d.precision = a.precision;
  d.uncertain = a.uncertain;
  d.approximate = a.approximate;
  d.first = a.first;
  d.last = a.last;
  if (slash != std::string_view::npos) {
    const detail::Part b = detail::part(text.substr(slash + 1), c);
    if (b.last < a.first) throw std::invalid_argument("the interval ends before it begins");
    d.interval = true;
    d.last = b.last;
    d.uncertain = d.uncertain || b.uncertain;
    d.approximate = d.approximate || b.approximate;
  }
  return d;
}

// Parses a reader's date in the historians' convention: Gregorian if it begins on or after the
// reform, Julian otherwise; an interval is read in the calendar of its beginning.
inline Date parse_historical(std::string_view text) {
  const Date g = parse_date(text, Calendar::gregorian);
  return g.first >= GREGORIAN_REFORM ? g : parse_date(text, Calendar::julian);
}

}  // namespace history
