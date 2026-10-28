// Calendar dates and times of day in a time scale, read and written as ISO 8601.
//
// A date and time in a scale is a Julian date; in UTC a day that ends with a leap second has
// 86 401 seconds, so 23:59:60 exists on it and the fraction of the day is scaled to its length
// (SOFA's iauDtf2d and iauD2dtf). Text is ISO 8601 (ISO 8601-1:2019), extended format:
//
//   2026-10-28T12:34:56.789Z           UTC (Z)
//   2026-10-28T12:34:56.789 TT         another scale, named after a space
//   2026-10-28 12:00 TDB               a space for the T, seconds optional
//   2026-10-28                          0h
//   JD 2451545.0 TT,  MJD 60000.5 UTC  Julian dates
//
// A time without a scale is refused unless the caller names a default: a time value always
// states its scale (CLAUDE.md §4.1).

#pragma once

#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include <metrology/time/julian.hpp>
#include <metrology/time/scales.hpp>

namespace metrology::time {

struct DateTime {
  int year = 2000, month = 1, day = 1, hour = 0, minute = 0;
  double second = 0;
};

// The Julian date of a calendar date and time in a scale (SOFA's iauDtf2d).
inline JulianDate from_datetime(const DateTime& d, Scale scale, const TimeScales& ts) {
  const JulianDate z = calendar_to_jd(d.year, d.month, d.day);
  double day = DAY, seclim = 60;
  if (scale == Scale::UTC) {
    const LeapSeconds& leap = ts.leap_seconds();
    const double dat0 = leap.tai_minus_utc(d.year, d.month, d.day, 0.0);
    const double dat12 = leap.tai_minus_utc(d.year, d.month, d.day, 0.5);
    const CalendarDate next = jd_to_calendar({z.d1 + z.d2, 1.5});
    const double dat24 = leap.tai_minus_utc(next.year, next.month, next.day, 0.0);
    const double dleap = dat24 - (2 * dat12 - dat0);
    day += dleap;
    if (d.hour == 23 && d.minute == 59) seclim += dleap;
  }
  if (d.hour < 0 || d.hour > 23) throw std::invalid_argument("metrology::time: hour out of range");
  if (d.minute < 0 || d.minute > 59) throw std::invalid_argument("metrology::time: minute out of range");
  if (d.second < 0 || d.second >= seclim) throw std::invalid_argument("metrology::time: second out of range for that minute");
  const double time = (60.0 * double(60 * d.hour + d.minute) + d.second) / day;
  return {z.d1 + z.d2, time};
}

// The calendar date and time of a Julian date in a scale, rounded to `decimals` of a second
// (SOFA's iauD2dtf; 23:59:60 on a UTC leap second).
inline DateTime to_datetime(const JulianDate& t, Scale scale, const TimeScales& ts, int decimals = 3) {
  if (decimals < 0 || decimals > 9) throw std::invalid_argument("metrology::time: 0 to 9 decimals of a second");
  CalendarDate c = jd_to_calendar(t);
  bool leap = false;
  if (scale == Scale::UTC) {
    const LeapSeconds& ls = ts.leap_seconds();
    const double dat0 = ls.tai_minus_utc(c.year, c.month, c.day, 0.0);
    const double dat12 = ls.tai_minus_utc(c.year, c.month, c.day, 0.5);
    const CalendarDate next = jd_to_calendar({t.d1 + 1.5, t.d2 - c.fraction});
    const double dat24 = ls.tai_minus_utc(next.year, next.month, next.day, 0.0);
    const double dleap = dat24 - (2 * dat12 - dat0);
    leap = std::fabs(dleap) > 0.5;
    if (leap) c.fraction += c.fraction * dleap / DAY;
  }
  // The time of day, rounded to the resolution (SOFA's iauD2tf).
  const double rs = std::pow(10.0, decimals);
  double a = dnint(rs * c.fraction * DAY);
  const double rm = rs * 60, rh = rm * 60;
  double ah = std::trunc(a / rh);
  a -= ah * rh;
  double am = std::trunc(a / rm);
  a -= am * rm;
  double as = a / rs;
  DateTime d{c.year, c.month, c.day, int(ah), int(am), as};
  if (d.hour > 23) {
    const CalendarDate next = jd_to_calendar({t.d1 + 1.5, t.d2 - c.fraction});
    if (!leap || d.second >= 1) {
      d = {next.year, next.month, next.day, 0, 0, 0.0};
    } else {
      d.hour = 23, d.minute = 59, d.second = 60 + d.second;  // within the leap second itself
    }
  }
  return d;
}

// ISO 8601 text of a Julian date in a scale: "…Z" for UTC, " TT" and the like otherwise.
inline std::string format_iso(const JulianDate& t, Scale scale, const TimeScales& ts, int decimals = 3) {
  const DateTime d = to_datetime(t, scale, ts, decimals);
  char seconds[32];
  if (decimals == 0)
    std::snprintf(seconds, sizeof seconds, "%02d", int(std::lround(d.second)));
  else
    std::snprintf(seconds, sizeof seconds, "%0*.*f", decimals + 3, decimals, d.second);
  char buf[96];
  std::snprintf(buf, sizeof buf, "%04d-%02d-%02dT%02d:%02d:%s", d.year, d.month, d.day, d.hour, d.minute, seconds);
  std::string s(buf);
  return scale == Scale::UTC ? s + "Z" : s + " " + std::string(name(scale));
}

// A date and time with its scale, from ISO 8601 or a Julian date; a missing scale is the
// default if one is given, and an error if not.
inline std::pair<JulianDate, Scale> parse_time(std::string_view text, const TimeScales& ts, std::optional<Scale> default_scale = std::nullopt) {
  std::string s(text);
  auto trim = [](std::string x) {
    const auto a = x.find_first_not_of(" \t"), b = x.find_last_not_of(" \t");
    return a == std::string::npos ? std::string() : x.substr(a, b - a + 1);
  };
  s = trim(s);
  auto scale_or_default = [&](const std::string& rest) {
    const std::string r = trim(rest);
    if (r.empty()) {
      if (!default_scale) throw std::invalid_argument("metrology::time: the time \"" + std::string(text) + "\" names no scale");
      return *default_scale;
    }
    return parse_scale(r);
  };
  // Julian dates.
  for (const std::string prefix : {"JD", "MJD"}) {
    if (s.rfind(prefix + " ", 0) == 0) {
      std::size_t used = 0;
      const std::string rest = s.substr(prefix.size() + 1);
      const double v = std::stod(rest, &used);
      const Scale sc = scale_or_default(rest.substr(used));
      return {prefix == "JD" ? JulianDate{v, 0.0} : JulianDate::from_mjd(v), sc};
    }
  }
  // ISO 8601: YYYY-MM-DD[(T| )hh:mm[:ss[.fff]]][Z| scale]
  DateTime d;
  int consumed = 0;
  if (std::sscanf(s.c_str(), "%d-%d-%d%n", &d.year, &d.month, &d.day, &consumed) != 3)
    throw std::invalid_argument("metrology::time: not an ISO 8601 date: " + std::string(text));
  std::string rest = s.substr(consumed);
  if (!rest.empty() && (rest[0] == 'T' || rest[0] == ' ') && rest.size() > 1 && std::isdigit(static_cast<unsigned char>(rest[1]))) {
    int used = 0;
    if (std::sscanf(rest.c_str() + 1, "%d:%d%n", &d.hour, &d.minute, &used) != 2)
      throw std::invalid_argument("metrology::time: not an ISO 8601 time: " + std::string(text));
    rest = rest.substr(1 + used);
    if (!rest.empty() && rest[0] == ':') {
      std::size_t n = 0;
      d.second = std::stod(rest.substr(1), &n);
      rest = rest.substr(1 + n);
    }
  }
  Scale sc;
  if (!rest.empty() && rest[0] == 'Z') {
    sc = Scale::UTC;
    if (!trim(rest.substr(1)).empty()) throw std::invalid_argument("metrology::time: text after Z: " + std::string(text));
  } else {
    sc = scale_or_default(rest);
  }
  return {from_datetime(d, sc, ts), sc};
}

}  // namespace metrology::time
