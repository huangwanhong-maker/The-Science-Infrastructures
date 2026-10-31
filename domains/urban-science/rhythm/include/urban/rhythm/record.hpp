// A record of counts in time (theory/street-oceanography.md §3, §4): values at equal steps of UTC from an
// instant, with gaps, and the zone's offsets that turn UTC into civil time.
//
// A street keeps the civil clock or the sun (Representation 3.12), so both times are kept: the record is laid
// out in UTC, the only time without repeated or missing hours; civil(t) adds the offset in force at t.
// offsets are the zone's transitions [(instant, hours)], sorted, the first at or before the record's start.
//
// spans() finds the stretches without gaps longer than max_gap steps; series() extracts one, filling short
// gaps by linear interpolation (stated with every result that uses it: a filled gap is not data).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <urban/rhythm/calendar.hpp>

namespace urban::rhythm {

// Days from 1970-01-01 of a proleptic Gregorian date (Hinnant's days_from_civil, chrono's algorithm).
inline long long days_from_civil(long long y, unsigned m, unsigned d) {
  y -= m <= 2;
  const long long era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + static_cast<long long>(doe) - 719468;
}

struct CivilDate { int year = 0; unsigned month = 0, day = 0; };

inline CivilDate civil_from_days(long long z) {
  z += 719468;
  const long long era = (z >= 0 ? z : z - 146096) / 146097;
  const unsigned doe = static_cast<unsigned>(z - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const long long y = static_cast<long long>(yoe) + era * 400;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  const unsigned d = doy - (153 * mp + 2) / 5 + 1;
  const unsigned m = mp < 10 ? mp + 3 : mp - 9;
  return {static_cast<int>(y + (m <= 2)), m, d};
}

// "YYYY-MM-DD", "YYYY-MM-DDTHH:MM[:SS]Z" (or without Z, read as UTC) to Unix seconds.
inline double parse_utc(const std::string& s) {
  int y = 0, mo = 0, d = 0, h = 0, mi = 0, se = 0;
  const int n = std::sscanf(s.c_str(), "%d-%d-%dT%d:%d:%d", &y, &mo, &d, &h, &mi, &se);
  if (n < 3 || mo < 1 || mo > 12 || d < 1 || d > 31) throw std::invalid_argument("urban::rhythm::parse_utc: not a date: " + s);
  return static_cast<double>(days_from_civil(y, static_cast<unsigned>(mo), static_cast<unsigned>(d))) * day_s + h * 3600.0 + mi * 60.0 + se;
}

inline std::string format_date(double unix_s) {
  const auto c = civil_from_days(static_cast<long long>(std::floor(unix_s / day_s)));
  char buf[16];
  std::snprintf(buf, sizeof buf, "%04d-%02u-%02u", c.year, c.month, c.day);
  return buf;
}

struct Record {
  std::string id;
  double t0 = 0;                                     // UTC of the first value, Unix s
  double step = 3600;                                // s
  std::vector<double> y;                             // NaN where missing
  std::vector<std::pair<double, double>> offsets;    // (UTC instant, offset in hours) from that instant on

  std::size_t size() const { return y.size(); }
  double time(std::size_t i) const { return t0 + step * static_cast<double>(i); }
  double offset(double t) const {
    double o = offsets.empty() ? 0.0 : offsets.front().second;
    for (const auto& [at, hours] : offsets) {
      if (at > t) break;
      o = hours;
    }
    return o;
  }
  double civil(double t) const { return t + 3600.0 * offset(t); }
  // Days since 1970-01-01 of the civil date at t, and its weekday (0 = Monday).
  long long civil_day(double t) const { return static_cast<long long>(std::floor(civil(t) / day_s)); }
  int weekday(double t) const { return static_cast<int>(((civil_day(t) - 4) % 7 + 7) % 7); }
  double civil_hour(double t) const { const double c = civil(t) / 3600.0; return c - 24.0 * std::floor(c / 24.0); }
};

struct Span { std::size_t begin = 0, end = 0; std::size_t filled = 0; };   // [begin, end), gaps filled inside

// Stretches of at least min_length steps whose gaps are each at most max_gap steps, longest first.
inline std::vector<Span> spans(const Record& r, std::size_t min_length, std::size_t max_gap = 0) {
  std::vector<Span> out;
  const std::size_t n = r.size();
  std::size_t i = 0;
  while (i < n) {
    while (i < n && !std::isfinite(r.y[i])) i++;
    if (i >= n) break;
    Span s{i, i, 0};
    std::size_t last = i, j = i + 1;
    while (j < n) {
      if (std::isfinite(r.y[j])) { last = j++; continue; }
      std::size_t k = j;
      while (k < n && !std::isfinite(r.y[k])) k++;
      if (k >= n || k - j > max_gap) break;
      s.filled += k - j;
      j = k;
    }
    s.end = last + 1;
    if (s.end - s.begin >= min_length) out.push_back(s);
    i = last + 1;
  }
  std::sort(out.begin(), out.end(), [](const Span& a, const Span& b) { return a.end - a.begin > b.end - b.begin; });
  return out;
}

// The values of a span with its short gaps filled linearly; transform ln(1 + y) when log1p.
inline std::vector<double> series(const Record& r, const Span& s, bool log1p = false) {
  std::vector<double> out(r.y.begin() + static_cast<std::ptrdiff_t>(s.begin), r.y.begin() + static_cast<std::ptrdiff_t>(s.end));
  std::size_t i = 0;
  while (i < out.size()) {
    if (std::isfinite(out[i])) { i++; continue; }
    std::size_t k = i;
    while (k < out.size() && !std::isfinite(out[k])) k++;
    const double a = i > 0 ? out[i - 1] : (k < out.size() ? out[k] : 0.0);
    const double b = k < out.size() ? out[k] : a;
    for (std::size_t j = i; j < k; j++) out[j] = a + (b - a) * static_cast<double>(j - i + 1) / static_cast<double>(k - i + 1);
    i = k;
  }
  if (log1p)
    for (double& v : out) v = std::log1p(std::max(0.0, v));
  return out;
}

// The finite values of the record (or of [begin, end)), their UTC times in hours since t_ref, transformed.
struct Points { std::vector<double> t, y; };
inline Points points(const Record& r, std::size_t begin, std::size_t end, double t_ref, bool log1p = false) {
  Points p;
  for (std::size_t i = begin; i < std::min(end, r.size()); i++) {
    if (!std::isfinite(r.y[i])) continue;
    p.t.push_back((r.time(i) - t_ref) / 3600.0);
    p.y.push_back(log1p ? std::log1p(std::max(0.0, r.y[i])) : r.y[i]);
  }
  return p;
}

}  // namespace urban::rhythm
