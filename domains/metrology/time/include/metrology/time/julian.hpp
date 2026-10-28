// Julian dates in two parts, and the Gregorian calendar.
//
// An instant of a time scale is held as a Julian date in two parts, jd = d1 + d2: a double
// alone resolves about 40 µs near the present (its 53 bits spread over 2.46 × 10⁶ days), the
// pair to well under a nanosecond. The split is free (2451545.0 + 0.25, or 2451545.25 + 0);
// the conversions keep the larger part as it is and carry the arithmetic in the smaller, as
// SOFA does (IAU SOFA Board, SOFA Time Scale and Calendar Tools, 2021, §2).
//
// The calendar is the proleptic Gregorian, from −4799 (Fliegel and Van Flandern 1968, Comm.
// ACM 11, 657; Explanatory Supplement to the Astronomical Almanac, 3rd ed., 2013, §15.11).
// A day of the calendar runs from 0h, which is a Julian date ending in .5; the fraction of a
// day is measured from there. The history domain holds the Julian calendar and the others.

#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <string>

namespace metrology::time {

inline constexpr double J2000 = 2451545.0;     // JD of J2000.0 (TT): 2000 January 1, 12h
inline constexpr double MJD_ZERO = 2400000.5;  // JD of MJD 0: 1858 November 17, 0h
inline constexpr double DAY = 86400.0;         // seconds per day
inline constexpr double JULIAN_CENTURY = 36525.0;
inline constexpr double JULIAN_MILLENNIUM = 365250.0;

struct JulianDate {
  double d1 = 0, d2 = 0;  // jd = d1 + d2

  // From an MJD: the whole days go into d1 with the MJD zero point (exactly, as a half-integer
  // Julian date), the fraction into d2, so that the fraction keeps its full precision.
  static JulianDate from_mjd(double mjd, double fraction = 0) {
    const double whole = std::floor(mjd);
    return {MJD_ZERO + whole, (mjd - whole) + fraction};
  }
  // One double: adequate for display and where 40 µs do not matter.
  double value() const { return d1 + d2; }
  double mjd() const { return (d1 - MJD_ZERO) + d2; }
  // Days since J2000.0, with the precision of the two parts.
  double since_j2000() const { return (d1 - J2000) + d2; }
  double centuries_since_j2000() const { return since_j2000() / JULIAN_CENTURY; }
  // The difference of two dates, in days, without forming either sum.
  double minus(const JulianDate& o) const { return (d1 - o.d1) + (d2 - o.d2); }
  // This date moved by an interval of seconds (added to the smaller part).
  JulianDate plus_seconds(double s) const {
    return std::fabs(d1) >= std::fabs(d2) ? JulianDate{d1, d2 + s / DAY} : JulianDate{d1 + s / DAY, d2};
  }
  JulianDate plus_days(double days) const {
    return std::fabs(d1) >= std::fabs(d2) ? JulianDate{d1, d2 + days} : JulianDate{d1 + days, d2};
  }
  // The same date with the day number in d1 (ending in .5, a calendar day's 0h) and the
  // fraction of that day in d2 ∈ [0, 1).
  JulianDate normalised() const {
    const double big = std::fabs(d1) >= std::fabs(d2) ? d1 : d2, small = std::fabs(d1) >= std::fabs(d2) ? d2 : d1;
    const double day = std::floor(big - 0.5) + 0.5;
    double f = (big - day) + small;
    const double carry = std::floor(f);
    return {day + carry, f - carry};
  }
};

struct CalendarDate {
  int year = 2000, month = 1, day = 1;
  double fraction = 0;  // of the day from 0h, in [0, 1)
};

// Rounding to the nearest integer, halves away from zero (SOFA's dnint).
inline double dnint(double a) { return a < 0 ? std::ceil(a - 0.5) : std::floor(a + 0.5); }

inline bool is_leap_year(int y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; }
inline int days_in_month(int y, int m) {
  static const int days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (m < 1 || m > 12) throw std::invalid_argument("metrology::time: month " + std::to_string(m) + " is not 1 to 12");
  return m == 2 && is_leap_year(y) ? 29 : days[m - 1];
}

// The Julian date of 0h of a Gregorian date (SOFA's iauCal2jd), whole in d1 (a half-integer,
// exact), so that a time of day added to d2 keeps its full precision. mjd() gives SOFA's MJD.
inline JulianDate calendar_to_jd(int iy, int im, int id) {
  if (iy < -4799) throw std::domain_error("metrology::time: years before -4799 are outside the Gregorian algorithm");
  if (id < 1 || id > days_in_month(iy, im)) throw std::invalid_argument("metrology::time: no day " + std::to_string(id) + " in that month");
  const long my = (im - 14) / 12;
  const long iypmy = long(iy) + my;
  const double mjd = double((1461L * (iypmy + 4800L)) / 4L + (367L * long(im - 2 - 12 * my)) / 12L -
                            (3L * ((iypmy + 4900L) / 100L)) / 4L + long(id) - 2432076L);
  return {MJD_ZERO + mjd, 0.0};
}

// The Gregorian date and fraction of a day of a two-part Julian date; SOFA's iauJd2cal, with
// its compensated summation of the fractions so that the fraction is right to rounding.
inline CalendarDate jd_to_calendar(const JulianDate& t) {
  const double dj = t.d1 + t.d2;
  if (dj < -68569.5 || dj > 1e9) throw std::domain_error("metrology::time: the Julian date is outside the calendar's range");
  double d = dnint(t.d1);
  const double f1 = t.d1 - d;
  long jd = long(d);
  d = dnint(t.d2);
  const double f2 = t.d2 - d;
  jd += long(d);
  double s = 0.5, cs = 0;
  for (const double x : {f1, f2}) {
    const double tt = s + x;
    cs += std::fabs(s) >= std::fabs(x) ? (s - tt) + x : (x - tt) + s;
    s = tt;
    if (s >= 1) {
      jd++;
      s -= 1;
    }
  }
  double f = s + cs;
  cs = f - s;
  if (f < 0) {
    f = s + 1;
    cs += (1 - f) + s;
    s = f;
    f = s + cs;
    cs = f - s;
    jd--;
  }
  if (f - 1 >= -DBL_EPSILON / 4) {
    const double tt = s - 1;
    cs += (s - tt) - 1;
    s = tt;
    f = s + cs;
    if (-DBL_EPSILON / 2 < f) {
      jd++;
      f = std::max(f, 0.0);
    }
  }
  long l = jd + 68569L;
  const long n = (4L * l) / 146097L;
  l -= (146097L * n + 3L) / 4L;
  const long i = (4000L * (l + 1L)) / 1461001L;
  l -= (1461L * i) / 4L - 31L;
  const long k = (80L * l) / 2447L;
  CalendarDate c;
  c.day = int(l - (2447L * k) / 80L);
  l = k / 11L;
  c.month = int(k + 2L - 12L * l);
  c.year = int(100L * (n - 49L) + i + l);
  c.fraction = f;
  return c;
}

}  // namespace metrology::time
