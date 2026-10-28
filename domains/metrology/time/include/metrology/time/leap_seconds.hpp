// TAI − UTC: the leap seconds since 1972, and the drifting offsets of 1960–1971.
//
// From 1972, UTC differs from TAI by a whole number of seconds, changed by the IERS in its
// Bulletin C at the end of a June or December. From 1960 to 1971 UTC ran at a rate offset from
// TAI's, with steps, and TAI − UTC was a linear function of the date (Explanatory Supplement
// to the Astronomical Almanac, 3rd ed., 2013, Table 3.1; the values below are SOFA's iauDat,
// release 2023, as ERFA carries them). A table carries the date up to which it is known to be
// complete: no leap second can have been decided before that date without its being in the
// table. After it, TAI − UTC is the last value, and a caller that cares is told so.
//
// The built-in table ends with the leap second of 2017 January 1 (TAI − UTC = 37 s); IANA's
// copy of the IERS list, retrieved 2026-10-28, declares itself valid until 2027 June 28 with
// no leap second since. A newer list read with `parse_iana` replaces it.

#pragma once

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <metrology/time/julian.hpp>

namespace metrology::time {

class LeapSeconds {
 public:
  struct Change {
    int year, month;   // the value holds from 0h UTC on the first of this month
    double offset;     // TAI − UTC, s (for 1960–1971, at the reference date)
    double drift_mjd;  // 1960–1971: TAI − UTC += (MJD − drift_mjd) × rate
    double rate;       // s per day; 0 from 1972
  };

  // SOFA's iauDat (release 2023) and the IERS list valid until 2027-06-28 (MJD 61584).
  static LeapSeconds built_in() {
    LeapSeconds t;
    t.changes_ = {
        {1960, 1, 1.4178180, 37300, 0.0012960}, {1961, 1, 1.4228180, 37300, 0.0012960}, {1961, 8, 1.3728180, 37300, 0.0012960},
        {1962, 1, 1.8458580, 37665, 0.0011232}, {1963, 11, 1.9458580, 37665, 0.0011232}, {1964, 1, 3.2401300, 38761, 0.0012960},
        {1964, 4, 3.3401300, 38761, 0.0012960}, {1964, 9, 3.4401300, 38761, 0.0012960}, {1965, 1, 3.5401300, 38761, 0.0012960},
        {1965, 3, 3.6401300, 38761, 0.0012960}, {1965, 7, 3.7401300, 38761, 0.0012960}, {1965, 9, 3.8401300, 38761, 0.0012960},
        {1966, 1, 4.3131700, 39126, 0.0025920}, {1968, 2, 4.2131700, 39126, 0.0025920},
    };
    const int steps[][3] = {{1972, 1, 10}, {1972, 7, 11}, {1973, 1, 12}, {1974, 1, 13}, {1975, 1, 14}, {1976, 1, 15}, {1977, 1, 16},
                            {1978, 1, 17}, {1979, 1, 18}, {1980, 1, 19}, {1981, 7, 20}, {1982, 7, 21}, {1983, 7, 22}, {1985, 7, 23},
                            {1988, 1, 24}, {1990, 1, 25}, {1991, 1, 26}, {1992, 7, 27}, {1993, 7, 28}, {1994, 7, 29}, {1996, 1, 30},
                            {1997, 7, 31}, {1999, 1, 32}, {2006, 1, 33}, {2009, 1, 34}, {2012, 7, 35}, {2015, 7, 36}, {2017, 1, 37}};
    for (const auto& s : steps) t.changes_.push_back({s[0], s[1], double(s[2]), 0, 0});
    t.valid_until_mjd_ = 61584;  // 2027-06-28
    t.source_ = "built in: SOFA iauDat (2023) and the IERS list via IANA, valid until 2027-06-28";
    return t;
  }

  // IANA's leap-seconds.list (https://data.iana.org/time-zones/data/leap-seconds.list): lines of
  // NTP seconds (from 1900-01-01) and TAI − UTC; "#@" gives the expiry. The pre-1972 drifting
  // offsets, which the list does not hold, are kept from the built-in table.
  static LeapSeconds parse_iana(std::string_view text) {
    LeapSeconds t = built_in();
    std::vector<Change> drifting;
    for (const Change& c : t.changes_)
      if (c.year < 1972) drifting.push_back(c);
    std::vector<Change> steps;
    double expiry = 0;
    std::istringstream in{std::string(text)};
    std::string line;
    while (std::getline(in, line)) {
      if (line.rfind("#@", 0) == 0) {
        expiry = std::stod(line.substr(2));
        continue;
      }
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      double ntp = 0, offset = 0;
      if (!(f >> ntp >> offset)) continue;
      const double mjd = 15020.0 + ntp / DAY;  // NTP epoch 1900-01-01 = MJD 15020
      const CalendarDate c = jd_to_calendar(JulianDate::from_mjd(mjd));
      if (c.day != 1) throw std::invalid_argument("leap-seconds.list: a change not on the first of a month");
      steps.push_back({c.year, c.month, offset, 0, 0});
    }
    if (steps.empty() || expiry == 0) throw std::invalid_argument("leap-seconds.list: no entries or no expiry date (#@)");
    t.changes_ = drifting;
    t.changes_.insert(t.changes_.end(), steps.begin(), steps.end());
    t.valid_until_mjd_ = 15020.0 + expiry / DAY;
    t.source_ = "IANA leap-seconds.list, valid until MJD " + std::to_string(int(t.valid_until_mjd_));
    return t;
  }

  // TAI − UTC (s) at a UTC calendar date and fraction of its day; SOFA's iauDat. Throws before
  // 1960, when UTC did not exist.
  double tai_minus_utc(int year, int month, int day, double fraction = 0) const {
    if (year < changes_.front().year) throw std::domain_error("metrology::time: UTC is not defined before 1960");
    const int m = 12 * year + month;
    int i = int(changes_.size()) - 1;
    while (i >= 0 && m < 12 * changes_[i].year + changes_[i].month) i--;
    if (i < 0) throw std::domain_error("metrology::time: UTC is not defined before 1960 January 1");
    const Change& c = changes_[i];
    double da = c.offset;
    if (c.rate != 0) da += (calendar_to_jd(year, month, day).mjd() + fraction - c.drift_mjd) * c.rate;
    return da;
  }
  // The MJD (UTC) up to which the table is known to be complete.
  double valid_until_mjd() const { return valid_until_mjd_; }
  const std::string& source() const { return source_; }
  const std::vector<Change>& changes() const { return changes_; }

 private:
  std::vector<Change> changes_;
  double valid_until_mjd_ = 0;
  std::string source_;
};

}  // namespace metrology::time
