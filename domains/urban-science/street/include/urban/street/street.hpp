// The street system (theory/street-oceanography.md §1; Huang, Street Oceanography, Definition 5.1): a
// segment of length L with its holdings, each at an arc length s along the street's axis, with a use from a
// finite set U and the hours its occupant declares, as intervals of the civil week.
//
// The week is the circle W = ℝ/168ℤ of hours, Monday 00:00 at 0. An interval [open, close) has
// 0 ≤ open < 168 and open < close ≤ open + 168; a close past 168 wraps into the next Monday (a bar open
// Sunday 22:00 to Monday 02:00 is [166, 170)). Hours are civil: they are what the occupant declares, in the
// zone of the street, daylight saving included.
//
// A holding without declared hours is kept (it is terrain) and is left out of every quantity of time.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace urban::street {

inline constexpr double week_hours = 168.0;

// An hour of the week reduced to [0, 168).
inline double week_angle(double hour) {
  double h = std::fmod(hour, week_hours);
  return h < 0 ? h + week_hours : h;
}

struct Holding {
  std::string id;                               // the source's identifier (an OSM element, a CLUE property)
  std::string name;
  double s = 0;                                 // arc length along the axis, m
  double offset = 0;                            // signed distance from the axis, m (left positive)
  std::size_t use = 0;                          // index into Street::uses
  std::vector<std::array<double, 2>> hours;     // [open, close) in hours of the week; empty: none declared
  bool declared = false;                        // whether hours were declared (open 24/7 is declared and full)
};

struct Street {
  std::string id, name;
  double length = 0;                            // m
  std::vector<std::string> uses;                // U, by index
  std::vector<Holding> holdings;
};

// Whether a holding is open at hour θ of the week (half-open intervals; wrapped intervals honoured).
inline bool open_at(const Holding& h, double theta) {
  const double t = week_angle(theta);
  for (const auto& [a, b] : h.hours) {
    if (t >= a && t < b) return true;
    if (b > week_hours && t + week_hours >= a && t + week_hours < b) return true;
  }
  return false;
}

// The hours a holding is open in a week.
inline double open_hours(const Holding& h) {
  double total = 0;
  for (const auto& [a, b] : h.hours) total += b - a;
  return total;
}

// Validate a street: positions within [0, L] (to a metre), uses in range, intervals well formed and
// disjoint after wrapping.
inline void validate(const Street& st) {
  if (!(st.length > 0)) throw std::invalid_argument("urban::street::validate: a street has positive length");
  for (const auto& h : st.holdings) {
    if (h.use >= st.uses.size()) throw std::invalid_argument("urban::street::validate: the use of " + h.id + " is not one of the street's");
    if (h.s < -1 || h.s > st.length + 1) throw std::invalid_argument("urban::street::validate: " + h.id + " lies off the street");
    for (const auto& [a, b] : h.hours)
      if (!(a >= 0 && a < week_hours && b > a && b <= a + week_hours))
        throw std::invalid_argument("urban::street::validate: an interval of " + h.id + " is not within the week");
  }
}

// Holdings converted from the day-and-minute form of the data library ([day, open_min, close_min], day 0 =
// Monday, minutes from that day's 00:00, close past 1440 for past midnight).
inline std::vector<std::array<double, 2>> week_intervals(const std::vector<std::array<double, 3>>& day_minutes) {
  std::vector<std::array<double, 2>> out;
  for (const auto& [day, open, close] : day_minutes) {
    if (day < 0 || day > 6 || open < 0 || close <= open) continue;
    const double a = 24 * day + open / 60.0, b = 24 * day + close / 60.0;
    if (b - a >= week_hours) { out.push_back({0, week_hours}); continue; }
    out.push_back({a, b});
  }
  return out;
}

}  // namespace urban::street
