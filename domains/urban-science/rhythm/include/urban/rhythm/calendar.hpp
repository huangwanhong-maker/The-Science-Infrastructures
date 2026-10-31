// Calendars as angles, and the lattice of their combinations (theory/street-oceanography.md §3.1–3.3).
//
// Three calendars drive a street periodically, each an angle advancing at its frequency in civil time:
//   - the civil week, ω_W = 1/168 h⁻¹, from Monday 00:00 civil; its 7th harmonic is the civil day (Remark 3.3:
//     the day and the week are one calendar, so the day is no independent generator);
//   - the synodic month, ω_M = 1/(29.530588861 d), from the mean new moon of k = 0 in Meeus (1998,
//     Astronomical Algorithms, eq. 49.1), JDE 2451550.09766 (2000 January 6, 14:20 TT; TT is taken as UTC,
//     64 s being 2.5 × 10⁻⁵ of a cycle);
//   - the tropical year, ω_Y = 1/(365.24219 d) (Meeus 1998, ch. 27), from J2000.0 (2000 January 1, 12:00).
//
// A combination k = (k_W, k_M, k_Y) ∈ ℤ³ has the frequency k·ω. The week's index is read in days and weeks:
// k_W = 7q + r with |r| ≤ 3 (unique, 7 being odd), named D<q>±W<|r|> — W8 is "D1+W1", the day's upper weekly
// sideband — and its order is |q| + |r|, the number of steps of a day or a week. The order of k is that plus
// |k_M| + |k_Y|. Names join the parts with signs: (7, 0, −1) is "D1−Y1", the year's lower sideband of the day.
//
// Proposition 3.9: with two or more incommensurate generators the lattice image is dense, so a frequency is
// "explained" only up to an order d. chance_rate() is the fraction of [0, f_max] within δ = R/T of some
// combination of order ≤ d with positive frequency: the probability that a frequency drawn at random would
// be explained by chance — the discriminator's false-explanation rate. Exact, as a union of intervals.
//
// Times are Unix seconds (UTC); civil time is UTC plus the zone's offset at the instant (record.hpp).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace urban::rhythm {

inline constexpr double hour_s = 3600.0;
inline constexpr double day_s = 86400.0;
inline constexpr double week_h = 168.0;
inline constexpr double synodic_days = 29.530588861;                          // Meeus 1998, eq. 49.1
inline constexpr double tropical_days = 365.24219;                            // mean tropical year
inline constexpr double new_moon_unix = (2451550.09766 - 2440587.5) * day_s;  // mean new moon, k = 0
inline constexpr double j2000_unix = (2451545.0 - 2440587.5) * day_s;
inline constexpr double monday_unix = 4 * day_s;                              // 1970-01-05 00:00, a Monday

// The generators, in cycles per hour: week, synodic month, tropical year.
inline std::array<double, 3> generators() { return {1.0 / week_h, 1.0 / (synodic_days * 24), 1.0 / (tropical_days * 24)}; }

inline double fraction(double x) { return x - std::floor(x); }
inline double week_angle(double civil_s) { return fraction((civil_s - monday_unix) / (week_h * hour_s)); }
inline double month_angle(double utc_s) { return fraction((utc_s - new_moon_unix) / (synodic_days * day_s)); }
inline double year_angle(double utc_s) { return fraction((utc_s - j2000_unix) / (tropical_days * day_s)); }

struct Combination {
  std::array<int, 3> k{0, 0, 0};
  double frequency = 0;            // cycles per hour
  int order = 0;                   // |q| + |r| + |k_M| + |k_Y|
  std::string name;
};

inline double frequency(const std::array<int, 3>& k) {
  const auto g = generators();
  return k[0] * g[0] + k[1] * g[1] + k[2] * g[2];
}

// k_W = 7q + r with r ∈ [−3, 3].
inline std::pair<int, int> days_and_weeks(int kw) {
  int q = static_cast<int>(std::lround(kw / 7.0));
  return {q, kw - 7 * q};
}

inline int order_of(const std::array<int, 3>& k) {
  const auto [q, r] = days_and_weeks(k[0]);
  return std::abs(q) + std::abs(r) + std::abs(k[1]) + std::abs(k[2]);
}

inline std::string name_of(const std::array<int, 3>& k) {
  std::string s;
  auto part = [&](int c, const char* letter) {
    if (c == 0) return;
    if (!s.empty()) s += c > 0 ? "+" : "−";
    else if (c < 0) s += "−";
    s += letter + std::to_string(std::abs(c));
  };
  const auto [q, r] = days_and_weeks(k[0]);
  part(q, "D");
  part(r, "W");
  part(k[1], "M");
  part(k[2], "Y");
  return s.empty() ? "Z0" : s;
}

inline Combination combination(std::array<int, 3> k) {
  if (frequency(k) < 0) k = {-k[0], -k[1], -k[2]};
  return {k, frequency(k), order_of(k), name_of(k)};
}

// Every combination of the generators in the mask (week, month, year) with positive frequency and order
// 1…max_order, at most max_day days (|q|), sorted by frequency.
inline std::vector<Combination> lattice(std::array<bool, 3> mask, int max_order, int max_day) {
  std::vector<Combination> out;
  const int qd = mask[0] ? max_day : 0, rd = mask[0] ? 3 : 0, bd = mask[1] ? max_order : 0, cd = mask[2] ? max_order : 0;
  for (int q = -qd; q <= qd; q++)
    for (int r = -rd; r <= rd; r++)
      for (int b = -bd; b <= bd; b++)
        for (int c = -cd; c <= cd; c++) {
          const int order = std::abs(q) + std::abs(r) + std::abs(b) + std::abs(c);
          if (order == 0 || order > max_order) continue;
          const std::array<int, 3> k{7 * q + r, b, c};
          if (frequency(k) <= 0) continue;
          out.push_back(combination(k));
        }
  std::sort(out.begin(), out.end(), [](const Combination& x, const Combination& y) { return x.frequency < y.frequency; });
  return out;
}

// The fraction of [0, f_max] within δ of a combination of order 1…d (at most max_day days) of the
// generators in the mask. Exact, by merging intervals.
inline double chance_rate(std::array<bool, 3> mask, int d, double delta, double f_max, int max_day = 1000) {
  std::vector<std::pair<double, double>> iv;
  for (const auto& c : lattice(mask, d, std::min(d, max_day))) {
    if (c.frequency - delta > f_max) continue;
    iv.emplace_back(std::max(0.0, c.frequency - delta), std::min(f_max, c.frequency + delta));
  }
  std::sort(iv.begin(), iv.end());
  double covered = 0, lo = 0, hi = 0;
  bool open = false;
  for (const auto& [a, b] : iv) {
    if (!open || a > hi) {
      if (open) covered += hi - lo;
      lo = a, hi = b, open = true;
    } else {
      hi = std::max(hi, b);
    }
  }
  if (open) covered += hi - lo;
  return covered / f_max;
}

struct Explanation {
  bool found = false;
  Combination combination;
  double distance = 0;            // |f − k·ω|, cycles per hour
};

// The combination of least order (then nearest) within δ of f, among those of order ≤ d of the generators
// in the mask.
inline Explanation explain(double f, std::array<bool, 3> mask, int d, double delta, int max_day = 1000) {
  Explanation best;
  best.distance = std::numeric_limits<double>::infinity();
  for (const auto& c : lattice(mask, d, std::min(d, max_day))) {
    const double dist = std::fabs(c.frequency - f);
    if (dist > delta) continue;
    if (!best.found || c.order < best.combination.order || (c.order == best.combination.order && dist < best.distance)) {
      best.found = true;
      best.combination = c;
      best.distance = dist;
    }
  }
  return best;
}

// ---- the week free -------------------------------------------------------------------------------------
//
// A forcing that is any function on the week's circle (hours that differ by weekday, a daily cycle of any
// shape) has every line m/168 h⁻¹ in its support (Remark 3.3): all the week's lines are forced, whatever
// their order in days and weeks. Explaining a line then costs only the steps across the incommensurate
// calendars: its cross order is |k_M| + |k_Y|, and a line of the week alone has cross order 0. These are the
// versions the spectrum uses when the week is a declared driving.

inline int cross_order(const Combination& c) { return std::abs(c.k[1]) + std::abs(c.k[2]); }

// Every combination of positive frequency up to f_max whose cross order is at most d, the week's index free.
inline std::vector<Combination> lattice_free_week(std::array<bool, 3> mask, int d, double f_max) {
  std::vector<Combination> out;
  const auto g = generators();
  const int bd = mask[1] ? d : 0, cd = mask[2] ? d : 0;
  const int wmax = mask[0] ? static_cast<int>(std::ceil(f_max / g[0])) + 2 : 0;
  for (int b = -bd; b <= bd; b++)
    for (int c = -cd; c <= cd; c++) {
      if (std::abs(b) + std::abs(c) > d) continue;
      for (int a = -wmax; a <= wmax; a++) {
        if (a == 0 && b == 0 && c == 0) continue;
        const std::array<int, 3> k{a, b, c};
        const double f = frequency(k);
        if (f <= 0 || f > f_max + 1e-12) continue;
        out.push_back(combination(k));
      }
    }
  std::sort(out.begin(), out.end(), [](const Combination& x, const Combination& y) { return x.frequency < y.frequency; });
  return out;
}

inline double chance_rate_free_week(std::array<bool, 3> mask, int d, double delta, double f_max) {
  std::vector<std::pair<double, double>> iv;
  for (const auto& c : lattice_free_week(mask, d, f_max + delta)) iv.emplace_back(std::max(0.0, c.frequency - delta), std::min(f_max, c.frequency + delta));
  std::sort(iv.begin(), iv.end());
  double covered = 0, lo = 0, hi = 0;
  bool open = false;
  for (const auto& [a, b] : iv) {
    if (a >= b) continue;
    if (!open || a > hi) {
      if (open) covered += hi - lo;
      lo = a, hi = b, open = true;
    } else {
      hi = std::max(hi, b);
    }
  }
  if (open) covered += hi - lo;
  return covered / f_max;
}

// The combination of least cross order (then nearest) within δ of f.
inline Explanation explain_free_week(double f, std::array<bool, 3> mask, int d, double delta) {
  Explanation best;
  best.distance = std::numeric_limits<double>::infinity();
  for (const auto& c : lattice_free_week(mask, d, f + delta)) {
    const double dist = std::fabs(c.frequency - f);
    if (dist > delta) continue;
    const int o = cross_order(c), ob = best.found ? cross_order(best.combination) : 1 << 20;
    if (!best.found || o < ob || (o == ob && dist < best.distance)) {
      best.found = true;
      best.combination = c;
      best.distance = dist;
    }
  }
  return best;
}

// The driving a line belongs to (Definition 4.2): the month if it involves the month, else the year if it
// involves the year, else the day for a multiple of the day, else the week.
inline std::string driving(const Combination& c) {
  if (c.k[1] != 0) return "month";
  if (c.k[2] != 0) return "year";
  return days_and_weeks(c.k[0]).second == 0 ? "day" : "week";
}

}  // namespace urban::rhythm
