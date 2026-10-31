// The tide of a street (theory/street-oceanography.md §4.6): a record of counts read as a tide gauge's.
//
// The prediction. A port's tide is predicted from its harmonic constants (Foreman 1977, Manual for Tidal
// Heights Analysis and Prediction, Pacific Marine Science Report 77-10; Pugh and Woodworth 2014, Sea-Level
// Science, ch. 4). A street's count is predicted in the same way, from its constituents on the calendars'
// lattice (calendar.hpp; the street constituents of analysis.hpp): least squares on the count itself, people
// an hour, without a trend. Each hourly count is placed at the middle of its hour. Time is civil, in hours
// from the week's origin (Monday 1970-01-05 00:00 civil), so the day and the week keep the civil clock, as
// the measured streets do (the civil share, analysis.hpp). The year's constituents are fitted on the same
// time. Over a year an hour of offset moves them by 1.1 × 10⁻⁴ of a cycle.
//
// The prediction is the calendars' alone. What they do not drive (holidays, the weather, events) is left in
// the residual, observed − predicted. This is the street's surge, as a storm surge is the sea's departure from
// its astronomical tide (Pugh and Woodworth 2014, ch. 6). A prediction below zero is reported as zero.
//
// Readings of the predicted mean week, the week's constituents alone (k_M = k_Y = 0), which the year's
// constituents modulate:
//   - high water and low water of the mean weekday (the pointwise mean of the curves of Monday to Friday)
//     and of the mean weekend day, and the civil hours at which they fall. The hour of high water plays the
//     part of a port's establishment, the interval at which its high water follows the Moon's transit;
//   - the range, high water − low water;
//   - the form. Tides are classified by the ratio of their diurnal species to their semidiurnal one,
//     F = (K1 + O1)/(M2 + S2) (Courtier 1938, The classification of tides, Hydrographic Review 15, 151).
//     Courtier's classes are F < 0.25 semidiurnal; 0.25 to 1.5 mixed, mainly semidiurnal; 1.5 to 3 mixed,
//     mainly diurnal; above 3 diurnal. A combination D<q>±W<r> of the lattice has q days, as Doodson's
//     first index gives a tidal constituent's species. So F = Σ_{q=1} A / Σ_{q=2} A over the week's
//     constituents. A street with two crests a day (a commuters' street) is semidiurnal; a street with one
//     crest, in the middle of the day or at night, is diurnal. The thresholds are borrowed from the sea: a
//     convention, not a law of streets;
//   - the weekend's level against the weekday's (their mean counts), and the year's swing (Y1's amplitude).
//
// Checks: tests/check_rhythm.cpp (a constructed street of known constituents: its prediction, high water,
// form and surge recovered).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <tyche/time_series/harmonic.hpp>

#include <urban/rhythm/analysis.hpp>
#include <urban/rhythm/calendar.hpp>
#include <urban/rhythm/record.hpp>

namespace urban::rhythm {

// Civil hours since the week's origin, Monday 1970-01-05 00:00 civil, of the UTC instant t of a record.
inline double week_hours(const Record& r, double t) { return (r.civil(t) - monday_unix) / hour_s; }

// Civil hours since the week's origin at the start of a civil day (days since 1970-01-01).
inline double day_start_hours(long long civil_day) { return (static_cast<double>(civil_day) * day_s - monday_unix) / hour_s; }

struct Tide {
  std::vector<Combination> constituents;
  ts::HarmonicFit fit;                    // t: week_hours at the middle of each hour; fitted and residual cleared
  double first = 0, last = 0;             // UTC of the first and last values fitted
  std::size_t observations = 0;
};

// The tide of a record over [from, to) (UTC, Unix s): least squares at the street constituents that a record of
// that length resolves (a month's too, if asked). At least four weeks of values.
inline Tide tide(const Record& r, double from, double to, bool month = false) {
  Points p;
  Tide out;
  out.first = std::numeric_limits<double>::infinity(), out.last = -out.first;
  for (std::size_t i = 0; i < r.size(); i++) {
    const double t = r.time(i);
    if (t < from || t >= to || !std::isfinite(r.y[i])) continue;
    p.t.push_back(week_hours(r, t) + 0.5 * r.step / hour_s);
    p.y.push_back(r.y[i]);
    out.first = std::min(out.first, t), out.last = std::max(out.last, t);
  }
  if (p.t.size() < 28 * 24) throw std::invalid_argument("urban::rhythm::tide: at least four weeks of values");
  const double hours = p.t.back() - p.t.front();
  out.constituents = street_constituents(hours, month);
  std::vector<ts::Constituent> cs;
  for (const auto& c : out.constituents) cs.push_back({c.name, c.frequency});
  ts::HarmonicOptions o;
  o.t0 = 0, o.coloured = false;
  out.fit = ts::harmonic_fit(p.t, p.y, cs, o);
  out.fit.fitted.clear(), out.fit.residual.clear();
  out.fit.fitted.shrink_to_fit(), out.fit.residual.shrink_to_fit();
  out.observations = p.t.size();
  return out;
}

// The predicted level (people an hour) at civil hours since the week's origin, from the constituents whose
// k satisfies `keep` (all by default); below zero, zero.
template <class Keep>
inline double level(const Tide& t, double hours, Keep keep) {
  double v = t.fit.mean;
  for (std::size_t k = 0; k < t.fit.terms.size(); k++) {
    if (!keep(t.constituents[k])) continue;
    const auto& term = t.fit.terms[k];
    const double arg = 2 * std::numbers::pi * term.frequency * hours;
    v += term.a * std::cos(arg) + term.b * std::sin(arg);
  }
  return std::max(0.0, v);
}
inline double level(const Tide& t, double hours) { return level(t, hours, [](const Combination&) { return true; }); }

inline bool of_the_week(const Combination& c) { return c.k[1] == 0 && c.k[2] == 0; }

// The predicted level through a civil day at per_hour points an hour, from 00:00 to 24:00 inclusive.
inline std::vector<double> day_curve(const Tide& t, long long civil_day, std::size_t per_hour = 4) {
  std::vector<double> out;
  const double h0 = day_start_hours(civil_day);
  for (std::size_t k = 0; k <= 24 * per_hour; k++) out.push_back(level(t, h0 + static_cast<double>(k) / static_cast<double>(per_hour)));
  return out;
}

// The record's counts of the 24 civil hours of a civil day (NaN where missing; the hour a change of offset
// skips is missing, and of an hour it repeats the first is kept).
inline std::array<double, 24> observed_day(const Record& r, long long civil_day) {
  std::array<double, 24> out;
  out.fill(std::nan(""));
  if (r.size() == 0) return out;
  const double guess = static_cast<double>(civil_day) * day_s - 14 * hour_s;
  const long long i0 = static_cast<long long>(std::floor((guess - r.t0) / r.step));
  for (long long i = std::max(0LL, i0); i < std::min(static_cast<long long>(r.size()), i0 + 40); i++) {
    const double t = r.time(static_cast<std::size_t>(i));
    if (r.civil_day(t) != civil_day) continue;
    const auto h = static_cast<std::size_t>(std::floor(r.civil_hour(t) + 1e-9));
    if (h < 24 && !std::isfinite(out[h])) out[h] = r.y[static_cast<std::size_t>(i)];
  }
  return out;
}

struct TideDay {
  double high = 0, high_at = 0;   // people an hour, civil hour of the day
  double low = 0, low_at = 0;
  double mean = 0;
  double range() const { return high - low; }
};

struct TideReadings {
  TideDay weekday, weekend;
  double mean = 0;                 // Z0, people an hour
  double diurnal = 0, semidiurnal = 0;   // Σ A over the week's constituents of species 1 and of species 2
  double form = 0;                 // F = diurnal / semidiurnal
  std::string form_class;          // Courtier's classes
  double annual = 0;               // A of Y1, people an hour
  double variance_explained = 0;
};

inline std::string form_class(double F) {
  if (!std::isfinite(F)) return "diurnal";
  if (F < 0.25) return "semidiurnal";
  if (F < 1.5) return "mixed, mainly semidiurnal";
  if (F <= 3.0) return "mixed, mainly diurnal";
  return "diurnal";
}

inline TideReadings readings(const Tide& t) {
  TideReadings out;
  out.mean = t.fit.mean;
  out.variance_explained = t.fit.variance_explained;
  for (std::size_t k = 0; k < t.fit.terms.size(); k++) {
    const auto& c = t.constituents[k];
    const double A = t.fit.terms[k].amplitude;
    if (c.k == std::array<int, 3>{0, 0, 1}) out.annual = A;
    if (!of_the_week(c)) continue;
    const int q = days_and_weeks(c.k[0]).first;
    if (q == 1) out.diurnal += A;
    if (q == 2) out.semidiurnal += A;
  }
  out.form = out.semidiurnal > 0 ? out.diurnal / out.semidiurnal : std::numeric_limits<double>::infinity();
  out.form_class = form_class(out.form);
  // The mean weekday and weekend day of the mean week, minute by minute.
  auto day = [&](std::initializer_list<int> days) {
    TideDay d;
    d.high = -1, d.low = std::numeric_limits<double>::infinity();
    double sum = 0;
    for (int m = 0; m < 1440; m++) {
      const double tau = m / 60.0;
      double v = 0;
      for (int w : days) v += level(t, 24.0 * w + tau, of_the_week);
      v /= static_cast<double>(days.size());
      sum += v;
      if (v > d.high) d.high = v, d.high_at = tau;
      if (v < d.low) d.low = v, d.low_at = tau;
    }
    d.mean = sum / 1440;
    return d;
  };
  out.weekday = day({0, 1, 2, 3, 4});
  out.weekend = day({5, 6});
  return out;
}

}  // namespace urban::rhythm
