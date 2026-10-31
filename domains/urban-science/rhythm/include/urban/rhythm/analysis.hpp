// The analysis of a street's rhythm (theory/street-oceanography.md §3–4): its spectrum against the standing
// null of integration, its lines read on the calendars' lattice, its harmonic constants, the selectivity of
// its response, the civil share of its clock, and the bicoherence of what the calendars leave.
//
// Spectrum (Definition 3.11). Thomson's multitaper estimate (calx `signal`) of the series (ln(1 + counts),
// hourly), its χ² interval with the estimate's degrees of freedom, and the red-noise background with its
// level (tyche `time_series`): a peak is a local maximum above the level, peaks closer than the bandwidth
// 2W merged.
//
// Lines on the lattice (Propositions 3.7–3.9, Remark 3.8). A peak is explained by the combination of least
// order ≤ d within δ = max(R/T, W) of the declared drivings (calendar.hpp); when the week is declared, its
// lines are all forced whatever their order (a forcing on the week's circle has them all), and the order
// counts only the steps across the month and the year (the cross order). It is
//   - "forced" when that combination is a line of one driving alone (a harmonic of the week, of the year,
//     of the month);
//   - a "combination" when it mixes drivings: disproportionate response when the declared forcing is
//     separable, and "undecided" otherwise (Remark 3.8: a linear response to a non-separable forcing has
//     these lines too);
//   - "own" when nothing of order ≤ d explains it, reported with the chance rate χ(d, T) that a frequency
//     drawn at random would have been explained.
//
// Harmonic constants. tyche's harmonic_fit at the calendar's combinations, keeping, in order of increasing
// order, only those the record resolves from every one already kept (the Rayleigh criterion: T_TIDE's
// choice; Foreman 1977).
//
// Selectivity (Definition 4.2). The fitted lines' variances A²/2 summed by driving (calendar.hpp, driving()),
// shares π_D, S = (|D| Σ π_D² − 1)/(|D| − 1); the forced fraction is their total over the series' variance.
//
// The civil share (Representation 3.12). At each change of the zone's offset, the UTC time of the daily
// fundamental's maximum is fitted over the same weekdays in windows before and after; ψ = −Δt / Δoffset.
//
// The constructed record of the paper's Figure 1 (Street Oceanography, §29.1), and a fourth case: a linear
// response to a non-separable forcing, which the paper's discriminator reads as disproportionate and the
// amended one does not (Remark 3.8).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/signal/bispectrum.hpp>
#include <calx/signal/spectrum.hpp>
#include <tyche/time_series/harmonic.hpp>
#include <tyche/time_series/spectral.hpp>
#include <urban/rhythm/calendar.hpp>
#include <urban/rhythm/record.hpp>

namespace urban::rhythm {

namespace ts = tyche::time_series;

// ---- the spectrum and its lines ---------------------------------------------------------------------------

struct Peak {
  double frequency = 0;          // cycles per hour
  double density = 0, level = 0; // the estimate and the red-noise level there
  double ratio = 0;              // density / level
  Explanation explanation;       // by the declared drivings, order ≤ d
  std::string verdict;           // forced | combination | undecided | own
};

struct SpectrumReport {
  std::vector<double> frequency, density, lower, upper, background, level;   // cycles per hour; one-sided densities
  double rho = 0, variance = 0;  // the red-noise null's ρ and σ²
  double nw = 0, bandwidth = 0;  // W in cycles per hour
  std::size_t tapers = 0;
  double record_hours = 0, resolution = 0;   // T, δ = R/T
  int order = 0;                 // d (the cross order when the week is declared)
  double chance = 0;             // χ(d, T) over [0, f_max], with δ = max(R/T, W)
  std::vector<Peak> peaks;       // by ratio, highest first
};

// One-sided multitaper spectrum of an hourly series with the red-noise null and its lines classified.
// mask: the declared drivings (week, month, year); separable: whether the declared forcing is separable
// across them; f_max: the highest frequency reported (cycles per hour).
inline SpectrumReport spectrum(const std::vector<double>& x, double nw, std::size_t tapers, std::array<bool, 3> mask, bool separable, int order,
                               double f_max = 0.5, double rayleigh = 1.0, std::size_t max_peaks = 60) {
  if (x.size() < 64) throw std::invalid_argument("urban::rhythm::spectrum: at least 64 values");
  const double fs = 1.0;  // per hour
  const auto mt = calx::signal::multitaper(x, fs, nw, tapers, true);
  const auto red = ts::red_noise(x);
  SpectrumReport out;
  out.rho = red.rho, out.variance = red.variance, out.nw = nw, out.tapers = tapers, out.bandwidth = mt.bandwidth;
  out.record_hours = static_cast<double>(x.size());
  out.resolution = rayleigh / out.record_hours;
  out.order = order;
  out.chance = mask[0] ? chance_rate_free_week(mask, order, std::max(out.resolution, mt.bandwidth), f_max)
                       : chance_rate(mask, order, std::max(out.resolution, mt.bandwidth), f_max);
  for (std::size_t j = 0; j < mt.frequency.size(); j++) {
    if (mt.frequency[j] > f_max) break;
    const double dof = j < mt.dof.size() ? mt.dof[j] : 2.0 * static_cast<double>(tapers);
    const auto ci = ts::spectral_interval(mt.density[j], dof);
    out.frequency.push_back(mt.frequency[j]);
    out.density.push_back(mt.density[j]);
    out.lower.push_back(ci.lower);
    out.upper.push_back(ci.upper);
    out.background.push_back(red.density(mt.frequency[j], fs));
    out.level.push_back(ts::red_noise_level(red, mt.frequency[j], fs, dof));
  }
  // Local maxima above the level, merged within 2W.
  const double df = out.frequency.size() > 1 ? out.frequency[1] - out.frequency[0] : 1;
  const std::size_t w = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(out.bandwidth / df)));
  std::vector<Peak> found;
  for (std::size_t j = 1; j + 1 < out.frequency.size(); j++) {
    if (out.density[j] <= out.level[j]) continue;
    bool top = true;
    for (std::size_t i = (j > w ? j - w : 0); i <= std::min(out.frequency.size() - 1, j + w) && top; i++) top = out.density[i] <= out.density[j];
    if (!top) continue;
    Peak p;
    p.frequency = out.frequency[j], p.density = out.density[j], p.level = out.level[j], p.ratio = p.density / p.level;
    if (!found.empty() && p.frequency - found.back().frequency < 2 * out.bandwidth) {
      if (p.density > found.back().density) found.back() = p;
      continue;
    }
    found.push_back(p);
  }
  const double delta = std::max(out.resolution, out.bandwidth);
  for (auto& p : found) {
    p.explanation = mask[0] ? explain_free_week(p.frequency, mask, order, delta) : explain(p.frequency, mask, order, delta);
    if (!p.explanation.found) p.verdict = "own";
    else {
      const auto& k = p.explanation.combination.k;
      const int parts = (k[0] != 0) + (k[1] != 0) + (k[2] != 0);
      p.verdict = parts <= 1 ? "forced" : (separable ? "combination" : "undecided");
    }
  }
  std::sort(found.begin(), found.end(), [](const Peak& a, const Peak& b) { return a.ratio > b.ratio; });
  if (found.size() > max_peaks) found.resize(max_peaks);
  out.peaks = std::move(found);
  return out;
}

// ---- harmonic constants -------------------------------------------------------------------------------------

// The calendar's combinations of the drivings in the mask up to the order and days given, ordered by
// increasing order (then frequency), keeping those resolved from every one kept before.
inline std::vector<Combination> resolvable_constituents(std::array<bool, 3> mask, int order, int max_day, double record_hours, double rayleigh = 1.0,
                                                        double f_max = 0.5) {
  auto all = lattice(mask, order, max_day);
  std::stable_sort(all.begin(), all.end(), [](const Combination& a, const Combination& b) { return a.order < b.order; });
  std::vector<Combination> kept;
  for (const auto& c : all) {
    if (c.frequency > f_max || c.frequency * record_hours < rayleigh) continue;   // separable from the mean
    bool ok = true;
    for (const auto& k : kept) ok = ok && ts::resolvable(c.frequency, k.frequency, record_hours, rayleigh);
    if (ok) kept.push_back(c);
  }
  std::sort(kept.begin(), kept.end(), [](const Combination& a, const Combination& b) { return a.frequency < b.frequency; });
  return kept;
}

// The street constituents: one list for every street, as tidal analysis has its standard constituents, so
// that streets are compared on the same lines. Every line of the week to six cycles a day (k_W = 1…45: the
// day's harmonics D1–D6 and their weekly sidebands to ±W3), the year's first three harmonics, the annual
// sidebands of D1, D2 and W1, and, when asked, the synodic month and its sidebands of D1. Kept in order of
// increasing order when the record resolves them from the lines kept before (Rayleigh, R = 1), and from
// the mean.
inline std::vector<Combination> street_constituents(double record_hours, bool month = false, double rayleigh = 1.0) {
  std::vector<std::array<int, 3>> ks;
  for (int w = 1; w <= 45; w++) ks.push_back({w, 0, 0});
  for (int y = 1; y <= 3; y++) ks.push_back({0, 0, y});
  for (int w : {7, 14, 1})
    for (int y : {-1, 1}) ks.push_back({w, 0, y});
  if (month) ks.push_back({0, 1, 0}), ks.push_back({7, -1, 0}), ks.push_back({7, 1, 0});
  std::vector<Combination> all;
  for (const auto& k : ks) all.push_back(combination(k));
  std::stable_sort(all.begin(), all.end(), [](const Combination& a, const Combination& b) { return a.order < b.order; });
  std::vector<Combination> kept;
  for (const auto& c : all) {
    if (c.frequency * record_hours < rayleigh) continue;
    bool ok = true;
    for (const auto& k : kept) ok = ok && ts::resolvable(c.frequency, k.frequency, record_hours, rayleigh);
    if (ok) kept.push_back(c);
  }
  std::sort(kept.begin(), kept.end(), [](const Combination& a, const Combination& b) { return a.frequency < b.frequency; });
  return kept;
}

struct Harmonics {
  std::vector<Combination> constituents;
  ts::HarmonicFit fit;
  double t_ref = 0;              // UTC of the phases' origin, Unix s
};

// Points (t in hours since t_ref, y) fitted at the given combinations.
inline Harmonics harmonics(const Points& p, const std::vector<Combination>& constituents, double t_ref, bool coloured = true) {
  Harmonics h;
  h.constituents = constituents;
  h.t_ref = t_ref;
  std::vector<ts::Constituent> cs;
  for (const auto& c : constituents) cs.push_back({c.name, c.frequency});
  ts::HarmonicOptions o;
  o.t0 = 0;
  o.coloured = coloured;
  h.fit = ts::harmonic_fit(p.t, p.y, cs, o);
  return h;
}

struct Selectivity {
  std::vector<std::string> drivings;
  std::vector<double> share;     // π_D
  double selectivity = 0;        // S
  double forced = 0;             // Σ A²/2 / var(y)
};

inline Selectivity selectivity(const Harmonics& h, double variance) {
  std::map<std::string, double> by;
  for (const auto& c : h.constituents) by[driving(c)] = 0;
  double total = 0;
  for (std::size_t i = 0; i < h.fit.terms.size() && i < h.constituents.size(); i++) {
    const double v = 0.5 * h.fit.terms[i].amplitude * h.fit.terms[i].amplitude;
    by[driving(h.constituents[i])] += v;
    total += v;
  }
  Selectivity s;
  double ipr = 0;
  for (const auto& [d, v] : by) {
    s.drivings.push_back(d);
    s.share.push_back(total > 0 ? v / total : 0);
    ipr += total > 0 ? (v / total) * (v / total) : 0;
  }
  const double n = static_cast<double>(by.size());
  s.selectivity = n > 1 ? (n * ipr - 1) / (n - 1) : 1;
  s.forced = variance > 0 ? total / variance : 0;
  return s;
}

// ---- the civil share of the clock ----------------------------------------------------------------------------

struct ClockShift {
  double at = 0;                 // UTC of the change, Unix s
  double change = 0;             // new offset − old, hours
  double before = 0, after = 0;  // UTC hour of the day of the daily fundamental's maximum
  double shift = 0, se = 0;      // after − before, hours, wrapped to (−12, 12]
  double civil_share = 0;        // ψ = −shift / change
  std::size_t days_before = 0, days_after = 0;
};

// For every change of offset inside the record: weekdays (Monday to Friday, civil) in `window` days before
// and after, the transition's day and the days in `skip` (holidays: civil days since 1970) left out.
inline std::vector<ClockShift> clock_shifts(const Record& r, int window = 21, const std::vector<long long>& skip = {}) {
  std::vector<ClockShift> out;
  for (std::size_t k = 1; k < r.offsets.size(); k++) {
    const double at = r.offsets[k].first, change = r.offsets[k].second - r.offsets[k - 1].second;
    if (change == 0 || at < r.t0 || at > r.time(r.size() - 1)) continue;
    auto fit_side = [&](double from, double to, std::size_t& days) -> std::pair<double, double> {
      Points p;
      std::vector<long long> seen;
      for (std::size_t i = 0; i < r.size(); i++) {
        const double t = r.time(i);
        if (t < from || t >= to || !std::isfinite(r.y[i])) continue;
        const int wd = r.weekday(t);
        const long long day = r.civil_day(t);
        if (wd > 4 || std::find(skip.begin(), skip.end(), day) != skip.end()) continue;
        p.t.push_back((t - at) / 3600.0);
        p.y.push_back(std::log1p(std::max(0.0, r.y[i])));
        if (seen.empty() || seen.back() != day) seen.push_back(day);
      }
      days = seen.size();
      if (p.t.size() < 48) return {std::nan(""), std::nan("")};
      const std::vector<ts::Constituent> cs{{"D1", 1.0 / 24}, {"D2", 1.0 / 12}};
      ts::HarmonicOptions o;
      o.coloured = false;
      const auto fit = ts::harmonic_fit(p.t, p.y, cs, o);
      // A cos(2π f t − φ) peaks at t = φ / (2π f): hours after the change, then as an hour of the UTC day.
      const double peak = fit.terms[0].phase / (2 * std::numbers::pi) * 24.0;
      const double hour_of_day = std::fmod(std::fmod(at / 3600.0 + peak, 24.0) + 24.0, 24.0);
      return {hour_of_day, fit.terms[0].phase_se / (2 * std::numbers::pi) * 24.0};
    };
    ClockShift c;
    c.at = at, c.change = change;
    const double day = 86400.0;
    const auto [b, sb] = fit_side(at - window * day, at - day, c.days_before);
    const auto [a, sa] = fit_side(at + day, at + (window + 1) * day, c.days_after);
    if (!std::isfinite(a) || !std::isfinite(b)) continue;
    c.before = b, c.after = a;
    double d = a - b;
    while (d > 12) d -= 24;
    while (d <= -12) d += 24;
    c.shift = d, c.se = std::sqrt(sa * sa + sb * sb);
    c.civil_share = -d / change;
    out.push_back(c);
  }
  return out;
}

// ---- the bicoherence of what the calendars leave ---------------------------------------------------------------

struct Bicoherence {
  std::vector<double> frequency;       // cycles per hour, bins
  std::size_t bins = 0, segments = 0;
  std::vector<double> b2;              // row-major bins × bins (principal domain)
  double level = 0;                    // tyche's level at α = 0.05
  struct Cell { std::size_t j, k; double b2; };
  std::vector<Cell> above;             // cells above the level, highest first (at most 40)
  double share_above = 0;              // of the principal domain's cells (expected ≈ α under the null)
};

inline Bicoherence bicoherence(const std::vector<double>& residual, std::size_t segment, std::size_t max_bin, double overlap = 0.5) {
  const auto b = calx::signal::bispectrum(residual, 1.0, segment, overlap, true, max_bin);
  Bicoherence out;
  out.frequency = b.frequency, out.bins = b.bins, out.segments = b.segments, out.b2 = b.bicoherence;
  out.level = ts::bicoherence_level(b.segments, 0.05);
  std::size_t cells = 0;
  for (std::size_t j = 1; j < b.bins; j++)
    for (std::size_t k = 1; k <= j && j + k < b.bins; k++) {
      cells++;
      const double v = b.b2(j, k);
      if (v > out.level) out.above.push_back({j, k, v});
    }
  out.share_above = cells ? static_cast<double>(out.above.size()) / static_cast<double>(cells) : 0;
  std::sort(out.above.begin(), out.above.end(), [](const Bicoherence::Cell& x, const Bicoherence::Cell& y) { return x.b2 > y.b2; });
  if (out.above.size() > 40) out.above.resize(40);
  return out;
}

// ---- the paper's constructed record, and the case it misreads ---------------------------------------------------

// Generators of any frequencies (no days and weeks): the combination of least |k|₁ ≤ d within δ.
struct GenericExplanation { bool found = false; std::vector<int> k; double frequency = 0, distance = 0; int order = 0; };
inline GenericExplanation explain_generic(double f, const std::vector<double>& g, int d, double delta) {
  GenericExplanation best;
  std::vector<int> k(g.size(), -d);
  while (true) {
    int order = 0;
    double fk = 0;
    for (std::size_t i = 0; i < g.size(); i++) order += std::abs(k[i]), fk += k[i] * g[i];
    if (order > 0 && order <= d && fk > 0) {
      const double dist = std::fabs(fk - f);
      if (dist <= delta && (!best.found || order < best.order || (order == best.order && dist < best.distance)))
        best = {true, k, fk, dist, order};
    }
    std::size_t i = 0;
    while (i < k.size() && ++k[i] > d) k[i++] = -d;
    if (i == k.size()) break;
  }
  return best;
}

inline double chance_rate_generic(const std::vector<double>& g, int d, double delta, double f_max) {
  std::vector<std::pair<double, double>> iv;
  std::vector<int> k(g.size(), -d);
  while (true) {
    int order = 0;
    double fk = 0;
    for (std::size_t i = 0; i < g.size(); i++) order += std::abs(k[i]), fk += k[i] * g[i];
    if (order > 0 && order <= d && fk > 0 && fk - delta < f_max) iv.emplace_back(std::max(0.0, fk - delta), std::min(f_max, fk + delta));
    std::size_t i = 0;
    while (i < k.size() && ++k[i] > d) k[i++] = -d;
    if (i == k.size()) break;
  }
  std::sort(iv.begin(), iv.end());
  double covered = 0, lo = 0, hi = 0;
  bool open = false;
  for (const auto& [a, b] : iv) {
    if (!open || a > hi) { if (open) covered += hi - lo; lo = a, hi = b, open = true; }
    else hi = std::max(hi, b);
  }
  if (open) covered += hi - lo;
  return covered / f_max;
}

struct ConstructedCase {
  std::string name, about;
  std::vector<double> frequency, power;   // cycles per day; periodogram normalised to the strongest driving's line
  std::map<std::string, double> at;       // normalised power at 2f₁, f₁ + f₂, f₁ − f₂, 1/11.3 d
  std::string paper_verdict;              // Discriminator 29.1 as written: power at sums, differences, multiples
  std::string amended_verdict;            // after removing the linear response to the declared forcing
  std::vector<std::string> residual_lines;   // lines of the residual above its level, with their explanation
};

// The record of Figure 1 (§29.1): drivings of 7 and 29.5 days, 24 samples a day for 512 days, the same
// white noise in every case; and case 4, a linear response to the non-separable forcing (1 + m cos θ₂) cos θ₁.
inline std::vector<ConstructedCase> constructed_cases(std::uint64_t seed = 1, double noise = 0.5) {
  const double f1 = 1.0 / 7.0, f2 = 1.0 / 29.5;   // cycles per day
  const std::size_t n = 512 * 24;
  calx::random::Stream rng(seed, 0x51, 0x29, 0);
  std::vector<double> e(n);
  for (auto& v : e) v = noise * rng.normal();
  const double two_pi = 2 * std::numbers::pi;
  std::vector<ConstructedCase> out;
  for (int c = 1; c <= 4; c++) {
    std::vector<double> y(n), t(n);
    for (std::size_t i = 0; i < n; i++) {
      t[i] = static_cast<double>(i) / 24.0;
      const double x1 = std::cos(two_pi * f1 * t[i]), x2 = 0.6 * std::cos(two_pi * f2 * t[i]);
      const double x = x1 + x2;
      double v = x;
      if (c == 2) v = x + 0.5 * x * x;
      if (c == 3) v = x + 0.8 * std::cos(two_pi * t[i] / 11.3 + 1.0);
      if (c == 4) v = (1 + 0.6 * std::cos(two_pi * f2 * t[i])) * x1 + x2;
      y[i] = v + e[i];
    }
    ConstructedCase k;
    k.name = c == 1 ? "only the drivings" : c == 2 ? "sums, differences and multiples" : c == 3 ? "a component of its own" : "a linear response to a non-separable forcing";
    k.about = c == 1 ? "y = x₁ + x₂ + noise" : c == 2 ? "y = x + ½x² + noise (a quadratic response)" : c == 3 ? "y = x + 0.8 cos(2πt/11.3 d + 1) + noise"
                     : "y = (1 + 0.6 cos 2πf₂t) cos 2πf₁t + x₂ + noise: linear in a forcing whose weekly amplitude follows the month";
    const auto pg = calx::signal::periodogram(y, 24.0);   // cycles per day
    double strongest = 0;
    const auto power_at = [&](double f) {
      std::size_t j = static_cast<std::size_t>(std::lround(f / (pg.frequency[1] - pg.frequency[0])));
      double best = 0;
      for (std::size_t i = (j > 2 ? j - 2 : 0); i <= std::min(pg.frequency.size() - 1, j + 2); i++) best = std::max(best, pg.density[i]);
      return best;
    };
    strongest = power_at(f1);
    for (std::size_t j = 0; j < pg.frequency.size() && pg.frequency[j] <= 0.35; j++) {
      k.frequency.push_back(pg.frequency[j]);
      k.power.push_back(pg.density[j] / strongest);
    }
    k.at["2f1"] = power_at(2 * f1) / strongest;
    k.at["f1+f2"] = power_at(f1 + f2) / strongest;
    k.at["f1-f2"] = power_at(f1 - f2) / strongest;
    k.at["11.3d"] = power_at(1 / 11.3) / strongest;
    // The paper: sums, differences or multiples well above the noise → disproportionate.
    double floor = 0;
    std::size_t cnt = 0;
    for (std::size_t j = 0; j < pg.frequency.size(); j++)
      if (pg.frequency[j] > 0.25 && pg.frequency[j] < 1.0) floor += pg.density[j] / strongest, cnt++;
    floor /= static_cast<double>(std::max<std::size_t>(cnt, 1));
    const bool combos = k.at["2f1"] > 100 * floor || k.at["f1+f2"] > 100 * floor || k.at["f1-f2"] > 100 * floor;
    k.paper_verdict = combos ? "disproportionate response" : (k.at["11.3d"] > 100 * floor ? "a process of its own" : "proportional response");
    // Amended: remove the linear response to the declared forcing (its Fourier support), then read the residual.
    std::vector<ts::Constituent> declared{{"f1", f1}, {"f2", f2}};
    if (c == 4) declared.push_back({"f1+f2", f1 + f2}), declared.push_back({"f1-f2", f1 - f2});
    ts::HarmonicOptions o;
    o.coloured = false;
    const auto fit = ts::harmonic_fit(t, y, declared, o);
    const auto res = calx::signal::periodogram(fit.residual, 24.0);
    // Residual lines: bins above 30 times the residual's median density (white noise: P(ratio > 30) ≈ e^{−30 ln 2}).
    std::vector<double> sorted(res.density.begin() + 1, res.density.end());
    std::nth_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(sorted.size() / 2), sorted.end());
    const double median = sorted[sorted.size() / 2];
    bool nonlinear = false, own = false;
    const double delta = 1.0 / 512.0;
    for (std::size_t j = 1; j + 1 < res.frequency.size() && res.frequency[j] <= 1.0; j++) {
      if (res.density[j] < 30 * median || res.density[j] < res.density[j - 1] || res.density[j] < res.density[j + 1]) continue;
      const auto ex = explain_generic(res.frequency[j], {f1, f2}, 3, 2 * delta);
      std::string line = "f = " + std::to_string(res.frequency[j]).substr(0, 6) + " d⁻¹";
      if (ex.found) {
        line += ": k = (" + std::to_string(ex.k[0]) + ", " + std::to_string(ex.k[1]) + "), order " + std::to_string(ex.order);
        nonlinear = true;
      } else {
        line += ": no combination of order ≤ 3";
        own = true;
      }
      k.residual_lines.push_back(line);
    }
    k.amended_verdict = own ? "a process of its own (no combination of order ≤ 3)" : nonlinear ? "disproportionate response" : "proportional response";
    out.push_back(std::move(k));
  }
  return out;
}

}  // namespace urban::rhythm
