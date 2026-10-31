// The spacing of clusters of a use (theory/street-oceanography.md, Representation 5.3; Huang, Street
// Oceanography, §11.3, fourth mechanism): the spacing a reinforcement with a longer-reaching suppression
// predicts, beside the spacing a street shows.
//
// Predicted. About a uniform state with kinetics' Jacobian J and diffusivities d_u (the reinforcing use) and
// d_v (the suppressing influence), the fastest-growing squared wavenumber k² (physis `chemistry`,
// turing_two_species) sets a wavelength 2π/k in the units of √(d/rate) — metres when the diffusivities are in
// m² per unit time and the rates per that unit. The condition that can fail: d_v/d_u above the critical ratio.
//
// Observed. The indicator of the use at the holdings' positions s_h, unevenly spaced, has a Lomb–Scargle
// periodogram (tyche, floating mean) over wavelengths from twice the median spacing of holdings to the
// street's length; its highest peak's wavelength, with Baluev's false-alarm probability.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include <physis/chemistry/reaction_diffusion.hpp>
#include <tyche/time_series/lomb_scargle.hpp>

namespace urban::affinity {

struct PredictedSpacing {
  physis::chemistry::TwoSpeciesTuring turing;
  double wavelength = 0;          // 2π/k at the fastest mode; NaN if not Turing-unstable
};

inline PredictedSpacing predicted_spacing(double j11, double j12, double j21, double j22, double du, double dv) {
  PredictedSpacing p;
  p.turing = physis::chemistry::turing_two_species(j11, j12, j21, j22, du, dv);
  p.wavelength = p.turing.unstable && p.turing.k2_fastest > 0 ? 2 * std::numbers::pi / std::sqrt(p.turing.k2_fastest) : std::nan("");
  return p;
}

struct ObservedSpacing {
  std::vector<double> wavelength, power;   // m; Lomb–Scargle power
  double best = 0;                         // the peak's wavelength, m
  double best_power = 0;
  double false_alarm = 1;                  // Baluev 2008
  std::size_t holdings = 0, of_use = 0;
};

inline ObservedSpacing observed_spacing(const std::vector<double>& s, const std::vector<std::size_t>& use, std::size_t which, double length,
                                        std::size_t grid = 400) {
  ObservedSpacing o;
  o.holdings = s.size();
  std::vector<double> y;
  for (auto u : use) y.push_back(u == which ? 1.0 : 0.0), o.of_use += u == which;
  if (s.size() < 8 || o.of_use < 2 || o.of_use == s.size()) return o;
  std::vector<double> sorted = s;
  std::sort(sorted.begin(), sorted.end());
  std::vector<double> gaps;
  for (std::size_t i = 1; i < sorted.size(); i++) gaps.push_back(sorted[i] - sorted[i - 1]);
  std::nth_element(gaps.begin(), gaps.begin() + static_cast<std::ptrdiff_t>(gaps.size() / 2), gaps.end());
  const double median = std::max(1.0, gaps[gaps.size() / 2]);
  const double f_lo = 1.0 / length, f_hi = 1.0 / (2 * median);
  if (!(f_hi > f_lo)) return o;
  std::vector<double> f(grid);
  for (std::size_t i = 0; i < grid; i++) f[i] = f_lo + (f_hi - f_lo) * static_cast<double>(i) / static_cast<double>(grid - 1);
  const auto ls = tyche::time_series::lomb_scargle(s, y, f, true);
  std::size_t best = 0;
  for (std::size_t i = 0; i < grid; i++) {
    o.wavelength.push_back(1.0 / f[i]);
    o.power.push_back(ls.power[i]);
    if (ls.power[i] > ls.power[best]) best = i;
  }
  o.best = 1.0 / f[best];
  o.best_power = ls.power[best];
  o.false_alarm = tyche::time_series::false_alarm_baluev(o.best_power, s, y, f_hi);
  return o;
}

}  // namespace urban::affinity
