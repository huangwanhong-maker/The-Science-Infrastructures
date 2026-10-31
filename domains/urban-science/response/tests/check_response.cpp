// Checks of the response module (theory/street-oceanography.md §4):
//
// - A constructed street whose log counts follow the hour of the week, a season and two lags of rain
//   (−0.20 and −0.10 per mm/h) and a holiday (+0.30): respond() recovers each slope within three of its
//   clustered standard errors, and the rain's susceptibility −0.30 as their sum, with the covariance of the lags.
// - A forcing of no effect is found consistent with zero.
// - The paired comparison: days of an alteration with totals 30% higher than their neighbours of the same
//   weekday give the effect ln 1.3 and a small randomisation p; without the alteration, p is not small.

#include <urban/response/susceptibility.hpp>
#include <urban/rhythm/record.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace urban::response;
using urban::rhythm::parse_utc;

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(5, 5, 5, 5);
  const double two_pi = 2 * std::numbers::pi;

  Record r;
  r.t0 = parse_utc("2023-01-02T00:00:00Z");
  r.offsets = {{r.t0 - 86400, 10.0}};
  const std::size_t n = 2 * 8760;
  std::vector<double> rain(n, 0), holiday(n, 0), placebo(n, 0);
  // Rain falls in spells; holidays are whole civil days; a placebo of the same shape as rain does nothing.
  for (std::size_t i = 0; i < n; i++) {
    if (rng.uniform() < 0.03) for (std::size_t k = i; k < std::min(n, i + 5); k++) rain[k] = rng.uniform(0, 4);
    placebo[i] = rng.uniform() < 0.1 ? rng.uniform(0, 3) : 0;
  }
  std::vector<long long> holidays;
  for (int k = 0; k < 24; k++) holidays.push_back(r.civil_day(r.t0) + 15 + 29 * k);
  for (std::size_t i = 0; i < n; i++) {
    const double t = r.time(i);
    const long long day = r.civil_day(t);
    if (std::find(holidays.begin(), holidays.end(), day) != holidays.end()) holiday[i] = 1;
    const double hour = r.civil(t) / 3600.0;
    const double base = 4 + 1.2 * std::cos(two_pi * (hour - 13) / 24) + 0.3 * std::cos(two_pi * (hour - 60) / 168) +
                        0.15 * std::cos(two_pi * urban::rhythm::year_angle(t));
    const double v = base - 0.20 * rain[i] - 0.10 * (i ? rain[i - 1] : 0) + 0.30 * holiday[i] + 0.15 * rng.normal();
    r.y.push_back(std::expm1(v));
  }
  const auto fit = respond(r, {{"rain", "mm/h", rain, {0, 1}}, {"holiday", "", holiday, {0}}, {"placebo", "", placebo, {0}}});
  const auto& R = fit.forcings[0];
  s.check("response: the rain's slope at lag 0, −0.20, within 3 standard errors", (R.beta[0] + 0.20) / (3 * R.se[0]), 1.0);
  s.check("response: the rain's slope at lag 1, −0.10, within 3 standard errors", (R.beta[1] + 0.10) / (3 * R.se[1]), 1.0);
  s.check("response: the rain's susceptibility −0.30, within 3 standard errors", (R.total + 0.30) / (3 * R.total_se), 1.0);
  s.check("response: the holiday's +0.30, within 3 standard errors", (fit.forcings[1].total - 0.30) / (3 * fit.forcings[1].total_se), 1.0);
  s.check("response: a forcing of no effect is within 3.5 standard errors of zero", fit.forcings[2].z / 3.5, 1.0);
  s.require("response: standard errors clustered by day (730 days)", fit.days >= 729 && fit.days <= 731, std::to_string(fit.days));

  // The paired comparison.
  {
    Record q;
    q.t0 = parse_utc("2023-01-02T00:00:00Z");
    q.offsets = {{q.t0 - 86400, 0.0}};
    std::vector<long long> a_days;
    for (int k = 0; k < 40; k++) a_days.push_back(q.civil_day(q.t0) + 20 + 9 * k);
    for (std::size_t i = 0; i < 400 * 24; i++) {
      const double t = q.time(i);
      const bool a = std::find(a_days.begin(), a_days.end(), q.civil_day(t)) != a_days.end();
      q.y.push_back(100 * (1 + 0.5 * std::cos(two_pi * static_cast<double>(i % 24) / 24)) * (a ? 1.3 : 1.0) * std::exp(0.05 * rng.normal()));
    }
    const auto e = paired_effect(q, a_days);
    s.check("paired: the effect of days 30% higher is ln 1.3", e.effect - std::log(1.3), 0.03);
    s.require("paired: its randomisation p is small", e.p < 0.001 && e.pairs >= 35, std::to_string(e.p) + ", " + std::to_string(e.pairs) + " pairs");
    std::vector<long long> none;
    for (int k = 0; k < 40; k++) none.push_back(q.civil_day(q.t0) + 24 + 9 * k);
    std::vector<long long> exclude = a_days;
    const auto f = paired_effect(q, none, exclude);
    s.require("paired: days without an alteration show no effect (p > 0.01)", f.p > 0.01, std::to_string(f.p));
  }
  return s.done();
}
