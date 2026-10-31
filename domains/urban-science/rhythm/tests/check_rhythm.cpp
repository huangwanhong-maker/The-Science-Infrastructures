// Checks of the rhythm module (theory/street-oceanography.md §3–4):
//
// - The calendar's combinations: names and orders read in days and weeks (W8 = D1+W1, order 2), the day's
//   frequency exactly 1/24 h⁻¹; the chance rate against a brute-force count of non-overlapping intervals and
//   its growth with the order; an annual sideband explained.
// - Dates: J2000.0 at its Unix second, the round trip of days and civil dates over four centuries, the weekday
//   of a Monday.
// - Spans of a record with gaps, short gaps filled and long ones splitting.
// - Harmonic constants of a constructed street (log counts with a day, a week and a year in civil time):
//   amplitudes and the day's phase recovered, and the selectivity of the shares in closed form.
// - The spectrum: the day forced, a component of 11.3 days found as "own".
// - The civil share: activity keeping the civil clock gives ψ ≈ 1 at each change of offset; keeping UTC, ψ ≈ 0.
// - The paper's Figure 1: the three cases read as the paper reads them, and the fourth (a linear response to
//   a non-separable forcing) read as disproportionate by the paper and as proportional by the amendment.
// - The tide: a constructed street on the lattice in civil time (Melbourne's offsets) predicted exactly; its
//   high water on the mean weekday at the brute force's minute; its form 2.5, mixed, mainly diurnal; a day's
//   dip of 300 an hour read back as its surge; the hour a change of offset skips is missing from its day.
// - Bicoherence: quadratic phase coupling of random phases found at its pair of bins; independent noise exceeds
//   the level in about α of the cells.

#include <urban/rhythm/analysis.hpp>
#include <urban/rhythm/calendar.hpp>
#include <urban/rhythm/record.hpp>
#include <urban/rhythm/tide.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace urban::rhythm;

int main() {
  calx::verification::Suite s;
  const double two_pi = 2 * std::numbers::pi;

  // ---- the calendar ----
  {
    s.require("calendar: W8 is D1+W1, of order 2", name_of({8, 0, 0}) == "D1+W1" && order_of({8, 0, 0}) == 2, name_of({8, 0, 0}));
    s.require("calendar: W6 is D1−W1 and W4 is D1−W3", name_of({6, 0, 0}) == "D1−W1" && name_of({4, 0, 0}) == "D1−W3");
    s.require("calendar: (7, 0, −1) is D1−Y1", name_of({7, 0, -1}) == "D1−Y1");
    s.check("calendar: the day's frequency is 1/24 h⁻¹", frequency({7, 0, 0}) - 1.0 / 24, 1e-17);
    const double delta = 1e-5, fmax = 0.1;
    std::size_t lines = 0;
    for (const auto& c : lattice({true, false, false}, 2, 2))
      if (c.frequency + delta <= fmax) lines++;
    s.check("chance rate: the week alone, order 2: 2δ times its lines", chance_rate({true, false, false}, 2, delta, fmax) - 2 * delta * lines / fmax, 1e-12);
    const double r2 = chance_rate({true, true, true}, 2, 1.0 / 17520, 0.5), r4 = chance_rate({true, true, true}, 4, 1.0 / 17520, 0.5);
    s.require("chance rate: grows with the order (Proposition 3.9)", r4 > r2 && r2 > 0, std::to_string(r2) + " → " + std::to_string(r4));
    const auto e = explain(1.0 / 24 + 1.0 / (tropical_days * 24), {true, false, true}, 3, 1e-6);
    s.require("explain: the day's upper annual sideband is D1+Y1", e.found && e.combination.name == "D1+Y1", e.found ? e.combination.name : "none");
  }

  // ---- dates and records ----
  {
    s.check("dates: J2000.0 is Unix 946728000", parse_utc("2000-01-01T12:00:00Z") - 946728000, 0);
    bool round = true;
    for (long long z = -50000; z < 100000; z += 37) {
      const auto c = civil_from_days(z);
      round = round && days_from_civil(c.year, c.month, c.day) == z;
    }
    s.require("dates: days ↔ civil dates round trip over four centuries", round);
    Record r;
    r.t0 = parse_utc("1970-01-05T00:00:00Z");
    r.y.assign(48, 1.0);
    s.require("dates: 1970-01-05 is a Monday", r.weekday(r.t0) == 0 && r.weekday(r.t0 + 86400) == 1);
    for (int i = 10; i < 12; i++) r.y[static_cast<std::size_t>(i)] = std::nan("");
    for (int i = 30; i < 40; i++) r.y[static_cast<std::size_t>(i)] = std::nan("");
    const auto sp = spans(r, 5, 3);
    s.require("spans: a gap of 2 is filled and a gap of 10 splits", sp.size() == 2 && sp[0].begin == 0 && sp[0].end == 30 && sp[0].filled == 2 && sp[1].begin == 40,
              std::to_string(sp.size()));
    const auto x = series(r, sp[0]);
    s.check("spans: the filled values are interpolated", x[10] - 1.0, 1e-15);
  }

  // ---- a constructed street: harmonic constants, selectivity, spectrum ----
  calx::random::Stream rng(7, 7, 7, 7);
  Record street;
  street.id = "constructed";
  street.t0 = parse_utc("2024-01-01T00:00:00Z");   // a Monday
  street.offsets = {{street.t0 - 86400, 10.0}};
  const std::size_t n = 2 * 8760;
  for (std::size_t i = 0; i < n; i++) {
    const double t = street.time(i), civ = street.civil(t) / 3600.0;
    const double v = 3 + 1.0 * std::cos(two_pi * (civ - 13) / 24) + 0.4 * std::cos(two_pi * (civ - 50) / 168) +
                     0.2 * std::cos(two_pi * year_angle(t)) + 0.3 * std::cos(two_pi * (t - street.t0) / (11.3 * 86400)) + 0.1 * rng.normal();
    street.y.push_back(std::expm1(v));
  }
  {
    const auto p = points(street, 0, n, street.t0, true);
    const auto cs = resolvable_constituents({true, false, true}, 2, 2, static_cast<double>(n));
    const auto h = harmonics(p, cs, street.t0, false);
    auto term = [&](const std::string& name) {
      for (std::size_t i = 0; i < h.constituents.size(); i++)
        if (h.constituents[i].name == name) return h.fit.terms[i];
      return tyche::time_series::HarmonicTerm{};
    };
    s.check("harmonics: the day's amplitude 1.0", term("D1").amplitude - 1.0, 0.01);
    s.check("harmonics: the week's amplitude 0.4", term("W1").amplitude - 0.4, 0.01);
    s.check("harmonics: the year's amplitude 0.2", term("Y1").amplitude - 0.2, 0.02);
    // cos(2π(t + 10 − 13)/24) = cos(2πt/24 − π/4): the phase π/4 with t from the start (civil 10:00).
    s.check("harmonics: the day's phase π/4 (13:00 civil, from a start at 10:00 civil)", term("D1").phase - std::numbers::pi / 4, 0.01, "rad");
    // Selectivity: the shares of variance 0.5, 0.08, 0.02 of day, week and year (the 11.3-day component left out
    // of the constituents), S = (3 Σπ² − 1)/2.
    double var = 0, mean = 0;
    for (double v : p.y) mean += v;
    mean /= static_cast<double>(p.y.size());
    for (double v : p.y) var += (v - mean) * (v - mean);
    var /= static_cast<double>(p.y.size());
    const auto sel = selectivity(h, var);
    const double pd = 0.5 / 0.6, pw = 0.08 / 0.6, py = 0.02 / 0.6;
    s.check("selectivity: S of shares 0.83, 0.13, 0.03 over three drivings", sel.selectivity - (3 * (pd * pd + pw * pw + py * py) - 1) / 2, 0.02);
  }
  {
    std::vector<double> x;
    for (double v : street.y) x.push_back(std::log1p(v));
    const auto rep = spectrum(x, 4, 7, {true, false, true}, false, 3, 0.25);
    bool day_forced = false, own = false;
    for (const auto& pk : rep.peaks) {
      if (pk.explanation.found && pk.explanation.combination.name == "D1" && pk.verdict == "forced") day_forced = true;
      if (std::fabs(pk.frequency - 1.0 / (11.3 * 24)) < 2 * rep.bandwidth && pk.verdict == "own") own = true;
    }
    s.require("spectrum: the day's line is forced", day_forced);
    const auto d4 = explain_free_week(4.0 / 24 + 1e-5, {true, false, true}, 0, 1e-4);
    s.require("spectrum: with the week free, the day's fourth harmonic is a line of the week (cross order 0)", d4.found && d4.combination.name == "D4" && cross_order(d4.combination) == 0,
              d4.found ? d4.combination.name : "none");
    // A band ending between the 41st and 42nd lines of the week (41/168 = 0.2440, 42/168 = 0.25): 41 whole intervals.
    const double week_only = chance_rate_free_week({true, false, false}, 0, 1e-4, 0.247);
    s.check("chance rate: the week's lines alone cover 41 intervals of 2δ in a band of 0.247 h⁻¹", week_only - 41 * 2e-4 / 0.247, 1e-12);
    s.require("spectrum: the component of 11.3 days is the street's own (order ≤ 3)", own);
    s.require("spectrum: the red-noise null has ρ in (0, 1)", rep.rho > 0 && rep.rho < 1, std::to_string(rep.rho));
  }

  // ---- the civil share of the clock ----
  {
    for (int civil = 0; civil < 2; civil++) {
      Record r;
      r.t0 = parse_utc("2023-08-01T00:00:00Z");
      r.offsets = {{r.t0 - 86400, 10.0}, {parse_utc("2023-09-30T16:00:00Z"), 11.0}, {parse_utc("2024-04-06T16:00:00Z"), 10.0}};
      calx::random::Stream z(3, 3, 3, static_cast<std::uint32_t>(civil));
      for (std::size_t i = 0; i < 300 * 24; i++) {
        const double t = r.time(i);
        const double hour = civil ? r.civil(t) / 3600.0 : t / 3600.0;
        r.y.push_back(std::expm1(4 + std::cos(two_pi * (hour - 13) / 24) + 0.05 * z.normal()));
      }
      const auto shifts = clock_shifts(r, 21);
      bool ok = shifts.size() == 2;
      for (const auto& c : shifts) ok = ok && std::fabs(c.civil_share - (civil ? 1.0 : 0.0)) < 0.05;
      s.require(civil ? "civil share: a street keeping the civil clock has ψ ≈ 1 at both changes" : "civil share: a street keeping UTC has ψ ≈ 0 at both changes", ok,
                shifts.empty() ? "none" : std::to_string(shifts[0].civil_share));
    }
  }

  // ---- the paper's Figure 1, and the case it misreads ----
  {
    const auto cases = constructed_cases(1);
    s.require("Figure 1, case 1: proportional, by the paper and by the amendment", cases[0].paper_verdict == "proportional response" && cases[0].amended_verdict == "proportional response",
              cases[0].paper_verdict + " / " + cases[0].amended_verdict);
    s.require("Figure 1, case 2: disproportionate, by both", cases[1].paper_verdict == "disproportionate response" && cases[1].amended_verdict == "disproportionate response",
              cases[1].paper_verdict + " / " + cases[1].amended_verdict);
    s.require("Figure 1, case 3: a process of its own, by both", cases[2].paper_verdict == "a process of its own" && cases[2].amended_verdict.rfind("a process of its own", 0) == 0,
              cases[2].paper_verdict + " / " + cases[2].amended_verdict);
    s.require("Figure 1, case 4: a linear response to a non-separable forcing, disproportionate by the paper, proportional by the amendment",
              cases[3].paper_verdict == "disproportionate response" && cases[3].amended_verdict == "proportional response", cases[3].paper_verdict + " / " + cases[3].amended_verdict);
    s.require("Figure 1, case 1: the sum's power at the noise (paper: 4.7 × 10⁻⁵)", cases[0].at.at("f1+f2") < 1e-3, std::to_string(cases[0].at.at("f1+f2")));
  }

  // ---- the tide of a street ----
  {
    Record r;
    r.t0 = parse_utc("2024-10-01T00:00:00Z");
    r.offsets = {{r.t0 - 86400, 10.0}, {parse_utc("2024-10-05T16:00:00Z"), 11.0}, {parse_utc("2025-04-05T16:00:00Z"), 10.0},
                 {parse_utc("2025-10-04T16:00:00Z"), 11.0}, {parse_utc("2026-04-04T16:00:00Z"), 10.0}};
    const double year = frequency({0, 0, 1});
    auto week = [&](double h) {   // civil hours since the week's origin
      return 900 + 400 * std::cos(two_pi * (h - 13) / 24) + 160 * std::cos(two_pi * 2 * (h - 8.5) / 24) + 120 * std::cos(two_pi * (h - 100) / 168);
    };
    auto truth = [&](double h) { return week(h) + 80 * std::cos(two_pi * year * h); };
    for (std::size_t i = 0; i < 2 * 365 * 24; i++) r.y.push_back(truth(week_hours(r, r.time(i)) + 0.5));
    const auto t = tide(r, -1e18, 1e18);
    double worst = 0;
    for (std::size_t i = 0; i < r.size(); i += 37) {
      const double h = week_hours(r, r.time(i)) + 0.5;
      worst = std::max(worst, std::fabs(level(t, h) - truth(h)));
    }
    s.check("tide: the prediction of a street on the lattice is exact (people an hour)", worst, 1e-6);
    const auto R = readings(t);
    double best = -1, best_at = 0;
    for (int mnt = 0; mnt < 1440; mnt++) {
      double v = 0;
      for (int d = 0; d < 5; d++) v += week(24.0 * d + mnt / 60.0) / 5;
      if (v > best) best = v, best_at = mnt / 60.0;
    }
    s.check("tide: high water on the mean weekday, at the brute force's minute (h)", R.weekday.high_at - best_at, 1e-9);
    s.check("tide: and its level", R.weekday.high - best, 1e-6);
    s.check("tide: the form F = 400/160", R.form - 2.5, 1e-9);
    s.require("tide: mixed, mainly diurnal (Courtier)", R.form_class == "mixed, mainly diurnal", R.form_class);
    s.check("tide: the year's swing", R.annual - 80, 1e-6);
    // A holiday: 300 an hour fewer on one Wednesday. The calendars cannot predict it, so it is the surge.
    Record q = r;
    const long long hol = days_from_civil(2025, 7, 16), other = days_from_civil(2025, 7, 9);
    for (std::size_t i = 0; i < q.size(); i++)
      if (q.civil_day(q.time(i)) == hol) q.y[i] -= 300;
    const auto tq = tide(q, -1e18, 1e18);
    auto surge = [&](long long day) {
      const auto obs = observed_day(q, day);
      double m = 0;
      for (int h = 0; h < 24; h++) m += (obs[static_cast<std::size_t>(h)] - level(tq, day_start_hours(day) + h + 0.5)) / 24;
      return m;
    };
    s.check("tide: the holiday's surge is its dip, −300 an hour (to 3%)", (surge(hol) + 300) / 300, 0.03);
    // The fit takes the dip's projection onto the constituents (a Wednesday among 104 moves every Wednesday by
    // about 300/104): what it spreads to another day is a few per cent of the dip.
    s.check("tide: a week before, the surge is under 3% of the dip", surge(other) / 300, 0.03);
    const auto skipped = observed_day(r, days_from_civil(2024, 10, 6));
    s.require("tide: the hour daylight saving skips is missing from its day", !std::isfinite(skipped[2]) && std::isfinite(skipped[1]) && std::isfinite(skipped[3]));
  }

  // ---- bicoherence ----
  {
    const std::size_t L = 256, S = 64;
    std::vector<double> coupled, noise;
    calx::random::Stream z(11, 0, 0, 0);
    for (std::size_t seg = 0; seg < S; seg++) {
      const double p1 = z.uniform(0, two_pi), p2 = z.uniform(0, two_pi);
      for (std::size_t t = 0; t < L; t++) {
        const double tt = static_cast<double>(t);
        coupled.push_back(std::cos(two_pi * 20 * tt / L + p1) + std::cos(two_pi * 12 * tt / L + p2) + std::cos(two_pi * 32 * tt / L + p1 + p2) + 0.3 * z.normal());
        noise.push_back(z.normal());
      }
    }
    // The phases are drawn segment by segment, so the estimate's segments are taken without overlap.
    const auto bc = bicoherence(coupled, L, 40, 0.0);
    s.require("bicoherence: quadratic phase coupling at (20, 12) has b² > 0.8", bc.b2[20 * bc.bins + 12] > 0.8, std::to_string(bc.b2[20 * bc.bins + 12]));
    const auto bn = bicoherence(noise, L, 40, 0.0);
    s.require("bicoherence: independent noise exceeds the level in about α of the cells", bn.share_above < 0.12, std::to_string(bn.share_above));
  }

  return s.done();
}
