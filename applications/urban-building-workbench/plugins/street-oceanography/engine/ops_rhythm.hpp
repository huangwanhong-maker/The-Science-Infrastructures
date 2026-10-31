// The engine's operations on records: a sensor's rhythm on the calendars' lattice, the civil share of its
// clock, the bicoherence of what the calendars leave, the paper's Figure 1, the response to forcing, and the
// classification of every sensor in the plane.

#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#include <urban/classification/plane.hpp>
#include <urban/response/susceptibility.hpp>
#include <urban/rhythm/analysis.hpp>
#include <urban/street/terrain.hpp>

#include "data.hpp"
#include "ops_street.hpp"

namespace so {

namespace ur = urban::rhythm;

// The part of a record in an era: "recent" from 2024-10-01, "history" before 2023-01-01, "all".
inline ur::Record era(const ur::Record& r, const std::string& which) {
  const double from = which == "recent" ? ur::parse_utc("2024-10-01") : -1e18;
  const double to = which == "history" ? ur::parse_utc("2023-01-01") : 1e18;
  ur::Record out;
  out.id = r.id, out.step = r.step, out.offsets = r.offsets;
  std::size_t first = r.size(), last = 0;
  for (std::size_t i = 0; i < r.size(); i++)
    if (r.time(i) >= from && r.time(i) < to) first = std::min(first, i), last = std::max(last, i);
  if (first > last) throw std::invalid_argument("the record of sensor " + r.id + " has no " + which + " part");
  out.t0 = r.time(first);
  out.y.assign(r.y.begin() + static_cast<std::ptrdiff_t>(first), r.y.begin() + static_cast<std::ptrdiff_t>(last + 1));
  return out;
}

// The longest span with gaps of at most six hours, of at least 60 days.
inline ur::Span longest(const ur::Record& r) {
  const auto sp = ur::spans(r, 60 * 24, 6);
  if (sp.empty()) throw std::invalid_argument("no stretch of 60 days with gaps of at most six hours in the record of sensor " + r.id);
  return sp.front();
}

// A spectrum's arrays thinned to at most `most` points, keeping the maximum of each group (lines survive).
inline Json thinned(const std::vector<double>& f, const std::vector<std::vector<double>>& ys, std::size_t most, double to_cpd = 24) {
  const std::size_t n = f.size(), g = std::max<std::size_t>(1, (n + most - 1) / most);
  std::vector<double> fx;
  std::vector<std::vector<double>> out(ys.size());
  for (std::size_t i = 0; i < n; i += g) {
    std::size_t best = i;
    for (std::size_t k = i; k < std::min(n, i + g); k++)
      if (ys[0][k] > ys[0][best]) best = k;
    fx.push_back(f[best] * to_cpd);
    for (std::size_t s = 0; s < ys.size(); s++) out[s].push_back(ys[s][best]);
  }
  Json o = Json::object({{"frequency", numbers(fx, 6)}});
  const char* names[] = {"density", "lower", "upper", "background", "level"};
  for (std::size_t s = 0; s < ys.size() && s < 5; s++) o.set(names[s], numbers(out[s], 9));
  return o;
}

inline double variance_of(const std::vector<double>& x) {
  double m = 0, v = 0;
  for (double a : x) m += a;
  m /= static_cast<double>(x.size());
  for (double a : x) v += (a - m) * (a - m);
  return v / static_cast<double>(x.size());
}

inline Json harmonics_json(const ur::Harmonics& h, std::size_t top = 40) {
  std::vector<std::size_t> idx(h.fit.terms.size());
  for (std::size_t i = 0; i < idx.size(); i++) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return h.fit.terms[a].amplitude > h.fit.terms[b].amplitude; });
  Json terms = Json::array();
  for (std::size_t q = 0; q < std::min(top, idx.size()); q++) {
    const auto& t = h.fit.terms[idx[q]];
    const auto& c = h.constituents[idx[q]];
    terms.push(Json::object({{"name", c.name}, {"driving", ur::driving(c)}, {"order", c.order}, {"cpd", c.frequency * 24}, {"amplitude", t.amplitude},
                             {"amplitudeSe", t.amplitude_se}, {"phase", t.phase}, {"phaseSe", t.phase_se}, {"snr", t.snr}}));
  }
  return Json::object({{"terms", terms}, {"constituents", static_cast<double>(h.constituents.size())}, {"varianceExplained", h.fit.variance_explained},
                       {"observations", static_cast<double>(h.fit.observations)}, {"recordHours", h.fit.record_length}});
}

// ---- the rhythm of a sensor ----------------------------------------------------------------------------------

inline Json op_rhythm(const Json& r) {
  const long long id = static_cast<long long>(num(r, "sensor", 1));
  const std::string which = r["era"].text("recent");
  const double nw = num(r, "nw", 4);
  const int order = static_cast<int>(num(r, "order", 3));
  const bool month = r["month"].is_bool() && r["month"].as_bool();
  const auto& full = record(id);
  const auto rec = era(full, which);
  const auto span = longest(rec);
  const auto x = ur::series(rec, span, true);
  const auto tapers = static_cast<std::size_t>(std::max(1.0, 2 * nw - 1));
  const auto rep = ur::spectrum(x, nw, tapers, {true, month, true}, false, order, 0.3);
  Json peaks = Json::array();
  for (const auto& p : rep.peaks)
    peaks.push(Json::object({{"cpd", p.frequency * 24}, {"ratio", p.ratio}, {"density", p.density}, {"verdict", p.verdict},
                             {"name", p.explanation.found ? p.explanation.combination.name : ""}, {"order", p.explanation.found ? ur::cross_order(p.explanation.combination) : 0},
                             {"driving", p.explanation.found ? ur::driving(p.explanation.combination) : ""}}));
  Json chances = Json::array();
  for (int d = 0; d <= 5; d++) chances.push(ur::chance_rate_free_week({true, month, true}, d, std::max(rep.resolution, rep.bandwidth), 0.3));
  // Harmonic constants on the street constituents, and the selectivity.
  const double hours = static_cast<double>(span.end - span.begin);
  const auto cs = ur::street_constituents(hours, month);
  const auto pts = ur::points(rec, span.begin, span.end, rec.time(span.begin), true);
  const auto h = ur::harmonics(pts, cs, rec.time(span.begin), true);
  const auto sel = ur::selectivity(h, variance_of(pts.y));
  // The weekly profile (civil hours) of ln(1 + y).
  std::vector<double> profile(168, 0), count(168, 0);
  for (std::size_t i = span.begin; i < span.end; i++) {
    if (!std::isfinite(rec.y[i])) continue;
    const auto k = static_cast<std::size_t>(std::floor(ur::week_angle(rec.civil(rec.time(i))) * 168 + 1e-9)) % 168;
    profile[k] += std::log1p(std::max(0.0, rec.y[i])), count[k] += 1;
  }
  for (std::size_t k = 0; k < 168; k++) profile[k] = count[k] > 0 ? profile[k] / count[k] : std::nan("");
  // Daily totals over the whole record, for the overview.
  const auto totals = urban::response::daily_totals(full);
  Json days = Json::array(), values = Json::array();
  for (const auto& [d, v] : totals) days.push(static_cast<double>(d)), values.push(std::round(v));
  Json drivings = Json::array();
  for (std::size_t k = 0; k < sel.drivings.size(); k++) drivings.push(Json::object({{"driving", sel.drivings[k]}, {"share", sel.share[k]}}));
  const auto& s = sensor(id);
  return Json::object({
      {"sensor", static_cast<double>(id)}, {"name", s.name}, {"era", which},
      {"span", Json::object({{"from", ur::format_date(rec.time(span.begin))}, {"to", ur::format_date(rec.time(span.end - 1))},
                             {"hours", hours}, {"filled", static_cast<double>(span.filled)}})},
      {"spectrum", thinned(rep.frequency, {rep.density, rep.lower, rep.upper, rep.background, rep.level}, 3000)},
      {"null", Json::object({{"rho", rep.rho}, {"variance", rep.variance}})},
      {"nw", rep.nw}, {"tapers", static_cast<double>(rep.tapers)}, {"bandwidthCpd", rep.bandwidth * 24}, {"resolutionCpd", rep.resolution * 24},
      {"order", order}, {"chance", rep.chance}, {"chanceByOrder", chances}, {"peaks", peaks},
      {"harmonics", harmonics_json(h)}, {"selectivity", Json::object({{"S", sel.selectivity}, {"forced", sel.forced}, {"drivings", drivings}})},
      {"profile", numbers(profile, 4)}, {"daily", Json::object({{"day", days}, {"total", values}})}});
}

inline Json op_clock(const Json& r) {
  const long long id = static_cast<long long>(num(r, "sensor", 1));
  const auto& rec = record(id);
  std::vector<long long> skip;
  try { skip = holiday_days(); } catch (...) {}
  const auto shifts = ur::clock_shifts(rec, static_cast<int>(num(r, "window", 21)), skip);
  Json list = Json::array();
  double wsum = 0, w = 0;
  for (const auto& c : shifts) {
    list.push(Json::object({{"date", ur::format_date(c.at)}, {"change", c.change}, {"before", c.before}, {"after", c.after}, {"shift", c.shift},
                            {"se", c.se}, {"civilShare", c.civil_share}, {"daysBefore", static_cast<double>(c.days_before)}, {"daysAfter", static_cast<double>(c.days_after)}}));
    if (c.se > 0 && c.days_before >= 10 && c.days_after >= 10) {
      const double q = 1 / (c.se * c.se);
      wsum += q * c.civil_share, w += q;
    }
  }
  return Json::object({{"sensor", static_cast<double>(id)}, {"name", sensor(id).name}, {"shifts", list},
                       {"civilShare", w > 0 ? wsum / w : std::nan("")}, {"civilShareSe", w > 0 ? 1 / std::sqrt(w) : std::nan("")}});
}

inline Json op_bicoherence(const Json& r) {
  const long long id = static_cast<long long>(num(r, "sensor", 1));
  const auto rec = era(record(id), r["era"].text("recent"));
  const auto span = longest(rec);
  const double hours = static_cast<double>(span.end - span.begin);
  const auto pts = ur::points(rec, span.begin, span.end, rec.time(span.begin), true);
  const auto h = ur::harmonics(pts, ur::street_constituents(hours), rec.time(span.begin), false);
  // The residual on the span's hours (points skip the gaps; they are few and short: filled by zero residual).
  std::vector<double> residual(span.end - span.begin, 0);
  for (std::size_t k = 0; k < pts.t.size(); k++) {
    const auto i = static_cast<std::size_t>(std::llround(pts.t[k]));
    if (i < residual.size()) residual[i] = h.fit.residual[k];
  }
  const std::size_t L = static_cast<std::size_t>(num(r, "segment", 336)), max_bin = static_cast<std::size_t>(num(r, "maxBin", 56));
  const auto b = ur::bicoherence(residual, L, max_bin, 0.5);
  Json grid = Json::array();
  for (std::size_t j = 0; j < b.bins; j++) {
    Json row = Json::array();
    for (std::size_t k = 0; k < b.bins; k++) row.push(k <= j && j + k < b.bins ? Json(std::round(b.b2[j * b.bins + k] * 1e4) / 1e4) : Json());
    grid.push(row);
  }
  Json above = Json::array();
  for (const auto& c : b.above) above.push(Json::array({b.frequency[c.j] * 24, b.frequency[c.k] * 24, c.b2}));
  return Json::object({{"sensor", static_cast<double>(id)}, {"segments", static_cast<double>(b.segments)}, {"segmentHours", static_cast<double>(L)},
                       {"cpd", numbers([&] { std::vector<double> f; for (double v : b.frequency) f.push_back(v * 24); return f; }(), 4)},
                       {"b2", grid}, {"level", b.level}, {"shareAbove", b.share_above}, {"above", above}});
}

inline Json op_constructed() {
  Json cases = Json::array();
  for (const auto& c : ur::constructed_cases(1)) {
    Json at = Json::object();
    for (const auto& [k, v] : c.at) at.set(k, v);
    std::vector<double> f, p;
    for (std::size_t i = 0; i < c.frequency.size(); i++) f.push_back(c.frequency[i]), p.push_back(c.power[i]);
    cases.push(Json::object({{"name", c.name}, {"about", c.about}, {"frequency", numbers(f, 5)}, {"power", numbers(p, 10)}, {"at", at},
                             {"paper", c.paper_verdict}, {"amended", c.amended_verdict}, {"residualLines", strings(c.residual_lines)}}));
  }
  return cases;
}

// ---- every sensor's week ------------------------------------------------------------------------------------

// The mean count at each civil hour of the week over an era, for every outdoor sensor with at least four weeks.
inline Json op_profiles(const Json& r) {
  const std::string which = r["era"].text("recent");
  static std::map<std::string, Json> cache;
  if (cache.count(which)) return cache.at(which);
  Json list = Json::array();
  for (const auto& s : sensors()) {
    if (s.location_type != "Outdoor") continue;
    try {
      const auto rec = era(record(s.id), which);
      std::vector<double> sum(168, 0), n(168, 0);
      for (std::size_t i = 0; i < rec.size(); i++) {
        if (!std::isfinite(rec.y[i])) continue;
        const auto k = static_cast<std::size_t>(std::floor(ur::week_angle(rec.civil(rec.time(i))) * 168 + 1e-9)) % 168;
        sum[k] += rec.y[i], n[k] += 1;
      }
      double weeks = 1e9;
      for (std::size_t k = 0; k < 168; k++) weeks = std::min(weeks, n[k]), sum[k] = n[k] > 0 ? sum[k] / n[k] : std::nan("");
      if (weeks < 4) continue;
      list.push(Json::object({{"id", static_cast<double>(s.id)}, {"name", s.name}, {"lat", s.lat}, {"lon", s.lon}, {"weeks", weeks}, {"profile", numbers(sum, 1)}}));
    } catch (...) {
    }
  }
  cache[which] = list;
  return list;
}

// ---- the response to forcing -------------------------------------------------------------------------------

struct ResponseSummary {
  urban::response::ResponseFit fit;
  urban::response::PairedEffect holidays;
};

inline ResponseSummary response_of(long long id, const std::string& which) {
  const auto rec = era(record(id), which);
  const auto& w = weather();
  auto rain = aligned(w.precipitation, w, rec);
  auto temp = aligned(w.temperature, w, rec);
  // The temperature's departure from its mean for the civil month and hour over the record.
  std::vector<double> sum(12 * 24, 0), cnt(12 * 24, 0);
  auto slot = [&](std::size_t i) {
    const auto c = ur::civil_from_days(rec.civil_day(rec.time(i)));
    return (c.month - 1) * 24 + static_cast<unsigned>(std::floor(rec.civil_hour(rec.time(i))));
  };
  for (std::size_t i = 0; i < rec.size(); i++)
    if (std::isfinite(temp[i])) sum[slot(i)] += temp[i], cnt[slot(i)] += 1;
  for (std::size_t i = 0; i < rec.size(); i++)
    if (std::isfinite(temp[i])) temp[i] -= sum[slot(i)] / cnt[slot(i)];
  std::vector<double> holiday(rec.size(), 0);
  std::vector<long long> hd;
  try { hd = holiday_days(); } catch (...) {}
  for (std::size_t i = 0; i < rec.size(); i++)
    if (std::binary_search(hd.begin(), hd.end(), rec.civil_day(rec.time(i)))) holiday[i] = 1;
  ResponseSummary s;
  s.fit = urban::response::respond(rec, {{"rain", "mm/h", rain, {0, 1, 2, 3}}, {"temperature departure", "°C", temp, {0}}, {"public holiday", "day", holiday, {0}}});
  std::vector<long long> weekday_holidays;
  for (long long d : hd) {
    const int wd = static_cast<int>(((d - 4) % 7 + 7) % 7);
    if (wd < 5) weekday_holidays.push_back(d);
  }
  s.holidays = urban::response::paired_effect(rec, weekday_holidays, hd, 4999, 1);
  return s;
}

inline Json susceptibility_json(const urban::response::Susceptibility& f) {
  return Json::object({{"name", f.name}, {"unit", f.unit}, {"lags", numbers(std::vector<double>(f.lags.begin(), f.lags.end()), 0)}, {"beta", numbers(f.beta)},
                       {"se", numbers(f.se)}, {"total", f.total}, {"totalSe", f.total_se}, {"z", f.z}, {"percent", f.percent}});
}

inline Json op_response(const Json& r) {
  const long long id = static_cast<long long>(num(r, "sensor", 1));
  const auto s = response_of(id, r["era"].text("recent"));
  Json forcings = Json::array();
  for (const auto& f : s.fit.forcings) forcings.push(susceptibility_json(f));
  Json pairs = Json::array();
  for (std::size_t k = 0; k < s.holidays.differences.size(); k++)
    pairs.push(Json::array({ur::format_date(static_cast<double>(s.holidays.a_days[k]) * 86400), ur::format_date(static_cast<double>(s.holidays.b_days[k]) * 86400),
                            s.holidays.differences[k]}));
  return Json::object({{"sensor", static_cast<double>(id)}, {"name", sensor(id).name}, {"forcings", forcings}, {"observations", static_cast<double>(s.fit.observations)},
                       {"days", static_cast<double>(s.fit.days)}, {"r2Within", s.fit.r2_within}, {"sigma", s.fit.sigma}, {"trend", s.fit.trend},
                       {"holidays", Json::object({{"pairs", static_cast<double>(s.holidays.pairs)}, {"effect", s.holidays.effect}, {"se", s.holidays.se},
                                                  {"percent", s.holidays.percent}, {"p", s.holidays.p}, {"list", pairs}})}});
}

// ---- every sensor: susceptibilities and the plane ------------------------------------------------------------

struct SensorSummary {
  long long id = 0;
  std::string name;
  double rain = std::nan(""), rain_se = std::nan(""), holiday = std::nan(""), holiday_p = std::nan("");
  double selectivity = std::nan(""), forced = std::nan(""), compression = std::nan(""), length = std::nan("");
  std::string error;
};

inline const std::vector<SensorSummary>& summaries() {
  static std::vector<SensorSummary> list;
  static bool done = false;
  static std::mutex lock;
  std::lock_guard<std::mutex> g(lock);
  if (done) return list;
  const Terrain* t = nullptr;
  try { t = &terrain(); } catch (...) {}
  for (const auto& s : sensors()) {
    if (s.location_type != "Outdoor") continue;
    SensorSummary m;
    m.id = s.id, m.name = s.name;
    try {
      const auto rec = era(record(s.id), "recent");
      const auto span = longest(rec);
      if (span.end - span.begin < 300 * 24) throw std::invalid_argument("less than 300 days");
      const double hours = static_cast<double>(span.end - span.begin);
      const auto pts = ur::points(rec, span.begin, span.end, rec.time(span.begin), true);
      const auto h = ur::harmonics(pts, ur::street_constituents(hours), rec.time(span.begin), false);
      const auto sel = ur::selectivity(h, variance_of(pts.y));
      m.selectivity = sel.selectivity, m.forced = sel.forced;
      const auto resp = response_of(s.id, "recent");
      m.rain = resp.fit.forcings[0].total, m.rain_se = resp.fit.forcings[0].total_se;
      m.holiday = resp.holidays.effect, m.holiday_p = resp.holidays.p;
    } catch (const std::exception& e) {
      m.error = e.what();
    }
    if (t && t->by_sensor.count(s.id)) {
      std::vector<double> x, y;
      std::vector<std::size_t> u;
      for (const auto& e : t->by_sensor.at(s.id))
        if (e.year == latest_year(*t)) x.push_back(e.x), y.push_back(e.y), u.push_back(t->use_index(e.use));
      if (x.size() >= 20) {
        const auto c = urban::street::use_correlation(x, y, u, 10, 200);
        m.length = c.length;
        m.compression = urban::street::compression(400, c);
      }
    }
    list.push_back(m);
  }
  done = true;
  return list;
}

inline Json op_sensors_summary() {
  Json list = Json::array();
  for (const auto& m : summaries())
    list.push(Json::object({{"id", static_cast<double>(m.id)}, {"name", m.name}, {"rain", m.rain}, {"rainSe", m.rain_se}, {"holiday", m.holiday},
                            {"holidayP", m.holiday_p}, {"selectivity", m.selectivity}, {"forced", m.forced}, {"compression", m.compression},
                            {"length", m.length}, {"error", m.error}}));
  return list;
}

inline Json op_plane() {
  std::vector<urban::classification::Segment> segs;
  for (const auto& m : summaries()) segs.push_back({std::to_string(m.id), m.name, m.compression, m.selectivity});
  const auto p = urban::classification::plane(segs);
  Json pts = Json::array();
  for (std::size_t i = 0; i < p.segments.size(); i++)
    pts.push(Json::object({{"id", p.segments[i].id}, {"name", p.segments[i].name}, {"compression", p.segments[i].compression},
                           {"selectivity", p.segments[i].selectivity}, {"x", p.x[i]}, {"y", p.y[i]}}));
  return Json::object({{"points", pts}, {"ratio", p.test.ratio}, {"z", p.test.z}, {"p", p.test.p}, {"reading", p.reading},
                       {"observed", p.test.observed}, {"expected", p.test.expected}, {"walk", 400.0}});
}

}  // namespace so
