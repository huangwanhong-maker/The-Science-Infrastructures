// The shared problem "A district's economy" made machine-readable (README.md): the district, its frontages and their
// relations from the domains; the people's rhythm calibrated from the City of Melbourne's counts by the urban-science
// domain; the affinities of the kinds from CLUE's establishments; and the problem written to results/problem.json.
//
//   calibrate <results/problem-inputs.json> <results directory>
//
// From the repository's root the inputs' paths to the data library's cache resolve; run it there.
//
// For each counter in the Hoddle Grid's box covering at least 80% of the reference year (52 weeks from Monday
// 27 October 2025, civil):
//   - its tide (urban::rhythm::tide): the street constituents fitted by least squares to the counts, people an hour,
//     over the reference year with Victoria's public holidays left out; its mean week, the week's constituents alone
//     (k_M = k_Y = 0), at the middle of each of the 168 hours; its readings (high and low water, Courtier's form);
//   - its spectrum (urban::rhythm::spectrum): Thomson's multitaper estimate (NW = 4, seven tapers; calx `signal`) of
//     ln(1 + counts) over the longest stretch of the year with no gap above six hours, the red-noise level (tyche),
//     and each peak's line on the calendars' lattice (week and year, cross order 2);
//   - its response (urban::response::respond): ln(1 + counts) on the hour of the week's 168 levels, the year's first
//     two harmonics, rain (mm in the hour and the hour before), the air's temperature and the public holidays,
//     with errors clustered by civil day; each forcing's susceptibility as a percentage change of the counts.
// The kinds' affinities (urban::affinity::affinity, on tyche's co-location quotients): CLUE's 2024 establishments in
// the Hoddle Grid's box, marked by kind, neighbours within 50 m, 199 permutations (free, and within 200 m cells).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <urban/affinity/colocation.hpp>
#include <urban/response/susceptibility.hpp>
#include <urban/rhythm/analysis.hpp>
#include <urban/rhythm/tide.hpp>

#include "economy.hpp"

namespace {

using economy::Json;
using economy::numbers;
namespace ur = urban::rhythm;

std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}
void write_file(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  f << text;
  if (!f) throw std::runtime_error("cannot write " + path);
}
double r(double x, double d) { return std::isfinite(x) ? std::round(x * d) / d : std::nan(""); }
Json num(double x, double d = 1e4) { return std::isfinite(x) ? Json(r(x, d)) : Json(); }

// A sensor's record from the data library's street-counts product: UTC hours, NaN where missing, civil offsets.
ur::Record record_of(const Json& c) {
  ur::Record rec;
  rec.id = std::to_string(static_cast<long long>(c["id"].as_number()));
  rec.t0 = ur::parse_utc(c["t0_utc"].text());
  rec.step = c["step_s"].as_number();
  for (const auto& v : c["counts"].as_array()) rec.y.push_back(v.is_null() ? std::nan("") : v.as_number());
  for (const auto& o : c["offsets"].as_array()) rec.offsets.push_back({ur::parse_utc(o[0].text()), o[1].as_number()});
  return rec;
}

// The record cut to [from, from + hours) of UTC.
ur::Record window(const ur::Record& rec, double from, std::size_t hours) {
  ur::Record w;
  w.id = rec.id, w.t0 = from, w.step = rec.step, w.offsets = rec.offsets;
  w.y.assign(hours, std::nan(""));
  const long long i0 = static_cast<long long>(std::llround((from - rec.t0) / rec.step));
  for (std::size_t i = 0; i < hours; i++) {
    const long long j = i0 + static_cast<long long>(i);
    if (j >= 0 && j < static_cast<long long>(rec.size())) w.y[i] = rec.y[static_cast<std::size_t>(j)];
  }
  return w;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: calibrate <results/problem-inputs.json> <results directory>\n";
    return 2;
  }
  const auto T0 = std::chrono::steady_clock::now();
  const std::string out = argv[2];
  const Json in = Json::parse(read_file(argv[1]));
  economy::Problem P(in);
  std::printf("district: %zu frontages (%zu front a street), %zu candidate links, %zu party walls among the parcels, %zu with the context; %zu of %zu buildings massed; %zu segments dropped\n",
              P.size(), static_cast<std::size_t>(std::count_if(P.frontages.begin(), P.frontages.end(), [](const auto& f) { return f.fronts; })), P.links.size(), P.party.size(),
              P.party_with_context, P.massed, P.size() + P.context_buildings, P.dropped.size());

  // ---- the people: the counters' rhythm and response ---------------------------------------------------------------
  const double ref0 = ur::parse_utc(in["reference"]["start_utc"].text());
  const std::size_t hours = static_cast<std::size_t>(in["reference"]["weeks"].as_number()) * 168;
  std::set<long long> holidays;
  {
    const Json H = Json::parse(read_file(in["files"]["holidays"].text()));
    for (const auto& h : H["holidays"].as_array()) holidays.insert(static_cast<long long>(std::floor(ur::parse_utc(h["date"].text()) / ur::day_s)));
  }
  const Json W = Json::parse(read_file(in["files"]["weather"].text()));
  const double w0 = ur::parse_utc(W["t0_utc"].text());
  const auto& rain = W["precipitation_mm"].as_array();
  const auto& temp = W["temperature_c"].as_array();
  const auto weather_at = [&](const Json::Array& a, double t) {
    const long long i = static_cast<long long>(std::llround((t - w0) / 3600.0));
    return i >= 0 && i < static_cast<long long>(a.size()) && !a[static_cast<std::size_t>(i)].is_null() ? a[static_cast<std::size_t>(i)].as_number() : std::nan("");
  };
  const Json& box = in["district"]["district_box"];
  Json sensors = Json::array(), spectra = Json::array();
  std::size_t used = 0;
  for (const auto& s : in["sensors"].as_array()) {
    const double coverage = s["coverage"].as_number();
    if (coverage < 0.8) continue;
    const ur::Record full = record_of(Json::parse(read_file(s["series"].text())));
    const ur::Record rec = window(full, ref0, hours);
    ur::Record typical = rec;  // the public holidays left out of the tide
    std::size_t holiday_hours = 0;
    for (std::size_t i = 0; i < typical.size(); i++)
      if (holidays.count(typical.civil_day(typical.time(i)))) typical.y[i] = std::nan(""), holiday_hours++;
    const auto tide = ur::tide(typical, ref0, ref0 + static_cast<double>(hours) * 3600.0);
    const auto read = ur::readings(tide);
    std::vector<double> mean_week(168);
    for (std::size_t h = 0; h < 168; h++) mean_week[h] = ur::level(tide, static_cast<double>(h) + 0.5, ur::of_the_week);
    // The largest constituents.
    std::vector<std::size_t> order(tide.fit.terms.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return tide.fit.terms[a].amplitude > tide.fit.terms[b].amplitude; });
    Json cons = Json::array();
    for (std::size_t q = 0; q < std::min<std::size_t>(order.size(), 16); q++) {
      const auto& t = tide.fit.terms[order[q]];
      cons.push(Json::object({{"name", tide.constituents[order[q]].name}, {"cycles_per_day", r(t.frequency * 24, 1e6)}, {"amplitude", r(t.amplitude, 100)},
                              {"phase_deg", r(t.phase * 180 / std::numbers::pi, 100)}, {"amplitude_se", r(t.amplitude_se, 100)}}));
    }
    // The spectrum of ln(1 + counts) over the year's longest stretch without a gap above six hours.
    Json peaks = Json::array();
    Json spec;
    const auto sp = ur::spans(rec, 28 * 24, 6);
    double span_hours = 0;
    if (!sp.empty()) {
      const auto x = ur::series(rec, sp[0], true);
      span_hours = static_cast<double>(x.size());
      const auto S = ur::spectrum(x, 4, 7, {true, false, true}, true, 2, 0.5, 1.0, 40);
      for (std::size_t q = 0; q < std::min<std::size_t>(S.peaks.size(), 12); q++) {
        const auto& pk = S.peaks[q];
        peaks.push(Json::object({{"cycles_per_day", r(pk.frequency * 24, 1e5)}, {"ratio", r(pk.ratio, 100)}, {"verdict", pk.verdict},
                                 {"line", pk.explanation.found ? Json(pk.explanation.combination.name) : Json()}}));
      }
      const double x0 = s["x"].as_number(), y0 = s["y"].as_number();
      if (x0 >= box["x0"].as_number() && x0 <= box["x1"].as_number() && y0 >= box["y0"].as_number() && y0 <= box["y1"].as_number()) {
        // Kept for drawing, every fourth frequency below 4 cycles a day.
        std::vector<double> f, d, lv;
        for (std::size_t j = 0; j < S.frequency.size() && S.frequency[j] * 24 <= 4.2; j += 4) f.push_back(S.frequency[j] * 24), d.push_back(S.density[j]), lv.push_back(S.level[j]);
        spec = Json::object({{"id", s["id"]}, {"name", s["name"]}, {"cycles_per_day", numbers(f, 5)}, {"density", numbers(d, 8)}, {"level", numbers(lv, 8)},
                             {"bandwidth_cycles_per_day", r(S.bandwidth * 24, 1e5)}, {"hours", span_hours}});
        spectra.push(spec);
      }
    }
    // The response to rain, temperature and holidays.
    std::vector<urban::response::Forcing> forcings(3);
    forcings[0] = {"rain", "mm in the hour", {}, {0, 1}};
    forcings[1] = {"temperature", "°C", {}, {0}};
    forcings[2] = {"holiday", "a public holiday in Victoria", {}, {0}};
    for (std::size_t i = 0; i < rec.size(); i++) {
      const double t = rec.time(i);
      forcings[0].value.push_back(weather_at(rain, t));
      forcings[1].value.push_back(weather_at(temp, t));
      forcings[2].value.push_back(holidays.count(rec.civil_day(t)) ? 1.0 : 0.0);
    }
    Json response = Json::object();
    try {
      const auto fit = urban::response::respond(rec, forcings, true, false);
      for (const auto& f : fit.forcings) response.set(f.name, Json::object({{"percent", r(f.percent, 100)}, {"total", r(f.total, 1e5)}, {"se", r(f.total_se, 1e5)}, {"z", r(f.z, 100)}, {"unit", f.unit}}));
      response.set("r2_within", r(fit.r2_within, 1e4));
      response.set("days", static_cast<double>(fit.days));
    } catch (const std::exception& e) {
      response.set("error", std::string(e.what()));
    }
    double total = 0;
    for (double v : mean_week) total += v;
    sensors.push(Json::object({{"id", s["id"]}, {"code", s["code"]}, {"name", s["name"]}, {"x", s["x"]}, {"y", s["y"]}, {"way", s["way"]}, {"coverage", coverage},
                               {"holiday_hours_left_out", static_cast<double>(holiday_hours)}, {"observations", static_cast<double>(tide.observations)},
                               {"mean_week", numbers(mean_week, 2)}, {"weekly_total", r(total, 1)},
                               {"tide", Json::object({{"mean", r(read.mean, 100)}, {"variance_explained", r(read.variance_explained, 1e4)}, {"form", r(read.form, 1e4)}, {"form_class", read.form_class},
                                                      {"weekday", Json::object({{"high", r(read.weekday.high, 10)}, {"high_at", r(read.weekday.high_at, 100)}, {"low", r(read.weekday.low, 10)}, {"low_at", r(read.weekday.low_at, 100)}, {"mean", r(read.weekday.mean, 10)}})},
                                                      {"weekend", Json::object({{"high", r(read.weekend.high, 10)}, {"high_at", r(read.weekend.high_at, 100)}, {"low", r(read.weekend.low, 10)}, {"low_at", r(read.weekend.low_at, 100)}, {"mean", r(read.weekend.mean, 10)}})},
                                                      {"annual_amplitude", r(read.annual, 100)}, {"constituents_fitted", static_cast<double>(tide.constituents.size())}})},
                               {"constituents", cons}, {"spectrum", Json::object({{"hours", span_hours}, {"peaks", peaks}})}, {"response", response}}));
    used++;
  }
  std::printf("%zu counters calibrated\n", used);

  // ---- the kinds' affinities -------------------------------------------------------------------------------------------
  std::vector<double> ex, ey;
  std::vector<std::size_t> ek;
  for (const auto& e : in["establishments"]["rows"].as_array()) ex.push_back(e[0].as_number()), ey.push_back(e[1].as_number()), ek.push_back(static_cast<std::size_t>(e[2].as_number()));
  std::vector<std::string> names;
  for (std::size_t k = 0; k < economy::kinds; k++) names.push_back(economy::kind_name(static_cast<std::int64_t>(k)));
  const auto A = urban::affinity::affinity(ex, ey, ek, names, 200, 50, 199, 7);
  Json a = Json::array(), clq = Json::array(), pu = Json::array(), pl = Json::array(), pus = Json::array();
  for (std::size_t i = 0; i < economy::kinds; i++) {
    Json ra = Json::array(), rc = Json::array(), ru = Json::array(), rl = Json::array(), rs = Json::array();
    for (std::size_t j = 0; j < economy::kinds; j++) {
      const double c = A.free.clq[i * economy::kinds + j];
      const double v = std::isfinite(c) && c > 0 ? std::clamp(std::log(c), -1.0, 1.0) : 0.0;
      P.P.affinity[i][j] = v;
      ra.push(r(v, 1e5)), rc.push(num(c, 1e5)), ru.push(num(A.free.p_upper[i * economy::kinds + j], 1e4)), rl.push(num(A.free.p_lower[i * economy::kinds + j], 1e4));
      rs.push(num(A.stratified.p_upper[i * economy::kinds + j], 1e4));
    }
    a.push(ra), clq.push(rc), pu.push(ru), pl.push(rl), pus.push(rs);
  }
  std::printf("affinities from %zu establishments\n", ex.size());

  // ---- the problem ----------------------------------------------------------------------------------------------------
  const auto& PP = P.P;
  const auto asm_ = [](const char* what) { return std::string("assumption: ") + what; };
  Json kinds_j = Json::array();
  for (std::size_t k = 0; k < economy::kinds; k++)
    kinds_j.push(Json::object({{"index", static_cast<double>(k)}, {"name", names[k]}, {"staff", PP.staff[k]}, {"propensity", PP.propensity[k]},
                               {"period_day_evening_night", numbers(PP.period[k], 3)}, {"spend_aud", PP.spend[k]}, {"margin", PP.margin[k]}, {"outgoings_aud_m2_week", PP.outgoings[k]},
                               {"hours_as_built", in["typical_hours"][names[k]]}}));
  Json params = Json::object({
      {"wage_aud_hour", Json::object({{"value", PP.wage}, {"source", asm_("a staff member's hourly cost with on-costs, near the adult rates of the retail and hospitality awards with superannuation")}})},
      {"staff", Json::object({{"source", asm_("staff on duty while open: retail 2, food and drink 3, services 2, culture 2; lobbies are the building's, not counted")}})},
      {"propensity", Json::object({{"source", asm_("a passer's chance to stop at an open frontage by day, in thousandths: retail 6, food and drink 10, services 2.5, culture 1.5")}})},
      {"period", Json::object({{"bands", "day 07-18, evening 18-24, night 00-07 (civil)"}, {"source", asm_("the multipliers of the propensity by band, after the purposes of passers at those hours")}})},
      {"spend_aud", Json::object({{"source", asm_("spending a visit: retail 45, food and drink 18, services 55, culture 30 (A$ 2026)")}})},
      {"margin", Json::object({{"source", asm_("gross margin, the share of spending not paid for goods: retail 0.40, food and drink 0.65, services 0.70, culture 0.75")}})},
      {"outgoings_aud_m2_week", Json::object({{"source", asm_("utilities, cleaning and insurance of a trading floor; rent is a transfer within the district and left out")}})},
      {"refit_aud_m2", Json::object({{"value", PP.refit}, {"years", PP.refit_years}, {"rate", PP.rate}, {"weekly_per_m2", r(PP.refit_weekly(1), 1e4)}, {"source", asm_("the fit-out of a change of kind (to any but vacant), amortised")}})},
      {"link_aud_m", Json::object({{"value", PP.link}, {"years", PP.link_years}, {"rate", PP.rate}, {"weekly_per_m", r(PP.link_weekly(1), 1e4)}, {"source", asm_("an arcade through a building, amortised")}})},
      {"depth_m", Json::object({{"value", PP.depth}, {"source", asm_("the trading floor is at most the frontage's length times this depth, and at most the footprint")}})},
      {"catalysis", Json::object({{"gamma", PP.catalysis}, {"coupling_length_m", r(PP.coupling, 100)}, {"clamp", numbers(std::vector<double>{PP.catalysis_floor, PP.catalysis_ceiling}, 2)},
                                  {"source", "gamma and the clamp: assumptions; the coupling length: the urban-science domain's catalysis (urban::affinity::walk_coupling_length at D = 1 m2/s, leave rate 1/1800 per s, spacing 10 m), itself a representation; a(i->j): ln CLQ of CLUE 2024, clamped to [-1, 1]"}})},
      {"access_cap_m", Json::object({{"value", PP.access_cap}, {"source", asm_("the walk counted when no frontage of a need is open in the area")}})},
  });
  Json fr = Json::array();
  for (std::size_t f = 0; f < P.size(); f++) {
    const auto& x = P.frontages[f];
    Json seg = Json::array();
    for (auto g : P.on_segment[f]) seg.push(Json(static_cast<double>(g)));
    fr.push(Json::object({{"index", static_cast<double>(f)}, {"parcel", x.parcel}, {"name", x.name}, {"clue_property", x.clue_property.empty() ? Json() : Json(x.clue_property)},
                          {"street", x.street}, {"way", x.way}, {"segment", std::to_string(x.segment)}, {"fronts", x.fronts}, {"length_m", r(x.length, 10)}, {"facade_m", r(x.facade, 10)},
                          {"footprint_m2", r(x.footprint_area, 10)}, {"trading_area_m2", r(x.trading_area, 10)}, {"height_m", x.height}, {"block", static_cast<double>(x.block)},
                          {"kind_as_built", economy::kind_name(x.kind0)}, {"kind_read", x.read0},
                          {"hours_as_built", Json::array({static_cast<double>(x.hours0[0]), static_cast<double>(x.hours0[1]), static_cast<double>(x.hours0[2]), static_cast<double>(x.hours0[3])})},
                          {"same_segment", seg}}));
  }
  Json party = Json::array();
  for (std::size_t q = 0; q < P.party.size(); q++) party.push(Json::array({static_cast<double>(P.party[q][0]), static_cast<double>(P.party[q][1]), r(P.party_area[q], 10)}));
  const auto D0 = P.distances(economy::hv::Bits(P.links.size(), 0));
  Json reach = Json::array();
  std::size_t reach_edges = 0;
  for (std::size_t f = 0; f < P.size(); f++)
    for (const auto& [g, w] : D0->near[f])
      if (g > f) reach.push(Json::array({static_cast<double>(f), static_cast<double>(g), r(D0->at(f, g), 10)})), reach_edges++;
  Json links = Json::array();
  for (const auto& l : P.links) links.push(Json::object({{"kind", l.kind}, {"length_m", r(l.length, 10)}, {"joins", Json::array({l.from, l.to})}, {"weekly_cost_aud", r(PP.link_weekly(l.length), 10)}}));
  const auto plan0 = P.as_built();
  Json doc = Json::object({
      {"problem", "a-district-economy"},
      {"about", "The shared problem of three case studies of planning a district's economy (research/experiments/2026-11-05-a-district-economy-the-problem/README.md): the ground-floor frontages of the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street (Melbourne's Hoddle Grid), their kinds and hours, and eleven links. A study of methods on a real place, not a planning proposal."},
      {"units", Json::object({{"value", "A$ per week of the mean week"}, {"time", "the civil week of 168 hours, Monday 00:00 first"}, {"distance", "m"}, {"passers", "people an hour passing a frontage"}})},
      {"kinds", kinds_j},
      {"parameters", params},
      {"encoding", Json::object({{"kind", "per frontage, an index of kinds"},
                                 {"hours", "per frontage, four integers: weekday start (0-23), weekday length (0-24), weekend start, weekend length; closed when the length is 0; open over [start, min(24, start + length))"},
                                 {"links", "per candidate link, 0 or 1"},
                                 {"relaxation", "u in [0,1]^(F*10 + L): per frontage six kind scores (argmax, first on ties) and four hour numbers (start floor(24u) capped 23, length floor(25u) capped 24); per link u >= 0.5"},
                                 {"dimension", static_cast<double>(P.dimension())}})},
      {"value", "V = sum over frontages and hours open of P_f(h) p_f(h) c_f(h) s_k m_k minus staffing (staff x wage x (open hours + one per day opened)), outgoings (per m2 of trading floor), refits (a change of kind to any but vacant) and the links' costs; p_f(h) = propensity x period multiplier x clamp(1 + gamma sum_g w_fg a(k_g->k_f), 0.5, 2) over open non-vacant neighbours within 3 coupling lengths; c_f(h) = (1 - prod(1 - p_g)) / sum(p_g) over the open frontages of f's kind on f's segment"},
      {"objectives", Json::object({{"value", "maximised"}, {"spread", "maximised: the mean line's share of the power of the spending's 168 hours (Parseval), 1/(1 + CV^2)"},
                                   {"mix", "maximised: H(K | block) of the non-vacant frontages, weighted by length (nats)"}, {"access", "minimised: the people-weighted walk to the nearest open retail, food and drink, and services frontage, averaged over the needs and the week's hours (m)"},
                                   {"descriptors", "evening share (18-24 h) and weekend share of the spending"}})},
      {"frontages", fr},
      {"relations", Json::object({{"party_walls", party}, {"party_walls_with_context", static_cast<double>(P.party_with_context)}, {"reach_within_3_coupling_lengths", reach},
                                  {"segments_fronted", static_cast<double>(P.segments.size())}, {"blocks", static_cast<double>(P.place_count)}})},
      {"links", links},
      {"as_built", economy::plan_json(plan0)},
      {"affinity", Json::object({{"a", a}, {"clq", clq}, {"p_upper_free", pu}, {"p_lower_free", pl}, {"p_upper_stratified", pus}, {"radius_m", 50.0}, {"cell_m", 200.0}, {"permutations", 199.0},
                                 {"establishments", static_cast<double>(ex.size())}, {"year", in["establishments"]["year"]}, {"rows_are", "a[i][j] = clamp(ln CLQ_{i->j}, -1, 1): j's over-representation among i's neighbours"}})},
      {"people", Json::object({{"reference", in["reference"]}, {"counters", static_cast<double>(used)}, {"file", "results/sensors.json"}})},
      {"report", "spatial (maps of the district), temporal (hourly and weekly profiles), frequency (spectra and constituents of flows and spending), people flow, and the arrangement of buildings and their relations"},
  });
  write_file(out + "/problem.json", doc.dump());
  write_file(out + "/sensors.json", Json::object({{"reference", in["reference"]}, {"sensors", sensors}, {"spectra", spectra},
                                                  {"credit", "City of Melbourne, Pedestrian Counting System (the City's open data, CC BY 4.0); weather: Open-Meteo (CC BY 4.0) from ERA5; holidays: Nager.Date"}})
                                        .dump());
  // Drawing geometry (extracts of OpenStreetMap, ODbL: not committed).
  {
    Json g = Json::object();
    Json fs = Json::array();
    for (const auto& x : P.frontages) {
      Json ring = Json::array();
      for (const auto& q : x.footprint) ring.push(Json::array({r(q[0], 100), r(q[1], 100)}));
      fs.push(Json::object({{"parcel", x.parcel}, {"footprint", ring}, {"anchor", numbers(std::vector<double>{x.anchor[0], x.anchor[1]}, 2)}, {"centre", numbers(std::vector<double>{x.centre[0], x.centre[1]}, 2)}}));
    }
    Json segs = Json::array();
    for (const auto id : P.base.of(1)) {
      Json line = Json::array();
      for (const auto& q : P.base.line.at(id)) line.push(Json::array({r(q[0], 100), r(q[1], 100)}));
      segs.push(Json::object({{"id", std::to_string(id)}, {"type", P.base.cell(id).type}, {"name", P.base.name.count(id) ? P.base.name.at(id) : ""}, {"line", line}}));
    }
    Json ctx = Json::array();
    for (const auto& j : in["district"]["context"].as_array()) ctx.push(Json::object({{"footprint", j["footprint"]}, {"height_m", j["height_m"]}}));
    Json lk = Json::array();
    for (const auto& l : P.links) lk.push(Json::object({{"kind", l.kind}, {"a", numbers(std::vector<double>{l.a[0], l.a[1]}, 2)}, {"b", numbers(std::vector<double>{l.b[0], l.b[1]}, 2)}}));
    Json blocks = Json::array();
    for (const auto id : P.base.of(2)) {
      Json ring = Json::array();
      for (const auto& q : P.base.ring(id)) ring.push(Json::array({r(q[0], 100), r(q[1], 100)}));
      blocks.push(ring);
    }
    g.set("frontages", fs), g.set("segments", segs), g.set("context", ctx), g.set("links", lk), g.set("blocks", blocks), g.set("quad", in["district"]["redesign_quad"]);
    g.set("open_spaces", in["district"]["open_spaces"]);
    g.set("credit", "(c) OpenStreetMap contributors, ODbL 1.0");
    write_file(out + "/district-geometry.json", g.dump());
  }
  std::printf("reach: %zu pairs within %.1f m; done in %.1f s\n", reach_edges, 3 * PP.coupling, std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count());
  return 0;
}
