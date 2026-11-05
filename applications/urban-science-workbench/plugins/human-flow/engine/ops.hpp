// The operations of the plug-in "Human flow": the urban-science domain's walkers (urban::movement) on the district the
// page sends, calibrated to the City's counters; runs, their frames, comparisons, permits, POMCP's advice and the
// simulated record. Every number is computed here; the page draws.

#pragma once

#include <chrono>
#include <numbers>

#include "runs.hpp"

namespace hf {

inline Json rnd(double x, int digits = 3) {
  if (!std::isfinite(x)) return Json();
  const double p = std::pow(10.0, digits);
  return Json(std::round(x * p) / p);
}
template <class C>
inline Json nums(const C& v, int digits = 2) {
  Json a = Json::array();
  for (const auto x : v) a.push(rnd(static_cast<double>(x), digits));
  return a;
}
inline Json pt(const um::Point& p) { return Json::array({rnd(p[0], 2), rnd(p[1], 2)}); }
inline Json ints(const std::vector<std::size_t>& v) {
  Json a = Json::array();
  for (auto x : v) a.push(static_cast<double>(x));
  return a;
}

// ---- the district ---------------------------------------------------------------------------------------------------

inline Json op_prepare(const Json& r) {
  const auto D = district_of(r);
  const auto& W = *D->W;
  Json segs = Json::array();
  for (std::size_t e = 0; e < W.edges(); e++) {
    const auto& s = W.segments[e];
    Json line = Json::array();
    for (const auto& p : s.line) line.push(pt(p));
    segs.push(Json::object({{"u", static_cast<double>(s.u)}, {"v", static_cast<double>(s.v)}, {"line", line}, {"name", s.name}, {"kind", s.kind},
                            {"width", rnd(s.width, 1)}, {"length", rnd(s.length, 1)}, {"main", static_cast<bool>(D->main[s.u])}}));
  }
  Json junctions = Json::array(), signals = Json::array();
  for (std::size_t j = 0; j < W.junctions(); j++) {
    junctions.push(pt(W.at[j]));
    if (W.signalled[j]) signals.push(static_cast<double>(j));
  }
  Json streets = Json::array();
  for (const auto& s : D->streets)
    streets.push(Json::object({{"name", s.name}, {"edges", ints(s.edges)}, {"ends", ints(s.ends)}, {"length", rnd(s.length, 0)}}));
  Json counters = Json::array();
  for (std::size_t k = 0; k < D->counters.size(); k++) {
    const auto& c = D->info[k];
    counters.push(Json::object({{"id", static_cast<double>(c.id)}, {"name", c.name}, {"street", c.street}, {"at", pt(c.at)}, {"edge", static_cast<double>(c.edge)},
                                {"along", rnd(c.along, 1)}, {"toward", Json::array({c.toward[0], c.toward[1]})}, {"share", D->counters[k].share},
                                {"hand", rnd(D->counters[k].hand[0], 3)}, {"r2", rnd(D->fit.counter_r2[k], 3)},
                                {"way1", nums(c.way[0], 1)}, {"way2", nums(c.way[1], 1)}, {"fitted1", nums(D->fit.fitted[k][0], 1)}, {"fitted2", nums(D->fit.fitted[k][1], 1)}}));
  }
  Json left = Json::array();
  for (const auto& s : D->left_out) left.push(s);
  Json entrances = Json::array();
  for (const auto& E : D->entrances)
    entrances.push(Json::object({{"node", static_cast<double>(E.node)}, {"at", pt(W.at[E.node])}, {"street", E.street}, {"source", E.source},
                                 {"counter", static_cast<double>(E.counter)}, {"distance", rnd(E.distance, 0)}, {"mean", rnd(E.mean, 1)}}));
  // The trips an hour of each pattern over the week, as calibrated.
  Json patterns = Json::array();
  for (int m = 0; m < modes; m++) {
    std::array<double, H> t{};
    bool any = false;
    for (const auto& s : D->week)
      if (s.mode == m) {
        any = true;
        for (std::size_t h = 0; h < H; h++) t[h] += s.rate[h];
      }
    if (any) patterns.push(Json::object({{"mode", m}, {"name", mode_name(m)}, {"trips", nums(t, 1)}}));
  }
  double floor = 0;
  for (double f : D->floor) floor += f;
  Json kinds_j = Json::object();
  for (int k = 0; k < kinds; k++) {
    double n = 0;
    for (double a : D->attract[static_cast<std::size_t>(k)]) n += a;
    kinds_j.set(kind_name(k), n);
  }
  Json dropped = Json::array();
  for (const auto& d : D->dropped) dropped.push(d);
  return Json::object({{"key", D->key}, {"segments", segs}, {"junctions", junctions}, {"signals", signals}, {"streets", streets}, {"counters", counters},
                       {"leftOut", left}, {"entrances", entrances}, {"patterns", patterns}, {"floor", rnd(floor, 0)}, {"establishments", kinds_j},
                       {"calibration", Json::object({{"r2Way", rnd(D->fit.r2_way, 3)}, {"r2Gross", rnd(D->fit.r2_gross, 3)}, {"visitShare", rnd(D->fit.visit_share, 3)},
                                                     {"rows", static_cast<double>(D->fit.rows)}, {"theta", D->law.theta}, {"seconds", rnd(D->fit.seconds, 1)}})},
                       {"radius", rnd(D->radius, 3)}, {"dropped", dropped}, {"seconds", rnd(D->seconds, 1)}});
}

// ---- runs --------------------------------------------------------------------------------------------------------

inline Json status_of(const Run& R) {
  std::lock_guard<std::mutex> g(R.lock);
  return Json::object({{"key", R.spec.key()}, {"from", R.spec.from}, {"to", R.spec.to}, {"reached", rnd(R.reached.load(), 2)}, {"done", R.done.load()},
                       {"error", R.error}, {"quarters", static_cast<double>(R.frames.size())}, {"peak", static_cast<double>(R.peak)},
                       {"waiting", static_cast<double>(Worker::get().waiting())}, {"seconds", rnd(R.seconds, 1)}});
}

inline Json op_run(const Json& r) {
  const auto D = district_of(r);
  return status_of(*run_of(D, spec_of(*D, r)));
}

// The quarter of a run at a civil hour: the walkers every ten seconds for a minute, each [id, segment, pattern, status,
// x, y (dm), heading (°)]; the densities on the segments; the queues at their mouths.
inline Json op_frame(const Json& r) {
  const auto D = district_of(r);
  const auto R = run_of(D, spec_of(*D, r));
  const double hour = num(r, "hour", R->spec.from);
  const auto q = static_cast<long long>(std::floor((hour - R->spec.from) * 4 + 1e-9));
  std::lock_guard<std::mutex> g(R->lock);
  if (!R->error.empty()) throw std::runtime_error(R->error);
  if (q < 0 || q >= static_cast<long long>(R->quarters())) return Json::object({{"outside", true}, {"from", R->spec.from}, {"to", R->spec.to}});
  if (q >= static_cast<long long>(R->frames.size())) return Json::object({{"pending", true}, {"reached", rnd(R->reached.load(), 2)}});
  const auto& fr = R->frames[static_cast<std::size_t>(q)];
  Json frames = Json::array();
  for (const auto& f : fr) {
    Json a = Json::array();
    for (const auto& s : f) {
      a.push(static_cast<double>(s.id)), a.push(static_cast<double>(s.edge)), a.push(static_cast<double>(s.mode)), a.push(static_cast<double>(s.status));
      a.push(static_cast<double>(std::lround(s.x * 10))), a.push(static_cast<double>(std::lround(s.y * 10))), a.push(static_cast<double>(s.heading));
    }
    frames.push(a);
  }
  Json queues = Json::array();
  for (const auto& [s, n] : R->queues[static_cast<std::size_t>(q)]) queues.push(Json::array({static_cast<double>(s), static_cast<double>(n)}));
  return Json::object({{"hour", R->spec.from + 0.25 * static_cast<double>(q)}, {"step", frame_step}, {"frames", frames}, {"density", nums(R->density[static_cast<std::size_t>(q)], 3)},
                       {"queues", queues}, {"done", R->done.load()}, {"reached", rnd(R->reached.load(), 2)}});
}

struct Summary {
  double present = 0, waiting = 0, los_e = 0, los_f = 0, departed = 0, arrived = 0, failed = 0, travel = 0, delay = 0;
};
inline Summary summary_of(const um::Simulation& s) {
  const auto& t = s.totals();
  return {t.present / 3600, t.waiting / 3600, t.los_e / 3600, t.los_f / 3600, t.departed, t.arrived, t.failed, t.arrived > 0 ? t.travel / t.arrived / 60 : 0,
          t.arrived > 0 ? t.delay / t.arrived / 60 : 0};
}
inline Json summary_json(const Summary& s) {
  return Json::object({{"personHours", rnd(s.present, 1)}, {"waitingHours", rnd(s.waiting, 1)}, {"crowdedHours", rnd(s.los_e, 1)}, {"crushHours", rnd(s.los_f, 2)},
                       {"departed", s.departed}, {"arrived", s.arrived}, {"failed", s.failed}, {"travelMinutes", rnd(s.travel, 2)}, {"delayMinutes", rnd(s.delay, 3)}});
}

inline std::shared_ptr<Run> done_run(const District& D, const std::shared_ptr<const District>& Dp, const Json& r) {
  const auto R = run_of(Dp, spec_of(D, r));
  if (!R->done) return nullptr;
  if (!R->error.empty()) throw std::runtime_error(R->error);
  return R;
}

// A done run's tallies: its totals; each hour's walkers, trips by pattern, crowded walkers; each counter's reading beside
// the counts and the calibration's expectation; each segment's walkers each hour.
inline Json op_result(const Json& r) {
  const auto D = district_of(r);
  const auto R = done_run(*D, D, r);
  if (!R) return Json::object({{"pending", true}, {"status", status_of(*run_of(D, spec_of(*D, r)))}});
  const auto& sim = *R->sim;
  const auto& T = sim.tallies();
  const std::size_t E = T.edges, N = sim.streams().size(), Hn = T.hours;
  Json hours = Json::array();
  for (std::size_t h = 0; h < Hn; h++) {
    double walkers = 0, crowded = 0;
    for (std::size_t e = 0; e < E; e++) walkers += T.occupancy[h * E + e] / 3600, crowded += T.los_e[h * E + e] / 3600;
    std::array<double, modes> trips{};
    for (std::size_t s = 0; s < N; s++) trips[static_cast<std::size_t>(sim.streams()[s].mode)] += T.trips[h * N + s];
    hours.push(Json::object({{"hour", static_cast<double>(R->spec.from + static_cast<int>(h))}, {"walkers", rnd(walkers, 1)}, {"crowded", rnd(crowded, 1)}, {"trips", nums(trips, 0)}}));
  }
  Json counters = Json::array();
  const std::size_t h0 = static_cast<std::size_t>(24 * R->spec.weekday + R->spec.from);
  for (std::size_t k = 0; k < sim.counters().size(); k++) {
    Json obs = Json::array(), fit = Json::array(), simj = Json::array();
    for (int w = 0; w < 2; w++) {
      std::vector<double> o, f, s;
      for (std::size_t h = 0; h < Hn; h++) {
        o.push_back(D->info[k].way[static_cast<std::size_t>(w)][h0 + h]), f.push_back(D->fit.fitted[k][static_cast<std::size_t>(w)][h0 + h]);
        s.push_back(sim.reading(k, h)[static_cast<std::size_t>(w)]);
      }
      obs.push(nums(o, 1)), fit.push(nums(f, 1)), simj.push(nums(s, 1));
    }
    counters.push(Json::object({{"id", static_cast<double>(D->info[k].id)}, {"name", D->info[k].name}, {"observed", obs}, {"expected", fit}, {"simulated", simj}}));
  }
  Json seg = Json::array();
  for (std::size_t e = 0; e < E; e++) {
    std::vector<double> v;
    for (std::size_t h = 0; h < Hn; h++) v.push_back(T.occupancy[h * E + e] / 3600);
    seg.push(nums(v, 2));
  }
  Json present = Json::array();
  for (const auto n : R->present) present.push(static_cast<double>(n));
  return Json::object({{"status", status_of(*R)}, {"totals", summary_json(summary_of(sim))}, {"hours", hours}, {"counters", counters}, {"segments", seg},
                       {"present", present}, {"routesSeconds", rnd(R->routes_seconds, 1)}});
}

// Two done runs compared over [from, to): each segment's mean walkers (both ways), the totals, the trips.
inline Json op_compare(const Json& r) {
  const auto D = district_of(r);
  Json ra = r, rb = r;
  for (const char* k : {"weekday", "from", "to", "seed", "regime", "event"}) ra.set(k, r["a"][k]), rb.set(k, r["b"][k]);
  const auto A = done_run(*D, D, ra), B = done_run(*D, D, rb);
  if (!A || !B) return Json::object({{"pending", true}});
  if (A->spec.from != B->spec.from || A->spec.to != B->spec.to) throw std::invalid_argument("two runs of one window are compared");
  const double from = std::clamp(num(r, "over", A->spec.from + 1), static_cast<double>(A->spec.from), static_cast<double>(A->spec.to));
  const double to = std::clamp(num(r, "until", A->spec.to), from, static_cast<double>(A->spec.to));
  const auto& Ta = A->sim->tallies(), &Tb = B->sim->tallies();
  const std::size_t E = Ta.edges, h0 = static_cast<std::size_t>(from - A->spec.from), h1 = static_cast<std::size_t>(std::ceil(to - A->spec.from));
  std::vector<double> a(E, 0.0), b(E, 0.0), ea(E, 0.0), eb(E, 0.0);
  for (std::size_t h = h0; h < h1 && h < Ta.hours; h++)
    for (std::size_t e = 0; e < E; e++) {
      a[e] += Ta.occupancy[h * E + e] / 3600 / static_cast<double>(h1 - h0), b[e] += Tb.occupancy[h * E + e] / 3600 / static_cast<double>(h1 - h0);
      ea[e] += Ta.entries[h * Ta.states + 2 * e] + Ta.entries[h * Ta.states + 2 * e + 1], eb[e] += Tb.entries[h * Tb.states + 2 * e] + Tb.entries[h * Tb.states + 2 * e + 1];
    }
  // Person-hours over the hours compared, and above Fruin's levels E and F.
  const auto window = [&](const um::Tallies& T) {
    double p = 0, le = 0, lf = 0;
    for (std::size_t h = h0; h < h1 && h < T.hours; h++)
      for (std::size_t e = 0; e < E; e++) p += T.occupancy[h * E + e] / 3600, le += T.los_e[h * E + e] / 3600, lf += T.los_f[h * E + e] / 3600;
    return Json::object({{"personHours", rnd(p, 1)}, {"crowdedHours", rnd(le, 1)}, {"crushHours", rnd(lf, 2)}});
  };
  return Json::object({{"over", from}, {"until", to}, {"a", nums(a, 2)}, {"b", nums(b, 2)}, {"entriesA", nums(ea, 0)}, {"entriesB", nums(eb, 0)},
                       {"totalsA", summary_json(summary_of(*A->sim))}, {"totalsB", summary_json(summary_of(*B->sim))}, {"windowA", window(Ta)}, {"windowB", window(Tb)}});
}

// ---- the City's event permits on the district's streets ------------------------------------------------------------

// A date's counts at a counter against the same weekday in the four weeks either side (the median): the surge.
struct Surge {
  double day = 0, usual = 0;
  std::array<double, 24> observed{}, typical{};
  bool known = false;
};
inline Surge surge(const ur::Record& rec, const std::string& date) {
  Surge s;
  const long long day = static_cast<long long>(std::floor(ur::parse_utc(date + "T00:00:00Z") / ur::day_s));
  const auto hours_of = [&](long long d, std::array<double, 24>& out) {
    out.fill(std::nan(""));
    const double t = static_cast<double>(d) * ur::day_s - 11 * 3600.0;
    const long long i0 = static_cast<long long>(std::floor((t - rec.t0) / rec.step));
    for (long long i = std::max(0LL, i0 - 2); i < std::min(static_cast<long long>(rec.size()), i0 + 30); i++) {
      const double ti = rec.time(static_cast<std::size_t>(i));
      if (rec.civil_day(ti) != d) continue;
      const auto h = static_cast<std::size_t>(std::floor((rec.civil(ti) - static_cast<double>(d) * ur::day_s) / 3600 + 1e-9));
      if (h < 24) out[h] = rec.y[static_cast<std::size_t>(i)];
    }
  };
  hours_of(day, s.observed);
  std::array<std::vector<double>, 24> others;
  for (int k : {-4, -3, -2, -1, 1, 2, 3, 4}) {
    std::array<double, 24> o{};
    hours_of(day + 7 * k, o);
    for (std::size_t h = 0; h < 24; h++)
      if (std::isfinite(o[h])) others[h].push_back(o[h]);
  }
  std::size_t have = 0;
  for (std::size_t h = 0; h < 24; h++) {
    auto& v = others[h];
    std::sort(v.begin(), v.end());
    s.typical[h] = v.size() >= 3 ? v[v.size() / 2] : std::nan("");
    if (std::isfinite(s.observed[h]) && std::isfinite(s.typical[h])) s.day += s.observed[h], s.usual += s.typical[h], have++;
  }
  s.known = have >= 20;
  return s;
}

inline Json op_permits(const Json& r) {
  const auto D = district_of(r);
  std::set<std::string> names;
  for (const auto& s : D->streets) names.insert(s.name);
  // The counters on each street (inside the box, any that count), with their whole records.
  std::map<std::string, std::vector<std::pair<std::string, const Weeks*>>> on_street;
  for (const auto& s : D->placed) {
    if (!D->inside(s.at)) continue;
    const auto it = sensors().find(s.id);
    if (it == sensors().end() || it->second.name.find("Crossing") != std::string::npos) continue;
    on_street[s.way].push_back({it->second.name, &weeks(s.id)});
  }
  Json out = Json::array();
  for (const auto& p : permits()) {
    std::vector<std::string> here;
    for (const auto& s : p.streets)
      if (names.count(s)) here.push_back(s);
    if (here.empty()) continue;
    Json streets = Json::array(), counts = Json::array();
    for (const auto& s : here) {
      streets.push(s);
      for (const auto& [name, w] : on_street[s]) {
        const auto g = surge(w->record, p.from);
        if (!g.known) continue;
        counts.push(Json::object({{"counter", name}, {"street", s}, {"day", rnd(g.day, 0)}, {"usual", rnd(g.usual, 0)}, {"observed", nums(g.observed, 0)}, {"typical", nums(g.typical, 0)}}));
      }
    }
    out.push(Json::object({{"id", p.id}, {"kind", p.kind}, {"from", p.from}, {"to", p.to}, {"location", p.location}, {"streets", streets}, {"counts", counts}}));
  }
  return Json::object({{"permits", out}, {"source", "City of Melbourne, event permits 2014–2018 (CC BY 4.0), via the data library's record of the Hoddle Grid's events"}});
}

// ---- the regime advice ----------------------------------------------------------------------------------------------

// POMCP (tyche) on the walkers themselves as the generative model (movement's RegimeModel): at the event's start the
// street's regime is decided every epoch; the hidden state is the walkers' configuration and how large the crowd is
// (a scale of ½, 1 or 2 on what was declared, equally likely); the observation, the counters on or nearest the street,
// binned; the reward, minus the person-hours present and ten times those above Fruin's level F. The world is a run with
// the crowd at a true scale the request names. The planner's and the world's streams keep the pairs carrying 90% of
// each stream's weight (each particle is a copy of the world). Each fixed regime is run from the same instant with the
// same draws, for the value of hindsight.
inline Json op_advice(const Json& r) {
  static std::map<std::string, std::shared_ptr<Advice>> cache;
  static std::mutex lock;
  const auto D = district_of(r);
  Spec s = spec_of(*D, r);
  if (s.crowd <= 0 || s.event_edges.empty()) throw std::invalid_argument("advice needs an event's crowd on a street");
  if (s.regime_edges.empty()) s.regime_edges = s.event_edges, s.regime_street = s.event_street;
  const int epochs = static_cast<int>(std::clamp(num(r, "epochs", 3), 1.0, 8.0));
  const double epoch = std::clamp(num(r, "epochMinutes", 10), 2.0, 30.0) * 60;
  const int simulations = static_cast<int>(std::clamp(num(r, "simulations", 12), 2.0, 200.0));
  const double truth = std::clamp(num(r, "truth", 1.5), 0.25, 4.0);
  const std::string key = D->key + "|" + s.key() + "|" + std::to_string(epochs) + "/" + std::to_string(epoch) + "/" + std::to_string(simulations) + "/" + std::to_string(truth);
  std::shared_ptr<Advice> A;
  {
    std::lock_guard<std::mutex> g(lock);
    if (auto it = cache.find(key); it != cache.end()) A = it->second;
    else {
      A = std::make_shared<Advice>();
      cache[key] = A;
      Worker::get().submit([A, D, s, epochs, epoch, simulations, truth](calx::scheduler::ThreadPool& pool) {
        try {
          const auto t0 = std::chrono::steady_clock::now();
          auto targets = std::make_shared<um::Targets>(*D->T);
          auto streams = streams_of(*D, s, &targets);
          const std::size_t crowd_stream = streams.size() - 1;
          std::shared_ptr<const um::Values> V = D->V;
          if (targets->size() != D->T->size()) V = std::make_shared<const um::Values>(D->R, *targets, &pool);
          um::Street street{s.regime_street, s.regime_edges, s.entry == none ? D->W->segments[s.regime_edges.front()].u : s.entry};
          if (street.entry == none || std::find_if(s.regime_edges.begin(), s.regime_edges.end(), [&](std::size_t e) { return D->W->segments[e].u == street.entry || D->W->segments[e].v == street.entry; }) == s.regime_edges.end()) street.entry = D->W->segments[s.regime_edges.front()].u;
          {
            std::map<std::size_t, int> degree;
            for (auto e : s.regime_edges) degree[D->W->segments[e].u]++, degree[D->W->segments[e].v]++;
            for (const auto& [j, d] : degree)
              if (d == 1) { street.entry = j; break; }
          }
          auto RR = std::make_shared<const um::RegimeRoutes>(um::regime_routes(D->W, D->law, street, *targets, &pool));
          streams = pruned(split_for_regime(std::move(streams), *RR), 0.9);
          um::Options o;
          o.seed = s.seed, o.hours = static_cast<std::size_t>(s.to - s.from), o.light = true;
          um::Simulation world(std::shared_ptr<const um::Targets>(targets), {RR->values[0][0], RR->values[0][1]}, streams, D->counters, o);
          world.set_scale(crowd_stream, truth);
          const double start = (std::max<double>(s.crowd_from, s.from + 1) - s.from) * 3600;
          world.run(start, &pool);
          // The sensed counters: those on the street, else the two nearest it.
          std::vector<std::size_t> sensed;
          std::set<std::size_t> on(s.regime_edges.begin(), s.regime_edges.end());
          for (std::size_t k = 0; k < D->counters.size(); k++)
            if (on.count(D->counters[k].edge)) sensed.push_back(k);
          if (sensed.empty()) {
            const auto& g = D->W->segments[s.regime_edges[s.regime_edges.size() / 2]];
            const um::Point mid = D->W->point(2 * (s.regime_edges[s.regime_edges.size() / 2]), g.length / 2);
            std::vector<std::pair<double, std::size_t>> d;
            for (std::size_t k = 0; k < D->counters.size(); k++) d.push_back({std::hypot(D->info[k].at[0] - mid[0], D->info[k].at[1] - mid[1]), k});
            std::sort(d.begin(), d.end());
            for (std::size_t q = 0; q < std::min<std::size_t>(2, d.size()); q++) sensed.push_back(d[q].second);
          }
          if (sensed.size() > 2) sensed.resize(2);
          // Bins of a sensed counter's reading in an epoch about what the calibration expects then.
          const std::size_t hw = static_cast<std::size_t>(24 * s.weekday) + static_cast<std::size_t>(std::floor(s.from + start / 3600));
          double expect = 0;
          for (const auto k : sensed) expect += (D->fit.fitted[k][0][std::min<std::size_t>(hw, H - 1)] + D->fit.fitted[k][1][std::min<std::size_t>(hw, H - 1)]) * epoch / 3600 / static_cast<double>(sensed.size());
          um::RegimeModel::Setup setup;
          setup.routes = RR, setup.epoch = epoch, setup.epochs = static_cast<std::size_t>(epochs), setup.sensed = sensed;
          setup.bins = {0.9 * expect, 1.3 * expect}, setup.scaled = {crowd_stream}, setup.scales = {0.5, 1.0, 2.0}, setup.prior = {1.0 / 3, 1.0 / 3, 1.0 / 3};
          setup.discount = 0.95;
          const um::Simulation at_start = world;
          um::RegimeModel model(at_start, setup);
          const double span = std::max(10.0, static_cast<double>(world.active()) * epoch / 3600 * 1.5);
          tyche::decision::Pomcp<um::RegimeModel> planner(model, {.simulations = simulations, .reward_span = span, .depth = std::min(epochs, 3), .particles = static_cast<std::size_t>(std::max(6, simulations / 2)),
                                                                  .refill_attempts = simulations, .seed = s.seed + 7});
          planner.set_leaf([](const um::RegimeModel::State&) { return 0.0; });
          Json steps = Json::array();
          double total = 0;
          for (int e = 0; e < epochs; e++) {
            const auto p0 = std::chrono::steady_clock::now();
            const std::size_t a = planner.plan();
            const double planning = std::chrono::duration<double>(std::chrono::steady_clock::now() - p0).count();
            const auto [obs, reward] = model.act(world, a);
            total += reward;
            // The belief over the crowd's scale: the particles' scales.
            std::array<double, 3> belief{};
            for (const auto& x : planner.belief()) {
              const double sc = x.sim.scale(crowd_stream);
              belief[sc < 0.75 ? 0 : sc < 1.5 ? 1 : 2] += 1.0 / static_cast<double>(planner.belief().size());
            }
            const bool kept = planner.advance(a, obs);
            steps.push(Json::object({{"hour", rnd(s.from + (start + e * epoch) / 3600, 3)}, {"action", um::regime_name(static_cast<um::Regime>(a))}, {"observation", static_cast<double>(obs)},
                                     {"reward", rnd(reward, 1)}, {"planningSeconds", rnd(planning, 2)}, {"belief", nums(belief, 3)}, {"particles", kept}}));
            A->epoch = e + 1;
            if (!kept) break;
          }
          // Hindsight: each fixed regime from the same instant, with the same draws.
          Json fixed = Json::array();
          for (int a = 0; a < 3; a++) {
            um::Simulation w = at_start;
            double g = 0;
            for (int e = 0; e < epochs; e++) g += model.act(w, static_cast<std::size_t>(a)).second;
            fixed.push(Json::object({{"regime", um::regime_name(static_cast<um::Regime>(a))}, {"reward", rnd(g, 1)}}));
          }
          Json sensed_j = Json::array();
          for (auto k : sensed) sensed_j.push(D->info[k].name);
          std::lock_guard<std::mutex> g(A->lock);
          A->result = Json::object({{"steps", steps}, {"reward", rnd(total, 1)}, {"fixed", fixed}, {"sensed", sensed_j}, {"bins", nums(setup.bins, 0)}, {"epochMinutes", epoch / 60},
                                    {"simulations", simulations}, {"truth", truth}, {"street", s.regime_street}, {"seconds", rnd(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), 1)}});
        } catch (const std::exception& e) {
          A->error = e.what();
        }
        A->done = true;
      });
    }
  }
  std::lock_guard<std::mutex> g(A->lock);
  if (!A->done) return Json::object({{"pending", true}, {"epoch", A->epoch.load()}});
  if (!A->error.empty()) throw std::runtime_error(A->error);
  return A->result;
}

// ---- the run's record, simulated --------------------------------------------------------------------------------

inline Json op_record(const Json& r) {
  const auto D = district_of(r);
  const auto R = done_run(*D, D, r);
  if (!R) return Json::object({{"pending", true}});
  std::vector<um::RegimeSpan> spans;
  if (R->spec.regime && !R->spec.regime_edges.empty())
    spans.push_back({{R->spec.regime_street, R->spec.regime_edges, R->spec.entry}, static_cast<um::Regime>(R->spec.regime), (std::max<double>(R->spec.from, R->spec.regime_from) - R->spec.from) * 3600,
                     (std::min<double>(R->spec.to, R->spec.regime_to) - R->spec.from) * 3600});
  um::RecordOptions o;
  o.mjd0 = num(r, "mjd0", 0), o.counts = false;
  const auto rec = um::simulated_record(*R->sim, spans, o);
  const auto problems = urban::events::problems(rec);
  Json events = Json::array();
  std::map<std::string, double> kinds;
  bool simulated = true;
  const std::size_t E = D->W->edges();
  for (const auto& e : rec.events) {
    kinds[e.type + "/" + e.kind] += 1;
    simulated = simulated && e.provenance == "simulated";
    const std::size_t c = e.couplings.empty() ? 0 : e.couplings.front();
    const auto& cp = rec.couplings[c];
    events.push(Json::object({{"id", e.id}, {"type", e.type}, {"kind", e.kind}, {"from", rnd(R->spec.from + (e.when.from - o.mjd0) * 24, 3)}, {"to", rnd(R->spec.from + (e.when.to - o.mjd0) * 24, 3)},
                              {"where", c < E ? cp.name : std::string("a junction")}, {"at", pt({cp.at[0], cp.at[1]})}, {"source", e.source}, {"provenance", e.provenance}}));
  }
  Json k = Json::object();
  for (const auto& [n, v] : kinds) k.set(n, v);
  Json p = Json::array();
  for (const auto& x : problems) p.push(x);
  return Json::object({{"events", events}, {"kinds", k}, {"simulated", simulated}, {"problems", p}, {"couplings", static_cast<double>(rec.couplings.size())},
                       {"primitives", static_cast<double>(rec.primitives.size())}});
}

}  // namespace hf
