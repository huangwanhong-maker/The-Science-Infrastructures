// Runs of the district's walkers: a window of civil hours on a day of the week, as built or with a street's regime
// changed (closed to walkers passing through, or one-way: phased as crowds are run at events) and with an event's crowd
// bound for a street, each a Simulation of urban::movement (synchronous, the same on any number of threads; common
// random numbers, so that two runs differ only by what was changed).
//
// A run is a job: it is queued and run on the engine's worker thread, and the page reads it as it goes. While it runs it
// keeps, at every quarter hour, a minute of the walkers' places every ten seconds (the page moves each walker between
// them), the density on each segment and the queues at the segments' mouths; when it is done, its tallies.
//
// The first hour of a window is its warm-up: the network starts empty, and walks last minutes.

#pragma once

#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <numbers>
#include <numeric>
#include <thread>

#include <tyche/decision/online.hpp>
#include <urban/events/record.hpp>

#include "district.hpp"

namespace hf {

// What a run is: its day and window, a regime on a street over hours, an event's crowd on a street over hours.
struct Spec {
  int weekday = 2;                 // 0 Monday … 6 Sunday
  int from = 0, to = 24;           // civil hours of the window
  std::uint64_t seed = 20261105;
  std::string regime_street;       // a named street (or empty)
  std::vector<std::size_t> regime_edges;   // its segments (the whole named street, or a block picked)
  std::size_t entry = none;        // where a one-way street is entered from
  int regime = 0;                  // 0 open, 1 one-way, 2 closed to walkers passing through
  double regime_from = 0, regime_to = 24;
  std::string event_street;
  std::vector<std::size_t> event_edges;
  double crowd = 0, crowd_from = 0, crowd_to = 0, dwell = 2700;   // walkers an hour; civil hours; s
  std::string key() const {
    std::string k = std::to_string(weekday) + "/" + std::to_string(from) + "-" + std::to_string(to) + "/" + std::to_string(seed);
    if (!regime_edges.empty()) {
      k += "/r" + std::to_string(regime) + ":" + std::to_string(entry) + ":" + std::to_string(regime_from) + "-" + std::to_string(regime_to) + ":";
      for (auto e : regime_edges) k += std::to_string(e) + ",";
    }
    if (crowd > 0 && !event_edges.empty()) {
      k += "/c" + std::to_string(crowd) + ":" + std::to_string(crowd_from) + "-" + std::to_string(crowd_to) + ":" + std::to_string(dwell) + ":";
      for (auto e : event_edges) k += std::to_string(e) + ",";
    }
    return k;
  }
};

inline std::vector<std::size_t> edges_from(const Json& j) {
  std::vector<std::size_t> out;
  if (j.is_array())
    for (const auto& e : j.as_array()) out.push_back(static_cast<std::size_t>(e.as_number()));
  return out;
}

// A spec from a request: {weekday, from, to, seed, regime: {street | edges, entry, kind, from, to}, event: {street | edges,
// crowd, from, to, dwell}}.
inline Spec spec_of(const District& D, const Json& r) {
  Spec s;
  s.weekday = static_cast<int>(std::clamp(num(r, "weekday", 2), 0.0, 6.0));
  s.from = static_cast<int>(std::clamp(std::floor(num(r, "from", 0)), 0.0, 23.0));
  s.to = static_cast<int>(std::clamp(std::ceil(num(r, "to", 24)), static_cast<double>(s.from + 1), 24.0));
  s.seed = static_cast<std::uint64_t>(num(r, "seed", 20261105));
  const auto street_edges = [&](const Json& j, std::string& name) {
    std::vector<std::size_t> e = edges_from(j["edges"]);
    name = j["street"].text();
    if (e.empty() && !name.empty())
      for (const auto& st : D.streets)
        if (st.name == name) e = st.edges;
    for (const auto x : e)
      if (x >= D.W->edges()) throw std::invalid_argument("no segment " + std::to_string(x) + " in the district");
    std::sort(e.begin(), e.end());
    return e;
  };
  if (const Json& g = r["regime"]; g.is_object()) {
    const std::string kind = g["kind"].text("open");
    s.regime = kind == "closed" ? 2 : kind == "one-way" ? 1 : 0;
    s.regime_edges = street_edges(g, s.regime_street);
    // A named street's stubs where it leaves the district stay open: walkers enter the district there.
    if (!edges_from(g["edges"]).size() && s.regime_edges.size() > 1)
      s.regime_edges.erase(std::remove_if(s.regime_edges.begin(), s.regime_edges.end(), [&](std::size_t e) {
        return D.W->leaving[D.W->segments[e].u].size() == 1 || D.W->leaving[D.W->segments[e].v].size() == 1;
      }), s.regime_edges.end());
    s.regime_from = std::clamp(num(g, "from", s.from), 0.0, 24.0), s.regime_to = std::clamp(num(g, "to", s.to), 0.0, 24.0);
    if (s.regime && s.regime_edges.empty()) throw std::invalid_argument("a regime on no street");
    if (s.regime == 1) {
      s.entry = g["entry"].is_number() ? static_cast<std::size_t>(g["entry"].as_number()) : none;
      if (s.entry == none) {
        // The street's end farthest south, then west, if none is named.
        std::map<std::size_t, int> degree;
        for (auto e : s.regime_edges) degree[D.W->segments[e].u]++, degree[D.W->segments[e].v]++;
        for (const auto& [j, d] : degree)
          if (d == 1 && (s.entry == none || D.W->at[j][1] < D.W->at[s.entry][1] - 1 || (std::fabs(D.W->at[j][1] - D.W->at[s.entry][1]) <= 1 && D.W->at[j][0] < D.W->at[s.entry][0]))) s.entry = j;
        if (s.entry == none) s.entry = D.W->segments[s.regime_edges.front()].u;
      }
    }
  }
  if (const Json& c = r["event"]; c.is_object()) {
    s.event_edges = street_edges(c, s.event_street);
    s.crowd = std::clamp(num(c, "crowd", 0), 0.0, 40000.0);
    s.crowd_from = std::clamp(num(c, "from", s.from + 1), 0.0, 24.0), s.crowd_to = std::clamp(num(c, "to", s.to), 0.0, 24.0);
    s.dwell = std::clamp(num(c, "dwell", 45), 0.0, 240.0) * 60;
    if (s.crowd > 0 && s.event_edges.empty()) throw std::invalid_argument("an event on no street");
  }
  return s;
}

// The streams of a spec's window: the district's week sliced to the day's hours, and the event's crowd: from the visits'
// origins (half the entrances, half the buildings) to the street, dwelling there, then back where it came from (route
// class 1: under a closure it may still enter the street it is bound for).
inline std::vector<um::Stream> streams_of(const District& D, const Spec& s, std::shared_ptr<um::Targets>* targets = nullptr) {
  std::vector<um::Stream> out = D.week;
  const std::size_t h0 = static_cast<std::size_t>(24 * s.weekday + s.from), n = static_cast<std::size_t>(s.to - s.from);
  for (auto& st : out) st.rate = std::vector<double>(st.rate.begin() + static_cast<std::ptrdiff_t>(h0), st.rate.begin() + static_cast<std::ptrdiff_t>(h0 + n));
  if (s.crowd > 0 && !s.event_edges.empty()) {
    std::size_t target = none;
    for (const auto& st : D.streets)
      if (st.edges == s.event_edges) target = st.target;
    if (target == none) {
      if (!targets) throw std::invalid_argument("the event's street is not a target of the district");
      target = (*targets)->segments(s.event_edges);
    }
    um::Stream c;
    c.name = "the crowd for " + (s.event_street.empty() ? std::string("the street") : s.event_street), c.mode = crowd, c.route_class = 1, c.dwell = s.dwell;
    for (std::size_t q = 0; q < D.W->junctions(); q++)
      if (D.origin[q] > 0 && D.junction_target[q] != none) c.origin.push_back(q), c.target.push_back(target), c.back.push_back(D.junction_target[q]), c.weight.push_back(D.origin[q]);
    c.rate.assign(n, 0.0);
    for (std::size_t h = 0; h < n; h++) {
      const double a = static_cast<double>(s.from) + static_cast<double>(h), b = a + 1;
      c.rate[h] = s.crowd * std::max(0.0, std::min(b, s.crowd_to) - std::max(a, s.crowd_from));
    }
    out.push_back(std::move(c));
  }
  return out;
}

// A street's regime keeps its own walkers: the pairs of a stream (of route class 0) that have no way to their target with
// the street closed to walkers passing through (they begin or end on it, or on a lane reached only through it) are
// split off into a stream of route class 1, which may enter the street under every regime; the rest stay in class 0.
// Poisson splitting keeps the law of the whole: each part's rate is the stream's times its share of the weight. The
// parts split off are appended, so that the streams before them keep their draws (common random numbers between a
// street open and changed). The case study does the same (research/experiments/2026-11-05-walkers-of-the-hoddle-grid).
inline std::vector<um::Stream> split_for_regime(std::vector<um::Stream> streams, const um::RegimeRoutes& RR) {
  const um::Values& closed = *RR.values[2][0];
  const um::Routes& R = closed.routes();
  const std::size_t n = streams.size();
  for (std::size_t k = 0; k < n; k++) {
    if (streams[k].route_class != 0) continue;
    um::Stream keep = streams[k], away = streams[k];
    keep.origin.clear(), keep.target.clear(), keep.back.clear(), keep.weight.clear();
    away.origin.clear(), away.target.clear(), away.back.clear(), away.weight.clear();
    double wk = 0, wa = 0;
    for (std::size_t p = 0; p < streams[k].origin.size(); p++) {
      const bool reach = R.total(streams[k].origin[p], closed[streams[k].target[p]]) > 0;
      um::Stream& to = reach ? keep : away;
      to.origin.push_back(streams[k].origin[p]), to.target.push_back(streams[k].target[p]), to.weight.push_back(streams[k].weight[p]);
      if (!streams[k].back.empty()) to.back.push_back(streams[k].back[p]);
      (reach ? wk : wa) += streams[k].weight[p];
    }
    if (wa <= 0) continue;
    const double total = wk + wa;
    for (auto& r : keep.rate) r *= wk / total;
    for (auto& r : away.rate) r *= wa / total;
    away.route_class = 1, away.name += " (its own street)";
    if (wk > 0) streams[k] = std::move(keep);
    else { streams[k].route_class = 1; continue; }
    streams.push_back(std::move(away));
  }
  return streams;
}

struct Snap {
  std::uint32_t id;
  std::uint16_t edge;
  std::uint8_t mode, status;
  float x, y;
  std::int16_t heading;   // degrees, counter-clockwise from east
};

inline constexpr int frames_per_quarter = 7;   // a minute, every ten seconds
inline constexpr double frame_step = 10;

class Run {
 public:
  Spec spec;
  std::string key, label;
  std::shared_ptr<const District> D;
  std::atomic<double> reached{-1};   // civil hours run
  std::atomic<bool> done{false};
  std::string error;
  double seconds = 0, routes_seconds = 0;
  std::size_t peak = 0;
  // Per quarter of the window: its frames, the densities on the segments (m⁻², the frames' mean), the queues.
  std::vector<std::vector<std::vector<Snap>>> frames;
  std::vector<std::vector<float>> density;
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> queues;   // (state, walkers waiting)
  std::vector<std::size_t> present;   // walkers at each quarter
  std::unique_ptr<um::Simulation> sim;
  std::shared_ptr<const um::RegimeRoutes> RR;
  std::shared_ptr<um::Targets> targets;
  mutable std::mutex lock;

  std::size_t quarters() const { return static_cast<std::size_t>(4 * (spec.to - spec.from)); }

  void execute(calx::scheduler::ThreadPool& pool) {
    const auto t0 = std::chrono::steady_clock::now();
    targets = std::make_shared<um::Targets>(*D->T);
    auto streams = streams_of(*D, spec, &targets);
    std::shared_ptr<const um::Values> V = D->V;
    if (targets->size() != D->T->size()) V = std::make_shared<const um::Values>(D->R, *targets, &pool);
    std::vector<std::shared_ptr<const um::Values>> routes{V, V};
    if (!spec.regime_edges.empty()) {
      um::Street street{spec.regime_street, spec.regime_edges, spec.entry == none ? D->W->segments[spec.regime_edges.front()].u : spec.entry};
      RR = std::make_shared<const um::RegimeRoutes>(um::regime_routes(D->W, D->law, street, *targets, &pool));
      routes = {RR->values[0][0], RR->values[0][1]};
      streams = split_for_regime(std::move(streams), *RR);
    }
    routes_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    um::Options o;
    o.seed = spec.seed, o.hours = static_cast<std::size_t>(spec.to - spec.from), o.minute_bins = 60;
    sim = std::make_unique<um::Simulation>(std::shared_ptr<const um::Targets>(targets), routes, streams, D->counters, o);
    const auto at = [&](double civil) { return (civil - spec.from) * 3600.0; };
    std::vector<std::pair<double, int>> switches;
    if (RR && spec.regime) switches = {{at(std::max<double>(spec.from, spec.regime_from)), spec.regime}, {at(std::min<double>(spec.to, spec.regime_to)), 0}};
    std::size_t next_switch = 0;
    const auto run_to = [&](double t) {
      while (next_switch < switches.size() && switches[next_switch].first <= t) {
        sim->run(switches[next_switch].first, &pool);
        um::apply(*sim, *RR, static_cast<um::Regime>(switches[next_switch].second));
        next_switch++;
      }
      sim->run(t, &pool);
    };
    const std::size_t Q = quarters(), E = D->W->edges();
    for (std::size_t q = 0; q < Q; q++) {
      std::vector<std::vector<Snap>> fr;
      std::vector<float> dens(E, 0.0f);
      std::vector<std::pair<std::uint32_t, std::uint32_t>> qs;
      for (int k = 0; k < frames_per_quarter; k++) {
        run_to(900.0 * static_cast<double>(q) + frame_step * k);
        std::vector<Snap> snap;
        for (const auto& p : sim->positions()) {
          Snap s{p.id, static_cast<std::uint16_t>(p.state / 2), static_cast<std::uint8_t>(sim->streams()[p.stream].mode), static_cast<std::uint8_t>(p.status),
                 static_cast<float>(p.at[0]), static_cast<float>(p.at[1]), static_cast<std::int16_t>(std::lround(std::atan2(p.heading[1], p.heading[0]) * 180 / std::numbers::pi))};
          snap.push_back(s);
        }
        const auto& d = sim->density();
        for (std::size_t e = 0; e < E; e++) dens[e] += static_cast<float>(d[e] / frames_per_quarter);
        if (k == 0) {
          const auto w = sim->queues();
          for (std::size_t s = 0; s < w.size(); s++)
            if (w[s] > 0) qs.push_back({static_cast<std::uint32_t>(s), static_cast<std::uint32_t>(w[s])});
        }
        peak = std::max(peak, snap.size());
        fr.push_back(std::move(snap));
      }
      {
        std::lock_guard<std::mutex> g(lock);
        present.push_back(fr.front().size());
        frames.push_back(std::move(fr)), density.push_back(std::move(dens)), queues.push_back(std::move(qs));
      }
      reached = spec.from + 0.25 * static_cast<double>(q + 1);
    }
    run_to(at(spec.to));
    reached = spec.to;
    seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
};

// The engine's one worker: jobs (runs, advice) one at a time, on its own pool.
class Worker {
 public:
  static Worker& get() {
    static Worker w;
    return w;
  }
  void submit(std::function<void(calx::scheduler::ThreadPool&)> job) {
    {
      std::lock_guard<std::mutex> g(m_);
      jobs_.push_back(std::move(job));
    }
    cv_.notify_one();
  }
  std::size_t waiting() const {
    std::lock_guard<std::mutex> g(m_);
    return jobs_.size();
  }

 private:
  Worker() : pool_(std::clamp<std::size_t>(std::thread::hardware_concurrency(), 1, 6)), thread_([this] { loop(); }) {}
  ~Worker() { thread_.detach(); }
  void loop() {
    while (true) {
      std::function<void(calx::scheduler::ThreadPool&)> job;
      {
        std::unique_lock<std::mutex> g(m_);
        cv_.wait(g, [&] { return !jobs_.empty(); });
        job = std::move(jobs_.front());
        jobs_.pop_front();
      }
      job(pool_);
    }
  }
  calx::scheduler::ThreadPool pool_;
  mutable std::mutex m_;
  std::condition_variable cv_;
  std::deque<std::function<void(calx::scheduler::ThreadPool&)>> jobs_;
  std::thread thread_;
};

// The runs, by their district and spec (the last few kept).
inline std::shared_ptr<Run> run_of(std::shared_ptr<const District> D, const Spec& s, const std::string& label = "") {
  static std::map<std::string, std::shared_ptr<Run>> cache;
  static std::vector<std::string> order;
  static std::mutex lock;
  const std::string key = D->key + "|" + s.key();
  std::lock_guard<std::mutex> g(lock);
  if (auto it = cache.find(key); it != cache.end()) return it->second;
  auto R = std::make_shared<Run>();
  R->spec = s, R->key = key, R->D = std::move(D), R->label = label;
  cache[key] = R, order.push_back(key);
  // Keep six runs, never dropping one not yet done.
  while (order.size() > 6) {
    auto it = std::find_if(order.begin(), order.end(), [&](const std::string& k) { return cache[k]->done.load(); });
    if (it == order.end()) break;
    cache.erase(*it), order.erase(it);
  }
  Worker::get().submit([R](calx::scheduler::ThreadPool& pool) {
    try {
      R->execute(pool);
    } catch (const std::exception& e) {
      R->error = e.what();
    }
    R->done = true;
  });
  return R;
}

// ---- the regime advice: POMCP on the walkers ----------------------------------------------------------------------

struct Advice {
  std::atomic<bool> done{false};
  std::atomic<int> epoch{0};
  std::string error;
  Json result;
  mutable std::mutex lock;
};

// A stream's pairs pruned to those carrying `keep` of its weight (the planner copies its world for every particle).
inline std::vector<um::Stream> pruned(std::vector<um::Stream> streams, double keep) {
  for (auto& s : streams) {
    if (s.origin.size() < 64) continue;
    std::vector<std::size_t> idx(s.origin.size());
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return s.weight[a] > s.weight[b]; });
    double total = 0, run = 0;
    for (double w : s.weight) total += w;
    um::Stream t = s;
    t.origin.clear(), t.target.clear(), t.back.clear(), t.weight.clear();
    for (const std::size_t i : idx) {
      t.origin.push_back(s.origin[i]), t.target.push_back(s.target[i]), t.weight.push_back(s.weight[i]);
      if (!s.back.empty()) t.back.push_back(s.back[i]);
      run += s.weight[i];
      if (run >= keep * total) break;
    }
    s = std::move(t);
  }
  return streams;
}

}  // namespace hf
