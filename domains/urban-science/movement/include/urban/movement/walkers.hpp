// Walkers as agents (theory/movement.md §4): each walker an agent of nomos's engine, with an origin, a target and a way
// chosen as it goes by the logit over walks (route.hpp), walking the network at its own speed, slowed by the crowd,
// held at red signals and at full streets.
//
// THE DEMAND. Trips come in streams, each a law over pairs (an origin junction, a target) and a rate, trips an hour,
// for every hour of the run: the patterns of the district-economy flow study (an amplitude times a fixed law of origins
// and destinations), or a crowd bound for a street. Departures are a Poisson process whose rate is constant through each
// hour: in each minute of the run a stream sends Poisson(rate/60) walkers, their pairs drawn from its law by Walker's
// alias table (calx), their instants uniform in the minute. A stream may make its walkers dwell at the target (an
// exponential time of given mean, standing on the segment where they stop) and then walk back to their origin.
//
// THE STEP (Δt, 2 s by default). A walker on a segment walks at its free speed v₀ (drawn once, normal of mean 1.34 and
// standard deviation 0.26 m/s: Weidmann 1993, Transporttechnik der Fussgänger, ETH Zürich, IVT 90, §4.1), times Tobler's
// factor of the grade (network.hpp), times Greenshields's 1 − ρ/ρ_j, ρ the density of the others (both ways) in the
// space ahead of it (the next 10 m of its segment, or what is left of it, times the walkable width) at the
// step's start, and ρ_j = 5.4 m⁻² (Weidmann's jam density, as urban::flow's kinematic street takes it). In a crowd it
// still creeps at 0.05 m/s; it stands when the walkers going its way ahead of it are at 0.95 ρ_j. Walkers dwelling
// stand aside (at a stall, a shop window): they count in the density that slows the others, not among the walkers
// going their way, so a crowd that dwells thins a street's stream without stopping it. So a queue stands in
// space behind a bottleneck, and its head, with room ahead, walks off: the local reading of the kinematic wave, whose
// cells would discharge a queue at capacity (Lighthill and Whitham 1955; Daganzo 1994). At the end of a segment it
// reaches its target or chooses the next segment by the route law. Crossing a carriageway at a signalled junction it
// waits for its phase: two phases, along the junction's axis and across it, each green for its share of the cycle (120 s,
// half and half, as urban::flow::StreetDay's crossings; the junctions' offsets are fixed and arbitrary, as their timings
// are not in the data). It then asks to enter the segment. A segment admits its supply, Godunov's receiving flow
// (Daganzo 1994, "The cell transmission model", Transp. Res. B 28, 269–287): its capacity v₀ ρ_j w / 4 walkers a second
// (Greenshields's greatest flow, at ρ_j / 2, times the width w) while the walkers going its way in the first 10 m past
// its mouth are no denser than ρ_j / 2, else the flow their density carries, and none at 0.95 ρ_j (or with everyone
// there, dwellers too, at ρ_j); the walkers it does not admit wait at its mouth, standing as a queue behind it at a
// queue's density (2.5 m⁻², an assumption), so that a queue at red does not read as an infinite density.
//
// THE SCHEDULE. Synchronous (nomos's Schedule::synchronous): every walker reads the configuration at the step's start
// and writes only its own part, so that the step runs on a pool. Admission is then decided segment by segment, the
// requests ranked by the instant asked and the walker's identity. Every draw is keyed by what it is for and by the
// walker's identity, not its place in the population (seed, "walk", identity, step), so that a run is the same on any
// number of threads, and a walker draws the same in two runs that differ elsewhere (common random numbers for comparing
// regimes). Walkers are born at departure and die on arriving; identities are their order of departure.
//
// WHAT A RUN KEEPS. The entries into every directed segment by the hour of entry and by the hour of departure (the latter
// is exactly the flow's expectation, route.hpp); each counter's passes each way by the hour; each segment's
// person-seconds by the hour (its mean walkers and density), those of walkers walking or dwelling about whom (5 m
// either way) the density is at Fruin's levels of service E and F (above 1.08 and 2.0 m⁻²: Fruin 1971, Pedestrian Planning and Design,
// Metropolitan Association of Urban Designers and Environmental Planners, New York, ch. 4), and the queues at each
// segment's mouth; optionally each segment's person-seconds, all and at level E, minute by minute; the trips made, with their times and their delays over walking the same way alone at v₀; and
// running totals for the decision (regime.hpp). The walkers' places at any instant are read by positions().

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <nomos/agents/engine.hpp>

#include <urban/movement/network.hpp>
#include <urban/movement/route.hpp>

namespace urban::movement {

struct Stream {
  std::string name;
  int mode = 0;                          // the caller's label (a pattern of trips)
  int route_class = 0;                   // the routes its walkers take (0: everyone's; a regime may treat a class apart)
  std::vector<std::size_t> origin, target, back;   // per pair: a junction, a target id, and the target of the walk back
  std::vector<double> weight;            // the law over pairs
  std::vector<double> rate;              // trips an hour, for each hour of the run
  double dwell = 0;                      // s, the mean dwell at the target (0: none, the trip ends there)
};

struct Kinematics {
  double step = 2;                    // s
  double free_speed = 1.34, speed_sd = 0.26;   // m/s (Weidmann 1993)
  double min_speed = 0.05;            // m/s: a crowd's creep
  double window = 10;                 // m: the space ahead whose density slows a walker (an assumption)
  double queue_density = 2.5;         // m⁻²: walkers standing in a queue (about 0.4 m² each; an assumption, after Fruin's
                                      // levels of service for queues, 0.3 to 0.65 m² a person at levels D and E)
  double jam = 5.4;                   // m⁻² (Weidmann 1993)
  bool capacity = true;               // entries limited to the capacity, none above 0.95 of the jam
  bool signals = true;
  double cycle = 120, green = 0.5;    // s; the green share of phase 0 (phase 1 has the rest)
  bool grade = true;                  // Tobler's factor on the speed
  double los_e = 1.08, los_f = 2.0;   // m⁻² (Fruin 1971)
};

struct Options {
  std::uint64_t seed = 20261105;
  Kinematics kinematics;
  std::size_t hours = 24;          // the hours the tallies keep (walkers still walking after them are kept, not tallied)
  double minute_bins = 0;          // s; > 0: each segment's density kept in bins of this length
  bool light = false;              // true: running totals and counters only (for a planner's copies)
};

// A counting line on a segment: walkers passing `along` metres from the segment's tail, each way. A footpath's counter
// sees a share of the street's walkers (½ on a street of two footpaths; the flow study's σ), and of the two ways (1 ± ρ)
// of that share (its handedness): what it reads is share × (1 ± ρ) × the passes, in expectation.
struct Counter {
  std::string name;
  std::size_t edge = 0;
  double along = 0;
  double share = 0.5;
  std::array<double, 2> hand{0, 0};          // ρ for its first and second ways: factors 1 + hand[w]
  std::array<std::size_t, 2> way{0, 1};      // the directed states of its first and second ways (2e, 2e + 1 in some order)
};

enum class Status : std::uint8_t { pending, walking, asking, signal, dwelling, done };

struct Walker {
  std::uint32_t id = 0, stream = 0, pair = 0;
  std::uint32_t state = static_cast<std::uint32_t>(none);   // the segment it is on
  std::uint32_t next = static_cast<std::uint32_t>(none);    // the segment it asks to enter
  std::uint32_t target = 0, back = 0, origin = 0;
  std::uint16_t hour = 0;           // of departure
  Status status = Status::pending;
  std::uint8_t leg = 0;             // 0 out, 1 back
  std::uint8_t cls = 0;
  std::uint8_t dwelt = 0;           // the dwell is over
  std::uint8_t failed = 0;          // no way to its target
  float speed = 1.34f;              // free speed, m/s
  double s = 0, s0 = 0;             // m along `state`, now and at the step's start
  double stop = -1;                 // m along a target segment where it will dwell
  double depart = 0, until = 0, ask = 0, end = 0;   // s
  double walked = 0, alone = 0;     // m walked; s it would have taken alone at v₀
  double rest = 0;                  // s dwelt
};

// Running totals since the run began (a planner differences two of them for an epoch's reward).
struct Totals {
  double present = 0;         // person-seconds in the network (walking, waiting, dwelling)
  double waiting = 0;         // person-seconds queued or at red
  double los_e = 0, los_f = 0;    // person-seconds on segments above Fruin's E and F
  double departed = 0, arrived = 0, failed = 0;
  double travel = 0, delay = 0;   // s over completed trips: their times, and those times less walking alone
};

struct Tallies {
  std::size_t hours = 0, states = 0, edges = 0;
  std::vector<double> entries, entries_by_departure;   // [hour * S + state]
  std::vector<double> occupancy, los_e, los_f;         // person-seconds [hour * E + edge]
  std::vector<double> queue, queue_max;                // walker-seconds queued, and the longest queue [hour * S + state]
  std::vector<double> counts;                          // passes [(counter * hours + hour) * 2 + way]
  std::vector<double> trips, travel, delay;            // completed trips by departure hour [hour * streams + stream]
  std::vector<float> minutes, minutes_e;               // person-seconds [bin * E + edge], all and at Fruin's E, when kept
  std::vector<double> starts, ends;                    // walks begun and ended at each junction
  double bin = 0;
};

class Simulation {
 public:
  // The network, the targets, and for each route class its routes' values (index 0 is everyone's and the fall-back);
  // the streams; the counters; the options.
  Simulation(std::shared_ptr<const Targets> targets, std::vector<std::shared_ptr<const Values>> routes, std::vector<Stream> streams, std::vector<Counter> counters,
             Options options)
      : targets_(std::move(targets)), routes_(std::move(routes)), streams_(std::move(streams)), counters_(std::move(counters)), o_(options) {
    if (routes_.empty()) throw std::invalid_argument("urban::movement::Simulation: routes for class 0 at least");
    W_ = routes_[0]->routes().shared_walkway();
    const std::size_t S = W_->states(), E = W_->edges();
    base_ = routes_[0];
    alias_.reserve(streams_.size());
    for (const auto& st : streams_) {
      if (st.origin.size() != st.target.size() || st.origin.size() != st.weight.size()) throw std::invalid_argument("urban::movement::Simulation: a stream's pairs and weights differ in number (" + st.name + ")");
      if (st.dwell > 0 && st.back.size() != st.origin.size()) throw std::invalid_argument("urban::movement::Simulation: a stream that dwells needs a target back for each pair (" + st.name + ")");
      if (static_cast<std::size_t>(st.route_class) >= routes_.size()) throw std::invalid_argument("urban::movement::Simulation: a stream of a class with no routes (" + st.name + ")");
      double total = 0;
      for (const double w : st.weight) total += w;
      alias_.push_back(total > 0 ? calx::random::AliasTable(st.weight) : calx::random::AliasTable());
    }
    scale_.assign(streams_.size(), 1.0);
    credit_.assign(S, 0.0);
    for (auto& c : counters_) {
      if (c.edge >= E) throw std::invalid_argument("urban::movement::Simulation: a counter on no segment (" + c.name + ")");
      // A line at a segment's very end would miss the walkers starting from it: kept half a metre inside (or at the middle).
      const double L = W_->segments[c.edge].length, margin = std::min(0.5, L / 2);
      c.along = std::clamp(c.along, margin, L - margin);
    }
    counters_on_.assign(E, {});
    for (std::size_t k = 0; k < counters_.size(); k++) counters_on_[counters_[k].edge].push_back(k);
    T_.hours = o_.hours, T_.states = S, T_.edges = E;
    const std::size_t H = o_.hours, C = counters_.size(), N = streams_.size();
    T_.counts.assign(C * H * 2, 0.0);
    if (!o_.light) {
      T_.entries.assign(H * S, 0.0), T_.entries_by_departure.assign(H * S, 0.0);
      T_.occupancy.assign(H * E, 0.0), T_.los_e.assign(H * E, 0.0), T_.los_f.assign(H * E, 0.0);
      T_.queue.assign(H * S, 0.0), T_.queue_max.assign(H * S, 0.0);
      T_.trips.assign(H * N, 0.0), T_.travel.assign(H * N, 0.0), T_.delay.assign(H * N, 0.0);
      T_.starts.assign(W_->junctions(), 0.0), T_.ends.assign(W_->junctions(), 0.0);
      if (o_.minute_bins > 0) {
        T_.bin = o_.minute_bins;
        T_.minutes.assign(static_cast<std::size_t>(std::ceil(static_cast<double>(H) * 3600.0 / o_.minute_bins)) * E, 0.0f);
        T_.minutes_e = T_.minutes;
      }
    }
    density_.assign(E, 0.0);
    for (int d = 0; d < 2; d++) line_[static_cast<std::size_t>(d)].assign(E, {});
    waiting_.assign(S, 0u);
    aside_.assign(E, {});
    passes_.assign(2 * counters_.size(), 0.0);
  }

  // ---- what may change while it runs ------------------------------------------------------------------------------

  // A class's routes (a regime coming into force): walkers take them from their next junction on.
  void set_routes(std::size_t cls, std::shared_ptr<const Values> v) {
    if (cls >= routes_.size()) throw std::invalid_argument("urban::movement::Simulation::set_routes: no such class");
    routes_[cls] = std::move(v);
  }
  // The widths of segments (a street given to walkers): the capacity and the density follow.
  void set_width(std::size_t edge, double width) {
    if (width_.empty()) {
      width_.resize(W_->edges());
      for (std::size_t e = 0; e < W_->edges(); e++) width_[e] = W_->segments[e].width;
    }
    width_.at(edge) = width;
  }
  double width(std::size_t edge) const { return width_.empty() ? W_->segments[edge].width : width_[edge]; }
  // A stream's rates scaled (a crowd larger or smaller than expected).
  void set_scale(std::size_t stream, double factor) { scale_.at(stream) = factor; }
  double scale(std::size_t stream) const { return scale_.at(stream); }
  // New draws from now on (the walkers present keep their speeds): a planner's copies diverge so.
  void reseed(std::uint64_t seed) { o_.seed = seed; }
  // A copy for a planner: running totals and counters only, at its own step (a coarser model of the same walkers).
  Simulation lightened(double step) const {
    Simulation c = *this;
    c.o_.light = true, c.o_.kinematics.step = step, c.o_.minute_bins = 0;
    Tallies t;
    t.hours = T_.hours, t.states = T_.states, t.edges = T_.edges, t.counts = T_.counts;
    c.T_ = std::move(t);
    return c;
  }

  // ---- running ----------------------------------------------------------------------------------------------------------

  double time() const { return time_; }
  std::uint32_t step_index() const { return tick_; }

  void run(double until, calx::scheduler::ThreadPool* pool = nullptr) {
    const double dt = o_.kinematics.step;
    while (time_ + 0.5 * dt <= until) step(pool);
  }

  void step(calx::scheduler::ThreadPool* pool = nullptr) {
    const double dt = o_.kinematics.step, t0 = time_, t1 = time_ + dt;
    release(t0, t1);
    densities();
    // Phase A: every walker moves and asks, reading the configuration at t0 (nomos's synchronous schedule).
    Mover M{this, t0, t1};
    nomos::agents::tick(M, pop_, nomos::agents::Schedule::synchronous, o_.seed, tick_, nullptr, pool);
    // The counting lines passed and the person-seconds, before admission (on the segments walked this step).
    tally_present(t0, dt);
    // Phase B: admission, segment by segment.
    admit(t0, t1);
    // The arrived leave the population.
    finish();
    time_ = t1;
    tick_++;
  }

  // ---- reading ---------------------------------------------------------------------------------------------------------

  struct Position {
    std::uint32_t id = 0, stream = 0;
    int mode = 0;
    std::size_t state = none;
    double along = 0;
    Point at{}, heading{};
    Status status = Status::walking;
  };
  std::vector<Position> positions() const {
    std::vector<Position> out;
    std::fill(waiting_.begin(), waiting_.end(), 0u);
    for (const auto& w : pop_.walkers) {
      if (w.state == static_cast<std::uint32_t>(none) || w.status == Status::done || w.status == Status::pending) continue;
      Position p;
      p.id = w.id, p.stream = w.stream, p.mode = streams_[w.stream].mode, p.state = w.state, p.along = place(w), p.status = w.status;
      p.at = W_->point(w.state, p.along);
      const Point q = W_->point(w.state, std::min(W_->length(w.state), p.along + 1.0)), r = W_->point(w.state, std::max(0.0, p.along - 1.0));
      const double n = std::hypot(q[0] - r[0], q[1] - r[1]);
      p.heading = n > 0 ? Point{(q[0] - r[0]) / n, (q[1] - r[1]) / n} : W_->heading_out(w.state);
      out.push_back(p);
    }
    return out;
  }
  std::size_t active() const { return pop_.walkers.size(); }
  const std::vector<Walker>& walkers() const { return pop_.walkers; }
  const Tallies& tallies() const { return T_; }
  const Totals& totals() const { return totals_; }
  const std::vector<Counter>& counters() const { return counters_; }
  const std::vector<Stream>& streams() const { return streams_; }
  const Walkway& walkway() const { return *W_; }
  const Targets& targets() const { return *targets_; }
  const Values& routes(std::size_t cls) const { return *routes_.at(cls); }
  const Options& options() const { return o_; }
  // The density (m⁻²) on each segment at the last step's start.
  const std::vector<double>& density() const { return density_; }
  // What a counter reads each way in an hour, in expectation over its footpath (share × (1 ± ρ) × passes).
  std::array<double, 2> reading(std::size_t counter, std::size_t hour) const {
    const auto& c = counters_[counter];
    const std::size_t i = (counter * T_.hours + hour) * 2;
    return {c.share * (1 + c.hand[0]) * T_.counts[i], c.share * (1 + c.hand[1]) * T_.counts[i + 1]};
  }
  // The walkers waiting at each segment's mouth now.
  std::vector<std::size_t> queues() const {
    std::vector<std::size_t> q(W_->states(), 0);
    for (const auto& w : pop_.walkers)
      if (w.status == Status::asking && w.next != static_cast<std::uint32_t>(none) && w.ask < time_) q[w.next]++;
    return q;
  }
  std::uint64_t born() const { return born_; }
  // Each counter's passes each way since the run began, [2 counter + way] (kept in light runs too).
  const std::vector<double>& passes() const { return passes_; }

  // ODD+D (Grimm et al. 2020; Müller et al. 2013), kept beside the model.
  static nomos::agents::Description describe() {
    return {"How walkers distribute over a district's footpaths and crossings, and how street regimes (closures, one-way running, events) move them, "
            "calibrated to the counters by hour and direction.",
            "Walkers (an origin junction, a target, a free speed, a route chosen as they go); the directed segments of the walking network with "
            "their lengths, widths, grades and signals; counters on footpaths.",
            "A step of 2 s; the district's network in its local frame (m); a run of hours, the week's 168 or a day's 24.",
            "Departures by stream (Poisson within each hour); walking at the free speed slowed by the density; at a junction, the next segment by "
            "the logit over walks; waiting for the phase at signals and for room at full segments; dwelling at a target; arriving.",
            "Emergence: link loads, queues and densities from individual choices. Sensing: the density of the segment walked, the signal. "
            "Stochasticity: departures, speeds, choices of way.",
            "A discrete choice of the next segment, P ∝ M_ab z(b): the logit over every walk to the target (Fosgerau et al. 2013); no planning.",
            "No walker at the run's start (or a configuration carried from an earlier run).",
            "The streams' rates by hour (fitted to the counters each way by NNLS), the network, the regimes in force.",
            "route.hpp (the walk sums z and π, the laws of the first and next segments); walkers.hpp (Greenshields's speed, Weidmann's speeds "
            "and jam density, capacity, signals).",
            "Synchronous, so that a step runs on a pool and a run does not depend on the threads; admission ranked by instant and identity."};
  }

 private:
  struct Population {
    std::vector<Walker> walkers;
  };
  // The agents' rule for nomos's tick: agent i of the population moves from t0 to t1, reading only the step's start.
  struct Mover {
    using State = Population;
    const Simulation* sim;
    double t0, t1;
    std::size_t agents(const State& x) const { return x.walkers.size(); }
    void update(std::size_t i, const State& read, State& write, calx::random::Stream&, nomos::agents::Log*) const {
      Walker w = read.walkers[i];
      sim->advance(w, t0, t1);
      write.walkers[i] = w;
    }
  };

  static constexpr std::uint32_t walk_stream = 0x77616c6bu;    // "walk"
  static constexpr std::uint32_t init_stream = 0x696e6974u;    // "init"
  static constexpr std::uint32_t arrive_stream = 0x61727276u;  // "arrv"
  static constexpr std::uint32_t stop_stream = 0x73746f70u;    // "stop"
  static constexpr std::uint32_t none32 = static_cast<std::uint32_t>(none);

  const std::vector<double>& z(const Walker& w) const {
    const std::size_t t = w.leg == 0 ? w.target : w.back;
    return (*routes_[w.cls])[t];
  }
  std::size_t leg_target(const Walker& w) const { return w.leg == 0 ? w.target : w.back; }

  // Departures in [t0, t1): each stream's minute is drawn whole when the step reaches it.
  void release(double t0, double t1) {
    while (next_minute_ * 60.0 < t1) {
      const std::size_t minute = next_minute_++;
      const auto hour = static_cast<std::size_t>(minute / 60);
      std::vector<Walker> fresh;
      for (std::size_t k = 0; k < streams_.size(); k++) {
        const auto& st = streams_[k];
        if (hour >= st.rate.size() || !(st.rate[hour] > 0) || alias_[k].size() == 0) continue;
        calx::random::Stream rng(o_.seed, arrive_stream, static_cast<std::uint32_t>(k), static_cast<std::uint32_t>(minute));
        const std::int64_t n = rng.poisson(scale_[k] * st.rate[hour] / 60.0);
        for (std::int64_t q = 0; q < n; q++) {
          Walker w;
          w.stream = static_cast<std::uint32_t>(k);
          w.pair = static_cast<std::uint32_t>(alias_[k](rng));
          w.depart = 60.0 * static_cast<double>(minute) + 60.0 * rng.uniform();
          fresh.push_back(w);
        }
      }
      std::stable_sort(fresh.begin(), fresh.end(), [](const Walker& a, const Walker& b) { return a.depart < b.depart; });
      for (auto& w : fresh) {
        const auto& st = streams_[w.stream];
        w.id = static_cast<std::uint32_t>(born_++);
        w.origin = static_cast<std::uint32_t>(st.origin[w.pair]);
        w.target = static_cast<std::uint32_t>(st.target[w.pair]);
        w.back = st.dwell > 0 ? static_cast<std::uint32_t>(st.back[w.pair]) : 0;
        w.cls = static_cast<std::uint8_t>(st.route_class);
        w.hour = static_cast<std::uint16_t>(std::min<std::size_t>(hour, 65535));
        calx::random::Stream r(o_.seed, init_stream, w.id, 0);
        w.speed = static_cast<float>(std::clamp(r.normal(o_.kinematics.free_speed, o_.kinematics.speed_sd), 0.5, 2.5));
        w.status = Status::pending;
        pending_.push_back(w);
      }
    }
    // Those departing in this step join the population.
    std::size_t k = 0;
    while (k < pending_.size() && pending_[k].depart < t1) {
      pop_.walkers.push_back(pending_[k]);
      totals_.departed += 1;
      k++;
    }
    pending_.erase(pending_.begin(), pending_.begin() + static_cast<std::ptrdiff_t>(k));
    (void)t0;
  }

  // The segments' mean densities, and each segment's walkers by their place along it (from its tail): those walking
  // (or waiting to walk on) each way, and those dwelling, who stand aside.
  void densities() {
    std::fill(density_.begin(), density_.end(), 0.0);
    for (int d = 0; d < 2; d++)
      for (auto& l : line_[d]) l.clear();
    for (auto& l : aside_) l.clear();
    std::fill(waiting_.begin(), waiting_.end(), 0u);
    for (const auto& w : pop_.walkers) {
      if (w.state == none32 || w.status == Status::done || w.status == Status::pending) continue;
      const std::size_t e = w.state / 2, d = w.state % 2;
      density_[e] += 1;
      const double s = place(w), x = d ? W_->segments[e].length - s : s;
      (w.status == Status::dwelling ? aside_[e] : line_[d][e]).push_back(x);
    }
    for (std::size_t e = 0; e < density_.size(); e++) {
      density_[e] /= W_->segments[e].length * width(e);
      std::sort(line_[0][e].begin(), line_[0][e].end());
      std::sort(line_[1][e].begin(), line_[1][e].end());
      std::sort(aside_[e].begin(), aside_[e].end());
    }
  }
  // Where a walker stands for the densities: a walker waiting at a segment's end (at red, or for room) stands behind
  // those who came before it, at a standing queue's density, rather than all at one point.
  double place(const Walker& w) const {
    const bool waits = (w.status == Status::asking || w.status == Status::signal) && w.s >= W_->length(w.state) - 1e-6;
    if (!waits) return w.s;
    const std::uint32_t rank = waiting_[w.state]++;
    return std::max(0.0, W_->length(w.state) - (rank + 0.5) / (o_.kinematics.queue_density * width(w.state / 2)));
  }
  // Walkers on segment e going way d (2: dwelling, aside) within [lo, hi] of its length (from its tail).
  std::size_t count(std::size_t e, std::size_t d, double lo, double hi) const {
    const auto& l = d < 2 ? line_[d][e] : aside_[e];
    return static_cast<std::size_t>(std::upper_bound(l.begin(), l.end(), hi) - std::lower_bound(l.begin(), l.end(), lo));
  }
  // The densities (m⁻²) in [a, b] of state c (its own coordinate, from its tail): of every walker (walking either way,
  // and dwelling), and of those walking its way; `self` walkers are left out, among those walking its way, or among the
  // dwellers when `self_aside` (the walker reading it).
  std::array<double, 2> local(std::size_t c, double a, double b, std::size_t self, bool self_aside = false) const {
    const std::size_t e = c / 2, d = c % 2;
    const double L = W_->segments[e].length;
    const double lo = d ? L - b : a, hi = d ? L - a : b;
    std::size_t own = count(e, d, lo, hi), other = count(e, 1 - d, lo, hi), aside = count(e, 2, lo, hi);
    if (self_aside) aside = aside > self ? aside - self : 0;
    else own = own > self ? own - self : 0;
    const double area = std::max(b - a, 1e-9) * width(e);
    return {static_cast<double>(own + other + aside) / area, static_cast<double>(own) / area};
  }
  // The density ahead of a walker at s on state c: the next `window` metres, or what is left of the segment (its area
  // taken as at least half a metre's), the walker itself left out.
  std::array<double, 2> ahead(std::size_t c, double s, std::size_t self) const {
    const double L = W_->length(c), b = std::min(L, s + o_.kinematics.window);
    auto rho = local(c, s, std::max(b, s + 1e-9), self);
    const double shrink = std::max(b - s, 1e-9) / std::max(b - s, 0.5);
    return {rho[0] * shrink, rho[1] * shrink};
  }

  double speed_on(const Walker& w, std::size_t c, double s, std::size_t self) const {
    const auto& k = o_.kinematics;
    const auto rho = ahead(c, s, self);
    if (rho[1] >= 0.95 * k.jam) return 0;   // a full queue ahead: it stands
    const double v = w.speed * (k.grade ? tobler(W_->grade(c)) : 1.0) * (1 - std::min(1.0, rho[0] / k.jam));
    return std::max(k.min_speed, v);
  }

  // Green for the phase that lets a walker leave along b at instant t? and the start of its next green.
  double next_green(std::size_t b, double t) const {
    const auto& k = o_.kinematics;
    const int p = W_->phase(b);
    const double from = p == 0 ? 0.0 : k.green * k.cycle, to = p == 0 ? k.green * k.cycle : k.cycle;
    const double offset = std::fmod(17.0 * static_cast<double>(W_->tail(b)), k.cycle);  // a fixed offset per junction
    const double phase = std::fmod(t + offset, k.cycle);
    if (phase >= from && phase < to) return t;
    return phase < from ? t + (from - phase) : t + (k.cycle - phase) + from;
  }

  // At the end of segment c at instant t: arrive, dwell, turn for the walk back, or choose the next segment and ask.
  void at_end(Walker& w, std::size_t c, double t, calx::random::Stream& rng) const {
    if (targets_->contains(leg_target(w), c)) {
      if (w.leg == 0 && streams_[w.stream].dwell > 0) {
        // The walk back from the junction it has reached.
        w.leg = 1;
        const std::size_t x = W_->head(c);
        const auto law = routes_[w.cls]->routes().first_law(x, z(w));
        std::size_t b = Routes::draw(rng, law);
        if (b == none) {
          const auto fb = base_->routes().first_law(x, (*base_)[w.back]);
          b = Routes::draw(rng, fb);
        }
        if (b == none) { w.status = Status::done, w.failed = 1, w.end = t; return; }
        w.next = static_cast<std::uint32_t>(b), w.ask = t, w.status = Status::asking;
        return;
      }
      w.status = Status::done, w.end = t;
      return;
    }
    auto law = routes_[w.cls]->routes().next_law(c, z(w));
    if (law.empty()) law = base_->routes().next_law(c, (*base_)[leg_target(w)]);
    const std::size_t b = Routes::draw(rng, law);
    if (b == none) { w.status = Status::done, w.failed = 1, w.end = t; return; }
    w.next = static_cast<std::uint32_t>(b);
    double when = t;
    if (o_.kinematics.signals && W_->signalled[W_->head(c)]) {
      const double cross = W_->crossings(c, b);
      if (cross >= 1 || (cross > 0 && rng.uniform() < cross)) when = next_green(b, t);
    }
    w.ask = when;
    w.status = when > t ? Status::signal : Status::asking;
  }

  // Phase A for one walker.
  void advance(Walker& w, double t0, double t1) const {
    calx::random::Stream rng(o_.seed, walk_stream, w.id, tick_);
    w.s0 = w.s;
    double now = t0;
    for (int guard = 0; guard < 8 && now < t1; guard++) {
      switch (w.status) {
        case Status::pending: {
          now = std::max(now, w.depart);
          const auto law = routes_[w.cls]->routes().first_law(w.origin, z(w));
          const std::size_t b = Routes::draw(rng, law);
          if (b == none) { w.status = Status::done, w.failed = 1, w.end = now; return; }
          w.next = static_cast<std::uint32_t>(b), w.ask = now, w.status = Status::asking;
          return;
        }
        case Status::walking: {
          const std::size_t c = w.state;
          const double L = W_->length(c), v = speed_on(w, c, w.s, 1);
          if (!(v > 0)) return;   // standing in a queue
          const double limit = (w.stop >= 0 && !w.dwelt) ? w.stop : L;
          const double reach = w.s + v * (t1 - now);
          if (reach < limit) { w.walked += reach - w.s, w.s = reach; return; }
          now += (limit - w.s) / v;
          w.walked += limit - w.s, w.s = limit;
          if (limit < L) {
            // Its dwell begins.
            calx::random::Stream d(o_.seed, stop_stream, w.id, 1);
            w.rest = streams_[w.stream].dwell * d.exponential();
            w.until = now + w.rest;
            w.dwelt = 1, w.status = Status::dwelling;
            continue;
          }
          at_end(w, c, now, rng);
          if (w.status == Status::asking || w.status == Status::signal || w.status == Status::done) return;
          continue;
        }
        case Status::dwelling: {
          if (w.until >= t1) return;
          now = std::max(now, w.until);
          w.status = Status::walking;
          continue;
        }
        case Status::signal: {
          if (w.ask >= t1) return;
          w.status = Status::asking;
          return;
        }
        case Status::asking: {
          w.ask = std::max(w.ask, t0);
          return;
        }
        case Status::done:
          return;
      }
    }
  }

  void count_passes(std::size_t c, double from, double to, double t) {
    for (const std::size_t k : counters_on_[c / 2]) {
      const auto& K = counters_[k];
      const double at = c % 2 ? W_->length(c) - K.along : K.along;
      if (at > from && at <= to) {
        const int way = c == K.way[0] ? 0 : (c == K.way[1] ? 1 : -1);
        const auto hour = static_cast<std::size_t>(std::floor(t / 3600));
        if (way >= 0) passes_[2 * k + static_cast<std::size_t>(way)] += 1;
        if (way >= 0 && hour < T_.hours) T_.counts[(k * T_.hours + hour) * 2 + static_cast<std::size_t>(way)] += 1;
      }
    }
  }

  void tally_present(double t0, double dt) {
    const auto hour = static_cast<std::size_t>(std::floor(t0 / 3600));
    const bool keep = !o_.light && hour < T_.hours;
    const std::size_t E = W_->edges(), S = W_->states();
    const std::size_t bin = (keep && T_.bin > 0) ? static_cast<std::size_t>(std::floor(t0 / T_.bin)) : 0;
    const bool minutes = keep && T_.bin > 0 && (bin + 1) * E <= T_.minutes.size();
    std::vector<double> q(keep ? S : 0, 0.0);
    for (const auto& w : pop_.walkers) {
      if (w.status == Status::pending) continue;
      const bool waiting = w.status == Status::asking || w.status == Status::signal;
      if (w.status == Status::asking && w.next != none32 && w.ask < t0 + dt) {
        if (keep) T_.queue[hour * S + w.next] += dt, q[w.next] += 1;
      }
      if (w.state == none32) {
        // At its origin, asking for its first segment.
        if (w.status != Status::done) totals_.present += dt, totals_.waiting += dt;
        continue;
      }
      const std::size_t c = w.state, e = c / 2;
      if (w.s > w.s0) count_passes(c, w.s0, w.s, t0);
      if (w.status == Status::done) continue;
      totals_.present += dt;
      // The density about it (5 m either way), from the step's start: for walkers walking or dwelling. A walker waiting
      // at red or for room stands in a queue whose density is assumed (queue_density), not read; Fruin's levels for
      // walkways do not apply to it.
      const double L = W_->length(c);
      const double a = std::max(0.0, w.s0 - 0.5 * o_.kinematics.window), b = std::min(L, std::max(a + 2.0, w.s0 + 0.5 * o_.kinematics.window));
      const double rho = waiting ? 0.0 : local(c, a, b, 1, w.status == Status::dwelling)[0];
      if (rho >= o_.kinematics.los_e) totals_.los_e += dt;
      if (rho >= o_.kinematics.los_f) totals_.los_f += dt;
      if (waiting) totals_.waiting += dt;
      if (keep) {
        T_.occupancy[hour * E + e] += dt;
        if (rho >= o_.kinematics.los_e) T_.los_e[hour * E + e] += dt;
        if (rho >= o_.kinematics.los_f) T_.los_f[hour * E + e] += dt;
        if (minutes) {
          T_.minutes[bin * E + e] += static_cast<float>(dt);
          if (rho >= o_.kinematics.los_e) T_.minutes_e[bin * E + e] += static_cast<float>(dt);
        }
      }
    }
    // The longest queue at each mouth in the hour.
    if (keep)
      for (std::size_t s = 0; s < S; s++) T_.queue_max[hour * S + s] = std::max(T_.queue_max[hour * S + s], q[s]);
  }

  void admit(double t0, double t1) {
    const auto& k = o_.kinematics;
    struct Ask {
      std::size_t state;
      double at;
      std::uint32_t id;
      std::size_t index;
    };
    std::vector<Ask> asks;
    for (std::size_t i = 0; i < pop_.walkers.size(); i++) {
      const auto& w = pop_.walkers[i];
      if (w.status == Status::asking && w.ask < t1) asks.push_back({w.next, w.ask, w.id, i});
    }
    std::sort(asks.begin(), asks.end(), [](const Ask& a, const Ask& b) {
      if (a.state != b.state) return a.state < b.state;
      if (a.at != b.at) return a.at < b.at;
      return a.id < b.id;
    });
    const double dt = t1 - t0;
    std::size_t q = 0;
    while (q < asks.size()) {
      const std::size_t b = asks[q].state;
      std::size_t r = q;
      while (r < asks.size() && asks[r].state == b) r++;
      const std::size_t e = b / 2;
      // Godunov's supply (the receiving flow of the cell transmission model, Daganzo 1994, Transp. Res. B 28, 269–287;
      // Lebacque 1996): the capacity while the segment is at or below the critical density ρ_j / 2, else the flow its
      // density carries, so that a queue discharges at capacity and the segment does not jam behind it.
      // The density of the walkers going b's way in the first metres past its mouth (those coming the other way, among
      // them a queue at red on the corner, slow a walker but do not take its room).
      const double rho = local(b, 0, std::min(W_->length(b), k.window), 0)[1], cap = k.free_speed * k.jam * width(e) / 4;   // walkers a second
      const double supply = rho <= k.jam / 2 ? cap : k.free_speed * rho * (1 - rho / k.jam) * width(e);
      double& credit = credit_[b];
      credit = std::min(credit + supply * dt, supply * dt + 1);
      // Full: those walking its way at 0.95 of the jam, or everyone there, dwellers too, at the jam itself.
      const bool full = k.capacity && (rho >= 0.95 * k.jam || local(b, 0, std::min(W_->length(b), k.window), 0)[0] >= k.jam);
      for (std::size_t j = q; j < r; j++) {
        Walker& w = pop_.walkers[asks[j].index];
        if (k.capacity && (full || credit < 1)) continue;
        if (k.capacity) credit -= 1;
        enter(w, b, std::max(w.ask, t0), t1);
      }
      q = r;
    }
  }

  void enter(Walker& w, std::size_t b, double at, double t1) {
    const auto hour = static_cast<std::size_t>(std::floor(at / 3600));
    if (!o_.light && hour < T_.hours) T_.entries[hour * W_->states() + b] += 1;
    if (!o_.light && w.hour < T_.hours) T_.entries_by_departure[static_cast<std::size_t>(w.hour) * W_->states() + b] += 1;
    if (!o_.light && w.state == none32) T_.starts[W_->tail(b)] += 1;
    w.state = static_cast<std::uint32_t>(b), w.next = none32;
    w.stop = -1;
    const double L = W_->length(b);
    w.alone += L / w.speed;
    if (w.leg == 0 && streams_[w.stream].dwell > 0 && !w.dwelt && targets_->contains(w.target, b)) {
      calx::random::Stream d(o_.seed, stop_stream, w.id, 0);
      w.stop = L * d.uniform();
    }
    const double v = std::max(o_.kinematics.min_speed, speed_on(w, b, 0, 0));
    const double limit = w.stop >= 0 ? w.stop : L * (1 - 1e-9);
    w.s = std::min(limit, v * (t1 - at));
    w.walked += w.s;
    w.s0 = 0;
    count_passes(b, 0, w.s, at);
    w.status = Status::walking;
    if (w.stop >= 0 && w.s >= w.stop) {
      calx::random::Stream d(o_.seed, stop_stream, w.id, 1);
      w.rest = streams_[w.stream].dwell * d.exponential();
      w.until = at + w.s / v + w.rest;  // it reached its place within the step
      w.dwelt = 1, w.status = Status::dwelling;
    }
  }

  void finish() {
    std::size_t keep = 0;
    for (std::size_t i = 0; i < pop_.walkers.size(); i++) {
      const Walker& w = pop_.walkers[i];
      if (w.status != Status::done) {
        if (keep != i) pop_.walkers[keep] = w;
        keep++;
        continue;
      }
      if (!o_.light && w.state != none32) T_.ends[W_->head(w.state)] += 1;
      if (w.failed) { totals_.failed += 1; continue; }
      totals_.arrived += 1;
      const double travel = w.end - w.depart - w.rest, delay = travel - w.alone;
      totals_.travel += travel;
      totals_.delay += delay;
      if (!o_.light && w.hour < T_.hours) {
        const std::size_t i2 = static_cast<std::size_t>(w.hour) * streams_.size() + w.stream;
        T_.trips[i2] += 1, T_.travel[i2] += travel, T_.delay[i2] += delay;
      }
    }
    pop_.walkers.resize(keep);
  }

  std::shared_ptr<const Targets> targets_;
  std::vector<std::shared_ptr<const Values>> routes_;
  std::shared_ptr<const Values> base_;
  std::shared_ptr<const Walkway> W_;
  std::vector<Stream> streams_;
  std::vector<calx::random::AliasTable> alias_;
  std::vector<double> scale_;
  std::vector<Counter> counters_;
  std::vector<std::vector<std::size_t>> counters_on_;
  Options o_;
  Population pop_;
  std::vector<Walker> pending_;
  std::vector<double> density_, credit_, width_, passes_;
  std::array<std::vector<std::vector<double>>, 2> line_;   // [way][edge]: places along the edge, sorted
  mutable std::vector<std::uint32_t> waiting_;              // per state, walkers placed in its end's queue this step
  std::vector<std::vector<double>> aside_;                  // [edge]: dwellers' places along it, sorted
  Tallies T_;
  Totals totals_;
  double time_ = 0;
  std::uint32_t tick_ = 0;
  std::size_t next_minute_ = 0;
  std::uint64_t born_ = 0;
};

}  // namespace urban::movement
