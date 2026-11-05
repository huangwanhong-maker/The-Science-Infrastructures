// The regime of a street as a decision under partial observation (theory/movement.md §6; ADR 0071, the design §9).
//
// A street (a run of segments) is held in one of three regimes, decided every epoch:
//   open      everyone may walk it, both ways;
//   one-way   it runs one way, from its entry end (a phased street, as crowds are run at events): both classes;
//   closed    closed to the walkers passing through (class 0), open to the street's own walkers (class 1: those whose
//             trip starts or ends on it, a crowd it holds among them; the caller puts their streams in class 1).
// A regime is the routes of each class with the street's states opened or closed (route.hpp): walkers take them from
// their next junction on.
//
// The decision is a POMDP: its hidden state is the configuration of the walkers (the crowd on the links no sensor
// sees, and how large the event's crowd is), its observations the counts at the sensors in the epoch, binned, its
// actions the regimes, its reward −(person-hours walked and waited) − λ (person-hours on segments above Fruin's level
// F, 2 m⁻²): walkers' time and their safety, λ an assumption. It is solved online by tyche's POMCP (Silver and Veness
// 2010) on the agent-based model itself as the generative model: a state is a Simulation, a step runs it for an
// epoch under the action's regime, its draws reseeded from the planner's stream.
//
// The exactly solvable case (`Corridor`): two ways from one place to another, one with a sensor, one without; groups
// arrive at a rate that is normal or a surge, hidden and Markov; each chooses the unsensed way by a logit whose share
// the regime sets; groups on it leave each epoch with a chance; it holds a few groups at most. Its law is a finite
// POMDP exactly (Poisson thinning: the groups on the two ways are independent Poisson numbers), which tyche solves with
// bounds, and its agents (each group an agent with its own draws, its choice nomos's logit) are a generative model of
// the same law, on which POMCP is checked against the exact value.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <nomos/agents/choice.hpp>
#include <tyche/decision/online.hpp>
#include <tyche/decision/pomdp.hpp>

#include <urban/movement/route.hpp>
#include <urban/movement/walkers.hpp>

namespace urban::movement {

enum class Regime : std::uint8_t { open = 0, one_way = 1, closed = 2 };
inline const char* regime_name(Regime r) { return r == Regime::open ? "open" : r == Regime::one_way ? "one-way" : "closed"; }

struct Street {
  std::string name;
  std::vector<std::size_t> edges;   // its segments
  std::size_t entry = 0;            // the junction a one-way street is entered from
};

// Each of the street's states, oriented away from its entry (walking the street's segments outward from it): +1 for a
// state running away from the entry, −1 for its reverse.
inline std::vector<int> orientation(const Walkway& W, const Street& S) {
  std::vector<int> o(W.states(), 0);
  std::vector<char> seen(W.junctions(), 0), used(W.edges(), 0), on(W.edges(), 0);
  for (const std::size_t e : S.edges) on.at(e) = 1;
  std::vector<std::size_t> frontier{S.entry};
  seen.at(S.entry) = 1;
  while (!frontier.empty()) {
    const std::size_t x = frontier.back();
    frontier.pop_back();
    for (const std::size_t s : W.leaving[x]) {
      const std::size_t e = s / 2;
      if (!on[e] || used[e]) continue;
      used[e] = 1, o[s] = 1, o[Walkway::reverse(s)] = -1;
      if (!seen[W.head(s)]) seen[W.head(s)] = 1, frontier.push_back(W.head(s));
    }
  }
  for (const std::size_t e : S.edges)
    if (!used[e]) throw std::invalid_argument("urban::movement::orientation: the street " + S.name + " is not connected to its entry");
  return o;
}

// The states open to a class under a street's regime (on top of `base`, what is open otherwise).
inline std::vector<char> open_states(const Walkway& W, const Street& S, Regime r, int cls, std::vector<char> base = {}) {
  if (base.empty()) base.assign(W.states(), 1);
  if (r == Regime::open) return base;
  if (r == Regime::one_way) {
    const auto o = orientation(W, S);
    for (std::size_t s = 0; s < W.states(); s++)
      if (o[s] < 0) base[s] = 0;
    return base;
  }
  if (cls == 0)
    for (const std::size_t e : S.edges) base[2 * e] = base[2 * e + 1] = 0;
  return base;
}

// The routes' values of both classes under each regime of a street: [regime][class].
struct RegimeRoutes {
  Street street;
  std::array<std::array<std::shared_ptr<const Values>, 2>, 3> values;
};

inline RegimeRoutes regime_routes(std::shared_ptr<const Walkway> W, const RouteLaw& law, const Street& S, const Targets& T, calx::scheduler::ThreadPool* pool = nullptr) {
  RegimeRoutes R;
  R.street = S;
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 2; c++) {
      auto routes = std::make_shared<const Routes>(W, law, open_states(*W, S, static_cast<Regime>(r), c));
      R.values[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] = std::make_shared<const Values>(routes, T, pool);
    }
  return R;
}

inline void apply(Simulation& sim, const RegimeRoutes& R, Regime r) {
  for (std::size_t c = 0; c < 2; c++) sim.set_routes(c, R.values[static_cast<std::size_t>(r)][c]);
}

struct RegimeReward {
  double time = 1;     // per person-hour present (walking, waiting, dwelling)
  double crush = 10;   // per person-hour on segments above Fruin's level F (an assumption)
};

// The street-regime POMDP on the walkers (a generative model for tyche::decision::Pomcp).
class RegimeModel {
 public:
  struct State {
    Simulation sim;
    std::size_t epoch = 0;
  };

  struct Setup {
    std::shared_ptr<const RegimeRoutes> routes;
    double epoch = 1800;                       // s
    std::size_t epochs = 1;                    // the decisions left from the start
    std::vector<std::size_t> sensed;           // counters observed
    std::vector<double> bins;                  // upper edges of the bins of a sensed counter's reading in an epoch
    std::vector<std::size_t> scaled;           // streams whose scale is hidden
    std::vector<double> scales, prior;         // the scales and their prior law
    RegimeReward reward;
    double discount = 0.95;
  };

  RegimeModel(Simulation start, Setup setup) : start_(std::move(start)), S_(std::move(setup)) {
    if (S_.scales.size() != S_.prior.size() || S_.scales.empty()) throw std::invalid_argument("urban::movement::RegimeModel: a prior over the hidden scales");
  }

  std::size_t actions() const { return 3; }
  double discount() const { return S_.discount; }
  const Setup& setup() const { return S_; }

  State initial(calx::random::Stream& rng) const {
    State x{start_, 0};
    const double scale = S_.scales[calx::random::categorical(rng, S_.prior)];
    for (const std::size_t k : S_.scaled) x.sim.set_scale(k, scale);
    x.sim.reseed(rng.next64());
    return x;
  }

  // The key of the binned readings of the sensed counters in an epoch.
  std::uint64_t observe(const std::vector<double>& before, const Simulation& sim) const {
    std::uint64_t key = 0;
    const std::uint64_t B = S_.bins.size() + 1;
    for (const std::size_t k : S_.sensed) {
      const auto& c = sim.counters()[k];
      const double read = c.share * ((1 + c.hand[0]) * (sim.passes()[2 * k] - before[2 * k]) + (1 + c.hand[1]) * (sim.passes()[2 * k + 1] - before[2 * k + 1]));
      std::uint64_t b = 0;
      while (b < S_.bins.size() && read > S_.bins[b]) b++;
      key = key * B + b;
    }
    return key;
  }

  double reward(const Totals& before, const Totals& after) const {
    return -S_.reward.time * (after.present - before.present) / 3600 - S_.reward.crush * (after.los_f - before.los_f) / 3600;
  }

  // One epoch under action a; its observation and reward.
  tyche::decision::Step<State> step(const State& x, std::size_t a, calx::random::Stream& rng) const {
    tyche::decision::Step<State> out{x, 0, 0};
    if (x.epoch >= S_.epochs) return out;   // past the last decision: nothing more is counted
    Simulation& sim = out.next.sim;
    sim.reseed(rng.next64());
    apply(sim, *S_.routes, static_cast<Regime>(a));
    const Totals before = sim.totals();
    const std::vector<double> passes = sim.passes();
    sim.run(sim.time() + S_.epoch);
    out.next.epoch++;
    out.observation = observe(passes, sim);
    out.reward = reward(before, sim.totals());
    return out;
  }

  // Runs the true world for an epoch under a, returning its observation and reward (the world's own draws).
  std::pair<std::uint64_t, double> act(Simulation& world, std::size_t a) const {
    apply(world, *S_.routes, static_cast<Regime>(a));
    const Totals before = world.totals();
    const std::vector<double> passes = world.passes();
    world.run(world.time() + S_.epoch);
    return {observe(passes, world), reward(before, world.totals())};
  }

 private:
  Simulation start_;
  Setup S_;
};

// ---- the exactly solvable case -----------------------------------------------------------------------------------------

struct Corridor {
  std::array<double, 2> demand{2.0, 6.0};            // groups arriving an epoch: normal, surge (λ_k)
  double persist = 0.8;                             // P(the demand level stays)
  std::size_t capacity = 4;                         // groups the unsensed way holds
  double leave = 0.5;                               // q: a group's chance to leave it in an epoch
  std::array<double, 3> share{0.5, 0.25, 0.0};      // p_a: the unsensed way's share under open, one-way, closed
  std::array<double, 3> detour{0.0, 0.15, 0.4};     // cost per arriving group of the regime's detours
  double crush = 0.5;                               // cost per group² on the unsensed way
  std::size_t observations = 8;                     // the sensed way's count, capped
  double discount = 0.9;

  std::size_t states() const { return 2 * (capacity + 1); }
  std::size_t state(std::size_t level, std::size_t groups) const { return level * (capacity + 1) + groups; }
  double reward(std::size_t level, std::size_t groups, std::size_t a) const {
    return -crush * static_cast<double>(groups * groups) - detour[a] * demand[level];
  }
};

namespace detail {
inline double poisson_pmf(std::size_t k, double mean) {
  if (mean <= 0) return k == 0 ? 1.0 : 0.0;
  return std::exp(static_cast<double>(k) * std::log(mean) - mean - std::lgamma(static_cast<double>(k) + 1));
}
inline double binomial_pmf(std::size_t k, std::size_t n, double p) {
  if (k > n) return 0;
  if (p <= 0) return k == 0 ? 1.0 : 0.0;
  if (p >= 1) return k == n ? 1.0 : 0.0;
  const double nn = static_cast<double>(n), kk = static_cast<double>(k);
  return std::exp(std::lgamma(nn + 1) - std::lgamma(kk + 1) - std::lgamma(nn - kk + 1) + kk * std::log(p) + (nn - kk) * std::log(1 - p));
}
}  // namespace detail

// The corridor's law as a finite POMDP, exactly: (k, n) → k' by the persistence, n' = min(N, Bin(n, 1 − q) + Poisson(λ_k' p_a)),
// o = min(O − 1, Poisson(λ_k' (1 − p_a))).
inline tyche::decision::Pomdp corridor_pomdp(const Corridor& C) {
  const std::size_t N = C.capacity, S = C.states(), O = C.observations;
  tyche::decision::Pomdp P(S, 3, O, C.discount);
  for (std::size_t a = 0; a < 3; a++)
    for (std::size_t k = 0; k < 2; k++)
      for (std::size_t n = 0; n <= N; n++) {
        const std::size_t s = C.state(k, n);
        P.reward(s, a) = C.reward(k, n, a);
        for (std::size_t k2 = 0; k2 < 2; k2++) {
          const double pk = k2 == k ? C.persist : 1 - C.persist;
          const double mu = C.demand[k2] * C.share[a];
          for (std::size_t j = 0; j <= n; j++) {
            const double pj = detail::binomial_pmf(j, n, 1 - C.leave);
            double tail = 1;
            for (std::size_t u = 0; j + u < N; u++) {
              const double pu = detail::poisson_pmf(u, mu);
              P.transition(a, s, C.state(k2, j + u)) += pk * pj * pu;
              tail -= pu;
            }
            P.transition(a, s, C.state(k2, N)) += pk * pj * std::max(0.0, tail);
          }
        }
      }
  for (std::size_t a = 0; a < 3; a++)
    for (std::size_t k2 = 0; k2 < 2; k2++) {
      const double mu = C.demand[k2] * (1 - C.share[a]);
      double tail = 1;
      std::vector<double> z(O, 0.0);
      for (std::size_t o = 0; o + 1 < O; o++) tail -= (z[o] = detail::poisson_pmf(o, mu));
      z[O - 1] = std::max(0.0, tail);
      for (std::size_t n = 0; n <= N; n++)
        for (std::size_t o = 0; o < O; o++) P.observe(a, C.state(k2, n), o) = z[o];
    }
  P.initial.assign(S, 0.0);
  P.initial[C.state(0, 0)] = 1;
  return P;
}

// The corridor's agents: each group on the unsensed way an agent of its own identity, leaving by its own draw; each
// arriving group choosing its way by nomos's logit. A generative model of the same law as corridor_pomdp.
struct CorridorAgents {
  struct State {
    std::size_t level = 0;
    std::vector<std::uint32_t> groups;   // the identities of the groups on the unsensed way
    std::uint32_t next = 0;              // the next identity
  };
  Corridor C;

  std::size_t actions() const { return 3; }
  double discount() const { return C.discount; }
  State initial(calx::random::Stream&) const { return {}; }
  std::size_t index(const State& x) const { return C.state(x.level, x.groups.size()); }

  tyche::decision::Step<State> step(const State& x, std::size_t a, calx::random::Stream& rng) const {
    tyche::decision::Step<State> out;
    out.reward = C.reward(x.level, x.groups.size(), a);
    State& y = out.next;
    y.next = x.next;
    y.level = rng.bernoulli(C.persist) ? x.level : 1 - x.level;
    // Each group on the way leaves or stays by its own draw.
    const std::uint64_t key = rng.next64();
    for (const std::uint32_t g : x.groups) {
      calx::random::Stream own(key, 0x67727570u, g, 0);  // "grup"
      if (!own.bernoulli(C.leave)) y.groups.push_back(g);
    }
    // The arrivals choose their way.
    const std::int64_t arrivals = rng.poisson(C.demand[y.level]);
    std::size_t sensed = 0;
    const double p = C.share[a];
    for (std::int64_t q = 0; q < arrivals; q++) {
      const std::uint32_t g = y.next++;
      calx::random::Stream own(key, 0x63686f6fu, g, 1);  // "choo"
      bool unsensed = false;
      if (p >= 1) unsensed = true;
      else if (p > 0) {
        const double u[2] = {std::log(p), std::log(1 - p)};
        unsensed = nomos::agents::choose(own, std::span<const double>(u, 2), 1.0) == 0;
      }
      if (!unsensed) sensed++;
      else if (y.groups.size() < C.capacity) y.groups.push_back(g);   // a full way turns the group away
    }
    out.observation = std::min<std::size_t>(sensed, C.observations - 1);
    return out;
  }
};

}  // namespace urban::movement
