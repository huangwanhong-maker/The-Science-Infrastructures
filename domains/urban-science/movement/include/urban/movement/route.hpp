// The walkers' choice of way: a logit over every walk of the network (theory/movement.md §3).
//
// The law is the district-economy flow study's (research/experiments/2026-11-05-a-district-economy-as-a-flow, flow.hpp,
// README §2.2), so that the mean of the walkers is that flow: a logit over walks (Bell 1995, "Alternatives to Dial's
// logit assignment algorithm", Transp. Res. B 29, 287–295; Akamatsu 1996, "Cyclic flows, Markov process and stochastic
// traffic assignment", Transp. Res. B 30, 369–386) whose states are the directed segments, so that a walker turns back
// only at a dead end (the link-based form of Fosgerau, Frejinger and Karlström 2013, "A link based network route choice
// model with unrestricted choice set", Transp. Res. B 56, 70–80). Entering state b from state a has the weight
//
//   M_ab = ω_b τ_ab,   ω_b = exp(−θ ℓ_b) (w_b / w_ref)^(κ ℓ_b / ℓ_parent) exp(−θ γ ℓ_b (1/φ(g_b) − 1)),
//   τ_ab = exp(−β_x c_ab − β_s c_ab [x signalled] (1 − g)² C / 2),
//
// b leaving the junction x that a enters, b not a reversed unless x leaves no other way open. Here θ is the route's
// sensitivity to length, w_b the walkable width, κ the weight of width (1: between two parallel segments walkers divide
// in proportion to width, as a fluid's flux divides by conductance), φ Tobler's factor of speed on the grade g_b (γ = 0:
// grade ignored), c_ab the carriageways crossed (network.hpp), β_x their cost, β_s the cost of a second of the expected
// wait at a signal of cycle C and green share g, (1 − g)² C / 2 for a walker arriving at random (the red's share times
// its mean remainder). With β_x = β_s = γ = 0 and κ = 1 the law is the flow study's.
//
// A trip from junction o ends on first entering a target A, a set of states: those entering a junction d, or the
// states of a set of segments (a street a crowd is bound for). With z_A(a) the weight of every walk from state a to A,
//
//   z_A(a) = 1 (a ∈ A),   z_A(a) = Σ_b M_ab z_A(b) (a ∉ A),
//
// a walker leaves o along b with probability ω_b z_A(b) / W_oA (W_oA = Σ_{b out of o} ω_b z_A(b)) and, from a ∉ A, takes
// b with probability M_ab z_A(b) / z_A(a): the logit over walks as a Markov chain (Doob's h-transform of the walk by
// z_A). Its mean visits to a state c are those of the flow study,
//
//   f_c(o → A) = π_oA(c) z_A(c) / W_oA,   π_oA(c) = ω_c [c out of o] + Σ_{a ∉ A} π_oA(a) M_ac,
//
// π the weight of the walks from o that reach c before A. The sums converge when the spectral radius of M is below one
// (checked when the routes are made: the walk sums diverge for θ too small). z and π are solved by Gauss–Seidel sweeps,
// each state to a declared relative tolerance, which needs no inverse of the transitions and serves networks of any size.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <nomos/agents/choice.hpp>

#include <urban/movement/network.hpp>

namespace urban::movement {

struct RouteLaw {
  double theta = 0.02;            // 1/m (the flow study's calibration over a grid, its README §3)
  double reference_width = 12;    // w_ref, m
  double width_exponent = 1;      // κ
  double crossing = 0;            // β_x, per carriageway crossed
  double signal = 0;              // β_s, per second of the expected wait at a signal
  double grade = 0;               // γ, the weight of Tobler's time penalty
  double cycle = 120, green = 0.5;   // the signals' cycle (s) and the green share of a phase
};

// The destinations of trips: each a sorted set of states.
class Targets {
 public:
  std::size_t add(std::vector<std::size_t> states) {
    std::sort(states.begin(), states.end());
    states.erase(std::unique(states.begin(), states.end()), states.end());
    if (states.empty()) throw std::invalid_argument("urban::movement::Targets: a target of no state");
    const auto it = index_.find(states);
    if (it != index_.end()) return it->second;
    sets_.push_back(states);
    index_.emplace(std::move(states), sets_.size() - 1);
    return sets_.size() - 1;
  }
  // Arriving at junction d.
  std::size_t junction(const Walkway& W, std::size_t d) {
    if (const auto it = junction_.find(d); it != junction_.end()) return it->second;
    const std::size_t id = add(W.entering.at(d));
    junction_[d] = id;
    return id;
  }
  // Entering any of a set of segments, either way.
  std::size_t segments(const std::vector<std::size_t>& edges) {
    std::vector<std::size_t> s;
    for (const std::size_t e : edges) s.push_back(2 * e), s.push_back(2 * e + 1);
    return add(std::move(s));
  }
  std::size_t size() const { return sets_.size(); }
  const std::vector<std::size_t>& states(std::size_t id) const { return sets_.at(id); }
  bool contains(std::size_t id, std::size_t s) const { return std::binary_search(sets_[id].begin(), sets_[id].end(), s); }
  // The junction a target arrives at, or none.
  std::size_t junction_of(std::size_t id) const {
    for (const auto& [d, t] : junction_)
      if (t == id) return d;
    return none;
  }

 private:
  std::vector<std::vector<std::size_t>> sets_;
  std::map<std::vector<std::size_t>, std::size_t> index_;
  std::map<std::size_t, std::size_t> junction_;
};

class Routes {
 public:
  struct Arc {
    std::size_t to = 0;
    double w = 0;
  };

  // The routes of a law on a network, with the states open to entry (all, by default). A closed state can be left
  // (a walker on it walks on) but not entered.
  Routes(std::shared_ptr<const Walkway> network, RouteLaw law, std::vector<char> open = {}) : W_(std::move(network)), L_(law), open_(std::move(open)) {
    const Walkway& W = *W_;
    const std::size_t S = W.states();
    if (open_.empty()) open_.assign(S, 1);
    if (open_.size() != S) throw std::invalid_argument("urban::movement::Routes: one entry of `open` per state");
    omega_.assign(S, 0.0);
    for (std::size_t s = 0; s < S; s++) {
      const auto& e = W.segments[s / 2];
      double logw = -L_.theta * e.length + L_.width_exponent * (e.length / e.parent) * std::log(e.width / L_.reference_width);
      if (L_.grade != 0) logw -= L_.theta * L_.grade * e.length * (1 / tobler(W.grade(s)) - 1);
      omega_[s] = open_[s] ? std::exp(logw) : 0.0;
    }
    const double wait = (1 - L_.green) * (1 - L_.green) * L_.cycle / 2;
    out_.assign(S, {}), in_.assign(S, {});
    for (std::size_t a = 0; a < S; a++) {
      const std::size_t x = W.head(a);
      std::size_t others = 0;
      for (const std::size_t b : W.leaving[x])
        if (open_[b] && b != Walkway::reverse(a)) others++;
      for (const std::size_t b : W.leaving[x]) {
        if (!open_[b]) continue;
        if (b == Walkway::reverse(a) && others > 0) continue;
        const double c = W.crossings(a, b);
        const double t = std::exp(-L_.crossing * c - (W.signalled[x] ? L_.signal * c * wait : 0.0));
        out_[a].push_back({b, omega_[b] * t});
        in_[b].push_back({a, omega_[b] * t});
      }
    }
  }

  const Walkway& walkway() const { return *W_; }
  std::shared_ptr<const Walkway> shared_walkway() const { return W_; }
  const RouteLaw& law() const { return L_; }
  bool open(std::size_t s) const { return open_[s] != 0; }
  const std::vector<char>& open_states() const { return open_; }
  double weight(std::size_t s) const { return omega_[s]; }
  const std::vector<Arc>& out(std::size_t a) const { return out_[a]; }
  const std::vector<Arc>& in(std::size_t b) const { return in_[b]; }
  std::size_t states() const { return out_.size(); }

  // ‖M^k 1‖^(1/k) for large k: the spectral radius of the transitions (the walk sums converge below one).
  double spectral_radius(int iterations = 3000) const {
    const std::size_t S = states();
    std::vector<double> x(S, 1.0), y(S);
    double log_growth = 0;
    int counted = 0;
    for (int it = 0; it < iterations; it++) {
      double top = 0;
      for (std::size_t a = 0; a < S; a++) {
        double s = 0;
        for (const auto& arc : out_[a]) s += arc.w * x[arc.to];
        y[a] = s, top = std::max(top, s);
      }
      if (top == 0) return 0;
      if (it >= iterations / 2) log_growth += std::log(top), counted++;
      for (std::size_t a = 0; a < S; a++) x[a] = y[a] / top;
    }
    return std::exp(log_growth / counted);
  }

  // z_A: the weight of every walk from each state to the target A, by Gauss–Seidel sweeps until no state changes by
  // more than `tolerance` of its own value. The walk sums span many orders of magnitude (e^{−θ ℓ} over a district:
  // 10⁻¹³ at 1.5 km), and the walkers' choices read their ratios, so each state is held to its own relative tolerance,
  // not to the largest. A state not yet reached changes by all of itself when it is, so the sweeps cannot stop early.
  std::vector<double> value(const std::vector<std::size_t>& target, double tolerance = 1e-13, int max_sweeps = 200000) const {
    const std::size_t S = states();
    std::vector<char> in_target(S, 0);
    for (const std::size_t s : target) in_target.at(s) = 1;
    std::vector<double> z(S, 0.0);
    for (std::size_t s = 0; s < S; s++)
      if (in_target[s]) z[s] = 1;
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
      double change = 0;
      for (std::size_t a = 0; a < S; a++) {
        if (in_target[a]) continue;
        double v = 0;
        for (const auto& arc : out_[a]) v += arc.w * z[arc.to];
        if (v != z[a]) change = std::max(change, std::fabs(v - z[a]) / std::fabs(v));
        z[a] = v;
      }
      if (change <= tolerance) return z;
    }
    throw std::runtime_error("urban::movement::Routes::value: the walk sums did not converge (is the spectral radius below one?)");
  }

  // π: the weight of the walks reaching each state before the target, from a source on the states (Gauss–Seidel, each
  // state to its own relative tolerance, as value()).
  std::vector<double> forward(const std::vector<double>& source, const std::vector<char>& in_target, double tolerance = 1e-13, int max_sweeps = 200000) const {
    const std::size_t S = states();
    std::vector<double> pi(source);
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
      double change = 0;
      for (std::size_t b = 0; b < S; b++) {
        double v = source[b];
        for (const auto& arc : in_[b])
          if (!in_target[arc.to]) v += pi[arc.to] * arc.w;  // arc.to is the state a the arc comes from
        if (v != pi[b]) change = std::max(change, std::fabs(v - pi[b]) / std::fabs(v));
        pi[b] = v;
      }
      if (change <= tolerance) return pi;
    }
    throw std::runtime_error("urban::movement::Routes::forward: the walk sums did not converge");
  }

  // W_oA = Σ_{b out of o} ω_b z_A(b): zero when A cannot be reached from o.
  double total(std::size_t origin, const std::vector<double>& z) const {
    double W = 0;
    for (const std::size_t b : W_->leaving[origin]) W += omega_[b] * z[b];
    return W;
  }

  // The mean visits of one trip from junction `origin` to the target to every state, f(o → A).
  std::vector<double> visits(std::size_t origin, const std::vector<std::size_t>& target, const std::vector<double>& z) const {
    const std::size_t S = states();
    std::vector<char> in_target(S, 0);
    for (const std::size_t s : target) in_target[s] = 1;
    std::vector<double> src(S, 0.0);
    const double W = total(origin, z);
    if (!(W > 0)) return std::vector<double>(S, 0.0);
    for (const std::size_t b : W_->leaving[origin]) src[b] = omega_[b] / W;
    auto pi = forward(src, in_target);
    for (std::size_t c = 0; c < S; c++) pi[c] *= z[c];
    return pi;
  }

  // The law of the first state of a trip from `origin` (states with their probabilities); empty if A cannot be reached.
  std::vector<std::pair<std::size_t, double>> first_law(std::size_t origin, const std::vector<double>& z) const {
    std::vector<std::pair<std::size_t, double>> law;
    const double W = total(origin, z);
    if (!(W > 0)) return law;
    for (const std::size_t b : W_->leaving[origin])
      if (omega_[b] * z[b] > 0) law.push_back({b, omega_[b] * z[b] / W});
    return law;
  }
  // The law of the next state after a ∉ A.
  std::vector<std::pair<std::size_t, double>> next_law(std::size_t a, const std::vector<double>& z) const {
    std::vector<std::pair<std::size_t, double>> law;
    double total = 0;
    for (const auto& arc : out_[a]) total += arc.w * z[arc.to];
    if (!(total > 0)) return law;
    for (const auto& arc : out_[a])
      if (arc.w * z[arc.to] > 0) law.push_back({arc.to, arc.w * z[arc.to] / total});
    return law;
  }

  // A draw from a law of states by nomos's logit choice on the log-weights (β = 1): none if the law is empty.
  static std::size_t draw(calx::random::Stream& rng, const std::vector<std::pair<std::size_t, double>>& law) {
    if (law.empty()) return none;
    if (law.size() == 1) return law[0].first;
    double u[16];
    std::vector<double> big;
    double* util = u;
    if (law.size() > 16) big.resize(law.size()), util = big.data();
    for (std::size_t k = 0; k < law.size(); k++) util[k] = std::log(law[k].second);
    return law[nomos::agents::choose(rng, std::span<const double>(util, law.size()), 1.0)].first;
  }

 private:
  std::shared_ptr<const Walkway> W_;
  RouteLaw L_;
  std::vector<char> open_;
  std::vector<double> omega_;
  std::vector<std::vector<Arc>> out_, in_;
};

// z for every target of a registry under one set of routes, solved once (on a pool when given one).
class Values {
 public:
  Values() = default;
  Values(std::shared_ptr<const Routes> routes, const Targets& targets, calx::scheduler::ThreadPool* pool = nullptr) : routes_(std::move(routes)) {
    const std::size_t T = targets.size();
    z_.assign(T, {});
    const auto solve = [&](std::size_t begin, std::size_t end, std::size_t) {
      for (std::size_t t = begin; t < end; t++) z_[t] = routes_->value(targets.states(t));
    };
    if (pool && T > 1) calx::scheduler::parallel_for(*pool, T, 1, solve);
    else solve(0, T, 0);
  }
  const Routes& routes() const { return *routes_; }
  std::shared_ptr<const Routes> shared_routes() const { return routes_; }
  const std::vector<double>& operator[](std::size_t target) const { return z_.at(target); }
  std::size_t size() const { return z_.size(); }

 private:
  std::shared_ptr<const Routes> routes_;
  std::vector<std::vector<double>> z_;
};

// A trip of the macroscopic flow: from a junction to a target, at a rate (trips an hour).
struct Trip {
  std::size_t origin = 0, target = 0;
  double rate = 0;
};

// The mean traversals of every state by a set of trips (an hour's, when the rates are an hour's): Σ_k rate_k f(o_k → A_k),
// one forward solve for each target.
inline std::vector<double> expected_flow(const Values& V, const Targets& T, const std::vector<Trip>& trips) {
  const Routes& R = V.routes();
  const std::size_t S = R.states();
  std::map<std::size_t, std::vector<double>> sources;
  for (const auto& t : trips) {
    if (t.rate == 0) continue;
    const auto& z = V[t.target];
    const double W = R.total(t.origin, z);
    if (!(W > 0)) continue;
    auto& src = sources[t.target];
    if (src.empty()) src.assign(S, 0.0);
    for (const std::size_t b : R.walkway().leaving[t.origin]) src[b] += t.rate * R.weight(b) / W;
  }
  std::vector<double> flow(S, 0.0);
  for (const auto& [target, src] : sources) {
    std::vector<char> in_target(S, 0);
    for (const std::size_t s : T.states(target)) in_target[s] = 1;
    const auto pi = R.forward(src, in_target);
    const auto& z = V[target];
    for (std::size_t c = 0; c < S; c++) flow[c] += pi[c] * z[c];
  }
  return flow;
}

}  // namespace urban::movement
