// Bubbles in microchannels as bits: microfluidic bubble logic (Prakash and Gershenfeld 2007)
// (substrates/fluids/README.md, simulator 8). Provisional.
//
// The model. Channels are the edges of a Kirchhoff network (physis::discretisation), each of resistance R plus
// R_b for every bubble in it (a bubble's lubricating film and menisci resist the flow); the liquid is driven by
// fixed pressures at the inlets and outlets. A bubble travels at its channel's mean speed |Q|/area along the
// flow; arriving at a junction it enters the outgoing channel of the largest outflow, the rule by which a
// bubble follows the stronger stream. A bubble that enters a channel raises its resistance and so diverts the
// next one: the nonlinearity of bubble logic. Bubbles reaching an outlet leave and are counted there.
//
// The gate. Two inlets A and B join at a junction with two outlets: one of lower resistance (A + B) and one of
// higher (A·B), with R₁ < R₂ < R₁ + R_b. A lone bubble takes the lower; two together take both, the first
// diverting the second: OR and AND in one junction. Arrivals at one instant are taken in the order of the
// channels, a declared tie-break.
//
// Limits (declared): bubbles as points with a fixed resistance; no break-up, no coalescence, no capillary number
// dependence; the junction rule a model of the stronger stream. Checks: tests/check_bubble_logic.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include <physis/discretisation/kirchhoff.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct Channel {
  std::size_t a = 0, b = 0;
  double resistance = 1;  // Pa s/m³
  double length = 1e-3;   // m
  double area = 1e-8;     // m², the cross-section
};

class BubbleCircuit {
 public:
  BubbleCircuit(std::size_t nodes, std::vector<Channel> channels, std::vector<std::pair<std::size_t, double>> pressures, double bubble_resistance)
      : n_(nodes), ch_(std::move(channels)), p_(std::move(pressures)), rb_(bubble_resistance), exits_(nodes, 0) {
    for (const auto& c : ch_)
      if (c.a >= n_ || c.b >= n_ || !(c.resistance > 0)) throw std::invalid_argument("physical_computing::fluids::BubbleCircuit: channels between nodes with positive resistance");
    for (const auto& [v, p] : p_) terminal_.push_back(v);
  }

  // A bubble enters channel c at its end `from` (a node of c).
  void inject(std::size_t c, std::size_t from) {
    if (c >= ch_.size() || (from != ch_[c].a && from != ch_[c].b)) throw std::invalid_argument("physical_computing::fluids::BubbleCircuit::inject: a channel and one of its ends");
    bubbles_.push_back({c, from, 0.0});
    injected_++;
  }

  physis::discretisation::KirchhoffSolution flows() const {
    physis::discretisation::ConductanceNetwork net;
    net.vertices = n_;
    for (std::size_t c = 0; c < ch_.size(); c++) {
      std::size_t in = 0;
      for (const auto& b : bubbles_) in += b.channel == c;
      net.add(ch_[c].a, ch_[c].b, 1 / (ch_[c].resistance + static_cast<double>(in) * rb_));
    }
    return physis::discretisation::solve_kirchhoff(net, p_);
  }

  // One step: every bubble advances with the flows at the step's start; the bubbles that reach a node are then
  // routed one at a time, in the order of their channels (the declared tie-break), each seeing the flows that
  // the bubbles routed before it have changed.
  void step(double dt) {
    const auto s = flows();
    std::vector<Bubble> placed, arriving;
    for (auto b : bubbles_) {
      const auto& c = ch_[b.channel];
      const double q = s.flow[b.channel];  // positive a → b
      // A bubble moves along the flow; against it, it stays (the flow in its channel has reversed).
      if ((b.from == c.a) == (q > 0)) b.position += std::fabs(q) / c.area * dt;
      (b.position < c.length ? placed : arriving).push_back(b);
    }
    std::stable_sort(arriving.begin(), arriving.end(), [](const Bubble& x, const Bubble& y) { return x.channel < y.channel; });
    for (const auto& b : arriving) {
      const auto& c = ch_[b.channel];
      const std::size_t node = b.from == c.a ? c.b : c.a;
      if (is_terminal(node)) {
        exits_[node]++;
        continue;
      }
      bubbles_ = placed;
      const auto now = flows();
      // The outgoing channel of the largest outflow, the arriving channel excluded; ties by order.
      std::size_t best = ch_.size();
      double best_out = 0;
      for (std::size_t k = 0; k < ch_.size(); k++) {
        if (k == b.channel) continue;
        double out = 0;
        if (ch_[k].a == node) out = now.flow[k];
        else if (ch_[k].b == node) out = -now.flow[k];
        else continue;
        if (out > best_out) best_out = out, best = k;
      }
      if (best == ch_.size()) placed.push_back({b.channel, b.from, c.length});  // held at the node
      else placed.push_back({best, node, 0.0});
    }
    bubbles_ = std::move(placed);
    t_ += dt;
  }
  void run(double duration, double dt) {
    const auto n = static_cast<std::size_t>(std::llround(duration / dt));
    for (std::size_t k = 0; k < n; k++) step(dt);
  }

  std::size_t exits(std::size_t node) const { return exits_[node]; }
  std::size_t in_flight() const { return bubbles_.size(); }
  std::size_t injected() const { return injected_; }
  std::size_t exited() const {
    std::size_t e = 0;
    for (auto x : exits_) e += x;
    return e;
  }

  // The AND–OR gate: nodes A = 0, B = 1, junction J = 2, outlets OR = 3, AND = 4; channels A→J (0), B→J (1),
  // J→OR (2), J→AND (3).
  struct Gate {
    std::size_t and_out = 0, or_out = 0;
  };
  static BubbleCircuit and_or_gate(double inlet_pressure = 1000) {
    const double R = 1e12, Rb = 0.8e12;
    std::vector<Channel> ch{{0, 2, R, 1e-3, 1e-8}, {1, 2, R, 1e-3, 1e-8}, {2, 3, 1.0 * R, 1e-3, 1e-8}, {2, 4, 1.5 * R, 1e-3, 1e-8}};
    return BubbleCircuit(5, ch, {{0, inlet_pressure}, {1, inlet_pressure}, {3, 0.0}, {4, 0.0}}, Rb);
  }
  static Gate evaluate(bool a, bool b) {
    auto g = and_or_gate();
    if (a) g.inject(0, 0);
    if (b) g.inject(1, 1);
    g.run(10.0, 1e-3);
    return {g.exits(4), g.exits(3)};
  }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "bubble logic";
    d.accepts = {tasks::Kind::function};
    d.inputs = {{"bubbles at the inlets", "1"}};
    d.outputs = {{"bubbles at the outlets", "1"}};
    d.materials.materials = {{"water", "the carrier liquid"}, {"air", "the bubbles"}, {"silicone (PDMS)", "the channels"}};
    d.model = "channels as a Kirchhoff network, bubbles as added resistances following the stronger stream (provisional; Prakash and Gershenfeld 2007)";
    d.limits = {"bubbles as points of fixed resistance", "no break-up or coalescence", "the junction rule a model"};
    return d;
  }

 private:
  struct Bubble {
    std::size_t channel, from;
    double position;
  };
  bool is_terminal(std::size_t v) const {
    for (auto t : terminal_)
      if (t == v) return true;
    return false;
  }
  std::size_t n_;
  std::vector<Channel> ch_;
  std::vector<std::pair<std::size_t, double>> p_;
  std::vector<std::size_t> terminal_;
  double rb_;
  std::vector<Bubble> bubbles_;
  std::vector<std::size_t> exits_;
  std::size_t injected_ = 0;
  double t_ = 0;
};

}  // namespace physical_computing::fluids
