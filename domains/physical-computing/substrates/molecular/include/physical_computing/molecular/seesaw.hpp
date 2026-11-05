// Seesaw gates: digital circuits of DNA strand displacement (Qian and Winfree 2011, Science 332, 1196)
// (substrates/molecular/README.md, simulator 2). Provisional.
//
// A seesaw gate holds one strand and exchanges it for another by reversible displacement, w_in + G:out ⇌ G:in + w_out
// at a rate k_s each way; a threshold complex absorbs an input fast and irreversibly, w_in + Th → waste at k_f ≫ k_s;
// and a fuel strand in excess drives the exchange back, w_fuel + G:in ⇌ G:fuel + w_in, so that each input molecule
// left over the threshold releases many outputs: amplification with a threshold, the sigmoid of digital logic.
//
// Two-input logic (Qian and Winfree 2011): an integrating gate turns each input into the same signal s, so
// that s = x₁ + x₂; an amplifying gate with threshold θ and a gate of output releases y ≈ its full load when
// s > θ and ≈ 0 when s < θ. With inputs of 1 unit for "on", θ = 0.6 makes OR and θ = 1.2 makes AND.
//
// Limits (declared): signals as single species; a reporter reads the output directly; rates as declared
// constants, not the toeholds' energies. Checks: tests/check_seesaw.cpp.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <morphe/arithmetic/rational.hpp>
#include <physis/chemistry/kinetics.hpp>
#include <physis/dynamics/flows.hpp>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::molecular {

struct SeesawGate {
  physis::chemistry::Network network;
  std::vector<double> initial;
  std::size_t x1 = 0, x2 = 0, y = 0;
};

// A two-input gate of threshold θ (0.6: OR; 1.2: AND), inputs "on" at 1 unit; gate load 1, fuel 2, as in the paper.
inline SeesawGate two_input_gate(double theta, double ks = 1, double kf = 100) {
  using morphe::arithmetic::Integer;
  using morphe::arithmetic::Rational;
  auto R = [](double v) { return Rational(Integer(static_cast<long long>(v * 1000 + 0.5)), Integer(1000)); };
  SeesawGate g;
  auto& n = g.network;
  n.species = {"x1", "x2", "s", "Gi1", "Gi2", "Gi1:x1", "Gi2:x2", "Th", "Ga:y", "Ga:s", "y", "fuel", "Ga:fuel", "waste"};
  // The integrating gate: each input frees one s (irreversible: its gates hold s, and the inputs are absorbed).
  n.add("x1 frees s", {{"x1", 1}, {"Gi1", 1}}, {{"Gi1:x1", 1}, {"s", 1}}, R(ks));
  n.add("x2 frees s", {{"x2", 1}, {"Gi2", 1}}, {{"Gi2:x2", 1}, {"s", 1}}, R(ks));
  // The threshold absorbs s fast.
  n.add("threshold", {{"s", 1}, {"Th", 1}}, {{"waste", 1}}, R(kf));
  // The amplifying seesaw: s + Ga:y ⇌ Ga:s + y; fuel + Ga:s ⇌ Ga:fuel + s.
  n.add("s releases y", {{"s", 1}, {"Ga:y", 1}}, {{"Ga:s", 1}, {"y", 1}}, R(ks));
  n.add("y returns", {{"Ga:s", 1}, {"y", 1}}, {{"s", 1}, {"Ga:y", 1}}, R(ks));
  n.add("fuel frees s", {{"fuel", 1}, {"Ga:s", 1}}, {{"Ga:fuel", 1}, {"s", 1}}, R(ks));
  n.add("s takes fuel's place", {{"Ga:fuel", 1}, {"s", 1}}, {{"fuel", 1}, {"Ga:s", 1}}, R(ks));
  g.initial.assign(n.species.size(), 0.0);
  g.initial[n.find("Gi1")] = 2, g.initial[n.find("Gi2")] = 2;
  g.initial[n.find("Th")] = theta;
  g.initial[n.find("Ga:y")] = 1, g.initial[n.find("fuel")] = 2;
  g.x1 = n.find("x1"), g.x2 = n.find("x2"), g.y = n.find("y");
  return g;
}

// Runs the gate on inputs (each 0 or 1 unit) and returns the output's concentration at time t.
inline double run_gate(SeesawGate g, bool a, bool b, double t = 200, double h = 1e-3) {
  physis::dynamics::State c(g.initial.begin(), g.initial.end());
  c[g.x1] = a ? 1 : 0, c[g.x2] = b ? 1 : 0;
  auto f = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { physis::chemistry::mass_action(g.network, y, dy); };
  physis::dynamics::integrate(f, 0.0, t, c, h);
  return c[g.y];
}

inline runtime::TargetDescription describe_seesaw() {
  runtime::TargetDescription d;
  d.family = "molecular", d.substrate = "seesaw circuits";
  d.accepts = {tasks::Kind::function};
  d.inputs = {{"input strands", "units of 1 for on"}};
  d.outputs = {{"output strand (by a reporter)", "units"}};
  d.materials.materials = {{"DNA", "gates, thresholds, fuels"}, {"magnesium salts", "the buffer"}};
  d.materials.reuse = computer::Reuse::consumed_per_run;
  d.materials.end_of_life = computer::EndOfLife::biodegradable;
  d.model = "seesaw gates by mass action (Qian and Winfree 2011)";
  d.limits = {"signals as single species", "declared rates", "a direct reporter"};
  return d;
}

}  // namespace physical_computing::molecular
