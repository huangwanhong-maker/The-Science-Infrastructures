// Reaction networks as programs (substrates/chemical/README.md, simulator 3).
//
// Functions on counts. A network computes y = g(x) stably when, from any start with xᵢ molecules of the inputs
// and none of the rest, every sequence of reactions ends where no reaction can occur with exactly g(x) molecules of
// the output. The deterministic networks compute exactly the semilinear functions (Chen, Doty and Soloveichik 2014):
// sums, multiples, minima, maxima, floors of quotients by constants, and their compositions. Each program here is a
// physis::chemistry::Network; it is run exactly by its Markov chain's law of absorption (physis's absorption, in
// rationals), so that "computes" means probability one, not most runs.
//
// Flows. A polynomial ODE, a task of kind flow, is lowered to a network by physis's dual rails, its encoding
// E(x) = (x⁺, x⁻) and readout R(p, q) = p − q, a fixed linear map: a non-trivial square (theory note §1.1),
// checked by computer::commutes against the ODE integrated directly.
//
// Limits (declared): ideal networks, with no leaks; realised by DNA in the molecular family.
// Checks: tests/check_programs.cpp.
#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include <physis/chemistry/dual_rail.hpp>
#include <physis/chemistry/kinetics.hpp>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::chemical {

using physis::chemistry::Network;
using physis::chemistry::State;

struct Program {
  std::string name;
  Network network;
  std::vector<std::string> inputs;
  std::string output;
};

inline Program sum_program() {
  Program p{"x₁ + x₂", {}, {"X1", "X2"}, "Y"};
  p.network.species = {"X1", "X2", "Y"};
  p.network.add("X1 → Y", {{"X1", 1}}, {{"Y", 1}}, 1);
  p.network.add("X2 → Y", {{"X2", 1}}, {{"Y", 1}}, 1);
  return p;
}
inline Program double_program() {
  Program p{"2x", {}, {"X"}, "Y"};
  p.network.species = {"X", "Y"};
  p.network.add("X → 2Y", {{"X", 1}}, {{"Y", 2}}, 1);
  return p;
}
inline Program half_program() {
  Program p{"⌊x/2⌋", {}, {"X"}, "Y"};
  p.network.species = {"X", "Y"};
  p.network.add("2X → Y", {{"X", 2}}, {{"Y", 1}}, 1);
  return p;
}
inline Program min_program() {
  Program p{"min(x₁, x₂)", {}, {"X1", "X2"}, "Y"};
  p.network.species = {"X1", "X2", "Y"};
  p.network.add("X1 + X2 → Y", {{"X1", 1}, {"X2", 1}}, {{"Y", 1}}, 1);
  return p;
}
// max(x₁, x₂) = x₁ + x₂ − min(x₁, x₂) (Chen, Doty and Soloveichik 2014).
inline Program max_program() {
  Program p{"max(x₁, x₂)", {}, {"X1", "X2"}, "Y"};
  p.network.species = {"X1", "X2", "Y", "A", "B", "K"};
  p.network.add("X1 → Y + A", {{"X1", 1}}, {{"Y", 1}, {"A", 1}}, 1);
  p.network.add("X2 → Y + B", {{"X2", 1}}, {{"Y", 1}, {"B", 1}}, 1);
  p.network.add("A + B → K", {{"A", 1}, {"B", 1}}, {{"K", 1}}, 1);
  p.network.add("Y + K → ∅", {{"Y", 1}, {"K", 1}}, {}, 1);
  return p;
}

// The exact law of the output when the program's chain is absorbed, from the given input counts.
inline std::map<long long, morphe::arithmetic::Rational> run_exactly(const Program& p, const std::vector<long long>& counts) {
  State n(p.network.species.size(), 0);
  for (std::size_t k = 0; k < p.inputs.size(); k++) n[p.network.find(p.inputs[k])] = counts[k];
  const auto law = physis::chemistry::absorption(p.network, n);
  const std::size_t out = p.network.find(p.output);
  std::map<long long, morphe::arithmetic::Rational> result;
  for (const auto& [state, prob] : law) result[state[out]] = result[state[out]] + prob;
  return result;
}

// A flow lowered to a network by dual rails.
inline physis::chemistry::DualRailNetwork lower_flow(const physis::chemistry::PolynomialSystem& system, const morphe::arithmetic::Rational& annihilation = 1) {
  return physis::chemistry::dual_rail(system, annihilation);
}

inline runtime::TargetDescription describe_programs() {
  runtime::TargetDescription d;
  d.family = "chemical", d.substrate = "reaction networks";
  d.accepts = {tasks::Kind::function, tasks::Kind::flow};
  d.inputs = {{"input species' counts or concentrations", "1"}};
  d.outputs = {{"the output species", "1"}};
  d.model = "networks of mass action and their Markov chains (physis::chemistry); dual rails (Fages et al. 2017)";
  d.limits = {"ideal networks without leaks", "realised by the molecular family"};
  return d;
}

}  // namespace physical_computing::chemical
