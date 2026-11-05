// Any formal reaction network, realised by DNA strand displacement (Soloveichik, Seelig and Winfree 2010)
// (substrates/molecular/README.md, simulator 1). Provisional until molecular-science is founded (ADR 0076 §2.5).
//
// The compilation. Each formal reaction r of rate k becomes cascades of strand displacement among fuel complexes
// held at a large concentration C_max, the formal species at concentrations far below it:
//
//   unimolecular, X → P:       X + G_r → O_r            at q_r = k / C_max,
//                              O_r + T_r → P            at q_max;
//   bimolecular,  X + Y → P:   X + L_r ⇌ H_r + B_r      at q_r = k forward, q_max back,
//                              Y + H_r → O_r            at q_max,
//                              O_r + T_r → P            at q_max,
//
// with G_r, T_r, L_r and B_r starting at C_max (wastes, inert, are left out). While the fuels stay near C_max and
// q_max C_max is fast against the formal rates, the formal species obey the formal network's mass action: for the
// unimolecular cascade at once (q_r C_max = k), for the bimolecular through H_r's quasi-steady state
// H_r ≈ (q_r/q_max) X L_r/B_r ≈ (k/q_max) X (Soloveichik et al. 2010). The fuels are consumed: each firing of
// r uses one G_r (or L_r) and one T_r, so the realisation runs down as its fuel does, and its error grows with
// the formal species' share of C_max.
//
// A leak (declared): a unimolecular gate also releases its output without its input, G_r + T_r → P at a declared
// rate, the leading defect of real strand displacement.
//
// Units: concentrations in a unit of the caller's (nM, M), rates consistently. Checks: tests/check_strand_displacement.cpp.
#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/arithmetic/rational.hpp>
#include <physis/chemistry/kinetics.hpp>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::molecular {

using morphe::arithmetic::Rational;
using physis::chemistry::Network;
using physis::chemistry::Reaction;

struct StrandDisplacement {
  Network network;                        // the realisation, formal species first, then the fuels
  std::size_t formal_species = 0;
  std::vector<std::size_t> consumed_fuel; // for each formal reaction, the fuel it uses up (G_r or L_r)
  std::vector<double> initial;            // the fuels at C_max, the formal species at 0 (set them)
};

// Compiles a formal network of uni- and bimolecular reactions (each with its formal rate).
inline StrandDisplacement compile_strand_displacement(const Network& formal, const Rational& c_max, const Rational& q_max, const Rational& leak = 0) {
  StrandDisplacement d;
  d.network.species = formal.species;
  d.formal_species = formal.species.size();
  auto add_species = [&](const std::string& name) {
    d.network.species.push_back(name);
    return d.network.species.size() - 1;
  };
  for (std::size_t r = 0; r < formal.reactions.size(); r++) {
    const auto& R = formal.reactions[r];
    const auto order = R.order();
    const std::string tag = std::to_string(r);
    if (order == 1) {
      const std::size_t X = R.reactants[0].first;
      const std::size_t G = add_species("G" + tag), O = add_species("O" + tag), T = add_species("T" + tag);
      Reaction a;
      a.name = R.name + ": input binds its gate";
      a.reactants = {{X, 1}, {G, 1}}, a.products = {{O, 1}}, a.rate = R.rate / c_max;
      Reaction b;
      b.name = R.name + ": the output released";
      b.reactants = {{O, 1}, {T, 1}}, b.products = R.products, b.rate = q_max;
      d.network.reactions.push_back(a), d.network.reactions.push_back(b);
      if (Rational(0) < leak) {
        Reaction l;
        l.name = R.name + ": leak";
        l.reactants = {{G, 1}, {T, 1}}, l.products = R.products, l.rate = leak;
        d.network.reactions.push_back(l);
      }
      d.consumed_fuel.push_back(G);
    } else if (order == 2) {
      std::size_t X, Y;
      if (R.reactants.size() == 1) X = Y = R.reactants[0].first;  // 2X → P
      else X = R.reactants[0].first, Y = R.reactants[1].first;
      const std::size_t L = add_species("L" + tag), H = add_species("H" + tag), B = add_species("B" + tag), O = add_species("O" + tag), T = add_species("T" + tag);
      Reaction a;
      a.name = R.name + ": the first input displaces";
      a.reactants = {{X, 1}, {L, 1}}, a.products = {{H, 1}, {B, 1}}, a.rate = R.rate;
      Reaction a2;
      a2.name = R.name + ": and is displaced back";
      a2.reactants = {{H, 1}, {B, 1}}, a2.products = {{X, 1}, {L, 1}}, a2.rate = q_max;
      Reaction b;
      b.name = R.name + ": the second input commits";
      b.reactants = {{Y, 1}, {H, 1}}, b.products = {{O, 1}}, b.rate = q_max;
      Reaction c;
      c.name = R.name + ": the output released";
      c.reactants = {{O, 1}, {T, 1}}, c.products = R.products, c.rate = q_max;
      d.network.reactions.push_back(a), d.network.reactions.push_back(a2), d.network.reactions.push_back(b), d.network.reactions.push_back(c);
      d.consumed_fuel.push_back(L);
    } else {
      throw std::invalid_argument("physical_computing::molecular::compile_strand_displacement: reactions of order one or two");
    }
  }
  d.initial.assign(d.network.species.size(), 0.0);
  for (std::size_t s = d.formal_species; s < d.network.species.size(); s++) {
    const char k = d.network.species[s][0];
    if (k == 'G' || k == 'T' || k == 'L' || k == 'B') d.initial[s] = c_max.to_double();
  }
  return d;
}

inline runtime::TargetDescription describe_strand_displacement() {
  runtime::TargetDescription d;
  d.family = "molecular", d.substrate = "DNA strand displacement";
  d.accepts = {tasks::Kind::function, tasks::Kind::flow};
  d.inputs = {{"input strands", "nM"}};
  d.outputs = {{"output strands (by fluorescent reporters)", "nM"}};
  d.materials.materials = {{"DNA", "the strands and gates"}, {"magnesium salts", "the buffer"}, {"phosphate rock", "the phosphate of DNA's synthesis"}};
  d.materials.reuse = computer::Reuse::consumed_per_run;
  d.materials.end_of_life = computer::EndOfLife::biodegradable;
  d.model = "the compilation of Soloveichik, Seelig and Winfree (2010), by mass action (physis::chemistry)";
  d.limits = {"wastes left out", "leaks one declared rate", "no secondary structure or toehold kinetics"};
  return d;
}

}  // namespace physical_computing::molecular
