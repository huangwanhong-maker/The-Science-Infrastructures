// Cells of the BZ medium coupled on a graph: the programmable chemical computer of Parrilla-Gutierrez et al.
// (2020) (substrates/chemical/README.md, simulator 2).
//
// Each cell is a stirred Oregonator (oregonator.hpp); neighbouring cells exchange the activator through their
// walls at a declared coupling (the Laplacian of the graph of cells, physis::chemistry::SparseLaplacian), and
// each cell's light φ (or, in the 2020 computer, its stirring) is its own input. The cells' phases, read at the
// peaks of their catalyst, are the array's state: coupled oscillators that lock to one another or not.
//
// Limits (declared): a cell as one well-stirred Oregonator; coupling through the activator alone (the bromine of
// real droplets inhibits instead, Toiya, Vanag and Epstein 2008, and can favour anti-phase); inputs as light.
// Checks: tests/check_bz_array.cpp.
#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <physis/chemistry/field_kinetics.hpp>

#include <physical_computing/chemical/oregonator.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::chemical {

class BZArray {
 public:
  // The cells' kinetics (each its own f and φ), and their couplings (a, b, strength): physis's FieldKinetics on the
  // graph of cells, the activator crossing the walls (D_u = 1 times the couplings), the catalyst not (D_v = 0),
  // stepped implicitly in the kinetics (they are stiff in the refractory phase).
  BZArray(std::vector<Oregonator> cells, const std::vector<std::tuple<std::size_t, std::size_t, double>>& couplings)
      : cells_(std::move(cells)),
        field_(2, physis::chemistry::FieldKinetics::VertexKinetics([this](std::size_t v, const double* c, double* dc) { cells_[v].kinetics(c, dc); }),
               physis::chemistry::SparseLaplacian::graph(cells_.size(), couplings), {1.0, 0.0}) {
    for (std::size_t v = 0; v < cells_.size(); v++) {
      const double u = cells_[v].steady_u();
      field_.state()[2 * v] = u * (1.1 + 0.1 * static_cast<double>(v));  // every cell disturbed from its steady state
      field_.state()[2 * v + 1] = u;
    }
    field_.set_non_negative();  // concentrations
  }
  BZArray(const BZArray&) = delete;  // the kinetics refer to this array's cells
  BZArray& operator=(const BZArray&) = delete;

  std::size_t cells() const { return cells_.size(); }
  const std::vector<double>& state() const { return field_.state(); }
  std::vector<double>& state() { return field_.state(); }
  Oregonator& cell(std::size_t v) { return cells_[v]; }
  double time() const { return field_.time(); }

  void step(double h) { field_.step_split(h); }
  // Runs for `duration`, returning each cell's times of the catalyst's peaks.
  std::vector<std::vector<double>> run(double duration, double h) {
    std::vector<std::vector<double>> peaks(cells_.size());
    std::vector<double> prev2(cells_.size()), prev1(cells_.size());
    for (std::size_t v = 0; v < cells_.size(); v++) prev2[v] = prev1[v] = field_.state()[2 * v + 1];
    const auto steps = static_cast<std::size_t>(std::llround(duration / h));
    for (std::size_t k = 0; k < steps; k++) {
      step(h);
      for (std::size_t v = 0; v < cells_.size(); v++) {
        const double now = field_.state()[2 * v + 1];
        if (prev1[v] > prev2[v] && prev1[v] >= now) peaks[v].push_back(field_.time() - h);
        prev2[v] = prev1[v], prev1[v] = now;
      }
    }
    return peaks;
  }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "chemical", d.substrate = "BZ array";
    d.accepts = {tasks::Kind::phenomenon, tasks::Kind::filter, tasks::Kind::function};
    for (std::size_t v = 0; v < cells_.size(); v++) d.inputs.push_back({"light on cell " + std::to_string(v), "1"});
    for (std::size_t v = 0; v < cells_.size(); v++) d.outputs.push_back({"phase of cell " + std::to_string(v), "rad"});
    d.materials.materials = {{"malonic acid", "substrate"}, {"sodium bromate", "oxidant"}, {"sulfuric acid", "acid"}, {"ferroin", "catalyst"}};
    d.materials.reuse = computer::Reuse::consumed_per_run;
    d.model = "stirred Oregonators coupled through the activator on a graph (Parrilla-Gutierrez et al. 2020)";
    d.limits = {"one stirred Oregonator per cell", "coupling by the activator alone", "inputs as light"};
    return d;
  }

 private:
  std::vector<Oregonator> cells_;
  physis::chemistry::FieldKinetics field_;
};

}  // namespace physical_computing::chemical
