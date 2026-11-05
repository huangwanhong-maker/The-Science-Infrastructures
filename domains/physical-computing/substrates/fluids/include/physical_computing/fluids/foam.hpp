// A two-dimensional foam coarsening: cells that grow or shrink by their number of sides (von Neumann 1952)
// (substrates/fluids/README.md, simulator 7).
//
// A dry foam's cells exchange gas through their films at a rate set by the films' curvature; in two dimensions
// a cell of n sides changes its area at dA/dt = (π/3) M γ (n − 6) whatever its shape (von Neumann 1952;
// Mullins 1956): cells of fewer than six sides shrink and vanish, cells of more grow, and the mean area grows
// in proportion to time. The model is the large-q Potts model with neighbour proposals at a small temperature
// (physis::statistical_mechanics::PottsLattice; Glazier, Anderson and Grest 1990), each Potts state a cell,
// started from a periodic Voronoi tessellation of random seeds. Its time is Monte Carlo sweeps; its law of
// coarsening is statistical: the mean dA/dt of the cells of each number of sides is linear in n and crosses
// zero near n = 6.
//
// Limits (declared): a lattice model, with the lattice's anisotropy reduced by the Moore neighbourhood and a
// small temperature, not removed; films of zero thickness; no drainage. Checks: tests/check_foam.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

#include <calx/random/stream.hpp>
#include <physis/statistical_mechanics/potts.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct FoamOptions {
  std::size_t size = 128;       // sites per side
  std::size_t cells = 400;      // initial cells
  double temperature = 0.3;     // in units of the film energy J
  std::uint64_t seed = 1;
};

class Foam {
 public:
  explicit Foam(const FoamOptions& o)
      : o_(o), lattice_(o.size, o.size, o.cells, 1.0, physis::statistical_mechanics::Neighbourhood::moore), rng_(o.seed, 0xF0A3u, 1, 0) {
    if (o.cells < 4 || o.cells > o.size * o.size / 4) throw std::invalid_argument("physical_computing::fluids::Foam: between 4 cells and a quarter of the sites");
    // A periodic Voronoi tessellation of random seeds.
    std::vector<std::pair<double, double>> seeds(o.cells);
    for (auto& s : seeds) s = {rng_.uniform() * o.size, rng_.uniform() * o.size};
    const double L = static_cast<double>(o.size);
    for (std::size_t y = 0; y < o.size; y++)
      for (std::size_t x = 0; x < o.size; x++) {
        std::size_t best = 0;
        double bd = INFINITY;
        for (std::size_t c = 0; c < seeds.size(); c++) {
          double dx = std::fabs(x + 0.5 - seeds[c].first), dy = std::fabs(y + 0.5 - seeds[c].second);
          dx = std::min(dx, L - dx), dy = std::min(dy, L - dy);
          const double d2 = dx * dx + dy * dy;
          if (d2 < bd) bd = d2, best = c;
        }
        lattice_.set_state(lattice_.site(x, y), best);
      }
  }

  void coarsen(std::size_t sweeps) {
    for (std::size_t s = 0; s < sweeps; s++) lattice_.sweep(o_.temperature, physis::statistical_mechanics::PottsProposal::neighbour, rng_);
    sweeps_ += sweeps;
  }
  std::size_t sweeps() const { return sweeps_; }
  const physis::statistical_mechanics::PottsLattice& lattice() const { return lattice_; }

  // Each living cell (by its Potts state): its area and its number of sides, those of its largest piece.
  struct Cell {
    std::size_t area = 0, sides = 0;
  };
  std::map<std::size_t, Cell> cells() const {
    const auto d = lattice_.domains();
    std::map<std::size_t, Cell> out;
    std::map<std::size_t, std::size_t> largest;
    for (std::size_t c = 0; c < d.area.size(); c++) {
      auto& cell = out[d.state[c]];
      cell.area += d.area[c];
      if (!largest.count(d.state[c]) || d.area[c] > d.area[largest[d.state[c]]]) largest[d.state[c]] = c, cell.sides = d.sides[c];
    }
    return out;
  }
  double mean_area() const {
    const auto c = cells();
    return static_cast<double>(o_.size * o_.size) / static_cast<double>(c.size());
  }

  // von Neumann's law measured: the mean rate of change of area, in sites per sweep, of the cells of each number
  // of sides at the first snapshot. A cell gone by the second counts with area zero, a rate that understates how
  // fast it shrank but not its sign: leaving it out would drop exactly the cells that shrink fastest.
  static std::map<std::size_t, std::pair<double, std::size_t>> area_rates(const std::map<std::size_t, Cell>& before, const std::map<std::size_t, Cell>& after, double sweeps) {
    std::map<std::size_t, std::pair<double, std::size_t>> rate;
    for (const auto& [id, c] : before) {
      const auto it = after.find(id);
      const double end = it == after.end() ? 0.0 : static_cast<double>(it->second.area);
      auto& r = rate[c.sides];
      r.first += (end - static_cast<double>(c.area)) / sweeps;
      r.second++;
    }
    for (auto& [n, r] : rate) r.first /= static_cast<double>(r.second);
    return rate;
  }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "foam";
    d.accepts = {tasks::Kind::phenomenon};
    d.inputs = {{"the initial cells", "1"}};
    d.outputs = {{"cells' areas and sides", "1"}};
    d.materials.materials = {{"water", "the films"}, {"soap", "the surfactant"}, {"air", "the cells"}};
    d.materials.reuse = computer::Reuse::consumed_per_run;
    d.materials.end_of_life = computer::EndOfLife::biodegradable;
    d.model = "the large-q Potts model of coarsening (physis; Glazier, Anderson and Grest 1990); von Neumann's law";
    d.limits = {"a lattice model, anisotropy reduced not removed", "films of zero thickness", "no drainage"};
    return d;
  }

 private:
  FoamOptions o_;
  physis::statistical_mechanics::PottsLattice lattice_;
  calx::random::Stream rng_;
  std::size_t sweeps_ = 0;
};

}  // namespace physical_computing::fluids
