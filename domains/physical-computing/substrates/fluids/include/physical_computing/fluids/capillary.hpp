// Networks of capillaries: slow flow in tubes, and the solving of mazes by it (Fuerstman et al. 2003)
// (substrates/fluids/README.md, simulator 6).
//
// A network of tubes (radius r, length L) carries a slow flow, each tube a Hagen–Poiseuille conductance
// π r⁴/(8 μ L) (physis::continuum::poiseuille_conductance) in a Kirchhoff network (physis::discretisation),
// driven by a pressure difference between terminals; a liquid drawn in by its meniscus alone fills a tube in
// Washburn's time 2 μ L²/(γ r cos θ). In a perfect maze (a spanning tree of a grid of cells, every two cells
// joined by exactly one path) driven from its entrance to its exit, a dead end carries no flow, since a branch
// that leads nowhere cannot conserve a non-zero flow: the flowing tubes are the maze's solution, read off at the
// cost of one linear solve (Fuerstman et al. 2003 dyed the flowing channels of real microfluidic mazes).
//
// Limits (declared): steady, laminar, fully developed flow; no entrance effects, no gravity; menisci only in
// Washburn's filling. Checks: tests/check_capillary.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/random/stream.hpp>
#include <physis/continuum/capillarity.hpp>
#include <physis/discretisation/kirchhoff.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct Tube {
  std::size_t a = 0, b = 0;
  double radius = 1e-4, length = 1e-2;  // m
};

class CapillaryNetwork {
 public:
  CapillaryNetwork(std::size_t vertices, std::vector<Tube> tubes, Fluid fluid = Fluid::water()) : n_(vertices), tubes_(std::move(tubes)), fluid_(std::move(fluid)) {
    for (const auto& t : tubes_)
      if (t.a >= n_ || t.b >= n_) throw std::invalid_argument("physical_computing::fluids::CapillaryNetwork: tubes between vertices");
  }

  physis::discretisation::ConductanceNetwork conductances() const {
    physis::discretisation::ConductanceNetwork net;
    net.vertices = n_;
    // The viscosity is the dynamic one, μ = ρ ν.
    const double mu = fluid_.density * fluid_.kinematic_viscosity;
    for (const auto& t : tubes_) net.add(t.a, t.b, physis::continuum::poiseuille_conductance(t.radius, t.length, mu));
    return net;
  }
  // The volume flow (m³/s) in each tube, positive from a to b, with Δp (Pa) from inlet to outlet.
  physis::discretisation::KirchhoffSolution flow(std::size_t inlet, std::size_t outlet, double pressure) const {
    return physis::discretisation::solve_kirchhoff(conductances(), {{inlet, pressure}, {outlet, 0.0}});
  }
  // Washburn's time to fill each tube by capillarity alone, at contact angle θ.
  std::vector<double> filling_times(double contact_angle = 0) const {
    const double mu = fluid_.density * fluid_.kinematic_viscosity;
    std::vector<double> t;
    for (const auto& tube : tubes_) t.push_back(physis::continuum::washburn_time(tube.length, tube.radius, fluid_.surface_tension, contact_angle, mu));
    return t;
  }

  std::size_t vertices() const { return n_; }
  const std::vector<Tube>& tubes() const { return tubes_; }

  // A perfect maze on a width × height grid of cells, by a randomised depth-first search from cell 0 (calx's
  // keyed stream): a spanning tree of the grid, each passage a tube of the given radius and length.
  static CapillaryNetwork maze(std::size_t width, std::size_t height, std::uint64_t seed, double radius = 1e-4, double length = 2e-3) {
    if (width < 2 || height < 2) throw std::invalid_argument("physical_computing::fluids::CapillaryNetwork::maze: at least 2 × 2 cells");
    calx::random::Stream rng(seed, 0x3A2Eu, 0, 0);
    const std::size_t n = width * height;
    std::vector<char> seen(n, 0);
    std::vector<std::size_t> stack{0};
    seen[0] = 1;
    std::vector<Tube> tubes;
    while (!stack.empty()) {
      const std::size_t c = stack.back(), x = c % width, y = c / width;
      std::vector<std::size_t> next;
      if (x > 0 && !seen[c - 1]) next.push_back(c - 1);
      if (x + 1 < width && !seen[c + 1]) next.push_back(c + 1);
      if (y > 0 && !seen[c - width]) next.push_back(c - width);
      if (y + 1 < height && !seen[c + width]) next.push_back(c + width);
      if (next.empty()) {
        stack.pop_back();
        continue;
      }
      const std::size_t d = next[static_cast<std::size_t>(rng.below(next.size()))];
      seen[d] = 1;
      tubes.push_back({c, d, radius, length});
      stack.push_back(d);
    }
    return CapillaryNetwork(n, std::move(tubes));
  }

  // The tubes that carry flow, |Q| > tolerance × max |Q|: in a perfect maze, its solution.
  std::vector<std::size_t> flowing(const physis::discretisation::KirchhoffSolution& s, double tolerance = 1e-9) const {
    double peak = 0;
    for (double q : s.flow) peak = std::max(peak, std::fabs(q));
    std::vector<std::size_t> out;
    for (std::size_t k = 0; k < s.flow.size(); k++)
      if (std::fabs(s.flow[k]) > tolerance * peak) out.push_back(k);
    return out;
  }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "capillary network";
    d.accepts = {tasks::Kind::function};
    d.inputs = {{"pressure at the inlet", "Pa"}};
    d.outputs = {{"flow in each tube", "m³/s"}};
    d.materials.materials = {{fluid_.name, "the liquid"}, {"silicone (PDMS)", "the channels"}, {"glass", "the cover"}};
    d.materials.reuse = computer::Reuse::reusable;
    d.model = "Hagen–Poiseuille conductances in a Kirchhoff network (physis); Washburn's filling; Fuerstman et al. 2003";
    d.limits = {"steady laminar flow", "no entrance effects or gravity", "menisci only in filling"};
    return d;
  }

 private:
  std::size_t n_;
  std::vector<Tube> tubes_;
  Fluid fluid_;
};

}  // namespace physical_computing::fluids
