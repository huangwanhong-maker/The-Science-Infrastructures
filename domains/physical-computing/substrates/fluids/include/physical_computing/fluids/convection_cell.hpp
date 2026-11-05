// A convection cell: a layer of fluid heated from below, whose rolls carry Lorenz's chaos
// (substrates/fluids/README.md, simulator 2).
//
// The cell. A layer of depth d between plates (rigid, or free), heated by ΔT, has the Rayleigh number
// Ra = g α ΔT d³/(ν κ) and the Prandtl number Pr = ν/κ; it convects above the onset Ra_c of its plates (physis's
// shooting: 1707.762 between rigid plates), so at ΔT_c = Ra_c ν κ/(g α d³).
//
// Its truncation. Lorenz's three modes (physis::continuum::lorenz_parameters) with b = 4/(1 + a_L²) for the
// rolls' wavenumber and σ = Pr. Lorenz derived them between free plates; the cell's r is taken relative to its
// own plates' onset, r = Ra/Ra_c(plates), so that the truncation convects when the cell does (a declared
// choice: `scaling`). Its time is τ = π²(1 + a_L²) κ t / d². The input is ΔT (K), the control of the
// heating; the readings are the rolls' circulation X, the temperature modes Y and Z, and the Nusselt number
// 1 + 2Z/r. A phenomenon of the cell is the subcritical Hopf bifurcation of its convecting states at
// r_H = σ(σ + b + 3)/(σ − b − 1) (24.74 for σ = 10, b = 8/3), beyond which it is chaotic.
//
// Limits (declared): three modes; rolls in two dimensions; the flow at large r is three-dimensional turbulence,
// which the truncation does not describe; Lorenz's modes are those of free plates. Checks: tests/check_convection_cell.cpp.
#pragma once

#include <cmath>
#include <cstdint>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

#include <physis/continuum/convection.hpp>
#include <physis/dynamics/flows.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

enum class LorenzScaling { cell_onset, free_plates };

struct ConvectionCellOptions {
  double depth = 0.01;                 // d, m
  double temperature_difference = 0.5; // ΔT, K
  Fluid fluid = Fluid::water();
  physis::continuum::Plate bottom = physis::continuum::Plate::rigid, top = physis::continuum::Plate::rigid;
  double roll_wavenumber = std::numbers::pi / std::numbers::sqrt2;  // a (times d)
  LorenzScaling scaling = LorenzScaling::cell_onset;
  double step = 0.002;                 // dτ of the truncation
};

class ConvectionCell final : public runtime::Backend {
 public:
  explicit ConvectionCell(const ConvectionCellOptions& o) : o_(o), onset_(physis::continuum::critical_rayleigh(o.bottom, o.top)) {
    if (!(o.depth > 0) || !(o.step > 0)) throw std::invalid_argument("physical_computing::fluids::ConvectionCell: a positive depth and step");
    reset(0);
  }

  double rayleigh(double dT) const {
    const auto& f = o_.fluid;
    return standard_gravity * f.thermal_expansion * dT * o_.depth * o_.depth * o_.depth / (f.kinematic_viscosity * f.thermal_diffusivity);
  }
  double rayleigh() const { return rayleigh(dT_); }
  double prandtl() const { return o_.fluid.kinematic_viscosity / o_.fluid.thermal_diffusivity; }
  const physis::continuum::Onset& onset() const { return onset_; }
  double critical_temperature_difference() const { return onset_.rayleigh / rayleigh(1.0); }

  // Lorenz's parameters at the present heating.
  physis::continuum::LorenzParameters lorenz() const {
    auto p = physis::continuum::lorenz_parameters(rayleigh(), prandtl(), o_.roll_wavenumber, o_.fluid.thermal_diffusivity, o_.depth);
    if (o_.scaling == LorenzScaling::cell_onset) p.r = rayleigh() / onset_.rayleigh;
    return p;
  }
  // The Hopf point of the convecting states, r_H = σ(σ + b + 3)/(σ − b − 1), and its ΔT.
  double hopf_r() const {
    const auto p = lorenz();
    return p.sigma * (p.sigma + p.b + 3) / (p.sigma - p.b - 1);
  }
  double hopf_temperature_difference() const {
    const double r = hopf_r();
    return o_.scaling == LorenzScaling::cell_onset ? r * critical_temperature_difference()
                                                   : r * critical_temperature_difference() * physis::continuum::marginal_rayleigh_free(o_.roll_wavenumber) / onset_.rayleigh;
  }

  const physis::dynamics::State& state() const { return xyz_; }
  void set_state(const physis::dynamics::State& x) { xyz_ = x; }
  void set_heating(double dT) { dT_ = dT; }

  runtime::TargetDescription describe() const override {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "convection cell";
    d.accepts = {tasks::Kind::flow, tasks::Kind::phenomenon, tasks::Kind::filter};
    d.inputs = {{"heating ΔT", "K"}};
    d.outputs = {{"circulation X", "1"}, {"temperature mode Y", "1"}, {"temperature mode Z", "1"}, {"Nusselt number", "1"}};
    d.time_scale = 1 / lorenz().time_scale;
    d.materials.materials = {{o_.fluid.name, "the layer"}, {"copper", "the plates"}, {"glass", "the side walls"}};
    d.materials.end_of_life = computer::EndOfLife::recyclable;
    d.model = "Rayleigh–Bénard onset by shooting (physis::continuum::critical_rayleigh); Lorenz's truncation (Saltzman 1962; Lorenz 1963)";
    d.limits = {"three modes, two-dimensional rolls", "not the turbulent flow at large r", "Lorenz's modes are those of free plates; r scaled to the cell's onset"};
    return d;
  }
  void reset(std::uint64_t seed) override {
    seed_ = seed, t_ = 0, dT_ = o_.temperature_difference;
    // A small disturbance of the conducting state, the same for every seed's run but signed by the seed.
    xyz_ = {seed % 2 ? -1e-3 : 1e-3, 0, 0};
  }
  // Holds the heating ΔT for `duration` s of the cell's time.
  void apply(std::span<const double> inputs, double duration) override {
    if (inputs.size() != 1) throw std::invalid_argument("physical_computing::fluids::ConvectionCell::apply: one input, ΔT");
    dT_ = inputs[0];
    const auto p = lorenz();
    const physis::dynamics::Lorenz field{p.sigma, p.r, p.b};
    physis::dynamics::integrate(field, 0.0, duration * p.time_scale, xyz_, o_.step);
    t_ += duration;
  }
  std::vector<double> read() const override {
    const auto p = lorenz();
    return {xyz_[0], xyz_[1], xyz_[2], physis::continuum::lorenz_nusselt(xyz_[2], p.r)};
  }
  double time() const override { return t_; }
  runtime::Provenance provenance() const override {
    auto p = simulator_provenance("convection cell", "Lorenz's truncation of Rayleigh–Bénard convection", seed_);
    p.conditions = {{"ΔT (K)", dT_}, {"depth (m)", o_.depth}, {"Ra", rayleigh()}};
    return p;
  }

 private:
  ConvectionCellOptions o_;
  physis::continuum::Onset onset_;
  physis::dynamics::State xyz_{0, 0, 0};
  double dT_ = 0, t_ = 0;
  std::uint64_t seed_ = 0;
};

}  // namespace physical_computing::fluids
