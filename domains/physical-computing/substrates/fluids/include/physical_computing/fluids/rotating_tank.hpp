// A rotating tank of water heated at its outer wall and cooled at its inner: the general circulation in a
// laboratory (Fultz et al. 1959; Hide and Mason 1975) (substrates/fluids/README.md, simulator 3).
//
// The tank. An annulus of gap L, mean radius R and depth D turns at Ω with a temperature contrast ΔT across
// the gap. Its thermal Rossby number Θ = g α ΔT D/(Ω² L²) and Taylor number Ta = 4 Ω² L⁵/(ν² D) set its
// regime (Hide and Mason 1975): axisymmetric flow at large Θ, regular baroclinic waves below a transition,
// vacillation and irregular flow at large Ta.
//
// Its model. Two layers of depth D/2 (physis::continuum::TwoLayerQG), the Coriolis parameter f = 2Ω, the
// reduced gravity g′ = g α ΔT (the stratification taken equal to the contrast across the gap: sloping
// convection), so F = f²/(g′ D/2); the layers' difference of velocity is the thermal wind between their
// mid-depths, U₁ − U₂ = g α ΔT D/(2 f L); the lower layer is damped by its Ekman layer at
// r = √(ν Ω)/(D/2) (Pedlosky 1987, §4.5; Greenspan 1968). The annulus is unrolled into a periodic channel of
// length 2πR between free-slip walls a gap L apart, by the method of images: physis's periodic box of width 2L
// with the state kept odd in y (TwoLayerQG's `channel`), whose modes are the channel's sin(jπy/L). The walls'
// boundary layers and the annulus's curvature are left out. An azimuthal wavenumber m is k = m/R.
//
// What follows. The deformation wavenumber is k_d L = 4/√Θ exactly, so without friction the azimuthal mode m
// of the lowest cross-gap structure is unstable when (m/R)² + π²/L² < k_d², that is
//
//   Θ < Θ_m = 16 / (π² + (m L / R)²),
//
// the model's upper symmetric transition for m (16/π² = 1.621 for a long channel; 1.559 for m = 1 at
// L/R = 0.625). The laboratory's transition is reported near Θ ≈ 1.58 (Hide and Mason 1975; to be checked
// against their figure); friction lowers the model's. Below it the most unstable m sets the number of waves.
//
// Limits (declared): quasi-geostrophic two-layer dynamics; no sidewall boundary layers, no curvature; the
// stratification tied to ΔT; a declared hyperviscous damping at the grid's scale. Checks: tests/check_rotating_tank.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

#include <calx/transforms/fft.hpp>
#include <physis/continuum/two_layer.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct RotatingTankOptions {
  double rotation = 0.5;              // Ω, rad/s
  double temperature_difference = 2;  // ΔT across the gap, K
  double depth = 0.1;                 // D, m
  double gap = 0.05;                  // L, m
  double mean_radius = 0.08;          // R, m
  Fluid fluid = Fluid::water();
  std::size_t nx = 256, ny = 32;
  double step = 0;                    // s; 0: a tenth of the advective limit at the shear
  double hyperviscous_rate = 1.0;     // 1/s at the truncation scale, declared
  bool ekman_drag = true;
};

class RotatingTank {
 public:
  explicit RotatingTank(const RotatingTankOptions& o) : o_(o) {
    if (!(o.rotation > 0) || !(o.temperature_difference > 0) || !(o.depth > 0) || !(o.gap > 0) || !(o.mean_radius > o.gap / 2))
      throw std::invalid_argument("physical_computing::fluids::RotatingTank: positive Ω, ΔT, depth and gap, and a radius beyond half the gap");
  }

  double coriolis() const { return 2 * o_.rotation; }
  double reduced_gravity() const { return standard_gravity * o_.fluid.thermal_expansion * o_.temperature_difference; }
  double thermal_rossby() const { return reduced_gravity() * o_.depth / (o_.rotation * o_.rotation * o_.gap * o_.gap); }
  double taylor() const {
    const double nu = o_.fluid.kinematic_viscosity, L = o_.gap;
    return 4 * o_.rotation * o_.rotation * L * L * L * L * L / (nu * nu * o_.depth);
  }
  double F() const { const double f = coriolis(); return f * f / (reduced_gravity() * o_.depth / 2); }
  double deformation_wavenumber() const { return std::sqrt(2 * F()); }
  double shear() const { return reduced_gravity() * o_.depth / (2 * coriolis() * o_.gap); }
  double drag() const { return o_.ekman_drag ? std::sqrt(o_.fluid.kinematic_viscosity * o_.rotation) / (o_.depth / 2) : 0; }
  // The model's upper symmetric transition without friction: for a long channel, and for azimuthal wavenumber m.
  static double inviscid_critical_thermal_rossby() { return 16 / (std::numbers::pi * std::numbers::pi); }
  double inviscid_critical_thermal_rossby(int m) const {
    const double a = m * o_.gap / o_.mean_radius;
    return 16 / (std::numbers::pi * std::numbers::pi + a * a);
  }

  physis::continuum::TwoLayerQG::Parameters qg_parameters() const {
    physis::continuum::TwoLayerQG::Parameters p;
    p.nx = o_.nx, p.ny = o_.ny;
    p.lx = 2 * std::numbers::pi * o_.mean_radius, p.ly = 2 * o_.gap;
    p.F1 = p.F2 = F();
    p.U1 = shear() / 2, p.U2 = -shear() / 2;
    p.drag = drag();
    const double dx = std::min(p.lx / p.nx, p.ly / p.ny);
    // A tenth of the advective limit at the shear: the equilibrated eddies reach several times the shear.
    p.step = o_.step > 0 ? o_.step : 0.1 * dx / std::max(shear(), 1e-9);
    p.hyperviscous_rate = o_.hyperviscous_rate;
    p.channel = true;
    return p;
  }
  physis::continuum::TwoLayerQG model() const { return physis::continuum::TwoLayerQG(qg_parameters()); }

  // The linear growth rate (1/s) of azimuthal wavenumber m in the channel's lowest cross-gap mode.
  double growth(double m) const { return model_for_growth().linear_growth(m / o_.mean_radius, std::numbers::pi / o_.gap).rate; }
  // The most unstable whole azimuthal wavenumber, and its rate (0 and 0 when none grows).
  std::pair<int, double> most_unstable(int up_to = 60) const {
    int best = 0;
    double rate = 0;
    const auto m = model_for_growth();
    for (int k = 1; k <= up_to; k++) {
      const double g = m.linear_growth(k / o_.mean_radius, std::numbers::pi / o_.gap).rate;
      if (g > rate) rate = g, best = k;
    }
    return {best, rate};
  }

  // The azimuthal wavenumber carrying most of a layer's streamfunction variance in a model's present state.
  static int dominant_wavenumber(const physis::continuum::TwoLayerQG& qg, int layer = 1) {
    const auto psi = qg.streamfunction(layer);
    const std::size_t nx = qg.nx(), ny = qg.ny();
    std::vector<double> power(nx / 2 + 1, 0.0);
    std::vector<std::complex<double>> row(nx);
    for (std::size_t j = 0; j < ny; j++) {
      for (std::size_t i = 0; i < nx; i++) row[i] = psi[j * nx + i];
      calx::transforms::fft(row);
      for (std::size_t i = 1; i <= nx / 2; i++) power[i] += std::norm(row[i]);
    }
    std::size_t best = 1;
    for (std::size_t i = 2; i <= nx / 2; i++)
      if (power[i] > power[best]) best = i;
    return static_cast<int>(best);
  }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "rotating tank";
    d.accepts = {tasks::Kind::phenomenon, tasks::Kind::flow};
    d.inputs = {{"rotation Ω", "rad/s"}, {"contrast ΔT", "K"}};
    d.outputs = {{"upper layer's streamfunction", "m²/s"}, {"lower layer's streamfunction", "m²/s"}};
    d.time_scale = 1 / coriolis();
    d.materials.materials = {{o_.fluid.name, "the working fluid"}, {"acrylic (PMMA)", "the annulus"}, {"copper", "the heated outer wall"}, {"aluminium", "the turntable"}};
    d.materials.end_of_life = computer::EndOfLife::recyclable;
    d.model = "two-layer quasi-geostrophic flow (physis::continuum::TwoLayerQG; Phillips 1954); the annulus as a periodic channel (Hide and Mason 1975)";
    d.limits = {"quasi-geostrophic, two layers", "no sidewall boundary layers or curvature", "stratification tied to ΔT", "declared hyperviscosity"};
    return d;
  }
  const RotatingTankOptions& options() const { return o_; }

 private:
  physis::continuum::TwoLayerQG model_for_growth() const {
    auto p = qg_parameters();
    p.nx = p.ny = 8;  // the linear rates need no grid
    return physis::continuum::TwoLayerQG(p);
  }
  RotatingTankOptions o_;
};

}  // namespace physical_computing::fluids
