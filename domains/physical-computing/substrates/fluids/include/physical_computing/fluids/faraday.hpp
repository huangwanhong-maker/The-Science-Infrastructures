// A layer of liquid on a shaker: the Faraday instability (Faraday 1831; Benjamin and Ursell 1954)
// (substrates/fluids/README.md, simulator 4).
//
// A container shaken vertically with acceleration A cos Ωt feels the effective gravity g − A cos Ωt. Each
// surface mode of wavenumber k on a layer of depth h then obeys a damped Mathieu equation (Benjamin and Ursell
// 1954; viscous damping after Lamb):
//
//   ä + 2γ ȧ + k tanh(kh) (g + σk²/ρ − A cos Ωt) a = 0,   γ = 2νk²,
//
// whose Floquet multipliers (physis::dynamics::parametric_oscillator) give its stability. Below threshold the
// surface is flat; the first modes to grow are those near half the drive's frequency, ω₀(k) = Ω/2, and they
// respond at Ω/2: the subharmonic waves Faraday saw. For small damping the threshold at exact resonance is
// A_c = 4γω₀ / (k tanh kh) (the damped tongue's tip ε = 4γ/ω₀, with ε = A k tanh(kh)/ω₀²).
//
// Limits (declared): linear modes in an unbounded layer (no container's modes, no meniscus); bulk viscous
// damping only (Kumar and Tuckerman's 1994 viscous Floquet theory gives larger thresholds; a declared factor
// stands for the boundary layers); no pattern selection, which is nonlinear. Checks: tests/check_faraday.cpp.
#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

#include <calx/optimisation/minimise.hpp>
#include <calx/solvers/roots.hpp>
#include <physis/continuum/waves.hpp>
#include <physis/dynamics/floquet.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct FaradayOptions {
  Fluid fluid = Fluid::water();
  double depth = 0.01;                                  // h, m
  double drive_frequency = 2 * std::numbers::pi * 40;   // Ω, rad/s
  double damping_factor = 1;                            // multiplies 2νk²: walls and films, declared
  std::size_t steps_per_period = 1500;                  // of the Floquet integration
};

class FaradayLayer {
 public:
  explicit FaradayLayer(const FaradayOptions& o) : o_(o) {
    if (!(o.depth > 0) || !(o.drive_frequency > 0) || !(o.damping_factor > 0)) throw std::invalid_argument("physical_computing::fluids::FaradayLayer: positive depth, drive and damping");
  }

  double natural_frequency(double k) const {
    return physis::continuum::gravity_capillary_frequency(k, o_.depth, standard_gravity, o_.fluid.surface_tension, o_.fluid.density);
  }
  double damping(double k) const { return o_.damping_factor * physis::continuum::viscous_damping_rate(k, o_.fluid.kinematic_viscosity); }
  // The largest |multiplier| of mode k at the drive's acceleration A (m/s²).
  double multiplier(double k, double A) const {
    const double w0 = natural_frequency(k), th = std::tanh(k * o_.depth);
    return physis::dynamics::parametric_oscillator(damping(k), w0 * w0, A * k * th, o_.drive_frequency, o_.steps_per_period).spectral_radius;
  }
  // The acceleration at which mode k becomes unstable, by bisection on the multiplier.
  double threshold(double k, double upper = 0) const {
    double lo = 0, hi = upper > 0 ? upper : 4 * standard_gravity;
    if (multiplier(k, hi) <= 1) return INFINITY;
    for (int it = 0; it < 60; it++) {
      const double mid = 0.5 * (lo + hi);
      (multiplier(k, mid) > 1 ? hi : lo) = mid;
    }
    return 0.5 * (lo + hi);
  }
  // The wavenumber resonant with the subharmonic, ω₀(k) = Ω/2.
  double resonant_wavenumber() const {
    const double half = o_.drive_frequency / 2;
    auto f = [&](double k) { return natural_frequency(k) - half; };
    double hi = 1;
    while (f(hi) < 0) hi *= 2;
    return calx::solvers::brent(f, hi / 2 > 1e-6 ? hi / 2 : 1e-6, hi, 1e-12).x;
  }
  // The small-damping threshold at exact resonance, 4γω₀/(k tanh kh).
  double small_damping_threshold(double k) const { return 4 * damping(k) * natural_frequency(k) / (k * std::tanh(k * o_.depth)); }

  struct Onset {
    double acceleration = 0, wavenumber = 0, wavelength = 0;
  };
  // The least threshold over k near resonance.
  Onset onset() const {
    const double kr = resonant_wavenumber();
    const auto m = calx::optimisation::brent_minimise([&](double k) { return threshold(k); }, 0.8 * kr, 1.25 * kr, 1e-7);
    return {m.f, m.x, 2 * std::numbers::pi / m.x};
  }
  // Whether the growing response of mode k at A is subharmonic (a negative real multiplier).
  bool subharmonic(double k, double A) const {
    const double w0 = natural_frequency(k), th = std::tanh(k * o_.depth);
    const auto F = physis::dynamics::parametric_oscillator(damping(k), w0 * w0, A * k * th, o_.drive_frequency, o_.steps_per_period);
    const auto& mu = std::abs(F.multipliers[0]) >= std::abs(F.multipliers[1]) ? F.multipliers[0] : F.multipliers[1];
    return mu.real() < 0 && std::fabs(mu.imag()) < 1e-9 * std::abs(mu);
  }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "Faraday layer";
    d.accepts = {tasks::Kind::phenomenon};
    d.inputs = {{"drive's acceleration A", "m/s²"}, {"drive's frequency Ω", "rad/s"}};
    d.outputs = {{"surface modes' amplitudes", "m"}};
    d.time_scale = 2 * std::numbers::pi / o_.drive_frequency;
    d.materials.materials = {{o_.fluid.name, "the layer"}, {"aluminium", "the container"}, {"copper", "the shaker's coil"}, {"neodymium magnet", "the shaker"}};
    d.materials.hazards = {"neodymium magnets need light rare earth elements (critical)"};
    d.materials.end_of_life = computer::EndOfLife::recyclable;
    d.model = "a damped Mathieu equation per mode (Benjamin and Ursell 1954), its Floquet multipliers (physis::dynamics)";
    d.limits = {"linear modes in an unbounded layer", "bulk viscous damping, times a declared factor", "no pattern selection"};
    return d;
  }
  const FaradayOptions& options() const { return o_; }

 private:
  FaradayOptions o_;
};

}  // namespace physical_computing::fluids
