// A droplet walking on a vibrated bath, guided by the waves of its own past bounces: hydrodynamic pilot-wave
// dynamics (Couder and Fort 2006; Bush 2015) (substrates/fluids/README.md, simulator 5).
//
// The stroboscopic model (Oza, Rosales and Bush 2013, J. Fluid Mech. 737, 552), dimensionless:
//
//   κ ẍ + ẋ = β ∫_{−∞}^{t} J₁(|x(t) − x(s)|) (x(t) − x(s))/|x(t) − x(s)| e^{−(t − s)} ds + F_ext(x, ẋ),
//
// lengths in units of 1/k_F (k_F the Faraday wavenumber), time in units of the waves' memory; κ the droplet's
// inertia relative to its drag, β the strength of the wave force (both grow with the drive's memory).
// Steady walking x = u t needs u = β ∫₀^∞ J₁(u z) e^{−z} dz = β (1 − 1/√(1 + u²))/u (the Laplace transform of
// J₁), with the non-zero root
//
//   u = (1/√2) √(2β − 1 − √(1 + 4β))   for β > 2,
//
// a pitchfork from rest at β = 2 whatever κ. The integral is the trapezoidal rule over the stored path, cut
// where e^{−(t − s)} < e^{−memory}; the step is semi-implicit in the drag. An external force (a harmonic well,
// Coriolis in a rotating frame) is a callable, so that the model can carry the experiments of Perrard et al.
// (2014) and Harris and Bush (2014).
//
// Limits (declared): the stroboscopic model averages the bouncing; no vertical dynamics; the "quantum-like"
// statistics it shows are claims to be held to reproducible checks (Andersen et al. 2015 did not reproduce
// single-particle interference). Checks: tests/check_walkers.cpp.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <stdexcept>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

using Vec2 = std::array<double, 2>;

struct WalkerOptions {
  double kappa = 1.0;   // κ
  double beta = 4.0;    // β
  double step = 0.02;   // dt
  double memory = 25;   // the path kept, in units of the memory time
  std::function<Vec2(const Vec2& x, const Vec2& v)> external;  // F_ext, or none
};

class Walker {
 public:
  explicit Walker(const WalkerOptions& o, Vec2 x0 = {0, 0}, Vec2 v0 = {0.1, 0}) : o_(o), x_(x0), v_(v0) {
    if (!(o.kappa > 0) || !(o.step > 0) || !(o.memory > 0)) throw std::invalid_argument("physical_computing::fluids::Walker: positive κ, step and memory");
    path_.push_back(x0);
  }

  // u for β > 2, else 0.
  static double steady_speed(double beta) { return beta > 2 ? std::sqrt(2 * beta - 1 - std::sqrt(1 + 4 * beta)) / std::sqrt(2.0) : 0.0; }

  // The wave force at the present position from the stored path.
  Vec2 wave_force() const {
    Vec2 f{0, 0};
    const std::size_t n = path_.size();
    if (n < 2) return f;
    for (std::size_t q = 0; q + 1 < n; q++) {  // q = 0 the oldest; the present point (q = n − 1) adds nothing
      const double age = static_cast<double>(n - 1 - q) * o_.step;
      const double dx = x_[0] - path_[q][0], dy = x_[1] - path_[q][1], r = std::hypot(dx, dy);
      if (r < 1e-12) continue;
      double w = std::cyl_bessel_j(1.0, r) / r * std::exp(-age);
      if (q == 0) w *= 0.5;  // the trapezoid's end
      f[0] += w * dx, f[1] += w * dy;
    }
    return {o_.beta * o_.step * f[0], o_.beta * o_.step * f[1]};
  }

  void step() {
    Vec2 F = wave_force();
    if (o_.external) {
      const Vec2 e = o_.external(x_, v_);
      F[0] += e[0], F[1] += e[1];
    }
    // κ v̇ = F − v, implicit in the drag.
    const double h = o_.step;
    for (int c = 0; c < 2; c++) v_[c] = (v_[c] + h * F[c] / o_.kappa) / (1 + h / o_.kappa);
    for (int c = 0; c < 2; c++) x_[c] += h * v_[c];
    path_.push_back(x_);
    const std::size_t keep = static_cast<std::size_t>(o_.memory / o_.step) + 1;
    while (path_.size() > keep) path_.pop_front();
    t_ += h;
  }
  void run(double duration) {
    const auto n = static_cast<std::size_t>(std::llround(duration / o_.step));
    for (std::size_t k = 0; k < n; k++) step();
  }

  const Vec2& position() const { return x_; }
  const Vec2& velocity() const { return v_; }
  double speed() const { return std::hypot(v_[0], v_[1]); }
  double time() const { return t_; }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "walking droplet";
    d.accepts = {tasks::Kind::phenomenon};
    d.inputs = {{"memory (β)", "1"}, {"external force", "1"}};
    d.outputs = {{"position", "1/k_F"}, {"velocity", "1/(k_F T_M)"}};
    d.materials.materials = {{"silicone (PDMS)", "the bath and the droplet"}, {"copper", "the shaker's coil"}, {"neodymium magnet", "the shaker"}};
    d.materials.end_of_life = computer::EndOfLife::neither;
    d.model = "the stroboscopic pilot-wave model (Oza, Rosales and Bush 2013)";
    d.limits = {"bouncing averaged out", "no vertical dynamics", "quantum-like statistics are claims, not shown (Andersen et al. 2015)"};
    return d;
  }

 private:
  WalkerOptions o_;
  Vec2 x_, v_;
  std::deque<Vec2> path_;
  double t_ = 0;
};

}  // namespace physical_computing::fluids
