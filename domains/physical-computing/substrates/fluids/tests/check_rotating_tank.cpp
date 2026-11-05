// Checks of the rotating tank (rotating_tank.hpp).
//
// - The thermal Rossby number g α ΔT D/(Ω² L²), and the deformation wavenumber k_d L = 4/√Θ exactly.
// - Without friction, the azimuthal mode m = 1 is unstable just below Θ₁ = 16/(π² + (L/R)²) and every mode is
//   neutral just above it: the model's upper symmetric transition (16/π² for a long channel).
// - Ekman friction lowers the growth.
// - From small noise in the wave regime, a run's dominant azimuthal wavenumber after the linear phase is the
//   most unstable one (within 1), and the flow stays finite as it equilibrates.

#include <physical_computing/fluids/rotating_tank.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  fluids::RotatingTankOptions o;
  const fluids::RotatingTank tank(o);
  const auto w = fluids::Fluid::water();
  s.check("tank: Θ = g α ΔT D/(Ω² L²) (relative)", tank.thermal_rossby() / (9.80665 * w.thermal_expansion * 2 * 0.1 / (0.25 * 0.0025)) - 1, 1e-14);
  s.check("tank: k_d L = 4/√Θ (relative)", tank.deformation_wavenumber() * o.gap / (4 / std::sqrt(tank.thermal_rossby())) - 1, 1e-14);

  auto at_theta = [&](double theta, bool drag) {
    fluids::RotatingTankOptions q = o;
    q.ekman_drag = drag;
    q.temperature_difference = theta * q.rotation * q.rotation * q.gap * q.gap / (9.80665 * w.thermal_expansion * q.depth);
    return fluids::RotatingTank(q);
  };
  const double T1 = tank.inviscid_critical_thermal_rossby(1), Tc = fluids::RotatingTank::inviscid_critical_thermal_rossby();
  s.check("tank: the long channel's transition 16/π²", Tc - 16 / (std::numbers::pi * std::numbers::pi), 1e-15);
  s.require("tank: without friction, m = 1 is unstable just below Θ₁ = 16/(π² + (L/R)²)", at_theta(0.98 * T1, false).growth(1) > 0);
  s.require("tank: and every mode is neutral just above it", at_theta(1.02 * T1, false).most_unstable().second < 1e-12);
  s.require("tank: Ekman friction lowers the growth", at_theta(0.5 * Tc, true).most_unstable().second < at_theta(0.5 * Tc, false).most_unstable().second);

  {
    fluids::RotatingTankOptions q = o;
    q.nx = 128, q.ny = 16;
    const fluids::RotatingTank t(q);
    const auto [m, rate] = t.most_unstable();
    auto model = t.model();
    model.perturb(1e-9, 21);
    const double efold = 1 / rate;
    while (model.time() < 12 * efold) model.step();
    const int dominant = fluids::RotatingTank::dominant_wavenumber(model, 1);
    s.require("tank: the run's dominant wavenumber is the most unstable one (within 1)", std::abs(dominant - m) <= 1,
              "dominant " + std::to_string(dominant) + ", most unstable " + std::to_string(m));
    while (model.time() < 30 * efold) model.step();
    s.require("tank: the flow stays finite as it equilibrates", std::isfinite(model.energy()) && model.energy() > 0);
  }
  return s.done();
}
