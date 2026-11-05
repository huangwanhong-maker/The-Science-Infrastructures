// Checks of the convection cell (convection_cell.hpp).
//
// - 1 cm of water between rigid plates convects at ΔT_c = 1707.76 ν κ/(g α d³) = 0.1208 K; its Prandtl number
//   is ν/κ = 7.02.
// - Below the onset the rolls die away; at r = 5 the cell settles on its convecting state, X² = b(r − 1), and
//   carries 1 + 2(r − 1)/r = 2.6 times the heat of conduction.
// - The Hopf point of the convecting states, r_H = σ(σ + b + 3)/(σ − b − 1): the Jacobian at C₊ is stable just
//   below it and has a complex pair with positive real part just above (calx's eigenvalues).
// - Beyond it the cell is chaotic: its largest Lyapunov exponent is positive (physis).
// - The runtime drives it with the heating as input, and the provenance records ΔT.

#include <physical_computing/fluids/convection_cell.hpp>
#include <physical_computing/runtime/run.hpp>

#include <calx/linalg/eigenvalues.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  fluids::ConvectionCellOptions o;
  fluids::ConvectionCell cell(o);
  const auto w = fluids::Fluid::water();
  const double expected = 1707.762 * w.kinematic_viscosity * w.thermal_diffusivity / (9.80665 * w.thermal_expansion * 1e-6);
  s.check("cell: ΔT_c of 1 cm of water between rigid plates (relative)", cell.critical_temperature_difference() / expected - 1, 1e-5);
  s.check("cell: and its value (K)", cell.critical_temperature_difference() - 0.1208, 5e-4, "K");
  s.check("cell: the Prandtl number ν/κ", cell.prandtl() - 1.004e-6 / 1.43e-7, 1e-12);

  const double dTc = cell.critical_temperature_difference();
  {
    fluids::ConvectionCell c(o);
    c.reset(1);
    c.apply(std::vector<double>{0.5 * dTc}, 3000);
    s.check("cell: below the onset the rolls die away (|X|)", c.read()[0], 1e-6);
  }
  {
    fluids::ConvectionCell c(o);
    c.reset(2);
    c.apply(std::vector<double>{5 * dTc}, 3000);
    const auto r = c.read();
    const auto p = c.lorenz();
    s.check("cell: at r = 5, the convecting state X² = b(r − 1) (relative)", r[0] * r[0] / (p.b * (p.r - 1)) - 1, 1e-8);
    s.check("cell: carries 1 + 2(r − 1)/r of conduction's heat", r[3] - 2.6, 1e-8);
  }
  {
    const double rH = cell.hopf_r();
    auto leading = [&](double r) {
      fluids::ConvectionCell c(o);
      c.set_heating(r * dTc);
      const auto p = c.lorenz();
      const double q = std::sqrt(p.b * (p.r - 1));
      const physis::dynamics::Lorenz f{p.sigma, p.r, p.b};
      const auto J = physis::dynamics::jacobian(f, 0.0, physis::dynamics::State{q, q, p.r - 1});
      double best = -1e300, imag = 0;
      for (const auto& e : calx::linalg::eigenvalues(J))
        if (e.real() > best) best = e.real(), imag = e.imag();
      return std::pair{best, imag};
    };
    const auto below = leading(0.98 * rH), above = leading(1.02 * rH);
    s.require("cell: the convecting state is stable just below r_H", below.first < 0);
    s.require("cell: and loses it to a complex pair just above (a Hopf bifurcation)", above.first > 0 && std::fabs(above.second) > 1);
    fluids::ConvectionCell c(o);
    c.set_heating(1.3 * rH * dTc);
    const auto p = c.lorenz();
    physis::dynamics::LyapunovOptions lo;
    lo.step = 0.005, lo.transient = 20, lo.time = 300, lo.count = 1;
    const auto L = physis::dynamics::lyapunov_spectrum(physis::dynamics::Lorenz{p.sigma, p.r, p.b}, physis::dynamics::State{1, 1, p.r}, lo);
    s.require("cell: beyond r_H it is chaotic, λ₁ > 0", L.exponents[0] > 0.3, "λ₁ = " + std::to_string(L.exponents[0]));
  }
  {
    fluids::ConvectionCell c(o);
    schedule::Schedule sch;
    sch.hold = 10;
    const auto d = runtime::drive(c, {{2 * dTc}, {3 * dTc}, {4 * dTc}}, sch);
    s.require("cell: the runtime drives it, four readings a step", d.runs[0].readings.size() == 3 && d.runs[0].readings[0].size() == 4);
    s.require("cell: the provenance records the heating", !d.runs[0].provenance.conditions.empty() && std::fabs(d.runs[0].provenance.conditions[0].second - 4 * dTc) < 1e-15);
  }
  return s.done();
}
