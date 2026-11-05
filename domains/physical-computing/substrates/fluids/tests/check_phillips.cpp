// Checks of Phillips's machine (phillips.hpp).
//
// - Its steady income is the multiplier's, I/(1 − c), and one time constant A R/(1 − c) brings it 1 − 1/e of
//   the way.
// - Its exact hold agrees with integrating A dh/dt = c h/R + I − h/R by Runge–Kutta (relative 1e−10).
// - The conjugacy: income h/R follows the income model dY/dt = (c Y + I − Y)/(A R), integrated apart, exactly
//   (relative 1e−10).
// - The runtime drives it as any backend.

#include <physical_computing/fluids/phillips.hpp>
#include <physical_computing/runtime/run.hpp>

#include <physis/dynamics/flows.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  fluids::PhillipsOptions o;
  fluids::PhillipsMachine m(o);
  const double I = 2e-4, tau = m.time_constant();
  m.reset(0);
  m.apply(std::vector<double>{I}, 40 * tau);
  s.check("Phillips: the steady income is the multiplier's, I/(1 − c) (relative)", m.income() / (I / (1 - o.propensity)) - 1, 1e-12);
  m.reset(0);
  m.apply(std::vector<double>{I}, tau);
  s.check("Phillips: one time constant brings it 1 − 1/e of the way (relative)", m.income() / (m.steady_income(I) * (1 - std::exp(-1.0))) - 1, 1e-12);

  physis::dynamics::State h{0};
  auto tank = [&](double, const physis::dynamics::State& x, physis::dynamics::State& dx) { dx[0] = (o.propensity * x[0] / o.resistance + I - x[0] / o.resistance) / o.area; };
  physis::dynamics::integrate(tank, 0.0, 0.7 * tau, h, tau / 2000);
  m.reset(0);
  m.apply(std::vector<double>{I}, 0.7 * tau);
  s.check("Phillips: the exact hold equals Runge–Kutta's integration of the tank (relative)", m.level() / h[0] - 1, 1e-10);
  physis::dynamics::State Y{0};
  auto income = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { dy[0] = (o.propensity * y[0] + I - y[0]) / (o.area * o.resistance); };
  physis::dynamics::integrate(income, 0.0, 0.7 * tau, Y, tau / 2000);
  s.check("Phillips: income h/R follows the income model, integrated apart (relative)", m.income() / Y[0] - 1, 1e-10);

  schedule::Schedule sch;
  sch.hold = tau / 4;
  const auto d = runtime::drive(m, {{I}, {I}, {2 * I}, {0}}, sch);
  s.require("Phillips: the runtime drives it, three readings a step", d.runs[0].readings.size() == 4 && d.runs[0].readings[0].size() == 3);
  return s.done();
}
