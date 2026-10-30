// Checks of a joint's geared motor (actuation/motor.hpp): its torque–speed line within the driver's peak,
// none at the no-load speed, full braking against its motion, and its electrical power.

#include <robotics/actuation/motor.hpp>

#include <calx/verification/check.hpp>

#include <cmath>

int main() {
  calx::verification::Suite s;
  const robotics::actuation::GearedMotor m{5.0, 30, 8e-4, 2e-3, 0.02, 0.45, 0.35};
  s.check("motor: at rest it gives its peak, not its stall torque (N m)", m.available(100, 0) - m.peak, 1e-12);
  s.check("motor: and its peak backward (N m)", m.available(-100, 0) + m.peak, 1e-12);
  s.check("motor: at its no-load speed it cannot drive (N m)", m.available(100, m.no_load), 1e-12);
  s.check("motor: at 80 % of it, a fifth of its stall torque (N m)", m.available(100, 0.8 * m.no_load) - 0.2 * m.stall(), 1e-12);
  s.check("motor: against its motion it brakes at its full peak (N m)", m.available(-100, 0.8 * m.no_load) + m.peak, 1e-12);
  s.check("motor: its electrical power, τω plus the copper loss (τ/k_τ)² R (W)", m.electrical_power(2.0, 10.0) - (20.0 + (2.0 / 0.45) * (2.0 / 0.45) * 0.35), 1e-12);
  return s.done();
}
