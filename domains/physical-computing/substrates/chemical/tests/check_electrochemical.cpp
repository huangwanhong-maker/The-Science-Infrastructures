// Checks of the electrochemical oscillators (electrochemical.hpp).
//
// - The transition of an array of 500 electrodes of Lorentzian spread γ is at K_c = 2γ; at 2K_c its order
//   parameter is √(1 − 2γ/K) (within 0.04); at K_c/2 it stays incoherent (r < 0.12).
// - Its materials name the platinum group metals of its counter electrode as critical.

#include <physical_computing/chemical/electrochemical.hpp>

#include <calx/verification/check.hpp>

#include <cmath>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  chemical::ElectrodeArray a(500, 2 * 3.141592653589793 * 0.4, 0.05);
  s.check("electrodes: K_c = 2γ", a.critical_coupling() - 0.1, 1e-15);
  const double K = 2 * a.critical_coupling();
  s.check("electrodes: r = √(1 − 2γ/K) at K = 2K_c", a.run(K, 400, 100) - a.predicted_order(K), 0.04);
  chemical::ElectrodeArray b(500, 2 * 3.141592653589793 * 0.4, 0.05);
  s.require("electrodes: incoherent at K_c/2", b.run(0.5 * b.critical_coupling(), 400, 100) < 0.12);
  const auto crit = a.describe().materials.critical();
  s.require("electrodes: the counter electrode's platinum group metals are critical", crit.size() == 1 && crit[0] == "platinum group metals");
  return s.done();
}
