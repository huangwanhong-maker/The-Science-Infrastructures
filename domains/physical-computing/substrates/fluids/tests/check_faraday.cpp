// Checks of the Faraday layer (faraday.hpp), water 1 cm deep shaken at 40 Hz.
//
// - The resonant wavenumber has ω₀ = Ω/2.
// - At it, the Floquet threshold equals the small-damping 4γω₀/(k tanh kh) (within 2%: γ/ω₀ ≈ 0.004).
// - The onset lies at the resonant wavenumber (within 1%), its response subharmonic.
// - Doubling the damping doubles the threshold (within 3%).

#include <physical_computing/fluids/faraday.hpp>

#include <calx/verification/check.hpp>

#include <cmath>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  fluids::FaradayOptions o;
  const fluids::FaradayLayer layer(o);
  const double kr = layer.resonant_wavenumber();
  s.check("Faraday: the resonant mode has ω₀ = Ω/2 (relative)", layer.natural_frequency(kr) / (o.drive_frequency / 2) - 1, 1e-10);
  const double A = layer.threshold(kr);
  s.check("Faraday: the threshold at resonance is 4γω₀/(k tanh kh) (relative)", A / layer.small_damping_threshold(kr) - 1, 0.02);
  const auto on = layer.onset();
  s.check("Faraday: the onset is at the resonant wavenumber (relative)", on.wavenumber / kr - 1, 0.01);
  s.require("Faraday: its response is subharmonic", layer.subharmonic(on.wavenumber, 1.1 * on.acceleration));
  fluids::FaradayOptions d = o;
  d.damping_factor = 2;
  s.check("Faraday: doubling the damping doubles the threshold (relative)", fluids::FaradayLayer(d).threshold(kr) / A / 2 - 1, 0.03);
  return s.done();
}
