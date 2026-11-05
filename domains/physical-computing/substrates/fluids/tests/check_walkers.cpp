// Checks of the walking droplet (walkers.hpp; Oza, Rosales and Bush 2013).
//
// - The steady speed u = √(2β − 1 − √(1 + 4β))/√2 solves u = β(1 − 1/√(1 + u²))/u, and vanishes at β = 2.
// - Simulated at β = 4, κ = 1, the droplet settles at that speed (within 1.5%); at β = 1.5 it comes to rest.
// - Started along a line, it stays on it.

#include <physical_computing/fluids/walkers.hpp>

#include <calx/verification/check.hpp>

#include <cmath>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  const double u = fluids::Walker::steady_speed(4);
  s.check("walker: the steady speed solves u = β(1 − 1/√(1 + u²))/u", u - 4 * (1 - 1 / std::sqrt(1 + u * u)) / u, 1e-14);
  s.require("walker: no steady walking at β ≤ 2", fluids::Walker::steady_speed(2) == 0 && fluids::Walker::steady_speed(2.001) > 0);
  {
    fluids::WalkerOptions o;
    o.beta = 4, o.kappa = 1;
    fluids::Walker w(o);
    w.run(60);
    s.check("walker: at β = 4 it settles at the steady speed (relative)", w.speed() / u - 1, 0.015);
    s.check("walker: started along a line, it stays on it", w.position()[1], 1e-12);
  }
  {
    fluids::WalkerOptions o;
    o.beta = 1.5;
    fluids::Walker w(o);
    w.run(60);
    s.check("walker: at β = 1.5 it comes to rest (speed)", w.speed(), 1e-3);
  }
  return s.done();
}
