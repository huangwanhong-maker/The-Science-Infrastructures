// Checks of strand displacement (strand_displacement.hpp).
//
// - A formal network X → Y, X + Y → Z compiles to fuels and cascades; as C_max and q_max grow tenfold together (fuel
//   in excess, cascades fast), the formal species' trajectories approach the formal network's (their error falls
//   by more than five), integrated implicitly since the cascades are stiff.
// - The unimolecular cascade conserves its fuel exactly: C_max − G = Y's share released plus the gates in flight
//   (O), at every time.
// - A leak releases the output with no input at all.

#include <physical_computing/molecular/strand_displacement.hpp>

#include <physis/dynamics/flows.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

#include <physis/chemistry/field_kinetics.hpp>

using namespace physical_computing;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;
  physis::chemistry::Network formal;
  formal.species = {"X", "Y", "Z"};
  formal.add("X → Y", {{"X", 1}}, {{"Y", 1}}, 1);
  formal.add("X + Y → Z", {{"X", 1}, {"Y", 1}}, {{"Z", 1}}, 1);
  auto formal_run = [&](double t) {
    physis::dynamics::State c{1, 0, 0};
    auto f = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { physis::chemistry::mass_action(formal, y, dy); };
    physis::dynamics::integrate(f, 0.0, t, c, 1e-4);
    return c;
  };
  // The limit is C_max → ∞ and q_max → ∞ together (fuel in excess, cascades fast); the cascades are then stiff
  // (q_max C_max up to 10⁷ against formal rates of 1), so they are integrated implicitly (physis's integrate_stiff).
  auto dsd_error = [&](long long cmax, long long qmax) {
    const auto D = molecular::compile_strand_displacement(formal, Rational(cmax), Rational(qmax));
    std::vector<double> c(D.initial.begin(), D.initial.end());
    c[0] = 1;
    c = physis::chemistry::integrate_stiff(c.size(), [&](const double* y, double* dy) {
      std::vector<double> yy(y, y + D.network.species.size()), dd(yy.size());
      physis::chemistry::mass_action(D.network, yy, dd);
      for (std::size_t i = 0; i < dd.size(); i++) dy[i] = dd[i];
    }, c, 2.0, 1e-3);
    const auto ref = formal_run(2.0);
    double e = 0;
    for (int i = 0; i < 3; i++) e = std::max(e, std::fabs(c[i] - ref[i]));
    return e;
  };
  const double e100 = dsd_error(100, 1000), e1000 = dsd_error(1000, 10000);
  s.require("compiled: the formal species approach the formal network's as C_max and q_max grow tenfold (error falls fivefold)", e100 / e1000 > 5,
            "errors " + std::to_string(e100) + " and " + std::to_string(e1000));
  s.check("compiled: at C_max = 1000, q_max = 10⁴ the trajectories agree to 1e−3", e1000, 1e-3);

  {
    physis::chemistry::Network uni;
    uni.species = {"X", "Y"};
    uni.add("X → Y", {{"X", 1}}, {{"Y", 1}}, 1);
    const auto D = molecular::compile_strand_displacement(uni, 100, 100);
    physis::dynamics::State c(D.initial.begin(), D.initial.end());
    c[0] = 1;
    const std::size_t G = D.network.find("G0"), O = D.network.find("O0"), Y = D.network.find("Y");
    double worst = 0;
    auto f = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { physis::chemistry::mass_action(D.network, y, dy); };
    physis::dynamics::integrate(f, 0.0, 3.0, c, 1e-5, [&](double, const physis::dynamics::State& y) { worst = std::max(worst, std::fabs((100 - y[G]) - (y[Y] + y[O]))); });
    s.check("fuel: C_max − G = Y + O at every time", worst, 1e-10);
  }
  {
    physis::chemistry::Network uni;
    uni.species = {"X", "Y"};
    uni.add("X → Y", {{"X", 1}}, {{"Y", 1}}, 1);
    const auto D = molecular::compile_strand_displacement(uni, 100, 100, Rational(Integer(1), Integer(1000000)));
    physis::dynamics::State c(D.initial.begin(), D.initial.end());
    auto f = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { physis::chemistry::mass_action(D.network, y, dy); };
    physis::dynamics::integrate(f, 0.0, 1.0, c, 1e-4);
    s.require("leak: the output appears with no input", c[D.network.find("Y")] > 1e-3);
  }
  return s.done();
}
