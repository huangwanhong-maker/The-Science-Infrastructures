// Checks of reaction networks as programs (programs.hpp).
//
// - Functions on counts, exactly: from every input with counts up to 4, the program's chain is absorbed with the
//   output g(x) with probability exactly 1 (rationals), for x₁ + x₂, 2x, ⌊x/2⌋, min and max.
// - A flow lowered by dual rails is a computer of the flow: the commuting square (computer::commutes) holds for
//   Lorenz's system over ten starts to t = 0.5, encoding and readout linear, within 1e−6.

#include <physical_computing/chemical/programs.hpp>
#include <physical_computing/computer/square.hpp>

#include <physis/dynamics/flows.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

using namespace physical_computing;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

int main() {
  calx::verification::Suite s;
  auto exact = [&](const chemical::Program& p, std::function<long long(const std::vector<long long>&)> g) {
    bool ok = true;
    std::vector<long long> x(p.inputs.size(), 0);
    while (true) {
      const auto law = chemical::run_exactly(p, x);
      ok = ok && law.size() == 1 && law.begin()->first == g(x) && law.begin()->second == Rational(1);
      std::size_t j = 0;
      while (j < x.size() && x[j] == 4) x[j] = 0, j++;
      if (j == x.size()) break;
      x[j]++;
    }
    s.require("program " + p.name + ": the output is g(x) with probability exactly 1, inputs up to 4", ok);
  };
  exact(chemical::sum_program(), [](const auto& x) { return x[0] + x[1]; });
  exact(chemical::double_program(), [](const auto& x) { return 2 * x[0]; });
  exact(chemical::half_program(), [](const auto& x) { return x[0] / 2; });
  exact(chemical::min_program(), [](const auto& x) { return std::min(x[0], x[1]); });
  exact(chemical::max_program(), [](const auto& x) { return std::max(x[0], x[1]); });

  {
    const auto P = physis::chemistry::lorenz_polynomials(10, 28, Rational(Integer(8), Integer(3)));
    const auto D = chemical::lower_flow(P);
    std::vector<physis::dynamics::State> K;
    for (int k = 0; k < 10; k++) K.push_back({-8.0 + 1.7 * k, 5.0 - 1.1 * k, 15.0 + 1.3 * k});
    auto encode = [&](const physis::dynamics::State& x) {
      physis::dynamics::State c(6, 0.0);
      for (int i = 0; i < 3; i++) c[D.plus[i]] = std::max(x[i], 0.0), c[D.minus[i]] = std::max(-x[i], 0.0);
      return c;
    };
    auto evolve = [&](physis::dynamics::State c) {
      auto f = [&](double, const physis::dynamics::State& y, physis::dynamics::State& dy) { physis::chemistry::mass_action(D.network, y, dy); };
      physis::dynamics::integrate(f, 0.0, 0.5, c, 1e-4);
      return c;
    };
    auto read = [&](const physis::dynamics::State& c) { return physis::dynamics::State{c[0] - c[1], c[2] - c[3], c[4] - c[5]}; };
    auto model = [&](physis::dynamics::State x) {
      physis::dynamics::integrate(physis::dynamics::Lorenz{}, 0.0, 0.5, x, 1e-4);
      return x;
    };
    auto metric = [](const physis::dynamics::State& a, const physis::dynamics::State& b) {
      double d = 0;
      for (std::size_t i = 0; i < a.size(); i++) d = std::max(d, std::fabs(a[i] - b[i]));
      return d;
    };
    const auto v = computer::commutes(K, encode, evolve, read, model, metric, 1e-6, computer::Map::linear, computer::Map::linear);
    s.require("dual rails: the square commutes for Lorenz's flow over ten starts", v.holds && v.cases == 10, "largest error " + std::to_string(v.error));
  }
  return s.done();
}
