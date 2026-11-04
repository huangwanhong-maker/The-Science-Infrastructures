// Checks of gravity (gravity.hpp): a quarter of a great circle is R₁π/2; on flows that equal their
// model exactly, every method returns the coefficients they were made with; on flows whose mean is
// the model's but half of which are zero (each cell 0 or twice its mean, at random), PPML with
// exporter and importer effects recovers the distance elasticity within four standard errors,
// keeping the zeros, which least squares on logarithms must leave out.

#include <economics/trade/gravity.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;
  constexpr double pi = std::numbers::pi;
  s.check("a quarter of the equator is R₁π/2 = 10 007.5 km (relative)", trade::great_circle_km(0, 0, 0, pi / 2) / (trade::earth_mean_radius_km * pi / 2) - 1, 1e-15);

  calx::random::Stream z(20261104, 0x67726176u, 0, 0);  // "grav"
  const std::size_t n = 40;
  std::vector<double> lat(n), lon(n), alpha(n), gamma(n), size(n);
  for (std::size_t i = 0; i < n; i++) {
    lat[i] = z.uniform(-1.1, 1.1);
    lon[i] = z.uniform(-pi, pi);
    alpha[i] = z.normal();
    gamma[i] = z.normal();
    size[i] = std::exp(2 + z.normal());
  }
  auto mean_effects = [&](std::size_t i, std::size_t j) {
    return std::exp(-1.0 * std::log(trade::great_circle_km(lat[i], lon[i], lat[j], lon[j])) + alpha[i] + gamma[j] + 8);
  };
  auto mean_sizes = [&](std::size_t i, std::size_t j) {
    return std::exp(1.0 + 0.8 * std::log(size[i]) + 0.9 * std::log(size[j]) - 1.1 * std::log(trade::great_circle_km(lat[i], lon[i], lat[j], lon[j])));
  };
  trade::Rectangular exact(n, n), sized(n, n), noisy(n, n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      if (i == j) continue;
      exact(i, j) = mean_effects(i, j);
      sized(i, j) = mean_sizes(i, j);
      noisy(i, j) = z.uniform() < 0.5 ? 0.0 : 2 * mean_effects(i, j);
    }

  const auto fe = trade::gravity_effects(exact, lat, lon);
  s.check("PPML with effects, exact flows: the distance elasticity is −1", fe.distance_elasticity + 1, 1e-8);
  s.require("... over every ordered pair, none of them zero", fe.pairs == n * (n - 1) && fe.zeros == 0);
  double worst = 0;
  for (std::size_t k = 0; k < fe.pairs; k++) worst = std::max(worst, std::fabs(fe.fitted[k] / fe.flow[k] - 1));
  s.check("... and fits every flow (relative)", worst, 1e-7);

  const auto ps = trade::gravity_sizes(sized, lat, lon, size, trade::GravityMethod::ppml_sizes);
  s.check("PPML with sizes, exact flows: (1, 0.8, 0.9, −1.1)",
          std::fabs(ps.beta[0] - 1) + std::fabs(ps.beta[1] - 0.8) + std::fabs(ps.beta[2] - 0.9) + std::fabs(ps.beta[3] + 1.1), 1e-7);
  const auto os = trade::gravity_sizes(sized, lat, lon, size, trade::GravityMethod::ols_sizes);
  s.check("least squares on logarithms, exact flows: the same", std::fabs(os.beta[0] - 1) + std::fabs(os.beta[1] - 0.8) + std::fabs(os.beta[2] - 0.9) + std::fabs(os.beta[3] + 1.1), 1e-9);

  const auto fn = trade::gravity_effects(noisy, lat, lon);
  s.require("PPML with effects, half the flows zero: the elasticity −1 within four standard errors", std::fabs(fn.distance_elasticity + 1) < 4 * fn.se[0],
            "β = " + calx::verification::exponential(fn.distance_elasticity) + ", se " + calx::verification::exponential(fn.se[0]));
  s.require("... keeping the zeros among its pairs", fn.zeros > fn.pairs / 3 && fn.pairs == n * (n - 1));
  const auto on = trade::gravity_sizes(noisy, lat, lon, size, trade::GravityMethod::ols_sizes);
  s.require("... where least squares on logarithms fits only the positive flows", on.zeros == 0 && on.pairs == fn.pairs - fn.zeros);

  return s.done();
}
