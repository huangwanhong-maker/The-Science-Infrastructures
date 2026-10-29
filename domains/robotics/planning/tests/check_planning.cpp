// Checks of survey planning (planning/survey.hpp; theory/vehicles.md §5), on the Moon's sphere.
//
// - A lawnmower: two waypoints a line; each line its length, along its azimuth; the lines' midpoints their
//   spacing apart on the great circle across them, to the exponential map's curvature (a midpoint lies at
//   tan(ρ)/ρ ≈ 1 + ρ²/3 of its offset, ρ the half-line's angle: 5 × 10⁻⁹ here); alternate lines driven the
//   other way.
// - A ring: every candidate at its radius, the first straight ahead, the third of eight a quarter turn to the left.
// - The best candidate: the greatest value less λ times the cost, the first among equals.

#include <robotics/planning/survey.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <vector>

using namespace robotics;
using namespace physis::geometry::operators;
using vehicles::Vec3;
namespace pg = physis::geometry;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;
  const vehicles::Sphere moon{1737.4e3, 4.9028e12};
  const Vec3 centre = vehicles::SurfacePose::at(0.13, -1.03, 0.0).place(moon);
  const double az = 0.3, L = 400, spacing = 25;
  const int n = 5;
  const auto w = planning::lawnmower(centre, az, L, spacing, n, moon);
  s.require("lawnmower: two waypoints a line", w.size() == std::size_t(2 * n));
  double length = 0, gap = 0, direction = 0;
  const auto axes = planning::local_axes(centre);
  const Vec3 along = std::cos(az) * axes.north + std::sin(az) * axes.east;
  for (int j = 0; j < n; j++) {
    length = std::max(length, std::fabs(vehicles::distance(w[2 * j], w[2 * j + 1], moon) / L - 1));
    // Driven along +along on even lines, −along on odd ones: the bearing at the line's start.
    const double c = pg::dot(vehicles::bearing(w[2 * j], w[2 * j + 1]), along);
    direction = std::max(direction, std::fabs(c - (j % 2 == 0 ? 1.0 : -1.0)));
    if (j > 0) {
      const Vec3 m0 = pg::unit(w[2 * j - 2] + w[2 * j - 1]), m1 = pg::unit(w[2 * j] + w[2 * j + 1]);
      gap = std::max(gap, std::fabs(vehicles::distance(m0, m1, moon) / spacing - 1));
    }
  }
  s.check("lawnmower: each line its length (relative)", length, 1e-6);
  s.check("lawnmower: the lines' midpoints their spacing apart (relative)", gap, 1e-7);
  s.check("lawnmower: each line along the azimuth, alternately one way and the other (1 − |cos|)", direction, 1e-6);

  const auto pose = vehicles::SurfacePose::at(0.13, -1.03, 1.0);
  const auto r = planning::ring(pose, 30, 8, moon);
  double radius = 0;
  for (const auto& c : r) radius = std::max(radius, std::fabs(vehicles::distance(pose.place(moon), c, moon) - 30));
  s.check("ring: every candidate at its radius (m)", radius, 1e-9, "m");
  s.check("ring: the first straight ahead", 1 - pg::dot(vehicles::bearing(pose.place(moon), r[0]), pose.forward()), 1e-12);
  s.check("ring: the third a quarter turn to the left", 1 - pg::dot(vehicles::bearing(pose.place(moon), r[2]), pose.left()), 1e-12);

  const std::vector<Vec3> cand{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}};
  const std::vector<double> value{3, 5, 5, 4}, cost{1, 2, 2, 0};
  const auto index = [&](const Vec3& c) { for (std::size_t i = 0; i < cand.size(); i++) if (c == cand[i]) return i; return std::size_t(0); };
  const std::size_t b0 = planning::best(cand, [&](const Vec3& c) { return value[index(c)]; }, [&](const Vec3& c) { return cost[index(c)]; }, 0.0);
  const std::size_t b1 = planning::best(cand, [&](const Vec3& c) { return value[index(c)]; }, [&](const Vec3& c) { return cost[index(c)]; }, 1.0);
  s.require("best: without cost, the greatest value, the first among equals", b0 == 1);
  s.require("best: with λ = 1, value less cost", b1 == 3);
  (void)pi;
  return s.done();
}
