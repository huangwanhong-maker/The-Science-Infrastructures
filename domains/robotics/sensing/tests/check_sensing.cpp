// Checks of sensors on a vehicle (sensing/mount.hpp; theory/vehicles.md §6), on the Moon's sphere.
//
// - A mast 1.5 m up puts its sensor 1.5 m above the sphere, over the vehicle.
// - A sensor turned a quarter turn about the vehicle's up reads the vehicle's forward as its −y and its
//   left as its +x; a world vector to the sensor's axes and back is itself.
// - White noise: sample k is the same whenever it is drawn; a hundred thousand samples have mean 0 and
//   standard deviation sd, and two channels are uncorrelated.

#include <robotics/sensing/mount.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>

using namespace robotics;
using namespace physis::geometry::operators;
using vehicles::Rotation;
using vehicles::Vec3;
namespace pg = physis::geometry;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;
  const vehicles::Sphere moon{1737.4e3, 4.9028e12};
  const auto pose = vehicles::SurfacePose::at(0.13, -1.03, 0.7);

  const sensing::Mount mast{{0, 0, 1.5}, Rotation{}};
  const auto p = sensing::sensor_pose(pose, moon, mast);
  s.check("mount: a mast 1.5 m up puts the sensor 1.5 m above the sphere (m)", std::fabs(pg::norm(p.place) - moon.radius - 1.5), 1e-9, "m");
  s.check("mount: over the vehicle", pg::norm(pg::unit(p.place) - pose.up()), 1e-15);

  const sensing::Mount turned{{0.2, 0, 1.0}, Rotation::about(pg::Axis::z, pi / 2)};
  const auto q = sensing::sensor_pose(pose, moon, turned);
  s.check("mount: turned a quarter about up, the vehicle's forward reads as the sensor's −y", pg::norm(q.to_sensor(pose.forward()) - Vec3{0, -1, 0}), 1e-15);
  s.check("mount: and its left as the sensor's +x", pg::norm(q.to_sensor(pose.left()) - Vec3{1, 0, 0}), 1e-15);
  const Vec3 B{12.5, -3.25, 40.0};
  s.check("mount: a vector to the sensor's axes and back is itself", pg::norm(q.from_sensor(q.to_sensor(B)) - B), 1e-13);

  const sensing::WhiteNoise noise{0.5, 20261029, 7};
  s.require("noise: sample k is the same whenever it is drawn", noise.sample(12345) == noise.sample(12345) && noise.vector(3)[2] == noise.vector(3)[2]);
  const int n = 100000;
  double sum = 0, sq = 0, cross = 0;
  const sensing::WhiteNoise other{0.5, 20261029, 8};
  for (int k = 0; k < n; k++) {
    const double a = noise.sample(std::uint64_t(k)), b = other.sample(std::uint64_t(k));
    sum += a, sq += a * a, cross += a * b;
  }
  // Standard errors: sd/√n = 0.0016 for the mean, sd/√(2n) ≈ 0.0011 for the sd, sd²/√n ≈ 0.0008 for the cross term.
  s.check("noise: mean 0 over 10⁵ samples", std::fabs(sum / n), 0.007);
  s.check("noise: standard deviation sd (relative)", std::fabs(std::sqrt(sq / n) / 0.5 - 1), 0.01);
  s.check("noise: two channels uncorrelated (their mean product)", std::fabs(cross / n), 0.0035);
  return s.done();
}
