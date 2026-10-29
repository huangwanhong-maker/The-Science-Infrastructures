// Planning a survey on a sphere: a fixed pattern of lines, and the greedy choice of the next place to
// measure by what is not yet known there (theory/vehicles.md §5).
//
// A lawnmower (boustrophedon) pattern covers a rectangle with parallel lines, driven alternately one way and
// the other: laid out in the tangent plane at its centre (north, east), the lines along an azimuth, and
// carried onto the sphere by the exponential map at the centre (surface.hpp), which keeps distances from the
// centre and their directions exactly; over a rectangle a hundredth of the radius across, the lines' spacing
// changes by less than 10⁻⁴ of itself.
//
// An adaptive survey chooses each next place from candidates (a ring about the vehicle, say) by a value less
// a cost: with the value the posterior variance of the field to be mapped there, this is greedy maximum
// entropy sampling, which for a Gaussian field is uncertainty sampling (Shewry and Wynn 1987, J. Appl.
// Statist. 14, 165–170); the cost, the distance to drive, weighted by λ. The value is the caller's: the
// planner names no field.

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <robotics/vehicles/surface.hpp>

namespace robotics::planning {

using robotics::vehicles::Sphere;
using robotics::vehicles::SurfacePose;
using robotics::vehicles::Vec3;

// The tangent unit vectors north and east at a place (not at a pole).
struct LocalAxes {
  Vec3 north, east;
};
inline LocalAxes local_axes(const Vec3& place) {
  const Vec3 u = physis::geometry::unit(place);
  const double h = std::hypot(u[0], u[1]);
  if (h < 1e-12) throw std::domain_error("robotics::planning::local_axes: north and east are undefined at a pole");
  const Vec3 east{-u[1] / h, u[0] / h, 0};
  return {physis::geometry::cross(u, east), east};
}

// `lines` lines of `length` metres, `spacing` apart, centred at `centre`, along the azimuth (radians from
// north, clockwise), driven alternately: the waypoints at the lines' ends, in the order of driving.
inline std::vector<Vec3> lawnmower(const Vec3& centre, double azimuth, double length, double spacing, int lines, const Sphere& s) {
  using namespace physis::geometry::operators;
  if (lines < 1 || !(length > 0) || !(spacing >= 0)) throw std::invalid_argument("robotics::planning::lawnmower: one line or more, of positive length");
  const LocalAxes a = local_axes(centre);
  const Vec3 along = std::cos(azimuth) * a.north + std::sin(azimuth) * a.east;
  const Vec3 across = std::cos(azimuth + std::numbers::pi / 2) * a.north + std::sin(azimuth + std::numbers::pi / 2) * a.east;
  std::vector<Vec3> w;
  for (int j = 0; j < lines; j++) {
    const double offset = (j - 0.5 * (lines - 1)) * spacing;
    const double from = (j % 2 == 0 ? -0.5 : 0.5) * length;
    w.push_back(robotics::vehicles::exponential(centre, offset * across + from * along, s));
    w.push_back(robotics::vehicles::exponential(centre, offset * across - from * along, s));
  }
  return w;
}

// `count` places on a circle of `radius` metres about a pose, the first straight ahead, counter-clockwise.
inline std::vector<Vec3> ring(const SurfacePose& pose, double radius, int count, const Sphere& s) {
  using namespace physis::geometry::operators;
  std::vector<Vec3> c;
  for (int i = 0; i < count; i++) {
    const double a = 2 * std::numbers::pi * i / count;
    c.push_back(robotics::vehicles::exponential(pose.place(s), radius * (std::cos(a) * pose.forward() + std::sin(a) * pose.left()), s));
  }
  return c;
}

// The index of the candidate of greatest value(c) − λ cost(c), the first among equals.
template <class Value, class Cost>
std::size_t best(const std::vector<Vec3>& candidates, Value&& value, Cost&& cost, double lambda) {
  if (candidates.empty()) throw std::invalid_argument("robotics::planning::best: one candidate or more");
  std::size_t k = 0;
  double top = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < candidates.size(); i++) {
    const double score = value(candidates[i]) - lambda * cost(candidates[i]);
    if (score > top) top = score, k = i;
  }
  return k;
}

}  // namespace robotics::planning
