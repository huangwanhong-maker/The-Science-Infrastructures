// A vehicle's place on a sphere: its pose as one rotation, and the geometry of the sphere it moves on
// (theory/vehicles.md §1).
//
// The pose. A vehicle standing on a sphere has a place p (|p| = R) and a heading e, a unit vector tangent
// there. The pair is one rotation A ∈ SO(3), whose columns are the vehicle's axes in the sphere's frame:
// forward e, left n = û × e and up û = p/R (the body axes x forward, y left, z up of ROS's REP 103). SO(3)
// acts simply transitively on the unit tangent vectors of S², so every (place, heading) is one rotation and
// every rotation one (place, heading): A is the configuration, without a constraint to keep.
//
// Its motion. Driving forward at speed v while turning left at ω turns the axes by the body angular
// velocity w = (0, v/R, ω): Ȧ = A w^ (d(up)/dt = (v/R) forward, d(forward)/dt = ω left − (v/R) up). For v
// and ω constant over a step the motion is exactly A(t + h) = A(t) exp(h w^): a small circle on the sphere
// of geodesic curvature ω/v, a great circle when ω = 0.
//
// The sphere. Its radius R and its gravitational parameter GM (gravity GM/R² at its surface), in the frame of
// the body it is the figure of (for the Moon its mean-Earth frame); a terrain's heights are not yet carried.
// Distances along it are R times the angle between places; the bearing from p to q is the direction of the
// great circle from p to q, (q̂ − (q̂·p̂) p̂) normalised, and the exponential map exp_p(d) = R (cos(|d|/R) p̂ +
// sin(|d|/R) d̂) carries a tangent displacement d at p to the place it reaches along the great circle. From
// unit vectors a bearing between places d apart carries the rounding R ε/d (10⁻¹¹ rad at 50 m on the Moon).

#pragma once

#include <cmath>
#include <stdexcept>

#include <physis/geometry/rotation.hpp>

namespace robotics::vehicles {

using physis::geometry::Mat3;
using physis::geometry::Rotation;
using physis::geometry::Vec3;

struct Sphere {
  double radius = 0;  // R, m
  double gm = 0;      // GM, m³/s²
  double gravity() const { return gm / (radius * radius); }
};

// Where a vehicle is and which way it faces: body axes x forward, y left, z up, in the sphere's frame.
class SurfacePose {
 public:
  SurfacePose() = default;
  explicit SurfacePose(const Rotation& attitude) : A_(attitude) {}

  // The pose at a place (any vector along it), facing the tangent direction nearest `toward`.
  static SurfacePose at(const Vec3& place, const Vec3& toward) {
    using namespace physis::geometry::operators;
    const Vec3 up = physis::geometry::unit(place);
    const Vec3 t = toward - physis::geometry::dot(toward, up) * up;
    const Vec3 forward = physis::geometry::unit(t);
    const Vec3 left = physis::geometry::cross(up, forward);
    const Mat3 m{{{forward[0], left[0], up[0]}, {forward[1], left[1], up[1]}, {forward[2], left[2], up[2]}}};
    return SurfacePose(Rotation::from_matrix(m));
  }
  // The pose at a latitude and an east longitude (radians), heading at an azimuth from north, clockwise
  // towards east (radians).
  static SurfacePose at(double latitude, double longitude, double azimuth) {
    const double cp = std::cos(latitude), sp = std::sin(latitude), cl = std::cos(longitude), sl = std::sin(longitude);
    const Vec3 up{cp * cl, cp * sl, sp}, north{-sp * cl, -sp * sl, cp}, east{-sl, cl, 0};
    using namespace physis::geometry::operators;
    return at(up, std::cos(azimuth) * north + std::sin(azimuth) * east);
  }

  const Rotation& attitude() const { return A_; }
  Vec3 up() const { return A_.apply({0, 0, 1}); }
  Vec3 forward() const { return A_.apply({1, 0, 0}); }
  Vec3 left() const { return A_.apply({0, 1, 0}); }
  Vec3 place(const Sphere& s) const {
    using namespace physis::geometry::operators;
    return s.radius * up();
  }
  double latitude() const {
    const Vec3 u = up();
    return std::atan2(u[2], std::hypot(u[0], u[1]));
  }
  double longitude() const {
    const Vec3 u = up();
    return std::atan2(u[1], u[0]);
  }
  // The heading's azimuth from north, clockwise towards east, in (−π, π].
  double azimuth() const {
    const Vec3 u = up(), f = forward();
    const Vec3 east{-u[1], u[0], 0};
    const double e = std::hypot(east[0], east[1]);
    if (e < 1e-15) return 0;  // at a pole north is undefined
    const Vec3 E{east[0] / e, east[1] / e, 0};
    const Vec3 N = physis::geometry::cross(u, E);
    return std::atan2(physis::geometry::dot(f, E), physis::geometry::dot(f, N));
  }

  // Driving the distance d forward while turning by the angle ψ, uniformly: A exp((0, d/R, ψ)^).
  SurfacePose moved(double distance, double turn, const Sphere& s) const { return SurfacePose(A_ * Rotation::exp(Vec3{0, distance / s.radius, turn})); }

 private:
  Rotation A_;
};

// The distance along the sphere between two places (any vectors along them), m.
inline double distance(const Vec3& a, const Vec3& b, const Sphere& s) {
  const Vec3 c = physis::geometry::cross(a, b);
  return s.radius * std::atan2(physis::geometry::norm(c), physis::geometry::dot(a, b));
}

// The unit tangent at a from which the great circle leads to b (b ≠ ±a).
inline Vec3 bearing(const Vec3& a, const Vec3& b) {
  using namespace physis::geometry::operators;
  const Vec3 u = physis::geometry::unit(a), v = physis::geometry::unit(b);
  const Vec3 t = v - physis::geometry::dot(v, u) * u;
  return physis::geometry::unit(t);
}

// exp_p(d): the place reached from p along the great circle of the tangent displacement d (m).
inline Vec3 exponential(const Vec3& p, const Vec3& d, const Sphere& s) {
  using namespace physis::geometry::operators;
  const Vec3 u = physis::geometry::unit(p);
  const double length = physis::geometry::norm(d);
  if (length == 0) return s.radius * u;
  const double angle = length / s.radius;
  return s.radius * (std::cos(angle) * u + std::sin(angle) * ((1 / length) * d));
}

}  // namespace robotics::vehicles
