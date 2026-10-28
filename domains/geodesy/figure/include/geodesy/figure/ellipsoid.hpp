// The figure of the Earth as an ellipsoid of revolution, and geodetic coordinates on it.
//
// A reference ellipsoid is fixed by its equatorial radius a and its flattening f = (a − b)/a;
// the first eccentricity e² = f(2 − f) and the second e′² = e²/(1 − e²) follow (Torge and Müller
// 2012, Geodesy, 4th ed., §4.1). A point is given geodetically by its latitude φ (the angle of
// the ellipsoid's normal with the equator), longitude λ (east of Greenwich, positive east) and
// height h along the normal; its Earth-centred Cartesian coordinates (in the frame the
// ellipsoid is attached to, the ITRS for WGS 84) are
//   X = (N + h) cos φ cos λ,  Y = (N + h) cos φ sin λ,  Z = (N(1 − e²) + h) sin φ,
// N = a / √(1 − e² sin² φ) the radius of curvature in the prime vertical (Torge and Müller,
// eq. 4.24). The inverse has a closed form (Vermeille 2002, J. Geodesy 76, 451–454), here
// followed by one Newton step on the latitude, which takes it to rounding everywhere outside a
// few kilometres of the centre; the checks hold it to ERFA's gc2gd (Fukushima 2006).
//
// The ellipsoids carry their definitions with their sources: WGS 84 (NIMA TR8350.2, 2000, §3),
// GRS 80 (Moritz 2000, J. Geodesy 74, 128), the IERS 2010 conventional ellipsoid (IERS
// Conventions 2010, Table 1.2), and a sphere of the IUGG mean radius R₁ = (2a + b)/3 of WGS 84.

#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

#include <physis/geometry/rotation.hpp>

namespace geodesy::figure {

using physis::geometry::Vec3;

struct Geodetic {
  double latitude = 0, longitude = 0;  // radians; longitude east of the prime meridian
  double height = 0;                   // metres above the ellipsoid, along its normal
};

class Ellipsoid {
 public:
  Ellipsoid(std::string name, double a, double inverse_flattening)
      : name_(std::move(name)), a_(a), f_(inverse_flattening == 0 ? 0.0 : 1 / inverse_flattening) {
    if (!(a > 0) || f_ < 0 || f_ >= 1) throw std::invalid_argument("geodesy: an ellipsoid needs a > 0 and 0 ≤ f < 1");
    e2_ = f_ * (2 - f_);
  }
  static Ellipsoid sphere(std::string name, double radius) { return Ellipsoid(std::move(name), radius, 0); }

  // WGS 84: a = 6 378 137 m, 1/f = 298.257223563 (NIMA TR8350.2, Table 3.1).
  static Ellipsoid wgs84() { return {"WGS 84", 6378137.0, 298.257223563}; }
  // GRS 80: a = 6 378 137 m, 1/f = 298.257222101 (Moritz 2000).
  static Ellipsoid grs80() { return {"GRS 80", 6378137.0, 298.257222101}; }
  // The IERS 2010 conventional ellipsoid: a = 6 378 136.6 m, 1/f = 298.25642.
  static Ellipsoid iers2010() { return {"IERS 2010", 6378136.6, 298.25642}; }
  // The Earth as a sphere of WGS 84's mean radius R₁ = 6 371 008.8 m (Moritz 2000).
  static Ellipsoid earth_sphere() { return sphere("Earth, mean sphere", 6371008.8); }

  const std::string& name() const { return name_; }
  double a() const { return a_; }
  double b() const { return a_ * (1 - f_); }
  double flattening() const { return f_; }
  double e2() const { return e2_; }
  double second_e2() const { return e2_ / (1 - e2_); }
  double mean_radius() const { return (2 * a_ + b()) / 3; }
  bool is_sphere() const { return f_ == 0; }

  // The radii of curvature in the prime vertical (N) and the meridian (M) at a latitude.
  double prime_vertical_radius(double phi) const { return a_ / std::sqrt(1 - e2_ * std::sin(phi) * std::sin(phi)); }
  double meridian_radius(double phi) const {
    const double w2 = 1 - e2_ * std::sin(phi) * std::sin(phi);
    return a_ * (1 - e2_) / (w2 * std::sqrt(w2));
  }

  Vec3 cartesian(const Geodetic& p) const {
    const double sp = std::sin(p.latitude), cp = std::cos(p.latitude);
    const double n = a_ / std::sqrt(1 - e2_ * sp * sp);
    return {(n + p.height) * cp * std::cos(p.longitude), (n + p.height) * cp * std::sin(p.longitude), (n * (1 - e2_) + p.height) * sp};
  }

  // The geodetic coordinates of an Earth-centred point (Vermeille 2002, then a Newton step).
  Geodetic geodetic(const Vec3& x) const {
    const double rho2 = x[0] * x[0] + x[1] * x[1], rho = std::sqrt(rho2), z = x[2];
    Geodetic g;
    g.longitude = (rho2 > 0) ? std::atan2(x[1], x[0]) : 0.0;
    if (is_sphere()) {
      g.latitude = std::atan2(z, rho);
      g.height = std::hypot(rho, z) - a_;
      return g;
    }
    const double a2 = a_ * a_, e4 = e2_ * e2_;
    const double p = rho2 / a2, q = (1 - e2_) / a2 * z * z, r = (p + q - e4) / 6;
    double phi;
    if (r > 0 || q > 0) {
      // Vermeille (2002), valid outside the evolute of the meridian ellipse (all but within
      // about 45 km of the centre).
      const double s = e4 * p * q / (4 * r * r * r);
      const double t = std::cbrt(1 + s + std::sqrt(s * (2 + s)));
      const double u = r * (1 + t + 1 / t);
      const double v = std::sqrt(u * u + e4 * q);
      const double w = e2_ * (u + v - q) / (2 * v);
      const double k = std::sqrt(u + v + w * w) - w;
      const double d = k * rho / (k + e2_);
      phi = 2 * std::atan2(z, d + std::hypot(d, z));
    } else {
      phi = std::atan2(z, rho);
    }
    // One Newton step on f(φ) = ρ sin φ − z cos φ − e² N sin φ cos φ = 0, whose root is the
    // geodetic latitude (the point's offset along the normal is perpendicular to the meridian).
    for (int i = 0; i < 2; i++) {
      const double sp = std::sin(phi), cp = std::cos(phi);
      const double w2 = 1 - e2_ * sp * sp, n = a_ / std::sqrt(w2);
      const double f = rho * sp - z * cp - e2_ * n * sp * cp;
      const double df = rho * cp + z * sp - e2_ * n * (cp * cp - sp * sp) - e2_ * e2_ * n * sp * sp * cp * cp / w2;
      if (df == 0) break;
      phi -= f / df;
    }
    g.latitude = phi;
    const double sp = std::sin(phi), cp = std::cos(phi);
    const double n = a_ / std::sqrt(1 - e2_ * sp * sp);
    // Height along the normal, well conditioned at every latitude.
    g.height = rho * cp + z * sp - a_ * a_ / n;
    return g;
  }

  // The outward unit normal at a geodetic latitude and longitude (the local vertical, "up").
  static Vec3 normal(double phi, double lambda) {
    return {std::cos(phi) * std::cos(lambda), std::cos(phi) * std::sin(lambda), std::sin(phi)};
  }
  // The geocentric latitude of a point on the ellipsoid at geodetic latitude φ.
  double geocentric_latitude(double phi) const { return std::atan((1 - e2_) * std::tan(phi)); }

 private:
  std::string name_;
  double a_ = 0, f_ = 0, e2_ = 0;
};

}  // namespace geodesy::figure
