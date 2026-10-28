// The transverse Mercator projection of an ellipsoid of revolution, by Krüger's series in the third
// flattening n to sixth order (Krüger 1912; Karney 2011, "Transverse Mercator with an accuracy of a
// few nanometers", J. Geodesy 85, 475–485).
//
// The projection is conformal: it maps the ellipsoid onto the plane with the central meridian λ₀
// true to scale k₀. It is composed of two conformal maps (Karney 2011, §2–3):
//   1. the ellipsoid onto the sphere by the conformal latitude χ, τ′ = tan χ from τ = tan φ by
//        τ′ = τ √(1 + σ²) − σ √(1 + τ²),   σ = sinh(e atanh(e τ / √(1 + τ²)))           (eqs. 7–9);
//   2. the sphere's Gauss–Schreiber transverse Mercator, ξ′ = atan2(τ′, cos λ),
//        η′ = asinh(sin λ / √(τ′² + cos² λ))                                                (eq. 10),
//   then the series that carries (ξ′, η′) to the meridian distance's rectified coordinates,
//        ξ = ξ′ + Σⱼ αⱼ sin 2jξ′ cosh 2jη′,   η = η′ + Σⱼ αⱼ cos 2jξ′ sinh 2jη′           (eq. 11),
//   and x = k₀ A η, y = k₀ A ξ, A = a/(1 + n) (1 + n²/4 + n⁴/64 + n⁶/256) the rectifying radius
//   (eq. 14). The inverse runs the series βⱼ back (eq. 15) and recovers τ from τ′ by Newton's
//   method (eqs. 19–21). The meridian convergence γ (the bearing of grid north, clockwise from
//   true north) and the point scale k follow from the derivatives of the maps (eqs. 12–13, 16–17).
//
// The coefficients αⱼ, βⱼ are Karney's eqs. 35–36 to n⁶; the checks hold the projection to PROJ's
// (pyproj 3.8, PROJ 9.8: its extended transverse Mercator of Poder and Engsager) and its inverse
// to its forward. Within 3900 km of the central meridian the error of the series is below 5 nm
// (Karney 2011, Fig. 1); beyond, it grows, and the projection is not meant for it (UTM's zones are
// 6° wide).

#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include <geodesy/figure/ellipsoid.hpp>

namespace geodesy::projection {

using geodesy::figure::Ellipsoid;

// A point of the plane, with the projection's local properties there.
struct GridPoint {
  double easting = 0, northing = 0;  // m, before false origins
  double convergence = 0;            // rad: the bearing of grid north, clockwise from true north (Karney's γ, PROJ's convergence)
  double scale = 1;                  // the point scale k
};
struct GeoPoint {
  double latitude = 0, longitude = 0;  // rad
  double convergence = 0, scale = 1;
};

class TransverseMercator {
 public:
  TransverseMercator(const Ellipsoid& e, double k0) : a_(e.a()), f_(e.flattening()), k0_(k0) {
    if (!(k0 > 0)) throw std::invalid_argument("geodesy: a transverse Mercator needs k0 > 0");
    e2_ = f_ * (2 - f_);
    es_ = std::sqrt(e2_);
    const double n = f_ / (2 - f_), n2 = n * n, n3 = n2 * n, n4 = n3 * n, n5 = n4 * n, n6 = n5 * n;
    A_ = a_ / (1 + n) * (1 + n2 / 4 + n4 / 64 + n6 / 256);
    // Karney 2011, eq. 35: from the conformal sphere to the rectified plane.
    alpha_ = {n / 2 - 2 * n2 / 3 + 5 * n3 / 16 + 41 * n4 / 180 - 127 * n5 / 288 + 7891 * n6 / 37800,
              13 * n2 / 48 - 3 * n3 / 5 + 557 * n4 / 1440 + 281 * n5 / 630 - 1983433 * n6 / 1935360,
              61 * n3 / 240 - 103 * n4 / 140 + 15061 * n5 / 26880 + 167603 * n6 / 181440,
              49561 * n4 / 161280 - 179 * n5 / 168 + 6601661 * n6 / 7257600,
              34729 * n5 / 80640 - 3418889 * n6 / 1995840,
              212378941 * n6 / 319334400};
    // Karney 2011, eq. 36: back.
    beta_ = {n / 2 - 2 * n2 / 3 + 37 * n3 / 96 - n4 / 360 - 81 * n5 / 512 + 96199 * n6 / 604800,
             n2 / 48 + n3 / 15 - 437 * n4 / 1440 + 46 * n5 / 105 - 1118711 * n6 / 3870720,
             17 * n3 / 480 - 37 * n4 / 840 - 209 * n5 / 4480 + 5569 * n6 / 90720,
             4397 * n4 / 161280 - 11 * n5 / 504 - 830251 * n6 / 7257600,
             4583 * n5 / 161280 - 108847 * n6 / 3991680,
             20648693 * n6 / 638668800};
  }

  double k0() const { return k0_; }
  double rectifying_radius() const { return A_; }

  // The conformal latitude's tangent τ′ from the geodetic latitude's τ (Karney 2011, eqs. 7–9).
  double taup(double tau) const {
    const double tau1 = std::hypot(1.0, tau);
    const double sig = std::sinh(es_ * std::atanh(es_ * tau / tau1));
    return std::hypot(1.0, sig) * tau - sig * tau1;
  }
  // Its inverse by Newton's method from τ = τ′/(1 − e²) (eqs. 19–21); a few steps to rounding.
  double tau(double taup) const {
    double t = taup / (1 - e2_);
    const double stol = 0.1 * std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::fabs(taup));
    for (int i = 0; i < 6; i++) {
      const double tp = this->taup(t);
      const double dtau = (taup - tp) * (1 + (1 - e2_) * t * t) / ((1 - e2_) * std::hypot(1.0, t) * std::hypot(1.0, tp));
      t += dtau;
      if (std::fabs(dtau) < stol) break;
    }
    return t;
  }

  // The projection of a point at latitude φ and longitude λ from the central meridian (rad).
  GridPoint forward(double phi, double lambda) const {
    const double tau = std::tan(phi), tp = taup(tau);
    const double c = std::cos(lambda), s = std::sin(lambda);
    double xip, etap;
    if (std::fabs(phi) >= std::numbers::pi / 2 - 1e-15) {  // the poles
      xip = std::copysign(std::numbers::pi / 2, phi);
      etap = 0;
    } else {
      xip = std::atan2(tp, c);
      etap = std::asinh(s / std::hypot(tp, c));
    }
    double xi = xip, eta = etap, p = 1, q = 0;
    for (int j = 1; j <= 6; j++) {
      const double a = alpha_[j - 1], s2 = std::sin(2 * j * xip), c2 = std::cos(2 * j * xip), sh = std::sinh(2 * j * etap), ch = std::cosh(2 * j * etap);
      xi += a * s2 * ch;
      eta += a * c2 * sh;
      p += 2 * j * a * c2 * ch;
      q += 2 * j * a * s2 * sh;
    }
    GridPoint g;
    g.easting = k0_ * A_ * eta;
    g.northing = k0_ * A_ * xi;
    // Convergence and scale (eqs. 12–13): the sphere's and the series'.
    const double gamma1 = std::atan2(tp * s, std::hypot(1.0, tp) * c);  // atan(τ′ tan λ / √(1 + τ′²)), quadrant-safe
    const double gamma2 = std::atan2(q, p);
    const double k1 = std::sqrt(1 - e2_ * std::sin(phi) * std::sin(phi)) * std::hypot(1.0, tau) / std::hypot(tp, c);
    const double k2 = A_ / a_ * std::hypot(p, q);
    g.convergence = gamma1 + gamma2;
    g.scale = k0_ * k1 * k2;
    return g;
  }

  // The inverse: latitude and longitude from the central meridian (rad) of a point of the plane.
  GeoPoint reverse(double x, double y) const {
    const double xi = y / (k0_ * A_), eta = x / (k0_ * A_);
    double xip = xi, etap = eta, p = 1, q = 0;
    for (int j = 1; j <= 6; j++) {
      const double b = beta_[j - 1], s2 = std::sin(2 * j * xi), c2 = std::cos(2 * j * xi), sh = std::sinh(2 * j * eta), ch = std::cosh(2 * j * eta);
      xip -= b * s2 * ch;
      etap -= b * c2 * sh;
      p -= 2 * j * b * c2 * ch;
      q -= 2 * j * b * s2 * sh;
    }
    const double sxp = std::sin(xip), cxp = std::cos(xip), shp = std::sinh(etap);
    const double r = std::hypot(shp, cxp);
    GeoPoint g;
    if (r == 0) {  // a pole
      g.latitude = std::copysign(std::numbers::pi / 2, xip);
      g.longitude = 0;
      g.convergence = 0;
      g.scale = k0_ * A_ / a_ * std::hypot(p, q) * std::sqrt(1 - e2_) * std::exp(es_ * std::atanh(es_));
      return g;
    }
    const double tp = sxp / r, t = tau(tp);
    g.latitude = std::atan(t);
    g.longitude = std::atan2(shp, cxp);
    // Convergence and scale (eqs. 16–17).
    const double gamma1 = std::atan2(sxp * std::tanh(etap), cxp);
    const double gamma2 = std::atan2(q, p);
    g.convergence = gamma1 - gamma2;
    const double sphi = t / std::hypot(1.0, t);
    const double k1 = std::sqrt(1 - e2_ * sphi * sphi) * std::hypot(1.0, t) * r;
    const double k2 = A_ / a_ / std::hypot(p, q);
    g.scale = k0_ * k1 * k2;
    return g;
  }

 private:
  double a_, f_, k0_, e2_, es_, A_;
  std::array<double, 6> alpha_{}, beta_{};
};

}  // namespace geodesy::projection
