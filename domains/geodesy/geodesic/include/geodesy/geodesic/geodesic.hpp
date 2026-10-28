// Geodesics on an ellipsoid of revolution: the direct and the inverse problem, and points along a
// geodesic, to the round-off of double precision (Karney 2013).
//
// A geodesic is the shortest path on the surface between two of its points (or a locally shortest
// one, prolonged). On an ellipsoid of revolution with equatorial radius a and flattening f it is
// mapped, by Bessel's method, onto a great circle of an auxiliary sphere: the reduced latitude β
// (tan β = (1 − f) tan φ) and the spherical arc σ along the great circle stand for the latitude and
// the distance, and Clairaut's constant sin α₀ = sin α cos β is kept along the line. The distance
// and the longitude are then integrals over σ, expanded in the small parameter
//   ε = (√(1 + k²) − 1) / (√(1 + k²) + 1),   k = e′ cos α₀,
// and summed by Clenshaw's method (Karney 2013, eqs. 7–9, 15–18, 23–25, 40–42). The direct problem
// (from a point, an azimuth and a distance, the end point) is then explicit: σ from the distance by
// the reverted series of eq. 20, and the end point from the auxiliary sphere (eqs. 5–6, 8). The
// inverse problem (the geodesic between two points) solves λ₁₂(α₁) = λ₁₂ for the azimuth at the
// first point by Newton's method, whose derivative is the reduced length (eq. 46), from a starting
// value on the sphere or, for nearly antipodal points, from the astroid of eq. 55 (§5–6); the
// iteration is bracketed and falls back on bisection, so that it converges for every pair. The
// reduced length m₁₂, the geodesic scales M₁₂, M₂₁ (eqs. 38–39) and the area S₁₂ between the
// geodesic and the equator (eqs. 58–60) come with the solution. The series are carried to the
// sixth order in ε and n = f/(2 − f), which keeps the error within 15 nm on WGS 84.
//
// This is a port of GeographicLib's Geodesic and GeodesicLine (its Python transcription, version
// 2.1), by Charles F. F. Karney, kept close to it, line for line, with the same reductions of
// angles in degrees (exact remainders and quadrants before any radian appears), the same rounding
// of tiny angles (AngRound) and the same safeguards; its notice follows. One change: in the
// lengths needed for Newton's derivative alone, every coefficient of the combined series
// A₁C₁ₗ − A₂C₂ₗ is formed (l = 1..6, as in the C++ library), where the Python transcription forms
// l = 1..5; the root is unchanged, only the derivative is the exact one.
//
//   Copyright (c) Charles Karney (2011-2022) <karney@alum.mit.edu> and licensed under the MIT/X11
//   License. For more information, see https://geographiclib.sourceforge.io/
//
//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software
//   and associated documentation files (the "Software"), to deal in the Software without
//   restriction, including without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions: The above copyright
//   notice and this permission notice shall be included in all copies or substantial portions of
//   the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//   PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
//   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
//   IN THE SOFTWARE.
//
// Conventions (the infrastructure's, CLAUDE.md §4.1). Angles at the interface are radians:
// latitudes geodetic, longitudes east of the prime meridian, azimuths clockwise from north; inside,
// angles are reduced in degrees as GeographicLib does, for exactness. Distances and lengths in
// metres, areas in square metres. A longitude returned is reduced to [−π, π] unless it is asked to
// be unrolled (counted continuously along the line). A latitude within 10⁻⁹ degree beyond a pole
// (the rounding of ±π/2 converted to degrees) is taken to be the pole; one farther beyond is not a
// latitude and gives NaN.
//
// Reference: C. F. F. Karney, "Algorithms for geodesics", J. Geodesy 87, 43–55 (2013),
// doi:10.1007/s00190-012-0578-z, with its addenda (geographiclib.sourceforge.io/geod-addenda.html).

#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <geodesy/figure/ellipsoid.hpp>

namespace geodesy::geodesic {

// ---- the arithmetic of angles in degrees (GeographicLib's Math) -------------------------------------

namespace math {

inline constexpr double nan = std::numeric_limits<double>::quiet_NaN();
inline constexpr double deg_to_rad = std::numbers::pi / 180.0;
inline constexpr double rad_to_deg = 180.0 / std::numbers::pi;

inline double sq(double x) { return x * x; }

// The unit vector of (x, y).
inline std::pair<double, double> norm(double x, double y) {
  const double r = std::hypot(x, y);
  return {x / r, y / r};
}

// u + v = s + t exactly, s the rounded sum (Knuth's error-free transformation).
inline std::pair<double, double> sum(double u, double v) {
  const volatile double s = u + v;
  volatile double up = s - v;
  volatile double vpp = s - up;
  up = up - u;
  vpp = vpp - v;
  const double t = s == 0 ? double(s) : 0.0 - (up + vpp);
  return {s, t};
}

// Horner's evaluation of the polynomial of degree N with coefficients p[s..s+N], highest first.
inline double polyval(int N, const double* p, int s, double x) {
  double y = N < 0 ? 0.0 : p[s];
  while (N > 0) {
    --N;
    ++s;
    y = y * x + p[s];
  }
  return y;
}

// An angle rounded so that tiny values underflow to zero: the smallest gap is then 2⁻⁵⁷ degree,
// 0.7 pm on the Earth, which removes nearly singular cases.
inline double ang_round(double x) {
  constexpr double z = 1 / 16.0;
  volatile double y = std::fabs(x);
  if (y < z) {
    volatile double w = z - y;
    y = z - w;
  }
  return std::copysign(double(y), x);
}

inline double remainder(double x, double y) { return std::isfinite(x) ? std::remainder(x, y) : nan; }

// An angle in degrees reduced to [−180, 180].
inline double ang_normalize(double x) {
  const double y = remainder(x, 360);
  return std::fabs(y) == 180 ? std::copysign(180.0, x) : y;
}

inline double lat_fix(double x) { return std::fabs(x) > 90 ? nan : x; }

// y − x reduced to [−180, 180], with the rounding error of the reduction in the second member.
inline std::pair<double, double> ang_diff(double x, double y) {
  auto [d, t] = sum(remainder(-x, 360), remainder(y, 360));
  auto [d2, t2] = sum(remainder(d, 360), t);
  d = d2;
  t = t2;
  if (d == 0 || std::fabs(d) == 180) d = std::copysign(d, t == 0 ? y - x : -t);
  return {d, t};
}

namespace detail {
inline std::pair<double, double> quadrant(double s, double c, int q) {
  q = ((q % 4) + 4) % 4;
  switch (q) {
    case 1: return {c, -s};
    case 2: return {-s, -c};
    case 3: return {-c, s};
    default: return {s, c};
  }
}
}  // namespace detail

// The sine and cosine of x degrees, exact at the multiples of 90 degrees.
inline std::pair<double, double> sincosd(double x) {
  double r = std::isfinite(x) ? std::fmod(x, 360.0) : nan;
  const int q = std::isnan(r) ? 0 : int(std::nearbyint(r / 90));  // rounding half to even, as Python's round
  r -= 90 * q;
  r *= deg_to_rad;
  auto [s, c] = detail::quadrant(std::sin(r), std::cos(r), q);
  c = c + 0.0;
  if (s == 0) s = std::copysign(s, x);
  return {s, c};
}

// The sine and cosine of x + t degrees, x in [−180, 180] and t its small correction.
inline std::pair<double, double> sincosde(double x, double t) {
  const int q = std::isfinite(x) ? int(std::nearbyint(x / 90)) : 0;
  double r = x - 90 * q;
  r = ang_round(r + t) * deg_to_rad;
  auto [s, c] = detail::quadrant(std::sin(r), std::cos(r), q);
  c = c + 0.0;
  if (s == 0) s = std::copysign(s, x);
  return {s, c};
}

// atan2(y, x) in degrees, exact at the quadrants.
inline double atan2d(double y, double x) {
  int q = 0;
  if (std::fabs(y) > std::fabs(x)) {
    q = 2;
    std::swap(x, y);
  }
  if (x < 0) {
    q += 1;
    x = -x;
  }
  double ang = std::atan2(y, x) * rad_to_deg;
  if (q == 1) ang = std::copysign(180.0, y) - ang;
  else if (q == 2) ang = 90 - ang;
  else if (q == 3) ang = -90 + ang;
  return ang;
}

}  // namespace math

// ---- what the problems answer -------------------------------------------------------------------

// The inverse problem: the geodesic between two points.
struct Inverse {
  double s12 = math::nan;   // distance (m)
  double azi1 = math::nan;  // azimuth at the first point (rad, clockwise from north)
  double azi2 = math::nan;  // azimuth at the second point, forward along the line (rad)
  double a12 = math::nan;   // arc length on the auxiliary sphere (rad)
  double m12 = math::nan;   // reduced length (m)
  double M12 = math::nan;   // geodesic scale of the second point relative to the first
  double M21 = math::nan;   // geodesic scale of the first point relative to the second
  double S12 = math::nan;   // area between the geodesic and the equator (m²)
};

// The direct problem: the end of a geodesic from a point, an azimuth and a distance.
struct Direct {
  double lat2 = math::nan, lon2 = math::nan;  // the end point (rad)
  double azi2 = math::nan;                    // azimuth there (rad)
  double s12 = math::nan;                     // distance (m)
  double a12 = math::nan;                     // arc length on the auxiliary sphere (rad)
  double m12 = math::nan, M12 = math::nan, M21 = math::nan, S12 = math::nan;
};

class GeodesicLine;

// The capabilities of a computation, as GeographicLib's bit masks (internal).
namespace caps {
inline constexpr unsigned C1 = 1u << 0, C1p = 1u << 1, C2 = 1u << 2, C3 = 1u << 3, C4 = 1u << 4, ALL_CAPS = 0x1F;
inline constexpr unsigned OUT_ALL = 0x7F80, OUT_MASK = 0xFF80;
inline constexpr unsigned EMPTY = 0;
inline constexpr unsigned LATITUDE = 1u << 7;
inline constexpr unsigned LONGITUDE = (1u << 8) | C3;
inline constexpr unsigned AZIMUTH = 1u << 9;
inline constexpr unsigned DISTANCE = (1u << 10) | C1;
inline constexpr unsigned STANDARD = LATITUDE | LONGITUDE | AZIMUTH | DISTANCE;
inline constexpr unsigned DISTANCE_IN = (1u << 11) | C1 | C1p;
inline constexpr unsigned REDUCEDLENGTH = (1u << 12) | C1 | C2;
inline constexpr unsigned GEODESICSCALE = (1u << 13) | C1 | C2;
inline constexpr unsigned AREA = (1u << 14) | C4;
inline constexpr unsigned LONG_UNROLL = 1u << 15;
inline constexpr unsigned ALL = OUT_ALL | ALL_CAPS;
}  // namespace caps

// ---- the ellipsoid's geodesics ----------------------------------------------------------------------

class Geodesic {
 public:
  static constexpr int order = 6;  // of the series in ε and n
  static constexpr int nA1 = order, nC1 = order, nC1p = order, nA2 = order, nC2 = order, nA3 = order, nA3x = nA3;
  static constexpr int nC3 = order, nC3x = (nC3 * (nC3 - 1)) / 2, nC4 = order, nC4x = (nC4 * (nC4 + 1)) / 2;
  static constexpr int maxit1 = 20, maxit2 = maxit1 + DBL_MANT_DIG + 10;

  // An ellipsoid of equatorial radius a (m) and flattening f (negative for a prolate one).
  Geodesic(double a, double f) : a_(a), f_(f) {
    f1_ = 1 - f_;
    e2_ = f_ * (2 - f_);
    ep2_ = e2_ / math::sq(f1_);
    n_ = f_ / (2 - f_);
    b_ = a_ * f1_;
    // The authalic radius squared.
    c2_ = (math::sq(a_) + math::sq(b_) *
                              (e2_ == 0 ? 1.0
                                        : (e2_ > 0 ? std::atanh(std::sqrt(e2_)) : std::atan(std::sqrt(-e2_))) /
                                              std::sqrt(std::fabs(e2_)))) /
          2;
    // The threshold of arc below which a line counts as very short (GeographicLib's etol2).
    etol2_ = 0.1 * tol2() / std::sqrt(std::max(0.001, std::fabs(f_)) * std::min(1.0, 1 - f_ / 2) / 2);
    if (!(std::isfinite(a_) && a_ > 0)) throw std::invalid_argument("geodesy::geodesic: the equatorial radius is not positive");
    if (!(std::isfinite(b_) && b_ > 0)) throw std::invalid_argument("geodesy::geodesic: the polar semi-axis is not positive");
    A3coeff();
    C3coeff();
    C4coeff();
  }
  explicit Geodesic(const figure::Ellipsoid& e) : Geodesic(e.a(), e.flattening()) {}
  static Geodesic wgs84() { return Geodesic(figure::Ellipsoid::wgs84()); }

  double a() const { return a_; }
  double flattening() const { return f_; }
  // The area of the whole ellipsoid (m²): 4π c², c the authalic radius.
  double ellipsoid_area() const { return 4 * std::numbers::pi * c2_; }

  // The inverse problem between (lat1, lon1) and (lat2, lon2), radians.
  Inverse inverse(double lat1, double lon1, double lat2, double lon2) const {
    return inverse_degrees(to_lat_degrees(lat1), lon1 * math::rad_to_deg, to_lat_degrees(lat2), lon2 * math::rad_to_deg);
  }
  // The direct problem from (lat1, lon1) along the azimuth azi1 (radians) for s12 metres. With
  // `unroll`, the longitude is counted continuously along the line (it may leave [−π, π]).
  Direct direct(double lat1, double lon1, double azi1, double s12, bool unroll = false) const;
  // The direct problem with the arc length a12 (radians) on the auxiliary sphere given instead.
  Direct arc_direct(double lat1, double lon1, double azi1, double a12, bool unroll = false) const;

  // A geodesic from (lat1, lon1) along azi1 (radians), whose points can be had at any distance.
  GeodesicLine line(double lat1, double lon1, double azi1) const;
  // The geodesic between two points, as a line whose distance() is theirs.
  GeodesicLine inverse_line(double lat1, double lon1, double lat2, double lon2) const;

  // The same problems in degrees, as GeographicLib states them (inputs and outputs in degrees),
  // for exact reductions of angles given in degrees.
  Inverse inverse_degrees(double lat1, double lon1, double lat2, double lon2) const {
    const auto g = gen_inverse(lat1, lon1, lat2, lon2, caps::ALL);
    Inverse r;
    r.s12 = g.s12;
    r.a12 = g.a12 * math::deg_to_rad;
    r.azi1 = math::atan2d(g.salp1, g.calp1) * math::deg_to_rad;
    r.azi2 = math::atan2d(g.salp2, g.calp2) * math::deg_to_rad;
    r.m12 = g.m12;
    r.M12 = g.M12;
    r.M21 = g.M21;
    r.S12 = g.S12;
    return r;
  }
  // The direct problem in degrees: lat2, lon2, azi2 and a12 in degrees in the answer as well.
  Direct direct_degrees(double lat1, double lon1, double azi1, double s12_a12, bool arcmode = false, bool unroll = false) const;

 private:
  friend class GeodesicLine;

  static double tiny() { return std::sqrt(DBL_MIN); }
  static double tol0() { return DBL_EPSILON; }
  static double tol1() { return 200 * tol0(); }
  static double tol2() { return std::sqrt(tol0()); }
  static double tolb() { return tol0(); }
  static double xthresh() { return 1000 * tol2(); }

  // A latitude given in radians, in degrees; the rounding of ±π/2 is the pole.
  static double to_lat_degrees(double lat) {
    const double d = lat * math::rad_to_deg;
    return (std::fabs(d) > 90 && std::fabs(d) < 90 + 1e-9) ? std::copysign(90.0, d) : d;
  }

  // Clenshaw's sum of Σ c[l] sin 2lx (l = 1..n) when sinp, else Σ c[l] cos (2l + 1)x (l = 0..n−1);
  // len is the length of c (c[0] unused for the sine series).
  static double sin_cos_series(bool sinp, double sinx, double cosx, const double* c, int len) {
    int k = len;
    int n = k - (sinp ? 1 : 0);
    const double ar = 2 * (cosx - sinx) * (cosx + sinx);  // 2 cos 2x
    double y1 = 0, y0;
    if (n & 1) {
      --k;
      y0 = c[k];
    } else {
      y0 = 0;
    }
    n /= 2;
    while (n--) {
      --k;
      y1 = ar * y0 - y1 + c[k];
      --k;
      y0 = ar * y1 - y0 + c[k];
    }
    return sinp ? 2 * sinx * cosx * y0 : cosx * (y0 - y1);
  }

  // The positive root k of k⁴ + 2k³ − (x² + y² − 1)k² − 2y²k − y² = 0 (the astroid, eq. 55).
  static double astroid(double x, double y) {
    const double p = math::sq(x), q = math::sq(y), r = (p + q - 1) / 6;
    double k;
    if (!(q == 0 && r <= 0)) {
      const double S = p * q / 4, r2 = math::sq(r), r3 = r * r2;
      const double disc = S * (S + 2 * r3);
      double u = r;
      if (disc >= 0) {
        double T3 = S + r3;
        T3 += T3 < 0 ? -std::sqrt(disc) : std::sqrt(disc);
        const double T = std::cbrt(T3);
        u += T + (T != 0 ? r2 / T : 0);
      } else {
        const double ang = std::atan2(std::sqrt(-disc), -(S + r3));
        u += 2 * r * std::cos(ang / 3);
      }
      const double v = std::sqrt(math::sq(u) + q);
      const double uv = u < 0 ? q / (v - u) : u + v;
      const double w = (uv - q) / (2 * v);
      k = uv / (std::sqrt(uv + math::sq(w)) + w);
    } else {
      k = 0;
    }
    return k;
  }

  // A₁ − 1 (eq. 17).
  static double A1m1f(double eps) {
    static const double coeff[] = {1, 4, 64, 0, 256};
    const int m = nA1 / 2;
    const double t = math::polyval(m, coeff, 0, math::sq(eps)) / coeff[m + 1];
    return (t + eps) / (1 - eps);
  }
  // C₁ₗ (eq. 18).
  static void C1f(double eps, double* c) {
    static const double coeff[] = {-1, 6, -16, 32, -9, 64, -128, 2048, 9, -16, 768, 3, -5, 512, -7, 1280, -7, 2048};
    const double eps2 = math::sq(eps);
    double d = eps;
    int o = 0;
    for (int l = 1; l <= nC1; ++l) {
      const int m = (nC1 - l) / 2;
      c[l] = d * math::polyval(m, coeff, o, eps2) / coeff[o + m + 1];
      o += m + 2;
      d *= eps;
    }
  }
  // C′₁ₗ, the coefficients of the reverted series (eq. 21).
  static void C1pf(double eps, double* c) {
    static const double coeff[] = {205, -432, 768, 1536, 4005, -4736, 3840, 12288, -225, 116, 384,
                                   -7173, 2695, 7680, 3467, 7680, 38081, 61440};
    const double eps2 = math::sq(eps);
    double d = eps;
    int o = 0;
    for (int l = 1; l <= nC1p; ++l) {
      const int m = (nC1p - l) / 2;
      c[l] = d * math::polyval(m, coeff, o, eps2) / coeff[o + m + 1];
      o += m + 2;
      d *= eps;
    }
  }
  // A₂ − 1 (eq. 42).
  static double A2m1f(double eps) {
    static const double coeff[] = {-11, -28, -192, 0, 256};
    const int m = nA2 / 2;
    const double t = math::polyval(m, coeff, 0, math::sq(eps)) / coeff[m + 1];
    return (t - eps) / (1 + eps);
  }
  // C₂ₗ (eq. 43).
  static void C2f(double eps, double* c) {
    static const double coeff[] = {1, 2, 16, 32, 35, 64, 384, 2048, 15, 80, 768, 7, 35, 512, 63, 1280, 77, 2048};
    const double eps2 = math::sq(eps);
    double d = eps;
    int o = 0;
    for (int l = 1; l <= nC2; ++l) {
      const int m = (nC2 - l) / 2;
      c[l] = d * math::polyval(m, coeff, o, eps2) / coeff[o + m + 1];
      o += m + 2;
      d *= eps;
    }
  }

  // The coefficients of A₃ in ε, as polynomials in n (eq. 24).
  void A3coeff() {
    static const double coeff[] = {-3, 128, -2, -3, 64, -1, -3, -1, 16, 3, -1, -2, 8, 1, -1, 2, 1, 1};
    int o = 0, k = 0;
    for (int j = nA3 - 1; j >= 0; --j) {
      const int m = std::min(nA3 - j - 1, j);
      A3x_[k++] = math::polyval(m, coeff, o, n_) / coeff[o + m + 1];
      o += m + 2;
    }
  }
  // The coefficients of C₃ₗ (eq. 25).
  void C3coeff() {
    static const double coeff[] = {3, 128, 2, 5, 128, -1, 3, 3, 64, -1, 0, 1, 8, -1, 1, 4, 5, 256, 1, 3, 128, -3, -2, 3, 64,
                                   1, -3, 2, 32, 7, 512, -10, 9, 384, 5, -9, 5, 192, 7, 512, -14, 7, 512, 21, 2560};
    int o = 0, k = 0;
    for (int l = 1; l < nC3; ++l) {
      for (int j = nC3 - 1; j >= l; --j) {
        const int m = std::min(nC3 - j - 1, j);
        C3x_[k++] = math::polyval(m, coeff, o, n_) / coeff[o + m + 1];
        o += m + 2;
      }
    }
  }
  // The coefficients of C₄ₗ, for the area (eq. 59).
  void C4coeff() {
    static const double coeff[] = {
        97,    15015,  1088,  156,   45045,  -224,  -4784, 1573,   45045, -10656, 14144, -4576, -858,  45045,
        64,    624,    -4576, 6864,  -3003,  15015, 100,   208,    572,   3432,   -12012, 30030, 45045, 1,
        9009,  -2944,  468,   135135, 5792,  1040,  -1287, 135135, 5952,  -11648, 9152,  -2574, 135135, -64,
        -624,  4576,   -6864, 3003,  135135, 8,     10725, 1856,   -936,  225225, -8448, 4992,  -1144, 225225,
        -1440, 4160,   -4576, 1716,  225225, -136,  63063, 1024,   -208,  105105, 3584,  -3328, 1144,  315315,
        -128,  135135, -2560, 832,   405405, 128,   99099};
    int o = 0, k = 0;
    for (int l = 0; l < nC4; ++l) {
      for (int j = nC4 - 1; j >= l; --j) {
        const int m = nC4 - j - 1;
        C4x_[k++] = math::polyval(m, coeff, o, n_) / coeff[o + m + 1];
        o += m + 2;
      }
    }
  }

  double A3f(double eps) const { return math::polyval(nA3 - 1, A3x_, 0, eps); }
  void C3f(double eps, double* c) const {
    double mult = 1;
    int o = 0;
    for (int l = 1; l < nC3; ++l) {
      const int m = nC3 - l - 1;
      mult *= eps;
      c[l] = mult * math::polyval(m, C3x_, o, eps);
      o += m + 1;
    }
  }
  void C4f(double eps, double* c) const {
    double mult = 1;
    int o = 0;
    for (int l = 0; l < nC4; ++l) {
      const int m = nC4 - l - 1;
      c[l] = mult * math::polyval(m, C4x_, o, eps);
      o += m + 1;
      mult *= eps;
    }
  }

  struct Lengths {
    double s12b = math::nan, m12b = math::nan, m0 = math::nan, M12 = math::nan, M21 = math::nan;
  };
  // The distance and reduced length over b, the secular coefficient m₀ of the reduced length, and the
  // geodesic scales, along the arc σ₁₂ (eqs. 38–42).
  Lengths lengths(double eps, double sig12, double ssig1, double csig1, double dn1, double ssig2, double csig2, double dn2,
                  double cbet1, double cbet2, unsigned outmask, double* C1a, double* C2a) const {
    outmask &= caps::OUT_MASK;
    Lengths r;
    double m0x = math::nan, J12 = math::nan, A1 = 0, A2 = 0;
    if (outmask & (caps::DISTANCE | caps::REDUCEDLENGTH | caps::GEODESICSCALE)) {
      A1 = A1m1f(eps);
      C1f(eps, C1a);
      if (outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE)) {
        A2 = A2m1f(eps);
        C2f(eps, C2a);
        m0x = A1 - A2;
        A2 = 1 + A2;
      }
      A1 = 1 + A1;
    }
    if (outmask & caps::DISTANCE) {
      const double B1 = sin_cos_series(true, ssig2, csig2, C1a, nC1 + 1) - sin_cos_series(true, ssig1, csig1, C1a, nC1 + 1);
      r.s12b = A1 * (sig12 + B1);
      if (outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE)) {
        const double B2 = sin_cos_series(true, ssig2, csig2, C2a, nC2 + 1) - sin_cos_series(true, ssig1, csig1, C2a, nC2 + 1);
        J12 = m0x * sig12 + (A1 * B1 - A2 * B2);
      }
    } else if (outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE)) {
      for (int l = 1; l <= nC2; ++l) C2a[l] = A1 * C1a[l] - A2 * C2a[l];
      J12 = m0x * sig12 + (sin_cos_series(true, ssig2, csig2, C2a, nC2 + 1) - sin_cos_series(true, ssig1, csig1, C2a, nC2 + 1));
    }
    if (outmask & caps::REDUCEDLENGTH) {
      r.m0 = m0x;
      // The products in parentheses cancel accurately for coincident points.
      r.m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
    }
    if (outmask & caps::GEODESICSCALE) {
      const double csig12 = csig1 * csig2 + ssig1 * ssig2;
      const double t = ep2_ * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
      r.M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
      r.M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
    }
    return r;
  }

  struct Start {
    double sig12, salp1, calp1, salp2, calp2, dnm;
  };
  // A starting azimuth for Newton's method; for very short lines the solution itself (sig12 ≥ 0).
  Start inverse_start(double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double lam12,
                      double slam12, double clam12, double* C1a, double* C2a) const {
    Start r{-1, 0, 0, math::nan, math::nan, math::nan};
    const double sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
    const double cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
    volatile double sbet12a_v = sbet2 * cbet1;
    sbet12a_v = sbet12a_v + cbet2 * sbet1;
    const double sbet12a = sbet12a_v;
    const bool shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
    double somg12, comg12;
    if (shortline) {
      double sbetm2 = math::sq(sbet1 + sbet2);
      sbetm2 /= sbetm2 + math::sq(cbet1 + cbet2);
      r.dnm = std::sqrt(1 + ep2_ * sbetm2);
      const double omg12 = lam12 / (f1_ * r.dnm);
      somg12 = std::sin(omg12);
      comg12 = std::cos(omg12);
    } else {
      somg12 = slam12;
      comg12 = clam12;
    }
    double salp1 = cbet2 * somg12;
    double calp1 = comg12 >= 0 ? sbet12 + cbet2 * sbet1 * math::sq(somg12) / (1 + comg12)
                               : sbet12a - cbet2 * sbet1 * math::sq(somg12) / (1 - comg12);
    const double ssig12 = std::hypot(salp1, calp1);
    const double csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;
    if (shortline && ssig12 < etol2_) {
      // Very short lines.
      double salp2 = cbet1 * somg12;
      double calp2 = sbet12 - cbet1 * sbet2 * (comg12 >= 0 ? math::sq(somg12) / (1 + comg12) : 1 - comg12);
      std::tie(salp2, calp2) = math::norm(salp2, calp2);
      r.salp2 = salp2;
      r.calp2 = calp2;
      r.sig12 = std::atan2(ssig12, csig12);
    } else if (std::fabs(n_) >= 0.1 || csig12 >= 0 || ssig12 >= 6 * std::fabs(n_) * std::numbers::pi * math::sq(cbet1)) {
      // The zeroth-order spherical estimate will do.
    } else {
      // Nearly antipodal: scale λ₁₂ and β₂ so that the antipodal point is at the origin and the
      // singular point at (−1, 0), and estimate α₁ from the astroid (§6).
      const double lam12x = std::atan2(-slam12, -clam12);
      double x, y, lamscale, betscale;
      if (f_ >= 0) {
        const double k2 = math::sq(sbet1) * ep2_;
        const double eps = k2 / (2 * (1 + std::sqrt(1 + k2)) + k2);
        lamscale = f_ * cbet1 * A3f(eps) * std::numbers::pi;
        betscale = lamscale * cbet1;
        x = lam12x / lamscale;
        y = sbet12a / betscale;
      } else {
        const double cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
        const double bet12a = std::atan2(sbet12a, cbet12a);
        const Lengths l = lengths(n_, std::numbers::pi + bet12a, sbet1, -cbet1, dn1, sbet2, cbet2, dn2, cbet1, cbet2,
                                  caps::REDUCEDLENGTH, C1a, C2a);
        x = -1 + l.m12b / (cbet1 * cbet2 * l.m0 * std::numbers::pi);
        betscale = x < -0.01 ? sbet12a / x : -f_ * math::sq(cbet1) * std::numbers::pi;
        lamscale = betscale / cbet1;
        y = lam12x / lamscale;
      }
      if (y > -tol1() && x > -1 - xthresh()) {
        // Near the cut.
        if (f_ >= 0) {
          salp1 = std::min(1.0, -x);
          calp1 = -std::sqrt(1 - math::sq(salp1));
        } else {
          calp1 = std::max(x > -tol1() ? 0.0 : -1.0, x);
          salp1 = std::sqrt(1 - math::sq(calp1));
        }
      } else {
        // The astroid gives ω₁₂ (better than α₁ directly), and the sphere α₁ from it.
        const double k = astroid(x, y);
        const double omg12a = lamscale * (f_ >= 0 ? -x * k / (1 + k) : -y * (1 + k) / k);
        somg12 = std::sin(omg12a);
        comg12 = -std::cos(omg12a);
        salp1 = cbet2 * somg12;
        calp1 = sbet12a - cbet2 * sbet1 * math::sq(somg12) / (1 - comg12);
      }
    }
    // A sanity check on the start; the backward test lets a NaN through.
    if (!(salp1 <= 0)) {
      std::tie(salp1, calp1) = math::norm(salp1, calp1);
    } else {
      salp1 = 1;
      calp1 = 0;
    }
    r.salp1 = salp1;
    r.calp1 = calp1;
    return r;
  }

  struct Lambda {
    double lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12;
  };
  // The longitude difference λ₁₂ of the geodesic leaving the first point along α₁ when it reaches
  // the latitude of the second (the hybrid problem), less λ₁₂ given; with its derivative (eq. 46).
  Lambda lambda12(double sbet1, double cbet1, double dn1, double sbet2, double cbet2, double dn2, double salp1, double calp1,
                  double slam120, double clam120, bool diffp, double* C1a, double* C2a, double* C3a) const {
    Lambda r{};
    if (sbet1 == 0 && calp1 == 0) calp1 = -tiny();  // the degenerate equatorial line, handled already
    const double salp0 = salp1 * cbet1;
    const double calp0 = std::hypot(calp1, salp1 * sbet1);
    double ssig1 = sbet1;
    const double somg1 = salp0 * sbet1;
    double csig1 = calp1 * cbet1;
    const double comg1 = calp1 * cbet1;
    std::tie(ssig1, csig1) = math::norm(ssig1, csig1);
    r.salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
    r.calp2 = (cbet2 != cbet1 || std::fabs(sbet2) != -sbet1)
                  ? std::sqrt(math::sq(calp1 * cbet1) + (cbet1 < -sbet1 ? (cbet2 - cbet1) * (cbet1 + cbet2)
                                                                       : (sbet1 - sbet2) * (sbet1 + sbet2))) /
                        cbet2
                  : std::fabs(calp1);
    double ssig2 = sbet2;
    const double somg2 = salp0 * sbet2;
    double csig2 = r.calp2 * cbet2;
    const double comg2 = r.calp2 * cbet2;
    std::tie(ssig2, csig2) = math::norm(ssig2, csig2);
    r.sig12 = std::atan2(std::max(0.0, csig1 * ssig2 - ssig1 * csig2) + 0.0, csig1 * csig2 + ssig1 * ssig2);
    const double somg12 = std::max(0.0, comg1 * somg2 - somg1 * comg2) + 0.0;
    const double comg12 = comg1 * comg2 + somg1 * somg2;
    const double eta = std::atan2(somg12 * clam120 - comg12 * slam120, comg12 * clam120 + somg12 * slam120);
    const double k2 = math::sq(calp0) * ep2_;
    r.eps = k2 / (2 * (1 + std::sqrt(1 + k2)) + k2);
    C3f(r.eps, C3a);
    const double B312 = sin_cos_series(true, ssig2, csig2, C3a, nC3) - sin_cos_series(true, ssig1, csig1, C3a, nC3);
    r.domg12 = -f_ * A3f(r.eps) * salp0 * (r.sig12 + B312);
    r.lam12 = eta + r.domg12;
    if (diffp) {
      if (r.calp2 == 0) {
        r.dlam12 = -2 * f1_ * dn1 / sbet1;
      } else {
        const Lengths l = lengths(r.eps, r.sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, caps::REDUCEDLENGTH, C1a, C2a);
        r.dlam12 = l.m12b * f1_ / (r.calp2 * cbet2);
      }
    } else {
      r.dlam12 = math::nan;
    }
    r.ssig1 = ssig1;
    r.csig1 = csig1;
    r.ssig2 = ssig2;
    r.csig2 = csig2;
    return r;
  }

  struct General {
    double a12 = math::nan, s12 = math::nan, salp1 = math::nan, calp1 = math::nan, salp2 = math::nan, calp2 = math::nan;
    double m12 = math::nan, M12 = math::nan, M21 = math::nan, S12 = math::nan;
  };
  // The inverse problem in degrees (GeographicLib's GenInverse): a12 in degrees, the azimuths as
  // their sines and cosines.
  General gen_inverse(double lat1, double lon1, double lat2, double lon2, unsigned outmask) const {
    General g;
    outmask &= caps::OUT_MASK;
    // The longitude difference, made positive; −180 only for west-going geodesics.
    auto [lon12, lon12s] = math::ang_diff(lon1, lon2);
    double lonsign = std::copysign(1.0, lon12);
    lon12 = lonsign * lon12;
    lon12s = lonsign * lon12s;
    const double lam12 = lon12 * math::deg_to_rad;
    auto [slam12, clam12] = math::sincosde(lon12, lon12s);
    lon12s = (180 - lon12) - lon12s;  // the supplementary longitude difference

    // Near the equator, on the equator.
    lat1 = math::ang_round(math::lat_fix(lat1));
    lat2 = math::ang_round(math::lat_fix(lat2));
    // The point with the larger |latitude| first; a NaN latitude becomes the first.
    const double swapp = (std::fabs(lat1) < std::fabs(lat2) || std::isnan(lat2)) ? -1 : 1;
    if (swapp < 0) {
      lonsign *= -1;
      std::swap(lat1, lat2);
    }
    // lat1 ≤ 0.
    const double latsign = std::copysign(1.0, -lat1);
    lat1 *= latsign;
    lat2 *= latsign;
    // Now 0 ≤ lon12 ≤ 180, −90 ≤ lat1 ≤ 0, lat1 ≤ lat2 ≤ −lat1.

    auto [sbet1, cbet1] = math::sincosd(lat1);
    sbet1 *= f1_;
    std::tie(sbet1, cbet1) = math::norm(sbet1, cbet1);
    cbet1 = std::max(tiny(), cbet1);  // +ε at the poles
    auto [sbet2, cbet2] = math::sincosd(lat2);
    sbet2 *= f1_;
    std::tie(sbet2, cbet2) = math::norm(sbet2, cbet2);
    cbet2 = std::max(tiny(), cbet2);
    // Force β₂ = ±β₁ exactly where the measures of |β₁| − |β₂| vanish.
    if (cbet1 < -sbet1) {
      if (cbet2 == cbet1) sbet2 = std::copysign(sbet1, sbet2);
    } else {
      if (std::fabs(sbet2) == -sbet1) cbet2 = cbet1;
    }
    const double dn1 = std::sqrt(1 + ep2_ * math::sq(sbet1));
    const double dn2 = std::sqrt(1 + ep2_ * math::sq(sbet2));

    double C1a[nC1 + 1] = {0}, C2a[nC2 + 1] = {0}, C3a[nC3] = {0};
    double a12 = math::nan, sig12 = 0, calp1 = 0, salp1 = 0, calp2 = 0, salp2 = 0;
    double s12x = math::nan, m12x = math::nan, M12 = math::nan, M21 = math::nan;
    double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0, domg12 = 0;

    bool meridian = lat1 == -90 || slam12 == 0;
    if (meridian) {
      // Both points on one full meridian: the geodesic may lie along it.
      calp1 = clam12;
      salp1 = slam12;  // toward the target's longitude
      calp2 = 1.0;
      salp2 = 0.0;  // heading north at the target
      ssig1 = sbet1;
      csig1 = calp1 * cbet1;
      ssig2 = sbet2;
      csig2 = calp2 * cbet2;
      sig12 = std::atan2(std::max(0.0, csig1 * ssig2 - ssig1 * csig2) + 0.0, csig1 * csig2 + ssig1 * ssig2);
      const Lengths l = lengths(n_, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
                                outmask | caps::DISTANCE | caps::REDUCEDLENGTH, C1a, C2a);
      s12x = l.s12b;
      m12x = l.m12b;
      M12 = l.M12;
      M21 = l.M21;
      // The meridian is the shortest path unless m12 < 0 (a prolate ellipsoid, nearly antipodal).
      if (sig12 < tol2() || m12x >= 0) {
        if (sig12 < 3 * tiny() || (sig12 < tol0() && (s12x < 0 || m12x < 0))) sig12 = m12x = s12x = 0.0;
        m12x *= b_;
        s12x *= b_;
        a12 = sig12 * math::rad_to_deg;
      } else {
        meridian = false;
      }
    }

    double somg12 = 2.0, comg12 = 0.0, omg12 = 0.0;  // somg12 = 2 marks "to be computed"
    if (!meridian && sbet1 == 0 && (f_ <= 0 || lon12s >= f_ * 180)) {
      // Along the equator.
      calp1 = calp2 = 0.0;
      salp1 = salp2 = 1.0;
      s12x = a_ * lam12;
      sig12 = omg12 = lam12 / f1_;
      m12x = b_ * std::sin(sig12);
      if (outmask & caps::GEODESICSCALE) M12 = M21 = std::cos(sig12);
      a12 = lon12 / f1_;
    } else if (!meridian) {
      // Neither meridional nor equatorial, the points within a hemisphere bounded by a meridian.
      const Start st = inverse_start(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, C1a, C2a);
      sig12 = st.sig12;
      salp1 = st.salp1;
      calp1 = st.calp1;
      if (sig12 >= 0) {
        // Short lines: the start is the solution.
        salp2 = st.salp2;
        calp2 = st.calp2;
        s12x = sig12 * b_ * st.dnm;
        m12x = math::sq(st.dnm) * b_ * std::sin(sig12 / st.dnm);
        if (outmask & caps::GEODESICSCALE) M12 = M21 = std::cos(sig12 / st.dnm);
        a12 = sig12 * math::rad_to_deg;
        omg12 = lam12 / (f1_ * st.dnm);
      } else {
        // Newton's method on f(α₁) = λ₁₂(α₁) − λ₁₂, which has one root in (0, π) with a positive
        // slope there; a bracket (α₁a, α₁b) of the root shrinks with each evaluation, and the
        // iteration restarts from its middle when the slope is not positive or the step leaves (0, π).
        int numit = 0;
        bool tripn = false, tripb = false;
        double salp1a = tiny(), calp1a = 1.0, salp1b = tiny(), calp1b = -1.0;
        for (;;) {
          const Lambda lm = lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1, slam12, clam12, numit < maxit1, C1a, C2a, C3a);
          const double v = lm.lam12, dv = lm.dlam12;
          salp2 = lm.salp2;
          calp2 = lm.calp2;
          sig12 = lm.sig12;
          ssig1 = lm.ssig1;
          csig1 = lm.csig1;
          ssig2 = lm.ssig2;
          csig2 = lm.csig2;
          eps = lm.eps;
          domg12 = lm.domg12;
          // The reversed test lets NaNs escape.
          if (tripb || !(std::fabs(v) >= (tripn ? 8 : 1) * tol0()) || numit == maxit2) break;
          // Update the bracket.
          if (v > 0 && (numit > maxit1 || calp1 / salp1 > calp1b / salp1b)) {
            salp1b = salp1;
            calp1b = calp1;
          } else if (v < 0 && (numit > maxit1 || calp1 / salp1 < calp1a / salp1a)) {
            salp1a = salp1;
            calp1a = calp1;
          }
          ++numit;
          if (numit < maxit1 && dv > 0) {
            const double dalp1 = -v / dv;
            if (std::fabs(dalp1) < std::numbers::pi) {
              const double sdalp1 = std::sin(dalp1), cdalp1 = std::cos(dalp1);
              const double nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
              if (nsalp1 > 0) {
                calp1 = calp1 * cdalp1 - salp1 * sdalp1;
                salp1 = nsalp1;
                std::tie(salp1, calp1) = math::norm(salp1, calp1);
                // Where the slope tends to zero, convergence is not quadratic: test against ε.
                tripn = std::fabs(v) <= 16 * tol0();
                continue;
              }
            }
          }
          // Bisection of the bracket.
          salp1 = (salp1a + salp1b) / 2;
          calp1 = (calp1a + calp1b) / 2;
          std::tie(salp1, calp1) = math::norm(salp1, calp1);
          tripn = false;
          tripb = (std::fabs(salp1a - salp1) + (calp1a - calp1) < tolb() || std::fabs(salp1 - salp1b) + (calp1 - calp1b) < tolb());
        }
        const unsigned lengthmask = outmask | ((outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE)) ? caps::DISTANCE : caps::EMPTY);
        const Lengths l = lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2, lengthmask, C1a, C2a);
        s12x = l.s12b * b_;
        m12x = l.m12b * b_;
        M12 = l.M12;
        M21 = l.M21;
        a12 = sig12 * math::rad_to_deg;
        if (outmask & caps::AREA) {
          // ω₁₂ = λ₁₂ − Δω₁₂
          const double sdomg12 = std::sin(domg12), cdomg12 = std::cos(domg12);
          somg12 = slam12 * cdomg12 - clam12 * sdomg12;
          comg12 = clam12 * cdomg12 + slam12 * sdomg12;
        }
      }
    }

    if (outmask & caps::DISTANCE) g.s12 = 0.0 + s12x;
    if (outmask & caps::REDUCEDLENGTH) g.m12 = 0.0 + m12x;
    if (outmask & caps::AREA) {
      double S12;
      const double salp0 = salp1 * cbet1, calp0 = std::hypot(calp1, salp1 * sbet1);
      if (calp0 != 0 && salp0 != 0) {
        double ss1 = sbet1, cs1 = calp1 * cbet1, ss2 = sbet2, cs2 = calp2 * cbet2;
        const double k2 = math::sq(calp0) * ep2_;
        const double e = k2 / (2 * (1 + std::sqrt(1 + k2)) + k2);
        const double A4 = math::sq(a_) * calp0 * salp0 * e2_;  // a² e² cos α₀ sin α₀
        std::tie(ss1, cs1) = math::norm(ss1, cs1);
        std::tie(ss2, cs2) = math::norm(ss2, cs2);
        double C4a[nC4];
        C4f(e, C4a);
        const double B41 = sin_cos_series(false, ss1, cs1, C4a, nC4), B42 = sin_cos_series(false, ss2, cs2, C4a, nC4);
        S12 = A4 * (B42 - B41);
      } else {
        S12 = 0.0;  // σ₁, σ₂ indeterminate on the equator
      }
      if (!meridian && somg12 == 2.0) {
        somg12 = std::sin(omg12);
        comg12 = std::cos(omg12);
      }
      double alp12;
      if (!meridian && comg12 > -0.7071 && sbet2 - sbet1 < 1.75) {
        // tan(Γ/2) from tan(ω₁₂/2) and the half-angles of β₁ and β₂.
        const double domg = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
        alp12 = 2 * std::atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1), domg * (sbet1 * sbet2 + dbet1 * dbet2));
      } else {
        double salp12 = salp2 * calp1 - calp2 * salp1, calp12 = calp2 * calp1 + salp2 * salp1;
        if (salp12 == 0 && calp12 < 0) {
          salp12 = tiny() * calp1;
          calp12 = -1.0;
        }
        alp12 = std::atan2(salp12, calp12);
      }
      S12 += c2_ * alp12;
      S12 *= swapp * lonsign * latsign;
      g.S12 = S12 + 0.0;
    }
    // Undo the transformations to the canonical form.
    if (swapp < 0) {
      std::swap(salp1, salp2);
      std::swap(calp1, calp2);
      if (outmask & caps::GEODESICSCALE) std::swap(M12, M21);
    }
    salp1 *= swapp * lonsign;
    calp1 *= swapp * latsign;
    salp2 *= swapp * lonsign;
    calp2 *= swapp * latsign;
    g.a12 = a12;
    g.salp1 = salp1;
    g.calp1 = calp1;
    g.salp2 = salp2;
    g.calp2 = calp2;
    g.M12 = M12;
    g.M21 = M21;
    return g;
  }

  double a_, f_, f1_, e2_, ep2_, n_, b_, c2_, etol2_;
  double A3x_[nA3x], C3x_[nC3x], C4x_[nC4x];
};

// ---- points along a geodesic ----------------------------------------------------------------------

// A geodesic from a point along an azimuth: its auxiliary-sphere quantities computed once, so that
// each point along it costs one evaluation of the series (Karney 2013, §3).
class GeodesicLine {
 public:
  // In degrees, with the azimuth's sine and cosine if known (as GeographicLib's constructor).
  GeodesicLine(const Geodesic& g, double lat1_deg, double lon1_deg, double azi1_deg,
               unsigned capabilities = caps::ALL | caps::DISTANCE_IN, double salp1 = math::nan, double calp1 = math::nan)
      : a_(g.a_), f_(g.f_), b_(g.b_), c2_(g.c2_), f1_(g.f1_) {
    caps_ = capabilities | caps::LATITUDE | caps::AZIMUTH | caps::LONG_UNROLL;
    lat1_ = math::lat_fix(lat1_deg);
    lon1_ = lon1_deg;
    if (std::isnan(salp1) || std::isnan(calp1)) {
      azi1_ = math::ang_normalize(azi1_deg);
      std::tie(salp1_, calp1_) = math::sincosd(math::ang_round(azi1_deg));
    } else {
      azi1_ = azi1_deg;
      salp1_ = salp1;
      calp1_ = calp1;
    }
    auto [sbet1, cbet1] = math::sincosd(math::ang_round(lat1_));
    sbet1 *= f1_;
    std::tie(sbet1, cbet1) = math::norm(sbet1, cbet1);
    cbet1 = std::max(Geodesic::tiny(), cbet1);
    dn1_ = std::sqrt(1 + g.ep2_ * math::sq(sbet1));
    // sin α₀ = sin α₁ cos β₁ (Clairaut).
    salp0_ = salp1_ * cbet1;
    calp0_ = std::hypot(calp1_, salp1_ * sbet1);
    // tan β₁ = tan σ₁ cos α₁, tan ω₁ = sin α₀ tan σ₁; σ = 0 at the northward crossing of the equator.
    ssig1_ = sbet1;
    somg1_ = salp0_ * sbet1;
    csig1_ = comg1_ = (sbet1 != 0 || calp1_ != 0) ? cbet1 * calp1_ : 1;
    std::tie(ssig1_, csig1_) = math::norm(ssig1_, csig1_);
    k2_ = math::sq(calp0_) * g.ep2_;
    const double eps = k2_ / (2 * (1 + std::sqrt(1 + k2_)) + k2_);
    if (caps_ & caps::C1) {
      A1m1_ = Geodesic::A1m1f(eps);
      Geodesic::C1f(eps, C1a_);
      B11_ = Geodesic::sin_cos_series(true, ssig1_, csig1_, C1a_, Geodesic::nC1 + 1);
      const double s = std::sin(B11_), c = std::cos(B11_);
      stau1_ = ssig1_ * c + csig1_ * s;  // τ₁ = σ₁ + B₁₁
      ctau1_ = csig1_ * c - ssig1_ * s;
    }
    if (caps_ & caps::C1p) Geodesic::C1pf(eps, C1pa_);
    if (caps_ & caps::C2) {
      A2m1_ = Geodesic::A2m1f(eps);
      Geodesic::C2f(eps, C2a_);
      B21_ = Geodesic::sin_cos_series(true, ssig1_, csig1_, C2a_, Geodesic::nC2 + 1);
    }
    if (caps_ & caps::C3) {
      g.C3f(eps, C3a_);
      A3c_ = -f_ * salp0_ * g.A3f(eps);
      B31_ = Geodesic::sin_cos_series(true, ssig1_, csig1_, C3a_, Geodesic::nC3);
    }
    if (caps_ & caps::C4) {
      g.C4f(eps, C4a_);
      A4_ = math::sq(a_) * calp0_ * salp0_ * g.e2_;
      B41_ = Geodesic::sin_cos_series(false, ssig1_, csig1_, C4a_, Geodesic::nC4);
    }
  }

  // The point s12 metres along the line (radians). With `unroll`, the longitude is counted
  // continuously from the start.
  Direct position(double s12, bool unroll = false) const { return to_radians(gen_position(false, s12, caps::ALL | (unroll ? caps::LONG_UNROLL : 0))); }
  // The point a12 radians of the auxiliary sphere along the line.
  Direct arc_position(double a12, bool unroll = false) const {
    return to_radians(gen_position(true, a12 * math::rad_to_deg, caps::ALL | (unroll ? caps::LONG_UNROLL : 0)));
  }
  // The same in degrees (lat2, lon2, azi2, a12 of the answer in degrees).
  Direct position_degrees(double s12_a12, bool arcmode = false, bool unroll = false) const {
    return gen_position(arcmode, s12_a12, caps::ALL | (unroll ? caps::LONG_UNROLL : 0));
  }

  // The start (radians).
  double latitude() const { return lat1_ * math::deg_to_rad; }
  double longitude() const { return lon1_ * math::deg_to_rad; }
  double azimuth() const { return azi1_ * math::deg_to_rad; }
  // The distance (m) and arc (rad) to the line's end, when it was made between two points or to a
  // distance; NaN for a line without an end.
  double distance() const { return s13_; }
  double arc() const { return a13_ * math::deg_to_rad; }

  // n + 1 points evenly spaced from the start to the end of the line (radians), for drawing.
  std::vector<Direct> points(int n, bool unroll = false) const {
    if (n < 1) throw std::invalid_argument("geodesy::geodesic: a path needs at least one step");
    std::vector<Direct> out;
    out.reserve(n + 1);
    for (int i = 0; i <= n; ++i) out.push_back(position(s13_ * i / n, unroll));
    return out;
  }

  void set_distance(double s13) {
    s13_ = s13;
    a13_ = gen_position(false, s13_, 0).a12;
  }
  void set_arc_degrees(double a13) {
    a13_ = a13;
    s13_ = gen_position(true, a13_, caps::DISTANCE).s12;
  }

  // GeographicLib's GenPosition, in degrees.
  Direct gen_position(bool arcmode, double s12_a12, unsigned outmask) const {
    Direct r;
    outmask &= caps_ & caps::OUT_MASK;
    if (!(arcmode || (caps_ & (caps::OUT_MASK & caps::DISTANCE_IN)))) return r;  // no distance possible
    double B12 = 0.0, AB1 = 0.0, sig12, ssig12, csig12;
    if (arcmode) {
      sig12 = s12_a12 * math::deg_to_rad;
      std::tie(ssig12, csig12) = math::sincosd(s12_a12);
    } else {
      // σ₁₂ from the distance, by the reverted series (eq. 20–21).
      double tau12 = s12_a12 / (b_ * (1 + A1m1_));
      tau12 = std::isfinite(tau12) ? tau12 : math::nan;
      const double s = std::sin(tau12), c = std::cos(tau12);
      B12 = -Geodesic::sin_cos_series(true, stau1_ * c + ctau1_ * s, ctau1_ * c - stau1_ * s, C1pa_, Geodesic::nC1p + 1);
      sig12 = tau12 - (B12 - B11_);
      ssig12 = std::sin(sig12);
      csig12 = std::cos(sig12);
      if (std::fabs(f_) > 0.01) {
        // The reverted series is not accurate enough beyond |f| = 1/100: one Newton step.
        const double ssig2 = ssig1_ * csig12 + csig1_ * ssig12, csig2 = csig1_ * csig12 - ssig1_ * ssig12;
        B12 = Geodesic::sin_cos_series(true, ssig2, csig2, C1a_, Geodesic::nC1 + 1);
        const double serr = (1 + A1m1_) * (sig12 + (B12 - B11_)) - s12_a12 / b_;
        sig12 = sig12 - serr / std::sqrt(1 + k2_ * math::sq(ssig2));
        ssig12 = std::sin(sig12);
        csig12 = std::cos(sig12);
      }
    }
    // σ₂ = σ₁ + σ₁₂.
    const double ssig2 = ssig1_ * csig12 + csig1_ * ssig12;
    double csig2 = csig1_ * csig12 - ssig1_ * ssig12;
    const double dn2 = std::sqrt(1 + k2_ * math::sq(ssig2));
    if (outmask & (caps::DISTANCE | caps::REDUCEDLENGTH | caps::GEODESICSCALE)) {
      if (arcmode || std::fabs(f_) > 0.01) B12 = Geodesic::sin_cos_series(true, ssig2, csig2, C1a_, Geodesic::nC1 + 1);
      AB1 = (1 + A1m1_) * (B12 - B11_);
    }
    // sin β₂ = cos α₀ sin σ₂.
    const double sbet2 = calp0_ * ssig2;
    double cbet2 = std::hypot(salp0_, calp0_ * csig2);
    if (cbet2 == 0) cbet2 = csig2 = Geodesic::tiny();  // sin α₀ = 0 and cos σ₂ = 0: the degeneracy broken
    const double salp2 = salp0_, calp2 = calp0_ * csig2;  // tan α₀ = cos σ₂ tan α₂
    if (outmask & caps::DISTANCE) r.s12 = arcmode ? b_ * ((1 + A1m1_) * sig12 + AB1) : s12_a12;
    if (outmask & caps::LONGITUDE) {
      const double somg2 = salp0_ * ssig2, comg2 = csig2;  // tan ω₂ = sin α₀ tan σ₂
      const double E = std::copysign(1.0, salp0_);         // east- or west-going
      const double omg12 = (outmask & caps::LONG_UNROLL)
                               ? E * (sig12 - (std::atan2(ssig2, csig2) - std::atan2(ssig1_, csig1_)) +
                                      (std::atan2(E * somg2, comg2) - std::atan2(E * somg1_, comg1_)))
                               : std::atan2(somg2 * comg1_ - comg2 * somg1_, comg2 * comg1_ + somg2 * somg1_);
      const double lam12 = omg12 + A3c_ * (sig12 + (Geodesic::sin_cos_series(true, ssig2, csig2, C3a_, Geodesic::nC3) - B31_));
      const double lon12 = lam12 * math::rad_to_deg;
      r.lon2 = (outmask & caps::LONG_UNROLL) ? lon1_ + lon12 : math::ang_normalize(math::ang_normalize(lon1_) + math::ang_normalize(lon12));
    }
    if (outmask & caps::LATITUDE) r.lat2 = math::atan2d(sbet2, f1_ * cbet2);
    if (outmask & caps::AZIMUTH) r.azi2 = math::atan2d(salp2, calp2);
    if (outmask & (caps::REDUCEDLENGTH | caps::GEODESICSCALE)) {
      const double B22 = Geodesic::sin_cos_series(true, ssig2, csig2, C2a_, Geodesic::nC2 + 1);
      const double AB2 = (1 + A2m1_) * (B22 - B21_);
      const double J12 = (A1m1_ - A2m1_) * sig12 + (AB1 - AB2);
      if (outmask & caps::REDUCEDLENGTH) r.m12 = b_ * ((dn2 * (csig1_ * ssig2) - dn1_ * (ssig1_ * csig2)) - csig1_ * csig2 * J12);
      if (outmask & caps::GEODESICSCALE) {
        const double t = k2_ * (ssig2 - ssig1_) * (ssig2 + ssig1_) / (dn1_ + dn2);
        r.M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1_ / dn1_;
        r.M21 = csig12 - (t * ssig1_ - csig1_ * J12) * ssig2 / dn2;
      }
    }
    if (outmask & caps::AREA) {
      const double B42 = Geodesic::sin_cos_series(false, ssig2, csig2, C4a_, Geodesic::nC4);
      double salp12, calp12;
      if (calp0_ == 0 || salp0_ == 0) {
        salp12 = salp2 * calp1_ - calp2 * salp1_;
        calp12 = calp2 * calp1_ + salp2 * salp1_;
      } else {
        // tan(α₂ − α₁) from tan α = tan α₀ sec σ, written to avoid cancellation.
        salp12 = calp0_ * salp0_ *
                 (csig12 <= 0 ? csig1_ * (1 - csig12) + ssig12 * ssig1_ : ssig12 * (csig1_ * ssig12 / (1 + csig12) + ssig1_));
        calp12 = math::sq(salp0_) + math::sq(calp0_) * csig1_ * csig2;
      }
      r.S12 = c2_ * std::atan2(salp12, calp12) + A4_ * (B42 - B41_);
    }
    r.a12 = arcmode ? s12_a12 : sig12 * math::rad_to_deg;
    return r;
  }

 private:
  static Direct to_radians(Direct d) {
    d.lat2 *= math::deg_to_rad;
    d.lon2 *= math::deg_to_rad;
    d.azi2 *= math::deg_to_rad;
    d.a12 *= math::deg_to_rad;
    return d;
  }

  double a_, f_, b_, c2_, f1_;
  unsigned caps_;
  double lat1_, lon1_, azi1_, salp1_, calp1_;
  double dn1_, salp0_, calp0_, ssig1_, csig1_, somg1_, comg1_, k2_;
  double A1m1_ = 0, B11_ = 0, stau1_ = 0, ctau1_ = 0, A2m1_ = 0, B21_ = 0, A3c_ = 0, B31_ = 0, A4_ = 0, B41_ = 0;
  double C1a_[Geodesic::nC1 + 1] = {0}, C1pa_[Geodesic::nC1p + 1] = {0}, C2a_[Geodesic::nC2 + 1] = {0};
  double C3a_[Geodesic::nC3] = {0}, C4a_[Geodesic::nC4] = {0};
  double s13_ = math::nan, a13_ = math::nan;
};

// ---- the problems through lines -------------------------------------------------------------------

inline Direct Geodesic::direct_degrees(double lat1, double lon1, double azi1, double s12_a12, bool arcmode, bool unroll) const {
  // The distance as input needs the reverted series (DISTANCE_IN).
  const GeodesicLine l(*this, lat1, lon1, azi1, caps::ALL | (arcmode ? 0u : caps::DISTANCE_IN));
  return l.gen_position(arcmode, s12_a12, caps::ALL | (unroll ? caps::LONG_UNROLL : 0));
}

inline Direct Geodesic::direct(double lat1, double lon1, double azi1, double s12, bool unroll) const {
  Direct d = direct_degrees(to_lat_degrees(lat1), lon1 * math::rad_to_deg, azi1 * math::rad_to_deg, s12, false, unroll);
  d.lat2 *= math::deg_to_rad;
  d.lon2 *= math::deg_to_rad;
  d.azi2 *= math::deg_to_rad;
  d.a12 *= math::deg_to_rad;
  return d;
}

inline Direct Geodesic::arc_direct(double lat1, double lon1, double azi1, double a12, bool unroll) const {
  Direct d = direct_degrees(to_lat_degrees(lat1), lon1 * math::rad_to_deg, azi1 * math::rad_to_deg, a12 * math::rad_to_deg, true, unroll);
  d.lat2 *= math::deg_to_rad;
  d.lon2 *= math::deg_to_rad;
  d.azi2 *= math::deg_to_rad;
  d.a12 *= math::deg_to_rad;
  return d;
}

inline GeodesicLine Geodesic::line(double lat1, double lon1, double azi1) const {
  return GeodesicLine(*this, to_lat_degrees(lat1), lon1 * math::rad_to_deg, azi1 * math::rad_to_deg);
}

inline GeodesicLine Geodesic::inverse_line(double lat1, double lon1, double lat2, double lon2) const {
  const double lat1d = to_lat_degrees(lat1), lon1d = lon1 * math::rad_to_deg;
  const General g = gen_inverse(lat1d, lon1d, to_lat_degrees(lat2), lon2 * math::rad_to_deg, 0);
  const double azi1 = math::atan2d(g.salp1, g.calp1);
  GeodesicLine l(*this, lat1d, lon1d, azi1, caps::ALL | caps::DISTANCE_IN | caps::DISTANCE, g.salp1, g.calp1);
  l.set_arc_degrees(g.a12);
  return l;
}

}  // namespace geodesy::geodesic
