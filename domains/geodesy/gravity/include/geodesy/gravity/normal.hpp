// The normal gravity field of a level ellipsoid: its derived constants, Somigliana's formula on the
// ellipsoid, and the exact closed field at any point outside it, with its vector.
//
// The level ellipsoid. An ellipsoid of revolution (semi-axes a > b) that is a surface of constant
// gravity potential U = U₀ of the field it generates, rotating at ω and enclosing the mass M: by the
// theorem of Stokes and Poincaré its exterior potential is fixed by four constants, whatever the
// density inside (Moritz 1980/2000, "Geodetic Reference System 1980", J. Geodesy 74, 128, §2;
// Heiskanen and Moritz 1967, Physical Geodesy, §§2-7 to 2-10, whose equation numbers Moritz cites
// and which are repeated here as "HM"). WGS 84 takes a, f, GM, ω (NGA.STND.0036_1.0.0_WGS84, 2014,
// Table 3.1, which keeps TR8350.2's values); GRS 80 takes a, J₂, GM, ω (IUGG 1979, Resolution 7).
//
// Ellipsoidal-harmonic coordinates (u, β, λ): u is the semi-minor axis of the ellipsoid confocal
// with the reference one (linear eccentricity E = √(a² − b²)) through the point, β the reduced
// latitude on it, λ the longitude (NGA eqs. 4-7 to 4-9). With x = E/u,
//   q(x)  = ½[(1 + 3/x²) arctan x − 3/x],      q′(x) = 3(1 + 1/x²)(1 − arctan(x)/x) − 1
// (NGA eqs. 4-11 to 4-13; HM 2-58, 2-67), q₀ = q(e′) and q₀′ = q′(e′) on the reference ellipsoid,
// e′ = E/b. Both differences cancel badly for small x (q ~ 2x³/15, and x ≤ e′ = 0.082 outside the
// ellipsoid), so they are summed from the arctangent's series, as Moritz gives 2q₀ (§3):
//   q  = Σₙ≥₁ (−1)ⁿ⁺¹ 2n x²ⁿ⁺¹/((2n+1)(2n+3)),   q′ = Σₙ≥₁ (−1)ⁿ⁺¹ 6 x²ⁿ/((2n+1)(2n+3)).
//
// Derived constants (Moritz 2000, §3): m = ω² a² b/GM (HM 2-70); the normal potential on the
// ellipsoid U₀ = (GM/E) arctan e′ + ω² a²/3 (HM 2-61); normal gravity at the equator and the poles
//   γₑ = GM/(ab) (1 − m − (m/6) e′ q₀′/q₀),   γₚ = GM/a² (1 + (m/3) e′ q₀′/q₀)   (HM 2-73, 2-74);
// k = bγₚ/(aγₑ) − 1; the gravity flattening f* = (γₚ − γₑ)/γₑ; the zonal coefficients
//   J₂ = (e²/3)(1 − (2/15) m e′/q₀),   J₂ₙ = (−1)ⁿ⁺¹ 3e²ⁿ (1 − n + 5n J₂/e²)/((2n+1)(2n+3))
// (HM 2-90, 2-92, 2-92′), and, given J₂, e² by Moritz's iteration e² = 3J₂ + (4/15)(ω²a³/GM) e³/(2q₀);
// the mean of normal gravity over the ellipsoid's surface,
//   γ̄ = ∫₀^{π/2} γ(φ) cos φ (1 − e² sin² φ)⁻² dφ / ∫₀^{π/2} cos φ (1 − e² sin² φ)⁻² dφ,
// by calx's adaptive Gauss–Kronrod quadrature.
//
// Normal gravity. On the ellipsoid, Somigliana's closed formula γ = γₑ (1 + k sin² φ)/√(1 − e² sin² φ)
// (NGA eq. 4-1; Moritz §3) at geodetic latitude φ. Off it, the potential
//   U(u, β) = (GM/E) arctan(E/u) + ½ ω² a² (q/q₀)(sin² β − ⅓) + ½ ω² (u² + E²) cos² β
// (Hofmann-Wellenhof and Moritz 2006, Physical Geodesy, 2nd ed., eq. 2-126), and its gradient in the
// orthogonal ellipsoidal frame, with w = √((u² + E² sin² β)/(u² + E²)) (NGA eq. 4-10),
//   γᵤ = −(1/w)[GM/(u² + E²) + (ω² a² E/(u² + E²))(q′/q₀)(½ sin² β − ⅙)] + (1/w) ω² u cos² β   (4-5)
//   γᵦ = (1/w)(ω² a²/√(u² + E²))(q/q₀) sin β cos β − (1/w) ω² √(u² + E²) sin β cos β       (4-6)
// in closed form, exact at every height ("with unlimited numerical accuracy", NGA §4.3). The vector
// γ = γᵤ eᵤ + γᵦ e_β = ∇U points down; its magnitude is √(γᵤ² + γᵦ²) (NGA eq. 4-4, there called an
// approximation of the component along the normal, which it is to 1 μGal below 20 km). Its components
// along the downward geodetic normal (γₕ) and northward in the meridian (γ_φ) are NGA's eqs. 4-16 and
// 4-23, here computed by projecting the vector on the local frame of geodesy::figure (ecef_to_enu),
// which is the same rotation as NGA's R₂R₁ (eqs. 4-17 to 4-20). The direction of γ departs from the
// geodetic normal by less than 4″ below 20 km (NGA §4.3); the second-order series in height
//   γₕ ≈ γ [1 − (2/a)(1 + f + m − 2f sin² φ) h + (3/a²) h²]   (NGA eq. 4-3)
// is given for comparison.
//
// Validity. The field is that of the ellipsoid's exterior: harmonic (with the centrifugal term)
// outside and on it, and the reference of gravimetry there. At points below the ellipsoid (h < 0)
// the same closed formulas give its analytic continuation, which is the field of no mass
// distribution in particular but is continuous through the surface, as free-air reductions assume;
// at the focal disc (u = 0) they fail, and the functions throw. Normal gravity is the reference of
// the real field, not the real field: they differ by the gravity disturbance, tens of mGal
// (10⁻⁴ m s⁻²) and up to a few hundred over mountains and trenches.

#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

#include <calx/quadrature/gauss_kronrod.hpp>
#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>

namespace geodesy::gravity {

using figure::Ellipsoid;
using figure::Geodetic;
using physis::geometry::Vec3;

namespace detail {
// q(x) and q′(x) of the header, x = E/u > 0: by their series below x = 1/2 (alternating, each term about
// x² of the last), in closed form above, where the closed forms keep some 1e-14.
inline double q_of(double x) {
  if (x >= 0.5) return 0.5 * ((1 + 3 / (x * x)) * std::atan(x) - 3 / x);
  const double x2 = x * x;
  double sum = 0, power = x * x2, sign = 1;
  for (int n = 1; n < 60; n++, power *= x2, sign = -sign) {
    const double term = sign * 2.0 * n * power / ((2.0 * n + 1) * (2.0 * n + 3));
    sum += term;
    if (std::fabs(term) < 1e-18 * std::fabs(sum)) break;
  }
  return sum;
}
inline double q_prime_of(double x) {
  if (x >= 0.5) return 3 * (1 + 1 / (x * x)) * (1 - std::atan(x) / x) - 1;
  const double x2 = x * x;
  double sum = 0, power = x2, sign = 1;
  for (int n = 1; n < 60; n++, power *= x2, sign = -sign) {
    const double term = sign * 6.0 * power / ((2.0 * n + 1) * (2.0 * n + 3));
    sum += term;
    if (std::fabs(term) < 1e-18 * std::fabs(sum)) break;
  }
  return sum;
}
}  // namespace detail

// Ellipsoidal-harmonic coordinates of a point (NGA eqs. 4-8, 4-9).
struct EllipsoidalHarmonic {
  double u = 0;                 // semi-minor axis of the confocal ellipsoid through the point (m)
  double reduced_latitude = 0;  // β (rad)
  double longitude = 0;         // λ (rad), east
};

// Normal gravity at a point: its components and its vector (m s⁻²).
struct NormalGravity {
  double magnitude_m_s2 = 0;   // |γ| = √(γᵤ² + γᵦ²) (NGA eq. 4-4)
  double u_m_s2 = 0;           // γᵤ (eq. 4-5), negative: toward the ellipsoid
  double beta_m_s2 = 0;        // γᵦ (eq. 4-6), toward increasing reduced latitude (northward)
  double down_m_s2 = 0;        // γₕ, along the downward geodetic normal (eq. 4-16)
  double north_m_s2 = 0;       // γ_φ, northward in the meridian, perpendicular to the normal (eq. 4-23)
  double deflection_rad = 0;   // the angle of γ from the downward normal, positive when γ leans north
  Vec3 ecef_m_s2{};            // γ = ∇U in the Earth-fixed axes of the ellipsoid
  Vec3 enu_m_s2{};             // γ in the local east-north-up frame of the point
};

class LevelEllipsoid {
 public:
  // From a (m), 1/f, GM (m³ s⁻²) and ω (rad s⁻¹), as WGS 84 is defined.
  static LevelEllipsoid from_flattening(std::string name, double a, double inverse_flattening, double GM, double omega) {
    return LevelEllipsoid(Ellipsoid(name, a, inverse_flattening), GM, omega);
  }
  // From a (m), J₂, GM (m³ s⁻²) and ω (rad s⁻¹), as GRS 80 is defined: e² by Moritz's iteration (§3).
  static LevelEllipsoid from_J2(std::string name, double a, double J2, double GM, double omega) {
    if (!(J2 > 0)) throw std::invalid_argument("geodesy: a level ellipsoid from J2 needs J2 > 0");
    const double c = 4.0 / 15.0 * omega * omega * a * a * a / GM;
    double e2 = 3 * J2;
    for (int i = 0; i < 100; i++) {
      const double e = std::sqrt(e2), ep = e / std::sqrt(1 - e2);
      const double next = 3 * J2 + c * e2 * e / (2 * detail::q_of(ep));
      const bool done = std::fabs(next - e2) <= 1e-17;
      e2 = next;
      if (done) break;
    }
    return LevelEllipsoid(Ellipsoid(name, a, 1 / (1 - std::sqrt(1 - e2))), GM, omega);
  }
  // WGS 84: a = 6 378 137 m, 1/f = 298.257 223 563, GM = 3.986 004 418 × 10¹⁴ m³ s⁻²,
  // ω = 7.292 115 × 10⁻⁵ rad s⁻¹ (NGA.STND.0036_1.0.0_WGS84, Table 3.1).
  static LevelEllipsoid wgs84() { return from_flattening("WGS 84", 6378137.0, 298.257223563, 3.986004418e14, 7.292115e-5); }
  // GRS 80: a = 6 378 137 m, J₂ = 108 263 × 10⁻⁸, GM = 3 986 005 × 10⁸ m³ s⁻², ω = 7 292 115 × 10⁻¹¹
  // rad s⁻¹ (IUGG 1979, Resolution 7, as Moritz 2000 quotes it).
  static LevelEllipsoid grs80() { return from_J2("GRS 80", 6378137.0, 108263e-8, 3986005e8, 7292115e-11); }

  const std::string& name() const { return figure_.name(); }
  const Ellipsoid& figure() const { return figure_; }

  // ---- defining and derived geometric constants ----
  double a() const { return figure_.a(); }                    // m
  double b() const { return figure_.b(); }                    // m
  double flattening() const { return figure_.flattening(); }
  double GM() const { return GM_; }                           // m³ s⁻²
  double omega() const { return omega_; }                     // rad s⁻¹
  double e2() const { return figure_.e2(); }
  double second_e2() const { return figure_.second_e2(); }
  double second_eccentricity() const { return ep_; }          // e′ = E/b
  double linear_eccentricity() const { return E_; }           // E (m)
  double polar_radius_of_curvature() const { return a() * a() / b(); }  // c = a²/b (m)

  // ---- derived physical constants (Moritz 2000, §3) ----
  double m() const { return m_; }                             // ω² a² b/GM
  double q0() const { return q0_; }
  double q0_prime() const { return q0p_; }
  double U0() const { return GM_ / E_ * std::atan(ep_) + omega_ * omega_ * a() * a() / 3; }  // m² s⁻²
  double gamma_equator() const { return GM_ / (a() * b()) * (1 - m_ - m_ / 6 * ep_ * q0p_ / q0_); }  // m s⁻²
  double gamma_pole() const { return GM_ / (a() * a()) * (1 + m_ / 3 * ep_ * q0p_ / q0_); }           // m s⁻²
  double k() const { return b() * gamma_pole() / (a() * gamma_equator()) - 1; }
  double gravity_flattening() const { return (gamma_pole() - gamma_equator()) / gamma_equator(); }  // f*
  double J2() const { return e2() / 3 * (1 - 2.0 / 15.0 * m_ * ep_ / q0_); }
  // J₂ₙ for n ≥ 1 (n = 1 is J₂ itself; HM 2-92′).
  double J2n(int n) const {
    if (n < 1) throw std::invalid_argument("geodesy: J2n needs n >= 1");
    if (n == 1) return J2();
    return (n % 2 ? 1.0 : -1.0) * 3 * std::pow(e2(), n) * (1 - n + 5 * n * J2() / e2()) / ((2.0 * n + 1) * (2.0 * n + 3));
  }
  // The mean of normal gravity over the ellipsoid's surface (m s⁻²).
  double mean_gravity() const {
    const double e2_ = e2();
    auto area = [e2_](double phi) { const double w = 1 - e2_ * std::sin(phi) * std::sin(phi); return std::cos(phi) / (w * w); };
    const double half_pi = std::numbers::pi / 2;
    const auto num = calx::quadrature::integrate([&](double phi) { return somigliana(phi) * area(phi); }, 0.0, half_pi, 0.0, 1e-15);
    const auto den = calx::quadrature::integrate(area, 0.0, half_pi, 0.0, 1e-15);
    return num.value / den.value;
  }

  // ---- normal gravity ----
  // On the ellipsoid, at geodetic latitude φ (rad): Somigliana (NGA eq. 4-1), m s⁻².
  double somigliana(double latitude) const {
    const double s2 = std::sin(latitude) * std::sin(latitude);
    return gamma_equator() * (1 + k() * s2) / std::sqrt(1 - e2() * s2);
  }
  // The second-order series in height h (m) above the ellipsoid (NGA eq. 4-3), m s⁻².
  double series(double latitude, double height) const {
    const double s2 = std::sin(latitude) * std::sin(latitude), f = flattening();
    return somigliana(latitude) * (1 - 2 / a() * (1 + f + m_ - 2 * f * s2) * height + 3 / (a() * a()) * height * height);
  }

  // The ellipsoidal-harmonic coordinates of an Earth-centred point (NGA eqs. 4-8, 4-9).
  EllipsoidalHarmonic ellipsoidal_harmonic(const Vec3& x) const {
    const double rho2 = x[0] * x[0] + x[1] * x[1], r2 = rho2 + x[2] * x[2], d = r2 - E_ * E_;
    EllipsoidalHarmonic h;
    // u² = ½ d [1 + √(1 + 4E²z²/d²)], written without cancellation for d < 0 (deep inside) as well.
    const double root = std::sqrt(d * d + 4 * E_ * E_ * x[2] * x[2]);
    const double u2 = d >= 0 ? 0.5 * (d + root) : 2 * E_ * E_ * x[2] * x[2] / (root - d);
    if (!(u2 > 0)) throw std::domain_error("geodesy: normal gravity is undefined on the focal disc of the ellipsoid");
    h.u = std::sqrt(u2);
    h.reduced_latitude = std::atan2(x[2] * std::sqrt(u2 + E_ * E_), h.u * std::sqrt(rho2));
    h.longitude = rho2 > 0 ? std::atan2(x[1], x[0]) : 0.0;
    return h;
  }

  // The normal potential U (m² s⁻²) at a point given geodetically (HWM 2006 eq. 2-126).
  double potential(const Geodetic& p) const {
    const EllipsoidalHarmonic h = ellipsoidal_harmonic(figure_.cartesian(p));
    const double sb = std::sin(h.reduced_latitude), cb = std::cos(h.reduced_latitude), w2 = omega_ * omega_;
    return GM_ / E_ * std::atan(E_ / h.u) + 0.5 * w2 * a() * a() * detail::q_of(E_ / h.u) / q0_ * (sb * sb - 1.0 / 3.0) +
           0.5 * w2 * (h.u * h.u + E_ * E_) * cb * cb;
  }

  // Normal gravity at a point given geodetically: the exact closed field (NGA eqs. 4-4 to 4-13, 4-16,
  // 4-23), with its vector in Earth-fixed axes and in the point's east-north-up frame.
  NormalGravity at(const Geodetic& p) const {
    const EllipsoidalHarmonic h = ellipsoidal_harmonic(figure_.cartesian(p));
    const double u = h.u, sb = std::sin(h.reduced_latitude), cb = std::cos(h.reduced_latitude);
    const double sl = std::sin(h.longitude), cl = std::cos(h.longitude);
    const double s = u * u + E_ * E_, rs = std::sqrt(s), w = std::sqrt((u * u + E_ * E_ * sb * sb) / s), w2 = omega_ * omega_;
    const double x = E_ / u, q = detail::q_of(x), qp = detail::q_prime_of(x);
    NormalGravity g;
    g.u_m_s2 = -(GM_ / s + w2 * a() * a() * E_ / s * qp / q0_ * (0.5 * sb * sb - 1.0 / 6.0)) / w + w2 * u * cb * cb / w;
    g.beta_m_s2 = (w2 * a() * a() / rs * q / q0_ - w2 * rs) * sb * cb / w;
    g.magnitude_m_s2 = std::hypot(g.u_m_s2, g.beta_m_s2);
    // eᵤ and e_β in Earth-fixed axes (the columns of NGA's R₁, eq. 4-18).
    const Vec3 eu{u * cb * cl / (w * rs), u * cb * sl / (w * rs), sb / w};
    const Vec3 eb{-sb * cl / w, -sb * sl / w, u * cb / (w * rs)};
    for (int i = 0; i < 3; i++) g.ecef_m_s2[i] = g.u_m_s2 * eu[i] + g.beta_m_s2 * eb[i];
    using namespace physis::geometry::operators;
    g.enu_m_s2 = figure::ecef_to_enu(p.latitude, p.longitude) * g.ecef_m_s2;
    g.down_m_s2 = -g.enu_m_s2[2];
    g.north_m_s2 = g.enu_m_s2[1];
    g.deflection_rad = std::atan2(g.north_m_s2, g.down_m_s2);
    return g;
  }

 private:
  LevelEllipsoid(Ellipsoid figure, double GM, double omega) : figure_(std::move(figure)), GM_(GM), omega_(omega) {
    if (!(GM > 0) || !(omega >= 0) || figure_.is_sphere()) throw std::invalid_argument("geodesy: a level ellipsoid needs GM > 0, omega >= 0 and f > 0");
    E_ = figure_.a() * std::sqrt(figure_.e2());
    ep_ = E_ / figure_.b();
    q0_ = detail::q_of(ep_);
    q0p_ = detail::q_prime_of(ep_);
    m_ = omega_ * omega_ * figure_.a() * figure_.a() * figure_.b() / GM_;
  }

  Ellipsoid figure_;
  double GM_ = 0, omega_ = 0, E_ = 0, ep_ = 0, q0_ = 0, q0p_ = 0, m_ = 0;
};

}  // namespace geodesy::gravity
