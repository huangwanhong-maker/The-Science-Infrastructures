// The orientation of the Earth in space: the IAU 2006/2000A precession-nutation, the Earth
// rotation angle and sidereal time, and polar motion; the matrices between the celestial and
// terrestrial frames by the two routes the IERS defines.
//
// Sources: IERS Conventions 2010 (Petit and Luzum, IERS Technical Note 36), ch. 5; Capitaine,
// Wallace and Chapront 2003 (A&A 412, 567) and Wallace and Capitaine 2006 (A&A 459, 981) for
// the IAU 2006 precession in Fukushima–Williams angles; Mathews, Herring and Buffett 2002 for
// the IAU 2000A nutation (series_nutation.hpp), with the IAU 2006 adjustments for J2's secular
// change; Capitaine, Wallace and Chapront 2005 for the CIO locator s (series_cio.hpp); SOFA's
// implementations (IAU SOFA Board 2023: iauPfw06, iauFw2m, iauNut06a, iauS06, iauEra00,
// iauGmst06, iauEect00, iauEors, iauSp00, iauPom00, iauC2ixys), which these follow step by
// step and against whose test values the checks compare them.
//
// Matrices are written as SOFA writes them: r maps the coordinates of a vector in one frame to
// its coordinates in the next, v_next = r v; the elementary rotations are rotations of the
// frame (physis::geometry::frame_rotation, SOFA's R_x, R_y, R_z).
//
// The two routes from the Geocentric Celestial Reference System to the International
// Terrestrial Reference System (IERS Conventions 2010, §5.4):
//   CIO-based:     ITRS ← W(x_p, y_p, s′) · R_z(θ) · C(X, Y, s) ← GCRS, θ the Earth rotation
//                  angle of UT1 (the Celestial Intermediate Origin route; CIRS, TIRS between);
//   equinox-based: ITRS ← W · R_z(GAST) · N P B ← GCRS, GAST the apparent sidereal time
//                  (the true equator and equinox of date, TOD, between).
// GAST computed as GMST + the equation of the equinoxes (with its complementary terms) is an
// independent route to the same matrix; the two agree to about a microarcsecond, and the frame
// graph measures their difference as the holonomy of the loop they close.

#pragma once

#include <cmath>
#include <numbers>

#include <metrology/frames/series_cio.hpp>
#include <metrology/frames/series_nutation.hpp>
#include <metrology/time/julian.hpp>
#include <physis/geometry/rotation.hpp>

namespace metrology::frames {

using physis::geometry::Axis;
using physis::geometry::frame_rotation;
using physis::geometry::Mat3;
using physis::geometry::Vec3;
using namespace physis::geometry::operators;

inline constexpr double PI = std::numbers::pi;
inline constexpr double TWO_PI = 2 * std::numbers::pi;
inline constexpr double ARCSEC = PI / 648000.0;   // radians per arcsecond
inline constexpr double TURNAS = 1296000.0;       // arcseconds in a full circle

// The angle in [0, 2π) and in [−π, π).
inline double normalised_angle(double a) {
  const double w = std::fmod(a, TWO_PI);
  return w < 0 ? w + TWO_PI : w;
}
inline double signed_angle(double a) {
  double w = std::fmod(a, TWO_PI);
  if (std::fabs(w) >= PI) w -= std::copysign(TWO_PI, a);
  return w;
}

// Julian centuries of TT since J2000.0.
inline double centuries(const time::JulianDate& tt) { return tt.centuries_since_j2000(); }

// ---- fundamental arguments (IERS Conventions 2003, as SOFA's iauFal03 … iauFapa03) ------------
namespace fundamental {
inline double l(double t) { return std::fmod(485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470)))), TURNAS) * ARCSEC; }
inline double lp(double t) { return std::fmod(1287104.793048 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))), TURNAS) * ARCSEC; }
inline double f(double t) { return std::fmod(335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))), TURNAS) * ARCSEC; }
inline double d(double t) { return std::fmod(1072260.703692 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))), TURNAS) * ARCSEC; }
inline double om(double t) { return std::fmod(450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939)))), TURNAS) * ARCSEC; }
inline double me(double t) { return std::fmod(4.402608842 + 2608.7903141574 * t, TWO_PI); }
inline double ve(double t) { return std::fmod(3.176146697 + 1021.3285546211 * t, TWO_PI); }
inline double ea(double t) { return std::fmod(1.753470314 + 628.3075849991 * t, TWO_PI); }
inline double ma(double t) { return std::fmod(6.203480913 + 334.0612426700 * t, TWO_PI); }
inline double ju(double t) { return std::fmod(0.599546497 + 52.9690962641 * t, TWO_PI); }
inline double sa(double t) { return std::fmod(0.874016757 + 21.3299104960 * t, TWO_PI); }
inline double ur(double t) { return std::fmod(5.481293872 + 7.4781598567 * t, TWO_PI); }
inline double pa(double t) { return (0.024381750 + 0.00000538691 * t) * t; }
}  // namespace fundamental

// ---- precession -------------------------------------------------------------------------------

// The mean obliquity of the ecliptic, IAU 2006 (iauObl06), radians.
inline double mean_obliquity(const time::JulianDate& tt) {
  const double t = centuries(tt);
  return (84381.406 + t * (-46.836769 + t * (-0.0001831 + t * (0.00200340 + t * (-0.000000576 + t * (-0.0000000434)))))) * ARCSEC;
}

// The Fukushima–Williams angles of bias and precession, IAU 2006 (iauPfw06): γ̄, φ̄, ψ̄ and ε_A.
struct FukushimaWilliams {
  double gamb, phib, psib, epsa;
};
inline FukushimaWilliams precession_angles(const time::JulianDate& tt) {
  const double t = centuries(tt);
  return {(-0.052928 + t * (10.556378 + t * (0.4932044 + t * (-0.00031238 + t * (-0.000002788 + t * (0.0000000260)))))) * ARCSEC,
          (84381.412819 + t * (-46.811016 + t * (0.0511268 + t * (0.00053289 + t * (-0.000000440 + t * (-0.0000000176)))))) * ARCSEC,
          (-0.041775 + t * (5038.481484 + t * (1.5584175 + t * (-0.00018522 + t * (-0.000026452 + t * (-0.0000000148)))))) * ARCSEC,
          mean_obliquity(tt)};
}
// The matrix of Fukushima–Williams angles (iauFw2m): R_x(−ε) R_z(−ψ) R_x(φ) R_z(γ).
inline Mat3 fukushima_williams_matrix(double gamb, double phib, double psi, double eps) {
  return frame_rotation(Axis::x, -eps) * frame_rotation(Axis::z, -psi) * frame_rotation(Axis::x, phib) * frame_rotation(Axis::z, gamb);
}
// The bias-precession matrix: GCRS to the mean equator and equinox of date (iauPmat06).
inline Mat3 bias_precession(const time::JulianDate& tt) {
  const auto a = precession_angles(tt);
  return fukushima_williams_matrix(a.gamb, a.phib, a.psib, a.epsa);
}
// The frame bias: GCRS to the mean equator and equinox of J2000.0 (EME2000), the bias-precession
// matrix at J2000.0 (iauBp06's B).
inline Mat3 frame_bias() { return bias_precession({time::J2000, 0.0}); }

// ---- nutation -------------------------------------------------------------------------------

struct Nutation {
  double dpsi, deps;  // in longitude and obliquity, radians
};
// The IAU 2000A nutation (iauNut00a): 678 luni-solar and 687 planetary terms.
inline Nutation nutation_2000a(const time::JulianDate& tt) {
  const double t = centuries(tt);
  const double u2r = ARCSEC / 1e7;  // 0.1 µas to radians
  const double el = fundamental::l(t);
  // Mean anomaly of the Sun and the elongation of the Moon as MHB2000 has them.
  const double elp = std::fmod(1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))), TURNAS) * ARCSEC;
  const double f = fundamental::f(t);
  const double d = std::fmod(1072260.70369 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))), TURNAS) * ARCSEC;
  const double om = fundamental::om(t);
  double dp = 0, de = 0;
  for (int i = 677; i >= 0; i--) {
    const auto& x = series::luni_solar[i];
    const double arg = std::fmod(x.nl * el + x.nlp * elp + x.nf * f + x.nd * d + x.nom * om, TWO_PI);
    const double s = std::sin(arg), c = std::cos(arg);
    dp += (x.sp + x.spt * t) * s + x.cp * c;
    de += (x.ce + x.cet * t) * c + x.se * s;
  }
  const double dpsils = dp * u2r, depsls = de * u2r;
  // Planetary nutation, with MHB2000's arguments for the Moon and Neptune.
  const double al = std::fmod(2.35555598 + 8328.6914269554 * t, TWO_PI);
  const double af = std::fmod(1.627905234 + 8433.466158131 * t, TWO_PI);
  const double ad = std::fmod(5.198466741 + 7771.3771468121 * t, TWO_PI);
  const double aom = std::fmod(2.18243920 - 33.757045 * t, TWO_PI);
  const double apa = fundamental::pa(t);
  const double alme = fundamental::me(t), alve = fundamental::ve(t), alea = fundamental::ea(t), alma = fundamental::ma(t),
               alju = fundamental::ju(t), alsa = fundamental::sa(t), alur = fundamental::ur(t);
  const double alne = std::fmod(5.321159000 + 3.8127774000 * t, TWO_PI);
  dp = de = 0;
  for (int i = 686; i >= 0; i--) {
    const auto& x = series::planetary[i];
    const double arg = std::fmod(x.nl * al + x.nf * af + x.nd * ad + x.nom * aom + x.nme * alme + x.nve * alve + x.nea * alea +
                                     x.nma * alma + x.nju * alju + x.nsa * alsa + x.nur * alur + x.nne * alne + x.npa * apa,
                                 TWO_PI);
    const double s = std::sin(arg), c = std::cos(arg);
    dp += x.sp * s + x.cp * c;
    de += x.se * s + x.ce * c;
  }
  return {dpsils + dp * u2r, depsls + de * u2r};
}
// The IAU 2000A nutation with the IAU 2006 adjustments (iauNut06a; Wallace and Capitaine 2006,
// eq. 5): the secular change of J2 and the P03 scale.
inline Nutation nutation(const time::JulianDate& tt) {
  const double t = centuries(tt), fj2 = -2.7774e-6 * t;
  const Nutation n = nutation_2000a(tt);
  return {n.dpsi + n.dpsi * (0.4697e-6 + fj2), n.deps + n.deps * fj2};
}

// The classical bias-precession-nutation matrix: GCRS to the true equator and equinox of date
// (iauPnm06a).
inline Mat3 bias_precession_nutation(const time::JulianDate& tt, const Nutation& n) {
  const auto a = precession_angles(tt);
  return fukushima_williams_matrix(a.gamb, a.phib, a.psib + n.dpsi, a.epsa + n.deps);
}
inline Mat3 bias_precession_nutation(const time::JulianDate& tt) { return bias_precession_nutation(tt, nutation(tt)); }

// ---- the Celestial Intermediate Pole and Origin ------------------------------------------------

// The CIP's X, Y: the bottom row of the NPB matrix (iauBpn2xy).
inline std::array<double, 2> cip_xy(const Mat3& npb) { return {npb[2][0], npb[2][1]}; }

// The CIO locator s, given X and Y (iauS06).
inline double cio_locator(const time::JulianDate& tt, double x, double y) {
  const double t = centuries(tt);
  const double fa[8] = {fundamental::l(t), fundamental::lp(t), fundamental::f(t), fundamental::d(t),
                        fundamental::om(t), fundamental::ve(t), fundamental::ea(t), fundamental::pa(t)};
  auto sum = [&](const series::ArgumentTerm* terms, int n, double w) {
    for (int i = n - 1; i >= 0; i--) {
      double a = 0;
      for (int j = 0; j < 8; j++) a += terms[i].nfa[j] * fa[j];
      w += terms[i].s * std::sin(a) + terms[i].c * std::cos(a);
    }
    return w;
  };
  const double* sp = series::s06_polynomial;
  const double w0 = sum(series::s06_s0, 33, sp[0]), w1 = sum(series::s06_s1, 3, sp[1]), w2 = sum(series::s06_s2, 25, sp[2]),
               w3 = sum(series::s06_s3, 4, sp[3]), w4 = sum(series::s06_s4, 1, sp[4]), w5 = sp[5];
  return (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t) * ARCSEC - x * y / 2;
}

// The matrix of the celestial-to-intermediate transformation, GCRS to CIRS, from X, Y and s
// (iauC2ixys): R_z(−(E + s)) R_y(d) R_z(E).
inline Mat3 celestial_to_intermediate(double x, double y, double s) {
  const double r2 = x * x + y * y;
  const double e = r2 > 0 ? std::atan2(y, x) : 0.0;
  const double d = std::atan(std::sqrt(r2 / (1 - r2)));
  return frame_rotation(Axis::z, -(e + s)) * frame_rotation(Axis::y, d) * frame_rotation(Axis::z, e);
}

// The equation of the origins, from the NPB matrix and s (iauEors; Wallace and Capitaine 2006,
// eq. 16): the angle from the CIO to the equinox along the intermediate equator.
inline double equation_of_origins(const Mat3& npb, double s) {
  const double x = npb[2][0], ax = x / (1 + npb[2][2]);
  const double xs = 1 - ax * x, ys = -ax * npb[2][1], zs = -x;
  const double p = npb[0][0] * xs + npb[0][1] * ys + npb[0][2] * zs;
  const double q = npb[1][0] * xs + npb[1][1] * ys + npb[1][2] * zs;
  return (p != 0 || q != 0) ? s - std::atan2(q, p) : s;
}

// ---- the rotation of the Earth ------------------------------------------------------------------

// The Earth rotation angle θ of UT1, in [0, 2π) (iauEra00; IERS Conventions 2010, eq. 5.15).
inline double earth_rotation_angle(const time::JulianDate& ut1) {
  const double d1 = std::min(ut1.d1, ut1.d2), d2 = std::max(ut1.d1, ut1.d2);
  const double t = d1 + (d2 - time::J2000);
  const double f = std::fmod(d1, 1.0) + std::fmod(d2, 1.0);
  return normalised_angle(TWO_PI * (f + 0.7790572732640 + 0.00273781191135448 * t));
}
// dθ/dt in rad per SI second of UT1, and corrected for the excess length of day.
inline double earth_rotation_rate(double lod_seconds = 0) {
  return TWO_PI * 1.00273781191135448 / time::DAY * (1 - lod_seconds / time::DAY);
}
// Greenwich mean sidereal time, IAU 2006 (iauGmst06).
inline double gmst(const time::JulianDate& ut1, const time::JulianDate& tt) {
  const double t = centuries(tt);
  return normalised_angle(earth_rotation_angle(ut1) +
                          (0.014506 + t * (4612.156534 + t * (1.3915817 + t * (-0.00000044 + t * (-0.000029956 + t * (-0.0000000368)))))) * ARCSEC);
}
// The complementary terms of the equation of the equinoxes (iauEect00).
inline double equation_of_equinoxes_complementary(const time::JulianDate& tt) {
  const double t = centuries(tt);
  const double fa[8] = {fundamental::l(t), fundamental::lp(t), fundamental::f(t), fundamental::d(t),
                        fundamental::om(t), fundamental::ve(t), fundamental::ea(t), fundamental::pa(t)};
  double s0 = 0, s1 = 0;
  for (int i = 32; i >= 0; i--) {
    double a = 0;
    for (int j = 0; j < 8; j++) a += series::eect00_e0[i].nfa[j] * fa[j];
    s0 += series::eect00_e0[i].s * std::sin(a) + series::eect00_e0[i].c * std::cos(a);
  }
  {
    double a = 0;
    for (int j = 0; j < 8; j++) a += series::eect00_e1[0].nfa[j] * fa[j];
    s1 += series::eect00_e1[0].s * std::sin(a) + series::eect00_e1[0].c * std::cos(a);
  }
  return (s0 + s1 * t) * ARCSEC;
}
// The equation of the equinoxes, Δψ cos ε_A plus the complementary terms (iauEe00).
inline double equation_of_equinoxes(const time::JulianDate& tt, const Nutation& n) {
  return n.dpsi * std::cos(mean_obliquity(tt)) + equation_of_equinoxes_complementary(tt);
}
// Greenwich apparent sidereal time by the CIO route, θ − EO (iauGst06a).
inline double gast(const time::JulianDate& ut1, const time::JulianDate& tt) {
  const Mat3 npb = bias_precession_nutation(tt);
  const auto xy = cip_xy(npb);
  return normalised_angle(earth_rotation_angle(ut1) - equation_of_origins(npb, cio_locator(tt, xy[0], xy[1])));
}
// Greenwich apparent sidereal time by the equinox route, GMST + EE: independent of the CIO
// series, it agrees with `gast` to about a microarcsecond.
inline double gast_classical(const time::JulianDate& ut1, const time::JulianDate& tt) {
  return normalised_angle(gmst(ut1, tt) + equation_of_equinoxes(tt, nutation(tt)));
}

// ---- polar motion -----------------------------------------------------------------------------

// The TIO locator s′ (iauSp00): −47 µas per century.
inline double tio_locator(const time::JulianDate& tt) { return -47e-6 * centuries(tt) * ARCSEC; }
// The polar motion matrix W, TIRS to ITRS (iauPom00): R_x(−y_p) R_y(−x_p) R_z(s′), radians.
inline Mat3 polar_motion(double xp, double yp, double sp) {
  return frame_rotation(Axis::x, -yp) * frame_rotation(Axis::y, -xp) * frame_rotation(Axis::z, sp);
}

// ---- GCRS to ITRS ---------------------------------------------------------------------------------

// The intermediate quantities of the celestial-to-terrestrial transformation at an instant.
struct EarthOrientation {
  Mat3 npb;      // GCRS → true equator and equinox of date
  Mat3 c2i;      // GCRS → CIRS
  double x, y, s;  // the CIP and the CIO locator, radians (X, Y with the observed offsets dX, dY)
  double era;    // Earth rotation angle, radians
  double gast;   // apparent sidereal time, CIO route
  Mat3 w;        // TIRS → ITRS
  Mat3 c2t;      // GCRS → ITRS, CIO route
};
// From TT, UT1 and the IERS's x_p, y_p, dX, dY (arcseconds).
inline EarthOrientation orientation_of_earth(const time::JulianDate& tt, const time::JulianDate& ut1, double xp_arcsec = 0,
                                          double yp_arcsec = 0, double dx_arcsec = 0, double dy_arcsec = 0) {
  EarthOrientation e;
  e.npb = bias_precession_nutation(tt);
  const auto xy = cip_xy(e.npb);
  e.x = xy[0] + dx_arcsec * ARCSEC;
  e.y = xy[1] + dy_arcsec * ARCSEC;
  e.s = cio_locator(tt, e.x, e.y);
  e.c2i = celestial_to_intermediate(e.x, e.y, e.s);
  e.era = earth_rotation_angle(ut1);
  e.gast = normalised_angle(e.era - equation_of_origins(e.npb, e.s));
  e.w = polar_motion(xp_arcsec * ARCSEC, yp_arcsec * ARCSEC, tio_locator(tt));
  e.c2t = e.w * frame_rotation(Axis::z, e.era) * e.c2i;
  return e;
}

// The derivative of R_z(θ) (a frame rotation) with respect to θ.
inline Mat3 frame_rotation_z_derivative(double theta) {
  const double c = std::cos(theta), s = std::sin(theta);
  return {{{-s, c, 0}, {-c, -s, 0}, {0, 0, 0}}};
}

}  // namespace metrology::frames
