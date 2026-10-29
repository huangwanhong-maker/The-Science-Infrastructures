// Checks of the normal gravity of a level ellipsoid (geodesy/gravity/normal.hpp).
//
// References (tests/fixtures/gravity-reference.txt, made by tests/make_gravity_fixtures.py): the derived
// constants as NGA.STND.0036_1.0.0_WGS84 (Tables 3.5, 3.6) and Moritz's Geodetic Reference System 1980
// (sec. 4) print them, each held to half a unit of its last printed digit; boule's level ellipsoid, its
// normal gravity (Li and Goetze 2001) and potential at 666 points to 100 km, and the derivatives of its
// potential along the ellipsoidal and the geodetic coordinates. And by construction: Somigliana's
// formula against the closed field on the ellipsoid, the ellipsoid as a level surface (U = U0), Clairaut's
// theorem, the vector as the gradient of the potential, NGA's bounds on the direction of normal gravity
// below 20 km, and the second-order series in height.

#include <cmath>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <geodesy/gravity/field.hpp>
#include <geodesy/gravity/normal.hpp>
#include <physis/fields/analytic.hpp>

using namespace geodesy::gravity;
using calx::verification::relative;
using calx::verification::Suite;

namespace {

struct Printed {
  double value, half_unit;
};

// The level ellipsoid's value of a published constant named as the fixture names it.
double constant(const LevelEllipsoid& e, const std::string& name) {
  if (name == "f") return e.flattening();
  if (name == "inverse_f") return 1 / e.flattening();
  if (name == "b") return e.b();
  if (name == "e2") return e.e2();
  if (name == "second_e") return e.second_eccentricity();
  if (name == "second_e2") return e.second_e2();
  if (name == "E") return e.linear_eccentricity();
  if (name == "c") return e.polar_radius_of_curvature();
  if (name == "J2") return e.J2();
  if (name == "C20") return -e.J2() / std::sqrt(5.0);  // NGA eq. 3-1
  if (name == "J4") return e.J2n(2);
  if (name == "J6") return e.J2n(3);
  if (name == "J8") return e.J2n(4);
  if (name == "U0") return e.U0();
  if (name == "m") return e.m();
  if (name == "gamma_e") return e.gamma_equator();
  if (name == "gamma_p") return e.gamma_pole();
  if (name == "gamma_mean") return e.mean_gravity();
  if (name == "gamma_45") return e.somigliana(std::numbers::pi / 4);
  if (name == "k") return e.k();
  if (name == "f_star") return e.gravity_flattening();
  throw std::invalid_argument("unknown constant " + name);
}

}  // namespace

int main() {
  Suite s;
  const double deg = std::numbers::pi / 180;
  const LevelEllipsoid wgs84 = LevelEllipsoid::wgs84(), grs80 = LevelEllipsoid::grs80();
  // boule builds GRS 80 from Moritz's printed 1/f; the comparisons with boule use the same.
  const LevelEllipsoid grs80_printed_f = LevelEllipsoid::from_flattening("GRS 80 (printed 1/f)", 6378137.0, 298.257222101, 3986005e8, 7292115e-11);

  // ---- read the references ------------------------------------------------------------------------------
  std::vector<std::tuple<std::string, std::string, Printed>> published;
  std::map<std::pair<std::string, std::string>, double> boule;
  std::vector<std::pair<std::string, std::vector<double>>> field;
  {
    std::ifstream in(GEODESY_FIXTURES "/gravity-reference.txt");
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      std::string kind, ellipsoid, name;
      f >> kind >> ellipsoid;
      if (kind == "published") {
        Printed p;
        f >> name >> p.value >> p.half_unit;
        published.emplace_back(ellipsoid, name, p);
      } else if (kind == "boule") {
        double v;
        f >> name >> v;
        boule[{ellipsoid, name}] = v;
      } else {
        std::vector<double> row;
        double v;
        while (f >> v) row.push_back(v);
        field.emplace_back(ellipsoid, row);
      }
    }
  }
  s.require("the references were read (33 published constants, 16 of boule's, 666 points)", published.size() == 33 && boule.size() == 16 && field.size() == 666);

  // ---- the defining constants --------------------------------------------------------------------------------
  s.check("GRS 80 from J2: J2 comes back as 108 263e-8 (Moritz's iteration for e2), to rounding through 1/f", grs80.J2() - 108263e-8, 2e-17);
  s.check("WGS 84 from f: 1/f comes back as 298.257 223 563", 1 / wgs84.flattening() - 298.257223563, 1e-9);
  {
    bool threw = false;
    try { LevelEllipsoid::from_J2("no", 6378137.0, -1e-3, 3986005e8, 7292115e-11); } catch (const std::invalid_argument&) { threw = true; }
    s.require("a level ellipsoid with J2 <= 0 is refused", threw);
  }

  // ---- the derived constants as the documents print them ------------------------------------------------------
  for (const auto& [ellipsoid, name, p] : published) {
    const LevelEllipsoid& e = ellipsoid == "WGS84" ? wgs84 : grs80;
    const std::string source = ellipsoid == "WGS84" ? "WGS 84 (NGA.STND.0036, Tables 3.5-3.6)" : "GRS 80 (Moritz, sec. 4)";
    std::ostringstream label;
    label.precision(15);
    label << source << ": " << name << " = " << p.value << ", in half-units of its last printed digit";
    s.check(label.str(), std::fabs(constant(e, name) - p.value) / p.half_unit, 1.0);
  }

  // ---- q and q': the series against the closed forms in long double ---------------------------------------------
  {
    // The closed forms cancel (q ~ 2x^3/15 from terms near 3/x); in long double, with an 64-bit significand,
    // they keep twelve digits or more for x >= 0.06, which the series must meet. In double, as boule sums
    // them, they keep some 1e-11 at e' = 0.082, which is why boule's gamma_e and gamma_p differ below.
    double worst_q = 0, worst_qp = 0;
    for (double x = 0.06; x <= 0.49; x += 0.0037) {
      const long double X = x, a = std::atan(X);
      const long double q = 0.5L * ((1 + 3 / (X * X)) * a - 3 / X), qp = 3 * (1 + 1 / (X * X)) * (1 - a / X) - 1;
      worst_q = std::max(worst_q, std::fabs(static_cast<double>((detail::q_of(x) - q) / q)));
      worst_qp = std::max(worst_qp, std::fabs(static_cast<double>((detail::q_prime_of(x) - qp) / qp)));
    }
    s.check("q(x) by its series against the closed form in long double, x from 0.06 to 0.49 (relative)", worst_q, 1e-12);
    s.check("q'(x) by its series against the closed form in long double (relative)", worst_qp, 1e-12);
    s.check("q and q' are continuous where the series hands over to the closed form (x = 1/2)",
            std::max(relative(detail::q_of(0.5 - 1e-15), detail::q_of(0.5)), relative(detail::q_prime_of(0.5 - 1e-15), detail::q_prime_of(0.5))), 1e-12);
  }

  // ---- against boule's level ellipsoid ----------------------------------------------------------------------------
  for (const std::string ellipsoid : {"WGS84", "GRS80"}) {
    const LevelEllipsoid& e = ellipsoid == "WGS84" ? wgs84 : grs80_printed_f;
    double geometric = 0, gravity = 0;
    for (const std::string name : {"b", "E", "e2", "second_e2", "m", "U0"}) geometric = std::max(geometric, relative(constant(e, name), boule.at({ellipsoid, name})));
    for (const std::string name : {"gamma_e", "gamma_p"}) gravity = std::max(gravity, relative(constant(e, name), boule.at({ellipsoid, name})));
    s.check(e.name() + ": b, E, e2, e'2, m and U0 against boule (relative)", geometric, 1e-14);
    s.check(e.name() + ": gamma_e and gamma_p against boule (relative; boule's q0 and q0' cancel in double)", gravity, 5e-13);
  }

  // ---- normal gravity off the ellipsoid: against boule, at 666 points to 100 km ------------------------------------
  {
    // field: lat lon h gamma(=-gamma_u) U u beta(deg) gamma_beta down north
    double gu = 0, U = 0, u = 0, beta = 0, gb = 0, down = 0, north = 0, mag = 0, east = 0;
    for (const auto& [ellipsoid, r] : field) {
      const LevelEllipsoid& e = ellipsoid == "WGS84" ? wgs84 : grs80_printed_f;
      const Geodetic p{r[0] * deg, r[1] * deg, r[2]};
      const NormalGravity g = e.at(p);
      const EllipsoidalHarmonic h = e.ellipsoidal_harmonic(e.figure().cartesian(p));
      gu = std::max(gu, std::fabs(g.u_m_s2 + r[3]));
      U = std::max(U, relative(e.potential(p), r[4]));
      u = std::max(u, std::fabs(h.u - r[5]));
      beta = std::max(beta, std::fabs(h.reduced_latitude / deg - r[6]));
      gb = std::max(gb, std::fabs(g.beta_m_s2 - r[7]));
      down = std::max(down, std::fabs(g.down_m_s2 - r[8]));
      north = std::max(north, std::fabs(g.north_m_s2 - r[9]));
      mag = std::max(mag, std::fabs(g.magnitude_m_s2 - std::hypot(r[3], r[7])));
      east = std::max(east, std::fabs(g.enu_m_s2[0]));
    }
    s.check("gamma_u (NGA eq. 4-5) against boule's normal gravity (Li and Goetze), 0-100 km", gu, 1e-11, "m/s2");
    s.check("the normal potential U against boule's (HWM eq. 2-126), relative", U, 1e-12);
    s.check("the ellipsoidal coordinate u (NGA eq. 4-8) against boule's", u, 1e-8, "m");
    s.check("the reduced latitude beta (NGA eq. 4-9) against boule's", beta, 1e-10, "deg");
    s.check("gamma_beta (NGA eq. 4-6) against the derivative of boule's U along beta", gb, 5e-11, "m/s2");
    s.check("|gamma| = sqrt(gamma_u^2 + gamma_beta^2) (NGA eq. 4-4) against boule's two components", mag, 1e-11, "m/s2");
    s.check("the component down the geodetic normal (eq. 4-16) against -dU/dh of boule's U", down, 5e-9, "m/s2");
    s.check("the component north in the meridian (eq. 4-23) against dU/dphi/(M+h) of boule's U", north, 5e-9, "m/s2");
    s.check("normal gravity has no east component", east, 1e-12, "m/s2");
  }

  // ---- on the ellipsoid: Somigliana, and the ellipsoid as a level surface ------------------------------------------
  for (const LevelEllipsoid* e : {&wgs84, &grs80}) {
    double som = 0, level = 0;
    for (double lat = -90; lat <= 90; lat += 2.5) {
      const Geodetic p{lat * deg, 0.37 * lat * deg, 0.0};
      som = std::max(som, std::fabs(e->somigliana(p.latitude) - e->at(p).magnitude_m_s2));
      level = std::max(level, relative(e->potential(p), e->U0()));
    }
    s.check(e->name() + ": Somigliana's formula (NGA eq. 4-1) is the closed field at h = 0, every 2.5 deg", som, 1e-13, "m/s2");
    s.check(e->name() + ": the ellipsoid is a level surface, U = U0 on it (relative)", level, 1e-15);
    s.check(e->name() + ": gamma_e and gamma_p are the field at the equator and the pole",
            std::max(std::fabs(e->at({0, 0, 0}).magnitude_m_s2 - e->gamma_equator()), std::fabs(e->at({std::numbers::pi / 2, 0, 0}).magnitude_m_s2 - e->gamma_pole())), 1e-13, "m/s2");
    // Clairaut's theorem in closed form, f + f* = (omega^2 b/gamma_e)(1 + e' q0'/(2 q0)) (HM 2-75), an identity
    // of the formulas for gamma_e and gamma_p; f* is a difference of numbers near 1, so it holds to a few ulps of 1.
    const double clairaut = e->omega() * e->omega() * e->b() / e->gamma_equator() * (1 + e->second_eccentricity() * e->q0_prime() / (2 * e->q0()));
    s.check(e->name() + ": Clairaut's theorem f + f* = (omega^2 b/gamma_e)(1 + e' q0'/2q0) (HM 2-75)", e->flattening() + e->gravity_flattening() - clairaut, 1e-15);
  }

  // ---- the vector: the gradient of the potential, and its direction ----------------------------------------------
  {
    // Central differences of this module's own potential along the geodetic height and latitude.
    double grad = 0;
    const double d = 200;
    for (double lat = -88; lat <= 88; lat += 8)
      for (double h : {0.0, 5000.0, 40000.0, 100000.0}) {
        const Geodetic p{lat * deg, 1.1 * lat * deg, h};
        const NormalGravity g = wgs84.at(p);
        auto U = [&](double dphi, double dh) { return wgs84.potential({p.latitude + dphi, p.longitude, p.height + dh}); };
        const double M = wgs84.figure().meridian_radius(p.latitude) + h, dphi = d / M;
        const double dUdh = (-U(0, 2 * d) + 8 * U(0, d) - 8 * U(0, -d) + U(0, -2 * d)) / (12 * d);
        const double dUdphi = (-U(2 * dphi, 0) + 8 * U(dphi, 0) - 8 * U(-dphi, 0) + U(-2 * dphi, 0)) / (12 * dphi);
        grad = std::max({grad, std::fabs(g.down_m_s2 + dUdh), std::fabs(g.north_m_s2 - dUdphi / M)});
      }
    s.check("WGS 84: the vector is the gradient of the potential (differences of U along h and phi, to 100 km)", grad, 1e-9, "m/s2");
    // NGA sec. 4.3: below 20 km, |epsilon| < 4", |gamma_phi| < 0.002 % of gamma_h, and eq. 4-4 within 1 uGal of eq. 4-16.
    double eps = 0, ratio = 0, approx = 0;
    for (double lat = -90; lat <= 90; lat += 1)
      for (double h : {1000.0, 5000.0, 10000.0, 20000.0}) {
        const NormalGravity g = wgs84.at({lat * deg, 0.0, h});
        eps = std::max(eps, std::fabs(g.deflection_rad) / deg * 3600);
        ratio = std::max(ratio, std::fabs(g.north_m_s2) / g.down_m_s2);
        approx = std::max(approx, g.magnitude_m_s2 - g.down_m_s2);
      }
    s.check("below 20 km the direction departs from the geodetic normal by less than 4\" (NGA sec. 4.3)", eps, 4.0, "arcsec");
    s.check("below 20 km gamma_phi is less than 0.002 % of gamma_h (NGA sec. 4.3)", ratio, 2e-5);
    s.check("below 20 km |gamma| (eq. 4-4) exceeds gamma_h (eq. 4-16) by less than 1 uGal (NGA sec. 4.3)", approx, 1e-8, "m/s2");
    const NormalGravity n = wgs84.at({45 * deg, 0.0, 10000.0});
    s.require("at 45 deg N, 10 km up, normal gravity points down and leans south", n.down_m_s2 > 9.7 && n.north_m_s2 < 0 && n.deflection_rad < 0);
  }

  // ---- the second-order series in height (NGA eq. 4-3) --------------------------------------------------------------
  {
    auto worst = [&](double h) {
      double w = 0;
      for (double lat = -90; lat <= 90; lat += 1) w = std::max(w, std::fabs(wgs84.series(lat * deg, h) - wgs84.at({lat * deg, 0.0, h}).down_m_s2));
      return w;
    };
    s.check("the series in height (eq. 4-3) against the closed field at 1 km", worst(1000), 1.5e-7, "m/s2");
    s.check("and at 20 km (0.15 mGal)", worst(20000), 1.5e-6, "m/s2");
    s.check("and at 100 km, where it is no longer for use (20 mGal)", worst(100000), 2e-4, "m/s2");
    s.check("the series is exact on the ellipsoid", worst(0), 1e-13, "m/s2");
  }

  // ---- the bounds of the closed formulas ------------------------------------------------------------------------------
  {
    bool threw = false;
    try { wgs84.at({0.0, 0.0, -wgs84.a()}); } catch (const std::domain_error&) { threw = true; }
    s.require("at the centre, on the focal disc, normal gravity is refused", threw);
    const NormalGravity below = wgs84.at({30 * deg, 0.0, -500.0}), above = wgs84.at({30 * deg, 0.0, 500.0}), on = wgs84.at({30 * deg, 0.0, 0.0});
    // Through the surface the field continues smoothly: its second difference over +-500 m is eq. 4-3's
    // 6 gamma d^2/a^2 (3.6e-7 m s-2), to the O(f) the series neglects.
    const double second = above.down_m_s2 - 2 * on.down_m_s2 + below.down_m_s2, predicted = 6 * on.down_m_s2 * 500.0 * 500.0 / (wgs84.a() * wgs84.a());
    s.check("below the ellipsoid the field continues it smoothly (second difference over +-500 m against 6 gamma d^2/a^2)", second - predicted, 1e-8, "m/s2");
  }

  // ---- normal gravity as a field (physis fields) -----------------------------------------------------------------------
  {
    namespace pf = physis::fields;
    const NormalGravityField g(wgs84);
    const pf::Descriptor d = pf::described(g);
    s.require("field: gamma on the Earth-fixed axes, a vector in m/s2, analytic, of the model, static",
              d.name == "earth.gravity.normal" && d.frame == "ITRS" && d.time == "none" && d.fibre.dimension == physis::units::dim::acceleration &&
                  d.representation.kind == pf::RepresentationKind::analytic && d.provenance == pf::Provenance::model);
    // The field at a point is the closed field at its geodetic place, found again on the ellipsoid.
    double same = 0;
    for (double lat = -88; lat <= 88; lat += 8)
      for (double h : {0.0, 5000.0, 100000.0, 1e6}) {
        const Geodetic p{lat * deg, 0.7 * lat * deg, h};
        const Vec3 a = g.at({0, wgs84.figure().cartesian(p)}).value, b = wgs84.at(p).ecef_m_s2;
        for (int i = 0; i < 3; i++) same = std::max(same, std::fabs(a[i] - b[i]));
      }
    s.check("field: at 92 places to 1000 km, the closed field at the place", same, 1e-12, "m/s2");
    const auto below = g.at({0, wgs84.figure().cartesian({0.5, 1.0, -5000.0})}), deep = g.at({0, wgs84.figure().cartesian({0.5, 1.0, -12000.0})});
    s.require("field: 5 km below the ellipsoid, the continuation, with the note; 12 km below, outside, with the reason",
              below.ok() && below.note == "below the ellipsoid: its analytic continuation, as free-air reduction takes it" &&
                  deep.quality == pf::Quality::outside && deep.reason == "more than 11 km below the ellipsoid, beneath the field's continuation");
    // A gradient has no curl; outside the masses its divergence is the centrifugal potential's
    // Laplacian, 2 omega^2 = 1.0635e-8 s-2. By central differences of step h = 1e-4 r each derivative errs
    // by h^2/6 times the third derivative, 24 GM/r^5 for GM/r^2: the divergence by about
    // 12 (h/r)^2 gamma/r, 2e-13 s-2 at the ground.
    calx::random::Stream rnd(20261029, 5, 0, 0);
    double div = 0, cur = 0;
    for (int i = 0; i < 20; i++) {
      const Vec3 x = wgs84.figure().cartesian({std::asin(rnd.uniform(-0.99, 0.99)), rnd.uniform(-std::numbers::pi, std::numbers::pi), rnd.uniform(0, 1e6)});
      const double r = std::hypot(x[0], x[1], x[2]);
      const auto D = pf::central_differences(g, {0, x}, 1e-4 * r);
      div = std::max(div, std::fabs(pf::divergence(D.value) - 2 * wgs84.omega() * wgs84.omega()));
      cur = std::max(cur, std::hypot(pf::curl(D.value)[0], pf::curl(D.value)[1], pf::curl(D.value)[2]));
    }
    s.check("field: its divergence is 2 omega^2, by differences of step 1e-4 r (20 points to 1000 km)", div, 1e-12, "1/s2");
    s.check("field: its curl vanishes", cur, 1e-12, "1/s2");
  }
  return s.done();
}
