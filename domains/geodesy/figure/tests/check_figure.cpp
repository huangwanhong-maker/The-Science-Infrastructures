// Checks of the figure of the Earth: geodetic and Earth-centred coordinates on WGS 84 and GRS 80,
// and local level frames (geodesy/figure).
//
// References: ERFA's gd2gc and gc2gd through pyerfa (tests/fixtures/geodetic-reference.txt), at
// points from the poles to the equator and from 10 km below the ellipsoid to geostationary
// height; the defining constants of WGS 84 (NIMA TR8350.2) for its derived ones; closed forms
// for the local frames.

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>

using namespace geodesy::figure;
using calx::verification::Suite;

namespace {
constexpr double pi = std::numbers::pi;
double max_abs(const Vec3& v) { return std::max({std::fabs(v[0]), std::fabs(v[1]), std::fabs(v[2])}); }
}  // namespace

int main() {
  Suite s;
  const Ellipsoid wgs = Ellipsoid::wgs84(), grs = Ellipsoid::grs80();

  // ---- the derived constants of WGS 84 (NIMA TR8350.2, Table 3.3) ----------------------------------
  s.check("WGS 84: b = 6 356 752.3142 m", wgs.b() - 6356752.3142, 1e-4, "m");
  s.check("WGS 84: e² = 6.69437999014 × 10⁻³", wgs.e2() - 6.69437999014e-3, 1e-14);
  s.check("WGS 84: e′² = 6.73949674228 × 10⁻³", wgs.second_e2() - 6.73949674228e-3, 1e-14);
  s.check("WGS 84: R₁ = 6 371 008.7714 m", wgs.mean_radius() - 6371008.7714, 1e-4, "m");

  // ---- against ERFA ---------------------------------------------------------------------------------
  // Each point was given geodetically, carried to Earth-centred coordinates by ERFA's gd2gc, and
  // back by ERFA's gc2gd. The inverse here is held to the original point, which it recovers to
  // rounding; ERFA's own inverse (Fukushima 2006) drifts by up to 2.4 × 10⁻¹¹ rad (0.15 mm at the
  // surface, 1 mm at geostationary height) at satellite heights, and is held to within that.
  {
    std::ifstream in(GEODESY_FIXTURES "/geodetic-reference.txt");
    std::string line;
    double forward = 0, latitude = 0, height = 0, longitude = 0, erfa_latitude = 0;
    int rows = 0;
    Geodetic original{};
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      std::string name;
      int n;
      double v[6];
      f >> name >> n >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5];
      const Ellipsoid& e = n == 1 ? wgs : grs;
      rows++;
      if (name == "gd2gc") {
        original = {v[1], v[0], v[2]};
        const Vec3 x = e.cartesian(original);
        forward = std::max(forward, max_abs(x - Vec3{v[3], v[4], v[5]}));
      } else {
        const Geodetic g = e.geodetic({v[0], v[1], v[2]});
        latitude = std::max(latitude, std::fabs(g.latitude - original.latitude));
        if (std::fabs(original.latitude) < pi / 2 - 1e-8) longitude = std::max(longitude, std::fabs(std::remainder(g.longitude - original.longitude, 2 * pi)));
        height = std::max(height, std::fabs(g.height - original.height) / std::max(1.0, std::fabs(original.height) * 1e-7));
        erfa_latitude = std::max(erfa_latitude, std::fabs(g.latitude - v[4]));
      }
    }
    s.require("the reference file was read", rows > 100);
    s.check("geodetic → Earth-centred against ERFA (98 points, WGS 84 and GRS 80)", forward, 1e-8, "m");
    s.check("Earth-centred → geodetic recovers the latitude", latitude, 1e-15, "rad");
    s.check("Earth-centred → geodetic recovers the longitude", longitude, 1e-15, "rad");
    s.check("Earth-centred → geodetic recovers the height (to 10⁻⁸ m, or 10⁻¹⁵ of it beyond 100 km)", height, 1e-8, "m");
    s.check("ERFA's inverse agrees within its own error at satellite heights", erfa_latitude, 3e-11, "rad");
  }

  // ---- round trips, near the centre and the poles -------------------------------------------------
  {
    double worst = 0;
    for (const double lat : {-pi / 2, -1.2, -0.3, 0.0, 1e-12, 0.6, 1.5, pi / 2})
      for (const double h : {-6000e3, -50e3, -1.0, 0.0, 1.0, 8848.86, 400e3, 35786e3, 384400e3}) {
        const Geodetic g{lat, 0.7, h};
        const Geodetic back = wgs.geodetic(wgs.cartesian(g));
        worst = std::max(worst, physis::geometry::norm(wgs.cartesian(back) - wgs.cartesian(g)));
      }
    s.check("geodetic ↔ Earth-centred round trips from 6000 km deep to the Moon's distance", worst, 1e-7, "m");
    const Geodetic pole = wgs.geodetic({0, 0, 6356752.314245179 + 10});
    s.check("10 m above the north pole", std::fabs(pole.latitude - pi / 2) + std::fabs(pole.height - 10), 1e-9);
    const Geodetic centre = wgs.geodetic({1, 0, 0});
    s.require("the centre of the Earth has a finite geodetic height", std::isfinite(centre.height));
  }

  // ---- local level frames -------------------------------------------------------------------------
  {
    const Place here = Place::at(wgs, {0.6, 2.0, 120});
    const Mat3 r = here.to_enu;
    s.check("ECEF → ENU is a rotation", physis::geometry::orthogonality_defect(r) + std::fabs(physis::geometry::determinant(r) - 1), 1e-15);
    s.check("up is the ellipsoid's normal", max_abs(Vec3{r[2][0], r[2][1], r[2][2]} - Ellipsoid::normal(0.6, 2.0)), 1e-16);
    // A point straight up is at elevation 90°; a point 1 km due east on the tangent plane at azimuth 90°.
    const LookAngles up = here.look_at(here.ecef_of({0, 0, 1000}));
    s.check("straight up: elevation 90°", up.elevation - pi / 2, 1e-12, "rad");
    const LookAngles east = here.look_at(here.ecef_of({1000, 0, 0}));
    s.check("due east: azimuth 90°", east.azimuth - pi / 2, 1e-12, "rad");
    s.check("due east: range 1 km", east.range - 1000, 1e-9, "m");
    // A point on the ellipsoid 10 km north along the meridian lies below the horizon by about
    // d / 2R (the Earth's curvature): −0.045°.
    const Geodetic north_point{0.6 + 10000 / wgs.meridian_radius(0.6), 2.0, 120};
    const LookAngles n = here.look_at(wgs.cartesian(north_point));
    s.check("10 km north along the meridian: azimuth 0°", std::remainder(n.azimuth, 2 * pi), 1e-9, "rad");
    s.check("10 km north: below the horizon by d / 2R", n.elevation + 10000 / (2 * wgs.meridian_radius(0.6)), 2e-6, "rad");
    const auto t = here.ned_to_earth();
    s.check("NED's x is north: the same as ENU's y", max_abs(t.direction({1, 0, 0}) - here.enu_to_earth().direction({0, 1, 0})), 1e-16);
    s.check("NED's z is down", max_abs(t.direction({0, 0, 1}) + here.enu_to_earth().direction({0, 0, 1})), 1e-16);
    s.check("look angles and back", max_abs(enu_from_look(look_angles({300, -400, 50})) - Vec3{300, -400, 50}), 1e-12, "m");
  }
  return s.done();
}
