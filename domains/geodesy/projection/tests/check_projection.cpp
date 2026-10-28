// Checks of the map projections and grids (geodesy/projection): the transverse Mercator, UTM, UPS
// and MGRS on WGS 84.
//
// References: PROJ 9.8 through pyproj 3.8 for UTM and UPS coordinates, convergence and scale, and
// NGA's GEOTRANS through the mgrs package for MGRS references (tests/fixtures/projection-reference.txt,
// made by tests/make_projection_fixtures.py); closed forms for the projections' properties (the
// central meridian at scale k₀, the equator at northing 0, the pole of UPS at 2000 km); round trips.

#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>

#include <calx/verification/check.hpp>
#include <geodesy/projection/mgrs.hpp>
#include <geodesy/projection/transverse_mercator.hpp>
#include <geodesy/projection/utm.hpp>

using namespace geodesy::projection;
using calx::verification::Suite;

int main() {
  Suite s;
  const UniversalGrids grids;
  const Mgrs mgrs(grids);
  const TransverseMercator& tm = grids.transverse_mercator();

  // ---- properties in closed form ----------------------------------------------------------------------
  {
    const GridPoint p = tm.forward(0.7, 0);
    s.check("transverse Mercator: on the central meridian the scale is k0", p.scale - UTM_K0, 1e-15);
    s.check("transverse Mercator: on the central meridian the convergence vanishes", p.convergence, 1e-16, "rad");
    s.check("transverse Mercator: the equator maps to northing 0", tm.forward(0, 0.05).northing, 1e-9, "m");
    // The quarter meridian of WGS 84 is 10 001 965.7293 m (Karney 2011, eq. 14 with A = 6 367 449.1458 m).
    s.check("transverse Mercator: the pole's northing is k0 times the quarter meridian, 10 001 965.7293 m",
            tm.forward(std::numbers::pi / 2, 0).northing / UTM_K0 - 10001965.7293, 1e-4, "m");
    s.check("transverse Mercator: the rectifying radius A = 6 367 449.1458 m", tm.rectifying_radius() - 6367449.1458, 1e-4, "m");
  }

  // ---- against PROJ and GEOTRANS --------------------------------------------------------------------------
  std::ifstream in(GEODESY_FIXTURES "/projection-reference.txt");
  std::string line;
  int n_utm = 0, n_ups = 0, n_mgrs = 0, mgrs_wrong = 0;
  double utm_xy = 0, utm_conv = 0, utm_scale = 0, utm_back = 0, ups_xy = 0, ups_back = 0, mgrs_back = 0;
  std::string first_wrong;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    std::string kind;
    f >> kind;
    if (kind == "utm") {
      int zone, north;
      double lat, lon, e, n, conv, scale;
      f >> zone >> north >> lat >> lon >> e >> n >> conv >> scale;
      const UtmUps p = grids.forward(lat * DEG, lon * DEG, zone);
      utm_xy = std::max(utm_xy, std::hypot(p.easting - e, p.northing - n));
      utm_conv = std::max(utm_conv, std::fabs(p.convergence - conv * DEG));
      utm_scale = std::max(utm_scale, std::fabs(p.scale - scale));
      const auto g = grids.reverse(zone, north, e, n);
      utm_back = std::max(utm_back, std::hypot(g.latitude - lat * DEG, std::remainder(g.longitude - lon * DEG, 2 * std::numbers::pi) * std::cos(lat * DEG)) * 6.4e6);
      n_utm++;
    } else if (kind == "ups") {
      int north;
      double lat, lon, e, n;
      f >> north >> lat >> lon >> e >> n;
      const UtmUps p = grids.forward(lat * DEG, lon * DEG, 0);
      ups_xy = std::max(ups_xy, std::hypot(p.easting - e, p.northing - n));
      const auto g = grids.reverse(0, north, e, n);
      const double dlon = std::fabs(lat) > 89.99999 ? 0 : std::remainder(g.longitude - lon * DEG, 2 * std::numbers::pi) * std::cos(lat * DEG);
      ups_back = std::max(ups_back, std::hypot(g.latitude - lat * DEG, dlon) * 6.4e6);
      n_ups++;
    } else if (kind == "mgrs") {
      double lat, lon;
      std::string ref;
      f >> lat >> lon >> ref;
      const std::string got = mgrs.forward(lat * DEG, lon * DEG, 5);
      if (got != ref) {
        mgrs_wrong++;
        if (first_wrong.empty()) first_wrong = line + " → " + got;
      }
      // The reference names a 1 m square: its south-west corner lies within a metre and a half of the point.
      const auto back = mgrs.reverse(ref, true);
      const double d = std::hypot(back.geodetic.latitude - lat * DEG,
                                  std::remainder(back.geodetic.longitude - lon * DEG, 2 * std::numbers::pi) * std::cos(lat * DEG)) * 6.4e6;
      mgrs_back = std::max(mgrs_back, d);
      n_mgrs++;
    }
  }
  s.require("the reference file was read", n_utm > 300 && n_ups > 100 && n_mgrs > 400);
  s.check("UTM: easting and northing against PROJ, out to 4° from the central meridian", utm_xy, 1e-6, "m");
  s.check("UTM: the meridian convergence against PROJ (its factors are numerical, to about 1e-10)", utm_conv, 1e-10, "rad");
  s.check("UTM: the point scale against PROJ (its factors are numerical, to about 1e-10)", utm_scale, 3e-10);
  s.check("UTM: grid to geodetic recovers PROJ's points", utm_back, 1e-6, "m");
  s.check("UPS: easting and northing against PROJ", ups_xy, 1e-6, "m");
  s.check("UPS: grid to geodetic recovers PROJ's points", ups_back, 1e-6, "m");
  s.require("MGRS: every reference as GEOTRANS writes it (with UPS, Norway and Svalbard) — " + (first_wrong.empty() ? std::string("all equal") : first_wrong),
            mgrs_wrong == 0);
  s.check("MGRS: a 1 m reference's corner lies within 1.5 m of its point", mgrs_back, 1.5, "m");

  // ---- convergence and scale by differences of the projection itself ------------------------------------------
  // PROJ derives its factors by numerical differentiation (to about 10⁻¹⁰); the analytic ones here are
  // held to it above at that level, and to central differences of this projection along the meridian:
  // the scale is the grid length over the ellipsoid's, M Δφ with M = a(1 − e²)/(1 − e² sin² φ)^{3/2};
  // the meridian's image runs at the bearing −γ from grid north.
  {
    const auto& e = grids.ellipsoid();
    double worst_k = 0, worst_g = 0;
    for (const auto& [phi, lam] : {std::pair{0.3, 0.05}, {-1.1, -0.06}, {1.3, 0.07}, {0.01, -0.02}, {-0.6, 0.069}}) {
      const double h = 1e-6;
      const GridPoint a = tm.forward(phi - h, lam), b = tm.forward(phi + h, lam), p = tm.forward(phi, lam);
      const double dE = b.easting - a.easting, dN = b.northing - a.northing;
      const double sp = std::sin(phi), M = e.a() * (1 - e.e2()) / std::pow(1 - e.e2() * sp * sp, 1.5);
      worst_k = std::max(worst_k, std::fabs(std::hypot(dE, dN) / (M * 2 * h) - p.scale));
      worst_g = std::max(worst_g, std::fabs(-std::atan2(dE, dN) - p.convergence));
    }
    s.check("transverse Mercator: the point scale against differences of the projection", worst_k, 1e-9);
    s.check("transverse Mercator: the convergence against differences of the projection", worst_g, 1e-9, "rad");
  }

  // ---- the zones and the bands -------------------------------------------------------------------------
  s.require("zones: 32 is widened over south-west Norway", standard_zone(60 * DEG, 5 * DEG) == 32 && standard_zone(60 * DEG, 2.9 * DEG) == 31);
  s.require("zones: Svalbard has only 31, 33, 35 and 37", standard_zone(78 * DEG, 8 * DEG) == 31 && standard_zone(78 * DEG, 10 * DEG) == 33 &&
                                                           standard_zone(78 * DEG, 30 * DEG) == 35 && standard_zone(78 * DEG, 40 * DEG) == 37);
  s.require("zones: UPS beyond 84° N and 80° S", standard_zone(84.1 * DEG, 0) == 0 && standard_zone(-80.1 * DEG, 0) == 0 && standard_zone(83.9 * DEG, 0) == 31);
  s.require("bands: C at 80° S, X from 72° to 84° N", utm_band(-79.9 * DEG) == 'C' && utm_band(72.1 * DEG) == 'X' && utm_band(83.9 * DEG) == 'X' && utm_band(0.1 * DEG) == 'N');

  // ---- MGRS at every precision, and its round trip ---------------------------------------------------------
  {
    const double lat = 48.8566 * DEG, lon = 2.3522 * DEG;
    const std::string r5 = mgrs.forward(lat, lon, 5);
    bool truncates = true;
    for (int d = 0; d <= 4; d++) {
      const std::string r = mgrs.forward(lat, lon, d);
      truncates = truncates && r == r5.substr(0, 5) + r5.substr(5, d) + r5.substr(10, d);
    }
    s.require("MGRS: a coarser reference is the finer one truncated (" + r5 + ")", truncates);
    double worst = 0;
    for (int d = 1; d <= 5; d++) {
      const auto back = mgrs.reverse(mgrs.forward(lat, lon, d));
      const double unit = std::pow(10.0, 5 - d);
      worst = std::max(worst, std::hypot(back.geodetic.latitude - lat, (back.geodetic.longitude - lon) * std::cos(lat)) * 6.4e6 / unit);
    }
    s.check("MGRS: the centre of a reference's square is within its half diagonal of the point (in squares)", worst, 0.71);
    bool refused = false;
    try { mgrs.reverse("31UDQ123"); } catch (const std::invalid_argument&) { refused = true; }
    s.require("MGRS: an odd number of digits is refused", refused);
    const auto spaced = mgrs.reverse("31u dq 52482 11717");
    s.check("MGRS: spaces and small letters are read", std::hypot(spaced.geodetic.latitude - lat, (spaced.geodetic.longitude - lon) * std::cos(lat)) * 6.4e6, 1.5, "m");
  }
  return s.done();
}
