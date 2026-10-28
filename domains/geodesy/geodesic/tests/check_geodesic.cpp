// Checks of geodesics on the ellipsoid (geodesy/geodesic): the inverse and the direct problem, points
// along a line, the sphere, and the hard cases of GeographicLib's own tests.
//
// References: GeographicLib 2.1 (Karney; its Python package, of which this is a port) through
// tests/fixtures/geodesic-reference.txt (made by tests/make_geodesic_fixtures.py): 620 inverse
// problems on WGS 84, anywhere, nearly antipodal, at the poles, along the equator and meridians,
// short and very short, and 400 direct problems from 1 mm to 40 000 km; GeographicLib's documented
// test cases with their values (test_geodesic.py, the 20 lines of its "testcases" and the GeodSolve
// examples quoted below with the tolerances it states); and on a sphere the closed forms of great
// circles. The fixtures and the port differ only by the libraries' rounding (CPython's on Windows,
// glibc's here); the limits say how far that carries.

#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <geodesy/geodesic/geodesic.hpp>

using namespace geodesy::geodesic;
using calx::verification::Suite;

namespace {
constexpr double pi = std::numbers::pi;
constexpr double D = pi / 180;  // radians per degree
// The difference of two angles in degrees, reduced to [−180, 180].
double dang(double a, double b) { return std::remainder(a - b, 360.0); }
double nanv() { return std::numeric_limits<double>::quiet_NaN(); }
}  // namespace

int main() {
  Suite s;
  const Geodesic g = Geodesic::wgs84();

  // ---- against GeographicLib: the fixtures ----------------------------------------------------------
  {
    std::ifstream in(GEODESY_FIXTURES "/geodesic-reference.txt");
    std::string line;
    int ninv = 0, ndir = 0, ndoc = 0;
    double inv_s = 0, inv_azi = 0, inv_a = 0, inv_m = 0, inv_M = 0, inv_S = 0;
    double dir_pos = 0, dir_azi = 0, dir_a = 0, dir_m = 0, dir_M = 0, dir_S = 0;
    double rad_s = 0, rad_azi = 0, round_trip = 0, round_trip_ill = 0, inv_shift = 0, rad_shift = 0;
    double doc_azi = 0, doc_s = 0, doc_a = 0, doc_m = 0, doc_M = 0, doc_S = 0;
    double doc_pos = 0, doc_dazi = 0, doc_da = 0, doc_dm = 0, doc_dM = 0, doc_dS = 0;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      std::string type;
      f >> type;
      if (type == "doc") {
        double v[12];
        for (double& x : v) f >> x;
        const double lat1 = v[0], lon1 = v[1], azi1 = v[2], lat2 = v[3], lon2 = v[4], azi2 = v[5];
        const double s12 = v[6], a12 = v[7], m12 = v[8], M12 = v[9], M21 = v[10], S12 = v[11];
        ndoc++;
        const Inverse r = g.inverse_degrees(lat1, lon1, lat2, lon2);
        doc_azi = std::max({doc_azi, std::fabs(dang(r.azi1 / D, azi1)), std::fabs(dang(r.azi2 / D, azi2))});
        doc_s = std::max(doc_s, std::fabs(r.s12 - s12));
        doc_a = std::max(doc_a, std::fabs(r.a12 / D - a12));
        doc_m = std::max(doc_m, std::fabs(r.m12 - m12));
        doc_M = std::max({doc_M, std::fabs(r.M12 - M12), std::fabs(r.M21 - M21)});
        doc_S = std::max(doc_S, std::fabs(r.S12 - S12));
        const Direct d = g.direct_degrees(lat1, lon1, azi1, s12, false, true);
        doc_pos = std::max({doc_pos, std::fabs(d.lat2 - lat2), std::fabs(d.lon2 - lon2)});
        doc_dazi = std::max(doc_dazi, std::fabs(dang(d.azi2, azi2)));
        doc_da = std::max(doc_da, std::fabs(d.a12 - a12));
        doc_dm = std::max(doc_dm, std::fabs(d.m12 - m12));
        doc_dM = std::max({doc_dM, std::fabs(d.M12 - M12), std::fabs(d.M21 - M21)});
        doc_dS = std::max(doc_dS, std::fabs(d.S12 - S12));
        continue;
      }
      std::string kind, bar;
      double in4[4], out8[8];
      f >> kind >> in4[0] >> in4[1] >> in4[2] >> in4[3] >> bar;
      for (double& x : out8) f >> x;
      if (type == "inv") {
        ninv++;
        const double lat1 = in4[0], lon1 = in4[1], lat2 = in4[2], lon2 = in4[3];
        const Inverse r = g.inverse_degrees(lat1, lon1, lat2, lon2);
        inv_s = std::max(inv_s, std::fabs(r.s12 - out8[0]));
        // An azimuth's error moves the other end of the line by it times the reduced length m12; for
        // lines of millimetres, or nearly antipodal ones, that is the measure (the angle itself is
        // ill-conditioned there), and the angle itself elsewhere.
        const double eazi = std::max(std::fabs(dang(r.azi1 / D, out8[1])), std::fabs(dang(r.azi2 / D, out8[2]))) * D;
        inv_shift = std::max(inv_shift, eazi * std::fabs(out8[4]));
        if (kind != "antipodal" && out8[0] > 1000) inv_azi = std::max(inv_azi, eazi);
        inv_a = std::max(inv_a, std::fabs(r.a12 / D - out8[3]));
        inv_m = std::max(inv_m, std::fabs(r.m12 - out8[4]));
        inv_M = std::max({inv_M, std::fabs(r.M12 - out8[5]), std::fabs(r.M21 - out8[6])});
        inv_S = std::max(inv_S, std::fabs(r.S12 - out8[7]) / std::max(1.0, std::fabs(out8[7]) * 1e-15));
        // The same through the interface in radians: inputs rounded once more in the conversion.
        const Inverse q = g.inverse(lat1 * D, lon1 * D, lat2 * D, lon2 * D);
        rad_s = std::max(rad_s, std::fabs(q.s12 - out8[0]));
        const double eq = std::max(std::fabs(std::remainder(q.azi1 - out8[1] * D, 2 * pi)), std::fabs(std::remainder(q.azi2 - out8[2] * D, 2 * pi)));
        rad_shift = std::max(rad_shift, eq * std::fabs(out8[4]));
        if (kind != "antipodal" && out8[0] > 1000) rad_azi = std::max(rad_azi, eq);
        // Round trip: the direct problem from the first point along azi1 for s12 reaches the second.
        const Direct d = g.direct_degrees(lat1, lon1, r.azi1 / D, r.s12);
        const double dn = (d.lat2 - lat2) * D * 6.4e6, de = dang(d.lon2, lon2) * D * 6.4e6 * std::cos(lat2 * D);
        const double miss = std::fabs(lat2) > 89.9999999 ? std::fabs(dn) : std::hypot(dn, de);
        (kind == "antipodal" ? round_trip_ill : round_trip) = std::max(kind == "antipodal" ? round_trip_ill : round_trip, miss);
      } else if (type == "dir") {
        ndir++;
        const double lat1 = in4[0], lon1 = in4[1], azi1 = in4[2], s12 = in4[3];
        const Direct d = g.direct_degrees(lat1, lon1, azi1, s12, false, true);
        dir_pos = std::max({dir_pos, std::fabs(d.lat2 - out8[0]), std::fabs(d.lon2 - out8[1])});
        dir_azi = std::max(dir_azi, std::fabs(dang(d.azi2, out8[2])));
        dir_a = std::max(dir_a, std::fabs(d.a12 - out8[3]));
        dir_m = std::max(dir_m, std::fabs(d.m12 - out8[4]));
        dir_M = std::max({dir_M, std::fabs(d.M12 - out8[5]), std::fabs(d.M21 - out8[6])});
        dir_S = std::max(dir_S, std::fabs(d.S12 - out8[7]) / std::max(1.0, std::fabs(out8[7]) * 1e-15));
      }
    }
    s.require("the reference file was read (620 inverse, 400 direct, 20 documented)", ninv == 620 && ndir == 400 && ndoc == 20,
              std::to_string(ninv) + ", " + std::to_string(ndir) + ", " + std::to_string(ndoc));
    s.check("inverse against GeographicLib: distance", inv_s, 15e-9, "m");
    s.check("inverse against GeographicLib: azimuths of lines over 1 km, not nearly antipodal", inv_azi, 1e-13, "rad");
    s.check("inverse against GeographicLib: azimuths, as the shift |Δα| m12 of the line's other end", inv_shift, 15e-9, "m");
    s.check("inverse against GeographicLib: arc length", inv_a * D, 1e-13, "rad");
    s.check("inverse against GeographicLib: reduced length", inv_m, 15e-9, "m");
    s.check("inverse against GeographicLib: geodesic scales", inv_M, 1e-14);
    s.check("inverse against GeographicLib: area to the equator (in units of 10⁻¹⁵ of it, or m²)", inv_S, 10);
    s.check("inverse in radians: distance (inputs rounded in the conversion)", rad_s, 15e-9, "m");
    s.check("inverse in radians: azimuths of lines over 1 km, not nearly antipodal", rad_azi, 1e-13, "rad");
    s.check("inverse in radians: azimuths, as the shift |Δα| m12 of the line's other end", rad_shift, 15e-9, "m");
    s.check("direct against GeographicLib: latitude and unrolled longitude", dir_pos * D, 1e-13, "rad");
    s.check("direct against GeographicLib: azimuth", dir_azi * D, 1e-13, "rad");
    s.check("direct against GeographicLib: arc length", dir_a * D, 1e-13, "rad");
    s.check("direct against GeographicLib: reduced length", dir_m, 15e-9, "m");
    s.check("direct against GeographicLib: geodesic scales", dir_M, 1e-14);
    s.check("direct against GeographicLib: area (in units of 10⁻¹⁵ of it, or m²)", dir_S, 10);
    s.check("inverse then direct reaches the second point", round_trip, 15e-9, "m");
    s.check("inverse then direct reaches the second point, nearly antipodal", round_trip_ill, 15e-9, "m");
    // GeographicLib's own test cases, with the tolerances its tests use.
    s.check("GeographicLib's test cases, inverse: azimuths", doc_azi, 1e-13, "deg");
    s.check("GeographicLib's test cases, inverse: distance", doc_s, 1e-8, "m");
    s.check("GeographicLib's test cases, inverse: arc length", doc_a, 1e-13, "deg");
    s.check("GeographicLib's test cases, inverse: reduced length", doc_m, 1e-8, "m");
    s.check("GeographicLib's test cases, inverse: geodesic scales", doc_M, 1e-15);
    s.check("GeographicLib's test cases, inverse: area", doc_S, 0.1, "m²");
    s.check("GeographicLib's test cases, direct: latitude and longitude", doc_pos, 1e-13, "deg");
    s.check("GeographicLib's test cases, direct: azimuth", doc_dazi, 1e-13, "deg");
    s.check("GeographicLib's test cases, direct: arc length", doc_da, 1e-13, "deg");
    s.check("GeographicLib's test cases, direct: reduced length", doc_dm, 1e-8, "m");
    s.check("GeographicLib's test cases, direct: geodesic scales", doc_dM, 1e-15);
    s.check("GeographicLib's test cases, direct: area", doc_dS, 0.1, "m²");
  }

  // ---- GeographicLib's GeodSolve examples (test_geodesic.py), with their stated tolerances ----------------
  {
    auto inv = [&](const Geodesic& e, double a, double b, double c, double d) { return e.inverse_degrees(a, b, c, d); };
    Inverse r = inv(g, 40.6, -73.8, 49.01666667, 2.55);  // JFK to CDG
    s.check("GeodSolve0: JFK to CDG, 5 853 226 m at 53.47022° (GeographicLib)",
            std::max({std::fabs(r.azi1 / D - 53.47022) / 0.5e-5, std::fabs(r.azi2 / D - 111.59367) / 0.5e-5, std::fabs(r.s12 - 5853226) / 0.5}), 1);
    Direct d = g.direct_degrees(40.63972222, -73.77888889, 53.5, 5850e3);
    s.check("GeodSolve1: the direct problem from JFK",
            std::max({std::fabs(d.lat2 - 49.01467), std::fabs(d.lon2 - 2.56106), std::fabs(d.azi2 - 111.62947)}) / 0.5e-5, 1);
    const Geodesic prolate(6.4e6, -1 / 150.0);
    r = inv(prolate, 0.07476, 0, -0.07476, 180);
    Inverse r2 = inv(prolate, 0.1, 0, -0.1, 180);
    s.check("GeodSolve2: antipodal points of a prolate ellipsoid",
            std::max({std::fabs(r.azi1 / D - 90.00078) / 0.5e-5, std::fabs(r.azi2 / D - 90.00078) / 0.5e-5, std::fabs(r.s12 - 20106193) / 0.5,
                      std::fabs(r2.azi1 / D - 90.00105) / 0.5e-5, std::fabs(r2.azi2 / D - 90.00105) / 0.5e-5, std::fabs(r2.s12 - 20106193) / 0.5}),
            1);
    r = inv(g, 36.493349428792, 0, 36.49334942879201, 0.0000008);
    s.check("GeodSolve4: a very short line, 0.072 m", (r.s12 - 0.072) / 0.5e-3, 1);
    d = g.direct_degrees(0.01777745589997, 30, 0, 10e6);
    s.check("GeodSolve5: a line ending at the pole", std::max({std::fabs(d.lat2 - 90), d.lon2 < 0 ? std::fabs(d.lon2 + 150) : std::fabs(d.lon2 - 30)}) / 0.5e-5, 1);
    const double hard[][5] = {{88.202499451857, 0, -88.202499451857, 179.981022032992859592, 20003898.214},
                              {89.262080389218, 0, -89.262080389218, 179.992207982775375662, 20003925.854},
                              {89.333123580033, 0, -89.333123580032997687, 179.99295812360148422, 20003926.881},
                              {56.320923501171, 0, -56.320923501171, 179.664747671772880215, 19993558.287},
                              {52.784459512564, 0, -52.784459512563990912, 179.634407464943777557, 19991596.095},
                              {48.522876735459, 0, -48.52287673545898293, 179.599720456223079643, 19989144.774}};
    double worst = 0;
    for (const auto& h : hard) worst = std::max(worst, std::fabs(inv(g, h[0], h[1], h[2], h[3]).s12 - h[4]) / 0.5e-3);
    s.check("GeodSolve6, 9, 10, 11: nearly antipodal pairs that failed other implementations", worst, 1);
    const Geodesic odd(89.8, -1.83);
    r = inv(odd, 0, 0, -10, 160);
    s.check("GeodSolve12: an extremely prolate ellipsoid",
            std::max({std::fabs(r.azi1 / D - 120.27) / 1e-2, std::fabs(r.azi2 / D - 105.15) / 1e-2, std::fabs(r.s12 - 266.7) / 1e-1}), 1);
    r = inv(g, 0, 0, 1, nanv());
    s.require("GeodSolve14: a NaN longitude gives NaNs", std::isnan(r.azi1) && std::isnan(r.azi2) && std::isnan(r.s12));
    r = inv(g, nanv(), 0, 0, 90);
    s.require("GeodSolve55: a NaN latitude gives NaNs", std::isnan(r.azi1) && std::isnan(r.azi2) && std::isnan(r.s12));
    // GeodSolve33: along the equator to near the antipode, on WGS 84, a sphere and a prolate ellipsoid.
    const Geodesic sphere(6.4e6, 0), prolate300(6.4e6, -1 / 300.0);
    struct Case { const Geodesic* e; double lat2, lon2, azi1, azi2, s12; };
    const Case eq[] = {{&g, 0, 179, 90, 90, 19926189},       {&g, 0, 179.5, 55.96650, 124.03350, 19980862},
                       {&g, 0, 180, 0, 180, 20003931},        {&g, 1, 180, 0, 180, 19893357},
                       {&sphere, 0, 179, 90, 90, 19994492},   {&sphere, 0, 180, 0, 180, 20106193},
                       {&sphere, 1, 180, 0, 180, 19994492},   {&prolate300, 0, 179, 90, 90, 19994492},
                       {&prolate300, 0, 180, 90, 90, 20106193}, {&prolate300, 0.5, 180, 33.02493, 146.97364, 20082617},
                       {&prolate300, 1, 180, 0, 180, 20027270}};
    worst = 0;
    for (const auto& c : eq) {
      r = inv(*c.e, 0, 0, c.lat2, c.lon2);
      worst = std::max({worst, std::fabs(std::fabs(r.azi1 / D) - c.azi1) / 0.5e-5, std::fabs(std::fabs(r.azi2 / D) - c.azi2) / 0.5e-5, std::fabs(r.s12 - c.s12) / 0.5});
    }
    s.check("GeodSolve33: lines along the equator toward the antipode (WGS 84, a sphere, a prolate ellipsoid)", worst, 1);
    r = inv(g, 5, 0.00000000000001, 10, 180);
    s.check("GeodSolve59: longitudes nearly 180° apart",
            std::max({std::fabs(r.azi1 / D - 0.000000000000035) / 1.5e-14, std::fabs(r.azi2 / D - 179.99999999999996) / 1.5e-14, std::fabs(r.s12 - 18345191.174332713) / 5e-9}), 1);
    d = g.direct_degrees(45, 0, -0.000000000000000003, 1e7, false, true);
    s.check("GeodSolve61: a small negative azimuth goes west",
            std::max({std::fabs(d.lat2 - 45.30632), std::fabs(d.lon2 + 180), std::fabs(std::fabs(d.azi2) - 180)}) / 0.5e-5, 1);
    d = g.direct_degrees(90, 10, 180, -1e6);
    s.check("GeodSolve73: backwards from the pole", std::max({std::fabs(d.lat2 - 81.04623), std::fabs(d.lon2 + 170), std::fabs(d.azi2)}) / 0.5e-5, 1);
    s.require("GeodSolve73: the azimuth there is +0", std::copysign(1.0, d.azi2) > 0);
    r = inv(g, 54.1589, 15.3872, 54.1591, 15.3877);
    s.check("GeodSolve74: a short line with all its quantities",
            std::max({std::fabs(r.azi1 / D - 55.723110355), std::fabs(r.azi2 / D - 55.723515675), std::fabs(r.s12 - 39.527686385), std::fabs(r.a12 / D - 0.000355495),
                      std::fabs(r.m12 - 39.527686385), std::fabs(r.M12 - 0.999999995), std::fabs(r.M21 - 0.999999995)}) / 5e-9,
            1);
    s.check("GeodSolve74: the area of the short line", std::fabs(r.S12 - 286698586.30197), 5e-4, "m²");
    r = inv(g, -(41 + 19 / 60.0), 174 + 49 / 60.0, 40 + 58 / 60.0, -(5 + 30 / 60.0));
    s.check("GeodSolve76: Wellington to Salamanca (a classic failure of Vincenty)",
            std::max({std::fabs(r.azi1 / D - 160.39137649664) / 0.5e-11, std::fabs(r.azi2 / D - 19.50042925176) / 0.5e-11, std::fabs(r.s12 - 19960543.857179) / 0.5e-6}), 1);
    r = inv(g, 27.2, 0.0, -27.1, 179.5);
    s.check("GeodSolve78: where the NGS calculator fails to converge",
            std::max({std::fabs(r.azi1 / D - 45.82468716758) / 0.5e-11, std::fabs(r.azi2 / D - 134.22776532670) / 0.5e-11, std::fabs(r.s12 - 19974354.765767) / 0.5e-6}), 1);
    r = inv(g, 20.001, 0, 20.001, 0);
    s.require("GeodSolve80: a geodesic of zero length is +0 long, heading south (180°), with unit scales",
              r.s12 == 0 && std::copysign(1.0, r.s12) > 0 && r.a12 == 0 && std::fabs(r.azi1 / D - 180) < 1e-13 && std::fabs(r.azi2 / D - 180) < 1e-13 &&
                  r.m12 == 0 && std::fabs(r.M12 - 1) < 1e-15 && std::fabs(r.M21 - 1) < 1e-15 && std::fabs(r.S12) < 1e-10);
    r = inv(g, 0, 0, 0, 90);
    s.check("GeodSolve80: the geodesic scale a quarter of the equator away", std::max(std::fabs(r.M12 + 0.00528427534), std::fabs(r.M21 + 0.00528427534)), 0.5e-10);
    // GeodSolve65 and 66: lines nearly south across the antipode, with the line's end.
    GeodesicLine l65 = g.inverse_line(30 * D, -0.000000000000000001 * D, -31 * D, 180 * D);
    d = l65.position_degrees(1e7, false, true);
    s.check("GeodSolve65: a point along a line across the pole",
            std::max({std::fabs(d.lat2 + 60.23169) / 0.5e-5, std::fabs(d.lon2) / 0.5e-5, std::fabs(std::fabs(d.azi2) - 180) / 0.5e-5, std::fabs(d.a12 - 90.06544) / 0.5e-5,
                      std::fabs(d.m12 - 6363636) / 0.5, std::fabs(d.M21 - 0.0013749) / 0.5e-7, std::fabs(d.S12) / 0.5}),
            1);
    d = l65.position_degrees(2e7, false, true);
    s.check("GeodSolve65: the same line near the antipode",
            std::max({std::fabs(d.lat2 + 30.03547) / 0.5e-5, std::fabs(d.lon2 + 180) / 0.5e-5, std::fabs(d.azi2) / 0.5e-5, std::fabs(d.a12 - 179.96459) / 0.5e-5,
                      std::fabs(d.m12 - 54342) / 0.5, std::fabs(d.M21 + 0.9954339) / 0.5e-7, std::fabs(d.S12 - 127516405431022.0) / 0.5}),
            1);
    GeodesicLine l66 = g.inverse_line(-5 * D, -0.000000000000002 * D, -10 * D, 180 * D);
    d = l66.position_degrees(0.5 * l66.distance(), false, true);
    s.check("GeodSolve66: halfway along a line slightly west of south, through the pole",
            std::max({std::fabs(d.lat2 + 87.52461), std::fabs(d.lon2), std::fabs(d.azi2 + 180)}) / 0.5e-5, 1);
  }

  // ---- points along a line --------------------------------------------------------------------------
  {
    std::mt19937_64 rng(20261028);
    std::uniform_real_distribution<double> u(-1, 1);
    double worst = 0, end = 0, unrolled = 0;
    for (int k = 0; k < 60; k++) {
      const double lat1 = std::asin(u(rng)), lon1 = pi * u(rng), lat2 = std::asin(u(rng)), lon2 = pi * u(rng);
      const GeodesicLine l = g.inverse_line(lat1, lon1, lat2, lon2);
      const auto pts = l.points(16);
      for (int i = 0; i <= 16; i++) {
        const Direct d = g.direct(lat1, lon1, l.azimuth(), l.distance() * i / 16);
        worst = std::max({worst, std::fabs(pts[i].lat2 - d.lat2), std::fabs(std::remainder(pts[i].lon2 - d.lon2, 2 * pi)), std::fabs(std::remainder(pts[i].azi2 - d.azi2, 2 * pi))});
      }
      const Direct last = pts.back();
      end = std::max(end, std::hypot(last.lat2 - lat2, std::remainder(last.lon2 - lon2, 2 * pi) * std::cos(lat2)) * 6.4e6);
      // Unrolled longitudes differ from the reduced ones by whole turns only.
      const Direct ur = l.position(l.distance(), true);
      unrolled = std::max(unrolled, std::fabs(std::remainder(ur.lon2 - last.lon2, 2 * pi)));
    }
    s.check("a line's points are the direct problem's, at every step (to rounding)", worst, 4e-15, "rad");
    s.check("a line made between two points ends at the second", end, 15e-9, "m");
    s.check("an unrolled longitude is the reduced one, in whole turns (to rounding)", unrolled, 4e-15, "rad");
    const GeodesicLine paris_tokyo = g.inverse_line(48.8566 * D, 2.3522 * D, 35.6895 * D, 139.6917 * D);
    const Inverse pt = g.inverse(48.8566 * D, 2.3522 * D, 35.6895 * D, 139.6917 * D);
    s.check("a line knows its distance: that of the inverse problem", paris_tokyo.distance() - pt.s12, 1e-8, "m");
    s.check("and its arc length", paris_tokyo.arc() - pt.a12, 1e-15, "rad");
  }

  // ---- the sphere: great circles -------------------------------------------------------------------
  {
    const double R = 6371008.8;
    const Geodesic sphere(R, 0);
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-1, 1);
    double ds = 0, dazi = 0, ddir = 0;
    for (int k = 0; k < 300; k++) {
      const double p1 = std::asin(u(rng)), l1 = pi * u(rng), p2 = std::asin(u(rng)), l2 = pi * u(rng);
      const double x1[3] = {std::cos(p1) * std::cos(l1), std::cos(p1) * std::sin(l1), std::sin(p1)};
      const double x2[3] = {std::cos(p2) * std::cos(l2), std::cos(p2) * std::sin(l2), std::sin(p2)};
      const double cr[3] = {x1[1] * x2[2] - x1[2] * x2[1], x1[2] * x2[0] - x1[0] * x2[2], x1[0] * x2[1] - x1[1] * x2[0]};
      const double sigma = std::atan2(std::hypot(cr[0], cr[1], cr[2]), x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]);
      const double azi = std::atan2(std::sin(l2 - l1) * std::cos(p2), std::cos(p1) * std::sin(p2) - std::sin(p1) * std::cos(p2) * std::cos(l2 - l1));
      const Inverse r = sphere.inverse(p1, l1, p2, l2);
      ds = std::max(ds, std::fabs(r.s12 - R * sigma));
      if (sigma < 179 * D && sigma > 1e-3) dazi = std::max(dazi, std::fabs(std::remainder(r.azi1 - azi, 2 * pi)));
      // The direct problem on the sphere: the point at angle σ along azimuth α.
      const double a = 2 * pi * u(rng), sg = pi * std::fabs(u(rng));
      const Direct d = sphere.direct(p1, l1, a, R * sg);
      const double p = std::asin(std::sin(p1) * std::cos(sg) + std::cos(p1) * std::sin(sg) * std::cos(a));
      const double l = l1 + std::atan2(std::sin(a) * std::sin(sg) * std::cos(p1), std::cos(sg) - std::sin(p1) * std::sin(p));
      ddir = std::max({ddir, std::fabs(d.lat2 - p), std::fabs(std::remainder(d.lon2 - l, 2 * pi)) * std::cos(p)});
    }
    s.check("on a sphere, the distance is R times the central angle", ds, 5e-8, "m");
    s.check("on a sphere, the azimuth is the great circle's", dazi, 1e-13, "rad");
    s.check("on a sphere, the direct problem lands where the great circle does", ddir, 1e-14, "rad");
  }

  // ---- the interface in radians ---------------------------------------------------------------------
  {
    const Inverse a = g.inverse(pi / 2, 0.3, -pi / 2, -1.1);
    s.check("pole to pole is twice the meridian quadrant, 20 003 931.4586 m", a.s12 - 2 * 10001965.7293, 1e-3, "m");
    const Inverse b = g.inverse(std::nextafter(pi / 2, 1.0), 0, 0, 0);
    s.require("the rounding of π/2 beyond the pole is the pole", std::isfinite(b.s12) && std::fabs(b.s12 - 10001965.7293) < 1e-3);
    const Inverse c = g.inverse(1.6, 0, 0, 0);
    s.require("a latitude beyond the pole is not one", std::isnan(c.s12));
    s.check("the area of WGS 84 is 510 065 621.7 km² (4π c², c the authalic radius)", g.ellipsoid_area() / 1e6 - 510065621.724, 1e-2, "km²");
  }

  return s.done();
}
