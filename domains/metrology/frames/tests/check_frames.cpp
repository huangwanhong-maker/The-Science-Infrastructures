// Checks of the Earth's orientation: precession, nutation, the CIO locator, sidereal time, polar
// motion and the celestial-to-terrestrial matrices (metrology/frames), and the frame graph of the
// geocentric frames with the holonomy of the IERS's two routes.
//
// References: ERFA's own test values (t_erfa_c.c, the SOFA test suite) written in below, and
// tests/fixtures/erfa-reference.txt, ERFA through pyerfa at fifteen dates from 1850 to 2200.

#include <cmath>
#include <string>

#include <calx/verification/check.hpp>
#include <metrology/frames/earth_frames.hpp>
#include <metrology/frames/iau1982.hpp>
#include <metrology/frames/iau2006.hpp>
#include <metrology/time/scales.hpp>
#include <physis/geometry/frame_graph.hpp>

#include "../../tests/reference.hpp"

using namespace metrology;
using namespace metrology::frames;
using calx::verification::Suite;
using physis::geometry::max_difference;

namespace {

Mat3 matrix_from(const checks::Row& r, std::size_t first) {
  Mat3 m{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m[i][j] = r.number(first + 3 * i + j);
  return m;
}
double angle_difference(double a, double b) { return std::fabs(signed_angle(a - b)); }

}  // namespace

int main() {
  Suite s;
  auto ref = checks::read_reference(METROLOGY_FIXTURES "/erfa-reference.txt");
  const time::JulianDate t53736{2400000.5, 53736.0};

  // ---- SOFA's own test values ------------------------------------------------------------------
  {
    const Nutation n = nutation(t53736);
    s.check("nutation Δψ, IAU 2006/2000A (SOFA test)", n.dpsi + 0.9630912025820308797e-5, 1e-13, "rad");
    s.check("nutation Δε, IAU 2006/2000A (SOFA test)", n.deps - 0.4063238496887249798e-4, 1e-13, "rad");
    s.check("the CIO locator s (SOFA test)", cio_locator(t53736, 0.5791308486706011000e-3, 0.4020579816732961219e-4) + 0.1220032213076463117e-7, 1e-18, "rad");
    s.check("the Earth rotation angle (SOFA test)", earth_rotation_angle({2400000.5, 54388.0}) - 0.4022837240028158102, 1e-12, "rad");
    s.check("GMST, IAU 2006 (SOFA test)", gmst(t53736, t53736) - 1.754174971870091203, 1e-12, "rad");
    s.check("GAST, IAU 2006/2000A (SOFA test)", gast(t53736, t53736) - 1.754166137675019159, 1e-12, "rad");
    s.check("the complementary terms of the equation of the equinoxes (SOFA test)", equation_of_equinoxes_complementary(t53736) - 0.2046085004885125264e-8, 1e-20, "rad");
    s.check("s′ (SOFA test)", tio_locator({2400000.5, 52541.0}) + 0.6216698469981019309e-11, 1e-12, "rad");
    const Mat3 npb = bias_precession_nutation({2400000.5, 50123.9999});
    const Mat3 npb_sofa{{{0.9999995832794205484, 0.8372382772630962111e-3, 0.3639684771140623099e-3},
                         {-0.8372533744743683605e-3, 0.9999996486492861646, 0.4132905944611019498e-4},
                         {-0.3639337469629464969e-3, -0.4163377605910663999e-4, 0.9999999329094260057}}};
    s.check("the NPB matrix (SOFA test, iauPnm06a)", max_difference(npb, npb_sofa), 1e-14);
    const Mat3 w = polar_motion(2.55060238e-7, 1.860359247e-6, -0.1367174580728891460e-10);
    const Mat3 w_sofa{{{0.9999999999999674721, -0.1367174580728846989e-10, 0.2550602379999972345e-6},
                       {0.1414624947957029801e-10, 0.9999999999982695317, -0.1860359246998866389e-5},
                       {-0.2550602379741215021e-6, 0.1860359247002414021e-5, 0.9999999999982370039}}};
    s.check("the polar motion matrix (SOFA test, iauPom00)", max_difference(w, w_sofa), 1e-16);
    const EarthOrientation e = orientation_of_earth(t53736, t53736, 2.55060238e-7 / ARCSEC, 1.860359247e-6 / ARCSEC);
    const Mat3 c2t_sofa{{{-0.1810332128305897282, 0.9834769806938592296, 0.6555550962998436505e-4},
                         {-0.9834768134136214897, -0.1810332203649130832, 0.5749800844905594110e-3},
                         {0.5773474024748545878e-3, 0.3961816829632690581e-4, 0.9999998325501747785}}};
    s.check("GCRS to ITRS, the CIO route (SOFA test, iauC2t06a)", max_difference(e.c2t, c2t_sofa), 1e-12);
  }

  // ---- the checks of the domain's JavaScript (Meeus 1998, Astronomical Algorithms, 2nd ed.) ------
  {
    // Example 12.a: 1987 April 10, 0h UT: GMST 13h 10m 46.3668s (IAU 1982; IAU 2006 differs by ms).
    const time::JulianDate ut{2446895.5, 0.0};
    const double hms = (13 + 10 / 60.0 + 46.3668 / 3600) * 15 * PI / 180;
    s.check("GMST on 1987-04-10 0h UT (Meeus 12.a)", signed_angle(gmst(ut, ut.plus_seconds(55.184)) - hms) * 43200 / PI, 0.01, "s");
    // Example 21.b: θ Persei from the mean equator of J2000.0 to that of 2028 November 13.19 TD.
    const double ra0 = (2 + 44 / 60.0 + 12.975 / 3600) * 15 * PI / 180, de0 = (49 + 13 / 60.0 + 39.90 / 3600) * PI / 180;
    const Vec3 v0{std::cos(de0) * std::cos(ra0), std::cos(de0) * std::sin(ra0), std::sin(de0)};
    const Mat3 p = bias_precession({2462088.69, 0.0}) * physis::geometry::transpose(frame_bias());
    const Vec3 v = p * v0;
    const double ra = std::atan2(v[1], v[0]), de = std::asin(v[2]);
    const double ra_meeus = (2 + 46 / 60.0 + 11.331 / 3600) * 15 * PI / 180, de_meeus = (49 + 20 / 60.0 + 54.54 / 3600) * PI / 180;
    s.check("precessed right ascension (Meeus 21.b; IAU 2006 against IAU 1976)", signed_angle(ra - ra_meeus) * 43200 / PI, 0.01, "s");
    s.check("precessed declination (Meeus 21.b)", (de - de_meeus) / ARCSEC, 0.1, "arcsec");
  }

  // ---- ERFA at fifteen dates, 1850–2200 --------------------------------------------------------
  {
    double nut = 0, nut00 = 0, pfw = 0, obl = 0, pnm = 0, s06 = 0, ee = 0, era = 0, gm = 0, gm82 = 0, gs = 0, gc = 0, sp = 0, c2t = 0, holonomy_ref = 0;
    for (const auto& r : ref["nut06a"]) {
      const Nutation n = nutation({r.number(0), r.number(1)});
      nut = std::max({nut, std::fabs(n.dpsi - r.number(2)), std::fabs(n.deps - r.number(3))});
    }
    for (const auto& r : ref["nut00a"]) {
      const Nutation n = nutation_2000a({r.number(0), r.number(1)});
      nut00 = std::max({nut00, std::fabs(n.dpsi - r.number(2)), std::fabs(n.deps - r.number(3))});
    }
    for (const auto& r : ref["pfw06"]) {
      const auto a = precession_angles({r.number(0), r.number(1)});
      pfw = std::max({pfw, std::fabs(a.gamb - r.number(2)), std::fabs(a.phib - r.number(3)), std::fabs(a.psib - r.number(4)), std::fabs(a.epsa - r.number(5))});
    }
    for (const auto& r : ref["obl06"]) obl = std::max(obl, std::fabs(mean_obliquity({r.number(0), r.number(1)}) - r.number(2)));
    for (const auto& r : ref["pnm06a"]) pnm = std::max(pnm, max_difference(bias_precession_nutation({r.number(0), r.number(1)}), matrix_from(r, 2)));
    for (const auto& r : ref["s06"]) s06 = std::max(s06, std::fabs(cio_locator({r.number(0), r.number(1)}, r.number(2), r.number(3)) - r.number(4)));
    for (const auto& r : ref["eect00"]) ee = std::max(ee, std::fabs(equation_of_equinoxes_complementary({r.number(0), r.number(1)}) - r.number(2)));
    for (const auto& r : ref["era00"]) era = std::max(era, angle_difference(earth_rotation_angle({r.number(0), r.number(1)}), r.number(2)));
    for (const auto& r : ref["gmst06"]) gm = std::max(gm, angle_difference(gmst({r.number(0), r.number(1)}, {r.number(2), r.number(3)}), r.number(4)));
    for (const auto& r : ref["gst06a"]) gs = std::max(gs, angle_difference(gast({r.number(0), r.number(1)}, {r.number(2), r.number(3)}), r.number(4)));
    for (const auto& r : ref["gmst82"]) gm82 = std::max(gm82, angle_difference(gmst1982(metrology::time::JulianDate{r.number(0), r.number(1)}), r.number(2)));
    for (const auto& r : ref["gast_classical"]) gc = std::max(gc, angle_difference(gast_classical({r.number(0), r.number(1)}, {r.number(2), r.number(3)}), r.number(4)));
    for (const auto& r : ref["sp00"]) sp = std::max(sp, std::fabs(tio_locator({r.number(0), r.number(1)}) - r.number(2)));
    for (const auto& r : ref["c2t06a"]) {
      const EarthOrientation e = orientation_of_earth({r.number(0), r.number(1)}, {r.number(2), r.number(3)}, r.number(4) / ARCSEC, r.number(5) / ARCSEC);
      c2t = std::max(c2t, max_difference(e.c2t, matrix_from(r, 6)));
    }
    // ERFA's own two sidereal times differ by what the frame graph should find as a holonomy.
    for (std::size_t i = 0; i < ref["gst06a"].size() && i < ref["gast_classical"].size(); i++)
      holonomy_ref = std::max(holonomy_ref, angle_difference(ref["gst06a"][i].number(4), ref["gast_classical"][i].number(4)));
    const bool all = !ref["nut06a"].empty() && !ref["c2t06a"].empty();
    s.require("the reference file has the frames' rows", all);
    s.check("IAU 2006/2000A nutation against ERFA, 1850–2200", nut, 1e-15, "rad");
    s.check("IAU 2000A nutation against ERFA", nut00, 1e-15, "rad");
    s.check("Fukushima–Williams precession angles against ERFA", pfw, 1e-15, "rad");
    s.check("the mean obliquity against ERFA", obl, 1e-15, "rad");
    s.check("the NPB matrix against ERFA", pnm, 1e-15);
    s.check("the CIO locator s against ERFA", s06, 1e-18, "rad");
    s.check("the equation of the equinoxes' complementary terms against ERFA", ee, 1e-18, "rad");
    s.check("the Earth rotation angle against ERFA", era, 1e-15, "rad");
    s.check("GMST against ERFA", gm, 1e-15, "rad");
    s.require("the reference file has GMST 1982's rows", !ref["gmst82"].empty());
    // Both compute in seconds, where the term 8640184.812866 s T is 1.7e7 s at T = 2 (2200) and its
    // last bit 3.7 ns, 2.7e-13 rad: the formula's own rounding, ERFA's and ours, bounds the agreement.
    s.check("GMST of IAU 1982 against ERFA's gmst82, to the rounding of its seconds", gm82, 3e-13, "rad");
    s.check("GAST (CIO route) against ERFA", gs, 1e-15, "rad");
    s.check("GAST (GMST + EE) against ERFA", gc, 1e-15, "rad");
    s.check("s′ against ERFA", sp, 1e-20, "rad");
    s.check("GCRS → ITRS against ERFA, with polar motion", c2t, 1e-15);
    s.check("the two sidereal times of ERFA differ by under 1 µas over 1850–2200", holonomy_ref, 4.85e-12, "rad");
  }

  // ---- the frame graph of the geocentric frames --------------------------------------------------
  {
    physis::geometry::FrameGraph<time::Instant> g;
    const int gcrs = g.root("GCRS");
    const int itrs = add_earth_frames(g, gcrs);
    // An instant with polar motion and pole offsets as the IERS gives them (2026-10-01).
    earth_orientation::Values eop;
    eop.mjd = 61314;
    eop.xp = 0.174596, eop.yp = 0.325269, eop.dut1 = -0.0225232, eop.lod = 0.6645e-3, eop.dx = 0.365e-3, eop.dy = -0.104e-3;
    eop.has_nutation = true;
    const auto eop_series = earth_orientation::Series::from_rows({eop, [&] { auto e = eop; e.mjd += 1; return e; }()}, "a check");
    const time::TimeScales ts(time::LeapSeconds::built_in(), &eop_series, time::Ut1Policy::refuse);
    const time::Instant now = ts.instant(time::JulianDate::from_mjd(61314.3), time::Scale::UTC);

    // Through the graph against the direct product W R_z(θ) C with the same inputs.
    const EarthOrientation e = orientation_of_earth(now.tt, now.ut1, eop.xp, eop.yp, eop.dx, eop.dy);
    const auto t = g.transform(gcrs, itrs, now);
    s.check("GCRS → ITRS through CIRS and TIRS in the graph equals W R_z(θ) C", max_difference(t.r, e.c2t), 1e-15);
    // Velocity: a point fixed on the equator at the prime meridian moves at ω R eastward in the GCRS.
    const auto back = g.transform(itrs, gcrs, now);
    const auto state = back.state({6378137.0, 0, 0}, {0, 0, 0});
    const double speed = physis::geometry::norm(state[1]);
    s.check("a point on the equator moves at ω a (465.1 m/s) in the GCRS", speed - earth_rotation_rate(eop.lod) * 6378137.0, 1e-3, "m/s");
    const auto ahead = g.transform(itrs, gcrs, ts.instant(time::JulianDate::from_mjd(61314.3).plus_seconds(0.5), time::Scale::UTC));
    const auto behind = g.transform(itrs, gcrs, ts.instant(time::JulianDate::from_mjd(61314.3).plus_seconds(-0.5), time::Scale::UTC));
    const auto fd = (1.0 / 1.0) * (ahead.position({6378137.0, 0, 0}) - behind.position({6378137.0, 0, 0}));
    s.check("that velocity against a central difference of positions 1 s apart", physis::geometry::norm(fd - state[1]), 1e-3, "m/s");
    // The holonomy of the IERS's two routes: about a microarcsecond of model, and the pole offsets.
    double route = 0, pef = 0;
    for (std::size_t i = 0; i < g.relations().size(); i++) {
      const auto d = physis::geometry::deviation(g.holonomy(int(i), now));
      if (g.relations()[i].name == "ITRS by the equinox route") route = d.angle;
      if (g.relations()[i].name == "PEF is TIRS") pef = d.angle;
    }
    // The equinox route ignores dX, dY (0.38 mas here); its holonomy is their size.
    s.check("holonomy of the CIO and equinox routes to ITRS: the pole offsets dX, dY (0.38 mas)", route - std::hypot(eop.dx, eop.dy) * ARCSEC, 2e-11, "rad");
    s.check("holonomy of PEF against TIRS: the same 0.38 mas", pef - std::hypot(eop.dx, eop.dy) * ARCSEC, 2e-11, "rad");
    // Without pole offsets the two routes agree to the model's microarcsecond.
    earth_orientation::Values eop0 = eop;
    eop0.dx = eop0.dy = 0;
    const auto series0 = earth_orientation::Series::from_rows({eop0, [&] { auto e2 = eop0; e2.mjd += 1; return e2; }()}, "a check");
    const time::TimeScales ts0(time::LeapSeconds::built_in(), &series0, time::Ut1Policy::refuse);
    const time::Instant now0 = ts0.instant(time::JulianDate::from_mjd(61314.3), time::Scale::UTC);
    double model = 0;
    for (std::size_t i = 0; i < g.relations().size(); i++) model = std::max(model, physis::geometry::deviation(g.holonomy(int(i), now0)).angle);
    s.check("without dX, dY the two routes agree to under 1 µas", model, 4.85e-12, "rad");
    // EME2000 is GCRS with the frame bias: 23 mas.
    const auto bias = physis::geometry::deviation(g.transform("EME2000", "GCRS", now));
    s.check("the frame bias between EME2000 and GCRS is 23.1 mas", bias.angle / ARCSEC - 0.0231, 5e-4, "arcsec");
  }
  return s.done();
}
