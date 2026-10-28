// The geocentric frames of the Earth's orientation, as frames of a physis frame graph.
//
// Each frame here has its origin at the geocentre, so the edges are rotations; the frame graph
// (physis/geometry/frame_graph.hpp) takes a metrology Instant as its time, from which each edge
// reads the scale it needs (TT for precession and nutation, UT1 for the Earth's rotation) and
// the Earth orientation parameters. Each edge maps coordinates in the child to coordinates in
// the parent, with its time derivative (in rad per SI second, TT), so that velocities transform
// too: the Earth's rotation exactly, by dθ/dt corrected for the length of day; precession and
// nutation by central differences of an hour, whose error is below 10⁻¹⁴ rad/s.
//
//   GCRS     the Geocentric Celestial Reference System (the root given)
//   EME2000  the mean equator and equinox of J2000.0, "J2000" of older software: GCRS with the frame bias
//   MOD      the mean equator and equinox of date (IAU 2006 precession)
//   TOD      the true equator and equinox of date (IAU 2006/2000A nutation)
//   CIRS     the Celestial Intermediate Reference System
//   TIRS     the Terrestrial Intermediate Reference System: CIRS turned by the Earth rotation angle
//   ITRS     the International Terrestrial Reference System: TIRS with polar motion
//   PEF      TOD turned by the apparent sidereal time: the "pseudo Earth-fixed" frame, ITRS
//            without polar motion by the equinox route
//
// Two relations close loops whose holonomy measures the agreement of the IERS's two routes:
//   "ITRS by the equinox route"  GCRS → ITRS as W R_z(GMST + EE) N P B, against the tree's CIO route;
//   "PEF is TIRS"                 the frames are the same by the IERS's definitions; their
//                                 holonomy is the difference of the two apparent sidereal times.

#pragma once

#include <string>

#include <metrology/frames/iau2006.hpp>
#include <metrology/time/scales.hpp>
#include <physis/geometry/frame_graph.hpp>
#include <physis/geometry/kinematics.hpp>

namespace metrology::frames {

using physis::geometry::FrameGraph;
using physis::geometry::StateTransform;
using physis::geometry::transpose;

// An edge that depends on TT alone, slowly: its rotation, with the derivative by central
// differences of an hour.
template <class F>
StateTransform slow_rotation(const time::Instant& t, F&& matrix_of_tt) {
  const double h = 3600.0;
  const Mat3 r = matrix_of_tt(t.tt);
  const Mat3 plus = matrix_of_tt(t.tt.plus_seconds(h)), minus = matrix_of_tt(t.tt.plus_seconds(-h));
  return {r, (1 / (2 * h)) * (plus - minus), {0, 0, 0}, {0, 0, 0}};
}

// The observed polar motion and pole offsets at the instant, arcseconds (zero when UT1 was
// approximated and no measurement holds).
struct PoleValues {
  double xp, yp, dx, dy;
};
inline PoleValues pole(const time::Instant& t) { return {t.eop.xp, t.eop.yp, t.eop.has_nutation ? t.eop.dx : 0.0, t.eop.has_nutation ? t.eop.dy : 0.0}; }

// Adds the frames above to a graph, under the frame `gcrs`, and the two relations. Returns the
// index of ITRS.
inline int add_earth_frames(FrameGraph<time::Instant>& g, int gcrs) {
  using Edge = FrameGraph<time::Instant>::Edge;
  g.add("EME2000", gcrs, Edge([](const time::Instant&) { return StateTransform::fixed(transpose(frame_bias())); }),
        "the mean equator and equinox of J2000.0, geocentric (GCRS with the frame bias)");
  g.add("MOD", gcrs, Edge([](const time::Instant& t) {
          return slow_rotation(t, [](const time::JulianDate& tt) { return transpose(bias_precession(tt)); });
        }),
        "the mean equator and equinox of date, geocentric (IAU 2006 precession)");
  const int tod = g.add("TOD", gcrs, Edge([](const time::Instant& t) {
                          return slow_rotation(t, [](const time::JulianDate& tt) { return transpose(bias_precession_nutation(tt)); });
                        }),
                        "the true equator and equinox of date, geocentric (IAU 2006/2000A)");
  const int cirs = g.add("CIRS", gcrs, Edge([](const time::Instant& t) {
                           const PoleValues p = pole(t);
                           return slow_rotation(t, [&](const time::JulianDate& tt) {
                             const Mat3 npb = bias_precession_nutation(tt);
                             const auto xy = cip_xy(npb);
                             const double x = xy[0] + p.dx * ARCSEC, y = xy[1] + p.dy * ARCSEC;
                             return transpose(celestial_to_intermediate(x, y, cio_locator(tt, x, y)));
                           });
                         }),
                         "the Celestial Intermediate Reference System (IAU 2006/2000A, with the IERS's dX, dY)");
  const int tirs = g.add("TIRS", cirs, Edge([](const time::Instant& t) {
                           const double era = earth_rotation_angle(t.ut1);
                           const double rate = earth_rotation_rate(t.eop.lod);
                           // v_CIRS = R_z(θ)ᵀ v_TIRS, so d/dt = (dR_z/dθ)ᵀ θ̇.
                           return StateTransform{transpose(frame_rotation(Axis::z, era)), rate * transpose(frame_rotation_z_derivative(era)), {0, 0, 0}, {0, 0, 0}};
                         }),
                         "the Terrestrial Intermediate Reference System: CIRS turned by the Earth rotation angle of UT1");
  const int itrs = g.add("ITRS", tirs, Edge([](const time::Instant& t) {
                           const PoleValues p = pole(t);
                           return StateTransform::fixed(transpose(polar_motion(p.xp * ARCSEC, p.yp * ARCSEC, tio_locator(t.tt))));
                         }),
                         "the International Terrestrial Reference System: TIRS with polar motion");
  const int pef = g.add("PEF", tod, Edge([](const time::Instant& t) {
                          const double gst = gast_classical(t.ut1, t.tt);
                          const double rate = earth_rotation_rate(t.eop.lod);
                          return StateTransform{transpose(frame_rotation(Axis::z, gst)), rate * transpose(frame_rotation_z_derivative(gst)), {0, 0, 0}, {0, 0, 0}};
                        }),
                        "the pseudo Earth-fixed frame: TOD turned by GMST + the equation of the equinoxes");
  g.relate("ITRS by the equinox route", gcrs, itrs, Edge([](const time::Instant& t) {
             const PoleValues p = pole(t);
             const Mat3 npb = bias_precession_nutation(t.tt);
             const Mat3 w = polar_motion(p.xp * ARCSEC, p.yp * ARCSEC, tio_locator(t.tt));
             return StateTransform::fixed(w * frame_rotation(Axis::z, gast_classical(t.ut1, t.tt)) * npb);
           }),
           "GCRS to ITRS by W R_z(GMST + EE) N P B, without the CIO series; against the CIO route of the tree");
  g.relate("PEF is TIRS", pef, tirs, Edge([](const time::Instant&) { return StateTransform::fixed(physis::geometry::identity3()); }),
           "the equinox and CIO routes name the same terrestrial intermediate frame; the holonomy is GAST(CIO) − GAST(equinox) and dX, dY");
  return itrs;
}

}  // namespace metrology::frames
