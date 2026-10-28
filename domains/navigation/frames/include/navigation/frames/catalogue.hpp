// The frames of the Solar System as one frame graph: from the barycentre's ICRS to the Earth's
// terrestrial frame, the Moon's mean-Earth frame and a place on either.
//
// navigation is the bridge where astronomy's ephemerides and rotations, metrology's time scales
// and Earth orientation, and geodesy's figure of the Earth meet; the catalogue joins them in a
// physis frame graph (physis/geometry/frame_graph.hpp) whose time argument is a metrology Instant.
// Each frame is an origin and axes:
//
//   ICRS        the solar-system barycentre, ICRS axes (the root; BCRS)
//   GCRS        the geocentre, ICRS axes (kinematically non-rotating); beneath it metrology's
//               EME2000, MOD, TOD, CIRS, TIRS, ITRS and PEF (earth_frames.hpp)
//   TEME        SGP4's true equator and mean equinox of date, hung from the ITRS by GMST 1982 and
//               polar motion (orbits/teme.hpp); a second route to TOD by the equation of the
//               equinoxes closes a loop whose holonomy is GMST 1982's departure from IAU 2006's
//   SUN_ICRS, MOON_ICRS, <BODY>_ICRS
//               the centre of the Sun, the Moon, a planet (its system's barycentre beyond Mars),
//               ICRS axes
//   MOON_PA, MOON_ME
//               the Moon's DE440 principal-axes and mean-Earth frames (with NAIF's PCK), else
//   IAU_MOON    the IAU model (always present)
//   IAU_<BODY>  the IAU frames of the Sun and the planets
//   EARTH_MOON_ROTATING
//               the Earth–Moon barycentre (DE440's body 3), x from the Earth toward the Moon, z
//               along the Moon's orbital angular momentum about the Earth, y = z × x: the synodic
//               frame of the restricted three-body problem, turning with the Moon (physis's RTN
//               frame of the Moon's orbit, its rate from the ephemeris's acceleration)
//   ECLIPTIC_J2000
//               the barycentre, the ecliptic and equinox of J2000.0 by IAU 2006 (the mean
//               obliquity 84381.406″ on the GCRS's axes with the frame bias; SPICE's ECLIPJ2000
//               uses IAU 1976's 84381.448″ on EME2000, 0.04″ away)
//   GALACTIC    the barycentre, galactic axes (the Hipparcos realisation, as SOFA's iauIcrs2g)
//   GALACTOCENTRIC
//               the Galaxy's centre, x toward it from the Sun's projection on the midplane, z
//               toward the north galactic pole (astropy's Galactocentric; galactocentric.hpp)
//
// and places added with `add_place` (ENU or NED at a geodetic point on the WGS 84 ellipsoid, or on
// the Moon's sphere of 1737.4 km in MOON_ME). Origins are moved by the ephemeris's geometric
// positions, without light time or aberration: these are frames, not apparent directions. The
// relativistic difference between the GCRS and the translated BCRS (a scale of order 10⁻⁸ and
// the geodesic precession of 1.9″ per century) is not modelled; the deviation is stated, not hidden.

#pragma once

#include <cctype>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <astronomy/ephemerides/bodies.hpp>
#include <astronomy/ephemerides/pck.hpp>
#include <astronomy/ephemerides/spk.hpp>
#include <astronomy/rotation/iau.hpp>
#include <astronomy/rotation/moon.hpp>
#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>
#include <metrology/frames/earth_frames.hpp>
#include <navigation/frames/galactocentric.hpp>
#include <navigation/orbits/teme.hpp>
#include <metrology/frames/iau2006.hpp>
#include <metrology/time/scales.hpp>
#include <physis/geometry/frame_graph.hpp>
#include <physis/mechanics/two_body.hpp>

namespace navigation::frames {

using Instant = metrology::time::Instant;
using Graph = physis::geometry::FrameGraph<Instant>;
using Edge = Graph::Edge;
using physis::geometry::Mat3;
using physis::geometry::StateTransform;
using physis::geometry::Vec3;
using namespace physis::geometry::operators;

// The ICRS → galactic matrix (SOFA's iauIcrs2g, from the Hipparcos catalogue's definition).
inline Mat3 icrs_to_galactic() {
  return {{{-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116},
           {+0.494109427875583673525222371358, -0.444829629960011178146614061616, +0.746982244497218890527388004556},
           {-0.867666149019004701181616534570, -0.198076373431201528180486091412, +0.455983776175066922272100478348}}};
}
// The ICRS → ecliptic of J2000.0 matrix, IAU 2006 (SOFA's iauEcm06 at J2000.0: R_x(ε₀) B).
inline Mat3 icrs_to_ecliptic_j2000() {
  return metrology::frames::frame_rotation(metrology::frames::Axis::x, metrology::frames::mean_obliquity({metrology::time::J2000, 0.0})) *
         metrology::frames::frame_bias();
}

// The data a catalogue reads: the ephemeris (required), the Moon's orientation (optional); and the
// parameters of the galactocentric frame, a scenario's choice (astropy's by default).
struct Kernels {
  std::shared_ptr<const astronomy::ephemerides::Spk> spk;
  std::shared_ptr<const astronomy::ephemerides::Pck> moon_pck;  // may be null
  Galactocentric galaxy{};
};

class Catalogue {
 public:
  explicit Catalogue(Kernels k) : k_(std::move(k)) {
    if (!k_.spk) throw std::invalid_argument("navigation: the catalogue needs an ephemeris");
    build();
  }
  const Graph& graph() const { return g_; }
  Graph& graph() { return g_; }
  bool has_moon_pa() const { return k_.moon_pck != nullptr; }
  const Kernels& kernels() const { return k_; }

  // A place on the Earth (WGS 84, in ITRS) or the Moon (sphere of 1737.4 km, in MOON_ME, or
  // IAU_MOON without the PCK), with its ENU frame and, for the Earth, NED.
  int add_place(const std::string& name, const std::string& body, const geodesy::figure::Geodetic& where) {
    if (body == "earth") {
      const auto place = geodesy::figure::Place::at(geodesy::figure::Ellipsoid::wgs84(), where);
      const int enu = g_.add(name + "_ENU", "ITRS", Edge([place](const Instant&) { return place.enu_to_earth(); }),
                             "east-north-up at a place on the WGS 84 ellipsoid");
      g_.add(name + "_NED", "ITRS", Edge([place](const Instant&) { return place.ned_to_earth(); }), "north-east-down at the same place");
      return enu;
    }
    if (body == "moon") {
      const auto place = geodesy::figure::Place::at(moon_sphere(), where);
      return g_.add(name + "_ENU", has_moon_pa() ? "MOON_ME" : "IAU_MOON", Edge([place](const Instant&) { return place.enu_to_earth(); }),
                    "east-north-up at a place on the Moon's sphere of 1737.4 km");
    }
    throw std::invalid_argument("navigation: places on " + body + " are not yet supported");
  }
  static geodesy::figure::Ellipsoid moon_sphere() { return geodesy::figure::Ellipsoid::sphere("Moon, IAU mean sphere", 1737.4e3); }

  // The geometric state of a body's centre relative to another's, in ICRS axes.
  astronomy::ephemerides::State body_state(int target, int observer, const Instant& t) const {
    return k_.spk->state(target, observer, t.ephemeris_seconds());
  }

 private:
  // A frame centred on a body, with ICRS axes: the barycentric state as its translation.
  Edge centred(int naif) const {
    auto spk = k_.spk;
    return Edge([spk, naif](const Instant& t) {
      const auto s = spk->barycentric(naif, t.ephemeris_seconds());
      return StateTransform{physis::geometry::identity3(), physis::geometry::zero3(), {s.position[0], s.position[1], s.position[2]},
                            {s.velocity[0], s.velocity[1], s.velocity[2]}};
    });
  }
  void build() {
    const int icrs = g_.root("ICRS", "the solar-system barycentre with ICRS axes (BCRS)");
    const int gcrs = g_.add("GCRS", icrs, centred(399), "the geocentre with ICRS axes (the GCRS, its relativistic scale not modelled)");
    metrology::frames::add_earth_frames(g_, gcrs);
    // TEME, the frame of SGP4's states: to the ITRS by GMST 1982, the Earth's rate and polar motion.
    const int teme = g_.add("TEME", "ITRS", Edge([](const Instant& t) {
                              constexpr double arcsec = metrology::frames::ARCSEC;
                              return navigation::orbits::teme_to_itrs(t.ut1, t.eop.xp * arcsec, t.eop.yp * arcsec, t.eop.lod);
                            }),
                            "the true equator and mean equinox of date, SGP4's frame (Vallado et al. 2006)");
    // TEME's x axis is the mean equinox on the true equator: TOD's axes turned by the equation of the
    // equinoxes, EE = GAST − GMST (IAU 2006/2000A). The loop through the ITRS measures GMST 1982
    // against IAU 2006's sidereal time.
    g_.relate("TEME by the equation of the equinoxes", teme, g_.require("TOD"), Edge([](const Instant& t) {
                const double ee = metrology::frames::gast_classical(t.ut1, t.tt) - metrology::frames::gmst(t.ut1, t.tt);
                return StateTransform::fixed(metrology::frames::frame_rotation(metrology::frames::Axis::z, -ee));
              }),
              "TOD from TEME by R3(−EE); the holonomy is GMST 1982 − GMST (IAU 2006), tens of milliarcseconds");
    // The Moon hangs from the Earth, by the ephemeris's Moon-relative-to-Earth state: paths between
    // the two never pass through the barycentre, and keep the precision of their own size.
    auto spk = k_.spk;
    const int moon = g_.add("MOON_ICRS", gcrs, Edge([spk](const Instant& t) {
                              const auto s = spk->state(301, 399, t.ephemeris_seconds());
                              return StateTransform{physis::geometry::identity3(), physis::geometry::zero3(), {s.position[0], s.position[1], s.position[2]},
                                                    {s.velocity[0], s.velocity[1], s.velocity[2]}};
                            }),
                            "the Moon's centre with ICRS axes (from the Earth's)");
    g_.add("SUN_ICRS", icrs, centred(10), "the Sun's centre with ICRS axes");
    // The Earth–Moon rotating frame: the RTN frame of the Moon's orbit about the Earth, its rate from
    // the Moon's acceleration (central differences of the ephemeris's velocity over a minute: the
    // Sun's pull turns the orbit's plane, which the central force alone would not), its origin the
    // Earth–Moon barycentre.
    g_.add("EARTH_MOON_ROTATING", gcrs, Edge([spk](const Instant& t) {
             const double et = t.ephemeris_seconds(), h = 60;
             const auto m = spk->state(301, 399, et), mp = spk->state(301, 399, et + h), mm = spk->state(301, 399, et - h);
             const auto b = spk->state(3, 399, et);
             const Vec3 r{m.position[0], m.position[1], m.position[2]}, v{m.velocity[0], m.velocity[1], m.velocity[2]};
             const Vec3 a{(mp.velocity[0] - mm.velocity[0]) / (2 * h), (mp.velocity[1] - mm.velocity[1]) / (2 * h), (mp.velocity[2] - mm.velocity[2]) / (2 * h)};
             StateTransform x = physis::mechanics::orbit_frame(physis::mechanics::OrbitFrame::rtn, r, v, a);
             x.b = {b.position[0], b.position[1], b.position[2]};
             x.b_dot = {b.velocity[0], b.velocity[1], b.velocity[2]};
             return x;
           }),
           "the Earth-Moon barycentre, x toward the Moon, z along its orbital angular momentum (the synodic frame)");
    if (k_.moon_pck) {
      auto pck = k_.moon_pck;
      g_.add("MOON_PA", moon, Edge([pck](const Instant& t) { return astronomy::rotation::icrf_to_moon_pa(*pck, t.ephemeris_seconds()).inverse(); }),
             "the Moon's principal axes of DE440 (NAIF PCK)");
      g_.add("MOON_ME", moon, Edge([pck](const Instant& t) { return astronomy::rotation::icrf_to_moon_me(*pck, t.ephemeris_seconds()).inverse(); }),
             "the Moon's mean-Earth/polar-axis frame of DE440 (DE421-aligned, NAIF frame 31009)");
    }
    g_.add("IAU_MOON", moon, Edge([](const Instant& t) { return astronomy::rotation::icrf_to_iau_frame("moon", t.tdb.value()).inverse(); }),
           "the Moon by the IAU WGCCRE 2015 model");
    for (const auto& b : astronomy::ephemerides::BODIES) {
      if (b.key == "earth" || b.key == "moon") continue;
      const std::string up = upper(b.key);
      int centre;
      if (b.key == "sun") {
        centre = g_.require("SUN_ICRS");
      } else {
        if (!k_.spk->has(b.ephemeris)) continue;
        centre = g_.add(up + "_ICRS", icrs, centred(b.ephemeris), std::string("the centre of ") + std::string(b.name) +
                                                                     (b.ephemeris < 10 ? " (its system's barycentre in DE440)" : "") + ", ICRS axes");
      }
      const std::string key(b.key);
      g_.add("IAU_" + up, centre, Edge([key](const Instant& t) { return astronomy::rotation::icrf_to_iau_frame(key, t.tdb.value()).inverse(); }),
             std::string(b.name) + " by the IAU WGCCRE 2015 model");
    }
    g_.add("ECLIPTIC_J2000", icrs, Edge([](const Instant&) { return StateTransform::fixed(physis::geometry::transpose(icrs_to_ecliptic_j2000())); }),
           "the barycentre, the ecliptic and mean equinox of J2000.0 (IAU 2006)");
    g_.add("GALACTIC", icrs, Edge([](const Instant&) { return StateTransform::fixed(physis::geometry::transpose(icrs_to_galactic())); }),
           "the barycentre, galactic axes (Hipparcos)");
    g_.add("GALACTOCENTRIC", icrs, Edge([gc = k_.galaxy](const Instant&) { return gc.to_icrs(); }),
           "the Galaxy's centre (Sgr A*), x toward it from the Sun's projection on the midplane, z toward the north galactic pole (astropy's Galactocentric)");
    // A second route to the Moon's mean-Earth frame: the IAU model; its holonomy against DE440 is
    // the model's error, tens of metres on the surface.
    if (k_.moon_pck) {
      g_.relate("MOON_ME by the IAU model", g_.require("IAU_MOON"), g_.require("MOON_ME"),
                Edge([](const Instant&) { return StateTransform::fixed(physis::geometry::identity3()); }),
                "IAU_MOON and MOON_ME name the same frame; the holonomy is the IAU model's departure from DE440");
    }
  }
  static std::string upper(std::string_view s) {
    std::string u(s);
    for (char& c : u) c = char(std::toupper(static_cast<unsigned char>(c)));
    return u;
  }

  Kernels k_;
  Graph g_;
};

}  // namespace navigation::frames
