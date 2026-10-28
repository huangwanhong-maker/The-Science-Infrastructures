// Checks of the catalogue of frames (navigation/frames): the Solar System's frames joined in one
// graph, against the ephemeris, NAIF's worked example and closed forms of rotating frames.
//
// The kernels come through the data library (jpl-de440s, naif-lunar-frames) and the Earth
// orientation through iers-eop; without the ephemeris the checks say so and are skipped.

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>

#include <calx/verification/check.hpp>
#include <navigation/frames/catalogue.hpp>

using namespace navigation::frames;
namespace mt = metrology::time;
using calx::verification::Suite;

namespace {
std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(NAVIGATION_ROOT "/../../data/cache")) + "/" + relative;
}
std::string read_text(const std::string& path) {
  std::ifstream in(path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}
}  // namespace

int main() {
  Suite s;
  const std::string spk_path = cache_path("jpl-de440s/de440s.bsp"), pck_path = cache_path("naif-lunar-frames/moon_pa_de440_200625.bpc"),
                    eop_path = cache_path("iers-eop/latest/finals2000A.all");
  if (!std::filesystem::exists(spk_path)) {
    s.skip("the JPL ephemeris DE440s is not in the data library's cache (python data/acquisition/fetch.py jpl-de440s)");
    return s.done();
  }
  Kernels k;
  k.spk = std::make_shared<const astronomy::ephemerides::Spk>(astronomy::ephemerides::Spk::read(spk_path));
  if (std::filesystem::exists(pck_path)) k.moon_pck = std::make_shared<const astronomy::ephemerides::Pck>(astronomy::ephemerides::Pck::read(pck_path));
  std::unique_ptr<metrology::earth_orientation::Series> eop;
  if (std::filesystem::exists(eop_path)) eop = std::make_unique<metrology::earth_orientation::Series>(metrology::earth_orientation::Series::parse_finals2000a(read_text(eop_path)));
  else s.skip("the IERS Earth orientation is not in the cache: UT1 is approximated (python data/acquisition/fetch.py iers-eop)");
  const mt::TimeScales ts(mt::LeapSeconds::built_in(), eop.get(), mt::Ut1Policy::approximate);

  Catalogue c(k);
  const Graph& g = c.graph();
  for (const char* name : {"ICRS", "GCRS", "EME2000", "MOD", "TOD", "CIRS", "TIRS", "ITRS", "PEF", "MOON_ICRS", "IAU_MOON", "SUN_ICRS", "IAU_SUN", "IAU_MARS",
                           "ECLIPTIC_J2000", "GALACTIC", "EARTH_MOON_ROTATING", "TEME"})
    s.require(std::string("the catalogue has ") + name, g.find(name) >= 0);

  // ---- the Earth's centre: the ephemeris -------------------------------------------------------
  const mt::Instant t = ts.instant({2459852.5, 0.0}, mt::Scale::TDB);  // 2022 SEP 30 TDB, NAIF's example
  {
    const auto gcrs_to_icrs = g.transform("GCRS", "ICRS", t);
    const auto earth = k.spk->barycentric(399, t.ephemeris_seconds());
    s.check("the GCRS's origin in the ICRS is the Earth's barycentric position", physis::geometry::norm(gcrs_to_icrs.b - Vec3{earth.position[0], earth.position[1], earth.position[2]}), 1e-6, "m");
    const auto round = physis::geometry::deviation(g.transform("ITRS", "MOON_ICRS", t) * g.transform("MOON_ICRS", "ITRS", t));
    s.check("ITRS → MOON_ICRS → ITRS is the identity, to the rounding of rotations at 4 × 10⁸ m", std::max(round.angle, round.translation), 5e-7);
  }
  // ---- a geostationary point and a point on the Moon's equator ---------------------------------
  {
    const auto itrs_to_gcrs = g.transform("ITRS", "GCRS", t);
    const auto geo = itrs_to_gcrs.state({42164.17e3, 0, 0}, {0, 0, 0});
    s.check("a geostationary point moves at ω r ≈ 3074.7 m/s in the GCRS", physis::geometry::norm(geo[1]) - metrology::frames::earth_rotation_rate(t.eop.lod) * 42164.17e3, 0.01, "m/s");
    const auto moon_iau = g.transform("IAU_MOON", "MOON_ICRS", t).state({1737.4e3, 0, 0}, {0, 0, 0});
    s.check("a point on the Moon's equator moves at about 4.6 m/s about its centre", physis::geometry::norm(moon_iau[1]) - 2 * std::numbers::pi / (27.321661 * 86400) * 1737.4e3, 0.05, "m/s");
  }
  // ---- the Moon's frames and NAIF's example ------------------------------------------------------
  if (!c.has_moon_pa()) {
    s.skip("the Moon's DE440 orientation is not in the cache (python data/acquisition/fetch.py naif-lunar-frames)");
  } else {
    const auto earth_in_me = g.transform("GCRS", "MOON_ME", t).state({0, 0, 0}, {0, 0, 0});
    s.check("the Earth from the Moon in MOON_ME through the catalogue (NAIF's example)",
            physis::geometry::norm(earth_in_me[0] - Vec3{374008.654e3, -23435.970e3, 10141.334e3}), 1.0, "m");
    s.check("and its velocity", physis::geometry::norm(earth_in_me[1] - Vec3{-26.38134, 45.32757, 118.94164}), 0.01, "m/s");
    int relation = -1;
    for (std::size_t i = 0; i < g.relations().size(); i++)
      if (g.relations()[i].name == "MOON_ME by the IAU model") relation = int(i);
    const auto h = physis::geometry::deviation(g.holonomy(relation, t));
    s.check("the holonomy of IAU_MOON against MOON_ME: the IAU model's error (≲ 0.005°)", h.angle * 180 / std::numbers::pi, 0.005, "deg");
    s.require("which is not zero: the IAU model is an approximation", h.angle > 1e-7);
  }
  // ---- TEME ------------------------------------------------------------------------------------------
  {
    int relation = -1;
    for (std::size_t i = 0; i < g.relations().size(); i++)
      if (g.relations()[i].name == "TEME by the equation of the equinoxes") relation = int(i);
    s.require("the catalogue closes a loop through TEME", relation >= 0);
    if (relation >= 0) {
      const auto h = physis::geometry::deviation(g.holonomy(relation, t));
      // GMST 1982 against IAU 2006's GMST at this instant; s′ (tens of µas) and polar motion's order
      // in the two routes make the rest.
      const double expected = std::fabs(std::remainder(metrology::frames::gmst1982(t.ut1) - metrology::frames::gmst(t.ut1, t.tt), 2 * std::numbers::pi));
      s.check("the holonomy through TEME is GMST 1982 − GMST (IAU 2006), within 0.05 mas", (h.angle - expected) / metrology::frames::ARCSEC * 1e3, 0.05, "mas");
      s.require("and it is tens of milliarcseconds, not nothing", h.angle / metrology::frames::ARCSEC * 1e3 > 10);
    }
  }
  // ---- the Earth–Moon rotating frame ----------------------------------------------------------------
  {
    const auto x = g.transform("EARTH_MOON_ROTATING", "GCRS", t);
    const auto moon = k.spk->state(301, 399, t.ephemeris_seconds());
    const Vec3 rm{moon.position[0], moon.position[1], moon.position[2]};
    const double d = physis::geometry::norm(rm);
    s.check("the rotating frame's x points from the Earth to the Moon", physis::geometry::norm(x.direction({1, 0, 0}) - (1 / d) * rm), 1e-15);
    // DE440 places the Earth about the Earth–Moon barycentre by its mass ratio EMRAT = 81.3005682214972154.
    s.check("its origin is the barycentre, at |r_Moon| / (1 + EMRAT) from the Earth (DE440)", physis::geometry::norm(x.b) * (1 + 81.3005682214972154) / d - 1, 1e-9);
    // The Moon stands still on its x axis in the rotating frame, up to the radial part of its motion.
    const auto still = g.transform("MOON_ICRS", "EARTH_MOON_ROTATING", t).state({0, 0, 0}, {0, 0, 0});
    s.check("the Moon lies on the rotating frame's x axis", std::hypot(still[0][1], still[0][2]), 1e-6, "m");
    s.check("and moves only along it (its radial speed)", std::hypot(still[1][1], still[1][2]), 1e-6, "m/s");
    const double w = physis::geometry::norm(x.angular_velocity());
    s.require("the frame turns at about the sidereal month's rate, 2.66 × 10⁻⁶ rad/s (within the eccentricity's ±11 %)", w > 2.3e-6 && w < 3.0e-6);
  }
  // ---- places -------------------------------------------------------------------------------------
  {
    Catalogue c2(k);
    const int gw = c2.add_place("GREENWICH", "earth", {51.4769 * std::numbers::pi / 180, -0.0005 * std::numbers::pi / 180, 46});
    (void)gw;
    const auto up = c2.graph().transform("GREENWICH_ENU", "ITRS", t).direction({0, 0, 1});
    s.check("up at Greenwich is the ellipsoid's normal at 51.4769° N", std::fabs(std::asin(up[2]) * 180 / std::numbers::pi - 51.4769), 1e-9, "deg");
    const auto ned = c2.graph().transform("GREENWICH_NED", "GREENWICH_ENU", t);
    s.check("NED's north is ENU's north", physis::geometry::norm(ned.direction({1, 0, 0}) - Vec3{0, 1, 0}), 1e-15);
    const int crater = c2.add_place("SHACKLETON", "moon", {-89.67 * std::numbers::pi / 180, 129.78 * std::numbers::pi / 180, 0});
    (void)crater;
    const auto o = c2.graph().transform("SHACKLETON_ENU", "MOON_ICRS", t).position({0, 0, 0});
    s.check("a place on the Moon lies on its sphere", physis::geometry::norm(o) - 1737.4e3, 1e-6, "m");
  }
  return s.done();
}
