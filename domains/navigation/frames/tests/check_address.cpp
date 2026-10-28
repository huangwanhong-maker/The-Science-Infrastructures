// Checks of addresses (navigation/frames/address.hpp) and of the galactocentric frame.
//
// The language: canonical words read back to themselves (print ∘ parse = id on canonical words),
// any address prints and reads back to itself (parse ∘ print = id), malformed words fail at the
// column where they go wrong, saying what was expected. The meaning: every chart's coordinates
// resolve and are answered back (geodetic, planetocentric, UTM, MGRS, spherical, cylindrical,
// Cartesian), steps from the horizon are seen again from it, a direction is carried without the
// translations; the Galaxy's centre named in the ICRS is the galactocentric origin; and the frame
// carries 300 stars' positions and velocities as astropy does (tests/fixtures/
// galactocentric-reference.txt, made by tests/make_frames_fixtures.py with astropy 8.0.2). The
// meaning needs the ephemeris (jpl-de440s) and is skipped without it; the language does not.

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <navigation/frames/address.hpp>

using namespace navigation::frames;
namespace mt = metrology::time;
using calx::verification::Suite;

namespace {
std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(NAVIGATION_ROOT "/../../data/cache")) + "/" + relative;
}
// The column of a syntax error, or −1 when the text parses.
long error_column(const std::string& text) {
  try {
    parse_address(text);
  } catch (const SyntaxError& e) {
    return long(e.column());
  }
  return -1;
}
bool resolve_fails(const Catalogue& c, const std::string& text) {
  try {
    resolve(parse_address(text), c);
  } catch (const ResolveError&) {
    return true;
  }
  return false;
}
double angle_gap(double a, double b) {  // degrees
  double d = std::fmod(std::fabs(a - b), 360.0);
  return std::min(d, 360.0 - d);
}
}  // namespace

int main() {
  Suite s;
  constexpr double deg = std::numbers::pi / 180;

  // ---- the language -------------------------------------------------------------------------------
  {
    const std::vector<std::string> canonical = {
        "earth(48.8583701° N, 2.2944813° E, 330 m)",
        "ITRS:geodetic(48.8583701° N, 2.2944813° E, 330 m) @ 2026-10-28T12:00:00Z",
        "ITRS:mgrs(\"31UDQ4825211954\", 330 m)",
        "ITRS:utm(31 N, 448252 m, 5411954 m)",
        "moon(43.31° S, 11.36° W)",
        "earth(35.36° N, 138.73° E)/aer(45°, 30°, 10 km)",
        "earth(35.36° N, 138.73° E)/enu(10 m, -20 m, 3.5 m)/ned(1 m, 2 m, 3 m)",
        "ICRS:equatorial(266.4051°, -28.936175°, 8.122 kpc)",
        "GALACTOCENTRIC:cylindrical(8.2 kpc, 180°, 20.8 pc)",
        "mars",
        "galaxy",
        "ssb @ now",
        "ICRS:equatorial(6.75247 h, -16.7161°)",
        "ICRS(1 au, 0 m, -2.5e+20 m) @ JD 2461000.5 TDB",
        "GALACTIC:galactic(120°, 1.5°, 3.2 ly)",
    };
    std::size_t same = 0;
    std::string first_bad;
    for (const auto& w : canonical) {
      const std::string again = print_address(parse_address(w));
      if (again == w) same++;
      else if (first_bad.empty()) first_bad = w + " → " + again;
    }
    s.require("print ∘ parse is the identity on canonical words (15 of them)", same == canonical.size(), first_bad);

    // Words that are not canonical: they read, print canonically, and read back to the same address.
    const std::vector<std::string> loose = {
        "Earth( 48°51′30.13″ N , 2°17′40.13″E )",
        "itrs:GEODETIC(1 deg N, 2 deg E)",
        "icrs:equatorial(6 h 45 m 8.9 s, -16° 42′ 58″)",
        "Moon(11.36 W, 43.31 S)",
        "ecliptic_j2000:ecliptic(1.5 rad, -0.25 rad, 1 au)",
        "earth(0 N, 0 E)/AER(90 deg, 0 deg, 1000 km)",
        "GCRS:cartesian(+6378137, 0, 0)",
        "ITRS:geodetic(10 N, 20 E) @   2026-01-01T00:00:00 TT  ",
    };
    std::size_t back = 0;
    for (const auto& w : loose) {
      const Address a = parse_address(w);
      const std::string p = print_address(a);
      if (parse_address(p) == a && print_address(parse_address(p)) == p) back++;
    }
    s.require("parse ∘ print is the identity on addresses, and canonical printing is idempotent (8 loose words)", back == loose.size());
    const Address dms = parse_address("Earth( 48°51′30.13″ N , 2°17′40.13″E )");
    s.check("degrees, minutes and seconds of arc add up", std::fabs(dms.coords[0].value - (48 + 51.0 / 60 + 30.13 / 3600)), 1e-12, "deg");
    const Address hms = parse_address("ICRS:equatorial(6 h 45 m 8.9 s, -16° 42′ 58″)");
    s.check("hours, minutes and seconds of time add up, and a negative sign takes the whole angle",
            std::fabs(hms.coords[0].value - (6 + 45.0 / 60 + 8.9 / 3600)) + std::fabs(hms.coords[1].value + (16 + 42.0 / 60 + 58.0 / 3600)), 1e-12);
    s.require("an instant is the text after @, trimmed", parse_address("ITRS:geodetic(10 N, 20 E) @   2026-01-01T00:00:00 TT  ").instant == "2026-01-01T00:00:00 TT");
    s.require("\"2.35E\" is 2.35 and the hemisphere E, not an exponent", parse_address("earth(1N, 2.35E)").coords[1].hemisphere == 'E' &&
                                                                          parse_address("earth(1N, 2.35E)").coords[1].value == 2.35);

    // Malformed words fail where they go wrong (columns from 0, in bytes).
    s.require("a chart without coordinates fails at the number (column 14)", error_column("ITRS:geodetic 48 N") == 14);
    s.require("an unclosed list fails at the end (column 15)", error_column("earth(48 N, 2 E") == 15);
    s.require("a body with a chart fails at the colon (column 5)", error_column("earth:geodetic(1 N, 2 E)") == 5);
    s.require("an empty instant fails after the @ (column 33)", error_column("ICRS:equatorial(1°, 2°, 3 pc) @ ") == 33);
    s.require("text after the point fails at it (column 15)", error_column("moon(1 N, 2 E) junk") == 15);
    s.require("an unclosed quotation fails at it (column 10)", error_column("ITRS:mgrs(\"31UDQ") == 10);
    s.require("an unknown step fails at its name (column 11)", error_column("earth(1,2)/xyz(1,2,3)") == 11);
    bool says = false;
    try { parse_address("earth(48 N, 2 E"); } catch (const SyntaxError& e) { says = e.expected().find(")") != std::string::npos; }
    s.require("and says what was expected", says);
  }

  // ---- the meaning ---------------------------------------------------------------------------------------
  const std::string spk_path = cache_path("jpl-de440s/de440s.bsp"), pck_path = cache_path("naif-lunar-frames/moon_pa_de440_200625.bpc");
  if (!std::filesystem::exists(spk_path)) {
    s.skip("the JPL ephemeris DE440s is not in the data library's cache (python data/acquisition/fetch.py jpl-de440s): addresses are not resolved");
    return s.done();
  }
  Kernels k;
  k.spk = std::make_shared<const astronomy::ephemerides::Spk>(astronomy::ephemerides::Spk::read(spk_path));
  if (std::filesystem::exists(pck_path)) k.moon_pck = std::make_shared<const astronomy::ephemerides::Pck>(astronomy::ephemerides::Pck::read(pck_path));
  const mt::TimeScales ts(mt::LeapSeconds::built_in(), nullptr, mt::Ut1Policy::approximate);
  const Catalogue c(k);
  const mt::Instant t = ts.instant({2461342.0, 0.0}, mt::Scale::TT);  // 2026-10-28T12:00 TT

  {
    const Event e = resolve(parse_address("earth(48.8583701° N, 2.2944813° E, 330 m)"), c);
    s.require("a place on the Earth is in the ITRS, on the Earth", e.frame == "ITRS" && e.body == "earth");
    const auto x = express(e, "ITRS", "", c, t);
    s.require("the ITRS answers in geodetic coordinates by default", x.chart == "geodetic" && x.values.size() == 3);
    s.check("geodetic coordinates come back", std::fabs(x.values[0].value - 48.8583701) + std::fabs(x.values[1].value - 2.2944813), 1e-11, "deg");
    s.check("and the height", std::fabs(x.values[2].value - 330), 1e-7, "m");
    const auto u = express(e, "ITRS", "utm", c, t);
    const Event eu = resolve(u.address, c);
    s.check("UTM: answered and read back", physis::geometry::norm(eu.position - e.position), 1e-6, "m");
    const auto m = express(e, "ITRS", "mgrs", c, t);
    const Event em = resolve(m.address, c);
    s.require("the Eiffel Tower's MGRS reference begins 31UDQ4825", m.address.coords[0].word.rfind("31UDQ4825", 0) == 0, m.address.coords[0].word);
    s.check("MGRS at 1 m: read back within its square's diagonal", physis::geometry::norm(em.position - e.position), 1.5, "m");
    // The same place written in a frame of the sky and read back.
    const auto sky = express(e, "ICRS", "", c, t);
    const Event back = carry(resolve(sky.address, c), "ITRS", c, t);
    s.check("an Earth place through the ICRS (1 au away) and back", physis::geometry::norm(back.position - e.position), 1e-4, "m");
  }
  {
    const Event fuji = resolve(parse_address("earth(35.36° N, 138.73° E)"), c);
    const Event p = resolve(parse_address("earth(35.36° N, 138.73° E)/aer(45°, 30°, 10 km)"), c);
    const auto seen = express_from(p, fuji, "aer", c, t);
    s.check("a step by azimuth, elevation and range is seen again from the observer",
            angle_gap(seen.values[0].value, 45) + std::fabs(seen.values[1].value - 30) + std::fabs(seen.values[2].value - 1e4) * 1e-3, 1e-9);
    const Event a = resolve(parse_address("earth(35.36 N, 138.73 E)/enu(10 m, -20 m, 3.5 m)"), c);
    const Event b = resolve(parse_address("earth(35.36 N, 138.73 E)/ned(-20 m, 10 m, -3.5 m)"), c);
    s.check("enu and ned name the same point", physis::geometry::norm(a.position - b.position), 1e-9, "m");
    const auto again = parse_address(print_address(seen.address));
    s.check("the answer's address reads back to the point", physis::geometry::norm(resolve(again, c).position - p.position), 1e-6, "m");
  }
  {
    const Event tycho = resolve(parse_address("moon(43.31° S, 11.36° W)"), c);
    s.require("a place on the Moon is in its mean-Earth frame", tycho.frame == (c.has_moon_pa() ? "MOON_ME" : "IAU_MOON") && tycho.body == "moon");
    const auto x = express(tycho, tycho.frame, "", c, t);
    s.check("selenographic coordinates come back", std::fabs(x.values[0].value + 43.31) + std::fabs(x.values[1].value + 11.36), 1e-11, "deg");
    const auto via = carry(carry(tycho, "ITRS", c, t), tycho.frame, c, t);
    s.check("Tycho through the Earth's frame and back", physis::geometry::norm(via.position - tycho.position), 1e-4, "m");
    const Event mars = resolve(parse_address("mars(0 N, 0 E)"), c);
    s.check("a place on Mars is on its IAU mean sphere of 3389.5 km", std::fabs(physis::geometry::norm(mars.position) - 3389.5e3), 1e-6, "m");
    const Event centre = resolve(parse_address("mars"), c);
    s.require("Mars alone is its centre", centre.frame == "MARS_ICRS" && physis::geometry::norm(centre.position) == 0);
  }
  {
    // The Galaxy's centre, as the ICRS sees it, is the galactocentric origin.
    const Event gc = carry(resolve(parse_address("ICRS:equatorial(266.4051°, -28.936175°, 8.122 kpc)"), c), "GALACTOCENTRIC", c, t);
    s.check("Sgr A* at 8.122 kpc in the ICRS is the galactocentric origin (relative)", physis::geometry::norm(gc.position) / (8.122e3 * PARSEC_M), 1e-14);
    const Event g = carry(resolve(parse_address("galaxy"), c), "ICRS", c, t);
    s.check("the Galaxy's centre is 8.122 kpc from the barycentre", std::fabs(physis::geometry::norm(g.position) / PARSEC_M - 8122), 1e-9, "pc");
    const auto sun = express(resolve(parse_address("ssb"), c), "GALACTOCENTRIC", "cartesian", c, t);
    s.check("the barycentre sits 20.8 pc above the midplane", std::fabs(sun.values[2].value / PARSEC_M - 20.8), 1e-9, "pc");
    // astropy's 300 stars.
    std::ifstream in(NAVIGATION_FIXTURES "/galactocentric-reference.txt");
    std::string line;
    double dx = 0, dv = 0;
    int n = 0;
    const auto& graph = c.graph();
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      Vec3 x, v, X, V;
      f >> x[0] >> x[1] >> x[2] >> v[0] >> v[1] >> v[2] >> X[0] >> X[1] >> X[2] >> V[0] >> V[1] >> V[2];
      const auto st = graph.transform("ICRS", "GALACTOCENTRIC", t).state(x, v);
      dx = std::max(dx, physis::geometry::norm(st[0] - X) / physis::geometry::norm(X));
      dv = std::max(dv, physis::geometry::norm(st[1] - V));
      n++;
    }
    s.require("astropy's stars were read", n == 300);
    s.check("galactocentric positions against astropy (relative)", dx, 1e-13);
    s.check("galactocentric velocities against astropy", dv, 1e-6, "m/s");
  }
  {
    // A direction is carried by the rotations alone.
    const Event sirius = resolve(parse_address("ICRS:equatorial(6 h 45 m 8.9 s, -16° 42′ 58″)"), c);
    s.require("an equatorial position without a distance is a direction", sirius.direction);
    const auto in_gcrs = express(sirius, "GCRS", "", c, t);
    s.check("the direction is the same from the geocentre (no translation)", angle_gap(in_gcrs.values[0].value, (6 + 45.0 / 60 + 8.9 / 3600) * 15) +
                                                                            std::fabs(in_gcrs.values[1].value + (16 + 42.0 / 60 + 58.0 / 3600)), 1e-11, "deg");
    const auto gal = express(sirius, "GALACTIC", "", c, t);
    const Vec3 direct = icrs_to_galactic() * sirius.position;
    s.check("and in galactic axes it is the matrix's", angle_gap(gal.values[0].value, std::atan2(direct[1], direct[0]) / deg) +
                                                         std::fabs(gal.values[1].value - std::asin(direct[2]) / deg), 1e-11, "deg");
    s.require("a direction has no Cartesian coordinates", [&] { try { express(sirius, "ICRS", "cartesian", c, t); } catch (const ResolveError&) { return true; } return false; }());
  }
  {
    s.require("a chart the frame does not have is refused", resolve_fails(c, "ITRS:ecliptic(1°, 2°)"));
    s.require("a frame the catalogue does not have is refused", resolve_fails(c, "NOWHERE(1 m, 2 m, 3 m)"));
    s.require("a step from a body's centre is refused", resolve_fails(c, "mars/enu(1 m, 2 m, 3 m)"));
    s.require("a latitude beyond 90° is refused", resolve_fails(c, "earth(91 N, 0 E)"));
    s.require("a length given in degrees is refused", resolve_fails(c, "ICRS(1°, 2 m, 3 m)"));
    s.require("charts: the ITRS's first is geodetic, the ICRS's equatorial, GALACTOCENTRIC's cartesian",
              charts_of(c, "ITRS").front() == "geodetic" && charts_of(c, "ICRS").front() == "equatorial" &&
                  charts_of(c, "GALACTOCENTRIC").front() == "cartesian" && charts_of(c, "IAU_MARS").front() == "planetocentric");
  }
  return s.done();
}
