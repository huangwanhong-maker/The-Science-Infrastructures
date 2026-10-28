// Checks of the astronomy domain's C++: the JPL ephemeris DE440 read from its SPK kernel, the
// Moon's DE440 frames from NAIF's binary PCK, and the IAU rotation models.
//
// References: JPL Horizons' barycentric vectors (DE441; tests/fixtures/horizons.txt) and sub-Earth
// points (tests/fixtures/orientation.txt), the worked example of NAIF's lunar frame kernel
// moon_de440_250416.tf (the Earth seen from the Moon on 2022-09-30 TDB in MOON_ME and MOON_PA),
// and the domain's JavaScript (rotation.js) for the IAU models. The kernels come through the data
// library (jpl-de440s, naif-lunar-frames); without them these checks say so and are skipped.

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <astronomy/ephemerides/bodies.hpp>
#include <astronomy/ephemerides/pck.hpp>
#include <astronomy/ephemerides/spk.hpp>
#include <astronomy/rotation/iau.hpp>
#include <astronomy/rotation/moon.hpp>
#include <calx/verification/check.hpp>

using namespace astronomy;
using calx::verification::Suite;
using physis::geometry::Vec3;
using namespace physis::geometry::operators;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double LIGHT = 299792458.0;  // m/s

std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(ASTRONOMY_ROOT "/../../data/cache")) + "/" + relative;
}
std::vector<std::vector<std::string>> rows(const std::string& path, const std::string& kind) {
  std::vector<std::vector<std::string>> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    std::vector<std::string> fields;
    std::string x;
    while (f >> x) fields.push_back(x);
    if (!fields.empty() && fields[0] == kind) out.push_back(fields);
  }
  return out;
}
double num(const std::string& s) { return std::strtod(s.c_str(), nullptr); }
Vec3 to_vec(const ephemerides::Vec3& v) { return {v[0], v[1], v[2]}; }

}  // namespace

int main() {
  Suite s;
  const std::string spk_path = cache_path("jpl-de440s/de440s.bsp");
  const std::string pck_path = cache_path("naif-lunar-frames/moon_pa_de440_200625.bpc");
  const bool have_spk = std::filesystem::exists(spk_path), have_pck = std::filesystem::exists(pck_path);

  // ---- the IAU models, which need no data ----------------------------------------------------------
  {
    const auto earth = rotation::pole_and_meridian("earth", 2451545.0);
    s.check("the Earth's IAU prime meridian at J2000.0 is 190.147°", earth->w - 190.147, 1e-12, "deg");
    const auto moon = rotation::pole_and_meridian("moon", 2451545.0);
    s.require("the Moon's IAU model is there", moon.has_value());
    const auto m = rotation::icrf_to_iau_frame("mars", 2460000.5);
    s.check("an IAU frame is a rotation", physis::geometry::orthogonality_defect(m.r), 1e-15);
    // The prime meridian's rate: Ṙ against a central difference of a minute.
    const auto a = rotation::icrf_to_iau_frame("jupiter", 2460000.5 + 30 / 86400.0), b = rotation::icrf_to_iau_frame("jupiter", 2460000.5 - 30 / 86400.0);
    const auto j = rotation::icrf_to_iau_frame("jupiter", 2460000.5);
    // A central difference over 60 s of a frame turning at 870° a day errs by (h²/6) ω³ ≈ 8 × 10⁻¹⁰.
    s.check("the IAU frame's rotation rate against differences", physis::geometry::max_difference((1.0 / 60) * (a.r - b.r), j.r_dot), 2e-9, "1/s");
    s.require("no model for an unknown body", !rotation::pole_and_meridian("vulcan", 2451545.0).has_value());
  }

  if (!have_spk) {
    s.skip("the JPL ephemeris DE440s is not in the data library's cache (python data/acquisition/fetch.py jpl-de440s)");
    return s.done();
  }
  const ephemerides::Spk spk = ephemerides::Spk::read(spk_path);

  // ---- DE440 against Horizons ---------------------------------------------------------------------
  {
    const auto cov = spk.coverage();
    s.require("DE440s covers 1850 to 2149", cov.first < -4.7e9 && cov.second > 4.7e9);
    double worst_planet = 0, worst_velocity = 0, worst_moon = 0;
    int n = 0;
    for (const auto& r : rows(ASTRONOMY_ROOT "/tests/fixtures/horizons.txt", "ssb")) {
      const int id = std::atoi(r[1].c_str());
      const double et = (num(r[2]) - 2451545.0) * 86400.0;
      if (!spk.covers(et) || !spk.has(id)) continue;
      const auto st = spk.barycentric(id, et);
      const Vec3 p{num(r[3]) * 1e3, num(r[4]) * 1e3, num(r[5]) * 1e3}, v{num(r[6]) * 1e3, num(r[7]) * 1e3, num(r[8]) * 1e3};
      const double dp = physis::geometry::norm(to_vec(st.position) - p), dv = physis::geometry::norm(to_vec(st.velocity) - v);
      if (id == 301) worst_moon = std::max(worst_moon, dp);
      else worst_planet = std::max(worst_planet, dp);
      worst_velocity = std::max(worst_velocity, dv);
      n++;
    }
    s.require("Horizons' vectors were read", n > 20);
    // DE440 against Horizons' DE441: the two agree to metres over the span of DE440s.
    s.check("barycentric positions of the Sun and planets against Horizons (DE441)", worst_planet, 2e3, "m");
    s.check("the Moon's barycentric position against Horizons", worst_moon, 50, "m");
    s.check("barycentric velocities against Horizons", worst_velocity, 1e-3, "m/s");
    // The Earth-Moon distance stays between perigee and apogee.
    const auto em = spk.state(301, 399, 717768000.0);
    const double d = physis::geometry::norm(to_vec(em.position));
    s.require("the Moon is between 356 000 and 407 000 km from the Earth", d > 356e6 && d < 407e6);
  }

  // ---- the Moon's frames: NAIF's worked example ----------------------------------------------------
  if (!have_pck) {
    s.skip("the Moon's DE440 orientation is not in the cache (python data/acquisition/fetch.py naif-lunar-frames)");
  } else {
    const ephemerides::Pck pck = ephemerides::Pck::read(pck_path);
    const double et = 717768000.0;  // 2022 SEP 30 TDB
    const auto earth = spk.state(399, 301, et);
    const auto pa = rotation::icrf_to_moon_pa(pck, et).state(to_vec(earth.position), to_vec(earth.velocity));
    const auto me = rotation::icrf_to_moon_me(pck, et).state(to_vec(earth.position), to_vec(earth.velocity));
    // moon_de440_250416.tf, "Program output": the Earth relative to the Moon, km and km/s.
    const Vec3 pa_naif{373997.028e3, -23558.987e3, 10284.057e3}, me_naif{374008.654e3, -23435.970e3, 10141.334e3};
    const Vec3 vpa_naif{-0.02641181e3, 0.04533642e3, 0.11893151e3}, vme_naif{-0.02638134e3, 0.04532757e3, 0.11894164e3};
    s.check("the Earth from the Moon in MOON_PA (NAIF's example, printed to 1 m)", physis::geometry::norm(pa[0] - pa_naif), 1.0, "m");
    s.check("the Earth from the Moon in MOON_ME (NAIF's example)", physis::geometry::norm(me[0] - me_naif), 1.0, "m");
    s.check("its velocity in MOON_PA, with the frame's rotation", physis::geometry::norm(pa[1] - vpa_naif), 1e-5 * 1e3, "m/s");
    s.check("its velocity in MOON_ME", physis::geometry::norm(me[1] - vme_naif), 1e-5 * 1e3, "m/s");
    // The rotation of PA against central differences (the frame's derivative).
    const auto plus = rotation::icrf_to_moon_pa(pck, et + 10), minus = rotation::icrf_to_moon_pa(pck, et - 10);
    s.check("the PA frame's rotation rate against differences", physis::geometry::max_difference((1.0 / 20) * (plus.r - minus.r), rotation::icrf_to_moon_pa(pck, et).r_dot), 1e-13, "1/s");
    // The sub-Earth point in MOON_ME against Horizons (light time from the Earth).
    double worst = 0, iau_worst = 0;
    int n = 0;
    for (const auto& r : rows(ASTRONOMY_ROOT "/tests/fixtures/orientation.txt", "sub")) {
      if (r[1] != "moon") continue;
      const double jd_tt = num(r[2]);
      const double et_tt = (jd_tt - 2451545.0) * 86400.0;  // TDB − TT (2 ms) is negligible here
      double tau = 0;
      Vec3 g{};
      for (int i = 0; i < 3; i++) {
        const auto e = spk.barycentric(399, et_tt), m = spk.barycentric(301, et_tt - tau);
        g = to_vec(e.position) - to_vec(m.position);
        tau = physis::geometry::norm(g) / LIGHT;
      }
      const Vec3 v = rotation::icrf_to_moon_me(pck, et_tt - tau).r * g;
      const double lon = std::atan2(v[1], v[0]) * 180 / pi, lat = std::asin(v[2] / physis::geometry::norm(v)) * 180 / pi;
      worst = std::max({worst, std::fabs(std::remainder(lon - num(r[3]), 360.0)), std::fabs(lat - num(r[4]))});
      const Vec3 vi = rotation::icrf_to_iau_frame("moon", jd_tt - tau / 86400).r * g;
      const double loni = std::atan2(vi[1], vi[0]) * 180 / pi, lati = std::asin(vi[2] / physis::geometry::norm(vi)) * 180 / pi;
      iau_worst = std::max({iau_worst, std::fabs(std::remainder(loni - num(r[3]), 360.0)), std::fabs(lati - num(r[4]))});
      n++;
    }
    s.require("Horizons' sub-Earth points on the Moon were read", n >= 4);
    s.check("the sub-Earth point in MOON_ME (DE440) against Horizons", worst, 2e-4, "deg");
    s.check("the same by the IAU model, which approximates MOON_ME", iau_worst, 0.02, "deg");
  }
  return s.done();
}
