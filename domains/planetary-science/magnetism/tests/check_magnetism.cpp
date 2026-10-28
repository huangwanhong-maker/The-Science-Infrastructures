// Checks of the magnetism module (planetary-science/magnetism): the Moon's crustal field and the
// interplanetary field carried to it.
//
// References:
// - pyshtools 4.14 (Wieczorek and Meschede 2018) on T2015_449 at five mounts on and above the Moon,
//   tests/fixtures/lunar-field.txt, written by tests/make_fixtures.py from the data library's file;
// - closed forms: a dipole's field, the geometry of the Sun–Earth line, GSE's axes.
// The model comes from the data library's cache (moon-magnetic-t2015/T2015_449.sh, or under
// SCIENCE_CACHE); without it only the checks of the closed forms run.

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <planetary_science/magnetism/crust.hpp>
#include <planetary_science/magnetism/solar_wind.hpp>

using calx::verification::Suite;
namespace pm = planetary_science::magnetism;
namespace pa = physis::analysis;

namespace {

std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(PLANETARY_SCIENCE_ROOT "/../../data/cache")) + "/" + relative;
}

double norm(const pm::Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
double dot(const pm::Vec3& a, const pm::Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

}  // namespace

int main() {
  Suite s;
  const double deg = std::acos(-1.0) / 180;

  // ---- the carriage from L1 -------------------------------------------------------------------------------------------
  {
    const pm::Vec3 sun{1.496e11, 0, 0};
    const double v = 4e5, moon = 3.844e8;
    const pm::Carriage day = pm::carry_from_l1(sun, {moon, 0, 0}, v);
    const pm::Carriage night = pm::carry_from_l1(sun, {-moon, 0, 0}, v);
    const pm::Carriage beside = pm::carry_from_l1(sun, {-moon, 1.7e8, 0}, v);
    s.check("carriage: a body sunward of the Earth sees the field of (x_L1 − x∥)/v before", std::fabs(day.delay_s - (1.5e9 - moon) / v), 1e-9, "s");
    s.check("carriage: on the night side, a longer travel", std::fabs(night.delay_s - (1.5e9 + moon) / v), 1e-9, "s");
    s.require("carriage: on the anti-Sun line, inside the magnetotail; sunward, not", night.in_magnetotail && !day.in_magnetotail);
    s.require("carriage: on the night side 1.7e8 m from the line, beyond the tail's 25 Earth radii (1.59e8 m)",
              !beside.in_magnetotail && std::fabs(beside.off_axis_m - 1.7e8) < 1e-6);
  }

  // ---- GSE's axes -----------------------------------------------------------------------------------------------------
  {
    const double e = 23.4392911 * deg;
    const pm::Vec3 pole{0, -std::sin(e), std::cos(e)}, sun{0.3e11, 1.2e11, 0.5e11};
    const pm::GseAxes a = pm::gse_axes(sun, pole);
    const pm::Vec3 xy{a.x[1] * a.y[2] - a.x[2] * a.y[1], a.x[2] * a.y[0] - a.x[0] * a.y[2], a.x[0] * a.y[1] - a.x[1] * a.y[0]};
    const double ortho = std::max({std::fabs(dot(a.x, a.y)), std::fabs(dot(a.y, a.z)), std::fabs(dot(a.z, a.x)), std::fabs(norm(a.x) - 1),
                                   std::fabs(norm(a.y) - 1), std::fabs(norm(a.z) - 1), norm({xy[0] - a.z[0], xy[1] - a.z[1], xy[2] - a.z[2]})});
    s.check("GSE: orthonormal and right-handed (x × y = z)", ortho, 1e-15);
    s.require("GSE: x toward the Sun, y in the ecliptic, z on the pole's side",
              std::fabs(dot(a.x, sun) / norm(sun) - 1) < 1e-15 && std::fabs(dot(a.y, pole)) < 1e-15 && dot(a.z, pole) > 0);
    const pm::Vec3 b{3e-9, -2e-9, 5e-9}, f = pm::from_gse(a, b);
    s.check("GSE: a vector carried to the frame reads back by its components", norm({dot(f, a.x) - b[0], dot(f, a.y) - b[1], dot(f, a.z) - b[2]}), 1e-23, "T");
  }

  // ---- the nearest sample, and the local horizon ------------------------------------------------------------------------
  {
    const std::vector<std::array<double, 4>> series{{100, 1, 0, 0}, {160, 2, 0, 0}, {220, 3, 0, 0}};
    const auto a = pm::nearest_sample(series, 175, 300), b = pm::nearest_sample(series, 900, 300), c = pm::nearest_sample({}, 0, 300);
    s.require("series: the nearest sample within the gap; none beyond it, none of an empty series", a && (*a)[1] == 2 && !b && !c);
    const pm::Vec3 up = pm::north_east_down({1, 0, 0}, 0, 0), north = pm::north_east_down({0, 0, 1}, 0, 0), east = pm::north_east_down({0, 1, 0}, 0, 0);
    s.require("horizon: at (0, 0) the x axis is up, z north, y east",
              up[2] == -1 && std::fabs(up[0]) < 1e-16 && north[0] == 1 && east[1] == 1);
  }

  // ---- low ground: a dipole -----------------------------------------------------------------------------------------------
  {
    pa::HarmonicPotential p = pa::HarmonicPotential::zero(1, 1737400);
    p.g[pa::HarmonicPotential::index(1, 0)] = 1e-9;
    const pm::CrustReading low = pm::crust_at(p, 1737400 - 2000, 0.3, 1.1), at = pm::crust_at(p, 1737400, 0.3, 1.1);
    const pm::CrustReading high = pm::crust_at(p, 1737400 + 5000, 0.3, 1.1);
    const pa::FieldValue direct = pa::field_at(p, 1737400 + 5000, 0.3, 1.1);
    s.require("low ground: 2 km below the sphere the model is read at the sphere, and the depth is said",
              low.field.down == at.field.down && low.field.north == at.field.north && low.radius_read_m == 1737400 && std::fabs(low.below_sphere_m - 2000) < 1e-9);
    s.require("low ground: above the sphere the reading is the expansion's own", high.field.down == direct.down && high.below_sphere_m == 0);
  }

  // ---- T2015_449 ----------------------------------------------------------------------------------------------------------
  const std::string path = cache_path("moon-magnetic-t2015/T2015_449.sh");
  if (!std::filesystem::exists(path)) {
    s.skip("T2015_449 is not in the data library's cache (python data/acquisition/fetch.py moon-magnetic-t2015)");
    return s.done();
  }
  const auto f = pm::read_lunar_t2015(path);
  s.require("T2015_449: degree 449 on the sphere of 1737.4 km", f->model().degree == 449 && f->model().radius == 1737400);
  {
    std::ifstream in(PLANETARY_SCIENCE_FIXTURES "/lunar-field.txt");
    std::string line;
    double worst = 0;
    int n = 0;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream row(line);
      std::string name;
      double lat, lon, h, bn, be, bd;
      row >> name >> lat >> lon >> h >> bn >> be >> bd;
      const pm::CrustReading r = pm::crust_at(f->model(), f->model().radius + h, lat * deg, lon * deg);
      const double b = std::sqrt(bn * bn + be * be + bd * bd);
      worst = std::max({worst, std::fabs(r.field.north - bn) / b, std::fabs(r.field.east - be) / b, std::fabs(r.field.down - bd) / b});
      n++;
    }
    s.require("T2015_449: five mounts in the fixture", n == 5);
    s.check("T2015_449: the crust's field at five mounts on and above the Moon as pyshtools gives it (of |B|)", worst, 1e-10);
  }
  {
    const double depth = f->depth_m();
    s.require("the law: its power falls at high degree (ρ < 1), from sources 1 to 20 km deep in effect", f->law_bounded() && depth > 1e3 && depth < 20e3);
    s.check("the law: it fits the spectrum on degrees 30–250 within a quarter in its natural logarithm (rms)", f->law().rms_log_residual, 0.25);
    // The unresolved part above the model's degree carries the law's power: the mean of Rₙ/R̂ₙ over
    // degrees 450–720, each a χ² of 2n + 1 degrees of freedom over its mean, within four standard errors.
    const auto u = f->unresolved(720, 7);
    const std::vector<double> R = pa::lowes_spectrum(*u);
    double mean = 0, var = 0;
    for (int k = 450; k <= 720; k++) mean += R[k] / f->law()(k), var += 2.0 / (2 * k + 1);
    mean /= 271, var /= 271.0 * 271.0;
    s.check("the unresolved part: its power above degree 449 is the law's, on average", std::fabs(mean - 1), 4 * std::sqrt(var));
    // At degree 500 the model has no power, so each seed draws its own coefficient there.
    const std::size_t k500 = pa::HarmonicPotential::index(500, 3);
    s.require("the unresolved part: one draw for a degree and a seed, another for another seed",
              f->unresolved(720, 7) == u && u->g[k500] != 0 && f->unresolved(720, 8)->g[k500] != u->g[k500]);
  }
  return s.done();
}
