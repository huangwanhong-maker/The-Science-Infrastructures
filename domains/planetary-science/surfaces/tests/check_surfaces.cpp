// Checks of the Moon's ground as a robot's foot meets it (surfaces/regolith.hpp, ADR 0059).
//
// - The presumption: the friction's 95 % interval is tan 30° to tan 50° (the regolith's angles of internal
//   friction, Lunar Sourcebook ch. 9); the load at 1 mm of the median is the rabbit's contact as built.
// - The realisation: over 2000 places a kilometre apart, each parameter's z = log(value/median)/spread has
//   mean 0 and deviation 1 (within four standard errors); two feet 5 cm apart stand on nearly the same
//   ground and two 4 m apart on independent ground; a seed fixes it and another seed changes it.
// - As fields: the presumed field is the median everywhere with provenance model, the realised field is the
//   realisation (synthetic), and a point far above the ground is refused with the reason.

#include <planetary_science/surfaces/regolith.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>

namespace ps = planetary_science::surfaces;
namespace pf = physis::fields;

namespace {

ps::RegolithGround::Point on_moon(double lat, double lon) {
  const double r = ps::MOON_RADIUS_M;
  return {r * std::cos(lat) * std::cos(lon), r * std::cos(lat) * std::sin(lon), r * std::sin(lat)};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const ps::RegolithLaw law;
  const double deg = std::numbers::pi / 180;
  s.check("presumption: the friction's 95 % interval begins at tan 30° (relative)", law.friction.at(-1.959964) / std::tan(30 * deg) - 1, 1e-6);
  s.check("presumption: and ends at tan 50° (relative)", law.friction.at(1.959964) / std::tan(50 * deg) - 1, 1e-6);
  s.check("presumption: the median load at 1 mm is the rabbit's contact as built, k = 3 × 10⁵ N/m^{3/2} (relative)", law.median().stiffness() / 3e5 - 1, 2e-3);

  const auto ground = std::make_shared<const ps::RegolithGround>(law, 7);
  {
    // 2000 places on a line of latitude, 1 km apart: independent draws of each parameter's z.
    const int n = 2000;
    double s1[3] = {}, s2[3] = {};
    for (int i = 0; i < n; i++) {
      const auto x = on_moon(7.5 * deg, (-59.0 + i * 1000 / (ps::MOON_RADIUS_M * std::cos(7.5 * deg)) / deg) * deg);
      const auto c = ground->at(x);
      const double z[3] = {std::log(c.friction / law.friction.median) / law.friction.log_sd, std::log(c.load_1mm_n / law.load.median) / law.load.log_sd,
                           std::log(c.damping_s_m / law.damping.median) / law.damping.log_sd};
      for (int k = 0; k < 3; k++) s1[k] += z[k], s2[k] += z[k] * z[k];
    }
    double worst_mean = 0, worst_sd = 0;
    for (int k = 0; k < 3; k++) {
      const double m = s1[k] / n, sd = std::sqrt(s2[k] / n - m * m);
      worst_mean = std::max(worst_mean, std::fabs(m) * std::sqrt(double(n)));
      worst_sd = std::max(worst_sd, std::fabs(sd - 1) * std::sqrt(2.0 * n));
    }
    s.check("realisation: over 2000 places each parameter's z has mean 0 (in standard errors)", worst_mean, 4);
    s.check("realisation: and deviation 1 (in standard errors)", worst_sd, 4);
  }
  {
    const auto a = on_moon(7.5 * deg, -59 * deg);
    const double step = 1 / ps::MOON_RADIUS_M;  // radians of longitude per metre at the equator, nearly
    const auto near = on_moon(7.5 * deg, -59 * deg + 0.05 * step / std::cos(7.5 * deg));
    const auto far = on_moon(7.5 * deg, -59 * deg + 4.0 * step / std::cos(7.5 * deg));
    // Over 400 seeds: the correlation of z at 5 cm and at 4 m.
    double sa = 0, sn = 0, sf = 0, saa = 0, snn = 0, sff = 0, san = 0, saf = 0;
    const int n = 400;
    for (int seed = 1; seed <= n; seed++) {
      const ps::RegolithGround g(law, std::uint64_t(seed));
      const double za = g.z(a)[0], zn = g.z(near)[0], zf = g.z(far)[0];
      sa += za, sn += zn, sf += zf, saa += za * za, snn += zn * zn, sff += zf * zf, san += za * zn, saf += za * zf;
    }
    const auto corr = [n](double x, double y, double xx, double yy, double xy) {
      return (xy / n - x * y / (n * n)) / std::sqrt((xx / n - x * x / (n * n)) * (yy / n - y * y / (n * n)));
    };
    s.check("realisation: two feet 5 cm apart stand on nearly the same ground (the correlation of their z less e^{−0.00125})",
            corr(sa, sn, saa, snn, san) - std::exp(-0.05 * 0.05 / 2), 0.01);
    s.check("realisation: two places 4 m apart on independent ground (the correlation, e^{−8} in law; four standard errors 0.2)", corr(sa, sf, saa, sff, saf), 0.2);
  }
  {
    const auto x = on_moon(-8.973 * deg, 15.5 * deg);
    const ps::RegolithGround a(law, 11), b(law, 11), c(law, 12);
    s.require("realisation: a seed fixes the ground, another changes it", a.at(x).friction == b.at(x).friction && a.at(x).friction != c.at(x).friction);
  }
  {
    const ps::RegolithField presumed(ground, ps::RegolithField::Quantity::friction, ps::RegolithField::Kind::presumed);
    const ps::RegolithField realised(ground, ps::RegolithField::Quantity::friction, ps::RegolithField::Kind::realised);
    const auto x = on_moon(7.5 * deg, -59 * deg);
    const pf::Vec3 p{x[0], x[1], x[2]};
    const auto e = presumed.at({0, p}), r = realised.at({0, p});
    s.require("fields: the presumed friction is the median, of provenance model",
              e.ok() && e.value == law.friction.median && presumed.descriptor().provenance == pf::Provenance::model && presumed.descriptor().name == "moon.ground.friction.presumed");
    s.require("fields: the realised friction is the ground's, synthetic", r.ok() && r.value == ground->at(x).friction && realised.descriptor().provenance == pf::Provenance::synthetic);
    const auto up = realised.at({0, {p[0] * 1.1, p[1] * 1.1, p[2] * 1.1}});
    s.require("fields: a point 174 km above the ground is refused, with the reason", !up.ok() && up.quality == pf::Quality::outside && !up.reason.empty(), std::string(up.reason));
  }
  return s.done();
}
