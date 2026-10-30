// Checks of the ocean's water as a body meets it (water/presumed.hpp, ADR 0059).
//
// - The presumption: the density's 95 % interval is 1020–1030 kg/m³, the range of surface seawater.
// - The realisation: over 1000 places 50 km apart, the density's z and each component of the current have
//   mean 0 and the presumed deviation (within four standard errors); the current is horizontal; a seed fixes
//   it.
// - As fields: presumed (model) and realised (synthetic); a point a kilometre up is refused with the reason.

#include <oceanography/water/presumed.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>

namespace ow = oceanography::water;
namespace pf = physis::fields;

namespace {

ow::UpperOcean::Point at_sea(double lat, double lon, double h = 0) {
  const double a = 6378137.0, f = 1 / 298.257223563, e2 = f * (2 - f), s = std::sin(lat), n = a / std::sqrt(1 - e2 * s * s);
  return {(n + h) * std::cos(lat) * std::cos(lon), (n + h) * std::cos(lat) * std::sin(lon), (n * (1 - e2) + h) * s};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const ow::WaterLaw law;
  s.check("presumption: the density's 95 % interval begins at 1020 kg/m³ (kg/m³)", law.density_median * std::exp(-1.959964 * law.density_log_sd) - 1020, 0.05);
  s.check("presumption: and ends at 1030 kg/m³ (kg/m³)", law.density_median * std::exp(1.959964 * law.density_log_sd) - 1030, 0.05);
  const auto ocean = std::make_shared<const ow::UpperOcean>(law, 3);
  {
    const int n = 1000;
    double s1[3] = {}, s2[3] = {}, radial = 0;
    const double d = std::numbers::pi / 180;
    for (int i = 0; i < n; i++) {
      // Along the equator of the Pacific, 50 km apart: independent draws.
      const auto x = at_sea(0.0, (-170.0 + i * 50e3 / 6378137.0 / d) * d);
      const auto w = ocean->at(x);
      const double z[3] = {std::log(w.density / law.density_median) / law.density_log_sd, w.east / law.current_sd, w.north / law.current_sd};
      for (int k = 0; k < 3; k++) s1[k] += z[k], s2[k] += z[k] * z[k];
      const auto c = ocean->current_fixed(x, w);
      const double r = std::hypot(x[0], x[1], x[2]);
      radial = std::max(radial, std::fabs((c[0] * x[0] + c[1] * x[1] + c[2] * x[2]) / r));
    }
    double wm = 0, ws = 0;
    for (int k = 0; k < 3; k++) {
      const double m = s1[k] / n, sd = std::sqrt(s2[k] / n - m * m);
      wm = std::max(wm, std::fabs(m) * std::sqrt(double(n))), ws = std::max(ws, std::fabs(sd - 1) * std::sqrt(2.0 * n));
    }
    s.check("realisation: over 1000 places 50 km apart the density's z and the current's components have mean 0 (in standard errors)", wm, 4);
    s.check("realisation: and their presumed deviations (in standard errors)", ws, 4);
    s.check("realisation: the current is horizontal (its component along the vertical, m/s)", radial, 1e-15);
  }
  {
    const auto x = at_sea(-10 * std::numbers::pi / 180, -120 * std::numbers::pi / 180);
    const ow::UpperOcean a(law, 9), b(law, 9), c(law, 10);
    s.require("realisation: a seed fixes the water, another changes it", a.at(x).density == b.at(x).density && a.at(x).east != c.at(x).east);
  }
  {
    const ow::DensityField presumed(ocean, true), realised(ocean, false);
    const ow::CurrentField current(ocean, false);
    const auto x = at_sea(0.2, -2.9, -1.0);
    const pf::Vec3 p{x[0], x[1], x[2]};
    const auto e = presumed.at({0, p}), r = realised.at({0, p});
    const auto c = current.at({0, p});
    s.require("fields: the presumed density is the median, of provenance model", e.ok() && e.value == 1025 && presumed.descriptor().provenance == pf::Provenance::model);
    s.require("fields: the realised density and current are the water's, synthetic",
              r.ok() && r.value == ocean->at(x).density && c.ok() && realised.descriptor().provenance == pf::Provenance::synthetic);
    const auto up = realised.at({0, {p[0] * (1 + 1000 / 6.37e6), p[1] * (1 + 1000 / 6.37e6), p[2] * (1 + 1000 / 6.37e6)}});
    s.require("fields: a point a kilometre up is refused, with the reason", !up.ok() && up.quality == pf::Quality::outside, std::string(up.reason));
  }
  return s.done();
}
