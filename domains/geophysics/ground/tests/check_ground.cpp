// Checks of the Earth's ground as a foot meets it (ground/presumed.hpp, ADR 0059), on CRUST1.0 from the data
// library's cache (skipped, and said, when it is not there).
//
// - The class from CRUST1.0: ice on the Greenland ice sheet and in Antarctica, the sea in the equatorial
//   Pacific, sediment in the Paris Basin.
// - The presumption: each class's medians; the sea has none.
// - The realisation: over 1000 places of one class, each parameter's z of mean 0 and deviation 1 (within four
//   standard errors); a seed fixes it.
// - As fields: presumed (model) and realised (synthetic), refused at sea and high above the ground with the
//   reason.

#include <geophysics/ground/presumed.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <numbers>
#include <string>

namespace gg = geophysics::ground;
namespace gs = geophysics::structure;
namespace pf = physis::fields;

namespace {

std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(GEOPHYSICS_ROOT "/../../data/cache")) + "/" + relative;
}

gg::EarthGround::Point on_earth(double lat_deg, double lon_deg, double h = 0) {
  const double a = 6378137.0, f = 1 / 298.257223563, e2 = f * (2 - f), la = lat_deg * std::numbers::pi / 180, lo = lon_deg * std::numbers::pi / 180;
  const double n = a / std::sqrt(1 - e2 * std::sin(la) * std::sin(la));
  return {(n + h) * std::cos(la) * std::cos(lo), (n + h) * std::cos(la) * std::sin(lo), (n * (1 - e2) + h) * std::sin(la)};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const std::string product = cache_path("products/crust-crust1.bin");
  if (!std::filesystem::exists(product)) {
    s.skip("ground: crust-crust1 is not in the cache (python data/acquisition/fetch.py crust1; python data/products/crust_model.py)");
    return s.done();
  }
  const auto crust = std::make_shared<const gs::Crust1>(gs::Crust1::read_product(product));
  const auto ground = std::make_shared<const gg::EarthGround>(crust, 7);
  s.require("class: the Greenland ice sheet is ice (72.6° N, 38.5° W, Summit)", ground->class_at(72.6, -38.5) == gg::GroundClass::ice);
  s.require("class: Antarctica's plateau is ice (90° S)", ground->class_at(-89.5, 0.0) == gg::GroundClass::ice);
  s.require("class: the equatorial Pacific is the sea (0°, 150° W)", ground->class_at(0.0, -150.0) == gg::GroundClass::sea);
  s.require("class: the Paris Basin is sediment (48.86° N, 2.35° E)", ground->class_at(48.86, 2.35) == gg::GroundClass::sediment);
  {
    const auto ice = ground->presumed(72.6, -38.5), soil = ground->presumed(48.86, 2.35);
    s.require("presumption: ice's friction is far below sediment's (its median 0.15 against 0.6)", ice.friction == 0.15 && soil.friction == 0.6);
    s.require("presumption: the sea has no ground's law", ground->presumed(0.0, -150.0).friction == 0);
  }
  {
    // 1000 places on the Greenland ice sheet, along 72.6° N a few tens of metres apart: independent draws.
    const int n = 1000;
    double s1 = 0, s2 = 0;
    for (int i = 0; i < n; i++) {
      const auto x = on_earth(72.6, -38.5 + i * 50.0 / (6.378e6 * std::cos(72.6 * std::numbers::pi / 180)) * 180 / std::numbers::pi);
      const auto g = ground->at(x);
      const double z = std::log(g.friction / 0.15) / 0.45;
      s1 += z, s2 += z * z;
    }
    const double m = s1 / n, sd = std::sqrt(s2 / n - m * m);
    s.check("realisation: over 1000 places of ice its friction's z has mean 0 (in standard errors)", m * std::sqrt(double(n)), 4);
    s.check("realisation: and deviation 1 (in standard errors)", (sd - 1) * std::sqrt(2.0 * n), 4);
    const gg::EarthGround a(crust, 3), b(crust, 3), c(crust, 4);
    const auto x = on_earth(48.86, 2.35);
    s.require("realisation: a seed fixes the ground, another changes it", a.at(x).friction == b.at(x).friction && a.at(x).friction != c.at(x).friction);
  }
  {
    const gg::GroundField presumed(ground, gg::GroundField::Quantity::friction, true), realised(ground, gg::GroundField::Quantity::friction, false);
    const auto x = on_earth(48.86, 2.35, 35.0);
    const pf::Vec3 p{x[0], x[1], x[2]};
    const auto e = presumed.at({0, p}), r = realised.at({0, p});
    s.require("fields: the presumed friction in Paris is sediment's median, of provenance model",
              e.ok() && e.value == 0.6 && presumed.descriptor().provenance == pf::Provenance::model && presumed.descriptor().name == "earth.ground.friction.presumed");
    s.require("fields: the realised one is the ground's, synthetic", r.ok() && r.value == ground->at(x).friction && realised.descriptor().provenance == pf::Provenance::synthetic);
    const auto sea = on_earth(0.0, -150.0);
    const auto w = realised.at({0, {sea[0], sea[1], sea[2]}});
    s.require("fields: at sea it is refused, the sea named", !w.ok() && w.quality == pf::Quality::outside && std::string(w.reason).find("sea") != std::string::npos, std::string(w.reason));
    const auto up = on_earth(48.86, 2.35, 3.0e4);
    s.require("fields: 30 km up it is refused", !realised.at({0, {up[0], up[1], up[2]}}).ok());
  }
  return s.done();
}
