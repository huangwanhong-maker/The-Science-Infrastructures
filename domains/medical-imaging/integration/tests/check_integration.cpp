// Checks of integration across modalities (integration.hpp; ADR 0070).
//
// One body, an organ with an arterial tree and a venous tree, imaged twice by the modalities' forward models:
// - a portal-venous contrast CT (arteries 250 HU, veins 180 HU, organ 60, fat −100) at 0.7 × 0.7 × 1.0 mm;
// - a TOF-MRA at 0.6 mm, the patient turned by (5°, −3°, 6°) and moved by (2, −1.5, 3) mm, arterial blood
//   entering fresh and venous blood through a saturation band, with Rician noise.
// Neither alone tells arteries from veins in the CT's frame. Integrated, and compared within the organ searched:
// - the MRA registered to the CT within 1 mm over the CT's box;
// - all vessels (Dice at least 0.85 against both trees);
// - arteries and veins told apart (Dice at least 0.8 each against their trees).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <eikon/phantoms/phantoms.hpp>
#include <medical_imaging/integration/integration.hpp>
#include <medical_imaging/modalities/modalities.hpp>

using calx::verification::Suite;
namespace eg = eikon::grid;
namespace ep = eikon::phantoms;
namespace ef = eikon::fusion;
namespace er = eikon::registration;
namespace mm = medical_imaging::modalities;
namespace mr = medical_imaging::reconstruction;
namespace mi = medical_imaging::integration;

int main() {
  Suite s;
  calx::random::Stream rng(20261104, 51, 0, 0);
  constexpr double deg = std::numbers::pi / 180;
  eg::Geometry fine;
  fine.size = {200, 172, 144};
  fine.spacing = {0.35, 0.35, 0.35};
  fine.origin = {-34.825, -29.925, -25.025};
  const double A = 30, B = 25, C = 22;
  eg::Image<float> organ(fine);
  for (std::size_t v = 0; v < fine.count(); v++) {
    const auto x = fine.world(double(fine.unlinear(v)[0]), double(fine.unlinear(v)[1]), double(fine.unlinear(v)[2]));
    organ.data[v] = x[0] * x[0] / (A * A) + x[1] * x[1] / (B * B) + x[2] * x[2] / (C * C) <= 1 ? 1.0f : 0.0f;
  }
  ep::TreeOptions to;
  to.depth = 2;
  to.angle = 0.45;
  to.length_ratio = 0.9;
  const auto arteries = ep::tree({-11, 0, -19}, {-0.1, 0.05, 1}, 13, 2.2, to, rng);   // rising on the left
  const auto veins = ep::tree({11, 2, 19}, {0.1, -0.05, -1}, 13, 2.6, to, rng);       // descending on the right
  const auto art = ep::rasterise(fine, ep::shapes(arteries), 2), ven = ep::rasterise(fine, ep::shapes(veins), 2);

  // the CT, the reference
  mm::CtProtocol cp;
  cp.grid.size = {100, 86, 50};
  cp.grid.spacing = {0.7, 0.7, 1.0};
  cp.grid.origin = {-34.65, -29.75, -24.5};
  cp.background = -100;
  const auto ct = mm::ct({{organ, 60}, {ven.occupancy, 180}, {art.occupancy, 250}}, cp, rng);
  mr::Mask ct_region(cp.grid);
  {
    const auto o = eg::resample(organ, cp.grid);
    for (std::size_t v = 0; v < o.data.size(); v++) ct_region.data[v] = o.data[v] >= 0.5f;
  }
  // the MRA, of the body moved: its world y = T(x) for the CT's world x
  const auto T = er::about({5 * deg, -3 * deg, 6 * deg}, {2, -1.5, 3}, {0, 0, 0});
  const auto Ti = T.inverse();
  mm::TofProtocol tp;
  tp.grid.size = {110, 96, 84};
  tp.grid.spacing = {0.6, 0.6, 0.6};
  tp.grid.origin = {-32.7, -28.5, -24.9};
  tp.slab_entry = {0, 0, -25};
  tp.slab_axis = {0, 0, 1};
  tp.noise = 0.008;
  const auto in_mra = [&](const eg::Image<float>& im) { return eg::resample(im, tp.grid, [&](const eg::Vec3& y) { return Ti(y); }); };
  const auto organ_m = in_mra(organ), art_m = in_mra(art.occupancy), ven_m = in_mra(ven.occupancy);
  eg::Image<float> tissue(tp.grid);
  for (std::size_t v = 0; v < tissue.data.size(); v++) tissue.data[v] = organ_m.data[v] + 0.4f * (1 - organ_m.data[v]);
  const auto mra = mm::tof(tissue, {{art_m, 1.6, 350, false}, {ven_m, 1.6, 120, true}}, tp, rng);
  mr::Mask mra_region(tp.grid);
  for (std::size_t v = 0; v < organ_m.data.size(); v++) mra_region.data[v] = organ_m.data[v] >= 0.5f;

  auto tof_options = mr::tof_mra(0.8, 3.0);
  std::vector<mi::Acquisition> acq{{"contrast CT", ct, mr::contrast_ct(0.8, 3.5), mi::Shows::all, &ct_region},
                                   {"TOF-MRA", mra, tof_options, mi::Shows::arteries, &mra_region}};
  const auto I = mi::integrate(acq);

  double worst = 0;
  for (int c = 0; c < 8; c++) {
    const auto x = cp.grid.world(c & 1 ? 99 : 0, c & 2 ? 85 : 0, c & 4 ? 49 : 0);
    const auto u = I.to_acquisition[1](x), w = T(x);
    worst = std::max(worst, std::sqrt((u[0] - w[0]) * (u[0] - w[0]) + (u[1] - w[1]) * (u[1] - w[1]) + (u[2] - w[2]) * (u[2] - w[2])));
  }
  s.check("the MRA registered to the CT, over the CT's box", worst, 1.0, "mm");
  const auto truth = [&](const ep::Raster& r) {
    const auto o = eg::resample(r.occupancy, cp.grid);
    mr::Mask m(cp.grid);
    for (std::size_t v = 0; v < o.data.size(); v++) m.data[v] = o.data[v] >= 0.5f && ct_region.data[v];  // within the organ searched
    return m;
  };
  const auto ta = truth(art), tv = truth(ven);
  mr::Mask both(cp.grid), found(cp.grid), fa(cp.grid), fv(cp.grid);
  for (std::size_t v = 0; v < both.data.size(); v++) {
    both.data[v] = ta.data[v] || tv.data[v];
    found.data[v] = I.vessels.data[v] && ct_region.data[v], fa.data[v] = I.arteries.data[v] && ct_region.data[v], fv.data[v] = I.veins.data[v] && ct_region.data[v];
  }
  const double dv = ef::dice(found, both), da = ef::dice(fa, ta), dn = ef::dice(fv, tv);
  int na = 0, nv = 0;
  for (int c : I.branch_class) na += c == 1, nv += c == 2;
  std::printf("      (%zu branches: %d arterial, %d venous; NMI %.3f)\n", I.branch_class.size(), na, nv, I.similarity[1]);
  std::printf("      (each alone: the CT against both trees, Dice %.3f; the MRA against the arteries, Dice %.3f)\n", ef::dice(I.on_reference[0], both),
              ef::dice(I.on_reference[1], ta));
  s.require("all vessels: Dice at least 0.85 against both trees", dv >= 0.85, std::to_string(dv));
  s.require("arteries told from veins: the arteries' Dice at least 0.8", da >= 0.8, std::to_string(da));
  s.require("and the veins'", dn >= 0.8, std::to_string(dn));
  return s.done();
}
