// The reconstruction of vessels on real CT: four cases of the Medical Segmentation Decathlon's hepatic vessels
// (Task 08, portal-venous phase, 1.5 mm slices), against their expert labels (reconstruction.hpp; ADR 0070).
//
// The Decathlon gives no liver mask, so the region searched is the annotated territory: the labelled vessels and
// tumours dilated by 10 mm. Within it the contrast-CT preset (radii 1–6 mm) is compared with the vessel label,
// tumour voxels left out of both. Reported per case: Dice, clDice with its centreline precision and sensitivity,
// and the tree's size. The bars are just below the baselines measured when the pipeline was written (2026-11-04:
// mean Dice 0.32 and clDice 0.32; cases 172 and 324 near 0.5, 238 near 0.2, 166 near 0.04, whose parenchyma at
// 82 HU shows weak enhancement), so that a change that makes the pipeline worse on real images fails. The labels
// are incomplete in the smallest branches (the dataset's notes), so precision against them is a lower bound.
// The portal-venous preset's thresholds were tuned on these cases: a held-out evaluation is the open question.
//
// The check reads the data library's product (python data/products/msd_vessels.py after fetching
// msd-hepatic-vessels) and skips without it.

#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <eikon/fusion/fusion.hpp>
#include <eikon/grid/nifti.hpp>
#include <medical_imaging/reconstruction/reconstruction.hpp>

#ifndef MEDICAL_IMAGING_PRODUCTS
#define MEDICAL_IMAGING_PRODUCTS "data/cache/products"
#endif

using calx::verification::Suite;
namespace eg = eikon::grid;
namespace ef = eikon::fusion;
namespace mr = medical_imaging::reconstruction;

int main() {
  Suite s;
  const std::filesystem::path dir = std::filesystem::path(MEDICAL_IMAGING_PRODUCTS) / "msd-hepatic-vessels";
  const std::vector<std::string> cases{"hepaticvessel_166", "hepaticvessel_172", "hepaticvessel_238", "hepaticvessel_324"};
  if (!std::filesystem::exists(dir / (cases[0] + "-image.nii"))) {
    s.skip("the Decathlon's hepatic vessels: the product is not built (python data/acquisition/fetch.py msd-hepatic-vessels; python data/products/msd_vessels.py)");
    return s.done();
  }
  double sum_dice = 0, sum_cl = 0;
  for (const auto& c : cases) {
    const auto image = eg::read_nifti(dir / (c + "-image.nii"));
    const auto label = eg::read_nifti(dir / (c + "-label.nii"));
    mr::Mask vessel(image.geometry), annotated(image.geometry), tumour(image.geometry);
    for (std::size_t v = 0; v < vessel.data.size(); v++) {
      vessel.data[v] = label.data[v] > 0.5f && label.data[v] < 1.5f;
      tumour.data[v] = label.data[v] > 1.5f;
      annotated.data[v] = label.data[v] > 0.5f;
    }
    auto region = eikon::distance::dilate(annotated, 10.0);
    for (std::size_t v = 0; v < region.data.size(); v++) region.data[v] = region.data[v] && !tumour.data[v];
    const auto R = mr::reconstruct(image, mr::portal_venous_ct(1.0, 6.0), &region);
    auto found = mr::embed(R.mask, image.geometry);
    for (std::size_t v = 0; v < found.data.size(); v++) found.data[v] = found.data[v] && region.data[v];
    const double dice = eikon::fusion::dice(found, vessel);
    const auto cl = eikon::fusion::cldice(found, vessel);
    std::printf("      %s: grid %zu × %zu × %zu at %.2f × %.2f × %.2f mm; background %.0f HU, lumen %.0f HU; Dice %.3f, clDice %.3f "
                "(precision %.3f, sensitivity %.3f); %zu branches\n",
                c.c_str(), image.geometry.size[0], image.geometry.size[1], image.geometry.size[2], image.geometry.spacing[0], image.geometry.spacing[1],
                image.geometry.spacing[2], R.background, R.lumen, dice, cl.cldice, cl.precision, cl.sensitivity, R.graph.branches.size());
    sum_dice += dice, sum_cl += cl.cldice;
  }
  const double n = double(cases.size());
  s.require("four real CTs: the mean Dice against the experts' vessels at least the baseline 0.30", sum_dice / n >= 0.30, std::to_string(sum_dice / n));
  s.require("the mean clDice at least the baseline 0.30", sum_cl / n >= 0.30, std::to_string(sum_cl / n));
  return s.done();
}
