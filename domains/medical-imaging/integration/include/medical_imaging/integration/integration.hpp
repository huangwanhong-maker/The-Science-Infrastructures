// Vessels reconstructed from several images of one body together: each acquisition registered to a reference,
// its vessels carried onto the reference grid, the evidence fused, and the vessels told apart by what each
// modality shows (ADR 0070).
//
// One modality rarely suffices. A portal-venous contrast CT shows arteries and veins alike; a TOF-MRA with a
// venous saturation band shows arteries only (modalities.hpp); an MR venogram or a black-blood sequence the
// veins or the walls. Integration rests on three steps, each eikon's or this domain's:
//
// 1. **One frame.** Every acquisition is registered rigidly to the reference by normalised mutual information,
//    which asks only that the two images' values depend on each other (eikon's registration), giving the
//    transform from the reference's world to its own. Patients move between scans; scanners place them
//    differently.
// 2. **One set of vessels.** Each acquisition is reconstructed on its own grid (reconstruction.hpp), with the
//    choices of its modality, and its vessels resampled onto the reference grid through that transform. The
//    vessels are the fusion by STAPLE of the acquisitions that show all vessels (when there are two or more;
//    the one, when one), joined by the arteries that the arterial ones show.
// 3. **Arteries and veins.** The fused vessels are thinned to their graph; a branch is an artery when most of
//    its centreline lies within a voxel of a vessel seen by an acquisition that shows arteries only, a vein when
//    an acquisition that shows veins only sees it or, with an arterial acquisition present, none does. Every
//    vessel voxel takes the class of the nearest centreline voxel (the feature transform).
//
// Classification by branch, not by voxel, keeps a vessel one class along its length, as anatomy has it.
// Checks: tests/check_integration.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include <eikon/fusion/fusion.hpp>
#include <eikon/registration/registration.hpp>
#include <medical_imaging/reconstruction/reconstruction.hpp>

namespace medical_imaging::integration {

using eikon::grid::Geometry;
using eikon::grid::Image;
using Mask = Image<std::uint8_t>;
using eikon::registration::Rigid;

enum class Shows { all, arteries, veins };

struct Acquisition {
  std::string name;
  Image<float> image;
  reconstruction::Options options;
  Shows shows = Shows::all;
  const Mask* region = nullptr;  // on its own grid
};

struct Integrated {
  Geometry reference;
  std::vector<Rigid> to_acquisition;                 // the reference's world → each acquisition's
  std::vector<double> similarity;                    // the normalised mutual information reached
  std::vector<reconstruction::Reconstruction> each;  // on each acquisition's grid (cropped)
  std::vector<Mask> on_reference;                    // each acquisition's vessels on the reference grid
  Mask vessels, arteries, veins;
  eikon::skeleton::Graph graph;                      // of the fused vessels
  std::vector<int> branch_class;                     // per branch: 1 artery, 2 vein, 0 unknown
};

struct Options {
  std::size_t reference = 0;
  bool register_images = true;
  eikon::registration::MutualInformationOptions mutual_information;
};

inline Integrated integrate(const std::vector<Acquisition>& acquisitions, const Options& o = {}) {
  namespace eg = eikon::grid;
  Integrated I;
  const auto& ref = acquisitions.at(o.reference);
  I.reference = ref.image.geometry;
  const std::size_t K = acquisitions.size();
  I.to_acquisition.assign(K, Rigid{});
  I.similarity.assign(K, 0);
  // 1. one frame
  for (std::size_t k = 0; k < K; k++) {
    if (k == o.reference || !o.register_images) continue;
    const auto r = eikon::registration::register_mutual_information(ref.image, acquisitions[k].image, {}, o.mutual_information);
    I.to_acquisition[k] = r.transform, I.similarity[k] = r.value;
  }
  // 2. one set of vessels
  std::vector<Mask> all, arterial, venous;
  for (std::size_t k = 0; k < K; k++) {
    const auto& a = acquisitions[k];
    I.each.push_back(reconstruction::reconstruct(a.image, a.options, a.region));
    const auto full = reconstruction::embed(I.each.back().mask, a.image.geometry);
    const Rigid T = I.to_acquisition[k];
    auto on = k == o.reference ? full : eg::resample(full, I.reference, [&](const eg::Vec3& x) { return T(x); }, true, 0);
    I.on_reference.push_back(on);
    (a.shows == Shows::all ? all : a.shows == Shows::arteries ? arterial : venous).push_back(on);
  }
  if (all.size() >= 2) I.vessels = eikon::fusion::staple(all).fused();
  else if (all.size() == 1) I.vessels = all[0];
  else I.vessels = Mask(I.reference, 0);
  for (const auto* group : {&arterial, &venous})
    for (const auto& m : *group)
      for (std::size_t v = 0; v < m.data.size(); v++) I.vessels.data[v] = I.vessels.data[v] || m.data[v];
  // 3. arteries and veins, by branch
  const auto D = eikon::distance::distance_to(eikon::distance::complement(I.vessels));
  I.graph = eikon::skeleton::graph(eikon::skeleton::prune(eikon::skeleton::thin(I.vessels), D), D);
  const double tol = std::max({I.reference.spacing[0], I.reference.spacing[1], I.reference.spacing[2]});
  const auto near_any = [&](const std::vector<Mask>& ms) {
    Mask u(I.reference, 0);
    for (const auto& m : ms)
      for (std::size_t v = 0; v < m.data.size(); v++) u.data[v] = u.data[v] || m.data[v];
    return eikon::distance::dilate(u, tol);
  };
  const auto art = arterial.empty() ? Mask(I.reference, 0) : near_any(arterial);
  const auto ven = venous.empty() ? Mask(I.reference, 0) : near_any(venous);
  Image<std::int32_t> label(I.reference, 0);
  for (const auto& b : I.graph.branches) {
    std::size_t na = 0, nv = 0;
    for (auto v : b.voxels) na += art.data[v], nv += ven.data[v];
    const double fa = double(na) / double(b.voxels.size()), fv = double(nv) / double(b.voxels.size());
    int c = 0;
    if (!arterial.empty() && fa > 0.5) c = 1;
    else if (!venous.empty() && fv > 0.5) c = 2;
    else if (!arterial.empty()) c = 2;  // seen, but not by the arterial acquisitions
    I.branch_class.push_back(c);
    for (auto v : b.voxels) label.data[v] = std::max(label.data[v], c);
  }
  Mask labelled(I.reference, 0);
  for (std::size_t v = 0; v < label.data.size(); v++) labelled.data[v] = label.data[v] != 0;
  const auto F = eikon::distance::transform(labelled);
  I.arteries = Mask(I.reference, 0), I.veins = Mask(I.reference, 0);
  for (std::size_t v = 0; v < I.vessels.data.size(); v++) {
    if (!I.vessels.data[v] || F.nearest.data[v] < 0) continue;
    const int c = label.data[std::size_t(F.nearest.data[v])];
    if (c == 1) I.arteries.data[v] = 1;
    if (c == 2) I.veins.data[v] = 1;
  }
  return I;
}

}  // namespace medical_imaging::integration
