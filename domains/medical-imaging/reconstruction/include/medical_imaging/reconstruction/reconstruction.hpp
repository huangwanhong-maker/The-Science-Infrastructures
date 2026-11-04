// Vessels reconstructed from one image: detected by their shape, delineated by their intensity, thinned to
// centrelines, read as a vascular tree with radii, and bounded by a surface (ADR 0070).
//
// The pipeline composes eikon's methods; what is medical is the choice of each step for a modality.
//
// 1. **Region.** The search is confined to a region (an organ's mask, or the body), cropped with a margin, since
//    other bright structures (bone and kidneys in contrast CT) look like vessels too.
// 2. **Detection by shape.** Vesselness across the scales σ = r/√2 for the radii sought (eikon's features), and
//    hysteresis on its response, where the image's value lies in the modality's window for a lumen: the cores.
//    Frangi's c is taken within the region, as half the 99.5th percentile of the Hessian's norm: taken over the
//    whole crop of a real CT, the edges of bone and air set it, and no vessel answered (the real-CT check).
// 3. **Delineation by intensity.** A vesselness response falls off before the wall, so the lumen is grown from
//    the cores through voxels beyond their own half-maximum: halfway between the region's background (its
//    median) and the value at the nearest core (eikon's feature transform), within r_max of it. For a blurred
//    tube the half-maximum is where its wall was; taken locally, it holds where the lumen's brightness varies
//    along the vessel (TOF's inflow fades with depth; a global half-maximum fattened the brightest arteries
//    twofold in the checks). Components smaller than a volume are dropped.
// 4. **Centrelines and tree.** Distance-ordered thinning, spurs pruned, the skeleton's graph with the distance
//    to the outside as radius (eikon's skeleton), and the anatomy domain's vascular tree.
// 5. **Surface.** The mask closed by a voxel (joining parts that touch only at edges or corners), smoothed by
//    half a voxel, its 0.5 level by marching tetrahedra, Taubin's smoothing, and islands smaller than the least
//    volume dropped.
//
// Presets state each modality's choices: contrast CT (bright lumen, Frangi's measure, a window in HU), its
// portal-venous variant for the liver, and TOF-MRA (bright arteries, Jerman's measure, which holds at
// bifurcations). Their thresholds are defaults to be
// tuned on reference labels (Lamy et al. 2022 found no filter's defaults sufficient). Checks: tests/check_reconstruction.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include <anatomy/vessels/vessels.hpp>
#include <eikon/distance/distance.hpp>
#include <eikon/features/vesselness.hpp>
#include <eikon/segmentation/segmentation.hpp>
#include <eikon/skeleton/skeleton.hpp>
#include <eikon/surfaces/surfaces.hpp>

namespace medical_imaging::reconstruction {

using eikon::grid::Geometry;
using eikon::grid::Image;
using Mask = Image<std::uint8_t>;

struct Options {
  double r_min = 1.0, r_max = 6.0;         // mm, the radii sought
  int scales = 6;
  eikon::features::VesselnessOptions vesselness;
  double low = 0.05, high = 0.2;           // hysteresis on the response
  double value_min = -INFINITY;            // a vessel's values lie in [value_min, value_max]: the cores', and
  double value_max = INFINITY;             // those of the voxels its lumen grows through
  double c_quantile = 0.995;               // Frangi's c, when not given: half this quantile of ‖H‖ in the region
  double min_volume = 20;                  // mm³
  double prune = 1.5;
  int taubin = 10;
  bool surface = true;
};

inline Options contrast_ct(double r_min, double r_max) {
  Options o;
  o.r_min = r_min, o.r_max = r_max;
  o.vesselness.measure = eikon::features::Measure::frangi;
  o.value_min = 100;  // HU: an enhanced lumen ...
  o.value_max = 700;  // ... and not cortical bone, which exceeds about 700 HU
  return o;
}

// Portal-venous CT of the liver: the contrast-CT preset with lower hysteresis (0.02, 0.08), tuned on four cases
// of the Medical Segmentation Decathlon's hepatic vessels (the real-CT check), where the parenchyma is itself
// enhanced and the vessels stand only 50–100 HU above it.
inline Options portal_venous_ct(double r_min, double r_max) {
  Options o = contrast_ct(r_min, r_max);
  o.low = 0.02, o.high = 0.08;
  return o;
}

inline Options tof_mra(double r_min, double r_max) {
  Options o;
  o.r_min = r_min, o.r_max = r_max;
  o.vesselness.measure = eikon::features::Measure::jerman;
  o.low = 0.3, o.high = 0.8;
  return o;
}

struct Reconstruction {
  Geometry geometry;                  // the cropped grid the results lie on
  Image<float> response, scale;       // vesselness and its scale
  Mask cores, mask, skeleton;
  Image<float> distance;              // to the outside of the mask: the radius along the skeleton
  eikon::skeleton::Graph graph;
  eikon::surfaces::TriangleSurface surface;
  double background = 0, lumen = 0, threshold = 0;  // the background, the cores' median and their half-maximum (for reports)
};

namespace detail {
inline double median_where(const Image<float>& im, const Mask& m) {
  std::vector<float> v;
  for (std::size_t n = 0; n < im.data.size(); n++)
    if (m.data[n]) v.push_back(im.data[n]);
  if (v.empty()) return 0;
  std::nth_element(v.begin(), v.begin() + std::ptrdiff_t(v.size() / 2), v.end());
  return v[v.size() / 2];
}
}  // namespace detail

// The vessels of an image within a region (a mask on the image's grid; all of it when null).
inline Reconstruction reconstruct(const Image<float>& image, const Options& o, const Mask* region = nullptr) {
  namespace eg = eikon::grid;
  Reconstruction R;
  // 1. the region, cropped with a margin of 2 r_max
  eg::Index3 lo{0, 0, 0}, hi{image.geometry.size[0] - 1, image.geometry.size[1] - 1, image.geometry.size[2] - 1};
  if (region) {
    const double hmin = std::min({image.geometry.spacing[0], image.geometry.spacing[1], image.geometry.spacing[2]});
    if (!eg::bounding_box(*region, [](std::uint8_t v) { return v != 0; }, std::size_t(std::ceil(2 * o.r_max / hmin)), lo, hi)) return R;
  }
  const auto im = eg::crop(image, lo, hi);
  Mask reg(im.geometry, 1);
  if (region) reg = eg::crop(*region, lo, hi);
  R.geometry = im.geometry;
  const bool bright = o.vesselness.polarity == eikon::features::Polarity::bright;
  // 2. detection
  const auto scales = eikon::features::log_scales(o.r_min / std::sqrt(2.0), o.r_max / std::sqrt(2.0), o.scales);
  auto vo = o.vesselness;
  if (vo.measure == eikon::features::Measure::frangi && !(vo.c > 0)) vo.c = eikon::features::frangi_c(im, scales, &reg, o.c_quantile);
  const auto V = eikon::features::vesselness(im, scales, vo);
  R.response = V.response, R.scale = V.scale;
  Image<float> gated(im.geometry);
  for (std::size_t n = 0; n < gated.data.size(); n++)
    gated.data[n] = reg.data[n] && im.data[n] >= o.value_min && im.data[n] <= o.value_max ? V.response.data[n] : 0.0f;
  R.cores = eikon::segmentation::hysteresis(gated, o.low, o.high);
  // 3. delineation at the half-maximum, within r_max of a core. The cores are kept only where the image is
  // beyond the global half-maximum, inside a lumen: a measure that saturates (Jerman's) spills past the wall at
  // coarse scales, and its cores then outgrew the vessels (1.75 times a TOF tree's volume in the checks).
  R.background = detail::median_where(im, reg);
  R.lumen = detail::median_where(im, R.cores);
  R.threshold = (R.background + R.lumen) / 2;
  for (std::size_t n = 0; n < R.cores.data.size(); n++)
    if (R.cores.data[n] && !(bright ? im.data[n] >= R.threshold : im.data[n] <= R.threshold)) R.cores.data[n] = 0;
  // each voxel against its own half-maximum: between the background and the local peak, the most extreme value
  // of the lightly smoothed image among its nearest core and that core's neighbours (the axis lies among them)
  const auto nearest = eikon::distance::transform(R.cores);
  const double r2 = o.r_max * o.r_max;
  const auto S = eikon::filters::smooth(im, 0.5 * std::min({im.geometry.spacing[0], im.geometry.spacing[1], im.geometry.spacing[2]}));
  Image<float> peak(im.geometry, float(R.background));
  std::vector<std::size_t> seeds;
  for (std::size_t n = 0; n < R.cores.data.size(); n++) {
    if (!R.cores.data[n]) continue;
    seeds.push_back(n);
    const auto ijk = im.geometry.unlinear(n);
    float best = S.data[n];
    for (int dk = -1; dk <= 1; dk++)
      for (int dj = -1; dj <= 1; dj++)
        for (int di = -1; di <= 1; di++) {
          const long a = long(ijk[0]) + di, b = long(ijk[1]) + dj, c = long(ijk[2]) + dk;
          if (a < 0 || b < 0 || c < 0 || a >= long(im.geometry.size[0]) || b >= long(im.geometry.size[1]) || c >= long(im.geometry.size[2])) continue;
          const std::size_t u = im.geometry.linear(std::size_t(a), std::size_t(b), std::size_t(c));
          if (R.cores.data[u]) best = bright ? std::max(best, S.data[u]) : std::min(best, S.data[u]);
        }
    peak.data[n] = best;
  }
  Image<float> admit(im.geometry);
  for (std::size_t n = 0; n < admit.data.size(); n++) {
    if (!reg.data[n] || nearest.nearest.data[n] < 0 || nearest.squared.data[n] > r2) continue;
    const double half = (R.background + peak.data[std::size_t(nearest.nearest.data[n])]) / 2;
    admit.data[n] = (bright ? im.data[n] >= half : im.data[n] <= half) && im.data[n] <= o.value_max && im.data[n] >= (bright ? -INFINITY : o.value_min) ? 1.0f : 0.0f;
  }
  for (auto s : seeds) admit.data[s] = 1;  // the cores lie beyond the global half-maximum, and their own peak
  auto grown = eikon::segmentation::grow(admit, seeds, [](float v) { return v > 0; });
  const auto comps = eikon::segmentation::components(grown);
  const std::size_t min_voxels = std::size_t(std::ceil(o.min_volume / im.geometry.voxel_volume()));
  R.mask = eikon::segmentation::keep(comps, 0, std::max<std::size_t>(1, min_voxels));
  // 4. centrelines
  R.distance = eikon::distance::distance_to(eikon::distance::complement(R.mask));
  R.skeleton = eikon::skeleton::prune(eikon::skeleton::thin(R.mask), R.distance, o.prune);
  R.graph = eikon::skeleton::graph(R.skeleton, R.distance);
  // 5. surface
  if (o.surface) {
    // closed by a voxel first: parts the mask joins only by an edge or a corner (26-connectivity) would be
    // separate pieces of a level set
    const double h = std::min({im.geometry.spacing[0], im.geometry.spacing[1], im.geometry.spacing[2]});
    const auto closed = eikon::distance::close(R.mask, std::max({im.geometry.spacing[0], im.geometry.spacing[1], im.geometry.spacing[2]}));
    Image<float> m(im.geometry);
    for (std::size_t n = 0; n < m.data.size(); n++) m.data[n] = closed.data[n];
    R.surface = eikon::surfaces::without_small_pieces(eikon::surfaces::smoothed(eikon::surfaces::isosurface(eikon::filters::smooth(m, 0.5 * h), 0.5), o.taubin), o.min_volume);
  }
  return R;
}

// The skeleton's graph as the anatomy domain's vascular tree (positions and radii in the image's world).
inline anatomy::vessels::Tree tree(const eikon::skeleton::Graph& g) {
  anatomy::vessels::Tree t;
  for (const auto& n : g.nodes) t.nodes.push_back(n.position);
  for (std::size_t e = 0; e < g.branches.size(); e++) {
    anatomy::vessels::Segment s;
    s.from = g.graph.edges[e][0], s.to = g.graph.edges[e][1];
    s.points = g.branches[e].points;
    s.radii = g.branches[e].distance;
    t.segments.push_back(s);
  }
  return t;
}

// A mask on a cropped grid carried to a grid it was cropped from (or any grid sharing its world).
inline Mask embed(const Mask& m, const Geometry& target) { return eikon::grid::resample(m, target, nullptr, true, 0); }

}  // namespace medical_imaging::reconstruction
