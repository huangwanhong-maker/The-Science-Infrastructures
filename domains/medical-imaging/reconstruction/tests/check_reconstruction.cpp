// Checks of the reconstruction of vessels from one image (reconstruction.hpp; ADR 0070).
//
// A contrast CT of a known body, made by the modalities' forward model: an organ (60 HU) in fat (−100 HU) with a
// tree of veins (200 HU) of depth 3 obeying Murray's law, imaged at 0.7 × 0.7 × 1.0 mm with a PSF of 0.8/1.2 mm
// and 15 HU of noise. Reconstructed within the organ, against the truth on the CT's grid:
// - overlap (Dice) at least 0.85 and centreline overlap (clDice) at least 0.9; one component, no loop;
// - the tree's 7 junctions and 9 ends, within one;
// - each segment's radius within half a voxel (0.35 mm) of its tube's, on average;
// - Murray's exponent fitted on the reconstructed tree, within 0.75 of the truth's 3 (radii of 1–3 mm at 0.7 mm
//   voxels: the exponent is sensitive to them, and the check states how sensitive);
// - the surface a closed manifold enclosing the truth's volume within 10%.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <eikon/fusion/fusion.hpp>
#include <eikon/phantoms/phantoms.hpp>
#include <medical_imaging/modalities/modalities.hpp>
#include <medical_imaging/reconstruction/reconstruction.hpp>

using calx::verification::Suite;
namespace eg = eikon::grid;
namespace ep = eikon::phantoms;
namespace ef = eikon::fusion;
namespace mm = medical_imaging::modalities;
namespace mr = medical_imaging::reconstruction;
namespace av = anatomy::vessels;

int main() {
  Suite s;
  calx::random::Stream rng(20261104, 41, 0, 0);
  // the truth on a fine grid
  eg::Geometry fine;
  fine.size = {200, 172, 144};
  fine.spacing = {0.35, 0.35, 0.35};
  fine.origin = {-34.825, -29.925, -25.025};
  const double A = 30, B = 25, C = 22;  // the organ's semi-axes
  eg::Image<float> organ(fine);
  for (std::size_t v = 0; v < fine.count(); v++) {
    const auto ijk = fine.unlinear(v);
    const auto x = fine.world(double(ijk[0]), double(ijk[1]), double(ijk[2]));
    organ.data[v] = (x[0] * x[0]) / (A * A) + (x[1] * x[1]) / (B * B) + (x[2] * x[2]) / (C * C) <= 1 ? 1.0f : 0.0f;
  }
  ep::TreeOptions to;
  to.depth = 3;
  to.exponent = 3;
  to.length_ratio = 0.85;
  to.angle = 0.6;
  const auto tree = ep::tree({-2, 1, -18}, {0.05, -0.05, 1}, 13, 2.6, to, rng);
  const auto vessels = ep::rasterise(fine, ep::shapes(tree), 2);
  // the CT
  mm::CtProtocol p;
  p.grid.size = {100, 86, 50};
  p.grid.spacing = {0.7, 0.7, 1.0};
  p.grid.origin = {-34.65, -29.75, -24.5};
  p.fwhm_inplane = 0.8, p.fwhm_z = 1.2, p.noise = 15, p.background = -100;
  const auto ct = mm::ct({{organ, 60}, {vessels.occupancy, 200}}, p, rng);
  const auto organ_ct = eg::resample(organ, p.grid);
  mr::Mask region(p.grid);
  for (std::size_t v = 0; v < region.data.size(); v++) region.data[v] = organ_ct.data[v] >= 0.5f;
  const auto truth_occ = eg::resample(vessels.occupancy, p.grid);
  mr::Mask truth(p.grid);
  for (std::size_t v = 0; v < truth.data.size(); v++) truth.data[v] = truth_occ.data[v] >= 0.5f;

  const auto R = mr::reconstruct(ct, mr::contrast_ct(0.8, 3.5), &region);
  const auto mask = mr::embed(R.mask, p.grid);
  const double dice = ef::dice(mask, truth);
  const auto cl = ef::cldice(mask, truth);
  const auto betti = ef::betti(mask);
  std::printf("      (background %.1f HU, lumen %.1f HU, threshold %.1f HU)\n", R.background, R.lumen, R.threshold);
  s.require("contrast CT of a known tree: Dice at least 0.85", dice >= 0.85, std::to_string(dice));
  s.require("clDice at least 0.9", cl.cldice >= 0.9, std::to_string(cl.cldice));
  s.require("one component, no loop", betti == ef::Betti{1, 0, 0}, std::to_string(betti.b0) + ", " + std::to_string(betti.b1));
  const long J = long(R.graph.count(eikon::skeleton::NodeKind::junction)), E = long(R.graph.count(eikon::skeleton::NodeKind::end));
  s.require("7 junctions and 9 ends, within one", std::abs(J - 7) <= 1 && std::abs(E - 9) <= 1, std::to_string(J) + " junctions, " + std::to_string(E) + " ends");

  auto t = mr::tree(R.graph);
  // root at the end nearest the tree's inlet
  std::size_t root = 0;
  for (std::size_t n = 0; n < t.nodes.size(); n++)
    if (av::distance(t.nodes[n], tree.tubes[0].a) < av::distance(t.nodes[root], tree.tubes[0].a)) root = n;
  const auto rooted = av::root_at(t, root);
  // each segment against the tube its points lie nearest
  double err = 0;
  int matched = 0;
  for (const auto& seg : t.segments) {
    std::map<std::size_t, int> votes;
    for (const auto& q : seg.points) {
      std::size_t best = 0;
      for (std::size_t i = 1; i < tree.tubes.size(); i++)
        if (std::fabs(ep::distance(ep::Tube{tree.tubes[i].a, tree.tubes[i].b, 0, 0, 1}, q)) < std::fabs(ep::distance(ep::Tube{tree.tubes[best].a, tree.tubes[best].b, 0, 0, 1}, q))) best = i;
      votes[best]++;
    }
    const auto top = std::max_element(votes.begin(), votes.end(), [](auto& a, auto& b) { return a.second < b.second; })->first;
    if (av::length(seg) < 2) continue;
    err += std::fabs(av::radius(seg) - tree.tubes[top].ra), matched++;
  }
  s.check("each segment's radius against its tube's, on average", matched ? err / matched : INFINITY, 0.35, "mm");
  const double gamma = av::murray_exponent(av::bifurcations(t, rooted));
  s.check("Murray's exponent of the reconstructed tree against the truth's 3", gamma - 3, 0.75);
  long chi = 0;
  {
    std::map<std::pair<std::size_t, std::size_t>, int> directed;
    for (const auto& tr : R.surface.triangles)
      for (int a = 0; a < 3; a++) directed[{tr[std::size_t(a)], tr[std::size_t((a + 1) % 3)]}]++;
    bool closed = true;
    for (const auto& [e, c] : directed) closed = closed && c == 1 && directed.count({e.second, e.first});
    chi = long(R.surface.points.size()) - long(directed.size() / 2) + long(R.surface.triangles.size());
    double vt = 0;
    for (float x : vessels.occupancy.data) vt += x;
    vt *= fine.voxel_volume();
    const double vs = eikon::surfaces::enclosed_volume(R.surface);
    // its pieces, joined through their triangles
    std::vector<std::size_t> up(R.surface.points.size());
    for (std::size_t i = 0; i < up.size(); i++) up[i] = i;
    const auto find = [&](std::size_t x) { while (up[x] != x) x = up[x] = up[up[x]]; return x; };
    for (const auto& tr : R.surface.triangles) up[find(tr[0])] = find(tr[1]), up[find(tr[1])] = find(tr[2]);
    std::map<std::size_t, std::size_t> pieces;
    for (const auto& tr : R.surface.triangles) pieces[find(tr[0])]++;
    std::string sizes;
    for (const auto& [r, c] : pieces) sizes += std::to_string(c) + " ";
    s.require("the surface a closed manifold of genus 0", closed && chi == 2, "χ " + std::to_string(chi) + (closed ? ", closed" : ", open") + ", pieces of " + sizes + "triangles");
    s.check("enclosing the truth's volume, relatively", vs / vt - 1, 0.1);
  }
  return s.done();
}
