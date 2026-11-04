// Checks of the sections of the body (sections.hpp).
//
// On boxes placed in the LPS frame:
// - a transverse section through a box is one piece, its loop counter-clockwise in the image (positive
//   area), of the box's cross-section;
// - a box on the subject's left lies on the image's right in the radiological transverse display and on its
//   left in the neurological one; a posterior box lies on the right of the sagittal display;
// - a box the plane does not reach is not cut.
// - The level at a height: the vertebra whose interval holds it, nothing above or below the column.
// On BodyParts3D, when the product has been built: the vertebrae C1 to the sacrum descend in order; the
// transverse sections of the left kidney integrate to its volume (Cavalieri's principle) by the midpoint rule
// on 400 planes; and a transverse section at the tenth thoracic vertebra cuts the aorta, the inferior vena
// cava and the oesophagus, with the aorta on the subject's left of the vena cava, as at that level it lies.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include <anatomy/body/body.hpp>
#include <anatomy/measures/measures.hpp>
#include <anatomy/sections/sections.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace as = anatomy::sections;
namespace af = anatomy::frames;
namespace cg = physis::computational_geometry;

namespace {
cg::TriangleSurface box(cg::Point3 lo, cg::Point3 hi) {
  cg::TriangleSurface S;
  for (int k = 0; k < 8; k++) S.points.push_back({k & 1 ? hi[0] : lo[0], k & 2 ? hi[1] : lo[1], k & 4 ? hi[2] : lo[2]});
  const int quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
  for (const auto& q : quads) {
    S.triangles.push_back({std::size_t(q[0]), std::size_t(q[1]), std::size_t(q[2])});
    S.triangles.push_back({std::size_t(q[0]), std::size_t(q[2]), std::size_t(q[3])});
  }
  return S;
}
double mean_u(const as::Piece& p) {
  double s = 0;
  std::size_t n = 0;
  for (const auto& l : p.loops)
    for (const auto& q : l) s += q[0], n++;
  return n ? s / n : 0;
}
}  // namespace

int main() {
  Suite s;
  {
    const auto left = box({40, -10, 900}, {80, 30, 1000});      // on the subject's left, astride z = 950
    const auto back = box({-10, 60, 900}, {10, 90, 1000});      // posterior, on the midline
    const auto high = box({-10, -10, 1200}, {10, 10, 1300});    // above the plane
    std::vector<const cg::TriangleSurface*> S{&left, &back, &high};
    std::vector<cg::Box3> bounds{cg::bounds(left.points), cg::bounds(back.points), cg::bounds(high.points)};
    const auto R = as::section(S, bounds, af::view(af::Plane::transverse, af::Convention::radiological), 950);
    s.require("a transverse section at z = 950 mm cuts the two boxes that reach it, not the one above", R.size() == 2 && R[0].part == 0 && R[1].part == 1);
    s.check("the left box's piece: one loop, counter-clockwise, of its cross-section 40 × 40 mm", R[0].area - 1600, 1e-9, "mm²");
    s.require("in the radiological display the left box lies on the image's right", mean_u(R[0]) > 0);
    const auto N = as::section(S, bounds, af::view(af::Plane::transverse, af::Convention::neurological), 950);
    s.require("in the neurological display, on its left", mean_u(N[0]) < 0 && std::fabs(N[0].area - 1600) < 1e-9);
    const auto G = as::section(S, bounds, af::view(af::Plane::sagittal), 0);
    s.require("a sagittal section through the midline cuts the two midline boxes, the posterior one on the image's right",
              G.size() == 2 && G[0].part == 1 && G[1].part == 2 && mean_u(G[0]) > 0 && G[0].area > 0 && std::fabs(mean_u(G[1])) < 1e-9);
  }
  {
    const std::vector<as::Level> L{{"T9", 1250}, {"T10", 1225}, {"T11", 1200}};
    s.require("a height is named by the vertebra whose interval holds it",
              as::level_at(L, 1226) == "T10" && as::level_at(L, 1238) == "T9" && as::level_at(L, 1212) == "T11" && as::level_at(L, 1300).empty() && as::level_at(L, 1150).empty());
  }
#ifdef ANATOMY_PRODUCTS
  {
    const std::filesystem::path dir = ANATOMY_PRODUCTS;
    if (!std::filesystem::exists(dir / "body-bp3d-elements.tsv")) {
      s.skip("sections of BodyParts3D: the product is not built (python data/products/body_anatomy.py)");
    } else {
      const auto B = anatomy::body::read_body(dir);
      std::map<std::string, std::size_t> by_name;
      for (std::size_t e = 0; e < B.elements.size(); e++)
        if (B.elements[e].duplicate_of == anatomy::parts::none) by_name.emplace(B.elements[e].name, e);
      std::vector<as::Level> levels;
      for (const auto& v : as::VERTEBRAE) {
        const auto c = B.atlas.find(v.fma);
        for (std::size_t e = 0; e < B.elements.size(); e++)
          if (B.elements[e].concept_ == c && B.elements[e].duplicate_of == anatomy::parts::none) {
            levels.push_back({std::string(v.label), anatomy::measures::measure(B.surface(e)).centre[2]});
            break;
          }
      }
      bool descending = levels.size() == as::VERTEBRAE.size();
      for (std::size_t k = 1; k < levels.size(); k++) descending = descending && levels[k].height < levels[k - 1].height;
      s.require("the 24 vertebrae and the sacrum, C1 to S, descend in order", descending, std::to_string(levels.size()) + " found");

      const auto it = by_name.find("Left kidney");
      if (it != by_name.end()) {
        const auto K = B.surface(it->second);
        const auto M = anatomy::measures::measure(K);
        const auto W = cg::welded(K, anatomy::measures::weld_tolerance);
        const auto V = af::view(af::Plane::transverse);
        const int n = 400;
        const double lo = M.bounds.lower[2], hi = M.bounds.upper[2], dz = (hi - lo) / n;
        double sum = 0;
        for (int i = 0; i < n; i++) {
          const auto P = as::section({&W}, {M.bounds}, V, lo + (i + 0.5) * dz);
          for (const auto& p : P) sum += p.area * dz;
        }
        s.check("Cavalieri on the left kidney: its transverse sections integrate to its volume, relative", sum / M.volume - 1, 2e-4);
      }

      double t10 = 0;
      for (const auto& l : levels)
        if (l.label == "T10") t10 = l.height;
      // Every element of each name (the inferior vena cava is drawn in two pieces).
      std::vector<cg::TriangleSurface> keep;
      std::vector<std::string> names;
      for (std::size_t e = 0; e < B.elements.size(); e++) {
        const auto& n = B.elements[e].name;
        if (B.elements[e].duplicate_of == anatomy::parts::none && (n == "Descending thoracic aorta" || n == "Inferior vena cava" || n == "Esophagus"))
          keep.push_back(B.surface(e)), names.push_back(n);
      }
      std::vector<const cg::TriangleSurface*> S;
      std::vector<cg::Box3> bounds;
      for (const auto& k : keep) S.push_back(&k), bounds.push_back(cg::bounds(k.points));
      const auto P = as::section(S, bounds, af::view(af::Plane::transverse), t10);
      std::map<std::string, double> u;
      for (const auto& p : P) u[names[p.part]] = mean_u(p);
      s.require("at T10 the section cuts the aorta, the inferior vena cava and the oesophagus", u.size() == 3,
                "T10 at " + std::to_string(t10) + " mm; " + std::to_string(u.size()) + " of 3 cut");
      if (u.size() == 3) s.require("the aorta lies to the subject's left of the vena cava (the image's right)", u["Descending thoracic aorta"] > u["Inferior vena cava"]);
    }
  }
#endif
  return s.done();
}
