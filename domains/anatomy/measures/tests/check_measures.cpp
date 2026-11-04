// Checks of the measures of a part (measures.hpp).
//
// On closed forms:
// - a unit cube: a sphere (one closed component, genus 0, consistent), area 6, volume 1, centre (½, ½, ½);
//   the same cube oriented inward is read reversed, with the same volume;
// - a box 30 × 20 × 10 mm turned by a rotation: its extents along its principal axes are 30, 20 and 10, the
//   long axis first, along the turned x axis;
// - a torus of radii 40 and 12 mm: genus 1, its volume 2π²Rr² within the polygon's error;
// - an open square: a disc, not closed, no volume;
// - two cubes combined: volume 2, the centre between them; and the two touch within 0.6 mm, not 0.4.
// On BodyParts3D, when the product has been built (data/products/body_anatomy.py):
// - the volumes of the parts whose meshes bound a solid, against the volumes the database published in the
//   files' headers (tests/fixtures/bodyparts3d-volumes.tsv, made by tests/make_bp3d_fixtures.py from the
//   files themselves). The median ratio is 1 to 0.5 %, and the parts over 100 cm³ agree to 1 %; smaller and
//   thinner parts scatter more, to factors of 40 for some small vessels, as if the database had computed its
//   volumes on its full-resolution meshes, of which the distributed files keep 1 % of the polygons
//   (theory/anatomy.md §4.3). The volumes here are those of the meshes as distributed, exact for them.
// - sides: of the parts named left or right whose centres lie more than 100 mm from the midline (the limbs,
//   where a name's side is the body's), all lie on their named side in LPS but the database's one
//   mislabelled pair, the flexor pollicis brevis (FJ1469 named left on the right, FJ1469M the reverse).
//   Nearer the midline a name's side may be a chamber's or a branch's: the right ventricle and the right
//   coronary artery's posterior branches lie left of the midline.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <anatomy/body/body.hpp>
#include <anatomy/measures/measures.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace am = anatomy::measures;
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

cg::TriangleSurface torus(double R, double r, int m, int n) {
  cg::TriangleSurface S;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      const double u = 2 * std::numbers::pi * i / m, v = 2 * std::numbers::pi * j / n;
      S.points.push_back({(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)});
    }
  auto at = [&](int i, int j) { return std::size_t(((i + m) % m) * n + (j + n) % n); };
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      S.triangles.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1)});
      S.triangles.push_back({at(i, j), at(i + 1, j + 1), at(i, j + 1)});
    }
  return S;
}

}  // namespace

int main() {
  Suite s;
  {
    const auto M = am::measure(box({0, 0, 0}, {1, 1, 1}));
    s.require("a unit cube: one closed component of genus 0, consistently oriented",
              M.topology.components == 1 && M.topology.closed == 1 && M.topology.genus == 0 && M.topology.consistent && M.topology.euler == 2);
    s.check("its area is 6", M.area - 6, 1e-15);
    s.check("its volume is 1", M.volume - 1, 1e-15);
    s.check("its centre is (½, ½, ½)", std::hypot(M.centre[0] - 0.5, M.centre[1] - 0.5, M.centre[2] - 0.5), 1e-15);
    auto inward = box({0, 0, 0}, {1, 1, 1});
    for (auto& t : inward.triangles) std::swap(t[1], t[2]);
    const auto N = am::measure(inward);
    s.require("the cube oriented inward is read reversed", N.inward && !M.inward);
    s.check("with the same volume", N.volume - 1, 1e-15);
  }
  {
    auto B = box({-15, -10, -5}, {15, 10, 5});
    const double c = std::cos(0.6), sn = std::sin(0.6);
    for (auto& p : B.points) p = {c * p[0] - sn * p[1], sn * p[0] + c * p[1], p[2]};
    const auto M = am::measure(B);
    s.check("a box 30 × 20 × 10 mm turned: its length along its long axis", M.extents[0] - 30, 1e-12, "mm");
    s.check("its breadth", M.extents[1] - 20, 1e-12, "mm");
    s.check("its thickness", M.extents[2] - 10, 1e-12, "mm");
    s.check("the long axis is the turned x axis", 1 - std::fabs(M.axes[0][0] * c + M.axes[0][1] * sn), 1e-12);
  }
  {
    const double R = 40, r = 12;
    const auto M = am::measure(torus(R, r, 96, 48));
    s.require("a torus: one closed component of genus 1", M.topology.components == 1 && M.topology.closed == 1 && M.topology.genus == 1 && M.topology.consistent);
    s.check("its volume 2π²Rr² within the polygon's error, relative", M.volume / (2 * std::numbers::pi * std::numbers::pi * R * r * r) - 1, 5e-3);
  }
  {
    cg::TriangleSurface Q;
    Q.points = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    Q.triangles = {{0, 1, 2}, {0, 2, 3}};
    const auto M = am::measure(Q);
    s.require("an open square: a surface with one boundary circle, not closed, no volume",
              M.topology.surfaces == 1 && M.topology.closed == 0 && M.topology.boundary_circles == 1 && !M.solid());
  }
  {
    const auto A = box({0, 0, 0}, {1, 1, 1}), B = box({1.5, 0, 0}, {2.5, 1, 1});
    const auto MA = am::measure(A), MB = am::measure(B);
    const auto C = am::combine({&MA, &MB}, {&A, &B});
    s.check("two cubes combined: volume 2", C.volume - 2, 1e-15);
    s.check("the centre between them", std::hypot(C.centre[0] - 1.25, C.centre[1] - 0.5, C.centre[2] - 0.5), 1e-15);
    s.check("their extent along the long axis, 2.5", C.extents[0] - 2.5, 1e-12);
    const auto near = am::touching({&A, &B}, 0.6), far = am::touching({&A, &B}, 0.4);
    s.require("they touch within 0.6 mm, 0.5 apart, and not within 0.4", near.size() == 1 && std::fabs(near[0].distance - 0.5) < 1e-15 && far.empty());
  }
#ifdef ANATOMY_PRODUCTS
  {
    const std::filesystem::path dir = ANATOMY_PRODUCTS;
    const std::filesystem::path fixture = std::filesystem::path(ANATOMY_FIXTURES) / "bodyparts3d-volumes.tsv";
    if (!std::filesystem::exists(dir / "body-bp3d-elements.tsv")) {
      s.skip("the body's measures against BodyParts3D: the product is not built (python data/products/body_anatomy.py)");
    } else {
      const auto Bd = anatomy::body::read_body(dir);
      std::map<std::string, double> published;
      std::ifstream in(fixture);
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line.rfind("id\t", 0) == 0) continue;
        std::istringstream row(line);
        std::string id, volume;
        std::getline(row, id, '\t');
        std::getline(row, volume, '\t');
        published[id] = std::stod(volume);
      }
      std::size_t solids = 0, sided = 0;
      std::vector<double> ratios, large;
      std::vector<std::string> wrong;
      for (std::size_t e = 0; e < Bd.elements.size(); e++) {
        const auto& E = Bd.elements[e];
        const auto M = am::measure(Bd.surface(e));
        const auto side = anatomy::parts::side_by_name(E.name);
        if (side != anatomy::parts::Side::none && std::fabs(M.centre[0]) > 100) {
          sided++;
          if ((side == anatomy::parts::Side::left) != (M.centre[0] > 0)) wrong.push_back(E.id);
        }
        if (!M.solid()) continue;
        solids++;
        const auto it = published.find(E.id);
        if (it == published.end() || !(it->second > 0)) continue;
        const double r = M.volume / 1000 / it->second;
        ratios.push_back(r);
        if (it->second > 100) large.push_back(r);
      }
      s.require("the fixture names the database's volume of every element", published.size() == Bd.elements.size());
      std::sort(ratios.begin(), ratios.end());
      const double median = ratios.empty() ? 0 : ratios[ratios.size() / 2];
      s.check("the median ratio of the volumes to the database's", median - 1, 5e-3);
      const auto agree = std::count_if(large.begin(), large.end(), [](double r) { return std::fabs(r - 1) <= 0.01; });
      s.check("of the parts over 100 cm³, the share not within 1 % of the database's volume", 1 - double(agree) / double(large.size()), 0.1);
      std::sort(wrong.begin(), wrong.end());
      s.require("parts named left or right more than 100 mm from the midline lie on their named side (LPS: left x > 0), "
                "but the database's mislabelled flexor pollicis brevis",
                sided > 200 && wrong == std::vector<std::string>{"FJ1469", "FJ1469M"}, std::to_string(wrong.size()) + " of " + std::to_string(sided));
      std::printf("      %zu of %zu elements bound a solid; %zu volumes compared, %zu over 100 cm3; %zu sided parts far from the midline\n",
                  solids, Bd.elements.size(), ratios.size(), large.size(), sided);
    }
  }
#endif
  return s.done();
}
