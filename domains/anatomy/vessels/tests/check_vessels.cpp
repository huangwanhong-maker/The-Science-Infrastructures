// Checks of vascular trees (vessels.hpp; ADR 0070).
//
// - A full binary tree of depth 4 whose segments are stored in random directions: rooted, every segment points
//   away from the root (each child begins where its parent ends), generations count bifurcations, and the root's
//   Strahler order is 5.
// - Murray's exponent: radii obeying r^γ = r₁^γ + r₂^γ with γ = 3 give 3 to 10⁻⁶; with γ = 2.6, 2.6; with radii
//   perturbed by 5%, within 0.3.
// - A segment's radius ignores the swelling near its junctions; a helix's tortuosity is √(a² + c²)/c.
// - Naming: segments named by an atlas of the same centrelines displaced by 0.5 mm, all correctly; a segment far
//   from every atlas line left unnamed.

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <anatomy/vessels/vessels.hpp>

using calx::verification::Suite;
namespace av = anatomy::vessels;

namespace {

// A full binary tree: node 0 the root, segments with constant radius obeying the law with exponent γ.
av::Tree binary_tree(int depth, double gamma, calx::random::Stream& rng, double noise = 0) {
  av::Tree t;
  t.nodes.push_back({0, 0, 0});
  t.nodes.push_back({0, 0, 20});
  struct Item { std::size_t from, to; double r; int d; av::Vec3 dir; double L; };
  std::vector<Item> items{{0, 1, 3.0, 0, {0, 0, 1}, 20}};
  for (std::size_t q = 0; q < items.size(); q++) {
    const auto it = items[q];
    av::Segment s;
    s.from = it.from, s.to = it.to;
    const auto a = t.nodes[it.from], b = t.nodes[it.to];
    for (int k = 0; k <= 40; k++) {
      const double u = k / 40.0;
      s.points.push_back({a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1]), a[2] + u * (b[2] - a[2])});
      s.radii.push_back(it.r * (1 + noise * rng.normal()));
    }
    s.name = "s" + std::to_string(q);
    if (rng.uniform() < 0.5) {  // stored in a random direction
      std::swap(s.from, s.to);
      std::reverse(s.points.begin(), s.points.end());
      std::reverse(s.radii.begin(), s.radii.end());
    }
    t.segments.push_back(s);
    if (it.d >= depth) continue;
    const double a1 = 0.8, r1 = it.r * std::pow(1 + std::pow(a1, gamma), -1 / gamma), r2 = a1 * r1;
    for (int c = 0; c < 2; c++) {
      const double ang = (c ? -0.5 : 0.5) + 0.3 * it.d, L = it.L * 0.8;
      const av::Vec3 dir{std::sin(ang) * std::cos(1.3 * it.d), std::sin(ang) * std::sin(1.3 * it.d), std::cos(ang)};
      t.nodes.push_back({b[0] + L * dir[0], b[1] + L * dir[1], b[2] + L * dir[2]});
      items.push_back({it.to, t.nodes.size() - 1, c ? r2 : r1, it.d + 1, dir, L});
    }
  }
  return t;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261104, 31, 0, 0);
  {
    auto t = binary_tree(4, 3.0, rng);
    const auto R = av::root_at(t, 0);
    bool away = true;
    for (std::size_t e = 0; e < t.segments.size(); e++)
      if (R.parent[e] >= 0) away = away && t.segments[std::size_t(R.parent[e])].to == t.segments[e].from && av::distance(t.segments[std::size_t(R.parent[e])].points.back(), t.segments[e].points.front()) < 1e-12;
    s.require("rooted: every segment points away from the root, children beginning where parents end", away && t.segments.size() == 31);
    std::size_t root_seg = 0;
    for (std::size_t e = 0; e < t.segments.size(); e++)
      if (R.parent[e] < 0) root_seg = e;
    const int maxgen = *std::max_element(R.generation.begin(), R.generation.end());
    s.require("generations count bifurcations (4 at the leaves); the root's Strahler order is 5", maxgen == 4 && R.strahler[root_seg] == 5,
              std::to_string(maxgen) + ", " + std::to_string(R.strahler[root_seg]));
    const auto bs = av::bifurcations(t, R);
    s.require("15 bifurcations", bs.size() == 15);
    s.check("Murray's exponent of radii obeying γ = 3", av::murray_exponent(bs) - 3, 1e-6);
  }
  {
    auto t = binary_tree(4, 2.6, rng);
    const auto R = av::root_at(t, 0);
    s.check("and of radii obeying γ = 2.6", av::murray_exponent(av::bifurcations(t, R)) - 2.6, 1e-6);
    auto tn = binary_tree(5, 3.0, rng, 0.05);
    const auto Rn = av::root_at(tn, 0);
    s.check("radii perturbed by 5%: γ = 3 within 0.3", av::murray_exponent(av::bifurcations(tn, Rn)) - 3, 0.3);
  }
  {
    av::Segment seg;
    for (int k = 0; k <= 100; k++) {
      seg.points.push_back({0, 0, double(k) * 0.2});
      const double z = k * 0.2;
      seg.radii.push_back(z < 2 || z > 18 ? 4.0 : 2.0);  // swollen within 2 mm of either junction
    }
    s.check("a segment's radius ignores the swelling at its junctions", av::radius(seg) - 2, 1e-12, "mm");
    av::Segment helix;
    const double a = 3, c = 1.5;
    for (int k = 0; k <= 4000; k++) {
      const double th = 4 * std::numbers::pi * k / 4000.0;
      helix.points.push_back({a * std::cos(th), a * std::sin(th), c * th});
    }
    s.check("a helix's tortuosity √(a² + c²)/c, relatively", av::tortuosity(helix) / (std::sqrt(a * a + c * c) / c) - 1, 1e-5);
  }
  {
    auto t = binary_tree(3, 3.0, rng);
    std::vector<std::pair<std::string, std::vector<av::Vec3>>> atlas;
    for (const auto& seg : t.segments) {
      std::vector<av::Vec3> line;
      for (const auto& p : seg.points) line.push_back({p[0] + 0.5, p[1] - 0.2, p[2] + 0.1});
      atlas.push_back({seg.name, line});
    }
    av::Tree named = t;
    av::Segment far;
    far.points = {{100, 100, 100}, {101, 100, 100}};
    far.radii = {1, 1};
    named.segments.push_back(far);
    av::name_by_atlas(named, atlas, 3.0);
    bool all = true;
    for (std::size_t e = 0; e < t.segments.size(); e++) all = all && named.segments[e].name == t.segments[e].name;
    s.require("named by an atlas displaced by 0.5 mm: every segment correctly", all);
    s.require("and a segment far from every atlas line left unnamed", named.segments.back().name.empty());
  }
  return s.done();
}
