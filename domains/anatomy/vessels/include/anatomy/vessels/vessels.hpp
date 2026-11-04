// Vascular trees: vessels as segments between junctions and ends, with their centrelines and radii, rooted
// where blood enters or leaves, ordered, measured, held to the law of branching, and named by an atlas
// (ADR 0070).
//
// A tree is nodes (junctions and ends, in the patient frame, mm) and segments between them, each a centreline
// with a radius at every point: what a reconstruction from images gives (the medical-imaging domain), or what an
// atlas holds. Rooted at a node (the inlet of an arterial tree, the outlet of a venous one), every segment points
// away from the root; its generation counts the bifurcations above it, and its Strahler order rises by one where
// two segments of equal order meet (Strahler 1957, Quantitative analysis of watershed geomorphology, Trans. AGU
// 38, 913), an ordering of rivers that serves vascular trees alike.
//
// Measures of a segment: its length along the centreline, its chord, tortuosity = length/chord, and its radius,
// the median of the radii away from its ends (within one radius of a junction the inscribed radius swells).
// At a bifurcation, the angle between the children's initial directions.
//
// Murray's law (Murray 1926, The physiological principle of minimum work, PNAS 12, 207): minimising the power of
// Poiseuille flow plus a metabolic cost proportional to blood volume gives r_p³ = Σ r_c³ at every branching. The
// exponent is fitted rather than assumed: γ minimising Σ_b (Σ_c (r_c/r_p)^γ − 1)² over the bifurcations, by
// calx's Brent minimiser on [1, 5]; measured exponents of real trees range around 2 to 3, and the survey behind
// ADR 0070 recommends fitting it.
//
// Naming by an atlas: each segment takes the name of the atlas centreline that most of its points lie nearest,
// within a distance, after the atlas has been registered to the patient (eikon's ICP); unnamed otherwise.
// Checks: tests/check_vessels.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include <calx/optimisation/minimise.hpp>

namespace anatomy::vessels {

using Vec3 = std::array<double, 3>;

struct Segment {
  std::size_t from = 0, to = 0;   // nodes
  std::vector<Vec3> points;       // the centreline, from `from` to `to`
  std::vector<double> radii;      // at each point
  std::string name;
};

struct Tree {
  std::vector<Vec3> nodes;
  std::vector<Segment> segments;
};

inline double distance(const Vec3& a, const Vec3& b) { return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])); }

inline double length(const Segment& s) {
  double L = 0;
  for (std::size_t k = 1; k < s.points.size(); k++) L += distance(s.points[k - 1], s.points[k]);
  return L;
}
inline double chord(const Segment& s) { return s.points.empty() ? 0 : distance(s.points.front(), s.points.back()); }
inline double tortuosity(const Segment& s) { const double c = chord(s); return c > 0 ? length(s) / c : 1; }

// The median radius of the points farther than one (local) radius from either end; all points if none is.
inline double radius(const Segment& s) {
  if (s.radii.empty()) return 0;
  std::vector<double> r;
  double along = 0;
  const double L = length(s);
  for (std::size_t k = 0; k < s.points.size(); k++) {
    if (k) along += distance(s.points[k - 1], s.points[k]);
    if (along > s.radii[k] && L - along > s.radii[k]) r.push_back(s.radii[k]);
  }
  if (r.empty()) r = s.radii;
  std::nth_element(r.begin(), r.begin() + std::ptrdiff_t(r.size() / 2), r.end());
  return r[r.size() / 2];
}

struct Rooted {
  std::size_t root = 0;
  std::vector<int> parent;                       // per segment: the segment above it, −1 at the root
  std::vector<std::vector<std::size_t>> children;
  std::vector<int> generation, strahler;
  std::vector<char> reached;                     // per segment: connected to the root
};

// The tree rooted at a node: segments turned to point away from it, with their parents, generations and orders.
inline Rooted root_at(Tree& t, std::size_t root) {
  Rooted R;
  R.root = root;
  const std::size_t S = t.segments.size();
  R.parent.assign(S, -1), R.children.assign(S, {}), R.generation.assign(S, 0), R.strahler.assign(S, 0), R.reached.assign(S, 0);
  std::vector<std::vector<std::size_t>> at(t.nodes.size());
  for (std::size_t e = 0; e < S; e++) at[t.segments[e].from].push_back(e), at[t.segments[e].to].push_back(e);
  std::queue<std::pair<std::size_t, int>> q;  // (node, segment that reached it)
  std::vector<char> seen(t.nodes.size(), 0);
  q.push({root, -1});
  seen[root] = 1;
  std::vector<std::size_t> order;
  while (!q.empty()) {
    const auto [v, via] = q.front();
    q.pop();
    for (auto e : at[v]) {
      if (int(e) == via || R.reached[e]) continue;
      auto& s = t.segments[e];
      if (s.to == v && s.from != v) {  // turn it to point away
        std::swap(s.from, s.to);
        std::reverse(s.points.begin(), s.points.end());
        std::reverse(s.radii.begin(), s.radii.end());
      }
      R.reached[e] = 1;
      R.parent[e] = via;
      if (via >= 0) R.children[std::size_t(via)].push_back(e), R.generation[e] = R.generation[std::size_t(via)] + 1;
      order.push_back(e);
      if (!seen[s.to]) seen[s.to] = 1, q.push({s.to, int(e)});
    }
  }
  for (auto it = order.rbegin(); it != order.rend(); ++it) {  // Strahler from the leaves up
    const auto e = *it;
    if (R.children[e].empty()) { R.strahler[e] = 1; continue; }
    int top = 0, count = 0;
    for (auto c : R.children[e]) {
      if (R.strahler[c] > top) top = R.strahler[c], count = 1;
      else if (R.strahler[c] == top) count++;
    }
    R.strahler[e] = count >= 2 ? top + 1 : top;
  }
  return R;
}

struct Bifurcation {
  std::size_t parent = 0;
  std::vector<std::size_t> children;
  double r_parent = 0;
  std::vector<double> r_children;
  double angle = 0;  // between the first two children's initial directions (rad)
};

inline std::vector<Bifurcation> bifurcations(const Tree& t, const Rooted& R) {
  std::vector<Bifurcation> out;
  const auto direction = [&](const Segment& s) {
    const double r = radius(s);
    std::size_t k = 1;
    double along = 0;
    while (k + 1 < s.points.size() && along < 2 * r) along += distance(s.points[k - 1], s.points[k]), k++;
    Vec3 d{s.points[k][0] - s.points[0][0], s.points[k][1] - s.points[0][1], s.points[k][2] - s.points[0][2]};
    const double n = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    for (auto& x : d) x /= n > 0 ? n : 1;
    return d;
  };
  for (std::size_t e = 0; e < t.segments.size(); e++) {
    if (R.children[e].size() < 2) continue;
    Bifurcation b;
    b.parent = e;
    b.children = R.children[e];
    b.r_parent = radius(t.segments[e]);
    for (auto c : b.children) b.r_children.push_back(radius(t.segments[c]));
    const auto d0 = direction(t.segments[b.children[0]]), d1 = direction(t.segments[b.children[1]]);
    b.angle = std::acos(std::clamp(d0[0] * d1[0] + d0[1] * d1[1] + d0[2] * d1[2], -1.0, 1.0));
    out.push_back(b);
  }
  return out;
}

// The exponent γ of r_p^γ = Σ r_c^γ fitted over the bifurcations.
inline double murray_exponent(const std::vector<Bifurcation>& bs) {
  const auto loss = [&](double g) {
    double s = 0;
    for (const auto& b : bs) {
      if (!(b.r_parent > 0)) continue;
      double sum = 0;
      for (double r : b.r_children) sum += std::pow(r / b.r_parent, g);
      s += (sum - 1) * (sum - 1);
    }
    return s;
  };
  return calx::optimisation::brent_minimise(loss, 1.0, 5.0, 1e-8).x;
}

// Each segment named by the atlas centreline most of its points lie nearest, within `within` (mm).
inline void name_by_atlas(Tree& t, const std::vector<std::pair<std::string, std::vector<Vec3>>>& atlas, double within) {
  for (auto& s : t.segments) {
    std::map<std::string, int> votes;
    for (const auto& p : s.points) {
      double best = within;
      const std::string* name = nullptr;
      for (const auto& [n, line] : atlas)
        for (const auto& q : line)
          if (const double d = distance(p, q); d < best) best = d, name = &n;
      if (name) votes[*name]++;
    }
    s.name.clear();
    int top = 0;
    for (const auto& [n, c] : votes)
      if (2 * c > int(s.points.size()) && c > top) top = c, s.name = n;
  }
}

}  // namespace anatomy::vessels
