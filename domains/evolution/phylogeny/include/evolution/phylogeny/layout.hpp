// Coordinates for drawing a tree: rectangular (a phylogram or cladogram), circular, and the
// equal-angle layout of an unrooted tree.
//
// - Rectangular: a node's x is its distance from the root (or its depth in branches, for a
//   cladogram or a tree without lengths); a leaf's y is its rank from the top, 0 to n − 1, and
//   an internal node's y is halfway between its first and last child's (Felsenstein 2004,
//   Inferring Phylogenies, on drawing trees).
// - Circular: the rectangular layout wrapped round a circle, the angle 2π y / n and the radius
//   x.
// - Equal-angle: each node owns a wedge of the plane proportional to its number of leaves; a
//   child sits at the middle of its wedge, at its branch's length from its parent (the
//   equal-angle algorithm of Meacham's DRAWTREE, PHYLIP; Felsenstein 2004). Branches
//   never cross.
//
// The coordinates are those of the tree's numbering; branches are the segments from a node to
// its parent (rectangular: through the corner (x of the parent, y of the node)). O(n).

#pragma once

#include <evolution/phylogeny/tree.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace evolution::phylogeny {

struct Layout {
  std::vector<double> x, y;            // rectangular: distance from the root, rank of leaves
  std::vector<double> angle, radius;   // circular: radians from the positive x axis, radius
  std::vector<double> ux, uy;          // equal-angle, unrooted, the root at the origin
  std::vector<int> leaf_order;         // the leaves, top to bottom
  std::vector<std::size_t> leaves_below;
  double width = 0;                    // the largest x
};

// With `lengths` false (or a tree lacking them) every branch counts 1.
inline Layout layout(const Tree& t, bool lengths = true) {
  Layout l;
  const std::size_t n = t.size();
  l.x.assign(n, 0);
  l.y.assign(n, 0);
  l.angle.assign(n, 0);
  l.radius.assign(n, 0);
  l.ux.assign(n, 0);
  l.uy.assign(n, 0);
  l.leaves_below.assign(n, 0);
  if (t.root < 0) return l;
  const bool use = lengths && t.has_lengths();
  const auto len = [&](int v) { return use ? std::max(0.0, t[v].length) : 1.0; };
  const auto pre = t.preorder();
  for (int v : pre)
    if (v != t.root) l.x[static_cast<std::size_t>(v)] = l.x[static_cast<std::size_t>(t[v].parent)] + len(v);
  double rank = 0;
  for (int v : pre)
    if (t.is_leaf(v)) {
      l.y[static_cast<std::size_t>(v)] = rank++;
      l.leaf_order.push_back(v);
    }
  for (int v : t.postorder()) {
    const auto& c = t[v].children;
    if (c.empty()) {
      l.leaves_below[static_cast<std::size_t>(v)] = 1;
      continue;
    }
    l.y[static_cast<std::size_t>(v)] = (l.y[static_cast<std::size_t>(c.front())] + l.y[static_cast<std::size_t>(c.back())]) / 2;
    for (int k : c) l.leaves_below[static_cast<std::size_t>(v)] += l.leaves_below[static_cast<std::size_t>(k)];
  }
  const double leaves = std::max(1.0, rank);
  for (std::size_t v = 0; v < n; v++) {
    l.width = std::max(l.width, l.x[v]);
    l.angle[v] = 2 * std::numbers::pi * l.y[v] / leaves;
    l.radius[v] = l.x[v];
  }
  // Equal angle: the wedge [lo, lo + 2π · leaves below / leaves] of each node.
  std::vector<double> lo(n, 0);
  for (int v : pre) {
    double start = lo[static_cast<std::size_t>(v)];
    for (int c : t[v].children) {
      const double w = 2 * std::numbers::pi * static_cast<double>(l.leaves_below[static_cast<std::size_t>(c)]) / leaves;
      const double a = start + w / 2;
      l.ux[static_cast<std::size_t>(c)] = l.ux[static_cast<std::size_t>(v)] + len(c) * std::cos(a);
      l.uy[static_cast<std::size_t>(c)] = l.uy[static_cast<std::size_t>(v)] + len(c) * std::sin(a);
      lo[static_cast<std::size_t>(c)] = start;
      start += w;
    }
  }
  return l;
}

}  // namespace evolution::phylogeny
