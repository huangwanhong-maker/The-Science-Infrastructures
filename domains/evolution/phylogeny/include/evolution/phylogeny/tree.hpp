// Phylogenetic trees: rooted and unrooted, with lengths and supports on their branches; their
// traversals, and rerooting at a node, on an edge, and at the midpoint.
//
// A tree is stored hanging from a node, its root: each node has a parent, its children, and
// the label it carries. A branch is named by the node below it, which holds its length and its
// support (the branch's data travel with the branch, not with a node, when a tree is rerooted).
// An unrooted tree is the same structure with `rooted` false: its root is any internal node,
// conventionally one of degree three, and has no meaning. A rooted tree's root is a node of
// degree two, the point on a branch where the tree's history begins. Between the two, the
// operations are those of Felsenstein (2004), Inferring Phylogenies:
//
//   reroot(t, v, x)   the root on the branch above v, at distance x from v: the undirected tree
//                     hung from a new point, every node of degree two left by the old root
//                     suppressed (its two branches joined, lengths added);
//   root_at(t, v)     the tree hung from node v, unrooted;
//   unroot(t)         a rooted tree's root suppressed, its two branches joined;
//   midpoint_root(t)  the root at the midpoint of the longest path between two leaves (Farris
//                     1972, "Estimating phylogenetic trees from distance matrices", American
//                     Naturalist 106, 645–668), found by two searches for the farthest leaf, which
//                     is exact on a tree with lengths that are not negative.
//
// Node indices are ints, the root's parent is −1, and every operation that changes the shape
// returns a new tree whose nodes are numbered in preorder. Absent lengths and supports are NaN
// (has_length, has_support). All operations are O(n) in the number of nodes, the leaf
// distance matrix O(n · leaves).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace evolution::phylogeny {

inline constexpr double ABSENT = std::numeric_limits<double>::quiet_NaN();

struct Node {
  int parent = -1;
  std::vector<int> children;
  double length = ABSENT;   // of the branch to the parent
  double support = ABSENT;  // of that branch: the support of the split it makes
  std::string label;
};

struct Tree {
  std::vector<Node> nodes;
  int root = -1;
  bool rooted = false;

  std::size_t size() const { return nodes.size(); }
  bool is_leaf(int v) const { return nodes[static_cast<std::size_t>(v)].children.empty(); }
  const Node& operator[](int v) const { return nodes[static_cast<std::size_t>(v)]; }
  Node& operator[](int v) { return nodes[static_cast<std::size_t>(v)]; }
  bool has_length(int v) const { return !std::isnan((*this)[v].length); }
  bool has_support(int v) const { return !std::isnan((*this)[v].support); }
  double length_or(int v, double fallback) const { return has_length(v) ? (*this)[v].length : fallback; }

  // Adds a node below `parent` (the root when parent is −1 and the tree is empty).
  int add(int parent, double length = ABSENT, std::string label = {}) {
    const int v = static_cast<int>(nodes.size());
    Node n;
    n.parent = parent;
    n.length = length;
    n.label = std::move(label);
    nodes.push_back(std::move(n));
    if (parent >= 0) nodes[static_cast<std::size_t>(parent)].children.push_back(v);
    else if (root < 0) root = v;
    else throw std::invalid_argument("evolution::phylogeny::Tree::add: a tree has one root");
    return v;
  }

  // Nodes from the root, each before its children, children in order.
  std::vector<int> preorder() const {
    std::vector<int> out;
    if (root < 0) return out;
    out.reserve(nodes.size());
    std::vector<int> stack{root};
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      out.push_back(v);
      const auto& c = (*this)[v].children;
      for (std::size_t k = c.size(); k-- > 0;) stack.push_back(c[k]);
    }
    return out;
  }
  // Each node after its children: the preorder of the mirrored tree, reversed.
  std::vector<int> postorder() const {
    std::vector<int> out;
    if (root < 0) return out;
    out.reserve(nodes.size());
    std::vector<int> stack{root};
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      out.push_back(v);
      for (int c : (*this)[v].children) stack.push_back(c);
    }
    std::reverse(out.begin(), out.end());
    return out;
  }
  // The leaves, left to right.
  std::vector<int> leaves() const {
    std::vector<int> out;
    for (int v : preorder())
      if (is_leaf(v)) out.push_back(v);
    return out;
  }
  std::size_t leaf_count() const {
    std::size_t n = 0;
    for (const auto& x : nodes) n += x.children.empty();
    return n;
  }
  // The leaf with a label, or −1.
  int find_leaf(std::string_view label) const {
    for (std::size_t v = 0; v < nodes.size(); v++)
      if (nodes[v].children.empty() && nodes[v].label == label) return static_cast<int>(v);
    return -1;
  }
  // The labels of the leaves, left to right.
  std::vector<std::string> leaf_labels() const {
    std::vector<std::string> out;
    for (int v : leaves()) out.push_back((*this)[v].label);
    return out;
  }
  // Whether every branch has a length.
  bool has_lengths() const {
    for (std::size_t v = 0; v < nodes.size(); v++)
      if (static_cast<int>(v) != root && std::isnan(nodes[v].length)) return false;
    return true;
  }
  // The sum of the branch lengths (absent lengths count 0).
  double total_length() const {
    double s = 0;
    for (std::size_t v = 0; v < nodes.size(); v++)
      if (static_cast<int>(v) != root && !std::isnan(nodes[v].length)) s += nodes[v].length;
    return s;
  }
  // Whether every internal node has two children (rooted) or the root three and the others two.
  bool is_binary() const {
    for (std::size_t v = 0; v < nodes.size(); v++) {
      const std::size_t k = nodes[v].children.size();
      if (k == 0) continue;
      if (static_cast<int>(v) == root ? (k != 2 && k != 3) : k != 2) return false;
    }
    return true;
  }
};

// The distance of every node from the root (absent lengths count `missing`).
inline std::vector<double> depths(const Tree& t, double missing = 0) {
  std::vector<double> d(t.size(), 0);
  for (int v : t.preorder())
    if (v != t.root) d[static_cast<std::size_t>(v)] = d[static_cast<std::size_t>(t[v].parent)] + t.length_or(v, missing);
  return d;
}

namespace detail {

// The tree as an undirected graph: for each node its neighbours, and for each neighbour the
// branch between them, named by the node below it in the original tree.
struct Undirected {
  std::vector<std::vector<int>> next;
  std::vector<std::vector<int>> branch;  // the child end of the edge
};
inline Undirected undirected(const Tree& t) {
  Undirected g;
  g.next.resize(t.size());
  g.branch.resize(t.size());
  for (std::size_t v = 0; v < t.size(); v++) {
    const int p = t.nodes[v].parent;
    if (p < 0) continue;
    g.next[v].push_back(p);
    g.branch[v].push_back(static_cast<int>(v));
    g.next[static_cast<std::size_t>(p)].push_back(static_cast<int>(v));
    g.branch[static_cast<std::size_t>(p)].push_back(static_cast<int>(v));
  }
  return g;
}

inline double add_lengths(double a, double b) {
  if (std::isnan(a)) return b;
  if (std::isnan(b)) return a;
  return a + b;
}

// Hangs the undirected tree from `start` (a node), its first branches leading to the nodes
// in `first` with the data given, suppressing nodes left with degree two (other than the
// start). `from` is the neighbour of start not to be entered (−1 for none).
struct Hang {
  const Tree& t;
  const Undirected& g;
  Tree out;

  // Copies the subtree entered at node v from node `came`, below parent `under` in out.
  void copy(int v, int came, int under, double length, double support) {
    // A node of degree two (one way on) is suppressed: its two branches are one.
    for (;;) {
      const auto& nx = g.next[static_cast<std::size_t>(v)];
      if (nx.size() != 2) break;
      const int k = nx[0] == came ? 1 : 0;
      const int b = g.branch[static_cast<std::size_t>(v)][static_cast<std::size_t>(k)];
      length = add_lengths(length, t[b].length);
      if (std::isnan(support)) support = t[b].support;
      came = v;
      v = nx[static_cast<std::size_t>(k)];
    }
    const int w = out.add(under, length, t[v].label);
    out[w].support = support;
    const auto& nx = g.next[static_cast<std::size_t>(v)];
    for (std::size_t k = 0; k < nx.size(); k++) {
      if (nx[k] == came) continue;
      const int b = g.branch[static_cast<std::size_t>(v)][k];
      copy(nx[k], v, w, t[b].length, t[b].support);
    }
  }
};

}  // namespace detail

// The tree with its root on the branch above v, at distance x from v (0 ≤ x ≤ the branch's
// length; with an absent length the two new branches have none). The result is rooted, with a
// root of degree two.
inline Tree reroot(const Tree& t, int v, double x) {
  if (v < 0 || static_cast<std::size_t>(v) >= t.size() || v == t.root)
    throw std::invalid_argument("evolution::phylogeny::reroot: the root goes on a branch, named by the node below it");
  const double len = t[v].length;
  if (!std::isnan(len) && (x < 0 || x > len)) throw std::invalid_argument("evolution::phylogeny::reroot: 0 <= x <= the branch's length");
  const auto g = detail::undirected(t);
  detail::Hang h{t, g, {}};
  const int r = h.out.add(-1);
  const double sup = t[v].support;
  h.copy(v, t[v].parent, r, std::isnan(len) ? ABSENT : x, sup);
  h.copy(t[v].parent, v, r, std::isnan(len) ? ABSENT : len - x, sup);
  h.out.rooted = true;
  return h.out;
}

// The tree hung from node v (unrooted; a former root of degree two is suppressed).
inline Tree root_at(const Tree& t, int v) {
  if (v < 0 || static_cast<std::size_t>(v) >= t.size()) throw std::invalid_argument("evolution::phylogeny::root_at: no such node");
  const auto g = detail::undirected(t);
  detail::Hang h{t, g, {}};
  const int r = h.out.add(-1, ABSENT, t[v].label);
  const auto& nx = g.next[static_cast<std::size_t>(v)];
  for (std::size_t k = 0; k < nx.size(); k++) {
    const int b = g.branch[static_cast<std::size_t>(v)][k];
    h.copy(nx[k], v, r, t[b].length, t[b].support);
  }
  h.out.rooted = false;
  return h.out;
}

// The tree without its root: a root of degree two is suppressed, the tree hung from an internal
// neighbour. A tree of two leaves keeps its root.
inline Tree unroot(const Tree& t) {
  if (t.root < 0) return t;
  const auto& c = t[t.root].children;
  if (c.size() == 2) {
    for (int k : c)
      if (!t.is_leaf(k)) return root_at(t, k);
  }
  Tree out = root_at(t, t.root);
  out.rooted = false;
  return out;
}

// The tree in preorder numbering (unchanged in shape).
inline Tree compact(const Tree& t) {
  Tree out;
  if (t.root < 0) return out;
  std::vector<int> map(t.size(), -1);
  for (int v : t.preorder()) {
    const int p = t[v].parent;
    const int w = out.add(p < 0 ? -1 : map[static_cast<std::size_t>(p)], t[v].length, t[v].label);
    out[w].support = t[v].support;
    map[static_cast<std::size_t>(v)] = w;
  }
  out.rooted = t.rooted;
  return out;
}

// The distances from node s to every node (absent lengths count `missing`): O(n).
inline std::vector<double> distances_from(const Tree& t, int s, double missing = 0) {
  const auto g = detail::undirected(t);
  std::vector<double> d(t.size(), -1);
  std::vector<int> stack{s};
  d[static_cast<std::size_t>(s)] = 0;
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    const auto& nx = g.next[static_cast<std::size_t>(v)];
    for (std::size_t k = 0; k < nx.size(); k++) {
      const int w = nx[k];
      if (d[static_cast<std::size_t>(w)] >= 0) continue;
      d[static_cast<std::size_t>(w)] = d[static_cast<std::size_t>(v)] + t.length_or(g.branch[static_cast<std::size_t>(v)][k], missing);
      stack.push_back(w);
    }
  }
  return d;
}

// The patristic distances between the leaves labelled by `taxa`, as a flat n × n matrix.
inline std::vector<double> leaf_distances(const Tree& t, const std::vector<std::string>& taxa) {
  const std::size_t n = taxa.size();
  std::vector<int> leaf(n);
  for (std::size_t i = 0; i < n; i++) {
    leaf[i] = t.find_leaf(taxa[i]);
    if (leaf[i] < 0) throw std::invalid_argument("evolution::phylogeny::leaf_distances: no leaf " + taxa[i]);
  }
  std::vector<double> out(n * n, 0);
  for (std::size_t i = 0; i < n; i++) {
    const auto d = distances_from(t, leaf[i]);
    for (std::size_t j = 0; j < n; j++) out[i * n + j] = d[static_cast<std::size_t>(leaf[j])];
  }
  return out;
}

// The root at the midpoint of the longest leaf-to-leaf path.
inline Tree midpoint_root(const Tree& t) {
  const auto ls = t.leaves();
  if (ls.size() < 2) return t;
  const auto farthest = [&](const std::vector<double>& d) {
    int best = ls[0];
    for (int v : ls)
      if (d[static_cast<std::size_t>(v)] > d[static_cast<std::size_t>(best)]) best = v;
    return best;
  };
  const int a = farthest(distances_from(t, ls[0]));
  const auto da = distances_from(t, a);
  const int b = farthest(da);
  const double half = da[static_cast<std::size_t>(b)] / 2;
  // Walk from b towards a until the midpoint is passed: the path is b's ancestors up to the
  // meeting point, then down to a; with distances from a, the edge (u, w) straddling the
  // midpoint has da[u] ≥ half ≥ da[w].
  const auto g = detail::undirected(t);
  std::vector<int> prev(t.size(), -2), via(t.size(), -1);
  std::vector<int> stack{a};
  prev[static_cast<std::size_t>(a)] = -1;
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    const auto& nx = g.next[static_cast<std::size_t>(v)];
    for (std::size_t k = 0; k < nx.size(); k++) {
      const int w = nx[k];
      if (prev[static_cast<std::size_t>(w)] != -2) continue;
      prev[static_cast<std::size_t>(w)] = v;
      via[static_cast<std::size_t>(w)] = g.branch[static_cast<std::size_t>(v)][k];
      stack.push_back(w);
    }
  }
  int u = b;
  while (prev[static_cast<std::size_t>(u)] >= 0 && da[static_cast<std::size_t>(prev[static_cast<std::size_t>(u)])] >= half) u = prev[static_cast<std::size_t>(u)];
  // Now da[u] ≥ half > da[prev[u]] (or u = a when the path has no length).
  const int p = prev[static_cast<std::size_t>(u)];
  if (p < 0) return reroot(t, ls[0] == a ? ls[1] : ls[0], 0);
  const int br = via[static_cast<std::size_t>(u)];  // the branch between p and u, named by its lower node
  const double from_u = da[static_cast<std::size_t>(u)] - half;
  // reroot measures from the branch's lower node.
  const double x = br == u ? from_u : t.length_or(br, 0) - from_u;
  return reroot(t, br, std::clamp(x, 0.0, t.length_or(br, 0)));
}

// Children sorted by their number of leaves, the larger last (ascending) or first: a drawing
// convention, the shape unchanged.
inline Tree ladderize(const Tree& t, bool larger_last = true) {
  std::vector<std::size_t> count(t.size(), 0);
  for (int v : t.postorder()) {
    if (t.is_leaf(v)) count[static_cast<std::size_t>(v)] = 1;
    for (int c : t[v].children) count[static_cast<std::size_t>(v)] += count[static_cast<std::size_t>(c)];
  }
  Tree s = t;
  for (auto& n : s.nodes)
    std::stable_sort(n.children.begin(), n.children.end(), [&](int x, int y) {
      return larger_last ? count[static_cast<std::size_t>(x)] < count[static_cast<std::size_t>(y)] : count[static_cast<std::size_t>(x)] > count[static_cast<std::size_t>(y)];
    });
  return compact(s);
}

}  // namespace evolution::phylogeny
