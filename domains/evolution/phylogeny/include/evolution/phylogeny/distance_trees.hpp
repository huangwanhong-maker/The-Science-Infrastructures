// Trees from distances: UPGMA, neighbour joining and BIONJ.
//
// - UPGMA (Sokal and Michener 1958, "A statistical method for evaluating systematic
//   relationships", University of Kansas Science Bulletin 38, 1409–1438): join the two
//   closest clusters, at the height of half their distance; the distance of the union to
//   another cluster is the mean of its members' distances (average linkage; WPGMA, the mean of
//   the two clusters' distances, on request). The tree is rooted and ultrametric, and it is the
//   tree of the matrix exactly when the matrix is an ultrametric (the three-point condition).
//   Ties are broken by the first pair in the order of the rows. O(n³) time in the worst case,
//   O(n²) memory.
// - Neighbour joining (Saitou and Nei 1987, Mol. Biol. Evol. 4, 406–425, in the form of Studier
//   and Keppler 1988, Mol. Biol. Evol. 5, 729–731): with r clusters and R_i = Σ_k d(i, k), join
//   the pair minimising Q(i, j) = (r − 2) d(i, j) − R_i − R_j at a new node u, with branches
//   d(i, u) = d(i, j)/2 + (R_i − R_j)/(2(r − 2)) and d(j, u) = d(i, j) − d(i, u), and set
//   d(u, k) = (d(i, k) + d(j, k) − d(i, j))/2. It recovers a tree from its tree metric exactly
//   (Saitou and Nei 1987; Studier and Keppler 1988), and the topology from any
//   matrix whose entries are all within half the shortest branch of the tree metric (Atteson
//   1999, Algorithmica 25, 251–278: the optimal ℓ∞ radius). The tree is unrooted, its root the node of the last three clusters. O(n³) time,
//   O(n²) memory; branch lengths may be negative on a matrix that is not additive, and are
//   left so unless the options clamp them to zero.
// - BIONJ (Gascuel 1997, Mol. Biol. Evol. 14, 685–695): neighbour joining's choice of pair and
//   branch lengths, with the new distances d(u, k) = λ(d(i, k) − d(i, u)) + (1 − λ)(d(j, k) −
//   d(j, u)) at the λ that minimises their variance under a model where variances are
//   proportional to distances: λ = 1/2 + Σ_k (V(j, k) − V(i, k)) / (2(r − 2) V(i, j)), clamped to
//   [0, 1], and V(u, k) = λV(i, k) + (1 − λ)V(j, k) − λ(1 − λ)V(i, j) (Gascuel 1997,
//   with V = d at the start). O(n³).
//
// The matrices are read in the order of their rows; leaves are labelled by the names.

#pragma once

#include <evolution/phylogeny/distance.hpp>
#include <evolution/phylogeny/tree.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace evolution::phylogeny {

enum class Linkage { average, weighted };  // UPGMA, WPGMA

inline Tree upgma(const DistanceMatrix& m, Linkage linkage = Linkage::average) {
  const std::size_t n = m.n;
  if (n == 0) throw std::invalid_argument("evolution::phylogeny::upgma: no taxa");
  Tree t;
  std::vector<double> d = m.d;
  std::vector<int> node(n), size(n, 1);
  std::vector<double> height(n, 0);
  std::vector<bool> active(n, true);
  if (n == 1) {
    t.add(-1, ABSENT, m.names[0]);
    t.rooted = true;
    return t;
  }
  // Build bottom-up: nodes are created as clusters join, the tree assembled at the end.
  struct Join { int a, b; double ha, hb; };
  std::vector<Join> joins;
  std::vector<double> join_height;
  for (std::size_t i = 0; i < n; i++) node[i] = -static_cast<int>(i) - 1;  // a leaf: −(row + 1)
  // The nearest active neighbour of each active row (among the rows after it), for O(n²) scans.
  std::vector<std::size_t> best(n, SIZE_MAX);
  const auto refresh = [&](std::size_t i) {
    best[i] = SIZE_MAX;
    for (std::size_t j = i + 1; j < n; j++)
      if (active[j] && (best[i] == SIZE_MAX || d[i * n + j] < d[i * n + best[i]])) best[i] = j;
  };
  for (std::size_t i = 0; i < n; i++) refresh(i);
  for (std::size_t r = n; r > 1; r--) {
    std::size_t bi = SIZE_MAX;
    for (std::size_t i = 0; i < n; i++)
      if (active[i] && best[i] != SIZE_MAX && (bi == SIZE_MAX || d[i * n + best[i]] < d[bi * n + best[bi]])) bi = i;
    const std::size_t i = bi, j = best[bi];
    const double h = d[i * n + j] / 2;
    joins.push_back({node[i], node[j], height[i], height[j]});
    join_height.push_back(h);
    const int id = static_cast<int>(joins.size()) - 1;
    for (std::size_t k = 0; k < n; k++) {
      if (!active[k] || k == i || k == j) continue;
      const double v = linkage == Linkage::average
                           ? (size[i] * d[i * n + k] + size[j] * d[j * n + k]) / static_cast<double>(size[i] + size[j])
                           : (d[i * n + k] + d[j * n + k]) / 2;
      d[i * n + k] = d[k * n + i] = v;
    }
    active[j] = false;
    size[i] += size[j];
    height[i] = h;
    node[i] = id;
    // Rows whose nearest neighbour was i or j, and i itself, look again; a row before i may now
    // be nearer to the new cluster.
    for (std::size_t k = 0; k < n; k++) {
      if (!active[k]) continue;
      if (k == i || best[k] == i || best[k] == j) refresh(k);
      else if (k < i && d[k * n + i] < d[k * n + best[k]]) best[k] = i;
    }
  }
  // Assemble from the last join down.
  const int root = t.add(-1);
  std::vector<std::pair<int, int>> stack{{static_cast<int>(joins.size()) - 1, root}};
  while (!stack.empty()) {
    const auto [jn, at] = stack.back();
    stack.pop_back();
    const auto& J = joins[static_cast<std::size_t>(jn)];
    for (int side = 0; side < 2; side++) {
      const int c = side == 0 ? J.a : J.b;
      const double hc = side == 0 ? J.ha : J.hb;
      const double len = join_height[static_cast<std::size_t>(jn)] - hc;
      if (c < 0) {
        t.add(at, len, m.names[static_cast<std::size_t>(-c - 1)]);
      } else {
        const int w = t.add(at, len);
        stack.push_back({c, w});
      }
    }
  }
  t.rooted = true;
  return compact(t);
}

struct JoiningOptions {
  bool clamp_negative = false;  // branch lengths below zero set to zero
};

namespace detail {

// Neighbour joining and BIONJ share everything but the update of the distances.
inline Tree join_neighbours(const DistanceMatrix& m, bool bionj, const JoiningOptions& o) {
  const std::size_t n = m.n;
  if (n < 2) throw std::invalid_argument("evolution::phylogeny::neighbour_joining: at least two taxa");
  Tree t;
  const int root = t.add(-1);
  if (n == 2) {
    const double h = m(0, 1) / 2;
    t.add(root, h, m.names[0]);
    t.add(root, h, m.names[1]);
    t.rooted = false;
    return t;
  }
  std::vector<double> d = m.d, V;
  if (bionj) V = m.d;
  // The node of each slot; children are attached when their parent is made, so the nodes are
  // created as free subtrees and joined under new nodes.
  std::vector<int> node(n);
  std::vector<Node> pool;  // assembled into the tree at the end
  const auto make = [&](std::string label) {
    Node x;
    x.label = std::move(label);
    pool.push_back(std::move(x));
    return static_cast<int>(pool.size()) - 1;
  };
  for (std::size_t i = 0; i < n; i++) node[i] = make(m.names[i]);
  std::vector<std::size_t> act(n);
  for (std::size_t i = 0; i < n; i++) act[i] = i;
  std::vector<double> R(n, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < n; k++) R[i] += d[i * n + k];
  const auto clamp = [&](double x) { return o.clamp_negative ? std::max(0.0, x) : x; };
  while (act.size() > 3) {
    const double r = static_cast<double>(act.size());
    std::size_t bi = 0, bj = 1;
    double bq = std::numeric_limits<double>::infinity();
    for (std::size_t a = 0; a < act.size(); a++) {
      const std::size_t i = act[a];
      const double* di = &d[i * n];
      for (std::size_t b = a + 1; b < act.size(); b++) {
        const std::size_t j = act[b];
        const double q = (r - 2) * di[j] - R[i] - R[j];
        if (q < bq) { bq = q; bi = a; bj = b; }
      }
    }
    const std::size_t i = act[bi], j = act[bj];
    const double dij = d[i * n + j];
    const double li = dij / 2 + (R[i] - R[j]) / (2 * (r - 2)), lj = dij - li;
    double lambda = 0.5;
    if (bionj) {
      if (V[i * n + j] > 0) {
        double s = 0;
        for (std::size_t k : act)
          if (k != i && k != j) s += V[j * n + k] - V[i * n + k];
        lambda = std::clamp(0.5 + s / (2 * (r - 2) * V[i * n + j]), 0.0, 1.0);
      }
    }
    const int u = make({});
    pool[static_cast<std::size_t>(node[i])].parent = u;
    pool[static_cast<std::size_t>(node[i])].length = clamp(li);
    pool[static_cast<std::size_t>(node[j])].parent = u;
    pool[static_cast<std::size_t>(node[j])].length = clamp(lj);
    pool[static_cast<std::size_t>(u)].children = {node[i], node[j]};
    // The new cluster takes slot i.
    for (std::size_t k : act) {
      if (k == i || k == j) continue;
      const double dik = d[i * n + k], djk = d[j * n + k];
      const double duk = bionj ? lambda * (dik - li) + (1 - lambda) * (djk - lj) : (dik + djk - dij) / 2;
      R[k] += duk - dik - djk;
      d[i * n + k] = d[k * n + i] = duk;
      if (bionj) {
        const double v = lambda * V[i * n + k] + (1 - lambda) * V[j * n + k] - lambda * (1 - lambda) * V[i * n + j];
        V[i * n + k] = V[k * n + i] = v;
      }
    }
    R[i] = 0;
    for (std::size_t k : act)
      if (k != j && k != i) R[i] += d[i * n + k];
    node[i] = u;
    act.erase(act.begin() + static_cast<std::ptrdiff_t>(bj));
  }
  // The last three clusters meet at the root.
  const std::size_t a = act[0], b = act[1], c = act[2];
  const double la = (d[a * n + b] + d[a * n + c] - d[b * n + c]) / 2;
  const double lb = (d[a * n + b] + d[b * n + c] - d[a * n + c]) / 2;
  const double lc = (d[a * n + c] + d[b * n + c] - d[a * n + b]) / 2;
  const int top = make({});
  for (auto [k, l] : {std::pair{a, la}, std::pair{b, lb}, std::pair{c, lc}}) {
    pool[static_cast<std::size_t>(node[k])].parent = top;
    pool[static_cast<std::size_t>(node[k])].length = clamp(l);
    pool[static_cast<std::size_t>(top)].children.push_back(node[k]);
  }
  // Copy the pool into a tree hung from `top`.
  std::vector<std::pair<int, int>> stack{{top, root}};
  while (!stack.empty()) {
    const auto [src, dst] = stack.back();
    stack.pop_back();
    for (int c2 : pool[static_cast<std::size_t>(src)].children) {
      const auto& x = pool[static_cast<std::size_t>(c2)];
      const int w = t.add(dst, x.length, x.label);
      stack.push_back({c2, w});
    }
  }
  t.rooted = false;
  return compact(t);
}

}  // namespace detail

inline Tree neighbour_joining(const DistanceMatrix& m, const JoiningOptions& o = {}) { return detail::join_neighbours(m, false, o); }
inline Tree bionj(const DistanceMatrix& m, const JoiningOptions& o = {}) { return detail::join_neighbours(m, true, o); }

}  // namespace evolution::phylogeny
