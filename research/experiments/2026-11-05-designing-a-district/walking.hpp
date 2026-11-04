// Walking distances through a district's network: the study's own provisional code, not a foundation's.
//
// The urban-science domain holds the walking network as a regular cell complex (urban::landscape::District: junctions,
// segments with their polylines, blocks) and places a holding on the segment nearest to it (attach_store: the segment
// and the arc length along it). What it does not hold, nor does any foundation, is the least walking distance between
// two places on such a network: morphe's combinatorics gives distances in edges (breadth first), and its semirings
// (MinPlus) give the algebra but no algorithm on a graph with cycles. The study writes Dijkstra's algorithm here,
// marked as its own, and reports the gap (README §8): a single-source shortest distance on a weighted multigraph belongs
// beside morphe::combinatorics::distances.
//
// Dijkstra (1959), "A note on two problems in connexion with graphs", Numer. Math. 1, 269–271: with non-negative
// lengths, settling nodes in order of their tentative distance gives the least distance to each. Sources and targets
// lie on segments (at arc length s from the segment's tail): a source on segment (u, v) of length L starts u at s and v
// at L − s; a target on (u′, v′) at s′ is reached at min(d(u′) + s′, d(v′) + L′ − s′), or directly along its own
// segment from a source on it, |s − s′|.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace study {

struct WalkGraph {
  std::size_t nodes = 0;
  struct Edge { std::size_t u = 0, v = 0; double length = 0; };
  std::vector<Edge> edges;
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> incident;  // node ↦ (edge, other end)
  void finish() {
    incident.assign(nodes, {});
    for (std::size_t e = 0; e < edges.size(); e++) {
      incident[edges[e].u].push_back({e, edges[e].v});
      incident[edges[e].v].push_back({e, edges[e].u});
    }
  }
};

// A place on the network: a segment and the arc length along it from the segment's tail.
struct OnSegment { std::size_t edge = 0; double s = 0; };

// The least distance from any of the sources to each node.
inline std::vector<double> from_sources(const WalkGraph& G, const std::vector<OnSegment>& sources) {
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> d(G.nodes, inf);
  using Item = std::pair<double, std::size_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> q;
  const auto relax = [&](std::size_t x, double v) { if (v < d[x]) { d[x] = v; q.push({v, x}); } };
  for (const auto& s : sources) {
    const auto& e = G.edges[s.edge];
    relax(e.u, s.s), relax(e.v, std::max(0.0, e.length - s.s));
  }
  while (!q.empty()) {
    const auto [dx, x] = q.top();
    q.pop();
    if (dx > d[x]) continue;
    for (const auto& [e, y] : G.incident[x]) relax(y, dx + G.edges[e].length);
  }
  return d;
}

// The least distance from the sources to a target, given the nodes' distances from them.
inline double to_target(const WalkGraph& G, const std::vector<double>& d, const std::vector<OnSegment>& sources, const OnSegment& t) {
  const auto& e = G.edges[t.edge];
  double best = std::min(d[e.u] + t.s, d[e.v] + std::max(0.0, e.length - t.s));
  for (const auto& s : sources)
    if (s.edge == t.edge) best = std::min(best, std::fabs(s.s - t.s));
  return best;
}

}  // namespace study
