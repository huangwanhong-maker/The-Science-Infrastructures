// The refined complex K̂ of a history (theory note §2.4): the carrier on which signed ties live.
//
//   vertices: the presences (b, o);
//   edges: first the temporal edges, a being between consecutive occasions (in the order of the
//     history complex's edges), then the rungs, one for each tie at each meeting within its life,
//     from (a, o) to (b, o);
//   2-cells: first the lenses, squares a[p, q] · r(q) · b[p, q]⁻¹ · r(p)⁻¹ (in the order of the
//     history complex's lenses), then the triangles of rungs at an occasion, one for each three
//     ties joining three beings pairwise there.
//
// Collapsing the rungs of each occasion gives back the history complex K. A frustrated triangle
// of K̂ is curvature at an occasion that K, with its rungs collapsed, cannot carry: that is why
// the signs of ties are read on K̂ (theory note §6).

#pragma once

#include <history/construction/complex.hpp>
#include <morphe/topology/cw_complex.hpp>

#include <array>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace history {

struct RefinedTriangle {
  std::size_t occasion = 0;
  std::array<std::size_t, 3> beings{};  // x, y, z
  std::array<std::size_t, 3> ties{};    // x–y, y–z, x–z
  std::array<std::size_t, 3> rungs{};   // the index of each tie's rung at the occasion
};

struct RefinedComplex {
  std::vector<std::array<std::size_t, 2>> vertices;        // (being, occasion)
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> vertex_of;
  std::size_t temporal = 0;                                // edges [0, temporal) are temporal
  std::vector<std::vector<std::size_t>> rung_edge;         // tie -> rung index -> edge
  std::vector<std::array<std::size_t, 2>> rung_of_edge;    // edge - temporal -> (tie, rung index)
  std::vector<RefinedTriangle> triangles;                  // faces [lenses, lenses + triangles)
  morphe::topology::CWComplex cw;
};

inline RefinedComplex refine(const HistoryComplex& K) {
  RefinedComplex R;
  const auto& r = K.record;
  for (std::size_t o : K.order)
    for (const auto& p : r.occasions[o].present) {
      R.vertex_of[{p.being, o}] = R.vertices.size();
      R.vertices.push_back({p.being, o});
    }
  auto& X = R.cw;
  X.name = K.cw.name + " (refined)";
  X.vertices = R.vertices.size();
  for (const auto& s : K.segments) {
    X.edges.push_back({R.vertex_of.at({s.being, s.from}), R.vertex_of.at({s.being, s.to})});
    X.edge_names.push_back(r.beings[s.being].id);
  }
  R.temporal = X.edges.size();
  R.rung_edge.assign(r.ties.size(), {});
  for (std::size_t t = 0; t < r.ties.size(); t++) {
    const Tie& tie = r.ties[t];
    for (std::size_t k = 0; k < K.rungs[t].size(); k++) {
      const std::size_t o = K.rungs[t][k].occasion;
      R.rung_edge[t].push_back(X.edges.size());
      R.rung_of_edge.push_back({t, k});
      X.edges.push_back({R.vertex_of.at({tie.a, o}), R.vertex_of.at({tie.b, o})});
      X.edge_names.push_back(tie.id);
    }
  }
  // The lenses as squares.
  for (const auto& L : K.lenses) {
    const Tie& tie = r.ties[L.tie];
    morphe::topology::Face f;
    f.base = R.vertex_of.at({tie.a, L.from});
    for (std::size_t e : detail::path(K, tie.a, L.from, L.to)) f.boundary.push_back({e, 1});
    f.boundary.push_back({R.rung_edge[L.tie][L.rung_to], 1});
    const auto back = detail::path(K, tie.b, L.from, L.to);
    for (auto it = back.rbegin(); it != back.rend(); ++it) f.boundary.push_back({*it, -1});
    f.boundary.push_back({R.rung_edge[L.tie][L.rung_from], -1});
    X.faces.push_back(std::move(f));
  }
  // The triangles of rungs at each occasion: every three ties joining three beings pairwise.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> at(r.occasions.size());  // occasion -> (tie, rung)
  for (std::size_t t = 0; t < r.ties.size(); t++)
    for (std::size_t k = 0; k < K.rungs[t].size(); k++) at[K.rungs[t][k].occasion].push_back({t, k});
  for (std::size_t o : K.order) {
    const auto& here = at[o];
    for (std::size_t i = 0; i < here.size(); i++)
      for (std::size_t j = i + 1; j < here.size(); j++)
        for (std::size_t k = j + 1; k < here.size(); k++) {
          const std::array<std::pair<std::size_t, std::size_t>, 3> three{here[i], here[j], here[k]};
          std::map<std::size_t, int> degree;
          for (const auto& [t, g] : three) {
            degree[r.ties[t].a]++;
            degree[r.ties[t].b]++;
          }
          if (degree.size() != 3) continue;
          bool cycle = true;
          for (const auto& d : degree) cycle = cycle && d.second == 2;
          if (!cycle) continue;
          // Order the beings x, y, z and the ties x–y, y–z, x–z.
          RefinedTriangle T;
          T.occasion = o;
          std::size_t q = 0;
          for (const auto& d : degree) T.beings[q++] = d.first;
          const auto joins = [&](std::size_t t, std::size_t u, std::size_t v) {
            return (r.ties[t].a == u && r.ties[t].b == v) || (r.ties[t].a == v && r.ties[t].b == u);
          };
          const std::array<std::array<std::size_t, 2>, 3> pairs{{{T.beings[0], T.beings[1]}, {T.beings[1], T.beings[2]}, {T.beings[0], T.beings[2]}}};
          for (std::size_t n = 0; n < 3; n++)
            for (const auto& [t, g] : three)
              if (joins(t, pairs[n][0], pairs[n][1])) { T.ties[n] = t; T.rungs[n] = g; }
          // The boundary x → y → z → x.
          morphe::topology::Face f;
          f.base = R.vertex_of.at({T.beings[0], o});
          const std::array<std::array<std::size_t, 2>, 3> walk{{{T.beings[0], T.beings[1]}, {T.beings[1], T.beings[2]}, {T.beings[2], T.beings[0]}}};
          const std::array<std::size_t, 3> tie_of_step{T.ties[0], T.ties[1], T.ties[2]}, rung_of_step{T.rungs[0], T.rungs[1], T.rungs[2]};
          for (std::size_t n = 0; n < 3; n++) {
            const std::size_t t = tie_of_step[n];
            f.boundary.push_back({R.rung_edge[t][rung_of_step[n]], r.ties[t].a == walk[n][0] ? 1 : -1});
          }
          X.faces.push_back(std::move(f));
          R.triangles.push_back(T);
        }
  }
  return R;
}

// The closed path of K̂ around a tie's meetings k < l: a's worldline from the one to the other, the
// rung at l, b's worldline back, the rung at k reversed. It bounds the tie's lenses from k to l, and
// its holonomy is u_k u_l in temporal gauge: the Wilson loop of that stretch of the tie's life.
inline morphe::topology::Path tie_loop(const HistoryComplex& K, const RefinedComplex& R, std::size_t t, std::size_t k, std::size_t l) {
  const Tie& tie = K.record.ties[t];
  const std::size_t p = K.rungs[t][k].occasion, q = K.rungs[t][l].occasion;
  morphe::topology::Path w;
  for (std::size_t e : detail::path(K, tie.a, p, q)) w.push_back({e, 1});
  w.push_back({R.rung_edge[t][l], 1});
  const auto back = detail::path(K, tie.b, p, q);
  for (auto it = back.rbegin(); it != back.rend(); ++it) w.push_back({*it, -1});
  w.push_back({R.rung_edge[t][k], -1});
  return w;
}

}  // namespace history
