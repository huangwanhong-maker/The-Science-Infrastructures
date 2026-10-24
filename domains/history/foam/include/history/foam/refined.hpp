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
// It is morphe's refined worldline complex of the record's encounters. Collapsing the rungs of
// each occasion gives back the history complex K. A frustrated triangle
// of K̂ is curvature at an occasion that K, with its rungs collapsed, cannot carry: that is why
// the signs of ties are read on K̂ (theory note §6).

#pragma once

#include <history/construction/complex.hpp>
#include <morphe/topology/cw_complex.hpp>
#include <morphe/topology/worldlines.hpp>

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
  const auto Fw = morphe::topology::refined_worldlines(K.encounters, K.worldlines);
  const auto& r = K.record;
  RefinedComplex R;
  for (const auto& [x, v] : Fw.presences) {
    R.vertex_of[{x, K.order[v]}] = R.vertices.size();
    R.vertices.push_back({x, K.order[v]});
  }
  R.temporal = Fw.temporal;
  R.rung_edge = Fw.rung_edge;
  R.rung_of_edge = Fw.rung_of_edge;
  for (const auto& T : Fw.triangles) R.triangles.push_back({K.order[T.event], T.participants, T.pairs, T.meetings});
  R.cw = Fw.cw;
  R.cw.name = K.cw.name + " (refined)";
  for (const auto& sg : K.segments) R.cw.edge_names.push_back(r.beings[sg.being].id);
  for (const auto& [t, k] : Fw.rung_of_edge) R.cw.edge_names.push_back(r.ties[t].id);
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
