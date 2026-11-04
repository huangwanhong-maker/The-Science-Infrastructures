// The event complex of a record and its refinement (the design, §2): morphe's worldline complexes of the record's
// encounters, in which the couplings are the participants and the primitives the pairs.
//
//   K̂ (refined_worldlines): vertices the presences (c, v), a coupling at an event; temporal edges the segments of a
//      coupling between two consecutive events it attends; rungs (p, v), a primitive at an event that meets it;
//      2-cells the lens squares of each primitive between two consecutive meetings, and the triangles of rungs at an
//      event (two buildings sharing a wall and fronting one street: the hinge where their wall meets the facade).
//   K  (worldline_complex): K̂ with the rungs of each event collapsed; vertices the events, edges the segments, faces
//      the lenses a[v, w] · b[v, w]⁻¹, whose cyclic order is time.
//
// Two boundary events are added, ∂₀ at the window's start and ∂₁ at its end, attended by every cell alive then; their
// links in K are the spin networks of the moment at the window's two ends (the foam with its boundaries coned off).
//
// Attendance (the design's A2): a coupling attends an event that touches it, begins or ends it, or meets, begins or
// ends one of its primitives. A primitive is met where an event names it, begins it or ends it; so an event at an
// interface is an event of both its sides.
//
// The cut at a time t (the spin network of the moment, §2.3): the couplings whose segment crosses t and the primitives
// whose lens crosses t. Homology: b₀ by union-find; b₂ mod 2 by morphe's cycles (the elimination of ∂₂, its rows the
// faces and its columns the edges, bit-packed) where the faces and edges are few enough, and b₁ from the Euler
// characteristic, which is always given.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/topology/cycles.hpp>
#include <morphe/topology/worldlines.hpp>
#include <urban/events/record.hpp>

namespace urban::events {

namespace mt = morphe::topology;

struct EventComplex {
  mt::Encounters encounters;          // participants = couplings, pairs = primitives
  mt::WorldlineComplex K;             // vertex v = event v in time
  mt::RefinedWorldlines refined;      // K̂
  std::vector<long> event_of;         // vertex ↦ the record's event, −1 for ∂₀, −2 for ∂₁
  std::vector<std::size_t> vertex_of; // the record's event ↦ its vertex
  std::vector<double> time;           // vertex ↦ the start of its label (MJD); ∂₀ the window's start, ∂₁ its end
  std::vector<std::array<std::size_t, 2>> life_c, life_p;  // each cell's first and last vertex
};

inline EventComplex construct(const Record& R) {
  if (const auto bad = problems(R); !bad.empty()) throw std::invalid_argument("urban::events::construct: " + bad.front());
  EventComplex X;
  const auto order = event_order(R);
  const auto L = lives(R, order);
  const std::size_t n = order.size(), V = n + 2;  // vertex 0 = ∂₀, 1..n the events, n + 1 = ∂₁
  X.event_of.assign(V, 0);
  X.event_of[0] = -1, X.event_of[V - 1] = -2;
  X.vertex_of.assign(R.events.size(), 0);
  X.time.assign(V, 0);
  X.time[0] = R.window.from, X.time[V - 1] = R.window.to;
  for (std::size_t k = 0; k < n; k++) {
    X.event_of[k + 1] = static_cast<long>(order[k]);
    X.vertex_of[order[k]] = k + 1;
    X.time[k + 1] = R.events[order[k]].when.from;
  }
  const auto vertex = [&](long position) { return static_cast<std::size_t>(position + 1); };  // −1 ↦ ∂₀, n ↦ ∂₁
  X.life_c.resize(R.couplings.size());
  X.life_p.resize(R.primitives.size());
  for (std::size_t c = 0; c < R.couplings.size(); c++) X.life_c[c] = {vertex(L.couplings[c][0]), vertex(L.couplings[c][1])};
  for (std::size_t p = 0; p < R.primitives.size(); p++) X.life_p[p] = {vertex(L.primitives[p][0]), vertex(L.primitives[p][1])};
  // Who attends each vertex, and where each primitive is met.
  std::vector<std::vector<char>> attends(V);
  std::vector<std::vector<std::size_t>> meetings(R.primitives.size());
  const auto attend = [&](std::size_t v, std::size_t c) {
    if (attends[v].empty()) attends[v].assign(R.couplings.size(), 0);
    attends[v][c] = 1;
  };
  for (std::size_t c = 0; c < R.couplings.size(); c++) attend(X.life_c[c][0], c), attend(X.life_c[c][1], c);
  for (std::size_t p = 0; p < R.primitives.size(); p++) {
    meetings[p].push_back(X.life_p[p][0]);
    meetings[p].push_back(X.life_p[p][1]);
  }
  for (std::size_t k = 0; k < n; k++) {
    const auto& e = R.events[order[k]];
    const std::size_t v = k + 1;
    for (auto c : e.couplings) attend(v, c);
    for (auto list : {&e.primitives, &e.begins_primitives, &e.ends_primitives})
      for (auto p : *list) meetings[p].push_back(v);
  }
  for (std::size_t p = 0; p < R.primitives.size(); p++) {
    auto& m = meetings[p];
    std::sort(m.begin(), m.end());
    m.erase(std::unique(m.begin(), m.end()), m.end());
    for (auto v : m) attend(v, R.primitives[p].a), attend(v, R.primitives[p].b);
  }
  auto& E = X.encounters;
  E.participants = R.couplings.size();
  E.present.assign(V, {});
  for (std::size_t v = 0; v < V; v++)
    for (std::size_t c = 0; c < attends[v].size(); c++)
      if (attends[v][c]) E.present[v].push_back(c);
  for (std::size_t p = 0; p < R.primitives.size(); p++) E.pairs.push_back({R.primitives[p].a, R.primitives[p].b, meetings[p]});
  X.K = mt::worldline_complex(E);
  X.refined = mt::refined_worldlines(E, X.K);
  return X;
}

// The spin network of the moment at t (not an event's time): the couplings alive at t, each by the segment that
// crosses it, and the primitives, each by the lens that crosses it.
struct Cut {
  double t = 0;
  std::vector<std::size_t> couplings, segments;   // a coupling and the segment of K crossing t
  std::vector<std::size_t> primitives, lenses;    // a primitive and the lens of K crossing t
};

inline Cut cut(const EventComplex& X, double t) {
  Cut S;
  S.t = t;
  for (std::size_t e = 0; e < X.K.segments.size(); e++) {
    const auto& s = X.K.segments[e];
    if (X.time[s.from] <= t && t < X.time[s.to]) S.couplings.push_back(s.participant), S.segments.push_back(e);
  }
  for (std::size_t f = 0; f < X.K.lenses.size(); f++) {
    const auto& l = X.K.lenses[f];
    if (X.time[l.from] <= t && t < X.time[l.to]) S.primitives.push_back(l.pair), S.lenses.push_back(f);
  }
  return S;
}

// The counts and homology of a complex: cells by dimension, the Euler characteristic, b₀ always, b₁ and b₂ over
// ℤ/2 when its faces number at most `faces_limit` and its edges at most `edges_limit` (the elimination of ∂₂ costs
// about F² E / 64 word operations), and nullopt otherwise.
struct Topology {
  std::array<std::size_t, 3> cells{};
  long euler = 0;
  std::size_t b0 = 0;
  std::optional<std::size_t> b1, b2;
};

inline Topology topology(const mt::CWComplex& C, std::size_t faces_limit = 6000, std::size_t edges_limit = 80000) {
  Topology T;
  T.cells = {C.vertices, C.edges.size(), C.faces.size()};
  T.euler = static_cast<long>(C.vertices) - static_cast<long>(C.edges.size()) + static_cast<long>(C.faces.size());
  std::vector<std::size_t> parent(C.vertices);
  std::iota(parent.begin(), parent.end(), 0);
  const auto root = [&](std::size_t x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; };
  T.b0 = C.vertices;
  for (const auto& e : C.edges) {
    const auto a = root(e[0]), b = root(e[1]);
    if (a != b) parent[a] = b, T.b0--;
  }
  if (C.faces.size() <= faces_limit && C.edges.size() <= edges_limit) {
    // b₂ = dim ker ∂₂ (no 3-cells); b₁ from the Euler characteristic, which holds over any field.
    T.b2 = mt::cycles_mod2(C, 2).basis.size();
    T.b1 = static_cast<std::size_t>(static_cast<long>(T.b0) + static_cast<long>(*T.b2) - T.euler);
  }
  return T;
}

}  // namespace urban::events
