// The construction of the history complex (theory note §2): from a record, with no judgement,
// the CW 2-complex K whose
//
//   vertices are the occasions, in the order of time;
//   edges are the beings between consecutive occasions they attend (their worldlines);
//   2-cells are the lenses: for each tie and each two consecutive meetings of its two beings
//   within its life, the disc bounded by the two worldlines between those meetings, with the
//   attaching word a[p, q] · b[p, q]⁻¹, pinched at the meetings.
//
// A tie lives from its first attestation to its last, or to its end (codebook §8). At each
// meeting of its beings within its life it has a rung: its sign there, attested, or carried
// from its last attestation (the latent value of the theory note §3, which the readings weigh
// by the persistence of the tie's colour).

#pragma once

#include <history/record/record.hpp>
#include <morphe/topology/cw_complex.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace history {

// An edge: a being between two consecutive occasions it attends.
struct Segment {
  std::size_t being = 0, from = 0, to = 0;  // occasions
};

// A tie's value at a meeting of its two beings within its life.
struct Rung {
  std::size_t occasion = 0;
  int sign = 1;                  // attested here, or carried from the last attestation
  bool attested = false;
  std::size_t since = 0;         // the occasion of the last attestation
  std::size_t lenses_since = 0;  // lenses crossed since it (0 when attested here)
};

// A 2-cell: a tie between two consecutive meetings of its beings.
struct Lens {
  std::size_t tie = 0, from = 0, to = 0;  // occasions
  std::size_t rung_from = 0, rung_to = 0; // indices into the tie's rungs
};

struct HistoryComplex {
  Record record;                                   // as read, with the bereaved derived (codebook §8)
  std::vector<std::size_t> order;                  // vertex v is the occasion order[v]
  std::vector<std::size_t> vertex;                 // the vertex of each occasion
  std::vector<std::vector<std::size_t>> worldline; // each being's occasions in time
  std::vector<Segment> segments;                   // edge e
  std::vector<std::vector<std::size_t>> segments_of;  // each being's edges in time
  std::vector<std::array<std::size_t, 2>> life;    // each tie's first and last occasion
  std::vector<std::vector<Rung>> rungs;            // each tie's meetings within its life
  std::vector<Lens> lenses;                        // face f
  std::vector<std::vector<std::size_t>> lenses_of; // each tie's faces in time
  morphe::topology::CWComplex cw;
  std::vector<Problem> problems;                   // the construction's own warnings

  bool flip(std::size_t f) const {
    const Lens& L = lenses[f];
    return rungs[L.tie][L.rung_from].sign != rungs[L.tie][L.rung_to].sign;
  }
};

namespace detail {

// The edges of a being's worldline from occasion p to occasion q (p before q).
inline std::vector<std::size_t> path(const HistoryComplex& K, std::size_t being, std::size_t p, std::size_t q) {
  const auto& w = K.worldline[being];
  std::size_t i = 0, j = 0;
  for (std::size_t k = 0; k < w.size(); k++) {
    if (w[k] == p) i = k;
    if (w[k] == q) j = k;
  }
  return {K.segments_of[being].begin() + static_cast<std::ptrdiff_t>(i), K.segments_of[being].begin() + static_cast<std::ptrdiff_t>(j)};
}

}  // namespace detail

// Builds the complex of a record. A record with errors is built as far as it can be; its
// problems travel with it.
inline HistoryComplex construct(const Record& record) {
  HistoryComplex K;
  K.record = record;
  Record& r = K.record;
  add_bereaved(r);
  const auto rank = r.ranks();

  // Vertices: the occasions in time.
  K.order = r.order;
  K.vertex.assign(r.occasions.size(), 0);
  for (std::size_t v = 0; v < K.order.size(); v++) K.vertex[K.order[v]] = v;

  // Edges: each being's consecutive occasions.
  K.worldline.assign(r.beings.size(), {});
  for (std::size_t o : K.order)
    for (const auto& p : r.occasions[o].present) K.worldline[p.being].push_back(o);
  K.segments_of.assign(r.beings.size(), {});
  for (std::size_t b = 0; b < r.beings.size(); b++) {
    const auto& w = K.worldline[b];
    for (std::size_t k = 0; k + 1 < w.size(); k++) {
      K.segments_of[b].push_back(K.segments.size());
      K.segments.push_back({b, w[k], w[k + 1]});
      // Two occasions of one being that overlap in time are ordered by the file (codebook §4).
      if (r.occasions[w[k + 1]].date.first <= r.occasions[w[k]].date.last)
        K.problems.push_back({false, r.beings[b].id,
                              "the order of \"" + r.occasions[w[k]].id + "\" and \"" + r.occasions[w[k + 1]].id + "\" rests on the file: their dates overlap"});
    }
  }

  // Each tie's life, its rungs at the meetings within it, and its lenses.
  K.life.assign(r.ties.size(), {0, 0});
  K.rungs.assign(r.ties.size(), {});
  K.lenses_of.assign(r.ties.size(), {});
  for (std::size_t t = 0; t < r.ties.size(); t++) {
    const Tie& tie = r.ties[t];
    if (tie.attested.empty()) continue;
    const std::size_t first = tie.attested.front().occasion;
    const std::size_t last = tie.ends && rank[*tie.ends] > rank[tie.attested.back().occasion] ? *tie.ends : tie.attested.back().occasion;
    K.life[t] = {first, last};
    std::size_t a = 0;  // the next attestation
    Rung carried;
    for (std::size_t o : K.worldline[tie.a]) {
      if (rank[o] < rank[first] || rank[o] > rank[last]) continue;
      if (!is_present(r.occasions[o], tie.b)) continue;
      Rung g;
      g.occasion = o;
      if (a < tie.attested.size() && tie.attested[a].occasion == o) {
        g.sign = tie.attested[a].sign;
        g.attested = true;
        g.since = o;
        g.lenses_since = 0;
        a++;
      } else {
        g.sign = carried.sign;
        g.since = carried.since;
        g.lenses_since = carried.lenses_since + 1;
      }
      carried = g;
      K.rungs[t].push_back(g);
    }
    for (std::size_t k = 0; k + 1 < K.rungs[t].size(); k++) {
      K.lenses_of[t].push_back(K.lenses.size());
      K.lenses.push_back({t, K.rungs[t][k].occasion, K.rungs[t][k + 1].occasion, k, k + 1});
    }
  }

  // The CW complex.
  auto& X = K.cw;
  X.name = r.canon.id.empty() ? "history" : r.canon.id;
  X.vertices = K.order.size();
  for (const auto& s : K.segments) {
    X.edges.push_back({K.vertex[s.from], K.vertex[s.to]});
    X.edge_names.push_back(r.beings[s.being].id);
  }
  for (const auto& L : K.lenses) {
    const Tie& tie = r.ties[L.tie];
    morphe::topology::Face f;
    f.base = K.vertex[L.from];
    for (std::size_t e : detail::path(K, tie.a, L.from, L.to)) f.boundary.push_back({e, 1});
    const auto back = detail::path(K, tie.b, L.from, L.to);
    for (auto it = back.rbegin(); it != back.rend(); ++it) f.boundary.push_back({*it, -1});
    X.faces.push_back(std::move(f));
  }
  return K;
}

}  // namespace history
