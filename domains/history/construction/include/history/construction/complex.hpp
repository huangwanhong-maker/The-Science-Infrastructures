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
//
// The complex itself is morphe's worldline complex of the record's encounters (the occasions in
// time, those present at each, and each tie's meetings within its life): the construction names
// no object of history, so it is written there, and the domain supplies what is the record's.

#pragma once

#include <history/record/record.hpp>
#include <morphe/topology/cw_complex.hpp>
#include <morphe/topology/worldlines.hpp>

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
  morphe::topology::Encounters encounters;       // the record as morphe reads it: vertices, presences, meetings
  morphe::topology::WorldlineComplex worldlines;  // morphe's construction, of which cw is the complex
  std::vector<Problem> problems;                   // the construction's own warnings

  bool flip(std::size_t f) const {
    const Lens& L = lenses[f];
    return rungs[L.tie][L.rung_from].sign != rungs[L.tie][L.rung_to].sign;
  }
};

namespace detail {

// The edges of a being's worldline from occasion p to occasion q (p before q).
inline std::vector<std::size_t> path(const HistoryComplex& K, std::size_t being, std::size_t p, std::size_t q) {
  return K.worldlines.path(being, K.vertex[p], K.vertex[q]);
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

  // The encounters: the occasions in time and those present at each. morphe's worldline
  // complex gives each being's worldline and its segments.
  K.encounters.participants = r.beings.size();
  for (std::size_t o : K.order) {
    std::vector<std::size_t> here;
    for (const auto& p : r.occasions[o].present) here.push_back(p.being);
    K.encounters.present.push_back(std::move(here));
  }
  K.worldlines = morphe::topology::worldline_complex(K.encounters);
  K.worldline.assign(r.beings.size(), {});
  for (std::size_t b = 0; b < r.beings.size(); b++)
    for (std::size_t v : K.worldlines.worldline[b]) K.worldline[b].push_back(K.order[v]);
  K.segments_of = K.worldlines.segments_of;
  for (const auto& s : K.worldlines.segments) {
    K.segments.push_back({s.participant, K.order[s.from], K.order[s.to]});
    // Two occasions of one being that overlap in time are ordered by the file (codebook §4).
    const auto& from = r.occasions[K.order[s.from]];
    const auto& to = r.occasions[K.order[s.to]];
    if (to.date.first <= from.date.last)
      K.problems.push_back({false, r.beings[s.participant].id, "the order of \"" + from.id + "\" and \"" + to.id + "\" rests on the file: their dates overlap"});
  }

  // Each tie's life, and its rungs at the meetings within it.
  K.life.assign(r.ties.size(), {0, 0});
  K.rungs.assign(r.ties.size(), {});
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
  }

  // The lenses and the CW complex: morphe's worldline complex with the ties' meetings.
  for (std::size_t t = 0; t < r.ties.size(); t++) {
    morphe::topology::Encounters::Pair pair;
    pair.a = r.ties[t].a;
    pair.b = r.ties[t].b;
    for (const auto& g : K.rungs[t]) pair.meetings.push_back(K.vertex[g.occasion]);
    K.encounters.pairs.push_back(std::move(pair));
  }
  K.worldlines = morphe::topology::worldline_complex(K.encounters);
  K.lenses_of = K.worldlines.lenses_of;
  for (const auto& L : K.worldlines.lenses) K.lenses.push_back({L.pair, K.order[L.from], K.order[L.to], L.k, L.k + 1});
  K.cw = K.worldlines.cw;
  K.cw.name = r.canon.id.empty() ? "history" : r.canon.id;
  for (const auto& sg : K.segments) K.cw.edge_names.push_back(r.beings[sg.being].id);
  return K;
}

}  // namespace history
