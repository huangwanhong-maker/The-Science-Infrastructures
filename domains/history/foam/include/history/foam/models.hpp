// Other complexes a record makes (theory note §7), beside the history complex K and its
// refinement K̂:
//
//   K³, the history complex with its closed shells filled: a 3-cell for each three mutually
//     tied beings carried between two meetings of all three, a collective episode held whole;
//   the Dowker complex of the presence relation (Dowker 1952): the simplicial complex on the
//     beings whose simplices are the sets of beings present together at an occasion, which is
//     Atkin's complex of a relation (Atkin 1974), with its Q-analysis; by Dowker's theorem it is
//     homotopy equivalent to the complex on the occasions whose simplices are the sets of
//     occasions a being attends, the duality of persons and groups (Breiger 1974) made
//     topological;
//   the history in time: K filtered by the order of its occasions, and its persistent homology
//     over ℤ/2: when a circle of encounters closes, and when a later meeting fills it.

#pragma once

#include <history/construction/complex.hpp>
#include <history/readings/signed.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/persistence.hpp>
#include <morphe/topology/q_analysis.hpp>
#include <morphe/topology/simplicial.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace history {

// K³: K with a 3-cell on each minimal closed shell. The boundary of the 3-cell is the sum of
// the shell's three lenses with the signs that make it a cycle (∂∂ = 0).
inline morphe::topology::CWComplex shells_filled(const HistoryComplex& K, const std::vector<Shell>& S) {
  morphe::topology::CWComplex X = K.cw;
  X.name = K.cw.name + " (shells filled)";
  std::vector<morphe::algebra::Vector> cells;
  for (const auto& s : S) {
    // ∂(lens of tie (a, b)) = P_a − P_b, P the being's path between the two meetings.
    for (int signs = 0; signs < 8; signs++) {
      std::map<std::size_t, int> sum;
      std::array<int, 3> c{};
      for (std::size_t q = 0; q < 3; q++) {
        c[q] = (signs >> q) & 1 ? -1 : 1;
        const Tie& t = K.record.ties[K.lenses[s.faces[q]].tie];
        sum[t.a] += c[q];
        sum[t.b] -= c[q];
      }
      if (std::all_of(sum.begin(), sum.end(), [](const auto& x) { return x.second == 0; })) {
        morphe::algebra::Vector v;
        std::array<std::size_t, 3> order{0, 1, 2};
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return s.faces[a] < s.faces[b]; });
        for (std::size_t q : order) v.push_back({s.faces[q], morphe::arithmetic::Integer(static_cast<long long>(c[q]))});
        cells.push_back(v);
        break;
      }
    }
  }
  X.higher.push_back(cells);
  return X;
}

// The Dowker complex on the beings: a simplex for the beings present at each occasion (derived
// presences left out: the bereaved did not meet).
inline morphe::topology::SimplicialComplex dowker_beings(const HistoryComplex& K) {
  std::vector<morphe::topology::Simplex> facets;
  for (const auto& o : K.record.occasions) {
    morphe::topology::Simplex s;
    for (const auto& p : o.present)
      if (!p.derived) s.push_back(p.being);
    std::sort(s.begin(), s.end());
    if (!s.empty()) facets.push_back(s);
  }
  return morphe::topology::SimplicialComplex::from_facets(facets, K.record.beings.size(), K.cw.name + " (Dowker, beings)");
}

// The Dowker complex on the occasions: a simplex for the occasions each being attends. Its
// simplices grow as 2^(occasions of a being): only for small records.
inline morphe::topology::SimplicialComplex dowker_occasions(const HistoryComplex& K) {
  std::vector<morphe::topology::Simplex> facets;
  std::vector<morphe::topology::Simplex> of(K.record.beings.size());
  for (std::size_t o = 0; o < K.record.occasions.size(); o++)
    for (const auto& p : K.record.occasions[o].present)
      if (!p.derived) of[p.being].push_back(o);
  for (auto& s : of)
    if (!s.empty()) { std::sort(s.begin(), s.end()); facets.push_back(s); }
  return morphe::topology::SimplicialComplex::from_facets(facets, K.record.occasions.size(), K.cw.name + " (Dowker, occasions)");
}

// Atkin's Q-analysis of the occasions as simplices on the beings present (the derived bereaved
// left out): morphe's q_analysis of that family (Atkin 1974, ch. 3). Q[q] for q = 0 .. the
// largest dimension.
inline std::vector<std::size_t> q_analysis(const HistoryComplex& K) {
  std::vector<morphe::topology::Simplex> family;
  for (const auto& o : K.record.occasions) {
    morphe::topology::Simplex s;
    for (const auto& p : o.present)
      if (!p.derived) s.push_back(p.being);
    if (!s.empty()) family.push_back(std::move(s));
  }
  return morphe::topology::q_analysis(std::move(family));
}

// The history in time: each cell of K enters at the rank of its last occasion (a vertex at its
// occasion, a segment when its being reaches the later occasion, a lens at the later meeting).
inline morphe::topology::Filtration time_filtration(const HistoryComplex& K) {
  const auto rank = K.record.ranks();
  morphe::topology::Filtration F;
  for (std::size_t v = 0; v < K.order.size(); v++) F.vertex.push_back(rank[K.order[v]]);
  for (const auto& s : K.segments) F.edge.push_back(rank[s.to]);
  for (const auto& L : K.lenses) F.face.push_back(rank[L.to]);
  return F;
}

// Persistence of the history in time; stages are the ranks of occasions, so a stage names the
// occasion at which a class is born or dies (K.order[stage]).
inline morphe::topology::PersistenceDiagram persistence_in_time(const HistoryComplex& K) {
  return morphe::topology::persistent_homology(K.cw, time_filtration(K));
}

}  // namespace history
