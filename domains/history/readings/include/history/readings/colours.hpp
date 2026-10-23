// The disciplines of a history as colours of its faces (codebook §6; theory note §2): the
// history of a discipline is the sub-complex of the lenses of its colour, and the integrated
// history the whole. Beings that bound lenses of several colours carry the disciplines into
// one another; occasions where lenses of several colours meet are interdisciplinary by
// construction, not by a coder's label.

#pragma once

#include <history/construction/complex.hpp>

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace history {

struct Colours {
  std::map<std::string, std::size_t> lenses;  // the lenses of each colour
  std::vector<std::pair<std::size_t, std::map<std::string, std::size_t>>> carriers;  // beings bounding several colours
  std::vector<std::pair<std::size_t, std::vector<std::string>>> meetings;            // occasions where several meet
  std::vector<std::string> of_being;  // each being's colour: that of most of the lenses it bounds, else of its ties
};

inline Colours colours_of(const HistoryComplex& K) {
  Colours C;
  const auto& r = K.record;
  std::vector<std::map<std::string, std::size_t>> bounds(r.beings.size());
  std::vector<std::set<std::string>> at(r.occasions.size());
  for (const auto& L : K.lenses) {
    const Tie& t = r.ties[L.tie];
    C.lenses[t.colour]++;
    bounds[t.a][t.colour]++;
    bounds[t.b][t.colour]++;
    at[L.from].insert(t.colour);
    at[L.to].insert(t.colour);
  }
  for (std::size_t b = 0; b < r.beings.size(); b++)
    if (bounds[b].size() > 1) C.carriers.push_back({b, bounds[b]});
  for (std::size_t o : K.order)
    if (at[o].size() > 1) C.meetings.push_back({o, {at[o].begin(), at[o].end()}});
  // A being with no lens takes the colour of its ties.
  std::vector<std::map<std::string, std::size_t>> ties_of(r.beings.size());
  for (const auto& t : r.ties) {
    ties_of[t.a][t.colour]++;
    ties_of[t.b][t.colour]++;
  }
  C.of_being.assign(r.beings.size(), "");
  for (std::size_t b = 0; b < r.beings.size(); b++) {
    const auto& m = bounds[b].empty() ? ties_of[b] : bounds[b];
    std::size_t best = 0;
    for (const auto& colour : colours())  // ties go to the codebook's earlier colour
      if (m.count(colour) && m.at(colour) > best) { best = m.at(colour); C.of_being[b] = colour; }
  }
  return C;
}

}  // namespace history
