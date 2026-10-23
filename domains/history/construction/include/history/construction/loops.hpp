// The loops of a history (theory note §4): each generator of H₁ of the complex, a 1-cycle,
// told as closed walks along the worldlines: the beings it follows, and the occasions where
// it passes from one being to another. A loop is a circle of encounters that no recorded tie
// fills: where the record is silent, and where a flat connection may wind.

#pragma once

#include <history/construction/complex.hpp>
#include <morphe/topology/invariants.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

namespace history {

struct LoopStep {
  std::size_t edge = 0;
  bool forward = true;  // with time, along the being's worldline
};

struct Loop {
  std::vector<LoopStep> steps;         // closed walks, one after another
  std::vector<std::size_t> beings;     // those it follows, in order of first passage
  std::vector<std::size_t> occasions;  // where it passes from one being to another
};

// A 1-cycle (edge, coefficient) as closed walks. Each edge with coefficient c is |c| steps in
// the direction of its sign; at each vertex the walk prefers to stay with the being it follows.
inline Loop loop_of(const HistoryComplex& K, const morphe::algebra::Chain& cycle) {
  struct Arrow { std::size_t edge, tail, head; bool forward; bool used = false; };
  std::vector<Arrow> arrows;
  for (const auto& entry : cycle) {
    const long long c = entry.value.to_int64();
    const auto& e = K.cw.edges[entry.index];
    for (long long k = 0; k < (c < 0 ? -c : c); k++)
      arrows.push_back(c > 0 ? Arrow{entry.index, e[0], e[1], true} : Arrow{entry.index, e[1], e[0], false});
  }
  std::multimap<std::size_t, std::size_t> out;  // tail vertex -> arrow
  for (std::size_t i = 0; i < arrows.size(); i++) out.insert({arrows[i].tail, i});
  Loop L;
  auto being_of = [&](std::size_t edge) { return K.segments[edge].being; };
  for (;;) {
    // Start each closed walk at the earliest vertex with an unused arrow.
    std::size_t start = arrows.size();
    for (std::size_t i = 0; i < arrows.size(); i++)
      if (!arrows[i].used && (start == arrows.size() || arrows[i].tail < arrows[start].tail)) start = i;
    if (start == arrows.size()) break;
    std::size_t at = start;
    const std::size_t home = arrows[start].tail;
    while (true) {
      Arrow& a = arrows[at];
      a.used = true;
      L.steps.push_back({a.edge, a.forward});
      if (a.head == home) {
        bool more = false;
        for (auto [it, end] = out.equal_range(home); it != end; ++it) more = more || !arrows[it->second].used;
        if (!more) break;
      }
      std::size_t next = arrows.size();
      for (auto [it, end] = out.equal_range(a.head); it != end; ++it) {
        if (arrows[it->second].used) continue;
        if (next == arrows.size() || being_of(arrows[it->second].edge) == being_of(a.edge)) next = it->second;
        if (being_of(arrows[next].edge) == being_of(a.edge)) break;
      }
      if (next == arrows.size()) break;  // not a cycle: stop where it ends
      at = next;
    }
  }
  // The beings in order of passage, and the occasions where the walk changes being.
  for (std::size_t k = 0; k < L.steps.size(); k++) {
    const auto& s = L.steps[k];
    const std::size_t b = being_of(s.edge);
    if (std::find(L.beings.begin(), L.beings.end(), b) == L.beings.end()) L.beings.push_back(b);
    const auto& next = L.steps[(k + 1) % L.steps.size()];
    if (being_of(next.edge) != b) {
      const Segment& seg = K.segments[s.edge];
      const std::size_t o = s.forward ? seg.to : seg.from;
      if (std::find(L.occasions.begin(), L.occasions.end(), o) == L.occasions.end()) L.occasions.push_back(o);
    }
  }
  return L;
}

// The loops of the generators of H₁.
inline std::vector<Loop> loops(const HistoryComplex& K, const std::vector<morphe::algebra::HomologyGroup>& H) {
  std::vector<Loop> out;
  if (H.size() < 2) return out;
  for (const auto& g : H[1].free_generators) out.push_back(loop_of(K, g));
  for (const auto& g : H[1].torsion_generators) out.push_back(loop_of(K, g));
  return out;
}

}  // namespace history
