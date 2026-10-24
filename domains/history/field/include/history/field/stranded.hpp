// The history complex as a Feynman graph (theory note §8.1).
//
// A lens of a tie t from meeting p to meeting q is bounded by a[p, q] · b[p, q]⁻¹. Read as a strand,
// it runs along a's segments from p to q, turns at q onto b's, runs back to p and turns again.
// So every segment of a being carries one strand for each lens along it: a field quantum
// φ(g₁, …, gₙ) with one argument for each tie of the being alive across it (its rank n, the n of
// the SU(2) kinematics). Every occasion is an interaction vertex whose strands pair up at its
// corners:
//
//   closing  (in_a, in_b)   where a lens of t ends: t's two beings meet here;
//   opening  (out_a, out_b) where the next lens of t begins;
//   passing  (in_x, out_x)  where a lens runs past an occasion only one of its beings attends.
//
// The graph whose nodes are an occasion's half-edges and whose edges are its corners is the link of
// the vertex in K: the bubble of the interaction. Following strands through propagators and
// corners closes them into exactly the lenses, so K is the 2-complex of this stranded graph.
//
// The Feynman rules of a group field theory on ℤ₂ (theory note §8.2) put a variable g on each
// strand end, a δ(g h_e g′) along each propagator (h_e the segment's own variable, the
// closure of the field), a δ at each passing and opening corner (with the rung u at an opening),
// and a Gibbs kernel e^{β χ(g u g′)} at each closing corner. Eliminating the δ's leaves, for each
// lens, one free variable and a kernel whose argument is the lens's holonomy on K̂: the temporal
// edges of both paths and the rungs at p and q. Hence the amplitude of the history as a Feynman
// graph is 2^{#lenses} 2^{#segments} Z, Z the weight of the spin foam in temporal gauge
// (Proposition 8.1). `feynman_forms` performs that elimination for any record and returns the
// arguments, so the identity is checked on every canon rather than assumed.

#pragma once

#include <history/construction/complex.hpp>
#include <history/foam/refined.hpp>
#include <morphe/combinatorics/canonical.hpp>
#include <physis/gauge/contraction.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace history {

struct HalfEdge {
  std::size_t segment = 0;
  bool out = false;  // the segment leaves this occasion (its `from` is here)
  bool operator==(const HalfEdge&) const = default;
};

struct Corner {
  enum Kind : unsigned char { passing, closing, opening };
  std::size_t lens = 0;
  Kind kind = passing;
  std::array<std::size_t, 2> ends{};  // indices into the bubble's half-edges
};

struct Bubble {
  std::size_t occasion = 0;
  std::vector<std::size_t> beings;     // present, in the occasion's order
  std::vector<HalfEdge> half_edges;
  std::vector<std::size_t> owner;      // half-edge -> its being
  std::vector<Corner> corners;
  std::size_t k_in = 0, k_new = 0, k_last = 0;  // beings arriving from before, first here, last here
  std::size_t meetings = 0, begun = 0, ended = 0;  // rungs here; ties whose first, last rung is here
  std::string type;                    // canonical form of the labelled bubble
  std::size_t automorphisms = 1;
  bool canonical = true;               // false where the search was cut short (then a certificate)
};

struct StrandedGraph {
  std::vector<std::vector<std::size_t>> strands;  // segment -> the lenses along it
  std::vector<Bubble> bubbles;                    // vertex v (the occasion K.order[v]) -> its bubble
  std::map<std::string, std::vector<std::size_t>> types;  // canonical type -> vertices

  std::size_t rank(std::size_t segment) const { return strands[segment].size(); }
};

inline StrandedGraph stranded(const HistoryComplex& K) {
  StrandedGraph G;
  const auto& r = K.record;
  const std::size_t V = K.order.size();
  G.strands.assign(K.segments.size(), {});
  G.bubbles.assign(V, {});
  // The half-edges of each occasion: for each present being, the segment arriving and the one leaving.
  std::vector<std::map<std::pair<std::size_t, bool>, std::size_t>> local(V);  // (segment, out) -> index
  for (std::size_t v = 0; v < V; v++) {
    Bubble& B = G.bubbles[v];
    B.occasion = K.order[v];
    for (const auto& p : r.occasions[B.occasion].present) B.beings.push_back(p.being);
  }
  for (std::size_t e = 0; e < K.segments.size(); e++) {
    const auto& s = K.segments[e];
    for (const bool out : {true, false}) {
      Bubble& B = G.bubbles[K.vertex[out ? s.from : s.to]];
      local[K.vertex[out ? s.from : s.to]][{e, out}] = B.half_edges.size();
      B.half_edges.push_back({e, out});
      B.owner.push_back(s.being);
    }
  }
  // The corners: each lens's boundary, step by step; between step i and step i + 1 it turns at a vertex.
  for (std::size_t f = 0; f < K.cw.faces.size(); f++) {
    const auto& steps = K.cw.faces[f].boundary;
    const std::size_t m = steps.size();
    for (const auto& step : steps) G.strands[step.edge].push_back(f);
    for (std::size_t i = 0; i < m; i++) {
      const auto& a = steps[i];
      const auto& b = steps[(i + 1) % m];
      // Step a ends at the head of its edge if traversed forwards (the segment arrives: an in
      // half-edge), at its tail if backwards (the segment leaves: an out half-edge); step b starts
      // at the tail of its edge if forwards (out), at the head if backwards (in).
      const HalfEdge ha{a.edge, a.sign < 0}, hb{b.edge, b.sign > 0};
      const auto& sa = K.segments[ha.segment];
      const std::size_t v = K.vertex[ha.out ? sa.from : sa.to];
      Corner c;
      c.lens = f;
      c.ends = {local[v].at({ha.segment, ha.out}), local[v].at({hb.segment, hb.out})};
      const bool same = K.segments[ha.segment].being == K.segments[hb.segment].being;
      c.kind = same ? Corner::passing : (ha.out ? Corner::opening : Corner::closing);
      G.bubbles[v].corners.push_back(c);
    }
  }
  // Labels, counts and canonical types.
  std::vector<std::size_t> first_rung(r.ties.size(), 0), last_rung(r.ties.size(), 0);
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> rungs_at(r.occasions.size());  // occasion -> (tie, rung index)
  for (std::size_t t = 0; t < r.ties.size(); t++)
    for (std::size_t k = 0; k < K.rungs[t].size(); k++) rungs_at[K.rungs[t][k].occasion].push_back({t, k});
  for (std::size_t v = 0; v < V; v++) {
    Bubble& B = G.bubbles[v];
    const std::size_t o = B.occasion, k = B.beings.size();
    std::map<std::size_t, std::size_t> index;
    for (std::size_t i = 0; i < k; i++) index[B.beings[i]] = i;
    std::vector<bool> in(k, false), out(k, false);
    for (std::size_t h = 0; h < B.half_edges.size(); h++) (B.half_edges[h].out ? out : in)[index.at(B.owner[h])] = true;
    std::vector<std::vector<std::string>> passing(k);
    for (const auto& c : B.corners)
      if (c.kind == Corner::passing) passing[index.at(B.owner[c.ends[0]])].push_back(r.ties[K.lenses[c.lens].tie].colour);
    std::vector<std::string> label(k);
    for (std::size_t i = 0; i < k; i++) {
      std::sort(passing[i].begin(), passing[i].end());
      label[i] = std::string(in[i] ? "i" : "") + (out[i] ? "o" : "") + "[";
      for (const auto& x : passing[i]) label[i] += x + " ";
      label[i] += "]";
      B.k_in += in[i];
      B.k_new += !in[i];
      B.k_last += !out[i];
    }
    // The bubble as a labelled multigraph: its beings labelled by (arrives, leaves, the colours of
    // the strands passing them), its meetings by (colour, a lens closes, a lens opens); its type
    // is morphe's canonical form, and the order of its automorphism group comes with it.
    morphe::combinatorics::LabelledGraph bubble;
    for (const auto& [t, g] : rungs_at[o]) {
      const bool closes = g > 0, opens = g + 1 < K.rungs[t].size();
      B.meetings++;
      B.begun += g == 0;
      B.ended += !opens;
      const std::string e = r.ties[t].colour + "/" + (closes ? "c" : "") + (opens ? "o" : "");
      bubble.edges.push_back({index.at(r.ties[t].a), index.at(r.ties[t].b), e});
    }
    bubble.vertex = label;
    const auto form = morphe::combinatorics::canonical_form(bubble);
    B.type = form.code;
    B.automorphisms = form.automorphisms;
    B.canonical = form.exact;
    G.types[B.type].push_back(v);
  }
  return G;
}

// Every strand end at every occasion lies on exactly one corner: the pairing of a vertex is
// complete, as an interaction vertex's must be.
inline bool pairing_complete(const HistoryComplex& K, const StrandedGraph& G) {
  for (const auto& B : G.bubbles) {
    std::map<std::pair<std::size_t, std::size_t>, int> used;  // (half-edge, lens) -> corners
    for (const auto& c : B.corners)
      for (std::size_t h : c.ends) used[{h, c.lens}]++;
    for (std::size_t h = 0; h < B.half_edges.size(); h++)
      for (std::size_t f : G.strands[B.half_edges[h].segment])
        if (used[{h, f}] != 1) return false;
    for (const auto& [key, n] : used)
      if (n != 1) return false;
  }
  (void)K;
  return true;
}

// The Feynman rules over ℤ₂ eliminated: the argument of each lens's kernel as a set of variables
// (segments h_e first, then the rungs of K̂ in its order), and the number of free strand variables
// left. The arguments must be the lenses' boundaries in K̂, and the free variables one per lens.
struct FeynmanForms {
  std::vector<std::vector<std::size_t>> argument;  // lens -> variables (sorted)
  std::size_t free_strands = 0;
  bool consistent = true;  // no δ constrained the segment and rung variables
};

inline FeynmanForms feynman_forms(const HistoryComplex& K, const RefinedComplex& R, const StrandedGraph& G) {
  using Bits = std::vector<std::uint64_t>;
  const std::size_t E = K.segments.size(), nvars = R.cw.edges.size();  // h_e for e < E; rungs after
  const std::size_t words = (nvars + 63) / 64;
  auto unit = [&](std::size_t i) {
    Bits b(words, 0);
    b[i / 64] |= std::uint64_t{1} << (i % 64);
    return b;
  };
  auto add = [](Bits& a, const Bits& b) {
    for (std::size_t w = 0; w < a.size(); w++) a[w] ^= b[w];
  };
  // Strand variables: (vertex, half-edge, lens). Union-find with the offset of each to its root.
  std::map<std::array<std::size_t, 3>, std::size_t> id;
  std::vector<std::size_t> parent;
  std::vector<Bits> offset;  // g = g_root + offset
  auto var = [&](std::size_t v, std::size_t h, std::size_t f) {
    const std::array<std::size_t, 3> key{v, h, f};
    const auto it = id.find(key);
    if (it != id.end()) return it->second;
    id[key] = parent.size();
    parent.push_back(parent.size());
    offset.push_back(Bits(words, 0));
    return parent.size() - 1;
  };
  std::function<std::size_t(std::size_t)> find = [&](std::size_t x) {
    if (parent[x] == x) return x;
    const std::size_t root = find(parent[x]);
    if (parent[x] != root) {
      add(offset[x], offset[parent[x]]);
      parent[x] = root;
    }
    return root;
  };
  FeynmanForms out;
  // g_x + c + g_y = 0
  auto relate = [&](std::size_t x, std::size_t y, const Bits& c) {
    const std::size_t rx = find(x), ry = find(y);
    Bits d = offset[x];
    add(d, offset[y]);
    add(d, c);  // g_rx + g_ry = d
    if (rx == ry) {
      if (std::any_of(d.begin(), d.end(), [](std::uint64_t w) { return w != 0; })) out.consistent = false;
      return;
    }
    parent[rx] = ry;
    offset[rx] = d;
  };
  auto local = [&](std::size_t v, const HalfEdge& he) {
    const auto& B = G.bubbles[v];
    for (std::size_t h = 0; h < B.half_edges.size(); h++)
      if (B.half_edges[h] == he) return h;
    return B.half_edges.size();
  };
  // Propagators: each strand along a segment joins its two ends through h_e.
  for (std::size_t e = 0; e < E; e++) {
    const auto& s = K.segments[e];
    const std::size_t vf = K.vertex[s.from], vt = K.vertex[s.to];
    const std::size_t hf = local(vf, {e, true}), ht = local(vt, {e, false});
    for (std::size_t f : G.strands[e]) relate(var(vf, hf, f), var(vt, ht, f), unit(e));
  }
  // The rung of K̂ at a meeting of a lens's tie.
  auto rung_var = [&](std::size_t f, std::size_t o) {
    const Lens& L = K.lenses[f];
    const std::size_t k = K.rungs[L.tie][L.rung_from].occasion == o ? L.rung_from : L.rung_to;
    return R.rung_edge[L.tie][k];
  };
  // Corners: δ at passing and opening corners; the kernel's argument at closing corners, read after.
  std::vector<std::array<std::size_t, 4>> closing;  // (x, y, rung variable, lens)
  for (std::size_t v = 0; v < G.bubbles.size(); v++)
    for (const auto& c : G.bubbles[v].corners) {
      const std::size_t x = var(v, c.ends[0], c.lens), y = var(v, c.ends[1], c.lens);
      if (c.kind == Corner::passing) relate(x, y, Bits(words, 0));
      else if (c.kind == Corner::opening) relate(x, y, unit(rung_var(c.lens, G.bubbles[v].occasion)));
      else closing.push_back({x, y, rung_var(c.lens, G.bubbles[v].occasion), c.lens});
    }
  out.argument.assign(K.lenses.size(), {});
  for (const auto& [x, y, u, f] : closing) {
    if (find(x) != find(y)) { out.consistent = false; continue; }
    Bits a = offset[x];
    add(a, offset[y]);
    add(a, unit(u));
    for (std::size_t i = 0; i < nvars; i++)
      if ((a[i / 64] >> (i % 64)) & 1) out.argument[f].push_back(i);
  }
  for (std::size_t x = 0; x < parent.size(); x++) out.free_strands += find(x) == x;
  return out;
}

// The boundary of each lens in K̂ as a set of edges (temporal edges, then rungs): the links its
// holonomy multiplies in physis's ℤ₂ theory of K̂, for comparison.
inline std::vector<std::vector<std::size_t>> lens_boundaries(const RefinedComplex& R, std::size_t lenses) {
  auto T = physis::gauge::Z2Theory::of(R.cw, std::vector<double>(R.cw.faces.size(), 0));
  T.plaquettes.resize(lenses);
  return T.plaquettes;
}

}  // namespace history
