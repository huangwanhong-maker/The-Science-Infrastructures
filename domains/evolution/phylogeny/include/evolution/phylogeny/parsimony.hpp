// Parsimony: the least number of changes (or least total cost) by which a tree explains the
// sequences at its leaves, and the search for the trees that need the fewest.
//
// - Fitch's algorithm (Fitch 1971, "Toward defining the course of evolution: minimum change
//   for a specific tree topology", Systematic Zoology 20, 406–416), in Hartigan's form for
//   nodes of any degree (Hartigan 1973, "Minimum mutation fits to a given tree", Biometrics 29,
//   53–65): from the leaves up, a node's set is the states found in the most children's sets,
//   and it costs the number of children whose sets lack them. For two children: their
//   intersection at no cost, or their union at a cost of one. The total is the least number of
//   changes over every assignment of states to the internal nodes, for each site; O(n · s).
// - Sankoff's algorithm (Sankoff 1975, "Minimal mutation trees of sequences", SIAM J. Appl.
//   Math. 28, 35–42; Sankoff and Rousseau 1975): S_v(x) = Σ_c min_y [c(x, y) + S_c(y)], the
//   least cost of v's subtree given state x at v, for any cost matrix; the score is
//   min_x S_root(x), and a traceback from the root gives one assignment that attains it.
//   O(n · s · 16).
//
//   Both are dynamic programmes over the tree in the tropical semiring (min, +): Felsenstein's
//   pruning with sums replaced by minima and products by sums (theory/phylogeny.md, §7).
//
// - The search: stepwise addition (each taxon in turn placed on the branch where it costs
//   least; Felsenstein 2004, Inferring Phylogenies, ch. 4), then nearest-neighbour interchanges
//   while one lowers the score; and an exact search by branch and bound over every order of
//   stepwise addition, pruning a partial tree as soon as its score reaches the best found,
//   since adding a taxon never lowers a score (Hendy and Penny 1982, "Branch and bound
//   algorithms to determine minimal evolutionary trees", Mathematical Biosciences 59,
//   277–290). The heuristic is O(n³ · s); the exact search is exponential, for a dozen taxa.
//
// Gaps and N are missing data: any state, at no cost. Trees are unrooted (the score does not
// depend on the root); leaves are matched to rows by label.

#pragma once

#include <evolution/phylogeny/alignment.hpp>
#include <evolution/phylogeny/tree.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace evolution::phylogeny {

namespace detail {

inline std::vector<std::size_t> rows_of_leaves(const Tree& t, const Patterns& p) {
  std::vector<std::size_t> row(t.size(), SIZE_MAX);
  for (std::size_t v = 0; v < t.size(); v++) {
    if (!t.nodes[v].children.empty()) continue;
    row[v] = p.index(t.nodes[v].label);
    if (row[v] == SIZE_MAX) throw std::invalid_argument("evolution::phylogeny: the leaf " + t.nodes[v].label + " has no sequence");
  }
  return row;
}

// Fitch–Hartigan with leaves given by rows; `sets` is scratch of nodes × patterns.
inline double fitch(const Tree& t, const Patterns& p, const std::vector<std::size_t>& row, std::vector<std::uint8_t>& sets) {
  const std::size_t P = p.count;
  sets.resize(t.size() * P);
  double score = 0;
  for (int v : t.postorder()) {
    std::uint8_t* s = &sets[static_cast<std::size_t>(v) * P];
    const auto& ch = t[v].children;
    if (ch.empty()) {
      const std::size_t r = row[static_cast<std::size_t>(v)];
      for (std::size_t q = 0; q < P; q++) {
        const Code c = p.states[q * p.taxa + r];
        s[q] = c < 4 ? static_cast<std::uint8_t>(1u << c) : std::uint8_t{15};
      }
      continue;
    }
    if (ch.size() == 2) {
      const std::uint8_t* a = &sets[static_cast<std::size_t>(ch[0]) * P];
      const std::uint8_t* b = &sets[static_cast<std::size_t>(ch[1]) * P];
      for (std::size_t q = 0; q < P; q++) {
        const std::uint8_t i = a[q] & b[q];
        if (i) s[q] = i;
        else {
          s[q] = a[q] | b[q];
          score += p.weight[q];
        }
      }
      continue;
    }
    for (std::size_t q = 0; q < P; q++) {
      int count[4] = {0, 0, 0, 0};
      for (int c : ch) {
        const std::uint8_t x = sets[static_cast<std::size_t>(c) * P + q];
        for (int k = 0; k < 4; k++) count[k] += (x >> k) & 1;
      }
      const int m = std::max(std::max(count[0], count[1]), std::max(count[2], count[3]));
      std::uint8_t out = 0;
      for (int k = 0; k < 4; k++)
        if (count[k] == m) out |= static_cast<std::uint8_t>(1u << k);
      s[q] = out;
      score += p.weight[q] * static_cast<double>(static_cast<int>(ch.size()) - m);
    }
  }
  return score;
}

}  // namespace detail

// The parsimony score of a tree: the least number of changes, weighted by the patterns.
inline double fitch_score(const Tree& t, const Patterns& p) {
  std::vector<std::uint8_t> sets;
  return detail::fitch(t, p, detail::rows_of_leaves(t, p), sets);
}

// The cost of a change from x to y, c[4x + y]; the unit cost (1 − δ) gives Fitch's score.
using CostMatrix = std::array<double, 16>;
inline CostMatrix unit_costs() {
  CostMatrix c{};
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++) c[static_cast<std::size_t>(x * 4 + y)] = x == y ? 0 : 1;
  return c;
}
// Transitions cost `transition`, transversions `transversion`.
inline CostMatrix transition_costs(double transition, double transversion) {
  CostMatrix c{};
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++) {
      const bool ts = (x ^ y) == 2;  // A↔G (0, 2) and C↔T (1, 3)
      c[static_cast<std::size_t>(x * 4 + y)] = x == y ? 0 : ts ? transition : transversion;
    }
  return c;
}

struct SankoffResult {
  double score = 0;
  std::vector<Code> states;  // nodes × patterns: one assignment attaining the score
  std::size_t patterns = 0;
  Code state(int node, std::size_t pattern) const { return states[static_cast<std::size_t>(node) * patterns + pattern]; }
};

inline SankoffResult sankoff(const Tree& t, const Patterns& p, const CostMatrix& cost) {
  const auto row = detail::rows_of_leaves(t, p);
  const std::size_t P = p.count;
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> S(t.size() * P * 4, 0);
  for (int v : t.postorder()) {
    double* s = &S[static_cast<std::size_t>(v) * P * 4];
    const auto& ch = t[v].children;
    if (ch.empty()) {
      for (std::size_t q = 0; q < P; q++) {
        const Code c = p.states[q * p.taxa + row[static_cast<std::size_t>(v)]];
        for (int x = 0; x < 4; x++) s[q * 4 + static_cast<std::size_t>(x)] = (c >= 4 || c == x) ? 0 : inf;
      }
      continue;
    }
    std::fill(s, s + P * 4, 0.0);
    for (int c : ch) {
      const double* sc = &S[static_cast<std::size_t>(c) * P * 4];
      for (std::size_t q = 0; q < P; q++)
        for (int x = 0; x < 4; x++) {
          double m = inf;
          for (int y = 0; y < 4; y++) m = std::min(m, cost[static_cast<std::size_t>(x * 4 + y)] + sc[q * 4 + static_cast<std::size_t>(y)]);
          s[q * 4 + static_cast<std::size_t>(x)] += m;
        }
    }
  }
  SankoffResult r;
  r.patterns = P;
  r.states.assign(t.size() * P, 0);
  const double* sr = &S[static_cast<std::size_t>(t.root) * P * 4];
  for (std::size_t q = 0; q < P; q++) {
    int best = 0;
    for (int x = 1; x < 4; x++)
      if (sr[q * 4 + static_cast<std::size_t>(x)] < sr[q * 4 + static_cast<std::size_t>(best)]) best = x;
    r.score += p.weight[q] * sr[q * 4 + static_cast<std::size_t>(best)];
    r.states[static_cast<std::size_t>(t.root) * P + q] = static_cast<Code>(best);
  }
  for (int v : t.preorder()) {
    if (v == t.root) continue;
    const double* sv = &S[static_cast<std::size_t>(v) * P * 4];
    for (std::size_t q = 0; q < P; q++) {
      const int x = r.states[static_cast<std::size_t>(t[v].parent) * P + q];
      int best = 0;
      double bv = inf;
      for (int y = 0; y < 4; y++) {
        const double c = cost[static_cast<std::size_t>(x * 4 + y)] + sv[q * 4 + static_cast<std::size_t>(y)];
        if (c < bv) { bv = c; best = y; }
      }
      r.states[static_cast<std::size_t>(v) * P + q] = static_cast<Code>(best);
    }
  }
  return r;
}

struct ParsimonyTree {
  Tree tree;
  double score = 0;
  std::size_t trees_scored = 0;
};

namespace detail {

// The tree with a new leaf on the branch above u.
inline Tree insert_leaf(const Tree& t, int u, const std::string& label) {
  Tree s = t;
  const int p = s[u].parent;
  const int w = static_cast<int>(s.nodes.size());
  Node mid;
  mid.parent = p;
  mid.children = {u};
  s.nodes.push_back(std::move(mid));
  for (int& c : s[p].children)
    if (c == u) c = w;
  s[u].parent = w;
  s.add(w, ABSENT, label);
  return s;
}

// Interchange across the branch above v: v's child at index `which` exchanged with v's first
// sibling.
inline bool interchange(Tree& t, int v, int which) {
  const int p = t[v].parent;
  if (p < 0 || t.is_leaf(v)) return false;
  int c = -1;
  for (int s : t[p].children)
    if (s != v) { c = s; break; }
  if (c < 0) return false;
  const int moved = t[v].children[static_cast<std::size_t>(which)];
  for (int& x : t[v].children)
    if (x == moved) x = c;
  for (int& x : t[p].children)
    if (x == c) x = moved;
  t[c].parent = v;
  t[moved].parent = p;
  return true;
}

inline Tree star3(const std::string& a, const std::string& b, const std::string& c) {
  Tree t;
  const int r = t.add(-1);
  t.add(r, ABSENT, a);
  t.add(r, ABSENT, b);
  t.add(r, ABSENT, c);
  return t;
}

}  // namespace detail

// Stepwise addition in the order of the rows, then nearest-neighbour interchanges while the
// score falls. The tree is unrooted and binary, without branch lengths.
inline ParsimonyTree parsimony_search(const Patterns& p, int rounds = 1000) {
  if (p.taxa < 3) throw std::invalid_argument("evolution::phylogeny::parsimony_search: at least three taxa");
  ParsimonyTree out;
  std::vector<std::uint8_t> sets;
  Tree t = detail::star3(p.names[0], p.names[1], p.names[2]);
  for (std::size_t k = 3; k < p.taxa; k++) {
    double best = std::numeric_limits<double>::infinity();
    Tree chosen;
    for (std::size_t u = 0; u < t.size(); u++) {
      if (static_cast<int>(u) == t.root) continue;
      Tree s = detail::insert_leaf(t, static_cast<int>(u), p.names[k]);
      const double sc = detail::fitch(s, p, detail::rows_of_leaves(s, p), sets);
      out.trees_scored++;
      if (sc < best) { best = sc; chosen = std::move(s); }
    }
    t = std::move(chosen);
  }
  double score = fitch_score(t, p);
  for (int r = 0; r < rounds; r++) {
    double best = score;
    Tree chosen;
    for (std::size_t v = 0; v < t.size(); v++) {
      if (static_cast<int>(v) == t.root || t.is_leaf(static_cast<int>(v))) continue;
      for (int which = 0; which < 2; which++) {
        Tree s = t;
        if (!detail::interchange(s, static_cast<int>(v), which)) continue;
        const double sc = detail::fitch(s, p, detail::rows_of_leaves(s, p), sets);
        out.trees_scored++;
        if (sc < best) { best = sc; chosen = std::move(s); }
      }
    }
    if (!(best < score)) break;
    t = std::move(chosen);
    score = best;
  }
  out.tree = compact(t);
  out.score = score;
  return out;
}

// The most parsimonious tree, exactly, by branch and bound (one of them if there are ties; the
// first found in the order of addition). For up to about a dozen taxa.
inline ParsimonyTree parsimony_exact(const Patterns& p) {
  if (p.taxa < 3) throw std::invalid_argument("evolution::phylogeny::parsimony_exact: at least three taxa");
  ParsimonyTree out = parsimony_search(p);
  double bound = out.score;
  std::vector<std::uint8_t> sets;
  // Depth-first over the positions of each taxon in turn.
  const auto search = [&](auto&& self, const Tree& t, std::size_t k) -> void {
    if (k == p.taxa) return;
    for (std::size_t u = 0; u < t.size(); u++) {
      if (static_cast<int>(u) == t.root) continue;
      Tree s = detail::insert_leaf(t, static_cast<int>(u), p.names[k]);
      const double sc = detail::fitch(s, p, detail::rows_of_leaves(s, p), sets);
      out.trees_scored++;
      // A partial tree that costs the best score already cannot lead to a better one.
      if (sc >= bound) continue;
      if (k + 1 == p.taxa) {
        bound = sc;
        out.tree = compact(s);
        out.score = sc;
        continue;
      }
      self(self, s, k + 1);
    }
  };
  search(search, detail::star3(p.names[0], p.names[1], p.names[2]), 3);
  return out;
}

}  // namespace evolution::phylogeny
