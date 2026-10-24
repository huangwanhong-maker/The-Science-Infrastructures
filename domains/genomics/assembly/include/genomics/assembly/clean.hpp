// Cleaning the compacted de Bruijn graph of the structures that sequencing errors and
// heterozygous sites make: tips, bubbles and low-coverage unitigs, each removal recorded with
// what was removed and why, and the graph compacted again after each step, so that an
// application can play the cleaning step by step.
//
// - A tip is a unitig with exactly one dead side (no neighbour there), shorter than
//   tip_length · k bases, and covered less than tip_ratio times the best competitor: a
//   unitig that leaves the same branching side as it does. An error in the last k bases of a
//   read makes one: at most k erroneous k-mers that leave the true path and stop (Zerbino and
//   Birney 2008, "Velvet: algorithms for de novo short read assembly using de Bruijn graphs",
//   Genome Research 18, 821–829, §"Tips", where tips shorter than 2k are clipped when their
//   multiplicity is below the alternative's). An island, a unitig shorter than tip_length · k
//   with no neighbour on either side, is removed in the same step: the errors that two reads
//   share, away from where they leave the genome, make them.
// - A bubble is two or more paths between the same two sides. An error inside a read, or a
//   heterozygous SNV, makes one whose arms are single unitigs of about 2k − 1 bases: k k-mers
//   that leave a branching side and rejoin. Velvet's Tour Bus (Zerbino and Birney 2008, §"Tour
//   bus") searches from each node, by Dijkstra's algorithm, for paths that rejoin, and merges
//   the less covered into the more covered when their sequences are alike. The simplified
//   form here pops bubbles whose arms are single unitigs: s → a → e where the arm a has s as
//   its only predecessor and e as its only successor, at most bubble_length · k bases long,
//   with lengths differing by at most bubble_indel among the arms; the best covered arm is
//   kept (ties to the lower number) and the others removed. Nested and multi-unitig bubbles
//   become simple ones as the rounds remove what is inside them, since each round compacts.
//   The sequences of the arms are not compared, and an arm is not merged into the one kept:
//   its k-mers are dropped, as the k-mers of a popped error should be. One structure that is
//   not a bubble of single unitigs is popped as well, the crossed bubble s → {a₁, a₂} →
//   {b₁, b₂} → e with every aᵢ joined to every bⱼ: two heterozygous sites exactly k bases
//   apart, which share no k-mer but one (k − 1)-mer; the better covered a and b are kept.
// - A low-coverage unitig is shorter than low_length · k bases and covered less than
//   low_ratio times its neighbours': the best coverage among the neighbours of each side,
//   the smaller of the two when both sides have neighbours. Such unitigs are the erroneous
//   connections that an error repeated by chance, or a chimeric read, puts between true paths.
//   *Caveat.* A short single-copy unitig between two copies of a high-copy repeat is removed
//   too; contiguity is already broken there.
//
// Each step finds its removals on the graph as it stands, all at once (so the result does not
// depend on an order of visit), removes them, and compacts (debruijn.hpp, compact). A round is
// tips, bubbles, low coverage; rounds repeat until one removes nothing, or `rounds` have run.
// Every step is O(n) in the unitigs of the graph, plus the compaction, O(total length).

#pragma once

#include <genomics/assembly/debruijn.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace genomics::assembly {

struct CleaningOptions {
  double tip_length = 2;          // tips shorter than tip_length · k bases
  double tip_ratio = 1;           // and covered less than tip_ratio · the best competitor
  double bubble_length = 3;       // arms at most bubble_length · k bases
  std::size_t bubble_indel = 10;  // whose lengths differ by at most this many bases
  double low_length = 2;          // low-coverage unitigs shorter than low_length · k bases
  double low_ratio = 0.2;         // and covered less than low_ratio · their neighbours'
  std::size_t rounds = 10;
  bool keep_graphs = false;       // keep the graph after every step
};

enum class StepKind { tips, bubbles, low_coverage };

inline const char* step_name(StepKind kind) {
  switch (kind) {
    case StepKind::tips: return "tips";
    case StepKind::bubbles: return "bubbles";
    case StepKind::low_coverage: return "low coverage";
  }
  return "";
}

// One unitig removed, numbered in the graph before its step, and why: the coverage it was
// compared with (the best competitor of a tip, the arm kept in a bubble, the neighbours of a
// low-coverage unitig).
struct Removal {
  std::uint32_t unitig = 0;
  std::size_t length = 0;
  double coverage = 0, reference = 0;
  std::uint32_t kept = NONE;      // bubbles: the arm kept
  std::uint32_t from = NONE, to = NONE;  // bubbles: the oriented unitigs the arms join
};

struct CleaningStep {
  StepKind kind = StepKind::tips;
  std::size_t round = 0;
  std::vector<Removal> removed;
  std::size_t unitigs_before = 0, unitigs_after = 0;
  // For each unitig after the step, the oriented unitigs before it that it is made of.
  std::vector<std::vector<std::uint32_t>> merged_from;
  Graph after;  // when CleaningOptions::keep_graphs
};

struct Cleaning {
  Graph graph;
  std::vector<CleaningStep> steps;
  std::size_t removed(StepKind kind) const {
    std::size_t n = 0;
    for (const auto& s : steps)
      if (s.kind == kind) n += s.removed.size();
    return n;
  }
};

inline std::vector<Removal> find_tips(const Graph& g, const CleaningOptions& o = {}) {
  std::vector<Removal> out;
  const double limit = o.tip_length * g.k;
  for (std::uint32_t u = 0; u < g.size(); u++) {
    const std::uint32_t plus = oriented(u, false);
    const bool in_dead = g.in_degree(plus) == 0, out_dead = g.out_degree(plus) == 0;
    if (static_cast<double>(g.length(u)) >= limit || (!in_dead && !out_dead)) continue;
    if (in_dead && out_dead) {
      // An island: short, and joined to nothing.
      out.push_back({u, g.length(u), g.coverage(u), 0});
      continue;
    }
    // Oriented so that its dead side is its end.
    const std::uint32_t t = out_dead ? plus : flip(plus);
    double best = 0;
    const Adjacency before = g.predecessors(t);  // held: a span of a temporary would dangle
    for (auto p : before.view())
      for (auto q : g.successors(p))
        if (element(q) != u) best = std::max(best, g.coverage(element(q)));
    if (best > 0 && g.coverage(u) < o.tip_ratio * best) out.push_back({u, g.length(u), g.coverage(u), best});
  }
  return out;
}

inline std::vector<Removal> find_bubbles(const Graph& g, const CleaningOptions& o = {}) {
  std::vector<Removal> out;
  const double limit = o.bubble_length * g.k;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> arms;  // (end, arm)
  for (std::uint32_t s = 0; s < 2 * g.size(); s++) {
    if (g.out_degree(s) < 2) continue;
    arms.clear();
    for (auto a : g.successors(s)) {
      if (element(a) == element(s) || g.in_degree(a) != 1 || g.out_degree(a) != 1) continue;
      if (static_cast<double>(g.length(element(a))) > limit) continue;
      const std::uint32_t e = g.successors(a)[0];
      if (element(e) == element(a) || element(e) == element(s)) continue;
      // Each bubble is seen from both of its sides; it is popped from the smaller.
      if (s > flip(e)) continue;
      arms.push_back({e, a});
    }
    std::sort(arms.begin(), arms.end());
    for (std::size_t i = 0; i < arms.size();) {
      std::size_t j = i;
      while (j < arms.size() && arms[j].first == arms[i].first) j++;
      if (j - i >= 2) {
        // The arm kept: the best covered, ties to the lower number.
        std::uint32_t keep = arms[i].second;
        for (std::size_t q = i + 1; q < j; q++) {
          const std::uint32_t a = arms[q].second;
          const double ca = g.coverage(element(a)), ck = g.coverage(element(keep));
          if (ca > ck || (ca == ck && element(a) < element(keep))) keep = a;
        }
        const std::size_t lk = g.length(element(keep));
        for (std::size_t q = i; q < j; q++) {
          const std::uint32_t a = arms[q].second;
          if (a == keep) continue;
          const std::size_t la = g.length(element(a));
          if ((la > lk ? la - lk : lk - la) > o.bubble_indel) continue;
          out.push_back({element(a), la, g.coverage(element(a)), g.coverage(element(keep)), element(keep), s, arms[q].first});
        }
      }
      i = j;
    }
    // A crossed bubble: s → a₁, a₂, each a → b₁, b₂, each b → e. Two heterozygous sites exactly k
    // apart make it, the k-mers of their alleles meeting in one (k − 1)-mer. The better covered
    // of each pair is kept.
    if (g.out_degree(s) != 2) continue;
    const std::uint32_t a1 = g.successors(s)[0], a2 = g.successors(s)[1];
    const auto pair_of = [&](std::uint32_t a) {
      std::array<std::uint32_t, 2> b{g.next[a].to[0], g.next[a].to[1]};
      std::sort(b.begin(), b.end());
      return b;
    };
    if (g.in_degree(a1) != 1 || g.in_degree(a2) != 1 || g.out_degree(a1) != 2 || g.out_degree(a2) != 2 || pair_of(a1) != pair_of(a2)) continue;
    const auto [b1, b2] = pair_of(a1);
    if (g.in_degree(b1) != 2 || g.in_degree(b2) != 2 || g.out_degree(b1) != 1 || g.out_degree(b2) != 1 || g.successors(b1)[0] != g.successors(b2)[0]) continue;
    const std::uint32_t e = g.successors(b1)[0];
    const std::array<std::uint32_t, 6> parts{element(s), element(a1), element(a2), element(b1), element(b2), element(e)};
    bool distinct = g.in_degree(e) == 2 && s <= flip(e);
    for (std::size_t x = 0; x < 6; x++)
      for (std::size_t y = x + 1; y < 6; y++) distinct = distinct && parts[x] != parts[y];
    for (auto v : {a1, a2, b1, b2}) distinct = distinct && static_cast<double>(g.length(element(v))) <= limit;
    if (!distinct) continue;
    const auto weaker = [&](std::uint32_t x, std::uint32_t y) {
      const double cx = g.coverage(element(x)), cy = g.coverage(element(y));
      return cx < cy || (cx == cy && element(x) > element(y)) ? std::pair(x, y) : std::pair(y, x);
    };
    for (const auto& [drop, keep] : {weaker(a1, a2), weaker(b1, b2)})
      out.push_back({element(drop), g.length(element(drop)), g.coverage(element(drop)), g.coverage(element(keep)), element(keep), s, e});
  }
  return out;
}

inline std::vector<Removal> find_low_coverage(const Graph& g, const CleaningOptions& o = {}) {
  std::vector<Removal> out;
  const double limit = o.low_length * g.k;
  for (std::uint32_t u = 0; u < g.size(); u++) {
    if (static_cast<double>(g.length(u)) >= limit) continue;
    const std::uint32_t plus = oriented(u, false);
    double left = 0, right = 0;
    bool has_left = false, has_right = false;
    const Adjacency before = g.predecessors(plus);
    for (auto p : before.view())
      if (element(p) != u) { has_left = true; left = std::max(left, g.coverage(element(p))); }
    for (auto q : g.successors(plus))
      if (element(q) != u) { has_right = true; right = std::max(right, g.coverage(element(q))); }
    if (!has_left && !has_right) continue;
    const double reference = has_left && has_right ? std::min(left, right) : std::max(left, right);
    if (g.coverage(u) < o.low_ratio * reference) out.push_back({u, g.length(u), g.coverage(u), reference});
  }
  return out;
}

// The graph cleaned of tips, bubbles and low-coverage unitigs, with every step recorded.
inline Cleaning clean(const Graph& g, const CleaningOptions& o = {}) {
  Cleaning c;
  c.graph = g;
  for (std::size_t round = 0; round < o.rounds; round++) {
    bool any = false;
    for (StepKind kind : {StepKind::tips, StepKind::bubbles, StepKind::low_coverage}) {
      std::vector<Removal> removed = kind == StepKind::tips ? find_tips(c.graph, o) : kind == StepKind::bubbles ? find_bubbles(c.graph, o) : find_low_coverage(c.graph, o);
      if (removed.empty()) continue;
      any = true;
      std::vector<char> flags(c.graph.size(), 0);
      for (const auto& r : removed) flags[r.unitig] = 1;
      CleaningStep step;
      step.kind = kind;
      step.round = round;
      step.unitigs_before = c.graph.size();
      Graph next = compact(c.graph, flags, &step.merged_from);
      step.unitigs_after = next.size();
      step.removed = std::move(removed);
      if (o.keep_graphs) step.after = next;
      c.graph = std::move(next);
      c.steps.push_back(std::move(step));
    }
    if (!any) break;
  }
  return c;
}

}  // namespace genomics::assembly
