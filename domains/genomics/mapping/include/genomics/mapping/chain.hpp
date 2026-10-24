// Chaining of seed anchors by dynamic programming with minimap2's gap cost, the backtracking of
// chains, their division into primary and secondary chains by overlap on the query, and the
// mapping quality of a primary chain (Li 2018, "Minimap2: pairwise alignment for nucleotide
// sequences", Bioinformatics 34, 3094–3100, §2.1).
//
// - Seeds. An anchor of the minimizer index (index/minimizer_index.hpp) is a k-mer shared by the
//   read and the genome. On the reverse strand it is turned into the coordinates of the read's
//   reverse complement, so that every seed is a forward match of `span` bases from query base
//   `query` (of the strand mapped) to target base `target`: query' = L − k − query for a read of
//   length L (an anchor on the reverse strand pairs query base j + i with target base
//   t + k − 1 − i, which is base L − k − j + (k − 1 − i) of the reverse complement). Seeds of
//   other spans (the SMEMs of mapper.hpp) are chained alike.
// - The chaining score (Li 2018, §2.1.2, eqs. 1–2). With the seeds of one strand and contig
//   sorted by the ends (target + span, query + span), write x_i for the end of seed i on the
//   target, y_i for its end on the query, and w_i for its span. Then
//
//     f(i) = max{ max_{j < i} [ f(j) + α(j, i) − β(j, i) ], w_i },
//     α(j, i) = min{ min{y_i − y_j, x_i − x_j}, w_i },
//     β(j, i) = ∞ when y_i ≤ y_j, x_i ≤ x_j, max{y_i − y_j, x_i − x_j} > G, or |l| > bandwidth;
//               γ(l) otherwise, with l = (y_i − y_j) − (x_i − x_j),
//     γ(l) = 0.01 · w̄ · |l| + 0.5 · log₂|l|  (l ≠ 0),  γ(0) = 0,
//
//   α being the bases the seed i adds to the chain, β the cost of the gap between them and w̄ the
//   mean span of the seeds. The coordinates are the seeds' ends, as minimap2 keeps an anchor by
//   its last base: for seeds of one span k the differences are those of their starts, but for
//   seeds of different spans only the ends make α the bases of i beyond j (with starts, a long
//   seed after a short one would be counted as overlapping it). Only the h nearest predecessors
//   of i on its strand and contig are tried, and the search back stops when x_i − x_j > G
//   (minimap2's -g; the predecessors are sorted by their ends on the target): O(N h) for N seeds.
//   The bandwidth (minimap2's -r) is a cut-off of minimap2's implementation, not of the paper's
//   equations.
//   *Change.* Scores are integers, exact and the same on every machine: γ(l) is computed as
//   ⌊(2 w̄ |l| + 100 ⌊log₂|l|⌋) / 200⌋, that is the paper's cost with log₂|l| replaced by the
//   bit length of |l| minus one and the sum rounded down, as minimap2's implementation rounds it
//   to an integer; w̄ is the mean span, rounded down (k for minimizers).
// - Backtracking (§2.1.3). Each seed keeps the predecessor that attained f. Seeds are visited
//   in decreasing f (ties by index); from each seed not yet in a chain the predecessors are
//   followed until a seed already in a chain, or none; the chain's score is f(i) minus f of the
//   seed where the walk stopped. Chains of fewer than `min_anchors` seeds or of a score below
//   `min_score` are dropped, their seeds staying used, as minimap2 does.
// - Primary and secondary chains (§2.1.4). In decreasing score, a chain is first dropped as a
//   fragment when an earlier chain kept lies on the same contig and strand, overlaps it on the
//   query, and has a diagonal (target − query at their first seeds) within the bandwidth of its
//   own: the backtracking leaves such chains of seeds not taken by a better chain of the same
//   placement, and they are not another place the read could come from. *Change:* the paper
//   states no such rule; without it a fragment would lower the mapping quality of a unique
//   read. Then a chain whose query interval overlaps that of an earlier primary chain by at
//   least mask_level (0.5) of the shorter of the two is secondary to the first such primary;
//   otherwise it is primary itself (a second primary is a supplementary part of the read, such
//   as either side of a chimera). A primary
//   keeps the best score of its secondaries, f₂, kept or not. Secondaries are kept when their
//   score is at least `secondary_ratio` (minimap2's -p, 0.8) of their primary's, at most
//   `best_n` of them (-N).
// - Mapping quality (§2.1.4): mapQ = 40 · (1 − f₂/f₁) · min{1, m/10} · ln f₁, with f₁ the
//   primary's score, m its number of seeds, and f₂ at least min_score (as minimap2's
//   implementation bounds the secondary score below); rounded down and held to [0, 60].
//   Secondaries have mapping quality 0. mapper.hpp refines f₂/f₁ once the chains are aligned.
//
// Complexity: O(N log N) to sort, O(N h) to chain, O(N log N) to backtrack, O(c²) for c chains
// to divide. Deterministic: every sort has a total order.

#pragma once

#include <genomics/index/minimizer_index.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace genomics::mapping {

// A forward match of `span` bases: query bases [query, query + span) of the strand mapped against
// target bases [target, target + span). Ordered by strand, contig and their ends.
struct Seed {
  std::uint32_t contig = 0;
  bool reverse = false;
  std::int64_t target = 0, query = 0;
  std::int32_t span = 0;
  std::int64_t target_end() const { return target + span; }
  std::int64_t query_end() const { return query + span; }
  friend bool operator==(const Seed&, const Seed&) = default;
  friend bool operator<(const Seed& a, const Seed& b) {
    if (a.reverse != b.reverse) return a.reverse < b.reverse;
    if (a.contig != b.contig) return a.contig < b.contig;
    if (a.target_end() != b.target_end()) return a.target_end() < b.target_end();
    if (a.query_end() != b.query_end()) return a.query_end() < b.query_end();
    return a.span < b.span;
  }
};

struct ChainOptions {
  std::int64_t max_gap = 5000;       // G: no predecessor further than this on either sequence (-g)
  std::int64_t bandwidth = 500;      // |l| at most this (-r)
  int max_predecessors = 50;         // h
  std::size_t min_anchors = 3;       // seeds a chain needs (-n)
  std::int64_t min_score = 40;       // score a chain needs (-m)
  double mask_level = 0.5;           // overlap on the query that makes a chain secondary
  double secondary_ratio = 0.8;      // secondaries kept down to this fraction of the primary's score (-p)
  std::size_t best_n = 5;            // secondaries kept at most (-N)
};

struct Chain {
  std::uint32_t contig = 0;
  bool reverse = false;
  std::int64_t score = 0;
  std::vector<Seed> seeds;  // in increasing target and query
  std::int64_t query_begin = 0, query_end = 0, target_begin = 0, target_end = 0;  // what the seeds span
  // From select_chains:
  bool primary = false;
  int parent = -1;              // the index of its primary chain (itself for a primary)
  std::int64_t sub_score = 0;   // f₂: the best score of its secondaries (primaries only)
  int mapq = 0;
};

// The seeds of a read of `read_length` bases from its anchors, with k the index's k.
inline std::vector<Seed> seeds_of(const std::vector<index::Anchor>& anchors, std::size_t read_length, int k) {
  std::vector<Seed> out;
  out.reserve(anchors.size());
  const auto L = static_cast<std::int64_t>(read_length);
  for (const auto& a : anchors) {
    const std::int64_t q = a.reverse ? L - k - static_cast<std::int64_t>(a.query) : static_cast<std::int64_t>(a.query);
    out.push_back({a.contig, a.reverse, static_cast<std::int64_t>(a.target), q, k});
  }
  std::sort(out.begin(), out.end());
  return out;
}

// γ(l) in integers: ⌊(2 w̄ |l| + 100 ⌊log₂|l|⌋) / 200⌋.
inline std::int64_t gap_cost(std::int64_t l, std::int64_t mean_span) {
  if (l == 0) return 0;
  const std::uint64_t a = static_cast<std::uint64_t>(l < 0 ? -l : l);
  const auto log2 = static_cast<std::int64_t>(std::bit_width(a)) - 1;
  return (2 * mean_span * static_cast<std::int64_t>(a) + 100 * log2) / 200;
}

// The DP of eqs. 1–2 over seeds sorted as Seed orders them: f(i) and the predecessor attaining it
// (−1 when f(i) = w_i).
struct ChainDp {
  std::vector<std::int64_t> f;
  std::vector<std::int32_t> predecessor;
};

// α − β of the step from seed j to seed i, or nullopt-like INT64_MIN when the step is not allowed.
inline std::int64_t chain_step(const Seed& j, const Seed& i, std::int64_t mean_span, const ChainOptions& o) {
  if (j.reverse != i.reverse || j.contig != i.contig) return std::numeric_limits<std::int64_t>::min();
  const std::int64_t dq = i.query_end() - j.query_end(), dt = i.target_end() - j.target_end();
  if (dq <= 0 || dt <= 0 || dq > o.max_gap || dt > o.max_gap) return std::numeric_limits<std::int64_t>::min();
  const std::int64_t l = dq - dt;
  if ((l < 0 ? -l : l) > o.bandwidth) return std::numeric_limits<std::int64_t>::min();
  const std::int64_t alpha = std::min<std::int64_t>(std::min(dq, dt), i.span);
  return alpha - gap_cost(l, mean_span);
}

inline ChainDp chain_dp(const std::vector<Seed>& seeds, const ChainOptions& o) {
  ChainDp dp;
  const std::size_t n = seeds.size();
  dp.f.resize(n);
  dp.predecessor.assign(n, -1);
  if (n == 0) return dp;
  std::int64_t span_sum = 0;
  for (const auto& s : seeds) span_sum += s.span;
  const std::int64_t mean_span = span_sum / static_cast<std::int64_t>(n);
  const auto h = static_cast<std::size_t>(std::max(0, o.max_predecessors));
  for (std::size_t i = 0; i < n; i++) {
    const Seed& si = seeds[i];
    std::int64_t best = si.span;
    std::int32_t from = -1;
    const std::size_t stop = i > h ? i - h : 0;
    for (std::size_t j = i; j-- > stop;) {
      const Seed& sj = seeds[j];
      if (sj.reverse != si.reverse || sj.contig != si.contig) break;
      if (si.target_end() - sj.target_end() > o.max_gap) break;
      const std::int64_t step = chain_step(sj, si, mean_span, o);
      if (step == std::numeric_limits<std::int64_t>::min()) continue;
      const std::int64_t v = dp.f[j] + step;
      if (v > best) { best = v; from = static_cast<std::int32_t>(j); }
    }
    dp.f[i] = best;
    dp.predecessor[i] = from;
  }
  return dp;
}

// The chains of seeds (sorted inside), in decreasing score.
inline std::vector<Chain> chain(std::vector<Seed> seeds, const ChainOptions& o = {}) {
  std::sort(seeds.begin(), seeds.end());
  const ChainDp dp = chain_dp(seeds, o);
  const std::size_t n = seeds.size();
  std::vector<std::uint32_t> order(n);
  std::iota(order.begin(), order.end(), 0u);
  std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) { return dp.f[a] != dp.f[b] ? dp.f[a] > dp.f[b] : a < b; });
  std::vector<char> used(n, 0);
  std::vector<Chain> chains;
  std::vector<std::uint32_t> path;
  for (const std::uint32_t i : order) {
    if (used[i]) continue;
    path.clear();
    std::int64_t k = i;
    while (k >= 0 && !used[static_cast<std::size_t>(k)]) {
      path.push_back(static_cast<std::uint32_t>(k));
      used[static_cast<std::size_t>(k)] = 1;
      k = dp.predecessor[static_cast<std::size_t>(k)];
    }
    const std::int64_t score = dp.f[i] - (k >= 0 ? dp.f[static_cast<std::size_t>(k)] : 0);
    if (path.size() < o.min_anchors || score < o.min_score) continue;
    Chain c;
    c.contig = seeds[i].contig;
    c.reverse = seeds[i].reverse;
    c.score = score;
    c.seeds.reserve(path.size());
    for (std::size_t p = path.size(); p-- > 0;) c.seeds.push_back(seeds[path[p]]);
    c.query_begin = c.seeds.front().query;
    c.target_begin = c.seeds.front().target;
    c.query_end = c.target_end = 0;
    for (const auto& s : c.seeds) {
      c.query_begin = std::min(c.query_begin, s.query);
      c.query_end = std::max(c.query_end, s.query + s.span);
      c.target_end = std::max(c.target_end, s.target + s.span);
    }
    chains.push_back(std::move(c));
  }
  std::stable_sort(chains.begin(), chains.end(), [](const Chain& a, const Chain& b) { return a.score > b.score; });
  return chains;
}

// Li's mapping quality from the primary's score f₁, its secondary score f₂ and its seed count m.
inline int mapping_quality(double f1, double f2, std::size_t m) {
  if (!(f1 > 1)) return 0;
  const double x = std::clamp(f2 / f1, 0.0, 1.0);
  const double q = 40.0 * (1 - x) * std::min(1.0, static_cast<double>(m) / 10.0) * std::log(f1);
  return std::clamp(static_cast<int>(std::floor(q)), 0, 60);
}

// The length of the overlap of [a, b) and [c, d).
inline std::int64_t overlap(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t d) { return std::max<std::int64_t>(0, std::min(b, d) - std::max(a, c)); }

// Divides chains (in decreasing score, as chain() returns them) into primaries and secondaries,
// sets each primary's f₂ and mapping quality, and drops the fragments and the secondaries not
// kept. The first chain
// left is the best primary; the others follow in decreasing score.
inline void select_chains(std::vector<Chain>& chains, const ChainOptions& o = {}) {
  std::vector<std::size_t> primaries, earlier;
  std::vector<char> fragment(chains.size(), 0);
  for (std::size_t c = 0; c < chains.size(); c++) {
    Chain& x = chains[c];
    x.primary = true;
    x.parent = static_cast<int>(c);
    const std::int64_t dx = x.seeds.front().target - x.seeds.front().query;
    for (const std::size_t p : earlier) {
      const Chain& y = chains[p];
      const std::int64_t dy = y.seeds.front().target - y.seeds.front().query;
      if (y.contig == x.contig && y.reverse == x.reverse && (dx > dy ? dx - dy : dy - dx) <= o.bandwidth &&
          overlap(x.query_begin, x.query_end, y.query_begin, y.query_end) > 0) {
        fragment[c] = 1;
        break;
      }
    }
    if (fragment[c]) continue;
    earlier.push_back(c);
    for (const std::size_t p : primaries) {
      Chain& y = chains[p];
      const std::int64_t ov = overlap(x.query_begin, x.query_end, y.query_begin, y.query_end);
      const std::int64_t shorter = std::min(x.query_end - x.query_begin, y.query_end - y.query_begin);
      if (static_cast<double>(ov) >= o.mask_level * static_cast<double>(shorter)) {
        x.primary = false;
        x.parent = static_cast<int>(p);
        y.sub_score = std::max(y.sub_score, x.score);
        break;
      }
    }
    if (x.primary) primaries.push_back(c);
  }
  for (const std::size_t p : primaries) {
    Chain& y = chains[p];
    y.mapq = mapping_quality(static_cast<double>(y.score), static_cast<double>(std::max(y.sub_score, o.min_score)), y.seeds.size());
  }
  // Keep the primaries and the best secondaries; parents renumbered.
  std::vector<int> renumber(chains.size(), -1);
  std::vector<Chain> kept;
  std::size_t secondaries = 0;
  for (std::size_t c = 0; c < chains.size(); c++) {
    Chain& x = chains[c];
    if (fragment[c]) continue;
    if (!x.primary) {
      const Chain& p = chains[static_cast<std::size_t>(x.parent)];
      if (static_cast<double>(x.score) < o.secondary_ratio * static_cast<double>(p.score) || secondaries >= o.best_n) continue;
      secondaries++;
      x.mapq = 0;
    }
    renumber[c] = static_cast<int>(kept.size());
    kept.push_back(std::move(x));
  }
  for (auto& x : kept) x.parent = renumber[static_cast<std::size_t>(x.parent)];
  chains = std::move(kept);
}

}  // namespace genomics::mapping
