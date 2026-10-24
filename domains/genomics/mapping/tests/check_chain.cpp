// Checks of chaining: the integer gap cost against the paper's real one; the DP of eqs. 1–2
// against brute force over every subset of a few seeds; the predecessor limit h; backtracking
// into disjoint chains, each increasing on both sequences and scored as f(end) − f(stop);
// primary and secondary chains by overlap on the query, fragments of one placement dropped;
// Li's mapping quality against its closed form; and seeds of the reverse strand placed on the
// read's reverse complement.

#include <genomics/index/minimizer_index.hpp>
#include <genomics/mapping/chain.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <vector>

using namespace genomics;
using mapping::Chain;
using mapping::ChainOptions;
using mapping::Seed;

namespace {

constexpr std::uint32_t CHECK = 0x4348414e;

// f(i) by brute force: the best chain ending at i over every subset of seeds before it, in the
// order of the sorted seeds, every step allowed.
std::vector<std::int64_t> brute_f(const std::vector<Seed>& seeds, const ChainOptions& o) {
  const std::size_t n = seeds.size();
  std::int64_t span = 0;
  for (const auto& s : seeds) span += s.span;
  span /= static_cast<std::int64_t>(n);
  std::vector<std::int64_t> f(n, std::numeric_limits<std::int64_t>::min());
  for (std::uint32_t mask = 1; mask < (1u << n); mask++) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < n; i++)
      if (mask >> i & 1) idx.push_back(i);
    std::int64_t v = seeds[idx[0]].span;
    bool ok = true;
    for (std::size_t t = 1; t < idx.size() && ok; t++) {
      const std::int64_t step = mapping::chain_step(seeds[idx[t - 1]], seeds[idx[t]], span, o);
      if (step == std::numeric_limits<std::int64_t>::min()) ok = false;
      else v += step;
    }
    if (ok) f[idx.back()] = std::max(f[idx.back()], v);
  }
  return f;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, CHECK, 0, 0);

  // ---- the gap cost -------------------------------------------------------------------------------------
  {
    // γ(l) = 0.01 w̄ |l| + 0.5 log₂|l| (Li 2018, eq. 2), against the integer form, which is never
    // above it and less than 1.5 below (⌊log₂⌋ loses < 1, halved; the floor < 1).
    double worst = 0;
    bool never_above = true;
    for (const std::int64_t w : {15, 21}) {
      for (std::int64_t l = 1; l <= 5000; l++) {
        const double real = 0.01 * static_cast<double>(w) * static_cast<double>(l) + 0.5 * std::log2(static_cast<double>(l));
        const double ours = static_cast<double>(mapping::gap_cost(l, w));
        never_above = never_above && ours <= real + 1e-12 && mapping::gap_cost(-l, w) == mapping::gap_cost(l, w);
        worst = std::max(worst, real - ours);
      }
    }
    s.require("the integer gap cost is symmetric and never above the paper's γ(l)", never_above);
    s.check("the integer gap cost is within 1.5 of the paper's γ(l) for |l| ≤ 5000", worst, 1.5);
    s.require("γ(0) = 0 and γ(10) = 3 at w̄ = 15 (1.5 + 0.5·3 rounded down)", mapping::gap_cost(0, 15) == 0 && mapping::gap_cost(10, 15) == 3,
              std::to_string(mapping::gap_cost(10, 15)));
  }

  // ---- the DP against brute force ------------------------------------------------------------------------
  {
    ChainOptions o;
    o.max_gap = 400;
    o.bandwidth = 60;
    o.max_predecessors = 100;
    bool same = true, limited = true;
    for (int trial = 0; trial < 300; trial++) {
      std::vector<Seed> seeds;
      const std::size_t n = 4 + r.below(9);
      for (std::size_t i = 0; i < n; i++) {
        const std::int64_t q = static_cast<std::int64_t>(r.below(600));
        const std::int64_t t = 1000 + q + static_cast<std::int64_t>(r.below(120)) - 60;
        seeds.push_back({0, false, t, q, 15});
      }
      if (trial % 3 == 0) seeds.push_back({0, true, 1200, 10, 15});  // another strand: never a predecessor
      std::sort(seeds.begin(), seeds.end());
      const auto dp = mapping::chain_dp(seeds, o);
      const auto f = brute_f(seeds, o);
      for (std::size_t i = 0; i < seeds.size(); i++) same = same && dp.f[i] == f[i];
      // With h = 1 only the next seed back may precede.
      ChainOptions one = o;
      one.max_predecessors = 1;
      const auto d1 = mapping::chain_dp(seeds, one);
      for (std::size_t i = 0; i < seeds.size(); i++) {
        std::int64_t expect = seeds[i].span;
        if (i > 0) {
          std::int64_t span = 0;
          for (const auto& x : seeds) span += x.span;
          span /= static_cast<std::int64_t>(seeds.size());
          const std::int64_t step = mapping::chain_step(seeds[i - 1], seeds[i], span, o);
          if (step != std::numeric_limits<std::int64_t>::min()) expect = std::max(expect, d1.f[i - 1] + step);
        }
        limited = limited && d1.f[i] == expect;
      }
    }
    s.require("f(i) of the DP equals the best chain over every subset of seeds (300 random sets of 4–12 seeds)", same);
    s.require("with h = 1 the DP tries only the nearest predecessor", limited);
  }

  // ---- backtracking --------------------------------------------------------------------------------------
  {
    ChainOptions o;
    o.min_anchors = 2;
    o.min_score = 20;
    bool disjoint = true, increasing = true, scored = true;
    for (int trial = 0; trial < 50; trial++) {
      std::vector<Seed> seeds;
      // Three diagonals and noise.
      for (int d = 0; d < 3; d++) {
        const std::int64_t diag = 5000 * d + static_cast<std::int64_t>(r.below(100));
        for (std::int64_t q = 0; q < 300; q += 7 + static_cast<std::int64_t>(r.below(10))) seeds.push_back({0, d == 2, diag + q, q, 15});
      }
      for (int x = 0; x < 20; x++) seeds.push_back({0, false, static_cast<std::int64_t>(r.below(20000)), static_cast<std::int64_t>(r.below(300)), 15});
      std::sort(seeds.begin(), seeds.end());
      seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());
      const auto chains = mapping::chain(seeds, o);
      std::set<std::tuple<bool, std::int64_t, std::int64_t>> seen;
      for (const auto& c : chains) {
        scored = scored && c.score >= o.min_score && c.seeds.size() >= o.min_anchors;
        for (std::size_t i = 0; i < c.seeds.size(); i++) {
          disjoint = disjoint && seen.insert({c.seeds[i].reverse, c.seeds[i].target, c.seeds[i].query}).second;
          if (i) increasing = increasing && c.seeds[i].target > c.seeds[i - 1].target && c.seeds[i].query > c.seeds[i - 1].query;
        }
      }
      for (std::size_t i = 1; i < chains.size(); i++) scored = scored && chains[i].score <= chains[i - 1].score;
    }
    s.require("backtracked chains share no seed", disjoint);
    s.require("every chain increases on both sequences", increasing);
    s.require("chains meet the minimum seeds and score, in decreasing score", scored);
  }

  // ---- primaries, secondaries, fragments, mapping quality ---------------------------------------------------
  {
    ChainOptions o;
    // Two copies of one placement 10 kb apart, the same seeds: equal chains, one primary, one secondary, MAPQ 0.
    std::vector<Seed> seeds;
    for (std::int64_t q = 0; q < 140; q += 6) {
      seeds.push_back({0, false, 20000 + q, q, 15});
      seeds.push_back({0, false, 30000 + q, q, 15});
    }
    auto chains = mapping::chain(seeds, o);
    mapping::select_chains(chains, o);
    s.require("two equal copies give one primary and one secondary chain", chains.size() == 2 && chains[0].primary && !chains[1].primary && chains[1].parent == 0);
    s.require("a read in two equal copies has mapping quality 0", !chains.empty() && chains[0].mapq == 0, std::to_string(chains.empty() ? -1 : chains[0].mapq));

    // One copy: high quality; Li's formula against its closed form.
    std::vector<Seed> one(seeds.begin(), seeds.end());
    one.erase(std::remove_if(one.begin(), one.end(), [](const Seed& x) { return x.target >= 30000; }), one.end());
    auto unique = mapping::chain(one, o);
    mapping::select_chains(unique, o);
    s.require("a unique placement is primary with mapping quality 60", unique.size() == 1 && unique[0].primary && unique[0].mapq == 60,
              unique.empty() ? "none" : std::to_string(unique[0].mapq));
    // 40 (1 − 80/100) min(1, 5/10) ln 100 = 18.42; 40 (1 − 40/150) ln 150 = 146.98, held to 60.
    s.require("mapQ = 40 (1 − f₂/f₁) min(1, m/10) ln f₁: 18 at (100, 80, 5), 60 at (150, 40, 20), 0 at f₂ = f₁",
              mapping::mapping_quality(100, 80, 5) == 18 && mapping::mapping_quality(150, 40, 20) == 60 && mapping::mapping_quality(120, 120, 30) == 0);

    // A fragment: seeds of the same placement left by the backtracking (a second, nearby diagonal).
    std::vector<Seed> frag = one;
    for (std::int64_t q = 3; q < 100; q += 6) frag.push_back({0, false, 20000 + q + 30, q, 15});
    auto fchains = mapping::chain(frag, o);
    const std::size_t before = fchains.size();
    mapping::select_chains(fchains, o);
    s.require("a chain on a diagonal within the band of a better one, overlapping it, is dropped as a fragment",
              before >= 2 && fchains.size() == 1 && fchains[0].mapq == 60, std::to_string(before) + " chains, " + std::to_string(fchains.size()) + " kept");

    // Two parts of a read on two places, not overlapping on the read: both primary (a supplementary part).
    std::vector<Seed> split;
    for (std::int64_t q = 0; q < 70; q += 6) split.push_back({0, false, 5000 + q, q, 15});
    for (std::int64_t q = 80; q < 150; q += 6) split.push_back({0, true, 90000 + q, q, 15});
    auto schains = mapping::chain(split, o);
    mapping::select_chains(schains, o);
    s.require("two parts of a read on two places are both primary", schains.size() == 2 && schains[0].primary && schains[1].primary);

    // A weak secondary (below 0.8 of the primary) is dropped but lowers the primary's quality.
    std::vector<Seed> weak = one;
    for (std::int64_t q = 0; q < 70; q += 6) weak.push_back({0, false, 50000 + q, q, 15});
    auto wchains = mapping::chain(weak, o);
    mapping::select_chains(wchains, o);
    const double f1 = static_cast<double>(wchains.empty() ? 0 : wchains[0].score);
    s.require("a secondary below 0.8 of the primary is dropped, its score kept as f₂", wchains.size() == 1 && wchains[0].sub_score > 0 && wchains[0].sub_score < 0.8 * f1);
  }

  // ---- seeds of the reverse strand -------------------------------------------------------------------------
  {
    calx::scheduler::ThreadPool pool(2);
    const sequence::Codes g = simulate::random_sequence(50000, 0.5, 7);
    index::MinimizerOptions mo;
    const index::MinimizerIndex idx(pool, {g}, mo);
    const sequence::Codes piece(g.begin() + 12000, g.begin() + 12600);
    const sequence::Codes read = sequence::reverse_complement(piece);
    const auto seeds = mapping::seeds_of(idx.anchors(read), read.size(), mo.k);
    const sequence::Codes rc = sequence::reverse_complement(read);
    bool exact = !seeds.empty();
    std::size_t reverse = 0;
    for (const auto& x : seeds) {
      if (!x.reverse) continue;
      reverse++;
      for (std::int32_t i = 0; i < x.span; i++) exact = exact && rc[static_cast<std::size_t>(x.query + i)] == g[static_cast<std::size_t>(x.target + i)];
    }
    s.require("a seed of the reverse strand is an exact forward match of the read's reverse complement", exact && reverse > 50, std::to_string(reverse) + " seeds");
    auto chains = mapping::chain(seeds, ChainOptions{});
    mapping::select_chains(chains);
    s.require("the reverse-complemented read chains to its place on the reverse strand",
              !chains.empty() && chains[0].reverse && chains[0].target_begin >= 12000 && chains[0].target_end <= 12600 && chains[0].mapq == 60);
  }

  return s.done();
}
