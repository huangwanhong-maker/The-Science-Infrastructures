// Trees as systems of splits: the bipartitions of the leaves that the branches make; the
// Robinson–Foulds distance; compatibility; the tree of a compatible system; majority-rule
// consensus; and the support of a tree's splits in a sample of trees.
//
// Removing a branch of an unrooted tree on leaves X cuts X in two, A | X∖A. A tree is
// determined by its set of splits, and a set of splits is the set of a tree exactly when its
// splits are pairwise compatible (Buneman 1971, "The recovery of trees from measures of
// dissimilarity", in Mathematics in the Archaeological and Historical Sciences, 387–395;
// Semple and Steel 2003, Phylogenetics, §3.1, the splits-equivalence theorem): A | A′
// and B | B′ are compatible when one of A ∩ B, A ∩ B′, A′ ∩ B, A′ ∩ B′ is empty. A split is held
// here as the side that does not contain the first taxon, a bit set over the taxa in the order
// given, so that compatibility is containment or disjointness of the sides held, and the
// splits of a tree are a laminar family (the clusters of the tree rooted at the first taxon).
//
// - The Robinson–Foulds distance (Robinson and Foulds 1981, "Comparison of phylogenetic trees",
//   Mathematical Biosciences 53, 131–147) is the size of the symmetric difference of the two
//   sets of non-trivial splits; at most 2(n − 3) for binary trees on n leaves.
// - Majority-rule consensus (Margush and McMorris 1981, "Consensus n-trees", Bulletin of
//   Mathematical Biology 43, 239–244): the splits found in more than half the trees are
//   pairwise compatible (two incompatible splits cannot both be in a majority, since no tree
//   holds both), so they make a tree; each split's frequency is its support. With a threshold
//   below one half the splits are added greedily by frequency, each kept when compatible with
//   those already kept (the greedy, or extended majority-rule, consensus).
//
// Complexity: the splits of a tree O(n² / 64) words; Robinson–Foulds O(n² / 64 + n log n);
// the tree of k compatible splits O(k² n / 64); consensus of m trees O(m n² / 64 + m n log(mn)).
// Leaves are matched by label; a tree whose leaves are not the taxa given is an invalid
// argument.

#pragma once

#include <evolution/phylogeny/tree.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace evolution::phylogeny {

using Bits = std::vector<std::uint64_t>;

struct Split {
  Bits side;                // the side without the first taxon
  std::size_t size = 0;     // its number of taxa
  int node = -1;            // the node below the branch that makes it, in the tree it came from
  double length = ABSENT;   // the branch's length
  double support = ABSENT;  // the branch's support, or a frequency in a sample
};

namespace detail {

inline bool test(const Bits& b, std::size_t i) { return (b[i >> 6] >> (i & 63)) & 1; }
inline void set(Bits& b, std::size_t i) { b[i >> 6] |= std::uint64_t{1} << (i & 63); }
inline std::size_t popcount(const Bits& b) {
  std::size_t n = 0;
  for (auto w : b) n += static_cast<std::size_t>(std::popcount(w));
  return n;
}
inline bool subset(const Bits& a, const Bits& b) {
  for (std::size_t k = 0; k < a.size(); k++)
    if (a[k] & ~b[k]) return false;
  return true;
}
inline bool disjoint(const Bits& a, const Bits& b) {
  for (std::size_t k = 0; k < a.size(); k++)
    if (a[k] & b[k]) return false;
  return true;
}

// The index of each leaf's taxon, by label.
inline std::vector<std::size_t> taxon_of_leaves(const Tree& t, const std::vector<std::string>& taxa) {
  std::vector<std::pair<std::string, std::size_t>> sorted;
  for (std::size_t i = 0; i < taxa.size(); i++) sorted.emplace_back(taxa[i], i);
  std::sort(sorted.begin(), sorted.end());
  for (std::size_t i = 1; i < sorted.size(); i++)
    if (sorted[i].first == sorted[i - 1].first) throw std::invalid_argument("evolution::phylogeny: the taxon " + sorted[i].first + " is named twice");
  std::vector<std::size_t> out(t.size(), SIZE_MAX);
  std::vector<bool> seen(taxa.size(), false);
  std::size_t leaves = 0;
  for (std::size_t v = 0; v < t.size(); v++) {
    if (!t.nodes[v].children.empty()) continue;
    leaves++;
    const auto it = std::lower_bound(sorted.begin(), sorted.end(), std::make_pair(t.nodes[v].label, std::size_t{0}));
    if (it == sorted.end() || it->first != t.nodes[v].label)
      throw std::invalid_argument("evolution::phylogeny: the leaf " + t.nodes[v].label + " is not among the taxa");
    if (seen[it->second]) throw std::invalid_argument("evolution::phylogeny: two leaves are labelled " + t.nodes[v].label);
    seen[it->second] = true;
    out[v] = it->second;
  }
  if (leaves != taxa.size()) throw std::invalid_argument("evolution::phylogeny: the tree's leaves are not the taxa");
  return out;
}

}  // namespace detail

// The labels of a tree's leaves, sorted: the default order of taxa.
inline std::vector<std::string> taxa_of(const Tree& t) {
  auto l = t.leaf_labels();
  std::sort(l.begin(), l.end());
  return l;
}

// Whether two splits (sides without the first taxon) are compatible.
inline bool compatible(const Bits& a, const Bits& b) {
  return detail::subset(a, b) || detail::subset(b, a) || detail::disjoint(a, b);
}

// The splits of a tree over `taxa`, sorted by side; the trivial ones (a single leaf against
// the rest) only when asked. The two branches at a root of degree two make one split, whose
// length is the sum of theirs.
inline std::vector<Split> splits(const Tree& t, const std::vector<std::string>& taxa, bool trivial = false) {
  const std::size_t n = taxa.size(), words = (n + 63) / 64;
  const auto taxon = detail::taxon_of_leaves(t, taxa);
  std::vector<Bits> below(t.size(), Bits(words, 0));
  std::vector<Split> out;
  for (int v : t.postorder()) {
    auto& b = below[static_cast<std::size_t>(v)];
    if (t.is_leaf(v)) detail::set(b, taxon[static_cast<std::size_t>(v)]);
    for (int c : t[v].children)
      for (std::size_t k = 0; k < words; k++) b[k] |= below[static_cast<std::size_t>(c)][k];
    if (v == t.root) continue;
    Split s;
    s.side = b;
    if (n > 0 && detail::test(s.side, 0)) {
      for (std::size_t k = 0; k < words; k++) s.side[k] = ~s.side[k];
      if (n % 64) s.side[words - 1] &= (std::uint64_t{1} << (n % 64)) - 1;
    }
    s.size = detail::popcount(s.side);
    if (s.size == 0) continue;
    if (!trivial && (s.size < 2 || s.size + 2 > n)) continue;
    s.node = v;
    s.length = t[v].length;
    s.support = t[v].support;
    out.push_back(std::move(s));
  }
  std::sort(out.begin(), out.end(), [](const Split& a, const Split& b) { return a.side < b.side; });
  // One split from the two branches at a root of degree two.
  std::vector<Split> merged;
  for (auto& s : out) {
    if (!merged.empty() && merged.back().side == s.side) {
      merged.back().length = detail::add_lengths(merged.back().length, s.length);
      if (std::isnan(merged.back().support)) merged.back().support = s.support;
      continue;
    }
    merged.push_back(std::move(s));
  }
  return merged;
}

// The number of non-trivial splits in one tree and not the other.
inline std::size_t robinson_foulds(const Tree& a, const Tree& b) {
  const auto taxa = taxa_of(a);
  const auto sa = splits(a, taxa), sb = splits(b, taxa);
  std::size_t common = 0, i = 0, j = 0;
  while (i < sa.size() && j < sb.size()) {
    if (sa[i].side == sb[j].side) { common++; i++; j++; }
    else if (sa[i].side < sb[j].side) i++;
    else j++;
  }
  return sa.size() + sb.size() - 2 * common;
}

// The Robinson–Foulds distance over its largest value for binary trees, 2(n − 3).
inline double normalised_robinson_foulds(const Tree& a, const Tree& b) {
  const double n = static_cast<double>(a.leaf_count());
  return n > 3 ? static_cast<double>(robinson_foulds(a, b)) / (2 * (n - 3)) : 0.0;
}

// The tree whose non-trivial splits are the given ones, which must be pairwise compatible
// (otherwise std::invalid_argument). Splits carry their lengths and supports onto their
// branches; trivial splits in the list give the leaves' branches their lengths. The result is
// unrooted, hung from the node next to the first taxon.
inline Tree tree_from_splits(const std::vector<std::string>& taxa, std::vector<Split> given) {
  const std::size_t n = taxa.size(), words = (n + 63) / 64;
  if (n < 2) throw std::invalid_argument("evolution::phylogeny::tree_from_splits: at least two taxa");
  std::vector<double> leaf_length(n, ABSENT);
  std::vector<Split> s;
  for (auto& x : given) {
    if (x.side.size() != words) throw std::invalid_argument("evolution::phylogeny::tree_from_splits: a split over other taxa");
    if (detail::test(x.side, 0)) {
      for (std::size_t k = 0; k < words; k++) x.side[k] = ~x.side[k];
      if (n % 64) x.side[words - 1] &= (std::uint64_t{1} << (n % 64)) - 1;
    }
    x.size = detail::popcount(x.side);
    if (x.size == 1) {
      for (std::size_t i = 0; i < n; i++)
        if (detail::test(x.side, i)) leaf_length[i] = x.length;
    } else if (x.size + 1 == n) {
      leaf_length[0] = x.length;
    } else if (x.size >= 2) {
      s.push_back(std::move(x));
    }
  }
  std::stable_sort(s.begin(), s.end(), [](const Split& a, const Split& b) { return a.size > b.size; });
  for (std::size_t i = 0; i < s.size(); i++)
    for (std::size_t j = i + 1; j < s.size(); j++)
      if (!compatible(s[i].side, s[j].side)) throw std::invalid_argument("evolution::phylogeny::tree_from_splits: the splits are not compatible");
  Tree t;
  const int root = t.add(-1);
  t.add(root, leaf_length[0], taxa[0]);
  std::vector<int> node(s.size());
  for (std::size_t i = 0; i < s.size(); i++) {
    int parent = root;
    for (std::size_t j = i; j-- > 0;)
      if (detail::subset(s[i].side, s[j].side)) {
        if (s[j].size == s[i].size) throw std::invalid_argument("evolution::phylogeny::tree_from_splits: a split is given twice");
        parent = node[j];
        break;
      }
    node[i] = t.add(parent, s[i].length);
    t[node[i]].support = s[i].support;
  }
  for (std::size_t x = 1; x < n; x++) {
    int parent = root;
    for (std::size_t j = s.size(); j-- > 0;)
      if (detail::test(s[j].side, x)) { parent = node[j]; break; }
    t.add(parent, leaf_length[x], taxa[x]);
  }
  t.rooted = false;
  return compact(t);
}

// The consensus of trees on the same taxa: the splits in more than `threshold` of them (at
// least one half gives the majority-rule tree), each with its frequency as support and the
// mean of its lengths where the trees have them.
inline Tree majority_consensus(const std::vector<Tree>& trees, double threshold = 0.5) {
  if (trees.empty()) throw std::invalid_argument("evolution::phylogeny::majority_consensus: no trees");
  const auto taxa = taxa_of(trees[0]);
  std::vector<Split> all;
  for (const auto& t : trees) {
    auto s = splits(t, taxa, true);
    for (auto& x : s) all.push_back(std::move(x));
  }
  std::sort(all.begin(), all.end(), [](const Split& a, const Split& b) { return a.side < b.side; });
  const double m = static_cast<double>(trees.size());
  std::vector<Split> counted;
  for (std::size_t i = 0; i < all.size();) {
    std::size_t j = i;
    double len = 0, lengths = 0;
    while (j < all.size() && all[j].side == all[i].side) {
      if (!std::isnan(all[j].length)) { len += all[j].length; lengths++; }
      j++;
    }
    Split s = all[i];
    s.support = static_cast<double>(j - i) / m;
    s.length = lengths > 0 ? len / lengths : ABSENT;
    s.node = -1;
    counted.push_back(std::move(s));
    i = j;
  }
  // Greedy by frequency (ties by side), each kept when above the threshold and compatible.
  std::stable_sort(counted.begin(), counted.end(), [](const Split& a, const Split& b) { return a.support > b.support; });
  std::vector<Split> kept;
  for (auto& s : counted) {
    if (!(s.support > threshold)) continue;
    bool ok = true;
    if (s.size >= 2 && s.size + 1 < taxa.size())
      for (const auto& k : kept)
        if (k.size >= 2 && k.size + 1 < taxa.size() && !compatible(s.side, k.side)) { ok = false; break; }
    if (ok) kept.push_back(std::move(s));
  }
  // A leaf's branch keeps its mean length whatever its frequency (trivial splits are in every tree).
  return tree_from_splits(taxa, std::move(kept));
}

// For each node of `reference`, the fraction of `trees` that hold the split of the branch above
// it (NaN for leaves, the root, and the second branch at a root of degree two).
inline std::vector<double> split_frequencies(const Tree& reference, const std::vector<Tree>& trees) {
  const auto taxa = taxa_of(reference);
  const auto ref = splits(reference, taxa);
  std::vector<std::size_t> count(ref.size(), 0);
  for (const auto& t : trees) {
    const auto s = splits(t, taxa);
    std::size_t i = 0, j = 0;
    while (i < ref.size() && j < s.size()) {
      if (ref[i].side == s[j].side) { count[i]++; i++; j++; }
      else if (ref[i].side < s[j].side) i++;
      else j++;
    }
  }
  std::vector<double> out(reference.size(), ABSENT);
  for (std::size_t i = 0; i < ref.size(); i++)
    out[static_cast<std::size_t>(ref[i].node)] = trees.empty() ? ABSENT : static_cast<double>(count[i]) / static_cast<double>(trees.size());
  return out;
}

}  // namespace evolution::phylogeny
