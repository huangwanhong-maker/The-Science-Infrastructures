// Sequences evolved along a tree, and random trees: the truths that tree building is checked
// against.
//
// - Sequences (as Seq-Gen does: Rambaut and Grassly 1997, "Seq-Gen: an application for the
//   Monte Carlo simulation of DNA sequence evolution along phylogenetic trees", CABIOS 13,
//   235–238): the root's sequence drawn from π; each site's rate drawn from the rate model
//   (a discrete Γ category or the invariant class, or a continuous Γ(α, α) rate); then along
//   each branch, each site's state drawn from the row of P(r t) for its parent's state, in
//   preorder. The leaves' sequences are the alignment, named by their labels (t<k> for an
//   unlabelled leaf k), in the order of the leaves. O(sites · nodes) after one 4 × 4
//   propagator per branch and category.
// - Random trees: Yule's pure-birth process (Yule 1925, Phil. Trans. R. Soc. B 213, 21–87):
//   from two lineages, each splits at rate λ, until n lineages, then a last waiting time
//   Exp(nλ) to the present; Kingman's coalescent (Kingman 1982, Stochastic Processes and their
//   Applications 13, 235–248): from n lineages, any two merge at rate 1, so with k lineages
//   the waiting time is Exp(k(k − 1)/2); and uniform unrooted topologies by inserting each new
//   leaf on a uniformly chosen branch, each of the (2n − 5)!! topologies equally likely
//   (Felsenstein 2004, Inferring Phylogenies, ch. 3), every branch's length independent, a
//   minimum plus an exponential.
//
// Every draw comes from a calx::random::Stream keyed by the seed, the purpose below and the
// object (the branch's node, the site), so a simulation is a function of its seed.

#pragma once

#include <evolution/phylogeny/alignment.hpp>
#include <evolution/phylogeny/models.hpp>
#include <evolution/phylogeny/tree.hpp>

#include <calx/random/stream.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace evolution::phylogeny {

inline constexpr std::uint32_t SIM_ROOT = 0x524f4f54, SIM_RATES = 0x52415445, SIM_BRANCH = 0x42524e43, SIM_TREE = 0x54524545;

struct SiteRates {
  Rates discrete;             // categories (used when continuous_alpha is 0)
  double continuous_alpha = 0;  // > 0: a Γ(α, α) rate for each site
  double invariant = 0;         // with a continuous Γ, the proportion of invariant sites
};

namespace detail {

inline Code draw_state(calx::random::Stream& r, const double* row) {
  const double u = r.uniform();
  double c = row[0];
  if (u < c) return 0;
  c += row[1];
  if (u < c) return 1;
  c += row[2];
  return u < c ? 2 : 3;
}

}  // namespace detail

// The leaves' sequences after `sites` sites evolve along t (branch lengths in substitutions per
// site; absent lengths count 0) under the model and the rates.
inline Alignment simulate_alignment(const Tree& t, const Model& model, const SiteRates& rates, std::size_t sites, std::uint64_t seed) {
  if (t.root < 0) throw std::invalid_argument("evolution::phylogeny::simulate_alignment: an empty tree");
  const bool continuous = rates.continuous_alpha > 0;
  const Rates& cats = rates.discrete;
  if (!continuous && cats.size() == 0) throw std::invalid_argument("evolution::phylogeny::simulate_alignment: no rate categories");
  // Each site's rate (and category).
  std::vector<double> rate(sites, 1);
  std::vector<std::uint32_t> cat(sites, 0);
  {
    calx::random::Stream r(seed, SIM_RATES, 0, 0);
    for (std::size_t s = 0; s < sites; s++) {
      if (continuous) {
        rate[s] = r.uniform() < rates.invariant ? 0.0 : r.gamma(rates.continuous_alpha) / rates.continuous_alpha / (1 - rates.invariant);
      } else {
        double u = r.uniform(), c = 0;
        std::uint32_t k = 0;
        for (; k + 1 < cats.size(); k++) {
          c += cats.weight[k];
          if (u < c) break;
        }
        cat[s] = k;
        rate[s] = cats.rate[k];
      }
    }
  }
  const std::size_t n = t.size();
  std::vector<Codes> seq(n);
  {
    calx::random::Stream r(seed, SIM_ROOT, 0, 0);
    seq[static_cast<std::size_t>(t.root)].resize(sites);
    for (auto& x : seq[static_cast<std::size_t>(t.root)]) x = detail::draw_state(r, model.pi.data());
  }
  std::vector<double> P(cats.size() * 16), Q(16);
  for (int v : t.preorder()) {
    if (v == t.root) continue;
    const double len = t.length_or(v, 0);
    if (!continuous)
      for (std::size_t k = 0; k < cats.size(); k++) model.transition(cats.rate[k] * len, &P[k * 16]);
    calx::random::Stream r(seed, SIM_BRANCH, static_cast<std::uint32_t>(v), 0);
    const Codes& parent = seq[static_cast<std::size_t>(t[v].parent)];
    Codes& x = seq[static_cast<std::size_t>(v)];
    x.resize(sites);
    for (std::size_t s = 0; s < sites; s++) {
      const double* M;
      if (continuous) {
        model.transition(rate[s] * len, Q.data());
        M = Q.data();
      } else {
        M = &P[cat[s] * 16];
      }
      x[s] = detail::draw_state(r, M + parent[s] * 4);
    }
    // A parent's sequence is no longer needed once its last child has it.
    const auto& sib = t[t[v].parent].children;
    if (sib.back() == v && t[v].parent != t.root) Codes().swap(seq[static_cast<std::size_t>(t[v].parent)]);
  }
  Alignment a;
  for (int v : t.leaves()) {
    a.names.push_back(t[v].label.empty() ? "t" + std::to_string(v) : t[v].label);
    a.rows.push_back(seq[static_cast<std::size_t>(v)]);
  }
  return a;
}

inline Alignment simulate_alignment(const Tree& t, const Model& model, const Rates& rates, std::size_t sites, std::uint64_t seed) {
  SiteRates r;
  r.discrete = rates;
  return simulate_alignment(t, model, r, sites, seed);
}

// A Yule tree on n leaves (t1 … tn), rooted and ultrametric, with birth rate λ.
inline Tree yule_tree(std::size_t n, std::uint64_t seed, double birth_rate = 1) {
  if (n < 2) throw std::invalid_argument("evolution::phylogeny::yule_tree: at least two leaves");
  calx::random::Stream r(seed, SIM_TREE, 1, 0);
  Tree t;
  const int root = t.add(-1);
  std::vector<int> tips{t.add(root, 0), t.add(root, 0)};
  std::vector<double> start(2 * n, 0);  // when each node's branch began
  double now = 0;
  while (tips.size() < n) {
    now += r.exponential() / (birth_rate * static_cast<double>(tips.size()));
    const std::size_t k = static_cast<std::size_t>(r.below(tips.size()));
    const int v = tips[k];
    t[v].length = now - start[static_cast<std::size_t>(v)];
    const int a = t.add(v, 0), b = t.add(v, 0);
    start[static_cast<std::size_t>(a)] = start[static_cast<std::size_t>(b)] = now;
    tips[k] = a;
    tips.push_back(b);
  }
  now += r.exponential() / (birth_rate * static_cast<double>(n));
  for (std::size_t k = 0; k < tips.size(); k++) {
    const int v = tips[k];
    t[v].length = now - start[static_cast<std::size_t>(v)];
    t[v].label = "t" + std::to_string(k + 1);
  }
  t.rooted = true;
  return compact(t);
}

// A Kingman coalescent tree on n leaves (t1 … tn), heights in units of `scale` (2N generations
// times the rate, for sequences).
inline Tree coalescent_tree(std::size_t n, std::uint64_t seed, double scale = 1) {
  if (n < 2) throw std::invalid_argument("evolution::phylogeny::coalescent_tree: at least two leaves");
  calx::random::Stream r(seed, SIM_TREE, 2, 0);
  // Build bottom-up into a pool, then hang from the last merger.
  struct Lineage { int node; double height; };
  std::vector<Node> pool(n);
  std::vector<double> height(n, 0);
  std::vector<Lineage> live;
  for (std::size_t i = 0; i < n; i++) {
    pool[i].label = "t" + std::to_string(i + 1);
    live.push_back({static_cast<int>(i), 0});
  }
  double now = 0;
  while (live.size() > 1) {
    const double k = static_cast<double>(live.size());
    now += r.exponential() / (k * (k - 1) / 2) * scale;
    const std::size_t i = static_cast<std::size_t>(r.below(live.size()));
    std::size_t j = static_cast<std::size_t>(r.below(live.size() - 1));
    if (j >= i) j++;
    Node m;
    m.children = {live[i].node, live[j].node};
    pool.push_back(m);
    height.push_back(now);
    const int id = static_cast<int>(pool.size()) - 1;
    for (int c : m.children) {
      pool[static_cast<std::size_t>(c)].parent = id;
      pool[static_cast<std::size_t>(c)].length = now - height[static_cast<std::size_t>(c)];
    }
    const std::size_t hi = std::max(i, j), lo = std::min(i, j);
    live.erase(live.begin() + static_cast<std::ptrdiff_t>(hi));
    live[lo] = {id, now};
  }
  Tree t;
  std::vector<std::pair<int, int>> stack{{live[0].node, t.add(-1)}};
  while (!stack.empty()) {
    const auto [src, dst] = stack.back();
    stack.pop_back();
    for (int c : pool[static_cast<std::size_t>(src)].children) {
      const int w = t.add(dst, pool[static_cast<std::size_t>(c)].length, pool[static_cast<std::size_t>(c)].label);
      stack.push_back({c, w});
    }
  }
  t.rooted = true;
  return compact(t);
}

// A uniformly random unrooted binary topology on n ≥ 3 leaves (t1 … tn), its branch lengths
// exponential with the given mean (plus `minimum`).
inline Tree random_unrooted_tree(std::size_t n, std::uint64_t seed, double mean_length = 0.1, double minimum = 0) {
  if (n < 3) throw std::invalid_argument("evolution::phylogeny::random_unrooted_tree: at least three leaves");
  calx::random::Stream r(seed, SIM_TREE, 3, 0);
  const auto length = [&] { return minimum + mean_length * r.exponential(); };
  Tree t;
  const int root = t.add(-1);
  for (int k = 1; k <= 3; k++) t.add(root, length(), "t" + std::to_string(k));
  for (std::size_t k = 4; k <= n; k++) {
    // A branch uniformly: the branch above a node other than the root.
    const int u = 1 + static_cast<int>(r.below(t.size() - 1));
    const int p = t[u].parent;
    const int w = static_cast<int>(t.nodes.size());
    Node mid;
    mid.parent = p;
    mid.children = {u};
    mid.length = length();
    t.nodes.push_back(std::move(mid));
    for (int& c : t[p].children)
      if (c == u) c = w;
    t[u].parent = w;
    t.add(w, length(), "t" + std::to_string(k));
  }
  t.rooted = false;
  return compact(t);
}

// Every branch length multiplied by a factor.
inline Tree scale_tree(Tree t, double factor) {
  for (auto& n : t.nodes)
    if (!std::isnan(n.length)) n.length *= factor;
  return t;
}

}  // namespace evolution::phylogeny
