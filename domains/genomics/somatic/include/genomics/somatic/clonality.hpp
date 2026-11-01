// Clonality: the cancer cell fraction of each somatic mutation, the clones its mutations cluster
// into, the trees of clones the sum rule allows, and the tumour mutational burden.
//
// - The cancer cell fraction (theory/the-somatic-genome.md, Proposition 2.2). A mutation present
//   in a fraction φ of tumour cells, at m copies in each, at a site of tumour copy number N = n_A +
//   n_B and purity ρ, is read by a fraction f = ρ φ m / (2(1 − ρ) + ρ N) of the reads. McGranahan et
//   al. 2015 (Sci. Transl. Med. 7, 283ra54, Supplementary Methods; the same relation in Dentro,
//   Wedge and Van Loo 2017, Cold Spring Harb. Perspect. Med. 7, a026625, Box 1): the mutation copy
//   number n_mut = (VAF/ρ)(ρ N + 2(1 − ρ)), the multiplicity m = max(1, round(n_mut)) held to the
//   major allele's copies, φ̂ = n_mut / m; and over φ ∈ {0.01, 0.02, …, 1} the posterior under a
//   uniform prior P(φ | a, d) ∝ Bin(a; d, ρ m φ / (2(1 − ρ) + ρ N)), its central 95% interval, the
//   mutation clonal when the interval reaches 1 (the paper's rule: "clonal if the upper band of the
//   95% confidence interval was ≥ 1").
// - Clones. The reads (a_i, d_i) of the mutations, with the scales c_i = ρ m_i / (2(1 − ρ) + ρ N_i),
//   are a mixture of binomials Bin(a_i; d_i, c_i φ_k) over clones k of CCF φ_k (DPClust's and
//   PyClone's model without copy-number uncertainty): tyche's EM fit for K = 1, …, K_max, K chosen
//   by BIC (tyche::multivariate::select_binomial_mixture).
// - Trees (clone_trees). Clones sorted by decreasing φ, the first the trunk (or, when its φ is
//   below 1 − tolerance, a trunk of φ = 1 without mutations of its own is put above it); each other
//   clone's parent is one of the clones before it, (K − 1)! trees, each kept when the sum rule
//   holds at every clone, Σ_children φ ≤ φ_parent + tolerance (El-Kebir et al. 2015,
//   Bioinformatics 31, i62–i70, the sum condition; Jiao et al. 2014, BMC Bioinformatics 15, 35).
//   The pigeonhole principle (Nik-Zainal et al. 2012, Cell 149, 994–1007) is its consequence for
//   two clones: φ_a + φ_b > 1 + tolerance forces one to descend from the other; the forced pairs
//   are reported. The tree chosen has the greatest least slack min_k (φ_k − Σ_children φ), the
//   one furthest from violating the rule. *Caveat.* One sample seldom determines the tree: a
//   chain and a branching often both satisfy the rule; the count of consistent trees says how
//   underdetermined it is, and several samples (regions, times) are what resolve it.
// - The mutational burden: passing calls per megabase of the callable territory (mutect.hpp),
//   substitutions and indels together and apart (Chalmers et al. 2017, Genome Medicine 9, 34,
//   counted over a territory, there an exome panel's).
//
// Complexity: O(n · 100) for the posteriors of n mutations; the mixture's EM; O((K − 1)! K) trees.

#pragma once

#include <genomics/somatic/copy_number.hpp>
#include <genomics/somatic/mutect.hpp>

#include <tyche/multivariate/binomial_mixture.hpp>

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace genomics::somatic {

struct MutationCcf {
  std::size_t call = 0;           // its index in the calls
  std::uint32_t contig = 0;
  std::size_t position = 0;
  std::size_t alt = 0, depth = 0;
  double vaf = 0;
  int total = 2, major = 1, minor = 1;
  double n_mut = 0;               // the mutation copy number
  int multiplicity = 1;
  double scale = 0;               // c = ρ m / (2(1 − ρ) + ρ N)
  double ccf = 0, ccf_low = 0, ccf_high = 0;
  bool clonal = false;
  int cluster = -1;
};

struct CcfOptions {
  double step = 0.01;            // the grid of φ
  double interval = 0.95;
  bool pass_only = true;
  bool indels = true;
};

namespace detail {

inline const CnSegment* segment_of(const std::vector<CnSegment>& s, std::uint32_t contig, std::size_t p) {
  auto it = std::upper_bound(s.begin(), s.end(), std::pair{contig, p}, [](const std::pair<std::uint32_t, std::size_t>& x, const CnSegment& g) {
    return x.first != g.contig ? x.first < g.contig : x.second < g.begin;
  });
  if (it == s.begin()) return nullptr;
  --it;
  return it->contig == contig && p < it->end ? &*it : nullptr;
}

}  // namespace detail

// The cancer cell fractions of the calls at the copy numbers of `segments` (ASCAT's integers) and
// purity ρ; calls outside every segment or at a total copy number of 0 are left out.
inline std::vector<MutationCcf> cancer_cell_fractions(const std::vector<SomaticCall>& calls, const std::vector<CnSegment>& segments, double purity, const CcfOptions& o = {}) {
  std::vector<MutationCcf> out;
  const double rho = purity;
  const std::size_t G = static_cast<std::size_t>(std::lround(1 / o.step));
  for (std::size_t i = 0; i < calls.size(); i++) {
    const SomaticCall& c = calls[i];
    if ((o.pass_only && !c.passed()) || (!o.indels && c.indel) || c.tumour_depth == 0) continue;
    const CnSegment* s = detail::segment_of(segments, c.contig, c.position);
    if (!s || s->n_a + s->n_b <= 0) continue;
    MutationCcf m;
    m.call = i;
    m.contig = c.contig;
    m.position = c.position;
    m.alt = c.tumour_alt;
    m.depth = c.tumour_depth;
    m.vaf = c.vaf();
    m.total = s->n_a + s->n_b;
    m.major = s->major();
    m.minor = s->minor();
    const double dna = 2 * (1 - rho) + rho * m.total;
    m.n_mut = m.vaf * dna / rho;
    m.multiplicity = std::clamp(static_cast<int>(std::lround(m.n_mut)), 1, std::max(1, m.major));
    m.scale = rho * m.multiplicity / dna;
    m.ccf = m.n_mut / m.multiplicity;
    // The posterior over the grid and its central interval.
    std::vector<double> lw(G);
    double top = -std::numeric_limits<double>::infinity();
    for (std::size_t g = 0; g < G; g++) {
      const double phi = static_cast<double>(g + 1) * o.step;
      const double p = std::min(1 - 1e-12, m.scale * phi);
      lw[g] = static_cast<double>(m.alt) * std::log(p) + static_cast<double>(m.depth - m.alt) * std::log1p(-p);
      top = std::max(top, lw[g]);
    }
    double total = 0;
    for (double& x : lw) { x = std::exp(x - top); total += x; }
    const double lo = (1 - o.interval) / 2, hi = 1 - lo;
    double cum = 0;
    bool low_set = false;
    m.ccf_high = 1;
    for (std::size_t g = 0; g < G; g++) {
      cum += lw[g] / total;
      if (!low_set && cum >= lo) { m.ccf_low = static_cast<double>(g + 1) * o.step; low_set = true; }
      if (cum >= hi) { m.ccf_high = static_cast<double>(g + 1) * o.step; break; }
    }
    m.clonal = m.ccf_high >= 1 - 1e-9;
    out.push_back(m);
  }
  return out;
}

struct Clone {
  double ccf = 0, weight = 0;
  std::size_t mutations = 0;
};

struct Clustering {
  std::vector<Clone> clones;           // by decreasing CCF
  std::vector<double> bic;             // for K = 1, 2, …
  std::size_t components = 0;
  double log_likelihood = 0;
};

// The clones of the mutations (their `cluster` set), by tyche's binomial mixture and BIC.
inline Clustering cluster_mutations(std::vector<MutationCcf>& muts, std::size_t max_clones = 6, tyche::multivariate::BinomialMixtureOptions mo = {}) {
  Clustering out;
  if (muts.empty()) return out;
  std::vector<double> k, n, c;
  for (const auto& m : muts) {
    k.push_back(static_cast<double>(m.alt));
    n.push_back(static_cast<double>(m.depth));
    c.push_back(m.scale);
  }
  mo.theta_max = 1;
  const auto sel = tyche::multivariate::select_binomial_mixture(k, n, c, max_clones, mo);
  for (const auto& f : sel.fits) out.bic.push_back(f.bic);
  const auto& f = sel.fits[sel.best];
  out.components = f.theta.size();
  out.log_likelihood = f.log_likelihood;
  out.clones.resize(f.theta.size());
  for (std::size_t j = 0; j < f.theta.size(); j++) out.clones[j] = {f.theta[j], f.weights[j], 0};
  for (std::size_t i = 0; i < muts.size(); i++) {
    muts[i].cluster = static_cast<int>(f.assignment[i]);
    out.clones[f.assignment[i]].mutations++;
  }
  return out;
}

struct CloneTree {
  std::vector<int> parent;   // over the nodes of TreeSearch (−1 for the root)
  double slack = 0;          // min_k (φ_k − Σ_children φ)
};

struct TreeSearch {
  std::vector<double> ccf;                     // the nodes' CCFs, by decreasing CCF; node 0 the trunk
  std::vector<int> clone;                      // each node's clone (−1: a trunk put above them)
  std::vector<CloneTree> trees;                // those satisfying the sum rule
  std::size_t examined = 0, best = 0;
  std::vector<std::pair<int, int>> forced;     // (ancestor, descendant) nodes by the pigeonhole principle
};

// The trees of clones (by CCF) that satisfy the sum rule within `tolerance`; at most 8 clones.
inline TreeSearch clone_trees(const std::vector<double>& clone_ccf, double tolerance = 0.05) {
  TreeSearch t;
  std::vector<int> order(clone_ccf.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return clone_ccf[static_cast<std::size_t>(a)] > clone_ccf[static_cast<std::size_t>(b)]; });
  if (order.empty() || clone_ccf[static_cast<std::size_t>(order[0])] < 1 - tolerance) {
    t.ccf.push_back(1);
    t.clone.push_back(-1);
  }
  for (int j : order) {
    t.ccf.push_back(clone_ccf[static_cast<std::size_t>(j)]);
    t.clone.push_back(j);
  }
  const std::size_t K = t.ccf.size();
  if (K > 9) throw std::invalid_argument("genomics::somatic::clone_trees: at most 8 clones and a trunk");
  for (std::size_t a = 1; a < K; a++)
    for (std::size_t b = a + 1; b < K; b++)
      if (t.ccf[a] + t.ccf[b] > 1 + tolerance) t.forced.push_back({static_cast<int>(a), static_cast<int>(b)});
  std::vector<int> parent(K, 0);
  parent[0] = -1;
  // Every parent vector with parent[k] ∈ [0, k), as a mixed-radix counter (parent[1] is always 0).
  for (;;) {
    t.examined++;
    std::vector<double> children(K, 0);
    for (std::size_t k = 1; k < K; k++) children[static_cast<std::size_t>(parent[k])] += t.ccf[k];
    double slack = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < K; k++) slack = std::min(slack, t.ccf[k] - children[k]);
    if (slack >= -tolerance) {
      t.trees.push_back({parent, slack});
      if (t.trees.back().slack > t.trees[t.best].slack) t.best = t.trees.size() - 1;
    }
    std::size_t k = K - 1;
    for (; k >= 2; k--) {
      if (++parent[k] < static_cast<int>(k)) break;
      parent[k] = 0;
    }
    if (K < 3 || k < 2) break;
  }
  return t;
}

// Whether node d descends from node a in a tree.
inline bool descends(const CloneTree& tree, int a, int d) {
  for (int x = d; x >= 0; x = tree.parent[static_cast<std::size_t>(x)])
    if (x == a) return true;
  return false;
}

struct Burden {
  std::size_t substitutions = 0, indels = 0, callable = 0;
  double per_megabase = 0, substitutions_per_megabase = 0, indels_per_megabase = 0;
};

inline Burden mutational_burden(const SomaticCalls& calls) {
  Burden b;
  b.callable = calls.callable;
  for (const auto& c : calls.calls)
    if (c.passed()) (c.indel ? b.indels : b.substitutions)++;
  if (b.callable) {
    const double mb = static_cast<double>(b.callable) / 1e6;
    b.per_megabase = static_cast<double>(b.substitutions + b.indels) / mb;
    b.substitutions_per_megabase = static_cast<double>(b.substitutions) / mb;
    b.indels_per_megabase = static_cast<double>(b.indels) / mb;
  }
  return b;
}

}  // namespace genomics::somatic
