// The guide tree of a progressive alignment: distances between sequences, without aligning
// them (k-mers) or from an alignment (corrected identities), and their clustering by UPGMA into
// a dendrogram; and the sequence weights of CLUSTAL W read from it.
//
// - The k-mer distance of MUSCLE's first stage (Edgar 2004, Nucleic Acids Research 32,
//   1792–1797, §Methods, "Stage 1"; Edgar 2004, BMC Bioinformatics 5, 113, §k-mer distance):
//   with n_x(τ) the number of occurrences of the k-mer τ in x,
//     F(x, y) = Σ_τ min(n_x(τ), n_y(τ)) / (min(|x|, |y|) − k + 1),   d = 1 − F,
//   the fraction of the shorter sequence's k-mers that the other shares. The letters may first be
//   compressed into classes (Edgar 2004, Nucleic Acids Research 32, 380–385): for protein,
//   Dayhoff's six groups AGPST, C, DENQ, HKR, ILMV, FWY (B and Z join DENQ). k-mers holding an
//   ambiguous letter (N; X, *) are skipped, and the denominator counts the k-mers kept. Each
//   sequence's k-mers are sorted once, O(L log L); a pair is a merge of two sorted lists, O(L);
//   all pairs O(n² L), in parallel by rows.
// - Corrected distances from an alignment (MUSCLE's second stage): D the fraction of differing
//   residues among the columns where both rows have an unambiguous one; for protein Kimura's
//   d = −ln(1 − D − D²/5) (Kimura 1983, The Neutral Theory of Molecular Evolution, Cambridge
//   University Press, p. 75), for DNA Jukes and Cantor's d = −(3/4) ln(1 − 4D/3) (Jukes
//   and Cantor 1969, in Munro (ed.), Mammalian Protein Metabolism III, 21–132). Both diverge as
//   the argument of the logarithm reaches 0; it is held at 0.05 or above (d ≤ 3.0 and 2.25),
//   where the corrections no longer measure anything (MUSCLE uses a table there instead).
// - UPGMA (Sokal and Michener 1958, University of Kansas Science Bulletin 38, 1409–1438): merge
//   the two nearest clusters, the distance of a cluster to the union being the size-weighted
//   mean d(k, a ∪ b) = (|a| d(k, a) + |b| d(k, b)) / (|a| + |b|) (the Lance–Williams update of
//   average linkage); the node of a merge at distance δ has height δ/2 in the ultrametric tree.
//   upgma() runs the nearest-neighbour chain (Benzécri 1982; Juan 1982; Murtagh 1983, The
//   Computer Journal 26, 354–359; Müllner 2011, arXiv:1109.2378, the algorithm NN-chain-linkage):
//   follow nearest neighbours until two clusters are each other's, and merge them. Average
//   linkage is reducible (d(a ∪ b, k) ≥ min(d(a, k), d(b, k))), so the merges found are those of
//   the greedy algorithm, in another order, and sorting them by distance gives its dendrogram.
//   O(n²) time and the n² distances in memory. On a tie the chain prefers its previous element,
//   which is what makes it terminate (Müllner 2011), and otherwise the smaller index.
//   upgma_naive() is the greedy algorithm itself, O(n³), the reference the checks hold it to:
//   the two agree whenever no two distances that matter are equal.
// - The dendrogram is SciPy's linkage: leaves 0 .. n − 1, the k-th merge creating node n + k, a
//   and b its children (a < b), in order of distance.
// - Sequence weights of CLUSTAL W (Thompson, Higgins and Gibson 1994, Nucleic Acids Research
//   22, 4673–4680, §Sequence weighting): a leaf's weight is the sum, over the branches from it
//   to the root, of the branch's length divided by the number of leaves below it, so that a
//   group of close relatives shares its common branches and an isolated sequence keeps its own.
//   Here normalised to mean 1.
//
// *Remark.* UPGMA, weights and corrected distances are general methods of phylogenetics and
// clustering; they are written here because the evolution domain stands on this one and a
// foundation for inference does not yet exist. They should move to that foundation's
// clustering when it does, and the evolution domain's phylogeny may then reuse them.

#pragma once

#include <genomics/msa/msa.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::msa {

struct DistanceMatrix {
  std::size_t n = 0;
  std::vector<double> d;  // row-major, n × n, symmetric, 0 on the diagonal

  DistanceMatrix() = default;
  explicit DistanceMatrix(std::size_t order) : n(order), d(order * order, 0.0) {}
  double operator()(std::size_t i, std::size_t j) const { return d[i * n + j]; }
  double& at(std::size_t i, std::size_t j) { return d[i * n + j]; }
};

// ---- compressed alphabets ------------------------------------------------------------------------------------
// The class of each code (−1: skipped), and the number of classes.
struct Compression {
  std::string name;
  int classes = 0;
  std::array<std::int8_t, 256> class_of{};
};

inline Compression dna_letters() {
  Compression c{"ACGT", 4, {}};
  c.class_of.fill(-1);
  for (int a = 0; a < 4; a++) c.class_of[static_cast<std::size_t>(a)] = static_cast<std::int8_t>(a);
  return c;
}
inline Compression protein_letters() {
  Compression c{"20 amino acids", 20, {}};
  c.class_of.fill(-1);
  for (int a = 0; a < 20; a++) c.class_of[static_cast<std::size_t>(a)] = static_cast<std::int8_t>(a);
  return c;
}
inline Compression dayhoff6() {
  Compression c{"Dayhoff(6)", 6, {}};
  c.class_of.fill(-1);
  const char* groups[6] = {"AGPST", "C", "DENQBZ", "HKR", "ILMV", "FWY"};
  for (int g = 0; g < 6; g++)
    for (const char* p = groups[g]; *p; p++) c.class_of[sequence::protein_code(*p)] = static_cast<std::int8_t>(g);
  return c;
}

// ---- k-mer distances ------------------------------------------------------------------------------------
// The k-mers of a sequence over the compressed alphabet, sorted and counted.
struct KmerProfile {
  std::vector<std::uint64_t> kmers;
  std::vector<std::uint32_t> counts;
  std::size_t total = 0;  // k-mers kept
};

inline KmerProfile kmer_profile(const Codes& x, int k, const Compression& c) {
  if (k < 1 || c.classes < 2 || std::pow(static_cast<double>(c.classes), k) > 1.8e19)
    throw std::invalid_argument("genomics::msa::kmer_profile: k >= 1 and classes^k below 2^64");
  std::vector<std::uint64_t> all;
  all.reserve(x.size());
  std::uint64_t top = 1;  // classes^(k − 1)
  for (int i = 1; i < k; i++) top *= static_cast<std::uint64_t>(c.classes);
  std::uint64_t v = 0;
  int valid = 0;
  for (Code a : x) {
    const int cl = a == GAP ? -1 : c.class_of[a];
    if (cl < 0) { valid = 0; v = 0; continue; }
    v = (valid >= k ? v % top : v) * static_cast<std::uint64_t>(c.classes) + static_cast<std::uint64_t>(cl);
    if (++valid >= k) all.push_back(v);
  }
  std::sort(all.begin(), all.end());
  KmerProfile p;
  p.total = all.size();
  for (std::size_t i = 0; i < all.size();) {
    std::size_t j = i;
    while (j < all.size() && all[j] == all[i]) j++;
    p.kmers.push_back(all[i]);
    p.counts.push_back(static_cast<std::uint32_t>(j - i));
    i = j;
  }
  return p;
}

// F(x, y): the shared k-mers over the k-mers of the shorter; 0 when either has none.
inline double kmer_similarity(const KmerProfile& x, const KmerProfile& y) {
  const std::size_t denominator = std::min(x.total, y.total);
  if (denominator == 0) return 0;
  std::size_t shared = 0, i = 0, j = 0;
  while (i < x.kmers.size() && j < y.kmers.size()) {
    if (x.kmers[i] < y.kmers[j]) i++;
    else if (y.kmers[j] < x.kmers[i]) j++;
    else { shared += std::min(x.counts[i], y.counts[j]); i++; j++; }
  }
  return std::min(1.0, static_cast<double>(shared) / static_cast<double>(denominator));
}

// d = 1 − F for every pair, in parallel by rows (each pair written by one task).
inline DistanceMatrix kmer_distances(const std::vector<Codes>& seqs, int k, const Compression& c, calx::scheduler::ThreadPool& pool) {
  const std::size_t n = seqs.size();
  std::vector<KmerProfile> profiles(n);
  calx::scheduler::parallel_for(pool, n, 4, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++) profiles[i] = kmer_profile(seqs[i], k, c);
  });
  DistanceMatrix d(n);
  calx::scheduler::parallel_for(pool, n, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++)
      for (std::size_t j = i + 1; j < n; j++) {
        const double x = 1.0 - kmer_similarity(profiles[i], profiles[j]);
        d.d[i * n + j] = x;
        d.d[j * n + i] = x;
      }
  });
  return d;
}

// ---- distances from an alignment ------------------------------------------------------------------------------------
// The fraction of differing residues among the columns where both rows hold one below
// `letters`; 1 when there is no such column.
inline double p_distance(const Codes& a, const Codes& b, int letters_) {
  std::size_t same = 0, total = 0;
  for (std::size_t c = 0; c < a.size(); c++)
    if (a[c] < letters_ && b[c] < letters_) {
      total++;
      same += a[c] == b[c];
    }
  return total ? 1.0 - static_cast<double>(same) / static_cast<double>(total) : 1.0;
}
inline double kimura_protein(double D) { return -std::log(std::max(0.05, 1 - D - 0.2 * D * D)); }
inline double jukes_cantor(double D) { return -0.75 * std::log(std::max(0.05, 1 - 4.0 * D / 3.0)); }

inline DistanceMatrix corrected_distances(const Msa& m, Alphabet alphabet, calx::scheduler::ThreadPool& pool) {
  require_rectangular(m, "genomics::msa::corrected_distances");
  const std::size_t n = m.size();
  DistanceMatrix d(n);
  calx::scheduler::parallel_for(pool, n, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++)
      for (std::size_t j = i + 1; j < n; j++) {
        const double D = p_distance(m.rows[i], m.rows[j], letters(alphabet));
        const double x = alphabet == Alphabet::protein ? kimura_protein(D) : jukes_cantor(D);
        d.d[i * n + j] = x;
        d.d[j * n + i] = x;
      }
  });
  return d;
}

// ---- UPGMA ------------------------------------------------------------------------------------
struct Merge {
  std::size_t a = 0, b = 0;  // the children, a < b
  double distance = 0;       // the average distance between them; the node's height is half of it
  std::size_t size = 0;      // leaves below the node
  friend bool operator==(const Merge&, const Merge&) = default;
};

struct Dendrogram {
  std::size_t leaves = 0;
  std::vector<Merge> merges;  // merge k creates node leaves + k

  std::size_t nodes() const { return leaves == 0 ? 0 : 2 * leaves - 1; }
  std::size_t root() const { return nodes() - 1; }
  bool is_leaf(std::size_t v) const { return v < leaves; }
  const Merge& merge_of(std::size_t v) const { return merges[v - leaves]; }
  double height(std::size_t v) const { return v < leaves ? 0.0 : merges[v - leaves].distance / 2; }

  // The parent of each node (SIZE_MAX at the root).
  std::vector<std::size_t> parents() const {
    std::vector<std::size_t> p(nodes(), SIZE_MAX);
    for (std::size_t k = 0; k < merges.size(); k++) p[merges[k].a] = p[merges[k].b] = leaves + k;
    return p;
  }
  // The leaves below v, left (a) before right (b).
  std::vector<std::size_t> leaves_under(std::size_t v) const {
    std::vector<std::size_t> out, stack{v};
    while (!stack.empty()) {
      const std::size_t u = stack.back();
      stack.pop_back();
      if (u < leaves) { out.push_back(u); continue; }
      stack.push_back(merges[u - leaves].b);
      stack.push_back(merges[u - leaves].a);
    }
    return out;
  }
  // The tree in Newick, branch lengths the differences of heights.
  std::string newick(const std::vector<std::string>& names) const {
    if (leaves == 0) return ";";
    const auto label = [&](std::size_t leaf) { return leaf < names.size() ? names[leaf] : std::to_string(leaf); };
    const auto number = [](double x) {
      char buf[32];
      std::snprintf(buf, sizeof buf, "%.6g", x);
      return std::string(buf);
    };
    const auto text = [&](auto&& self, std::size_t v) -> std::string {
      if (v < leaves) return label(v);
      const Merge& m = merges[v - leaves];
      return "(" + self(self, m.a) + ":" + number(height(v) - height(m.a)) + "," + self(self, m.b) + ":" + number(height(v) - height(m.b)) + ")";
    };
    return text(text, root()) + ";";
  }
};

namespace detail {

struct RawMerge {
  std::size_t x = 0, y = 0;  // slots: each the smallest leaf of its cluster at the time
  double distance = 0;
};

// From merges of slots to SciPy's linkage: sort by distance (stably), then name the clusters by
// union–find over the leaves.
inline Dendrogram label(std::size_t n, std::vector<RawMerge> raw) {
  std::stable_sort(raw.begin(), raw.end(), [](const RawMerge& p, const RawMerge& q) { return p.distance < q.distance; });
  Dendrogram t;
  t.leaves = n;
  std::vector<std::size_t> parent(2 * n, SIZE_MAX), size(2 * n, 1);
  const auto find = [&](std::size_t x) {
    std::size_t r = x;
    while (parent[r] != SIZE_MAX) r = parent[r];
    while (parent[x] != SIZE_MAX) { const std::size_t next = parent[x]; parent[x] = r; x = next; }
    return r;
  };
  for (std::size_t k = 0; k < raw.size(); k++) {
    const std::size_t a = find(raw[k].x), b = find(raw[k].y), id = n + k;
    parent[a] = parent[b] = id;
    size[id] = size[a] + size[b];
    t.merges.push_back({std::min(a, b), std::max(a, b), raw[k].distance, size[id]});
  }
  return t;
}

inline void validate(const DistanceMatrix& d, const char* who) {
  if (d.d.size() != d.n * d.n) throw std::invalid_argument(std::string(who) + ": the matrix is not n × n");
  for (double x : d.d)
    if (!std::isfinite(x)) throw std::invalid_argument(std::string(who) + ": a distance that is not finite");
}

}  // namespace detail

// UPGMA by the nearest-neighbour chain, O(n²).
inline Dendrogram upgma(const DistanceMatrix& input) {
  detail::validate(input, "genomics::msa::upgma");
  const std::size_t n = input.n;
  if (n == 0) return {};
  std::vector<double> D = input.d;
  std::vector<std::size_t> size(n, 1);
  std::vector<bool> alive(n, true);
  std::vector<std::size_t> chain;
  chain.reserve(n);
  std::vector<detail::RawMerge> raw;
  raw.reserve(n - 1);
  std::size_t remaining = n, first_alive = 0;
  while (remaining > 1) {
    if (chain.empty()) {
      while (!alive[first_alive]) first_alive++;
      chain.push_back(first_alive);
    }
    const std::size_t a = chain.back();
    const std::size_t prev = chain.size() >= 2 ? chain[chain.size() - 2] : SIZE_MAX;
    std::size_t best = prev;
    double bd = prev != SIZE_MAX ? D[a * n + prev] : std::numeric_limits<double>::infinity();
    const double* row = D.data() + a * n;
    for (std::size_t c = 0; c < n; c++)
      if (alive[c] && c != a && row[c] < bd) { bd = row[c]; best = c; }
    if (best != prev || prev == SIZE_MAX) {
      chain.push_back(best);
      continue;
    }
    chain.pop_back();
    chain.pop_back();
    const std::size_t lo = std::min(a, prev), hi = std::max(a, prev);
    raw.push_back({lo, hi, bd});
    const double wa = static_cast<double>(size[lo]), wb = static_cast<double>(size[hi]);
    for (std::size_t k = 0; k < n; k++) {
      if (!alive[k] || k == lo || k == hi) continue;
      const double x = (wa * D[lo * n + k] + wb * D[hi * n + k]) / (wa + wb);
      D[lo * n + k] = D[k * n + lo] = x;
    }
    size[lo] += size[hi];
    alive[hi] = false;
    remaining--;
  }
  return detail::label(n, std::move(raw));
}

// UPGMA by the greedy algorithm, O(n³): the reference.
inline Dendrogram upgma_naive(const DistanceMatrix& input) {
  detail::validate(input, "genomics::msa::upgma_naive");
  const std::size_t n = input.n;
  if (n == 0) return {};
  std::vector<double> D = input.d;
  std::vector<std::size_t> size(n, 1);
  std::vector<bool> alive(n, true);
  std::vector<detail::RawMerge> raw;
  for (std::size_t step = 0; step + 1 < n; step++) {
    std::size_t bi = SIZE_MAX, bj = SIZE_MAX;
    double bd = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; i++) {
      if (!alive[i]) continue;
      for (std::size_t j = i + 1; j < n; j++)
        if (alive[j] && D[i * n + j] < bd) { bd = D[i * n + j]; bi = i; bj = j; }
    }
    raw.push_back({bi, bj, bd});
    const double wa = static_cast<double>(size[bi]), wb = static_cast<double>(size[bj]);
    for (std::size_t k = 0; k < n; k++) {
      if (!alive[k] || k == bi || k == bj) continue;
      const double x = (wa * D[bi * n + k] + wb * D[bj * n + k]) / (wa + wb);
      D[bi * n + k] = D[k * n + bi] = x;
    }
    size[bi] += size[bj];
    alive[bj] = false;
  }
  return detail::label(n, std::move(raw));
}

// CLUSTAL W's weights from a rooted tree, normalised to mean 1 (all 1 when every branch is 0).
inline std::vector<double> clustal_weights(const Dendrogram& t) {
  const std::size_t n = t.leaves;
  std::vector<double> w(n, 0);
  if (n <= 1) return std::vector<double>(n, 1.0);
  const auto parent = t.parents();
  std::vector<std::size_t> below(t.nodes(), 1);
  for (std::size_t k = 0; k < t.merges.size(); k++) below[n + k] = t.merges[k].size;
  for (std::size_t leaf = 0; leaf < n; leaf++)
    for (std::size_t v = leaf; parent[v] != SIZE_MAX; v = parent[v])
      w[leaf] += (t.height(parent[v]) - t.height(v)) / static_cast<double>(below[v]);
  const double total = std::accumulate(w.begin(), w.end(), 0.0);
  if (!(total > 0)) return std::vector<double>(n, 1.0);
  for (double& x : w) x *= static_cast<double>(n) / total;
  return w;
}

}  // namespace genomics::msa
