// Taxonomic classification of reads by exact k-mers and the lowest common ancestor: Kraken's
// method (Wood and Salzberg 2014, "Kraken: ultrafast metagenomic sequence classification using
// exact alignments", Genome Biology 15, R46).
//
// - The database maps each canonical k-mer of the reference genomes to the LCA of the taxa of
//   the genomes that contain it: a k-mer of one genome maps to that genome's taxon, one shared by
//   two species of a genus to the genus, one shared by a bacterium and a phage to the root. It is
//   built in one pass, each k-mer's taxon met with the one already stored; the LCA being a
//   semilattice operation (taxonomy.hpp), the result does not depend on the order of the genomes
//   or of their k-mers. It is held in the kmers module's open-addressing table (`KmerMap`, 12
//   bytes a slot at a load of at most 5/8), O(1) expected per k-mer, O(depth) for each LCA.
// - Classification (ibid., §"Sequence classification algorithm"). Each k-mer of the read free of
//   N is looked up, and each one found adds one to its taxon's hits. The classification tree is
//   the set of hit taxa with all their ancestors; each of its root-to-leaf paths is weighted by
//   the hits of the taxa along it, and the read is given the leaf of the heaviest path, or the
//   LCA of the leaves of all the heaviest paths when several tie. A read with no k-mer in the
//   database is unclassified.
//   *Remark.* As the hits are positive, a hit taxon with a hit descendant always weighs less than
//   that descendant's path, so the heaviest paths are found among the hit taxa themselves, with
//   no need to build the classification tree (theory, Proposition 2.6): O(m² · depth) for m
//   distinct taxa hit, a few in practice.
// - Thresholds (ours). A read needs `min_hits` k-mers in the database to be classified (1 is
//   Kraken's rule). With `confidence` τ > 0, the label then climbs toward the root until the
//   k-mers whose taxa lie in its clade are at least τ of the read's k-mers free of N, and the read
//   is unclassified when even the root's clade, every hit, falls short.
//
// Limits. Exact k-mers make classification as specific as the database allows and no more: a
// read of an organism absent from it goes to its nearest relatives present, at the level of
// their shared k-mers, which may be a relative's own leaf (theory §2.4). A sequencing error
// removes up to k k-mers from a read, rarely adds a false one.

#pragma once

#include <genomics/classify/taxonomy.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::classify {

using sequence::Codes;
using sequence::Kmer;

class Database {
 public:
  explicit Database(Taxonomy taxonomy, int k = 31) : taxonomy_(std::move(taxonomy)), k_(k), map_(1 << 16) { sequence::check_k(k, "classify::Database"); }

  // Every canonical k-mer of a genome (both strands at once), under its taxon.
  void add(const Codes& genome, TaxonId taxon) {
    if (taxon >= taxonomy_.size()) throw std::invalid_argument("genomics::classify::Database::add: a taxon of the taxonomy");
    sequence::for_each_kmer(genome, k_, [&](std::size_t, Kmer f, Kmer r) {
      const Kmer x = f < r ? f : r;
      std::uint32_t& v = map_.at(x, kmers::mix64(x));
      v = v == 0 ? taxon + 1 : taxonomy_.lca(v - 1, taxon) + 1;  // 0 is a fresh slot
    });
    genomes_++;
  }

  // The taxon of a canonical k-mer, or NO_TAXON when no genome has it.
  TaxonId lookup(Kmer canonical) const {
    const std::uint32_t* v = map_.find(canonical, kmers::mix64(canonical));
    return v ? *v - 1 : NO_TAXON;
  }
  void prefetch(Kmer canonical) const { map_.prefetch(kmers::mix64(canonical)); }

  const Taxonomy& taxonomy() const { return taxonomy_; }
  int k() const { return k_; }
  std::size_t size() const { return map_.size(); }  // distinct canonical k-mers
  std::size_t genomes() const { return genomes_; }
  std::size_t memory_bytes() const { return map_.memory_bytes(); }
  // f(canonical k-mer, taxon) for every entry, in slot order.
  template <class F>
  void for_each(F&& f) const {
    map_.for_each([&](Kmer x, std::uint32_t v) { f(x, static_cast<TaxonId>(v - 1)); });
  }

 private:
  Taxonomy taxonomy_;
  int k_;
  kmers::KmerMap map_;
  std::size_t genomes_ = 0;
};

struct ClassifyOptions {
  std::uint32_t min_hits = 1;  // k-mers in the database a read needs
  double confidence = 0;       // τ: the fraction of the read's k-mers its label's clade must hold
};

struct Classification {
  TaxonId taxon = NO_TAXON;      // the label; NO_TAXON when unclassified
  std::uint32_t kmers = 0;       // k-mers of the read free of N
  std::uint32_t hits = 0;        // of them found in the database
  std::uint32_t clade_hits = 0;  // of them with taxa in the label's clade
  bool classified() const { return taxon != NO_TAXON; }
};

// The label of a read from the hits of its taxa, (taxon, count) with distinct taxa and positive
// counts: the leaf of the heaviest root-to-leaf path of the classification tree, or the LCA of
// the leaves of the heaviest paths; NO_TAXON when there are no hits.
inline TaxonId heaviest_path(const Taxonomy& tax, const std::vector<std::pair<TaxonId, std::uint32_t>>& hits) {
  std::uint64_t best = 0;
  TaxonId label = NO_TAXON;
  for (const auto& [t, c] : hits) {
    std::uint64_t w = 0;
    for (const auto& [s, d] : hits)
      if (tax.above(s, t)) w += d;
    if (w > best) {
      best = w;
      label = t;
    } else if (w == best && label != NO_TAXON) {
      label = tax.lca(label, t);
    }
  }
  return label;
}

// The hits in the clade of `label` (every taxon at or below it).
inline std::uint32_t clade_hits(const Taxonomy& tax, TaxonId label, const std::vector<std::pair<TaxonId, std::uint32_t>>& hits) {
  std::uint32_t n = 0;
  for (const auto& [s, d] : hits)
    if (tax.above(label, s)) n += d;
  return n;
}

namespace detail {

// The (taxon, count) pairs of a read's hit taxa, sorted by taxon.
inline void count_hits(std::vector<TaxonId>& taxa, std::vector<std::pair<TaxonId, std::uint32_t>>& out) {
  std::sort(taxa.begin(), taxa.end());
  out.clear();
  for (std::size_t i = 0; i < taxa.size();) {
    std::size_t j = i;
    while (j < taxa.size() && taxa[j] == taxa[i]) j++;
    out.push_back({taxa[i], static_cast<std::uint32_t>(j - i)});
    i = j;
  }
}

}  // namespace detail

// The classification of one read.
inline Classification classify(const Database& db, const Codes& read, const ClassifyOptions& o = {}) {
  Classification c;
  std::vector<TaxonId> taxa;
  taxa.reserve(read.size());
  sequence::for_each_kmer(read, db.k(), [&](std::size_t, Kmer f, Kmer r) {
    c.kmers++;
    const TaxonId t = db.lookup(f < r ? f : r);
    if (t != NO_TAXON) taxa.push_back(t);
  });
  c.hits = static_cast<std::uint32_t>(taxa.size());
  if (c.hits == 0 || c.hits < o.min_hits) return c;
  std::vector<std::pair<TaxonId, std::uint32_t>> hits;
  detail::count_hits(taxa, hits);
  const Taxonomy& tax = db.taxonomy();
  TaxonId label = heaviest_path(tax, hits);
  std::uint32_t inside = clade_hits(tax, label, hits);
  if (o.confidence > 0) {
    const double need = o.confidence * static_cast<double>(c.kmers);
    while (static_cast<double>(inside) < need && label != ROOT) {
      label = tax[label].parent;
      inside = clade_hits(tax, label, hits);
    }
    if (static_cast<double>(inside) < need) return c;
  }
  c.taxon = label;
  c.clade_hits = inside;
  return c;
}

// The classifications of many reads on the pool, each written to its own slot: the same on any
// number of threads.
inline std::vector<Classification> classify_all(const Database& db, const std::vector<Codes>& reads, calx::scheduler::ThreadPool& pool, const ClassifyOptions& o = {}) {
  std::vector<Classification> out(reads.size());
  calx::scheduler::parallel_for(pool, reads.size(), 256, [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t i = begin; i < end; i++) out[i] = classify(db, reads[i], o);
  });
  return out;
}

}  // namespace genomics::classify
