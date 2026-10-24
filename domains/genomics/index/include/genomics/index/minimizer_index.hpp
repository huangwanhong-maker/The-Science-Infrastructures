// The minimizer index of a genome: every (w, k)-minimizer of its contigs, with where it occurs,
// looked up by its hash; the seed index of a minimap2-like mapper (Li 2018, "Minimap2: pairwise
// alignment for nucleotide sequences", Bioinformatics 34, 3094–3100, §2.1).
//
// - What is indexed. The minimizers of genomics::sequence::minimizers (hash64 of the canonical
//   k-mer, every minimal position of each window kept, palindromic k-mers skipped; Roberts et al.
//   2004), each as (contig, position, strand): the strand says whether the canonical k-mer is
//   the reverse complement of the forward one.
// - Construction, in parallel and deterministic. The contigs are cut into chunks of a fixed
//   length (not a number of threads), each taken with a margin of w − 1 bases before it and
//   w + k − 2 after, so that every window holding a k-mer that starts in the chunk lies inside;
//   the minimizers of the extended chunk starting in the chunk are exactly the genome's there
//   (a window cut by the margin holds no k-mer of the chunk, and a run of fewer than w k-mers cut
//   by it starts outside). The chunks are computed on the pool, each into its own slot; the
//   entries are scattered into 2^b buckets by the top bits of their hashes, at offsets counted
//   in chunk order; each bucket is sorted on the pool by (hash, contig, position, strand). The
//   result is one flat array sorted by hash, whatever the number of threads. O(N/t + N log N / t)
//   for N bases on t threads.
// - Lookup. A table with open addressing (linear probing, load 0.7, the slot chosen by Fibonacci
//   hashing of the key, Knuth 1998, TAOCP 3, §6.4, reduced to the table's size by a multiply and a
//   shift, Lemire 2019, Software: Practice and Experience 49, 1–14) maps each distinct hash to its
//   range [offset, offset + count) of the flat array: O(1) expected per lookup.
// - Repetitive minimizers are ignored as minimap2's -f does (mm_idx_cal_max_occ): with c the
//   ⌊(1 − f) n⌋-th smallest of the n distinct minimizers' counts, a minimizer occurring more than
//   c + 1 times gives no anchors. f = 2·10⁻⁴ by default, minimap2's.
// - Anchors. Each minimizer of a query that is in the index and not ignored pairs with each of its
//   occurrences: (contig, target position, query position, reverse), positions being the first
//   base of the k-mer on each forward strand, and reverse when the two strands differ. An anchor
//   on the reverse strand places query base j + i against target base t + k − 1 − i. Sorted by
//   (reverse, contig, target, query), the order chaining reads them in.
//
// Space: 8 bytes an occurrence and 16 bytes a slot of the table, 23 bytes a distinct minimizer:
// about 5.6 bytes a base at (w, k) = (10, 15), where nearly every minimizer is distinct (26 MB
// for Escherichia coli). Exact: no probabilistic step.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::index {

using sequence::Code;
using sequence::Codes;

struct MinimizerOptions {
  int k = 15, w = 10;
  double ignore = 2e-4;            // the fraction of the most frequent distinct minimizers ignored
  std::size_t chunk = 1 << 20;     // bases a task (fixed, so that the result is the same on any pool)
};

// An occurrence of a minimizer in the genome.
struct Occurrence {
  std::uint32_t contig = 0, position = 0;
  bool reverse = false;
  friend bool operator==(const Occurrence&, const Occurrence&) = default;
};

// A query minimizer matched to an occurrence in the genome.
struct Anchor {
  std::uint32_t contig = 0, target = 0, query = 0;
  bool reverse = false;
  friend bool operator==(const Anchor&, const Anchor&) = default;
  friend bool operator<(const Anchor& a, const Anchor& b) {
    if (a.reverse != b.reverse) return a.reverse < b.reverse;
    if (a.contig != b.contig) return a.contig < b.contig;
    if (a.target != b.target) return a.target < b.target;
    return a.query < b.query;
  }
};

class MinimizerIndex {
 public:
  MinimizerIndex() = default;

  MinimizerIndex(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, MinimizerOptions options = {}) : options_(options) {
    sequence::check_k(options.k, "MinimizerIndex");
    if (options.w < 1 || options.chunk == 0) throw std::invalid_argument("genomics::index::MinimizerIndex: w >= 1 and a chunk of at least one base");
    for (const auto& c : contigs)
      if (c.size() >= (1ULL << 31)) throw std::invalid_argument("genomics::index::MinimizerIndex: contigs shorter than 2^31 bases");
    // The chunks: (contig, begin, end).
    struct Chunk { std::uint32_t contig; std::size_t begin, end; };
    std::vector<Chunk> chunks;
    for (std::uint32_t i = 0; i < contigs.size(); i++)
      for (std::size_t b = 0; b < contigs[i].size(); b += options.chunk) chunks.push_back({i, b, std::min(contigs[i].size(), b + options.chunk)});
    const std::size_t before = static_cast<std::size_t>(options.w - 1), after = static_cast<std::size_t>(options.w + options.k - 2);
    std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>> found(chunks.size());
    calx::scheduler::parallel_for(pool, chunks.size(), 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
      for (std::size_t c = lo; c < hi; c++) {
        const Chunk& ch = chunks[c];
        const Codes& seq = contigs[ch.contig];
        const std::size_t from = ch.begin >= before ? ch.begin - before : 0, to = std::min(seq.size(), ch.end + after);
        const Codes slice(seq.begin() + static_cast<std::ptrdiff_t>(from), seq.begin() + static_cast<std::ptrdiff_t>(to));
        for (const auto& m : sequence::minimizers(slice, options_.k, options_.w)) {
          const std::size_t p = from + m.position;
          if (p >= ch.begin && p < ch.end) found[c].push_back({m.hash, pack(ch.contig, static_cast<std::uint32_t>(p), m.reverse)});
        }
      }
    });
    // Buckets by the top bits of the hash (hashes lie below 4^k).
    const int key_bits = 2 * options.k;
    const int bucket_bits = std::min(key_bits, 10);
    const std::size_t buckets = std::size_t{1} << bucket_bits;
    const auto bucket = [&](std::uint64_t h) { return static_cast<std::size_t>(h >> (key_bits - bucket_bits)); };
    std::vector<std::size_t> counts(chunks.size() * buckets, 0);
    calx::scheduler::parallel_for(pool, chunks.size(), 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
      for (std::size_t c = lo; c < hi; c++)
        for (const auto& e : found[c]) counts[c * buckets + bucket(e.first)]++;
    });
    // Offsets in (bucket, chunk) order.
    std::vector<std::size_t> offset(chunks.size() * buckets), bucket_start(buckets + 1, 0);
    std::size_t total = 0;
    for (std::size_t b = 0; b < buckets; b++) {
      bucket_start[b] = total;
      for (std::size_t c = 0; c < chunks.size(); c++) {
        offset[c * buckets + b] = total;
        total += counts[c * buckets + b];
      }
    }
    bucket_start[buckets] = total;
    if (total >= std::numeric_limits<std::uint32_t>::max()) throw std::invalid_argument("genomics::index::MinimizerIndex: more than 2^32 occurrences");
    std::vector<std::pair<std::uint64_t, std::uint64_t>> all(total);
    calx::scheduler::parallel_for(pool, chunks.size(), 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
      for (std::size_t c = lo; c < hi; c++) {
        std::size_t* at = &offset[c * buckets];
        for (const auto& e : found[c]) all[at[bucket(e.first)]++] = e;
      }
    });
    found.clear();
    found.shrink_to_fit();
    calx::scheduler::parallel_for(pool, buckets, 4, [&](std::size_t lo, std::size_t hi, std::size_t) {
      for (std::size_t b = lo; b < hi; b++)
        std::sort(all.begin() + static_cast<std::ptrdiff_t>(bucket_start[b]), all.begin() + static_cast<std::ptrdiff_t>(bucket_start[b + 1]));
    });
    // The flat array and the distinct hashes with their ranges.
    entries_.resize(total);
    std::vector<std::uint64_t> keys;
    std::vector<std::uint32_t> starts;
    for (std::size_t i = 0; i < total; i++) {
      entries_[i] = all[i].second;
      if (i == 0 || all[i].first != all[i - 1].first) {
        keys.push_back(all[i].first);
        starts.push_back(static_cast<std::uint32_t>(i));
      }
    }
    starts.push_back(static_cast<std::uint32_t>(total));
    all.clear();
    all.shrink_to_fit();
    distinct_ = keys.size();
    // The table, of ⌈distinct / 0.7⌉ slots.
    table_.assign(std::max<std::size_t>(16, static_cast<std::size_t>(std::ceil(static_cast<double>(distinct_) / 0.7))), Slot{});
    std::vector<std::uint32_t> occurrence_counts;
    occurrence_counts.reserve(distinct_);
    for (std::size_t d = 0; d < distinct_; d++) {
      const std::uint32_t n = starts[d + 1] - starts[d];
      occurrence_counts.push_back(n);
      std::size_t s = slot(keys[d]);
      while (table_[s].count != 0) s = s + 1 == table_.size() ? 0 : s + 1;
      table_[s] = {keys[d], starts[d], n};
    }
    // The cut-off of repetitive minimizers.
    max_occ_ = std::numeric_limits<std::uint32_t>::max();
    if (options.ignore > 0 && !occurrence_counts.empty()) {
      const auto at = static_cast<std::size_t>(std::floor((1 - options.ignore) * static_cast<double>(occurrence_counts.size())));
      const std::size_t q = std::min(at, occurrence_counts.size() - 1);
      std::nth_element(occurrence_counts.begin(), occurrence_counts.begin() + static_cast<std::ptrdiff_t>(q), occurrence_counts.end());
      max_occ_ = occurrence_counts[q] + 1;
    }
  }

  int k() const { return options_.k; }
  int w() const { return options_.w; }
  const MinimizerOptions& options() const { return options_; }
  std::size_t distinct() const { return distinct_; }
  std::size_t occurrences() const { return entries_.size(); }
  // A minimizer occurring more often than this gives no anchors.
  std::uint32_t max_occurrences() const { return max_occ_; }

  // The occurrences of a hash, packed, in the order (contig, position, strand); empty if absent.
  std::span<const std::uint64_t> lookup(std::uint64_t hash) const {
    if (table_.empty()) return {};
    for (std::size_t s = slot(hash);; s = s + 1 == table_.size() ? 0 : s + 1) {
      const Slot& x = table_[s];
      if (x.count == 0) return {};
      if (x.hash == hash) return {entries_.data() + x.offset, x.count};
    }
  }
  static Occurrence decode(std::uint64_t packed) {
    return {static_cast<std::uint32_t>(packed >> 33), static_cast<std::uint32_t>((packed >> 1) & 0xFFFFFFFFULL), static_cast<bool>(packed & 1)};
  }

  // The anchors of a query, sorted by (reverse, contig, target, query).
  std::vector<Anchor> anchors(const Codes& query) const {
    std::vector<Anchor> out;
    for (const auto& m : sequence::minimizers(query, options_.k, options_.w)) {
      const auto hits = lookup(m.hash);
      if (hits.size() > max_occ_) continue;
      for (const std::uint64_t e : hits) {
        const Occurrence o = decode(e);
        out.push_back({o.contig, o.position, m.position, o.reverse != m.reverse});
      }
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  // The flat array of occurrences sorted by (hash, contig, position, strand), for whoever needs to
  // compare indexes or write them out.
  const std::vector<std::uint64_t>& entries() const { return entries_; }

  std::size_t bytes() const { return entries_.size() * 8 + table_.size() * sizeof(Slot); }

 private:
  struct Slot {
    std::uint64_t hash = 0;
    std::uint32_t offset = 0, count = 0;  // count 0: an empty slot
  };
  static std::uint64_t pack(std::uint32_t contig, std::uint32_t position, bool reverse) {
    return (static_cast<std::uint64_t>(contig) << 33) | (static_cast<std::uint64_t>(position) << 1) | static_cast<std::uint64_t>(reverse);
  }
  // Fibonacci hashing, reduced to [0, slots) by a multiply and a shift (Lemire 2019).
  std::size_t slot(std::uint64_t hash) const {
    return static_cast<std::size_t>((static_cast<unsigned __int128>(hash * 0x9e3779b97f4a7c15ULL) * table_.size()) >> 64);
  }

  MinimizerOptions options_;
  std::vector<std::uint64_t> entries_;
  std::vector<Slot> table_;
  std::size_t distinct_ = 0;
  std::uint32_t max_occ_ = std::numeric_limits<std::uint32_t>::max();
};

}  // namespace genomics::index
