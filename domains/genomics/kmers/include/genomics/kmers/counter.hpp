// Counting canonical k-mers of many sequences, in parallel and deterministically.
//
// A k-mer and its reverse complement are one piece of a double strand, so they are counted
// under one name, the canonical k-mer (the smaller of the two, sequence/kmer.hpp); k-mers
// containing N are skipped. k ≤ 32.
//
// The scheme is that of Jellyfish and KMC in spirit (Marçais and Kingsford 2011, "A fast,
// lock-free approach for efficient parallel counting of occurrences of k-mers",
// Bioinformatics 27, 764–770; Kokot, Długosz and Deorowicz 2017, "KMC 3: counting and
// manipulating k-mer statistics", Bioinformatics 33, 2759–2761): the k-mers are partitioned
// by the high bits of a hash, and each partition is counted by one task in a table of its
// own, so that no two threads ever write to one table and no locks are needed.
//
// - The hash is SplitMix64's finaliser (Steele, Lea and Flood 2014, "Fast splittable
//   pseudorandom number generators", OOPSLA, §3; the constants of Stafford's Mix13), a
//   bijection of 64-bit words: its high bits choose the partition and its low bits the slot,
//   so the two are independent. (hash64 of sequence/kmer.hpp is invertible on 2k bits but
//   weak below k = 16, see its remark; the partitions need all 64 bits mixed.)
// - A table is open addressing with linear probing (Knuth 1998, TAOCP 3, §6.4, Algorithm L):
//   a power-of-two capacity, keys and 32-bit counts in two arrays, the empty key ~0, which is
//   never a canonical k-mer (the canonical form of T^32 is A^32 = 0). It doubles when its load
//   passes 5/8, where an unsuccessful search probes (1 + 1/(1 − α)²)/2 ≤ 3.9 slots on average
//   (Knuth §6.4, eq. (45)). Counts saturate at 2^32 − 1.
// - Counting runs in batches of about batch_bases bases. A batch is cut into chunks of about
//   chunk_bases bases (a long sequence into pieces overlapping by k − 1, so that no k-mer is
//   lost or counted twice); each chunk, one task, rolls its k-mers twice, first counting how
//   many fall in each partition and then writing them to their place in the chunk's buffer
//   (a counting sort). Then each partition, one task, inserts its k-mers from every chunk in
//   chunk order. The chunks depend only on the input and the options, never on the threads,
//   so every table receives the same keys in the same order on any number of threads, and
//   ends in the same state: counts, iteration and memory are the same on 1 or 64 threads.
//
// Memory: 12 bytes a slot, at a load between 5/16 and 5/8, so 19 to 38 bytes for each
// distinct k-mer; and 8 bytes for each k-mer of a batch (128 MB at the default 16 Mb).
// memory_bytes() says what is held. Time: O(total bases) for the rolling, and an expected
// O(1) probes for each k-mer.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <span>
#include <vector>

namespace genomics::kmers {

using sequence::Code;
using sequence::Codes;
using sequence::Kmer;

// SplitMix64's finaliser: a bijection of 64-bit words with full avalanche.
constexpr std::uint64_t mix64(std::uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

struct KmerCount {
  Kmer kmer = 0;
  std::uint32_t count = 0;
  friend bool operator==(const KmerCount&, const KmerCount&) = default;
};

// A map from canonical k-mers to 32-bit values by open addressing with linear probing; the
// caller passes the k-mer's mix64 with it.
class KmerMap {
 public:
  static constexpr Kmer EMPTY = ~Kmer{0};

  explicit KmerMap(std::size_t capacity = 16) { reset(capacity); }

  std::size_t size() const { return size_; }
  std::size_t capacity() const { return keys_.size(); }
  std::size_t memory_bytes() const { return keys_.capacity() * sizeof(Kmer) + values_.capacity() * sizeof(std::uint32_t); }

  // The value of x, inserted as 0 if absent.
  std::uint32_t& at(Kmer x, std::uint64_t h) {
    if ((size_ + 1) * 8 > keys_.size() * 5) grow();
    std::size_t i = h & mask_;
    for (;;) {
      if (keys_[i] == x) return values_[i];
      if (keys_[i] == EMPTY) {
        keys_[i] = x;
        size_++;
        return values_[i];
      }
      i = (i + 1) & mask_;
    }
  }
  // The value of x, or nullptr.
  const std::uint32_t* find(Kmer x, std::uint64_t h) const {
    std::size_t i = h & mask_;
    for (;;) {
      if (keys_[i] == x) return &values_[i];
      if (keys_[i] == EMPTY) return nullptr;
      i = (i + 1) & mask_;
    }
  }
  void prefetch(std::uint64_t h) const { __builtin_prefetch(&keys_[h & mask_]); }

  // f(kmer, value) for every entry, in slot order.
  template <class F>
  void for_each(F&& f) const {
    for (std::size_t i = 0; i < keys_.size(); i++)
      if (keys_[i] != EMPTY) f(keys_[i], values_[i]);
  }

 private:
  void reset(std::size_t capacity) {
    std::size_t c = 16;
    while (c < capacity) c <<= 1;
    keys_.assign(c, EMPTY);
    values_.assign(c, 0);
    mask_ = c - 1;
    size_ = 0;
  }
  void grow() {
    std::vector<Kmer> keys = std::move(keys_);
    std::vector<std::uint32_t> values = std::move(values_);
    reset(keys.size() * 2);
    for (std::size_t j = 0; j < keys.size(); j++) {
      if (keys[j] == EMPTY) continue;
      std::size_t i = mix64(keys[j]) & mask_;
      while (keys_[i] != EMPTY) i = (i + 1) & mask_;
      keys_[i] = keys[j];
      values_[i] = values[j];
      size_++;
    }
  }

  std::vector<Kmer> keys_;
  std::vector<std::uint32_t> values_;
  std::size_t mask_ = 15, size_ = 0;
};

struct CounterOptions {
  int k = 31;
  unsigned partition_bits = 10;          // 2^bits partitions
  std::size_t chunk_bases = 1u << 18;    // bases of a chunk, the unit of the scatter
  std::size_t batch_bases = 1u << 24;    // bases scattered before they are counted
};

class KmerCounter {
 public:
  explicit KmerCounter(const CounterOptions& options = {}) : options_(options) {
    sequence::check_k(options.k, "kmers::KmerCounter");
    if (options.partition_bits > 16) throw std::invalid_argument("genomics::kmers::KmerCounter: at most 2^16 partitions");
    if (options.chunk_bases == 0 || options.batch_bases == 0) throw std::invalid_argument("genomics::kmers::KmerCounter: chunks and batches of at least one base");
    tables_.resize(std::size_t{1} << options.partition_bits);
  }

  int k() const { return options_.k; }
  const CounterOptions& options() const { return options_; }
  std::size_t partitions() const { return tables_.size(); }
  std::size_t partition_of(std::uint64_t h) const { return options_.partition_bits ? static_cast<std::size_t>(h >> (64 - options_.partition_bits)) : 0; }

  // Counts the canonical k-mers of sequences 0, …, n − 1, get(i) giving the i-th as a const
  // Codes&. May be called many times; the counts add up.
  template <class Get>
  void add(std::size_t n, Get&& get, calx::scheduler::ThreadPool& pool) {
    const std::size_t k = static_cast<std::size_t>(options_.k);
    std::size_t i = 0, offset = 0;  // the next base to cut: sequence i, from offset
    std::vector<Piece> pieces;
    std::vector<std::size_t> chunk_end;  // the end of each chunk in pieces
    while (i < n) {
      // Cut one batch into chunks.
      pieces.clear();
      chunk_end.clear();
      std::size_t batch = 0, chunk = 0;
      while (i < n && batch < options_.batch_bases) {
        const Codes& s = get(i);
        if (s.size() < k) { i++; offset = 0; continue; }
        const std::size_t room = options_.chunk_bases > chunk ? options_.chunk_bases - chunk : 1;
        // A piece [offset, end) of s: its k-mers are those starting in [offset, end − k + 1).
        std::size_t end = s.size();
        if (end - offset > room + k - 1) end = offset + room + k - 1;
        pieces.push_back({i, offset, end});
        const std::size_t bases = end - offset;
        batch += bases;
        chunk += bases;
        if (end == s.size()) { i++; offset = 0; }
        else offset = end - (k - 1);
        if (chunk >= options_.chunk_bases) { chunk_end.push_back(pieces.size()); chunk = 0; }
      }
      if (chunk_end.empty() || chunk_end.back() != pieces.size()) chunk_end.push_back(pieces.size());
      count_batch(pieces, chunk_end, get, pool);
    }
  }
  void add(const std::vector<Codes>& sequences, calx::scheduler::ThreadPool& pool) {
    add(sequences.size(), [&](std::size_t j) -> const Codes& { return sequences[j]; }, pool);
  }

  // The count of a k-mer, given in either orientation.
  std::uint32_t count(Kmer x) const {
    const Kmer c = sequence::canonical(x, options_.k);
    const std::uint64_t h = mix64(c);
    const std::uint32_t* v = tables_[partition_of(h)].find(c, h);
    return v ? *v : 0;
  }
  // Distinct canonical k-mers, and their occurrences.
  std::size_t distinct() const {
    std::size_t d = 0;
    for (const auto& t : tables_) d += t.size();
    return d;
  }
  std::uint64_t total() const { return total_; }
  // Bytes held: the tables and the buffers kept for the next batch.
  std::size_t memory_bytes() const {
    std::size_t b = 0;
    for (const auto& t : tables_) b += t.memory_bytes();
    for (const auto& x : buffers_) b += x.capacity() * sizeof(Kmer);
    for (const auto& x : offsets_) b += x.capacity() * sizeof(std::size_t);
    return b;
  }
  // The largest memory_bytes() seen.
  std::size_t peak_bytes() const { return std::max(peak_, memory_bytes()); }
  // Frees the scatter buffers (they are kept between calls of add).
  void release_buffers() {
    peak_ = peak_bytes();
    buffers_ = {};
    offsets_ = {};
  }

  const KmerMap& table(std::size_t p) const { return tables_.at(p); }

  // f(kmer, count) for every distinct k-mer: partitions in order, slots in order. The order
  // depends on the input and the options, never on the threads.
  template <class F>
  void for_each(F&& f) const {
    for (const auto& t : tables_) t.for_each(f);
  }

  // The k-mers counted at least min_count times, in increasing order of the k-mer: the
  // canonical order, independent of the hash. A bucket sort by the k-mer's leading bits: each
  // partition counts its entries by bucket, then writes them to the buckets' places in
  // partition order, and each bucket is sorted on its own task. O(n) to distribute, O(n log n)
  // in all, every step in parallel.
  std::vector<KmerCount> sorted(std::uint32_t min_count, calx::scheduler::ThreadPool& pool) const {
    const std::size_t P = tables_.size();
    const unsigned bits = static_cast<unsigned>(std::min(10, 2 * options_.k));
    const unsigned shift = static_cast<unsigned>(2 * options_.k) - bits;
    const std::size_t B = std::size_t{1} << bits;
    std::vector<std::vector<std::size_t>> at(P, std::vector<std::size_t>(B + 1, 0));
    calx::scheduler::parallel_for(pool, P, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t p = b; p < e; p++)
        tables_[p].for_each([&](Kmer x, std::uint32_t c) { if (c >= min_count) at[p][(x >> shift) + 1]++; });
    });
    // Bucket-major places: bucket j holds partition 0's entries of j, then partition 1's, …
    std::size_t total = 0;
    for (std::size_t j = 0; j < B; j++)
      for (std::size_t p = 0; p < P; p++) {
        const std::size_t n = at[p][j + 1];
        at[p][j] = total;
        total += n;
      }
    std::vector<std::size_t> bucket_begin(B + 1, total);
    for (std::size_t j = 0; j < B; j++) bucket_begin[j] = P ? at[0][j] : 0;
    std::vector<KmerCount> out(total);
    calx::scheduler::parallel_for(pool, P, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t p = b; p < e; p++)
        tables_[p].for_each([&](Kmer x, std::uint32_t c) { if (c >= min_count) out[at[p][x >> shift]++] = {x, c}; });
    });
    calx::scheduler::parallel_for(pool, B, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t j = b; j < e; j++)
        std::sort(out.begin() + static_cast<std::ptrdiff_t>(bucket_begin[j]), out.begin() + static_cast<std::ptrdiff_t>(bucket_begin[j + 1]),
                  [](const KmerCount& x, const KmerCount& y) { return x.kmer < y.kmer; });
    });
    return out;
  }

 private:
  struct Piece {
    std::size_t sequence = 0, begin = 0, end = 0;
  };

  // f(canonical) for every k-mer starting in [begin, end − k + 1) of s.
  template <class F>
  void roll(const Codes& s, std::size_t begin, std::size_t end, Codes& scratch, F&& f) const {
    (void)scratch;  // the stretch is read in place
    sequence::for_each_kmer(std::span<const Code>(s.data() + begin, end - begin), options_.k, [&](std::size_t, Kmer fw, Kmer rv) { f(fw < rv ? fw : rv); });
  }

  template <class Get>
  void count_batch(const std::vector<Piece>& pieces, const std::vector<std::size_t>& chunk_end, Get& get, calx::scheduler::ThreadPool& pool) {
    const std::size_t C = chunk_end.size(), P = tables_.size();
    if (buffers_.size() < C) {
      buffers_.resize(C);
      offsets_.resize(C);
    }
    std::vector<std::uint64_t> totals(C, 0);
    // A small batch is counted on the calling thread: waking the pool would cost more.
    std::size_t bases = 0;
    for (const auto& q : pieces) bases += q.end - q.begin;
    const std::size_t grain = bases < (std::size_t{1} << 16) ? std::max<std::size_t>(C, P) : 1;
    // Scatter: each chunk counts its k-mers by partition, then writes them in partition order.
    calx::scheduler::parallel_for(pool, C, grain, [&](std::size_t b, std::size_t e, std::size_t) {
      Codes scratch;
      for (std::size_t c = b; c < e; c++) {
        auto& off = offsets_[c];
        off.assign(P + 1, 0);
        const std::size_t first = c ? chunk_end[c - 1] : 0;
        for (std::size_t q = first; q < chunk_end[c]; q++)
          roll(get(pieces[q].sequence), pieces[q].begin, pieces[q].end, scratch, [&](Kmer x) { off[partition_of(mix64(x)) + 1]++; });
        for (std::size_t p = 0; p < P; p++) off[p + 1] += off[p];
        auto& buf = buffers_[c];
        buf.resize(off[P]);
        totals[c] = off[P];
        std::vector<std::size_t> at(off.begin(), off.end() - 1);
        for (std::size_t q = first; q < chunk_end[c]; q++)
          roll(get(pieces[q].sequence), pieces[q].begin, pieces[q].end, scratch, [&](Kmer x) { buf[at[partition_of(mix64(x))]++] = x; });
      }
    });
    // Count: each partition takes its k-mers from every chunk, in chunk order.
    calx::scheduler::parallel_for(pool, P, grain, [&](std::size_t b, std::size_t e, std::size_t) {
      constexpr std::size_t AHEAD = 8;
      for (std::size_t p = b; p < e; p++) {
        KmerMap& t = tables_[p];
        for (std::size_t c = 0; c < C; c++) {
          const Kmer* x = buffers_[c].data() + offsets_[c][p];
          const std::size_t m = offsets_[c][p + 1] - offsets_[c][p];
          for (std::size_t j = 0; j < m; j++) {
            if (j + AHEAD < m) t.prefetch(mix64(x[j + AHEAD]));
            std::uint32_t& v = t.at(x[j], mix64(x[j]));
            if (v != UINT32_MAX) v++;
          }
        }
      }
    });
    for (std::size_t c = 0; c < C; c++) total_ += totals[c];
    peak_ = std::max(peak_, memory_bytes());
  }

  CounterOptions options_;
  std::vector<KmerMap> tables_;
  std::vector<std::vector<Kmer>> buffers_;
  std::vector<std::vector<std::size_t>> offsets_;
  std::uint64_t total_ = 0;
  std::size_t peak_ = 0;
};

// The canonical k-mers of the sequences, counted.
inline KmerCounter count_kmers(const std::vector<Codes>& sequences, calx::scheduler::ThreadPool& pool, const CounterOptions& options = {}) {
  KmerCounter counter(options);
  counter.add(sequences, pool);
  return counter;
}

}  // namespace genomics::kmers
