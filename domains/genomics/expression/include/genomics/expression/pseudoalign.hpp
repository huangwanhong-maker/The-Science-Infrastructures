// Pseudoalignment of RNA-seq reads to a transcriptome by k-mer equivalence classes (Bray,
// Pimentel, Melsted and Pachter 2016, "Near-optimal probabilistic RNA-seq quantification",
// Nature Biotechnology 34, 525–527): which transcripts a read is compatible with, without where
// in them it lies.
//
// - The equivalence class of a k-mer is the set of transcripts that contain it (canonically:
//   the k-mer or its reverse complement, so reads of either strand are placed). The index maps
//   every canonical k-mer of the transcriptome (k = 31 by default, odd, so that no k-mer is its
//   own reverse complement) to its class, and holds each class once. It is built transcript by
//   transcript: a k-mer met first in t gets {t}; a k-mer of class C met in t gets C ∪ {t}, and
//   since t is greater than every member of C the union is C with t appended. A memo of the
//   classes made while reading t sends every k-mer of class C to the same C ∪ {t}, so a class
//   is made once: O(total k-mers) expected time, with the map of kmers::KmerMap (open addressing
//   on SplitMix64's mix, 12 bytes a slot) and the classes in one flat array.
// - The class of a read (or of a pair: both mates' k-mers together) is the intersection of the
//   classes of its k-mers that the index holds; k-mers it does not hold (a sequencing error, or
//   sequence of no transcript) are skipped, as kallisto skips them. A read none of whose k-mers
//   is held, or whose intersection is empty, is not pseudoaligned. The classes form a
//   ∩-semilattice of the subsets of transcripts, and the read's class is the meet of the
//   classes of its k-mers (theory note, §4). Consecutive k-mers of one class are met once, and an
//   intersection is a merge of two sorted lists; every k-mer is looked up, the read's slots all
//   prefetched before the first is read, since each lookup misses the cache. (kallisto also skips
//   along the unitigs of the transcriptome's de Bruijn graph, on which the class is constant;
//   that needs the graph, and is an optimisation of the same answer.) O(L) lookups for a read of
//   length L.
// - An intersection need not be a class of the index (two k-mers of classes {1, 2, 3} and
//   {2, 3, 4} give {2, 3}); such classes are added to the result's table as reads meet them, in
//   the order of the reads, so that the table and the counts are the same on any number of
//   threads: chunks of reads are pseudoaligned in parallel, each read writing its own slot, and
//   the new classes are numbered afterwards in read order.

#pragma once

#include <genomics/expression/transcriptome.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace genomics::expression {

using sequence::Kmer;

// Sets of transcripts, each sorted and held once, numbered in order of first insertion.
class ClassTable {
 public:
  static constexpr std::uint32_t NONE = UINT32_MAX;

  std::size_t size() const { return begin_.size() - 1; }
  std::span<const std::uint32_t> members(std::uint32_t c) const {
    return {members_.data() + begin_[c], static_cast<std::size_t>(begin_[c + 1] - begin_[c])};
  }
  // The class with exactly these members, or NONE.
  std::uint32_t find(std::span<const std::uint32_t> set) const {
    const auto it = head_.find(hash(set));
    if (it == head_.end()) return NONE;
    for (std::uint32_t c = it->second; c != NONE; c = next_[c])
      if (equal(members(c), set)) return c;
    return NONE;
  }
  // The class with these (sorted) members, made if new.
  std::uint32_t insert(std::span<const std::uint32_t> set) {
    const auto [it, fresh] = head_.try_emplace(hash(set), NONE);
    for (std::uint32_t c = it->second; c != NONE; c = next_[c])
      if (equal(members(c), set)) return c;
    const auto id = static_cast<std::uint32_t>(size());
    members_.insert(members_.end(), set.begin(), set.end());
    begin_.push_back(members_.size());
    next_.push_back(it->second);
    it->second = id;
    return id;
  }
  std::size_t memory_bytes() const {
    return begin_.capacity() * sizeof(std::uint64_t) + (members_.capacity() + next_.capacity()) * sizeof(std::uint32_t) + head_.size() * 32;
  }

  static std::uint64_t hash(std::span<const std::uint32_t> s) {
    std::uint64_t h = 0x9e3779b97f4a7c15ULL ^ s.size();
    for (std::uint32_t x : s) h = kmers::mix64(h + x);
    return h;
  }

 private:
  static bool equal(std::span<const std::uint32_t> a, std::span<const std::uint32_t> b) { return std::equal(a.begin(), a.end(), b.begin(), b.end()); }

  std::vector<std::uint64_t> begin_{0};
  std::vector<std::uint32_t> members_, next_;
  std::unordered_map<std::uint64_t, std::uint32_t> head_;
};

class Index {
 public:
  static constexpr std::uint32_t NONE = ClassTable::NONE;

  explicit Index(const Transcriptome& tx, int k = 31) : k_(k) {
    sequence::check_k(k, "expression::Index");
    lengths_ = tx.lengths();
    // The map's values are class + 1, so that 0 is a k-mer not yet seen.
    std::unordered_map<std::uint32_t, std::uint32_t> memo;
    std::vector<std::uint32_t> scratch;
    for (std::uint32_t t = 0; t < tx.size(); t++) {
      memo.clear();
      std::uint32_t alone = NONE;  // the class {t}
      sequence::for_each_kmer(tx.transcripts[t].sequence, k, [&](std::size_t, Kmer f, Kmer r) {
        const Kmer c = f < r ? f : r;
        std::uint32_t& v = map_.at(c, kmers::mix64(c));
        if (v == 0) {
          if (alone == NONE) alone = classes_.insert(std::span<const std::uint32_t>(&t, 1));
          v = alone + 1;
          return;
        }
        const std::uint32_t old = v - 1;
        if (classes_.members(old).back() == t) return;  // already met in t
        const auto hit = memo.find(old);
        if (hit != memo.end()) { v = hit->second + 1; return; }
        const auto m = classes_.members(old);
        scratch.assign(m.begin(), m.end());
        scratch.push_back(t);
        const std::uint32_t id = classes_.insert(scratch);
        memo.emplace(old, id);
        v = id + 1;
      });
    }
  }

  int k() const { return k_; }
  std::size_t transcripts() const { return lengths_.size(); }
  const std::vector<std::size_t>& lengths() const { return lengths_; }
  std::size_t kmers() const { return map_.size(); }
  const ClassTable& classes() const { return classes_; }
  std::size_t memory_bytes() const { return map_.memory_bytes() + classes_.memory_bytes(); }

  // The class of a canonical k-mer, or NONE.
  std::uint32_t class_of(Kmer canonical) const {
    const std::uint32_t* v = map_.find(canonical, kmers::mix64(canonical));
    return v ? *v - 1 : NONE;
  }

  struct Scratch {
    std::vector<std::uint32_t> current, next;
    std::vector<Kmer> kmers;
    std::vector<std::uint64_t> hashes;
  };
  // The class of a read, or of a pair when `mate` is given: a class of the index, or NONE with
  // scratch.current holding a set that is not one (empty when the read is not pseudoaligned).
  // The read's k-mers are hashed and their slots prefetched first, so that the lookups' cache
  // misses overlap rather than follow one another.
  std::uint32_t pseudoalign(const Codes& read, const Codes* mate, Scratch& s) const {
    s.current.clear();
    s.kmers.clear();
    s.hashes.clear();
    const auto collect = [&](std::size_t, Kmer f, Kmer r) {
      const Kmer c = f < r ? f : r;
      const std::uint64_t h = kmers::mix64(c);
      map_.prefetch(h);
      s.kmers.push_back(c);
      s.hashes.push_back(h);
    };
    if (read.size() >= static_cast<std::size_t>(k_)) sequence::for_each_kmer(read, k_, collect);
    if (mate && mate->size() >= static_cast<std::size_t>(k_)) sequence::for_each_kmer(*mate, k_, collect);
    std::uint32_t cls = NONE, last = NONE;
    bool started = false;
    for (std::size_t q = 0; q < s.kmers.size(); q++) {
      const std::uint32_t* v = map_.find(s.kmers[q], s.hashes[q]);
      if (!v) continue;
      const std::uint32_t e = *v - 1;
      if (e == last) continue;
      last = e;
      const auto m = classes_.members(e);
      if (!started) {
        started = true;
        cls = e;
        s.current.assign(m.begin(), m.end());
        continue;
      }
      if (e == cls) continue;
      s.next.clear();
      std::set_intersection(s.current.begin(), s.current.end(), m.begin(), m.end(), std::back_inserter(s.next));
      if (s.next.empty()) {
        s.current.clear();
        return NONE;
      }
      if (s.next.size() == m.size()) cls = e;
      else if (s.next.size() != s.current.size()) cls = NONE;
      std::swap(s.current, s.next);
    }
    if (!started) return NONE;
    if (cls == NONE) cls = classes_.find(s.current);
    return cls;
  }

 private:
  int k_ = 31;
  std::vector<std::size_t> lengths_;
  kmers::KmerMap map_;
  ClassTable classes_;
};

// The pseudoalignment of a sample: the classes met (the index's first, then the new ones in
// the order of the fragments), the number of fragments in each, and each fragment's class.
struct Pseudoalignment {
  ClassTable classes;
  std::vector<std::uint32_t> counts;          // fragments per class
  std::vector<std::uint32_t> fragment_class;  // per fragment, NONE when not pseudoaligned
  std::size_t fragments = 0, aligned = 0;
};

// Pseudoaligns reads (single-end), or pairs when `mates` is not empty (mates[i] the mate of reads[i]).
inline Pseudoalignment pseudoalign(const Index& index, const std::vector<Codes>& reads, const std::vector<Codes>& mates, calx::scheduler::ThreadPool& pool,
                                   std::size_t chunk = 2048) {
  if (!mates.empty() && mates.size() != reads.size()) throw std::invalid_argument("genomics::expression::pseudoalign: one mate per read");
  constexpr std::uint32_t NEW = 0x80000000u;
  const std::size_t n = reads.size();
  chunk = std::max<std::size_t>(1, chunk);
  const std::size_t chunks = (n + chunk - 1) / chunk;
  Pseudoalignment out;
  out.fragments = n;
  out.fragment_class.assign(n, ClassTable::NONE);
  // The sets of each chunk that are not classes of the index, flat.
  std::vector<std::vector<std::uint32_t>> fresh(chunks), fresh_begin(chunks);
  std::vector<Index::Scratch> scratch(pool.size());
  calx::scheduler::parallel_for(pool, chunks, 1, [&](std::size_t b, std::size_t e, std::size_t thread) {
    Index::Scratch& s = scratch[thread];
    for (std::size_t c = b; c < e; c++) {
      fresh[c].clear();
      fresh_begin[c].assign(1, 0);
      for (std::size_t i = c * chunk; i < std::min(n, (c + 1) * chunk); i++) {
        const std::uint32_t cls = index.pseudoalign(reads[i], mates.empty() ? nullptr : &mates[i], s);
        if (cls != ClassTable::NONE) { out.fragment_class[i] = cls; continue; }
        if (s.current.empty()) continue;
        out.fragment_class[i] = NEW | static_cast<std::uint32_t>(fresh_begin[c].size() - 1);
        fresh[c].insert(fresh[c].end(), s.current.begin(), s.current.end());
        fresh_begin[c].push_back(static_cast<std::uint32_t>(fresh[c].size()));
      }
    }
  });
  out.classes = index.classes();
  for (std::size_t c = 0; c < chunks; c++)
    for (std::size_t i = c * chunk; i < std::min(n, (c + 1) * chunk); i++) {
      std::uint32_t& cls = out.fragment_class[i];
      if (cls == ClassTable::NONE || !(cls & NEW)) continue;
      const std::uint32_t local = cls & ~NEW;
      const std::uint32_t* p = fresh[c].data() + fresh_begin[c][local];
      cls = out.classes.insert(std::span<const std::uint32_t>(p, fresh_begin[c][local + 1] - fresh_begin[c][local]));
    }
  out.counts.assign(out.classes.size(), 0);
  for (std::uint32_t cls : out.fragment_class)
    if (cls != ClassTable::NONE) {
      out.counts[cls]++;
      out.aligned++;
    }
  return out;
}

// The pseudoalignment of simulated fragments (pairs when they have second reads).
inline Pseudoalignment pseudoalign(const Index& index, const std::vector<Fragment>& fragments, calx::scheduler::ThreadPool& pool) {
  std::vector<Codes> r1, r2;
  r1.reserve(fragments.size());
  const bool paired = !fragments.empty() && !fragments[0].read2.empty();
  for (const auto& f : fragments) {
    r1.push_back(f.read1);
    if (paired) r2.push_back(f.read2);
  }
  return pseudoalign(index, r1, r2, pool);
}

}  // namespace genomics::expression
