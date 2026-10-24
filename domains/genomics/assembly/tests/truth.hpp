// What the assembly checks compare with: a simulated genome, indexed by its k-mers on both
// strands, so that a contig can be placed on it exactly; and reads simulated in shards.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace truth {

using genomics::sequence::Codes;
using genomics::sequence::Kmer;

// Where a sequence lies on the genome: forward coordinates [begin, end), and the strand.
struct Place {
  std::size_t begin = 0, end = 0;
  bool reverse = false;
};

class GenomeIndex {
 public:
  GenomeIndex(const Codes& genome, int k) : k_(k), forward_(genome), reverse_(genomics::sequence::reverse_complement(genome)) {
    index(forward_, forward_at_);
    index(reverse_, reverse_at_);
  }
  // Every place where s occurs exactly, on either strand.
  std::vector<Place> places(const Codes& s) const {
    std::vector<Place> out;
    if (s.size() < static_cast<std::size_t>(k_)) return out;
    const Kmer x = genomics::sequence::encode_kmer(s.data(), k_);
    const auto scan = [&](const Codes& g, const std::unordered_map<Kmer, std::vector<std::uint32_t>>& at, bool rev) {
      const auto it = at.find(x);
      if (it == at.end()) return;
      for (auto p : it->second)
        if (p + s.size() <= g.size() && std::equal(s.begin(), s.end(), g.begin() + p)) {
          const std::size_t b = rev ? g.size() - p - s.size() : p;
          out.push_back({b, b + s.size(), rev});
        }
    };
    scan(forward_, forward_at_, false);
    scan(reverse_, reverse_at_, true);
    return out;
  }
  bool contains(const Codes& s) const { return !places(s).empty(); }
  // Whether the canonical k-mer x occurs in the genome.
  bool has(Kmer x) const {
    return forward_at_.count(x) || forward_at_.count(genomics::sequence::reverse_complement(x, k_));
  }
  // The fraction of the k-mers of s that occur in the genome.
  double kmer_fraction(const Codes& s) const {
    std::size_t in = 0, all = 0;
    genomics::sequence::for_each_kmer(s, k_, [&](std::size_t, Kmer f, Kmer) { all++; in += forward_at_.count(f) || reverse_at_.count(f); });
    return all ? static_cast<double>(in) / static_cast<double>(all) : 1;
  }
  const Codes& genome() const { return forward_; }

 private:
  void index(const Codes& g, std::unordered_map<Kmer, std::vector<std::uint32_t>>& at) {
    genomics::sequence::for_each_kmer(g, k_, [&](std::size_t p, Kmer f, Kmer) { at[f].push_back(static_cast<std::uint32_t>(p)); });
  }
  int k_;
  Codes forward_, reverse_;
  std::unordered_map<Kmer, std::vector<std::uint32_t>> forward_at_, reverse_at_;
};

// The reads of a sample, their bases only, simulated in shards on the pool (each read has its
// own stream, so the sample is the one short_reads makes alone).
inline std::vector<Codes> reads(const genomics::simulate::Source& source, std::size_t n, const genomics::simulate::ShortReadModel& model, std::uint64_t seed,
                                calx::scheduler::ThreadPool& pool) {
  const std::size_t shards = std::max<std::size_t>(1, n / 2000);
  std::vector<std::vector<Codes>> part(shards);
  calx::scheduler::parallel_for(pool, shards, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t q = b; q < e; q++) {
      const auto r = calx::scheduler::part(n, shards, q);
      for (auto& read : genomics::simulate::short_reads(source, r.size(), model, seed, r.begin)) part[q].push_back(std::move(read.bases));
    }
  });
  std::vector<Codes> out;
  out.reserve(n);
  for (auto& p : part)
    for (auto& r : p) out.push_back(std::move(r));
  return out;
}

}  // namespace truth
