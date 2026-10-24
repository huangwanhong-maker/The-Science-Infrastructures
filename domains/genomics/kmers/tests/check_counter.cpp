// Checks of kmers/counter.hpp: canonical k-mer counts of simulated reads (with errors and Ns)
// against a std::unordered_map filled one k-mer at a time, for k = 1, 15, 21, 31 and 32; the
// same tables, in the same order and the same memory, on 1, 2 and 8 threads; long sequences
// cut into chunks that overlap by k − 1, against the same reference; any number of partitions
// and any cut into calls of add giving the same counts; the canonical order of sorted(); and
// the memory against its statement.

#include <genomics/kmers/counter.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

using namespace genomics;
using kmers::Kmer;
using sequence::Codes;

namespace {

std::unordered_map<Kmer, std::uint32_t> reference_counts(const std::vector<Codes>& seqs, int k) {
  std::unordered_map<Kmer, std::uint32_t> m;
  for (const auto& s : seqs)
    for (std::size_t i = 0; i + static_cast<std::size_t>(k) <= s.size(); i++) {
      bool clean = true;
      for (std::size_t j = i; j < i + static_cast<std::size_t>(k); j++) clean = clean && s[j] < 4;
      if (!clean) continue;
      m[sequence::canonical(sequence::encode_kmer(&s[i], k), k)]++;
    }
  return m;
}

bool same_as(const kmers::KmerCounter& c, const std::unordered_map<Kmer, std::uint32_t>& ref) {
  if (c.distinct() != ref.size()) return false;
  std::uint64_t total = 0;
  for (const auto& [x, n] : ref) {
    total += n;
    if (c.count(x) != n) return false;
  }
  bool ok = total == c.total();
  c.for_each([&](Kmer x, std::uint32_t n) { auto it = ref.find(x); ok = ok && it != ref.end() && it->second == n; });
  return ok;
}

// The slot order of every table, as one sequence.
std::vector<kmers::KmerCount> slot_order(const kmers::KmerCounter& c) {
  std::vector<kmers::KmerCount> v;
  c.for_each([&](Kmer x, std::uint32_t n) { v.push_back({x, n}); });
  return v;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // A genome with repeats, and reads with errors from it, a few Ns sprinkled in.
  simulate::GenomeSpec spec;
  spec.length = 60'000;
  spec.repeats = {3, 4, 400, 0.01};
  const auto genome = simulate::make_genome(spec, 7);
  const auto source = simulate::Source::haploid({genome.sequence});
  auto reads = simulate::short_reads(source, 4000, {}, 11);
  calx::random::Stream r(5, 0x4e4e4e4e, 0, 0);
  for (auto& read : reads)
    for (auto& b : read.bases)
      if (r.uniform() < 0.002) b = sequence::N;
  std::vector<Codes> seqs;
  for (const auto& read : reads) seqs.push_back(read.bases);

  calx::scheduler::ThreadPool one(1), two(2), eight(8);

  // Against the reference map, for several k.
  for (int k : {1, 15, 21, 31, 32}) {
    kmers::CounterOptions o;
    o.k = k;
    o.chunk_bases = 20'000;
    o.batch_bases = 150'000;
    const auto c = kmers::count_kmers(seqs, eight, o);
    const auto ref = reference_counts(seqs, k);
    s.require("counter: canonical counts equal a hash map's, k = " + std::to_string(k) + ", 4000 reads with errors and Ns, 8 threads",
              same_as(c, ref), std::to_string(c.distinct()) + " distinct, " + std::to_string(c.total()) + " in all");
  }

  // The same tables on 1, 2 and 8 threads.
  {
    kmers::CounterOptions o;
    o.chunk_bases = 30'000;
    o.batch_bases = 200'000;
    const auto c1 = kmers::count_kmers(seqs, one, o), c2 = kmers::count_kmers(seqs, two, o), c8 = kmers::count_kmers(seqs, eight, o);
    const auto o1 = slot_order(c1);
    s.require("counter: the same counts, in the same slot order, on 1, 2 and 8 threads", o1 == slot_order(c2) && o1 == slot_order(c8), std::to_string(o1.size()) + " entries");
    s.require("counter: the same memory on 1, 2 and 8 threads", c1.memory_bytes() == c2.memory_bytes() && c1.memory_bytes() == c8.memory_bytes());
    s.require("counter: sorted() the same on 1 and 8 threads", c1.sorted(1, one) == c8.sorted(1, eight));
  }

  // Long sequences cut into chunks overlapping by k − 1: the genome itself, with Ns, and its pieces.
  {
    Codes g = genome.sequence;
    for (std::size_t i = 1000; i < g.size(); i += 7919) g[i] = sequence::N;
    const std::vector<Codes> contigs = {g, Codes(g.begin(), g.begin() + 5000), Codes(10, sequence::A)};
    bool ok = true;
    for (std::size_t chunk : {1u, 31u, 97u, 1000u, 100000u}) {
      kmers::CounterOptions o;
      o.chunk_bases = chunk;
      o.batch_bases = std::max<std::size_t>(7 * chunk + 3, 4000);
      ok = ok && same_as(kmers::count_kmers(contigs, eight, o), reference_counts(contigs, 31));
    }
    s.require("counter: a 60 kb sequence with Ns cut into chunks of 1 to 10^5 bases overlapping by k − 1: counts unchanged", ok);
  }

  // Partitions and calls of add.
  {
    const auto want = kmers::count_kmers(seqs, eight).sorted(1, eight);
    bool ok = true;
    for (unsigned bits : {0u, 3u, 8u, 12u}) {
      kmers::CounterOptions o;
      o.partition_bits = bits;
      ok = ok && kmers::count_kmers(seqs, eight, o).sorted(1, eight) == want;
    }
    s.require("counter: 1, 8, 256 and 4096 partitions give the sorted counts of the default 1024", ok);
    kmers::KmerCounter c;
    const std::size_t half = seqs.size() / 2;
    c.add(half, [&](std::size_t i) -> const Codes& { return seqs[i]; }, two);
    c.add(seqs.size() - half, [&](std::size_t i) -> const Codes& { return seqs[half + i]; }, eight);
    s.require("counter: counts added in two calls equal counts added in one", c.sorted(1, eight) == want);
    bool increasing = true;
    for (std::size_t i = 1; i < want.size(); i++) increasing = increasing && want[i - 1].kmer < want[i].kmer;
    bool canonical = true;
    for (const auto& e : want) canonical = canonical && e.kmer == sequence::canonical(e.kmer, 31) && c.count(sequence::reverse_complement(e.kmer, 31)) == e.count;
    s.require("counter: sorted() in increasing order of canonical k-mers, each found in either orientation", increasing && canonical);
    std::size_t solid = 0;
    for (const auto& e : want) solid += e.count >= 3;
    s.require("counter: sorted(3) keeps exactly the k-mers counted three times or more", c.sorted(3, eight).size() == solid, std::to_string(solid) + " solid");
    {
      const auto ref = reference_counts(seqs, 31);
      std::size_t absent = 0, wrong = 0;
      for (std::uint32_t t = 0; absent < 1000; t++) {
        calx::random::Stream q(9, 0x41425345, t, 0);
        const Kmer x = sequence::canonical(q.next64() & sequence::kmer_mask(31), 31);
        if (ref.count(x)) continue;
        absent++;
        wrong += c.count(x) != 0;
      }
      s.require("counter: 1000 random k-mers absent from the reads count 0", wrong == 0);
    }
  }

  // Memory: 12 bytes a slot at a load between 5/16 and 5/8, plus the scatter buffers.
  {
    auto c = kmers::count_kmers(seqs, eight);
    std::size_t slots = 0;
    for (std::size_t p = 0; p < c.partitions(); p++) slots += c.table(p).capacity();
    c.release_buffers();
    const double per = static_cast<double>(c.memory_bytes()) / static_cast<double>(c.distinct());
    s.require("counter: the tables hold 12 bytes a slot, and every table below a load of 5/8",
              c.memory_bytes() == 12 * slots && [&] {
                for (std::size_t p = 0; p < c.partitions(); p++)
                  if (c.table(p).size() * 8 > c.table(p).capacity() * 5) return false;
                return true;
              }(),
              calx::verification::exponential(per) + " bytes for each distinct k-mer");
  }
  return s.done();
}
