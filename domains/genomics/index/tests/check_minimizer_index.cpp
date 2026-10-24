// Checks of the minimizer index: its occurrences, built from chunks with margins, equal the
// minimizers of each whole contig (contigs shorter than a window, runs of N, chunks of 997 bases);
// the index is the same on one and eight threads; lookups return exactly each hash's occurrences;
// the cut-off of repetitive minimizers is minimap2's quantile plus one; anchors equal their
// definition by brute force; reads drawn from the genome on either strand anchor on their true
// diagonal; and, when the data library has it, the index of Escherichia coli with its time and size.

#include <genomics/index/minimizer_index.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/sequence/kmer.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace genomics;
using index::Anchor;
using index::MinimizerIndex;
using sequence::Codes;

namespace {

constexpr std::uint32_t CHECK = 0x4d494e49;

std::string read_file(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// hash → occurrences (contig, position, strand), from the minimizers of each whole contig.
std::map<std::uint64_t, std::vector<std::tuple<std::uint32_t, std::uint32_t, bool>>> brute_index(const std::vector<Codes>& contigs, int k, int w) {
  std::map<std::uint64_t, std::vector<std::tuple<std::uint32_t, std::uint32_t, bool>>> m;
  for (std::uint32_t i = 0; i < contigs.size(); i++)
    for (const auto& x : sequence::minimizers(contigs[i], k, w)) m[x.hash].push_back({i, x.position, x.reverse});
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, CHECK, 0, 0);
  calx::scheduler::ThreadPool one(1), eight(8);

  // A genome with a family of repeats and runs of N, cut into contigs of very different sizes.
  simulate::GenomeSpec spec;
  spec.length = 200000;
  spec.repeats = {4, 40, 400, 0.01};
  const Codes whole = simulate::make_genome(spec, 3).sequence;
  std::vector<Codes> contigs;
  for (auto [b, e] : {std::pair<std::size_t, std::size_t>{0, 120000}, {120000, 120005}, {120005, 120030}, {120030, 170000}, {170000, 200000}})
    contigs.emplace_back(whole.begin() + static_cast<std::ptrdiff_t>(b), whole.begin() + static_cast<std::ptrdiff_t>(e));
  for (std::size_t p = 5000; p < 5020; p++) contigs[0][p] = sequence::N;
  for (std::size_t p = 999; p < 1001; p++) contigs[0][p] = sequence::N;
  for (std::size_t p = 0; p < 50000; p += 1009) contigs[3][p] = sequence::N;
  contigs[2][12] = sequence::N;

  // ---- construction ----------------------------------------------------------------------------------
  {
    bool chunked = true, threads = true;
    for (auto [k, w] : {std::pair{15, 10}, std::pair{21, 11}, std::pair{5, 3}, std::pair{19, 1}}) {
      index::MinimizerOptions o;
      o.k = k;
      o.w = w;
      o.chunk = 997;
      const MinimizerIndex a(one, contigs, o), b(eight, contigs, o);
      const auto brute = brute_index(contigs, k, w);
      std::size_t total = 0;
      bool same = a.distinct() == brute.size();
      for (const auto& [h, occ] : brute) {
        const auto got = a.lookup(h);
        total += occ.size();
        same = same && got.size() == occ.size();
        for (std::size_t i = 0; same && i < occ.size(); i++) {
          const auto o2 = MinimizerIndex::decode(got[i]);
          same = same && std::tuple{o2.contig, o2.position, o2.reverse} == occ[i];
        }
      }
      chunked = chunked && same && a.occurrences() == total;
      threads = threads && a.entries() == b.entries() && a.distinct() == b.distinct() && a.max_occurrences() == b.max_occurrences();
      o.chunk = 1 << 20;
      const MinimizerIndex c(eight, contigs, o);
      threads = threads && c.entries() == a.entries();
    }
    s.require("the index built from chunks of 997 bases with margins holds exactly each whole contig's minimizers, and each lookup its occurrences", chunked);
    s.require("the index is the same on one and eight threads, and for any chunk length", threads);
    bool absent = true;
    index::MinimizerOptions o;
    const MinimizerIndex a(eight, contigs, o);
    const auto brute = brute_index(contigs, o.k, o.w);
    for (int q = 0; q < 10000; q++) {
      const std::uint64_t h = r.below(1ULL << 30);
      if (!brute.count(h)) absent = absent && a.lookup(h).empty();
    }
    s.require("a hash that is no minimizer of the genome finds nothing", absent);

    // The cut-off: the ⌊(1 − f) n⌋-th smallest count plus one (minimap2's mm_idx_cal_max_occ).
    for (double f : {2e-4, 0.01, 0.1}) {
      o.ignore = f;
      const MinimizerIndex x(eight, contigs, o);
      std::vector<std::uint32_t> counts;
      for (const auto& [h, occ] : brute) counts.push_back(static_cast<std::uint32_t>(occ.size()));
      std::sort(counts.begin(), counts.end());
      const auto q = static_cast<std::size_t>((1 - f) * static_cast<double>(counts.size()));
      const std::uint32_t want = counts[std::min(q, counts.size() - 1)] + 1;
      std::size_t ignored = 0;
      for (auto c : counts) ignored += c > want;
      s.require("the cut-off at f = " + std::to_string(f) + " is the (1 − f) quantile of counts plus one", x.max_occurrences() == want,
                "at most " + std::to_string(want) + " occurrences; " + std::to_string(ignored) + " of " + std::to_string(counts.size()) + " minimizers ignored");
    }
  }

  // ---- anchors ----------------------------------------------------------------------------------
  {
    index::MinimizerOptions o;
    o.ignore = 0.01;
    const MinimizerIndex idx(eight, contigs, o);
    const auto brute = brute_index(contigs, o.k, o.w);
    bool definition = true, diagonal = true;
    std::size_t reads = 0, on_diagonal = 0;
    for (int q = 0; q < 400; q++) {
      const std::uint32_t c = q % 2 ? 0 : 4;
      const Codes& g = contigs[c];
      const std::size_t len = 300, b = r.below(g.size() - len);
      Codes read(g.begin() + static_cast<std::ptrdiff_t>(b), g.begin() + static_cast<std::ptrdiff_t>(b + len));
      const bool rev = q % 4 >= 2;
      if (rev) read = sequence::reverse_complement(read);
      if (q % 5 == 0)
        for (auto& x : read)
          if (r.uniform() < 0.02) x = static_cast<sequence::Code>((x + 1 + r.below(3)) % 4);
      const auto got = idx.anchors(read);
      std::vector<Anchor> want;
      for (const auto& m : sequence::minimizers(read, o.k, o.w)) {
        const auto it = brute.find(m.hash);
        if (it == brute.end() || it->second.size() > idx.max_occurrences()) continue;
        for (const auto& [ct, pos, strand] : it->second) want.push_back({ct, pos, m.position, strand != m.reverse});
      }
      std::sort(want.begin(), want.end());
      definition = definition && got == want && std::is_sorted(got.begin(), got.end());
      // On the true diagonal: forward, t − q = b; reverse, t + q + k = b + len.
      if (q % 5 == 0) continue;
      reads++;
      std::size_t on = 0;
      for (const auto& a : got)
        if (a.contig == c && a.reverse == rev &&
            (rev ? static_cast<std::size_t>(a.target) + a.query + static_cast<std::size_t>(o.k) == b + len : static_cast<std::size_t>(a.target) == b + a.query))
          on++;
      on_diagonal += on > 0;
      diagonal = diagonal && on > 0;
    }
    s.require("anchors equal their definition: every query minimizer against every occurrence, the repetitive ignored", definition);
    s.require("every exact read of 300 bases, on either strand, has anchors on its true diagonal", diagonal, std::to_string(on_diagonal) + " of " + std::to_string(reads));
  }

  // ---- Escherichia coli ----------------------------------------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
    const auto file = cache / "ncbi-genomes" / "NC_000913.3.gb";
    if (!std::filesystem::exists(file)) {
      s.skip("Escherichia coli NC_000913.3 is not in the data cache (data/cache/ncbi-genomes)");
    } else {
      const std::vector<Codes> ecoli{sequence::encode_dna(sequence::read_genbank(read_file(file)).at(0).sequence)};
      const auto time = [&](calx::scheduler::ThreadPool& pool, MinimizerIndex& out) {
        const auto t0 = std::chrono::steady_clock::now();
        out = MinimizerIndex(pool, ecoli, {15, 10, 2e-4, 1 << 18});
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      };
      MinimizerIndex a, b;
      const double t1 = time(one, a), t8 = time(eight, b);
      char buf[200];
      std::snprintf(buf, sizeof buf, "%zu occurrences of %zu minimizers, %.1f MB; %.3f s on one thread, %.3f s on eight; cut-off %u", a.occurrences(), a.distinct(),
                    static_cast<double>(a.bytes()) / 1e6, t1, t8, a.max_occurrences());
      s.require("Escherichia coli's (10, 15)-minimizer index, the same on one and eight threads", a.entries() == b.entries(), buf);
      const double density = static_cast<double>(a.occurrences()) / static_cast<double>(ecoli[0].size());
      s.check("its density of minimizers against 2/(w + 1), within the excess of hash64 at k = 15 (kmer.hpp)", density / (2.0 / 11) - 1, 0.05);
    }
  }

  return s.done();
}
