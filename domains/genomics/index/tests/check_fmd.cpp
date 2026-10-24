// Checks of the FMD-index: its BWT and occurrence counts against a BWT built independently (prefix
// doubling on the same text); counts and positions of patterns on both strands against brute
// force, across contigs and runs of N; bidirectional extension (forwards, backwards and mixed)
// against backward search of the pattern and of its reverse complement; super-maximal exact
// matches against their definition by brute force; SMEMs of many queries the same on one and four
// threads; and, when the data library has fetched it, Escherichia coli K-12 MG1655 (NC_000913.3):
// the time and memory of the build, every sampled 32-mer found at its true place, and the SMEMs
// of simulated reads.

#include <genomics/index/fmd.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/simulate/genome.hpp>

#include <morphe/formal/suffix_array.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace genomics;
using index::BiInterval;
using index::FmdIndex;
using index::Hit;
using index::Smem;
using sequence::Codes;

namespace {

constexpr std::uint32_t CHECK = 0x464d4443;  // the purpose of the checks' own draws

// The contigs as the index sees them: each N replaced, in order, by a base drawn from the
// stream (seed, FMD_N, contig, 0), as fmd.hpp documents.
std::vector<Codes> as_indexed(const std::vector<Codes>& contigs, std::uint64_t seed) {
  std::vector<Codes> out = contigs;
  for (std::uint32_t i = 0; i < out.size(); i++) {
    calx::random::Stream r(seed, index::FMD_N, i, 0);
    for (auto& c : out[i])
      if (c > 3) c = static_cast<sequence::Code>(r.below(4));
  }
  return out;
}

// Every occurrence of p on both strands, as hits.
std::vector<Hit> brute_hits(const std::vector<Codes>& genome, const Codes& p) {
  std::vector<Hit> out;
  if (p.empty()) return out;
  for (std::uint32_t i = 0; i < genome.size(); i++) {
    const Codes& f = genome[i];
    const Codes rc = sequence::reverse_complement(f);
    if (p.size() > f.size()) continue;
    for (std::size_t o = 0; o + p.size() <= f.size(); o++) {
      if (std::equal(p.begin(), p.end(), f.begin() + static_cast<std::ptrdiff_t>(o))) out.push_back({i, o, false});
      if (std::equal(p.begin(), p.end(), rc.begin() + static_cast<std::ptrdiff_t>(o))) out.push_back({i, f.size() - o - p.size(), true});
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

// The longest match in the genome (either strand) starting at each query position.
std::vector<std::size_t> longest_from(const std::vector<Codes>& genome, const Codes& q, std::uint64_t min_occ) {
  std::vector<Codes> strands;
  for (const auto& c : genome) {
    strands.push_back(c);
    strands.push_back(sequence::reverse_complement(c));
  }
  std::vector<std::size_t> best(q.size(), 0);
  for (std::size_t b = 0; b < q.size(); b++) {
    if (q[b] > 3) continue;
    // The lengths of the match at every text position; the min_occ-th largest of them.
    std::vector<std::size_t> top(min_occ, 0);
    for (const auto& s : strands)
      for (std::size_t t = 0; t < s.size(); t++) {
        std::size_t h = 0;
        while (b + h < q.size() && t + h < s.size() && q[b + h] < 4 && q[b + h] == s[t + h]) h++;
        if (h > top.back()) {
          top.back() = h;
          for (std::size_t j = top.size() - 1; j > 0 && top[j] > top[j - 1]; j--) std::swap(top[j], top[j - 1]);
        }
      }
    best[b] = top.back();
  }
  return best;
}

// The SMEMs of a query by definition: maximal matches [b, b + L(b)) not extendable to the left
// (L(b − 1) ≤ L(b)), and contained in no other.
std::vector<std::pair<std::uint32_t, std::uint32_t>> brute_smems(const std::vector<Codes>& genome, const Codes& q, std::uint64_t min_occ) {
  const auto L = longest_from(genome, q, min_occ);
  std::vector<std::pair<std::uint32_t, std::uint32_t>> mems, out;
  for (std::size_t b = 0; b < q.size(); b++)
    if (L[b] > 0 && (b == 0 || L[b - 1] <= L[b])) mems.push_back({static_cast<std::uint32_t>(b), static_cast<std::uint32_t>(b + L[b])});
  for (const auto& m : mems) {
    bool contained = false;
    for (const auto& o : mems)
      if (o != m && o.first <= m.first && m.second <= o.second) contained = true;
    if (!contained) out.push_back(m);
  }
  return out;
}

std::string read_file(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

double peak_megabytes() {
#ifdef __linux__
  std::ifstream in("/proc/self/status");
  std::string line;
  while (std::getline(in, line))
    if (line.rfind("VmHWM:", 0) == 0) return std::atof(line.c_str() + 6) / 1024;
#endif
  return 0;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, CHECK, 0, 0);

  // A genome of five contigs with repeats and runs of N.
  simulate::GenomeSpec spec;
  spec.length = 24000;
  spec.repeats = {3, 6, 300, 0.02};
  spec.tandem = {2, 7, 12, 0.0};
  const Codes whole = simulate::make_genome(spec, 7).sequence;
  std::vector<Codes> contigs;
  for (std::size_t b : {0, 9000, 15000, 15040, 21000}) {
    const std::size_t e = b == 0 ? 9000 : b == 9000 ? 15000 : b == 15000 ? 15040 : b == 15040 ? 21000 : whole.size();
    contigs.emplace_back(whole.begin() + static_cast<std::ptrdiff_t>(b), whole.begin() + static_cast<std::ptrdiff_t>(e));
  }
  for (std::size_t p = 100; p < 130; p++) contigs[0][p] = sequence::N;
  for (std::size_t p = 4000; p < 4003; p++) contigs[3][p] = sequence::N;
  contigs[4][0] = sequence::N;
  const index::FmdOptions options{8, 11};
  const FmdIndex idx(contigs, options);
  const auto genome = as_indexed(contigs, options.seed);

  // ---- the BWT against one built independently ----------------------------------------------------------------------
  {
    std::vector<std::uint32_t> text;
    const std::uint32_t seqs = static_cast<std::uint32_t>(2 * genome.size());
    for (std::uint32_t i = 0; i < genome.size(); i++) {
      for (auto c : genome[i]) text.push_back(seqs + c);
      text.push_back(2 * i);
      for (auto c : sequence::reverse_complement(genome[i])) text.push_back(seqs + c);
      text.push_back(2 * i + 1);
    }
    const auto sa = morphe::formal::suffix_array_doubling(text, seqs + 4);
    bool symbols = true, occ = true, occ4 = true;
    std::array<std::uint64_t, 4> count{};
    for (std::uint64_t i = 0; i <= text.size(); i++) {
      const auto o4 = idx.occ4(i);
      for (int c = 0; c < 4; c++) {
        occ = occ && idx.occ(static_cast<sequence::Code>(c), i) == count[static_cast<std::size_t>(c)];
        occ4 = occ4 && o4[static_cast<std::size_t>(c)] == count[static_cast<std::size_t>(c)];
      }
      if (i == text.size()) break;
      const std::uint32_t prior = text[sa[i] == 0 ? text.size() - 1 : sa[i] - 1];
      const sequence::Code want = prior < seqs ? FmdIndex::DOLLAR : static_cast<sequence::Code>(prior - seqs);
      symbols = symbols && idx.bwt(i) == want;
      if (want < 4) count[want]++;
    }
    s.require("the FMD-index's BWT over {$, A, C, G, T} equals one from an independent suffix array (prefix doubling)", symbols && idx.size() == text.size(),
              std::to_string(text.size()) + " rows");
    s.require("Occ(c, i) from the interleaved blocks equals counting, at every row", occ);
    s.require("Occ of the four bases at once equals counting", occ4);
    s.require("the text is closed under reverse complementation: Occ(A) = Occ(T) and Occ(C) = Occ(G) at its end",
              count[0] == count[3] && count[1] == count[2]);
    std::vector<index::AmbiguousRun> runs{{0, 100, 130}, {3, 4000, 4003}, {4, 0, 1}};
    s.require("the runs of N are recorded where they are", idx.ambiguous() == runs);
  }

  // ---- count and locate ----------------------------------------------------------------------------------
  {
    bool counts = true, hits = true, with_n = true, flagged = true;
    std::size_t patterns = 0, occurrences = 0;
    for (int q = 0; q < 3000; q++) {
      Codes p;
      const std::uint32_t c = static_cast<std::uint32_t>(r.below(genome.size()));
      const Codes& g = genome[c];
      if (q % 3 != 2 && g.size() > 50) {
        const std::size_t len = 1 + r.below(40), b = r.below(g.size() - len);
        p.assign(g.begin() + static_cast<std::ptrdiff_t>(b), g.begin() + static_cast<std::ptrdiff_t>(b + len));
        if (q % 3 == 1) p = sequence::reverse_complement(p);
      } else {
        p.resize(1 + r.below(10));
        for (auto& x : p) x = static_cast<sequence::Code>(r.below(4));
      }
      const auto want = brute_hits(genome, p);
      counts = counts && idx.count(p) == want.size();
      hits = hits && idx.locate(p) == want;
      patterns++;
      occurrences += want.size();
      if (q % 50 == 0 && p.size() > 2) {
        Codes pn = p;
        pn[p.size() / 2] = sequence::N;
        with_n = with_n && idx.count(pn) == 0;
      }
    }
    // A pattern across the run of N in contig 0 is found on the drawn bases, and flagged.
    const Codes across(genome[0].begin() + 95, genome[0].begin() + 135);
    const auto h = idx.locate(across);
    flagged = !h.empty() && std::all_of(h.begin(), h.end(), [&](const Hit& x) { return idx.overlaps_ambiguous(x, across.size()); }) &&
              !idx.overlaps_ambiguous({0, 200, false}, 40) && idx.overlaps_ambiguous({0, 129, false}, 1);
    s.require("counts of patterns on both strands equal brute force", counts, std::to_string(patterns) + " patterns, " + std::to_string(occurrences) + " occurrences");
    s.require("positions (contig, position, strand) by the sampled suffix array equal brute force", hits);
    s.require("a pattern holding N occurs nowhere", with_n);
    s.require("a match on the bases drawn for N is found and flagged as overlapping the run", flagged);
  }

  // ---- bidirectional extension ----------------------------------------------------------------------------------
  {
    bool forward = true, backward = true, mixed = true;
    for (int q = 0; q < 2000; q++) {
      const Codes& g = genome[r.below(genome.size())];
      if (g.size() < 60) continue;
      const std::size_t len = 1 + r.below(50), b = r.below(g.size() - len);
      Codes p(g.begin() + static_cast<std::ptrdiff_t>(b), g.begin() + static_cast<std::ptrdiff_t>(b + len));
      if (q % 2) p = sequence::reverse_complement(p);
      const BiInterval want = idx.search(p);
      const BiInterval rc = idx.search(sequence::reverse_complement(p));
      backward = backward && want.s > 0 && want.l == rc.k && rc.l == want.k && rc.s == want.s;
      BiInterval x = idx.init(p[0]);
      for (std::size_t i = 1; i < p.size(); i++) x = idx.forward(x, p[i]);
      forward = forward && x == want;
      // From a middle base, forwards to the end, then backwards to the start.
      const std::size_t m = r.below(p.size());
      BiInterval y = idx.init(p[m]);
      for (std::size_t i = m + 1; i < p.size(); i++) y = idx.forward(y, p[i]);
      for (std::size_t i = m; i-- > 0;) y = idx.backward(y, p[i]);
      mixed = mixed && y == want;
    }
    s.require("the bi-interval of P holds the rows of P̄: l(P) = k(P̄) and k(P) = l(P̄)", backward);
    s.require("forward extension letter by letter gives the bi-interval of backward search", forward);
    s.require("extension forwards then backwards from any base gives the same bi-interval", mixed);
  }

  // ---- super-maximal exact matches ----------------------------------------------------------------------------------
  {
    bool exact = true, counted = true, filtered = true, repeated = true;
    std::size_t smems = 0;
    std::vector<Codes> queries;
    for (int q = 0; q < 120; q++) {
      Codes x;
      const Codes& g = genome[r.below(genome.size())];
      if (q % 6 == 5 || g.size() < 200) {
        x.resize(80);
        for (auto& c : x) c = static_cast<sequence::Code>(r.below(4));
      } else {
        const std::size_t b = r.below(g.size() - 150);
        x.assign(g.begin() + static_cast<std::ptrdiff_t>(b), g.begin() + static_cast<std::ptrdiff_t>(b + 150));
        for (auto& c : x)
          if (r.uniform() < 0.03) c = static_cast<sequence::Code>((c + 1 + r.below(3)) % 4);
        if (q % 4 == 0) x[r.below(x.size())] = sequence::N;
        if (q % 2) x = sequence::reverse_complement(x);
      }
      queries.push_back(x);
      const auto got = idx.smems(x);
      std::vector<std::pair<std::uint32_t, std::uint32_t>> g2;
      for (const auto& m : got) {
        g2.push_back({m.begin, m.end});
        const Codes f(x.begin() + m.begin, x.begin() + m.end);
        counted = counted && m.interval.s == brute_hits(genome, f).size();
      }
      exact = exact && g2 == brute_smems(genome, x, 1);
      smems += got.size();
      const auto long_ones = idx.smems(x, 19);
      std::vector<Smem> want;
      for (const auto& m : got)
        if (m.length() >= 19) want.push_back(m);
      filtered = filtered && long_ones == want;
      if (q < 30) {
        std::vector<std::pair<std::uint32_t, std::uint32_t>> g3;
        for (const auto& m : idx.smems(x, 1, 2)) g3.push_back({m.begin, m.end});
        repeated = repeated && g3 == brute_smems(genome, x, 2);
      }
    }
    s.require("SMEMs (Li 2012, Algorithm 5) equal their definition by brute force on mutated reads, both strands, N", exact,
              std::to_string(smems) + " SMEMs in " + std::to_string(queries.size()) + " queries");
    s.require("each SMEM's bi-interval counts its occurrences on both strands", counted);
    s.require("a minimum length keeps exactly the SMEMs that long", filtered);
    s.require("SMEMs occurring at least twice equal their definition by brute force", repeated);
    calx::scheduler::ThreadPool one(1), four(4);
    s.require("the SMEMs of many queries are the same on one and four threads", idx.smems(one, queries, 10) == idx.smems(four, queries, 10));
  }

  // ---- Escherichia coli ----------------------------------------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
    const auto file = cache / "ncbi-genomes" / "NC_000913.3.gb";
    if (!std::filesystem::exists(file)) {
      s.skip("Escherichia coli NC_000913.3 is not in the data cache (data/cache/ncbi-genomes)");
    } else {
      const auto recs = sequence::read_genbank(read_file(file));
      const std::vector<Codes> ecoli{sequence::encode_dna(recs.at(0).sequence)};
      const auto t0 = std::chrono::steady_clock::now();
      const FmdIndex e(ecoli);
      const double build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      char buf[200];
      std::snprintf(buf, sizeof buf, "%zu bases, %.2f s, index %.1f MB, peak resident %.0f MB", ecoli[0].size(), build, static_cast<double>(e.bytes()) / 1e6,
                    peak_megabytes());
      s.require("Escherichia coli K-12 MG1655 is indexed (both strands)", e.size() == 2 * (ecoli[0].size() + 1) && build < 30, buf);

      bool found = true;
      const auto t1 = std::chrono::steady_clock::now();
      for (int q = 0; q < 1000; q++) {
        const std::size_t b = r.below(ecoli[0].size() - 32);
        Codes p(ecoli[0].begin() + static_cast<std::ptrdiff_t>(b), ecoli[0].begin() + static_cast<std::ptrdiff_t>(b + 32));
        const bool rev = q % 2;
        if (rev) p = sequence::reverse_complement(p);
        const auto h = e.locate(p);
        found = found && std::find(h.begin(), h.end(), Hit{0, b, rev}) != h.end();
      }
      const double locate = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
      std::snprintf(buf, sizeof buf, "%.1f µs a 32-mer located", locate / 1000 * 1e6);
      s.require("1000 32-mers sampled from both strands are located at their true places", found, buf);

      std::vector<Codes> reads;
      std::vector<std::pair<std::size_t, bool>> truth;
      for (int q = 0; q < 4000; q++) {
        const std::size_t b = r.below(ecoli[0].size() - 150);
        Codes x(ecoli[0].begin() + static_cast<std::ptrdiff_t>(b), ecoli[0].begin() + static_cast<std::ptrdiff_t>(b + 150));
        for (auto& c : x)
          if (r.uniform() < 0.01) c = static_cast<sequence::Code>((c + 1 + r.below(3)) % 4);
        const bool rev = q % 2;
        if (rev) x = sequence::reverse_complement(x);
        reads.push_back(x);
        truth.push_back({b, rev});
      }
      calx::scheduler::ThreadPool pool(8);
      const auto t2 = std::chrono::steady_clock::now();
      const auto seeds = e.smems(pool, reads, 19);
      const double smem = std::chrono::duration<double>(std::chrono::steady_clock::now() - t2).count();
      std::size_t placed = 0, total = 0;
      for (std::size_t q = 0; q < reads.size(); q++) {
        total += seeds[q].size();
        // The longest SMEM placed within a few bases of the read's true diagonal.
        const Smem* best = nullptr;
        for (const auto& m : seeds[q])
          if (!best || m.length() > best->length()) best = &m;
        if (!best || best->interval.s > 20) continue;
        for (const Hit& h : e.locate(best->interval, best->length())) {
          const auto [b, rev] = truth[q];
          const std::int64_t start = rev ? static_cast<std::int64_t>(h.position) - (150 - best->end) : static_cast<std::int64_t>(h.position) - best->begin;
          if (h.reverse == rev && std::llabs(start - static_cast<std::int64_t>(b)) <= 0) { placed++; break; }
        }
      }
      std::snprintf(buf, sizeof buf, "%zu of %zu reads; %.1f SMEMs a read; %.1f µs a read on 8 threads", placed, reads.size(),
                    static_cast<double>(total) / static_cast<double>(reads.size()), smem / static_cast<double>(reads.size()) * 1e6);
      s.require("the longest SMEM of 150-base reads with 1 per cent substitutions lies on the true diagonal for at least 99 per cent", placed >= 3960, buf);
    }
  }

  return s.done();
}
