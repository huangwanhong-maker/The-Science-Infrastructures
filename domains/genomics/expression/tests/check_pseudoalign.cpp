// Checks of pseudoalignment: a hand-built transcriptome of three transcripts on three exons (AB,
// AC, BC) whose six classes are known; every k-mer's class against an independent construction
// (all (k-mer, transcript) pairs sorted); every read's class against brute force (the
// intersection, by sorted lists, of the sets of the read's k-mers that some transcript holds),
// single and paired, with sequencing errors; on a small transcriptome, against the definition
// itself (each k-mer searched for in every transcript and its reverse complement); the true
// transcript of an error-free read always in its class; the result the same on 1 and 8 threads;
// and the speed of the index and of pseudoalignment, reported.

#include <genomics/expression/pseudoalign.hpp>
#include <genomics/expression/transcriptome.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace genomics::expression;
namespace seq = genomics::sequence;

namespace {

Codes random_dna(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(4));
  return x;
}

// The independent construction: every (canonical k-mer, transcript) pair, sorted and unique.
struct Brute {
  std::vector<std::pair<Kmer, std::uint32_t>> pairs;
  int k;
  Brute(const Transcriptome& tx, int k_) : k(k_) {
    for (std::uint32_t t = 0; t < tx.size(); t++)
      seq::for_each_kmer(tx.transcripts[t].sequence, k, [&](std::size_t, Kmer f, Kmer r) { pairs.push_back({std::min(f, r), t}); });
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  }
  std::vector<std::uint32_t> set_of(Kmer c) const {
    std::vector<std::uint32_t> out;
    auto it = std::lower_bound(pairs.begin(), pairs.end(), std::pair<Kmer, std::uint32_t>{c, 0});
    for (; it != pairs.end() && it->first == c; ++it) out.push_back(it->second);
    return out;
  }
  // The read's class by definition: the intersection over its k-mers held by some transcript.
  std::vector<std::uint32_t> read_class(const Codes& a, const Codes* b) const {
    std::vector<std::uint32_t> cur;
    bool started = false;
    const auto visit = [&](std::size_t, Kmer f, Kmer r) {
      const auto s = set_of(std::min(f, r));
      if (s.empty()) return;
      if (!started) { cur = s; started = true; return; }
      std::vector<std::uint32_t> next;
      std::set_intersection(cur.begin(), cur.end(), s.begin(), s.end(), std::back_inserter(next));
      cur.swap(next);
    };
    if (a.size() >= static_cast<std::size_t>(k)) seq::for_each_kmer(a, k, visit);
    if (b && b->size() >= static_cast<std::size_t>(k)) seq::for_each_kmer(*b, k, visit);
    return cur;
  }
};

std::vector<std::uint32_t> as_vector(std::span<const std::uint32_t> s) { return {s.begin(), s.end()}; }

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x50534555, 0, 0);
  calx::scheduler::ThreadPool one(1), eight(8);

  // ---- three transcripts on three exons ------------------------------------------------------------------------------------
  {
    const Codes A = random_dna(r, 60), B = random_dna(r, 60), C = random_dna(r, 60);
    const auto cat = [](const Codes& x, const Codes& y) { Codes z = x; z.insert(z.end(), y.begin(), y.end()); return z; };
    const Transcriptome tx = make_transcriptome({{"t0", "g", cat(A, B)}, {"t1", "g", cat(A, C)}, {"t2", "g", cat(B, C)}});
    const Index index(tx, 15);
    const auto& cl = index.classes();
    std::vector<std::vector<std::uint32_t>> got;
    for (std::uint32_t c = 0; c < cl.size(); c++) got.push_back(as_vector(cl.members(c)));
    std::sort(got.begin(), got.end());
    const std::vector<std::vector<std::uint32_t>> want = {{0}, {0, 1}, {0, 2}, {1}, {1, 2}, {2}};
    s.require("transcripts AB, AC, BC have the six classes {0}, {1}, {2} (the junctions) and {0, 1}, {0, 2}, {1, 2} (the exons)", got == want);
    Index::Scratch sc;
    const Codes in_a(A.begin() + 10, A.begin() + 50), across(A.begin() + 40, A.end());
    Codes junction = across;
    junction.insert(junction.end(), B.begin(), B.begin() + 20);
    const auto ca = index.pseudoalign(in_a, nullptr, sc), cj = index.pseudoalign(junction, nullptr, sc);
    s.require("a read inside exon A is {0, 1}; a read across the junction of A and B is {0}",
              ca != Index::NONE && as_vector(cl.members(ca)) == std::vector<std::uint32_t>{0, 1} && cj != Index::NONE &&
                  as_vector(cl.members(cj)) == std::vector<std::uint32_t>{0});
    const Codes rc = seq::reverse_complement(junction);
    s.require("its reverse complement has the same class", index.pseudoalign(rc, nullptr, sc) == cj);
    // A pair with one mate in A and one in C: {0, 1} ∩ {1, 2} = {1}.
    const Codes in_c(C.begin() + 5, C.begin() + 45);
    const auto cp = index.pseudoalign(in_a, &in_c, sc);
    s.require("a pair with a mate in A and a mate in C is {1}", cp != Index::NONE && as_vector(cl.members(cp)) == std::vector<std::uint32_t>{1});
  }

  // ---- a simulated transcriptome against the independent construction ------------------------------------------------------------------------------------
  {
    TranscriptomeSpec spec;
    spec.genes = 120;
    spec.paralogs = 0.3;
    spec.paralog_divergence = 0.02;
    const auto sim = simulate_transcriptome(spec, 7);
    const Transcriptome& tx = sim.transcriptome;
    const int k = 31;
    const Index index(tx, k);
    const Brute brute(tx, k);
    bool kmers = true;
    std::size_t distinct = 0;
    for (std::size_t i = 0; i < brute.pairs.size();) {
      std::size_t j = i;
      std::vector<std::uint32_t> set;
      while (j < brute.pairs.size() && brute.pairs[j].first == brute.pairs[i].first) set.push_back(brute.pairs[j++].second);
      const auto c = index.class_of(brute.pairs[i].first);
      kmers = kmers && c != Index::NONE && as_vector(index.classes().members(c)) == set;
      distinct++;
      i = j;
    }
    s.require("every k-mer of " + std::to_string(tx.size()) + " transcripts (isoforms and paralogs) has the class of the sorted (k-mer, transcript) pairs",
              kmers && distinct == index.kmers(), std::to_string(index.kmers()) + " k-mers, " + std::to_string(index.classes().size()) + " classes");
    bool unique = true;
    for (std::uint32_t c = 0; c < index.classes().size(); c++) unique = unique && index.classes().find(index.classes().members(c)) == c;
    s.require("each class is held once", unique);

    std::vector<double> rho = simulate_abundances(tx.size(), 1.0, 0.1, 11);
    const auto d = FragmentLengths::normal(200, 25);
    for (const bool paired : {false, true}) {
      ReadModel model;
      model.length = 75;
      model.paired = paired;
      model.error = 0.003;
      const auto frags = simulate_fragments(tx, rho, d, model, 6000, 13);
      const Pseudoalignment p = pseudoalign(index, frags, eight);
      bool same = true;
      std::size_t aligned = 0;
      for (std::size_t i = 0; i < frags.size(); i++) {
        const auto want = brute.read_class(frags[i].read1, paired ? &frags[i].read2 : nullptr);
        const std::uint32_t c = p.fragment_class[i];
        if (want.empty()) same = same && c == ClassTable::NONE;
        else { same = same && c != ClassTable::NONE && as_vector(p.classes.members(c)) == want; aligned++; }
      }
      s.require(std::string("every ") + (paired ? "pair's" : "read's") + " class (6000 fragments, errors 0.003 per base) equals the brute-force intersection", same && aligned == p.aligned,
                std::to_string(p.aligned) + " pseudoaligned, " + std::to_string(p.classes.size() - index.classes().size()) + " classes new to the index");
    }
    ReadModel clean;
    const auto frags = simulate_fragments(tx, rho, d, clean, 6000, 17);
    const Pseudoalignment p = pseudoalign(index, frags, one);
    bool truth = true;
    for (std::size_t i = 0; i < frags.size(); i++) {
      const auto c = p.fragment_class[i];
      truth = truth && c != ClassTable::NONE;
      if (c != ClassTable::NONE) {
        const auto m = p.classes.members(c);
        truth = truth && std::binary_search(m.begin(), m.end(), frags[i].transcript);
      }
    }
    s.require("the class of every error-free read holds the transcript it came from", truth);
    const Pseudoalignment p8 = pseudoalign(index, frags, eight);
    bool threads = p8.counts == p.counts && p8.fragment_class == p.fragment_class && p8.classes.size() == p.classes.size();
    for (std::uint32_t c = 0; threads && c < p.classes.size(); c++) threads = as_vector(p.classes.members(c)) == as_vector(p8.classes.members(c));
    s.require("the classes, their numbering and the counts are the same on 1 and 8 threads", threads);
  }

  // ---- the definition itself ------------------------------------------------------------------------------------
  {
    std::vector<Transcript> ts;
    const Codes shared = random_dna(r, 80);
    for (int t = 0; t < 6; t++) {
      Codes x = random_dna(r, 100 + r.below(100));
      if (t % 2) x.insert(x.begin() + static_cast<std::ptrdiff_t>(r.below(x.size())), shared.begin(), shared.end());
      if (t == 4) x = seq::reverse_complement(x);
      ts.push_back({"t" + std::to_string(t), "g" + std::to_string(t / 2), x});
    }
    const Transcriptome tx = make_transcriptome(ts);
    const int k = 15;
    const Index index(tx, k);
    Index::Scratch sc;
    bool ok = true;
    for (int i = 0; i < 300 && ok; i++) {
      // A read from a random transcript, on a random strand, with a few errors.
      const auto& src = tx.transcripts[r.below(tx.size())].sequence;
      const std::size_t len = 30 + r.below(40), b = r.below(src.size() - len);
      Codes read(src.begin() + static_cast<std::ptrdiff_t>(b), src.begin() + static_cast<std::ptrdiff_t>(b + len));
      if (r.bernoulli(0.5)) read = seq::reverse_complement(read);
      if (r.bernoulli(0.3)) read[r.below(len)] = static_cast<Code>(r.below(4));
      std::vector<std::uint32_t> cur;
      bool started = false;
      for (std::size_t p = 0; p + k <= read.size(); p++) {
        const Codes w(read.begin() + static_cast<std::ptrdiff_t>(p), read.begin() + static_cast<std::ptrdiff_t>(p + k));
        const Codes wr = seq::reverse_complement(w);
        std::vector<std::uint32_t> in;
        for (std::uint32_t t = 0; t < tx.size(); t++) {
          const auto& x = tx.transcripts[t].sequence;
          if (std::search(x.begin(), x.end(), w.begin(), w.end()) != x.end() || std::search(x.begin(), x.end(), wr.begin(), wr.end()) != x.end()) in.push_back(t);
        }
        if (in.empty()) continue;
        if (!started) { cur = in; started = true; continue; }
        std::vector<std::uint32_t> next;
        std::set_intersection(cur.begin(), cur.end(), in.begin(), in.end(), std::back_inserter(next));
        cur.swap(next);
      }
      const auto c = index.pseudoalign(read, nullptr, sc);
      const std::vector<std::uint32_t> got = c != Index::NONE ? as_vector(index.classes().members(c)) : sc.current;
      ok = got == cur;
    }
    s.require("300 reads with errors, on either strand, have the class that searching every transcript and its reverse complement for each k-mer gives", ok);
  }

  // ---- speed ------------------------------------------------------------------------------------
  {
    TranscriptomeSpec spec;
    spec.genes = 2000;
    const auto sim = simulate_transcriptome(spec, 23);
    std::size_t bases = 0;
    for (const auto& t : sim.transcriptome.transcripts) bases += t.sequence.size();
    const auto t0 = std::chrono::steady_clock::now();
    const Index index(sim.transcriptome, 31);
    const double build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    ReadModel model;
    model.error = 0.002;
    const auto frags = simulate_fragments(sim.transcriptome, simulate_abundances(sim.transcriptome.size(), 1.0, 0.1, 3), FragmentLengths::normal(200, 25), model, 200000, 5);
    std::vector<Codes> reads;
    for (const auto& f : frags) reads.push_back(f.read1);
    const auto t1 = std::chrono::steady_clock::now();
    const Pseudoalignment p1 = pseudoalign(index, reads, {}, one);
    const double a1 = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    const auto t8 = std::chrono::steady_clock::now();
    const Pseudoalignment p8 = pseudoalign(index, reads, {}, eight);
    const double a8 = std::chrono::duration<double>(std::chrono::steady_clock::now() - t8).count();
    std::printf("index of %zu transcripts (%.1f Mb, %zu k-mers, %zu classes, %.1f MB) in %.2f s; 200 000 reads of 75 bases: %.2f M reads/s on 1 thread, %.2f on 8\n",
                sim.transcriptome.size(), static_cast<double>(bases) / 1e6, index.kmers(), index.classes().size(), static_cast<double>(index.memory_bytes()) / 1e6, build,
                0.2 / a1, 0.2 / a8);
    s.require("200 000 reads pseudoaligned the same on 1 and 8 threads", p1.counts == p8.counts && p1.fragment_class == p8.fragment_class,
              std::to_string(p1.aligned) + " aligned");
  }
  return s.done();
}
