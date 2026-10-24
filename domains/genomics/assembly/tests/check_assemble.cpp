// Checks of assembly/assemble.hpp: the statistics N50, L50, NG50 and LG50 on a worked example;
// the number of contigs of error-free reads at low coverage against Lander and Waterman's
// N e^{−cσ} and the exact expectation for uniform reads; and the assembly of 30× short reads
// with errors of a 200 kb genome with repeats longer than k (700 bases) and shorter (20
// bases): contigs that lie exactly on the genome, covering nearly all of it, an N50 above the
// spacing of the long repeats, the short repeats inside contigs, the genome size from the
// spectrum, a cleaned graph without tips or bubbles, and the same contigs on 1 and 8 threads.

#include "truth.hpp"

#include <genomics/assembly/assemble.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);
  const int k = 31;

  // The statistics of contigs of 100, 80, 50, 30, 20 and 10 bases (290 in all).
  {
    const auto st = assembly::statistics(std::vector<std::size_t>{30, 100, 10, 50, 80, 20}, 400);
    s.require("statistics: total 290, largest 100, N50 80 and L50 2; against a genome of 400 bases NG50 50 and LG50 3",
              st.total == 290 && st.largest == 100 && st.n50 == 80 && st.l50 == 2 && st.ng50 == 50 && st.lg50 == 3 && st.contigs == 6);
    const auto big = assembly::statistics(std::vector<std::size_t>{30, 100, 10, 50, 80, 20}, 1000);
    s.require("statistics: NG50 is 0 when the contigs hold less than half the genome", big.ng50 == 0 && big.lg50 == 0);
  }

  // Lander–Waterman: error-free reads of 100 bases at coverage 2, 3 and 4 on 200 kb without
  // repeats; contigs break where consecutive reads overlap by fewer than k − 1 bases.
  {
    const std::size_t G = 200'000, L = 100;
    const Codes genome = simulate::random_sequence(G, 0.5, 51);
    const auto source = simulate::Source::haploid({genome});
    for (std::size_t c : {2u, 3u, 4u}) {
      const std::size_t n = c * G / L;
      const auto reads = truth::reads(source, n, {L, 0, 0, 0}, 52 + c, eight);
      assembly::AssemblyOptions o;
      o.min_count = 1;
      o.clean = false;
      const auto a = assembly::assemble(reads, eight, o);
      // Starts uniform on about G − L positions: a contig ends where the next start is more than
      // L − k + 1 bases on, (1 − g/G′)^n for each of the n − 1 spacings.
      const double g = static_cast<double>(L - k + 1), span = static_cast<double>(G - L);
      const double exact = 1 + (static_cast<double>(n) - 1) * std::pow(1 - g / span, static_cast<double>(n));
      const double lw = assembly::lander_waterman_contigs(static_cast<double>(n), static_cast<double>(L), static_cast<double>(G), static_cast<double>(k - 1));
      const double observed = static_cast<double>(a.contigs.size());
      s.check("Lander–Waterman: contigs of error-free reads at " + std::to_string(c) + "× against N e^(−cσ), σ = 1 − (k − 1)/L (deviation in √E)",
              (observed - lw) / std::sqrt(lw), 4, "");
      s.check("Lander–Waterman: the same against the exact expectation for uniform starts (deviation in √E)", (observed - exact) / std::sqrt(exact), 4, "");
    }
  }

  // 30× reads with errors of a genome with repeats longer and shorter than k.
  simulate::GenomeSpec spec;
  spec.length = 200'000;
  spec.repeats = {6, 3, 700, 0.0};
  auto genome = simulate::make_genome(spec, 60);
  struct Copy {
    std::size_t begin, end;
  };
  std::vector<Copy> short_copies;
  {
    calx::random::Stream r(61, 0x53484f52, 0, 0);
    for (int e = 0; e < 5; e++) {
      const Codes unit = simulate::random_sequence(20, 0.5, 62, static_cast<std::uint32_t>(e));
      for (int c = 0; c < 8; c++) {
        const std::size_t at = r.below(spec.length - 20);
        std::copy(unit.begin(), unit.end(), genome.sequence.begin() + static_cast<std::ptrdiff_t>(at));
        short_copies.push_back({at, at + 20});
      }
    }
  }
  const std::size_t G = spec.length;
  const simulate::ShortReadModel model;
  const auto reads = truth::reads(simulate::Source::haploid({genome.sequence}), 30 * G / model.length, model, 63, eight);
  const assembly::AssemblyOptions o;
  const auto a = assembly::assemble(reads, eight, o, G);
  const truth::GenomeIndex index(genome.sequence, k);

  s.check("assembly: the genome size from the spectrum (relative)", calx::verification::relative(a.fit.genome_size, static_cast<double>(G)), 0.03);
  {
    std::size_t placed = 0, placed_bases = 0, bases = 0;
    double kmers_in = 0, kmers_all = 0;
    std::vector<truth::Place> places;  // of the contigs placed once
    for (const auto& c : a.contigs) {
      const auto p = index.places(c);
      bases += c.size();
      if (!p.empty()) {
        placed++;
        placed_bases += c.size();
        if (p.size() == 1) places.push_back(p[0]);
      }
      const double km = static_cast<double>(c.size() - static_cast<std::size_t>(k) + 1);
      kmers_in += index.kmer_fraction(c) * km;
      kmers_all += km;
    }
    s.require("assembly: contigs and their statistics", !a.contigs.empty(),
              std::to_string(a.stats.contigs) + " contigs, " + std::to_string(a.stats.total) + " bases, largest " + std::to_string(a.stats.largest) + ", N50 " +
                  std::to_string(a.stats.n50) + ", NG50 " + std::to_string(a.stats.ng50) + ", threshold " + std::to_string(a.threshold));
    s.check("assembly: contigs not found exactly on the genome or its reverse complement (fraction)", 1 - static_cast<double>(placed) / static_cast<double>(a.contigs.size()), 0.01);
    s.check("assembly: contig bases in contigs not found exactly (fraction)", 1 - static_cast<double>(placed_bases) / static_cast<double>(bases), 1e-3);
    s.check("assembly: contig k-mers absent from the genome (fraction)", 1 - kmers_in / kmers_all, 1e-4);
    // The genome's k-mers that the contigs hold.
    kmers::KmerCounter held = kmers::count_kmers(a.contigs, eight);
    std::size_t genome_kmers = 0, covered = 0;
    kmers::count_kmers({genome.sequence}, eight).for_each([&](sequence::Kmer x, std::uint32_t) {
      genome_kmers++;
      covered += held.count(x) > 0;
    });
    s.check("assembly: genome k-mers missing from the contigs (fraction)", 1 - static_cast<double>(covered) / static_cast<double>(genome_kmers), 5e-3);
    // The spacing of the long repeats: the genome less their copies, over the gaps between them.
    std::size_t copies = 0, repeat_bases = 0;
    for (const auto& r : genome.repeats) {
      copies++;
      repeat_bases += r.end - r.begin;
    }
    const double spacing = static_cast<double>(G - repeat_bases) / static_cast<double>(copies + 1);
    s.require("assembly: N50 above the mean spacing of the repeats longer than k", static_cast<double>(a.stats.n50) > spacing,
              "N50 " + std::to_string(a.stats.n50) + ", spacing " + calx::verification::exponential(spacing) + " over " + std::to_string(copies) + " copies");
    // Short repeats (20 < k − 1 bases) do not break contigs: each copy away from the long repeats
    // and the ends lies inside a contig placed once, with k − 1 bases to spare on each side.
    std::size_t tested = 0, inside = 0;
    for (const auto& c : short_copies) {
      bool near = c.begin < 2 * static_cast<std::size_t>(k) || c.end + 2 * static_cast<std::size_t>(k) > G;
      for (const auto& r : genome.repeats) near = near || (c.begin < r.end + 2 * static_cast<std::size_t>(k) && r.begin < c.end + 2 * static_cast<std::size_t>(k));
      if (near) continue;
      tested++;
      for (const auto& p : places)
        if (p.begin + static_cast<std::size_t>(k) <= c.begin && c.end + static_cast<std::size_t>(k) <= p.end) { inside++; break; }
    }
    s.check("assembly: copies of 20-base repeats not inside a contig (fraction)", 1 - static_cast<double>(inside) / static_cast<double>(tested), 0.1);
    // At the spectrum's threshold few errors reach the graph; what does is cleaned (check_clean
    // stresses the cleaning at a threshold of 2).
    s.require("assembly: the cleaned graph has no tip and no bubble", assembly::find_tips(a.cleaning.graph).empty() && assembly::find_bubbles(a.cleaning.graph).empty(),
              std::to_string(a.raw.size()) + " unitigs raw, " + std::to_string(a.cleaning.graph.size()) + " cleaned, " + std::to_string(a.cleaning.steps.size()) + " steps");
  }
  {
    const auto b = assembly::assemble(reads, one, o, G);
    s.require("assembly: the same contigs on 1 and 8 threads", b.contigs == a.contigs);
  }
  return s.done();
}
