// Checks of kmers/spectrum.hpp on simulated reads of known genomes: the peak and the fitted
// coverage at the k-mer coverage c (L − k + 1)/L of error-free reads, and at that coverage
// thinned by the chance of an error-free k-mer for reads with errors; the error threshold
// separating the reads' errors from the genome's k-mers; the genome size within a few per
// cent, repeats included; the heterozygosity of a diploid genome against the truth of its
// mutations, and nearly none for a haploid one; the spectrum's sums against the counter.

#include <genomics/kmers/counter.hpp>
#include <genomics/kmers/spectrum.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

// Reads of a source in shards on the pool (each read its own stream, so the sample is the same).
std::vector<simulate::Read> reads_of(const simulate::Source& source, std::size_t n, const simulate::ShortReadModel& model, std::uint64_t seed,
                                     calx::scheduler::ThreadPool& pool) {
  const std::size_t shards = 64;
  std::vector<std::vector<simulate::Read>> part(shards);
  calx::scheduler::parallel_for(pool, shards, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t q = b; q < e; q++) {
      const auto r = calx::scheduler::part(n, shards, q);
      part[q] = simulate::short_reads(source, r.size(), model, seed, r.begin);
    }
  });
  std::vector<simulate::Read> out;
  for (auto& p : part)
    for (auto& r : p) out.push_back(std::move(r));
  return out;
}

kmers::KmerCounter count_reads(const std::vector<simulate::Read>& reads, calx::scheduler::ThreadPool& pool, int k = 31) {
  kmers::CounterOptions o;
  o.k = k;
  kmers::KmerCounter c(o);
  c.add(reads.size(), [&](std::size_t i) -> const Codes& { return reads[i].bases; }, pool);
  return c;
}

// The mean over the k-mers of a read of the chance that none of its bases is an error, for the
// simulator's model: substitutions rising linearly along the read, and indels.
double error_free_kmers(const simulate::ShortReadModel& m, int k) {
  const std::size_t L = m.length;
  double sum = 0;
  for (std::size_t j = 0; j + static_cast<std::size_t>(k) <= L; j++) {
    double p = 1;
    for (std::size_t i = j; i < j + static_cast<std::size_t>(k); i++) {
      const double x = static_cast<double>(i) / static_cast<double>(L - 1);
      p *= (1 - m.error_first - (m.error_last - m.error_first) * x) * (1 - m.indel);
    }
    sum += p;
  }
  return sum / static_cast<double>(L - static_cast<std::size_t>(k) + 1);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(8);
  const int k = 31;

  // A. Error-free reads of a random haploid genome at 30×.
  {
    const std::size_t G = 400'000, L = 150;
    const auto genome = simulate::random_sequence(G, 0.45, 101);
    simulate::ShortReadModel model{L, 0, 0, 0};
    const std::size_t n = 30 * G / L;
    const auto reads = reads_of(simulate::Source::haploid({genome}), n, model, 102, pool);
    const auto counts = count_reads(reads, pool);
    const auto sp = kmers::spectrum(counts, pool);
    const auto fit = kmers::fit_spectrum(sp);
    // Reads start uniformly on [0, G − L − 8]: the k-mer coverage is n (L − k + 1) / (G − L − 7).
    const double ck = static_cast<double>(n) * static_cast<double>(L - k + 1) / static_cast<double>(G - L - 7);
    s.require("spectrum: error-free reads at 30× have their peak at c (L − k + 1)/L = 24", std::abs(static_cast<double>(fit.peak) - ck) <= 1,
              "peak " + std::to_string(fit.peak) + ", expected " + calx::verification::exponential(ck));
    s.check("spectrum: the fitted coverage of error-free reads against n (L − k + 1)/(G − L − 7) (relative)", calx::verification::relative(fit.coverage, ck), 0.01);
    s.check("spectrum: the genome size of error-free reads (relative)", calx::verification::relative(fit.genome_size, static_cast<double>(G)), 0.01);
    std::uint64_t distinct = 0, total = 0;
    for (std::size_t c = 1; c < sp.histogram.size(); c++) distinct += sp.histogram[c];
    total = sp.total;
    s.require("spectrum: the histogram sums to the distinct k-mers, and the total to the k-mers of the reads",
              distinct == counts.distinct() && total == n * (L - k + 1));
  }

  // B. Reads with errors of a haploid genome with repeats: four copies each of ten elements.
  {
    simulate::GenomeSpec spec;
    spec.length = 500'000;
    spec.repeats = {10, 4, 2000, 0.0};
    const auto genome = simulate::make_genome(spec, 201);
    const simulate::ShortReadModel model;  // 150 bases, errors 0.1 % to 1 %, indels 0.02 %
    const std::size_t L = model.length, G = spec.length, n = 30 * G / L;
    const auto reads = reads_of(simulate::Source::haploid({genome.sequence}), n, model, 202, pool);
    const auto counts = count_reads(reads, pool);
    const auto sp = kmers::spectrum(counts, pool);
    const auto fit = kmers::fit_spectrum(sp);
    const double ck = static_cast<double>(n) * static_cast<double>(L - k + 1) / static_cast<double>(G - L - 7) * error_free_kmers(model, k);
    s.check("spectrum: the fitted coverage of reads with errors against c (L − k + 1)/L thinned by error-free k-mers (relative)", calx::verification::relative(fit.coverage, ck), 0.02);
    s.require("spectrum: the error threshold lies between 2 and λ/3", fit.threshold >= 2 && fit.threshold <= ck / 3,
              "threshold " + std::to_string(fit.threshold) + ", λ " + calx::verification::exponential(fit.coverage));
    s.check("spectrum: the genome size of reads with errors, repeats counted with multiplicity (relative)", calx::verification::relative(fit.genome_size, static_cast<double>(G)), 0.02);
    // The genome's k-mers against the threshold.
    kmers::CounterOptions o;
    const auto truth = kmers::count_kmers({genome.sequence}, pool, o);
    std::size_t genomic_solid = 0, error_below = 0, below = 0;
    truth.for_each([&](sequence::Kmer x, std::uint32_t) { genomic_solid += counts.count(x) >= fit.threshold; });
    counts.for_each([&](sequence::Kmer x, std::uint32_t c) {
      if (c >= fit.threshold) return;
      below++;
      error_below += truth.count(x) == 0;
    });
    s.check("spectrum: genome k-mers below the threshold (fraction)", 1 - static_cast<double>(genomic_solid) / static_cast<double>(truth.distinct()), 1e-3);
    s.check("spectrum: k-mers below the threshold that are in the genome (fraction)", 1 - static_cast<double>(error_below) / static_cast<double>(below), 1e-3);
    s.check("spectrum: the 4-copy repeats, 4 × 10 × 2000 bases, found in the 4λ component (relative)",
            calx::verification::relative(fit.kmers[4], 10 * (2000 - k + 1)), 0.15);
    s.check("spectrum: a haploid genome's heterozygosity", fit.heterozygosity, 2e-4);
  }

  // C. A diploid genome, every variant heterozygous, at 40× over both haplotypes.
  {
    const std::size_t G = 400'000;
    const auto ref = simulate::random_sequence(G, 0.5, 301);
    simulate::MutationRates rates;
    rates.snv = 2e-3;
    rates.indel = 0;
    rates.heterozygous = 1;
    const auto d = simulate::mutate_diploid(ref, rates, 302);
    const simulate::ShortReadModel model;
    const std::size_t L = model.length, n = 40 * G / L;
    const auto reads = reads_of(simulate::Source::diploid({d}), n, model, 303, pool);
    const auto counts = count_reads(reads, pool);
    const auto sp = kmers::spectrum(counts, pool);
    const auto fit = kmers::fit_spectrum(sp);
    // The truth: the rate of heterozygous sites, and the fraction of windows of k bases holding one.
    const double r = static_cast<double>(d.truth.size()) / static_cast<double>(G);
    std::vector<char> het(G, 0);
    for (const auto& v : d.truth) het[v.position] = 1;
    std::size_t windows = 0, het_windows = 0, run = 0;
    for (std::size_t i = 0; i < G; i++) {
      run = het[i] ? 0 : run + 1;
      // The window ending at i holds a site when the last site lies within k − 1 bases of i.
      if (i + 1 >= static_cast<std::size_t>(k)) {
        windows++;
        het_windows += run < static_cast<std::size_t>(k);
      }
    }
    const double f = static_cast<double>(het_windows) / static_cast<double>(windows);
    s.check("spectrum: the fraction of heterozygous windows against the truth of the mutations (relative)", calx::verification::relative(fit.heterozygous_windows, f), 0.05);
    s.check("spectrum: the heterozygosity against the rate of heterozygous SNVs, 2e-3 (relative)", calx::verification::relative(fit.heterozygosity, r), 0.08);
    s.check("spectrum: the haploid genome size of a diploid sample (relative)", calx::verification::relative(fit.genome_size, static_cast<double>(G)), 0.02);
    const double ck = static_cast<double>(n) * static_cast<double>(L - k + 1) / static_cast<double>(G - L - 7) * error_free_kmers(model, k);
    s.check("spectrum: a diploid sample's coverage λ is the homozygous peak, against the simulation (relative)", calx::verification::relative(fit.coverage, ck), 0.02);
  }

  // D. A high heterozygosity puts the mode at λ/2; the coverage given resolves it.
  {
    const std::size_t G = 300'000;
    const auto ref = simulate::random_sequence(G, 0.5, 401);
    simulate::MutationRates rates;
    rates.snv = 2e-2;
    rates.indel = 0;
    rates.heterozygous = 1;
    const auto d = simulate::mutate_diploid(ref, rates, 402);
    const simulate::ShortReadModel model{150, 0, 0, 0};
    const std::size_t n = 50 * G / 150;
    const auto reads = reads_of(simulate::Source::diploid({d}), n, model, 403, pool);
    const auto counts = count_reads(reads, pool);
    const auto sp = kmers::spectrum(counts, pool);
    const double ck = static_cast<double>(n) * 120.0 / static_cast<double>(G - 157);
    kmers::FitOptions fo;
    fo.coverage = ck * 1.1;
    const auto fit = kmers::fit_spectrum(sp, fo);
    const auto blind = kmers::fit_spectrum(sp);
    s.require("spectrum: at 2 % heterozygosity the mode is the half-coverage peak", std::abs(static_cast<double>(blind.peak) - ck / 2) <= 2,
              "peak " + std::to_string(blind.peak) + ", λ/2 " + calx::verification::exponential(ck / 2));
    s.check("spectrum: with the coverage given, the fitted λ (relative)", calx::verification::relative(fit.coverage, ck), 0.02);
    s.check("spectrum: with the coverage given, the heterozygosity against the truth (relative)",
            calx::verification::relative(fit.heterozygosity, static_cast<double>(d.truth.size()) / static_cast<double>(G)), 0.08);
  }
  return s.done();
}
