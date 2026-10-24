// The speed of counting and assembly: a simulated genome of 5 Mb (GENOMICS_SPEED_MB to change
// it) with repeats, 30× of 150-base reads with errors, counted and assembled on 1 and 8
// threads. The times and the memory are reported; the checks are that the two runs give the
// same contigs, that the genome size is read from the spectrum within 3 %, and that the
// contigs' k-mers are the genome's.

#include "truth.hpp"

#include <genomics/assembly/assemble.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

// The resident peak of this process in MB, where the system says it (Linux), else 0.
double peak_resident_mb() {
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line))
    if (line.rfind("VmHWM:", 0) == 0) return std::strtod(line.c_str() + 6, nullptr) / 1024.0;
  return 0;
}

std::string fixed(double x, int digits = 2) {
  char b[64];
  std::snprintf(b, sizeof b, "%.*f", digits, x);
  return b;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const char* mb = std::getenv("GENOMICS_SPEED_MB");
  const std::size_t G = static_cast<std::size_t>((mb ? std::strtod(mb, nullptr) : 5.0) * 1e6);
  simulate::GenomeSpec spec;
  spec.length = G;
  spec.repeats = {G / 250'000, 4, 1500, 0.01};
  const auto genome = simulate::make_genome(spec, 80);
  const simulate::ShortReadModel model;
  const std::size_t n = 30 * G / model.length;
  std::vector<Codes> reads;
  {
    calx::scheduler::ThreadPool pool(8);
    const auto t = std::chrono::steady_clock::now();
    reads = truth::reads(simulate::Source::haploid({genome.sequence}), n, model, 81, pool);
    s.require("speed: the sample", true,
              fixed(static_cast<double>(G) / 1e6, 1) + " Mb at 30×, " + std::to_string(reads.size()) + " reads of 150 bases, simulated in " +
                  fixed(std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count()) + " s on 8 threads");
  }
  std::vector<std::vector<Codes>> contigs;
  for (std::size_t threads : {1u, 8u}) {
    calx::scheduler::ThreadPool pool(threads);
    const auto a = assembly::assemble(reads, pool, {}, G);
    s.require("speed: " + std::to_string(threads) + (threads == 1 ? " thread" : " threads"), true,
              "count " + fixed(a.seconds.count) + " s, spectrum " + fixed(a.seconds.spectrum) + " s, graph " + fixed(a.seconds.graph) + " s, clean " + fixed(a.seconds.clean) +
                  " s, in all " + fixed(a.seconds.total()) + " s; counter " + fixed(static_cast<double>(a.count_bytes) / 1048576.0, 0) + " MB, graph " +
                  fixed(static_cast<double>(a.graph_bytes) / 1048576.0, 0) + " MB; " + std::to_string(a.spectrum.distinct) + " distinct k-mers, " +
                  std::to_string(a.solid_kmers) + " solid; " + std::to_string(a.stats.contigs) + " contigs, N50 " + std::to_string(a.stats.n50) + ", NG50 " +
                  std::to_string(a.stats.ng50));
    if (threads == 8) {
      s.check("speed: the genome size from the spectrum (relative)", calx::verification::relative(a.fit.genome_size, static_cast<double>(G)), 0.03);
      const auto truth = kmers::count_kmers({genome.sequence}, pool);
      std::size_t all = 0, in = 0;
      for (const auto& c : a.contigs)
        sequence::for_each_kmer(c, 31, [&](std::size_t, sequence::Kmer f, sequence::Kmer) {
          all++;
          in += truth.count(f) > 0;
        });
      s.check("speed: contig k-mers absent from the genome (fraction)", 1 - static_cast<double>(in) / static_cast<double>(all), 1e-4);
    }
    contigs.push_back(a.contigs);
  }
  s.require("speed: the same contigs on 1 and 8 threads", contigs[0] == contigs[1]);
  s.require("speed: the resident peak of the process", true, fixed(peak_resident_mb(), 0) + " MB, the reads " + fixed(static_cast<double>(n) * 174 / 1048576.0, 0) + " MB of it");
  return s.done();
}
