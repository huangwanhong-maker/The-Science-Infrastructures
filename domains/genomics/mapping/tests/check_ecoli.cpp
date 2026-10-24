// The mapper on Escherichia coli K-12 MG1655 (NC_000913.3, from the data library's cache; the
// check is skipped when it is absent): short reads (seeded by minimizers and by SMEMs) and long
// reads simulated from it map to their true places, with the time of the index (the minimizer
// index on eight threads, the FMD-index on one) and of the mapping on one and eight threads,
// and the accuracy curve by mapping quality. Its seven rRNA operons (about 5 kb each, nearly identical)
// and its insertion sequences are the repeats where reads become ambiguous.

#include <genomics/mapping/evaluate.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

}  // namespace

int main() {
  calx::verification::Suite s;
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
  const auto file = cache / "ncbi-genomes" / "NC_000913.3.gb";
  if (!std::filesystem::exists(file)) {
    s.skip("Escherichia coli NC_000913.3 is not in the data cache (data/cache/ncbi-genomes)");
    return s.done();
  }
  const std::vector<Codes> genome{sequence::encode_dna(sequence::read_genbank(read_file(file)).at(0).sequence)};
  s.require("Escherichia coli K-12 MG1655 is read: 4 641 652 bases", genome[0].size() == 4641652, std::to_string(genome[0].size()));
  calx::scheduler::ThreadPool one(1), eight(8);
  const auto source = simulate::Source::haploid(genome);
  char buf[320];

  for (const int mode : {0, 1, 2}) {
    const bool long_reads = mode == 2;
    auto t0 = std::chrono::steady_clock::now();
    const mapping::Mapper mapper(eight, genome,
                                 mode == 0 ? mapping::MapperOptions::short_reads() : mode == 1 ? mapping::MapperOptions::short_reads_smem() : mapping::MapperOptions::long_reads());
    const double build = seconds_since(t0);
    std::vector<simulate::Read> reads;
    if (long_reads) {
      simulate::LongReadModel model;
      model.mean_length = 8000;
      model.error = 0.08;
      reads = simulate::long_reads(source, 400, model, 31);
    } else {
      simulate::ShortReadModel model;
      model.error_first = 0.005;
      model.error_last = 0.015;
      reads = simulate::short_reads(source, 40000, model, 30);
    }
    std::size_t bases = 0;
    for (const auto& r : reads) bases += r.bases.size();
    t0 = std::chrono::steady_clock::now();
    const auto on8 = mapper.map_all(eight, reads, long_reads ? 1 : 16);
    const double s8 = seconds_since(t0);
    t0 = std::chrono::steady_clock::now();
    const auto on1 = mapper.map_all(one, reads, long_reads ? 1 : 16);
    const double s1 = seconds_since(t0);
    const auto acc = mapping::accuracy_of(on8, reads);
    const char* what = mode == 2 ? "long reads (8 kb, 8% errors)" : mode == 1 ? "short reads (150 b, 1% substitutions) seeded by SMEMs" : "short reads (150 b, 1% substitutions)";
    std::snprintf(buf, sizeof buf, "%.4f mapped, %.4f correct of %zu; at MAPQ >= 20: %.4f of reads, error rate %.1e", acc.fraction_mapped(), acc.fraction_correct(), acc.reads,
                  acc.curve[4].fraction_mapped, acc.curve[4].error_rate);
    s.require(std::string("Escherichia coli: ") + what + " map to their true places (above 0.98), wrong at MAPQ ≥ 20 less than once in a thousand",
              acc.fraction_correct() > 0.98 && acc.curve[4].error_rate < 1e-3, buf);
    bool same = on1.size() == on8.size();
    for (std::size_t i = 0; same && i < on1.size(); i++) {
      same = on1[i].size() == on8[i].size();
      for (std::size_t j = 0; same && j < on1[i].size(); j++)
        same = on1[i][j].target_begin == on8[i][j].target_begin && on1[i][j].alignment.cigar == on8[i][j].alignment.cigar && on1[i][j].mapq == on8[i][j].mapq;
    }
    s.require(std::string("Escherichia coli: ") + what + ", the same mappings on one and eight threads", same);
    const mapping::Throughput t1{reads.size(), bases, s1}, t8{reads.size(), bases, s8};
    std::snprintf(buf, sizeof buf, "index %.2f s; %.0f reads/s (%.2f Mb/s) on 1 thread, %.0f reads/s (%.2f Mb/s) on 8", build, t1.reads_per_second(),
                  t1.bases_per_second() / 1e6, t8.reads_per_second(), t8.bases_per_second() / 1e6);
    s.require(std::string("Escherichia coli: ") + what + ", throughput", true, buf);
  }
  return s.done();
}
