// The toolbox's operations on the user's reads (src/tools/data.js): a FASTQ file checked; reads
// mapped to a genome, with their mapping qualities, coverage and SAM; the variants of a diploid
// called from short reads; a genome assembled from its reads; and transcripts quantified from an
// RNA-seq sample.
//
// Files of reads may be large: they are read a record at a time (inputs.hpp), checked or counted as
// they stream, mapped in batches, and held only where the method needs them all (the caller's
// pileups, the pseudoalignment); each tool's bounds say how many reads it takes, and the result
// what it read. Every computation is the domain's (genomics::mapping, genomics::variant,
// genomics::kmers, genomics::assembly, genomics::expression): the engine gathers and reports.

#pragma once

#include "genomes.hpp"
#include "inputs.hpp"
#include "shared.hpp"
#include "shared_assembly.hpp"
#include "shared_expression.hpp"
#include "shared_mapping.hpp"

#include <genomics/kmers/counter.hpp>
#include <genomics/mapping/evaluate.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/mapping/sam.hpp>
#include <genomics/sequence/composition.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/vcf.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace workbench {

namespace gvar = genomics::variant;

namespace detail {

inline double seconds_since(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }

// The name of a sample: its file's name without the extensions of reads.
inline std::string sample_name(const std::string& source) {
  std::string s = source;
  for (const char* ext : {".gz", ".fastq", ".fq", ".fasta", ".fa", ".txt"})
    if (s.size() > std::string(ext).size() && s.compare(s.size() - std::string(ext).size(), std::string::npos, ext) == 0) s.resize(s.size() - std::string(ext).size());
  return s.empty() ? "sample" : s;
}

inline void require_reads(const Json& request) {
  if (str(request, "reads_file").empty() && !request["text"].is_string()) throw Fault("choose a file of reads");
}

}  // namespace detail

// A file of reads read through: their number and lengths, GC and N, and the mean Phred quality at
// each position along the read (in bins of positions for long reads).
BENCH_OP("tool.reads.qc") {
  detail::require_reads(request);
  detail::Bounds b;
  b.bytes = std::size_t{8} << 30;
  b.records = 200'000'000;
  b.whole = false;
  constexpr std::size_t POSITIONS = 100000;
  std::size_t reads = 0, bases = 0, acgt = 0, gc = 0, n = 0, shortest = SIZE_MAX, longest = 0, scored = 0, q30 = 0;
  double qtotal = 0;
  std::vector<double> qsum;
  std::vector<std::size_t> qcount;
  std::map<std::size_t, std::size_t> lengths;
  const auto reading = each_read(request, b, [&](ReadRecord&& r) {
    const std::size_t L = r.bases.size();
    reads++;
    bases += L;
    shortest = std::min(shortest, L);
    longest = std::max(longest, L);
    lengths[L]++;
    const auto k = seq::base_counts(r.bases);
    acgt += k[0] + k[1] + k[2] + k[3];
    gc += k[1] + k[2];
    n += k[4];
    if (r.quality.size() == L) {
      const std::size_t P = std::min(L, POSITIONS);
      if (qsum.size() < P) { qsum.resize(P, 0); qcount.resize(P, 0); }
      for (std::size_t i = 0; i < L; i++) {
        const int q = seq::phred_of(r.quality[i]);
        qtotal += q;
        q30 += q >= 30;
        if (i < P) { qsum[i] += q; qcount[i]++; }
      }
      scored += L;
    }
  });
  if (reads == 0) throw Fault("no reads could be read from " + reading.source);
  // Quality along the read, in bins of positions when the reads are long.
  const std::size_t qbin = std::max<std::size_t>(1, (qsum.size() + 499) / 500);
  std::vector<double> quality;
  for (std::size_t i = 0; i < qsum.size(); i += qbin) {
    double s = 0, c = 0;
    for (std::size_t j = i; j < std::min(qsum.size(), i + qbin); j++) { s += qsum[j]; c += static_cast<double>(qcount[j]); }
    quality.push_back(c > 0 ? s / c : 0);
  }
  // Lengths, exactly when they span at most a thousand, else in about 200 bins.
  const std::size_t lbin = longest - shortest <= 1000 ? 1 : (longest + 200) / 200;
  std::map<std::size_t, std::size_t> binned;
  for (const auto& [len, c] : lengths) binned[lbin == 1 ? len : (len / lbin) * lbin + lbin / 2] += c;
  Json lens = Json::array();
  for (const auto& [len, c] : binned) lens.push(Json::array({len, c}));
  return Json::object({{"reads", reads}, {"bases", bases}, {"mean_length", static_cast<double>(bases) / static_cast<double>(reads)}, {"shortest", shortest}, {"longest", longest},
                       {"gc", acgt ? static_cast<double>(gc) / static_cast<double>(acgt) : 0.0}, {"n_rate", bases ? static_cast<double>(n) / static_cast<double>(bases) : 0.0},
                       {"qualities", scored > 0}, {"mean_quality", scored ? Json(detail::round_to(qtotal / static_cast<double>(scored), 100)) : Json()},
                       {"q30", scored ? Json(static_cast<double>(q30) / static_cast<double>(scored)) : Json()},
                       {"quality", numbers(quality, 2)}, {"quality_bin", qbin}, {"lengths", lens}, {"length_bin", lbin}, {"reading", reading.json()}});
}

// Reads placed on a genome by the domain's minimizer mapper (anchors, chains, base-level
// alignment; after Li 2018), in batches on the worker's threads: how many are mapped and at what
// quality, the depth along the genome, and the first records as SAM.
BENCH_OP("tool.reads.map") {
  detail::require_reads(request);
  const auto g = worker.genome(request);
  const bool long_reads = detail::long_reads_of(request);
  const auto max_reads = static_cast<std::size_t>(bounded(request, "max_reads", long_reads ? 200'000 : 2'000'000, 1, 50'000'000));
  const std::size_t L = g->length();
  if (L == 0) throw Fault("the genome is empty");
  double index_seconds = 0;
  const auto mapper = detail::mapper_of(worker, *g, long_reads, &index_seconds);
  std::vector<std::string> names;
  std::vector<std::size_t> lengths, offsets;
  for (const auto& c : g->contigs) {
    offsets.push_back(std::accumulate(lengths.begin(), lengths.end(), std::size_t{0}));
    names.push_back(c.name);
    lengths.push_back(c.codes.size());
  }
  const std::size_t B = std::min<std::size_t>(1000, L);
  std::vector<double> covered(B, 0);
  const auto bin_of = [&](std::size_t x) { return std::min(B - 1, x * B / L); };
  const auto edge = [&](std::size_t k) { return k * L / B; };
  std::size_t reads = 0, bases = 0, mapped = 0, mapq20 = 0, aligned = 0, sam_records = 0;
  constexpr std::size_t SAM_CAP = 200;
  std::ostringstream sam;
  sam << gmap::sam_header(names, lengths, "genome-workbench", "tool.reads.map");
  double map_seconds = 0;
  std::vector<ReadRecord> batch;
  const std::size_t per_batch = long_reads ? 4000 : 100000;
  const auto flush = [&] {
    if (batch.empty()) return;
    const auto t = std::chrono::steady_clock::now();
    const auto ms = mapper->map_all(worker.pool, batch, long_reads ? 1 : 64);
    map_seconds += detail::seconds_since(t);
    for (std::size_t i = 0; i < batch.size(); i++) {
      reads++;
      bases += batch[i].bases.size();
      const gmap::Mapping* m = gmap::primary_of(ms[i]);
      if (m) {
        mapped++;
        mapq20 += m->mapq >= 20;
        aligned += m->target_end - m->target_begin;
        // The bases it covers, added to the bins they fall in.
        const std::size_t a = offsets[m->contig] + m->target_begin, e = offsets[m->contig] + m->target_end;
        if (e > a)
          for (std::size_t k = bin_of(a) > 0 ? bin_of(a) - 1 : 0; k <= std::min(B - 1, bin_of(e - 1) + 1); k++) {
            const std::size_t lo = std::max(a, edge(k)), hi = std::min(e, edge(k + 1));
            if (hi > lo) covered[k] += static_cast<double>(hi - lo);
          }
      }
      if (sam_records < SAM_CAP) {
        sam << gmap::sam_line(m ? gmap::sam_record(*m, batch[i].name, batch[i].bases, batch[i].quality, names) : gmap::sam_unmapped(batch[i].name, batch[i].bases, batch[i].quality));
        sam_records++;
      }
    }
    batch.clear();
  };
  detail::Bounds b;
  b.bytes = std::size_t{8} << 30;
  b.records = max_reads;
  b.whole = false;
  const auto reading = each_read(request, b, [&](ReadRecord&& r) {
    batch.push_back(std::move(r));
    if (batch.size() >= per_batch) flush();
  });
  flush();
  if (reads == 0) throw Fault("no reads could be read from " + reading.source);
  Json coverage = Json::array();
  for (std::size_t k = 0; k < B; k++) coverage.push(detail::round_to(covered[k] / static_cast<double>(std::max<std::size_t>(1, edge(k + 1) - edge(k))), 100));
  return Json::object({{"reads", reads}, {"bases", bases}, {"mapped", mapped}, {"mapq20", mapq20}, {"unmapped", reads - mapped}, {"length", L},
                       {"mean_depth", static_cast<double>(aligned) / static_cast<double>(L)}, {"coverage", coverage}, {"contigs", names},
                       {"kind", long_reads ? "long" : "short"}, {"seconds", map_seconds}, {"index_seconds", index_seconds}, {"threads", worker.pool.size()},
                       {"sam", sam.str()}, {"sam_records", sam_records}, {"sam_truncated", reads > sam_records}, {"reading", reading.json()}});
}

// Short reads mapped to a genome and the variants of a diploid called from their pileups:
// genotype likelihoods (Li 2011), priors, the posterior genotype and its quality, filters; the
// calls as VCF. The reads are held, since the caller reads every base of every pileup.
BENCH_OP("tool.reads.call") {
  detail::require_reads(request);
  const auto g = worker.genome(request);
  const std::size_t L = g->length();
  if (L > 200'000'000) throw Fault("calling here takes genomes of up to 200 Mb");
  const auto max_reads = static_cast<std::size_t>(bounded(request, "max_reads", 2'000'000, 1, 20'000'000));
  detail::Bounds b;
  b.bytes = std::size_t{8} << 30;
  b.records = max_reads;
  b.whole = false;
  std::vector<ReadRecord> reads;
  const auto reading = each_read(request, b, [&](ReadRecord&& r) { reads.push_back(std::move(r)); });
  if (reads.empty()) throw Fault("no reads could be read from " + reading.source);
  double index_seconds = 0;
  const auto mapper = detail::mapper_of(worker, *g, false, &index_seconds);
  auto t = std::chrono::steady_clock::now();
  const auto mappings = mapper->map_all(worker.pool, reads, 64);
  const double map_seconds = detail::seconds_since(t);
  t = std::chrono::steady_clock::now();
  const auto calls = gvar::call_variants(worker.pool, detail::contigs_of(*g), mappings, gvar::read_views(reads));
  const double call_seconds = detail::seconds_since(t);
  std::vector<std::string> names;
  std::vector<std::size_t> lengths;
  for (const auto& c : g->contigs) {
    names.push_back(c.name);
    lengths.push_back(c.codes.size());
  }
  std::size_t mapped = 0, aligned = 0, passed = 0, snvs = 0, indels = 0;
  for (const auto& ms : mappings)
    if (const auto* m = gmap::primary_of(ms)) {
      mapped++;
      aligned += m->target_end - m->target_begin;
    }
  Json list = Json::array();
  std::vector<gvar::VcfRecord> records;
  for (const auto& c : calls) {
    passed += c.passed();
    (c.indel ? indels : snvs)++;
    if (records.size() < 20000) records.push_back(gvar::vcf_record(c, names));
    if (list.size() >= 5000) continue;
    std::string alt, filter;
    for (std::size_t i = 0; i < c.alts.size(); i++) alt += (i ? "," : "") + c.alts[i];
    for (std::size_t i = 0; i < c.filters.size(); i++) filter += (i ? ";" : "") + c.filters[i];
    list.push(Json::object({{"contig", names[c.contig]}, {"position", c.position}, {"ref", c.ref}, {"alt", alt},
                            {"genotype", std::to_string(c.genotype[0]) + "/" + std::to_string(c.genotype[1])}, {"gq", c.gq}, {"qual", detail::round_to(c.qual, 10)},
                            {"depth", c.depth}, {"filter", c.passed() ? std::string("PASS") : filter}, {"indel", c.indel}}));
  }
  const std::string vcf = gvar::write_vcf(gvar::caller_header(names, lengths, detail::sample_name(reading.source), "genome-workbench tool.reads.call"), records);
  return Json::object({{"calls", list}, {"called", calls.size()}, {"listed", list.size()}, {"passed", passed}, {"snvs", snvs}, {"indels", indels},
                       {"reads", reads.size()}, {"mapped", mapped}, {"length", L}, {"mean_depth", static_cast<double>(aligned) / static_cast<double>(std::max<std::size_t>(1, L))},
                       {"vcf", vcf}, {"vcf_records", records.size()},
                       {"seconds", Json::object({{"index", index_seconds}, {"map", map_seconds}, {"call", call_seconds}})}, {"threads", worker.pool.size()},
                       {"reading", reading.json()}});
}

// A genome assembled from its reads: their canonical k-mers counted a batch of reads at a time,
// the spectrum fitted (its error threshold, its peak, the genome's size), the solid k-mers joined
// in a de Bruijn graph and cleaned of tips, bubbles and weak paths; the contigs that remain.
BENCH_OP("tool.assembly.reads") {
  detail::require_reads(request);
  const int k = static_cast<int>(bounded(request, "k", 31, 11, 31));
  if (k % 2 == 0) throw Fault("k is odd, so that no k-mer is its own reverse complement");
  constexpr std::size_t DISTINCT = 150'000'000, BATCH = std::size_t{16} << 20;
  genomics::kmers::CounterOptions co;
  co.k = k;
  genomics::kmers::KmerCounter counter(co);
  std::vector<seq::Codes> batch;
  std::size_t in_batch = 0, reads = 0, bases = 0;
  bool full = false;
  const auto t0 = std::chrono::steady_clock::now();
  const auto flush = [&] {
    if (batch.empty()) return;
    counter.add(batch, worker.pool);
    batch.clear();
    in_batch = 0;
    full = counter.distinct() > DISTINCT;
  };
  detail::Bounds b;
  b.bytes = std::size_t{4} << 30;
  b.records = 100'000'000;
  b.whole = false;
  auto reading = each_read(request, b, [&](ReadRecord&& r) {
    reads++;
    bases += r.bases.size();
    in_batch += r.bases.size();
    batch.push_back(std::move(r.bases));
    if (in_batch >= BATCH) flush();
    return !full;
  });
  flush();
  if (full) reading.stopped = "reading stopped at " + std::to_string(counter.distinct()) + " distinct k-mers, the most a worker holds here";
  if (reads == 0) throw Fault("no reads could be read from " + reading.source);
  const double counting = detail::seconds_since(t0);
  counter.release_buffers();
  asmb::AssemblyOptions o;
  o.k = k;
  o.counter.k = k;
  o.cleaning.keep_graphs = false;
  const auto a = asmb::assemble(counter, worker.pool, o, 0);
  Json out = Json::object({{"reads", reads}, {"bases", bases}, {"k", k}});
  detail::assembly_report(a, out);
  Json seconds = detail::assembly_seconds(a);
  seconds.set("read_and_count", counting);
  out.set("seconds", seconds);
  // The contigs, longest first, as FASTA, up to 16 MB of it.
  std::string fasta;
  std::size_t written = 0;
  for (std::size_t i = 0; i < a.contigs.size() && fasta.size() < (std::size_t{16} << 20); i++, written++) {
    fasta += ">contig_" + std::to_string(i + 1) + " length=" + std::to_string(a.contigs[i].size()) + "\n";
    const std::string s = seq::decode_dna(a.contigs[i]);
    for (std::size_t j = 0; j < s.size(); j += 70) fasta += s.substr(j, 70) + "\n";
  }
  out.set("contigs_fasta", fasta);
  out.set("contigs_written", written);
  out.set("threads", worker.pool.size());
  out.set("reading", reading.json());
  return out;
}

// Transcripts quantified from single-end RNA-seq reads: pseudoaligned to the transcriptome's index
// of k-mers by their equivalence classes, effective lengths from a normal law of fragment lengths
// (mean and deviation given, as kallisto takes them for single-end reads), and each transcript's
// abundance by EM (Bray, Pimentel, Melsted and Pachter 2016).
BENCH_OP("tool.expr.quant") {
  detail::require_reads(request);
  const std::string tpath = str(request, "transcripts_file");
  if (tpath.empty()) throw Fault("choose a FASTA file of transcripts");
  detail::Reading transcripts_reading;
  transcripts_reading.source = detail::file_name(tpath);
  transcripts_reading.format = "fasta";
  std::error_code ec;
  transcripts_reading.file_bytes = static_cast<std::size_t>(std::filesystem::file_size(tpath, ec));
  if (transcripts_reading.file_bytes > (std::size_t{1} << 30)) throw Fault(transcripts_reading.source + " has " + detail::megabytes(transcripts_reading.file_bytes) + "; a transcriptome here has up to 1 GB");
  const std::string text = detail::read_text(tpath);
  if (detail::sniff(text) != "fasta") throw Fault(transcripts_reading.source + " is not FASTA");
  transcripts_reading.bytes = text.size();
  const auto tx = gex::read_transcriptome(text, &transcripts_reading.problems);
  if (tx.size() == 0) throw Fault("no transcripts could be read from " + transcripts_reading.source);
  transcripts_reading.records = tx.size();
  for (const auto& t : tx.transcripts) transcripts_reading.letters += t.sequence.size();
  const double mean = num(request, "fragment", 200), sd = num(request, "fragment_sd", 20);
  if (!(mean >= 20 && mean <= 2000) || !(sd > 0 && sd <= 1000)) throw Fault("a mean fragment length of 20 to 2000 bases, and a deviation of up to 1000");
  const auto max_reads = static_cast<std::size_t>(bounded(request, "max_reads", 10'000'000, 1, 50'000'000));
  detail::Bounds b;
  b.bytes = std::size_t{8} << 30;
  b.records = max_reads;
  b.whole = false;
  std::vector<seq::Codes> reads;
  const auto t0 = std::chrono::steady_clock::now();
  const auto reading = each_read(request, b, [&](ReadRecord&& r) { reads.push_back(std::move(r.bases)); });
  const double read_seconds = detail::seconds_since(t0);
  if (reads.empty()) throw Fault("no reads could be read from " + reading.source);
  const auto d = gex::FragmentLengths::normal(mean, sd);
  const auto q = detail::quantify_sample(tx, d, worker.pool, [&](const gex::Index& index) { return gex::pseudoalign(index, reads, {}, worker.pool); });
  Json rows = Json::array();
  for (std::size_t i = 0; i < tx.size() && rows.size() < 100000; i++)
    rows.push(Json::object({{"name", tx.transcripts[i].name}, {"gene", tx.genes[tx.gene_of[i]]}, {"length", tx.transcripts[i].sequence.size()},
                            {"effective", detail::round_to(q.effective[i], 10)}, {"count", detail::round_to(q.abundances.counts[i], 100)},
                            {"tpm", detail::round_to(q.abundances.tpm[i], 1000)}}));
  return Json::object({{"transcripts", rows}, {"genes", tx.genes.size()}, {"fragments", reads.size()}, {"aligned", q.pa.aligned}, {"classes", q.pa.classes.size()},
                       {"rounds", q.abundances.rounds}, {"converged", q.abundances.converged}, {"kmers", q.kmers}, {"index_bytes", q.index_bytes},
                       {"fragment", Json::object({{"mean", mean}, {"sd", sd}})},
                       {"seconds", Json::object({{"read", read_seconds}, {"index", q.index_seconds}, {"pseudoalign", q.pseudoalign_seconds}, {"em", q.em_seconds}})},
                       {"threads", worker.pool.size()}, {"transcripts_reading", transcripts_reading.json()}, {"reading", reading.json()}});
}

}  // namespace workbench
