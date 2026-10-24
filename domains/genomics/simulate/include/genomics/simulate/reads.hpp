// Reads simulated from a genome or from a diploid individual, with their truth: where each
// came from on the reference, on which strand, from which haplotype, and how many errors the
// sequencer made in it.
//
// - Short reads (Illumina-like): a fixed length; substitution errors whose probability rises
//   linearly from the first base to the last (the decline of quality along a read, Dohm et
//   al. 2008, Nucleic Acids Research 36, e105), with a small fraction of insertions and
//   deletions of one base; each base's quality is the Phred value of its error probability,
//   with a little noise, so that qualities are informative as a sequencer's are.
// - Long reads (nanopore- or PacBio-like): lengths from a gamma law of the given mean and
//   shape, at least `min_length`; errors at a rate per base, divided among insertions,
//   deletions and substitutions (indels dominate, as in Weirather et al. 2017, F1000Research
//   6, 100); a flat quality.
//
// Read i is drawn from its own stream, keyed by the seed and i, so that reads a to b are the
// same whether made alone or as part of reads 0 to n: a shard of a simulated sample is two
// numbers. A genome of several contigs gives reads from each in proportion to its length;
// reads never cross the end of a contig.

#pragma once

#include <genomics/sequence/fasta.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace genomics::simulate {

struct Read {
  std::string name;
  Codes bases;
  std::string quality;          // Phred + 33, one per base
  std::size_t contig = 0;       // of the reference
  std::size_t begin = 0, end = 0;  // the reference interval it covers (truth)
  bool reverse = false;         // read from the reverse strand
  int haplotype = 0;
  std::size_t errors = 0;       // substitutions, insertions and deletions made
};

struct ShortReadModel {
  std::size_t length = 150;
  double error_first = 0.001, error_last = 0.01;  // substitution probability at the first and last base
  double indel = 0.0002;                          // per base, half insertions, half deletions
};

struct LongReadModel {
  double mean_length = 8000, shape = 2;  // the gamma law of lengths
  std::size_t min_length = 500;
  double error = 0.08;                   // per base
  double insertion = 0.35, deletion = 0.45;  // fractions of the errors; the rest substitutions
  int quality = 12;                      // the Phred value written for every base
};

// What reads are drawn from: one or two copies (haplotypes) of each contig, with their maps to
// the reference.
struct Source {
  std::vector<std::vector<Codes>> copies;     // copies[h][contig]
  std::vector<std::vector<Liftover>> maps;    // maps[h][contig]

  static Source haploid(const std::vector<Codes>& contigs) {
    Source s;
    s.copies.push_back(contigs);
    s.maps.emplace_back();
    for (const auto& c : contigs) s.maps.back().push_back(Liftover::identity(c.size()));
    return s;
  }
  // A diploid genome of several contigs, one Diploid per contig.
  static Source diploid(const std::vector<Diploid>& contigs) {
    Source s;
    s.copies.resize(2);
    s.maps.resize(2);
    for (int h = 0; h < 2; h++)
      for (const auto& d : contigs) {
        s.copies[h].push_back(d.haplotype[h]);
        s.maps[h].push_back(d.liftover[h]);
      }
    return s;
  }
};

namespace detail {

// A contig drawn in proportion to its length, then a start within it for a fragment of `len`.
inline bool place(calx::random::Stream& r, const std::vector<Codes>& contigs, std::size_t len, std::size_t& contig, std::size_t& start) {
  std::size_t total = 0;
  for (const auto& c : contigs) total += c.size() >= len ? c.size() - len + 1 : 0;
  if (total == 0) return false;
  std::size_t u = r.below(total);
  for (contig = 0; contig < contigs.size(); contig++) {
    const std::size_t room = contigs[contig].size() >= len ? contigs[contig].size() - len + 1 : 0;
    if (u < room) { start = u; return true; }
    u -= room;
  }
  return false;
}

}  // namespace detail

// Reads first, first + 1, ..., first + count − 1 of a short-read sample.
inline std::vector<Read> short_reads(const Source& source, std::size_t count, const ShortReadModel& model, std::uint64_t seed, std::size_t first = 0) {
  std::vector<Read> out;
  out.reserve(count);
  const std::size_t L = model.length;
  for (std::size_t i = first; i < first + count; i++) {
    calx::random::Stream r(seed, READS, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i >> 32));
    Read read;
    read.name = "read" + std::to_string(i);
    read.haplotype = static_cast<int>(r.below(source.copies.size()));
    const auto& contigs = source.copies[static_cast<std::size_t>(read.haplotype)];
    // A fragment a little longer than the read, for the deletions to draw on.
    const std::size_t fragment = L + 8;
    std::size_t start = 0;
    if (!detail::place(r, contigs, fragment, read.contig, start)) continue;
    const Codes& src = contigs[read.contig];
    read.reverse = r.bernoulli(0.5);
    // The template, read 5' to 3' on the strand sequenced.
    Codes t(src.begin() + static_cast<std::ptrdiff_t>(start), src.begin() + static_cast<std::ptrdiff_t>(start + fragment));
    if (read.reverse) t = sequence::reverse_complement(t);
    std::size_t j = 0;  // the next base of the template
    while (read.bases.size() < L && j < t.size()) {
      const double x = L > 1 ? static_cast<double>(read.bases.size()) / static_cast<double>(L - 1) : 0;
      const double perr = model.error_first + (model.error_last - model.error_first) * x;
      const double u = r.uniform();
      int q = sequence::phred_from_probability(perr) + static_cast<int>(std::lround(r.normal(0, 1.5)));
      q = std::clamp(q, 2, 41);
      if (u < model.indel / 2) {  // an inserted base
        read.bases.push_back(static_cast<Code>(r.below(4)));
        read.quality += sequence::phred_char(q);
        read.errors++;
        continue;
      }
      if (u < model.indel) {  // a deleted base
        j++;
        read.errors++;
        continue;
      }
      Code b = t[j++];
      if (r.uniform() < perr) { b = substitute(r, b); read.errors++; }
      read.bases.push_back(b);
      read.quality += sequence::phred_char(q);
    }
    // The bases of the template it used, back on the forward strand of the haplotype, then the reference.
    const std::size_t used = j;
    const std::size_t hb = read.reverse ? start + fragment - used : start, he = read.reverse ? start + fragment : start + used;
    const Liftover& map = source.maps[static_cast<std::size_t>(read.haplotype)][read.contig];
    read.begin = map.to_reference(hb);
    read.end = std::max(read.begin, map.to_reference(he - 1) + 1);
    out.push_back(std::move(read));
  }
  return out;
}

// Reads first, ..., first + count − 1 of a long-read sample.
inline std::vector<Read> long_reads(const Source& source, std::size_t count, const LongReadModel& model, std::uint64_t seed, std::size_t first = 0) {
  std::vector<Read> out;
  out.reserve(count);
  const double pin = model.error * model.insertion, pdel = model.error * model.deletion, psub = model.error - pin - pdel;
  for (std::size_t i = first; i < first + count; i++) {
    calx::random::Stream r(seed, LONG_READS, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i >> 32));
    Read read;
    read.name = "long" + std::to_string(i);
    read.haplotype = static_cast<int>(r.below(source.copies.size()));
    const auto& contigs = source.copies[static_cast<std::size_t>(read.haplotype)];
    std::size_t longest = 0;
    for (const auto& c : contigs) longest = std::max(longest, c.size());
    std::size_t len = static_cast<std::size_t>(r.gamma(model.shape, model.mean_length / model.shape));
    len = std::clamp(len, model.min_length, longest);
    std::size_t start = 0;
    if (!detail::place(r, contigs, len, read.contig, start)) continue;
    const Codes& src = contigs[read.contig];
    read.reverse = r.bernoulli(0.5);
    Codes t(src.begin() + static_cast<std::ptrdiff_t>(start), src.begin() + static_cast<std::ptrdiff_t>(start + len));
    if (read.reverse) t = sequence::reverse_complement(t);
    read.bases.reserve(len + len / 10);
    for (std::size_t j = 0; j < t.size();) {
      const double u = r.uniform();
      if (u < pin) { read.bases.push_back(static_cast<Code>(r.below(4))); read.errors++; continue; }
      if (u < pin + pdel) { j++; read.errors++; continue; }
      Code b = t[j++];
      if (u < pin + pdel + psub) { b = substitute(r, b); read.errors++; }
      read.bases.push_back(b);
    }
    read.quality.assign(read.bases.size(), sequence::phred_char(model.quality));
    const Liftover& map = source.maps[static_cast<std::size_t>(read.haplotype)][read.contig];
    read.begin = map.to_reference(start);
    read.end = std::max(read.begin, map.to_reference(start + len - 1) + 1);
    out.push_back(std::move(read));
  }
  return out;
}

// Reads as FASTQ records.
inline std::vector<sequence::Record> fastq_records(const std::vector<Read>& reads) {
  std::vector<sequence::Record> out;
  out.reserve(reads.size());
  for (const auto& r : reads) out.push_back({r.name, "", sequence::decode_dna(r.bases), r.quality});
  return out;
}

}  // namespace genomics::simulate
