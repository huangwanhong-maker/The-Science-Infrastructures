// The calling of single-nucleotide variants and short indels of a diploid from its mapped reads:
// candidates from the pileup, genotypes by their posterior (genotype.hpp), and filters on depth
// and strand bias.
//
// - SNVs. At each reference base the pileup's bases (quality ≥ 13) are counted; the most
//   frequent other base is the candidate alternate allele when it has at least `min_alt` reads
//   (2) and a fraction `min_alt_fraction` (0.1) of the depth. A second alternate with at least
//   `min_alt` reads and a fifth of the depth makes the site triallelic: the six genotypes over
//   {R, A₁, A₂} are weighed, and the site is reduced to the alleles of the genotype of largest
//   posterior: to {R, A₁} or {R, A₂} (and weighed again over its three genotypes) unless it is
//   A₁/A₂, which is reported as one multi-allelic record (ALT = A₁,A₂, GT 1/2).
// - Indels. At each anchor, the reads' indels (left-aligned in the pileup) are grouped by their
//   allele; the most frequent is the candidate, carried by k reads. An indel in a repeat can be
//   placed anywhere from its left-aligned to its right-aligned position, and only a read that
//   spans that whole interval of ambiguity (from the anchor to the base after the rightmost
//   placement) can show the reference allele there: a read ending inside the repeat looks the
//   same with the indel or without it. The reference support is the reads spanning it that do
//   not carry the indel, d = k + those, weighed by the counting model of genotype.hpp (ε = 0.01,
//   θ = 1.25 · 10⁻⁴). A second indel allele at the same anchor is not reported (a site carrying
//   two different indels is rare and is left to a haplotype caller). REF and ALT are written
//   anchored on the base before, and normalised (normalize.hpp).
// - A site is called when its genotype of largest posterior is not 0/0 and QUAL ≥ `min_qual`
//   (20). Filters (VCF's FILTER): LowDepth below `min_depth` reads (5); HighDepth above
//   `max_depth` (0: none); StrandBias when FS > 60 for an SNV, 200 for an indel (GATK's hard
//   filters; Van der Auwera et al. 2013, Current Protocols in Bioinformatics 43, 11.10).
// - Calls as VCF: caller_header (the contigs; INFO DP and FS; FILTER LowDepth, HighDepth,
//   StrandBias; FORMAT GT:GQ:DP:AD) and vcf_record, for vcf.hpp's writer.
//
// The genome is called in windows of `window` bases (32 768), each a pileup of the mappings that
// overlap it, on calx::scheduler::parallel_for; each window's calls go to its own slot, and the
// slots are joined in order, so that the calls are the same on any number of threads.
//
// Complexity: O(Σ ℓ_r) for the pileups (every aligned base once), O(d) per candidate site.
// Limits: a pileup caller, site by site: an indel and the SNVs near it are weighed apart, reads
// are not realigned to candidate haplotypes (as GATK's HaplotypeCaller and DeepVariant do), and
// indels longer than the reads hold are not seen.

#pragma once

#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/variant/genotype.hpp>
#include <genomics/variant/normalize.hpp>
#include <genomics/variant/pileup.hpp>
#include <genomics/variant/vcf.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::variant {

struct CallerOptions {
  PileupOptions pileup;
  double theta = 1e-3, theta_indel = 1.25e-4;
  double indel_error = 0.01;
  std::size_t min_alt = 2;
  double min_alt_fraction = 0.1;
  double second_allele_fraction = 0.2;
  double min_qual = 20;
  std::size_t min_depth = 5, max_depth = 0;
  double max_strand_bias = 60, max_strand_bias_indel = 200;
  std::size_t window = 1 << 15;
};

struct Call {
  std::uint32_t contig = 0;
  std::size_t position = 0;           // 0-based, of REF's first base
  std::string ref;
  std::vector<std::string> alts;
  double qual = 0;
  std::vector<std::string> filters;   // empty: PASS
  std::size_t depth = 0;
  std::vector<std::size_t> allele_depths;  // REF, then each ALT
  std::array<int, 2> genotype{0, 0};  // allele indices, a ≤ b
  int gq = 0;
  std::vector<int> pl;
  double strand_bias = 0;             // FS
  bool indel = false;
  bool passed() const { return filters.empty(); }
};

namespace detail {

inline void apply_filters(Call& c, const CallerOptions& o) {
  if (c.depth < o.min_depth) c.filters.push_back("LowDepth");
  if (o.max_depth && c.depth > o.max_depth) c.filters.push_back("HighDepth");
  if (c.strand_bias > (c.indel ? o.max_strand_bias_indel : o.max_strand_bias)) c.filters.push_back("StrandBias");
}

inline void snv_at(const Pileup& pu, const Codes& reference, std::size_t p, const CallerOptions& o, std::vector<BaseCall>& scratch, std::vector<Call>& out) {
  const sequence::Code ref = reference[p];
  if (ref > 3) return;
  const auto obs = pu.at(p);
  if (obs.size() < o.min_alt) return;
  std::array<std::size_t, 5> count{};
  std::array<std::size_t, 5> forward{};
  for (const auto& b : obs) {
    count[std::min<std::size_t>(b.base, 4)]++;
    forward[std::min<std::size_t>(b.base, 4)] += !b.reverse;
  }
  sequence::Code a1 = 4, a2 = 4;
  for (sequence::Code b = 0; b < 4; b++) {
    if (b == ref) continue;
    if (a1 == 4 || count[b] > count[a1]) { a2 = a1; a1 = b; }
    else if (a2 == 4 || count[b] > count[a2]) a2 = b;
  }
  const double depth = static_cast<double>(obs.size());
  if (a1 == 4 || count[a1] < o.min_alt || static_cast<double>(count[a1]) < o.min_alt_fraction * depth) return;
  scratch.clear();
  for (const auto& b : obs) scratch.push_back({b.base, b.quality});
  std::vector<sequence::Code> alleles{ref, a1};
  if (a2 != 4 && count[a2] >= o.min_alt && static_cast<double>(count[a2]) >= o.second_allele_fraction * depth) alleles.push_back(a2);
  GenotypeCall g = posterior(log10_likelihoods(scratch, alleles), genotype_prior(alleles.size(), o.theta));
  if (alleles.size() == 3 && !(g.allele_a == 1 && g.allele_b == 2)) {
    // Reduced to the alleles of the best genotype.
    const sequence::Code keep = (g.allele_a == 2 || g.allele_b == 2) ? a2 : a1;
    alleles = {ref, keep};
    g = posterior(log10_likelihoods(scratch, alleles), genotype_prior(2, o.theta));
  }
  if (g.best == 0 || g.qual < o.min_qual) return;
  Call c;
  c.contig = pu.contig;
  c.position = p;
  c.ref = std::string(1, sequence::dna_letter(ref));
  for (std::size_t a = 1; a < alleles.size(); a++) c.alts.push_back(std::string(1, sequence::dna_letter(alleles[a])));
  c.qual = g.qual;
  c.depth = obs.size();
  for (const auto a : alleles) c.allele_depths.push_back(count[a]);
  c.genotype = {static_cast<int>(g.allele_a), static_cast<int>(g.allele_b)};
  c.gq = g.gq;
  c.pl = g.pl;
  std::size_t alt_f = 0, alt_r = 0;
  for (std::size_t a = 1; a < alleles.size(); a++) {
    alt_f += forward[alleles[a]];
    alt_r += count[alleles[a]] - forward[alleles[a]];
  }
  c.strand_bias = strand_bias(forward[ref], count[ref] - forward[ref], alt_f, alt_r);
  apply_filters(c, o);
  out.push_back(std::move(c));
}

// The indel candidates at one anchor: indels[first, last) of the pileup.
inline void indel_at(const Pileup& pu, const Codes& reference, std::size_t first, std::size_t last, const CallerOptions& o, std::vector<Call>& out) {
  // The most frequent allele; ties to the first in the pileup's order (shorter, then by bases).
  std::size_t best = first, best_n = 0;
  for (std::size_t i = first; i < last;) {
    std::size_t j = i;
    while (j < last && pu.indels[j].length == pu.indels[i].length && pu.indels[j].bases == pu.indels[i].bases) j++;
    if (j - i > best_n) { best = i; best_n = j - i; }
    i = j;
  }
  const IndelObservation& e = pu.indels[best];
  std::size_t kf = 0, k = 0;
  std::vector<std::uint32_t> carriers;
  for (std::size_t i = best; i < last && pu.indels[i].length == e.length && pu.indels[i].bases == e.bases; i++) {
    k++;
    kf += !pu.indels[i].reverse;
    carriers.push_back(pu.indels[i].read);
  }
  std::sort(carriers.begin(), carriers.end());
  const std::size_t p = e.position, size = reference.size();
  // The interval a read must span to tell the indel from the reference: from the anchor to the
  // base after the indel's rightmost equivalent placement.
  std::size_t informative = 0;
  if (e.length < 0) {
    const auto L = static_cast<std::size_t>(-e.length);
    std::size_t a = p + 1;
    while (a + L < size && reference[a] == reference[a + L]) a++;
    informative = a + L + 1;
  } else {
    std::size_t shift = 0;
    while (p + 1 + shift < size && sequence::dna_code(e.bases[shift % e.bases.size()]) == reference[p + 1 + shift]) shift++;
    informative = p + shift + 2;
  }
  if (informative > size) return;
  std::size_t ref_f = 0, ref_r = 0;
  pu.for_each_spanning(p, informative, [&](const ReadSpan& s) {
    if (std::binary_search(carriers.begin(), carriers.end(), s.read)) return;
    (s.reverse ? ref_r : ref_f)++;
  });
  const std::size_t d = ref_f + ref_r + k;
  if (k < o.min_alt || static_cast<double>(k) < o.min_alt_fraction * static_cast<double>(d)) return;
  const auto l = log10_likelihoods(d - k, k, o.indel_error);
  const GenotypeCall g = posterior(l, genotype_prior(2, o.theta_indel));
  if (g.best == 0 || g.qual < o.min_qual) return;
  std::string ref(1, sequence::dna_letter(reference[p])), alt = ref;
  if (e.length > 0) alt += e.bases;
  else ref += e.bases;
  Normalized n = normalize(p, ref, alt, reference);
  Call c;
  c.contig = pu.contig;
  c.position = n.position;
  c.ref = std::move(n.ref);
  c.alts = std::move(n.alts);
  c.indel = true;
  c.qual = g.qual;
  c.depth = d;
  c.allele_depths = {d - k, k};
  c.genotype = {static_cast<int>(g.allele_a), static_cast<int>(g.allele_b)};
  c.gq = g.gq;
  c.pl = g.pl;
  c.strand_bias = strand_bias(ref_f, ref_r, kf, k - kf);
  apply_filters(c, o);
  out.push_back(std::move(c));
}

}  // namespace detail

// The calls in [begin, end) of one contig from the mappings (any order) of `reads`.
inline std::vector<Call> call_region(const Codes& reference, std::uint32_t contig, std::size_t begin, std::size_t end, std::span<const mapping::Mapping> mappings,
                                     std::span<const ReadView> reads, const CallerOptions& o = {}) {
  const Pileup pu = pileup(reference, contig, begin, end, mappings, reads, o.pileup);
  std::vector<Call> out;
  std::vector<BaseCall> scratch;
  std::size_t next = 0;  // the first indel not yet visited
  for (std::size_t p = begin; p < end; p++) {
    detail::snv_at(pu, reference, p, o, scratch, out);
    std::size_t last = next;
    while (last < pu.indels.size() && pu.indels[last].position == p) last++;
    if (last > next) detail::indel_at(pu, reference, next, last, o, out);
    next = last;
  }
  // Normalisation may move an indel before an SNV called at an earlier anchor.
  std::stable_sort(out.begin(), out.end(), [](const Call& a, const Call& b) { return a.position < b.position; });
  return out;
}

// The calls of a whole genome from every read's mappings (Mapper::map_all's output), in order of
// contig and position, the same on any number of threads.
inline std::vector<Call> call_variants(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const std::vector<std::vector<mapping::Mapping>>& mappings,
                                       std::span<const ReadView> reads, const CallerOptions& o = {}) {
  if (o.window == 0) throw std::invalid_argument("genomics::variant::call_variants: a window of at least one base");
  std::vector<mapping::Mapping> flat;
  for (const auto& ms : mappings)
    for (const auto& m : ms)
      if (m.primary && m.mapq >= o.pileup.min_mapping_quality) flat.push_back(m);
  std::sort(flat.begin(), flat.end(), [](const mapping::Mapping& a, const mapping::Mapping& b) {
    if (a.contig != b.contig) return a.contig < b.contig;
    if (a.target_begin != b.target_begin) return a.target_begin < b.target_begin;
    return a.read < b.read;
  });
  std::size_t longest = 0;
  for (const auto& m : flat) longest = std::max(longest, m.target_end - m.target_begin);
  struct Window { std::uint32_t contig; std::size_t begin, end; };
  std::vector<Window> windows;
  for (std::uint32_t c = 0; c < contigs.size(); c++)
    for (std::size_t b = 0; b < contigs[c].size(); b += o.window) windows.push_back({c, b, std::min(contigs[c].size(), b + o.window)});
  std::vector<std::vector<Call>> slots(windows.size());
  calx::scheduler::parallel_for(pool, windows.size(), 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
    for (std::size_t w = lo; w < hi; w++) {
      const Window& x = windows[w];
      const auto key = [&](std::uint32_t contig, std::size_t pos) {
        return std::lower_bound(flat.begin(), flat.end(), std::pair{contig, pos}, [](const mapping::Mapping& m, const std::pair<std::uint32_t, std::size_t>& k) {
          return m.contig != k.first ? m.contig < k.first : m.target_begin < k.second;
        });
      };
      const auto a = key(x.contig, x.begin > longest ? x.begin - longest : 0), b = key(x.contig, x.end);
      const std::span<const mapping::Mapping> part(flat.data() + (a - flat.begin()), static_cast<std::size_t>(b - a));
      slots[w] = call_region(contigs[x.contig], x.contig, x.begin, x.end, part, reads, o);
    }
  });
  std::vector<Call> out;
  for (auto& s : slots)
    for (auto& c : s) out.push_back(std::move(c));
  std::stable_sort(out.begin(), out.end(), [](const Call& a, const Call& b) { return a.contig != b.contig ? a.contig < b.contig : a.position < b.position; });
  return out;
}

// ---- calls as VCF ---------------------------------------------------------------------------------
// The header of the caller's VCF: its contigs, INFO, FILTER and FORMAT definitions, one sample.
inline VcfHeader caller_header(const std::vector<std::string>& names, const std::vector<std::size_t>& lengths, std::string sample, std::string source = "genomics::variant") {
  if (names.size() != lengths.size()) throw std::invalid_argument("genomics::variant::caller_header: one length for each contig");
  VcfHeader h;
  h.meta.push_back("source=" + source);
  for (std::size_t i = 0; i < names.size(); i++) {
    h.meta.push_back("contig=<ID=" + names[i] + ",length=" + std::to_string(lengths[i]) + ">");
    h.contigs.push_back({names[i], lengths[i]});
  }
  h.meta.push_back("INFO=<ID=DP,Number=1,Type=Integer,Description=\"Reads at the site after the base and mapping quality filters\">");
  h.meta.push_back("INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value of Fisher's exact test of strand bias\">");
  h.meta.push_back("FILTER=<ID=LowDepth,Description=\"Fewer reads than the minimum depth\">");
  h.meta.push_back("FILTER=<ID=HighDepth,Description=\"More reads than the maximum depth\">");
  h.meta.push_back("FILTER=<ID=StrandBias,Description=\"FS above its limit (60 for SNVs, 200 for indels)\">");
  h.meta.push_back("FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  h.meta.push_back("FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality, -10 log10 P(genotype wrong)\">");
  h.meta.push_back("FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">");
  h.meta.push_back("FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth of each allele\">");
  h.samples.push_back(std::move(sample));
  return h;
}

// The record of a call, its contig named from `names`.
inline VcfRecord vcf_record(const Call& c, const std::vector<std::string>& names) {
  if (c.contig >= names.size()) throw std::invalid_argument("genomics::variant::vcf_record: a contig without a name");
  VcfRecord r;
  r.chrom = names[c.contig];
  r.position = c.position;
  r.ref = c.ref;
  r.alts = c.alts;
  r.qual = std::round(c.qual * 100) / 100;
  r.filters = c.passed() ? std::vector<std::string>{"PASS"} : c.filters;
  r.info = {{"DP", std::to_string(c.depth)}, {"FS", detail::qual_string(c.strand_bias)}};
  r.format = {"GT", "GQ", "DP", "AD"};
  std::string ad;
  for (std::size_t i = 0; i < c.allele_depths.size(); i++) ad += (i ? "," : "") + std::to_string(c.allele_depths[i]);
  r.samples = {{genotype_string({{c.genotype[0], c.genotype[1]}, false}), std::to_string(c.gq), std::to_string(c.depth), ad}};
  return r;
}

}  // namespace genomics::variant
