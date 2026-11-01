// A tumour derived from a normal diploid, with the truth of everything it carries: its clones and
// their cancer cell fractions, its somatic single-base substitutions and short indels drawn from
// mutational signatures, its allele-specific copy-number profile, its microsatellites that
// slipped, its purity; and reads of the tumour sample (tumour cells mixed with normal cells) and
// of the matched normal.
//
// - Clones. The tumour's cells descend from one cell, the trunk (clone 0); every other clone has a
//   parent of smaller index and a cancer cell fraction φ, the fraction of tumour cells descending
//   from it (itself included), so φ(trunk) = 1 and, since the descendants of a clone's children
//   are disjoint sets of its own descendants, the sum rule Σ_children φ ≤ φ(parent) holds (El-Kebir
//   et al. 2015, Bioinformatics 31, i62–i70, the sum condition; checked on entry). The fraction of
//   tumour cells that are of clone j and of none of its descendants is φ_j − Σ_children φ.
// - Copy number. Every tumour cell carries one allele-specific profile (the copy-number changes
//   are clonal, made in the trunk): a step function of the reference giving the copies of each
//   parental haplotype (simulate::CopySegment of structural.hpp), 1 + 1 where nothing changed;
//   gains, losses, copy-neutral loss of heterozygosity (2 + 0) and amplifications, and with a
//   whole-genome doubling every state doubled (random_copy_number). Normal cells are 1 + 1.
// - Mutations and their multiplicity. A mutation of the trunk is early (made in the diploid
//   ancestor, before the copy-number changes; a fraction `early` of them) or late. An early
//   mutation lies on a haplotype h chosen with probability 1/2 and is carried by every copy of h
//   the profile has there (multiplicity m = n_h; lost when n_h = 0); a late one, and every
//   mutation of another clone, hits one copy of the genome, on h with probability n_h/(n_A + n_B),
//   and has m = 1 (McGranahan et al. 2015, Sci. Transl. Med. 7, 283ra54: early mutations are
//   clonal and at the copies of an amplified or doubled allele, late ones at one copy). The
//   expected fraction of reads carrying it in the tumour sample is then
//     f = ρ φ m / (2(1 − ρ) + ρ (n_A + n_B))
//   for purity ρ (theory/the-somatic-genome.md, Proposition 2.2).
// - Signatures. A clone's substitutions are drawn from a mixture of mutational signatures, each a
//   profile over the 96 pyrimidine-centred classes (sequence/contexts.hpp): a signature by the
//   clone's exposures, a class by the signature's profile, then a position uniformly among those
//   whose trinucleotide (on either strand) is the class's, by rejection. The catalogue of the
//   truth therefore follows the mixture whatever the genome's trinucleotide composition, which is
//   how signatures are defined (Alexandrov et al. 2013, Nature 500, 415–421). COSMIC's published
//   profiles are not fetched: their terms forbid redistribution and display on free platforms
//   (somatic/README.md, §3), so `synthetic_signatures` gives profiles constructed here from the
//   mechanisms the literature describes, labelled synthetic.
// - Indels: uniform positions, half deletions and half insertions, lengths geometric of mean
//   `indel_mean` at most `indel_max`, anchored as VCF anchors them (mutate.hpp). Microsatellite
//   instability: each site of `msi.sites` (found in the reference by the caller, e.g.
//   somatic::find_microsatellites) slips in the trunk with probability `msi.rate`, losing one unit
//   (probability `msi.contraction`) or gaining one, written as a deletion or insertion of the unit
//   anchored on the base before the repeat. The sites are kept free of other somatic events, and
//   so are germline variants and every other variant, by `spacing` bases.
// - Reads. A tumour read (tumour_reads) starts at a reference position x drawn uniformly and kept
//   with probability T(x)/max T, T(x) = 2(1 − ρ) + ρ (n_A(x) + n_B(x)) the DNA there, so that depth
//   follows copy number; its molecule is a normal cell's haplotype (weight 1 − ρ each), or copy c
//   of haplotype h of a tumour cell of clone j (weight ρ · (fraction of j) · [c < n_h(x)]); the
//   molecule's bases are the germline haplotype with the somatic variants that molecule carries
//   (clone j descends from the variant's clone, same haplotype, c < m); the read is sequenced from
//   it by reads.hpp's short-read model (detail::sequence_short, one order of draws). Read i comes
//   from its own stream, keyed by the seed, TUMOUR_READS and i, so a shard is two numbers. The
//   normal is mutate_diploid's diploid, read by short_reads (normal_source).
//
// Every draw comes from calx::random::Stream keyed by the seed and a purpose (SOMATIC for the
// mutations, COPY_PROFILE for random profiles, TUMOUR_READS for reads), so a tumour is the same on
// any machine. Limits: copy-number changes are clonal (no subclonal copy number, as Battenberg
// would fit); a mutation's multiplicity is n_h or 1, never a value between (a mutation between two
// successive gains); no structural rearrangements (simulate/structural.hpp has them); no
// doublet substitutions; reads are single-end, without GC bias.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/contexts.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>
#include <genomics/simulate/structural.hpp>

#include <calx/random/multinomial.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace genomics::simulate {

// The purposes of this file's streams.
inline constexpr std::uint32_t SOMATIC = 0x534f4d41, COPY_PROFILE = 0x434e5052, TUMOUR_READS = 0x54554d52;  // "SOMA", "CNPR", "TUMR"

struct CloneSpec {
  int parent = -1;                // −1: the trunk
  double ccf = 1;                 // the fraction of tumour cells descending from it
  std::size_t snvs = 0, indels = 0;
  std::vector<double> exposures;  // weights over TumourSpec::signatures; empty: equal weights
};

// A microsatellite of the reference: `copies` units of `unit` bases from `begin`.
struct RepeatSite {
  std::uint32_t contig = 0;
  std::size_t begin = 0, unit = 1, copies = 0;
  std::size_t end() const { return begin + unit * copies; }
};

struct MsiSpec {
  std::vector<RepeatSite> sites;
  double rate = 0;            // the probability that a site slips in the trunk
  double contraction = 0.75;  // the probability that a slip loses a unit
};

struct TumourSpec {
  double purity = 0.7;
  std::vector<CloneSpec> clones{CloneSpec{}};
  std::vector<CopySegment> copy_number;          // the profile; bases outside every segment have 1 + 1
  bool whole_genome_doubling = false;            // recorded with the truth: the profile given includes it
  double early = 0.5;                            // the fraction of the trunk's mutations made before its copy-number changes
  std::vector<std::vector<double>> signatures;   // SBS-96 profiles; empty: one uniform profile
  double indel_mean = 1.5;
  std::size_t indel_max = 10;
  MsiSpec msi;
  std::size_t spacing = 5;                       // bases kept free around every variant
};

struct SomaticVariant {
  std::uint32_t contig = 0;
  std::size_t position = 0;  // 0-based, of REF's first base
  std::string ref, alt;      // letters, indels anchored on the base before
  int haplotype = 0;         // the parental haplotype it lies on
  int clone = 0;
  int multiplicity = 1;      // copies of its haplotype carrying it in a cell of its clone's descendants
  int signature = -1;        // the signature that made it (substitutions)
  int sbs96 = -1;            // its class (substitutions)
  bool early = false;
  bool microsatellite = false;  // a slip of a microsatellite (MSI); `site` its index in MsiSpec::sites
  std::size_t site = 0;
  double ccf = 1;            // its clone's
  double vaf = 0;            // the expected fraction of tumour-sample reads carrying it
  bool snv() const { return ref.size() == 1 && alt.size() == 1; }
  bool insertion() const { return alt.size() > ref.size(); }
  bool deletion() const { return ref.size() > alt.size(); }
};

struct Tumour {
  std::vector<Codes> reference;
  std::vector<Diploid> germline;           // the normal, one per contig
  double purity = 0;
  std::vector<CloneSpec> clones;
  std::vector<double> fraction;            // of tumour cells in each clone and none of its descendants
  std::vector<std::uint64_t> lineage;      // bit c of lineage[j]: clone c is j or an ancestor of j
  std::vector<CopySegment> copy_number;    // covering every contig, by contig and position
  bool whole_genome_doubling = false;
  std::vector<SomaticVariant> truth;       // by contig and position
  std::vector<std::vector<double>> signatures;
  std::vector<std::size_t> slipped;        // indices of the MSI sites that slipped
  std::size_t msi_sites = 0;               // the sites considered

  // The segment of the profile holding base p of a contig.
  const CopySegment& segment_at(std::uint32_t contig, std::size_t p) const {
    auto it = std::upper_bound(copy_number.begin(), copy_number.end(), std::pair{contig, p}, [](const std::pair<std::uint32_t, std::size_t>& x, const CopySegment& s) {
      return x.first != s.contig ? x.first < s.contig : x.second < s.begin;
    });
    if (it == copy_number.begin()) throw std::out_of_range("genomics::simulate::Tumour::segment_at: before the profile");
    return *std::prev(it);
  }
  // The mean tumour copy number over the reference (ASCAT's ψ when logR is normalised by length).
  double ploidy() const {
    double s = 0, n = 0;
    for (const auto& c : copy_number) {
      s += static_cast<double>(c.total()) * static_cast<double>(c.end - c.begin);
      n += static_cast<double>(c.end - c.begin);
    }
    return n > 0 ? s / n : 2;
  }
  // The expected fraction of reads carrying a variant of multiplicity m and CCF φ at base p.
  double expected_vaf(std::uint32_t contig, std::size_t p, double ccf, int m) const {
    const CopySegment& s = segment_at(contig, p);
    return purity * ccf * m / (2 * (1 - purity) + purity * s.total());
  }
};

namespace detail {

// The profile completed to cover every contig, segments sorted and checked.
inline std::vector<CopySegment> complete_profile(const std::vector<Codes>& reference, std::vector<CopySegment> given) {
  std::sort(given.begin(), given.end(), [](const CopySegment& a, const CopySegment& b) { return a.contig != b.contig ? a.contig < b.contig : a.begin < b.begin; });
  std::vector<CopySegment> out;
  std::size_t g = 0;
  for (std::uint32_t c = 0; c < reference.size(); c++) {
    std::size_t at = 0;
    for (; g < given.size() && given[g].contig == c; g++) {
      const CopySegment& s = given[g];
      if (s.begin < at || s.end > reference[c].size() || s.begin >= s.end || s.copies[0] < 0 || s.copies[1] < 0 || s.total() == 0)
        throw std::invalid_argument("genomics::simulate::make_tumour: copy-number segments inside their contig, disjoint, with copies, not both zero");
      if (s.begin > at) out.push_back({c, at, s.begin, {1, 1}});
      out.push_back(s);
      at = s.end;
    }
    if (at < reference[c].size()) out.push_back({c, at, reference[c].size(), {1, 1}});
  }
  if (g < given.size()) throw std::invalid_argument("genomics::simulate::make_tumour: a copy-number segment on a contig not given");
  return out;
}

// A position drawn uniformly over the contigs: (contig, position).
inline std::pair<std::uint32_t, std::size_t> uniform_position(calx::random::Stream& r, const std::vector<Codes>& reference, std::size_t total) {
  std::size_t u = r.below(total);
  for (std::uint32_t c = 0; c < reference.size(); c++) {
    if (u < reference[c].size()) return {c, u};
    u -= reference[c].size();
  }
  return {0, 0};
}

// Marks [a, b) of a bitmap, clamped.
inline void mark(std::vector<char>& taken, std::size_t a, std::size_t b) {
  for (std::size_t i = a; i < std::min(b, taken.size()); i++) taken[i] = 1;
}
inline bool free_span(const std::vector<char>& taken, const Codes& x, std::size_t a, std::size_t b) {
  if (b > x.size()) return false;
  for (std::size_t i = a; i < b; i++)
    if (taken[i] || x[i] > 3) return false;
  return true;
}

}  // namespace detail

// Mechanism-motivated SBS-96 profiles, constructed here (not COSMIC's, whose terms forbid their
// redistribution): each a background of weight 1 on every class with its process's motif raised.
//   0 "clock-like CpG (synthetic)": C>T at NpCpG ×40, the deamination of 5-methylcytosine;
//   1 "APOBEC-like (synthetic)": C>T and C>G at TpCpA and TpCpT ×30, APOBEC3's TCW motif;
//   2 "UV-like (synthetic)": C>T after a pyrimidine (CpC, TpC) ×25, and CC>TT's single-base part;
//   3 "tobacco-like (synthetic)": C>A everywhere ×6, ×12 at CpCpN, the guanine adducts of
//     benzo[a]pyrene read on the pyrimidine strand;
//   4 "flat (synthetic)": T>C ×2 and C>T ×1.5 over the background, a flat clock-like profile.
struct NamedSignature {
  std::string name;
  std::vector<double> profile;
};
inline std::vector<NamedSignature> synthetic_signatures() {
  using namespace sequence;
  std::vector<NamedSignature> out(5);
  const char* names[5] = {"clock-like CpG (synthetic)", "APOBEC-like (synthetic)", "UV-like (synthetic)", "tobacco-like (synthetic)", "flat (synthetic)"};
  for (int s = 0; s < 5; s++) {
    out[static_cast<std::size_t>(s)].name = names[s];
    auto& p = out[static_cast<std::size_t>(s)].profile;
    p.assign(SBS96, 1.0);
    for (std::size_t k = 0; k < SBS96; k++) {
      const SubstitutionContext c = sbs96_context(k);
      const bool ct = c.ref == C && c.alt == T, cg = c.ref == C && c.alt == G, ca = c.ref == C && c.alt == A, tc = c.ref == T && c.alt == C;
      if (s == 0 && ct && c.right == G) p[k] *= 40;
      if (s == 1 && (ct || cg) && c.left == T && (c.right == A || c.right == T)) p[k] *= 30;
      if (s == 2 && ct && (c.left == C || c.left == T)) p[k] *= 25;
      if (s == 3 && ca) p[k] *= c.left == C ? 12 : 6;
      if (s == 4 && tc) p[k] *= 2;
      if (s == 4 && ct) p[k] *= 1.5;
    }
    double t = 0;
    for (double x : p) t += x;
    for (double& x : p) x /= t;
  }
  return out;
}

struct CopyNumberSpec {
  std::size_t segments = 8;        // per contig
  std::size_t min_length = 20000;  // bases
  // The relative weights of the states (copies of the two haplotypes before doubling, the
  // haplotypes swapped with probability 1/2): 1+1, 2+1, 1+0, 2+0 (copy-neutral LOH), 3+1, 2+2.
  std::array<double, 6> weights{4, 2, 1.5, 1, 1, 0.5};
  bool whole_genome_doubling = false;
};

// A random allele-specific profile: each contig cut at `segments − 1` points drawn uniformly (at
// least min_length apart, else fewer), each piece a state drawn by the weights.
inline std::vector<CopySegment> random_copy_number(const std::vector<std::size_t>& lengths, const CopyNumberSpec& spec, std::uint64_t seed) {
  static constexpr std::array<std::array<int, 2>, 6> STATES{{{1, 1}, {2, 1}, {1, 0}, {2, 0}, {3, 1}, {2, 2}}};
  std::vector<CopySegment> out;
  for (std::uint32_t c = 0; c < lengths.size(); c++) {
    calx::random::Stream r(seed, COPY_PROFILE, c, 0);
    const std::size_t n = lengths[c];
    std::vector<std::size_t> cuts;
    for (std::size_t k = 1; k < spec.segments; k++) cuts.push_back(r.below(n));
    std::sort(cuts.begin(), cuts.end());
    std::vector<std::size_t> kept{0};
    for (std::size_t x : cuts)
      if (x >= kept.back() + spec.min_length && x + spec.min_length <= n) kept.push_back(x);
    kept.push_back(n);
    double total = 0;
    for (double w : spec.weights) total += w;
    for (std::size_t k = 0; k + 1 < kept.size(); k++) {
      double u = r.uniform() * total;
      std::size_t s = 0;
      while (s + 1 < STATES.size() && u >= spec.weights[s]) { u -= spec.weights[s]; s++; }
      std::array<int, 2> copies = STATES[s];
      if (r.bernoulli(0.5)) std::swap(copies[0], copies[1]);
      if (spec.whole_genome_doubling) { copies[0] *= 2; copies[1] *= 2; }
      out.push_back({c, kept[k], kept[k + 1], copies});
    }
  }
  return out;
}

// The tumour of a normal (`germline`, one diploid per contig of `reference`).
inline Tumour make_tumour(const std::vector<Codes>& reference, std::vector<Diploid> germline, const TumourSpec& spec, std::uint64_t seed) {
  if (germline.size() != reference.size()) throw std::invalid_argument("genomics::simulate::make_tumour: one diploid for each contig");
  if (!(spec.purity > 0 && spec.purity <= 1)) throw std::invalid_argument("genomics::simulate::make_tumour: a purity in (0, 1]");
  const std::size_t K = spec.clones.size();
  if (K == 0 || K > 64 || spec.clones[0].parent != -1 || spec.clones[0].ccf != 1)
    throw std::invalid_argument("genomics::simulate::make_tumour: 1 to 64 clones, the first the trunk (parent −1, CCF 1)");
  Tumour t;
  t.reference = reference;
  t.germline = std::move(germline);
  t.purity = spec.purity;
  t.clones = spec.clones;
  t.whole_genome_doubling = spec.whole_genome_doubling;
  t.fraction.assign(K, 0);
  t.lineage.assign(K, 0);
  std::vector<double> children(K, 0);
  for (std::size_t j = 0; j < K; j++) {
    const auto& c = spec.clones[j];
    if (j > 0 && (c.parent < 0 || static_cast<std::size_t>(c.parent) >= j || !(c.ccf > 0 && c.ccf <= 1)))
      throw std::invalid_argument("genomics::simulate::make_tumour: each clone's parent before it, its CCF in (0, 1]");
    t.lineage[j] = (j > 0 ? t.lineage[static_cast<std::size_t>(c.parent)] : 0) | (std::uint64_t{1} << j);
    if (j > 0) children[static_cast<std::size_t>(c.parent)] += c.ccf;
  }
  for (std::size_t j = 0; j < K; j++) {
    t.fraction[j] = spec.clones[j].ccf - children[j];
    if (t.fraction[j] < -1e-12) throw std::invalid_argument("genomics::simulate::make_tumour: the sum rule (the children's CCFs at most their parent's)");
    t.fraction[j] = std::max(0.0, t.fraction[j]);
  }
  t.copy_number = detail::complete_profile(reference, spec.copy_number);
  t.signatures = spec.signatures.empty() ? std::vector<std::vector<double>>{std::vector<double>(sequence::SBS96, 1.0 / sequence::SBS96)} : spec.signatures;
  for (const auto& s : t.signatures)
    if (s.size() != sequence::SBS96) throw std::invalid_argument("genomics::simulate::make_tumour: signatures of 96 classes");

  // What is taken: germline variants and their surroundings, then each somatic event's.
  std::vector<std::vector<char>> taken(reference.size());
  std::size_t total = 0;
  for (std::uint32_t c = 0; c < reference.size(); c++) {
    taken[c].assign(reference[c].size(), 0);
    total += reference[c].size();
    for (const auto& v : t.germline[c].truth) {
      const std::size_t a = v.position > spec.spacing ? v.position - spec.spacing : 0;
      detail::mark(taken[c], a, v.position + v.ref.size() + spec.spacing);
    }
    detail::mark(taken[c], 0, 1 + spec.spacing);
    const std::size_t n = reference[c].size();
    if (n > spec.spacing + 1) detail::mark(taken[c], n - spec.spacing - 1, n);
  }
  if (total == 0) throw std::invalid_argument("genomics::simulate::make_tumour: a reference with bases");
  calx::random::Stream r(seed, SOMATIC, 0, 0);

  // The haplotype and multiplicity of a mutation of clone j at base p: false when it is lost.
  const auto place = [&](SomaticVariant& v, std::size_t j) -> bool {
    const CopySegment& s = t.segment_at(v.contig, v.position);
    const int na = s.copies[0], nb = s.copies[1];
    v.clone = static_cast<int>(j);
    v.ccf = spec.clones[j].ccf;
    v.early = j == 0 && r.uniform() < spec.early;
    if (v.early) {
      v.haplotype = r.bernoulli(0.5) ? 1 : 0;
      v.multiplicity = s.copies[static_cast<std::size_t>(v.haplotype)];
      if (v.multiplicity == 0) return false;
    } else {
      v.haplotype = static_cast<int>(r.below(static_cast<std::uint64_t>(na + nb))) < na ? 0 : 1;
      v.multiplicity = 1;
    }
    v.vaf = t.expected_vaf(v.contig, v.position, v.ccf, v.multiplicity);
    return true;
  };

  // The microsatellites: all kept free; those that slip, slipped.
  t.msi_sites = 0;
  for (std::size_t i = 0; i < spec.msi.sites.size(); i++) {
    const RepeatSite& site = spec.msi.sites[i];
    if (site.contig >= reference.size() || site.begin < 1 + spec.spacing || site.end() + spec.spacing + 1 > reference[site.contig].size() || site.unit == 0) continue;
    const Codes& x = reference[site.contig];
    const std::size_t a = site.begin - 1 - spec.spacing, b = site.end() + spec.spacing;
    if (!detail::free_span(taken[site.contig], x, a, b)) continue;
    detail::mark(taken[site.contig], a, b);
    t.msi_sites++;
    if (!(r.uniform() < spec.msi.rate)) continue;
    SomaticVariant v;
    v.contig = site.contig;
    v.position = site.begin - 1;
    v.microsatellite = true;
    v.site = i;
    std::string unit;
    for (std::size_t k = 0; k < site.unit; k++) unit += sequence::dna_letter(x[site.begin + k]);
    const std::string anchor(1, sequence::dna_letter(x[site.begin - 1]));
    if (r.uniform() < spec.msi.contraction) {
      v.ref = anchor + unit;
      v.alt = anchor;
    } else {
      v.ref = anchor;
      v.alt = anchor + unit;
    }
    if (!place(v, 0)) continue;
    t.slipped.push_back(i);
    t.truth.push_back(std::move(v));
  }

  // Substitutions and indels, clone by clone.
  std::vector<double> cumulative;
  for (std::size_t j = 0; j < K; j++) {
    const CloneSpec& clone = spec.clones[j];
    std::vector<double> e = clone.exposures;
    if (e.empty()) e.assign(t.signatures.size(), 1.0);
    if (e.size() != t.signatures.size()) throw std::invalid_argument("genomics::simulate::make_tumour: one exposure for each signature");
    double et = 0;
    for (double x : e) et += std::max(0.0, x);
    if (!(et > 0)) throw std::invalid_argument("genomics::simulate::make_tumour: positive exposures");
    for (std::size_t m = 0; m < clone.snvs; m++) {
      double u = r.uniform() * et;
      std::size_t s = 0;
      while (s + 1 < e.size() && u >= std::max(0.0, e[s])) { u -= std::max(0.0, e[s]); s++; }
      const auto& prof = t.signatures[s];
      double pt = 0;
      for (double x : prof) pt += x;
      double w = r.uniform() * pt;
      std::size_t k = 0;
      while (k + 1 < prof.size() && w >= prof[k]) { w -= prof[k]; k++; }
      const sequence::SubstitutionContext want = sequence::sbs96_context(k);
      // A position of the class's trinucleotide on either strand, by rejection.
      for (int attempt = 0; attempt < 20000; attempt++) {
        const auto [c, p] = detail::uniform_position(r, reference, total);
        const Codes& x = reference[c];
        if (p == 0 || p + 1 >= x.size() || !detail::free_span(taken[c], x, p - 1, p + 2)) continue;
        sequence::Code alt;
        if (x[p - 1] == want.left && x[p] == want.ref && x[p + 1] == want.right) alt = want.alt;
        else if (x[p - 1] == sequence::complement(want.right) && x[p] == sequence::complement(want.ref) && x[p + 1] == sequence::complement(want.left)) alt = sequence::complement(want.alt);
        else continue;
        SomaticVariant v;
        v.contig = c;
        v.position = p;
        v.ref = std::string(1, sequence::dna_letter(x[p]));
        v.alt = std::string(1, sequence::dna_letter(alt));
        v.signature = static_cast<int>(s);
        v.sbs96 = static_cast<int>(k);
        if (!place(v, j)) break;
        detail::mark(taken[c], p > spec.spacing ? p - spec.spacing : 0, p + 1 + spec.spacing);
        t.truth.push_back(std::move(v));
        break;
      }
    }
    for (std::size_t m = 0; m < clone.indels; m++) {
      const double q = 1 / std::max(1.0, spec.indel_mean);
      const std::size_t len = std::min<std::size_t>(spec.indel_max, 1 + r.geometric(q));
      const bool deletion = r.bernoulli(0.5);
      for (int attempt = 0; attempt < 1000; attempt++) {
        const auto [c, p] = detail::uniform_position(r, reference, total);
        const Codes& x = reference[c];
        const std::size_t span = deletion ? len + 1 : 1;
        if (!detail::free_span(taken[c], x, p, p + span)) continue;
        SomaticVariant v;
        v.contig = c;
        v.position = p;
        for (std::size_t i = p; i < p + span; i++) v.ref += sequence::dna_letter(x[i]);
        if (deletion) {
          v.alt = v.ref.substr(0, 1);
        } else {
          v.alt = v.ref;
          for (std::size_t i = 0; i < len; i++) v.alt += sequence::DNA_LETTERS[r.below(4)];
        }
        if (!place(v, j)) break;
        detail::mark(taken[c], p > spec.spacing ? p - spec.spacing : 0, p + span + spec.spacing);
        t.truth.push_back(std::move(v));
        break;
      }
    }
  }
  std::sort(t.truth.begin(), t.truth.end(), [](const SomaticVariant& a, const SomaticVariant& b) { return a.contig != b.contig ? a.contig < b.contig : a.position < b.position; });
  return t;
}

// The normal sample's source: the germline diploid.
inline Source normal_source(const Tumour& t) { return Source::diploid(t.germline); }

// The catalogues of a cohort drawn from the model of signature analysis: tumour i has
// N_ij ~ Poisson(exposures[i][j]) mutations of signature j, distributed over the classes as
// Multinomial(N_ij, signature j) (calx's exact samplers), so that its catalogue has Poisson counts of
// means Σ_j exposures[i][j] · signature j (the Poisson model behind NMF's Kullback–Leibler objective).
// Tumour i draws from a stream keyed by the seed, CATALOGUES and i.
inline constexpr std::uint32_t CATALOGUES = 0x43415441;  // "CATA"
inline std::vector<std::vector<double>> simulate_catalogues(const std::vector<std::vector<double>>& signatures, const std::vector<std::vector<double>>& exposures, std::uint64_t seed) {
  std::vector<std::vector<double>> out;
  const std::size_t m = signatures.empty() ? 0 : signatures[0].size();
  for (std::size_t i = 0; i < exposures.size(); i++) {
    if (exposures[i].size() != signatures.size()) throw std::invalid_argument("genomics::simulate::simulate_catalogues: one exposure per signature");
    calx::random::Stream r(seed, CATALOGUES, static_cast<std::uint32_t>(i), 0);
    std::vector<double> c(m, 0);
    for (std::size_t j = 0; j < signatures.size(); j++) {
      const auto n = r.poisson(std::max(0.0, exposures[i][j]));
      if (n == 0) continue;
      const auto k = calx::random::multinomial(r, n, signatures[j]);
      for (std::size_t q = 0; q < m; q++) c[q] += static_cast<double>(k[q]);
    }
    out.push_back(std::move(c));
  }
  return out;
}

// The number of reads giving a mean depth over the reference (normal: 2 copies everywhere; the
// tumour's depth then varies with copy number about this mean).
inline std::size_t reads_for_depth(const Tumour& t, double depth, std::size_t read_length) {
  std::size_t n = 0;
  for (const auto& c : t.reference) n += c.size();
  return static_cast<std::size_t>(std::llround(depth * static_cast<double>(n) / static_cast<double>(read_length)));
}

// Reads first, ..., first + count − 1 of the tumour sample.
inline std::vector<Read> tumour_reads(const Tumour& t, std::size_t count, const ShortReadModel& model, std::uint64_t seed, std::size_t first = 0) {
  std::vector<Read> out;
  out.reserve(count);
  const double rho = t.purity;
  double tmax = 0;
  for (const auto& s : t.copy_number) tmax = std::max(tmax, 2 * (1 - rho) + rho * s.total());
  std::size_t total = 0;
  for (const auto& c : t.reference) total += c.size();
  const std::size_t F = model.length + 8, W = F + 64;
  std::vector<double> population(t.fraction);  // cumulative weights of the clones
  for (std::size_t j = 1; j < population.size(); j++) population[j] += population[j - 1];
  Codes bases;
  std::vector<std::size_t> ref_of;
  std::vector<char> inserted;
  for (std::size_t i = first; i < first + count; i++) {
    calx::random::Stream r(seed, TUMOUR_READS, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i >> 32));
    Read read;
    read.name = "tumour" + std::to_string(i);
    // A start, kept in proportion to the DNA there.
    std::uint32_t c = 0;
    std::size_t x = 0;
    const CopySegment* seg = nullptr;
    bool placed = false;
    for (int attempt = 0; attempt < 1000 && !placed; attempt++) {
      std::tie(c, x) = detail::uniform_position(r, t.reference, total);
      if (x + F > t.reference[c].size()) continue;
      seg = &t.segment_at(c, x);
      placed = r.uniform() * tmax < 2 * (1 - rho) + rho * seg->total();
    }
    if (!placed) continue;
    // The molecule.
    const double v = r.uniform() * (2 * (1 - rho) + rho * seg->total());
    int h = 0;
    long long clone = -1, copy = 0;
    if (v < 2 * (1 - rho)) {
      h = v < 1 - rho ? 0 : 1;
    } else {
      const double w = v - 2 * (1 - rho);
      h = w < rho * seg->copies[0] ? 0 : 1;
      copy = static_cast<long long>(r.below(static_cast<std::uint64_t>(seg->copies[static_cast<std::size_t>(h)])));
      const double u = r.uniform() * population.back();
      clone = 0;
      while (static_cast<std::size_t>(clone) + 1 < population.size() && u >= population[static_cast<std::size_t>(clone)]) clone++;
    }
    read.contig = c;
    read.haplotype = h;
    // The window of the germline haplotype from the start, with what each base is on the reference.
    const Diploid& g = t.germline[c];
    const Codes& hap = g.haplotype[static_cast<std::size_t>(h)];
    const Liftover& lift = g.liftover[static_cast<std::size_t>(h)];
    const std::size_t hs = detail::to_haplotype(lift, x);
    const std::size_t he = std::min(hap.size(), hs + W);
    bases.assign(hap.begin() + static_cast<std::ptrdiff_t>(hs), hap.begin() + static_cast<std::ptrdiff_t>(he));
    ref_of.resize(bases.size());
    inserted.assign(bases.size(), 0);
    for (std::size_t k = 0; k < bases.size(); k++) {
      const Locus l = lift.locate(hs + k);
      ref_of[k] = l.position;
      inserted[k] = l.inserted;
    }
    // The somatic variants this molecule carries.
    if (clone >= 0 && !bases.empty()) {
      const std::size_t lo = ref_of.front() > 32 ? ref_of.front() - 32 : 0, hi = ref_of.back();
      auto it = std::lower_bound(t.truth.begin(), t.truth.end(), std::pair{c, lo}, [](const SomaticVariant& a, const std::pair<std::uint32_t, std::size_t>& k) {
        return a.contig != k.first ? a.contig < k.first : a.position < k.second;
      });
      // Right to left, so that the indices of the bases before each variant stay put.
      std::vector<const SomaticVariant*> carried;
      for (; it != t.truth.end() && it->contig == c && it->position <= hi; ++it)
        if (it->haplotype == h && ((t.lineage[static_cast<std::size_t>(clone)] >> it->clone) & 1) && copy < it->multiplicity) carried.push_back(&*it);
      for (auto p = carried.rbegin(); p != carried.rend(); ++p) {
        const SomaticVariant& s = **p;
        std::size_t k = 0;
        while (k < bases.size() && (inserted[k] || ref_of[k] != s.position)) k++;
        if (s.deletion()) {
          // Remove the deleted reference bases [position + 1, position + |REF|) found in the window.
          std::size_t w = 0;
          for (std::size_t q = 0; q < bases.size(); q++) {
            const bool gone = !inserted[q] && ref_of[q] > s.position && ref_of[q] < s.position + s.ref.size();
            if (gone) continue;
            bases[w] = bases[q];
            ref_of[w] = ref_of[q];
            inserted[w] = inserted[q];
            w++;
          }
          bases.resize(w);
          ref_of.resize(w);
          inserted.resize(w);
        } else if (k < bases.size()) {
          if (s.snv()) {
            bases[k] = sequence::dna_code(s.alt[0]);
          } else {
            Codes add;
            for (std::size_t q = 1; q < s.alt.size(); q++) add.push_back(sequence::dna_code(s.alt[q]));
            bases.insert(bases.begin() + static_cast<std::ptrdiff_t>(k + 1), add.begin(), add.end());
            ref_of.insert(ref_of.begin() + static_cast<std::ptrdiff_t>(k + 1), add.size(), s.position + 1);
            inserted.insert(inserted.begin() + static_cast<std::ptrdiff_t>(k + 1), add.size(), 1);
          }
        }
      }
    }
    if (bases.size() < F) continue;
    read.reverse = r.bernoulli(0.5);
    Codes tpl(bases.begin(), bases.begin() + static_cast<std::ptrdiff_t>(F));
    if (read.reverse) tpl = sequence::reverse_complement(tpl);
    const std::size_t used = detail::sequence_short(r, tpl, model, read.bases, read.quality, read.errors);
    if (used == 0) continue;
    const std::size_t a = read.reverse ? F - used : 0, b = read.reverse ? F : used;
    read.begin = ref_of[a];
    read.end = std::max(read.begin, ref_of[b - 1] + 1);
    out.push_back(std::move(read));
  }
  return out;
}

}  // namespace genomics::simulate
