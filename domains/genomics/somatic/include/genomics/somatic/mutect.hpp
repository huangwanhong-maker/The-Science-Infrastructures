// Somatic substitutions and indels called from a tumour and its matched normal, by MuTect's two
// likelihood-ratio tests and its filters (Cibulskis et al. 2013, "Sensitive detection of somatic
// point mutations in impure and heterogeneous cancer samples", Nature Biotechnology 31, 213–219,
// Online Methods and Table 1).
//
// - The tumour test. At a site with reference base r, the candidate m is the most frequent other
//   base of the tumour's pileup. Under M_f^m a fraction f of the molecules carry m; a base b read
//   with error probability e (its Phred quality) has
//     P(b | e, r, m, f) = f P(b | e, m) + (1 − f) P(b | e, r),   P(b | e, A) = 1 − e if b = A, e/3 otherwise,
//   the reads independent; M₀ is f = 0. LOD_T = log₁₀ L(M_{f̂}^m) / L(M₀) at f̂ = (reads of m)/depth
//   (MuTect does not maximise over f). The site is a candidate when LOD_T ≥ θ_T = log₁₀ δ_T −
//   log₁₀ P(m, f) with δ_T = 2 and P(m, f) = P(m) P(f) = (3 · 10⁻⁶ / 3) · 1: θ_T = 6.3.
//   *Remark.* At most one read of m can never pass: its LOD is below log₁₀(3/e) ≤ 4.6 for e ≥
//   10^{−4.1} (Phred 41), so sites with fewer than two reads of m are skipped unweighed.
// - The normal test. LOD_N = log₁₀ L(M₀) / L(M_{0.5}^m) on the normal's pileup: the evidence that the
//   normal is not heterozygous for m. Somatic when LOD_N ≥ θ_N = log₁₀ δ_N + log₁₀ [P(germline) /
//   P(somatic)], δ_N = 10, P(somatic) = 3 · 10⁻⁶: θ_N = 2.2 with P(germline) = 5 · 10⁻⁵ (a site not
//   in dbSNP), 5.5 with 0.095 (a dbSNP site; `known_germline` lists such sites).
// - Filters (Table 1), as VCF FILTER names:
//     proximal_gap      ≥ 3 tumour reads with insertions, or ≥ 3 with deletions, anchored in the
//                       11 bases centred on the site;
//     poor_mapping      ≥ 50% of the reads over the site in tumour and normal together of mapping
//                       quality 0, or no read of m of mapping quality ≥ 20;
//     triallelic_site   the normal heterozygous for another allele B ≠ m (≥ 2 reads and ≥ 20% of
//                       its depth);
//     strand_bias       the LOD of the reads of one strand below 2.0 where the power to reach 2.0
//                       on that strand is ≥ 0.9;
//     clustered_position the reads of m at a consistent distance from their alignment's start or
//                       end: median ≤ 10 and median absolute deviation ≤ 3, at either end;
//     observed_in_control the normal with ≥ 2 reads of m or ≥ 3% of its depth, and their
//                       qualities summing above 20.
//   The panel of normals (sites called in two or more normals) needs other normals; not here.
// - Power (the paper's sensitivity calculation): for n reads of error e and allele fraction f, the
//   least k for which k reads of m and n − k of r give LOD ≥ θ, and P(X ≥ k) for X ~ Binomial(n,
//   p_m), p_m = f(1 − e) + (1 − f) e/3 the model's probability of reading m.
// - Indels (not MuTect's, which calls substitutions only): the same two tests on counts, every
//   read spanning the indel's interval of ambiguity (the caller's rule, variant/caller.hpp) showing
//   it with probability f(1 − ε) + (1 − f) ε, ε = `indel_error`; the strand, control and gap
//   filters as for substitutions, the gap filter counting other reads' indels only.
// - The callable territory: the bases where the tumour has at least `callable_tumour` reads and
//   the normal `callable_normal`: 14, the least tumour depth at which a mutation at f = 0.3 is
//   called with power ≥ 0.8 at Q30, and 8, the least normal depth at which reference reads alone
//   (Q30) reach LOD_N ≥ 2.2 (both computed in the checks); the denominator of the mutational burden.
//
// The genome is called in windows on calx::scheduler::parallel_for, each window's calls in its
// own slot, joined in order: the same calls on any number of threads.
//
// Complexity: O(Σ ℓ) for the two pileups, O(d) per site with two reads of m, plus O(d ℓ) at the
// candidates for the reads' positions. Limits: no local realignment or assembly (Mutect2 and
// Strelka2 realign reads to candidate haplotypes; here a site's bases are weighed as aligned), one
// alternate allele per site, no estimate of contamination (MuTect's ContEst), and no panel of
// normals.

#pragma once

#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/somatic/coverage.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/normalize.hpp>
#include <genomics/variant/pileup.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace genomics::somatic {

using sequence::Code;
using sequence::Codes;

struct MutectOptions {
  variant::PileupOptions pileup;
  double theta_tumour = 6.3, theta_normal = 2.2, theta_normal_known = 5.5;
  std::size_t min_alt = 2;
  // Filters.
  std::size_t gap_reads = 3, gap_half_window = 5;
  double poor_mapping_fraction = 0.5;
  int alt_mapq = 20;
  std::size_t triallelic_reads = 2;
  double triallelic_fraction = 0.2;
  double strand_lod = 2.0, strand_power = 0.9;
  double cluster_median = 10, cluster_mad = 3;
  std::size_t control_reads = 2;
  double control_fraction = 0.03, control_quality = 20;
  // Indels.
  bool indels = true;
  double indel_error = 0.005;
  // Territory.
  std::size_t callable_tumour = 14, callable_normal = 8;
  std::size_t window = 1 << 15;
  // Sites known to vary in the germline (dbSNP's), by contig and position: θ_N = 5.5 there.
  std::vector<std::pair<std::uint32_t, std::size_t>> known_germline;
};

struct SomaticCall {
  std::uint32_t contig = 0;
  std::size_t position = 0;     // 0-based, of REF's first base
  std::string ref, alt;
  bool indel = false;
  double lod_tumour = 0, lod_normal = 0;
  std::size_t tumour_depth = 0, tumour_alt = 0, normal_depth = 0, normal_alt = 0;
  std::size_t alt_forward = 0, alt_reverse = 0;
  double lod_forward = 0, lod_reverse = 0, power_forward = 0, power_reverse = 0;
  double median_start = 0, mad_start = 0, median_end = 0, mad_end = 0;  // the reads of m: distance from their alignment's ends
  std::vector<std::string> filters;  // empty: PASS
  bool passed() const { return filters.empty(); }
  double vaf() const { return tumour_depth ? static_cast<double>(tumour_alt) / static_cast<double>(tumour_depth) : 0; }
};

struct SomaticCalls {
  std::vector<SomaticCall> calls;   // by contig and position
  std::size_t callable = 0;         // bases of the callable territory
  std::size_t sites_weighed = 0;    // sites whose LOD_T was computed
};

// ---- the likelihoods ---------------------------------------------------------------------------------

// log₁₀ L(M_f^m) of bases (base, Phred quality) at a site with reference r.
inline double log10_likelihood(std::span<const variant::BaseCall> bases, Code r, Code m, double f) {
  double l = 0;
  for (const auto& b : bases) {
    const double e = std::pow(10.0, -static_cast<double>(b.quality) / 10.0);
    double p;
    if (b.base == m) p = f * (1 - e) + (1 - f) * e / 3;
    else if (b.base == r) p = f * e / 3 + (1 - f) * (1 - e);
    else p = e / 3;
    l += std::log10(p);
  }
  return l;
}

// LOD_T = log₁₀ L(M_{f̂}^m)/L(M₀), f̂ the fraction of the bases that are m.
inline double tumour_lod(std::span<const variant::BaseCall> bases, Code r, Code m) {
  if (bases.empty()) return 0;
  std::size_t k = 0;
  for (const auto& b : bases) k += b.base == m;
  const double f = static_cast<double>(k) / static_cast<double>(bases.size());
  return log10_likelihood(bases, r, m, f) - log10_likelihood(bases, r, m, 0);
}

// LOD_N = log₁₀ L(M₀)/L(M_{0.5}^m).
inline double normal_lod(std::span<const variant::BaseCall> bases, Code r, Code m) {
  return log10_likelihood(bases, r, m, 0) - log10_likelihood(bases, r, m, 0.5);
}

// The counting form (indels): k of d reads show the event, each misreported with probability ε.
inline double log10_likelihood_counts(std::size_t k, std::size_t d, double error, double f) {
  const double show = f * (1 - error) + (1 - f) * error;
  return static_cast<double>(k) * std::log10(show) + static_cast<double>(d - k) * std::log10(1 - show);
}
inline double tumour_lod_counts(std::size_t k, std::size_t d, double error) {
  if (d == 0) return 0;
  const double f = static_cast<double>(k) / static_cast<double>(d);
  return log10_likelihood_counts(k, d, error, f) - log10_likelihood_counts(k, d, error, 0);
}
inline double normal_lod_counts(std::size_t k, std::size_t d, double error) {
  return log10_likelihood_counts(k, d, error, 0) - log10_likelihood_counts(k, d, error, 0.5);
}

// The LOD of k reads of m and n − k of r at one error e (the power's model).
inline double lod_at(std::size_t k, std::size_t n, double e) {
  if (n == 0) return 0;
  const double f = static_cast<double>(k) / static_cast<double>(n);
  const double pm1 = f * (1 - e) + (1 - f) * e / 3, pr1 = f * e / 3 + (1 - f) * (1 - e);
  return static_cast<double>(k) * (std::log10(pm1) - std::log10(e / 3)) + static_cast<double>(n - k) * (std::log10(pr1) - std::log10(1 - e));
}

// MuTect's sensitivity: P(LOD ≥ θ) for n reads of error e from molecules a fraction f of which carry m.
inline double power(std::size_t n, double f, double e, double theta) {
  if (n == 0) return 0;
  std::size_t k = 1;
  while (k <= n && lod_at(k, n, e) < theta) k++;
  if (k > n) return 0;
  const double pm = f * (1 - e) + (1 - f) * e / 3;
  // P(X ≥ k) = 1 − F(k − 1).
  return k == 0 ? 1 : 1 - calx::special_functions::binomial_cdf(static_cast<double>(k - 1), static_cast<double>(n), pm);
}

// The median and the median absolute deviation (unscaled) of a sample.
inline std::pair<double, double> median_mad(std::vector<double> x) {
  if (x.empty()) return {0, 0};
  const auto med = [](std::vector<double>& v) {
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
  };
  const double m = med(x);
  for (double& v : x) v = std::fabs(v - m);
  return {m, med(x)};
}

// The read's base aligned to reference base p (on the strand mapped), or −1.
inline long long query_offset(const mapping::Mapping& m, std::size_t p) {
  std::size_t i = m.alignment.query_begin, j = m.alignment.target_begin;
  for (const auto& op : m.alignment.cigar) {
    const bool q = align::consumes_query(op.op), t = align::consumes_target(op.op);
    if (q && t) {
      if (p >= j && p < j + op.length) return static_cast<long long>(i + (p - j));
    } else if (t && p >= j && p < j + op.length) {
      return -1;
    }
    if (q) i += op.length;
    if (t) j += op.length;
  }
  return -1;
}

namespace detail {

struct Window {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
};

struct SiteMappings {
  std::size_t total = 0, mapq0 = 0;
};
inline SiteMappings mapping_qualities(std::span<const mapping::Mapping> ms, std::size_t p) {
  SiteMappings s;
  for (const auto& m : ms)
    if (m.target_begin <= p && p < m.target_end) {
      s.total++;
      s.mapq0 += m.mapq == 0;
    }
  return s;
}

// The distances of the reads of m from their alignment's start and end at p.
inline void positions_in_reads(std::span<const mapping::Mapping> ms, std::span<const std::uint32_t> reads, std::size_t p, SomaticCall& c) {
  std::vector<double> from_start, from_end;
  for (const auto& m : ms) {
    if (!(m.target_begin <= p && p < m.target_end)) continue;
    if (!std::binary_search(reads.begin(), reads.end(), m.read)) continue;
    const long long i = query_offset(m, p);
    if (i < 0) continue;
    from_start.push_back(static_cast<double>(static_cast<std::size_t>(i) - m.alignment.query_begin));
    from_end.push_back(static_cast<double>(m.alignment.query_end - 1 - static_cast<std::size_t>(i)));
  }
  std::tie(c.median_start, c.mad_start) = median_mad(from_start);
  std::tie(c.median_end, c.mad_end) = median_mad(from_end);
}

inline double mean_error(std::span<const variant::BaseCall> b) {
  if (b.empty()) return 1e-3;
  double s = 0;
  for (const auto& x : b) s += std::pow(10.0, -static_cast<double>(x.quality) / 10.0);
  return s / static_cast<double>(b.size());
}

// Reads with insertions and with deletions anchored in [p − w, p + w], other than those in `skip`.
inline std::array<std::size_t, 2> gaps_near(const variant::Pileup& pu, std::size_t p, std::size_t w, std::span<const std::uint32_t> skip = {}) {
  std::vector<std::uint32_t> ins, del;
  const std::size_t lo = p > w ? p - w : 0;
  auto it = std::lower_bound(pu.indels.begin(), pu.indels.end(), lo, [](const variant::IndelObservation& e, std::size_t x) { return e.position < x; });
  for (; it != pu.indels.end() && it->position <= p + w; ++it) {
    if (std::binary_search(skip.begin(), skip.end(), it->read)) continue;
    (it->length > 0 ? ins : del).push_back(it->read);
  }
  const auto distinct = [](std::vector<std::uint32_t>& v) {
    std::sort(v.begin(), v.end());
    return static_cast<std::size_t>(std::unique(v.begin(), v.end()) - v.begin());
  };
  return {distinct(ins), distinct(del)};
}

struct Context {
  const Codes& reference;
  const variant::Pileup& tumour;
  const variant::Pileup& normal;
  std::span<const mapping::Mapping> tumour_maps, normal_maps;
  const MutectOptions& o;
  bool known(std::uint32_t contig, std::size_t p) const { return std::binary_search(o.known_germline.begin(), o.known_germline.end(), std::pair{contig, p}); }
};

inline void substitution_at(const Context& x, std::size_t p, std::vector<variant::BaseCall>& tb, std::vector<variant::BaseCall>& nb, SomaticCalls& out) {
  const Code r = x.reference[p];
  if (r > 3) return;
  const auto tobs = x.tumour.at(p);
  std::array<std::size_t, 5> count{};
  for (const auto& b : tobs) count[std::min<std::size_t>(b.base, 4)]++;
  Code m = 4;
  for (Code b = 0; b < 4; b++)
    if (b != r && (m == 4 || count[b] > count[m])) m = b;
  if (m == 4 || count[m] < x.o.min_alt) return;
  tb.clear();
  for (const auto& b : tobs) tb.push_back({b.base, b.quality});
  out.sites_weighed++;
  const double lod_t = tumour_lod(tb, r, m);
  if (lod_t < x.o.theta_tumour) return;
  const auto nobs = x.normal.at(p);
  nb.clear();
  for (const auto& b : nobs) nb.push_back({b.base, b.quality});
  const double lod_n = normal_lod(nb, r, m);
  const double theta_n = x.known(x.tumour.contig, p) ? x.o.theta_normal_known : x.o.theta_normal;
  if (lod_n < theta_n) return;
  SomaticCall c;
  c.contig = x.tumour.contig;
  c.position = p;
  c.ref = std::string(1, sequence::dna_letter(r));
  c.alt = std::string(1, sequence::dna_letter(m));
  c.lod_tumour = lod_t;
  c.lod_normal = lod_n;
  c.tumour_depth = tb.size();
  c.tumour_alt = count[m];
  c.normal_depth = nb.size();
  std::array<std::size_t, 5> ncount{};
  double nquality = 0;
  for (const auto& b : nb) {
    ncount[std::min<std::size_t>(b.base, 4)]++;
    if (b.base == m) nquality += b.quality;
  }
  c.normal_alt = ncount[m];
  // Strands.
  std::vector<variant::BaseCall> fw, rv;
  std::vector<std::uint32_t> alt_reads;
  int best_alt_mapq = 0;
  for (const auto& b : tobs) {
    (b.reverse ? rv : fw).push_back({b.base, b.quality});
    if (b.base == m) {
      (b.reverse ? c.alt_reverse : c.alt_forward)++;
      alt_reads.push_back(b.read);
      best_alt_mapq = std::max<int>(best_alt_mapq, b.mapq);
    }
  }
  std::sort(alt_reads.begin(), alt_reads.end());
  const double f = c.vaf();
  c.lod_forward = tumour_lod(fw, r, m);
  c.lod_reverse = tumour_lod(rv, r, m);
  c.power_forward = power(fw.size(), f, mean_error(fw), x.o.strand_lod);
  c.power_reverse = power(rv.size(), f, mean_error(rv), x.o.strand_lod);
  positions_in_reads(x.tumour_maps, alt_reads, p, c);
  // Filters.
  const auto gaps = gaps_near(x.tumour, p, x.o.gap_half_window);
  if (gaps[0] >= x.o.gap_reads || gaps[1] >= x.o.gap_reads) c.filters.push_back("proximal_gap");
  const auto mt = mapping_qualities(x.tumour_maps, p), mn = mapping_qualities(x.normal_maps, p);
  if ((mt.total + mn.total > 0 && static_cast<double>(mt.mapq0 + mn.mapq0) >= x.o.poor_mapping_fraction * static_cast<double>(mt.total + mn.total)) || best_alt_mapq < x.o.alt_mapq)
    c.filters.push_back("poor_mapping");
  for (Code b = 0; b < 4; b++)
    if (b != r && b != m && ncount[b] >= x.o.triallelic_reads && static_cast<double>(ncount[b]) >= x.o.triallelic_fraction * static_cast<double>(nb.size())) {
      c.filters.push_back("triallelic_site");
      break;
    }
  if ((c.power_forward >= x.o.strand_power && c.lod_forward < x.o.strand_lod) || (c.power_reverse >= x.o.strand_power && c.lod_reverse < x.o.strand_lod))
    c.filters.push_back("strand_bias");
  if ((c.median_start <= x.o.cluster_median && c.mad_start <= x.o.cluster_mad) || (c.median_end <= x.o.cluster_median && c.mad_end <= x.o.cluster_mad))
    c.filters.push_back("clustered_position");
  if ((c.normal_alt >= x.o.control_reads || static_cast<double>(c.normal_alt) >= x.o.control_fraction * static_cast<double>(nb.size())) && nquality > x.o.control_quality)
    c.filters.push_back("observed_in_control");
  out.calls.push_back(std::move(c));
}

// The interval a read must span to tell an indel at anchor p from the reference (caller.hpp's rule).
inline std::size_t informative_end(const Codes& reference, const variant::IndelObservation& e) {
  const std::size_t p = e.position, size = reference.size();
  if (e.length < 0) {
    const auto L = static_cast<std::size_t>(-e.length);
    std::size_t a = p + 1;
    while (a + L < size && reference[a] == reference[a + L]) a++;
    return a + L + 1;
  }
  std::size_t shift = 0;
  while (p + 1 + shift < size && sequence::dna_code(e.bases[shift % e.bases.size()]) == reference[p + 1 + shift]) shift++;
  return p + shift + 2;
}

// Carriers of an allele in a pileup's indels at one anchor, and the reads spanning its interval.
struct IndelSupport {
  std::size_t carriers = 0, carriers_forward = 0, depth = 0, forward = 0;
  std::vector<std::uint32_t> reads;  // the carriers, sorted
};
inline IndelSupport indel_support(const variant::Pileup& pu, std::size_t p, std::int32_t length, const std::string& bases, std::size_t informative) {
  IndelSupport s;
  auto it = std::lower_bound(pu.indels.begin(), pu.indels.end(), p, [](const variant::IndelObservation& e, std::size_t x) { return e.position < x; });
  for (; it != pu.indels.end() && it->position == p; ++it)
    if (it->length == length && it->bases == bases) {
      s.carriers++;
      s.carriers_forward += !it->reverse;
      s.reads.push_back(it->read);
    }
  std::sort(s.reads.begin(), s.reads.end());
  std::size_t ref_f = 0, ref_r = 0;
  pu.for_each_spanning(p, informative, [&](const variant::ReadSpan& sp) {
    if (std::binary_search(s.reads.begin(), s.reads.end(), sp.read)) return;
    (sp.reverse ? ref_r : ref_f)++;
  });
  s.depth = s.carriers + ref_f + ref_r;
  s.forward = s.carriers_forward + ref_f;
  return s;
}

inline void indel_at(const Context& x, std::size_t first, std::size_t last, SomaticCalls& out) {
  const auto& ind = x.tumour.indels;
  std::size_t best = first, best_n = 0;
  for (std::size_t i = first; i < last;) {
    std::size_t j = i;
    while (j < last && ind[j].length == ind[i].length && ind[j].bases == ind[i].bases) j++;
    if (j - i > best_n) { best = i; best_n = j - i; }
    i = j;
  }
  if (best_n < x.o.min_alt) return;
  const variant::IndelObservation e = ind[best];
  const std::size_t informative = informative_end(x.reference, e);
  if (informative > x.reference.size()) return;
  const IndelSupport t = indel_support(x.tumour, e.position, e.length, e.bases, informative);
  out.sites_weighed++;
  const double lod_t = tumour_lod_counts(t.carriers, t.depth, x.o.indel_error);
  if (lod_t < x.o.theta_tumour) return;
  const IndelSupport n = indel_support(x.normal, e.position, e.length, e.bases, informative);
  const double lod_n = normal_lod_counts(n.carriers, n.depth, x.o.indel_error);
  if (lod_n < x.o.theta_normal) return;
  std::string ref(1, sequence::dna_letter(x.reference[e.position])), alt = ref;
  if (e.length > 0) alt += e.bases;
  else ref += e.bases;
  variant::Normalized nz = variant::normalize(e.position, ref, alt, x.reference);
  SomaticCall c;
  c.contig = x.tumour.contig;
  c.position = nz.position;
  c.ref = std::move(nz.ref);
  c.alt = std::move(nz.alts[0]);
  c.indel = true;
  c.lod_tumour = lod_t;
  c.lod_normal = lod_n;
  c.tumour_depth = t.depth;
  c.tumour_alt = t.carriers;
  c.normal_depth = n.depth;
  c.normal_alt = n.carriers;
  c.alt_forward = t.carriers_forward;
  c.alt_reverse = t.carriers - t.carriers_forward;
  const std::size_t fw = t.forward, rv = t.depth - t.forward;
  const double f = c.vaf(), eps = x.o.indel_error;
  c.lod_forward = tumour_lod_counts(c.alt_forward, fw, eps);
  c.lod_reverse = tumour_lod_counts(c.alt_reverse, rv, eps);
  // The power on a strand, by the counting model: the least k reaching the strand LOD among n.
  const auto strand_power = [&](std::size_t nreads) {
    if (nreads == 0) return 0.0;
    std::size_t k = 1;
    while (k <= nreads && tumour_lod_counts(k, nreads, eps) < x.o.strand_lod) k++;
    if (k > nreads) return 0.0;
    const double ps = f * (1 - eps) + (1 - f) * eps;
    return 1 - calx::special_functions::binomial_cdf(static_cast<double>(k - 1), static_cast<double>(nreads), ps);
  };
  c.power_forward = strand_power(fw);
  c.power_reverse = strand_power(rv);
  const auto gaps = gaps_near(x.tumour, e.position, x.o.gap_half_window, t.reads);
  if (gaps[0] >= x.o.gap_reads || gaps[1] >= x.o.gap_reads) c.filters.push_back("proximal_gap");
  if ((c.power_forward >= x.o.strand_power && c.lod_forward < x.o.strand_lod) || (c.power_reverse >= x.o.strand_power && c.lod_reverse < x.o.strand_lod))
    c.filters.push_back("strand_bias");
  if (n.carriers >= x.o.control_reads || static_cast<double>(n.carriers) >= x.o.control_fraction * static_cast<double>(n.depth)) c.filters.push_back("observed_in_control");
  out.calls.push_back(std::move(c));
}

}  // namespace detail

// The somatic calls of [begin, end) of one contig.
inline SomaticCalls call_somatic_region(const Codes& reference, std::uint32_t contig, std::size_t begin, std::size_t end, std::span<const mapping::Mapping> tumour_maps,
                                        std::span<const variant::ReadView> tumour_reads, std::span<const mapping::Mapping> normal_maps,
                                        std::span<const variant::ReadView> normal_reads, const MutectOptions& o = {}) {
  const variant::Pileup tp = variant::pileup(reference, contig, begin, end, tumour_maps, tumour_reads, o.pileup);
  const variant::Pileup np = variant::pileup(reference, contig, begin, end, normal_maps, normal_reads, o.pileup);
  const detail::Context x{reference, tp, np, tumour_maps, normal_maps, o};
  SomaticCalls out;
  std::vector<variant::BaseCall> tb, nb;
  std::size_t next = 0;
  for (std::size_t p = begin; p < end; p++) {
    if (tp.depth(p) >= o.callable_tumour && np.depth(p) >= o.callable_normal) out.callable++;
    detail::substitution_at(x, p, tb, nb, out);
    std::size_t last = next;
    while (last < tp.indels.size() && tp.indels[last].position == p) last++;
    if (o.indels && last > next) detail::indel_at(x, next, last, out);
    next = last;
  }
  std::stable_sort(out.calls.begin(), out.calls.end(), [](const SomaticCall& a, const SomaticCall& b) { return a.position < b.position; });
  return out;
}

// The somatic calls of a genome from the tumour's and the normal's mapped reads, the same on any
// number of threads.
inline SomaticCalls call_somatic(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const MappedReads& tumour, const MappedReads& normal, MutectOptions o = {}) {
  if (o.window == 0) throw std::invalid_argument("genomics::somatic::call_somatic: a window of at least one base");
  std::sort(o.known_germline.begin(), o.known_germline.end());
  std::vector<detail::Window> windows;
  for (std::uint32_t c = 0; c < contigs.size(); c++)
    for (std::size_t b = 0; b < contigs[c].size(); b += o.window) windows.push_back({c, b, std::min(contigs[c].size(), b + o.window)});
  std::vector<SomaticCalls> slots(windows.size());
  calx::scheduler::parallel_for(pool, windows.size(), 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
    for (std::size_t w = lo; w < hi; w++) {
      const auto& x = windows[w];
      slots[w] = call_somatic_region(contigs[x.contig], x.contig, x.begin, x.end, tumour.overlapping(x.contig, x.begin, x.end), tumour.reads(),
                                     normal.overlapping(x.contig, x.begin, x.end), normal.reads(), o);
    }
  });
  SomaticCalls out;
  for (auto& s : slots) {
    out.callable += s.callable;
    out.sites_weighed += s.sites_weighed;
    for (auto& c : s.calls) out.calls.push_back(std::move(c));
  }
  std::stable_sort(out.calls.begin(), out.calls.end(), [](const SomaticCall& a, const SomaticCall& b) { return a.contig != b.contig ? a.contig < b.contig : a.position < b.position; });
  return out;
}

// A call as the variant module's Call (genotype 0/1), for its evaluation and its VCF.
inline variant::Call as_call(const SomaticCall& s) {
  variant::Call c;
  c.contig = s.contig;
  c.position = s.position;
  c.ref = s.ref;
  c.alts = {s.alt};
  c.qual = std::min(999.0, 10 * s.lod_tumour);
  c.filters = s.filters;
  c.depth = s.tumour_depth;
  c.allele_depths = {s.tumour_depth - s.tumour_alt, s.tumour_alt};
  c.genotype = {0, 1};
  c.indel = s.indel;
  return c;
}

}  // namespace genomics::somatic
