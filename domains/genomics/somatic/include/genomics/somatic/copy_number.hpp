// Purity, ploidy and allele-specific copy number of a tumour from its reads and its normal's:
// logR in bins and at the normal's heterozygous sites, the B-allele fraction at those sites, their
// joint segmentation (ASPCF), and ASCAT's fit of purity ρ and ploidy ψ with integer copy numbers of
// both alleles in every segment (Van Loo et al. 2010, "Allele-specific copy number analysis of
// tumors", PNAS 107, 16910–16915, eqs. 1–4 and SI; the R package ASCAT, VanLoo-lab/ascat,
// ascat.runAscat.R and ascat.aspcf.R, whose rules are followed where the paper is not explicit).
//
// - The model (eqs. 1–2). A sample of tumour cells (fraction ρ) of copy numbers n_A, n_B at a
//   site and normal cells (1 + 1) has DNA 2(1 − ρ) + ρ(n_A + n_B) there, and on average
//   2(1 − ρ) + ρψ, ψ the tumour's mean copy number; so
//     r = γ log₂ [(2(1 − ρ) + ρ(n_A + n_B)) / (2(1 − ρ) + ρψ)],   b = (1 − ρ + ρ n_B) / (2 − 2ρ + ρ(n_A + n_B)),
//   r the logR (γ = 1 for sequencing; 0.55 for Illumina's 109K arrays) and b the B-allele fraction
//   at a germline heterozygous site. Inverted (eqs. 3–4; the code's form), with T = 2(1 − ρ) + ρψ,
//     n_A = (ρ − 1 − (b − 1) 2^{r/γ} T) / ρ,   n_B = (ρ − 1 + b 2^{r/γ} T) / ρ.
// - logR here: the tumour's depth over the normal's in bins of `bin` bases, each relative to its
//   mean over the genome, log₂; the normal has two copies everywhere, so the ratio is the tumour
//   sample's DNA over its mean, eq. 1 with ψ the length-weighted mean copy number. A site carries
//   the mean logR of the bins between the midpoints to its neighbouring sites (ascat.aspcf.R's
//   logRaveraged). The B allele is the site's ALT; BAF = ALT reads / (REF + ALT reads).
// - Heterozygous sites: the normal's SNV calls of genotype 0/1 that pass (variant/caller.hpp),
//   with at least `min_depth` reads of the two alleles in each sample.
// - ASPCF (Nilsen et al. 2012, BMC Genomics 13, 591; ascat.aspcf.R): one segmentation of the two
//   series at the sites of a contig, logR and the mirrored BAF min(b, 1 − b), minimising
//     Σ_segments [Σ (r − r̄)²/σ_r² + Σ (b̃ − b̃̄)²/σ_b²] + γ_pcf · (segments − 1),
//   segments of at least `min_sites` sites (6) and penalty γ_pcf = 70 (ASCAT's defaults), by
//   tyche's PELT on its SumCost of two mean costs, exact (ASCAT runs the same recursion in
//   overlapping windows of 1000 sites, an approximation). σ² by tyche's robust_variance (the MAD
//   of first differences; ASCAT's getMad takes the MAD of deviations from a running median of 25).
//   A segment's BAF is 0.5 + μ, μ the mean of |b − 0.5|, set to 0 (a balanced segment) when
//   √(σ_b² + μ²) < 2σ_b (fastAspcf's rule). ASCAT's winsorisation of logR is not done.
// - The fit. For ψ on 1.00, 1.05, …, 6.00 and ρ on 0.10, 0.11, …, 1.05 (min_ploidy − 0.5 to
//   max_ploidy + 0.5 and min_purity to max_purity), the distance
//     d(ψ, ρ) = Σ_segments |n_m − max(round(n_m), 0)|² · s · w,
//   n_m the minor allele (of n_A, n_B, the one of smaller sum over segments), s the segment's
//   sites, w = 0.05 for balanced segments (b = 0.5) and 1 otherwise (create_distance_matrix); the
//   goodness of fit 100 (1 − d / d_max), d_max = Σ 0.25 s w, the distance when every segment is
//   half a copy from an integer. The "sunrise" is d over the grid. Candidates are strict local
//   minima in 7 × 7 windows of the grid; ploidy recomputed as Σ (n_A + n_B) s / Σ s; in turn (the
//   first pass that finds any): (1) ploidy in (min, max), ρ ≥ 0.2, goodness > 80, the fraction of
//   sites with an allele at 0 above 0.02, the sample aberrant; (2) without the zeros' condition but
//   ploidy in (1.7, 2.3) and the zeros among unbalanced segments above 0.1; (3) ρ ≤ 1 included at
//   the border, any of the zeros' and parity conditions; (4) ploidy in (1.7, 2.3) only. The least
//   distance among the candidates is the fit, ρ held to 1. Not aberrant: unbalanced segments over
//   at most 3% of the sites and none over 0.5% (MINABB, MINABBREGION).
// - Integer copy numbers of a segment at the fit: n_A, n_B from eqs. 3–4, negative values moved to
//   the other allele (or both 0), rounded half to even (R's round); in a balanced segment an odd
//   total more than half a copy above or below the rounded one adds or removes a copy (limitround).
//
// Complexity: O(m) for the bins of m mappings, O(d ℓ) per site, PELT's expected O(n) per contig,
// O(|ψ| |ρ| S) for the grid on S segments (101 × 96 × S). Limits: clonal copy number only (a
// subclonal segment fits half-integers badly; Battenberg would fit two states), no GC or
// replication-timing correction of logR, autosomes (no sex chromosomes), no germline homozygous
// stretches inferred.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/somatic/coverage.hpp>
#include <genomics/variant/caller.hpp>

#include <tyche/time_series/changepoint.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::somatic {

struct CopyNumberOptions {
  std::size_t bin = 1000;
  int min_mapq = 20;
  double gamma = 1;              // ASCAT's γ: 1 for sequencing
  std::size_t min_depth = 10;    // reads of the two alleles at a site, in each sample
  double penalty = 70;           // ASPCF's
  std::size_t min_sites = 6;     // ASPCF's kmin
  double min_ploidy = 1.5, max_ploidy = 5.5, min_purity = 0.1, max_purity = 1.05;
  variant::PileupOptions pileup;
};

struct HetSite {
  std::uint32_t contig = 0;
  std::size_t position = 0;
  Code ref = 0, alt = 0;
  std::uint32_t tumour_ref = 0, tumour_alt = 0, normal_ref = 0, normal_alt = 0;
  double baf = 0.5;          // the tumour's
  double normal_baf = 0.5;
  double logr = 0;
  std::size_t segment = 0;
};

struct LogRBin {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  double tumour = 0, normal = 0;  // mean depths
  double logr = std::numeric_limits<double>::quiet_NaN();  // NaN where the normal has no reads
  std::size_t segment = SIZE_MAX;
};

struct CnSegment {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;          // bases: from the midpoints to the neighbouring segments' sites
  std::size_t first = 0, last = 0;         // its sites [first, last) in the profile's list
  double logr = 0, baf = 0.5;              // r̄, 0.5 + μ
  double n_a_raw = 0, n_b_raw = 0;         // eqs. 3–4 at the fit
  int n_a = 1, n_b = 1;                    // integers
  std::size_t sites() const { return last - first; }
  int major() const { return std::max(n_a, n_b); }
  int minor() const { return std::min(n_a, n_b); }
  bool balanced() const { return baf == 0.5; }
};

struct AscatFit {
  bool found = false;
  double purity = 0, psi = 0;      // the grid point
  double ploidy = 0;               // Σ (n_A + n_B) s / Σ s at the grid point (raw)
  double goodness = 0;             // %
  double distance = 0;
  int pass = 0;                    // which of the four passes found it
  bool nonaberrant = false;
  std::vector<double> psi_grid, rho_grid;
  std::vector<double> distances;   // d(ψ_i, ρ_j) at i · |ρ| + j
  struct Optimum {
    double psi = 0, rho = 0, ploidy = 0, goodness = 0, distance = 0;
  };
  std::vector<Optimum> optima;     // the candidates of the pass that found the fit
};

struct CopyNumberProfile {
  std::vector<LogRBin> bins;
  std::vector<HetSite> sites;
  std::vector<CnSegment> segments;
  AscatFit fit;
  double ploidy = 0;               // the integer copy numbers' mean over the segments' bases
};

// ---- ASCAT's equations -------------------------------------------------------------------------------

// Eq. 1: the logR of copy numbers n_A, n_B at purity ρ and ploidy ψ.
inline double ascat_logr(double na, double nb, double rho, double psi, double gamma = 1) {
  return gamma * std::log2((2 * (1 - rho) + rho * (na + nb)) / (2 * (1 - rho) + rho * psi));
}
// Eq. 2: the B-allele fraction.
inline double ascat_baf(double na, double nb, double rho) { return (1 - rho + rho * nb) / (2 - 2 * rho + rho * (na + nb)); }
// Eqs. 3–4: the copy numbers of a logR and a BAF.
inline double ascat_na(double r, double b, double rho, double psi, double gamma = 1) {
  return (rho - 1 - (b - 1) * std::pow(2.0, r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho;
}
inline double ascat_nb(double r, double b, double rho, double psi, double gamma = 1) {
  return (rho - 1 + b * std::pow(2.0, r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho;
}

// R's round: to the nearest integer, halves to even.
inline double round_half_even(double x) {
  const double f = std::floor(x), d = x - f;
  if (d > 0.5) return f + 1;
  if (d < 0.5) return f;
  return std::fmod(f, 2.0) == 0 ? f : f + 1;
}

// ---- sites, bins and logR ------------------------------------------------------------------------------

// The normal's heterozygous SNVs (PASS, genotype 0/1), read in both samples.
inline std::vector<HetSite> heterozygous_sites(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const std::vector<variant::Call>& normal_calls,
                                               const MappedReads& tumour, const MappedReads& normal, const CopyNumberOptions& o = {}) {
  std::vector<Site> where;
  std::vector<HetSite> cand;
  for (const auto& c : normal_calls) {
    if (!c.passed() || c.indel || c.ref.size() != 1 || c.alts.size() != 1 || c.alts[0].size() != 1 || c.genotype[0] != 0 || c.genotype[1] != 1) continue;
    HetSite h;
    h.contig = c.contig;
    h.position = c.position;
    h.ref = sequence::dna_code(c.ref[0]);
    h.alt = sequence::dna_code(c.alts[0][0]);
    if (h.ref > 3 || h.alt > 3) continue;
    cand.push_back(h);
    where.push_back({c.contig, c.position});
  }
  const auto tc = allele_counts(pool, contigs, where, tumour, o.pileup), nc = allele_counts(pool, contigs, where, normal, o.pileup);
  std::vector<HetSite> out;
  for (std::size_t i = 0; i < cand.size(); i++) {
    HetSite h = cand[i];
    h.tumour_ref = tc[i][h.ref];
    h.tumour_alt = tc[i][h.alt];
    h.normal_ref = nc[i][h.ref];
    h.normal_alt = nc[i][h.alt];
    const std::uint32_t td = h.tumour_ref + h.tumour_alt, nd = h.normal_ref + h.normal_alt;
    if (td < o.min_depth || nd < o.min_depth) continue;
    h.baf = static_cast<double>(h.tumour_alt) / td;
    h.normal_baf = static_cast<double>(h.normal_alt) / nd;
    out.push_back(h);
  }
  std::sort(out.begin(), out.end(), [](const HetSite& a, const HetSite& b) { return a.contig != b.contig ? a.contig < b.contig : a.position < b.position; });
  return out;
}

// The bins with the tumour's and the normal's depths and their logR.
inline std::vector<LogRBin> logr_bins(const std::vector<Codes>& contigs, const MappedReads& tumour, const MappedReads& normal, const CopyNumberOptions& o = {}) {
  const auto t = depth_bins(contigs, tumour, o.bin, o.min_mapq), n = depth_bins(contigs, normal, o.bin, o.min_mapq);
  std::vector<LogRBin> out(t.size());
  double ts = 0, ns = 0, len = 0;
  for (std::size_t i = 0; i < t.size(); i++) {
    out[i] = {t[i].contig, t[i].begin, t[i].end, t[i].depth, n[i].depth};
    if (n[i].depth > 0) {
      const double w = static_cast<double>(t[i].end - t[i].begin);
      ts += t[i].depth * w;
      ns += n[i].depth * w;
      len += w;
    }
  }
  if (!(ts > 0 && ns > 0)) return out;
  const double tm = ts / len, nm = ns / len;
  for (auto& b : out)
    if (b.normal > 0 && b.tumour > 0) b.logr = std::log2((b.tumour / tm) / (b.normal / nm));
  return out;
}

// Each site's logR: the mean of the bins whose centres lie between the midpoints to its neighbours.
inline void site_logr(std::vector<HetSite>& sites, const std::vector<LogRBin>& bins) {
  std::size_t lo = 0;
  for (std::size_t i = 0; i < sites.size(); i++) {
    const auto& s = sites[i];
    const bool first = i == 0 || sites[i - 1].contig != s.contig, last = i + 1 == sites.size() || sites[i + 1].contig != s.contig;
    const double a = first ? 0 : (static_cast<double>(sites[i - 1].position) + static_cast<double>(s.position)) / 2;
    const double b = last ? std::numeric_limits<double>::infinity() : (static_cast<double>(s.position) + static_cast<double>(sites[i + 1].position)) / 2;
    while (lo < bins.size() && (bins[lo].contig < s.contig || (bins[lo].contig == s.contig && static_cast<double>(bins[lo].end) <= a))) lo++;
    double sum = 0;
    std::size_t k = 0, own = SIZE_MAX;
    for (std::size_t j = lo; j < bins.size() && bins[j].contig == s.contig && static_cast<double>(bins[j].begin) < b; j++) {
      const double centre = (static_cast<double>(bins[j].begin) + static_cast<double>(bins[j].end)) / 2;
      if (bins[j].begin <= s.position && s.position < bins[j].end) own = j;
      if (centre >= a && centre < b && std::isfinite(bins[j].logr)) {
        sum += bins[j].logr;
        k++;
      }
    }
    sites[i].logr = k ? sum / static_cast<double>(k) : (own != SIZE_MAX && std::isfinite(bins[own].logr) ? bins[own].logr : 0);
  }
}

// ---- ASPCF ----------------------------------------------------------------------------------------------

// The joint segmentation of the sites' logR and mirrored BAF, contig by contig.
inline std::vector<CnSegment> aspcf(std::vector<HetSite>& sites, const std::vector<std::size_t>& contig_lengths, const CopyNumberOptions& o = {}) {
  namespace ts = tyche::time_series;
  std::vector<CnSegment> out;
  for (std::size_t a = 0; a < sites.size();) {
    std::size_t b = a;
    while (b < sites.size() && sites[b].contig == sites[a].contig) b++;
    const std::size_t n = b - a;
    std::vector<double> r(n), m(n);
    for (std::size_t i = 0; i < n; i++) {
      r[i] = sites[a + i].logr;
      m[i] = std::min(sites[a + i].baf, 1 - sites[a + i].baf);
    }
    const double vr = ts::robust_variance(r), vb = ts::robust_variance(m);
    std::vector<std::size_t> cuts;  // first sites of segments, relative
    if (n >= 2 * o.min_sites) {
      ts::CostOptions cr, cb;
      cr.variance = vr;
      cb.variance = vb;
      const ts::SumCost cost({ts::SegmentCost(r, cr), ts::SegmentCost(m, cb)});
      const auto seg = ts::pelt(cost, ts::penalty_terms(ts::Penalty::manual, n, 2, o.penalty), o.min_sites);
      cuts = seg.changepoints;
    }
    cuts.insert(cuts.begin(), 0);
    cuts.push_back(n);
    const double sd2 = std::sqrt(vb);
    for (std::size_t k = 0; k + 1 < cuts.size(); k++) {
      CnSegment s;
      s.contig = sites[a].contig;
      s.first = a + cuts[k];
      s.last = a + cuts[k + 1];
      double sr = 0, mu = 0;
      for (std::size_t i = s.first; i < s.last; i++) {
        sr += sites[i].logr;
        mu += std::fabs(sites[i].baf - 0.5);
        sites[i].segment = out.size();
      }
      const double len = static_cast<double>(s.last - s.first);
      s.logr = sr / len;
      mu /= len;
      if (std::sqrt(sd2 * sd2 + mu * mu) < 2 * sd2) mu = 0;
      s.baf = 0.5 + mu;
      const bool first = k == 0, last = k + 2 == cuts.size();
      s.begin = first ? 0 : (sites[s.first - 1].position + sites[s.first].position) / 2 + 1;
      s.end = last ? (s.contig < contig_lengths.size() ? contig_lengths[s.contig] : sites[s.last - 1].position + 1) : (sites[s.last - 1].position + sites[s.last].position) / 2 + 1;
      out.push_back(s);
    }
    a = b;
  }
  return out;
}

// ---- the fit ----------------------------------------------------------------------------------------

namespace detail {

struct Grid {
  std::vector<double> psi, rho;
};
inline Grid ascat_grid(const CopyNumberOptions& o) {
  Grid g;
  const auto steps = [](double a, double b, double h, std::vector<double>& v) {
    for (long k = 0;; k++) {
      const double x = std::round((a + static_cast<double>(k) * h) * 100) / 100;
      if (x > b + 1e-9) break;
      v.push_back(x);
    }
  };
  steps(o.min_ploidy - 0.5, o.max_ploidy + 0.5, 0.05, g.psi);
  steps(std::round(o.min_purity * 100) / 100, std::round(o.max_purity * 100) / 100, 0.01, g.rho);
  return g;
}

struct Evaluation {
  double ploidy = 0, zero = 0, zero_unbalanced = 0, odd_even = 0;
};
inline Evaluation evaluate(const std::vector<CnSegment>& s, double rho, double psi, double gamma) {
  Evaluation e;
  double len = 0, len_unbal = 0;
  for (const auto& x : s) {
    const double na = ascat_na(x.logr, x.baf, rho, psi, gamma), nb = ascat_nb(x.logr, x.baf, rho, psi, gamma);
    const double w = static_cast<double>(x.sites()), wu = x.balanced() ? 0 : w;
    const double ra = round_half_even(na), rb = round_half_even(nb);
    e.ploidy += (na + nb) * w;
    e.zero += ((ra == 0) + (rb == 0)) * w;
    e.zero_unbalanced += ((ra == 0) + (rb == 0)) * wu;
    const bool ae = std::fmod(std::fabs(ra), 2.0) == 0, be = std::fmod(std::fabs(rb), 2.0) == 0;
    e.odd_even += (ae != be) * w;
    len += w;
    len_unbal += wu;
  }
  e.ploidy /= len;
  e.zero /= len;
  e.odd_even /= len;
  e.zero_unbalanced = len_unbal > 0 ? e.zero_unbalanced / len_unbal : 0;
  return e;
}

}  // namespace detail

// ASCAT's distance at one (ψ, ρ).
inline double ascat_distance(const std::vector<CnSegment>& s, double rho, double psi, double gamma = 1) {
  double sa = 0, sb = 0;
  std::vector<double> na(s.size()), nb(s.size());
  for (std::size_t i = 0; i < s.size(); i++) {
    na[i] = ascat_na(s[i].logr, s[i].baf, rho, psi, gamma);
    nb[i] = ascat_nb(s[i].logr, s[i].baf, rho, psi, gamma);
    sa += na[i];
    sb += nb[i];
  }
  const auto& minor = sa < sb ? na : nb;
  double d = 0;
  for (std::size_t i = 0; i < s.size(); i++) {
    const double e = minor[i] - std::max(round_half_even(minor[i]), 0.0);
    d += e * e * static_cast<double>(s[i].sites()) * (s[i].balanced() ? 0.05 : 1);
  }
  return d;
}

// The fit of ρ and ψ to the segments, and their integer copy numbers (written into them).
inline AscatFit ascat(std::vector<CnSegment>& segments, const CopyNumberOptions& o = {}) {
  AscatFit f;
  if (segments.empty()) return f;
  const auto g = detail::ascat_grid(o);
  f.psi_grid = g.psi;
  f.rho_grid = g.rho;
  const std::size_t I = g.psi.size(), J = g.rho.size();
  f.distances.assign(I * J, 0);
  for (std::size_t i = 0; i < I; i++)
    for (std::size_t j = 0; j < J; j++) f.distances[i * J + j] = ascat_distance(segments, g.rho[j], g.psi[i], o.gamma);
  double dmax = 0, all = 0, unbal = 0, biggest = 0;
  for (const auto& s : segments) {
    const double w = static_cast<double>(s.sites());
    dmax += 0.25 * w * (s.balanced() ? 0.05 : 1);
    all += w;
    if (!s.balanced()) {
      unbal += w;
      biggest = std::max(biggest, w);
    }
  }
  f.nonaberrant = unbal / all <= 0.03 && biggest / all <= 0.005;
  constexpr double MINRHO = 0.2, MINGOODNESS = 80, MINZERO = 0.02, MINZEROABB = 0.1, MINODDEVEN = 0.05, STRICT_LO = 1.7, STRICT_HI = 2.3;
  std::vector<double> d = f.distances;
  const auto search = [&](int pass) {
    std::vector<AscatFit::Optimum> found;
    if (I < 7 || J < 7) return found;
    for (std::size_t i = 3; i + 3 < I; i++)
      for (std::size_t j = 3; j + 3 < J; j++) {
        const double m = d[i * J + j];
        bool strict = true;
        for (std::size_t a = i - 3; a <= i + 3 && strict; a++)
          for (std::size_t b = j - 3; b <= j + 3; b++)
            if ((a != i || b != j) && !(d[a * J + b] > m)) { strict = false; break; }
        if (!strict) continue;
        const double psi = g.psi[i], rho = g.rho[j];
        const detail::Evaluation e = detail::evaluate(segments, rho, psi, o.gamma);
        const double good = (1 - m / dmax) * 100;
        bool keep = false;
        if (pass == 1) keep = !f.nonaberrant && e.ploidy > o.min_ploidy && e.ploidy < o.max_ploidy && rho >= MINRHO && good > MINGOODNESS && e.zero > MINZERO;
        if (pass == 2) keep = e.ploidy > STRICT_LO && e.ploidy < STRICT_HI && rho >= MINRHO && good > MINGOODNESS && e.zero_unbalanced > MINZEROABB;
        if (pass == 3)
          keep = !f.nonaberrant && e.ploidy > o.min_ploidy && e.ploidy < o.max_ploidy && rho >= MINRHO && good > MINGOODNESS &&
                 (e.zero_unbalanced > MINZEROABB || e.zero > MINZERO || e.odd_even > MINODDEVEN);
        if (pass == 4) keep = e.ploidy > STRICT_LO && e.ploidy < STRICT_HI && rho >= MINRHO && good > MINGOODNESS;
        if (keep) found.push_back({psi, rho, e.ploidy, good, m});
      }
    return found;
  };
  const bool strict_allowed = o.min_ploidy < STRICT_HI && o.max_ploidy > STRICT_LO;
  std::vector<AscatFit::Optimum> found = search(1);
  f.pass = 1;
  if (found.empty() && strict_allowed) { found = search(2); f.pass = 2; }
  if (found.empty()) {
    for (std::size_t j = 0; j < J; j++)
      if (g.rho[j] > 1)
        for (std::size_t i = 0; i < I; i++) d[i * J + j] = 1e20;
    found = search(3);
    f.pass = 3;
  }
  if (found.empty() && strict_allowed) { found = search(4); f.pass = 4; }
  f.optima = found;
  if (found.empty()) { f.pass = 0; return f; }
  const auto best = *std::min_element(found.begin(), found.end(), [](const AscatFit::Optimum& a, const AscatFit::Optimum& b) { return a.distance < b.distance; });
  f.found = true;
  f.psi = best.psi;
  f.purity = std::min(1.0, best.rho);
  f.ploidy = best.ploidy;
  f.goodness = best.goodness;
  f.distance = best.distance;
  // Integer copy numbers (runAscat's segment loop).
  const double rho = f.purity, psi = f.psi;
  for (auto& s : segments) {
    double na = ascat_na(s.logr, s.baf, rho, psi, o.gamma), nb = ascat_nb(s.logr, s.baf, rho, psi, o.gamma);
    if (na + nb < 0) { na = 0; nb = 0; }
    else if (na < 0) { nb = na + nb; na = 0; }
    else if (nb < 0) { na = na + nb; nb = 0; }
    s.n_a_raw = na;
    s.n_b_raw = nb;
    double ra = round_half_even(na), rb = round_half_even(nb);
    if (s.balanced()) {
      if (na + nb > ra + rb + 0.5) ra += 1;
      else if (na + nb < ra + rb - 0.5) rb -= 1;
    }
    s.n_a = static_cast<int>(ra);
    s.n_b = static_cast<int>(rb);
  }
  return f;
}

// The whole profile: bins, sites, segments, fit; each bin and site labelled with its segment.
inline CopyNumberProfile copy_number(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const MappedReads& tumour, const MappedReads& normal,
                                     const std::vector<variant::Call>& normal_calls, const CopyNumberOptions& o = {}) {
  CopyNumberProfile p;
  p.bins = logr_bins(contigs, tumour, normal, o);
  p.sites = heterozygous_sites(pool, contigs, normal_calls, tumour, normal, o);
  site_logr(p.sites, p.bins);
  std::vector<std::size_t> lengths;
  for (const auto& c : contigs) lengths.push_back(c.size());
  p.segments = aspcf(p.sites, lengths, o);
  p.fit = ascat(p.segments, o);
  std::size_t s = 0;
  double cn = 0, len = 0;
  for (auto& b : p.bins) {
    const std::size_t centre = (b.begin + b.end) / 2;
    while (s < p.segments.size() && (p.segments[s].contig < b.contig || (p.segments[s].contig == b.contig && p.segments[s].end <= centre))) s++;
    if (s < p.segments.size() && p.segments[s].contig == b.contig && p.segments[s].begin <= centre) b.segment = s;
  }
  for (const auto& g : p.segments) {
    cn += static_cast<double>(g.n_a + g.n_b) * static_cast<double>(g.end - g.begin);
    len += static_cast<double>(g.end - g.begin);
  }
  p.ploidy = len > 0 ? cn / len : 0;
  return p;
}

}  // namespace genomics::somatic
