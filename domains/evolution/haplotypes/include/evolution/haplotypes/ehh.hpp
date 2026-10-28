// Haplotype scans for selection: the extended haplotype homozygosity around a core site, by
// allele; its integral over the genetic map; and the integrated haplotype score iHS,
// standardised within bins of the derived allele's frequency.
//
// EHH (Sabeti et al. 2002, Nature 419, 832–837). Let C_a be the haplotypes carrying allele a
// at the core site c, n_a = |C_a| ≥ 2. For a site j, the alleles of C_a at the sites between c
// and j (both included) partition C_a into classes of identical extended haplotypes, and
//
//   EHH_a(j) = Σ_classes C(n_h, 2) / C(n_a, 2),
//
// the probability that two distinct carriers drawn at random are identical over the whole
// interval. EHH_a(c) = 1, and EHH is non-increasing away from the core, since each site can only
// refine the partition. A carrier missing an allele (−1) is counted as different from every
// other from that site on. The walk refines the partition one site at a time and drops the
// classes of one, which never contribute again: O(n_a) per site, less as the classes shrink.
//
// iHH. iHH_a is the integral of EHH_a over the genetic map, on both sides of the core, by the
// trapezoid rule over the sites walked: Σ_i |g_i − g_{i−1}| (EHH_{i−1} + EHH_i)/2 outward from
// g_0 = g_c with EHH_0 = 1. The walk on each side stops at the first site where EHH has fallen
// to the cut-off or below, that last interval included; Voight et al.'s cut-off, 0.05, is the
// default. A side whose walk reaches the end of the data first is truncated, and a truncated
// core is left out of the scan unless it is asked for (its integral would be cut short).
//
// iHS (Voight, Kudaravalli, Wen and Pritchard 2006, PLoS Biology 4, e72). The unstandardised
// score of a core site is ln(iHH_A / iHH_D), A the ancestral allele (0) and D the derived (1).
// Its distribution depends on the derived allele's frequency p, so it is standardised within
// bins of p:
//
//   iHS = (ln(iHH_A / iHH_D) − m_b) / s_b,
//
// with m_b and s_b the mean and standard deviation of the unstandardised scores of the sites in
// the same bin b: those of the scan itself, or a standardisation fitted on others (a genome's
// background) and applied. With s_b's divisor the bin's count, each bin's scores have mean 0
// and variance 1 exactly, by construction. Negative scores are unusually long haplotypes
// carrying the derived allele, the signature of an allele that has risen fast and not yet fixed.
//
// Complexity: EHH and iHH cost O(n L) for n haplotypes and the L sites walked; a scan, the sum
// over its core sites, in parallel over the sites (calx's dynamic loop), the same on any number
// of threads. Standardisation is O(S) for S sites.
//
// Limits: a large gap between sites is integrated over like any interval (implementations of
// the scan often penalise or refuse long gaps; none is applied here); the polarisation of the
// alleles is the caller's (the coalescent's 0 is ancestral); bins are equal widths of p on [0, 1].

#pragma once

#include <evolution/population/genotypes.hpp>

#include <calx/scheduler/parallel_for.hpp>
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

namespace evolution::haplotypes {

using population::Haplotypes;

namespace detail {

[[noreturn]] inline void ehh_fail(const char* where, const std::string& why) {
  throw std::invalid_argument(std::string("evolution::haplotypes::") + where + ": " + why);
}

// Walks from the core in one direction (+1 or −1), calling visit(j, ehh) at each site after the
// core until it returns false or the sites end; returns whether the walk was stopped.
template <class Visit>
bool ehh_walk(const Haplotypes& h, std::size_t core, std::int8_t allele, int direction, Visit&& visit) {
  const std::size_t n = h.sequences;
  std::vector<std::uint32_t> active;
  for (std::size_t i = 0; i < n; i++)
    if (h.alleles[core * n + i] == allele) active.push_back(static_cast<std::uint32_t>(i));
  const double na = static_cast<double>(active.size());
  if (active.size() < 2) return true;
  const double pairs = na * (na - 1) / 2;
  std::vector<std::uint32_t> cls(n, 0), size, next_active;
  std::vector<std::int32_t> id;
  std::size_t classes = 1;
  for (std::size_t step = 1;; step++) {
    if (direction < 0 ? step > core : core + step >= h.sites) return false;
    const std::size_t j = direction < 0 ? core - step : core + step;
    const std::int8_t* row = h.alleles.data() + j * n;
    id.assign(2 * classes, -1);
    size.clear();
    for (const std::uint32_t i : active) {
      const std::int8_t a = row[i];
      if (a < 0) {
        cls[i] = std::numeric_limits<std::uint32_t>::max();
        continue;
      }
      std::int32_t& k = id[2 * cls[i] + static_cast<std::size_t>(a)];
      if (k < 0) {
        k = static_cast<std::int32_t>(size.size());
        size.push_back(0);
      }
      cls[i] = static_cast<std::uint32_t>(k);
      size[static_cast<std::size_t>(k)]++;
    }
    // Keep the carriers in classes of two or more, renumbered.
    std::vector<std::int32_t> renumber(size.size(), -1);
    double homozygous = 0;
    std::size_t kept = 0;
    for (std::size_t k = 0; k < size.size(); k++)
      if (size[k] >= 2) {
        renumber[k] = static_cast<std::int32_t>(kept++);
        homozygous += static_cast<double>(size[k]) * static_cast<double>(size[k] - 1) / 2;
      }
    next_active.clear();
    for (const std::uint32_t i : active)
      if (cls[i] != std::numeric_limits<std::uint32_t>::max() && renumber[cls[i]] >= 0) {
        cls[i] = static_cast<std::uint32_t>(renumber[cls[i]]);
        next_active.push_back(i);
      }
    active.swap(next_active);
    classes = std::max<std::size_t>(kept, 1);
    if (!visit(j, homozygous / pairs)) return true;
  }
}

inline void check_core(const Haplotypes& h, std::size_t core, std::int8_t allele, const char* where) {
  if (h.alleles.size() != h.sites * h.sequences) ehh_fail(where, "the alleles are sites × sequences");
  if (core >= h.sites) ehh_fail(where, "the core is a site of the haplotypes");
  if (allele != 0 && allele != 1) ehh_fail(where, "the core allele is 0 or 1");
}

}  // namespace detail

// The number of haplotypes carrying `allele` at the core.
inline std::size_t carriers(const Haplotypes& h, std::size_t core, std::int8_t allele) {
  detail::check_core(h, core, allele, "carriers");
  return static_cast<std::size_t>(std::count(h.alleles.begin() + static_cast<std::ptrdiff_t>(core * h.sequences),
                                             h.alleles.begin() + static_cast<std::ptrdiff_t>((core + 1) * h.sequences), allele));
}

// EHH_a at one site (1 at the core; NaN with fewer than two carriers).
inline double ehh_at(const Haplotypes& h, std::size_t core, std::int8_t allele, std::size_t site) {
  detail::check_core(h, core, allele, "ehh_at");
  if (site >= h.sites) detail::ehh_fail("ehh_at", "a site of the haplotypes");
  if (carriers(h, core, allele) < 2) return std::numeric_limits<double>::quiet_NaN();
  if (site == core) return 1;
  double e = 0;
  detail::ehh_walk(h, core, allele, site < core ? -1 : 1, [&](std::size_t j, double x) {
    e = x;
    return j != site && x > 0;
  });
  return e;
}

// EHH_a at the sites around the core, in increasing order of site, each side walked until EHH
// has fallen to the cut-off or below (that site included) or the sites end; cut-off 0 walks
// until no pair is left identical.
struct EhhCurve {
  std::size_t core = 0, carriers = 0;
  std::int8_t allele = 0;
  std::vector<std::size_t> sites;
  std::vector<double> ehh;
};

inline EhhCurve ehh(const Haplotypes& h, std::size_t core, std::int8_t allele, double cutoff = 0) {
  detail::check_core(h, core, allele, "ehh");
  EhhCurve c;
  c.core = core;
  c.allele = allele;
  c.carriers = carriers(h, core, allele);
  if (c.carriers < 2) return c;
  std::vector<std::size_t> left_sites;
  std::vector<double> left;
  detail::ehh_walk(h, core, allele, -1, [&](std::size_t j, double x) {
    left_sites.push_back(j);
    left.push_back(x);
    return x > cutoff;
  });
  c.sites.assign(left_sites.rbegin(), left_sites.rend());
  c.ehh.assign(left.rbegin(), left.rend());
  c.sites.push_back(core);
  c.ehh.push_back(1);
  detail::ehh_walk(h, core, allele, 1, [&](std::size_t j, double x) {
    c.sites.push_back(j);
    c.ehh.push_back(x);
    return x > cutoff;
  });
  return c;
}

// ---- iHH ------------------------------------------------------------------------------------------------

struct Ihh {
  double upstream = 0, downstream = 0;  // the integrals towards smaller and larger positions
  bool truncated = false;               // a walk reached the end of the data above the cut-off
  std::size_t carriers = 0;
  double total() const { return upstream + downstream; }
};

// iHH_a by the trapezoid rule over the genetic map (the haplotypes' positions unless given).
inline Ihh ihh(const Haplotypes& h, std::size_t core, std::int8_t allele, std::span<const double> map = {}, double cutoff = 0.05) {
  detail::check_core(h, core, allele, "ihh");
  if (map.empty()) map = h.positions;
  if (map.size() != h.sites) detail::ehh_fail("ihh", "one genetic position for each site");
  Ihh out;
  out.carriers = carriers(h, core, allele);
  if (out.carriers < 2) return out;
  for (const int direction : {-1, 1}) {
    double area = 0, last_g = map[core], last_e = 1;
    const bool stopped = detail::ehh_walk(h, core, allele, direction, [&](std::size_t j, double x) {
      area += std::fabs(map[j] - last_g) * (last_e + x) / 2;
      last_g = map[j];
      last_e = x;
      return x > cutoff;
    });
    if (!stopped) out.truncated = true;
    (direction < 0 ? out.upstream : out.downstream) = area;
  }
  return out;
}

// ---- iHS ------------------------------------------------------------------------------------------------

struct IhsOptions {
  double cutoff = 0.05;         // EHH at which each walk stops
  double min_maf = 0.05;        // core sites with a rarer minor allele are left out
  bool keep_truncated = false;  // keep cores whose walk reached the end of the data
};

struct IhsSite {
  std::size_t site = 0;
  double position = 0;
  double frequency = 0;  // of the derived allele (1) among the haplotypes typed at the core
  double ihh_ancestral = 0, ihh_derived = 0;
  double unstandardised = std::numeric_limits<double>::quiet_NaN();  // ln(iHH_A / iHH_D)
  double ihs = std::numeric_limits<double>::quiet_NaN();             // after standardisation
  bool truncated = false;
};

// ln(iHH_A / iHH_D) at every core site that passes the options, in order of site.
inline std::vector<IhsSite> ihs_scan(const Haplotypes& h, std::span<const double> map = {}, const IhsOptions& o = {},
                                     calx::scheduler::ThreadPool* pool = nullptr) {
  if (map.empty()) map = h.positions;
  if (map.size() != h.sites) detail::ehh_fail("ihs_scan", "one genetic position for each site");
  std::vector<IhsSite> all(h.sites);
  std::vector<std::uint8_t> keep(h.sites, 0);
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t j = begin; j < end; j++) {
      std::size_t typed = 0, derived = 0;
      for (const std::int8_t a : h.site(j))
        if (a >= 0) {
          typed++;
          derived += static_cast<std::size_t>(a);
        }
      if (typed == 0) continue;
      const double p = static_cast<double>(derived) / static_cast<double>(typed);
      if (std::min(p, 1 - p) < o.min_maf || derived < 2 || typed - derived < 2) continue;
      const Ihh a = ihh(h, j, 0, map, o.cutoff), d = ihh(h, j, 1, map, o.cutoff);
      IhsSite s;
      s.site = j;
      s.position = j < h.positions.size() ? h.positions[j] : map[j];
      s.frequency = p;
      s.ihh_ancestral = a.total();
      s.ihh_derived = d.total();
      s.truncated = a.truncated || d.truncated;
      if (s.ihh_ancestral > 0 && s.ihh_derived > 0) s.unstandardised = std::log(s.ihh_ancestral / s.ihh_derived);
      if (s.truncated && !o.keep_truncated) continue;
      all[j] = s;
      keep[j] = 1;
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, h.sites, 16, body);
  else body(0, h.sites, 0);
  std::vector<IhsSite> out;
  for (std::size_t j = 0; j < h.sites; j++)
    if (keep[j]) out.push_back(all[j]);
  return out;
}

// The mean and standard deviation (divisor n) of the unstandardised scores in each of `bins`
// equal bins of the derived frequency on [0, 1].
struct IhsStandardisation {
  std::size_t bins = 20;
  std::vector<double> mean, sd;
  std::vector<std::size_t> count;

  std::size_t bin(double p) const { return std::min(bins - 1, static_cast<std::size_t>(std::max(0.0, p) * static_cast<double>(bins))); }
  // NaN in a bin with fewer than two sites or no spread.
  double apply(double p, double unstandardised) const {
    const std::size_t b = bin(p);
    if (count[b] < 2 || !(sd[b] > 0)) return std::numeric_limits<double>::quiet_NaN();
    return (unstandardised - mean[b]) / sd[b];
  }
};

inline IhsStandardisation fit_ihs_standardisation(std::span<const IhsSite> sites, std::size_t bins = 20) {
  if (bins == 0) detail::ehh_fail("fit_ihs_standardisation", "at least one bin");
  IhsStandardisation z;
  z.bins = bins;
  z.mean.assign(bins, 0);
  z.sd.assign(bins, 0);
  z.count.assign(bins, 0);
  for (const IhsSite& s : sites)
    if (std::isfinite(s.unstandardised)) {
      const std::size_t b = z.bin(s.frequency);
      z.count[b]++;
      z.mean[b] += s.unstandardised;
    }
  for (std::size_t b = 0; b < bins; b++)
    if (z.count[b]) z.mean[b] /= static_cast<double>(z.count[b]);
  for (const IhsSite& s : sites)
    if (std::isfinite(s.unstandardised)) {
      const std::size_t b = z.bin(s.frequency);
      const double d = s.unstandardised - z.mean[b];
      z.sd[b] += d * d;
    }
  for (std::size_t b = 0; b < bins; b++) z.sd[b] = z.count[b] ? std::sqrt(z.sd[b] / static_cast<double>(z.count[b])) : 0;
  return z;
}

// iHS of each site under a fitted standardisation.
inline void standardise(std::span<IhsSite> sites, const IhsStandardisation& z) {
  for (IhsSite& s : sites) s.ihs = std::isfinite(s.unstandardised) ? z.apply(s.frequency, s.unstandardised) : std::numeric_limits<double>::quiet_NaN();
}

// Standardised on the sites themselves; returns the standardisation.
inline IhsStandardisation standardise(std::span<IhsSite> sites, std::size_t bins = 20) {
  const IhsStandardisation z = fit_ihs_standardisation(sites, bins);
  standardise(sites, z);
  return z;
}

}  // namespace evolution::haplotypes
