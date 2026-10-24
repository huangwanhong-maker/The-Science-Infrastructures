// Genotype matrices of biallelic sites: what every statistic of population genetics reads.
//
// A sample of n diploid individuals typed at m biallelic sites is an n × m matrix G of
// alternate-allele counts (dosages) g_ij ∈ {0, 1, 2}, with −1 for a missing call; each site
// has a chromosome and a 0-based position; when the calls are phased, each individual's two
// haplotypes h_{i,0}, h_{i,1} ∈ {0, 1} are kept too, with g_ij = h_{i,0} + h_{i,1}. Haplotypes
// are the same thing for haploid sequences: n sequences at S sites, 0 the reference (or
// ancestral) allele and 1 the alternate (or derived), as the coalescent produces them, with
// real-valued positions (base pairs, or [0, 1) on the coalescent's continuous chromosome).
//
// Both are stored site by site (sites × samples, row-major), as VCF lists them and as every
// statistic here runs: allele counts, Hardy–Weinberg tests, linkage between two sites, the
// spectra and F_ST all read one site's row at a time, and principal components accumulate one
// site's outer product at a time. The matrix in the mathematical sense, samples × sites, is
// G(i, j) = dosage[j · n + i].
//
// Reading VCF is the genomics domain's variant module's (written later): it is to produce
// the dosages site by site, the phased alleles when the GT fields are phased (|), the positions
// (VCF's POS − 1) and the chromosome of each record, and to split or drop multiallelic records.
// The constructors here take those raw arrays and check them.
//
// Allele counts at each site (copies of the alternate allele, and the number of chromosomes
// observed there) are the sufficient statistic of every single-site estimator in diversity.hpp;
// filters by minor-allele frequency and missingness select sites before analysis.
//
// Complexity: every function is one pass over the data it reads, O(n m).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace evolution::population {

inline constexpr std::int8_t MISSING = -1;

// n sequences at S biallelic sites, site by site; positions in increasing order.
struct Haplotypes {
  std::size_t sequences = 0, sites = 0;
  std::vector<std::int8_t> alleles;  // sites × sequences: 0, 1, or −1 missing
  std::vector<double> positions;     // of each site

  std::int8_t operator()(std::size_t sequence, std::size_t site) const { return alleles[site * sequences + sequence]; }
  std::span<const std::int8_t> site(std::size_t j) const { return {alleles.data() + j * sequences, sequences}; }
};

// n diploid individuals at m biallelic sites on one or more chromosomes, site by site, sorted
// by chromosome and then position.
struct Genotypes {
  std::size_t samples = 0, sites = 0;
  std::vector<std::int8_t> dosage;        // sites × samples: copies of the alternate allele, −1 missing
  std::vector<std::int8_t> haplotypes;    // sites × 2·samples when phased (copies 2i, 2i + 1 of sample i), else empty
  std::vector<std::uint64_t> positions;   // 0-based, of each site
  std::vector<std::uint32_t> chromosome;  // of each site: an index into chromosomes
  std::vector<std::string> chromosomes;   // their names
  std::vector<std::string> names;         // of the samples, or empty

  bool phased() const { return !haplotypes.empty(); }
  // G(i, j): the dosage of sample i at site j.
  std::int8_t operator()(std::size_t sample, std::size_t site) const { return dosage[site * samples + sample]; }
  std::span<const std::int8_t> site(std::size_t j) const { return {dosage.data() + j * samples, samples}; }
  // The 2n phased alleles of site j.
  std::span<const std::int8_t> haplotypes_at(std::size_t j) const {
    if (!phased()) throw std::invalid_argument("evolution::population::Genotypes::haplotypes_at: the genotypes are not phased");
    return {haplotypes.data() + j * 2 * samples, 2 * samples};
  }
  // The sites [first, last) of one chromosome.
  std::pair<std::size_t, std::size_t> chromosome_sites(std::uint32_t c) const {
    const auto lo = std::lower_bound(chromosome.begin(), chromosome.end(), c);
    const auto hi = std::upper_bound(chromosome.begin(), chromosome.end(), c);
    return {static_cast<std::size_t>(lo - chromosome.begin()), static_cast<std::size_t>(hi - chromosome.begin())};
  }
};

namespace detail {

inline void check_layout(Genotypes& g) {
  const auto fail = [](const char* why) { throw std::invalid_argument(std::string("evolution::population::Genotypes: ") + why); };
  if (g.dosage.size() != g.samples * g.sites) fail("the dosages are samples × sites");
  if (g.positions.size() != g.sites) fail("one position for each site");
  if (g.chromosome.empty()) g.chromosome.assign(g.sites, 0);
  if (g.chromosome.size() != g.sites) fail("one chromosome for each site");
  if (g.chromosomes.empty()) {
    const std::uint32_t top = g.sites ? *std::max_element(g.chromosome.begin(), g.chromosome.end()) : 0;
    for (std::uint32_t c = 0; c <= top && g.sites; c++) g.chromosomes.push_back(std::to_string(c + 1));
  }
  for (std::size_t j = 0; j < g.sites; j++) {
    if (g.chromosome[j] >= g.chromosomes.size()) fail("a chromosome index names a chromosome");
    if (j && (g.chromosome[j] < g.chromosome[j - 1] || (g.chromosome[j] == g.chromosome[j - 1] && g.positions[j] < g.positions[j - 1])))
      fail("the sites are sorted by chromosome and position");
  }
  for (const std::int8_t x : g.dosage)
    if (x < -1 || x > 2) fail("a dosage is 0, 1, 2 or −1 (missing)");
  if (!g.names.empty() && g.names.size() != g.samples) fail("one name for each sample, or none");
}

}  // namespace detail

// From dosages listed site by site (sites × samples), as VCF lists them. An empty chromosome
// vector puts every site on one chromosome; empty names number the chromosomes from 1.
inline Genotypes genotypes_from_dosages(std::size_t samples, std::vector<std::int8_t> by_site, std::vector<std::uint64_t> positions,
                                        std::vector<std::uint32_t> chromosome = {}, std::vector<std::string> chromosomes = {}) {
  Genotypes g;
  g.samples = samples;
  g.sites = samples ? by_site.size() / samples : positions.size();
  g.dosage = std::move(by_site);
  g.positions = std::move(positions);
  g.chromosome = std::move(chromosome);
  g.chromosomes = std::move(chromosomes);
  detail::check_layout(g);
  return g;
}

// From the matrix in its mathematical layout, samples × sites (row i is sample i).
inline Genotypes genotypes_from_matrix(std::size_t samples, std::size_t sites, std::span<const std::int8_t> by_sample,
                                       std::vector<std::uint64_t> positions, std::vector<std::uint32_t> chromosome = {},
                                       std::vector<std::string> chromosomes = {}) {
  if (by_sample.size() != samples * sites) throw std::invalid_argument("evolution::population::genotypes_from_matrix: samples × sites values");
  std::vector<std::int8_t> by_site(samples * sites);
  for (std::size_t i = 0; i < samples; i++)
    for (std::size_t j = 0; j < sites; j++) by_site[j * samples + i] = by_sample[i * sites + j];
  return genotypes_from_dosages(samples, std::move(by_site), std::move(positions), std::move(chromosome), std::move(chromosomes));
}

// From phased alleles listed site by site (sites × 2·samples: copies 2i and 2i + 1 are sample
// i's); a dosage is missing when either copy is.
inline Genotypes genotypes_from_haplotypes(std::size_t samples, std::vector<std::int8_t> by_site, std::vector<std::uint64_t> positions,
                                           std::vector<std::uint32_t> chromosome = {}, std::vector<std::string> chromosomes = {}) {
  if (samples == 0 || by_site.size() % (2 * samples) != 0)
    throw std::invalid_argument("evolution::population::genotypes_from_haplotypes: sites × 2·samples alleles");
  const std::size_t sites = by_site.size() / (2 * samples);
  std::vector<std::int8_t> dosage(sites * samples);
  for (std::size_t j = 0; j < sites; j++)
    for (std::size_t i = 0; i < samples; i++) {
      const std::int8_t a = by_site[j * 2 * samples + 2 * i], b = by_site[j * 2 * samples + 2 * i + 1];
      if (a < -1 || a > 1 || b < -1 || b > 1) throw std::invalid_argument("evolution::population::genotypes_from_haplotypes: an allele is 0, 1 or −1");
      dosage[j * samples + i] = (a < 0 || b < 0) ? MISSING : static_cast<std::int8_t>(a + b);
    }
  Genotypes g = genotypes_from_dosages(samples, std::move(dosage), std::move(positions), std::move(chromosome), std::move(chromosomes));
  g.haplotypes = std::move(by_site);
  return g;
}

// Diploid individuals from haploid sequences, sequences 2i and 2i + 1 making individual i;
// positions x on [0, 1) (the coalescent's) become ⌊x · length⌋ base pairs on one chromosome.
inline Genotypes diploids(const Haplotypes& h, std::uint64_t length, std::string chromosome_name = "1") {
  if (h.sequences % 2) throw std::invalid_argument("evolution::population::diploids: an even number of sequences");
  std::vector<std::uint64_t> pos(h.sites);
  for (std::size_t j = 0; j < h.sites; j++) {
    const double x = h.positions[j] * static_cast<double>(length);
    pos[j] = static_cast<std::uint64_t>(std::clamp(x, 0.0, static_cast<double>(length ? length - 1 : 0)));
  }
  return genotypes_from_haplotypes(h.sequences / 2, h.alleles, std::move(pos), {}, {std::move(chromosome_name)});
}

// The haplotypes of phased genotypes, as 2n sequences.
inline Haplotypes haplotypes_of(const Genotypes& g) {
  if (!g.phased()) throw std::invalid_argument("evolution::population::haplotypes_of: the genotypes are not phased");
  Haplotypes h;
  h.sequences = 2 * g.samples;
  h.sites = g.sites;
  h.alleles = g.haplotypes;
  h.positions.assign(g.positions.begin(), g.positions.end());
  return h;
}

// ---- allele counts ------------------------------------------------------------------------------------

// A view of the counts of a set of sites: copies of the alternate (or derived) allele, and
// chromosomes observed, at each.
struct CountsView {
  std::span<const std::uint32_t> alt, total;
  std::size_t sites() const { return alt.size(); }
  CountsView slice(std::size_t first, std::size_t last) const { return {alt.subspan(first, last - first), total.subspan(first, last - first)}; }
};

struct AlleleCounts {
  std::vector<std::uint32_t> alt, total;
  std::size_t sites() const { return alt.size(); }
  CountsView view() const { return {alt, total}; }
  CountsView slice(std::size_t first, std::size_t last) const { return view().slice(first, last); }
  operator CountsView() const { return view(); }
};

// The counts of one site's dosages, over the samples listed (all when the list is empty).
inline std::pair<std::uint32_t, std::uint32_t> count_site(std::span<const std::int8_t> dosages, std::span<const std::size_t> which = {}) {
  std::uint32_t alt = 0, total = 0;
  const auto add = [&](std::int8_t x) {
    if (x < 0) return;
    alt += static_cast<std::uint32_t>(x);
    total += 2;
  };
  if (which.empty())
    for (const std::int8_t x : dosages) add(x);
  else
    for (const std::size_t i : which) add(dosages[i]);
  return {alt, total};
}

inline AlleleCounts allele_counts(const Genotypes& g, std::span<const std::size_t> which = {}) {
  AlleleCounts c;
  c.alt.resize(g.sites);
  c.total.resize(g.sites);
  for (std::size_t j = 0; j < g.sites; j++) std::tie(c.alt[j], c.total[j]) = count_site(g.site(j), which);
  return c;
}

inline AlleleCounts allele_counts(const Haplotypes& h, std::span<const std::size_t> which = {}) {
  AlleleCounts c;
  c.alt.resize(h.sites);
  c.total.resize(h.sites);
  for (std::size_t j = 0; j < h.sites; j++) {
    std::uint32_t alt = 0, total = 0;
    const auto row = h.site(j);
    const auto add = [&](std::int8_t x) {
      if (x < 0) return;
      alt += static_cast<std::uint32_t>(x);
      total++;
    };
    if (which.empty())
      for (const std::int8_t x : row) add(x);
    else
      for (const std::size_t i : which) add(row[i]);
    c.alt[j] = alt;
    c.total[j] = total;
  }
  return c;
}

// The samples of each population, from a label per sample (labels 0 … K − 1).
inline std::vector<std::vector<std::size_t>> members(std::span<const std::uint32_t> population) {
  std::vector<std::vector<std::size_t>> m;
  for (std::size_t i = 0; i < population.size(); i++) {
    if (population[i] >= m.size()) m.resize(population[i] + 1);
    m[population[i]].push_back(i);
  }
  return m;
}

// The alternate allele's frequency at each site; NaN where nothing was observed.
inline std::vector<double> frequencies(CountsView c) {
  std::vector<double> f(c.sites());
  for (std::size_t j = 0; j < c.sites(); j++)
    f[j] = c.total[j] ? static_cast<double>(c.alt[j]) / c.total[j] : std::numeric_limits<double>::quiet_NaN();
  return f;
}

// The minor allele's frequency at each site, min(p, 1 − p); NaN where nothing was observed.
inline std::vector<double> minor_allele_frequencies(CountsView c) {
  std::vector<double> f = frequencies(c);
  for (double& x : f) x = std::min(x, 1 - x);
  return f;
}

// The fraction of samples missing at each site, and of sites missing in each sample.
inline std::vector<double> missing_by_site(const Genotypes& g) {
  std::vector<double> f(g.sites);
  for (std::size_t j = 0; j < g.sites; j++) {
    std::size_t k = 0;
    for (const std::int8_t x : g.site(j)) k += x < 0;
    f[j] = g.samples ? static_cast<double>(k) / static_cast<double>(g.samples) : 0;
  }
  return f;
}
inline std::vector<double> missing_by_sample(const Genotypes& g) {
  std::vector<std::size_t> k(g.samples, 0);
  for (std::size_t j = 0; j < g.sites; j++) {
    const auto row = g.site(j);
    for (std::size_t i = 0; i < g.samples; i++) k[i] += row[i] < 0;
  }
  std::vector<double> f(g.samples);
  for (std::size_t i = 0; i < g.samples; i++) f[i] = g.sites ? static_cast<double>(k[i]) / static_cast<double>(g.sites) : 0;
  return f;
}

// The sites whose minor-allele frequency is at least min_maf and whose missing fraction is at
// most max_missing, in order. A monomorphic site has frequency 0, so min_maf > 0 drops it.
inline std::vector<std::size_t> passing_sites(const Genotypes& g, double min_maf, double max_missing = 1) {
  std::vector<std::size_t> keep;
  for (std::size_t j = 0; j < g.sites; j++) {
    const auto [alt, total] = count_site(g.site(j));
    if (total == 0) continue;
    const double p = static_cast<double>(alt) / total, maf = std::min(p, 1 - p);
    const double missing = 1 - static_cast<double>(total) / (2.0 * static_cast<double>(g.samples));
    if (maf >= min_maf && missing <= max_missing) keep.push_back(j);
  }
  return keep;
}

// The genotypes of the sites listed (in increasing order), and of the samples listed.
inline Genotypes select_sites(const Genotypes& g, std::span<const std::size_t> sites) {
  Genotypes r;
  r.samples = g.samples;
  r.sites = sites.size();
  r.chromosomes = g.chromosomes;
  r.names = g.names;
  r.dosage.reserve(r.sites * g.samples);
  if (g.phased()) r.haplotypes.reserve(r.sites * 2 * g.samples);
  for (const std::size_t j : sites) {
    if (j >= g.sites) throw std::invalid_argument("evolution::population::select_sites: a site of the matrix");
    const auto row = g.site(j);
    r.dosage.insert(r.dosage.end(), row.begin(), row.end());
    if (g.phased()) {
      const auto h = g.haplotypes_at(j);
      r.haplotypes.insert(r.haplotypes.end(), h.begin(), h.end());
    }
    r.positions.push_back(g.positions[j]);
    r.chromosome.push_back(g.chromosome[j]);
  }
  detail::check_layout(r);
  return r;
}

inline Genotypes select_samples(const Genotypes& g, std::span<const std::size_t> samples) {
  Genotypes r;
  r.samples = samples.size();
  r.sites = g.sites;
  r.positions = g.positions;
  r.chromosome = g.chromosome;
  r.chromosomes = g.chromosomes;
  for (const std::size_t i : samples) {
    if (i >= g.samples) throw std::invalid_argument("evolution::population::select_samples: a sample of the matrix");
    if (!g.names.empty()) r.names.push_back(g.names[i]);
  }
  r.dosage.resize(r.sites * r.samples);
  if (g.phased()) r.haplotypes.resize(r.sites * 2 * r.samples);
  for (std::size_t j = 0; j < g.sites; j++)
    for (std::size_t k = 0; k < samples.size(); k++) {
      r.dosage[j * r.samples + k] = g.dosage[j * g.samples + samples[k]];
      if (g.phased()) {
        r.haplotypes[j * 2 * r.samples + 2 * k] = g.haplotypes[j * 2 * g.samples + 2 * samples[k]];
        r.haplotypes[j * 2 * r.samples + 2 * k + 1] = g.haplotypes[j * 2 * g.samples + 2 * samples[k] + 1];
      }
    }
  return r;
}

}  // namespace evolution::population
