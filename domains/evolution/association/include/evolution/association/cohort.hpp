// Cohorts with truth: genotypes with the linkage of real populations, a trait of a chosen
// heritability (quantitative, and case–control by a liability threshold), and the expression of genes
// under cis-acting variants, some of which act on the trait. Simulated, and labelled so: every
// association method of this module is measured against what the simulation put in.
//
// THE GENOTYPES. Each region (a stretch of chromosome, simulated independently of the others, so
// that regions play chromosomes) is a founder population's haplotypes from the coalescent with
// recombination (population/coalescent.hpp: Hudson's algorithm, θ = 4N₀μ and ρ = 4N₀r for the region),
// whose linkage disequilibrium decays with distance as it does in a panmictic population. The founders
// are split into K groups, one for each population, and each individual of population k carries two
// haplotypes that are mosaics of its population's founders: the copying model of Li and Stephens
// (haplotypes/copying.hpp: the chain whose jumps go to a founder drawn uniformly, at rate `switches`
// per haplotype over the region, without miscopying), as HAPGEN2 makes new haplotypes from a reference
// panel (Su, Marchini and Donnelly 2011, Bioinformatics 27, 2304–2305). The mosaic adds the recent
// recombination of a population descended from its founders; drawing each population's founders
// apart gives it its own allele frequencies (F_ST ≈ 1/(2F) for F founder haplotypes, a founder
// effect), which principal components then find; an admixed individual copies from every population's
// founders. Sites keep the coalescent's infinite-sites alleles (0 ancestral, 1 derived), at
// ⌊x · length⌋ base pairs (a site landing on an occupied base is dropped), and those whose minor allele
// frequency in the cohort is below min_maf are dropped.
//
// THE TRAIT (the conventions of GCTA's --simu-qt; Yang et al. 2011, AJHG 88, 76–82). Direct causal
// variants are drawn uniformly among the sites, with effects b_j ~ N(0, 1) on standardised genotypes
// z_ij = (g_ij − 2p_j)/√(2p_j(1 − p_j)), so that rarer variants have larger effects per allele (the
// α = −1 model); the direct genetic value Σ_j b_j z_ij is scaled to variance h² − h²_med in the cohort.
// Mediating genes act on the trait through their expression E_g (unit variance), with effects ±a scaled
// so that their cis-genetic part has variance h²_med; the residual is Gaussian with the variance that
// makes the trait's core variance 1, so that the genetic value G = direct + mediated-genetic has
// variance h². On the core the cohort adds covariates (sex, age) and, as confounding, a mean shift per
// population. The case–control trait is the indicator that the liability (the trait, standardised)
// exceeds its upper `prevalence` quantile in the cohort (Falconer 1965's threshold model).
//
// THE EXPRESSION. Each gene has a transcription start site (TSS) in a region and causal variants among
// the common sites (minor allele frequency ≥ 0.05) within cis_window of it: none, one, or between two
// and max_causal, in the proportions asked; their effects, N(0, 1) on standardised genotypes, are scaled
// to the gene's cis-heritability (uniform between cis_h2_low and cis_h2_high). Hidden factors (shared
// by every gene, as batch and cell-type composition are) take factor_share of each gene's variance, sex
// a little, and Gaussian noise the rest, so each gene has unit variance in the cohort. The first
// `reference` individuals are the expression panel (as GTEx is), the others the association sample.
//
// DETERMINISM. Each draw is a calx stream keyed by (seed, purpose, object): a region's coalescent by
// its index, each haplotype's mosaic by (region, haplotype), the architecture by purpose; regions are
// simulated in parallel and assembled in order, so a cohort is a function of its seed alone.
//
// COMPLEXITY. The coalescent of K·F founders per region (milliseconds for a few hundred with ρ in the
// hundreds), then O(2n S) for the mosaics of n individuals at S sites, and O(n m) for the traits.

#pragma once

#include <evolution/haplotypes/copying.hpp>
#include <evolution/population/coalescent.hpp>
#include <evolution/population/genotypes.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace evolution::association {

namespace pop = evolution::population;

// The purposes of this module's streams.
inline constexpr std::uint32_t COHORT = 0x434f4852, MOSAIC = 0x4d4f5341, ARCHITECTURE = 0x41524348, PHENOTYPE = 0x50484e4f;

struct CohortSpec {
  std::size_t individuals = 3000;
  std::size_t regions = 4;                // independent stretches of chromosome, one chromosome each
  std::uint64_t region_length = 1'000'000;  // base pairs
  double theta = 300;                     // 4N₀μ over a region, of the founders' coalescent
  double rho = 300;                       // 4N₀r over a region
  std::size_t populations = 2;            // K
  std::size_t founders = 150;             // founder haplotypes of each population
  double admixed = 0;                     // the fraction of individuals copying from every population's founders
  double switches = 20;                   // expected switches of a mosaic over a region
  double min_maf = 0.01;
};

struct StudySpec {
  CohortSpec cohort;
  // The trait.
  double heritability = 0.4;              // h²: the genetic value's share of the trait's core variance
  std::size_t causal = 20;                // direct causal variants
  double mediated = 0.05;                 // h²_med: the share through the cis-genetic expression of the mediating genes
  std::size_t mediating = 4;              // genes whose expression acts on the trait
  double population_effect = 0.25;        // mean shift per population index, in units of the core sd (confounding)
  double sex_effect = 0.2, age_effect = 0.01;  // per unit of sex (0/1) and per year of age
  double prevalence = 0.2;                // of the case–control trait
  // The expression.
  std::size_t genes = 40;
  double none = 0.25, single = 0.45;      // fractions of genes with no and with one causal variant (the rest several)
  std::size_t max_causal = 4;
  double cis_h2_low = 0.05, cis_h2_high = 0.4;
  std::uint64_t cis_window = 100'000;     // causal variants lie within this of the TSS
  std::size_t factors = 3;                // hidden factors shared by the genes
  double factor_share = 0.15;             // their share of each gene's variance
  double expression_sex = 0.1;            // the effect of sex on expression, per unit
  std::size_t reference = 500;            // individuals whose expression is measured
};

struct Gene {
  std::string name;
  std::uint32_t chromosome = 0;           // the region
  std::uint64_t tss = 0;                  // 0-based base pair
  std::vector<std::size_t> causal;        // sites of its cis-eQTLs (indices into the cohort's genotypes)
  std::vector<double> effects;            // per allele (dosage scale)
  double cis_h2 = 0;                      // realised: the variance of its cis-genetic part
  bool mediating = false;
  double trait_effect = 0;                // a_g, the trait's change per unit of expression
};

struct Study {
  bool simulated = true;                  // always: these are simulated individuals, never real ones
  std::uint64_t seed = 0;
  pop::Genotypes genotypes;               // phased; one chromosome per region
  std::size_t individuals = 0, populations = 0;
  std::vector<std::uint32_t> population;  // per individual; K for an admixed one
  std::vector<double> sex, age;           // covariates
  std::vector<std::size_t> founder_sites; // per region: sites segregating among the founders
  std::vector<Gene> genes;
  std::vector<double> expression;         // genes × individuals: the measured value (the truth, for every individual)
  std::vector<double> expression_genetic; // genes × individuals: its cis-genetic part
  std::vector<double> factors;            // individuals × factors: the hidden factors
  std::vector<double> trait;              // quantitative, every individual
  std::vector<double> liability;          // the standardised trait
  std::vector<double> cases;              // 0/1, the case–control trait
  std::vector<double> genetic;            // G: the trait's genetic value
  std::vector<std::size_t> causal;        // the direct causal sites, increasing
  std::vector<double> effects;            // their per-allele effects
  std::vector<std::size_t> trait_causal;  // every site acting on the trait: the direct ones and the mediating genes' eQTLs, increasing
  double heritability = 0;                // realised var G / var core
  std::vector<std::size_t> reference, sample;  // the expression panel and the association sample (indices)

  // The dosages of site j over the individuals listed, as doubles.
  std::vector<double> dosages(std::size_t j, std::span<const std::size_t> who) const {
    std::vector<double> d(who.size());
    const auto row = genotypes.site(j);
    for (std::size_t k = 0; k < who.size(); k++) d[k] = row[who[k]];
    return d;
  }
  std::span<const double> gene_expression(std::size_t g) const { return {expression.data() + g * individuals, individuals}; }
};

namespace detail {

inline double mean_of(std::span<const double> x) {
  double s = 0;
  for (const double v : x) s += v;
  return x.empty() ? 0 : s / static_cast<double>(x.size());
}

inline double variance_of(std::span<const double> x) {
  const double m = mean_of(x);
  double s = 0;
  for (const double v : x) s += (v - m) * (v - m);
  return x.size() > 1 ? s / static_cast<double>(x.size() - 1) : 0;
}

// One region: the founders' haplotypes and the individuals' mosaics, as haplotypes sites × 2n.
struct Region {
  std::vector<std::uint64_t> positions;
  std::vector<std::int8_t> alleles;  // sites × 2n
  std::size_t founder_sites = 0;
};

inline Region simulate_region(const CohortSpec& c, const std::vector<std::uint32_t>& population, std::uint64_t seed, std::uint32_t region) {
  const std::size_t K = c.populations, F = c.founders, n = population.size();
  pop::CoalescentSpec cs;
  cs.samples = K * F;
  cs.theta = c.theta;
  cs.rho = c.rho;
  const auto genealogy = pop::coalescent(cs, seed, region);
  const auto h = pop::haplotypes(genealogy);
  Region out;
  out.founder_sites = h.sites;
  if (h.sites == 0) return out;
  // Sites at distinct base pairs.
  std::vector<std::size_t> keep;
  std::uint64_t last = 0;
  for (std::size_t j = 0; j < h.sites; j++) {
    const auto bp = static_cast<std::uint64_t>(std::min(h.positions[j] * static_cast<double>(c.region_length), static_cast<double>(c.region_length - 1)));
    if (!keep.empty() && bp <= last) continue;
    keep.push_back(j);
    out.positions.push_back(bp);
    last = bp;
  }
  const std::size_t S = keep.size();
  // The jumps of the copying chain over these sites: copying_chain's r_j = 1 − e^{−ρ d_j/K'} on a panel of
  // K' haplotypes, with ρ = switches · K', is 1 − e^{−switches·d_j} whatever the panel, so one chain of
  // all the founders serves every population's (ε = 0: no miscopying).
  std::vector<double> jump;
  {
    pop::Haplotypes panel;
    panel.sequences = h.sequences;
    panel.sites = S;
    panel.alleles.resize(S * h.sequences);
    panel.positions.resize(S);
    for (std::size_t j = 0; j < S; j++) {
      const auto row = h.site(keep[j]);
      std::copy(row.begin(), row.end(), panel.alleles.begin() + static_cast<std::ptrdiff_t>(j * h.sequences));
      panel.positions[j] = h.positions[keep[j]];
    }
    haplotypes::CopyingModel model;
    model.rho = c.switches * static_cast<double>(h.sequences);
    model.epsilon = 0;
    jump = haplotypes::copying_chain(panel, model).jump;
  }
  out.alleles.assign(S * 2 * n, 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t copy = 0; copy < 2; copy++) {
      const std::size_t hap = 2 * i + copy;
      calx::random::Stream r(seed, MOSAIC, region, static_cast<std::uint32_t>(hap));
      const bool admixed = population[i] >= K;
      const std::size_t base = admixed ? 0 : population[i] * F, width = admixed ? K * F : F;
      std::size_t z = 0;
      for (std::size_t j = 0; j < S; j++) {
        if (j == 0 || r.bernoulli(jump[j])) z = base + r.below(width);
        out.alleles[j * 2 * n + hap] = h(z, keep[j]);
      }
    }
  return out;
}

}  // namespace detail

// A cohort with its genotypes, trait and expression, and the truth of each.
inline Study simulate_study(const StudySpec& spec, std::uint64_t seed, calx::scheduler::ThreadPool* pool = nullptr) {
  const CohortSpec& c = spec.cohort;
  if (c.individuals < 10 || c.regions == 0 || c.populations == 0 || c.founders < 2 || !(c.admixed >= 0 && c.admixed <= 1) || !(c.switches >= 0) ||
      c.region_length < 1000)
    throw std::invalid_argument("evolution::association::simulate_study: at least 10 individuals, a region, a population of two founders");
  if (!(spec.heritability >= 0 && spec.heritability < 1) || !(spec.mediated >= 0 && spec.mediated <= spec.heritability) || !(spec.prevalence > 0 && spec.prevalence < 1))
    throw std::invalid_argument("evolution::association::simulate_study: 0 ≤ h²_med ≤ h² < 1 and a prevalence in (0, 1)");
  if (spec.reference >= c.individuals) throw std::invalid_argument("evolution::association::simulate_study: the expression panel is part of the cohort");
  const std::size_t n = c.individuals, K = c.populations;
  Study st;
  st.seed = seed;
  st.individuals = n;
  st.populations = K;
  // Individuals: their population (in blocks), sex and age.
  {
    calx::random::Stream r(seed, COHORT, 0, 0);
    st.population.resize(n);
    st.sex.resize(n);
    st.age.resize(n);
    const auto admixed = static_cast<std::size_t>(std::llround(c.admixed * static_cast<double>(n)));
    for (std::size_t i = 0; i < n; i++) {
      st.population[i] = i < n - admixed ? static_cast<std::uint32_t>(i * K / std::max<std::size_t>(1, n - admixed)) : static_cast<std::uint32_t>(K);
      st.sex[i] = r.bernoulli(0.5) ? 1 : 0;
      st.age[i] = std::round(r.uniform(30, 70));
    }
    // The panel and the sample: a seeded permutation, so that both hold every population.
    std::vector<std::size_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    for (std::size_t i = n; i > 1; i--) std::swap(perm[i - 1], perm[r.below(i)]);
    st.reference.assign(perm.begin(), perm.begin() + static_cast<std::ptrdiff_t>(spec.reference));
    st.sample.assign(perm.begin() + static_cast<std::ptrdiff_t>(spec.reference), perm.end());
    std::sort(st.reference.begin(), st.reference.end());
    std::sort(st.sample.begin(), st.sample.end());
  }
  // The regions, in parallel, assembled in order.
  std::vector<detail::Region> regions(c.regions);
  const auto one = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t r = begin; r < end; r++) regions[r] = detail::simulate_region(c, st.population, seed, static_cast<std::uint32_t>(r));
  };
  if (pool) calx::scheduler::parallel_for(*pool, c.regions, 1, one);
  else one(0, c.regions, 0);
  {
    std::vector<std::int8_t> alleles;
    std::vector<std::uint64_t> positions;
    std::vector<std::uint32_t> chromosome;
    std::vector<std::string> names;
    for (std::size_t r = 0; r < c.regions; r++) {
      names.push_back(std::to_string(r + 1));
      st.founder_sites.push_back(regions[r].founder_sites);
      const std::size_t S = regions[r].positions.size();
      for (std::size_t j = 0; j < S; j++) {
        const std::int8_t* row = &regions[r].alleles[j * 2 * n];
        std::size_t alt = 0;
        for (std::size_t k = 0; k < 2 * n; k++) alt += static_cast<std::size_t>(row[k]);
        const double f = static_cast<double>(alt) / (2.0 * static_cast<double>(n));
        if (std::min(f, 1 - f) < c.min_maf) continue;
        alleles.insert(alleles.end(), row, row + 2 * n);
        positions.push_back(regions[r].positions[j]);
        chromosome.push_back(static_cast<std::uint32_t>(r));
      }
      regions[r] = {};
    }
    st.genotypes = pop::genotypes_from_haplotypes(n, std::move(alleles), std::move(positions), std::move(chromosome), std::move(names));
  }
  const std::size_t m = st.genotypes.sites;
  if (m == 0) throw std::invalid_argument("evolution::association::simulate_study: no site passes the frequency filter");
  // Standardised genotypes of a site.
  std::vector<double> freq(m);
  for (std::size_t j = 0; j < m; j++) {
    const auto [alt, total] = pop::count_site(st.genotypes.site(j));
    freq[j] = static_cast<double>(alt) / static_cast<double>(total);
  }
  const auto add_standardised = [&](std::size_t j, double b, std::vector<double>& out) {
    const double p = freq[j], sd = std::sqrt(2 * p * (1 - p));
    const auto row = st.genotypes.site(j);
    for (std::size_t i = 0; i < n; i++) out[i] += b * (row[i] - 2 * p) / sd;
  };
  calx::random::Stream arch(seed, ARCHITECTURE, 0, 0);
  // The genes: TSS, causal variants and effects.
  st.genes.resize(spec.genes);
  st.expression_genetic.assign(spec.genes * n, 0.0);
  for (std::size_t g = 0; g < spec.genes; g++) {
    Gene& gene = st.genes[g];
    gene.name = "GENE" + std::to_string(g + 1);
    gene.chromosome = static_cast<std::uint32_t>(g % c.regions);
    const auto [lo, hi] = st.genotypes.chromosome_sites(gene.chromosome);
    const double margin = std::min(0.2, static_cast<double>(spec.cis_window) / static_cast<double>(c.region_length));
    gene.tss = static_cast<std::uint64_t>(arch.uniform(margin, 1 - margin) * static_cast<double>(c.region_length));
    const double u = arch.uniform();
    std::size_t k = u < spec.none ? 0 : u < spec.none + spec.single ? 1 : 2 + arch.below(std::max<std::size_t>(1, spec.max_causal - 1));
    std::vector<std::size_t> candidates;
    for (std::size_t j = lo; j < hi; j++) {
      const auto d = st.genotypes.positions[j] > gene.tss ? st.genotypes.positions[j] - gene.tss : gene.tss - st.genotypes.positions[j];
      if (d <= spec.cis_window && std::min(freq[j], 1 - freq[j]) >= 0.05) candidates.push_back(j);
    }
    k = std::min(k, candidates.size());
    for (std::size_t t = 0; t < k; t++) {
      const std::size_t pick = t + arch.below(candidates.size() - t);
      std::swap(candidates[t], candidates[pick]);
      gene.causal.push_back(candidates[t]);
    }
    std::sort(gene.causal.begin(), gene.causal.end());
    if (gene.causal.empty()) continue;
    const double h2 = arch.uniform(spec.cis_h2_low, spec.cis_h2_high);
    std::vector<double> b(gene.causal.size());
    std::vector<double> cis(n, 0.0);
    for (std::size_t t = 0; t < b.size(); t++) {
      b[t] = arch.normal();
      add_standardised(gene.causal[t], b[t], cis);
    }
    const double v = detail::variance_of(cis), scale = v > 0 ? std::sqrt(h2 / v) : 0;
    for (std::size_t t = 0; t < b.size(); t++) {
      const double p = freq[gene.causal[t]];
      gene.effects.push_back(b[t] * scale / std::sqrt(2 * p * (1 - p)));
    }
    for (std::size_t i = 0; i < n; i++) st.expression_genetic[g * n + i] = cis[i] * scale;
    gene.cis_h2 = detail::variance_of(std::span<const double>(&st.expression_genetic[g * n], n));
  }
  // The hidden factors and the measured expression, of unit variance.
  st.factors.assign(n * spec.factors, 0.0);
  {
    calx::random::Stream r(seed, PHENOTYPE, 1, 0);
    for (double& f : st.factors) f = r.normal();
  }
  st.expression.assign(spec.genes * n, 0.0);
  for (std::size_t g = 0; g < spec.genes; g++) {
    calx::random::Stream r(seed, PHENOTYPE, 2, static_cast<std::uint32_t>(g));
    std::vector<double> loading(spec.factors);
    for (double& l : loading) l = r.normal();
    std::vector<double> hidden(n, 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t f = 0; f < spec.factors; f++) hidden[i] += loading[f] * st.factors[i * spec.factors + f];
    const double hv = detail::variance_of(hidden), hs = hv > 0 ? std::sqrt(spec.factor_share / hv) : 0;
    const double sexv = spec.expression_sex * spec.expression_sex * 0.25;
    const double rest = std::max(0.02, 1 - st.genes[g].cis_h2 - spec.factor_share - sexv);
    for (std::size_t i = 0; i < n; i++)
      st.expression[g * n + i] = st.expression_genetic[g * n + i] + hs * hidden[i] + spec.expression_sex * (st.sex[i] - 0.5) + std::sqrt(rest) * r.normal();
  }
  // The trait: direct causal variants, mediating genes, noise, covariates and confounding.
  std::vector<double> direct(n, 0.0), mediated(n, 0.0), mediated_genetic(n, 0.0);
  {
    std::vector<std::size_t> all(m);
    std::iota(all.begin(), all.end(), 0);
    const std::size_t M = std::min(spec.causal, m);
    for (std::size_t t = 0; t < M; t++) std::swap(all[t], all[t + arch.below(m - t)]);
    st.causal.assign(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(M));
    std::sort(st.causal.begin(), st.causal.end());
    std::vector<double> b(M);
    for (std::size_t t = 0; t < M; t++) add_standardised(st.causal[t], b[t] = arch.normal(), direct);
    const double target = spec.heritability - spec.mediated, v = detail::variance_of(direct);
    const double scale = v > 0 && target > 0 ? std::sqrt(target / v) : 0;
    for (double& x : direct) x *= scale;
    for (std::size_t t = 0; t < M; t++) {
      const double p = freq[st.causal[t]];
      st.effects.push_back(b[t] * scale / std::sqrt(2 * p * (1 - p)));
    }
    // Mediating genes: among those with a cis-eQTL.
    std::vector<std::size_t> eligible;
    for (std::size_t g = 0; g < spec.genes; g++)
      if (!st.genes[g].causal.empty()) eligible.push_back(g);
    const std::size_t G = std::min(spec.mediating, eligible.size());
    for (std::size_t t = 0; t < G; t++) std::swap(eligible[t], eligible[t + arch.below(eligible.size() - t)]);
    for (std::size_t t = 0; t < G; t++) {
      Gene& gene = st.genes[eligible[t]];
      gene.mediating = true;
      gene.trait_effect = arch.bernoulli(0.5) ? 1 : -1;
      for (std::size_t i = 0; i < n; i++) {
        mediated_genetic[i] += gene.trait_effect * st.expression_genetic[eligible[t] * n + i];
        mediated[i] += gene.trait_effect * st.expression[eligible[t] * n + i];
      }
    }
    const double mv = detail::variance_of(mediated_genetic), ms = mv > 0 && spec.mediated > 0 ? std::sqrt(spec.mediated / mv) : 0;
    for (std::size_t i = 0; i < n; i++) {
      mediated_genetic[i] *= ms;
      mediated[i] *= ms;
    }
    for (Gene& gene : st.genes)
      if (gene.mediating) gene.trait_effect *= ms;
  }
  st.genetic.resize(n);
  std::vector<double> core(n);
  {
    std::vector<double> part(n);
    for (std::size_t i = 0; i < n; i++) part[i] = direct[i] + mediated[i];
    const double rest = std::max(0.05, 1 - detail::variance_of(part));
    calx::random::Stream r(seed, PHENOTYPE, 3, 0);
    for (std::size_t i = 0; i < n; i++) {
      st.genetic[i] = direct[i] + mediated_genetic[i];
      core[i] = part[i] + std::sqrt(rest) * r.normal();
    }
  }
  st.heritability = detail::variance_of(st.genetic) / detail::variance_of(core);
  st.trait_causal = st.causal;
  for (const Gene& gene : st.genes)
    if (gene.mediating) st.trait_causal.insert(st.trait_causal.end(), gene.causal.begin(), gene.causal.end());
  std::sort(st.trait_causal.begin(), st.trait_causal.end());
  st.trait_causal.erase(std::unique(st.trait_causal.begin(), st.trait_causal.end()), st.trait_causal.end());
  st.trait.resize(n);
  for (std::size_t i = 0; i < n; i++) {
    const double shift = st.population[i] >= K ? 0.5 * static_cast<double>(K - 1) : static_cast<double>(st.population[i]);
    st.trait[i] = core[i] + spec.population_effect * shift + spec.sex_effect * st.sex[i] + spec.age_effect * (st.age[i] - 50);
  }
  {
    const double mu = detail::mean_of(st.trait), sd = std::sqrt(detail::variance_of(st.trait));
    st.liability.resize(n);
    for (std::size_t i = 0; i < n; i++) st.liability[i] = (st.trait[i] - mu) / sd;
    std::vector<double> sorted = st.liability;
    const auto cut = static_cast<std::size_t>(std::floor((1 - spec.prevalence) * static_cast<double>(n)));
    std::nth_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(cut), sorted.end());
    const double threshold = sorted[cut];
    st.cases.resize(n);
    for (std::size_t i = 0; i < n; i++) st.cases[i] = st.liability[i] >= threshold ? 1 : 0;
  }
  return st;
}

}  // namespace evolution::association
