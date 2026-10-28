# population

*Status: draft. A module of the evolution domain.* Population genetics: the variation of a sample at biallelic sites and what it says about the population it came from. Genotype matrices and their allele counts; tests of Hardy–Weinberg equilibrium; linkage disequilibrium and its decay; diversity, site frequency spectra and F_ST; the coalescent, which generates samples backwards in time with mutation, changing size and recombination; Wright–Fisher populations forward in time with selection, mutation and migration; a selective sweep forward in time, whose haplotypes the haplotype scans are checked on; and population structure by principal components and by ancestry proportions. Serves use (tested estimators), learning (every estimator checked against the process that generates its data) and research. The theory is [theory/population.md](../theory/population.md).

## 1. Concepts

- **Two directions of time.** A sample's genealogy is generated backwards by the coalescent (`coalescent.hpp`), and a population's allele frequency forwards by Wright–Fisher (`wright_fisher.hpp`). The coalescent is the Wright–Fisher model's genealogy in the limit of large N, so each checks the other's theory: the estimators of θ against the coalescent's expectations, the fixation probability against Kimura's diffusion.
- **Units.** Coalescent time is in units of 2N₀ generations (a pair coalesces at rate 1); ms's are 4N₀, half as many (`ms_time`). θ = 4N₀μ and ρ = 4N₀r are for the whole simulated chromosome, which is [0, 1).
- **Allele counts are the sufficient statistic** of every single-site estimator: θ_W, π, θ_H, the spectra and Hudson's F_ST read only the counts of the alternate (or derived) allele and of the chromosomes typed. Weir and Cockerham's θ also reads heterozygosity, so it takes genotypes.
- **Missing data** make the sample size vary by site. θ_W and π sum each site's own term; the spectra, and Tajima's D and Fay and Wu's H read from them, are projected to a common sample size by hypergeometric sampling.
- **Ratios of averages.** F_ST over many sites is Σ numerator / Σ denominator, never a mean of per-site ratios (Bhatia et al. 2013); so are windows.
- **Determinism.** Every draw is a calx stream keyed by its purpose and its object (replicate, site, individual); parallel work is in fixed shards merged in order. A genealogy, a trajectory, a structured sample, an ancestry estimate and every summary are the same on any number of threads.

## 2. Content

| Header | Content |
|---|---|
| `genotypes.hpp` | `Genotypes` (samples × sites dosages 0/1/2/−1, stored site by site; positions and chromosomes; phased haplotypes optional) from raw arrays (`genotypes_from_dosages`, `genotypes_from_matrix`, `genotypes_from_haplotypes`); `Haplotypes` (n sequences, as the coalescent makes them), `diploids`, `haplotypes_of`; `AlleleCounts` and `CountsView`, `allele_counts` (over all samples or a subset), `members`, `frequencies`, `minor_allele_frequencies`, `missing_by_site`, `missing_by_sample`, `passing_sites`, `select_sites`, `select_samples` |
| `hardy_weinberg.hpp` | `GenotypeCounts`, `genotype_counts`; `hwe_chi_square` (χ², its P value and F̂); `heterozygote_law` (the exact conditional law of the heterozygote count by Wigginton's recurrences); `hwe_exact` (P_HWE, P_low, P_high) for a site or a whole matrix |
| `linkage.hpp` | `Linkage` (p_A, p_B, p_AB, D, D′, r²); `linkage_phased`; `PackedHaplotypes` (64 haplotypes a word, a pair in popcounts); `GenotypeTable`, `linkage_unphased` (the EM of haplotype frequencies, from three starts), `haplotype_log_likelihood`; `ld_decay` (mean r² in bins of distance, phased or unphased, in parallel); `ld_matrix` for drawing |
| `diversity.hpp` | `watterson_a`, `watterson_b`; `segregating_sites`, `theta_watterson`, `nucleotide_diversity`, `theta_h`; `tajima_constants`, `tajima_d`; `fay_wu_h`; `site_frequency_spectrum` (unfolded, projected), `fold`, `summarise`, `polarise`; `hudson_site`, `fst_hudson`; `weir_cockerham_site`, `fst_weir_cockerham`; `windows`, `window_statistics` |
| `coalescent.hpp` | `Demography` (epochs of constant or exponential size; `exponential`, `bottleneck`; Λ and its inverse); `CoalescentSpec`; `coalescent` (Hudson's algorithm: a `Genealogy` of nodes with times, edges with intervals, and infinite-sites mutations); `for_each_tree` (the sweep of marginal trees), `tree_count`, `tmrca_at`; `haplotypes`; `Tree`, `marginal_trees`, `tree_at` (trees for drawing: nodes, times, children, mutations on branches); `coalescent_replicates` (in parallel); `expected_tmrca`, `expected_tree_length`, `ms_time` |
| `wright_fisher.hpp` | `WrightFisherSpec` (N, demes, s, h, μ, ν, m, p₀); `wright_fisher_trajectory` (for drawing); `wright_fisher_replicates` (a calx task graph of shards; fixation, settling time, mean p, H_S, H_T, F_ST, the trace); `kimura_fixation`, `fixation_probability` (the diffusion with dominance, by quadrature), `heterozygosity_decay`, `island_fst`, `island_fst_exact` |
| `sweep.hpp` | `SweepSpec` (N, n, θ, ρ, s, position, stop frequency), `sweep(spec, seed, replicate)` → `Sweep` (the sample's haplotypes at their polymorphic sites, the selected site among them, the beneficial allele's trajectory, generations, attempts): a hard sweep in a Wright–Fisher population of 2N chromosomes whose standing variation is a coalescent sample, with genic selection and crossovers, attempts that lose the allele restarted; no new mutation during the sweep (theory: [haplotypes.md §5](../theory/haplotypes.md)) |
| `structure.hpp` | `kinship` (M Mᵀ/m on vector lanes, in parallel); `pca` (Patterson's components, eigenvalues on his scale, Tracy–Widom statistics, effective markers); `admixture` (the EM of ancestry, SQUAREM-accelerated, restarts), `admixture_log_likelihood`; `simulate_structure` (Balding–Nichols populations, unadmixed and Dirichlet-admixed individuals) |
| `vcf.hpp` | `genotypes_from_vcf(records, samples, report)`, `genotypes_from_vcf(vcf, report)` → `Genotypes` from the genomics domain's VCF records: biallelic records as sites (others skipped and counted), dosages from GT (missing alleles, other ploidies, unknown alleles and absent calls missing), positions POS − 1, chromosomes by first appearance and sites sorted, haplotypes kept when every site is phased; `VcfConversion{records, sites, not_biallelic, other_ploidy, unknown_allele, absent, phased_records, phased}` |

## 3. Dependencies

- **Foundations:** calx: `random` (streams: exponential, Poisson, binomial, gamma, beta), `linalg` (the symmetric eigensolver), `special_functions` (χ² tails, log Γ), `quadrature` (the diffusion's integral), `scheduler` (thread pool, dynamic loops, task graphs), `cpu` (vector lanes); `verification` for the checks.
- **Domains:** genomics (`variant`, for VCF).

## 4. Checks

Each check program runs in seconds: `node domains/evolution/tests/run-tests.mjs` builds and runs them all.

| Check | What it verifies |
|---|---|
| genotypes | The three constructors agree; malformed input refused; counts, frequencies, missingness and filters on a matrix worked by hand; selecting sites and samples keeps the haplotypes; diploids from the coalescent's sequences |
| hardy_weinberg | Wigginton et al. (2005) Table 1 to six decimals (N = 100, 21 minor copies: probabilities, P_HWE, P_high, P_low); the paper's rejection regions and true sizes (0.010 exact, 0.070 χ²); the law against enumeration of every allele arrangement (N ≤ 10) and against eq. 1 for a million individuals; the exact test never above its size; χ² = N F̂²; size kept in one simulated population, the Wahlund effect in two pooled |
| linkage | D, D′, r² by hand, repulsion, missing calls; packed equals plain; the EM exact without double heterozygotes, at the likelihood's maximum on a grid, and recovering the phased r² of coalescent diploids; the decay the same on 1 and 8 threads and from phased genotypes as from haplotypes |
| diversity | S, θ_W, π, θ_H, H by hand and π by pairs; Tajima's D tending to Tajima's D_min and D_max; spectra, folding, projection against enumeration, polarisation; Hudson's and Weir–Cockerham's components on scikit-allel's worked example; both F_ST estimating Balding–Nichols' F; windows partitioning the sites, their sums adding up |
| coalescent | E[T_MRCA] = 2(1 − 1/n), Var[T_MRCA], E[L] = 2a_n; E[S] = θa_n, Var[S] = θa_n + θ²b_n; E[π] = θ, Var[π] = b₁θ + b₂θ², Var[π − S/a₁] = c₁θ + c₂θ²; E[ξ_i] = θ/i; E[θ_H] = θ; Tajima's D near 0; no four gametes on one tree; P(T₂ > t) = e^{−Λ(t)} under growth, a bottleneck and three epochs; D negative after an expansion; 1 and 8 threads alike (each expectation within four standard errors) |
| recombination | The marginal trees partition the chromosome, each Kingman's at a fixed position; trees grow in number with ρ; Cov[T_A, T_B] = (ρ + 18)/(ρ² + 13ρ + 18) for two sequences at two loci; four gametes; r² falling with distance; the drawn trees binary, ordered in time, with their mutations |
| wright_fisher | Neutral fixation p and absorption time −4N(p ln p + q ln q); Kimura's u(p) at 4Ns = −1, 1, 2, 4; the diffusion with a recessive advantage; heterozygosity (1 − 1/(2N))^t; island F_ST ≈ 1/(1 + 4Nm) at 4Nm = 1 and 10; the mutation balance ν/(μ + ν); trajectories seeded; 1 and 8 threads alike |
| sweep | The probability of fixing from one copy against Kimura's u(1/(2N)) = (1 − e^{−2s})/(1 − e^{−4Ns}) (300 sweeps); one generation's law, Binomial(2N, ψ(p)) with ψ(p) = p(1 + s)/(1 + sp), by its residuals where conditioning on fixation does not reach (mean 0, variance 1 within 1%); a sweep that fixes shorter than the deterministic logistic from one copy; hitchhiking: π within 0.01 of the selected site at 14% of what drift alone leaves of the same standing variation, (1 − 1/(2N))^G, and within 10% of it far away; a partial sweep stopped at 0.7, its selected site segregating in the sample at the population's frequency; replicates reproduced from their seeds |
| structure | Balding–Nichols variance F p(1 − p); the kinship matrix against a naive M Mᵀ/m, 1 and 8 threads alike; eigenpairs, Patterson's scale; K − 1 axes standing out with Tracy–Widom statistics above 3, none without structure, and the null statistic distributed about as TW₁; the populations separated by two components; ancestry recovered (mean absolute error about 0.012), likelihood above the truth's, SQUAREM faster than plain EM, 1 and 8 threads alike |
| vcf | The VCF 4.2 specification's example (§1.1): its two biallelic sites with their dosages, positions, names, unphased, three records skipped; a phased VCF written from known genotypes on two chromosomes read back to the same dosages, haplotypes, positions, chromosomes and names, a multi-allelic record skipped; unsorted records sorted; haploid calls, unknown alleles and missing columns made missing and counted |

## 5. Performance

Measured on the development machine (i7-11370H, WSL, `-O3 -march=native`):

- a genealogy of 1000 samples with θ = 100 and its haplotype matrix: about 1 ms;
- 100 samples with ρ = 100 and θ = 100 (about 1100 recombinations, 450 marginal trees): about 2 ms;
- r² of 300,000 pairs among 400 haplotypes, packed, on 8 threads: about 3 ms (some 10⁸ pairs a second);
- 1000 Wright–Fisher replicates of 300 generations in 5 demes: about 80 ms on 8 threads;
- the kinship matrix and principal components of 120 samples × 5000 sites: about 30 ms.

## 6. Open Questions

1. **VCF.** Genotypes are built from raw arrays. The genomics domain's `variant` module is to read VCF into them: dosages site by site from GT, phased alleles when every GT of a record is phased, POS − 1, the chromosome's index, biallelic records only (multiallelic records split or dropped, as the caller chooses).
2. **Normalised H and other tests.** Fay and Wu's H is unnormalised, as they defined it; Zeng et al.'s (2006) normalisation, Fu and Li's tests and the Hudson–Kaplan bound R_m would come next.
3. **Tracy–Widom P values** need the TW₁ distribution function, which belongs in calx's special functions; and a hypergeometric law, used here for projecting spectra, belongs in calx too.
4. **Structured coalescent.** Migration among demes in the coalescent (ms's -I) is not simulated yet, nor the structured coalescent of a sweep (Kaplan, Hudson and Langley 1989); Wright–Fisher carries migration and selection forward in time, and `sweep.hpp` a sweep with its chromosomes, without new mutation during it (drift lowers the neutral diversity by (1 − 1/(2N))^G meanwhile). Soft sweeps from standing variation would come next.
5. **An accelerator.** SQUAREM is written here for the EM of ancestry; it is a general accelerator of fixed-point iterations and may move to calx's solvers.

## 7. References

- Alexander, D. H., Novembre, J. and Lange, K. (2009), Fast model-based estimation of ancestry in unrelated individuals, *Genome Research* 19, 1655–1664.
- Balding, D. J. and Nichols, R. A. (1995), A method for quantifying differentiation between populations at multi-allelic loci, *Genetica* 96, 3–12.
- Bhatia, G., Patterson, N., Sankararaman, S. and Price, A. L. (2013), Estimating and interpreting F_ST: the impact of rare variants, *Genome Research* 23, 1514–1521.
- Excoffier, L. and Slatkin, M. (1995), Maximum-likelihood estimation of molecular haplotype frequencies in a diploid population, *Molecular Biology and Evolution* 12, 921–927.
- Fay, J. C. and Wu, C.-I. (2000), Hitchhiking under positive Darwinian selection, *Genetics* 155, 1405–1413.
- Fu, Y.-X. (1995), Statistical properties of segregating sites, *Theoretical Population Biology* 48, 172–197.
- Griffiths, R. C. (1981), Neutral two-locus multiple allele models with recombination, *Theoretical Population Biology* 19, 169–186.
- Griffiths, R. C. and Tavaré, S. (1994), Sampling theory for neutral alleles in a varying environment, *Phil. Trans. R. Soc. B* 344, 403–410.
- Hill, W. G. (1974), Estimation of linkage disequilibrium in randomly mating populations, *Heredity* 33, 229–239.
- Hill, W. G. and Robertson, A. (1968), Linkage disequilibrium in finite populations, *Theoretical and Applied Genetics* 38, 226–231.
- Hudson, R. R. (1983), Properties of a neutral allele model with intragenic recombination, *Theoretical Population Biology* 23, 183–201.
- Hudson, R. R. (1990), Gene genealogies and the coalescent process, *Oxford Surveys in Evolutionary Biology* 7, 1–44.
- Hudson, R. R. (2002), Generating samples under a Wright–Fisher neutral model of genetic variation, *Bioinformatics* 18, 337–338.
- Hudson, R. R., Slatkin, M. and Maddison, W. P. (1992), Estimation of levels of gene flow from DNA sequence data, *Genetics* 132, 583–589.
- Kaplan, N. L., Hudson, R. R. and Langley, C. H. (1989), The "hitchhiking effect" revisited, *Genetics* 123, 887–899.
- Kelleher, J., Etheridge, A. M. and McVean, G. (2016), Efficient coalescent simulation and genealogical analysis for large sample sizes, *PLoS Computational Biology* 12, e1004842.
- Kimura, M. (1962), On the probability of fixation of mutant genes in a population, *Genetics* 47, 713–719.
- Kimura, M. and Ohta, T. (1969), The average number of generations until fixation of a mutant gene in a finite population, *Genetics* 61, 763–771.
- Kingman, J. F. C. (1982), The coalescent, *Stochastic Processes and their Applications* 13, 235–248.
- Lewontin, R. C. (1964), The interaction of selection and linkage. I. General considerations; heterotic models, *Genetics* 49, 49–67.
- Maynard Smith, J. and Haigh, J. (1974), The hitch-hiking effect of a favourable gene, *Genetical Research* 23, 23–35.
- Nei, M. and Li, W.-H. (1979), Mathematical model for studying genetic variation in terms of restriction endonucleases, *PNAS* 76, 5269–5273.
- Patterson, N., Price, A. L. and Reich, D. (2006), Population structure and eigenanalysis, *PLoS Genetics* 2, e190.
- Slatkin, M. and Hudson, R. R. (1991), Pairwise comparisons of mitochondrial DNA sequences in stable and exponentially growing populations, *Genetics* 129, 555–562.
- Tajima, F. (1983), Evolutionary relationship of DNA sequences in finite populations, *Genetics* 105, 437–460.
- Tajima, F. (1989), Statistical method for testing the neutral mutation hypothesis by DNA polymorphism, *Genetics* 123, 585–595.
- Tang, H., Peng, J., Wang, P. and Risch, N. J. (2005), Estimation of individual admixture: analytical and study design considerations, *Genetic Epidemiology* 28, 289–301.
- Varadhan, R. and Roland, C. (2008), Simple and globally convergent methods for accelerating the convergence of any EM algorithm, *Scandinavian Journal of Statistics* 35, 335–353.
- Watterson, G. A. (1975), On the number of segregating sites in genetical models without recombination, *Theoretical Population Biology* 7, 256–276.
- Weir, B. S. and Cockerham, C. C. (1984), Estimating F-statistics for the analysis of population structure, *Evolution* 38, 1358–1370.
- Wigginton, J. E., Cutler, D. J. and Abecasis, G. R. (2005), A note on exact tests of Hardy–Weinberg equilibrium, *American Journal of Human Genetics* 76, 887–893.
- Wright, S. (1931), Evolution in Mendelian populations, *Genetics* 16, 97–159.
