# association

*Status: draft. A module of the evolution domain (ADR 0064).* Genotype to trait and genotype to expression: association scans and genomic control, linkage disequilibrium and clumping, cis-eQTL mapping with permutations and q-values, fine-mapping (ABF, SuSiE, SuSiE-RSS), colocalisation, expression predicted from genotypes and its association with a trait (PrediXcan, S-PrediXcan, the TWAS), SNP heritability (GREML, Haseman–Elston, the liability scale) and polygenic scores (C+T, LDpred-inf); and a simulator of cohorts with truth to measure them all. Serves use (tested methods on real-shaped data), learning (each held to its definition and to the simulated truth) and research. The theory is [theory/association.md](../theory/association.md).

## 1. Concepts

- **Truth by simulation, clearly labelled.** `cohort.hpp` makes cohorts whose linkage is a population's: founders from the coalescent with recombination, individuals as Li–Stephens mosaics of their population's founders (as HAPGEN2 makes them), so that principal components find the populations and LD decays as it does. Traits and expression are put on top with known causal variants, heritability, cis-heritability, mediation and confounding. Every `Study` says `simulated = true`.
- **Statistics are tyche's; genetics is here.** SuSiE, REML, Haseman–Elston, q-values, the beta approximation of permutation minima and the elastic net name no object in nature and live in tyche (its decisions 0006 and 0007). This module gives them variants, genes, cis windows, LD from genotypes and credible sets of variants.
- **One projection a study.** Covariates are projected out once (Frisch–Waugh–Lovell), so a scan costs one inner product a site and equals the full regression exactly.
- **Summary statistics suffice** where they should: SuSiE-RSS from z-scores and LD, S-PrediXcan from z-scores and LD, coloc from estimates and standard errors; the checks hold each to its individual-level counterpart.

## 2. Content

| Header | Content |
|---|---|
| `cohort.hpp` | `CohortSpec` (individuals, regions, θ, ρ, populations, founders, admixed, switches, min_maf), `StudySpec` (h², causal, mediated, mediating, confounding, covariates' effects, prevalence; genes, eQTL architecture, cis window, hidden factors, the expression panel's size), `Gene`, `Study` (phased genotypes, populations, covariates, genes, expression and its genetic part, factors, the trait, liability and cases, the genetic value, the causal variants and effects, the realised h², the panel and the sample), `simulate_study(spec, seed, pool)` |
| `scan.hpp` | `Residualizer` (the projection off covariates), `sample_dosages`; `Scan`; `linear_scan` (vectorised, exact), `logistic_scan` (the score test, Wald refits of the strongest), `genomic_control`, `qq_points` (with the beta band) |
| `ld.hpp` | `standardised_dosages`, `ld_correlation`, `sites_between`, `ld_scores` (bias-corrected), `ClumpOptions`, `Clump`, `clump` (PLINK's), `Components`, `principal_components` (through the m′ × m′ covariance, equal to population::pca) |
| `eqtl.hpp` | `GeneLocus`, `CisOptions`, `CisGene`, `CisMapping`, `cis_eqtl` (nominal pass, permutations, effective df, the beta approximation, q-values, eGenes and their nominal thresholds), `expression_components` |
| `finemap.hpp` | `log_abf`, `SingleCausal`, `finemap_abf`; `finemap_susie` (individual data, covariates projected out), `finemap_rss`; `FinemapTruth`, `against_truth` |
| `coloc.hpp` | `ColocPriors`, `Coloc`, `coloc_from_lbf`, `coloc_abf` (coloc.abf), `ColocPair`, `coloc_susie` |
| `predict.hpp` | `PredictOptions`, `GeneModel`, `train_expression_model` (the elastic net with nested cross-validation, PredictDB's criteria), `train_expression_models`, `predict_expression`, `Association`, `associate` (the TWAS of individual data), `s_predixcan` |
| `heritability.hpp` | `genetic_relationship` (the GRM, from population's kinship), `Heritability`, `snp_heritability` (GREML and Haseman–Elston), `liability_scale` |
| `prs.hpp` | `polygenic_score`, `score_r2`, `ThresholdScore`, `clump_and_threshold`, `LdpredInf`, `ldpred_inf` |

## 3. Dependencies

- **Within evolution:** `population` (genotypes, the coalescent, kinship, PCA for the checks, F_ST), `haplotypes` (the copying chain's jumps).
- **Foundations:** tyche `regression` (OLS, GLMs, SuSiE, variance components, the elastic net), `testing` (q-values, permutation minima), `probability`; calx `linalg`, `random`, `scheduler`, `cpu`; `verification` for the checks.

## 4. Checks

| Check | What it verifies |
|---|---|
| scan | The cohort: a function of its seed, the same on 1 and 8 threads; frequency filter; LD decaying with distance; F_ST of the founder groups; realised h², cis-h² and case fraction. Principal components through the m′ × m′ covariance equal population::pca's axes and eigenvalues (to 10⁻⁹) and separate the populations. The linear scan equals tyche's ols on [C, g] (β, se, p to 10⁻⁹); λ_GC near 1 under a null trait, raised by confounding without components and restored with them; causal variants found. The logistic score χ² equals its definition from tyche's null glm; Wald refits equal glm. The Q–Q band; LD scores against direct correlations; clumps in LD with their index |
| studies | ABF in closed form; coloc.abf against R's coloc 5.2.3 (PP.H0–H4 and SNP.PP.H4 to 10⁻¹², three scenarios) and H₃, H₄ against the enumeration of configurations; SuSiE's credible sets in associated regions of a simulated cohort covering causal variants, SuSiE-RSS on z-scores and LD finding the same sets; cis-eQTL mapping finding eGenes among genes with causal eQTLs with few false ones, leads causal or in LD; PrediXcan's nested R² within (0, cis-h²], predictions tracking true genetic expression, no models for genes without eQTLs; S-PrediXcan's z equal to the individual TWAS z; the mediating genes found; GREML and Haseman–Elston recovering h²; the liability scale; LDpred-inf equal to ridge regression, and on held-out samples LDpred-inf and C+T predicting below h², C+T ahead on the sparse simulated architecture |

The statistics they rest on are checked in tyche against susieR 0.14.2, R's regress, qvalue and SciPy (`check_tyche_susie`, `check_tyche_variance_components`, `check_tyche_multiple`, `check_tyche_permutation`) and against glmnet's and scikit-learn's fits (`check_tyche_penalized`).

## 5. Performance

The engine's op `association.study` (3000 individuals, 4 regions of 1 Mb, 30 genes, 300 permutations) runs the whole study in seconds on 8 threads; the costliest steps are the permutations of the eQTL pass, the nested cross-validation of PrediXcan and the O(n³) REML (on at most 1000 samples there).

## 6. Open Questions

1. **LD score regression** (Bulik-Sullivan et al. 2015): the LD scores are here; the weighted regression and its block jackknife (a general method, for tyche's planned `estimation`) are planned.
2. **Real panels.** A region of the 1000 Genomes Project's phase 3 (free to use) through the data library would replace the simulated founders for real LD; the methods take any `Genotypes`.
3. **Mixed-model scans** (BOLT-LMM, REGENIE's first step) for related samples, using tyche's variance components.
4. **Rank-based inverse normal transforms** of expression (GTEx's), a general transform for tyche's `descriptive`.
5. **Conditional analysis and COJO** from summary statistics; **coloc.susie** in the engine (written here, not yet drawn).

## 7. References

- Barbeira, A. N. et al. (2018), Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics, *Nature Communications* 9, 1825.
- Bulik-Sullivan, B. K. et al. (2015), LD score regression distinguishes confounding from polygenicity in genome-wide association studies, *Nature Genetics* 47, 291–295.
- Choi, S. W., Mak, T. S.-H. and O'Reilly, P. F. (2020), Tutorial: a guide to performing polygenic risk score analyses, *Nature Protocols* 15, 2759–2772.
- Devlin, B. and Roeder, K. (1999), Genomic control for association studies, *Biometrics* 55, 997–1004.
- Gamazon, E. R. et al. (2015), A gene-based association method for mapping traits using reference transcriptome data, *Nature Genetics* 47, 1091–1098.
- Giambartolomei, C. et al. (2014), Bayesian test for colocalisation between pairs of genetic association studies using summary statistics, *PLoS Genetics* 10, e1004383.
- GTEx Consortium (2017), Genetic effects on gene expression across human tissues, *Nature* 550, 204–213.
- Gusev, A. et al. (2016), Integrative approaches for large-scale transcriptome-wide association studies, *Nature Genetics* 48, 245–252.
- Haseman, J. K. and Elston, R. C. (1972), The investigation of linkage between a quantitative trait and a marker locus, *Behavior Genetics* 2, 3–19.
- Lee, S. H., Wray, N. R., Goddard, M. E. and Visscher, P. M. (2011), Estimating missing heritability for disease from genome-wide association studies, *AJHG* 88, 294–305.
- Ongen, H. et al. (2016), Fast and efficient QTL mapper for thousands of molecular phenotypes, *Bioinformatics* 32, 1479–1485.
- Price, A. L. et al. (2006), Principal components analysis corrects for stratification in genome-wide association studies, *Nature Genetics* 38, 904–909.
- Purcell, S. et al. (2007), PLINK: a tool set for whole-genome association and population-based linkage analyses, *AJHG* 81, 559–575.
- Shabalin, A. A. (2012), Matrix eQTL: ultra fast eQTL analysis via large matrix operations, *Bioinformatics* 28, 1353–1358.
- Storey, J. D. and Tibshirani, R. (2003), Statistical significance for genomewide studies, *PNAS* 100, 9440–9445.
- Su, Z., Marchini, J. and Donnelly, P. (2011), HAPGEN2: simulation of multiple disease SNPs, *Bioinformatics* 27, 2304–2305.
- Vilhjálmsson, B. J. et al. (2015), Modeling linkage disequilibrium increases accuracy of polygenic risk scores, *AJHG* 97, 576–592.
- Wakefield, J. (2009), Bayes factors for genome-wide association studies: comparison with P-values, *Genetic Epidemiology* 33, 79–86.
- Wallace, C. (2021), A more accurate method for colocalisation analysis allowing for multiple causal variants, *PLoS Genetics* 17, e1009440.
- Wang, G., Sarkar, A., Carbonetto, P. and Stephens, M. (2020), A simple new approach to variable selection in regression, with application to genetic fine mapping, *JRSS B* 82, 1273–1300.
- Yang, J., Lee, S. H., Goddard, M. E. and Visscher, P. M. (2011), GCTA: a tool for genome-wide complex trait analysis, *AJHG* 88, 76–82.
- Zou, Y., Carbonetto, P., Wang, G. and Stephens, M. (2022), Fine-mapping from summary data with the "Sum of Single Effects" model, *PLoS Genetics* 18, e1010299.
