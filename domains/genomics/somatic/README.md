# somatic: The Somatic Genome of a Tumour

*Status: draft (2026-11-01).* A tumour read against its matched normal: the somatic substitutions and indels it carries (MuTect's tumour and normal tests and its filters), its purity, ploidy and allele-specific copy number (logR and BAF segmented jointly, ASCAT's fit), the cancer cell fraction of each mutation and the clones they cluster into with the trees the sum rule allows, its mutational catalogue and the signatures that explain it (de novo by NMF over a cohort, refitted by NNLS, bootstrapped), its microsatellite instability (MSIsensor), and its mutational burden. Every method is checked against the simulator's truth (`simulate/tumour.hpp`) and, where a reference exists, against it. Serves use, learning and research (ADR 0064); the theory is [the-somatic-genome.md](../theory/the-somatic-genome.md).

## 1. Concepts

- **A tumour sample is a mixture.** Normal cells (1 − ρ) and the cells of a tree of clones, each carrying its ancestors' mutations and all one copy-number profile. A mutation of CCF φ and multiplicity m at tumour copy number N is read by a fraction ρφm/(2(1 − ρ) + ρN) of the reads (theory, Proposition 2.2); every method here inverts some part of that relation.
- **Two likelihood ratios decide a somatic call.** The tumour's evidence for an allele at fraction f against errors, and the normal's evidence against a heterozygous germline variant, each against a threshold derived from priors (6.3 and 2.2; Proposition 4.2), then filters for the artefacts that imitate mutations.
- **Copy number is identified up to a dilation.** logR and BAF determine (ρ, ψ, n) only up to (t, ψ, n) ↦ (λt, λψ, λn), t = (1 − ρ)/ρ (Proposition 3.3); integrality chooses, which is what ASCAT's distance measures.
- **One sample does not determine the tree.** The sum rule leaves several trees; the module enumerates them and says how many.
- **Signatures are a property of a cohort.** NMF of one catalogue is meaningless; extraction needs many tumours, refitting needs a known set.

## 2. Content

| Header | Content |
|---|---|
| [`coverage.hpp`](include/genomics/somatic/coverage.hpp) | `MappedReads` (a sample's primary mappings sorted, `overlapping(contig, b, e)`); `depth_bins(contigs, sample, width)`; `Site`, `allele_counts(pool, contigs, sites, sample)` (the variant module's pileup at each site) |
| [`mutect.hpp`](include/genomics/somatic/mutect.hpp) | `MutectOptions` (θ_T 6.3, θ_N 2.2 and 5.5 at known germline sites, the filters' parameters, indels, the callable depths 14 and 8); `log10_likelihood`, `tumour_lod`, `normal_lod` (Definition 4.1); the counting forms for indels; `lod_at`, `power` (MuTect's sensitivity); `median_mad`; `query_offset(mapping, p)`; `SomaticCall` (LODs, depths, strands, strand LODs and powers, position-in-read statistics, filters, `vaf()`); `call_somatic_region`, `call_somatic(pool, contigs, tumour, normal, options)` → `SomaticCalls{calls, callable, sites_weighed}`; `as_call` (a `variant::Call` for evaluation and VCF). Filters: proximal_gap, poor_mapping, triallelic_site, strand_bias, clustered_position, observed_in_control |
| [`copy_number.hpp`](include/genomics/somatic/copy_number.hpp) | ASCAT's eqs. 1–4 (`ascat_logr`, `ascat_baf`, `ascat_na`, `ascat_nb`), `round_half_even`; `heterozygous_sites` (the normal's 0/1 SNVs read in both samples), `logr_bins`, `site_logr`; `aspcf` (joint PELT of logR and mirrored BAF on tyche's SumCost, penalty 70, 6 sites); `ascat_distance`, `ascat(segments)` → `AscatFit` (the grid, the "sunrise" distances, the four passes, the optimum, goodness of fit, integer copy numbers written into the segments); `copy_number(...)` → `CopyNumberProfile{bins, sites, segments, fit, ploidy}` |
| [`clonality.hpp`](include/genomics/somatic/clonality.hpp) | `cancer_cell_fractions(calls, segments, purity)` → `MutationCcf` (n_mut, multiplicity, CCF, the 95% grid interval, clonal by McGranahan's rule); `cluster_mutations` (tyche's binomial mixture, K by BIC); `clone_trees(ccf, tolerance)` (every tree of (K − 1)! satisfying the sum rule, pigeonhole pairs, the tree of greatest least slack); `descends`; `mutational_burden` |
| [`catalogue.hpp`](include/genomics/somatic/catalogue.hpp) | `catalogue(calls, contigs)` → `Catalogue{sbs96, dbs78, …}` (adjacent substitutions as doublets, runs of three in neither); `read_profiles(text)` → `ProfileTable` (COSMIC's or SigProfiler's tab-separated matrices by class names, in any order, problems with their lines) |
| [`signatures.hpp`](include/genomics/somatic/signatures.hpp) | `extract_signatures(catalogues, rank_min, rank_max, stability options)` (tyche's NMF stability; the largest stable rank); `match_signatures` (Hungarian assignment by cosine, calx); `refit(catalogue, signatures)` (tyche's NNLS, reconstruction cosine); `bootstrap_exposures` (multinomial resampling, calx; means, sd, 2.5–97.5% quantiles) |
| [`msi.hpp`](include/genomics/somatic/msi.hpp) | `find_microsatellites(contigs, options)` (maximal runs of primitive units of 1–5 bases, MSIsensor's bounds, flanks of 5); `copies_in_read`; `read_lengths`; `chi_square_two` (Pearson's χ², MSIsensor's X2BetweenTwo); `msisensor_fdr` (its rank rule); `msi(pool, sites, tumour, normal)` → `MsiResult{sites, tested, somatic, score, high}` |

The simulator is `simulate/tumour.hpp`; the classes of substitutions are `sequence/contexts.hpp`.

## 3. What Is Not Implemented, and Why

- **COSMIC's signatures are not fetched.** COSMIC's terms (https://www.cosmickb.org/terms/, read 2026-11-01) allow free use only for internal academic research by not-for-profit organisations, forbid redistribution and display on free-access platforms, and cover the mutational signatures; the data library therefore does not register them. The simulator uses five synthetic, mechanism-motivated profiles (`synthetic_signatures`, labelled synthetic everywhere), and `read_profiles` reads COSMIC's matrix format from a file a licensed user supplies.
- **ID-83 indel catalogues** (Alexandrov et al. 2020): their 83 classes are in the fixture from SigProfilerMatrixGenerator, but its classifier (repeat units and microhomology read from the reference) is embedded in that program's file-reading code and cannot be run on its own; a classifier written from the published definitions and checked against it is planned.
- **Mutect2 and Strelka2.** Local reassembly of haplotypes and realignment of reads (Mutect2's), and Strelka2's tiered haplotype model, are not written; the calls weigh bases as aligned, and indels are tested by MuTect's model on counts (an extension). MuTect's estimate of cross-sample contamination (ContEst) and its panel of normals are not written.
- **Subclonal copy number** (Battenberg), **ABSOLUTE's karyotype prior**, and copy-number states in the clustering (PyClone's) are not written: copy number is clonal here, and so it is in the simulator.
- **Variational Dirichlet-process clustering** (PyClone-VI) is not written; the binomial mixture with BIC is.
- **SigProfilerExtractor's choice of rank** by a test of reconstruction errors among stable ranks is not written; the largest stable rank is chosen. Its GMM normalisation is not written (tyche's NMF has log2 and 100× normalisations).
- **MSIsensor-pro's** tumour-only model of polymerase slippage, and MSIsensor's coverage normalisation (off by default there), are not written.

## 4. Dependencies

`genomics::variant` (pileups, the normal's caller, normalisation, evaluation), `genomics::mapping`, `genomics::sequence`; tyche (`time_series`: PELT and its costs; `regression`: NNLS; `multivariate`: NMF and its stability, the binomial mixture; `descriptive`: quantiles); calx (`special_functions`: log Γ, the χ² and binomial tails; `solvers`: Brent; `random`: streams, multinomial; `optimisation`: the Hungarian method; `scheduler`). The checks also use `genomics::simulate` and `calx::verification`.

## 5. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time (`check_genomics_somatic_<topic>`).

| Check | What it verifies |
|---|---|
| `mutect` | LOD_T and LOD_N of a small pileup against the model written term by term; θ_T, θ_N and θ_N at dbSNP from their priors (6.30, 2.22, 5.50); one read of m never passes; MuTect's sensitivity against the sum over every count (to 4 · 10⁻¹⁴); 14 tumour reads detect f = 0.3 at Q30 with power 0.839 and 13 with 0.797; 8 reference reads give LOD_N 2.407 and 7 give 2.106; the indel forms on counts; median and MAD; offsets through a CIGAR |
| `copy_number` | Eqs. 3–4 invert eqs. 1–2 (2000 draws, to 6 · 10⁻¹⁵); R's rounding; the dilation of Proposition 3.3 (to 2 · 10⁻¹⁶); the fit of a noise-free profile exact (ρ 0.65, ψ 2.40, goodness 100%, every integer); twenty noisy profiles at purities 0.3 to 0.87: every integer right in all twenty, purity within 0.05 (measured 0.04), the grid's ψ within 0.25 (measured 0.2: the distance weighs the minor allele only); a doubled genome fitted at ψ 4.80; ASPCF on 500 noisy sites of six states, every breakpoint within two sites, then ASCAT within 0.02 in purity and 0.01 in ψ with every integer |
| `clonality` | n_mut, multiplicity and CCF by hand; the posterior interval by brute force; a subclonal mutation; two clones recovered from 400 mutations; 6 of 6 trees examined for four clones, 5 satisfying the sum rule as counted directly; pigeonhole nesting; the trunk put above clones that fall short of 1; the burden |
| `catalogue` | SBS-96 and DBS-78 of a small reference against a hand count (purine-centred calls read on the other strand, adjacent calls as a doublet, a run of three in neither, filtered calls left out); tables by class names in any order; their problems; a DBS-78 table |
| `signatures` | NNLS refit of a noise-free mixture exact (to 10⁻¹²); a cohort of 40 simulated tumours of three synthetic signatures: ranks 2–3 stable (0.998, 0.994), 4–5 not (0.52, 0.26), rank 3 chosen; extracted signatures against the truth by cosine (0.999, 0.998, 0.994; limit 0.97), one-to-one matching; refitted exposures of ten tumours within 14% (limit 15%); bootstrap 95% intervals covering the true exposure 19 of 22 times |
| `msi` | The scan against a brute-force count of maximal primitive runs; planted sites found and short ones not; copies between flanks; χ² against scipy on 40 tables (statistics to 2 · 10⁻¹⁶, p-values to 2 · 10⁻¹⁴, relative); MSIsensor's rank rule against a hand example where it differs from Benjamini–Hochberg |
| `tumour` | End to end on a 500 kb reference with 215 planted microsatellites, a normal diploid, a tumour of purity 0.75 with clones of CCF 1, 0.55 and 0.3 (1756 somatic variants, three synthetic signatures), six copy-number states (ploidy 2.68) and 56 slipped microsatellites, 100× and 40× of 150-base reads, mapped: somatic SNVs precision ≥ 0.97 (measured 1.000) and recall 0.954 overall, 0.985 at expected VAF ≥ 0.15 (limit 0.95), 0.877 below; indels precision 1.000; the same calls on one thread and eight; purity 0.75 within 0.05 (exact on the grid), ploidy within 0.1 (−0.003), integer copy numbers right over 99.7% of the bases (limit 95%); three clones of CCF 0.973, 0.531, 0.323 (limits 0.05, 0.07, 0.08) and the true chain among the 2 trees the sum rule allows; clonal or not right for 86% (limit 85%); the trunk's multiplicities 99.8%; the catalogue's exposures to the three signatures within 4.3% (limit 12%), reconstruction cosine 0.993; MSI: 207 sites tested, precision 1.000, recall 0.87 (limit 0.8), MSI-high (score 22.7%); the callable territory 99.98% |

## 6. Performance

Measured on the development machine (i7-11370H, WSL, Release, `-march=native`), 2026-11-01, eight threads, with other builds running: the 500 kb tumour at 100× and its normal at 40× (333 333 and 133 333 reads) are simulated and mapped in 8.4 s; MuTect's tests weigh 24 732 candidate sites and call in 0.6 s; heterozygous sites, logR, ASPCF and ASCAT's grid (101 × 96 points) take 0.1 s; the binomial mixture for K = 1 to 6 and the MSI test take well under a second.

## 7. References

- Cibulskis, K. et al. (2013), Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples, *Nature Biotechnology* 31, 213–219.
- Kim, S. et al. (2018), Strelka2: fast and accurate calling of germline and somatic variants, *Nature Methods* 15, 591–594.
- Van Loo, P. et al. (2010), Allele-specific copy number analysis of tumors, *PNAS* 107, 16910–16915; the R package ASCAT (VanLoo-lab/ascat), ascat.runAscat.R and ascat.aspcf.R, read 2026-11-01.
- Nilsen, G. et al. (2012), Copynumber: efficient algorithms for single- and multi-track copy number segmentation, *BMC Genomics* 13, 591.
- Carter, S. L. et al. (2012), Absolute quantification of somatic DNA alterations in human cancer, *Nature Biotechnology* 30, 413–421.
- McGranahan, N. et al. (2015), Clonal status of actionable driver events and the timing of mutational processes in cancer evolution, *Science Translational Medicine* 7, 283ra54.
- Dentro, S. C., Wedge, D. C. and Van Loo, P. (2017), Principles of reconstructing the subclonal architecture of cancers, *Cold Spring Harbor Perspectives in Medicine* 7, a026625.
- Nik-Zainal, S. et al. (2012), The life history of 21 breast cancers, *Cell* 149, 994–1007.
- El-Kebir, M. et al. (2015), Reconstruction of clonal trees and tumor composition from multi-sample sequencing data, *Bioinformatics* 31, i62–i70.
- Gillis, S. and Roth, A. (2020), PyClone-VI, *BMC Bioinformatics* 21, 571.
- Alexandrov, L. B. et al. (2013), *Nature* 500, 415–421 and *Cell Reports* 3, 246–259; (2020), *Nature* 578, 94–101.
- Islam, S. M. A. et al. (2022), SigProfilerExtractor, *Cell Genomics* 2, 100179.
- Bergstrom, E. N. et al. (2019), SigProfilerMatrixGenerator, *BMC Genomics* 20, 685.
- Niu, B. et al. (2014), MSIsensor, *Bioinformatics* 30, 1015–1016; the program ding-lab/msisensor (chi.cpp, homo.cpp, sample.cpp), read 2026-11-01.
- Jia, P. et al. (2020), MSIsensor-pro, *Genomics, Proteomics & Bioinformatics* 18, 65–71.
- Chalmers, Z. R. et al. (2017), Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden, *Genome Medicine* 9, 34.

## 8. Open Questions

1. **Haplotype-aware calling.** Realigning reads to candidate haplotypes (Mutect2's local assembly, or the pair HMM of `align/pairhmm.hpp`) would weigh an indel and the substitutions near it together, and phase adjacent substitutions before they are counted as doublets.
2. **Subclonal copy number and copy-number-aware clustering.** A segment present in a fraction of the tumour cells fits half-integers badly under ASCAT; Battenberg fits two states; PyClone sums over copy-number states in each mutation's likelihood.
3. **Several samples of one tumour** would make the sum rule decide the tree (theory §1, Remark).
4. **Where the binomial mixture belongs.** It is tyche's (`multivariate`), as a mixture with known scales; a Dirichlet-process alternative belongs beside it.
