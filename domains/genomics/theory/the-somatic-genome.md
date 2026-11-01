# The Somatic Genome

*Status: draft. Serves learning, use and research.* A tumour is read as a mixture: normal cells and the cells of a tree of clones, each clone carrying the mutations of its ancestors, all of them carrying one allele-specific copy-number profile. Everything the module [`somatic`](../somatic/) computes is an inversion of that mixture from reads: which bases changed (MuTect's two tests), how much DNA each allele has (ASCAT's logR and BAF), in what fraction of cells each mutation lies (the cancer cell fraction), which mutations travel together (clones and their tree), which processes made them (signatures), and whether the microsatellites slipped (MSIsensor). This note states the model, proves the relations the code relies on, and says where the inversion is not unique. The simulator that makes the truth is [`simulate/tumour.hpp`](../simulate/include/genomics/simulate/tumour.hpp); the reads and their pileups are those of [mapping and variants](mapping-and-variants.md).

## 1. The Tumour as a Mixture of Cell Populations

**Definition 1.1 (sample, clones, cancer cell fraction).** A *tumour sample* is a finite multiset of cells, a fraction ρ ∈ (0, 1] of them tumour cells (the *purity*) and the rest normal. The tumour cells descend from one cell, the *trunk*; a *clone* is a set of tumour cells consisting of one cell's descendants, and the clones a sample is described by form a rooted tree T under inclusion, the trunk at the root. The *cancer cell fraction* of a clone k is φ_k, the fraction of tumour cells in k (k's own cells and its descendants'). A mutation of clone k is carried by exactly the cells of k.

**Proposition 1.2 (the sum rule).** For every clone k, Σ_{c child of k} φ_c ≤ φ_k, and φ_trunk = 1.

*Proof.* The children of k are clones contained in k; two of them are disjoint, since a cell has one ancestor of each generation, so the cell sets of distinct children (each the descendants of a different cell of the same generation below k) do not meet. Their measures add up to at most k's. The trunk contains every tumour cell. ∎

**Corollary 1.3 (the pigeonhole principle; Nik-Zainal et al. 2012).** If φ_a + φ_b > 1 for two clones, one contains the other.

*Proof.* Two clones of a tree are nested or disjoint. Disjoint clones have φ_a + φ_b ≤ 1. ∎

The fraction of tumour cells that are of clone k and of none of its descendants is φ_k − Σ_children φ_c ≥ 0, and these fractions over all clones sum to φ_trunk = 1; the simulator draws a tumour molecule's cell by them.

*Remark.* One sample gives the CCFs, not the tree: many trees satisfy Proposition 1.2 for one vector φ. With K clones sorted by decreasing φ, each non-trunk clone's parent is one of the clones before it, so there are (K − 1)! candidate trees, which `clone_trees` enumerates and filters; the count of survivors measures how underdetermined the tree is. Corollary 1.3 is the only nesting a single sample can force. Several samples of one tumour (regions, times) give one vector φ per sample and the sum rule in each; their intersection is what pins a tree down (El-Kebir et al. 2015).

## 2. Allele Fractions and the Cancer Cell Fraction

**Definition 2.1.** At a site of the reference, a tumour cell has n_A copies of the first parental haplotype and n_B of the second, N = n_A + n_B; a normal cell has 1 + 1. A mutation on haplotype h has *multiplicity* m ∈ {1, …, n_h}: the number of copies of h carrying it in a cell that carries it.

**Proposition 2.2 (the expected allele fraction).** If reads are drawn uniformly from the molecules of the sample, a mutation of CCF φ and multiplicity m at a site of tumour copy number N is read by a fraction

f = ρ φ m / (2(1 − ρ) + ρ N)

of the reads at that site.

*Proof.* Per cell of the sample, the site contributes 2 molecules from a normal cell (weight 1 − ρ) and N from a tumour cell (weight ρ): 2(1 − ρ) + ρN in all. The mutated molecules are m per tumour cell carrying the mutation, a fraction φ of the tumour cells: ρφm. A uniformly drawn molecule is mutated with the ratio of the two. ∎

The simulator's reads are drawn exactly so (a start kept with probability proportional to the DNA there, then a molecule by weight), and its check holds the reads carrying the substitutions to Σ f · depth within four standard deviations.

**Corollary 2.3 (McGranahan et al. 2015).** Given an observed allele fraction v, purity ρ and copy number N, the *mutation copy number* n_mut = (v/ρ)(ρN + 2(1 − ρ)) estimates φm, and φ̂ = n_mut/m for the multiplicity m.

*Remark (what is not identifiable).* Proposition 2.2 depends on φ and m only through their product. A mutation at two copies in half the cells and one at one copy in all the cells are read alike. The convention m = max(1, round(n_mut)), held to the major allele's copies, assumes that a mutation present at more than one copy is clonal; early mutations in gained segments (made before the gain) are the reason for multiplicity above one, and the simulator makes them exactly so (early: m = n_h; late or subclonal: m = 1). The checks measure how often the convention recovers the truth: 99.8% of the trunk's multiplicities on the simulated tumour.

For one mutation with a alternate reads of d, the code places a uniform prior on φ ∈ {0.01, …, 1} and computes the posterior ∝ Bin(a; d, c φ) with c = ρm/(2(1 − ρ) + ρN); the mutation is *clonal* when the upper end of the central 95% interval reaches 1 (McGranahan et al.'s rule, deliberately generous to clonality).

## 3. Copy Number from logR and BAF

**Proposition 3.1 (ASCAT's eqs. 1–2; Van Loo et al. 2010).** Let ψ be the tumour's mean copy number over the genome, and let the logR of a site be r = γ log₂(D/D̄), D the sample's DNA there and D̄ its mean. With γ = 1 (sequencing),

r = log₂ [(2(1 − ρ) + ρN) / (2(1 − ρ) + ρψ)],   b = (1 − ρ + ρ n_B) / (2(1 − ρ) + ρN),

b the expected fraction of reads carrying the B allele at a germline heterozygous site (normal 1 + 1, the B allele on the haplotype with n_B tumour copies).

*Proof.* D = 2(1 − ρ) + ρN as in Proposition 2.2, and its mean over the genome is 2(1 − ρ) + ρψ by linearity. The B allele is on one normal copy (weight 1 − ρ) and n_B tumour copies (weight ρ). ∎

The code's logR is the tumour's depth over the normal's in bins, each relative to its genome-wide mean, log₂: the normal's depth is proportional to 2 everywhere, so the ratio is D/D̄, and ψ is the length-weighted mean copy number.

**Proposition 3.2 (eqs. 3–4).** With T = 2(1 − ρ) + ρψ and R = 2^{r/γ} T,

n_A = (ρ − 1 − (b − 1) R)/ρ,   n_B = (ρ − 1 + bR)/ρ.

*Proof.* R = 2(1 − ρ) + ρN by eq. 1, and bR = 1 − ρ + ρn_B by eq. 2, which gives n_B; then ρn_A = R − 2(1 − ρ) − ρn_B = (1 − b)R + ρ − 1. ∎

The checks hold the two directions to rounding over 2000 random (ρ, ψ, γ, n_A, n_B).

**Proposition 3.3 (the dilation).** Write t = (1 − ρ)/ρ. The pair (r, b) of every segment depends on (ρ, ψ, n_A, n_B) only through (t, ψ, n_A, n_B) and is invariant under

(t, ψ, n_A, n_B) ↦ (λt, λψ, λn_A, λn_B),   λ > 0.

*Proof.* Divide numerators and denominators by ρ: r = log₂[(2t + N)/(2t + ψ)] and b = (t + n_B)/(2t + N), both homogeneous of degree 0 in (t, ψ, n_A, n_B). ∎

So logR and BAF determine the tumour only up to a one-parameter family; what selects a member is *integrality*: the copy numbers of every segment must be whole. If (ρ, ψ, n) fits with integers, so does (ρ_λ, λψ, λn) for every λ that keeps λn integral, ρ_λ = 1/(1 + λt): λ = 2 is the *genome-doubled* reading of the same data, with ψ doubled and purity lowered. ASCAT's distance (below) is zero at both; the grid, which contains ρ but in general not ρ_λ, and the bounds on ploidy decide between them. This is the classical ambiguity of ASCAT and ABSOLUTE (Carter et al. 2012 resolve it with a prior over karyotypes); the check of a doubled genome accepts either the truth or its halving.

**Definition 3.4 (ASCAT's distance and goodness of fit).** For segments s with logR r_s, BAF b_s ≥ 1/2 and S_s sites, at (ψ, ρ):

d(ψ, ρ) = Σ_s |ν_s − max(round(ν_s), 0)|² S_s w_s,

ν_s the minor allele's copy number by Proposition 3.2 (the allele of smaller sum over the segments), w_s = 0.05 for balanced segments (b_s = 1/2) and 1 otherwise; goodness 100(1 − d/d_max) with d_max = Σ_s 0.25 S_s w_s, the distance when every segment sits half a copy from an integer (ascat.runAscat.R).

*Remark.* d measures the minor allele alone. A segment's minor copy number ν = ((1 − b) R − (1 − ρ))/ρ moves with ψ only through R = 2^{r} T, by the factor T, so a change of ψ is partly absorbed by a change of ρ along the valleys of d: ψ is held more loosely than ρ. On profiles with logR noise of 0.02 the grid's ψ was off by up to 0.2 (ρ by up to 0.04) while every segment's integers, and hence the ploidy they give, were right (check).

The minimum of d over a grid of ψ ∈ [1, 6] and ρ ∈ [0.1, 1.05] is sought among strict local minima in 7 × 7 neighbourhoods, filtered by plausibility (ploidy, purity at least 0.2, goodness above 80%, some allele at zero), in ASCAT's four passes (copy_number.hpp). Its surface over (ψ, ρ) is the "sunrise" picture. On noise-free profiles the code finds the truth with d = 0 and goodness 100%; on noisy ones, purity to 0.03 and ψ to 0.1, with every segment's integers right.

**Segmentation (ASPCF).** The sites of a contig are cut into segments of constant logR and constant mirrored BAF min(b, 1 − b) by minimising

Σ_segments [Σ (r_i − r̄)²/σ_r² + Σ (b̃_i − b̄̃)²/σ_b²] + γ_pcf (segments − 1)

over segmentations with segments of at least 6 sites (Nilsen et al. 2012; ASCAT's γ_pcf = 70). Each term is the cost of tyche's Gaussian mean model, and PELT with a minimum segment length is exact for it (tyche/time_series/changepoint.hpp, Killick et al. 2012, Theorem 3.1, with the departure stated there for ℓ_min > 1). A segment's BAF is 1/2 + mean |b_i − 1/2|, set to 1/2 when √(σ_b² + μ²) < 2σ_b: the mean of |b − 1/2| of a balanced segment is σ_b √(2/π), not 0, and the test keeps such segments balanced.

## 4. Somatic Calling: MuTect's Tests

**Definition 4.1 (the models).** At a site with reference base r, candidate m and reads b_1, …, b_d with error probabilities e_i, model M_f^m says a fraction f of the molecules carry m:

P(b | e, r, m, f) = f P(b | e, m) + (1 − f) P(b | e, r),   P(b | e, A) = 1 − e (b = A), e/3 (b ≠ A),

the reads independent (Cibulskis et al. 2013, Online Methods). M₀ is f = 0. LOD_T = log₁₀ L(M_{f̂}^m)/L(M₀) with f̂ the fraction of reads showing m; LOD_N = log₁₀ L(M₀)/L(M_{1/2}^m) in the normal.

**Proposition 4.2 (the thresholds as posterior odds).** Calling a site when its posterior odds of M_f^m against M₀ exceed δ_T, with prior P(m, f) = P(m)P(f), P(f) uniform and P(m) = μ/3 for a somatic rate μ per base, is calling when LOD_T ≥ log₁₀ δ_T − log₁₀(μ/3) (to first order in μ). With δ_T = 2 and μ = 3 · 10⁻⁶: θ_T = 6.3. Calling a site somatic rather than germline when the normal's odds exceed δ_N is LOD_N ≥ log₁₀ δ_N + log₁₀[P(germline)/P(somatic)]: with δ_N = 10, P(somatic) = 3 · 10⁻⁶ and P(germline) = 5 · 10⁻⁵ off dbSNP or 0.095 at a dbSNP site, θ_N = 2.2 or 5.5.

*Proof.* Posterior odds are likelihood ratio times prior odds; take log₁₀ and use P/(1 − P) ≈ P for small P. The numbers are the paper's. ∎

The check computes the three thresholds from their priors, and LOD_T and LOD_N of small pileups term by term.

**Proposition 4.3 (the counts that pass).** Fix one error probability e for all reads, n reads, k of them showing m and n − k the reference, and f = k/n. Then LOD(k, n) = n g(k/n) with

g(f) = f log₁₀[(e/3 + cf)/(e/3)] + (1 − f) log₁₀[(1 − e − cf)/(1 − e)],   c = 1 − 4e/3,

g strictly convex on [0, 1] with g(0) = 0. Hence for every θ > 0 the counts with LOD ≥ θ are exactly k ≥ k*, and MuTect's sensitivity is P(X ≥ k*) for X ~ Bin(n, f(1 − e) + (1 − f)e/3).

*Proof.* With p_m = e/3 + cf and p_r = 1 − e − cf (the probabilities of reading m and r at fraction f), g′(f) = log(p_m/(e/3)) − log(p_r/(1 − e)) + (e/3)(1/p_r − 1/p_m) up to the factor 1/ln 10, using f c/p_m = 1 − (e/3)/p_m and (1 − f)c/p_r = 1 − (e/3)/p_r (as 1 − e − c = e/3). Differentiating again, g″(f) ∝ c/p_m + c/p_r + (e/3)c/p_r² + (e/3)c/p_m² > 0. The sublevel set {g < θ/n} of a convex function is an interval, and it contains 0; its complement in [0, 1] is an interval [f*, 1]. ∎

*Remark.* g′(0) < 0: a single read of m among many is evidence *against* a variant, the reads of m expected from errors alone being about de/3. The plug-in f̂ = k/n is not the maximiser of the likelihood; it is MuTect's choice for speed, and Proposition 4.3 shows it costs nothing in the shape of the decision.

The callable territory uses it: 14 tumour reads detect f = 0.3 at Q30 with power 0.84 (13 give 0.797), and 8 reference reads at Q30 give LOD_N = 2.41 ≥ 2.2 (7 give 2.11); the check computes both and holds the power to the sum over every k.

The filters are MuTect's Table 1 (mutect.hpp). The strand filter uses Proposition 4.3 on each strand: a strand LOD below 2 is evidence of bias only where the strand had the power (≥ 0.9) to reach 2.

*Caveat.* The model weighs bases as aligned. Mutect2 and Strelka2 realign reads to candidate haplotypes (local assembly; Strelka2's tiered haplotype model, Kim et al. 2018), which matters most for indels and for substitutions near them; this module applies the same two tests to indel counts (Definition 4.1 with the read's indel status in place of its base and a fixed ε), which is an extension, not MuTect.

## 5. Clones as a Binomial Mixture

The mutations' reads (a_i, d_i) with scales c_i = ρm_i/(2(1 − ρ) + ρN_i) are a finite mixture

P(a_i) = Σ_k π_k Bin(a_i; d_i, c_i φ_k),

the clones' CCFs φ_k and their shares π_k of the mutations (DPClust and PyClone without copy-number uncertainty). It is fitted by EM (tyche/multivariate/binomial_mixture.hpp): the E-step's responsibilities, then the shares and, for each clone, the root of a strictly concave score, so that the log-likelihood never decreases (Dempster, Laird and Rubin 1977); K by BIC. The winner's curse biases the smallest clone upward: mutations of a subclone at low VAF are called when their reads happen to be many. On the simulated tumour (purity 0.75, 100×) the clones of CCF 1, 0.55 and 0.3 come back as 0.973, 0.531 and 0.323, the smallest biased upward as the winner's curse says, and BIC chooses three.

## 6. Catalogues and Signatures

**Proposition 6.1 (SBS-96 and DBS-78 as orbits).** Let C₂ act on directed substitutions by reverse complementation. It acts freely on the 192 directed trinucleotide substitutions, which therefore fall into 96 orbits, each with exactly one member whose mutated base is a pyrimidine. On the 144 directed doublet substitutions (both bases changed) it has 78 orbits.

*Proof.* A trinucleotide xyz is fixed only if y equals its complement, which no base does; so every orbit has two members, and complementation exchanges purines and pyrimidines at the middle base. For doublets, Burnside's lemma: the 144 directed doublets are 16 references times 9 alternates; a fixed one has a palindromic reference (AT, CG, GC, TA: 4) and a palindromic alternate a₀a₁ = rc(a₀a₁) with a₀ ≠ r₀ and a₁ ≠ r₁; for a palindromic reference r₀ r̄₀ that is a₀ ∈ {the three bases ≠ r₀}, a₁ = ā₀, and a₁ ≠ r̄₀ holds automatically: 3 each, 12 in all. Orbits: (144 + 12)/2 = 78. ∎

The check verifies both by enumeration, and the 78 representatives against SigProfilerMatrixGenerator's list (which member names an orbit is a convention).

**The model of signature analysis.** A cohort's catalogues V (96 × n) are Poisson counts with means WH, W's columns the signatures (distributions over the classes) and H the exposures. Then −log L(W, H) = D(V‖WH) + (terms in V), the generalised Kullback–Leibler divergence, which Lee and Seung's multiplicative updates do not increase; the simulator's cohorts are drawn from exactly this model (Poisson numbers per signature, multinomial classes).

**Proposition 6.2 (non-identifiability).** If V = WH with W, H ≥ 0, then V = (WQ)(Q⁻¹H) for every invertible Q with WQ ≥ 0 and Q⁻¹H ≥ 0; among them, every Q = PD with P a permutation and D a positive diagonal. Normalising W's columns to sum 1 removes D, not P, and not the other Q.

*Proof.* Direct. A non-monomial Q exists with both products non-negative whenever W's columns are not extreme rays of the cone they span tightly enough (Donoho and Stodden 2004; Huang, Sidiropoulos and Swami 2014 give the "sufficiently scattered" condition under which only PD remains). ∎

So extraction is judged up to permutation (the Hungarian assignment of extracted to reference signatures by cosine), and its reproducibility over resampled catalogues (tyche's stability, SigProfilerExtractor's silhouettes) is the evidence that the data determine the rank. On a cohort of 40 simulated tumours with three synthetic signatures, ranks 2 and 3 are stable (average stabilities 0.998 and 0.994) and ranks 4 and 5 are not (0.52 and 0.26, with clusters of negative silhouette), the rank chosen is 3, and the extracted signatures match their truths with cosines 0.999, 0.998 and 0.994 (check).

**Refitting and its bootstrap.** For a known set S, the exposures of one catalogue c are the non-negative least squares argmin_{e ≥ 0} ‖Se − c‖₂, characterised by the Karush–Kuhn–Tucker conditions (tyche's NNLS). Resampling the n mutations of c with replacement gives Multinomial(n, c/n) catalogues exactly; refitting each gives the sampling distribution of the exposures given the catalogue, from which the bootstrap interval comes. It does not measure the model's error (signatures missing from S, or wrong).

## 7. Microsatellite Instability

At a microsatellite of the reference (a unit of 1–5 bases repeated), each read spanning it with both flanks gives a length: the number of whole units between the flanks. The normal's and the tumour's lengths form a 2 × K table, tested by Pearson's χ² with K − 1 degrees of freedom (MSIsensor; Niu et al. 2014). A slip on one haplotype of a tumour of purity ρ moves a fraction ρ/2 of the tumour's reads by one unit, so the test's power grows with ρ and depth.

**Proposition 7.1 (MSIsensor's rule and Benjamini–Hochberg).** Let p₍₁₎ ≤ … ≤ p₍N₎ be the sites' p-values. MSIsensor calls site (k) somatic when p₍k₎N/k ≤ α; Benjamini and Hochberg's step-up rejects (1), …, (k̂) with k̂ = max{k : p₍k₎ ≤ kα/N}. The first set is contained in the second.

*Proof.* A site (k) with p₍k₎ ≤ kα/N has k ≤ k̂ by the definition of k̂, so it is among BH's rejections. ∎

The inclusion can be strict (p = 0.01, 0.04, 0.045 at α = 0.05: MSIsensor takes the first and the third, BH all three; checked), so MSIsensor's rule does not control the false discovery rate in BH's sense; it is followed because the score, the percentage of somatic sites among the tested ones with MSI-high at 3.5%, is MSIsensor's.

## 8. What Is Checked

| Claim | Check |
|---|---|
| Proposition 2.2 for simulated reads | `check_genomics_simulate_tumour`: reads carrying substitutions against Σ f · depth |
| Propositions 3.1–3.3 | `check_genomics_somatic_copy_number`: inversion to 10⁻¹², the dilation to 10⁻¹⁴, the fit on noise-free, noisy and doubled profiles |
| Propositions 4.2, 4.3 | `check_genomics_somatic_mutect`: the thresholds from their priors, the power against the sum over every count, the depths 14 and 8 |
| Proposition 6.1 | `check_genomics_contexts`: enumeration of the 192 and 144 directed substitutions, the orbit sizes, SigProfilerMatrixGenerator's names |
| Proposition 6.2, extraction | `check_genomics_somatic_signatures`: rank and cosine on a simulated cohort |
| Proposition 7.1 | `check_genomics_somatic_msi`: the example where the sets differ; χ² against scipy |
| The whole inversion | `check_genomics_somatic_tumour`: calls, purity, ploidy, integers, CCFs, tree, exposures, MSI against the simulator's truth |

## References

- Alexandrov, L. B. et al. (2013), Signatures of mutational processes in human cancer, *Nature* 500, 415–421; Deciphering signatures of mutational processes operative in human cancer, *Cell Reports* 3, 246–259.
- Alexandrov, L. B. et al. (2020), The repertoire of mutational signatures in human cancer, *Nature* 578, 94–101.
- Carter, S. L. et al. (2012), Absolute quantification of somatic DNA alterations in human cancer, *Nature Biotechnology* 30, 413–421.
- Cibulskis, K. et al. (2013), Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples, *Nature Biotechnology* 31, 213–219.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *JRSS B* 39, 1–38.
- Dentro, S. C., Wedge, D. C. and Van Loo, P. (2017), Principles of reconstructing the subclonal architecture of cancers, *Cold Spring Harbor Perspectives in Medicine* 7, a026625.
- Donoho, D. and Stodden, V. (2004), When does non-negative matrix factorization give a correct decomposition into parts?, *NIPS* 16.
- El-Kebir, M., Oesper, L., Acheson-Field, H. and Raphael, B. J. (2015), Reconstruction of clonal trees and tumor composition from multi-sample sequencing data, *Bioinformatics* 31, i62–i70.
- Gillis, S. and Roth, A. (2020), PyClone-VI: scalable inference of clonal population structures using whole genome data, *BMC Bioinformatics* 21, 571.
- Huang, K., Sidiropoulos, N. D. and Swami, A. (2014), Non-negative matrix factorization revisited: uniqueness and algorithm for symmetric decomposition, *IEEE Trans. Signal Process.* 62, 211–224.
- Islam, S. M. A. et al. (2022), Uncovering novel mutational signatures by de novo extraction with SigProfilerExtractor, *Cell Genomics* 2, 100179.
- Kim, S. et al. (2018), Strelka2: fast and accurate calling of germline and somatic variants, *Nature Methods* 15, 591–594.
- Killick, R., Fearnhead, P. and Eckley, I. A. (2012), Optimal detection of changepoints with a linear computational cost, *JASA* 107, 1590–1598.
- McGranahan, N. et al. (2015), Clonal status of actionable driver events and the timing of mutational processes in cancer evolution, *Science Translational Medicine* 7, 283ra54.
- Nik-Zainal, S. et al. (2012), The life history of 21 breast cancers, *Cell* 149, 994–1007.
- Nilsen, G. et al. (2012), Copynumber: efficient algorithms for single- and multi-track copy number segmentation, *BMC Genomics* 13, 591.
- Niu, B. et al. (2014), MSIsensor: microsatellite instability detection using paired tumor-normal sequence data, *Bioinformatics* 30, 1015–1016.
- Van Loo, P. et al. (2010), Allele-specific copy number analysis of tumors, *PNAS* 107, 16910–16915.
