# expression: Transcript Abundance and Differential Expression

*Status: draft (2026-10-24).* From RNA-seq reads to what changed: transcriptomes and simulated samples with their truth, pseudoalignment of reads to k-mer equivalence classes, transcript abundances by EM with effective lengths, TPM and the bootstrap, and differential expression by negative binomial models with dispersions shrunk toward a trend, Wald tests and the false discovery rate.

## 1. Purpose

An RNA-seq sample is a few million fragments of the transcripts of a tissue, each copied in proportion to its abundance. Two questions follow, and this module answers both:

- **How much of each transcript?** kallisto's insight (Bray, Pimentel, Melsted and Pachter 2016) is that the likelihood of the abundances needs only which transcripts each read is *compatible* with, not where it aligns. A read's compatible set is the intersection of the sets of transcripts holding its k-mers (its *equivalence class*), found by hashing k-mers alone; the likelihood of the abundances is then a product over the classes, and the EM algorithm finds its maximum, sharing each class's reads among its transcripts by their current abundances.
- **What changed between conditions?** Counts of a gene across samples are overdispersed relative to Poisson: the negative binomial, a gamma–Poisson mixture, has variance μ + αμ². DESeq2 (Love, Huber and Anders 2014) fits a GLM per gene with a log link, estimates α by Cox–Reid adjusted likelihood, shrinks it toward a trend in the mean (few samples estimate α poorly; many genes estimate the trend well), and tests fold changes by Wald's statistic, the p-values then adjusted for the false discovery rate (Benjamini and Hochberg 1995).

The theory, with the EM as a fixed point on the simplex and the negative binomial as a mixture, is [theory/msa-and-expression.md](../theory/msa-and-expression.md) §4–7.

*Remark.* The EM, the bootstrap, the negative binomial GLM, the dispersion's shrinkage, the tests and the FDR are general statistics. They are written here because genomics needs them first and the inference foundation is planned but absent; they move there when it exists (CLAUDE.md §8, open decision 3), and this module then keeps only the transcriptome, the pseudoalignment and the biological model of a sample.

## 2. Content

| Header | Content |
|---|---|
| `transcriptome.hpp` | `Transcript{name, gene, sequence}`, `Transcriptome{transcripts, genes, gene_of}` with `lengths()`; `make_transcriptome`, `read_transcriptome(fasta)` (the gene from `gene=` or `gene:` fields, as GENCODE and Ensembl write them); `FragmentLengths{p}` with `normal(mean, sd, max)`, `mean()`; `effective_lengths(lengths, d)` (kallisto's l − μ_l + 1); `TranscriptomeSpec`, `simulate_transcriptome(spec, seed)` → `SimulatedTranscriptome{transcriptome, exons}` (genes of exons, isoforms sharing them, paralogs); `simulate_abundances(T, sigma, zero, seed)`; `ReadModel{length, paired, error}`, `Fragment{transcript, begin, length, reverse, read1, read2}`, `fragment_fractions(rho, effective)`, `simulate_fragments(tx, rho, d, model, count, seed, first)` |
| `pseudoalign.hpp` | `ClassTable` (sorted sets of transcripts, each once, with `members`, `find`, `insert`); `Index(transcriptome, k = 31)` with `class_of(canonical k-mer)`, `classes()`, `kmers()`, `memory_bytes()`, and `pseudoalign(read, mate, scratch)`; `Pseudoalignment{classes, counts, fragment_class, fragments, aligned}`, `pseudoalign(index, reads, mates, pool)` and `pseudoalign(index, fragments, pool)` |
| `quantify.hpp` | `EmOptions` (kallisto's rule: 50 to 10 000 rounds, stop when every count above 0.01 changes by under 1 %), `Abundances{counts, tpm, effective_length, fragments, rounds, converged, log_likelihood}`, `quantify(classes, counts, effective, options, pool)` and `quantify(pseudoalignment, …)`; `log_likelihood`, `tpm`, `equal_split` (the baseline); `resample` (multinomial by conditional binomials), `bootstrap(…, replicates, seed, …)` → `Bootstrap{replicates, mean, sd}` |
| `differential.hpp` | `CountMatrix{genes, samples, counts}`; `median`; `size_factors` (the median of ratios); `Design{samples, p, x}` with `two_groups(group)`; `nb_log_likelihood`, `nb_draw` (gamma–Poisson); `fit_glm(y, s, X, α)` → `NbFit{beta, se, mu, log_likelihood, iterations, converged}` by IRLS; `dispersion_objective`, `dispersion_score` (in log α, Cox–Reid adjusted, with an optional normal prior), `estimate_dispersion`, `moments_dispersion`; `benjamini_hochberg`; `DifferentialOptions`, `differential(counts, design, options, pool)` → `Differential{size_factors, genes, trend_a0, trend_a1, trend_parametric, log_dispersion_variance, prior_variance}` with `GeneResult{base_mean, log2_fold_change, lfc_se, statistic, p_value, padj, dispersion_gene, dispersion_trend, dispersion, dispersion_outlier, converged}`; `simulate_counts(base, dispersion, log_fold, group, size, seed)` |

## 3. Conventions

- **k-mers are canonical**, k = 31 by default (odd, so no k-mer is its own reverse complement); reads of either strand pseudoalign alike (unstranded libraries).
- **Abundances.** α is the fraction of fragments, ρ the fraction of molecules (TPM = 10⁶ ρ), related by α_t ∝ ρ_t l̃_t. Estimated counts are N α, N the fragments pseudoaligned.
- **Fold changes** are log₂ of group 1 over group 0 in `Design::two_groups`; internally natural logarithms.
- **Randomness** (simulation, bootstrap) is a `calx::random::Stream` per transcript, fragment, gene or replicate, with this module's purposes (`TRANSCRIPTOME`, `ABUNDANCE`, `FRAGMENTS`, `BOOTSTRAP`, `COUNTS`); every parallel computation writes its own slots and is the same on any number of threads (checked).

## 4. Dependencies

`genomics::sequence` (codes, canonical k-mers, FASTA); `genomics::kmers` (`KmerMap` and `mix64` of `counter.hpp`: the index's open-addressing map from k-mers to classes); `calx::special_functions` (log Γ, ψ, ψ′, the incomplete gamma function for the χ²₁ tail), `calx::linalg` (the solves of IRLS and the Cox–Reid determinant), `calx::solvers` (Brent's root of the dispersion's score), `calx::random`, `calx::scheduler`. The checks also use `genomics::simulate`.

## 5. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time from the build (`check_genomics_expression_<topic>`).

| Check | What it verifies |
|---|---|
| `pseudoalign` | Transcripts AB, AC, BC on three exons with their six classes; a read inside an exon, across a junction, its reverse complement, a pair; every k-mer of 285 simulated transcripts (isoforms, paralogs) with the class of all (k-mer, transcript) pairs sorted, each class once; every read's and every pair's class (6000 each, errors 0.003 per base) equal to the brute-force intersection; on six transcripts, 300 reads against the definition itself (each k-mer searched for in every transcript and its reverse complement); the true transcript of every error-free read in its class; classes, numbering and counts the same on 1 and 8 threads; the speed (§6) |
| `quantify` | Effective lengths against a point mass (l − f + 1) and against Σ p(f)(l − f + 1)/P(F ≤ l); the EM at the closed form N n₀/(n₀ + n₁) of two transcripts sharing a class at equal lengths, and at the maximum of log L found by golden-section search at unequal ones (10⁻⁸); 400 000 simulated fragments of 936 transcripts of 400 genes (15 % paralogs): correlation with the true counts 0.996 (equal split 0.883), median relative error 0.074, mean 0.31 against 1.77 for the equal split, 0.038 on transcripts with 10 % of their reads their own, 0.003 on genes' totals; log L never decreasing over 811 rounds; TPM summing to 10⁶, genes' log TPM against the truth (1 − r = 0.003); the same on 1 and 8 threads; the bootstrap the same on 1 and 8 threads, its ±2 deviations covering the truth for 226 of 237 identifiable transcripts, and wider where reads are shared |
| `differential` | BH on a hand-worked list, with a NaN, and on Benjamini and Hochberg's fifteen p-values (four rejections at 0.05, q₍₄₎ = 0.035625, q₍₅₎ = 0.0603); size factors from exactly scaled counts (10⁻¹²) and from NB counts of 2000 genes (0.9 %); the NB's probabilities summing to 1 and equal to C(y + 2, 2) p³ (1 − p)^y at α = 1/3; gamma–Poisson draws of mean 20 and variance 120; IRLS at the closed form of a two-group design, and its score equations at zero with unequal size factors; at α = 0.1 (5 + 5 samples, mean 100) the median Cox–Reid estimate of 400 genes 0.092 (maximum likelihood 0.072), each a zero of its score and a local maximum; under the null of 4000 genes the Wald p-values' Kolmogorov–Smirnov distance from uniform 0.017 (1 % critical value 0.026; limit 0.04) and 6.4 % below 0.05; the trend's coefficients (0.031, 3.77 against 0.03, 3 scattered); 400 two-fold changes among 4000 genes found with power 0.87 at padj < 0.1, observed FDR 0.111; the same on 1 and 8 threads |

## 6. Performance

Measured 2026-10-24 on the i7-11370H in WSL (GCC 11.5, `-O3 -march=native`), with other builds running (load average about 2.5). The index of 4567 simulated transcripts (2000 genes, 5.5 Mb, 2.80 million distinct 31-mers, 8470 classes) is built in 0.42 s on one thread and takes 101 MB (`KmerMap` at a load between 5/16 and 5/8, 12 bytes a slot, and the classes). Pseudoalignment of 200 000 single-end reads of 75 bases runs at 0.59 million reads a second on one thread and 2.30 million on eight: a read's 45 lookups miss the cache, so their slots are prefetched together before any is read (2.5 times faster than looking them up one after another). The EM over 2221 classes converges in 811 rounds in 9 ms; below 2^18 incidences a round runs on the calling thread, since waking the pool 1622 times cost 1.8 s in WSL. The differential analysis of 4000 genes in 10 samples takes 0.72 s on one thread and 0.18 s on eight.

## 7. Deviations from kallisto and DESeq2

- **No skipping along unitigs.** kallisto jumps along the unitigs of the transcriptome's de Bruijn graph, on which the class is constant; here every k-mer is looked up (the same answer, slower by the length of the jumps).
- **The fragment-length law is given**, not estimated from pairs.
- **The likelihood leaves out P(F ≤ l)**, as kallisto's does; the simulator draws fragments by the same model, so the truth matches the likelihood.
- **DESeq2's outlier replacement (Cook's distances), independent filtering, fold-change shrinkage and local trend fit** are not done; the dispersion's gene-wise estimate takes two passes of fit and estimate (DESeq2 one).

## 8. References

- Anders, S. and Huber, W. (2010), Differential expression analysis for sequence count data, *Genome Biology* 11, R106.
- Benjamini, Y. and Hochberg, Y. (1995), Controlling the false discovery rate: a practical and powerful approach to multiple testing, *J. R. Statist. Soc. B* 57, 289–300.
- Bray, N. L., Pimentel, H., Melsted, P. and Pachter, L. (2016), Near-optimal probabilistic RNA-seq quantification, *Nature Biotechnology* 34, 525–527.
- Cox, D. R. and Reid, N. (1987), Parameter orthogonality and approximate conditional inference, *J. R. Statist. Soc. B* 49, 1–39.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *J. R. Statist. Soc. B* 39, 1–38.
- Li, B. and Dewey, C. N. (2011), RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome, *BMC Bioinformatics* 12, 323.
- Love, M. I., Huber, W. and Anders, S. (2014), Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, *Genome Biology* 15, 550.
- McCarthy, D. J., Chen, Y. and Smyth, G. K. (2012), Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation, *Nucleic Acids Research* 40, 4288–4297.
- McCullagh, P. and Nelder, J. A. (1989), *Generalized Linear Models*, 2nd ed., Chapman and Hall.
- Neuhaus, K. L. et al. (1992), Improved thrombolysis in acute myocardial infarction with front-loaded administration of alteplase, *J. Am. Coll. Cardiol.* 19, 885–891 (the p-values of Benjamini and Hochberg's example).
- Pachter, L. (2011), Models for transcript quantification from RNA-Seq, arXiv:1104.3889.
- Patro, R., Duggal, G., Love, M. I., Irizarry, R. A. and Kingsford, C. (2017), Salmon provides fast and bias-aware quantification of transcript expression, *Nature Methods* 14, 417–419.

## 9. Open Questions

1. **Identifiability.** An isoform whose every k-mer lies in other isoforms (all its exons and junctions shared) is not identifiable from classes; the EM then picks a point on a ridge of the likelihood (§5, the mean error of 0.31 against the median 0.074). Positions within transcripts (kallisto's pseudobam, Salmon's selective alignment) or fragment lengths of pairs would separate some of them; reporting the ridge is better than hiding it.
2. **Real data.** A GENCODE transcriptome and ENA reads by range requests, through the data library.
3. **The inference foundation** takes the statistics when it exists (§1).
