# sketch: Sketches of k-mer Sets

*Status: draft (2026-10-24).* Bottom-s MinHash and FracMinHash sketches of the sets of canonical k-mers of genomes; the Jaccard index, the Mash distance and its p-value, containment and ANI; HyperLogLog cardinalities, unions and Jaccard; all-against-all distance matrices as task graphs.

## 1. Purpose

A genome of millions of k-mers is compared with thousands of others through a sketch of a few thousand hashes. Each sketch is a homomorphism of some structure on sets, and that structure decides what it can estimate:

- **Bottom-s MinHash** preserves unions, S(A ∪ B) = bottom_s(S(A) ∪ S(B)), and estimates the Jaccard index with an error √(J(1 − J)/s) independent of the genomes' sizes (Broder 1997; Ondov et al. 2016, Mash).
- **FracMinHash** keeps every hash below a threshold, so it preserves both unions and intersections. It estimates containment, and from containment the ANI, for sets of very different sizes (Irber et al. 2022; Hera, Pierce-Ward and Koslicki 2023).
- **HyperLogLog** preserves unions as the register-wise maximum, in 2^p bytes, and estimates cardinalities with a relative error of 1.04/√m (Flajolet et al. 2007; Ertl 2017). It estimates Jaccard by inclusion and exclusion (Baker and Langmead 2019, Dashing).

The trees of genomes in the Genome Atlas are drawn from these distances. The theory, with the variances, is [theory/indexes.md §7](../theory/indexes.md).

## 2. Content

| Header | Content |
|---|---|
| `hash.hpp` | `mix64` (SplitMix64's finaliser, a bijection of 64-bit words), `kmer_hash(canonical, seed)`, `for_each_kmer_hash(sequence or contigs, k, seed, f)` |
| `minhash.hpp` | `MinHash{k, size, seed, hashes, kmers}` with `cardinality()`; `minhash(contigs, k, size, seed)` and `minhash_all(pool, genomes, …)`; `compare(a, b)` → `MashComparison{shared, considered, jaccard, distance, p_value}` (Ondov et al. eqs. 3, 4; 1, 5, 8); `mash_distance`, `mash_p_value` (with the k-mer space as an argument), `canonical_kmer_space`, `random_kmer_probability`, `binomial_upper_tail`. `FracMinHash{k, scaled, seed, hashes}` with `max_hash()`, `cardinality()`; `frac_minhash`; `intersection`, `containment` (Hera et al. eq. 2), `containment_debiased` (eq. 3), `frac_jaccard`, `ani` (C^{1/k}), `containment_variance` (Theorem 3) |
| `hyperloglog.hpp` | `HyperLogLog(p, seed)` with `insert(hash)`, `add(sequence or contigs, k)`, `merge`, `estimate()` (Ertl's improved raw estimator, eq. 32 and Algorithm 6), `estimate_classic()` (Flajolet et al., with linear counting; 64-bit hashes, after Heule et al.), `standard_error()`; `union_of`, `jaccard` |
| `distance.hpp` | `distance_matrix(pool, sketches, distance, block)` → `DistanceMatrix{n, values, block, tasks, trace}`: the blocks on and above the diagonal as tasks of a `calx::scheduler::TaskGraph`, each cell written once, deterministic; `mash_distances(pool, sketches, block)` |

## 3. Dependencies

`genomics::sequence` (canonical k-mers); `calx::scheduler` (the pool, parallel loops, task graphs); `calx::special_functions` (the regularised incomplete beta function, for the binomial tail of the p-value). The checks also use `genomics::simulate`.

## 4. Checks

| Check | What it verifies |
|---|---|
| `minhash` | On a simulated megabase, the sketch is the s smallest hashes of the exact k-mer set, and the same on the reverse strand; the cardinality from the sketch is within 4/√s. For pairs diverged at p from 0.001 to 0.2: the Jaccard estimate within 3.5 binomial standard errors of the exact Jaccard; the Mash distance within Ondov et al.'s Table 1 bounds of the distance of the exact sets; and that distance tracking d under Mash's Poisson model and −ln(1 − p) under independent substitutions. The binomial tail against 638/1024 and a direct sum. Unrelated genomes at k = 12 sharing about r·s hashes, r from the canonical k-mer space, and not significant (with Mash's 4^k they would be); related ones significant. FracMinHash: its size within 4σ of |A|·s, the containment of a subset 1 exactly, containment within 4σ of Theorem 3, ANI = Ĉ^{1/k} recovering 1 − p, Jaccard within 4σ |
| `hyperloglog` | Over 300 sets of 10 to 3·10⁵ elements at p = 10 and 14, the root mean square relative error of Ertl's estimate within 1.15 × 1.04/√m and its mean within 4 standard errors; the classical estimate above 6m likewise. Re-insertion changes nothing. The sketch of a union is the register-wise maximum. Jaccard by inclusion and exclusion within 4√3 (1 + J) ε. The distinct 21-mers of a 2 Mb genome within 4σ, and the reverse strand the same |
| `distance` | 40 genomes of a star phylogeny, sketched on the pool. The matrix symmetric with a zero diagonal, each entry the Mash distance of its sketches, the same to the bit on one and eight threads for blocks of 1, 7 and 32, with one task per block in the trace. The star's distances additive within 0.02, and each leaf's distance from the root within 0.01 of its rate |

*Caveat (Mash's null model).* Mash's p-value takes the k-mer space as 4^k (Ondov et al. eq. 1), though it sketches canonical k-mers, of which there are about half as many. Its p-values for canonical sketches are then too small: two unrelated random megabases at k = 12 share 63 hashes of 1000, where Mash's model expects 28 and the canonical one 56. `compare` uses the canonical space.

## 5. References

- Broder, A. Z. (1997), On the resemblance and containment of documents, *Compression and Complexity of Sequences*, 21–29.
- Ondov, B. D., Treangen, T. J., Melsted, P., Mallonee, A. B., Bergman, N. H., Koren, S. and Phillippy, A. M. (2016), Mash: fast genome and metagenome distance estimation using MinHash, *Genome Biology* 17, 132.
- Irber, L. et al. (2022), Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers, bioRxiv 2022.01.11.475838.
- Hera, M. R., Pierce-Ward, N. T. and Koslicki, D. (2023), Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash, *Genome Research* 33, 1061–1068.
- Flajolet, P., Fusy, É., Gandouet, O. and Meunier, F. (2007), HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm, *AofA*, DMTCS Proceedings AH, 127–146.
- Heule, S., Nunkesser, M. and Hall, A. (2013), HyperLogLog in practice, *EDBT*, 683–692.
- Ertl, O. (2017), New cardinality estimation algorithms for HyperLogLog sketches, arXiv:1702.01284.
- Baker, D. N. and Langmead, B. (2019), Dashing: fast and accurate genomic distances with HyperLogLog, *Genome Biology* 20, 265.
- Steele, G. L., Lea, D. and Flood, C. H. (2014), Fast splittable pseudorandom number generators, *OOPSLA*, 453–472.
- Topcuoglu, H., Hariri, S. and Wu, M.-Y. (2002), Performance-effective and low-complexity task scheduling for heterogeneous computing, *IEEE TPDS* 13, 260–274.

## 6. Open Questions

1. **Compatibility with Mash and sourmash.** Their sketches hash k-mer strings with MurmurHash3, so their hash values differ from ours, and sketches made by the two cannot be compared. Reading their sketch files would need their hash.
2. **Ertl's maximum-likelihood estimators** of unions and intersections (the "joint MLE" Dashing uses) would sharpen HyperLogLog's Jaccard for small J.
3. **The confidence intervals of Hera et al. (Theorems 7–9)** need the moments of the number of mutated k-mers (Blanca et al. 2022), which belong with a statistics module.
