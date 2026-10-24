# msa: Multiple Sequence Alignment

*Status: draft (2026-10-24).* The alignment of many sequences at once: the alignment and its formats, guide trees from k-mer distances by UPGMA, progressive alignment of profiles with MUSCLE's re-estimated tree and refinement, the centre star, alignment by partial orders, and the scores that judge them: the sum of pairs against a scoring, and Q and TC against a reference.

## 1. Purpose

A multiple alignment says which residues of many sequences descend from one residue of their ancestor. Its natural objective, the sum over all pairs of the scores of their projected pairwise alignments, is NP-hard to optimise (Wang and Jiang 1994), so every practical method is a heuristic with a reason:

- **Progressive alignment** (Feng and Doolittle 1987; Higgins and Sharp 1988, CLUSTAL) aligns along a guide tree, closest relatives first, each node aligning the profiles of its two children: a greedy order on the tree, "once a gap, always a gap". MUSCLE (Edgar 2004) makes the tree cheap (k-mer distances), re-estimates it from the first alignment, and then refines: each edge of the tree splits the rows in two, and the two sides are re-aligned and kept when the SP score rises.
- **The centre star** (Gusfield 1993) aligns every sequence to the one closest to all the others and merges the pairwise alignments; under a metric cost it is within a factor 2 − 2/k of the optimum.
- **Partial orders** (Lee, Grasso and Sharlow 2002) align each sequence to a graph of all the previous ones rather than to a profile, so that shared indels are branches.

The theory, with the proofs of the bounds, is [theory/msa-and-expression.md](../theory/msa-and-expression.md) §1–3. The Genome Atlas draws the alignments in its Trees view; the evolution domain's phylogeny reads them.

## 2. Content

| Header | Content |
|---|---|
| `msa.hpp` | `Msa{names, rows}` (codes, `sequence::GAP` in the gaps) with `size()`, `columns()`, `rectangular()`; `Alphabet{dna, protein}`, `letters(a)`; `ungapped`, `sequences`, `single`; `column`, `residues_in_column`, `gap_fraction(s)`, `slice`, `without_gap_columns`, `sub_alignment(m, rows)`; `consensus(m, letters, max_gap)`; `column_entropy(m, letters, count_gaps)` and `column_information` (Shannon, bits); `clustal_conservation(m, alphabet)` (ClustalW's `*`, `:`, `.` by its strong and weak groups); `read_aligned_fasta` / `write_aligned_fasta`, `read_clustal` / `write_clustal`, readers reporting `Problem{line, message}` |
| `guide_tree.hpp` | `DistanceMatrix{n, d}`; `Compression` with `dna_letters()`, `protein_letters()`, `dayhoff6()`; `KmerProfile`, `kmer_profile(x, k, compression)`, `kmer_similarity` (MUSCLE's fractional common k-mer count F), `kmer_distances(seqs, k, compression, pool)` (1 − F, in parallel by rows); `p_distance`, `jukes_cantor`, `kimura_protein`, `corrected_distances(msa, alphabet, pool)`; `Merge{a, b, distance, size}`, `Dendrogram{leaves, merges}` (SciPy's linkage) with `height`, `parents`, `leaves_under`, `newick`; `upgma(d)` by the nearest-neighbour chain, O(n²), and `upgma_naive(d)`, O(n³); `clustal_weights(tree)` (Thompson, Higgins and Gibson 1994) |
| `progressive.hpp` | `Profile` (weighted residue frequencies, sparse, and occupancies) and `make_profile(msa, weights, sigma)`; `align_profiles(x, y, scoring, terminal_open)` → `ProfileAlignment{score, path}`: Gotoh's recurrences with sum-of-pairs column scores and gap opens scaled by the occupancies on both sides (MUSCLE's open and close halves); `merge(x, y, path)`; `align_along(seqs, names, tree, weights, scoring, terminal_open, pool)` (the tree's nodes as a `calx::scheduler::TaskGraph`); `refine(msa, tree, weights, …, rounds, pool, report)` (tree-dependent restricted partitioning, `Refinement{rounds, trials, accepted, gain}`); `ProgressiveOptions`, `default_scoring(alphabet)`, `default_options(alphabet)`, `progressive(seqs, names, alphabet, options, pool)` → `ProgressiveResult{msa, tree, weights, sp_first, sp_second, sp, second_kept, refinement}` |
| `centre_star.hpp` | `centre_star(seqs, names, scoring, pool)` → `CentreStar{msa, centre, scores}`: all pairwise global scores in parallel, the centre, the star merged so that it induces each alignment to the centre exactly |
| `partial_order.hpp` | `partial_order_alignment(seqs, names, scoring, order, graph)`: the sequences added to an `align::PoaGraph` in a given order, and the alignment the graph holds |
| `scores.hpp` | `pair_score(a, b, scoring, terminal_open)` (the projection, natural affine gap costs), `sp_score(msa, scoring, terminal_open, pool)`, `sp_across(msa, left, right, …)`; `Accuracy{q, tc, precision, pairs_*, columns_*}` and `accuracy(test, reference)`: BAliBASE's SP (Q) and TC scores, and the precision of the pairs |

## 3. Conventions

- **Rows are in the input's order** in every result, whatever order the tree or the graph visited them in.
- **Scorings** are `align::Scoring`: a substitution matrix on the codes and affine gaps, a gap of length L costing open + L · extend. The defaults (`default_scoring`) are DNA +2/−1 with open 4, extend 1, and protein BLOSUM62 with open 6, extend 1, chosen by the grid of §6.2.
- **Terminal gaps** pay their extension but not their opening by default (`terminal_open = false`), in the profile alignment and in the SP objective alike; `sp_score(…, true)` charges them in full.
- **Parallel and deterministic.** Distances by rows, the tree's nodes as a task graph (siblings at once), SP scores by rows; every task writes its own slot, so the alignment is the same on any number of threads (checked).

## 4. Dependencies

`genomics::sequence`; `genomics::align` (`types.hpp`, `scoring.hpp`, `gotoh.hpp` for the centre star, `poa.hpp` for partial orders); `calx::scheduler` (pool, parallel loops, task graphs). The checks also use `genomics::simulate` (families with their true alignment), `align/statistics.hpp` (BLOSUM62's λ for the protein families) and `calx::random`.

*Remark.* UPGMA, CLUSTAL W's weights and the corrected distances are clustering and phylogenetics, written here because the evolution domain stands on genomics and no inference foundation exists yet. They belong in that foundation's clustering when it does; the evolution domain's phylogeny can then reuse them.

## 5. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time from the build (`check_genomics_msa_<topic>`).

| Check | What it verifies |
|---|---|
| `msa` | Aligned FASTA and Clustal written and read back exactly (40 random DNA and protein alignments); a ragged FASTA padded and reported at its lines; a file in Clustal Omega's form (human α and β globin, UniProt P69905 and P68871, aligned by hand) read whole; a Clustal file without its header read and reported; ClustalW's marks on columns built for its strong and weak groups; columns, gap fractions, slices, sub-alignments; the consensus; column entropies 0, 2, 1 bits and log₂ 20, and information as their complement |
| `scores` | Pair scores of hand-worked projections (two runs, terminal runs with and without their openings, a row of gaps only); the pair score of 500 random row pairs equal to `align::rescore` of their projection; the SP score the sum of pair scores on 1 and 8 threads, and its part across a split plus the two sides' SP scores; Q = 6/8, TC = 3/4, precision 1 on a hand-worked test; an alignment against itself 1, 1, 1; Q(A, B) = precision(B, A) on 100 random pairs; different sequences refused |
| `guide_tree` | k-mer profiles against counts by brute force; F = 2/4 on ACGTACGT and ACGTTT; Dayhoff(6) making residues of a group one class; the distance matrix symmetric, in [0, 1] and the same on 1 and 8 threads; Jukes–Cantor and Kimura at p = 0.3 against their closed forms; UPGMA of the 5S rRNA distances (Erdmann and Wolters 1986, the worked example of Wikipedia's UPGMA article) merging at heights 8.5, 11, 14, 16.5, with its Newick string; the nearest-neighbour chain equal to the greedy algorithm on 200 random matrices and 200 planar point sets; ultrametrics giving back their trees; CLUSTAL W's weights 6/7, 6/7, 9/7 on ((A:1, B:1):2, C:3) |
| `progressive` | The profile alignment of two sequences equal to Gotoh's global optimum, score and path (300 pairs); eleven DNA families and seven protein families (table, §6.1): every row of every method's alignment without its gaps its sequence; on the families of at most 16 members at substitution ≤ 0.05 and indels ≤ 0.01 per branch, Q ≥ 0.95 and TC ≥ 0.75; protein at substitution ≤ 0.08, Q ≥ 0.95; the progressive Q never below the centre star's (less 0.005) wherever the family is not saturated (mean identity ≥ 0.5 for DNA, ≥ 0.25 for protein); the SP score reported by the refinement that of the alignment, never below the first stage's; the same alignment on 1 and 8 threads (33 sequences), and the centre star too; the centre star inducing each optimal alignment to the centre exactly, and its SP cost under edit distance within (2 − 2/k) Σ D(i, j), the lower bound of the optimum (30 families; largest ratio 1.42) |

## 6. Accuracy and Performance

### 6.1 Accuracy on Simulated Families

`simulate::evolve_family` evolves a root along a balanced tree with the given substitution and indel rates *per branch*, so that the leaves of n members are up to 2 log₂ n branches apart; the identities are those of the true alignment's pairs (mean / least). Q / TC against the truth; dSP is the SP score (terminal openings free) minus the truth's, per pair of rows: positive means the method found an alignment the scoring prefers to the truth.

DNA (+2/−1, open 4, extend 1), root lengths 250–300:

| n | sub | indel | identity | progressive | first stage | centre star | partial order | dSP prog / cs / poa |
|---|---|---|---|---|---|---|---|---|
| 8 | 0.01 | 0.003 | 0.959 / 0.928 | 0.986 / 0.945 | 0.985 / 0.948 | 0.982 / 0.942 | 0.987 / 0.955 | +1.3 / −1.7 / −0.6 |
| 8 | 0.03 | 0.010 | 0.865 / 0.809 | 0.972 / 0.846 | 0.950 / 0.787 | 0.882 / 0.627 | 0.958 / 0.846 | +4.3 / −21.6 / −2.6 |
| 8 | 0.05 | 0.010 | 0.781 / 0.722 | 0.964 / 0.869 | 0.962 / 0.869 | 0.936 / 0.790 | 0.934 / 0.803 | +14.2 / +5.0 / +1.4 |
| 16 | 0.03 | 0.005 | 0.817 / 0.732 | 0.971 / 0.830 | 0.973 / 0.861 | 0.966 / 0.737 | 0.948 / 0.799 | −0.3 / −2.5 / −12.1 |
| 32 | 0.01 | 0.003 | 0.912 / 0.851 | 0.991 / 0.892 | 0.988 / 0.892 | 0.979 / 0.763 | 0.983 / 0.805 | +1.7 / −2.9 / −0.5 |
| 32 | 0.03 | 0.010 | 0.792 / 0.687 | 0.895 / 0.268 | 0.864 / 0.365 | 0.795 / 0.134 | 0.794 / 0.213 | +5.5 / −36.1 / −41.7 |
| 32 | 0.06 | 0.003 | 0.629 / 0.509 | 0.945 / 0.569 | 0.937 / 0.691 | 0.880 / 0.353 | 0.823 / 0.369 | +6.2 / −20.6 / −34.5 |
| 40 | 0.06 | 0.010 | 0.612 / 0.452 | 0.792 / 0.121 | 0.721 / 0.156 | 0.662 / 0.049 | 0.634 / 0.066 | +10.0 / −23.7 / −60.3 |
| 8 | 0.10 | 0.020 | 0.619 / 0.492 | 0.805 / 0.479 | 0.782 / 0.461 | 0.606 / 0.200 | 0.601 / 0.255 | +28.7 / −20.8 / −22.8 |
| 16 | 0.20 | 0.030 | 0.364 / 0.226 | 0.117 / 0.000 | 0.144 / 0.000 | 0.084 / 0.000 | 0.147 / 0.000 | saturated |
| 8 | 0.30 | 0.050 | 0.331 / 0.204 | 0.087 / 0.000 | 0.066 / 0.000 | 0.077 / 0.000 | 0.101 / 0.000 | saturated |

Protein (BLOSUM62, open 6, extend 1; substitutions drawn from BLOSUM62's target frequencies), root length 200:

| n | sub | indel | identity | progressive | centre star | partial order |
|---|---|---|---|---|---|---|
| 8 | 0.03 | 0.010 | 0.835 / 0.758 | 0.980 / 0.914 | 0.969 / 0.861 | 0.983 / 0.909 |
| 8 | 0.08 | 0.010 | 0.677 / 0.552 | 0.977 / 0.899 | 0.955 / 0.825 | 0.963 / 0.880 |
| 8 | 0.15 | 0.010 | 0.446 / 0.320 | 0.963 / 0.868 | 0.846 / 0.621 | 0.955 / 0.840 |
| 32 | 0.03 | 0.010 | 0.768 / 0.653 | 0.965 / 0.725 | 0.934 / 0.511 | 0.966 / 0.714 |
| 32 | 0.08 | 0.003 | 0.534 / 0.381 | 0.982 / 0.808 | 0.946 / 0.571 | 0.974 / 0.848 |
| 32 | 0.15 | 0.003 | 0.303 / 0.169 | 0.932 / 0.592 | 0.673 / 0.085 | 0.909 / 0.552 |
| 24 | 0.30 | 0.020 | 0.136 / 0.030 | 0.102 / 0.000 | 0.097 / 0.000 | 0.118 / 0.016 |

*Remark.* The per-branch rates the task named (substitution 0.05 to 0.3, indels 0.01 to 0.05) saturate families of 16 or more members: at 0.2 per branch over eight branches DNA keeps 36 % identity, near the 25 % of unrelated sequence, and no aligner recovers the truth there (Q about 0.1 for all four methods). The table therefore spans the alignable range and keeps two saturated families to show the limit. Where the families are alignable the full progressive method has a higher Q than the centre star and partial orders on every DNA family but the first (where partial orders lead by 0.001) and on four of the six protein ones (partial orders lead by 0.001 and 0.003 on the other two); the second stage and the refinement add up to 0.07 to the first stage's Q, most on the large indel-rich families, and once cost 0.002. Its SP score is at or above the truth's (at worst 0.3 below, per pair), so its remaining errors belong to the scoring, not to the search.

### 6.2 The Default Scorings

A grid over DNA scorings (+5/−4, +2/−3, +1/−2, +2/−1, +1/−1, with opens 2 to 12 and extensions 1 and 2) on twelve families of 8 and 32 members gave mean Q from 0.944 to 0.948 among the best five; +2/−1, open 4, extend 1 was best (0.9481). For protein, BLOSUM62 with (open, extend) in (11, 1), (8, 1), (6, 1), (4, 1), (10, 2), (6, 2) gave mean Q 0.949, 0.957, 0.961, 0.950, 0.950, 0.957 on twelve families; (6, 1) is the default.

### 6.3 Speed

Measured 2026-10-24 on the i7-11370H in WSL (GCC 11.5, `-O3 -march=native`), with other builds running (load average about 2). The profile alignment runs at 7.0 ns a cell for two single sequences of 2 kb and 8.9 ns a cell for two profiles of 20 rows (1265 × 1265 columns), where `align::gotoh` runs at 7.3 ns a cell: its pair scores and its gap openings into x are computed a whole row at a time on contiguous vectors, which the compiler puts on lanes, and only the recurrences of D and H run along the row. A progressive alignment of 8 sequences of 300 bases, with both later stages, takes 25 ms on one thread; of 40, 315 ms on one thread and 208 ms on eight, two thirds of it in the refinement's 154 re-alignments of two profiles. The centre star of 40 × 300 takes 53 ms on eight threads.

*Remark.* A parallel loop costs a wake-up of the pool's threads, slow in WSL; the SP scores therefore run on the calling thread below 2^20 column visits, where the wake-up would cost more than the loop.

## 7. References

- Altschul, S. F. (1989), Gap costs for multiple sequence alignment, *J. Theor. Biol.* 138, 297–309.
- Carrillo, H. and Lipman, D. (1988), The multiple sequence alignment problem in biology, *SIAM J. Appl. Math.* 48, 1073–1082.
- Edgar, R. C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, *Nucleic Acids Research* 32, 1792–1797; MUSCLE: a multiple sequence alignment method with reduced time and space complexity, *BMC Bioinformatics* 5, 113; Local homology recognition and distance measures in linear time using compressed amino acid alphabets, *Nucleic Acids Research* 32, 380–385.
- Feng, D.-F. and Doolittle, R. F. (1987), Progressive sequence alignment as a prerequisite to correct phylogenetic trees, *J. Mol. Evol.* 25, 351–360.
- Gusfield, D. (1993), Efficient methods for multiple sequence alignment with guaranteed error bounds, *Bull. Math. Biol.* 55, 141–154; (1997), *Algorithms on Strings, Trees, and Sequences*, Cambridge University Press, §14.6.
- Higgins, D. G. and Sharp, P. M. (1988), CLUSTAL: a package for performing multiple sequence alignment on a microcomputer, *Gene* 73, 237–244.
- Jukes, T. H. and Cantor, C. R. (1969), Evolution of protein molecules, in Munro, H. N. (ed.), *Mammalian Protein Metabolism* III, Academic Press, 21–132; Kimura, M. (1983), *The Neutral Theory of Molecular Evolution*, Cambridge University Press.
- Kececioglu, J. and Starrett, D. (2004), Aligning alignments exactly, *Proceedings of RECOMB 2004*.
- Lee, C., Grasso, C. and Sharlow, M. F. (2002), Multiple sequence alignment using partial order graphs, *Bioinformatics* 18, 452–464.
- Müllner, D. (2011), Modern hierarchical, agglomerative clustering algorithms, arXiv:1109.2378; Murtagh, F. (1983), A survey of recent advances in hierarchical clustering algorithms, *The Computer Journal* 26, 354–359.
- Sauder, J. M., Arthur, J. W. and Dunbrack, R. L. (2000), Large-scale comparison of protein sequence alignment algorithms with structure alignments, *Proteins* 40, 6–22.
- Sokal, R. R. and Michener, C. D. (1958), A statistical method for evaluating systematic relationships, *University of Kansas Science Bulletin* 38, 1409–1438.
- Thompson, J. D., Higgins, D. G. and Gibson, T. J. (1994), CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice, *Nucleic Acids Research* 22, 4673–4680.
- Thompson, J. D., Plewniak, F. and Poch, O. (1999), A comprehensive comparison of multiple sequence alignment programs, *Nucleic Acids Research* 27, 2682–2690.
- Wang, L. and Jiang, T. (1994), On the complexity of multiple sequence alignment, *J. Comput. Biol.* 1, 337–348.

## 8. Open Questions

1. **The log-expectation score** of MUSCLE (Edgar 2004) and consistency (T-Coffee's library, ProbCons' posterior probabilities from the pair HMM of `align/pairhmm.hpp`) score columns better than sum-of-pairs; the second would join two modules already here.
2. **Re-aligning only what changed** in the second stage, as MUSCLE does, would halve its cost; here the whole tree is re-aligned.
3. **Real benchmarks.** BAliBASE's reference alignments (core blocks) would hold the aligners to structure-based truth rather than simulated truth.
4. **A foundation for clustering.** UPGMA and the weights move to the inference foundation when it exists (§4).
