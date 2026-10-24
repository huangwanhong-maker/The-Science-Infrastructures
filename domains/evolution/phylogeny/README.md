# phylogeny: Trees of Descent from Sequences

*Status: draft. Module of the evolution domain.* The trees that relate aligned sequences, and the methods that infer them: trees as split systems and tree metrics, distances under models of substitution, trees from distances, the likelihood of a tree under a continuous-time Markov chain with its branch lengths, topology and parameters fitted, parsimony as its tropical shadow, the bootstrap, simulation along trees, and the molecular clock read from dated tips. The theory, with the reasoning behind each algorithm, is [theory/phylogeny.md](../theory/phylogeny.md). Serves use, learning and research.

## 1. Concepts

- **A tree is its splits.** A branch cuts the leaves in two; a tree is determined by its splits, and a set of splits is the set of a tree exactly when they are pairwise compatible (Buneman). Robinson–Foulds distance and majority-rule consensus are operations on split systems.
- **A tree with lengths is a metric.** Leaf-to-leaf path lengths satisfy the four-point condition; neighbour joining recovers a tree from its tree metric exactly, and the topology from any matrix within half its shortest branch (Atteson).
- **A tree with a Markov chain is a probability model.** The likelihood of a column is a sum over hidden states of products along branches; Felsenstein's pruning evaluates it by the distributive law. With a reversible model the root may be anywhere (the pulley principle), so each branch is fitted as if it held the root.
- **Parsimony is the same computation in the tropical semiring** (min, +): Fitch's and Sankoff's algorithms are pruning with sums replaced by minima and products by sums.

## 2. Content

| Header | Content |
|---|---|
| `tree.hpp` | `Node`, `Tree` (parent, children, length, support, label; `rooted`); preorder and postorder; `reroot(t, v, x)`, `root_at`, `unroot`, `midpoint_root` (Farris 1972), `compact`, `ladderize`, `depths`, `distances_from`, `leaf_distances` |
| `newick.hpp` | `read_newick(text, NewickOptions)` → `NewickTrees {trees, problems}` (Felsenstein's grammar: quoted labels, nested comments, `[&R]`/`[&U]`, internal labels as supports; problems with line and column, never thrown); `write_newick(t, NewickWriting)` |
| `splits.hpp` | `Split` (the side without the first taxon, a bit set); `splits`, `compatible`, `robinson_foulds`, `normalised_robinson_foulds`, `tree_from_splits`, `majority_consensus` (Margush–McMorris; greedy below one half), `split_frequencies` |
| `layout.hpp` | `layout(t)` → rectangular (depth, leaf rank), circular, and equal-angle unrooted coordinates for drawing |
| `alignment.hpp` | `Alignment {names, rows}`, `alignment_from_letters`; `Patterns` (distinct columns with weights) by `compress`; `base_frequencies` |
| `distance.hpp` | `distances(patterns, DistanceOptions)` → `DistanceMatrix`: p, JC69, K80, F81, TN93 with Γ forms, LogDet (paralinear); pairwise deletion; saturation flagged and capped, never NaN; the single-pair formulas `jc69_distance`, `k80_distance`, `f81_distance`, `tn93_distance`, `logdet_distance` |
| `distance_trees.hpp` | `upgma(m, Linkage)` (UPGMA, WPGMA), `neighbour_joining(m)` (Saitou–Nei in Studier–Keppler's O(n³) form), `bionj(m)` (Gascuel 1997) |
| `models.hpp` | `Model` (reversible, π and exchangeabilities, normalised Q, spectral form; `transition(t, P[, dP, d²P])`); `jc69`, `k80`, `f81`, `hky85`, `tn93`, `gtr`; `discrete_gamma(α, K[, median])` (Yang 1994), `gamma_rates(α, K, p_inv)`, `Rates`, `RateSpec` |
| `likelihood.hpp` | `Likelihood(patterns, model, rates)`: `log_likelihood(t)` (pruning, rescaled), `optimise_branch_lengths(t)` (Newton on each branch, both-sided partials), `nni_search(t)`; `maximum_likelihood(patterns, model, RateSpec[, start], SearchOptions)` → `MLTree {tree, log_likelihood, model, rates, interchanges}` |
| `parsimony.hpp` | `fitch_score` (Fitch; Hartigan at polytomies), `sankoff(t, patterns, CostMatrix)` with a traceback, `unit_costs`, `transition_costs`; `parsimony_search` (stepwise addition, then NNI), `parsimony_exact` (branch and bound, Hendy–Penny) |
| `bootstrap.hpp` | `bootstrap_weights(patterns, seed, b)`; `bootstrap(patterns, reference, B, seed, pool, build)` on a `calx::scheduler::TaskGraph`, deterministic on any number of threads; `bootstrap_neighbour_joining` |
| `simulate.hpp` | `simulate_alignment(t, model, rates | SiteRates, sites, seed)` (Seq-Gen's method, discrete or continuous Γ, +I); `yule_tree`, `coalescent_tree`, `random_unrooted_tree` (uniform topologies), `scale_tree` |
| `dating.hpp` | `root_to_tip(t, labels, dates, ClockObjective)` → `ClockFit` (TempEst: the best root by residual mean square, R² or correlation, in closed form on each branch; rate, root date, residuals); `root_to_tip_at_root` |

Conventions: branch lengths in expected substitutions per site; node indices `int`, preorder after any change of shape; absent lengths and supports NaN; supports as fractions in [0, 1]; gaps and N missing data; randomness only by `calx::random::Stream(seed, purpose, object, 0)` with the purposes in `simulate.hpp` and `bootstrap.hpp`.

## 3. Dependencies

- **Domains:** genomics (`sequence`: codes, `Problem`).
- **Foundations:** calx — `linalg` (`symmetric_eigen` for the spectral decomposition, `determinant` for LogDet; `expm` in the checks), `optimisation` (Brent's method for the model parameters), `special_functions` (`gamma_p`, `gamma_p_inverse` for Γ rates), `random` (streams), `scheduler` (thread pool, task graphs, parallel loops), `cpu` (lanes for the pruning kernels), `verification` (checks).

## 4. Checks

`node domains/evolution/tests/run-tests.mjs`, or one program at a time from the domain's build (`ctest -R evolution_phylogeny`).

| Check | What it verifies |
|---|---|
| `tree` | Newick read and written back unchanged, quoted labels, nested comments, `[&R]`/`[&U]`, supports; faults found on their lines and columns; rerooting on every branch keeps every leaf distance and the splits; the midpoint root halves the longest path; n − 3 splits that give the tree back; Robinson–Foulds worked by hand; majority-rule consensus with frequencies and mean lengths; layouts; random topologies uniform over the 15 of five leaves (χ²); coalescent mean height 2(1 − 1/n), Yule mean height |
| `distance` | Kimura's worked example (K = 0.348) and his codon-position values; the Γ JC69 closed form (p = 0.3, α = 0.5 gives 2/3) and Jin and Nei's A8 and A11; TN93 → K80, F81 → JC69; LogDet = −tr(Q) d/4; pairwise deletion; saturation capped and flagged; JC69, K80, TN93, Γ-JC69 on simulated pairs within four standard errors; NJ on Wikipedia's example and on Saitou and Nei's Table 1, exact; NJ and BIONJ exact on random additive trees; NJ within Atteson's radius; UPGMA on Wikipedia's example (8.5, 11, 14, 16.5), ultrametric, exact on ultrametrics |
| `models` | For JC69, K80, F81, HKY85, TN93, GTR: rate 1, P(0) = I, rows sum to 1, P(s + t) = P(s)P(t), detailed balance, P → Π, P = expm(Qt) (calx), dP/dt = QP, d²P/dt² = Q²P; JC69 and K80 closed forms; Yang's discrete Γ rates for α = 0.5 (0.0334, 0.2519, 0.8203, 2.8944), exact at α = 1, medians; Γ + I |
| `likelihood` | Two sequences under JC69 against the closed form; the ML distance equal to the JC69 distance; pruning against the sum over all internal states (GTR + Γ + I, missing data); the pulley principle on every branch and node; rescaling on 700 leaves; fitted lengths at least as likely as the true ones; ML (HKY85 + Γ4) and NJ recovering a simulated 12-taxon tree, κ and α near the truth; NNI from a random start to a local maximum |
| `parsimony` | Fitch and Hartigan against enumeration of internal states; Sankoff against enumeration, with unit costs equal to Fitch, its traceback attaining its score; branch and bound against all 10,395 trees of 8 taxa; the heuristic recovering a 14-taxon tree |
| `bootstrap` | Replicate weights (L draws, a function of seed and replicate, mean the counts); supports near 1 for well-supported splits; a split with no branch behind it poorly supported; supports and every replicate tree identical on 1 and 8 threads |
| `dating` | On a strict clock with tips sampled over 30 years: the true root by each objective, the rate and the root's date exact, R² = 1; from ML lengths of 20,000 simulated sites, the rate within 5 % and the root's date within 5 years |

## 5. Performance

Measured on the development machine (i7-11370H, AVX-512, WSL, g++ 11.5 `-O3 -march=native`), best of repeated runs while other builds shared the machine:

| Operation | Size | Time |
|---|---|---|
| Pruning (log-likelihood) | 50 taxa, 4,141 patterns, Γ4 | 4.0 ms, about 5 ns per internal node, pattern and category |
| One sweep of all branch lengths | the same | 47 ms |
| Neighbour joining / BIONJ / UPGMA | 1,000 taxa | 0.24 s / 0.25 s / 0.02 s |
| TN93 distance matrix | 50 taxa, 5,000 sites | 9 ms |
| NJ bootstrap, 200 replicates | 50 taxa, 5,000 sites | 3.0 s on 1 thread, 0.60 s on 8 (5.0 times) |
| Maximum likelihood, HKY85 + Γ4 with κ and α, from the true tree | 12 taxa, 3,000 sites (1,383 patterns) | 1.8 s |
| NNI search from a random tree to the true one (13 interchanges) | the same | 0.7 s |

## 6. References

- Felsenstein, J. (1981), Evolutionary trees from DNA sequences: a maximum likelihood approach, *J. Mol. Evol.* 17, 368–376.
- Felsenstein, J. (1985), Confidence limits on phylogenies: an approach using the bootstrap, *Evolution* 39, 783–791.
- Felsenstein, J. (2004), *Inferring Phylogenies*, Sinauer.
- Saitou, N. and Nei, M. (1987), The neighbor-joining method, *Mol. Biol. Evol.* 4, 406–425; Studier, J. A. and Keppler, K. J. (1988), *Mol. Biol. Evol.* 5, 729–731; Gascuel, O. (1997), BIONJ, *Mol. Biol. Evol.* 14, 685–695; Atteson, K. (1999), *Algorithmica* 25, 251–278.
- Sokal, R. R. and Michener, C. D. (1958), *University of Kansas Science Bulletin* 38, 1409–1438.
- Jukes, T. H. and Cantor, C. R. (1969); Kimura, M. (1980), *J. Mol. Evol.* 16, 111–120; Hasegawa, M., Kishino, H. and Yano, T. (1985), *J. Mol. Evol.* 22, 160–174; Tamura, K. and Nei, M. (1993), *Mol. Biol. Evol.* 10, 512–526; Tavaré, S. (1986), *Lectures on Mathematics in the Life Sciences* 17, 57–86; Jin, L. and Nei, M. (1990), *Mol. Biol. Evol.* 7, 82–102; Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994), *Mol. Biol. Evol.* 11, 605–612; Lake, J. A. (1994), *PNAS* 91, 1455–1459.
- Yang, Z. (1994), Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites, *J. Mol. Evol.* 39, 306–314; Yang, Z. (2014), *Molecular Evolution: A Statistical Approach*, Oxford.
- Guindon, S. and Gascuel, O. (2003), A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood, *Syst. Biol.* 52, 696–704.
- Fitch, W. M. (1971), *Systematic Zoology* 20, 406–416; Hartigan, J. A. (1973), *Biometrics* 29, 53–65; Sankoff, D. (1975), *SIAM J. Appl. Math.* 28, 35–42; Hendy, M. D. and Penny, D. (1982), *Mathematical Biosciences* 59, 277–290.
- Robinson, D. F. and Foulds, L. R. (1981), *Mathematical Biosciences* 53, 131–147; Margush, T. and McMorris, F. R. (1981), *Bulletin of Mathematical Biology* 43, 239–244; Buneman, P. (1971), in *Mathematics in the Archaeological and Historical Sciences*, 387–395; Semple, C. and Steel, M. (2003), *Phylogenetics*, Oxford.
- Rambaut, A. and Grassly, N. C. (1997), Seq-Gen, *CABIOS* 13, 235–238; Rambaut, A., Lam, T. T., Max Carvalho, L. and Pybus, O. G. (2016), TempEst, *Virus Evolution* 2, vew007.

## 7. Open Questions

1. **UPGMA's place.** The genomics domain's multiple aligner needs a guide tree and cannot depend on this domain (height 1). UPGMA and neighbour joining name no object in nature; they belong in a foundation (morphe, beside its combinatorics of trees, or the planned inference foundation), with `DistanceMatrix` and the `Tree` structure, so that genomics::msa and evolution::phylogeny both use them.
2. **More states.** The models and the likelihood are for four states; codon models (Goldman–Yang, Muse–Gaut) and amino-acid models (WAG, LG) need the same machinery over 61 and 20 states.
3. **Larger trees.** SPR moves, lazy recomputation of partials after an interchange, and lanes across patterns (eight patterns to an AVX-512 vector) would bring hundreds of taxa within reach.
