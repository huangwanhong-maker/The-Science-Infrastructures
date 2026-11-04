# Rugged Landscapes: NK Landscapes of Twenty Bits, Counted Whole

*Status: draft (2026-11-05). Research: a case study of heuresis's `landscapes` module and of six of its searches ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).*

## 1. The Question

A search can beat random search only on problems whose structure it exploits (Wolpert and Macready 1997; heuresis's theory note, §2). Kauffman's NK landscapes vary one kind of structure by one integer. Each of n bits contributes a value that depends on itself and on k other bits. At k = 0 the landscape is separable and has one optimum. At k = n − 1 every contribution depends on every bit, and the values of neighbouring points are uncorrelated.

With n = 20 the whole space, 2²⁰ = 1 048 576 points, can be counted. The study asks:

1. how the measures of ruggedness that heuresis computes (local optima, basins, autocorrelation, correlation length, fitness–distance correlation) change with k, measured exactly rather than estimated;
2. how seven searches, each with the same fixed budget, fare as k grows, judged against the enumerated optimum.

## 2. The Method

**The landscapes.** `heuresis::benchmarks::NK(20, k, 2026, instance)` draws, for each bit, k neighbours uniformly among the other bits and a table of 2^(k+1) values uniform in [0, 1). The fitness F(x) is the mean of the 20 table entries that x selects, and the objective, minimised, is 1 − F(x). There are five instances for each k ∈ {0, 1, …, 10, 19}, sixty landscapes in all.

**Exact structure.** Each landscape is tabulated whole by `landscapes::Hypercube::tabulate`. Then, with one-bit flips as the neighbourhood:
- the global optimum is the least value of the table;
- the local optima are the points that no neighbour improves (`local_minima`);
- their basins are found by steepest descent, ties to the lower bit (`basins`);
- the one-step correlation is ρ(1) = 1 − E[(f(x) − f(y))²]/(2 Var f) over every point x and neighbour y (`one_step_correlation`);
- the correlation length is −1/ln ρ(1);
- the fitness–distance correlation is taken over all 2²⁰ points, with the Hamming distance to the global optimum (`fitness_distance_correlation`).

A random walk of 2 × 10⁵ one-bit steps (`random_walk`, `autocorrelation`) estimates ρ(s) for s ≤ 20, as one would on a space too large to count.

**The searches.** Each has the same budget of evaluations, 5 000 (0.48% of the space) and, separately, 50 000. There are ten seeds per instance, so fifty runs per k, search and budget:

| Search | heuresis | Settings |
|---|---|---|
| random search | `search::RandomSearch<Bits>` | uniform strings |
| (1+1) EA | `genetic::one_plus_one_ea` | standard bit mutation at 1/n |
| GA | `genetic::GeneticAlgorithm<Bits>` | μ = 50, uniform crossover at 0.9, bit flips at 1/n, binary tournaments, generational with one elite |
| simulated annealing | `annealing::SimulatedAnnealing<Bits>` | one-bit flips; the modified Lam schedule from T₀ = 0.05 over the budget |
| tabu search | `local::TabuSearch<Bits>` | the one-flip neighbourhood, tenure 5, a random start |
| compact GA | `distributions::DistributionSearch` with `Compact` | K = 50 |
| UMDA | `distributions::DistributionSearch` with `Umda` | λ = 100, μ = 50, borders 1/n |

None restarts. Each run is driven by `search::run`, and its best value is compared with the enumerated optimum:
- the **gap** F* − F(best);
- the **rank** of the best value among all 2²⁰ points (0 for the optimum);
- **success**, whether the run found the optimum.

## 3. The Code

- heuresis `6b3cc88`; calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`; the infrastructure's root `4153cb9`, in which heuresis is a submodule (`results/versions.json`). The runs were made from the working tree that became that commit.
- `rugged.cpp`, built by `CMakeLists.txt` as a project of its own (C++20, Release, g++ 11.5 in WSL):

```
cmake -S . -B ~/.cache/science-build/exp-nk -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-nk -j4
~/.cache/science-build/exp-nk/rugged results
```

- `figures.py` draws the figures from `results/` in the house style (`../style/search_style.py`).

The program was run twice, the second time with a per-instance table added. Every other table came out identical to the byte, as heuresis's keyed draws promise.

## 4. The Parameters

All of them are in `results/reference.json`: n = 20; k ∈ {0, …, 10, 19}; five instances (NK seed 2026, instances 0–4); ten seeds per instance; budgets of 5 000 and 50 000; a walk of 2 × 10⁵ steps. The run took 312 s on 6 threads of a machine shared with other work.

## 5. The Data

None from outside. The landscapes are heuresis's benchmark, and every number comes from the run.

## 6. The Results

### 6.1 Ruggedness Grows Smoothly with k, and Weinberger's Law Holds Once the Tables Are Counted

Means over the five instances of each k (ranges in brackets):

| k | Local optima | ρ(1), exact | Correlation length | Fitness–distance correlation | Global optimum's basin |
|---|---|---|---|---|---|
| 0 | 1 | 0.900 | 9.49 | 0.81 | the whole space |
| 1 | 18 [2–44] | 0.863 | 6.81 | 0.58 | 44.6% |
| 2 | 51 [19–85] | 0.826 | 5.27 | 0.36 | 16.4% |
| 3 | 167 [136–194] | 0.788 | 4.20 | 0.20 | 7.2% |
| 4 | 425 [275–539] | 0.743 | 3.36 | 0.10 | 3.4% |
| 6 | 1 582 [1 506–1 629] | 0.647 | 2.30 | 0.04 | 1.3% |
| 8 | 3 811 [3 740–3 902] | 0.547 | 1.66 | 0.02 | 0.40% |
| 10 | 7 540 [7 351–7 829] | 0.450 | 1.25 | 0.01 | 0.17% |
| 19 | 49 940 [49 830–50 119] | 0.0001 | — | 0.000 | 0.004% |

At k = 19 the landscape is uncorrelated. A point is then a local optimum with probability 1/(n + 1), since it must be the best of itself and its n neighbours, so the expected count is 2²⁰/21 = 49 932. The five instances give 49 830 to 50 119.

Weinberger's ρ(1) = 1 − (k + 1)/n (1991) assumes tables of infinitely many entries. Each table here has only 2^(k+1) entries, and their empirical variance has expectation σ²(1 − 2^−(k+1)), while a flip that changes an entry changes it by a mean square of 2σ². The deficit (k + 1)/n therefore grows by the factor 2^(k+1)/(2^(k+1) − 1). With that correction the prediction meets the exact mean within 0.004 at every k: 0.867 against 0.863 at k = 1, 0.829 against 0.826 at k = 2, and 0.900 at k = 0, which is the elementary value 1 − 2/n of heuresis's check. Uncorrected, it is off by 0.05 at k = 0 and by 0.04 at k = 1. The walk's estimate agrees with the exact ρ(1) within 0.005 everywhere.

The autocorrelation along the walk decays more slowly than ρ(1)ˢ for 0 < k < 19, because an NK landscape is a mixture of Walsh orders 1 to k + 1 and not an elementary landscape. At k = 19, ρ(2) ≈ 0.05 = 1/n, the chance that a walk's second step undoes its first.

From k = 4 on, the local optima lie at a mean Hamming distance of 9.9 to 10.0 from the global optimum, that of random points (n/2), and the fitness–distance correlation is at most 0.10: there is no big valley left to follow. Within a landscape, though, fitter local optima have larger basins: the correlation of F with the logarithm of the basin size, averaged over the instances, is 0.84–0.90 at every k from 1 to 10, as Ochoa, Tomassini, Vérel and Darabos (2008) found for NK landscapes. The global optimum has the largest basin in 19 of the 60 landscapes, including all five at k = 0, four of five at k = 1, and none above k = 8.

### 6.2 The Searches Part Ways from k = 2

The share of 50 runs (five instances × ten seeds) that found the optimum:

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 19 |
|---|---|---|---|---|---|---|---|---|---|---|
| **5 000 evaluations** | | | | | | | | | | |
| random search | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| (1+1) EA | 1 | 0.86 | 0.36 | 0.32 | 0.10 | 0.08 | 0.02 | 0 | 0.02 | 0 |
| GA | 1 | 1 | 0.72 | 0.38 | 0.14 | 0.28 | 0.12 | 0.02 | 0 | 0.02 |
| simulated annealing | 1 | 0.98 | 0.92 | 0.74 | 0.44 | 0.30 | 0.26 | 0.10 | 0.08 | 0 |
| tabu search | 1 | 0.90 | 0.54 | 0.40 | 0.54 | 0.70 | 0.52 | 0.28 | 0.18 | 0 |
| compact GA | 1 | 0.90 | 0.36 | 0.26 | 0.04 | 0.08 | 0.06 | 0 | 0 | 0 |
| UMDA | 1 | 0.98 | 0.58 | 0.22 | 0.08 | 0.04 | 0.08 | 0 | 0.02 | 0 |
| **50 000 evaluations** | | | | | | | | | | |
| random search | 0.10 | 0.04 | 0.06 | 0.02 | 0.02 | 0.08 | 0.10 | 0.04 | 0.06 | 0.06 |
| (1+1) EA | 1 | 0.94 | 0.74 | 0.28 | 0.22 | 0.14 | 0.06 | 0 | 0 | 0 |
| GA | 1 | 0.98 | 0.82 | 0.36 | 0.28 | 0.36 | 0.28 | 0.10 | 0.10 | 0.10 |
| simulated annealing | 1 | 1 | 1 | 1 | 1 | 0.98 | 0.90 | 0.64 | 0.26 | 0.02 |
| tabu search | 1 | 0.86 | 0.54 | 0.46 | 0.40 | 0.62 | 0.54 | 0.24 | 0.48 | 0 |
| compact GA | 1 | 0.92 | 0.54 | 0.20 | 0.12 | 0.12 | 0.10 | 0.04 | 0.02 | 0 |
| UMDA | 1 | 1 | 0.48 | 0.24 | 0.14 | 0.16 | 0.04 | 0.02 | 0 | 0.02 |

- **At k = 0** every search but random search finds the optimum within 5 000 evaluations.
- **The searches that follow the landscape's neighbourhood** do best from k = 2 on. These are simulated annealing, which walks it with the Metropolis rule, and tabu search, which walks it deterministically. At 5 000 evaluations tabu search leads from k = 4 to 10, and its median best point is the optimum (rank 0) at k = 4, 5 and 6. At 50 000, annealing on the Lam schedule finds the optimum in every run up to k = 4 and in 64% of runs at k = 8.
- **Tabu search does not profit from more evaluations.** Its success at 50 000 is no higher than at 5 000 up to k = 6: 0.54 and 0.54 at k = 2, 0.54 and 0.40 at k = 4, 0.52 and 0.54 at k = 6 (the two budgets use different seeds). Per instance (`results/instances.csv`), success varies from seed to seed on the same landscape, for example 4, 5, 7, 2 and 9 of 10 on the five instances at k = 2 with 5 000 evaluations. The random start decides it: a run with tenure 5 and no restarts that is drawn into a cycle of moves stays in it.
- **The (1+1) EA, the compact GA and UMDA** fall off fastest. Without restarts they settle in a local optimum, and from k = 4 on ten times the evaluations gains them little (at k = 4, from 0.10 to 0.22, 0.04 to 0.12 and 0.08 to 0.14). UMDA's and the compact GA's product distributions model no dependence between bits, which is exactly what k adds.
- **The GA** keeps a population. Its median gap at 5 000 evaluations is smaller than the (1+1) EA's at every k ≥ 2, and larger than either local walk's at every k from 3 to 10.
- **At k = 19**, where the landscape is uncorrelated, no search does clearly better than random search. At 5 000 evaluations the median rank of random search's best point is 124, against an expected median of 2²⁰(1 − 2^(−1/5000)) ≈ 145. The GA's is 129, annealing's 234 and the (1+1) EA's 863. At 50 000 random search finds the optimum in 3 of 50 runs, near the 4.8% expected, and the median ranks are 16 for random search, 13 for the GA and 18 for annealing, differences within the spread of fifty runs. On an uncorrelated landscape a neighbour's value says nothing about a point's, so following neighbours gains nothing, and searches that revisit points lose evaluations: (1 − 1/n)ⁿ ≈ 36% of the (1+1) EA's offspring are copies of their parent, evaluated again. This is the situation of the no-free-lunch theorem, which holds for classes of functions closed under the permutations of the space.

## 7. The Figures

- `search-27-ruggedness.jpg`: local optima, ρ(1), correlation length, fitness–distance correlation, the global optimum's basin and the walk's autocorrelation against k.
- `search-28-performance-against-k.jpg`: the gap, the rank and the success of the seven searches against k, at both budgets.
- `search-29-local-optima-and-basins.jpg`: every local optimum of three landscapes, by fitness and basin size.

Their captions are in [`gallery.md`](gallery.md).

## 8. Caveats

- **One parameter setting each.** The tenure (5), the compact GA's K (50), UMDA's λ and μ, the GA's μ and the annealing schedule's T₀ were set once from the literature's defaults, not tuned per k. A tuned tabu search with restarts, or an annealing schedule fitted to the budget, would do better. The comparison is of these settings.
- **No restarts.** Restarts would help the searches that settle (the (1+1) EA, the compact GA, UMDA, tabu search) more than those that keep moving. heuresis has restart wrappers for CMA-ES only; a general restart wrapper for any `Optimiser` is not yet written.
- **Five instances per k.** The spread between instances is large at small k (2 to 44 local optima at k = 1), so the means at k ≤ 3 are rough.
- **Basins are those of steepest descent.** A first-improvement climber, or the (1+1) EA, has other basins. The local optima network's edges (the probabilities of escaping from one basin to another) are not computed: heuresis's `landscapes` lists them as an open question, so search-29 shows the optima and their basins without the network.

## 9. References

- Jones, T. and Forrest, S. (1995), Fitness distance correlation as a measure of problem difficulty for genetic algorithms, *ICGA 6*, 184–192.
- Kauffman, S. A. and Levin, S. (1987), Towards a general theory of adaptive walks on rugged landscapes, *J. Theor. Biol.* 128, 11–45.
- Kauffman, S. A. (1993), *The Origins of Order*, Oxford University Press.
- Ochoa, G., Tomassini, M., Vérel, S. and Darabos, C. (2008), A study of NK landscapes' basins and local optima networks, *GECCO*, 555–562.
- Stadler, P. F. (1996), Landscapes and their correlation functions, *J. Math. Chem.* 20, 1–45.
- Weinberger, E. D. (1990), Correlated and uncorrelated fitness landscapes and how to tell the difference, *Biol. Cybern.* 63, 325–336.
- Weinberger, E. D. (1991), Local properties of Kauffman's N-k model: a tunably rugged energy landscape, *Phys. Rev. A* 44, 6399–6413.
- Wolpert, D. H. and Macready, W. G. (1997), No free lunch theorems for optimization, *IEEE Trans. Evol. Comput.* 1, 67–82.
