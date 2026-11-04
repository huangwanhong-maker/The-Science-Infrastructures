# The Genetic Algorithm as a Population

*Status: draft (2026-11-05).* Learning and research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) beside the evolution domain's population genetics. Purely mathematical models: no organism, no sequence.

## 1. The Question

Evolutionary computation borrowed its words from population genetics: population, selection, drift, mutation, fixation. The borrowing is exact in one case. A generational genetic algorithm without crossover, whose parents are drawn by fitness-proportional selection with replacement and whose children copy them, is the haploid Wright–Fisher model with selection (Fisher 1930; Wright 1931). If k of its M individuals carry a variant of fitness 1 + s and the rest have fitness 1, each child's parent is a carrier with probability

  ψ = k(1 + s) / (k(1 + s) + M − k),

independently of the other children, so the carriers among the children are Binomial(M, ψ): one generation of Wright–Fisher with genic selection. The theory of the one then holds for the other: Kimura's (1962) diffusion gives the probability that one copy of the variant fixes,

  u = (1 − e^{−2s}) / (1 − e^{−2Ms}),

and without selection the expected heterozygosity 2p(1 − p) of any label decays by 1 − 1/M a generation.

The bridge is useful in both directions. It tells the designer of a search how much of its population's loss of diversity is drift rather than selection, and how that depends on the scheme. It lets population genetics be learnt and checked with the operators of a search library. This study shows the bridge by simulation with heuresis's operators, against the exact Markov chain and the evolution domain's Wright–Fisher replicates and diffusion formulas, and follows it to where it ends: other selection schemes are not Wright–Fisher, and their drift has a rate of its own.

## 2. The Method

Every GA here is built from heuresis's parts: parents by `selection::roulette` (one draw proportional to nonnegative weights), `selection::universal_sampling` (Baker's stochastic universal sampling), `selection::tournament` (the best of k drawn with replacement, ties to the earlier index) or `selection::truncation` (the best μ, then a parent uniform among them); children copied; for OneMax, `variation::standard_bit_mutation` and `benchmarks::onemax`, and for tournaments on OneMax the library's own `genetic::GeneticAlgorithm`. Every draw comes from a stream keyed by what it is for (`search::Draws`: a purpose, the generation, the replicate). Replicates run in parallel in fixed chunks on a calx thread pool, summed in order, so every number is the same on any number of threads.

**One generation's law.** M = 40, k = 10 carriers of fitness 1.5, 200,000 independent generations for each of roulette, universal sampling and a tournament of two. The theory: roulette Binomial(40, 1/3); a tournament of two Binomial(40, 1 − (1 − k/M)²) = Binomial(40, 7/16), since a child's parent is a carrier unless both contestants are not; universal sampling ⌊Mψ⌋ or ⌈Mψ⌉ copies, the latter with probability Mψ − ⌊Mψ⌋.

**Fixation.** One copy of the variant among M ∈ {20, 50, 100, 200}, at scaled selection coefficients S = 2Ms ∈ {−2, −1, 0, 1, 2, 4, 8, 16}, run under roulette until it is lost or fixed: 20,000 runs a point, 40,000 for S < 0. Compared with:
- the exact fixation probability of the haploid Wright–Fisher chain, u_i = Σ_j P_ij u_j on the transient states with P_ij = Binomial(M, ψ_i)(j), solved by calx's dense Gaussian elimination;
- Kimura's diffusion from the evolution domain, `evolution::population::kimura_fixation(N, s, p)` with N = M/2 diploids (2N = M gametes), the allele's advantage s and p = 1/M, which is the formula above; and `fixation_probability` (the diffusion with dominance, by quadrature) at h = ½, which agrees with it to the digits written;
- the domain's own Wright–Fisher replicates, `wright_fisher_replicates`, with N = M/2 diploids, genic selection (fitnesses 1 + 2s, 1 + s, 1) and p₀ = 1/M, the same number of runs.

**Drift.** A neutral label carried by half of M = 50 individuals, at positions drawn uniformly, under equal fitness, 4000 runs of 300 generations for each of five schemes. For a scheme whose children choose their parents independently with probabilities q_j (given the population's order), and labels arranged exchangeably, the expected heterozygosity obeys

  E[2p′(1 − p′)] = (1 − c) · 2p(1 − p),  c = Σ_j q_j²,

the probability that two children share a parent. (Proof: given the arrangement, p′ is a mean of M independent Bernoulli variables of mean p_q = Σ q_j g_j, so E[p′(1 − p′)] = (1 − 1/M) p_q(1 − p_q); averaging p_q over the exchangeable arrangement of the labels gives (1 − Σ q_j²)·p(1 − p).) So:
- roulette with equal weights: q_j = 1/M, c = 1/M, Wright–Fisher's rate;
- universal sampling with equal weights: every individual exactly one copy, c = 0, no drift at all;
- a tournament of two with ties to the earlier index (the library's rule): under equal fitness the earlier contestant always wins, so q_j = (2(M − j) − 1)/M² for j = 0 … M − 1, and c = (4M² − 1)/(3M³) ≈ 4/(3M);
- the same tournament with ties broken by keys drawn afresh for each tournament (the caller's order): q_j = 1/M, c = 1/M;
- truncation to the best M/2 (by index, under equal fitness), a parent uniform among them: c = 2/M.

The domain's neutral replicates (N = 25 diploids, p₀ = ½) are run beside them.

**Takeover.** One best individual (fitness 1 + s) among M = 500 others of fitness 1, under tournaments of two and four, truncation to M/2 and M/4, and roulette at s = 1 and s = 4, 1000 runs each, without mutation or elites. Each scheme's runs at s = 1 and s = 4 use the same streams. The time to fill the population, among the runs in which the best fixed, is compared with the deterministic recursion of its share p from 1/M (tournament: 1 − (1 − p)^k; truncation to τM: min(1, p/τ); roulette: p(1 + s)/(1 + sp); its takeover the first generation with fewer than half an individual left), and with Goldberg and Deb's (1991) closed forms: (ln M + ln ln M)/ln k for tournaments, ln M/ln(1/τ) for truncation, and, from odds 1/(M − 1) to M − 1, 2 ln(M − 1)/ln(1 + s) for roulette between two values.

**Mutation–selection balance on OneMax.** Strings of n = 100 bits from uniform random strings, standard bit mutation at rate χ/n with χ ∈ {¼, ½, 1, 2, 4, 8}, generational replacement without elites or crossover, M ∈ {20, 100, 500}, four seeds, 2000 generations; the equilibrium is the mean share of ones over generations 1000–2000. Two schemes:
- roulette with weights the number of ones, built from heuresis's parts;
- the library's `genetic::GeneticAlgorithm` with tournaments of two, crossover rate 0, no elites, minimising OneMax's zeros.

They are compared with the infinite population on the classes 0 … n of the number of ones. Selection (roulette: x_k ← k x_k/Σ_j j x_j; a tournament of two: the winner's class at least c with probability 1 − (1 − G(c))², G(c) = Σ_{j≥c} x_j) is followed by mutation (from class k, a ~ Bin(k, μ) ones lost and b ~ Bin(n − k, μ) gained), iterated from the binomial law of random strings to its fixed point. Because fitness and mutation depend on a string only through its class, this is exact for an infinite asexual population: Eigen's (1971) quasispecies on OneMax's classes. For roulette, a closed form follows if the loci stay independent, so that the variance of the number of ones is binomial: Fisher's selection differential Var/k̄ = n(1 − f) balances mutation's loss χ(2f − 1) when f = (1 + χ)/(1 + 2χ).

## 3. The Code

- `ga_population.cpp`, built by `CMakeLists.txt` against calx, morphe, physis, tyche, heuresis, the genomics domain and the evolution domain by relative path. Its modes `law`, `drift`, `takeover`, `balance` and `fixation` write `results/`.
- `figures.py` draws the three figures from `results/` alone.
- Versions: heuresis 6b3cc88, The Science Infrastructures 4153cb9 (the evolution domain is part of it). The runs were made on 2026-11-05 with g++ 11.5 (`-O3` by CMake's Release) in WSL on the development machine, eight threads, shared at the time with other work: law 2 s, drift 63 s, takeover 10 s, balance 121 s, fixation 246 s.

```
cmake -S . -B ~/.cache/science-build/exp-ga-population -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-ga-population --target ga_population -j4
~/.cache/science-build/exp-ga-population/ga_population results all
../../../.venv/Scripts/python.exe figures.py
```

## 4. The Parameters

| Part | Parameters |
|---|---|
| One generation | M 40, k 10, s 0.5; 200,000 generations a scheme; draws seed 11 |
| Fixation | M 20, 50, 100, 200; S = 2Ms from −2 to 16; 20,000 runs (40,000 for S < 0); draws seed 1000 + 100·(M's index) + (S's index); the domain's replicates seed 5000 + the same |
| Drift | M 50, p₀ ½ at uniform positions, 300 generations, 4000 runs a scheme; seeds 21 (arrangement) and 22 (selection); the domain's replicates seed 23 |
| Takeover | M 500, one best, s 1 and 4, 1000 runs, at most 200 generations; seed 31 |
| Balance | n 100, χ ¼ to 8, M 20, 100, 500, 4 seeds, 2000 generations, equilibrium over 1000–2000 |

## 5. The Data

None. Every number comes from the models above.

## 6. The Results

### 6.1 One Generation Is Wright–Fisher's

| Scheme | Mean (exact) | Variance (exact) | Total variation from the exact law |
|---|---|---|---|
| roulette | 13.3234 (13.3333) | 8.8460 (8.8889) | 0.0047 |
| tournament of two | 17.5034 (17.5000) | 9.8101 (9.8438) | 0.0033 |
| universal sampling | 13.3349 (13.3333) | 0.2228 (0.2222) | 0.0016 |

Roulette's carriers are Binomial(M, ψ), the tournament's Binomial(M, 1 − (1 − p)²): both are Wright–Fisher generations, the second with a selection of its own. Universal sampling is not: it gives 13 or 14 copies and nothing else, its variance 40 times smaller.

### 6.2 Fixation

The GA's fixation probabilities agree with the exact chain within two standard errors at 30 of the 32 points, and the sum of their squared standardised deviations is 28.2 on 32 points, as sampling error alone would give. The two outside are M = 50, S = 16 (−2.05) and M = 200, S = 16 (−2.57).

| M | S | GA | exact chain | Kimura's diffusion | the domain's replicates |
|---|---|---|---|---|---|
| 20 | 0 | 0.0505 | 0.0500 | 0.0500 | 0.0518 |
| 20 | 8 | 0.3168 | 0.3113 | 0.3298 | 0.3083 |
| 20 | 16 | 0.5060 | 0.5066 | 0.5507 | 0.4925 |
| 100 | −2 | 0.00332 | 0.00314 | 0.00316 | 0.00343 |
| 100 | 4 | 0.0388 | 0.0397 | 0.0399 | 0.0378 |
| 100 | 16 | 0.1408 | 0.1442 | 0.1479 | 0.1451 |
| 200 | 0 | 0.00510 | 0.00500 | 0.00500 | 0.00460 |
| 200 | 16 | 0.0711 | 0.0759 | 0.0769 | 0.0744 |

(`results/fixation.csv` holds all 32 points.) Kimura's diffusion is the limit s → 0 at fixed S, and its error shows where s is not small: it exceeds the exact chain by 8.7% at M = 20, S = 16 (s = 0.4), by 4.7% at M = 50, S = 16, and by at most 1.4% wherever s ≤ 0.04. In the scaled form M·u against S, the points of every M collapse onto the diffusion's limit S/(1 − e^{−S}) as M grows (search-44). The domain's diploid replicates agree with the haploid exact chain within 2.8 standard errors at every point but M = 20, S = 16 (−4.0), where s = 0.4 and diploid genic selection is no longer the haploid's.

### 6.3 Drift, by Scheme

| Scheme | c in theory | c measured over 100 generations | E[H] at 100: measured (theory) |
|---|---|---|---|
| universal sampling | 0 | 0 | 0.5000 (0.5000) |
| roulette | 1/M = 0.0200 | 0.0198 | 0.0679 (0.0663) |
| tournament of two, fresh ties | 1/M = 0.0200 | 0.0200 | 0.0666 (0.0663) |
| tournament of two, ties to the earlier index | (4M² − 1)/3M³ = 0.0267 | 0.0274 | 0.0311 (0.0335) |
| truncation to M/2 | 2/M = 0.0400 | 0.0389 | 0.0095 (0.0084) |
| the domain's Wright–Fisher, 2N = 50 | 1/(2N) = 0.0200 | 0.0194 | 0.0704 (0.0663) |

Roulette and the domain's Wright–Fisher drift at the same rate, 1 − 1/M. Universal sampling under equal weights does not drift at all. A tournament with the library's ties drifts a third faster than Wright–Fisher, an effective size of 3M/4: under equal fitness the order of the population becomes a ranking, and the earlier individuals have more children. Drawing fresh keys for each tournament restores Wright–Fisher's rate exactly. By the same formula, a random ranking drawn once a generation would not: its c is the index's.

### 6.4 Takeover

| Scheme | Fixed of 1000 | Takeover (generations), mean ± sd | Deterministic recursion | Goldberg and Deb |
|---|---|---|---|---|
| truncation to M/4 | 983 | 5.13 ± 0.53 | 5 | 4.48 |
| tournament of four | 974 | 6.46 ± 0.60 | 6 | 5.80 |
| roulette, s = 4 | 993 | 8.61 ± 0.90 | 9 | 7.72 |
| truncation to M/2 | 800 | 9.67 ± 1.48 | 9 | 8.97 |
| tournament of two | 809 | 12.39 ± 1.53 | 12 | 11.60 |
| roulette, s = 1 | 780 | 18.79 ± 2.38 | 19 | 17.93 |

The deterministic recursion predicts each mean within a generation; Goldberg and Deb's closed forms lie 0.6 to 0.9 generations below. Without elites the best is lost by drift in 0.7% (roulette, s = 4) to 22% (roulette, s = 1) of runs. Tournaments and truncation read fitness only through its order, so their runs at s = 1 and s = 4 are identical, count for count (`takeover_invariance.json`), as the theory note's Proposition 3.3 requires; roulette's time falls from 18.8 to 8.6 generations between s = 1 and s = 4.

### 6.5 Mutation–Selection Balance on OneMax

| χ | Roulette: M 20, 100, 500 | infinite | (1 + χ)/(1 + 2χ) | Tournament of two: M 20, 100, 500 | infinite |
|---|---|---|---|---|---|
| ¼ | 0.618, 0.709, 0.755 | 0.825 | 0.833 | 0.996, 0.997, 0.997 | 0.997 |
| ½ | 0.602, 0.664, 0.696 | 0.738 | 0.750 | 0.978, 0.990, 0.990 | 0.990 |
| 1 | 0.578, 0.614, 0.639 | 0.655 | 0.667 | 0.874, 0.912, 0.925 | 0.930 |
| 2 | 0.555, 0.581, 0.586 | 0.592 | 0.600 | 0.765, 0.795, 0.803 | 0.808 |
| 4 | 0.539, 0.546, 0.550 | 0.550 | 0.556 | 0.673, 0.690, 0.695 | 0.697 |
| 8 | 0.522, 0.524, 0.524 | 0.525 | 0.529 | 0.602, 0.610, 0.611 | 0.612 |

Roulette on OneMax selects weakly, since the numbers of ones differ little relative to their mean, and settles far from the optimum: 0.655 ones at χ = 1 even in an infinite population, 1.7% below the independent-loci formula's 2/3. A tournament of two, which reads only the order, holds 93% ones at χ = 1. Finite populations fall short of the infinite one, by more where selection is weaker: under roulette at χ = ¼, M = 500 reaches 0.755 against the infinite population's 0.825, and M = 20 only 0.618. Selection acts through the population's variance (Fisher's theorem), which drift erodes in proportion to 1/M; under tournaments the shortfall is at most 0.056 (M = 20, χ = 1).

## 7. The Figures

| Figure | What it shows |
|---|---|
| [search-44](../../../assets/gallery/search-44-ga-as-wright-fisher.jpg) | One generation's law under three schemes; the fixation probability of one copy, M·u against S, against the exact chain and the diffusion's limit |
| [search-45](../../../assets/gallery/search-45-drift-and-takeover.jpg) | The heterozygosity of a neutral label under five schemes and the domain's Wright–Fisher, against (1 − c)^t; takeover from one copy among 500 |
| [search-46](../../../assets/gallery/search-46-mutation-selection-balance.jpg) | The equilibrium share of ones on OneMax against χ, for roulette and tournaments at three sizes, against the infinite population; trajectories at M = 100 |

`gallery.md` holds them with their captions.

## 8. The Caveats

- **Haploid and diploid.** The GA is a haploid model; the domain's Wright–Fisher is diploid with genic selection, compared at 2N = M gametes. The two agree to first order in s and differ at s = 0.4.
- **What the bridge needs.** It is exact for generational replacement without elites, crossover or mutation, with roulette on fitnesses read as relative fitnesses. Elites, steady-state replacement, crossover and rank-based schemes each change the process. The drift rates of §6.3 hold for any scheme whose children choose parents independently, given an exchangeable arrangement of the labels.
- **Ties.** The faster drift of tournaments under equal fitness comes from the library's documented rule (ties to the earlier index) and the population's order; in the library's own GA, survivors are sorted best first, so within a class of equal fitness the earlier-born are favoured. On real problems ties are common on plateaus.
- **Takeover without elites** loses the best in up to 22% of runs, so the takeover times are conditional on fixation; the deterministic recursion ignores that loss.
- **Equilibria.** At χ ≤ ½ under roulette the population was still creeping upward (at M = 100, seed 1: 0.712 over generations 1000–1500, 0.716 over 1500–2000), so those equilibria are slight underestimates.
- **Cost.** `selection::roulette` builds the cumulative weights afresh at each draw, so a generation of M parents costs O(M²). That limited the fixation study to M ≤ 200.

## 9. References

- Baker, J. E. (1987), Reducing bias and inefficiency in the selection algorithm, *ICGA 2*, 14–21.
- Eigen, M. (1971), Selforganization of matter and the evolution of biological macromolecules, *Naturwissenschaften* 58, 465–523.
- Ewens, W. J. (2004), *Mathematical Population Genetics I: Theoretical Introduction*, 2nd ed., Springer.
- Fisher, R. A. (1930), *The Genetical Theory of Natural Selection*, Clarendon Press.
- Goldberg, D. E. and Deb, K. (1991), A comparative analysis of selection schemes used in genetic algorithms, *FOGA 1*, 69–93.
- Kimura, M. (1962), On the probability of fixation of mutant genes in a population, *Genetics* 47, 713–719.
- Wright, S. (1931), Evolution in Mendelian populations, *Genetics* 16, 97–159.
