# Benchmarking Continuous Search

*Status: draft (2026-11-05).* Research and use: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).

## 1. The Question

Which of heuresis's continuous searches should a domain reach for, on what kind of problem, in what dimension? And where does each fail?

A comparison of search algorithms means something only on a stated class of problems, under stated budgets and targets: by the No Free Lunch theorems no algorithm is better than another averaged over all problems, so a result is a statement about the structures a class has and an algorithm exploits (heuresis's theory note, §2). The COCO platform made such comparisons reproducible for continuous black-box search: functions with known optima, instances that move and turn them so that no algorithm profits from an optimum at the centre or from the coordinate axes, many targets per run, and runtimes summarised as empirical distributions and expected running times (Hansen, Auger, Ros, Mersmann, Tušar and Brockhoff 2021). This case study applies that method, on a smaller suite, to ten algorithms of heuresis at their defaults. It gives the domains a first map of which search to use, and it checks heuresis's algorithms against what the literature leads one to expect of them.

## 2. The Method

**Functions.** Eight functions of `heuresis::benchmarks`, each minimised on the box [−5, 5]ⁿ:

| Function | Instance | Structure |
|---|---|---|
| sphere | f(R(x − s)) | unimodal, isotropic |
| ellipsoid | f(x − s) | unimodal, separable, condition 10⁶ along the axes |
| rotated ellipsoid | f(R(x − s)) | the same, its axes turned |
| Rosenbrock | f(R(x − s) + 1) | a curved valley (its minimiser (1, …, 1) moved to s) |
| Rastrigin | f(R(x − s)) | 10ⁿ-fold multimodal on a global bowl, rotated |
| Ackley | f(R(x − s)) | multimodal, a deep funnel on a nearly flat plateau |
| Griewank | f(100 R(x − s)) | multimodal, at its classical scale [−500, 500]ⁿ |
| Schwefel | f(100 σ ⊙ x) + 10⁴ Σ max(0, \|xᵢ\| − 5)² | separable and deceptive: the second-best minimum of each coordinate lies far from the best |

The shift s is uniform in [−4, 4]ⁿ and the rotation R is drawn from the Haar measure on O(n) (`benchmarks::Rotation`, the Q of a Gaussian matrix's QR factorisation, signed by R's diagonal), both keyed by (function, instance), so that an instance is the same wherever it is made. The separable ellipsoid is shifted but not rotated, so that the pair of ellipsoids measures what an algorithm gains from separability. Schwefel's function keeps its separability, which is its structure: its coordinates are scaled by 100 into the box and each sign flipped with probability ½, so its minimiser is ±4.2096874… in each coordinate. Outside [−500, 500] Schwefel's function falls below f* without bound, so, as in COCO's f20, a penalty 10⁴ Σ max(0, |xᵢ| − 5)² is added outside the box; it is 0 inside and exceeds anything the function can gain outside (§8 says what happened without it). Its value at the minimiser is n times the minimum per coordinate, found by golden sections; for the other functions f* = 0. Before any run the program checks that every instance's minimiser has f(x*) − f* below 10⁻⁹.

**Dimensions, instances, budget, targets.** n ∈ {2, 5, 10, 20}; 15 instances per function and dimension; a budget of 10⁴·n evaluations per run; eleven targets f − f* ∈ {10², 10¹, …, 10⁻⁸}. A run ends at f − f* ≤ 10⁻⁸, at the budget, or when the algorithm stops by its own criteria. Each algorithm runs once on each instance: 8 × 4 × 15 × 10 = 4800 runs.

**Driving and measuring.** Every algorithm is driven by hand through the ask–tell protocol, the best point and its improvements kept by `search::Elite` as `search::run` keeps them. A run's progress is read by `heuresis::analysis`:
- its runtime to a target, the evaluations at its first improvement at or below it (`analysis::runtime`);
- the empirical distribution of runtimes over (run, target) pairs, COCO's anytime measure (`analysis::runtime_distribution`): for each budget, the fraction of pairs solved;
- the expected running time to a target, ERT = (evaluations spent over all runs until success, or in all when a run failed) / (successes), the expected cost of independent restarts until one succeeds (`analysis::expected_running_time`), with the success rate and Clopper and Pearson's 95% interval (`analysis::clopper_pearson`).

The runs are spread over eight threads; each is a pure function of its seed, so the threads change nothing. Each run's cost is measured as its thread's processor time, which does not depend on what else the machine runs.

## 3. The Code

- heuresis `6b3cc88` (`git -C foundations/heuresis rev-parse --short HEAD`), attached as a submodule of the root; calx, morphe, physis and tyche at the commits the root records.
- The root repository at `4153cb9`.
- [`benchmark.cpp`](benchmark.cpp), built by [`CMakeLists.txt`](CMakeLists.txt) as a project of its own on the foundations (C++20, Release, g++ 11.5 in WSL), and [`figures.py`](figures.py).

```
cmake -S . -B ~/.cache/science-build/exp-benchmarking-continuous-search -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-benchmarking-continuous-search -j4
~/.cache/science-build/exp-benchmarking-continuous-search/benchmark results
.venv/Scripts/python.exe research/experiments/2026-11-05-benchmarking-continuous-search/figures.py
```

## 4. The Parameters

Every algorithm runs at heuresis's defaults unless stated. The (1+1)-ES, CMA-ES and the cross-entropy method start at a point uniform in [−4, 4]ⁿ (the same point for the three on a given instance) with σ₀ = 2, and search ℝⁿ without bounds; IPOP and BIPOP draw each restart's mean uniformly in [−4, 4]ⁿ; the population methods start uniform in [−5, 5]ⁿ and stay in it.

| Algorithm | Construction |
|---|---|
| random | `search::RandomSearch`, uniform in [−5, 5]ⁿ |
| (1+1)-ES | `strategies::OnePlusOneEs`, the smooth one-fifth rule; stopped once σ < 10⁻¹² (`min_sigma`; the default, 0, never stops) |
| CMA-ES | `strategies::Cma`, the tutorial's parameters (λ = 4 + ⌊3 ln n⌋), active weights, its own stopping criteria (TolX, TolFun, ConditionCov, NoEffectAxis, NoEffectCoord, Stagnation, TolUpSigma) |
| IPOP-CMA-ES | `strategies::CmaRestarts`, IPOP: λ doubled at each restart, up to 30 runs (the budget ends every run first) |
| BIPOP-CMA-ES | `strategies::CmaRestarts`, BIPOP: IPOP's regime interleaved with short runs of small populations and step sizes σ₀·10^(−2U) |
| DE | `strategies::DifferentialEvolution`, rand/1/bin, F = 0.5, CR = 0.9, NP = 10n, midpoint repair at the bounds |
| SHADE | `strategies::Shade`, NP = 100, a memory of H = 100, p ∈ [2/NP, 0.2], an external archive |
| PSO | `swarms::ParticleSwarm`, the inertia rule (w = 0.7298, c₁ = c₂ = 1.49618, Clerc and Kennedy's constriction) on a ring of 40 particles |
| SPSO 2011 | `swarms::ParticleSwarm` at its defaults: the Standard PSO 2011, 40 particles, adaptive random informants |
| CE | `distributions::DistributionSearch` with `Gaussian` and `CrossEntropy`: λ = max(50, 10n), elite max(10, λ/5), α = 0.7, dynamic smoothing β = 0.8, q = 5 |

Two particle swarms are compared because heuresis's `swarms` README predicts that they differ: SPSO 2011's move is invariant under rotations, but its particles' second moments diverge from three dimensions on, while the inertia rule converges and depends on the coordinates.

## 5. The Data

No external data. Everything in `results/` is written by `benchmark.cpp`:

| File | Content |
|---|---|
| `runs.csv` | one line per run (4800): evaluations spent, the gap f − f* reached, how it ended (CMA-ES with its stopping criterion), the runtime to each of the eleven targets |
| `ecdf.csv` | the runtime distributions at budgets 10^(k/20)·n, k = 0, …, 80: per dimension over all functions (`all`) and per function |
| `ert.csv` | ERT, successes, success rate and its 95% interval, per function, dimension, algorithm and target |
| `timing.csv` | processor seconds and evaluations spent per algorithm and dimension |

## 6. The Results

All numbers are for the final target f − f* ≤ 10⁻⁸ unless stated; "solved" means an instance reached it within 10⁴·n evaluations. ERT is given per dimension, ERT/n.

**Overall.** Of the 120 instances in each dimension (8 functions × 15), the algorithms solved:

| Algorithm | n = 2 | n = 5 | n = 10 | n = 20 | ECDF at the budget, n = 2 / 5 / 10 / 20 |
|---|---|---|---|---|---|
| random | 0 | 0 | 0 | 0 | 0.317 / 0.102 / 0.068 / 0.035 |
| (1+1)-ES | 42 | 22 | 18 | 17 | 0.457 / 0.327 / 0.297 / 0.252 |
| CMA-ES | 77 | 76 | 81 | 87 | 0.714 / 0.696 / 0.718 / 0.745 |
| IPOP-CMA-ES | 113 | 106 | 104 | 95 | 0.958 / 0.891 / 0.874 / 0.829 |
| BIPOP-CMA-ES | 110 | 101 | 93 | 92 | 0.951 / 0.856 / 0.800 / 0.789 |
| DE | 106 | 89 | 82 | 31 | 0.924 / 0.828 / 0.783 / 0.501 |
| SHADE | 109 | 90 | 90 | 102 | 0.985 / 0.832 / 0.831 / 0.867 |
| PSO (inertia, ring) | 100 | 50 | 45 | 53 | 0.914 / 0.548 / 0.468 / 0.486 |
| SPSO 2011 | 89 | 25 | 0 | 0 | 0.870 / 0.383 / 0.167 / 0.064 |
| CE | 23 | 1 | 0 | 0 | 0.648 / 0.343 / 0.219 / 0.129 |

The ECDF column is the fraction of the 1320 (run, target) pairs solved at the budget. Which algorithm leads depends on the budget (search-06):
- up to 100·n evaluations, the (1+1)-ES, which takes the easy targets fastest (at 100·n it has solved 0.21 to 0.32 of the pairs, against 0.21 to 0.31 for CMA-ES);
- from about 300·n to 3000·n, the CMA-ES family in five, ten and twenty dimensions (IPOP 0.76 at 10³·n in twenty), and DE in two (0.91 at 10³·n);
- at the full budget, IPOP-CMA-ES in five and ten dimensions and SHADE in two and twenty. Where both succeed, plain CMA-ES needs 2.2 to 21 times fewer evaluations than SHADE (ERT/n on the sphere, the two ellipsoids, Rosenbrock's and Ackley's functions); SHADE's lead at the budget comes mostly from Schwefel's function.

**Unimodal functions.**
- The (1+1)-ES is the fastest on the sphere at every dimension (ERT/n 92, 83, 81, 79), faster than CMA-ES (131, 149, 152, 142). On either ellipsoid it fails in every dimension: its isotropic step must shrink to the scale of the narrow axes, so it hardly moves along the wide ones (median gaps 26.5, 32.8, 94 and 266 for n = 2, 5, 10, 20 on the separable ellipsoid).
- CMA-ES solves both ellipsoids at the same cost, as its invariance under rotations requires: ERT/n 255, 311, 428, 679 on the separable one and 253, 318, 417, 670 on the rotated one. Its restarts cost nothing here (IPOP: 267, 309, 426, 666).
- Algorithms whose operators work coordinate by coordinate exploit separability and lose it under a rotation. DE (CR = 0.9) solves the separable ellipsoid at every dimension (ERT/n 517 to 8085) but the rotated one only to n = 10, at 2.4 times the cost there (7380 against 3029). The inertia PSO solves the separable ellipsoid at every dimension and the rotated one in 2 of 15 instances at n = 2 and none beyond (median gaps 14, 297 and 1850 at n = 5, 10, 20). SHADE solves both at every dimension, the rotated one at 1.6 to 4.4 times the cost.
- Rosenbrock's function is solved by the CMA-ES family in every dimension (plain CMA-ES 14 of 15 at n = 20, ERT/n 944) and by SHADE in 12 of 15 at n = 20; DE solves 6 and 7 of 15 at n = 5 and 10 and none at 20.

**Multimodal functions.**
- Rotated Rastrigin separates the restart strategies. IPOP-CMA-ES solves 15, 15, 14 and 5 of 15 instances (ERT/n 1545, 4060, 7446, 29028); BIPOP, which spends about half its budget on small populations, 15, 12, 3 and 2; plain CMA-ES 1, 1, 0, 0. SHADE solves it only at n = 2, with median gaps 0.14, 3.0 and 10.7 at n = 5, 10, 20.
- Ackley's function is solved in every dimension by the CMA-ES family and SHADE, and by DE up to n = 10.
- Griewank's function at its classical scale becomes easier with dimension: plain CMA-ES solves 0, 0, 6 and 13 of 15, IPOP 12, 13, 15 and 15, SHADE 4, 0, 0 and 15.
- Schwefel's separable, deceptive function is SHADE's: it solves all 15 instances in every dimension (ERT/n 3747 to 4591), and DE all 15 at n = 5 and 10 and none at 20. The searches by rotation-invariant Gaussians fail from n = 5 on: IPOP solves 11 and 3 of 15 at n = 2 and 5, plain CMA-ES 1 at n = 2, and their median gaps at n = 20 are 2170 and 3530. The structure of Schwefel's function is its separability, and the algorithms that cross coordinates over exploit it.

**The two swarms.** The inertia rule on a ring solves 100, 50, 45 and 53 instances; SPSO 2011 solves 89 and 25 at n = 2 and 5 and none from n = 10 on, not even the sphere (median gaps 5.8 × 10⁻³ at n = 10, 2.5 at n = 20). That is what heuresis's `swarms` module predicts: SPSO 2011's particles' second moments grow from three dimensions on, so the swarm progresses only through its personal bests.

**The cross-entropy method** converges on the sphere, but slowly: 13 of 15 instances at n = 2, 1 at n = 5, with median gaps 2.2 × 10⁻⁷ at n = 5 and 7.0 × 10⁻³ at n = 10. Dynamic smoothing makes its covariance decay polynomially, so its convergence is sublinear by construction: it trades speed for not stalling.

**Endings.** Plain CMA-ES stopped by its own criteria in 159 runs (155 by TolFun, 4 by Stagnation), the (1+1)-ES in 213 (its step below 10⁻¹²). No run of any algorithm raised an error.

**Cost.** The 4800 runs took 869 processor-seconds, 397 s of wall time on eight threads of a machine shared with other work. Per evaluation, from 0.6 µs (random search at n = 2, mostly the function) to 10.4 µs (CMA-ES at n = 20, whose update is O(n²) per sample and whose eigendecomposition, at the tutorial's interval, runs every generation there).

## 7. The Figures

| Image | Content |
|---|---|
| [`search-06-ecdf-by-dimension.jpg`](../../../assets/gallery/search-06-ecdf-by-dimension.jpg) | runtime distributions over all eight functions, one panel per dimension |
| [`search-07-ecdf-by-function.jpg`](../../../assets/gallery/search-07-ecdf-by-function.jpg) | runtime distributions on each function at n = 10 |
| [`search-08-ert.jpg`](../../../assets/gallery/search-08-ert.jpg) | ERT/n to f − f* ≤ 10⁻⁸ against dimension, one panel per function; filled points where every instance was solved |
| [`search-09-solved.jpg`](../../../assets/gallery/search-09-solved.jpg) | for each algorithm, function and dimension, the instances solved, or the median gap reached where none was |

## 8. The Caveats

- **A small suite, not COCO's.** Eight functions with shifts and rotations, without BBOB's asymmetry and oscillation transformations (T_asy, T_osz) or its 24 functions; the instances are heuresis's, not COCO's, so the numbers are not comparable with COCO's archive. heuresis's `benchmarks` README lists BBOB's functions as an open question.
- **Fifteen runs per cell.** Each algorithm runs once on each of 15 instances. A success count of 10 of 15 has a 95% interval of [0.38, 0.88] (`ert.csv` gives each); ERTs have no interval, since bootstrapped intervals of ERT are planned in `analysis`, not written.
- **Defaults, not tuning.** Every algorithm runs at heuresis's defaults: DE's CR = 0.9 does not favour separable functions, as a CR near 0.1 would; SHADE's population of 100 is costly in two dimensions; the cross-entropy method runs at the parameters of the library's own check.
- **The budget.** 10⁴·n evaluations is small by COCO's standards. IPOP's success on Rastrigin's function at n = 20 (5 of 15) would grow with a larger budget; the plain strategies that stop early (CMA-ES, the (1+1)-ES) have ERTs that already count their restarts.
- **Bounds.** The (1+1)-ES, CMA-ES, IPOP, BIPOP and the cross-entropy method search ℝⁿ; the population methods are confined to [−5, 5]ⁿ, which contains every optimum. A first run of this study gave Schwefel's function without the penalty outside the box: the searches of ℝⁿ then "solved" it in 12 to 15 of 15 instances by leaving the box, to gaps as low as −3660, where Schwefel's function falls below f* without bound. The penalty, as in COCO's f20, removes those points; the record keeps only the corrected run.
- **Timing.** Processor time per thread does not depend on the other work on the machine, but shared caches and memory bandwidth do; the costs per evaluation are indicative.
