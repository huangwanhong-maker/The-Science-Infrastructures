# Search on a Landscape

*Status: draft (2026-11-05).* Learning and research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).

## 1. The Question

Six algorithms of heuresis search the same two-dimensional landscape: how does each move over it, and what does each find within a given number of evaluations?

The algorithms rest on different ideas of where to look next. An evolution strategy keeps one point and a step size; CMA-ES keeps a Gaussian law and learns its shape; differential evolution moves a population by differences of its own members; a particle swarm flies particles towards remembered bests; the cross-entropy method refits a Gaussian to its elite; Bayesian optimisation models the whole landscape and evaluates where the model expects most. In two dimensions their states can be drawn on the landscape itself, so the ideas become visible: a covariance turning to follow a valley, a population contracting onto a basin, a law collapsing where the optimum is not. The case study is for learning first. Its anytime curves and its record of every state are also a baseline for the later comparisons and for an interactive viewer.

## 2. The Method

Three landscapes of `heuresis::benchmarks` with known minima:

| Function | Box | f* | Minimisers | Character |
|---|---|---|---|---|
| Rosenbrock | [−2, 2] × [−1, 3] | 0 | (1, 1) | a narrow curved valley |
| Rastrigin | [−5.12, 5.12]² | 0 | (0, 0) | a lattice of 120 local minima on a bowl |
| Branin | [−5, 10] × [0, 15] | 0.397887357729738 | (−π, 12.275), (π, 2.275), (3π, 2.475) | three equal minima |

Each algorithm is driven by hand through the ask–tell protocol (`search/protocol.hpp`): ask a generation, evaluate it, tell it. Before each ask the program reads the state the generation is drawn from (the parent and σ, CMA-ES's m, σ and C, the cross-entropy method's m and C); after each tell it reads the state the generation left (DE's population after selection, PSO's velocities, the Bayesian surrogate). The best point and its improvements are kept by `search::Elite`, as `search::run` keeps them, so that `heuresis::analysis` reads these runs as it reads any outcome. Every random draw is keyed by its purpose through `search::Draws(seed)`, so each run is a pure function of its seed and the 651 runs are spread over eight threads without changing any of them.

Each (function, algorithm) pair runs with seeds 0 to 30. A run ends when f − f* ≤ 10⁻¹², when its budget is spent, or when the algorithm stops by its own criteria. From the runs the program writes:
- the anytime curves: the best f − f* after b evaluations, its median and quartiles over the 31 seeds (`analysis::anytime_quantile`, quantiles of type 7);
- the spread of each search against its distance from the nearest minimiser, quartiles over the seeds, read from the state after the last generation completed within b evaluations;
- the whole trajectory of seed 0 of every algorithm on every function, for the small multiples and a viewer.

The spread is the root mean square distance of a candidate from the centre of the search: σ√n for the (1+1)-ES, σ√tr C for CMA-ES (its law is N(m, σ²C)), √tr C for the cross-entropy method, and the r.m.s. distance from the centroid of DE's population and of PSO's particles. Bayesian optimisation has no spread in this sense and is left out of that panel.

## 3. The Code

- heuresis `6b3cc88` (`git -C foundations/heuresis rev-parse --short HEAD`), attached as a submodule of the root; calx, morphe, physis and tyche at the commits the root records.
- The root repository at `4153cb9`.
- [`landscape.cpp`](landscape.cpp), built by [`CMakeLists.txt`](CMakeLists.txt) as a project of its own on the foundations (C++20, Release, g++ 11.5 in WSL), and [`figures.py`](figures.py).

```
cmake -S . -B ~/.cache/science-build/exp-search-on-a-landscape -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-search-on-a-landscape -j4
~/.cache/science-build/exp-search-on-a-landscape/landscape results
.venv/Scripts/python.exe research/experiments/2026-11-05-search-on-a-landscape/figures.py
```

The 651 runs took 96 s on eight threads (207 s when the machine was busier with other work), most of it Bayesian optimisation's refits. The results do not depend on the threads: a second run wrote the same numbers.

## 4. The Parameters

Every algorithm uses heuresis's defaults unless stated.

| Algorithm | Construction | Start |
|---|---|---|
| (1+1)-ES | `strategies::OnePlusOneEs`, the one-fifth rule in Kern et al.'s smooth form (target success 1/5, damping 1 + n/2) | the start point, σ₀ |
| CMA-ES | `strategies::Cma`, the tutorial's parameters (λ = 6, μ = 3 in two dimensions), active weights, no bounds | the start point, σ₀ |
| DE | `strategies::DifferentialEvolution`, rand/1/bin, F = 0.5, CR = 0.9, NP = 10n = 20, midpoint repair | uniform in the box |
| PSO | `swarms::ParticleSwarm`, the Standard PSO 2011: 40 particles, adaptive random informants (K = 3), w = 1/(2 ln 2), c = ½ + ln 2 | uniform in the box |
| CE | `distributions::DistributionSearch` with `Gaussian` and `CrossEntropy{elite 10, α = 0.7, β = 0.8, q = 5}`, λ = 50: Kroese, Porotsky and Rubinstein's dynamic smoothing, the values of the library's own check | N(start, σ₀² I) |
| BO | `models::BayesianOptimisation`: Matérn-5/2 kernel with one length scale per coordinate, expected improvement, an initial design of 2n + 2 = 6 uniform points, one point per ask | uniform in the box |
| random | `search::RandomSearch`, uniform in the box | — |

The start points and initial step sizes of the (1+1)-ES, CMA-ES and the cross-entropy method are fixed per function, the same for every seed: (−1.5, −0.5) with σ₀ = 0.5 on Rosenbrock's function, (3.5, −3.0) with σ₀ = 2 on Rastrigin's, (2.5, 12.0) with σ₀ = 3 on Branin's. These three search ℝ² without bounds; the others search the box. Budgets: 2000 evaluations, and 100 for Bayesian optimisation, whose evaluations are meant to be expensive (each of its asks refits the surrogate's hyperparameters by maximum likelihood).

## 5. The Data

No external data: the functions are those of `heuresis::benchmarks` (`continuous.hpp`, with Branin's constants from Surjanovic and Bingham's Virtual Library of Simulation Experiments). Everything in `results/` is written by `landscape.cpp`:

| File | Content |
|---|---|
| `runs.csv` | one line per run: evaluations spent, final f − f*, how it ended (CMA-ES with its stopping criterion), runtimes to 10⁻², 10⁻⁴ and 10⁻⁸ |
| `anytime.csv` | the quartiles of the best f − f* after 1, 2, …, 2000 evaluations (ten budgets to a decade) |
| `spread.csv` | the quartiles of the distance to the nearest minimiser and of the spread at the same budgets |
| `landscapes.json` | f − f* of each function on a 141 × 141 grid of its box (four significant digits) |
| `trajectories.json` | seed 0 of every algorithm on every function, generation by generation (1.5 MB) |

`trajectories.json` is meant for a viewer. Every value in it is a gap f − f*, rounded to four significant digits; for each function it holds the box, f*, the minimisers, the start and σ₀, and for each algorithm its ending and its generations. A generation holds `e` (evaluations after it), `gap` (the best gap so far), `pts` and `f` (its candidates and their gaps), and as the algorithm has them `m`, `s` and `C` (the state it was drawn from: the parent and σ of the (1+1)-ES, the mean, σ and C of CMA-ES, whose law is N(m, s²C), the mean and covariance of the cross-entropy method; C as its entries (1,1), (1,2), (2,2)), `v` (PSO's velocities after the generation: from the second generation on, the steps that brought the particles to `pts`), and `pop` and `pf` (DE's population after selection). Bayesian optimisation adds `gp`, its surrogate's mean and standard deviation on a 50 × 50 grid of the box after 6, 15, 30, 60 and 100 evaluations.

## 6. The Results

Medians over the 31 seeds; quartiles in brackets.

**Rosenbrock's valley.**
- CMA-ES reaches f − f* ≤ 10⁻⁸ in every run, after a median of 595 evaluations. Its covariance turns to lie along the valley and its mean follows the parabola round to (1, 1) (search-01).
- DE reaches 10⁻⁸ in 27 of 31 runs (median 964 evaluations); its population contracts onto the valley floor and then along it.
- Bayesian optimisation is the most economical at low precision: 30 of 31 runs reach 10⁻² (median 46 evaluations) and 27 reach 10⁻⁴ (median 77) within its 100 evaluations. It ends at 2.1 × 10⁻⁵ [7.6 × 10⁻⁶, 5.5 × 10⁻⁵]; no run reaches 10⁻⁸.
- The (1+1)-ES ends at 3.1 × 10⁻³ [1.4 × 10⁻³, 4.6 × 10⁻³] after 2000 evaluations, only six times better than random search's 1.9 × 10⁻² [4.4 × 10⁻³, 3.5 × 10⁻²]. Its isotropic step cannot lie along the valley, so its step size shrinks to the valley's width (its spread σ√n is 8 × 10⁻⁴ at the end) and it creeps.
- The Standard PSO 2011 ends at 5.3 × 10⁻⁴; its particles still fly at a spread of 0.67 while its best lies 0.038 from the minimiser.
- The cross-entropy method stalls: it ends at f − f* = 1.11 [0.84, 1.5], its best point 1.45 from the minimiser (with seed 0, near (0, 0) on the valley floor), with a spread of 0.014. Its covariance shrinks faster than its mean travels along the curved valley, which dynamic smoothing slows but does not prevent at these parameters.

**Rastrigin's egg crate.**
- DE finds the global minimum in every run: 31 of 31 reach 10⁻⁸, after a median of 1359 evaluations. Its population of 20, spread over the box, sees the bowl beneath the lattice.
- CMA-ES with its default population of 6, started at (3.5, −3.0), finds the global minimum in 4 of 31 runs; the other 27 stop by TolFun at a local minimum after 564 to 828 evaluations. Its median final gap, 1.23, is no better than random search's 1.22 after 2000 evaluations. In 7 of the 27 stopped runs the best point seen is not the local minimum the mean converged to but a lucky early sample in a better basin: with seed 21 a sample at (−0.009, −0.040) has f = 0.331, in the global basin, while the mean settled at (−0.995, 0.995).
- The (1+1)-ES is trapped soonest: 2 of 31 runs find the global minimum, and the median ends at 4.97 with a spread σ√n of 6 × 10⁻¹⁶, its point 2.2 from the origin. Without a stopping rule (heuresis's default `min_sigma` is 0) it spends the rest of its budget there.
- The cross-entropy method ends at 1.0, a local minimum next to the global one; Bayesian optimisation at 1.01 [0.25, 1.9] after 100 evaluations, its surrogate a field of small basins round its samples; the swarm at 9.2 × 10⁻³, still flying (spread 1.65).

**Branin's three minima.**
- The (1+1)-ES, CMA-ES and DE reach 10⁻⁸ in every run (medians 210, 305 and 1021 evaluations), each settling in one of the three minima.
- Bayesian optimisation reaches 10⁻⁴ in every run after a median of 32 evaluations, and ends at 1.5 × 10⁻⁶ after 100. On this smooth function it is the most economical: the (1+1)-ES, the fastest of the others, needs a median of 116 evaluations for 10⁻⁴, CMA-ES 176.
- The cross-entropy method ends at 7.3 × 10⁻⁶ and the swarm at 6.6 × 10⁻⁵ after 2000.

**Across the three.** The spread of a search tells whether it is still searching. Where an algorithm succeeds its spread falls with its distance to the optimum, at about the same rate (CMA-ES on Rosenbrock's function: distance 6.2 × 10⁻⁷, spread 4.8 × 10⁻⁶). Where it fails, the spread falls while the distance stays near 1: CMA-ES and the (1+1)-ES on Rastrigin's function, the cross-entropy method on Rosenbrock's. The Standard PSO 2011 is the opposite case: its spread barely falls at all. heuresis's `swarms` README states the cause, the spectral radius 0.98 of its particles' second moments in two dimensions, so that they settle only slowly about their attractors.

## 7. The Figures

| Image | Content |
|---|---|
| [`search-01-rosenbrock.jpg`](../../../assets/gallery/search-01-rosenbrock.jpg) | the six searches on Rosenbrock's function after about 50, 300 and 1500 evaluations (Bayesian optimisation: 6, 30 and 100) |
| [`search-02-rastrigin.jpg`](../../../assets/gallery/search-02-rastrigin.jpg) | the same on Rastrigin's function |
| [`search-03-branin.jpg`](../../../assets/gallery/search-03-branin.jpg) | the same on Branin's function |
| [`search-04-anytime.jpg`](../../../assets/gallery/search-04-anytime.jpg) | the best f − f* against evaluations, median and quartiles over 31 seeds, seven algorithms including random search |
| [`search-05-spread.jpg`](../../../assets/gallery/search-05-spread.jpg) | the distance of the best point to the nearest minimiser, and the spread of each search, against evaluations |

## 8. The Caveats

- **Two dimensions.** The landscapes were chosen to be drawn. Behaviour in two dimensions does not predict behaviour in twenty: the Standard PSO 2011, slow here, solves not even the sphere from ten dimensions on ([the benchmark](../2026-11-05-benchmarking-continuous-search/README.md)), and CMA-ES's default population grows only as 4 + ⌊3 ln n⌋.
- **One start.** The (1+1)-ES, CMA-ES and the cross-entropy method start from one fixed point per function, which is what makes their trajectories comparable; their success rates on Rastrigin's function depend on it. Population methods start uniform in the box and see the whole of it at once, an advantage on Rastrigin's function that the fixed start does not give the others.
- **Defaults, not tuning.** Every algorithm runs at heuresis's defaults, the cross-entropy method at the values of the library's own check. Other settings change the picture: a larger population for CMA-ES on Rastrigin's function (as IPOP's restarts do in the benchmark), slower smoothing for the cross-entropy method.
- **Unequal budgets.** Bayesian optimisation is given 100 evaluations and the others 2000. It spends far more computation per evaluation (seconds against microseconds) and is meant for objectives that cost far more than that.
- **Rounding.** The trajectories are rounded to four significant digits, so positions near an optimum lose resolution there (a point within 10⁻⁵ of (1, 1) is written as (1, 1)); `spread.csv` and `anytime.csv` are computed at full precision.
