# Fronts in Three Objectives: Five Searches on DTLZ1, DTLZ2, DTLZ7, ZDT1 and ZDT3

*Status: draft (2026-11-05). Research: a case study of heuresis's `multiobjective` module ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).*

## 1. The Question

With several objectives, a search returns not a point but a finite set that should lie on the Pareto front and spread over it. Two questions about such a set are separate:
- **convergence**: how close it is to the front;
- **spread**: how it covers the front, including fronts made of several disconnected pieces.

The four algorithms of heuresis answer the second differently:
- NSGA-II cuts its last front by crowding distance;
- NSGA-III spreads the population among Das and Dennis's reference directions;
- MOEA/D solves one scalar problem per weight, by Tchebycheff's function or by PBI;
- SMS-EMOA removes, one at a time, the point that contributes least hypervolume.

The study asks how they compare on five standard problems whose fronts are known exactly. Each is measured by the exact hypervolume and by IGD against evaluations. The questions for the comparison are three: who converges on the multimodal DTLZ1, who covers the disconnected fronts of DTLZ7 and ZDT3, and where each puts its points.

## 2. The Method

### 2.1 The Problems

All five are heuresis's benchmarks (`benchmarks/multiobjective.hpp`), minimised on [0, 1]ⁿ:

| Problem | m | n | Front | Reference point r | HV of the front | Budget |
|---|---|---|---|---|---|---|
| DTLZ1 | 3 | 7 | the simplex Σfᵢ = ½, behind 11⁵ − 1 local fronts | (1, 1, 1) | 1 − 2⁻ᵐ/m! = 0.979167 | 40 000 |
| DTLZ2 | 3 | 12 | the positive octant of the unit sphere | (1, 1, 1) | 1 − π/6 = 0.476401 | 40 000 |
| DTLZ7 | 3 | 22 | four disconnected patches | (1, 1, 7) | 2.457669 (§2.2) | 40 000 |
| ZDT1 | 2 | 30 | f₂ = 1 − √f₁ | (1, 1) | 2/3 | 25 000 |
| ZDT3 | 2 | 30 | five disconnected segments | (1, 1) | 1.044426 (§2.2) | 25 000 |

### 2.2 The Fronts of ZDT3 and DTLZ7

heuresis's benchmarks give these two problems but not their fronts or hypervolumes, so they are derived here.

**ZDT3.** At g = 1 the problem is f = (t, z(t)) with z(t) = 1 − √t − t sin 10πt. The point at t is non-dominated among the curve's points exactly when z(t) < z(s) for every s < t. The front is therefore the set of strict running records of z from below.

**DTLZ7 (m = 3).** At g = 1 the problem is f = (t₁, t₂, 6 − φ(t₁) − φ(t₂)) with φ(t) = t(1 + sin 3πt). Because f₃ is separable, (t₁, t₂) is non-dominated exactly when each tᵢ is a strict running record of φ from above.
- If a′ ≤ a with φ(a′) ≥ φ(a), and the two are not both equal, then (a′, b) dominates (a, b).
- Conversely, any point dominating (a, b) makes a′ or b′ such a witness.

The front is the product T* × T* of the record set with itself. T* is two intervals, so the front is four patches.

**Hypervolumes.** Above each (y₁) or (y₁, y₂), the dominated region rises from the front's lower envelope to r. So, with Ψ(y) = min over t ≤ y of z(t) and Φ(y) = max over t ≤ y of φ(t):

  HV(ZDT3) = ∫₀¹ (1 − Ψ(y)) dy,    HV(DTLZ7) = 1 + 2 ∫₀¹ Φ(y) dy.

Each integral is computed by the midpoint rule on 4 × 10⁶ cells, and agrees with 10⁵ and 10⁶ cells to 10⁻⁹. Three checks are in `results/reference.json`:
- non-dominated sorting of a 61 × 61 grid of DTLZ7's g = 1 surface keeps exactly the 961 = 31² points of the product of the records;
- of 2 001 points of ZDT3's curve it keeps exactly the 536 records;
- the exact hypervolumes of sampled fronts converge to the integrals from below: for DTLZ7, 2.449182, 2.455549 and 2.456822 at grids of 201, 801 and 2 001, against 2.457669; for ZDT3, 1.043518, 1.044337 and 1.044417 at 1 001, 10 001 and 100 001, against 1.044426.

### 2.3 The Searches

| Search | heuresis | Population | Variation |
|---|---|---|---|
| NSGA-II | `Nsga2<Real>` | 92 (ZDT: 100) | SBX η_c = 20 at 0.9, polynomial mutation η_m = 20 at 1/n |
| NSGA-III | `Nsga3<Real>`, Das–Dennis directions p = 12 (91; ZDT: p = 99, 100) | 92 (ZDT: 100) | SBX η_c = 30 at 1.0, η_m = 20 (Deb and Jain 2014, §VI) |
| MOEA/D, Tchebycheff | `Moead<Real>`, the same weights, T = 20, n_r = T | 91 (ZDT: 100) | as NSGA-II |
| MOEA/D, PBI | the same, θ = 5 | 91 (ZDT: 100) | as NSGA-II |
| SMS-EMOA | `SmsEmoa<Real>`, one offspring per step | 92 (ZDT: 100) | as NSGA-II |

Each algorithm is driven by hand through the ask–tell protocol, with 15 seeds (keyed draws, seeds 1–15). At every 500 evaluations its current non-dominated set is measured, that is, what it held after its last completed update; for the generational algorithms, the last generation told whole. The measures are:
- the **hypervolume** by heuresis's exact algorithms (the sweep in two objectives, Beume et al.'s in three), divided by the front's;
- **IGD**, the mean distance from a dense reference set to the nearest point. The reference sets are 5 050 Das–Dennis points (p = 99) placed on the front for DTLZ1 and DTLZ2, 5 476 points of DTLZ7's front (the records of a grid of 151, squared), 1 000 points of ZDT1's and 1 066 records of ZDT3's.

At the end:
- a point counts as **on the front** when it lies within 0.01 of it, in Euclidean distance;
- the front's **components** (DTLZ7's patches, ZDT3's segments) are counted by the points on each;
- for DTLZ1 and DTLZ2, a reference direction is **occupied** when it is the nearest, by heuresis's NSGA-III association (`associate`), to a point on the front. Points of DTLZ1 are first divided by its intercept ½.

To set the hypervolumes of 91 or 100 points against something attainable, the same measures are taken of a natural set of μ points lying on the front, recorded in `results/reference.json`: the 91 Das–Dennis directions where they meet the fronts of DTLZ1 and DTLZ2, and 100 points of ZDT1's front evenly spaced in f₁. A finite set covers less than the whole front, so its hypervolume ratio is below 1.

## 3. The Code

- heuresis `6b3cc88`; calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`; the infrastructure's root `4153cb9`, in which heuresis is a submodule (`results/versions.json`). The runs were made from the working tree that became that commit.
- `fronts.cpp`, built by `CMakeLists.txt` as a project of its own (C++20, Release, g++ 11.5 in WSL):

```
cmake -S . -B ~/.cache/science-build/exp-fronts -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-fronts -j4
~/.cache/science-build/exp-fronts/fronts results
```

- `figures.py` draws the figures from `results/` in the house style (`../style/search_style.py`).

The program was run twice, the second time with the reference sets of §2.3 added. Every table but the seconds per run came out identical, as heuresis's keyed draws promise.

## 4. The Parameters

All of them are in `results/reference.json` and §2: the budgets, populations, directions, reference points and reference sets; 15 seeds; checkpoints every 500 evaluations; the tolerance 0.01 for "on the front". The second run took 1 268 s on 8 threads of a machine shared with other work.

## 5. The Data

None from outside. The problems are heuresis's benchmarks, and every number comes from the runs.

## 6. The Results

Medians over 15 seeds, with the interquartile range or the extremes where they matter. The tables are `results/final.csv` (each run's end) and `results/curves.csv` (the quartiles at each checkpoint).

| Problem | NSGA-II | NSGA-III | MOEA/D, Tchebycheff | MOEA/D, PBI | SMS-EMOA | μ points on the front |
|---|---|---|---|---|---|---|
| **HV / HV of the front** | | | | | | |
| DTLZ1 | 0.9898 | 0.9940 | 0.9877 | 0.9941 | 0.9945 | 0.9944 |
| DTLZ2 | 0.7714 | 0.8675 | 0.7753 | 0.8685 | 0.8916 | 0.8687 |
| DTLZ7 | 0.9254 | 0.9305 | 0.8517 | 0.9037 | 0.9676 | — |
| ZDT1 | 0.9889 | 0.9905 | 0.9898 | 0.9775 | 0.9929 | 0.9921 |
| ZDT3 | 0.9955 | 0.9941 | 0.9937 | 0.9723 | 0.9976 | — |
| **IGD** | | | | | | |
| DTLZ1 | 0.0305 | 0.0207 | 0.0310 | 0.0206 | 0.0205 | 0.0205 |
| DTLZ2 | 0.0731 | 0.0543 | 0.0757 | 0.0543 | 0.0773 | 0.0543 |
| DTLZ7 | 0.0811 | 0.0773 | 0.2036 | 0.1498 | 0.0775 | — |
| ZDT1 | 0.0049 | 0.0041 | 0.0045 | 0.0105 | 0.0036 | 0.0037 |
| ZDT3 | 0.0054 | 0.0061 | 0.0109 | 0.0163 | 0.0043 | — |

### 6.1 Every Search Converges on DTLZ1, NSGA-II Not Always

DTLZ1's distance function is Rastrigin's in five variables, with 11⁵ − 1 local fronts parallel to the true one.
- NSGA-III, both MOEA/Ds and SMS-EMOA end on the front in all 15 runs, with ratios from 0.9873 to 0.9945. The 91 Das–Dennis points on the front give 0.9944, which SMS-EMOA passes and NSGA-III and PBI meet (0.9940, 0.9941; their IGD equals the set's, 0.0205).
- NSGA-II passes 0.98 in 13 runs. In one (seed 3) no point of its final population lies within 0.01 of the front (HV ratio 0.776), and another ends at 0.975.
- In the median curves, MOEA/D with Tchebycheff reaches 90% of the front's hypervolume first (9 000 evaluations), then MOEA/D with PBI (10 500), SMS-EMOA (11 500), NSGA-III (13 000) and NSGA-II (17 500).

Crowding gives NSGA-II no pressure towards the front within a front, and its tournaments on crowding keep diversity at the cost of speed.

### 6.2 The Hypervolume of μ Points Has a Ceiling, and SMS-EMOA Passes the Evenly Spread One

On DTLZ2 no search comes near the front's whole hypervolume, because 91 or 92 points on a sphere cover only part of it. The 91 Das–Dennis directions where they meet the sphere have a ratio of 0.8687.
- MOEA/D with PBI (0.8685) and NSGA-III (0.8675) stand at that set, within 0.0012: they converge onto those directions.
- SMS-EMOA exceeds it, at 0.8916, because it maximises hypervolume itself.
- NSGA-II (0.7714) and MOEA/D with Tchebycheff (0.7753) fall below it.

IGD ranks them otherwise. SMS-EMOA's IGD on DTLZ2 (0.0773) is the worst of the five, and the evenly spread sets of NSGA-III and PBI are the best (0.0543, the IGD of the Das–Dennis set itself). A hypervolume-maximising set is not an even one. It puts many points on the front's edges, where its contributions are large because SMS-EMOA's own reference point is the worst value plus one in each objective, not (1, 1, 1). Of the 91 directions it occupies only 56 on DTLZ2 (search-23).

### 6.3 Where the Points Go

The occupied directions, out of 91, in the median run, with the extremes over the 15 runs in brackets:

| | NSGA-II | NSGA-III | MOEA/D, Tchebycheff | MOEA/D, PBI | SMS-EMOA |
|---|---|---|---|---|---|
| DTLZ1 | 50 [0–57] | 91 [79–91] | 49 [48–50] | 91 [91–91] | 84 [81–86] |
| DTLZ2 | 41 [29–47] | 91 [88–91] | 49 [47–50] | 91 [91–91] | 56 [54–58] |

- **NSGA-III** places its points on the reference directions by construction.
- **MOEA/D with PBI** places them there too. With θ = 5 the optimum of each subproblem lies on its own weight's ray.
- **MOEA/D with Tchebycheff** does not. On a front such as DTLZ1's simplex, the minimiser of maxᵢ wᵢ|fᵢ − zᵢ*| lies where every wᵢ(fᵢ − zᵢ*) is equal, on the ray through (1/w₁, 1/w₂, 1/w₃), not through w. Its points therefore crowd towards the middle and away from the edges, in rings (search-20, search-23). The weights with a zero component, which heuresis replaces by 10⁻⁶, send many subproblems to the same edge points. Of its 91 solutions only 83 (DTLZ1) and 85 (DTLZ2) are non-dominated, and they occupy 49 directions.
- **NSGA-II** spreads its points by crowding distance in objective space, which is uneven in three objectives.

### 6.4 Disconnected Fronts: SMS-EMOA Covers Them Best, MOEA/D Leaves Them Uneven

On DTLZ7, every search covers all four patches in at least 14 of 15 runs:
- SMS-EMOA has the highest ratio by far, 0.9676, and the most even share of points among the patches (17, 23, 23 and 30 in the median run). In one run (seed 15) it lost the two patches with f₁ in the upper interval (HV ratio 0.847). The single failures of NSGA-II and MOEA/D with Tchebycheff likewise lose two patches at once, a whole interval of one coordinate.
- NSGA-III covers all four in every run (0.9305), and NSGA-II in 14 (0.9254).
- MOEA/D covers them unevenly. Its weights that point at the gaps between the patches have their optima on the patches' edges, so with PBI the median run puts 2, 7, 7 and 32 points on the four patches (ratio 0.9037, IGD 0.150). With Tchebycheff the median ratio is 0.8517 and IGD 0.204, the worst.

On ZDT3's five segments:
- SMS-EMOA is again best (0.9976, all five segments in 14 of 15 runs);
- then NSGA-II (0.9955, 13), NSGA-III (0.9941, 15) and MOEA/D with Tchebycheff (0.9937, 14);
- MOEA/D with PBI is last (0.9723, 12), its points off the front at the joins between segments.

### 6.5 ZDT1 and the Cost of the Exact Contributions

On ZDT1 every search but MOEA/D with PBI ends within 1.1% of the front's hypervolume. SMS-EMOA comes closest (0.9929), above 100 points evenly spaced in f₁ (0.9921), and with a lower IGD than theirs (0.0036 against 0.0037). PBI's penalty, with θ = 5, slows it on this convex front (0.9775, and 0.938 in its worst run).

SMS-EMOA's lead costs time. In three objectives heuresis computes each point's exclusive contribution as a difference of hypervolumes, about 93 three-dimensional sweeps for every evaluation. In the committed run its median times were 183 s on DTLZ1, 271 s on DTLZ2 and 195 s on DTLZ7, against 0.5 to 3 s for every other search on every problem. heuresis's `multiobjective` README names the remedy, Guerreiro and Fonseca's incremental HV3D⁺.

## 7. The Figures

- `search-20-fronts-dtlz1-dtlz2.jpg`: the final populations of the median runs on DTLZ1 and DTLZ2 in three dimensions, against the true fronts.
- `search-21-disconnected-fronts.jpg`: the same on DTLZ7's four patches and ZDT3's five segments.
- `search-22-convergence.jpg`: hypervolume and IGD against evaluations.
- `search-23-spread-on-the-simplex.jpg`: the final points projected on the simplex among the 91 reference directions.

Their captions are in [`gallery.md`](gallery.md).

## 8. Caveats

- **One setting each.** The populations, η's, T, n_r and θ are the published defaults, not tuned per problem. MOEA/D's n_r = T (every improved neighbour replaced) is the original's choice. Li and Zhang's n_r = 2 keeps more diversity and might cover DTLZ7 more evenly.
- **heuresis's MOEA/D is the batch form.** A generation of 91 children is bred from the solutions as they stood at its start, then applied in subproblem order (decision 0003). The original updates after each child, which spreads information faster.
- **The reference point decides the hypervolume.** With r = (1, 1, 1), points on DTLZ2's boundary (some fᵢ = 1) contribute nothing, and SMS-EMOA's internal reference (worst plus one) rewards them. Another r would reorder SMS-EMOA and the evenly spread sets.
- **IGD depends on its reference set.** The Das–Dennis sample of DTLZ2's front is denser towards the middle of the octant than near its edges, which favours the evenly spread sets.
- **The median run** shown in the pictures is the eighth of fifteen by final hypervolume, one run, not an average.
- **Times are indicative.** The machine ran other work during both runs (load 20–50 on 8 cores), so the seconds compare searches only roughly.

## 9. References

- Beume, N., Fonseca, C. M., López-Ibáñez, M., Paquete, L. and Vahrenhold, J. (2009), On the complexity of computing the hypervolume indicator, *IEEE Trans. Evol. Comput.* 13, 1075–1082.
- Beume, N., Naujoks, B. and Emmerich, M. (2007), SMS-EMOA: multiobjective selection based on dominated hypervolume, *Eur. J. Oper. Res.* 181, 1653–1669.
- Das, I. and Dennis, J. E. (1998), Normal-boundary intersection, *SIAM J. Optim.* 8, 631–657.
- Deb, K. and Jain, H. (2014), An evolutionary many-objective optimization algorithm using reference-point-based nondominated sorting approach, part I, *IEEE Trans. Evol. Comput.* 18, 577–601.
- Deb, K., Pratap, A., Agarwal, S. and Meyarivan, T. (2002), A fast and elitist multiobjective genetic algorithm: NSGA-II, *IEEE Trans. Evol. Comput.* 6, 182–197.
- Deb, K., Thiele, L., Laumanns, M. and Zitzler, E. (2005), Scalable test problems for evolutionary multiobjective optimization, in *Evolutionary Multiobjective Optimization*, Springer, 105–145.
- Guerreiro, A. P. and Fonseca, C. M. (2018), Computing and updating hypervolume contributions in up to four dimensions, *IEEE Trans. Evol. Comput.* 22, 449–463.
- Li, H. and Zhang, Q. (2009), Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II, *IEEE Trans. Evol. Comput.* 13, 284–302.
- Zhang, Q. and Li, H. (2007), MOEA/D: a multiobjective evolutionary algorithm based on decomposition, *IEEE Trans. Evol. Comput.* 11, 712–731.
- Zitzler, E., Deb, K. and Thiele, L. (2000), Comparison of multiobjective evolutionary algorithms: empirical results, *Evol. Comput.* 8, 173–195.
