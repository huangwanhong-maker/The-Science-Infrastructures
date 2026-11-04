# Repertoires of an Arm: Quality and Diversity on a Planar Arm of Twelve Joints

*Status: draft (2026-11-05). Research: a case study of heuresis's `diversity` module ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).*

## 1. The Question

A search for quality and diversity does not ask for the best solution but for the best solution of every kind. For an arm, the kinds are the places its tip can reach, and the quality of an arm for a place is how evenly it bends to get there. The answer is a repertoire: one arm for each place, which is what a robot that must adapt after damage reaches for (Cully, Clune, Tarapore and Mouret 2015).

The study asks three things of heuresis's searches on its benchmark arm of twelve joints:

1. how much of the reachable disc each search fills, and how fast;
2. how good the arms it keeps are;
3. whether the iso+line variation of Vassiliades and Mouret (2018), which steps along the line between two elites, does better than Gaussian mutation at any of three step sizes.

A single-objective search (CMA-ES minimising the bend alone) is the contrast: it answers "what is best?", not "what is possible here?".

## 2. The Method

**The arm.** `heuresis::benchmarks::planar_arm` has d = 12 links of length 1/d from a base at the origin. Its genotype is the vector of joint angles θ ∈ [−π, π]¹², each angle measured from the previous link, so the j-th link points at θ₁ + … + θⱼ. Its descriptor is the tip, mapped from [−1, 1]² to [0, 1]², so the reachable set is the disc of radius ½ about (½, ½). Its objective, minimised, is the variance of the twelve angles. An arm with every angle equal has variance 0, so the best arms form a one-parameter family: for a common angle a the tip is (1/12) Σⱼ (cos ja, sin ja), a curve from the straight arm at the rim (a = 0) to the closed regular dodecagon at the base (a = 30°).

**The searches**, every one of heuresis, each given 100 000 evaluations and run with 15 seeds:

| Search | heuresis | Variation |
|---|---|---|
| MAP-Elites, iso+line | `diversity::MapElites<Real, Grid>`, 1 000 random samples first | iso+line, σ_iso = 0.05, σ_line = 0.3, then reflection into the box |
| MAP-Elites, Gaussian 0.03, 0.1, 0.3 | the same | Gaussian mutation of every angle at σ (rad), then reflection |
| CVT-MAP-Elites, iso+line | `diversity::MapElites<Real, Voronoi>`, 1 000 sites | iso+line, as above |
| novelty search | `diversity::NoveltySearch<Real>`: population 100, k = 15, initial threshold 0.02 | Gaussian, σ = 0.1 |
| random search | `search::RandomSearch<Real>`, uniform on the box | none |
| CMA-ES, IPOP | `strategies::CmaRestarts`, σ₀ = 1, the box by Hansen's transformation, up to 64 runs | its own |

MAP-Elites breeds one candidate per ask (the classical steady state). The CVT sites are 1 000 points drawn uniformly in the reachable disc from one keyed stream. heuresis asks its caller for the sites, because k-means (which would make them centroidal) belongs to tyche and is not yet written.

**One measure for all.** Every evaluation of every search is also offered to a passive archive on a grid of 50 × 50 cells over [0, 1]². All searches are measured on it:

- **coverage**, the share of the 2 040 cells that meet the reachable disc which hold an elite;
- **the QD-score**, Σ (π² − var) over the filled cells (heuresis's `Archive::qd_score` with offset π², the largest variance of angles in [−π, π]), divided by π² × 2 040, so that it lies in [0, 1];
- **the elites' mean variance**, which says how good the kept arms are, apart from how many there are.

For grid MAP-Elites the passive archive is its own archive. For CVT-MAP-Elites its own coverage of its 1 000 cells is recorded as well.

## 3. The Code

- heuresis `6b3cc88`; calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`; the infrastructure's root `4153cb9`, in which heuresis is a submodule (`results/versions.json`). The runs were made from the working tree that became that commit.
- `repertoires.cpp`, built by `CMakeLists.txt` as a project of its own (C++20, Release, g++ 11.5 in WSL):

```
cmake -S . -B ~/.cache/science-build/exp-arm -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-arm -j4
~/.cache/science-build/exp-arm/repertoires results
```

- `figures.py` draws the figures from `results/` in the house style (`../style/search_style.py`).

## 4. The Parameters

All of them are in `results/reference.json`: 12 joints in [−π, π]; 100 000 evaluations; 15 seeds (keyed draws, seeds 1–15); a grid of 50 × 50; 1 000 CVT sites; QD offset π². The run took 233 s on 8 threads of a machine shared with other work.

## 5. The Data

None from outside. The arm is heuresis's benchmark, and every number comes from the run.

## 6. The Results

### 6.1 The Repertoire Searches Fill the Disc, the Others Half of It

After 100 000 evaluations (medians over 15 seeds; ranges in brackets):

| Search | Coverage at 10 000 | Coverage at 100 000 | QD-score | Mean variance of the elites |
|---|---|---|---|---|
| MAP-Elites, iso+line | 82.8% | 96.9% [95.8, 98.4] | 0.965 | 0.039 |
| MAP-Elites, Gaussian 0.1 | 59.6% | 96.6% [95.5, 98.0] | 0.961 | 0.048 |
| MAP-Elites, Gaussian 0.3 | 67.9% | 95.5% [95.2, 96.0] | 0.946 | 0.096 |
| MAP-Elites, Gaussian 0.03 | 42.7% | 80.9% [74.1, 85.5] | 0.787 | 0.318 |
| CVT-MAP-Elites, iso+line | 85.6% | 95.7% [94.6, 97.0] | 0.954 | 0.035 |
| novelty search | 80.9% | 96.9% [95.7, 97.5] | 0.940 | 0.285 |
| random search | 36.3% | 51.4% [50.9, 52.8] | 0.452 | 1.213 |
| CMA-ES, IPOP | 32.1% | 50.1% [43.9, 64.4] | 0.465 | 0.758 |

Random search reaches only half the cells, because a uniformly random arm rarely ends near the rim. With every angle uniform, the absolute direction of each link is uniform and independent of the others, so the tip is the mean of twelve independent random unit steps, whose root-mean-square distance from the base is 1/√12 ≈ 0.29 of the reach. CMA-ES does no better. Each of its eight runs (λ from 11 to 1 408) converges onto the curve of equal-angle arms, which it finds exactly (best variance 1.1 × 10⁻¹⁷), and its evaluations on the way fill the cells they pass through.

### 6.2 Iso+line Gets There Two to Three Times Sooner

The median MAP-Elites run with iso+line covers 80% of the cells after 9 000 evaluations, 90% after 18 000 and 95% after 42 000. With the best Gaussian step (σ = 0.1) the same takes 26 000, 43 000 and 66 000. A larger Gaussian step (0.3) is faster at first but keeps worse arms (mean variance 0.096), and a smaller one (0.03) has not reached 90% after 100 000. At the end, iso+line and Gaussian 0.1 cover the same share (96.9% and 96.6%), but iso+line's elites are better (mean variance 0.039 against 0.048, interquartile ranges [0.038, 0.041] and [0.046, 0.051], which do not overlap).

*Rationale.* The elites of neighbouring cells are similar arms, so the elites lie near a low-dimensional set in the space of genotypes, the "elite hypervolume" of Vassiliades and Mouret. A step along the line between two elites stays near that set, and its length scales with the distance between them. An isotropic Gaussian step does neither: it must be small to stay near the set and large to travel along it.

### 6.3 Voronoi Cells Fill Their Own Tessellation Almost Whole

CVT-MAP-Elites fills 98.6% of its 1 000 cells (median over seeds; range 97.6–99.7%). It is the fastest to 80% and 90% of the common grid (8 000 and 15 000 evaluations), and its elites are the evenest (mean variance 0.035). It keeps one elite per Voronoi cell, about half as many as the grid has cells, and on the common grid it ends at 95.7%, below grid MAP-Elites.

### 6.4 Novelty Search Finds Every Place but Not the Best Arm for It

Novelty search selects for being unlike what it has seen and ignores the objective. It reaches as many cells as MAP-Elites (96.9%), and as fast at first (80.9% at 10 000), but the best arm it ever evaluated for a cell has seven times MAP-Elites' variance (0.285 against 0.039). Its archive of behaviours ends with a median of 1 088 points. Computing sparseness by a scan over the population and that archive made it the slowest search by far, 89–105 s a run against 0.3–7 s for the MAP-Elites variants, on the shared machine.

### 6.5 The Maps Show the Equal-Angle Curve

In the maps of MAP-Elites and CMA-ES (search-24), the darkest cells, those with the evenest arms, trace the curve of equal-angle arms described in §2: it runs from the straight arm at the right of the rim, round the top (a > 0) and the bottom (a < 0), into curls about the base. The rest of the repertoire grows outwards from it. In the arms drawn (search-25), the elites reaching near the base are coiled into nearly closed polygons, and those reaching the rim are nearly straight.

## 7. The Figures

- `search-24-repertoire-maps.jpg`: the elites of the median run of six searches over the descriptor space, coloured by the standard deviation of the angles.
- `search-25-arms.jpg`: arms of the median MAP-Elites (iso+line) run for places on rings about the base.
- `search-26-coverage-and-qd.jpg`: coverage, QD-score and the elites' spread against evaluations.

Their captions are in [`gallery.md`](gallery.md).

## 8. Caveats

- **The QD-score follows coverage.** With the offset π² ≈ 9.87 and elites' variances below 0.4, every filled cell adds nearly the same amount, so the normalised QD-score is coverage less a few per cent. The elites' mean variance is reported beside it for that reason. A smaller offset would weigh quality more, but it would no longer bound every variance.
- **Coverage cannot reach 1.** A cell on the rim that meets the disc only at its corner is reached only by a nearly straight arm pointing exactly there. The grid cells meeting the disc are counted whatever the size of the overlap.
- **The CVT sites are uniform, not centroidal.** Uniform random sites give Voronoi cells of unequal area, so the coverage of CVT-MAP-Elites' own tessellation weighs places unevenly. Centroidal sites need tyche's k-means, which is not yet written.
- **One arm, one budget.** The angle box [−π, π], the 12 joints and the budget follow heuresis's check of the module. Vassiliades and Mouret (2018) use other bounds and far more joints, and their comparison of iso+line with Gaussian variation concerns those.
- **Times are indicative.** The machine ran other work during the run, with a load near 20 on 8 cores, so the seconds per run compare searches only roughly.

## 9. References

- Cully, A., Clune, J., Tarapore, D. and Mouret, J.-B. (2015), Robots that can adapt like animals, *Nature* 521, 503–507.
- Lehman, J. and Stanley, K. O. (2011), Abandoning objectives: evolution through the search for novelty alone, *Evol. Comput.* 19, 189–223.
- Mouret, J.-B. and Clune, J. (2015), Illuminating search spaces by mapping elites, arXiv:1504.04909.
- Hansen, N. and Ostermeier, A. (2001), Completely derandomized self-adaptation in evolution strategies, *Evol. Comput.* 9, 159–195; Auger, A. and Hansen, N. (2005), A restart CMA evolution strategy with increasing population size, *CEC*, 1769–1776.
- Vassiliades, V., Chatzilygeroudis, K. and Mouret, J.-B. (2018), Using centroidal Voronoi tessellations to scale up the multidimensional archive of phenotypic elites algorithm, *IEEE Trans. Evol. Comput.* 22, 623–630.
- Vassiliades, V. and Mouret, J.-B. (2018), Discovering the elite hypervolume by leveraging interspecies correlation, *GECCO*, 149–156.
