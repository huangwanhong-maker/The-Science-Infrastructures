# Tuning an Image Pipeline on Phantoms

*Status: draft (2026-11-05).* Research and use: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) on the foundation eikon's reconstruction of tubes from images. Synthetic phantoms only.

## 1. The Question

A chain that recovers vessels from an image has many parameters: the scales of its vesselness filter, the filter's constants, the thresholds that cut its response, the level at which a vessel's boundary is drawn. Their defaults are known to be insufficient (Lamy et al. 2022), so they are tuned, usually by hand on a few images. Tuning is a search over a space without gradients whose every point costs a reconstruction, and it has two difficulties of its own. Its objectives disagree: overlap of volumes (Dice) rewards the thick vessels that hold most voxels, overlap of centrelines (clDice; Shit et al. 2021) rewards every branch equally. And a chain tuned on a few images may not carry to the next.

On phantoms of exact geometry both can be measured. This study tunes eikon's chain on two phantoms with NSGA-II for Dice against clDice, and with CMA-ES on one score, each against random search with the same budget, and then holds the chains it found to three phantoms not seen in tuning.

## 2. The Method

**The phantoms** (`eikon::phantoms`). A tree of round cones built as eikon's skeleton check builds it: depth 3 (15 tubes, 7 junctions, 8 ends), from (0, 0, −16) mm along (0.1, 0.05, 1), the root 10 mm long and of radius 2 mm, children's radii keeping Murray's law r³ = r₁³ + r₂³ with asymmetry 0.8 (radii 2.0 down to 0.68 mm), branching angle 0.7 rad, length ratio 0.9, each branching plane turned at random. Four balls of radius 2.5 to 3.5 mm are added, their centres uniform in ±14 mm and kept at least r + 2 mm from every tube and ball: blobs of the vessels' intensity that are not vessels, which a vesselness measure is meant to reject. Everything is rasterised with partial volumes (4³ sub-samples where a surface crosses a voxel) on a grid of 72³ voxels of 0.5 mm, as intensity 0 to 1, with Gaussian noise of deviation 0.5. The truth is the set of voxels whose centres lie inside the tubes (4291 voxels on phantom 1, 4137 on phantom 4). Each phantom is keyed by its seed: 1 and 2 for tuning, 3, 4 and 5 held out.

**The chain** (all eikon):
1. Frangi's vesselness (`features::frangi`, Frangi et al. 1998, eq. 13), the largest over the scales σ_lo … σ_hi of a ladder σ_k = 0.35 mm · 2^{k/4}, k = 0 … 12 (0.35 to 2.8 mm). The γ-normalised Hessian's eigenvalues at every scale of the ladder (`filters::hessian`, `features::hessian_eigenvalues`) are computed once per phantom, so an evaluation applies Frangi's measure to cached eigenvalues.
2. Hysteresis on the response (`segmentation::hysteresis`, 26-connected): the voxels above `low` joined to one above `high`.
3. Growth (`segmentation::grow`) from those voxels into the image smoothed at 0.5 mm, over voxels at or above `level`.

**The eight parameters**, from the unit cube: lo = round(12x₀) and hi = lo + round((12 − lo)x₁) on the ladder; α = 0.1 + 1.4x₂ and β = 0.1 + 1.4x₃; c = 10^{−3+3x₄}; high = 0.02 + 0.88x₅ and low = (0.05 + 0.95x₆)·high; level = 0.1 + 0.8x₇.

**The scores** (`fusion::dice`, `fusion::cldice`), each the mean over the two training phantoms. For NSGA-II, two objectives minimised: 1 − Dice and 1 − clDice. For CMA-ES, one score: 1 − (Dice + clDice)/2.

**The searches** (heuresis):
- NSGA-II (`multiobjective::Nsga2`), μ = 32, SBX (η 20, probability 0.9) and polynomial mutation (η 20, rate 1/8) on the box, 25 generations, 800 evaluations;
- random search (`search::RandomSearch`), uniform in the box, 800 evaluations on the same two objectives;
- CMA-ES (`strategies::Cma`) on the score, mean at the centre, σ₀ = 0.3, λ = 10, the box transformation, 320 evaluations, seeds 1 to 3;
- random search on the score, 320 evaluations, seeds 1 to 3.

Every generation is evaluated on a pool of eight threads (`search::Pool`); eikon's filters use threads of their own inside.

**The hand-set chain**, for reference: scales 0.5 to 2 mm (k = 2 … 10), Frangi's α = β = 0.5 and c half the largest Hessian norm over those scales (`features::frangi_c`, 0.274 on training phantom 1), hysteresis 0.1 to 0.3, growth at the half-maximum 0.5.

**Held out.** Every chain chosen (the hand-set chain; NSGA-II's best Dice, best clDice and best mean; CMA-ES's best of each seed; random search's best mean) and every non-dominated chain of NSGA-II's and random search's runs is applied to the three held-out phantoms. They are measured by Dice, clDice, Betti numbers (`fusion::betti`), surface distances (`fusion::surface_distances`) and the share of the balls' voxels kept.

## 3. The Code

- `tune_pipeline.cpp`, built by `CMakeLists.txt` against calx, morphe, physis, tyche, heuresis and eikon by relative path; `tune_pipeline results all` writes `results/`.
- `figures.py` draws the three figures from `results/` alone.
- Versions: heuresis 6b3cc88, eikon 10b108a, The Science Infrastructures 4153cb9. The runs were made on 2026-11-05 with g++ 11.5 (`-O3` by CMake's Release) in WSL on the development machine, eight threads, shared at the time with other work: 439 s in all, of which NSGA-II 86 s, random search 90 s, each CMA-ES run 31–36 s, and the five phantoms with their Hessians 13 s for each two.

```
cmake -S . -B ~/.cache/science-build/exp-tune-pipeline -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-tune-pipeline --target tune_pipeline -j4
~/.cache/science-build/exp-tune-pipeline/tune_pipeline results all
../../../.venv/Scripts/python.exe figures.py
```

## 4. The Parameters

| What | Value |
|---|---|
| Grid | 72³ voxels of 0.5 mm, centred |
| Tree | depth 3, root 10 mm × 2 mm radius, γ = 3, asymmetry 0.8, angle 0.7 rad, length ratio 0.9 |
| Balls | 4, radius 2.5–3.5 mm, at least r + 2 mm from the tree and each other |
| Image | intensity 0 (background) to 1, Gaussian noise σ 0.5, smoothed at 0.5 mm for growth |
| Phantoms | seeds 1, 2 (tuning); 3, 4, 5 (held out) |
| Scale ladder | 0.35 mm · 2^{k/4}, k = 0 … 12 |
| NSGA-II | μ 32, 25 generations, 800 evaluations; draws seed 41 |
| Random search (two objectives) | 800 evaluations; seed 42 |
| CMA-ES | λ 10, σ₀ 0.3, 320 evaluations; seeds 51–53 |
| Random search (score) | 320 evaluations; seeds 61–63 |

## 5. The Data

None. The phantoms are eikon's, made by this program from their seeds; no patient image is used.

## 6. The Results

### 6.1 Tuning Pays, and Random Search Is a Strong Rival

| Chain | Training Dice | Training clDice | Held-out Dice | Held-out clDice | Balls kept (held out) |
|---|---|---|---|---|---|
| hand-set | 0.608 | 0.950 | 0.665 | 0.950 | 89% |
| NSGA-II, best Dice | 0.932 | 0.988 | 0.933 | 0.989 | 0% |
| NSGA-II, best clDice | 0.864 | 0.994 | 0.577 | 0.660 | 0% |
| NSGA-II, best mean | 0.931 | 0.990 | 0.619 | 0.659 | 0% |
| CMA-ES, seed 1 | 0.933 | 0.990 | 0.620 | 0.659 | 0% |
| CMA-ES, seed 2 | 0.932 | 0.988 | 0.932 | 0.988 | 0% |
| CMA-ES, seed 3 | 0.931 | 0.990 | 0.886 | 0.985 | 12% |
| random search, best mean | 0.930 | 0.990 | 0.932 | 0.988 | 0% |

(Means over the phantoms of each set; `results/held_out.csv` holds each phantom.) The hand-set chain finds the whole tree (clDice 0.95) but grows into the balls: it keeps 88–93% of their voxels and falls into 5 or 6 components. Every tuned chain rejects the balls on the training phantoms, raising Dice from 0.61 to 0.93.

On the single score, CMA-ES's three runs reached 0.0388, 0.0398 and 0.0395 in 320 evaluations, random search's 0.0414, 0.0441 and 0.0398. On the two objectives, NSGA-II's set and random search's have hypervolumes (to Dice = clDice = 0) of 0.9262 and 0.9246 after 800 evaluations each. The searches beat random search by little because the landscape is three plateaus. Of random search's 1760 evaluations, 46% keep nothing (score 1), 19% lie on the hand-set chain's plateau, the tree with the balls (scores 0.19–0.26), and 0.5% reach below 0.05. CMA-ES's evaluations straddle the plateaus (22% nothing, 47% with the balls, 1.8% below 0.05). NSGA-II concentrates on the good one (38% of its evaluations below 0.05).

### 6.2 The Front Is a Corner

The trade-off between the objectives is small once the balls are rejected: NSGA-II's non-dominated set spans Dice 0.864 to 0.932 at clDice 0.994 to 0.988. A wider growth (level 0.29 against 0.39 to 0.43) keeps more of every branch, at the price of a thicker reconstruction. The 89 non-dominated chains of NSGA-II's run share only 7 distinct pairs of scores: the scores are step functions of the parameters, and the population filled the plateaus of the front with many different chains.

### 6.3 What Does Not Carry

Held to the three phantoms not seen, 60 of NSGA-II's 89 non-dominated chains keep nothing at all on held-out phantom 4 (Dice 0). Fifty-eight of them reconstruct phantoms 3 and 5 as well as the training ones (Dice 0.86 to 0.93); the other two keep nothing on phantom 5 either. Three of the seven chosen chains fail so: NSGA-II's best clDice and best mean, and CMA-ES's seed 1. All of NSGA-II's front uses only the coarsest scales, σ from 1.98 or 2.35 mm upwards, where the balls answer the vesselness filter too, and a high threshold of 0.098 to 0.170. Random search's 5 non-dominated chains use finer scales (σ from 0.59 to 1.98 mm) and higher thresholds (0.21 to 0.72), and all of them carry: their worst held-out Dice is 0.86 to 0.93.

The searches found what the two training phantoms reward, a window of thresholds at coarse scales between the trees' response and the balls'. Phantom 4's tree evidently falls outside that window, since nothing of it is kept. The more a search exploits its training set, the more closely it fits such a window. Random search, which exploits nothing, carried every chain it chose. Its chains are no better on the training phantoms; they are less specific to them.

## 7. The Figures

| Figure | What it shows |
|---|---|
| [search-47](../../../assets/gallery/search-47-vessel-chain-front.jpg) | Every evaluation of NSGA-II and random search in (Dice, clDice), their non-dominated chains, the hand-set chain and CMA-ES's bests; each non-dominated chain's training Dice against its worst held-out Dice |
| [search-48](../../../assets/gallery/search-48-vessel-reconstructions.jpg) | Held-out phantom 4: a slice, its maximum projection, the truth, and four chains' reconstructions (kept, kept wrongly, missed) |
| [search-49](../../../assets/gallery/search-49-vessel-chain-score.jpg) | CMA-ES against random search on the score; how each algorithm's evaluations fall among the plateaus |

`gallery.md` holds them with their captions.

## 8. The Caveats

- **Two training phantoms.** The failure to carry is a property of tuning on two images; more phantoms, or an objective over the worst phantom, or a margin on the thresholds, would be the remedies to test next. Three held-out phantoms measure it coarsely: one phantom decides it.
- **Phantoms are not patients.** The tubes are perfect round cones with sharp walls, the noise white and Gaussian, the contrast uniform, the distractors perfect balls. Real images have partial voluming at different resolutions, inhomogeneous contrast, touching structures and noise that is neither white nor Gaussian. The numbers here are of this phantom family only.
- **Quantised parameters.** The scales are chosen from a ladder of 13 (σ rounded to it), so the scores have plateaus in two of the eight coordinates by construction, besides those that thresholds make. CMA-ES and SBX see the plateaus as flat regions.
- **One run of NSGA-II and of random search** on the two objectives, three of CMA-ES. The hypervolumes and best scores are single draws.
- **Growth's smoothing** (0.5 mm) and the 26-connectivity are fixed, not tuned.
- **clDice thins both masks** at every evaluation (`fusion::cldice`), which is most of an evaluation's cost.

## 9. References

- Frangi, A. F., Niessen, W. J., Vincken, K. L. and Viergever, M. A. (1998), Multiscale vessel enhancement filtering, *MICCAI*, LNCS 1496, 130–137.
- Deb, K., Pratap, A., Agarwal, S. and Meyarivan, T. (2002), A fast and elitist multiobjective genetic algorithm: NSGA-II, *IEEE Trans. Evol. Comput.* 6, 182–197.
- Hansen, N. (2016), The CMA evolution strategy: a tutorial, arXiv:1604.00772.
- Lamy, J., Merveille, O., Kerautret, B. and Passat, N. (2022), A benchmark framework for multiregion analysis of vesselness filters, *IEEE Trans. Med. Imaging* 41, 3649–3662.
- Murray, C. D. (1926), The physiological principle of minimum work: I. The vascular system and the cost of blood volume, *PNAS* 12, 207–214.
- Shit, S., Paetzold, J. C., Sekuboyina, A. et al. (2021), clDice: a novel topology-preserving loss function for tubular structure segmentation, *CVPR*, 16560–16569.
