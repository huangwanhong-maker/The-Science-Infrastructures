---
title: Tropical Cyclones Found in a Field of Pressure
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/cyclones.js
  - applications/geography/earth-atlas/src/layers/tracks.js
references: [UllrichZarzycki2017, ZarzyckiUllrich2017, Holland1980, Emanuel2003]
---

# Tropical Cyclones Found in a Field of Pressure

## 1. Recognition as a Definition

A forecast model does not announce a typhoon. It produces fields: pressure, wind and temperature on a grid, every few hours. A tropical cyclone is a *feature* of these fields, and to compare the cyclones of several models, or of a model and the world, one must first say what counts as one.

A tracker is a map from a sequence of fields to a finite set of curves in space and time. Its definition has thresholds, and the thresholds are conventions: Zarzycki and Ullrich (2017) show that the number of cyclones a tracker finds in a climate simulation changes by tens of per cent between reasonable choices. A track is therefore reported together with the definition that found it. The definition used here follows TempestExtremes (Ullrich and Zarzycki 2017), with one substitution (§3.2).

## 2. The Pressure Deficit and the Wind

*Hypotheses.* An axisymmetric vortex in the steady state, near enough to its centre that the Coriolis force is small beside the centripetal acceleration (cyclostrophic balance), with air of constant density ρ:

  v² / r = (1/ρ) ∂p/∂r.  (1)

**Proposition 1 (modified Rankine vortex).** Let v = v_m r/r_m for r ≤ r_m, and v = v_m (r_m/r)^x for r > r_m, with 0 < x ≤ 1. Then the deficit of central pressure is

  Δp = p(∞) − p(0) = ρ v_m² (1 + x) / (2x).  (2)

*Proof.* Integrate (1) from 0 to ∞. Inside, ∫₀^{r_m} v_m² r / r_m² dr = v_m²/2. Outside, ∫_{r_m}^∞ v_m² r_m^{2x} r^{−2x−1} dr = v_m²/(2x). ∎

For a solid core and a potential vortex outside (x = 1), Δp = ρ v_m². For the flatter profiles that are observed (x ≈ 0.5), Δp = 1.5 ρ v_m², so a deficit of 50 hPa with ρ = 1.15 kg m⁻³ goes with a greatest wind of about 54 m s⁻¹. Holland (1980) gives the same relation in a smoother family of profiles, v_m = √(B Δp / (ρ e)) with B between 1 and 2.5.

*Remark.* Pressure is the quantity to detect, because it is smooth where the wind is not. By (2), half of the deficit (x = 1) or a third of it (x = 0.5) lies inside the radius of greatest wind, typically 20 to 50 km. A grid of 1° (about 110 km) cannot resolve that core, so a model's, or a product's, central pressure is higher than the cyclone's. The intensities in §6 are those of fields at 1°, and are compared with one another, not with the observed intensity.

## 3. Detection at One Time

Let p be the sea-level pressure on a regular longitude–latitude grid G, let d be the great-circle angle between places, and let γ_θ(s) be the place reached from x after an arc s at bearing θ (the direct problem on the sphere).

**Definition 1 (candidate).** A grid point x is a *candidate* if

1. p(x) ≤ p_max and |φ(x)| ≤ φ_max;
2. p(x) ≤ p(y) for every grid point y with d(x, y) ≤ r_m, ties going to the first point in grid order;
3. (*closed contour by rays*) for each bearing θ_k = 2πk/K, k = 0, …, K − 1, there is an s ∈ (0, R] with p(γ_{θ_k}(s)) ≥ p(x) + Δ.

The values are r_m = 6°, R = 5.5°, Δ = 200 Pa, K = 16, p_max = 1005 hPa and φ_max = 50°. Condition 2 makes any two candidates at least r_m apart.

### 3.1 Rays and the Flood Fill

TempestExtremes states the closed contour differently: the connected component C(x) of {y : p(y) < p(x) + Δ} that contains x lies within the disc of radius R about x. Call this *closed by flood fill*.

**Proposition 2.** Closed by flood fill implies closed by rays. The converse is false.

*Proof.* Each ray from x is a connected path from x out to distance R. Since C(x) lies within the open disc of radius R, the ray leaves C(x) at some s ≤ R, and where it leaves, p ≥ p(x) + Δ. For the converse, let a narrow channel of low pressure run outward between two adjacent rays. Every ray then meets a rise, while the component escapes the disc. ∎

The ray test is the weaker, and so admits more candidates. With K = 16, adjacent rays are 2π sin R / K ≈ 2.2° apart at the outer circle, two cells of a 1° grid, so a channel narrower than that can slip through. In practice such a channel is a trough joining two systems, and condition 2 together with the filters of §4 removes most of what it admits. Along each ray the field is read every half cell, at the nearest grid point.

### 3.2 The Warm Core, Replaced

TempestExtremes also asks for a warm core: the thickness between 300 and 500 hPa must fall away from the centre. Not every product the atlas shows carries both levels. The requirement is replaced by p ≤ p_max for a candidate, and by filters on the track's genesis and depth (§4). These keep out the summer heat lows over land and the lows of the westerlies. They do not keep out a monsoon depression over the sea, some of which a warm-core test would reject.

### 3.3 The Centre Between Grid Points

The minimum is placed between grid points by fitting, in each direction, the parabola through the point and its two neighbours. With f_x = (p₊ − p₋)/2 and f_xx = p₊ − 2p₀ + p₋ in grid units, the offset is −f_x / f_xx, limited to half a cell, and the central value is p₀ − f_x² / (2 f_xx), summed over the two directions.

**Proposition 3.** The fit is exact for any p that is a separable quadratic, p = a + b(ξ − ξ₀)² + c(η − η₀)² with b, c > 0 and |ξ₀|, |η₀| ≤ ½ in grid units about the point.

*Proof.* The neighbours lie at ξ = ±1. For a quadratic the centred differences are exact, f_x = −2bξ₀ and f_xx = 2b, and the terms in η cancel from them. Hence the offset is ξ₀, and likewise η₀; and p₀ − f_x²/(2f_xx) − f_y²/(2f_yy) = (a + bξ₀² + cη₀²) − bξ₀² − cη₀² = a. ∎

A Gaussian low is not a quadratic. On a 1° grid, the test of a 933 hPa low of radius 2.5° places it within 0.02° and 0.5 hPa.

## 4. Linking in Time

**Definition 2 (track).** Given the candidates N_i at times t_i, a *track* is a sequence n_i, n_{i+1}, …, n_j with n_k ∈ N_k, each candidate in at most one track, and

  d(n_k, n_{k+1}) ≤ v_s (t_{k+1} − t_k),  v_s = 8° per 6 h.  (3)

The bound is 41 m s⁻¹. Tropical cyclones usually move at 5 to 10 m s⁻¹, and faster when they recurve into the westerlies.

At each time, the open tracks and the new candidates form a bipartite graph whose edges satisfy (3). Edges are taken in increasing order of length, skipping those whose ends are already taken. The result is a maximal matching, but not necessarily the one of least total length that the Hungarian algorithm would find. Since candidates are at least r_m = 6° apart, two tracks rarely compete for one candidate. A track not continued ends. A candidate not taken begins a track.

A track is kept if

1. it lasts at least 24 h;
2. it deepens to 1000 hPa or below;
3. it begins within 25° of the equator, as tropical cyclones form over warm water away from the equator (Emanuel 2003);
4. it begins over the sea, judged by the relief (elevation below zero in ETOPO 2022), when a land–sea test is given.

## 5. Positions Between Times, and Names

A track is drawn at the time shown by linear interpolation, in time, of longitude (the short way across the date line) and latitude between the two positions around it. This is not the great circle. For a step of 3° along a parallel, the great circle's midpoint lies 0.005° to 0.010° poleward of the linear one at latitudes 15° to 45° (0.5 to 1.1 km), and for 6° four times as much. Both are far below the uncertainty of a forecast position.

A model's track takes the name of an official track (GDACS, relaying the warning centres) if their positions at some common time are within 5°.

## 6. Verification

| Test | Agreement |
|---|---|
| A Gaussian low of 933 hPa on a 1° grid: position and central pressure | 0.02°; 0.5 hPa |
| One low gives one candidate | Exact |
| A long trough, open at its ends | No candidate |
| A cyclone moving 3° in 6 h, over two days | One track of 9 positions |
| A low born at 35° N; a low that never deepens to 1000 hPa; a low born over land | Each left out |
| Interpolation across the date line | Exact |

These are tests of the definition, not of the world. The comparison with observed cyclones is in the Earth Atlas, where model tracks are drawn beside the observed track and the official forecast.

*Caveat.*
- The thresholds are TempestExtremes' values for grids of about 1°. On finer grids, more and smaller lows pass the contour test.
- Without the warm-core test, a strong monsoon depression or a subtropical low may be taken for a tropical cyclone.
- Intensities at 1° are not the cyclone's (§2).
