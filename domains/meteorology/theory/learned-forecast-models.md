---
title: Forecasting by a Learned Map: Pangu-Weather Run in the Laboratory
status: draft
implemented-by:
  - data/products/ai_forecast.py
references: [Bi2023, Hersbach2020, Liu2021Swin]
---

# Forecasting by a Learned Map: Pangu-Weather Run in the Laboratory

## 1. The Forecast as a Map

Let X be the space of states of the atmosphere on a grid. For Pangu-Weather (Bi et al. 2023), a state is five variables (geopotential, specific humidity, temperature and the two components of the wind) at 13 pressure levels, and four at the surface (pressure at sea level, the wind at 10 m and temperature at 2 m), on a grid of 0.25°:

  X = ℝ^{5 × 13 × 721 × 1440} × ℝ^{4 × 721 × 1440},  dim X = 71 638 560.  (1)

The equations of motion define a flow, Φ_τ : X → X, the state τ hours later. A numerical model approximates Φ_τ by discretising the equations, as the laboratory's [barotropic](barotropic-vorticity-model.md) and [shallow-water](shallow-water-model.md) models do. A learned model instead fits a map F_τ ≈ Φ_τ to pairs of states (x(t), x(t + τ)) drawn from a reanalysis, here ERA5 (Hersbach et al. 2020), 1979–2017.

The two approaches fail differently.
- A numerical model inherits the conservation laws and the dispersion relations of its equations, and fails where its resolution or its physics does.
- A learned model has no equations, so nothing guarantees that it conserves mass or energy. It is as good as its training data allow, and it inherits their faults: ERA5's tropical cyclones are too weak, and so are Pangu-Weather's (Bi et al. 2023).

## 2. The Semigroup and Hierarchical Temporal Aggregation

The true flow is a semigroup: Φ_{s+t} = Φ_s ∘ Φ_t, and Φ_0 = id. A learned family need not be. F_24 and F_6⁴ are fitted separately and differ, and each application of a learned map adds error that the next application carries forward. A forecast of T hours should therefore use as few applications as possible.

Pangu-Weather provides F_1, F_3, F_6 and F_24, and composes T greedily: as many 24 h steps as fit, then 6 h, 3 h and 1 h steps (Bi et al. 2023, Methods).

**Proposition 1.** For the step set {1, 3, 6, 24}, the greedy decomposition of any T ∈ ℕ uses the fewest steps.

*Proof.* Let T = c₁ + 3c₃ + 6c₆ + 24c₂₄ use the fewest steps. Then c₁ ≤ 2, c₃ ≤ 1 and c₆ ≤ 3, for otherwise three 1 h steps, two 3 h steps or four 6 h steps could be replaced by a single step. Hence c₁ + 3c₃ + 6c₆ ≤ 2 + 3 + 18 = 23 < 24, so c₂₄ = ⌊T/24⌋. In the same way, c₆ = ⌊(T mod 24)/6⌋, and so on. These bounds make (c₁, c₃, c₆, c₂₄) the digits of T in the mixed radix (3, 2, 4), which are unique, and they are what the greedy procedure produces. ∎

*Remark.* The proof uses only the fact that each step divides the next, so that the step set is a chain under divisibility. For a set like {1, 4, 6}, greedy fails: 8 = 6 + 1 + 1 by greedy, three steps, but 4 + 4 is two.

The laboratory uses F_24 and F_6, and forecasts every 6 hours: a lead of T = 24a + 6b (0 ≤ b ≤ 3) takes a + b steps. Each state at a multiple of 24 h is kept, and the leads between it and the next are reached from it by F_6. A forecast to 168 h with outputs every 6 h takes 7 + 7 × 3 = 28 steps.

## 3. The Architecture, Briefly

Pangu-Weather is a transformer of the Swin family (Liu et al. 2021), arranged in three dimensions: pressure level, latitude and longitude (Bi et al. 2023).

- The upper-air state is cut into patches of 2 × 4 × 4 points (levels × latitude × longitude) and the surface into 4 × 4, and each patch becomes a token.
- Attention is computed within windows of tokens and the windows shift between layers, so that information crosses their edges. The encoder has 2 + 6 layers, with one step of down-sampling between them, and the decoder mirrors it.
- The grid is not translation-invariant on the sphere. A cell near the pole is smaller than one at the equator, and a level near the ground behaves differently from one in the stratosphere. Swin's relative positional bias is therefore replaced by an *Earth-specific* bias, indexed by the absolute latitude and level of the window. The grid is periodic in longitude only, and the windows wrap there.

The model has about 256 million parameters and was trained on 192 GPUs. Each forecast step costs about 1.4 s on one such GPU, or 85 s on the eight CPU threads of the laboratory's computer.

## 4. An Initial State from Another System

The model was trained on ERA5 states, but the laboratory starts it from the GFS analysis (NOAA), at 0.25° and the same 13 levels. The two analyses differ in their models, their assimilation and their biases. A learned map applied to states slightly unlike those it was fitted on degrades by an amount nobody can derive. The laboratory therefore measures it directly (§5).

GFS gives geopotential height, which is converted to geopotential by Φ = g₀ Z with g₀ = 9.80665 m s⁻² (WMO's standard gravity, which defines geopotential height).

## 5. A Forecast Verified

### 5.1 The Height of the 500 hPa Surface

The forecast from 2026-10-19 00:00 UTC was compared with the GFS analyses at the leads that have since passed. It took 12 steps and 17 minutes. The measure is the root-mean-square error of the 500 hPa height over 20–80° N on a 1° grid, as for the laboratory's other models. At +72 h the reference is the GFS's forecast, as no analysis yet existed when the data were fetched.

| Lead | Pangu-Weather | Shallow water | Barotropic | Persistence |
|---|---|---|---|---|
| +6 h | 3.6 m | — | — | 29 m |
| +24 h | 6.6 m | 43 m | 58 m | 81 m |
| +48 h | 12.7 m | 74 m | 85 m | 108 m |
| +72 h | 18.3 m | 100 m | 108 m | 126 m |

At three days the learned model's error is a fifth of the shallow-water model's and a seventh of persistence's. This is the gap between one layer of fluid and a model of the whole column trained on 39 years of the atmosphere.

### 5.2 A Typhoon's Track

Typhoon Champi was then crossing the Philippine Sea. The cyclone tracker ([its note](tropical-cyclone-tracking.md)) finds it in the hindcast's sea-level pressure. Its positions are compared with the positions observed by the warning centres, as GDACS relays them.

| Lead | Distance from the observed centre | Central pressure at 1° |
|---|---|---|
| +24 h | 49 km | 973 hPa |
| +48 h | 58 km | 971 hPa |
| +60 h | 79 km | 968 hPa |

The track is good: within 80 km at two and a half days. The intensity is not. At 21 October 18:00 UTC the GFS, sampled on the same 1° grid, gives 921 hPa at the centre, and the hindcast about 967 hPa. This agrees with what Pangu-Weather's authors report: its tracks are good, and its cyclones too weak, as ERA5's are (Bi et al. 2023).

*Caveat.*
- The reference is the GFS's own analysis, which is also the initial state, so the measure favours models close to the GFS. ERA5 or another centre's analysis would be a fairer judge.
- One start is an anecdote. Skill varies from day to day, and only an average over many starts measures it.
- Errors at 500 hPa say little about rainfall, cyclone intensity or the surface, where a learned model's smoothing is greatest.

## 6. Use and Licence

The model's parameters are licensed CC BY-NC-SA 4.0. The laboratory uses them for learning and research, never commercially. They are fetched into its cache and not redistributed, and every forecast product says so in its source.
