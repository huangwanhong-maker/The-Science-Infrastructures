---
title: The Shallow-Water Equations on the Sphere, and the Initialisation of a Forecast
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/shallowwater.js
  - applications/geography/earth-atlas/src/core/spectral.js
  - applications/geography/earth-atlas/src/workers/labmodel.worker.js
references: [HackJakob1992, Williamson1992, LynchHuang1992, Lamb1932, Cressman1958, Durran2010]
---

# The Shallow-Water Equations on the Sphere, and the Initialisation of a Forecast

## 1. One Step Beyond the Barotropic Model

The barotropic vorticity equation (see [the barotropic note](barotropic-vorticity-model.md)) keeps only the advection of vorticity by a non-divergent flow. Its gravest fault is in the longest waves: their Rossby phase speed, −2Ω/(n(n + 1)) per radian of longitude, is large for small n, so planetary waves run westward much faster than the atmosphere's do. Cressman (1958) traced this to the missing divergence.

The shallow-water equations restore it. A single layer of fluid, of depth h and geopotential Φ = g h, may converge and diverge. Its surface rises and falls as it does, so the model carries gravity waves as well as Rossby waves. And its Rossby waves feel a finite radius of deformation, L = √Φ / f, beyond which the free surface, not the vorticity, takes up the motion. For the 500 hPa surface, √Φ is about 230 m s⁻¹ and L about 2300 km at mid-latitudes. Waves longer than L are slowed, which is the correction the barotropic model lacks.

## 2. The Equations

*Hypotheses.*
- A single layer of homogeneous, incompressible fluid, in hydrostatic balance, on a rotating sphere of radius a.
- Its horizontal velocity V is independent of height.
- There is no friction, heating or topography.

**Definition 1.** With η = ζ + f the absolute vorticity, δ = ∇ · V the divergence, and V = **k** × ∇ψ + ∇χ (so that ζ = ∇²ψ and δ = ∇²χ), the shallow-water equations in vorticity–divergence form are (Hack and Jakob 1992, eqs. 2.10–2.12)

  ∂ζ/∂t = −∇ · (η V),  (1)
  ∂δ/∂t = **k** · ∇ × (η V) − ∇²(Φ + |V|²/2),  (2)
  ∂Φ/∂t = −∇ · (Φ V).  (3)

*Derivation.* The momentum equation in vector-invariant form is ∂V/∂t = −η **k** × V − ∇(Φ + |V|²/2). Taking **k** · ∇× gives (1), since **k** · ∇ × (η **k** × V) = ∇ · (η V). Taking ∇· gives (2), since ∇ · (**k** × W) = −**k** · ∇ × W. Equation (3) is the conservation of mass. ∎

**Proposition 1.** Over the sphere, (1)–(3) conserve the mass ⟨Φ⟩ and the total energy ⟨Φ|V|²/2 + Φ²/2⟩ (Williamson et al. 1992, §2).

**Proposition 2 (gravity waves; Lamb 1932, §307).** On a non-rotating sphere, linearise about rest and a uniform depth H. A height perturbation proportional to a spherical harmonic of degree n then oscillates with frequency

  ω_n = √(g H n(n + 1)) / a.  (4)

*Proof.* Linearised, (2) and (3) read ∂δ/∂t = −∇²Φ′ and ∂Φ′/∂t = −gH δ. Hence ∂²Φ′/∂t² = gH ∇²Φ′, and ∇² Y_n = −n(n + 1)/a² Y_n. ∎

## 3. The Numerical Method

The model holds ζ, δ and Φ as spherical-harmonic coefficients at triangular truncation T42. It evaluates the products η V, Φ V and |V|²/2 on a Gaussian grid of 128 × 64 points, which is fine enough that quadratic products alias nothing. It transforms them back by the divergence and curl operators of the spectral transform method.

The transforms in longitude are fast Fourier transforms of 128 points (Cooley and Tukey 1965). This halves the cost of a step against direct sums.

Steps are of fourth-order Runge–Kutta, 15 minutes long. On the imaginary axis, the method is stable for |ω Δt| < 2.83. The fastest wave the model holds is the gravity wave of degree 42, with ω ≈ 1.5 × 10⁻³ s⁻¹ for √Φ ≈ 230 m s⁻¹. A step of 900 s keeps |ω Δt| below 1.4.

The smallest scales are damped by ∇⁴, as in the barotropic model.

## 4. Initialisation by a Digital Filter

An analysis is not in exact balance. Started from one, the model rings with gravity waves of periods from minutes to hours. Their amplitude in height is small, but their tendencies are not. Lynch and Huang (1992) proposed removing them with a filter in time, applied to a short integration of the model itself.

*Definition 2 (the filter).* Let N steps of length Δt span half of the interval T_s. The weights of a low-pass filter with cut-off period T_c are

  h_k = σ_k · sin(k θ_c) / (k π), with h_0 = θ_c/π, θ_c = 2π Δt / T_c and k = −N, …, N.

The Lanczos factors σ_k = sin(kπ/(N + 1)) / (kπ/(N + 1)) damp the ringing of the truncated series (the Gibbs phenomenon). The weights are normalised so that Σ h_k = 1.

*The procedure.*
1. Integrate the model backward N steps, without damping, from the analysis.
2. Integrate forward 2N steps, accumulating Σ h_k x_k over the states x_k.
3. The sum is valid at the time of the analysis. It keeps the slow, balanced motion and removes oscillations faster than T_c.

The atlas uses T_s = T_c = 6 h.

## 5. Verification

The tests run the model at T21 (`tests/run-tests.mjs`).

| Test | Result |
|---|---|
| Williamson et al. (1992), case 2: steady nonlinear geostrophic flow, after one day | Height error 7 × 10⁻¹⁵ (relative, l₂) |
| The same about an axis tilted by 45°, so that the flow crosses the grid's poles | 1 × 10⁻¹⁴ |
| Mass, over that day | Kept to 6 × 10⁻¹⁶ |
| Total energy, over a day of unbalanced flow | Kept to 6 × 10⁻⁹ |
| Gravity wave of degree 5, H = 1 km (period 20.5 h by (4)) | Reversed after half a period and restored after a whole one, to 3 × 10⁻⁶ |
| Digital filter: gravity waves of degree 20 (1.7 h) and degree 2 (46 h) | Amplitudes 0.003 and 0.99 after filtering |

Case 2 is exact because its height and stream function are polynomials of degree two on the sphere, which any truncation holds exactly. The error is that of rounding.

## 6. A Forecast Verified

Both laboratory models were started from the GFS analysis of 2026-10-19 00:00 UTC. They are compared with the GFS's own analyses (forecast beyond the present) over 20–80° N.

| Lead | Barotropic | Shallow water | Persistence |
|---|---|---|---|
| +24 h | 58 m | 43 m | 81 m |
| +48 h | 85 m | 74 m | 108 m |
| +72 h | 108 m | 100 m | 127 m |

The shallow-water model is better at every lead. Two reasons are likely. First, it starts from the analysed height itself, rather than from heights diagnosed by linear balance. Second, its radius of deformation slows the longest waves.

*Caveat.*
- One case is an anecdote, not a statistic. A fair comparison would average many starts.
- Both models remain far from operational skill: about 10 m at a day and 30 m at three days. The difference is what baroclinic dynamics, moisture, topography and data assimilation add.
