---
title: The Barotropic Vorticity Equation, and the First Numerical Forecast
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/spectral.js
  - applications/geography/earth-atlas/src/core/barotropic.js
  - applications/geography/earth-atlas/src/layers/labmodel.js
references: [Charney1950, Haurwitz1940, Williamson1992, Durran2010, HoltonHakim2013]
---

# The Barotropic Vorticity Equation, and the First Numerical Forecast

## 1. The Model of 1950

In March 1950, Charney, Fjørtoft and von Neumann ran the first numerical weather forecasts on the ENIAC. They forecast the 500 hPa flow over North America for 24 hours, and each forecast took about a day of computing. Their model was the barotropic vorticity equation. It keeps, of the dynamics of the atmosphere, only the advection of absolute vorticity by a non-divergent flow.

The atlas runs the same equation from today's analyses, in a few seconds, and verifies it the way they did: against what happened, and against persistence (the forecast that nothing changes).

## 2. The Equation

*Hypotheses.*
- The flow at the level of non-divergence, near 500 hPa, is horizontal and non-divergent.
- The flow is barotropic: it has no vertical shear to release energy.
- Friction and heating are neglected.

**Definition 1.** With V = **k** × ∇ψ the non-divergent velocity, ζ = ∇²ψ its relative vorticity, and f = 2Ω sin φ the Coriolis parameter, the barotropic vorticity equation is

  ∂ζ/∂t = −V · ∇(ζ + f) = −∇ · ((ζ + f) V).  (1)

The second form follows from ∇ · V = 0.

**Proposition 1.** Over the whole sphere, (1) conserves the kinetic energy E = ½ ⟨|∇ψ|²⟩ and the enstrophy Z = ½ ⟨ζ²⟩, where ⟨·⟩ is the mean over the sphere.

*Proof (sketch).* Write q = ζ + f and J(a, b) = **k** · (∇a × ∇b), so that V · ∇q = J(ψ, q) and (1) reads ∂ζ/∂t = −J(ψ, q). On a closed surface the integral of any Jacobian vanishes, and J(F(a), b) = F′(a) J(a, b).
- *Energy.* Because E = −½⟨ψζ⟩, we have dE/dt = −⟨ψ ∂ζ/∂t⟩ = ⟨ψ J(ψ, q)⟩ = ⟨J(½ψ², q)⟩ = 0.
- *Enstrophy.* ⟨q ∂q/∂t⟩ = −⟨q J(ψ, q)⟩ = −⟨J(ψ, ½q²)⟩ = 0, so the absolute enstrophy ½⟨q²⟩ is conserved. It differs from Z by ⟨fζ⟩ + ½⟨f²⟩. Since f ∝ μ is a spherical harmonic of degree 1, ⟨fζ⟩ = ⟨ψ∇²f⟩ = −(2/a²)⟨fψ⟩, which is proportional to the axial angular momentum of the flow, which (1) also conserves (Durran 2010, §6.4.2). ∎

*Remark.* Energy and enstrophy cascade in opposite directions (Fjørtoft 1953): energy to large scales, enstrophy to small ones. A truncated model must remove the enstrophy that reaches its smallest scales, or it accumulates there. The atlas uses a damping proportional to (n(n + 1))², which acts only near the truncation.

## 3. The Spectral Transform Method

- **Representation.** ζ is held as coefficients of spherical harmonics Y_n^m = P̄_n^m(μ) e^{imλ}, with μ = sin φ and triangular truncation n ≤ M = 42 (T42).
- **Laplacian.** It is diagonal: ∇² Y_n^m = −n(n + 1)/a² · Y_n^m. Inverting ζ = ∇²ψ is therefore exact.
- **The product (ζ + f)V.** It is formed on a Gaussian grid of 128 × 64 points. The grid is chosen so that Gaussian quadrature integrates products of two truncated fields exactly ("quadratic" grid, nlat ≥ (3M + 1)/2), which leaves the products free of aliasing.
- **Back to coefficients.** The divergence of the product is transformed back by quadrature. A term in ∂/∂μ is moved onto the basis function by integration by parts (Hack and Jakob 1992, §3).
- **Time stepping.** The classical fourth-order Runge–Kutta method, with steps of 30 minutes.

*Verification* (tests/run-tests.mjs):
- Transforms return their coefficients to 1 × 10⁻¹⁴.
- The curl of the wind of ψ is ∇²ψ, and its divergence is 0, to 6 × 10⁻¹⁵.

## 4. An Exact Solution for Testing

**Proposition 2** (Haurwitz 1940). The stream function

  ψ = −a² ω μ + a² K (1 − μ²)^{R/2} μ cos Rλ  (2)

solves (1) exactly. The pattern turns eastward without change of shape at the angular velocity

  ν = [R(3 + R) ω − 2Ω] / [(1 + R)(2 + R)].  (3)

For R = 4 and ω = K = 7.848 × 10⁻⁶ s⁻¹ (Williamson et al. 1992, test 6), ν is 12.2° a day.

*Verification.* After one day, the model's vorticity differs from (2) by 7 × 10⁻¹⁰ of its size, and its phase from (3) by 7 × 10⁻¹⁰ rad. Over a day of random flow, energy changes by 2 × 10⁻⁶ and enstrophy by 3 × 10⁻⁵; both losses come from the Runge–Kutta step at the smallest scales.

## 5. From Stream Function to Height

Forecasters read the 500 hPa flow as heights. The linear balance equation (Holton and Hakim 2013, eq. 6.62),

  g ∇²z = ∇ · (f ∇ψ),  (4)

gives the height consistent with a stream function.

- **Check.** For the zonal flow u = u₀ cos φ, (4) gives z = −(Ω a u₀ / g)(μ² − ⅓) + constant, which the tests reproduce to 5 × 10⁻¹³ m.
- **Initial balance.** The analysed heights are not exactly balanced. The atlas therefore adds the *change* of the balanced height to the analysed height at the start, so a forecast begins exactly at the analysis.

## 6. A Forecast Verified

The atlas starts the model from the GFS analysis of the 500 hPa wind and compares its heights with the GFS's own analyses at later times, over 20–80° N. Persistence is the baseline.

A run from 2026-10-19 00:00 UTC gave the following root-mean-square errors:

| Lead | Model | Persistence |
|---|---|---|
| +24 h | 58 m | 81 m |
| +48 h | 85 m | 108 m |
| +72 h | 108 m | 127 m |

These are errors of the right size for a barotropic model. Operational models today make errors of about 10 m at a day and 30 m at three days. The difference measures what baroclinic dynamics, moisture and data assimilation add.

The shallow-water model, started from the same analysis, does better at every lead: see [its note](shallow-water-model.md).

*Caveat.*
- Over the tropics, where f is small and the flow divergent, the model has little skill. Linear balance there says little about heights.
- The model also lacks the "equivalent barotropic" divergence term that later versions added (Cressman 1958), so very long waves drift westward too fast.
