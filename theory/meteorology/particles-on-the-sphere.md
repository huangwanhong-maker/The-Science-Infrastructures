---
title: Particles Carried by a Flow on the Sphere
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/flow.js
references: [Hairer2006, Hairer1993]
---

# Particles Carried by a Flow on the Sphere

## 1. The Equation of Motion

*Hypotheses.* The Earth is the sphere S² of radius R. A horizontal velocity field is given by its eastward and northward components u(λ, φ, t) and v(λ, φ, t), in m s⁻¹.

**Definition 1.** At a point **p** ∈ S² (a unit vector, Earth-fixed) with longitude λ and latitude φ, the local unit vectors are

  **e** = (−sin λ, cos λ, 0),  **n** = (−sin φ cos λ, −sin φ sin λ, cos φ).  (1)

Together with **p** they form an orthonormal basis. The velocity as a vector tangent to the sphere is **w** = u **e** + v **n**.

**Definition 2.** A particle carried by the flow is a solution of

  d**p**/dt = **w**(**p**, t) / R,  **p** ∈ S².  (2)

*Remark.* Writing (2) in λ and φ gives dλ/dt = u / (R cos φ), which is singular at the poles although the flow is not. The form (2), in the ambient space ℝ³, has no such singularity. **e** and **n** are undefined at the poles, but the vector **w** is defined there, as a limit, wherever the field is continuous.

## 2. A Step That Stays on the Sphere

The midpoint rule for d**p**/dt = **f**(**p**) is

  **m** = **p** + (h/2) **f**(**p**),  **p**′ = **p** + h **f**(**m**),  (3)

which is second-order accurate (Hairer, Nørsett and Wanner 1993, §II.1). Rule (3) leaves the sphere: |**p**′| ≠ 1. The atlas uses a variant that stays on it.

1. Take a half step and project it back: **m** = normalise(**p** + (h / 2R) **w**(**p**)).
2. Sample the velocity at **m** and carry it back to **p** by parallel transport along the great circle from **m** to **p**:

  **t** = **w**(**m**) − ((**w**(**m**) · **p**) / (1 + **m** · **p**)) (**m** + **p**).  (4)

3. Turn **p** about the axis **p** × **t̂** by the angle θ = |**t**| h / R:

  **p**′ = cos θ **p** + sin θ **t̂**,  **t̂** = **t** / |**t**|.  (5)

**Proposition 1.** For a unit vector **w** with **w** · **m** = 0, the vector (4) is tangent at **p** and has the length of **w**.

*Proof.* Write c = **m** · **p**. Using **w** · **m** = 0,

  **t** · **p** = **w** · **p** − (**w** · **p**)(c + 1)/(1 + c) = 0.

For the length,

  |**t**|² = |**w**|² − 2 (**w**·**p**)/(1+c) · (**w**·**m** + **w**·**p**) + (**w**·**p**)²/(1+c)² · |**m** + **p**|².

Substituting **w** · **m** = 0 and |**m** + **p**|² = 2 + 2c, the last two terms cancel. ∎

**Proposition 2.** The step (4)–(5) is second-order accurate, and it is exact for a steady flow of constant speed along a great circle.

*Sketch.* The transported vector **t** differs from **w**(**m**) by a vector of order h along **p**. In (5) this changes **p**′ only at order h³. The rest of the step agrees with (3) to order h², so the local error is O(h³). For flow along a great circle, transport along that circle carries the velocity at **m** exactly onto the velocity at **p**, and (5) is then the exact rotation. ∎

*Verification* (tests/run-tests.mjs):
- a uniform eastward flow along the equator is followed exactly;
- in solid-body rotation, 2000 steps keep the latitude to within 8 × 10⁻⁶ °.

The simpler rule, which moves along the great circle of the initial bearing, drifts by 0.9° in the same test.

## 3. Showing a Flow

- **Time step.** The step h is chosen so that a particle moving at 1 m s⁻¹ moves 0.16 pixels per frame, whatever the zoom. The legend states the resulting speed-up (of order 10⁵).
- **Life and rebirth.** Particles live for 40 to 110 steps, and are reborn uniformly in area within the visible part of the Earth.
- **Trails.** Each draws a trail of its last 15 steps, fading with age.

*Caveat.* Particles show streamlines of the field at each instant: the field is interpolated in time, but the particles are not air parcels followed through the forecast. Trajectories would integrate (2) at the real time step, from one step of the forecast to the next.
