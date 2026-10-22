---
title: The Tidal Spectrum as a Lattice
status: seed
date: 2026-10-21
related: [core/domains/tides, lib/morphe/algebra, core/inference/estimation]
---

# The Tidal Spectrum as a Lattice

## 1. Idea

- **Doodson's expansion.** It writes the tide-generating potential as a sum of constituents.
  - Each constituent's argument is k·θ(t).
  - Here θ = (τ, s, h, p, N′, pₛ) are six fundamental astronomical angles.
  - k ∈ ℤ⁶ is the constituent's Doodson number, which is conventionally written with an offset.
- **The frequencies form a lattice image.** They are the image of the rank-6 lattice ℤ⁶ under the linear map k ↦ k·ω.
- **The tide lives on a torus.** As a function of time, the tide factors through the 6-torus T⁶. Each constituent is a character of T⁶, restricted to a line that winds through the torus.

## 2. Consequences to Explore

- **Resolution in harmonic analysis.** A record of length T can separate two constituents only if their frequencies differ by more than about 1/T (the Rayleigh criterion). So which constituents can be resolved is a question about lattice differences k with small k·ω: a problem in Diophantine approximation.
- **Nodal modulation.** The lunar nodal cycle of 18.6 years, the N′ coordinate, is usually handled by correction factors. In the torus picture it is simply one more coordinate.
- **Shallow-water constituents.** In shallow water, nonlinear interactions produce sums and differences of constituents. Examples are M4 = 2·M2 and MS4 = M2 + S2. These are additions in the lattice.

## 3. Question

Does the torus picture give a cleaner basis for tidal analysis than the classical lists of constituents, and one whose completeness can be proved?
