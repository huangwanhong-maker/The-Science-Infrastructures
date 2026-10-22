---
title: How the Energy of Ocean Waves Travels
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/derived.js
  - applications/geography/earth-atlas/src/layers/flow.js
references: [Holthuijsen2007, Phillips1977]
---

# How the Energy of Ocean Waves Travels

## 1. Dispersion in Deep Water

*Hypotheses.* The waves are small-amplitude surface gravity waves on water of depth d. Surface tension is neglected, and the water is inviscid and irrotational.

Linear theory gives the dispersion relation between the angular frequency ω and the wavenumber k:

  ω² = g k tanh(k d).  (1)

**Definition 1.** Water is *deep* for a wave when d > L/2, where L = 2π/k is the wavelength. Then tanh(k d) > 0.996, and

  ω² = g k.  (2)

## 2. Phase and Group Speed

**Proposition 1.** In deep water, for a wave of period T = 2π/ω:

  L = g T² / 2π,  c = ω/k = g T / 2π,  c_g = dω/dk = c / 2 = g T / 4π.  (3)

*Proof.* From (2), k = ω²/g, so c = ω/k = g/ω = g T / 2π and L = c T. Differentiating ω = (g k)^{1/2} gives dω/dk = ½ (g/k)^{1/2} = c/2. ∎

*Remark.* The energy of a group of waves, and with it the wave field forecast by a model, travels at the group speed c_g, not at the speed of the crests. A swell of period 14 s travels at 10.9 m s⁻¹, about 940 km a day. Swell born in a Southern Ocean storm therefore reaches Hawaii or California a week later.

## 3. In the Atlas

- **Velocity.** WaveWatch III gives the peak period T_p and the peak direction θ_p. The direction is the one the waves come *from*, in degrees clockwise from north. The atlas draws particles moving with the velocity

  (u, v) = c_g(T_p) · (−sin θ_p, −cos θ_p),  (4)

  that is, at the group speed, toward the direction the waves travel.
- **Directions.** Peak directions are interpolated through their unit vectors (mean of 350° and 10° is 0°, not 180°).

*Caveat.* Equation (3) holds only in deep water. Where the depth is less than half a wavelength (on shelves, for long swell), the waves slow down, and (4) overstates their speed. The peak describes only part of the spectrum: a sea with swell and wind waves from different directions has energy travelling in more than one direction.
