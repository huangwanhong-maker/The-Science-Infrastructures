---
title: Relations Between Galaxies: Supergalactic Frame, Tidal Index, Hubble Flow
status: draft
implemented-by:
  - applications/astronomy/celestial-atlas/src/core/frames.js
  - applications/astronomy/celestial-atlas/src/core/galaxies.js
  - applications/astronomy/celestial-atlas/src/views/galaxy-diagrams.js
references: [deVaucouleurs1976, Karachentsev1999, Karachentsev2013, Tully2015]
---

# Relations Between Galaxies: Supergalactic Frame, Tidal Index, Hubble Flow

## 1. The Supergalactic Frame

Nearby galaxies crowd towards a plane, the plane of the Local Supercluster. The supergalactic coordinate system (de Vaucouleurs, de Vaucouleurs and Corwin 1976) is defined in galactic coordinates by two directions:

- its north pole, at (l, b) = (47.37°, +6.32°);
- the zero of supergalactic longitude, at (l, b) = (137.37°, 0°).

**Proposition 1.** The two directions are orthogonal.

*Proof.* With unit vectors **p** for the pole and **o** for the origin,

  **p** · **o** = cos 6.32° · cos 0° · cos(47.37° − 137.37°) + sin 6.32° · sin 0° = cos 6.32° · cos(−90°) = 0. ∎

The rows **o**, **p** × **o** and **p** therefore form an orthonormal, right-handed basis. The matrix with these rows maps galactic to supergalactic coordinates:

  R_SG = [ **o** ; **p** × **o** ; **p** ].  (1)

Composed with the ICRS-to-galactic matrix of the Hipparcos catalogue, it maps ICRS directly:

  R = R_SG · R_G.  (2)

*Verification.* The implementation reproduces the supergalactic angles that Tully (2015) lists for 43 036 galaxies. The median difference is 0.00012°, which is the catalogue's rounding.

## 2. The Tidal Index

Let galaxy *i* have neighbours *k* of masses M_k (in solar masses), at distances D_ik (in Mpc) in space. Its tidal index (Karachentsev and Makarov 1999) is

  Θ₁(i) = max_k log₁₀( M_k / D_ik³ ) + C,   C = −10.96.  (3)

The neighbour that attains the maximum is the *main disturber* of *i*.

The masses are estimated from K-band luminosities with M/L_K = 1. The constant C is chosen so that Θ₁ = 0 when the Keplerian period of *i* around its main disturber equals the age of the universe. So:

- Θ₁ > 0 marks a galaxy bound in a group;
- Θ₁ < 0 marks a galaxy in the field.

The index Θ₅ replaces the maximum by the sum of the five largest terms.

*Remark.* The tidal force per unit mass exerted by k scales as M_k / D_ik³. Equation (3) is the logarithm of the strongest such force, normalised.

### 2.1 Inferring a Missing Mass

The catalogue omits the K-band luminosity of the Milky Way. Equation (3) can be inverted for each galaxy *i* whose published main disturber is *k*:

  log₁₀ M_k = Θ₁(i) − C + 3 log₁₀ D_ik.  (4)

The median of (4) over the galaxies disturbed by *k* estimates its mass. For the Milky Way, over 38 satellites, it gives log₁₀ M = 10.55 (range 10.34 to 10.81).

*Verification.* With this mass, recomputing (3) for the 770 galaxies whose published main disturber is in the catalogue gives:

- the same main disturber for 99.2 % of them. The rest are ties at the published precision of 0.1.
- Θ₁ to a median of 0.027.

### 2.2 Tidal Groups

**Definition.** Join *i* to its main disturber whenever Θ₁(i) > 0. A *tidal group* is a connected component of the resulting graph. Its leader is its most luminous member.

The largest tidal groups of the Local Volume are those of NGC 3368 (Leo I), M 81, M 31, the Milky Way, NGC 5128 (Centaurus A), NGC 4258 and M 94.

## 3. The Hubble Flow and Peculiar Velocities

A galaxy at distance d moves away from us at about v = H₀ d, the Hubble flow, plus a peculiar velocity that gravity induces. With velocities in the Local Group frame:

  v_pec = v_LG − H₀ d.  (5)

The Hubble diagram plots v_LG against d. Within about 1 Mpc, the Local Group's members do not follow the flow. Beyond it, galaxies scatter about the line, and the members of groups scatter most. The atlas draws the line with H₀ = 73 km s⁻¹ Mpc⁻¹.

*Caveat.* In the group catalogue of Tully (2015), distances are derived from velocities, in units of h⁻¹ Mpc (H₀ = 100 h). They cannot be used to study peculiar velocities. Doing so needs distances measured independently, as in Cosmicflows-4.

## 4. Mass and Light of Groups

For a group with line-of-sight velocity dispersion σ and size R, the virial theorem gives M ∝ σ² R / G. Plotted against the group's K-band luminosity, virial masses show that the mass-to-light ratio M/L_K increases with mass, from galaxy pairs to rich clusters (Tully 2015). Lines of constant M/L are straight lines of unit slope in the log–log plane.

## 5. References

- de Vaucouleurs, G., de Vaucouleurs, A., & Corwin, H. G. (1976). *Second Reference Catalogue of Bright Galaxies*. University of Texas Press.
- Karachentsev, I. D., & Makarov, D. I. (1999). In *Galaxy Interactions at Low and High Redshift*, IAU Symp. 186, 109.
- Karachentsev, I. D., Makarov, D. I., & Kaisina, E. I. (2013). Updated Nearby Galaxy Catalog. *AJ* 145, 101.
- Tully, R. B. (2015). Galaxy Groups: A 2MASS Catalog. *AJ* 149, 171.
