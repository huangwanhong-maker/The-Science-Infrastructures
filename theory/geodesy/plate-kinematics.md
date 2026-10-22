---
title: The Motion of Rigid Plates on a Sphere
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/plates.js
  - applications/geography/earth-atlas/src/layers/plates.js
  - applications/geography/earth-atlas/src/layers/stations.js
references: [ArgusGordonDeMets2011, Altamimi2023, Bird2003, ArgusGordon1991, CoxHart1986, Blewitt2016]
---

# The Motion of Rigid Plates on a Sphere

## 1. Plates as Rotations

The Earth's outer shell, the lithosphere, is broken into plates that move as nearly rigid caps over the sphere. A rigid motion of the sphere that keeps its centre fixed is an element of SO(3).

**Proposition 1 (Euler 1776).** Every element of SO(3) other than the identity is a rotation about a unique axis through the centre.

*Proof.* A matrix R ∈ SO(3) has eigenvalues of modulus 1 whose product is det R = 1. The non-real ones come in conjugate pairs. So 1 is an eigenvalue, and its eigenvector is the axis. In the plane orthogonal to the axis, R acts as an element of SO(2), a rotation. The axis is unique unless R = I. ∎

A plate moving through time is a path R(t) in SO(3). Its velocity, Ṙ R⁻¹, is a skew-symmetric matrix, an element of the Lie algebra so(3). The isomorphism ω ↦ [ω]×, which sends a vector to the matrix of its cross product, identifies so(3) with ℝ³. So the instantaneous motion of a plate is one vector, its angular velocity ω, and a point r of the plate moves at

  v = ω × r.  (1)

The direction of ω is the plate's *Euler pole*, and |ω| is its rate of rotation.

**Proposition 2 (relative motion).** The velocity of plate B relative to plate A at any point r is (ω_B − ω_A) × r.

*Proof.* v_B − v_A = ω_B × r − ω_A × r, by the bilinearity of the cross product. ∎

So angular velocities subtract as vectors, although finite rotations do not commute. Over a finite time the composition of rotations involves the bracket [ω_A, ω_B] = ω_A × ω_B (the Baker–Campbell–Hausdorff series): this matters for reconstructing plates over millions of years, and not for their velocities today.

## 2. The Frame: No Net Rotation

Equation (1) needs a frame in which the ω are measured. Plate motions fix only the *relative* angular velocities ω_B − ω_A. The no-net-rotation frame chooses the absolute ones so that the surface as a whole does not turn.

**Definition 1.** The *net rotation* of a model is the angular velocity whose rigid rotation best fits the velocity field of the whole surface S:

  ω_net = (3 / 8π R⁴) ∫_S r × v dA.  (2)

A model has no net rotation when ω_net = 0 (Argus and Gordon 1991).

**Proposition 3.** For plates P with angular velocities ω_P,

  ω_net = (3 / 8π) Σ_P Q_P ω_P,  where Q_P = ∫_P (I − r̂ r̂ᵀ) dA,  (3)

with the integral over the unit sphere. If all the plates turned together (ω_P = ω), then ω_net = ω.

*Proof.* With v = ω_P × r on plate P and |r| = R, r × (ω × r) = ω R² − r (r · ω). Dividing by R⁴, (2) becomes (3). Over the whole unit sphere ∫ r̂ r̂ᵀ dA = (4π/3) I by symmetry, so Σ_P Q_P = 4π I − (4π/3) I = (8π/3) I. ∎

The tensors Q_P depend on the plates' shapes and act like moments of inertia. The tests evaluate (3) on a grid of a quarter degree, with NNR-MORVEL56's angular velocities and PB2002's outlines. The result is 0.0008°/Myr, about a thousandth of the Pacific's 0.651°/Myr. What remains comes from the four MORVEL plates that PB2002 does not outline.

## 3. Two Models of the Plates' Motion

**NNR-MORVEL56** (Argus, Gordon and DeMets 2011) gives 56 angular velocities in the no-net-rotation frame. They are estimated from spreading rates at mid-ocean ridges (over 3 million years), the azimuths of transform faults and earthquake slip vectors, and, for Bird's 31 small plates, his estimates. It describes the plates over geological time, and the atlas uses it to show the plates and their boundaries.

**ITRF2020's plate motion model** (Altamimi et al. 2023) gives angular velocities for 13 plates. They are estimated from GNSS velocities at stations in the plates' stable interiors, in the frame of ITRF2020 itself. The frame differs from the NNR frame by a small rotation, and its origin drifts by a small *origin rate bias* (0.37, 0.35, 0.74) mm/yr, which is added to ω × r. The model describes the present motion in the frame of the GPS velocities. The atlas therefore subtracts it from a station's velocity, wherever it covers the station's plate.

The choice matters. Relative to NNR-MORVEL56, GPS stations deep in four stable plates keep median velocities of 1.7 to 3.2 mm/yr, mostly the difference of frames. Relative to ITRF2020's model they keep 0.3 to 0.6 mm/yr. Those are the stations' errors and the plates' slight non-rigidity.

## 4. Across a Boundary

At a step of a plate boundary heading at azimuth α, let t = (sin α, cos α) be the step's direction and n = (cos α, −sin α) its normal to the right, both in (east, north). Following Bird (2003), let v be the velocity of the right-hand plate relative to the left. Then

  opening = v · n,  right-lateral = −v · t.  (4)

The opening is positive when the plates separate and negative when they converge. The right-lateral component is positive when, seen from either side, the other plate moves to the right; this sense is the same from both plates.

Bird classes each step by its setting and its motion: spreading ridge (OSR), oceanic transform (OTF) and oceanic convergent boundary (OCB) at sea; rift (CRB), continental transform (CTF) and continental convergent boundary (CCB) on land; and subduction zone (SUB), where a Benioff zone or a volcanic arc shows that one plate sinks. His identifier XX/YY means the right-hand plate YY sinks beneath XX; XX\\YY, the opposite.

## 5. GPS Velocities

A GNSS station's position, measured daily, drifts with its plate, jumps in earthquakes and when its antenna is changed, and swings with the seasons. MIDAS (Blewitt et al. 2016) estimates its velocity as the median of the slopes between pairs of positions about one year apart. Pairs a year apart cancel the seasonal swing. The median ignores the few pairs that straddle a jump. Its uncertainty is estimated from the spread of the slopes about the median.

## 6. Verification

| Test | Agreement |
|---|---|
| NNR-MORVEL56 velocities at four places against the EarthScope Plate Motion Calculator | 0.005 mm/yr |
| ITRF2020's model, the Pacific plate, against the same | 0.003 mm/yr |
| PB2002's 52 plates tile the sphere | No cell unassigned; areas sum to 4π within 2 × 10⁻⁹ |
| Net rotation of NNR-MORVEL56 on PB2002, by (3) | 0.0008°/Myr |
| Ridge and rift steps that open; subduction steps that close | 2,342 of 2,349; 1,122 of 1,127 |
| Transform steps that slip in the sense Bird gives (where his rate exceeds 2 mm/yr) | 1,563 of 1,572 |
| GPS stations in four stable interiors, relative to ITRF2020's model (median) | 0.3 to 0.6 mm/yr |
| A station's 20-year daily record against its MIDAS velocity (P123) | 0.3 mm/yr |

*Caveat.*
- Plates are not perfectly rigid. Wide zones (the Himalaya and Tibet, the western United States, the Aegean) deform across hundreds of kilometres, and their GPS velocities are not a plate's.
- PB2002 does not outline MORVEL's Capricorn, Lwandle, Macquarie and Sur plates; there the containing plate's motion is used.
