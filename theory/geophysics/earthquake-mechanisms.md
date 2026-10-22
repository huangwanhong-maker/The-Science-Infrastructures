---
title: Earthquake Mechanisms: Moment Tensors, Beach Balls and Slabs
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/moment.js
  - applications/geography/earth-atlas/src/render/beachballs.js
  - applications/geography/earth-atlas/src/ui/stereonet.js
  - applications/geography/earth-atlas/src/core/section.js
  - data/products/mechanism_catalog.py
references: [AkiRichards2002, JostHerrmann1989, Dziewonski1981, Ekstrom2012, HanksKanamori1979, Frohlich1992]
---

# Earthquake Mechanisms: Moment Tensors, Beach Balls and Slabs

## 1. The Moment Tensor

Seen from far away, an earthquake is a point source of forces acting in balanced pairs (couples). Its strength and geometry are a symmetric tensor M, the *moment tensor*: M_ij is the moment of the couple of forces in direction i separated in direction j. The Global CMT catalogue gives it in the basis (r up, t south, p east), and the atlas works in (north, east, down): N = −t, E = p, D = −r.

For slip on a fault of area A, with average slip D and rigidity μ, the source is a *double couple*,

  M = M₀ (n dᵀ + d nᵀ),  M₀ = μ A D,  (1)

with n the unit normal of the fault and d the unit direction of slip (n · d = 0). Its size is given by the *moment magnitude*

  Mw = (2/3)(log₁₀ M₀ − 16.1),  M₀ in dyne·cm,  (2)

(Hanks and Kanamori 1979). Each unit of Mw is a factor of 10^1.5 ≈ 32 in moment. The Tōhoku earthquake of 2011 has M₀ = 5.3 × 10²⁹ dyne·cm (5.3 × 10²² N·m) and Mw 9.1.

**Proposition 1.** The double couple (1) has eigenvalues M₀, 0 and −M₀, with eigenvectors

  T = (n + d)/√2,  B = n × d,  P = (n − d)/√2.  (3)

*Proof.* With n · d = 0 and |n| = |d| = 1: M(n + d) = M₀(n (d · (n + d)) + d (n · (n + d))) = M₀(n + d). Likewise M(n − d) = M₀(−n + d) = −M₀(n − d). And M(n × d) = 0, since n × d is orthogonal to both. ∎

**Corollary (the fault-plane ambiguity).** Exchanging n and d leaves M unchanged. So the fault (normal n, slip d) and the *auxiliary plane* (normal d, slip n) radiate identical waves. Seismograms alone cannot tell which plane broke; the aftershocks, the surface rupture or the tectonic setting decide.

*Remark (the symmetry group).* M is also unchanged when n and d both change sign. The two operations, exchange and joint reversal, generate a Klein four-group V₄ ≅ ℤ/2 × ℤ/2 acting on pairs (n, d). A mechanism is an orbit of this action. Its four elements are the two planes, each with two directions of the normal. Strike, dip and rake fix the normal's direction by convention (into the hanging wall), and so name the two planes.

Conversely, (3) gives the planes from a measured tensor: n = (T + P)/√2 and d = (T − P)/√2, or these exchanged. A real tensor is not exactly a double couple: its middle eigenvalue is not zero. The planes are those of the nearest double couple, from the eigenvectors of largest and smallest eigenvalue. The atlas computes them by Jacobi rotations. For the 70,451 events of the catalogue, 99 % agree with the catalogue's planes within 0.76°, which is the catalogue's own rounding to whole degrees.

## 2. Strike, Dip and Rake

A plane is named by its *strike* φ (the azimuth of a horizontal line in it, with the plane dipping to the right), its *dip* δ, and the *rake* λ of the slip within it. Rake is measured from the strike direction: positive when the hanging wall rises (reverse faulting), negative when it falls (normal), and near 0° or ±180° for strike-slip. In (north, east, down) (Aki and Richards 2002, eqs. 4.88 and 4.89):

  n = (−sin δ sin φ, sin δ cos φ, −cos δ),
  d = (cos λ cos φ + sin λ cos δ sin φ, cos λ sin φ − sin λ cos δ cos φ, −sin λ sin δ).  (4)

The atlas inverts (4) once n points upward, into the hanging wall. The tests build double couples from four triples (φ, δ, λ) by (1) and (4), and recover each to about 10⁻⁶ degrees.

## 3. The Beach Ball

The first motion of a P wave leaving the source in the direction x is a push (compression) where xᵀMx > 0, and a pull where it is negative. The two nodal planes divide the directions into four quadrants: compressional around T, dilatational around P.

A *beach ball* draws the lower hemisphere of directions, seen from above, by the Lambert azimuthal equal-area projection (the Schmidt net). Compressional quadrants are coloured.

**Proposition 2.** The map from a direction at angle θ from the downward vertical and azimuth α to the point (r sin α, r cos α), with r = √2 sin(θ/2), multiplies areas by a constant.

*Proof.* dr = (√2/2) cos(θ/2) dθ, so r dr dα = sin(θ/2) cos(θ/2) dθ dα = ½ sin θ dθ dα: one half of the element of area on the unit sphere. ∎

So equal areas of the ball stand for equal solid angles of rays. In the stereographic projection (the Wulff net), r = tan(θ/2); it preserves angles instead, and is the one structural geologists use for measuring them.

On the map, each ball is a disc whose shader evaluates xᵀMx for the direction each pixel stands for. The disc is turned so that its north is the map's north at the earthquake. The stereonet in the place card adds the nodal planes as great circles, and the P and T axes.

## 4. Kinds of Faulting

**Proposition 3.** For the orthonormal axes T, B, P with plunges δ_T, δ_B, δ_P,

  sin² δ_T + sin² δ_B + sin² δ_P = 1.  (5)

*Proof.* sin δ is an axis's vertical component. The three components of the vertical unit vector in an orthonormal basis have squares summing to 1. ∎

So every mechanism is a point of a triangle (a ternary diagram), and Frohlich (1992) divides the triangle into four classes:

| Class | Condition | Where it happens |
|---|---|---|
| Thrust | T plunges more than 50° | where plates converge: trenches, mountain fronts |
| Normal | P plunges more than 60° | where the crust is stretched: rifts, outer rises, ridges |
| Strike-slip | B plunges more than 60° | transform faults |
| Oblique | none of these | combined motions |

The beach balls on the map are coloured by class.

## 5. Sections Through Subduction Zones

Where an oceanic plate sinks beneath another, earthquakes follow it into the mantle. In a section across the trench, they lie along an inclined band, the *Wadati–Benioff zone*: thrusts at the shallow plate interface, then events within the bending and sinking slab.

The deepest earthquakes, near 660 to 700 km, stop where the mantle's olivine has transformed through wadsleyite (near 410 km) and ringwoodite into bridgmanite and ferropericlase (near 660 km). The section draws both discontinuities.

The atlas's sections run across the nearest plate boundary, toward the overriding plate, so that a slab always dips to the right. They show the Global CMT centroids within a band about the line. The along-track and cross-track distances on the sphere are

  across = R asin(sin δ₁₃ sin(θ₁₃ − θ₁₂)),  along = R atan2(sin δ₁₃ cos(θ₁₃ − θ₁₂), cos δ₁₃),  (6)

for a line leaving point 1 at azimuth θ₁₂ and a place 3 at angular distance δ₁₃ and azimuth θ₁₃. The second form keeps its sign for places behind the line's start.

## 6. Verification

| Test | Agreement |
|---|---|
| Four double couples (thrust, normal, strike-slip, oblique) built by (1) and (4), recovered | 10⁻⁶ ° |
| Their classes by the plunges of the axes | Exact |
| Eigenvalues of a unit double couple: 1, 0, −1 | 9 × 10⁻¹⁶ |
| A push along T and a pull along P | Exact |
| The equal-area projection, there and back | 4 × 10⁻¹⁶ |
| Nodal planes of all 70,451 Global CMT events against the catalogue's | 99 % within 0.76° |
| Tōhoku 2011: a thrust on a plane striking 203° and dipping 10° | Within 1° |
| A place 80 km beside a line and 400 km along it, by (6) | 10⁻¹¹ km |

*Caveat.*
- Centroids are the centres of the moment release. For great earthquakes they can lie tens of kilometres from where rupture began. Their depths are less precise than those of relocated hypocentres.
- The catalogue is complete only above about magnitude 5.5, so smaller events are missing from the sections.
