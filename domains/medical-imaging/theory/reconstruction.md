# The Reconstruction of Vessels from Images

*Status: draft (2026-11-04). The theory behind the domain's `reconstruction` and `modalities` ([ADR 0070](../../../development/decisions/0070-structures-reconstructed-from-images.md)).*

## 1. The Half-Maximum Locates a Blurred Wall

**Setting.** An image of a vessel is the true map f, constant c₁ inside the lumen and c₀ outside, convolved with a point-spread function h, which is symmetric and integrates to 1: g = f ∗ h.

**Proposition 1.** For a planar wall, f = c₀ + (c₁ − c₀)·1{x > 0}, the level (c₀ + c₁)/2 of g lies exactly at x = 0, whatever the symmetric h.

*Proof.* g(0) = c₀ + (c₁ − c₀) ∫_{y<0} h(y) dy. A symmetric h puts half its mass on each side, so g(0) = c₀ + (c₁ − c₀)/2. ∎

**Remark (curved walls).** For a cylinder of radius R, blurred by a Gaussian of deviation σ, the half-maximum contour moves inward by a term of order σ²/R: the wall's curvature puts less of the lumen within the blur on the inside. With σ ≈ 0.35 mm (a PSF of 0.8 mm FWHM) and R = 1.5 mm, the shift is a fraction of a voxel. This agrees with the checks' radius error of 0.10 mm on average.

**Why the maximum must be local.** The proposition needs c₁, the lumen's own value. Where the lumen's brightness varies along the vessel, a single global c₁ misplaces the wall. Contrast fades in small vessels, and TOF's inflow signal falls with depth into the slab. The pipeline therefore takes c₁ as the local peak: the largest lightly smoothed value among the nearest core voxel and its neighbours. The checks found a global level fattened TOF's brightest arteries twofold.

## 2. The Scale of a Tube

**Proposition 2.** For a cylinder of radius R with a sharp wall, the γ-normalised second derivative across the axis, σ²|∂²(f ∗ G_σ)|, is largest at σ = R/√2 (eikon's `features`).

*Proof.* In the cross-section, take the disc's indicator and the 2D Gaussian G. The divergence theorem gives the Laplacian on the axis: Δ(1_disc ∗ G)(0) = ∮ ∂_n G = 2πR · G′(R) = −(R²/σ⁴) e^{−R²/2σ²}. By symmetry the Laplacian splits equally between the two cross-sectional directions, so σ² |L_xx| = u e^{−u} with u = R²/2σ². This is largest at u = 1. ∎

Vesselness therefore detects; it does not delineate. The response falls off before the wall, and its scale reads the radius only to within the scale step (2^{1/4} in the presets). The pipeline uses it to find cores, and the half-maximum to find walls.

## 3. Time-of-Flight Contrast

**Setting.** In a spoiled gradient echo with repetition time TR and flip angle α, the longitudinal magnetisation just before pulse n+1 obeys M_{n+1} = M_n cos α · E₁ + M₀(1 − E₁), with E₁ = e^{−TR/T₁} (the Bloch equations between pulses, the transverse part spoiled).

**Proposition 3.** M_n = M_ss + (M_z(0) − M_ss)(E₁ cos α)^n, with M_ss = M₀(1 − E₁)/(1 − E₁ cos α).

*Proof.* M_ss is the fixed point of the affine recursion. Its deviation from the fixed point is multiplied by E₁ cos α at each step. ∎

**Consequences.**

- Stationary tissue sits at M_ss.
- Blood that entered the slab n pulses ago, fresh (M_z(0) = M₀), is brighter by (M₀ − M_ss)(E₁ cos α)^n.
- Blood that crossed a saturation band (M_z(0) = 0) is darker than its steady state. So arteries flowing into the slab shine and veins entering through the band do not: this is the contrast `integration` uses to tell them apart.
- The Ernst angle, cos α_E = E₁, maximises the steady signal sin α · M_ss (checked).

## 4. Telling Arteries from Veins

Let the vessels found be V (fused from the acquisitions that show all vessels) and the arteries seen be A (from those that show arteries only).

A branch b of V's skeleton is classed as an artery when more than half its centreline lies within one voxel of A; otherwise it is a vein. Classing by branch rather than by voxel is a choice, not a theorem: a vessel keeps one class along its length, as anatomy has it.

Errors come from two sources:

- **Registration.** A misregistration by δ moves A against V. Any branch whose radius exceeds δ minus one voxel keeps most of its centreline within reach. In the checks δ = 0.25 mm, and every branch found was classed rightly: arteries at Dice 0.87, veins at 0.94.
- **Junctions.** Where an artery and a vein touch, they merge in V and share a junction. The branch rule then splits them at the junction.
