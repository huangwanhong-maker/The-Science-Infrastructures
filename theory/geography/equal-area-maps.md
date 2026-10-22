---
title: Equal-Area Maps and Data per Place
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/projections.js
  - applications/geography/earth-atlas/src/core/raster.js
  - applications/geography/earth-atlas/src/core/join.js
references: [Savric2018, Snyder1987]
---

# Equal-Area Maps and Data per Place

## 1. Why Area Matters on a Thematic Map

A reader judges a coloured map by area: the more of the map a colour covers, the more weight it seems to carry. A map of a quantity over the Earth should therefore keep areas in proportion.

**Definition 1.** A map projection P: (λ, φ) ↦ (x, y) is *equal-area* when its area scale

  s(λ, φ) = |∂(x, y)/∂(λ, φ)| / cos φ  (1)

is the same at every point. The Jacobian determinant in (1) is divided by the area element cos φ dλ dφ of the unit sphere.

*Remark.* Mercator's projection is conformal: it keeps angles, so its area scale is sec² φ. It doubles areas at 45° and quadruples them at 60°. It suits navigation, and misleads on maps of quantities.

## 2. Equal Earth

Equal Earth (Šavrič, Patterson and Jenny 2018) is a pseudocylindrical projection: its parallels are straight lines. It is defined with a parametric latitude θ by

  sin θ = (√3 / 2) sin φ,
  x = 2√3 λ cos θ / [3 (9 A₄ θ⁸ + 7 A₃ θ⁶ + 3 A₂ θ² + A₁)],
  y = A₄ θ⁹ + A₃ θ⁷ + A₂ θ³ + A₁ θ,  (2)

with A₁ = 1.340264, A₂ = −0.081106, A₃ = 0.000893 and A₄ = 0.003796.

**Proposition 1.** The map (2) is equal-area, with s = 1.

*Proof.* Write y′ for dy/dθ = 9A₄θ⁸ + 7A₃θ⁶ + 3A₂θ² + A₁, so that x = 2√3 λ cos θ / (3 y′).
- Because x is linear in λ and y does not depend on λ, the Jacobian is (∂x/∂λ)(∂y/∂φ).
- ∂x/∂λ = 2√3 cos θ / (3 y′).
- ∂y/∂φ = y′ dθ/dφ, and dθ/dφ = (√3/2) cos φ / cos θ by differentiating sin θ = (√3/2) sin φ.
- The product is (2√3/3)(√3/2) cos φ = cos φ, so s = 1. ∎

*Verification.* The atlas evaluates (1) numerically over the whole sphere and finds it constant to 2 × 10⁻¹⁰. The ratio of the map's axes is 2.0546, which the paper gives as "about 2.05".

**The inverse** solves y(θ) = y₀ by Newton's method, starting at θ = y₀ (as in the paper, §3). It converges to 10⁻¹² in a few steps.

## 3. Joining Values to Places

An indicator gives a value per country. The base map gives a polygon per country. A choropleth needs a *join*: a map from records to polygons. The join is not trivial, because the two sources use different identifier systems:

| Country | Natural Earth's codes | World Bank's code |
|---|---|---|
| Kosovo | KOS (own), XK (ISO alpha-2), KSV (World Bank column) | XKX |
| France, Norway | −99 in the ISO alpha-3 column (Natural Earth's ISO_A3_EH column gives FRA, NOR) | FRA, NOR |

The atlas tries the identifiers of each country in a stated order: World Bank, corrected ISO, ISO, Natural Earth, then the two-letter codes. It lists the records that matched no polygon (Gibraltar, at 1:50 million). A gap in the map is then explained, not silent.

The countries are rasterised by the even-odd rule at pixel centres (Foley et al. 1990, §3.6). Neighbouring polygons that share an edge never both claim a pixel. Area tests compare the raster with the polygons' exact area: for the land of the Earth they agree to within 8 × 10⁻⁶ of the sphere.
