---
title: A Primitive-Equation Model, the Next Laboratory Model
status: seed
date: 2026-10-22
related: [core/domains/meteorology, lib/physis/continuum, lib/calx, theory/meteorology/shallow-water-model.md]
---

# A Primitive-Equation Model, the Next Laboratory Model

## 1. Idea

The laboratory's models now retrace the history of numerical weather prediction:
- the barotropic vorticity equation (1950);
- the shallow-water equations (ADR 0019).

The next step is the one operational models took in the 1960s and 1970s: the hydrostatic primitive equations on several levels. With them comes baroclinic instability, the growth of storms from the vertical shear of the wind, which neither single-layer model can represent.

## 2. Sketch

- **The equations.** Vorticity, divergence, temperature and surface pressure, in σ coordinates. The spectral transform in the horizontal, as now; finite differences in the vertical (Simmons and Burridge 1981); and a semi-implicit step, since gravity waves of the external mode are fast.
- **Tests.**
  - The baroclinic wave of Jablonowski and Williamson (2006): a balanced zonal jet, perturbed, grows a wave whose evolution has published reference solutions.
  - The Held and Suarez (1994) forcing, to test the model's climate.
- **Data.** GFS on its isobaric levels, through THREDDS (the `ncss` kind already reads them). Ten to twenty levels would be needed, not three.
- **Where.**
  - A C++ dynamical core in `core/domains/meteorology`, on calx's spectral transforms, is the right home (ADR 0003).
  - The JavaScript stand-ins have shown the design, and their tests could be carried over.

## 3. Open Questions

- Moisture: dry first, as the test cases are.
- Topography, which the shallow-water model also lacks.
- Whether a learned correction (as in NOAA's AIGFS) could be studied against our own dynamics, to separate what the physics gives from what the data teach.
