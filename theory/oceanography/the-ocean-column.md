---
title: The Ocean Beneath a Place: Currents and Temperature by Depth
status: draft
implemented-by:
  - applications/geography/earth-atlas/src/core/ocean.js
  - applications/geography/earth-atlas/src/ui/column.js
  - data/products/ocean_column.py
references: [Gill1982, deBoyerMontegut2004, Barton2021, Bleck2002, RDI2011]
---

# The Ocean Beneath a Place: Currents and Temperature by Depth

## 1. Two Ways of Seeing a Column

A map shows the ocean's surface. Beneath any point lies a column several kilometres deep, in which currents change speed and direction with depth and the water grows colder and denser. Oceanographers see such a column in two ways.

**By measuring it.** An acoustic Doppler current profiler (ADCP), on a ship's hull or a mooring, sends pulses of sound down through the water. Particles and plankton drifting with the water scatter some of the sound back. The Doppler shift of the echo gives the speed of the scatterers along the beam: for a pulse of frequency f₀ and sound speed c, the shift is Δf = 2 f₀ v / c. Three or four beams, inclined about 20° to 30° from the vertical, give the three components of velocity. Gating the echo in time gives them in depth cells, because an echo from depth z returns after 2z / c (RDI 2011). The strength of the echo is recorded too. It shows the scatterers themselves: the layers of zooplankton that rise towards the surface at dusk and sink at dawn, the largest daily migration of animals on Earth.

**By computing it.** An ocean model integrates the equations of motion of a rotating, stratified fluid and assimilates observations. The laboratory reads the US Navy's global model, the ocean component of ESPC-D-V02 (HYCOM on a 1/12° grid, with 40 output levels from the surface to 5000 m; Barton et al. 2021). The model gives currents, temperature and salinity at every point, now and a week ahead. It has no plankton, so it cannot show the echoes.

The Earth Atlas draws the model's column as a profiler would record it: depth downward, time across.

## 2. Inertial Oscillations

The currents in such a record sway with a period of about a day, and the period depends on latitude. The cause is the Earth's rotation.

*Hypotheses.* Horizontal motion of a layer of water, with no pressure gradient and no friction, on an f-plane: the Coriolis parameter f = 2Ω sin φ is taken constant near latitude φ.

**Proposition 1 (inertial oscillation; Gill 1982, §8.3).** The velocity (u, v) satisfies

  du/dt = f v,  dv/dt = −f u,  (1)

and so turns clockwise in the northern hemisphere (anticlockwise in the southern) at constant speed, with period

  T_i = 2π / |f| = π / (Ω |sin φ|).  (2)

*Proof.* Let w = u + i v. Then (1) reads dw/dt = −i f w, so w(t) = w(0) e^{−i f t}: rotation at angular rate f, clockwise for f > 0. ∎

At 30° latitude, 2Ω sin φ = Ω, so T_i is one sidereal day, 23.93 h. At 60° it is 13.8 h. Towards the equator f vanishes and T_i grows without bound; within a few degrees of the equator there are no inertial circles. A gust of wind sets the surface layer moving, and the layer then circles inertially for days, so the currents near the surface swing round at T_i. The figure's caption gives T_i for the place, and the tests check (2) at 30°.

## 3. The Mixed Layer

Near the surface, wind and cooling stir the water until its temperature hardly changes with depth. Below this mixed layer, temperature falls through the thermocline.

**Definition 1 (de Boyer Montégut et al. 2004).** The mixed-layer depth is the shallowest depth below 10 m at which the temperature differs from its value at 10 m by 0.2 °C.

The reference depth of 10 m avoids the daily warming of the top few metres. The threshold of 0.2 °C is small enough to follow the base of the stirred layer, and large enough to ignore noise. The atlas interpolates linearly between the model's levels. For a profile that is uniform to 50 m and then cools by 0.01 °C per metre, the depth is exactly 70 m; the tests check this.

The mixed layer is shallow in summer and in the tropics (tens of metres) and deep in winter at high latitude (hundreds of metres). It sets how much heat the ocean exchanges with the atmosphere over a season.

## 4. What the Record Shows

| Panel | Quantity | Colours |
|---|---|---|
| Eastward current | u (m s⁻¹) | red east, blue west, white still; symmetric about zero |
| Northward current | v (m s⁻¹) | red north, blue south |
| Temperature | T (°C) | thermal, with isotherms and the mixed layer dashed |

Above the panels, a strip shows the Sun's height at the place: pale for day and plum for night, with twilight between. Time runs in UTC. A dashed line marks the model's latest run, after which values are forecast. Beneath the deepest level with water, the record is grey, and the depth of the sea floor from the relief (ETOPO 2022) is dashed. Currents are named by the direction they flow *toward*, the oceanographer's convention, which is opposite to the meteorologist's for wind.

## 5. Data on Demand

The model's grid has 19 million points, about 12 million of them at sea, so the laboratory cannot fetch every column ahead of time. When a place is chosen, the local server asks the data library for that column. The library reads the model's coordinates, finds the grid point nearest the place, and requests the column from the model's OPeNDAP server, two days back and four ahead every 6 hours, in short runs of times, four at a time. That takes about four minutes the first time; afterwards the column is in the cache (ADR 0023).

## 6. Verification

| Test | Agreement |
|---|---|
| Inertial period at 30°, against one sidereal day (86 164.1 s) | 6 × 10⁻⁴ s |
| Inertial period at 60° S against (2) | Exact |
| Direction toward which a current flows: E, N, W, S, NE | Exact |
| Mixed-layer depth of the constructed profile (70 m) | 1 × 10⁻¹³ m |
| A column mixed to its deepest value | Reported as reaching the bottom |
| Interpolation in time and depth of f = a + bt + cz + dtz, on unequal steps | 9 × 10⁻¹⁶ |
| Values beneath the sea floor | None |

*Caveat.*
- The model's grid point may lie a few kilometres from the place chosen; the figure says how far.
- The model's analyses lag about two days behind the present (the latest run was two days old on 22 October 2026), so "now" is already a forecast.
- The model's collection on the server is assembled from many runs, and its record of which run each time came from is not always consistent. The product records it as given.
- A model cannot show the scattering layers of plankton, which only an echo sounder sees.
