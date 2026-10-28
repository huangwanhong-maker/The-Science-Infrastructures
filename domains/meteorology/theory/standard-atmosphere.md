---
title: The U.S. Standard Atmosphere 1976 to 86 km
status: draft
implemented-by:
  - domains/meteorology/atmosphere/include/meteorology/atmosphere/standard.hpp
checked-by:
  - domains/meteorology/atmosphere/tests/check_standard.cpp
references: [USStandardAtmosphere1976, ICAO1993]
---

# The U.S. Standard Atmosphere 1976 to 86 km

Equation, table and page numbers are those of the document: NOAA, NASA and USAF (1976), *U.S. Standard Atmosphere, 1976*, NOAA-S/T 76-1562 (NASA-TM-X-74335, NTRS 19770009539), Part 1.

## 1. The Purpose of a Standard Atmosphere

A standard atmosphere is a convention, not a measurement. It is one steady, horizontally uniform column of dry air, idealised from the mean state of middle latitudes, against which instruments are calibrated and measurements compared: an altimeter reads pressure and reports the height at which the standard has that pressure, an aircraft's performance is quoted in it, and an instrument mounted at a height can be given the air it would meet there on an average day. Its value lies in being fixed and shared, so its equations are definitions: a computation is right when it reproduces them, whatever the weather does.

The 1976 standard shares its temperature profile with the ICAO Standard Atmosphere (ICAO 1993) to 80 km′, which is ICAO's top; ICAO's constants differ from it in the sixth figure, and its pressures by up to 10⁻⁵ (§9). Below 86 km the 1976 is a mixed ideal gas in hydrostatic equilibrium; above 86 km it is a different model (§8), not implemented here.

## 2. Geopotential Height

*Hypotheses.* Gravity falls with geometric height $Z$ by the inverse square,

$$ g(Z) = g_0 \left(\frac{r_0}{r_0 + Z}\right)^2, \qquad g_0 = 9.806\,65\ \mathrm{m\,s^{-2}},\quad r_0 = 6\,356\,766\ \mathrm{m} \tag{eq. 17} $$

where $r_0$ is the effective radius of the Earth at the latitude where sea-level gravity is $g_0$ (p. 4; it absorbs the centrifugal acceleration there). It is the 1976's own radius, not that of any reference ellipsoid.

**Definition 1 (eq. 15).** The geopotential height of a point is $H = \frac{1}{g_0'} \int_0^Z g\, \mathrm{d}Z$, in standard geopotential metres (m′), where $g_0' = 9.806\,65$ m² s⁻² m′⁻¹: one m′ is the rise through which a kilogram gains 9.806 65 J.

**Proposition 1 (eqs. 18, 19).** With $\Gamma = g_0/g_0' = 1$ m′ m⁻¹,

$$ H = \Gamma\, \frac{r_0 Z}{r_0 + Z}, \qquad Z = \frac{r_0 H}{\Gamma r_0 - H}. $$

*Proof.* $\int_0^Z r_0^2 (r_0 + Z')^{-2}\, \mathrm{d}Z' = r_0^2 \left(\tfrac{1}{r_0} - \tfrac{1}{r_0 + Z}\right) = r_0 Z/(r_0 + Z)$; the second form solves the first for $Z$. ∎

*Rationale.* In $H$ the hydrostatic equation has constant gravity, $g\,\mathrm{d}Z = g_0'\,\mathrm{d}H$ (eq. 16), and integrates in closed form. At 86 km the two heights differ by 1148 m (84 852 m′ against 86 000 m).

## 3. Molecular-Scale Temperature and the Seven Layers

**Definition 2 (eq. 22).** The molecular-scale temperature is $T_M = T\, M_0/M$, with $T$ the kinetic temperature, $M$ the mean molecular weight and $M_0 = 28.9644$ kg kmol⁻¹ its sea-level value (eq. 21, from the composition of Table 3).

*Rationale.* The ideal-gas law $\rho = P M/(R^* T)$ involves $T$ and $M$ only through $T/M = T_M/M_0$, so one variable carries both. Below 80 km the air is fully mixed and $M = M_0$, so $T_M = T$.

**Definition 3 (eq. 23, Table 4).** $T_M$ is linear in $H$ in seven layers $b = 0, \dots, 6$,

$$ T_M(H) = T_{M,b} + L_{M,b}\, (H - H_b), \qquad H_b \le H < H_{b+1}, $$

from $T_{M,0} = T_0 = 288.15$ K, with each $T_{M,b}$ the value of the layer below at its base:

| $b$ | $H_b$ (km′) | $L_{M,b}$ (K/km′) | $T_{M,b}$ (K) | $P_b$ (Pa) | Name here |
|---|---|---|---|---|---|
| 0 | 0 | −6.5 | 288.15 | 101 325 | troposphere |
| 1 | 11 | 0.0 | 216.65 | 22 632.06 | tropopause |
| 2 | 20 | +1.0 | 216.65 | 5 474.889 | stratosphere |
| 3 | 32 | +2.8 | 228.65 | 868.0187 | stratosphere |
| 4 | 47 | 0.0 | 270.65 | 110.9063 | stratopause |
| 5 | 51 | −2.8 | 270.65 | 66.938 87 | mesosphere |
| 6 | 71 | −2.0 | 214.65 | 3.956 420 | mesosphere |
| (7) | 84.852 | | 186.946 | 0.373 384 | |

The base pressures are computed by Proposition 2 from $P_0$, not tabulated by the document as defining constants (p. 12). The document numbers its layers and does not name them; the names are those of the thermal structure each idealises (falling, constant or rising $T_M$). The model extends layer 0 down to −5 km′ (p. 12) and no further.

## 4. Pressure

*Hypotheses.* The air is an ideal gas, $P = \rho R^* T/M$ with $R^* = 8.314\,32$ J mol⁻¹ K⁻¹ (Table 2), in hydrostatic equilibrium, $\mathrm{d}P/P = -g M\, \mathrm{d}Z/(R^* T)$ (eq. 5).

**Proposition 2 (eqs. 33a, 33b).** In layer $b$,

$$ P = P_b \left[\frac{T_{M,b}}{T_{M,b} + L_{M,b}(H - H_b)}\right]^{g_0' M_0/(R^* L_{M,b})} \quad (L_{M,b} \neq 0), \qquad P = P_b \exp\left[\frac{-g_0' M_0 (H - H_b)}{R^* T_{M,b}}\right] \quad (L_{M,b} = 0). $$

*Proof.* By eq. 16 and Definition 2, $\mathrm{d}\ln P = -g_0' M_0\, \mathrm{d}H/(R^* T_M(H))$. With $T_M$ linear, $\mathrm{d}H = \mathrm{d}T_M/L_{M,b}$ and $\ln(P/P_b) = -\frac{g_0' M_0}{R^* L_{M,b}} \ln(T_M/T_{M,b})$; with $T_M$ constant, the exponential. ∎

**Proposition 3.** $T_M$, $P$ and $\mathrm{d}P/\mathrm{d}H$ are continuous at every base $H_b$; the second derivative of $P$ jumps where $L_{M,b}$ does.

*Proof.* $T_{M,b}$ and $P_b$ are defined as the lower layer's values at $H_b$, and $\mathrm{d}P/\mathrm{d}H = -g_0' M_0 P/(R^* T_M)$ is continuous with them. ∎

## 5. From 80 to 86 km: the Mean Molecular Weight

Above about 80 km dissociation of oxygen and the onset of diffusive separation lower $M$. The 1976 lets $M/M_0$ fall smoothly from 1 at 80 km to 0.999 579 at 86 km, where $M = 28.9522$ kg kmol⁻¹ joins the diffusive model above, and tabulates the ratio at every 0.5 geometric kilometre (Table 8, p. 9). Then $T = T_M\, M/M_0$, and at $Z_7 = 86$ km, $T_7 = 186.8673$ K (eq. 25), 0.0787 K below $T_M$.

The module interpolates Table 8 linearly in $Z$ between its nodes (the document names no interpolation), and computes every property from the kinetic $T$ where the document's equation uses $T$. Pressure, density, the speed of sound, the mean particle speed and the pressure scale height depend only on $T_M$ and $M_0$ (eqs. 33, 42, 44, 46, 50) and are unchanged.

*Caveat.* The document's own tables below 86 km were computed without Table 8 ("these values of $M/M_0$ were not included in the program used for computing the tables", p. 9), so they show a step of up to 0.04 per cent between 85.5 and 86 km. The document directs the correction: $T$, $M$ and the mean free path are multiplied by $M/M_0$, the number density and collision frequency divided by it, and the viscosities and conductivity recomputed. The module does this; its values between 80 and 86 km therefore differ from the printed tables there by those factors.

## 6. The Properties of the Air

With $N_A = 6.022\,169 \times 10^{23}$ mol⁻¹, $\sigma = 3.65 \times 10^{-10}$ m, $\gamma = 1.40$, $\beta = 1.458 \times 10^{-6}$ kg s⁻¹ m⁻¹ K⁻¹ᐟ² and $S = 110.4$ K (Table 2):

| Quantity | Formula | Eq. |
|---|---|---|
| gravity | $g = g_0 (r_0/(r_0 + Z))^2$ | 17 |
| number density | $N = N_A P/(R^* T)$ | 41 |
| mass density | $\rho = P M_0/(R^* T_M)$ | 42 |
| mole volume | $v_m = R^* T/P$ | 43 |
| pressure scale height | $H_P = R^* T_M/(g M_0)$ | 44 |
| mean particle speed | $V = (8 R^* T_M/(\pi M_0))^{1/2}$ | 46 |
| mean free path | $L = \sqrt{2}\, R^* T/(2\pi N_A \sigma^2 P)$ | 47 |
| collision frequency | $\nu = V/L$ | 48, 49 |
| speed of sound | $C_s = (\gamma R^* T_M/M_0)^{1/2}$ | 50 |
| dynamic viscosity | $\mu = \beta T^{3/2}/(T + S)$ | 51 |
| kinematic viscosity | $\eta = \mu/\rho$ | 52 |
| thermal conductivity | $k_t = 2.646\,38 \times 10^{-3}\, T^{3/2}/(T + 245.4 \times 10^{-12/T})$ | 53 |

*Remark.* Table 2 and the discussion on p. 4 print Sutherland's constant as "110 K"; the text of eq. 51 (p. 19) gives 110.4 K, and the tables are computed with it (1.7894 × 10⁻⁵ N s m⁻² at sea level in Table III). The module uses 110.4 K.

*Remark.* $R^*$ is the 1962 value and "is not exactly consistent with the cited values of $k$ and $N_A$" (p. 3): $k N_A/R^* - 1 = 2.3 \times 10^{-6}$. The module computes $N$ from $N_A$ and $R^*$, as eq. 41 does; $P = N k T$ then holds to that inconsistency.

*Remark.* Table III prints the thermal conductivity at sea level as 2.5362 × 10⁻⁵ in units of J m⁻¹ s⁻¹ K⁻¹; eq. 53 gives 2.5326 × 10⁻² W m⁻¹ K⁻¹. The table's values are 0.14 per cent above eq. 53 throughout and three powers of ten small, as if in kJ; the module follows eq. 53, and its conductivity is checked against an implementation of eq. 53 (fluids), not against Table III.

## 7. The Pressure Altitude

Proposition 2 is solved for $H$ in the layer whose base pressures bracket $P$:

$$ H = H_b + \frac{T_{M,b}}{L_{M,b}}\left[\left(\frac{P}{P_b}\right)^{-R^* L_{M,b}/(g_0' M_0)} - 1\right], \qquad H = H_b - \frac{R^* T_{M,b}}{g_0' M_0} \ln\frac{P}{P_b}. $$

This is the height an altimeter set to the standard reports, and the standard height of a radiosonde's pressure level.

## 8. Above 86 km

Above $Z_7 = 86$ km the 1976 no longer defines $T_M$: the kinetic temperature is given directly in geometric height in four segments (isothermal to 91 km, an elliptical arc to 110 km, linear to 120 km, and an exponential approach to 1000 K above; eqs. 24–31), the number density of each of N₂, O, O₂, Ar, He and H follows its own diffusion equation with eddy and molecular diffusion (eqs. 6–8, 35–40), and the pressure is the sum of the partial pressures (eq. 33c). The module does not implement this model, and returns no value above 86 km rather than an extrapolation of layer 6, which would be wrong in temperature by hundreds of kelvin within 50 km.

## 9. Verification

The checks hold the module to:

1. the document's printed Tables I–III, transcribed at thirty-eight heights from −1 to 86 km. The rounded columns agree to half a unit of their last printed digit on Table I's geopotential-altitude pages, and to one unit on the geometric-altitude pages, where some derived columns (density at 2 and 85 km, number density and collision frequency from 50 to 80 km) depart from the document's own equations by that much. Table I's pressure in millibars is *truncated*, not rounded: the equations' value less the printed one lies between 0 and one unit at all thirty-six heights (226.9996 mb is printed 226.99 at 11 km; 845.560 mb is printed 845.55 at 1500 m′ beside a $P/P_0$ of 0.834 50, which is 845.560 mb);
2. fluids' implementation of the 1976 (Bell; MIT licence), to 10⁻¹³ in $T_M$, $P$, $\rho$, $C_s$ and $g$ to 86 km;
3. ambiance's implementation of the ICAO 1993 atmosphere (Dettmann; Apache 2.0) to 80 km: to rounding in $T$, $\mu$ and $g$; within 1.5 × 10⁻⁵ in pressure and density, which ICAO computes from base pressures rounded to six figures; within 10⁻⁶ in the speeds and scale height, from ICAO's $R = 287.052\,87$ J kg⁻¹ K⁻¹; within 8 × 10⁻⁵ in $N$, $L$ and $\nu$, from ICAO's $N_A = 6.022\,57 \times 10^{23}$ mol⁻¹;
4. the hydrostatic equation $\mathrm{d}P/\mathrm{d}Z = -\rho g$ by central differences at 2000 heights (to 4 × 10⁻¹⁰), continuity at the bases of the layers, Table 8 and eq. 25.

The domain's JavaScript standard atmosphere (`js/settling.js`, ICAO's constants, to 20 km) is held to the same fluids references: to rounding in temperature and viscosity and within 3 × 10⁻⁶ in pressure and density, so it agrees with the C++ to that.

## References

- NOAA, NASA and USAF (1976). *U.S. Standard Atmosphere, 1976*. NOAA-S/T 76-1562, U.S. Government Printing Office, Washington, D.C. NASA-TM-X-74335, https://ntrs.nasa.gov/citations/19770009539.
- ICAO (1993). *Manual of the ICAO Standard Atmosphere (extended to 80 kilometres)*, Doc 7488/3, 3rd ed. International Civil Aviation Organization, Montreal.
