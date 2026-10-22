---
title: Volcanic Ash Carried by the Wind and Settling Through the Air
status: draft
implemented-by:
  - domains/meteorology/js/settling.js
  - applications/earth-atlas/plugins/volcanic-ash.js
references: [ICAO1993, USStandardAtmosphere1976, Stokes1851, SchillerNaumann1933, SeinfeldPandis2016]
---

# Volcanic Ash Carried by the Wind and Settling Through the Air

## 1. A Falling Parcel

An eruption column lifts ash to a height of several kilometres, and the wind then carries it while it falls. The simplest description that keeps both is a *parcel*: a point that moves horizontally with the wind at its height and falls relative to the air at the terminal velocity of its particles. A parcel has no mass or concentration; its path says where ash of one size goes, not how much arrives.

## 2. The Standard Atmosphere

*Hypotheses.* Dry air, an ideal gas with $R = 287.05287$ J kg⁻¹ K⁻¹, in hydrostatic balance under standard gravity $g_0 = 9.80665$ m s⁻²; the temperature falls linearly with geopotential height $z$ at $L = 6.5$ K/km from $T_0 = 288.15$ K at $p_0 = 101\,325$ Pa, to 11 km, and is constant (216.65 K) from 11 to 20 km (ICAO 1993; the same as the U.S. Standard Atmosphere 1976 below 32 km).

**Proposition 1.** The pressure is

$$ p(z) = p_0 \left(\frac{T_0 - L z}{T_0}\right)^{g_0/(R L)} \;\; (z \le 11\ \mathrm{km}), \qquad p(z) = p_{11}\, e^{-g_0 (z - 11\,\mathrm{km})/(R\, T_{11})} \;\; (11 < z \le 20\ \mathrm{km}), \tag{1} $$

with $g_0/(RL) = 5.25588$ and $p_{11} = 22\,632$ Pa.

*Proof.* Hydrostatic balance, $\mathrm{d}p/\mathrm{d}z = -\rho g_0 = -p g_0/(R T)$, gives $\mathrm{d}\ln p = -g_0\,\mathrm{d}z/(R T(z))$. With $T = T_0 - Lz$, $\mathrm{d}z = -\mathrm{d}T/L$, so $\ln(p/p_0) = (g_0/RL)\ln(T/T_0)$; with $T$ constant, the exponential. ∎

Inverting (1) gives the *pressure altitude* of a level: 850 hPa at 1457 m, 500 hPa at 5574 m, 250 hPa at 10 363 m. The viscosity of air follows Sutherland's law,

$$ \mu(T) = \frac{\beta\, T^{3/2}}{T + S}, \qquad \beta = 1.458 \times 10^{-6}\ \mathrm{kg\,m^{-1}\,s^{-1}\,K^{-1/2}},\; S = 110.4\ \mathrm{K}, \tag{2} $$

(U.S. Standard Atmosphere 1976, eq. 51): $1.7894 \times 10^{-5}$ Pa s at sea level, $1.4216 \times 10^{-5}$ Pa s above 11 km.

## 3. Stokes' Law

*Hypotheses.* A rigid sphere of diameter $d$ and density $\rho_p$ falls steadily through still air of density $\rho_a$ and viscosity $\mu$, slowly enough that inertia is negligible beside viscosity: the Reynolds number $Re = \rho_a w d / \mu \ll 1$.

**Proposition 2 (Stokes 1851).** The drag on the sphere is $3\pi\mu d\, w$, and the terminal velocity is

$$ w = \frac{(\rho_p - \rho_a)\, g\, d^2}{18\, \mu}. \tag{3} $$

*Proof (of (3) from the drag).* At terminal velocity the drag balances the weight less the buoyancy, $3\pi\mu d\, w = (\rho_p - \rho_a)\, g\, \pi d^3/6$. The drag itself is the solution of the Stokes equations about a sphere (Stokes 1851; Seinfeld and Pandis 2016, §9.3). ∎

**Validity.** The drag grows beyond Stokes' as $Re$ grows; Schiller and Naumann's correction multiplies it by $1 + 0.15\, Re^{0.687}$, so (3) overstates $w$ by 3% at $Re = 0.1$ and 15% at $Re = 1$. For volcanic glass ($\rho_p = 2500$ kg m⁻³), $Re = 0.1$ at $d = 27$ µm near the ground and 34 µm at 10 km. The slip of air at a small particle (Cunningham's correction, $1 + 2.52\, \lambda/d$ with $\lambda$ the mean free path, 0.23 µm at 12 km) adds less than 3% for $d > 25$ µm below 12 km and is neglected. Real ash is irregular, and falls more slowly than a sphere of the same volume.

*Example.* $d = 30$ µm, $\rho_p = 2500$ kg m⁻³ at sea level: $w = 2498.775 \times 9.80665 \times (30 \times 10^{-6})^2 / (18 \times 1.7894 \times 10^{-5}) = 0.06847$ m/s, with $Re = 0.14$. At 10 km the thinner air gives 0.0841 m/s.

## 4. The Wind at a Parcel's Height

The forecast gives the wind at 10 m and at 850, 500 and 250 hPa. The wind at height $z$ is linear in $z$ between the levels placed at their pressure altitudes (10, 1457, 5574 and 10 363 m), the 10 m wind below 10 m and the 250 hPa wind above 10 363 m.

*Caveat.* A level stands at its own geopotential height, which departs from the standard atmosphere's by a few hundred metres in middle latitudes and more in the tropics and in winter. The wind between levels is not linear, and the jet stream, near 250 hPa, is often stronger above it than at it. The air's own rising and sinking, a few cm/s on the large scale and much more in storms, are of the size of the settling speed and are left out.

## 5. The Path

**Definition 1.** A parcel at longitude $\lambda$, latitude $\varphi$ (radians) and height $z$ moves by

$$ \frac{\mathrm{d}\lambda}{\mathrm{d}t} = \frac{u}{a \cos\varphi}, \qquad \frac{\mathrm{d}\varphi}{\mathrm{d}t} = \frac{v}{a}, \qquad \frac{\mathrm{d}z}{\mathrm{d}t} = -w(z), \tag{4} $$

with $(u, v)$ the wind at $(\lambda, \varphi, z, t)$ and $a$ the Earth's mean radius, integrated by the classical fourth-order Runge–Kutta method, until $z$ falls to the ground (the library's relief, or sea level). The landing is placed within the last step, where the height above the ground, linear in the step, reaches zero.

**Proposition 3.** In still air, a parcel released at $z_0$ lands where it was released, after the time $\int_0^{z_0} \mathrm{d}z / w(z)$. In a uniform eastward wind $U$ along the equator, it moves $U t$ in time $t$.

*Proof.* In still air (4) leaves $\lambda, \varphi$ fixed and separates in $z$. Along the equator, $\cos\varphi = 1$ and $v = 0$, so $\varphi$ stays 0 and $a\,\mathrm{d}\lambda/\mathrm{d}t = U$. ∎

The Runge–Kutta step is exact for the constant rates of the second case, and of fourth order in the first.

## 6. Verification

| Test | Agreement |
|---|---|
| Standard atmosphere at sea level: 1.225 kg m⁻³, $1.7894 \times 10^{-5}$ Pa s | $10^{-5}$ (relative) |
| At 11 km: 22 632 Pa, 216.65 K | $2 \times 10^{-6}$ |
| Pressure altitudes of 850, 500, 250 hPa: 1457.3, 5574.4, 10 362.9 m | 0.04 m |
| Stokes velocity, 30 µm, 2500 kg m⁻³, sea level: 0.06847 m/s (by hand) | $2 \times 10^{-6}$ m/s |
| Still air: lands at the start, after $\int \mathrm{d}z/w$ from 5 km | Exact; $1.4 \times 10^{-7}$ (relative) |
| Uniform wind along the equator: displacement $U t$ | $2 \times 10^{-7}$ m |
| Plug-in: still air, three release heights to 3 km: landings at the vent, after $\int \mathrm{d}z / w$ from 1 km | Exact; $1.6 \times 10^{-7}$ |
| Plug-in: uniform 10 m/s: landing $U t$ from the vent | $5 \times 10^{-12}$ (relative) |
| Plug-in on GFS from Sakurajima, 12 parcels of 25 µm to 12 km: the low ones land on Kyushu and Shikoku in 6–17 h, the highest 3000 km east after 62 h; 0.15 s in the page | — |

*Caveat.* The plug-in shows where the wind carries ash of one size from a column, not the concentration of ash in the air, which needs the mass erupted, its spread of sizes, turbulence, the aggregation of fine ash into clusters that fall sooner, and rain that washes it out. The Volcanic Ash Advisory Centres model all of these.
