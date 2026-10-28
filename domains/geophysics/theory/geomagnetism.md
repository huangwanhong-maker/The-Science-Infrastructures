---
title: The Geomagnetic Field of the World Magnetic Model
status: draft
implemented-by:
  - domains/geophysics/geomagnetism/include/geophysics/geomagnetism/wmm.hpp
  - foundations/calx/numerics/special_functions/include/calx/special_functions/legendre.hpp (the Schmidt functions of §2; geomagnetism/schmidt.hpp gives them the module's names)
references: [Chulliat2025, NCEI2024WMM, Langel1987, AbramowitzStegun1972, HeiskanenMoritz1967, NIMA2000, McIlwain1961, CODATA2022]
---

# The Geomagnetic Field of the World Magnetic Model

*Serves use and learning.* The module [`geomagnetism/`](../geomagnetism/) evaluates the World Magnetic Model 2025 (WMM2025) of NOAA's National Centers for Environmental Information and the British Geological Survey: the Earth's main magnetic field and its secular variation at a geodetic place and an instant, the field on the Earth-fixed axes, the centred dipole with its McIlwain $L$, and the model's own statement of where and when it holds. It follows the model's technical report (Chulliat et al. 2025, cited below as *the report*) equation by equation; this note derives what the report states, says where the module departs from the letter of the report and why, and lists the checks (§9).

**Conventions.** SI throughout the code: tesla, radians, metres; rates per year of the decimal-year scale (§5.3). The coefficients and the report are in nanotesla, and the coefficients are kept so. A place is geodetic on WGS 84: latitude $\varphi$, longitude $\lambda$ east of Greenwich, height $h$ along the ellipsoid's normal. North, east and down are the place's own; azimuths, and the declination, run from north through east.

| Symbol | Meaning | Unit | Code |
|---|---|---|---|
| $a$ | geomagnetic reference radius, 6 371 200 m | m | `Coefficients::reference_radius_m` |
| $g_n^m, h_n^m$; $\dot g_n^m, \dot h_n^m$ | Gauss coefficients at $t_0$; their secular variation | nT; nT/yr | `g_nT`, `h_nT`; `g_nT_per_year`, `h_nT_per_year` |
| $t_0$ | the model's epoch, 2025.0 | yr | `epoch_year` |
| $\breve P_n^m$ | Schmidt semi-normalised associated Legendre function | 1 | `SchmidtFunctions::p` |
| $\lambda, \varphi', r$ | geocentric longitude, latitude, radius | rad, rad, m | `GeocentricPosition` |
| $\psi = \varphi' - \varphi$ | geocentric minus geodetic latitude | rad | (in `field`) |
| $X', Y', Z'$ | north, east, down on the sphere | T | `GeocentricField::north_T` … |
| $X, Y, Z$ | north, east, down on the ellipsoid | T | `MagneticElements::north_T` … |
| $H, F$ | horizontal and total intensity | T | `horizontal_T`, `total_T` |
| $I, D$ | inclination (down positive), declination (east positive) | rad | `inclination_rad`, `declination_rad` |
| $GV$ | grid variation | rad | `grid_variation_rad` |
| $\mathbf m$, $M$ | the centred dipole's moment and its strength | A m² | `Dipole::moment_vector_A_m2()`, `moment_A_m2` |
| $\lambda_m$, $L$ | geomagnetic latitude; McIlwain's shell parameter of the dipole | rad; 1 | `geomagnetic_latitude_rad`, `mcilwain_l` |

## 1. The Main Field as a Potential Field

The field measured at the Earth's surface is the sum of a field from the liquid outer core, by far the largest, a field from magnetised rocks of the crust and upper mantle, and fields from currents in the ionosphere and magnetosphere with those they induce in the Earth and the oceans. The WMM represents the first two at wavelengths longer than about 3 200 km, the *main field*, and leaves out the external fields by its choice of data (report §1.1).

Between the Earth's surface and the ionosphere the air carries no appreciable current, so there $\nabla \times \mathbf B = 0$, and always $\nabla \cdot \mathbf B = 0$. A curl-free field in a simply connected shell is a gradient, $\mathbf B = -\nabla V$, and then $\nabla^2 V = 0$.

**Definition 1.1 (the internal expansion).** The main field is $\mathbf B_m = -\nabla V$ (report eq. 3) with

$$ V(\lambda, \varphi', r, t) = a \sum_{n=1}^{N} \Big(\frac{a}{r}\Big)^{n+1} \sum_{m=0}^{n} \big(g_n^m(t)\cos m\lambda + h_n^m(t)\sin m\lambda\big)\, \breve P_n^m(\sin\varphi'), \tag{1.1} $$

with $N = 12$ (report eq. 4), in geocentric spherical coordinates.

**Proposition 1.2.** Each term of (1.1) is harmonic in $r > 0$, and decays at infinity.

*Proof.* In spherical coordinates, $\nabla^2 \big(r^{-(n+1)} Y\big) = r^{-(n+3)} \big(n(n+1) Y + \Delta_{S^2} Y\big)$ for a function $Y$ on the sphere, and $Y = \breve P_n^m(\sin\varphi')\cos m\lambda$ (or $\sin m\lambda$) is a spherical harmonic of degree $n$: $\Delta_{S^2} Y = -n(n+1) Y$, which is the associated Legendre equation for $P_{n,m}$ after separating $\lambda$ (Heiskanen and Moritz 1967, ch. 1). ∎

*Remark.* The sum starts at $n = 1$: a term of degree 0, $a^2 g_0^0/r$, would be a magnetic monopole, whose flux through a sphere $4\pi a^2 g_0^0$ would contradict $\nabla\cdot\mathbf B = 0$. Terms growing as $r^n$, of external origin, are left out by construction.

## 2. Schmidt Semi-Normalised Functions

Write $\mu = \sin\varphi'$, $\nu = \cos\varphi'$, and let $P_{n,m}(\mu) = (1 - \mu^2)^{m/2}\, d^m P_n(\mu)/d\mu^m$, the associated Legendre functions without the Condon–Shortley phase, $P_{n,m} = (-1)^m P_n^m$ in the notation of Abramowitz and Stegun (1972, ch. 8).

**Definition 2.1 (report eq. 5; Langel 1987, eq. 8 as the report cites it).**

$$ \breve P_n^m(\mu) = \sqrt{\frac{2\,(n-m)!}{(n+m)!}}\; P_{n,m}(\mu) \quad (m > 0), \qquad \breve P_n^0 = P_{n,0} = P_n. \tag{2.1} $$

**Proposition 2.2 (normalisation).** $\displaystyle\int_{-1}^{1} \breve P_n^m(\mu)^2\, d\mu = \frac{2\,(2 - \delta_{m0})}{2n + 1}$.

*Proof.* $\int_{-1}^{1} P_{n,m}^2\, d\mu = \frac{2}{2n+1}\frac{(n+m)!}{(n-m)!}$ (Abramowitz and Stegun 1972, §8.14); multiply by the square of the factor in (2.1). ∎

*Rationale.* The mean over the sphere of $\big(\breve P_n^m(\sin\varphi')\cos m\lambda\big)^2$ is then $1/(2n+1)$ for every order $m$, $m = 0$ included, since the mean of $\cos^2 m\lambda$ is $\tfrac12$ for $m > 0$ and 1 for $m = 0$. Coefficients of one degree are thereby on one scale, and $g_n^m$, $h_n^m$ measure directly how much of the field each harmonic carries. This is Schmidt's semi-normalisation, the convention of geomagnetism; the full normalisation of geodesy would make the mean 1.

**Proposition 2.3 (recurrences).** With $\breve P_{m-1}^m := 0$,

$$ \breve P_0^0 = 1,\quad \breve P_1^1 = \nu,\quad \breve P_m^m = \sqrt{\tfrac{2m-1}{2m}}\;\nu\,\breve P_{m-1}^{m-1}\ (m \ge 2), \tag{2.2} $$

$$ \breve P_n^m = \frac{(2n-1)\,\mu\,\breve P_{n-1}^m - \sqrt{(n-1)^2 - m^2}\;\breve P_{n-2}^m}{\sqrt{n^2 - m^2}} \quad (n > m). \tag{2.3} $$

*Proof.* (2.2): $P_{m,m} = (2m-1)!!\,\nu^m$, so $\breve P_m^m = \sqrt{2/(2m)!}\,(2m-1)!!\,\nu^m$ for $m \ge 1$, whence $\breve P_1^1 = \nu$, and for $m \ge 2$ the ratio $\breve P_m^m/\breve P_{m-1}^{m-1} = (2m-1)\sqrt{(2m-2)!/(2m)!}\;\nu = \sqrt{(2m-1)/(2m)}\;\nu$. (2.3): at fixed $m$ the functions satisfy $(n - m)P_{n,m} = (2n-1)\mu P_{n-1,m} - (n+m-1)P_{n-2,m}$ (Abramowitz and Stegun 1972, §8.5; the phase $(-1)^m$ is common to all three terms). Put $P_{k,m} = c_k \breve P_k^m$ with $c_k = \sqrt{(k+m)!/(2(k-m)!)}$ (or 1 when $m = 0$). Then $c_{n-1}/c_n = \sqrt{(n-m)/(n+m)}$ and $c_{n-2}/c_n = \sqrt{(n-m)(n-m-1)/((n+m)(n+m-1))}$, and dividing the recurrence by $(n - m)c_n$ gives (2.3), since $(n+m-1)(n-m-1) = (n-1)^2 - m^2$. ∎

**Corollary 2.4 (derivatives in latitude).** Differentiating (2.2) and (2.3) in $\varphi'$, with $d\mu/d\varphi' = \nu$ and $d\nu/d\varphi' = -\mu$, gives recurrences for $d\breve P_n^m/d\varphi'$ that divide by nothing, and so hold at the poles:

$$ \frac{d\breve P_m^m}{d\varphi'} = k_m\Big(\nu\,\frac{d\breve P_{m-1}^{m-1}}{d\varphi'} - \mu\,\breve P_{m-1}^{m-1}\Big), \qquad \frac{d\breve P_n^m}{d\varphi'} = \frac{(2n-1)\big(\nu\breve P_{n-1}^m + \mu\, d\breve P_{n-1}^m/d\varphi'\big) - \sqrt{(n-1)^2-m^2}\; d\breve P_{n-2}^m/d\varphi'}{\sqrt{n^2-m^2}}, \tag{2.4} $$

with $k_1 = 1$ and $k_m = \sqrt{(2m-1)/(2m)}$.

**Proposition 2.5 (the report's eq. 16).** $\dfrac{d\breve P_n^m}{d\varphi'} = (n+1)\tan\varphi'\;\breve P_n^m - \sqrt{(n+1)^2 - m^2}\;\sec\varphi'\;\breve P_{n+1}^m$ for $|\varphi'| < \pi/2$.

*Proof.* $(1 - \mu^2)\,dP_{n,m}/d\mu = (n+1)\mu P_{n,m} - (n - m + 1)P_{n+1,m}$ (Abramowitz and Stegun 1972, §8.5). Since $dP/d\varphi' = \nu\, dP/d\mu$, divide by $\nu$, and rescale with $c_{n+1}/c_n = \sqrt{(n+1+m)/(n+1-m)}$: $(n-m+1)\,c_{n+1}/c_n = \sqrt{(n+1)^2 - m^2}$. ∎

The module uses (2.4), not Proposition 2.5, which needs the degree $N + 1$ and divides by $\cos\varphi'$; the checks hold the two together.

**Proposition 2.6 (the east component's functions).** For $m \ge 1$, $\breve P_n^m/\nu$ is a polynomial in $\mu$ times $\nu^{m-1}$; it obeys (2.2) and (2.3) started from $\breve P_1^1/\nu = 1$; and at the poles

$$ \lim_{\varphi' \to \pm\pi/2} \frac{\breve P_n^1(\sin\varphi')}{\cos\varphi'} = (\pm 1)^{n+1}\sqrt{\frac{n(n+1)}{2}}, \qquad \lim_{\varphi' \to \pm\pi/2} \frac{\breve P_n^m(\sin\varphi')}{\cos\varphi'} = 0 \ (m \ge 2). \tag{2.5} $$

*Proof.* $P_{n,m} = \nu^m\, d^m P_n/d\mu^m$, and (2.3) is linear with coefficients free of $\nu$, so dividing every term by $\nu$ commutes with it. For $m = 1$, $P_{n,1}/\nu = P_n'(\mu)$, and $P_n'(1) = n(n+1)/2$, $P_n'(-1) = (-1)^{n+1} n(n+1)/2$; the factor of (2.1) is $\sqrt{2/(n(n+1))}$. For $m \ge 2$ a factor $\nu^{m-1}$ remains. ∎

At the north pole these are the limits $1, \sqrt3, \sqrt6, \sqrt{10}, \ldots, \sqrt{78}$ that the report tabulates (§1.4).

## 3. From the Ellipsoid to the Sphere

The coefficients belong to the sphere; a place is given on the ellipsoid. With $A = 6\,378\,137$ m, $1/f = 298.257223563$ and $e^2 = f(2-f)$ (WGS 84: NIMA 2000, Table 3.1; report eq. 8), and $R_c = A/\sqrt{1 - e^2\sin^2\varphi}$ the radius of curvature in the prime vertical,

$$ p = (R_c + h)\cos\varphi, \quad z = \big(R_c(1-e^2) + h\big)\sin\varphi, \quad r = \sqrt{p^2 + z^2}, \quad \varphi' = \arcsin(z/r) \tag{3.1} $$

(report eq. 7): $p$ and $z$ are the place's distance from the rotation axis and from the equatorial plane, so $(\lambda, \varphi', r)$ are its spherical coordinates. The module takes $\varphi' = \operatorname{atan2}(z, p)$, the same angle, which keeps full precision where $z/r$ is near 1.

*Remark.* NOAA's C library (`GeomagnetismLibrary.c`, `MAG_SetDefaults`) takes $e^2 = 1 - b^2/A^2$ with $b$ rounded to 6 356.752 3142 km, which changes $e^2$ by $1.4 \times 10^{-11}$. The report's numerical example (its Table 3b) is met to $10^{-11}$ rad in $\varphi'$ with $e^2 = f(2-f)$, as its eq. 8 states, and the module follows it, by geodesy's ellipsoid.

## 4. The Components and the Elements

### 4.1 The Gradient on the Sphere

In spherical coordinates the gradient has the components $\big(\partial_r V,\ r^{-1}\partial_{\varphi'} V,\ (r\cos\varphi')^{-1}\partial_\lambda V\big)$ along up, north and east. So $\mathbf B = -\nabla V$ has north, east and down components

$$ X' = -\frac1r\frac{\partial V}{\partial\varphi'}, \qquad Y' = -\frac{1}{r\cos\varphi'}\frac{\partial V}{\partial\lambda}, \qquad Z' = \frac{\partial V}{\partial r}, \tag{4.1} $$

and differentiating (1.1) term by term,

$$ X' = -\sum_n \Big(\frac ar\Big)^{n+2} \sum_m \big(g_n^m\cos m\lambda + h_n^m \sin m\lambda\big)\frac{d\breve P_n^m}{d\varphi'}, \quad Y' = \frac{1}{\cos\varphi'}\sum_n \Big(\frac ar\Big)^{n+2}\sum_m m\big(g_n^m\sin m\lambda - h_n^m\cos m\lambda\big)\breve P_n^m, $$
$$ Z' = -\sum_n (n+1)\Big(\frac ar\Big)^{n+2}\sum_m \big(g_n^m\cos m\lambda + h_n^m\sin m\lambda\big)\breve P_n^m \tag{4.2} $$

(report eqs. 10–12). In $Y'$ the module sums $m\,(\ldots)\,\breve P_n^m/\cos\varphi'$ with the functions of Proposition 2.6, which never divide by $\cos\varphi'$.

### 4.2 The Rotation to the Ellipsoid

The geodetic and geocentric frames of a place share the east axis, and their north and down axes differ by a turn through $\psi = \varphi' - \varphi$ about it.

**Proposition 4.1 (report eq. 17).** $X = X'\cos\psi - Z'\sin\psi$, $Y = Y'$, $Z = X'\sin\psi + Z'\cos\psi$.

*Proof.* On the Earth-fixed axes, north at latitude $\beta$ is $\hat{\mathbf n}(\beta) = (-\sin\beta\cos\lambda, -\sin\beta\sin\lambda, \cos\beta)$ and down is $\hat{\mathbf d}(\beta) = -(\cos\beta\cos\lambda, \cos\beta\sin\lambda, \sin\beta)$. Then $\hat{\mathbf n}(\varphi)\cdot\hat{\mathbf n}(\varphi') = \cos\psi$, $\hat{\mathbf n}(\varphi)\cdot\hat{\mathbf d}(\varphi') = -\sin\psi$, $\hat{\mathbf d}(\varphi)\cdot\hat{\mathbf n}(\varphi') = \sin\psi$, $\hat{\mathbf d}(\varphi)\cdot\hat{\mathbf d}(\varphi') = \cos\psi$; project $\mathbf B = X'\hat{\mathbf n}(\varphi') + Y'\hat{\mathbf e} + Z'\hat{\mathbf d}(\varphi')$. ∎

The same projection, read backwards, gives the field on the Earth-fixed axes, $\mathbf B = X\hat{\mathbf n}(\varphi) + Y\hat{\mathbf e} + Z\hat{\mathbf d}(\varphi)$, with the basis of geodesy's local frame (`ecef_to_ned`); the axes are those of WGS 84, which coincide with the ITRS's far within this model's accuracy. A vector on those axes can be carried to any other frame by the infrastructure's frame graph.

### 4.3 The Elements

$$ H = \sqrt{X^2 + Y^2}, \qquad F = \sqrt{H^2 + Z^2}, \qquad I = \operatorname{atan2}(Z, H), \qquad D = \operatorname{atan2}(Y, X) \tag{4.3} $$

(report eq. 19), so that $I \in [-\pi/2, \pi/2]$, positive where the field dips down, and $D \in (-\pi, \pi]$, positive east of true north. At a dip pole, $H = 0$ and $D$ is undefined.

Near the geographic poles $D$ turns with longitude and ceases to guide; the military specification (MIL-PRF-89500B, as the report gives it, §1.1.2) uses instead the grid variation

$$ GV_N = D - \lambda \ \ (\varphi > 55^\circ), \qquad GV_S = D + \lambda \ \ (\varphi < -55^\circ), \tag{4.4} $$

(report eq. 1), the angle from the grid north of the polar stereographic grid to the horizontal field, reduced to $(-\pi, \pi]$ as the report's Table 6 gives it. In general $GV = D - C$, with $C$ the convergence of the meridians of a chosen grid (report eq. 2); the module computes (4.4) only, taking $|\varphi| \ge 55^\circ$ as NOAA's software does, and returns no value elsewhere.

## 5. Secular Variation

### 5.1 The Linear Model

**Definition 5.1 (report eq. 9).** $g_n^m(t) = g_n^m(t_0) + (t - t_0)\,\dot g_n^m(t_0)$, and likewise $h_n^m$, with $t_0 = 2025.0$.

Since (4.2) is linear in the coefficients, at a fixed place $X'$, $Y'$, $Z'$ are affine in $t$, and their rates are (4.2) with $\dot g$, $\dot h$ in place of $g$, $h$ (report eqs. 13–15); Proposition 4.1 carries them to $\dot X, \dot Y, \dot Z$ unchanged (report eq. 18), $\psi$ not depending on time.

*Caveat.* The report's eq. 15 prints $\dot g_n^m\cos m\lambda - \dot h_n^m\sin m\lambda$ in $\dot Z'$. The derivative of its eq. 12 has a plus, as its eqs. 13 and 14 do in their places, and as NOAA's software has (`MAG_SecVarSummation`). The module uses the plus; its checks hold $\dot Z$ to the derivative of $Z$ in time and to NOAA's test values.

### 5.2 The Rates of the Elements

Differentiating (4.3), with $d\operatorname{atan2}(y, x) = (x\,dy - y\,dx)/(x^2 + y^2)$:

$$ \dot H = \frac{X\dot X + Y\dot Y}{H}, \quad \dot F = \frac{X\dot X + Y\dot Y + Z\dot Z}{F}, \quad \dot I = \frac{H\dot Z - Z\dot H}{F^2}, \quad \dot D = \frac{X\dot Y - Y\dot X}{H^2}, \quad \dot{GV} = \dot D \tag{5.1} $$

(report eq. 20). These are exact derivatives of the elements, which, unlike the components, are not affine in time.

### 5.3 The Decimal Year

The model's time is the decimal year. The module follows NOAA's software (`MAG_DateToYear`): a Gregorian date with a fraction $f$ of its day is

$$ t = Y + \frac{d - 1 + f}{365 \text{ or } 366}, \tag{5.2} $$

with $d$ the day of the year. So a "year" of the rates is a calendar year, of 365 or 366 days.

*Remark.* The date is the civil date, UTC, but the scale does not matter: $TT - UTC = 69.184$ s is $2.2 \times 10^{-6}$ year, over which the field changes by less than $5 \times 10^{-4}$ nT, since the secular variation stays below 230 nT per year at the surface (its greatest, 223 nT per year, is off the coast of Brazil; §9). The convention does matter: decimal years counted in Julian years of 365.25 days from J2000.0 stray from (5.2) by up to half a day in 2025–2030, and the field by up to 0.3 nT, more than the 0.1 nT the report allows an implementation (§1.10).

## 6. The Geographic Poles

At a pole the longitude is arbitrary, and with it the north and east of the frame; the field is regular there, its components not. The report resolves this by continuity (§1.4): a pole given the longitude $\lambda$ has the north of a point approaching it along the meridian $\lambda$. With Proposition 2.6, (4.2) then holds at the pole itself, and the module evaluates it there with no special case.

*Caveat.* The report gives, at $r = a$ on 1 January 2025 with the pole given $\lambda = 0$, $(X', Y', Z') = (1709.5, 417.8, 56517.5)$ nT at the north pole and $(14185.8, -8732.7, -51366.5)$ nT at the south pole, which the module meets to 0.05 nT, and says that the field lies in the half-planes of the meridians 166.27° W and 31.62° W. In its own frame, the north of the north pole along $\lambda = 0$ is $-\hat{\mathbf x}$ and its east $+\hat{\mathbf y}$, so the horizontal field $X'(-\hat{\mathbf x}) + Y'\hat{\mathbf y}$ points to the meridian $\operatorname{atan2}(417.8, -1709.5) = 166.27^\circ$ *east*; at the south pole, whose north is $+\hat{\mathbf x}$, $\operatorname{atan2}(-8732.7, 14185.8) = 31.62^\circ$ west, as printed. The checks hold the module to the report's components and its south meridian, and to 166.27° on the east for the north.

## 7. The Centred Dipole and $L$

### 7.1 The Dipole of the First Degree

Since $\breve P_1^0 = \sin\varphi'$ and $\breve P_1^1 = \cos\varphi'$, the degree-1 part of (1.1) is

$$ V_1 = \frac{a^3}{r^2}\big(g_1^0\sin\varphi' + (g_1^1\cos\lambda + h_1^1\sin\lambda)\cos\varphi'\big) = \frac{a^3\,(\mathbf g\cdot\hat{\mathbf r})}{r^2}, \qquad \mathbf g = (g_1^1, h_1^1, g_1^0) \tag{7.1} $$

on the Earth-fixed axes. The scalar potential of a point dipole $\mathbf m$ at the centre, with $\mathbf B = -\nabla V$, is $\mu_0(\mathbf m\cdot\hat{\mathbf r})/(4\pi r^2)$.

**Proposition 7.1.** The degree-1 field is that of the centred dipole $\mathbf m = 4\pi a^3\mathbf g/\mu_0$:

$$ \mathbf B_1 = \Big(\frac ar\Big)^3\big(3(\mathbf g\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf g\big) = \frac{\mu_0}{4\pi r^3}\big(3(\mathbf m\cdot\hat{\mathbf r})\hat{\mathbf r} - \mathbf m\big). \tag{7.2} $$

*Proof.* Compare the potentials; then $-\nabla(\mathbf g\cdot\mathbf r\, r^{-3}) = -\mathbf g\,r^{-3} + 3(\mathbf g\cdot\mathbf r)\mathbf r\, r^{-5}$. ∎

The field's strength on the dipole's equator at $r = a$ is $B_0 = |\mathbf g|$, and $M = |\mathbf m| = 4\pi a^3 B_0/\mu_0$, with $\mu_0$ the CODATA 2022 value (it has been measured, not fixed, since the SI of 2019, and differs from $4\pi\times10^{-7}$ N A⁻² by $5\times10^{-10}$). For WMM2025 at 2025.0, $B_0 = 29\,735.2$ nT and $M = 7.690\times10^{22}$ A m². Since $g_1^0 < 0$, $\mathbf m$ points south; the axis meets the sphere in the northern hemisphere in the direction $-\mathbf g$, at the *geomagnetic north pole*, with geocentric latitude $\arcsin(-g_1^0/B_0)$ and longitude $\operatorname{atan2}(-h_1^1, -g_1^1)$, and the angle $\arccos(-g_1^0/B_0)$ between the axis and the rotation axis. The point of the ellipsoid on the axis has the geodetic latitude $\varphi$ with $\tan\varphi = \tan\varphi'/(1 - e^2)$. The report gives 72.76° W, 80.79° N geocentric (80.85° N geodetic), and a tilt of 9.21° (§1.6, Table 4); the south geomagnetic pole is the antipode.

*Remark.* The geomagnetic poles are not the dip poles, where the whole field is vertical (in 2025 near 85.8° N, 139.3° E and 63.9° S, 135.1° E, report Table 4): the higher degrees move the dip poles, not the dipole's axis. The eccentric dipole, from the first eight coefficients (report §1.6), is not implemented.

### 7.2 Field Lines and McIlwain's $L$

Let $\lambda_m$ be the latitude about the dipole's axis, $\sin\lambda_m = \hat{\mathbf r}\cdot(-\hat{\mathbf g})$. By (7.2) the dipole field has the radial and northward components $B_r = -2B_0(a/r)^3\sin\lambda_m$ and $B_\lambda = B_0(a/r)^3\cos\lambda_m$.

**Proposition 7.2.** The field lines of the centred dipole are the curves $r = L\,a\cos^2\lambda_m$ in the meridian planes of the dipole, one for each $L > 0$.

*Proof.* Along a field line $dr/(r\,d\lambda_m) = B_r/B_\lambda = -2\tan\lambda_m$, so $\ln r = 2\ln\cos\lambda_m + \text{const}$. ∎

So $L = r/(a\cos^2\lambda_m)$ is constant along each field line: it is the distance, in units of $a$, at which the line crosses the dipole's equator. McIlwain (1961) defined $L$ for the real field from the adiabatic invariants of particles trapped in it, and showed that it reduces to this in a dipole; the module computes the dipole's $L$ only, which departs from McIlwain's by several per cent where the field is far from a dipole, as over the South Atlantic.

## 8. The Model's Validity

The report states the domain in which the model holds; the module encodes it, and says in words why a place and instant fall outside it.

- **Time.** "Valid from 2025.0 to 2030.0" (§1.1), "to be used from January 1, 2025 to December 31, 2029" (§1.9): the module takes $t_0 \le t \le t_0 + 5$, as NOAA's software accepts dates up to $t_0 + 5$.
- **Height.** The WMM gives the field "from 1 km below the World Geodetic System 1984 (WGS 84) ellipsoid surface to at least 850 km above it" (§1.1, after MIL-PRF-89500B). In quiet geomagnetic conditions the report finds it within the specification higher, to 1 900 km for $Z$ and further for the other components (§3.5, Table 17), and at no height for some components in strong storms; the module keeps the specified range.
- **Compass zones.** Where $H$ is weak a compass, and the model's $D$, cannot be trusted: the *blackout zone* where $H < 2000$ nT, and the *caution zone* around it where $2000 \le H < 6000$ nT (§1.8). The report draws them on the ellipsoid; NOAA's calculators warn by $H$ at the point asked for, and so does the module, which agrees with the drawing at $h = 0$.
- **Uncertainty.** The model's error model (§3.4) gives constant uncertainties of $X$, $Y$, $Z$ (137, 89, 141 nT; eqs. 27–29), of $H$, $F$ (133, 138 nT) and $I$ (0.20°; eqs. 37–39), and of $D$

  $$ \delta D = \sqrt{0.26^2 + (5417/H)^2}\ \text{degrees}, \quad H \text{ in nT} \tag{8.1} $$

  (eq. 43), which grows without bound towards the dip poles and is 0.29° where $H$ is greatest at the surface, 41 875 nT. They estimate the model's departure from the true field, the crustal and disturbance fields included, not the computation's.

*Caveat.* The report prints the third of eqs. 27–29 as $\delta X$; it is the uncertainty of $Z$, as eqs. 30–32 for WMMHR2025 show.

## 9. Verification

The checks are `geomagnetism/tests/check_wmm.cpp`, on calx's verification, with the references written by `tests/make_wmm_fixtures.py` from the files the data library fetches (`python data/acquisition/fetch.py noaa-wmm`); every number below is a check's. The first three rows, on the Schmidt functions, are now calx's checks (`special_functions/tests/check_legendre.cpp`, which also holds them to the standard library's `std::assoc_legendre`), since the functions moved there.

| Test | Agreement |
|---|---|
| $\int\breve P^2\,d\mu$ against Proposition 2.2, to degree 12 (calx's Gauss–Kronrod) | $9\times10^{-16}$ relative |
| $\breve P_3^m$ against the report's eq. 6; $d\breve P/d\varphi'$ against eq. 16 (Proposition 2.5) | $2\times10^{-16}$; $2\times10^{-14}$ |
| $\breve P/\cos\varphi'$ at the poles against (2.5) | $9\times10^{-15}$ |
| The report's example (Tables 3a, 3b): $\varphi'$, $r$; $X', Y', Z'$ and their rates; $X, Y, Z, H, F$ and rates; $D, I$ and rates | $1.3\times10^{-11}$ rad, $10^{-9}$ m; $5\times10^{-11}$ nT; $4\times10^{-11}$ nT; $3.5\times10^{-11}$ rad |
| The example's grid variation, south grid, and $D$ (Table 3d) | $3\times10^{-11}$° |
| NOAA's 100 test points: $D$, $I$ (printed to 0.01°); $Y$ and the rates of $Y$, $Z$, $D$, $I$ (printed to $10^{-6}$) | within the printing; $5\times10^{-7}$ |
| The same: $X$, $Z$, $H$, $F$; the rates of $X$, $H$, $F$ | $7.2\times10^{-4}$ nT; $1.5\times10^{-6}$ nT/yr (see the caveat) |
| The report's Table 6 (12 points, with $GV$), printed to 0.1 nT and 0.01° | within the printing |
| The field above the poles at $r = a$ (§1.4) | 0.05 nT, within the printing |
| $-\nabla V$ by differences of the potential (20 places from −1 to 850 km); the rates as derivatives in time | $10^{-5}$ nT; $2\times10^{-8}$ nT/yr |
| The Earth-fixed vector back to north, east, down, and from the sphere's basis | $1.4\times10^{-11}$ nT |
| The geomagnetic poles and tilt at 2025.0 (§1.6) | 0.003°, within the printing |
| The degree-1 model as the field (7.2) of $\mathbf m$; $\mathbf B\cdot\nabla L = 0$; $L = 2$ at $r = 2a$ on the dipole's equator | $3\times10^{-11}$ nT; $2\times10^{-10}$; $2\times10^{-16}$ |
| The greatest $H$ on the ellipsoid at 2025.0 against the report's 41 875 nT; $\delta D = 0.29°$ there | 0.3 nT; $4\times10^{-4}$° |
| The compass zones of the 100 points (2 in the blackout zone, 7 in the caution zone) by their printed $H$ | all |
| The greatest $\lvert\dot{\mathbf B}\rvert$ on the ellipsoid (a 1° grid) below 230 nT/yr; decimal years of Julian years within half a day of (5.2), 2025–2030 | 223 nT/yr; 0.5 day |

*Caveat.* NOAA's file of 100 test points departs from the report's equations by up to $7.2\times10^{-4}$ nT in $X$, and so in $H$, $F$, and through $\psi$ slightly in $Z$, while the module meets the report's own high-precision example (Table 3b) to $10^{-10}$ nT. NOAA's C library departs from the file alike: `GeomagnetismLibrary.c` of the WMM2025 Linux software (`wmm2025_Linux.zip` from NCEI), compiled with gcc 11.5 and called (`MAG_GeodeticToSpherical`, `MAG_TimelyModifyMagneticModel`, `MAG_Geomag`) on the 100 points on 2026-10-28, gives $X$ $7.2\times10^{-4}$ nT from the file's and within $5\times10^{-7}$ nT of the module's, the rounding of its printing. The file was evidently made by other code. All of this lies far within the 0.1 nT the report allows an implementation (§1.10), and further still within the model's uncertainty.

## 10. References

- Abramowitz, M. and Stegun, I. A. (eds.) (1972). *Handbook of Mathematical Functions*, 10th printing, ch. 8 (Legendre functions). National Bureau of Standards.
- Chulliat, A., Brown, W., Nair, M., Gómez Pérez, N., Young, L.-Y., Watson, C., Boneh, N., Beggan, C., Meyer, B. and Paniccia, M. (2025). *The US/UK World Magnetic Model for 2025–2030: Technical Report*. NOAA National Centers for Environmental Information. doi:10.25923/prbc-s316.
- Heiskanen, W. A. and Moritz, H. (1967). *Physical Geodesy*, ch. 1. W. H. Freeman.
- Langel, R. A. (1987). The main field. In J. A. Jacobs (ed.), *Geomagnetism*, vol. 1, pp. 249–512. Academic Press. (Cited here as the report cites it, for the definition of $\breve P_n^m$.)
- McIlwain, C. E. (1961). Coordinates for mapping the distribution of magnetically trapped particles. *Journal of Geophysical Research* 66, 3681–3691.
- Mohr, P. J., Newell, D. B., Taylor, B. N. and Tiesinga, E. (2025). CODATA recommended values of the fundamental physical constants: 2022. *Reviews of Modern Physics* 97, 025002. (The value of $\mu_0$, as scipy.constants 1.18 gives it.)
- National Imagery and Mapping Agency (2000). *Department of Defense World Geodetic System 1984*, TR8350.2, 3rd ed., Table 3.1.
- NOAA NCEI Geomagnetic Modeling Team and British Geological Survey (2024). *World Magnetic Model 2025*. doi:10.25921/aqfd-sd83. The coefficients, test values and software, public domain.
