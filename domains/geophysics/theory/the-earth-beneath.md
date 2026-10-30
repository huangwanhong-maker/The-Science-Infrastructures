# The Earth beneath a Place

*Status: draft (2026-10-30). Learning and use.* The notes behind the domain's `structure` in C++ ([ADR 0055](../../../development/decisions/0055-the-earth-beneath.md)): the column of CRUST1.0, the radial Earth of PREM and what follows from it, the column they compose, and sections along geodesics. Code names are in `code`.

## 1. The Crust of CRUST1.0

**Definition 1.1.** CRUST1.0 (Laske, Masters, Ma and Pasyanos 2013) assigns to each cell $c$ of a grid of $1^\circ$ nine heights $t_0(c) \ge t_1(c) \ge \dots \ge t_8(c)$ above sea level, the tops of water, ice, three sediments, three crystalline layers and the mantle, with a speed of P waves $v_P$, of S waves $v_S$ and a density $\rho$ for each. Layer $k$ is *present* when $t_k > t_{k+1}$. Its base is $t_{k+1}$.

The derived quantities (`CrustColumn`):
- the solid surface $t_2$;
- the Moho's depth $-t_8$;
- the crust's thickness $t_2 - t_8$;
- the sediments' thickness $t_2 - t_5$;
- the water's $t_0 - t_1$ and the ice's $t_1 - t_2$.

A column is *oceanic* when it lies under more than 500 m of water and has less than 15 km of rock.

*Remark.* A cell's values are its average, and the model asks to be read so ("please use the midpoint of this cell"). A place therefore takes its cell's column (`CrustGrid::cell`); nothing is interpolated between cells.

**The time through the crust.** A P wave rising vertically from the Moho reaches the solid surface after
$$T_P = \sum_{k \ge 2} \frac{t_k - t_{k+1}}{v_{P,k}},$$
the sum over the sediments and the crystalline layers present (the column's `crust_p_time_s`).

## 2. The Radial Earth of PREM

**Definition 2.1.** PREM (Dziewonski and Anderson 1981, Table I) is a sphere of radius $a = 6371$ km in thirteen shells $r_{k-1} < r \le r_k$. In each, with $x = r/a$,
$$\rho(r) = \sum_{j=0}^{3} c^{\rho}_{kj}\, x^j, \qquad v_P(r) = \sum_j c^{P}_{kj}\, x^j, \qquad v_S(r) = \sum_j c^{S}_{kj}\, x^j .$$
The model is isotropic, at the reference period of 1 s (`PREM_SHELLS`). Its discontinuities are those of its polynomials: the sea floor (3 km), 15 km, the Moho (24.4 km), 220, 400 and 670 km, the core–mantle boundary (2891 km) and the inner core boundary (5149.5 km). At 771 and 2741 km its polynomials change and its values do not jump.

**Proposition 2.2** (mass). The mass within radius $r$, in the shell $k$, is
$$M(r) = 4\pi a^3 \Big[ \sum_{i<k} \sum_j c^{\rho}_{ij}\, \frac{x_i^{\,j+3} - x_{i-1}^{\,j+3}}{j+3} + \sum_j c^{\rho}_{kj}\, \frac{x^{j+3} - x_{k-1}^{\,j+3}}{j+3} \Big],$$
with $x_i = r_i / a$. The moment of inertia about an axis is $I = \tfrac{8\pi}{3} \int_0^a \rho\, r^4\, dr$, in closed form in the same way.

*Proof.* $M(r) = \int_0^r 4\pi \rho(s)\, s^2\, ds$; substitute $s = a\xi$ and integrate each monomial.

The code is `prem_mass_within` and `prem_moment_of_inertia`. Its checks find $M(a)$ within $1.7 \times 10^{-4}$ of $GM_\oplus / G$ (WGS 84's $GM$; CODATA 2018's $G$) and $I = 0.3308\, M a^2$.

**Proposition 2.3** (gravity). In a spherically symmetric Earth, gravity at radius $r$ is $g(r) = G M(r)/r^2$, directed to the centre (`prem_gravity`). It satisfies
$$\frac{dg}{dr} = 4\pi G \Big( \rho(r) - \tfrac{2}{3}\, \bar\rho(r) \Big), \qquad \bar\rho(r) = \frac{3 M(r)}{4\pi r^3},$$
so $g$ grows outward exactly where the local density exceeds two thirds of the mean density within.

*Proof.* Differentiate $G M(r) r^{-2}$, using $M'(r) = 4\pi r^2 \rho(r)$.

*Remark.* In the outer core $\rho$ (9.9 to 12.2 g/cm³) exceeds $\tfrac{2}{3}\bar\rho$, and $g$ rises outward. Just above the core–mantle boundary the mantle's 5.6 g/cm³ falls below $\tfrac{2}{3}\bar\rho(r_{\mathrm{CMB}})$, about 7.3 g/cm³, and $g$ falls. Gravity is therefore greatest at the core–mantle boundary, 10.69 m/s² in PREM. The check finds the maximum there, by search on a grid of 1 km.

## 3. The Column: CRUST1.0 over PREM

**Definition 3.1** (`interior_at`, `interior_profile`). At a place with CRUST1.0's column $C$ and a height $z$ above sea level, the interior is:
- nothing, above the column's top $t_0$ (in the air);
- the layer of $C$ that contains $z$, down to the Moho;
- below the Moho, PREM at depth $d = -z$. Where $d$ is less than PREM's Moho (24.4 km), the lid's polynomial is continued upward.

Depth is taken below mean sea level along the ellipsoid's normal, and PREM's spherical depth is identified with it. This omits the ellipticity correction, an error of the order of the flattening, at most some 21 km near the centre. As fields (`InteriorField`), the three quantities are refused, with the reason, in the air and beyond the centre.

## 4. Sections along Geodesics

**Definition 4.1** (`SectionLine`). The section through a place $P_0$ at azimuth $\alpha_0$ is the geodesic $\gamma$ of the ellipsoid with $\gamma(0) = P_0$, $\gamma'(0)$ at azimuth $\alpha_0$, parametrised by arc length $s$. It is sampled at $s = k\,\Delta s$, $|s| \le h$, by the direct problem (Karney 2013, geodesy's `GeodesicLine`).

**Definition 4.2.** For a place $Q$ near the line:
- its *foot* is $s^\* = \arg\min_s \operatorname{dist}(\gamma(s), Q)$;
- its distance *along* the line is $s^\*$;
- its distance *across* is $\pm \operatorname{dist}(\gamma(s^\*), Q)$, positive when $Q$ lies to the right of the direction of travel, that is, when the azimuth from the foot to $Q$ exceeds the line's azimuth there by an angle whose sine is positive.

At the foot the geodesic to $Q$ meets $\gamma$ at a right angle (the first variation of length).

*Computation.* Brent's minimisation (calx) finds $s^\*$ to within $\sqrt{\varepsilon}\,|s|$, some millimetres at 1000 km. Two steps of the right spherical triangle at the current foot then take it to rounding: with $d$ the distance to $Q$ and $\Delta\alpha$ the angle between the line and the geodesic to $Q$,
$$\tan\frac{\delta s}{R} = \tan\frac{d}{R}\, \cos\Delta\alpha$$
(Napier's rules). Each step leaves an error of the order of the flattening times the step before.

*Check.* A place built by the direct problem 80 km to the right of the line at a right angle, 400 km along it, is found at 400 km along and 80 km across within $10^{-9}$ m. So is one 120 km to the left, 1200 km along.

## References

- Laske, G., Masters, G., Ma, Z. and Pasyanos, M. (2013). Update on CRUST1.0: a 1-degree global model of Earth's crust. *Geophys. Res. Abstracts* 15, EGU2013-2658.
- Dziewonski, A. M. and Anderson, D. L. (1981). Preliminary reference Earth model. *Phys. Earth Planet. Inter.* 25, 297–356.
- Karney, C. F. F. (2013). Algorithms for geodesics. *J. Geodesy* 87, 43–55.
- Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*, ch. 5. Prentice-Hall.
