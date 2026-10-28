---
title: The Magnetism of a Body without a Dynamo
status: draft
date: 2026-10-29
implemented-by: [magnetism/include/planetary_science/magnetism/crust.hpp, magnetism/include/planetary_science/magnetism/solar_wind.hpp]
---

# The Magnetism of a Body without a Dynamo

*Status: draft, 2026-10-29.* This note states the model of the Moon's magnetic field that the module `magnetism` implements, and the reasons for each of its parts. The general mathematics, potential fields outside a sphere and random fields with a given spectrum, is physis's ([analysis](https://github.com/huangwanhong-maker/physis/tree/main/analysis), decision 0010); this note names the body.

## 1. The Crustal Field

**Definition 1.1.** The Moon's field at and above its surface, from its crust, is $B=-\nabla V$ with the internal potential
$$
V(r,\psi,\lambda)=a\sum_{n=1}^{449}\Big(\frac ar\Big)^{n+1}\sum_{m=0}^{n}\big(g_n^m\cos m\lambda+h_n^m\sin m\lambda\big)\,\breve P_n^m(\sin\psi),\qquad r\ge a=1737.4\ \text{km},
$$
in the Moon's mean-Earth frame, with the coefficients of T2015_449 in tesla.

T2015_449 is Wieczorek's expansion of Tsunakawa et al. (2015)'s surface vector mapping of Kaguya's and Lunar Prospector's magnetometers: an assimilated model, the estimate of the potential to its degree. It resolves scales of about $\pi a/449\approx12$ km.

*Remark.* The Moon has no dynamo today; its crust was magnetised while it had one, and impacts and their ejecta have made the anomalies we see. The field is static over any human time: there is no history of it to learn from, only an inverse problem from orbit.

## 2. The Law of Its Spectrum, and What It Leaves Unresolved

**Definition 2.1 (the spectral law).** The Lowes–Mauersberger spectrum at $a$ is $R_n=(n+1)\sum_m[(g_n^m)^2+(h_n^m)^2]$. Its law on the band $[30,250]$ is $\hat R_n=A\,n^\beta\rho^{2n}$, fitted by least squares on $\log R_n$.

For sources at a depth $d$ the factor $\rho^{2n}$ is the attenuation $((a-d)/a)^{2n}$ of their field at the sphere, so $d=a(1-\rho)$ (Voorhies, Sabaka and Purucker 2002). T2015_449 gives $\beta=2.04$, $\rho=0.99621$ and $d=6.6$ km.

*Caveat.* This $d$ is an effective depth, of the sources and of the model's regularisation together, not a measured one. The band begins above the largest scales, whose power is not the crust's statistics alone, and ends before the degrees the regularisation smooths.

**Definition 2.2 (the unresolved part).** Let the prior of the coefficients be isotropic and Gaussian with the law's power. The model is the posterior mean, and by the law of total variance the residual's power is, on average, the law's less the model's. The unresolved part is drawn independently for each coefficient, $u_n^m\sim\mathcal N(0,\sigma_n^2)$, with
$$
(n+1)(2n+1)\,\sigma_n^2=\max\big(0,\hat R_n-R_n^{\text{model}}\big),\qquad R_n^{\text{model}}=0\ \text{for}\ n>449,
$$
from a counter-based stream keyed by the session's seed and the Moon.

**Proposition 2.3.** $\mathbb E[R_n(u)]=\max(0,\hat R_n-R_n^{\text{model}})$. Above degree 449 the draw carries the law's whole power: the mean of $R_n(u)/\hat R_n$ over $n\in[450,720]$ is 1, with a standard error of $\big(\sum_n 2/(2n+1)\big)^{1/2}/271$.

*Proof.* Degree $n$ has $2n+1$ coefficients of variance $\sigma_n^2$, weighted by $n+1$, so $\mathbb E[R_n(u)]=(n+1)(2n+1)\sigma_n^2$. $R_n(u)/\mathbb E R_n(u)$ is a $\chi^2_{2n+1}/(2n+1)$, of variance $2/(2n+1)$, independently for each $n$. $\square$

*Caveat.* The Moon's crust is neither stationary nor isotropic: its anomalies cluster about the South Pole–Aitken basin and at the swirls. The drawn part says how much finer structure there is, and of what spectrum, not where it is.

## 3. Low Ground

The expansion holds outside its sources. Continued downward to $r<a$, its degree $n$ grows as $(a/r)^{n+2}$: at 2 km below the sphere by $e^{451\cdot0.00115}\approx1.7$ at degree 449, and without bound for the unresolved scales. The maria lie kilometres below the sphere, the floor of South Pole–Aitken some 9 km (Smith et al. 2010). A mount on low ground therefore reads the model at the sphere above it, and its reading says how far below the sphere it is.

## 4. The Field of the Solar Wind

**Definition 4.1 (the carriage from L1).** Let $\hat x$ be the direction of the Sun from the Earth and $x_\parallel=\hat x\cdot x_{\text{body}}$ the body's distance sunward of the Earth. By Taylor's frozen-in hypothesis the interplanetary field is advected unchanged at the solar wind's speed $v$. So the body sees the field measured at L1, at $x_{L1}=1.5\times10^9$ m, a travel time
$$
\tau=\frac{x_{L1}-x_\parallel}{v}
$$
before. For the Moon at 400 km/s, $\tau$ ranges from 2800 s at new Moon to 4700 s at full Moon.

The Earth's magnetotail is taken as a cylinder of 25 Earth radii about the anti-Sun line; a body with $x_\parallel<0$ within it is in the tail, where there is no solar wind, and the part is left out.

**GSE.** The field is measured in GSE axes: $x=\hat x$, $z$ toward the ecliptic's pole and orthogonal to $x$, $y=z\times x$. The module builds them from the Sun's direction and the ecliptic's pole in any frame, and carries a measured vector into that frame.

*Caveat.* Ballistic carriage ignores the tilt of the field's phase fronts, which moves the arrival by minutes (Mailyan, Munteanu and Haaland 2008). The Moon's wake and the plasma's interaction with its strongest anomalies, which stand off the solar wind above Reiner Gamma (Lue et al. 2011), are not modelled. The tail's shape is a cylinder only to first order, and its lobes' field is not given.

## 5. References

See the [README](../README.md), §6.
