---
title: Daily Weather from a Forecast
status: draft
implemented-by:
  - domains/meteorology/js/fireweather.js
  - domains/meteorology/js/heat.js
  - applications/earth-atlas/plugins/fire-weather.js
  - applications/earth-atlas/plugins/heat.js
  - applications/earth-atlas/plugins/heavy-rain.js
references: [Wilks2011, VanWagner1987]
---

# Daily Weather from a Forecast

## 1. The Day Is Local, the Forecast Universal

The measures of hazard are defined on days: the weather at noon and the rain of the day before (fire weather), the greatest heat of the day (heat stress), the rain of a day or three (heavy rain). A day is local: noon in Tokyo is midnight in London. A forecast is universal: it gives its fields at fixed instants of UTC, every few hours, for the whole Earth at once. This note says how the one is read from the other, and what is lost.

All instants are UTC, in milliseconds in the code. Longitudes $\lambda$ are in degrees east.

## 2. Local Days

**Definition 1.** The *local mean solar time* at longitude $\lambda$ and instant $t$ is $t + (\lambda / 15)$ h. The *local noon* of the UTC calendar day beginning at $d$ is

$$ n(d, \lambda) = d + 12\,\mathrm{h} - (\lambda / 15)\,\mathrm{h}, \tag{1} $$

and its *local day* is the interval $[\, d - (\lambda/15)\,\mathrm{h},\; d + 24\,\mathrm{h} - (\lambda/15)\,\mathrm{h} \,)$.

Mean solar noon differs from apparent noon by the equation of time, at most 16 minutes, and from noon of the standard time zone by up to an hour or two where zones are wide. Over a day of 24 hours these are small; against samples 6 hours apart they are negligible.

*Remark.* At 180° E the local day begins at $d - 12$ h and at 180° W it ends at $d + 36$ h: the "same" day spans 48 hours of UTC over the Earth. A forecast must cover all of it for a map of that day to be whole.

## 3. Samples at Fixed Instants

Let the forecast's times be $t_0 < t_1 < \dots < t_N$, equally spaced by $\Delta$.

**Definition 2.** The *noon sample* at $(d, \lambda)$ is $t_k$ with $k = \arg\min_j |t_j - n(d, \lambda)|$, if $|t_k - n(d, \lambda)| \le \Delta/2$.

**Proposition 1.** Within the forecast's span, the noon sample exists and lies within $\Delta/2$ of noon. With $\Delta = 6$ h, the weather used is that of a local time between 09:00 and 15:00.

*Proof.* Consecutive samples are $\Delta$ apart, so every instant in $[t_0, t_N]$ is within $\Delta/2$ of one; and $\Delta/2 = 3$ h. ∎

**Proposition 2 (whole days).** Call the day $d$ *whole* for noon weather if at every longitude the rain of the 24 hours to noon lies within the samples' intervals (§4) and a noon sample exists; this holds if and only if $d - 24\,\mathrm{h} \ge t_0 - \Delta$ and $d + 24\,\mathrm{h} \le t_N$. Call it whole for daily maxima if at every longitude each instant $t_0 + j\Delta$ ($j \in \mathbb{Z}$) within the local day is a sample; this holds if and only if $d - 12\,\mathrm{h} > t_0 - \Delta$ and $d + 36\,\mathrm{h} \le t_N + \Delta$.

*Proof.* Noon is earliest at 180° E, $n(d, 180) = d$, and latest at 180° W, $n(d, -180) = d + 24$ h. The intervals of the samples cover $(t_0 - \Delta, t_N]$, so the rain windows $(n - 24\,\mathrm{h}, n]$ are covered for all $n \in [d, d + 24\,\mathrm{h}]$ if and only if $d - 24\,\mathrm{h} \ge t_0 - \Delta$ and $d + 24\,\mathrm{h} \le t_N$; then every noon lies in $[t_0, t_N]$, and Proposition 1 gives its sample. For the local day $[s, s + 24\,\mathrm{h})$: if $s \le t_0 - \Delta$, the instant $t_0 - \Delta$ lies in it and is not a sample; if $s > t_0 - \Delta$, the first instant of the sequence at or after $s$ is at least $t_0$. Likewise the last instant before $s + 24$ h is at most $t_N$ if and only if $s + 24\,\mathrm{h} \le t_N + \Delta$. Over all longitudes $s$ runs through $[d - 12\,\mathrm{h}, d + 12\,\mathrm{h}]$. ∎

For GFS every 6 h from 20 October 18 UTC to 26 October 12 UTC, both rules give the days 22–25 October (checked).

## 4. Amounts from Rates

A forecast of precipitation gives a rate $r$ (kg m⁻² s⁻¹; 1 kg m⁻² of water is 1 mm). GFS's rate at a forecast hour is the mean over the hours of the forecast before it (3 or 6 h).

**Definition 3.** With the sample interval $I_k = (t_k - \Delta_k, t_k]$, $\Delta_k = t_k - t_{k-1}$ ($\Delta_0 = \Delta_1$), the amount in a window $(a, b]$ is

$$ A(a, b) = \sum_k r_k\, \bigl| I_k \cap (a, b] \bigr|. \tag{2} $$

**Proposition 3.** If the rate is constant on each $I_k$, (2) is the exact amount. $A$ is additive over adjacent windows, $A(a, c) = A(a, b) + A(b, c)$.

*Proof.* The integral of a step function is the sum of its values times the lengths of its steps within the window; additivity is that of the measure $|\cdot|$. ∎

*Caveat.* Sampling every 6 h from a product that alternates 3 h and 6 h means, a 3 h mean stands for 6 h: a shower in the other 3 h is missed, or one in these 3 h is counted twice. Rain from rates is a sum of samples, not the model's own accumulation.

## 5. Maxima from Samples

**Proposition 4.** The greatest sample within a local day is at most the day's greatest value. If the diurnal cycle is $T(s) = \bar{T} + a \cos(2\pi (s - s^*)/24\,\mathrm{h})$, the deficit is at most $a\,[1 - \cos(\pi \Delta / 24\,\mathrm{h})]$: $0.29\,a$ for $\Delta = 6$ h.

*Proof.* The first claim is that of any subset. Some sample lies within $\Delta/2$ of the peak $s^*$, where $T \ge \bar T + a \cos(\pi \Delta / 24\,\mathrm{h})$. ∎

With a diurnal range $2a$ of 10 °C, the sampled maximum is low by up to 1.5 °C; with 3-hourly samples, by up to 0.4 °C.

## 6. Rain Against the Normal

**Definition 4.** The *fraction of the month* of an amount $A$ over a window ending in month $m$ is $A / \bar{P}_m$, where $\bar{P}_m$ is the normal precipitation of that month (TerraClimate, 1991–2020).

The fraction compares the amount with an ordinary month. It is not a measure of rarity, which depends on how the month's rain is distributed among its days.

**Proposition 5.** Let daily amounts be independent, a day being wet with probability $q$ and a wet day's amount exponential with mean $\mu$, so that a month of 30 days has mean $\bar P = 30\, q \mu$ (Wilks 2011, §4.4). The probability that a day brings more than the fraction $\rho$ of the month is

$$ \Pr(A > \rho \bar{P}) = q\, e^{-30 \rho q}. \tag{3} $$

*Proof.* $\Pr(A > x) = q\, e^{-x/\mu}$, and $x/\mu = 30 \rho q$. ∎

For $\rho = 0.5$: where rain comes in rare heavy storms ($q = 0.1$), $\Pr = 0.022$, a day in 45; where it falls steadily ($q = 0.6$), $\Pr = 7.4 \times 10^{-5}$, a day in 37 years. The same fraction differs in rarity by a factor of 300. A return period needs the distribution of the heaviest days at the place, a climatology of extremes, which the library does not hold; so the plug-in reports the fraction and says what it does and does not mean. Floods depend further on the ground's wetness, the steepness of the land and the rivers that gather the rain, none of which the fraction sees.

Where $\bar P_m < 1$ mm (by default), the fraction is left out: it divides by almost nothing.

## 7. Verification

| Test | Agreement |
|---|---|
| Local noon at 90° E is 06:00 UTC | Exact |
| A rate of $10^{-3}$ kg m⁻² s⁻¹ for 24 h is 86.4 mm | Exact |
| A window over half of a 6 h interval takes half its rain | Exact |
| GFS every 6 h, 20 Oct 18 UTC to 26 Oct 12 UTC: whole days 22–25 October, by both rules | Exact |
| Heavy rain plug-in: 0.36 mm/h gives 8.64 mm in 24 h and 25.92 mm in 72 h; fractions of 100 mm | $10^{-6}$ |
| Heavy rain plug-in: no 72 h sum before 72 h are covered; no fraction where the normal is below 1 mm | Exact |
| Heavy rain plug-in on GFS at Ho Chi Minh City: the 72 h sum to 23 October 12 UTC (41.6 mm) is the sum of its three 24 h windows (23.9 + 11.6 + 6.0) | 0.1 mm |
