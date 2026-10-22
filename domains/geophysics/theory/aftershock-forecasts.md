---
title: Aftershock Forecasts
status: draft
implemented-by:
  - domains/geophysics/js/aftershocks.js
  - applications/earth-atlas/plugins/aftershocks.js
references: [ReasenbergJones1989, WellsCoppersmith1994, Page2016, Utsu1995]
---

# Aftershock Forecasts

## 1. Aftershocks as a Point Process

A large earthquake changes the stress on the faults around it, and for days to years afterwards earthquakes there are many times more frequent. Two empirical laws describe them.

- **Gutenberg–Richter.** The number of earthquakes of magnitude $\ge M$ falls by a factor $10^{b}$ per unit of magnitude, with $b$ near 1.
- **The modified Omori law** (Utsu 1995). The rate of aftershocks decays with time $t$ since the mainshock as $(t + c)^{-p}$, with $p$ near 1 and a small delay $c$.

Reasenberg and Jones (1989) joined them into a nonhomogeneous Poisson process in time and magnitude, and so into probabilities.

## 2. The Rate and the Expected Number

**Definition 1.** After a mainshock of magnitude $M_m$, the rate (per day) of aftershocks of magnitude $\ge M$ at time $t$ (days) is

$$ \lambda(t, M) = 10^{\,a + b\,(M_m - M)}\,(t + c)^{-p}. \tag{1} $$

Reasenberg and Jones's generic parameters, fitted to 62 Californian sequences of mainshocks of magnitude 5 or more from 1933 to 1987, are $a = -1.67$, $b = 0.91$, $p = 1.08$, $c = 0.05$ day.

**Proposition 1.** The expected number of aftershocks of magnitude $\ge M$ in $[t_1, t_2]$ is, with $K = 10^{a + b(M_m - M)}$,

$$ N = \int_{t_1}^{t_2} \lambda\,\mathrm{d}t = \begin{cases} \dfrac{K}{p - 1}\left[(t_1 + c)^{1-p} - (t_2 + c)^{1-p}\right], & p \ne 1, \\[1.5ex] K \ln\dfrac{t_2 + c}{t_1 + c}, & p = 1, \end{cases} \tag{2} $$

and the first is continuous in $p$ at $p = 1$.

*Proof.* The antiderivative of $(t + c)^{-p}$ is $(t + c)^{1-p}/(1 - p)$ for $p \ne 1$ and $\ln(t + c)$ for $p = 1$. As $p \to 1$, $[x^{1-p} - y^{1-p}]/(p-1) = [e^{(1-p)\ln x} - e^{(1-p)\ln y}]/(p - 1) \to \ln y - \ln x$. ∎

*Remark.* For $p > 1$ the total over all time is finite, $K c^{1-p}/(p - 1)$. With the generic parameters the expected number of aftershocks as large as the mainshock, ever, is 0.34, so the chance of a later, larger earthquake is about 29%; in the first week it is 10%.

## 3. The Probability

**Proposition 2.** The probability of at least one aftershock of magnitude $\ge M$ in $[t_1, t_2]$ is

$$ P = 1 - e^{-N}. \tag{3} $$

For a window of fixed length $\tau$, $P$ decreases as its start $t_1$ moves away from the mainshock.

*Proof.* (3) is the probability of a nonzero count of a Poisson variable of mean $N$. $\partial N / \partial t_1 = \lambda(t_1 + \tau, M) - \lambda(t_1, M) < 0$, since $\lambda$ decreases in $t$, and $P$ increases with $N$. ∎

*Example (by hand).* After a mainshock of $M_m = 7$, the expected number of aftershocks of $M \ge 5$ in the first day is $10^{0.15}\,[0.05^{-0.08} - 1.05^{-0.08}]/0.08 = 1.412538 \times 0.274711 / 0.08 = 4.8505$, and $P = 1 - e^{-4.8505} = 0.992$.

## 4. The Zone of Aftershocks

Aftershocks fall mostly within about one rupture length of the rupture. Wells and Coppersmith (1994, Table 2A) relate the subsurface rupture length to moment magnitude, for all kinds of slip:

$$ \log_{10} L\,(\mathrm{km}) = -2.44 + 0.59\, M, \qquad \sigma(\log_{10} L) = 0.16, \qquad 4.8 \le M \le 8.1. \tag{4} $$

An earthquake of $M$ 7 breaks about 49 km, one of $M$ 7.7 about 127 km. The epicentre, where the rupture began, may lie anywhere along it; the plug-in draws the zone as a circle of radius $L$ about the epicentre, which holds the rupture wherever it began.

## 5. Observed and Expected

The catalogue's earthquakes within the zone after the mainshock may be compared with the model's expectation for the same magnitude threshold and time. After the $M_w$ 7.7 earthquake off Panama on 9 October 2026, USGS ComCat holds 77 earthquakes of $M \ge 4.5$ within 127 km to 21 October; the generic Californian model expects 98.5. The catalogue holds every earthquake there, not only aftershocks, misses many in the first hours, when the waves of large aftershocks overlap, and mixes magnitude types.

## 6. Verification

| Test | Agreement |
|---|---|
| (2) against Simpson's rule on a logarithmic variable, four windows, $p = 1.08$, 1 and 0.9 | $5 \times 10^{-14}$ (relative) |
| One unit of magnitude more: $10^{-0.91}$ times as many | $10^{-17}$ |
| The first day after $M$ 7, $M \ge 5$: 4.8505 (by hand) | $10^{-5}$ |
| $P$ over the next day, $M \ge 5$ after $M$ 7, from 0, 1, 3, 7, 30, 100 days | 0.992, 0.600, 0.304, 0.147, 0.035, 0.010: decreasing |
| A day, a week, a month: increasing | Yes |
| Wells and Coppersmith: $M$ 7, 48.98 km | 0.002 km |
| Plug-in: the mainshock is the largest within 300 km of the place; without a place, the largest of the month | Yes |
| Plug-in: the zone's circle at the rupture length | $10^{-11}$ km |
| Plug-in on USGS ComCat, Panama, $M$ 7.7: observed 77, expected 98.5; 58 ms in the page | — |

*Caveat.*
- The generic parameters are California's. Productivity and decay differ by tectonic region (subduction zones, mid-ocean ridges, stable interiors) and from one sequence to the next by up to a factor of ten in rate (Page et al. 2016), whose regional parameters may be set in the plug-in.
- The model and (4) are fitted to moment magnitude; a mainshock given in $m_b$ or $M_L$ may differ from its $M_w$ by several tenths.
- For $M$ at or above the mainshock's, $P$ is the chance of a larger earthquake to come, which the model counts among the aftershocks.
- These are probabilities of a model, not warnings. Operational aftershock forecasts are the USGS's and national surveys'.
