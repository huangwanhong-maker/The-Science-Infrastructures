---
title: The Canadian Forest Fire Weather Index System
status: draft
implemented-by:
  - domains/meteorology/js/fireweather.js
  - applications/earth-atlas/plugins/fire-weather.js
references: [VanWagner1987, VanWagnerPickett1985, LawsonArmitage2008, Wang2017cffdrs]
---

# The Canadian Forest Fire Weather Index System

## 1. Fuel Moisture as Bookkeeping

Whether a forest burns, and how fiercely, depends on how dry its fuels are, and the dryness of a fuel is the memory of the weather it has seen. The Canadian Forest Fire Weather Index (FWI) System (Van Wagner 1987) keeps that memory by daily bookkeeping. It follows three layers of a standard forest floor, each with its own time lag:

| Code | Layer | Nominal depth | Time lag |
|---|---|---|---|
| FFMC, $F$ | litter and fine fuels | 1.2 cm | 2/3 day |
| DMC, $P$ | loosely compacted organic layers (duff) | 7 cm | 12 days |
| DC, $D$ | deep, compact organic layers | 18 cm | 52 days |

Formally the system is a discrete dynamical system. With the codes $c_n = (F_n, P_n, D_n)$ on day $n$ and the day's weather $w_n = (T, H, W, r_o)$ (noon temperature in °C, relative humidity in %, wind speed at 10 m in km/h, and the rain of the 24 hours to noon in mm),

$$ c_{n+1} = \Phi(c_n, w_n; \text{month}, \varphi), \tag{0} $$

where $\Phi$ acts on each code separately: the three layers do not exchange water. Three indices are then functions of the state alone: the Initial Spread Index $R(F, W)$, the Buildup Index $U(P, D)$, and the Fire Weather Index $S(R, U)$. Each code is a monotone transform of a moisture content, chosen so that larger means drier.

Equation numbers below are Van Wagner's (1987), so that the code (`fireweather.js`) and the report can be read side by side.

## 2. The Fine Fuel Moisture Code

**Definition 1.** The FFMC $F \in [0, 101]$ and the moisture content $m$ (% of dry weight) of fine fuel are related by

$$ m = \frac{147.2\,(101 - F)}{59.5 + F}, \tag{1} \qquad F = \frac{59.5\,(250 - m)}{147.2 + m}. \tag{10} $$

*Remark.* (1) and (10) are not exact inverses. Solving (10) for $m$ gives $m = 147.2\,(101.05 - F)/(59.5 + F)$, since $59.5 \cdot 250 / 147.2 = 101.05$; (1) rounds this to 101. A round trip $F \mapsto m \mapsto F$ therefore moves $F$ by up to 0.05, which is the published rounding (§8).

**Rain.** If $r_o > 0.5$ mm, the first 0.5 mm is intercepted by the canopy, $r_f = r_o - 0.5$ (2), and

$$ m_r = m_o + 42.5\, r_f\, e^{-100/(251 - m_o)} \left(1 - e^{-6.93/r_f}\right) \;\bigl[\, + 0.0015\,(m_o - 150)^2 r_f^{1/2} \text{ if } m_o > 150 \bigr], \tag{3a, 3b} $$

with $m_r \le 250$.

**Drying and wetting.** The equilibrium moisture contents for drying and for wetting are

$$ E_d = 0.942 H^{0.679} + 11\, e^{(H-100)/10} + 0.18\,(21.1 - T)(1 - e^{-0.115 H}), \tag{4} $$
$$ E_w = 0.618 H^{0.753} + 10\, e^{(H-100)/10} + 0.18\,(21.1 - T)(1 - e^{-0.115 H}). \tag{5} $$

If $m_o > E_d$ the fuel dries, with log drying rate $k_d = k_o \cdot 0.581\, e^{0.0365 T}$, $k_o = 0.424\,[1 - (H/100)^{1.7}] + 0.0694\, W^{1/2}\,[1 - (H/100)^8]$ (6a, 6b); if $m_o < E_w$ it wets, with $k_w$ from (7a, 7b), the same with $H$ replaced by $100 - H$; between the two, $m = m_o$. Then

$$ m = E_d + (m_o - E_d)\, 10^{-k_d}, \tag{8} \qquad m = E_w - (E_w - m_o)\, 10^{-k_w}. \tag{9} $$

**Proposition 1 (time lag).** Equation (8) is the value after one day of the solution of $\mathrm{d}m/\mathrm{d}t = -\kappa\,(m - E_d)$ with $\kappa = k_d \ln 10$ per day. The time lag $1/\kappa$ is 0.62 day for $k_d = 0.7$, a typical summer value.

*Proof.* The solution is $m(t) = E_d + (m_o - E_d)\, e^{-\kappa t}$, and $e^{-\kappa} = 10^{-k_d}$. ∎

## 3. The Duff Moisture Code

**Definition 2.** The DMC $P \ge 0$ and the duff's moisture content $M$ (%) are related by

$$ M = 20 + e^{5.6348 - P/43.43}, \tag{12} \qquad P = 244.72 - 43.43 \ln(M - 20), \tag{15} $$

which are exact inverses, since $244.72 / 43.43 = 5.6348$.

**Rain.** If $r_o > 1.5$ mm, the effective rain is $r_e = 0.92\, r_o - 1.27$ (11), and $M_r = M_o + 1000\, r_e / (48.77 + b\, r_e)$ (14), with $b$ by (13a–c) according as $P_o \le 33$, $33 < P_o \le 65$ or $P_o > 65$; then $P_r$ by (15), not below 0.

**Drying.** With the effective day length $L_e$ (hours, §5) and $T \ge -1.1$ °C,

$$ K = 1.894\,(T + 1.1)(100 - H)\, L_e \times 10^{-6}, \tag{16} \qquad P = P_r + 100 K. \tag{17} $$

## 4. The Drought Code

**Definition 3.** The DC $D \ge 0$ is defined through a moisture equivalent $Q$ (units of 0.254 mm), $Q = 800\, e^{-D/400}$ (19), so that $D = 400 \ln(800 / Q)$ (21).

**Rain.** If $r_o > 2.8$ mm, $r_d = 0.83\, r_o - 1.27$ (18) and $Q_r = Q_o + 3.937\, r_d$ (20): 3.937 is $1/0.254$, rain in units of 0.254 mm.

**Drying.** With the day-length adjustment $L_f$ (§5) and $T \ge -2.8$ °C, the potential evapotranspiration is $V = 0.36\,(T + 2.8) + L_f \ge 0$ (22), and $D = D_r + 0.5\, V$ (23).

**Proposition 2 (memory of the start).** Through days without effective rain ($r_o \le 1.5$ mm for the DMC, $r_o \le 2.8$ mm for the DC), two runs started from different DMC or DC keep the same difference. Rain reduces it: for the DC, a difference $\delta$ in $D_o$ becomes, after rain $r_d$,

$$ \delta' = 400 \ln \frac{Q_o + 3.937\, r_d}{Q_o\, e^{-\delta/400} + 3.937\, r_d} < \delta. $$

*Proof.* Without rain, (17) and (23) add to $P$ and $D$ an amount that does not depend on them. With rain, let $Q_o$ be the first run's moisture equivalent and $\delta > 0$; the second's is $Q_o e^{-\delta/400}$, and (20, 21) give $\delta' = 400\, g(3.937\, r_d)$ with $g(a) = \ln[(Q_o + a)/(Q_o e^{-\delta/400} + a)]$, so long as neither $D_r$ is clamped at 0. Now $g(0) = \delta/400$ and $g'(a) = 1/(Q_o + a) - 1/(Q_o e^{-\delta/400} + a) < 0$, so $\delta' < \delta$ for $r_d > 0$. The DMC is alike, through (12), (14) and (15). ∎

*Remark.* The FFMC forgets its start geometrically (Proposition 1); the DMC and DC forget it only through rain. In the checks, 20 mm of rain halves a DMC difference of 20 and reduces a DC difference of 200 to 178. Hence the spin-up caveat of §7.

## 5. Day Length by Latitude

The factors $L_e$ and $L_f$ express the seasonal length of the drying day. Van Wagner's tables are for 46° N. Lawson and Armitage (2008) give factors for other bands, as the cffdrs package carries them (Wang et al. 2017):

| Band | $L_e$ (DMC) | $L_f$ (DC) |
|---|---|---|
| north of 30° N | 46° N table: 6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0 | north of 20° N: −1.6, −1.6, −1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, −1.6, −1.6 |
| 10° N to 30° N | 7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8 | 20° S to 20° N: 1.4 every month |
| 10° S to 10° N | 9 every month | |
| 30° S to 10° S | 10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2 | south of 20° S: the northern values six months on |
| south of 30° S | 11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10.0, 11.2, 11.8 | |

## 6. The Indices

**Initial Spread Index.** With $m$ from (1),

$$ f(W) = e^{0.05039 W}, \tag{24} \quad f(F) = 91.9\, e^{-0.1386 m} \left[1 + \frac{m^{5.31}}{4.93 \times 10^7}\right], \tag{25} \quad R = 0.208\, f(W)\, f(F). \tag{26} $$

**Buildup Index.**

$$ U = \begin{cases} \dfrac{0.8\, P D}{P + 0.4 D}, & P \le 0.4 D, \\[1ex] P - \left[1 - \dfrac{0.8 D}{P + 0.4 D}\right]\left[0.92 + (0.0114 P)^{1.7}\right], & P > 0.4 D. \end{cases} \tag{27a, 27b} $$

**Proposition 3.** $U$ is continuous across $P = 0.4 D$, where $U = P$.

*Proof.* At $P = 0.4 D$, (27a) gives $0.8 \cdot 0.4 D \cdot D / (0.8 D) = 0.4 D = P$, and the first bracket of (27b) is $1 - 0.8D/0.8D = 0$. ∎

**Fire Weather Index.**

$$ f(D) = \begin{cases} 0.626\, U^{0.809} + 2, & U \le 80, \\ 1000 / (25 + 108.64\, e^{-0.023 U}), & U > 80, \end{cases} \tag{28a, 28b} \qquad B = 0.1\, R\, f(D), \tag{29} $$

$$ \ln S = 2.72\, (0.434 \ln B)^{0.647} \;\; (B > 1), \qquad S = B \;\; (B \le 1). \tag{30a, 30b} $$

The Daily Severity Rating is $0.0272\, S^{1.77}$ (31). $S$ is continuous at $B = 1$ (both branches give 1); $f(D)$ jumps by 0.02 at $U = 80$ (23.69 below, 23.67 above), an artefact of the fitted forms.

## 7. Days of a Forecast

The system wants noon weather and the rain of the day before; a forecast gives fields at fixed instants. The plug-in takes, at each longitude, the forecast time nearest to local mean solar noon, and the rain of the 24 hours before it from the precipitation rate ([Daily weather from a forecast](daily-weather-from-a-forecast.md)). A day is computed only where every longitude has a time within half the spacing of its noon and 24 hours of rain before it; for GFS every 6 h from 20 October 18 UTC to 26 October 12 UTC, these are 22–25 October.

*Caveat (spin-up).* The codes must start from something. The standard start, FFMC 85, DMC 6, DC 15, is that of spring after the snow. By Proposition 2, a forecast of a few days started there describes those days, not the season's drought: in a dry autumn the DMC and DC, and so the BUI and FWI, are far too low. Starting codes from an agency's current analysis carry the season in; the plug-in accepts them as parameters.

## 8. Verification

| Test | Agreement |
|---|---|
| Van Wagner and Pickett (1985), standard test data: FFMC 85, DMC 6, DC 15 on 13 April at 46° N; five days, rain of 2.4 mm on the second; all six codes and indices | 0.048 at most (published to 0.1) |
| Day one: FFMC 87.7, DMC 8.5, DC 19.0, ISI 10.9, BUI 8.5, FWI 10.1 | 0.046 at most |
| (1) and (10) as inverses | 0.053 (§2, Remark) |
| Southern-hemisphere $L_f$ is the northern six months on; $L_e = 9$ h at the equator | Exact |
| Plug-in: the test day's weather at local noon only, cold and still otherwise: day one at six longitudes and three latitudes | 0.046 at most |
| Plug-in: 2.4 mm a day from a constant rate, against the domain's function | $6 \times 10^{-6}$ |
| Plug-in on GFS, 22–25 October 2026, 1° grid, in the page | 0.2 s |

*Caveat.* The system was calibrated for jack and lodgepole pine stands in Canada. Elsewhere its codes are indices of dryness, not the moisture of the local fuels, and its classes of danger must be set for the region. It says nothing over ice sheets, deserts without fuel, or the sea, where the plug-in leaves it out. The forecast's noon humidity and wind at 1° are not those of a forest clearing.
