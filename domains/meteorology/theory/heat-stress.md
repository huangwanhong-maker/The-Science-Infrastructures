---
title: Heat Stress: the Wet-Bulb Temperature and the Heat Index
status: draft
implemented-by:
  - domains/meteorology/js/heat.js
  - applications/earth-atlas/plugins/heat.js
references: [Stull2011, Rothfusz1990, Steadman1979, Bolton1980, SherwoodHuber2010, Vecellio2022, Raymond2020, NWSHeatIndex]
---

# Heat Stress: the Wet-Bulb Temperature and the Heat Index

## 1. Two Measures of Heat

A body at rest makes about 100 W of heat, and must lose it. In hot air it loses it mainly by evaporating sweat, and evaporation is limited by the humidity of the air. Two measures follow from this.

- The **wet-bulb temperature** $T_w$ is the temperature to which air is cooled by evaporating water into it. It bounds from below the temperature of any wet surface in that air, skin included. It is physics, with no model of a person.
- The **heat index** is the temperature of dry air that would feel as hot as the actual air to a model person: Steadman's (1979) *apparent temperature*, a model of a clothed adult walking in the shade, fitted by regression.

## 2. The Wet-Bulb Temperature

**Definition 1 (isobaric wet-bulb temperature).** At pressure $p$, air of temperature $T$ and mixing ratio $r$ has wet-bulb temperature $T_w$ such that

$$ c_p\,(T - T_w) = L_v(T_w)\,\bigl(r_s(T_w, p) - r\bigr), \tag{1} $$

where $r_s$ is the saturation mixing ratio, $L_v$ the latent heat of vaporisation and $c_p$ the specific heat of moist air: the heat drawn from the air evaporates the water that saturates it.

**Proposition 1.** $T_w \le T$, with equality if and only if the air is saturated.

*Proof.* The right side of (1) increases with $T_w$, and the left side decreases. At $T_w = T$ the left side is 0 and the right side is $L_v (r_s(T) - r) \ge 0$, with equality only at saturation; so the root, where the two sides meet, is at most $T$. ∎

**Stull's fit.** Stull (2011, eq. 1) fitted $T_w$ at $p = 101.325$ kPa, as a function of $T$ (°C) and $RH$ (%):

$$ T_w = T \arctan\!\bigl[0.151977\,(RH + 8.313659)^{1/2}\bigr] + \arctan(T + RH) - \arctan(RH - 1.676331) + 0.00391838\, RH^{3/2} \arctan(0.023101\, RH) - 4.686035. \tag{2} $$

It is valid for $5\% \le RH \le 99\%$ and $-20$ °C $\le T \le 50$ °C, except where it is both cold and dry, with errors between $-1$ °C and $+0.65$ °C and a mean absolute error below 0.3 °C. Stull's example: $T = 20$ °C, $RH = 50\%$ gives 13.7 °C.

*Remark.* Against our own solution of (1) at 101.325 kPa (Bolton's 1980 saturation pressure, $L_v = 2.501 \times 10^6 - 2370\, T$ J/kg), the fit differs by up to 1.17 °C, 0.37 °C on average, over its range. Stull's stated bounds are against his own reference, and the two references differ in their constants; ours gives 13.79 °C in his example.

The plug-in leaves $T_w$ out beyond the range, and takes $RH$ above 99% as 99%, which changes $T_w$ by less than 0.2 °C. On high ground the pressure is lower and (2) less accurate: at 850 hPa by about 1 °C in warm, humid air.

## 3. The Heat Index

The National Weather Service computes the heat index $HI$ (°F) from $T$ (°F) and $RH$ (%) in three steps (NWS Weather Prediction Center; Rothfusz 1990).

1. The simple formula, fitted to Steadman's values:

   $$ HI_s = 0.5\,\{T + 61.0 + 1.2\,(T - 68.0) + 0.094\, RH\}. \tag{3} $$

   If $\tfrac12 (HI_s + T) < 80$, then $HI = HI_s$.
2. Otherwise Rothfusz's regression:

   $$ \begin{aligned} HI ={}& -42.379 + 2.04901523\, T + 10.14333127\, RH - 0.22475541\, T\,RH - 0.00683783\, T^2 \\ & - 0.05481717\, RH^2 + 0.00122874\, T^2 RH + 0.00085282\, T\,RH^2 - 0.00000199\, T^2 RH^2, \end{aligned} \tag{4} $$

   which reproduces Steadman's table within ±1.3 °F.
3. Adjustments: where $RH < 13\%$ and $80 \le T \le 112$ °F, subtract $\tfrac{13 - RH}{4}\sqrt{\tfrac{17 - |T - 95|}{17}}$; where $RH > 85\%$ and $80 \le T \le 87$ °F, add $\tfrac{RH - 85}{10} \cdot \tfrac{87 - T}{5}$.

*Remark.* The rule is not continuous: along the curve where it changes from (3) to (4), the two differ by up to 2.3 °F (near 79 °F and 92%). In the ten cells of the NWS chart checked (§6), away from the adjustments, the chart is (4) rounded to whole degrees.

*Caveat.* The heat index describes a person in the shade, in light wind. Full sun can add up to 15 °F (8 °C). It is not defined for cold air, where (3) returns roughly $T$.

## 4. Thresholds

- **Wet-bulb temperature.** A sustained $T_w$ of 35 °C is the theoretical bound of human heat loss: skin at 35 °C can no longer shed metabolic heat (Sherwood and Huber 2010). Measurements in young, healthy people find their limits lower, near 31 °C in humid heat (Vecellio et al. 2022). Wet-bulb temperatures above 31 °C are rare, and above 35 °C have been recorded only for an hour or two at a few stations by the Persian Gulf and in South Asia (Raymond et al. 2020).
- **Heat index.** The NWS names four classes: caution, 80–90 °F (27–32 °C); extreme caution, 90–103 °F (32–39 °C); danger, 103–124 °F (39–51 °C); extreme danger, above 124 °F (51 °C).

## 5. The Excess over the Normal

The plug-in also compares the day's greatest forecast temperature with the normal daily maximum of the month (TerraClimate, 1991–2020): $\Delta T = \max_{\text{day}} T - \overline{T}_{\max, m}$. This is an anomaly, not a measure of stress; a day 5 °C above normal is a heatwave in Paris and comfortable in Delhi.

*Caveat.*
- The forecast's maximum is the greatest of its samples (every 6 h for GFS), which misses the afternoon peak by up to 1.5 °C ([Daily weather from a forecast](daily-weather-from-a-forecast.md), Proposition 4); the normal is a true daily maximum. The excess is therefore biased low by about as much.
- The normal is read at the forecast's 1° points from TerraClimate's finer grid; in mountains the two grids stand at different heights, and the excess carries the lapse rate of about 6.5 °C per km between them.

## 6. Verification

| Test | Agreement |
|---|---|
| Stull's example, 20 °C and 50% → 13.7 °C | 0.0007 °C |
| Stull against the isobaric solution of (1) over his range, worst and mean | 1.17 °C; 0.37 °C |
| $T_w \le T$ over the range (Proposition 1) | No exception |
| NWS chart, ten cells away from the adjustments: (80, 40) 80; (90, 40) 91; (100, 40) 109; (110, 40) 136; (90, 50) 95; (100, 50) 118; (96, 60) 116; (100, 65) 136; (88, 80) 106; (90, 90) 122 (°F, %) | 0.42 °F at most |
| The simple formula below 80 °F | Exact |
| Plug-in: 30 °C within 3 h of local noon, 20 °C otherwise, at 50%: daily maxima of $T_w$ and the heat index | $2 \times 10^{-5}$ °C |
| Plug-in: 5 °C above a normal of 25 °C, at five longitudes | $10^{-5}$ °C |
| Plug-in on GFS, 22–25 October 2026: Delhi's daily maxima 0.4–2.1 °C below the October normal, as the sampling predicts; 0.7 s in the page | — |
