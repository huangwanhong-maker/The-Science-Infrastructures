// Heat: the wet-bulb temperature, the heat index, and the daily maximum of a
// quantity from a forecast's samples.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).
//
// Wet-bulb temperature. Stull's (2011) empirical fit, eq. 1, of the
// psychrometric wet-bulb temperature at standard sea-level pressure
// (101.325 kPa) as a function of temperature T (°C) and relative humidity
// RH (%):
//   Tw = T atan[0.151977 (RH + 8.313659)^½] + atan(T + RH) - atan(RH - 1.676331)
//        + 0.00391838 RH^(3/2) atan(0.023101 RH) - 4.686035.
// Valid for 5 % ≤ RH ≤ 99 % and -20 °C ≤ T ≤ 50 °C, except where it is both
// cold and dry; within that range the error is between -1 °C and +0.65 °C, the
// mean absolute error less than 0.3 °C. At other pressures it is not valid
// (at 850 hPa it is off by about 1 °C in warm, humid air).
//
// Heat index. The National Weather Service's: Rothfusz's (1990) regression on
// Steadman's (1979) apparent temperature, in °F, with the NWS's two
// adjustments, and below 80 °F the simple formula of Steadman's values:
//   HI_simple = 0.5 {T + 61.0 + 1.2 (T - 68.0) + 0.094 RH};
// if the mean of HI_simple and T is 80 °F or more, the regression
//   HI = -42.379 + 2.04901523 T + 10.14333127 RH - 0.22475541 T RH - 0.00683783 T²
//        - 0.05481717 RH² + 0.00122874 T² RH + 0.00085282 T RH² - 0.00000199 T² RH²,
// less [(13 - RH)/4] √{[17 - |T - 95|]/17} where RH < 13 % and 80 ≤ T ≤ 112 °F,
// plus [(RH - 85)/10] [(87 - T)/5] where RH > 85 % and 80 ≤ T ≤ 87 °F.
// The regression reproduces Steadman's table within ±1.3 °F; it describes a
// person of given size and clothing in the shade, in light wind, and is not
// meant for temperatures below 80 °F (27 °C).
//
// References:
//   Stull, R. (2011), Wet-bulb temperature from relative humidity and air
//   temperature, J. Appl. Meteor. Climatol. 50, 2267-2269.
//   Rothfusz, L. P. (1990), The heat index "equation" (or, more than you ever
//   wanted to know about heat index), NWS Southern Region Technical
//   Attachment SR 90-23, Fort Worth.
//   Steadman, R. G. (1979), The assessment of sultriness. Part I, J. Appl.
//   Meteor. 18, 861-873.
//   NWS Weather Prediction Center, The Heat Index Equation,
//   https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml (retrieved 2026-10-23).

export const STULL_RANGE = { tMin: -20, tMax: 50, rhMin: 5, rhMax: 99 };

// Wet-bulb temperature (°C) from temperature (°C) and relative humidity (%), Stull (2011) eq. 1.
export function wetBulbStull(T, RH) {
  return T * Math.atan(0.151977 * Math.sqrt(RH + 8.313659)) + Math.atan(T + RH) - Math.atan(RH - 1.676331)
    + 0.00391838 * RH ** 1.5 * Math.atan(0.023101 * RH) - 4.686035;
}

// Whether (T °C, RH %) lies in the range Stull (2011) gives for eq. 1.
export const withinStullRange = (T, RH) => T >= STULL_RANGE.tMin && T <= STULL_RANGE.tMax && RH >= STULL_RANGE.rhMin && RH <= STULL_RANGE.rhMax;

// The NWS heat index (°F) from temperature (°F) and relative humidity (%).
export function heatIndexF(T, RH) {
  const simple = 0.5 * (T + 61.0 + (T - 68.0) * 1.2 + RH * 0.094);
  if ((simple + T) / 2 < 80) return simple;
  let hi = -42.379 + 2.04901523 * T + 10.14333127 * RH - 0.22475541 * T * RH - 0.00683783 * T * T
    - 0.05481717 * RH * RH + 0.00122874 * T * T * RH + 0.00085282 * T * RH * RH - 0.00000199 * T * T * RH * RH;
  if (RH < 13 && T >= 80 && T <= 112) hi -= ((13 - RH) / 4) * Math.sqrt((17 - Math.abs(T - 95)) / 17);
  else if (RH > 85 && T >= 80 && T <= 87) hi += ((RH - 85) / 10) * ((87 - T) / 5);
  return hi;
}

// The heat index in kelvin, from temperature in kelvin and relative humidity (%).
export const heatIndex = (T, RH) => ((heatIndexF((T - 273.15) * 1.8 + 32, RH) - 32) / 1.8) + 273.15;

// The greatest of a quantity's samples at `times` within [t0, t1): value(k) at
// times[k]. NaN if a sample within is missing, or if none lies within.
export function sampledMaximum(times, value, t0, t1) {
  let best = -Infinity;
  for (let k = 0; k < times.length; k++) {
    if (times[k] < t0 || times[k] >= t1) continue;
    const v = value(k);
    if (!(v === v)) return NaN;
    if (v > best) best = v;
  }
  return best === -Infinity ? NaN : best;
}

// The UTC days (ms at 00:00) whose local solar day, from 00:00 to 24:00 mean
// solar time, holds all its samples at every longitude: no sample of the
// forecast's spacing is missing from it at either end.
export function forecastLocalDays(times) {
  if (times.length < 2) return [];
  const H = 3600000, spacing = times[1] - times[0], last = times[times.length - 1], out = [];
  for (let d = Math.floor(times[0] / (24 * H)) * 24 * H; d - 12 * H <= last; d += 24 * H) {
    // The local day at 180° E begins at d - 12 h UTC; at 180° W it ends at d + 36 h.
    if (d - 12 * H > times[0] - spacing && d + 36 * H <= last + spacing) out.push(d);
  }
  return out;
}

// The local solar day at a longitude, as an interval of UTC: [start, end) in ms.
export function localDay(dayMs, lon) {
  const start = dayMs - (lon / 15) * 3600000;
  return [start, start + 86400000];
}
