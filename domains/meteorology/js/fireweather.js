// Fire weather: the Canadian Forest Fire Weather Index (FWI) System, and the
// daily weather it is computed from, taken from a forecast.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).
//
// The system tracks the moisture of three layers of a standard forest floor by
// daily bookkeeping, and combines them with the wind into indices of fire
// behaviour (Van Wagner 1987):
//   FFMC  Fine Fuel Moisture Code: litter and fine fuels, time lag about 2/3 day;
//   DMC   Duff Moisture Code: loosely compacted organic layers, about 12 days;
//   DC    Drought Code: deep compact organic layers, about 52 days;
//   ISI   Initial Spread Index: FFMC and wind, the expected rate of spread;
//   BUI   Buildup Index: DMC and DC, the fuel available for combustion;
//   FWI   Fire Weather Index: ISI and BUI, the intensity of a spreading fire;
//   DSR   Daily Severity Rating, from the FWI, for averaging over time.
// Each code is a function of the day's weather and the code of the day before,
// so the codes must be started somewhere: in spring, after the snow, from
// FFMC 85, DMC 6, DC 15 (the standard start, STANDARD_START).
//
// Inputs, by the system's definition: the weather at local noon (standard
// time) of the day, in the open: temperature T (°C), relative humidity H (%),
// wind speed at 10 m W (km/h), and the rain of the 24 hours before (mm).
// Equation numbers are Van Wagner's (1987).
//
// The day-length factors of the DMC and the DC are those of 46° N (Van Wagner
// 1987, Tables 1 and 2); for other latitude bands, those of Lawson and
// Armitage (2008, Tables 3 and 4), as the cffdrs package carries them.
//
// Range of validity: the system was calibrated on jack and lodgepole pine
// stands in Canada; elsewhere its codes are indices of fuel dryness, not
// moisture contents of the local fuels. The arithmetic reproduces the
// published test data to their rounding (0.1).
//
// References:
//   Van Wagner, C. E. (1987), Development and Structure of the Canadian Forest
//   Fire Weather Index System, Canadian Forestry Service, Forestry Technical
//   Report 35, Ottawa.
//   Van Wagner, C. E. and Pickett, T. L. (1985), Equations and FORTRAN Program
//   for the Canadian Forest Fire Weather Index System, Forestry Technical
//   Report 33, Ottawa.
//   Lawson, B. D. and Armitage, O. B. (2008), Weather Guide for the Canadian
//   Forest Fire Danger Rating System, Natural Resources Canada, Edmonton.
//   Wang, X., Wotton, B. M., Cantin, A. S., et al. (2017), cffdrs: an R package
//   for the Canadian Forest Fire Danger Rating System, Ecol. Process. 6, 5.

export const STANDARD_START = { ffmc: 85, dmc: 6, dc: 15 };

// The effective day length Le (hours) of the DMC, by month (January first):
// at 46° N (Van Wagner 1987, Table 1) and in other bands (Lawson and Armitage 2008).
const LE_46N = [6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0];
const LE_20N = [7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 8.6, 8.1, 7.8];
const LE_20S = [10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2];
const LE_40S = [11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10.0, 11.2, 11.8];
// The day-length adjustment Lf of the DC, by month: north of 20° N (Van Wagner
// 1987, Table 2), south of 20° S (the same, six months on), and 1.4 between.
const LF_NORTH = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6];
const LF_SOUTH = [6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8];

// Le for a month (1-12) and a latitude (degrees): 46° N's north of 30° N,
// the bands' to 10° N, 10° S and 30° S, and 9 h within 10° of the equator.
export function dmcDayLength(month, lat) {
  const m = month - 1;
  if (lat > 30) return LE_46N[m];
  if (lat > 10) return LE_20N[m];
  if (lat > -10) return 9;
  if (lat > -30) return LE_20S[m];
  return LE_40S[m];
}

// Lf for a month (1-12) and a latitude (degrees).
export function dcDayLength(month, lat) {
  if (lat > 20) return LF_NORTH[month - 1];
  if (lat <= -20) return LF_SOUTH[month - 1];
  return 1.4;
}

// ---- the moisture codes ----------------------------------------------------------------

// FFMC from yesterday's F0, noon temperature T (°C), humidity H (%), wind W (km/h)
// and 24 h rain ro (mm). Equations 1-10.
export function fineFuelMoistureCode(F0, T, H, W, ro) {
  let mo = (147.2 * (101 - F0)) / (59.5 + F0);                                        // (1)
  if (ro > 0.5) {
    const rf = ro - 0.5;                                                              // (2)
    let mr = mo + 42.5 * rf * Math.exp(-100 / (251 - mo)) * (1 - Math.exp(-6.93 / rf)); // (3a)
    if (mo > 150) mr += 0.0015 * (mo - 150) ** 2 * Math.sqrt(rf);                    // (3b)
    mo = Math.min(mr, 250);
  }
  const Ed = 0.942 * H ** 0.679 + 11 * Math.exp((H - 100) / 10) + 0.18 * (21.1 - T) * (1 - Math.exp(-0.115 * H)); // (4)
  let m = mo;
  if (mo > Ed) {
    const ko = 0.424 * (1 - (H / 100) ** 1.7) + 0.0694 * Math.sqrt(W) * (1 - (H / 100) ** 8);           // (6a)
    const kd = ko * 0.581 * Math.exp(0.0365 * T);                                                       // (6b)
    m = Ed + (mo - Ed) * 10 ** -kd;                                                                     // (8)
  } else {
    const Ew = 0.618 * H ** 0.753 + 10 * Math.exp((H - 100) / 10) + 0.18 * (21.1 - T) * (1 - Math.exp(-0.115 * H)); // (5)
    if (mo < Ew) {
      const kl = 0.424 * (1 - ((100 - H) / 100) ** 1.7) + 0.0694 * Math.sqrt(W) * (1 - ((100 - H) / 100) ** 8); // (7a)
      const kw = kl * 0.581 * Math.exp(0.0365 * T);                                                         // (7b)
      m = Ew - (Ew - mo) * 10 ** -kw;                                                                       // (9)
    }
  }
  return (59.5 * (250 - m)) / (147.2 + m);                                            // (10)
}

// DMC from yesterday's P0, T (°C), H (%), 24 h rain ro (mm), the month (1-12)
// and the latitude (degrees). Equations 11-17; T below -1.1 °C is taken as -1.1.
export function duffMoistureCode(P0, T, H, ro, month, lat = 46) {
  let Po = P0;
  if (ro > 1.5) {
    const re = 0.92 * ro - 1.27;                                                      // (11)
    const Mo = 20 + Math.exp(5.6348 - Po / 43.43);                                    // (12)
    const b = Po <= 33 ? 100 / (0.5 + 0.3 * Po)                                       // (13a)
      : Po <= 65 ? 14 - 1.3 * Math.log(Po)                                            // (13b)
        : 6.2 * Math.log(Po) - 17.2;                                                  // (13c)
    const Mr = Mo + (1000 * re) / (48.77 + b * re);                                   // (14)
    Po = Math.max(0, 244.72 - 43.43 * Math.log(Mr - 20));                             // (15)
  }
  const K = 1.894 * (Math.max(T, -1.1) + 1.1) * (100 - H) * dmcDayLength(month, lat) * 1e-6;   // (16)
  return Po + 100 * K;                                                                // (17)
}

// DC from yesterday's D0, T (°C), 24 h rain ro (mm), the month and the latitude.
// Equations 18-23; T below -2.8 °C is taken as -2.8, and V is not negative.
export function droughtCode(D0, T, ro, month, lat = 46) {
  let Do = D0;
  if (ro > 2.8) {
    const rd = 0.83 * ro - 1.27;                                                      // (18)
    const Qo = 800 * Math.exp(-Do / 400);                                             // (19)
    const Qr = Qo + 3.937 * rd;                                                       // (20)
    Do = Math.max(0, 400 * Math.log(800 / Qr));                                       // (21)
  }
  const V = Math.max(0, 0.36 * (Math.max(T, -2.8) + 2.8) + dcDayLength(month, lat)); // (22)
  return Do + 0.5 * V;                                                                // (23)
}

// ---- the fire behaviour indices -----------------------------------------------------------

// ISI from the FFMC and the wind W (km/h). Equations 24-26, with m from (1).
export function initialSpreadIndex(F, W) {
  const m = (147.2 * (101 - F)) / (59.5 + F);
  const fW = Math.exp(0.05039 * W);                                                   // (24)
  const fF = 91.9 * Math.exp(-0.1386 * m) * (1 + m ** 5.31 / 4.93e7);                 // (25)
  return 0.208 * fW * fF;                                                             // (26)
}

// BUI from the DMC (P) and the DC (D). Equations 27a and 27b.
export function buildupIndex(P, D) {
  if (P === 0 && D === 0) return 0;
  const U = P <= 0.4 * D ? (0.8 * P * D) / (P + 0.4 * D)                              // (27a)
    : P - (1 - (0.8 * D) / (P + 0.4 * D)) * (0.92 + (0.0114 * P) ** 1.7);             // (27b)
  return Math.max(0, U);
}

// FWI from the ISI (R) and the BUI (U). Equations 28-30.
export function fireWeatherIndex(R, U) {
  const fD = U <= 80 ? 0.626 * U ** 0.809 + 2 : 1000 / (25 + 108.64 * Math.exp(-0.023 * U)); // (28a, b)
  const B = 0.1 * R * fD;                                                             // (29)
  return B > 1 ? Math.exp(2.72 * (0.434 * Math.log(B)) ** 0.647) : B;                 // (30a, b)
}

// DSR from the FWI. Equation 31.
export const dailySeverityRating = (S) => 0.0272 * S ** 1.77;

// One day of the system: yesterday's codes { ffmc, dmc, dc } and today's noon
// weather { T (°C), H (%), W (km/h), rain (mm, 24 h) } to today's codes and indices.
export function fireWeatherDay(yesterday, { T, H, W, rain }, month, lat = 46) {
  const ffmc = fineFuelMoistureCode(yesterday.ffmc, T, H, W, rain);
  const dmc = duffMoistureCode(yesterday.dmc, T, H, rain, month, lat);
  const dc = droughtCode(yesterday.dc, T, rain, month, lat);
  const isi = initialSpreadIndex(ffmc, W);
  const bui = buildupIndex(dmc, dc);
  const fwi = fireWeatherIndex(isi, bui);
  return { ffmc, dmc, dc, isi, bui, fwi, dsr: dailySeverityRating(fwi) };
}

// ---- daily weather from a forecast ----------------------------------------------------------
//
// A forecast gives its fields at fixed instants (UTC), every few hours; the
// system wants the weather at each place's local noon, and the rain of the day
// before it. Local noon is taken as mean solar noon, 12 h - lon / 15 h UTC
// (within 16 minutes of apparent noon, by the equation of time, and within
// an hour or so of standard-time noon in most zones). Times are ms UTC.

const HOUR = 3600000;

// Mean solar noon (ms UTC) on the UTC calendar day starting at dayMs, at a longitude (degrees east).
export const localNoon = (dayMs, lon) => dayMs + 12 * HOUR - (lon / 15) * HOUR;

// The index of the time nearest to ms, or -1 if none lies within `within` ms.
export function nearestIndex(times, ms, within = Infinity) {
  let best = -1, gap = Infinity;
  for (let k = 0; k < times.length; k++) {
    const d = Math.abs(times[k] - ms);
    if (d < gap) { gap = d; best = k; }
  }
  return gap <= within ? best : -1;
}

// The amount fallen in the window (t0, t1], from a rate sampled at `times`:
// rate(k) is the rate at times[k] (kg m-2 s-1; 1 kg m-2 is 1 mm of water), taken
// as the mean over the interval since the time before (the first time's
// interval as long as the second's), as GFS's precipitation rates are means over
// the preceding hours. Returns { amount (kg m-2), covered (s of the window
// within the samples' intervals) }. Exact for a rate constant over each interval.
export function accumulated(times, rate, t0, t1) {
  let amount = 0, covered = 0;
  for (let k = 0; k < times.length; k++) {
    const span = k > 0 ? times[k] - times[k - 1] : (times[1] - times[0]) || 0;
    const a = Math.max(times[k] - span, t0), b = Math.min(times[k], t1);
    if (b <= a) continue;
    const r = rate(k);
    if (!(r === r)) return { amount: NaN, covered: 0 };
    amount += (Math.max(0, r) * (b - a)) / 1000;
    covered += (b - a) / 1000;
  }
  return { amount, covered };
}

// The UTC days (ms at 00:00) on which a forecast holds, at every longitude, a
// time within half its spacing of local noon, and the rain of the 24 hours
// before it: the days on which the system can be computed everywhere.
export function forecastNoonDays(times) {
  if (times.length < 2) return [];
  const spacing = times[1] - times[0], last = times[times.length - 1], out = [];
  for (let d = Math.floor(times[0] / (24 * HOUR)) * 24 * HOUR; localNoon(d, 180) <= last; d += 24 * HOUR) {
    // Noon is earliest at 180° E (00:00 UTC) and latest at 180° W (24:00 UTC);
    // the first sample's rain is that of the interval before it (accumulated()).
    if (localNoon(d, 180) - 24 * HOUR >= times[0] - spacing && localNoon(d, -180) <= last) out.push(d);
  }
  return out;
}
