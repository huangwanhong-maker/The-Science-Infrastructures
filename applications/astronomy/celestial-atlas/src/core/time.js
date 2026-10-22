// Time scales: UTC, UT1, TAI, TT and TDB, and the rotation angle of the Earth.
//
// Provisional stand-in for core/metrology/time (laboratory ADR 0011). When the
// C++ core exists, this module is replaced by it, and tests/run-tests.mjs is
// rerun against the same reference values.
//
// Sources:
//   IERS Bulletin C, leap seconds (TAI - UTC) from 1972.
//   IERS Conventions (2010), eq. 5.15: Earth rotation angle.
//   Capitaine, Wallace and Chapront (2005); SOFA iauGmst06: GMST (IAU 2006).
//   USNO Circular 179 (Kaplan 2005), eq. 2.6: TDB - TT.
//   Espenak and Meeus (2006), NASA/TP-2006-214141: Delta T before 1972.

export const J2000 = 2451545.0;           // JD of J2000.0 (TT)
export const DAY = 86400;                 // seconds per day
export const JULIAN_CENTURY = 36525;      // days
const JD_UNIX_EPOCH = 2440587.5;          // JD of 1970-01-01T00:00 UTC
const TT_MINUS_TAI = 32.184;              // seconds, exact

// [JD (UTC) at which the value takes effect, TAI - UTC in seconds]
const LEAP_SECONDS = [
  [2441317.5, 10], [2441499.5, 11], [2441683.5, 12], [2442048.5, 13],
  [2442413.5, 14], [2442778.5, 15], [2443144.5, 16], [2443509.5, 17],
  [2443874.5, 18], [2444239.5, 19], [2444786.5, 20], [2445151.5, 21],
  [2445516.5, 22], [2446247.5, 23], [2447161.5, 24], [2447892.5, 25],
  [2448257.5, 26], [2448804.5, 27], [2449169.5, 28], [2449534.5, 29],
  [2450083.5, 30], [2450630.5, 31], [2451179.5, 32], [2453736.5, 33],
  [2454832.5, 34], [2456109.5, 35], [2457204.5, 36], [2457754.5, 37],
];
const FIRST_LEAP = LEAP_SECONDS[0][0];

// TAI - UTC in seconds for a UTC Julian date on or after 1972-01-01.
// Dates after the last entry keep the last value: future leap seconds are unknown.
export function taiMinusUtc(jdUTC) {
  let value = LEAP_SECONDS[0][1];
  for (const [start, seconds] of LEAP_SECONDS) {
    if (jdUTC >= start) value = seconds;
    else break;
  }
  return value;
}

// Delta T = TT - UT1 in seconds, for years before 1972 (Espenak and Meeus 2006).
export function deltaTBefore1972(year) {
  if (year < 1860) {
    const t = year - 1800;
    return 13.72 - 0.332447 * t + 0.0068612 * t ** 2 + 0.0041116 * t ** 3
      - 0.00037436 * t ** 4 + 0.0000121272 * t ** 5 - 0.0000001699 * t ** 6
      + 0.000000000875 * t ** 7;
  }
  if (year < 1900) {
    const t = year - 1860;
    return 7.62 + 0.5737 * t - 0.251754 * t ** 2 + 0.01680668 * t ** 3
      - 0.0004473624 * t ** 4 + t ** 5 / 233174;
  }
  if (year < 1920) {
    const t = year - 1900;
    return -2.79 + 1.494119 * t - 0.0598939 * t ** 2 + 0.0061966 * t ** 3 - 0.000197 * t ** 4;
  }
  if (year < 1941) {
    const t = year - 1920;
    return 21.20 + 0.84493 * t - 0.076100 * t ** 2 + 0.0020936 * t ** 3;
  }
  if (year < 1961) {
    const t = year - 1950;
    return 29.07 + 0.407 * t - t ** 2 / 233 + t ** 3 / 2547;
  }
  const t = year - 1975;
  return 45.45 + 1.067 * t - t ** 2 / 260 - t ** 3 / 718;
}

// TDB - TT in seconds; amplitude 1.7 ms (USNO Circular 179, eq. 2.6).
export function tdbMinusTt(jdTT) {
  const g = (357.53 + 0.98560028 * (jdTT - J2000)) * Math.PI / 180;
  return 0.001657 * Math.sin(g) + 0.000014 * Math.sin(2 * g);
}

// All time scales for an instant given as Unix milliseconds of civil time.
//
// From 1972 the civil time is UTC: TT = UTC + (TAI - UTC) + 32.184 s, and
// UT1 is taken equal to UTC (|UT1 - UTC| < 0.9 s).
// Before 1972 the civil time is taken as UT1, and TT = UT1 + Delta T.
export function scalesFromUnixMs(unixMs) {
  const jdUTC = unixMs / 86400000 + JD_UNIX_EPOCH;
  let jdUT1, jdTT, deltaT;
  if (jdUTC >= FIRST_LEAP) {
    const ttMinusUtc = taiMinusUtc(jdUTC) + TT_MINUS_TAI;
    jdUT1 = jdUTC;
    jdTT = jdUTC + ttMinusUtc / DAY;
    deltaT = ttMinusUtc;
  } else {
    const year = 2000 + (jdUTC - J2000) / 365.25;
    deltaT = deltaTBefore1972(year);
    jdUT1 = jdUTC;
    jdTT = jdUT1 + deltaT / DAY;
  }
  const jdTDB = jdTT + tdbMinusTt(jdTT) / DAY;
  return {
    unixMs,
    jdUTC,
    jdUT1,
    jdTT,
    jdTDB,
    et: (jdTDB - J2000) * DAY,          // seconds past J2000 TDB, the SPK time argument
    deltaT,                            // TT - UT1, seconds
  };
}

// Inverse of the main relation: Unix milliseconds for a TDB Julian date (within ~2 ms).
export function unixMsFromJdTDB(jdTDB) {
  let ms = (jdTDB - JD_UNIX_EPOCH) * 86400000;
  for (let i = 0; i < 3; i++) {
    const s = scalesFromUnixMs(ms);
    ms -= (s.jdTDB - jdTDB) * 86400000;
  }
  return ms;
}

const TWO_PI = 2 * Math.PI;
const ARCSEC = Math.PI / 648000;

function normaliseAngle(a) {
  const r = a % TWO_PI;
  return r < 0 ? r + TWO_PI : r;
}

// Earth rotation angle, radians (IERS Conventions 2010, eq. 5.15).
export function earthRotationAngle(jdUT1) {
  const du = jdUT1 - J2000;
  const frac = (jdUT1 % 1) + (J2000 % 1);   // keep precision in the fractional day
  return normaliseAngle(TWO_PI * (frac + 0.7790572732640 + 0.00273781191135448 * du));
}

// Greenwich mean sidereal time, radians (IAU 2006; SOFA iauGmst06).
export function gmst(jdUT1, jdTT) {
  const t = (jdTT - J2000) / JULIAN_CENTURY;
  const poly = 0.014506 + t * (4612.156534 + t * (1.3915817 + t * (-0.00000044
    + t * (-0.000029956 + t * -0.0000000368))));
  return normaliseAngle(earthRotationAngle(jdUT1) + poly * ARCSEC);
}

// Calendar text, UTC, to the second: "2026-10-22 13:04:05".
export function formatUtc(unixMs) {
  const d = new Date(unixMs);
  if (Number.isNaN(d.getTime())) return '—';
  const iso = d.toISOString();
  return `${iso.slice(0, 10)} ${iso.slice(11, 19)}`;
}
