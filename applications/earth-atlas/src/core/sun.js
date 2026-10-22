// The Sun seen from the Earth: its apparent position, the subsolar point, the
// Sun's altitude anywhere, and the day, twilight and night that follow.
//
// Provisional stand-in for core/domains/astronomy (laboratory ADR 0014). The
// Celestial Atlas computes the Sun from DE440; here a low-precision theory is
// enough, because the terminator is drawn on a globe, not measured.
//
// Accuracy: 0.01 degree in right ascension and declination from 1950 to 2050
// (Astronomical Almanac); the subsolar point agrees with JPL Horizons to about
// 0.01 degree (about 1 km), tested in tests/run-tests.mjs. UT1 is taken equal to
// UTC (at most 0.9 s, or 0.004 degree of longitude), and GMST stands for GAST
// (the equation of the equinoxes is at most 1.2 s).
//
// Sources:
//   Astronomical Almanac (2019), section C5, "Low precision formulas for the Sun".
//   Meeus (1998), Astronomical Algorithms, eq. 12.4: mean sidereal time at Greenwich.
//   Twilight: civil, nautical and astronomical at solar altitudes of -6, -12 and
//   -18 degrees; sunrise and sunset at -0.833 degree (refraction 34' and the
//   solar semidiameter 16'), Astronomical Almanac, section A.

import { DEG, toVector, wrapLongitude } from './sphere.js';

export const J2000 = 2451545.0;
export const SUNRISE_ALTITUDE = -0.833;
export const TWILIGHT = { civil: -6, nautical: -12, astronomical: -18 };

export const julianDate = (ms) => ms / 86400000 + 2440587.5;

// Greenwich mean sidereal time, degrees, for a UT1 Julian date.
export function gmst(jd) {
  const t = (jd - J2000) / 36525;
  const theta = 280.46061837 + 360.98564736629 * (jd - J2000) + 0.000387933 * t * t - (t * t * t) / 38710000;
  return ((theta % 360) + 360) % 360;
}

// Apparent right ascension and declination (degrees, equinox of date), the
// ecliptic longitude and the distance (au) of the Sun at a Julian date.
export function sunPosition(jd) {
  const n = jd - J2000;
  const L = 280.460 + 0.9856474 * n;
  const g = (357.528 + 0.9856003 * n) * DEG;
  const lambda = (L + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g)) * DEG;
  const epsilon = (23.439 - 0.0000004 * n) * DEG;
  const ra = Math.atan2(Math.cos(epsilon) * Math.sin(lambda), Math.cos(lambda)) / DEG;
  const dec = Math.asin(Math.sin(epsilon) * Math.sin(lambda)) / DEG;
  return {
    ra: ((ra % 360) + 360) % 360,
    dec,
    longitude: ((lambda / DEG) % 360 + 360) % 360,
    distanceAu: 1.00014 - 0.01671 * Math.cos(g) - 0.00014 * Math.cos(2 * g),
  };
}

// The place where the Sun is at the zenith at a UTC instant (ms since 1970).
export function subsolarPoint(ms) {
  const jd = julianDate(ms);
  const { ra, dec } = sunPosition(jd);
  return { lon: wrapLongitude(ra - gmst(jd)), lat: dec };
}

// The Sun's geocentric altitude and azimuth (degrees; azimuth from north through
// east) at a place, given the subsolar point. Parallax (8.8") is neglected.
export function sunAltAz(lon, lat, sub) {
  const p = lat * DEG, d = sub.lat * DEG, h = (lon - sub.lon) * DEG;     // h: west hour angle
  const sinAlt = Math.sin(p) * Math.sin(d) + Math.cos(p) * Math.cos(d) * Math.cos(h);
  const alt = Math.asin(Math.max(-1, Math.min(1, sinAlt)));
  const az = Math.atan2(-Math.cos(d) * Math.sin(h), Math.sin(d) * Math.cos(p) - Math.cos(d) * Math.sin(p) * Math.cos(h));
  return { altitude: alt / DEG, azimuth: ((az / DEG) % 360 + 360) % 360 };
}

export function daylight(altitude) {
  if (altitude >= SUNRISE_ALTITUDE) return 'day';
  if (altitude >= TWILIGHT.civil) return 'civil twilight';
  if (altitude >= TWILIGHT.nautical) return 'nautical twilight';
  if (altitude >= TWILIGHT.astronomical) return 'astronomical twilight';
  return 'night';
}

// Unit vector toward the Sun in Earth-fixed coordinates.
export const sunVector = (sub) => toVector(sub.lon, sub.lat);

// Local mean solar time at a longitude, hours in [0, 24).
export const localMeanSolarTime = (ms, lon) => ((((ms / 3600000 + lon / 15) % 24) + 24) % 24);
