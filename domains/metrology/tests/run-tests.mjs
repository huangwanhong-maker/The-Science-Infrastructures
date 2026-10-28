// Checks of the metrology domain: time scales, precession and horizontal
// coordinates, against the IERS and Meeus (the provisional JavaScript); and its C++
// (time scales, leap seconds, Earth orientation, IAU 2006/2000A, the geocentric frames)
// against SOFA's test values and ERFA, built with CMake and run on calx's verification.
// Run from the repository root:  node domains/metrology/tests/run-tests.mjs

import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, cppChecks } from '../../../development/tools/checks.mjs';
import {
  scalesFromUnixMs, gmst, taiMinusUtc, deltaTBefore1972, DEG, fromSpherical, toSpherical, apply,
  precessionMatrix, horizontal,
} from '../js/index.js';

const { check, assert, skip, lines, done } = suite();
const hms = (h, m, s) => (h + m / 60 + s / 3600) * 15 * DEG;
const dms = (sign, d, m, s) => sign * (d + m / 60 + s / 3600) * DEG;

// Time scales.
const t2026 = scalesFromUnixMs(Date.UTC(2026, 9, 22, 12));
check('TT - UTC in 2026 is 69.184 s', t2026.deltaT - 69.184, 1e-9, 's');
check('TAI - UTC on 2016-12-31 is 36 s', taiMinusUtc(2457753.5) - 36, 0, 's');
check('Delta T in 1950 near 29.07 s', deltaTBefore1972(1950) - 29.07, 1e-9, 's');
{
  // Meeus, Astronomical Algorithms, example 12.a: 1987-04-10 0h UT, GMST 13h10m46.3668s.
  const t = scalesFromUnixMs(Date.UTC(1987, 3, 10));
  const diff = (gmst(t.jdUT1, t.jdTT) - hms(13, 10, 46.3668)) / (15 * DEG) * 3600;
  check('GMST 1987-04-10 0h UT (Meeus 12.a)', diff, 0.01, 's');
}

// Precession (Meeus example 21.b, theta Persei, mean place of date 2028-11-13.19 TD).
{
  const v0 = fromSpherical(hms(2, 44, 12.975), dms(1, 49, 13, 39.90));
  const v = apply(precessionMatrix(2462088.69), v0);
  const s = toSpherical(v);
  check('precessed RA (Meeus 21.b)', (s.lon - hms(2, 46, 11.331)) / (15 * DEG) * 3600, 0.01, 's');
  check('precessed Dec (Meeus 21.b)', (s.lat - dms(1, 49, 20, 54.54)) / DEG * 3600, 0.1, 'arcsec');
}

// Horizontal coordinates (Meeus example 13.b, Venus seen from the US Naval Observatory).
{
  const lst = hms(8, 34, 56.853) + dms(-1, 77, 3, 56);
  const h = horizontal(fromSpherical(hms(23, 9, 16.641), dms(-1, 6, 43, 11.61)), lst, dms(1, 38, 55, 17));
  check('altitude (Meeus 13.b)', h.alt / DEG - 15.1249, 0.001, 'deg');
  check('azimuth from north (Meeus 13.b)', h.az / DEG - (68.0337 + 180), 0.001, 'deg');
}

// The C++.
cppChecks({ assert, skip, lines }, { project: resolve(dirname(fileURLToPath(import.meta.url)), '..'), name: 'metrology' });

done();
