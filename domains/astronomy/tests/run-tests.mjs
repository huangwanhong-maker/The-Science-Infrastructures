// Checks of the astronomy domain: the DE440 ephemeris, apparent directions,
// orbits and the rotation of the bodies against JPL Horizons; the Sun seen from
// the Earth; and the relations between galaxies against Tully (2015) and
// Karachentsev et al. (2013).
// Run from the repository root:  node domains/astronomy/tests/run-tests.mjs
// The ephemeris and the galaxy catalogues come through the data library
// (python data/acquisition/fetch.py jpl-de440s | tully-2mass-groups | ungc);
// the checks that need them are skipped, and say so, when they are missing.
// The domain's C++ (the DE440 ephemeris and the Moon's DE440 frames, against Horizons and NAIF's
// worked example; the IAU rotation models) is built with CMake and run on calx's verification.

import { existsSync, readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, CACHE, cppChecks } from '../../../development/tools/checks.mjs';
import {
  SPK, GM_SUN, osculating, bodyFixedMatrix, supergalactic, tidalIndices, tidalGroups,
  inferDisturberLuminosities, displayName, subsolarPoint, sunAltAz, gmst, julianDate, daylight,
  localMeanSolarTime,
} from '../js/index.js';
import {
  DEG, AU_KM, LIGHT_KM_S, fromSpherical, apply, sub, norm, dot, cross, scalesFromUnixMs, ICRS_TO_SUPERGALACTIC,
} from '../../metrology/js/index.js';
import { wrapLongitude, centralAngle, eastNorth, fromVector } from '../../geodesy/js/index.js';

const { check, assert, skip, lines, done } = suite();
const here = dirname(fileURLToPath(import.meta.url));
const fixture = (name) => JSON.parse(readFileSync(resolve(here, 'fixtures', name), 'utf8'));

const etOf = (jdTdb) => (jdTdb - 2451545.0) * 86400;
const unixMsOf = (jd) => (jd - 2440587.5) * 86400000;
const separation = (a, b) => Math.atan2(norm(cross(a, b)), dot(a, b));

// ---- the ephemeris -----------------------------------------------------------
const KERNEL = resolve(CACHE, 'jpl-de440s/de440s.bsp');
if (existsSync(KERNEL)) {
  const kernel = readFileSync(KERNEL);
  const spk = new SPK(kernel.buffer.slice(kernel.byteOffset, kernel.byteOffset + kernel.byteLength));
  const horizons = fixture('horizons.json');

  // DE440 evaluator against Horizons barycentric vectors (DE441).
  for (const [target, entries] of Object.entries(horizons.vectors_ssb_icrf_km)) {
    let dp = 0, dv = 0;
    for (const e of entries) {
      const s = spk.stateSSB(Number(target), etOf(e.jd_tdb));
      dp = Math.max(dp, norm(sub(s.p, e.p)));
      dv = Math.max(dv, norm(sub(s.v, e.v)));
    }
    check(`position of ${target} vs Horizons, max over 1850-2150`, dp, 2, 'km');
    check(`velocity of ${target} vs Horizons`, dv, 1e-6, 'km/s');
  }

  // Geocentric astrometric directions (light-time corrected) against Horizons.
  const naifOf = { 10: 10, 301: 301, 499: 4, 599: 5, 699: 6 };
  for (const [target, entries] of Object.entries(horizons.astrometric_geocentric_icrf)) {
    let worst = 0;
    for (const e of entries) {
      const t = scalesFromUnixMs(unixMsOf(e.jd_utc));
      const earth = spk.stateSSB(399, t.et).p;
      let tau = 0, g;
      for (let i = 0; i < 3; i++) {
        g = sub(spk.stateSSB(naifOf[target], t.et - tau).p, earth);
        tau = norm(g) / LIGHT_KM_S;
      }
      const ref = fromSpherical(e.ra_deg * DEG, e.dec_deg * DEG);
      worst = Math.max(worst, separation(g, ref) / DEG * 3600);
    }
    check(`astrometric direction of ${target} vs Horizons`, worst, 0.5, 'arcsec');
  }

  // Osculating orbit of the Earth-Moon barycentre at J2000.
  {
    const emb = spk.stateSSB(3, 0), sun = spk.stateSSB(10, 0);
    const o = osculating(sub(emb.p, sun.p), sub(emb.v, sun.v), GM_SUN + 403503.235502);
    check('EMB semi-major axis near 1.0000 AU', o.a / AU_KM - 1.0000, 0.0005, 'AU');
    check('EMB eccentricity near 0.0167', o.e - 0.0167, 0.0005, '');
  }

  // Orientation of the bodies (IAU rotation models, Archinal et al. 2018) against the
  // sub-Earth points of JPL Horizons, 2000-2035: longitude of the prime meridian and
  // direction of the pole. Horizons gives planetodetic latitudes; ours are converted.
  {
    const orientation = fixture('orientation.json').entries;
    const NAIF = { sun: 10, mercury: 199, venus: 299, moon: 301, mars: 4, jupiter: 5, saturn: 6, uranus: 7, neptune: 8, pluto: 9 };
    const FLATTENING = { mars: 0.00589, jupiter: 0.06487, saturn: 0.09796, uranus: 0.02293, neptune: 0.01708 };
    const worst = {};
    for (const e of orientation) {
      const et = etOf(e.jd_tt);
      const earth = spk.stateSSB(399, et).p;
      let tau = 0, g;
      for (let i = 0; i < 3; i++) { g = sub(earth, spk.stateSSB(NAIF[e.body], et - tau).p); tau = norm(g) / LIGHT_KM_S; }
      const v = apply(bodyFixedMatrix(e.body, e.jd_tt - tau / 86400), g);
      let lon = Math.atan2(v[1], v[0]) / DEG;
      if (e.west_positive) lon = -lon;
      const f = FLATTENING[e.body] || 0;
      const lat = Math.atan(Math.tan(Math.asin(v[2] / norm(v))) / (1 - f) ** 2) / DEG;
      const dl = Math.abs(((lon - e.lon_deg + 540) % 360) - 180), db = Math.abs(lat - e.lat_deg);
      worst[e.body] = Math.max(worst[e.body] || 0, dl, db);
    }
    for (const [body, w] of Object.entries(worst)) check(`sub-Earth point of ${body} vs Horizons (IAU rotation)`, w, 0.02, 'deg');
  }
} else {
  skip('ephemeris checks: run python data/acquisition/fetch.py jpl-de440s');
}

// ---- the Sun seen from the Earth ----------------------------------------------
{
  const points = fixture('subsolar.json').points;
  let worstLon = 0, worstLat = 0, worstAngle = 0;
  for (const f of points) {
    const ms = (f.jd_ut - 2440587.5) * 86400000;
    const s = subsolarPoint(ms);
    const dl = Math.abs(wrapLongitude(s.lon - f.lon_deg));
    worstLon = Math.max(worstLon, dl);
    worstLat = Math.max(worstLat, Math.abs(s.lat - f.lat_deg));
    worstAngle = Math.max(worstAngle, centralAngle(s.lon, s.lat, f.lon_deg, f.lat_deg) / DEG);
  }
  check(`sun: subsolar longitude, ${points.length} instants 1955-2045 (Horizons)`, worstLon, 0.03, 'deg');
  check('sun: subsolar latitude (Horizons)', worstLat, 0.01, 'deg');
  check('sun: subsolar point, great-circle error (Horizons)', worstAngle, 0.03, 'deg');
  // Meeus example 12.a: 1987 April 10, 0h UT, GMST = 13h10m46.3668s.
  check('sun: GMST, Meeus example 12.a', (gmst(2446895.5) - (13 + 10 / 60 + 46.3668 / 3600) * 15) * 3600, 0.01, 'arcsec');
  check('sun: Julian date of the Unix epoch', julianDate(0) - 2440587.5, 0, 'd');
  // At the subsolar point the Sun is at the zenith; 90 degrees away, on the horizon.
  const sun = subsolarPoint(Date.UTC(2026, 9, 22, 6));
  check('sun: altitude at the subsolar point', sunAltAz(sun.lon, sun.lat, sun).altitude - 90, 1e-6, 'deg');
  const [e] = eastNorth(sun.lon, sun.lat);
  const away = fromVector(e);
  check('sun: altitude on the terminator', sunAltAz(away[0], away[1], sun).altitude, 1e-6, 'deg');
  // Due south of the subsolar point the Sun stands due north.
  const south = sunAltAz(sun.lon, sun.lat - 30, sun);
  check('sun: azimuth seen from the south', Math.min(south.azimuth, 360 - south.azimuth), 1e-6, 'deg');
  assert('sun: daylight classes', daylight(10) === 'day' && daylight(-3) === 'civil twilight' && daylight(-15) === 'astronomical twilight' && daylight(-30) === 'night');
  check('sun: local mean solar time at 90 E at 06:00 UTC', localMeanSolarTime(Date.UTC(2026, 0, 1, 6), 90) - 12, 1e-9, 'h');
}

// ---- galaxies -------------------------------------------------------------------

function readCsv(path) {
  const split = (line) => {
    const out = [];
    let cur = '', quoted = false;
    for (const ch of line) {
      if (ch === '"') quoted = !quoted;
      else if (ch === ',' && !quoted) { out.push(cur); cur = ''; } else cur += ch;
    }
    out.push(cur);
    return out;
  };
  const lines = readFileSync(path, 'utf8').trim().split(/\r?\n/);
  const head = split(lines[0]);
  return lines.slice(1).map((l) => Object.fromEntries(split(l).map((v, i) => [head[i], v])));
}
const num = (s) => (s === '' ? NaN : Number(s));
const percentile = (sorted, p) => sorted[Math.min(sorted.length - 1, Math.floor(p * sorted.length))];

// Supergalactic coordinates against Tully (2015), which lists for every
// galaxy its ICRS position, supergalactic angles and SGX, SGY, SGZ.
// Known anomaly: PGC 4665585's VizieR-computed position is the Galactic centre
// (a placeholder), inconsistent with its own supergalactic coordinates.
// Elsewhere the catalogue is self-consistent to about 0.01°, without pattern on
// the sky; a wrong transformation would show a systematic offset everywhere.
const TULLY = resolve(CACHE, 'tully-2mass-groups/all/result.csv');
if (existsSync(TULLY)) {
  const anomalies = new Set(['4665585']);
  const angles = [];
  let worstXYZ = 0;
  for (const r of readCsv(TULLY)) {
    if (anomalies.has(r.PGC)) continue;
    const u = apply(ICRS_TO_SUPERGALACTIC, fromSpherical(num(r._RA_icrs) * DEG, num(r._DE_icrs) * DEG));
    angles.push(separation(u, fromSpherical(num(r.SGLON) * DEG, num(r.SGLAT) * DEG)) / DEG);
    if (Number.isFinite(num(r.SGX))) {
      const p = supergalactic(num(r._RA_icrs), num(r._DE_icrs), num(r.Dist));
      worstXYZ = Math.max(worstXYZ, Math.abs(p[0] - num(r.SGX)), Math.abs(p[1] - num(r.SGY)), Math.abs(p[2] - num(r.SGZ)));
    }
  }
  angles.sort((a, b) => a - b);
  check(`supergalactic direction vs Tully (2015), median of ${angles.length}`, percentile(angles, 0.5), 0.0005, 'deg');
  check('supergalactic direction, 99.9th percentile', percentile(angles, 0.999), 0.005, 'deg');
  check('supergalactic direction, maximum', angles[angles.length - 1], 0.02, 'deg');
  check('supergalactic X, Y, Z vs Tully (2015), maximum', worstXYZ, 0.03, 'h^-1 Mpc');
} else {
  skip('supergalactic checks: run python data/acquisition/fetch.py tully-2mass-groups');
}

// Tidal indices and main disturbers against the Updated Nearby Galaxy
// Catalog (Karachentsev et al. 2013), recomputed from positions and
// luminosities; the Milky Way's missing luminosity is inferred from the
// published indices of its satellites.
const UNGC = resolve(CACHE, 'ungc/all/result.csv');
if (existsSync(UNGC)) {
  const rows = readCsv(UNGC);
  const index = new Map(rows.map((r, i) => [r.Name, i]));
  const positions = new Float64Array(3 * rows.length);
  rows.forEach((r, i) => positions.set(supergalactic(num(r.RAJ2000), num(r.DEJ2000), num(r.Dist)), 3 * i));
  const logL = Float64Array.from(rows, (r) => num(r.KLum));
  const published = Int32Array.from(rows, (r) => (index.has(r.MD) ? index.get(r.MD) : -1));
  const publishedTheta = Float64Array.from(rows, (r) => num(r.Ti1));
  const filled = inferDisturberLuminosities(positions, logL, published, publishedTheta);
  const mw = filled.inferred.find((x) => rows[x.index].Name === 'Milky Way');
  check('inferred log L_K of the Milky Way within 10.55 ± 0.3', mw ? mw.logL - 10.55 : 99, 0.3, 'dex');
  const ours = tidalIndices(positions, filled.logL);
  let comparable = 0, same = 0;
  const diffs = [];
  rows.forEach((r, i) => {
    if (published[i] < 0 || !Number.isFinite(publishedTheta[i])) return;
    comparable++;
    if (ours.disturber[i] === published[i]) {
      same++;
      diffs.push(Math.abs(ours.theta1[i] - publishedTheta[i]));
    }
  });
  diffs.sort((a, b) => a - b);
  check(`main disturber identical for ≥ 98 % of ${comparable} galaxies (miss rate)`, 1 - same / comparable, 0.02, '');
  check('median |Θ1 − published| (published values have 0.1 precision)', percentile(diffs, 0.5), 0.05, 'dex');
  check('90th percentile |Θ1 − published|', percentile(diffs, 0.9), 0.1, 'dex');
  const groups = tidalGroups(ours.disturber, ours.theta1, filled.logL);
  const leaders = groups.groups.slice(0, 6).map((g) => rows[g.leader].Name);
  check('M 81 and M 31 lead two of the six largest tidal groups',
    ['MESSIER081', 'MESSIER031'].every((n) => leaders.includes(n)) ? 0 : 1, 0, '');
} else {
  skip('tidal-index checks: run python data/acquisition/fetch.py ungc');
}

// Names.
{
  const cases = [['MESSIER081', 'M 81'], ['NGC0253', 'NGC 253'], ['ESO499-037', 'ESO 499-037'],
    ['UGCA086', 'UGCA 86'], ['KKH011', 'KKH 11'], ['Milky Way', 'Milky Way'], ['And XXX', 'And XXX']];
  check('display names', cases.filter(([a, b]) => displayName(a) !== b).length, 0, 'mismatches');
}

// The C++.
cppChecks({ assert, skip, lines }, { project: resolve(dirname(fileURLToPath(import.meta.url)), '..'), name: 'astronomy' });

done();
