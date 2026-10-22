// Tests for the Earth Atlas's provisional computational core.
// Run from the repository root:  node applications/geography/earth-atlas/tests/run-tests.mjs
// The checks on real data need the products of data/products/ (see the README);
// they are skipped, and say so, when a product is missing.

import { existsSync, readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { DEG, EARTH_RADIUS_M, toVector, fromVector, eastNorth, centralAngle, distanceM, move, wrapLongitude } from '../src/core/sphere.js';
import { subsolarPoint, sunAltAz, gmst, julianDate, daylight, localMeanSolarTime } from '../src/core/sun.js';
import { equalEarth, equirectangular, forward, inverse, extent, areaScale } from '../src/core/projections.js';
import { gridPosition, bilinear, timeBracket, FieldProduct, areaMean } from '../src/core/grid.js';
import { isolines, levels, extrema } from '../src/core/contour.js';
import { fillRings, indexRaster, maskRaster, areaFraction, featureAreaFraction, waterMask } from '../src/core/raster.js';
import { ParticleSystem, uniformSpawn } from '../src/core/flow.js';
import { windFromDirection, compassPoint, deepWaterGroupSpeed, waveEnergyVelocity, meanDirection } from '../src/core/derived.js';
import { displayUnit, toShown, toStored, format } from '../src/core/units.js';
import { joinCountries, valueForYear, latestWellCoveredYear } from '../src/core/join.js';
import { gaussLegendre, SphericalHarmonics, EARTH_RADIUS as RADIUS, OMEGA } from '../src/core/spectral.js';
import { BarotropicModel, GRAVITY } from '../src/core/barotropic.js';
import { ShallowWaterModel, lanczosWeights } from '../src/core/shallowwater.js';
import { FFT } from '../src/core/spectral.js';
import { detect, stitch, arcDeg, positionAt } from '../src/core/cyclones.js';
import { koppen, koppenQuantities, koppenGrid, KOPPEN_CODES, dryMonths, meanTemperature } from '../src/core/climate.js';
import { ClassRaster } from '../src/core/grid.js';
import { precipitationScale, climateDiagram } from '../src/ui/climate.js';
import { inertialPeriodHours, flowToward, mixedLayerDepth, columnSampler, symmetricRange, atDepth } from '../src/core/ocean.js';
import { NNR_MORVEL56, poleOf, omega, poleOfOmega, velocityAt, subtract, acrossStep, azimuth, morvelCode, itrfOmega, itrfVelocityAt } from '../src/core/plates.js';
import { mechanism, tensorNED, eigenSymmetric, radiation, fromEqualArea, toEqualArea, planeDifference } from '../src/core/moment.js';
import { cellIndex, column as crustColumn } from '../src/core/crust.js';
import { travel, relativeTo, sectionLine, R_KM } from '../src/core/section.js';

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(here, '../../../..');
const products = resolve(root, 'data/cache/products');

let failures = 0, checks = 0;
function check(name, value, limit, unit = '') {
  checks++;
  const ok = Number.isFinite(value) && Math.abs(value) <= limit;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'}  ${name}: ${Number.isFinite(value) ? value.toExponential(3) : value} ${unit} (limit ${limit} ${unit})`);
}
function assert(name, condition, detail = '') {
  checks++;
  if (!condition) failures++;
  console.log(`${condition ? 'ok  ' : 'FAIL'}  ${name}${detail ? `: ${detail}` : ''}`);
}
function seeded(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
}
const random = seeded(20261022);

// ---- the sphere ------------------------------------------------------------
{
  let worst = 0;
  for (let i = 0; i < 1000; i++) {
    const lon = random() * 360 - 180, lat = Math.asin(2 * random() - 1) / DEG;
    const [l2, p2] = fromVector(toVector(lon, lat));
    worst = Math.max(worst, Math.abs(wrapLongitude(l2 - lon)) * Math.cos(lat * DEG), Math.abs(p2 - lat));
  }
  check('sphere: vector round trip', worst, 1e-10, 'deg');
  const [e, n] = eastNorth(35, 40);
  const p = toVector(35, 40);
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  check('sphere: east, north and up are orthonormal', Math.max(Math.abs(dot(e, n)), Math.abs(dot(e, p)), Math.abs(dot(n, p)), Math.abs(dot(e, e) - 1), Math.abs(dot(n, n) - 1)), 1e-12);
  // One degree of a great circle is R pi / 180 = 111.195 km.
  check('sphere: one degree of latitude', distanceM(10, 0, 10, 1) - (EARTH_RADIUS_M * Math.PI) / 180, 1e-6, 'm');
  // Paris to Tokyo, 9 713 km on the sphere of mean radius (independently: 0.152 rad x R).
  const angle = centralAngle(2.3522, 48.8566, 139.6917, 35.6895);
  check('sphere: haversine equals the spherical law of cosines', angle - Math.acos(dot(toVector(2.3522, 48.8566), toVector(139.6917, 35.6895))), 1e-12, 'rad');
  // Moving along the equator by a quarter of the circumference.
  const q = move(0, 0, (Math.PI / 2) * EARTH_RADIUS_M, 0);
  check('sphere: a quarter turn east along the equator', Math.hypot(q[0] - 90, q[1]), 1e-9, 'deg');
  // Moving north through the pole comes down the opposite meridian.
  const r = move(30, 89, 0, (2 * Math.PI / 180) * EARTH_RADIUS_M);
  check('sphere: through the north pole', Math.hypot(wrapLongitude(r[0] - (-150)), r[1] - 89), 1e-9, 'deg');
}

// ---- the Sun ---------------------------------------------------------------
{
  const fixtures = JSON.parse(readFileSync(resolve(here, 'fixtures/subsolar.json'), 'utf8')).points;
  let worstLon = 0, worstLat = 0, worstAngle = 0;
  for (const f of fixtures) {
    const ms = (f.jd_ut - 2440587.5) * 86400000;
    const s = subsolarPoint(ms);
    const dl = Math.abs(wrapLongitude(s.lon - f.lon_deg));
    worstLon = Math.max(worstLon, dl);
    worstLat = Math.max(worstLat, Math.abs(s.lat - f.lat_deg));
    worstAngle = Math.max(worstAngle, centralAngle(s.lon, s.lat, f.lon_deg, f.lat_deg) / DEG);
  }
  check(`sun: subsolar longitude, ${fixtures.length} instants 1955-2045 (Horizons)`, worstLon, 0.03, 'deg');
  check('sun: subsolar latitude (Horizons)', worstLat, 0.01, 'deg');
  check('sun: subsolar point, great-circle error (Horizons)', worstAngle, 0.03, 'deg');
  // Meeus example 12.a: 1987 April 10, 0h UT, GMST = 13h10m46.3668s.
  check('sun: GMST, Meeus example 12.a', (gmst(2446895.5) - (13 + 10 / 60 + 46.3668 / 3600) * 15) * 3600, 0.01, 'arcsec');
  check('sun: Julian date of the Unix epoch', julianDate(0) - 2440587.5, 0, 'd');
  // At the subsolar point the Sun is at the zenith; 90 degrees away, on the horizon.
  const sub = subsolarPoint(Date.UTC(2026, 9, 22, 6));
  check('sun: altitude at the subsolar point', sunAltAz(sub.lon, sub.lat, sub).altitude - 90, 1e-6, 'deg');
  const [e] = eastNorth(sub.lon, sub.lat);
  const away = fromVector(e);
  check('sun: altitude on the terminator', sunAltAz(away[0], away[1], sub).altitude, 1e-6, 'deg');
  // Due south of the subsolar point the Sun stands due north.
  const south = sunAltAz(sub.lon, sub.lat - 30, sub);
  check('sun: azimuth seen from the south', Math.min(south.azimuth, 360 - south.azimuth), 1e-6, 'deg');
  assert('sun: daylight classes', daylight(10) === 'day' && daylight(-3) === 'civil twilight' && daylight(-15) === 'astronomical twilight' && daylight(-30) === 'night');
  check('sun: local mean solar time at 90 E at 06:00 UTC', localMeanSolarTime(Date.UTC(2026, 0, 1, 6), 90) - 12, 1e-9, 'h');
}

// ---- projections -----------------------------------------------------------
{
  for (const p of [equalEarth, equirectangular]) {
    let worst = 0;
    for (let i = 0; i < 2000; i++) {
      const lon = random() * 359.9 - 179.95, lat = random() * 179.8 - 89.9, lambda0 = random() * 360 - 180;
      const [x, y] = forward(p, lon, lat, lambda0);
      const back = inverse(p, x, y, lambda0);
      worst = Math.max(worst, Math.abs(wrapLongitude(back[0] - lon)) * Math.cos(lat * DEG), Math.abs(back[1] - lat));
    }
    check(`projection ${p.name}: forward then inverse`, worst, 1e-9, 'deg');
  }
  // Equal area: the Jacobian over cos(phi) is the same everywhere (Savric et al. 2018).
  let lo = Infinity, hi = -Infinity;
  for (let lat = -85; lat <= 85; lat += 5) for (let lon = -175; lon <= 175; lon += 35) {
    const s = areaScale(equalEarth, lon, lat);
    lo = Math.min(lo, s); hi = Math.max(hi, s);
  }
  check('projection Equal Earth: area scale constant within', (hi - lo) / lo, 1e-5);
  // The paper gives the ratio of the axes as about 2.05 : 1.
  const [xm, ym] = extent(equalEarth);
  check('projection Equal Earth: ratio of the axes, 2.05', xm / ym - 2.05, 0.006);
  check('projection Equal Earth: area of the map equals that of the sphere (4 pi)', areaScale(equalEarth, 0, 0) - 1, 1e-6);
  assert('projection: a point beyond the outline has no place', inverse(equalEarth, xm * 0.99, ym * 0.99) === null);
}

// ---- grids -----------------------------------------------------------------
{
  // A linear function of longitude and latitude is reproduced exactly by bilinear interpolation.
  const grid = { lon0: -180, dlon: 2, nlon: 180, lat0: 89, dlat: 2, nlat: 90, periodic: true };
  const values = new Float32Array(grid.nlon * grid.nlat);
  for (let r = 0; r < grid.nlat; r++) for (let c = 0; c < grid.nlon; c++) values[r * grid.nlon + c] = 3 * (grid.lat0 - r * grid.dlat) + 0.5 * c;
  let worst = 0;
  for (let i = 0; i < 500; i++) {
    const lon = -180 + random() * 356, lat = -88 + random() * 176;
    worst = Math.max(worst, Math.abs(bilinear(values, grid, lon, lat) - (3 * lat + 0.5 * ((lon + 180) / 2))));
  }
  check('grid: bilinear interpolation of a linear field', worst, 1e-3);
  // Periodic: between the last column (178 E) and the first (180 W), halfway.
  const wrap = bilinear(values, grid, 179, 1);
  check('grid: interpolation across the antimeridian', wrap - (3 * 1 + 0.5 * (179 + 0) / 2), 1e-4);
  const pos = gridPosition(grid, -181, 0);
  check('grid: longitude wraps', pos[0] - 179.5, 1e-9);
  // A missing corner: weights renormalised while they carry half the weight.
  const g2 = { lon0: 0, dlon: 1, nlon: 2, lat0: 1, dlat: 1, nlat: 2, periodic: false };
  const v2 = new Float32Array([1, NaN, 3, 5]);
  check('grid: renormalised weights near a missing value', bilinear(v2, g2, 0.25, 0.25) - ((0.75 * 0.25 * 1 + 0.75 * 0.75 * 3 + 0.25 * 0.75 * 5) / (1 - 0.25 * 0.25)), 1e-6);
  assert('grid: no value next to a missing value', Number.isNaN(bilinear(v2, g2, 0.9, 0.9)));
  const times = [0, 3600e3, 7200e3];
  const b = timeBracket(times, 5400e3);
  assert('grid: time bracket', b.i0 === 1 && b.i1 === 2 && Math.abs(b.w - 0.5) < 1e-12);
  const c = timeBracket(times, 9000e3);
  assert('grid: time after the last is clamped, with its offset', c.i0 === 2 && c.offset === 1800e3);
  assert('grid: nearest time for means', timeBracket(times, 2000e3, true).i0 === 1);
  // Area mean of the field cos^0 = 1 is 1; of sin(lat) is 0.
  const ones = new Float32Array(values.length).fill(1);
  check('grid: area mean of a constant', areaMean(ones, grid) - 1, 1e-6);
}

// ---- contours --------------------------------------------------------------
{
  // The isoline of f = distance from (0, 0) is a circle (level chosen off the grid values).
  const grid = { lon0: -5, dlon: 0.1, nlon: 101, lat0: 5, dlat: 0.1, nlat: 101, periodic: false };
  const values = new Float32Array(grid.nlon * grid.nlat);
  for (let r = 0; r < grid.nlat; r++) for (let c = 0; c < grid.nlon; c++) values[r * grid.nlon + c] = Math.hypot(grid.lon0 + c * grid.dlon, grid.lat0 - r * grid.dlat);
  const R0 = 2.95;
  const seg = isolines(values, grid, R0);
  let worst = 0, length = 0;
  for (let i = 0; i < seg.length; i += 4) {
    worst = Math.max(worst, Math.abs(Math.hypot(seg[i], seg[i + 1]) - R0), Math.abs(Math.hypot(seg[i + 2], seg[i + 3]) - R0));
    length += Math.hypot(seg[i + 2] - seg[i], seg[i + 3] - seg[i + 1]);
  }
  check('contour: points of the isoline lie on the level', worst, 2e-3, 'deg');
  check('contour: length of the circle 2 pi r', length / (2 * Math.PI * R0) - 1, 1e-3);
  // Every segment end is shared by exactly one other segment: the line is closed.
  const ends = new Map();
  for (let i = 0; i < seg.length; i += 2) {
    const k = `${seg[i].toFixed(5)},${seg[i + 1].toFixed(5)}`;
    ends.set(k, (ends.get(k) || 0) + 1);
  }
  assert('contour: the isoline is closed', [...ends.values()].every((n) => n === 2), `${ends.size} vertices`);
  assert('contour: levels at multiples of a step', JSON.stringify(levels([3.2, 9.9], 2)) === '[4,6,8]');
  // A single high and a single low.
  const g = { lon0: -180, dlon: 2, nlon: 180, lat0: 89, dlat: 2, nlat: 90, periodic: true };
  const p = new Float32Array(g.nlon * g.nlat);
  for (let r = 0; r < g.nlat; r++) for (let c = 0; c < g.nlon; c++) {
    const lon = g.lon0 + c * g.dlon, lat = g.lat0 - r * g.dlat;
    p[r * g.nlon + c] = 101300 + 2000 * Math.exp(-((lon - 40) ** 2 + (lat - 45) ** 2) / 200) - 2500 * Math.exp(-((lon + 30) ** 2 + (lat + 51) ** 2) / 150);
  }
  const ex = extrema(p, g, { radius: 6, prominence: 300 });
  const high = ex.find((e) => e.kind === 'high'), low = ex.find((e) => e.kind === 'low');
  assert('contour: one high and one low found', ex.length === 2 && high && low, JSON.stringify(ex.map((e) => [e.kind, e.lon, e.lat])));
  if (high && low) check('contour: positions of the high and the low', Math.hypot(high.lon - 40, high.lat - 45) + Math.hypot(low.lon + 30, low.lat + 51), 1e-9, 'deg');
}

// ---- rasters ---------------------------------------------------------------
{
  // A square 10 x 10 degrees with a 4 x 4 hole, on a 1-degree raster.
  const coords = new Float32Array([0, 0, 10, 0, 10, 10, 0, 10, 0, 0, 3, 3, 3, 7, 7, 7, 7, 3, 3, 3]);
  const target = new Uint8Array(360 * 180);
  const n = fillRings(coords, [[0, 5], [5, 5]], 360, 180, 1, target);
  assert('raster: a square with a hole fills 100 - 16 pixels', n === 84, `${n}`);
  // Two squares sharing an edge claim each of its pixels once.
  const two = indexRaster(new Float32Array([0, 0, 5, 0, 5, 5, 0, 5, 5, 0, 10, 0, 10, 5, 5, 5]), [{ rings: [[0, 4]] }, { rings: [[4, 4]] }], 360, 180);
  const counts = [0, 0, 0];
  two.forEach((v) => counts[v]++);
  assert('raster: shared edges are not filled twice', counts[1] === 25 && counts[2] === 25, `${counts}`);
}

// ---- particles -------------------------------------------------------------
{
  // A uniform eastward flow of 10 m/s moves a particle on the equator by 10 dt metres.
  const ps = new ParticleSystem(1, { minLife: 1000, maxLife: 1001, random: seeded(1) });
  ps.reset(() => [0, 0]);
  ps.lon[0] = 0; ps.lat[0] = 0; ps.age[0] = 0;
  for (let i = 0; i < 100; i++) ps.step(() => [10, 0], 3600, () => [0, 0]);
  check('flow: uniform flow along the equator, 100 h', ps.lon[0] - (10 * 3600 * 100) / EARTH_RADIUS_M / DEG, 1e-9, 'deg');
  // Solid-body rotation about the pole: a particle keeps its latitude (the midpoint rule).
  const omega = 1e-5;   // rad/s
  const ring = new ParticleSystem(1, { minLife: 100000, maxLife: 100001, random: seeded(2) });
  ring.reset(() => [0, 60]);
  ring.lon[0] = 0; ring.lat[0] = 60; ring.age[0] = 0;
  for (let i = 0; i < 2000; i++) ring.step((lon, lat) => [omega * EARTH_RADIUS_M * Math.cos(lat * DEG), 0], 600, () => [0, 60]);
  check('flow: solid-body rotation keeps the latitude (2000 steps)', ring.lat[0] - 60, 0.05, 'deg');
  // A particle dies where the field has no value.
  const dead = new ParticleSystem(1, { random: seeded(3) });
  dead.reset(() => [5, 5]);
  dead.step(() => null, 60, () => [7, 7]);
  assert('flow: rebirth where there is no flow', dead.alive[0] === 0 && dead.lon[0] === 7);
  // Uniform spawning: the mean of sin(latitude) is near 0.
  const spawn = uniformSpawn();
  let s = 0;
  const rnd = seeded(4);
  for (let i = 0; i < 20000; i++) s += Math.sin(spawn(rnd)[1] * DEG);
  check('flow: spawning uniform in area', s / 20000, 0.01);
}

// ---- derived quantities and units --------------------------------------------
{
  assert('derived: wind from the west is 270', Math.abs(windFromDirection(5, 0) - 270) < 1e-12);
  assert('derived: wind from the north is 0', Math.abs(windFromDirection(0, -5)) < 1e-12);
  assert('derived: compass points', compassPoint(225) === 'SW' && compassPoint(359) === 'N' && compassPoint(100) === 'E');
  check('derived: deep-water group speed g T / 4 pi, T = 10 s', deepWaterGroupSpeed(10) - 7.80388, 1e-5, 'm/s');
  const [u, v] = waveEnergyVelocity(270, 10);           // waves from the west travel east
  check('derived: waves from the west travel east', Math.hypot(u - deepWaterGroupSpeed(10), v), 1e-12);
  check('derived: mean of 350 and 10 degrees is 0', Math.min(meanDirection([350, 10]), 360 - meanDirection([350, 10])), 1e-9, 'deg');
  const k = displayUnit('K');
  check('units: 300 K in degrees Celsius', toShown(300, k) - 26.85, 1e-9);
  check('units: back to kelvin', toStored(26.85, k) - 300, 1e-9);
  check('units: 1e-3 kg m-2 s-1 is 3.6 mm/h', toShown(1e-3, displayUnit('kg m-2 s-1')) - 3.6, 1e-12);
  check('units: 10 m/s in knots', toShown(10, displayUnit('m s-1', 'kn')) - 19.4384, 1e-4);
  assert('units: numbers for reading', format(124516650) === '125 M' && format(3.14159) === '3.14' && format(1013.2) === '1013');
}

// ---- spherical harmonics and the barotropic model -------------------------------------
{
  const { mu, w } = gaussLegendre(64);
  check('spectral: Gauss-Legendre weights sum to 2', w.reduce((s, x) => s + x, 0) - 2, 1e-13);
  check('spectral: quadrature of mu^2 is 2/3', w.reduce((s, x, i) => s + x * mu[i] ** 2, 0) - 2 / 3, 1e-13);
  const sh = new SphericalHarmonics(42);
  let worst = 0;
  for (let k = 0; k < sh.size; k++) {
    let s = 0;
    for (let j = 0; j < sh.nlat; j++) s += sh.w[j] * sh.P[j * sh.size + k] ** 2;
    worst = Math.max(worst, Math.abs(s - 1));
  }
  check('spectral: associated Legendre functions are normalised (T42)', worst, 1e-12);
  const rnd = seeded(42);
  const r = sh.zeros();
  for (let k = 0; k < sh.size; k++) { r.re[k] = rnd() - 0.5; r.im[k] = sh.m[k] === 0 ? 0 : rnd() - 0.5; }
  const back = sh.analyse(sh.synthesise(r));
  let rt = 0;
  for (let k = 0; k < sh.size; k++) rt = Math.max(rt, Math.abs(back.re[k] - r.re[k]), Math.abs(back.im[k] - r.im[k]));
  check('spectral: synthesis then analysis returns the coefficients', rt, 1e-12);
  // Solid-body rotation u = u0 cos(phi): psi = -a u0 mu, vorticity 2 u0 mu / a.
  const u0 = 20, psi = sh.zeros();
  psi.re[sh.index(0, 1)] = (-RADIUS * u0) / Math.sqrt(1.5);
  const { U, V } = sh.windFromStreamFunction(psi);
  let eu = 0, ev = 0;
  for (let j = 0; j < sh.nlat; j++) for (let k = 0; k < sh.nlon; k++) {
    eu = Math.max(eu, Math.abs(U[j * sh.nlon + k] - u0 * (1 - sh.mu[j] ** 2)));
    ev = Math.max(ev, Math.abs(V[j * sh.nlon + k]));
  }
  check('spectral: wind of a solid-body rotation', Math.max(eu, ev), 1e-12, 'm/s');
  const zg = sh.synthesise(sh.curl(U, V));
  let ez = 0;
  for (let j = 0; j < sh.nlat; j++) ez = Math.max(ez, Math.abs(zg[j * sh.nlon] - (2 * u0 * sh.mu[j]) / RADIUS));
  check('spectral: curl of that wind is 2 u0 sin(phi) / a', ez / ((2 * u0) / RADIUS), 1e-12);
  // The curl of the wind of a random stream function is its Laplacian; its divergence is zero.
  const ps = sh.zeros();
  for (let k = 0; k < sh.size; k++) if (sh.n[k] > 0 && sh.n[k] < 30) { ps.re[k] = 1e7 * (rnd() - 0.5); ps.im[k] = sh.m[k] === 0 ? 0 : 1e7 * (rnd() - 0.5); }
  const wind = sh.windFromStreamFunction(ps);
  const curl = sh.curl(wind.U, wind.V), div = sh.divergence(wind.U, wind.V);
  const lap = sh.scale({ re: Float64Array.from(ps.re), im: Float64Array.from(ps.im) }, (n) => sh.laplacianEigen(n));
  let ec = 0, ed = 0, scale = 0;
  for (let k = 0; k < sh.size; k++) {
    ec = Math.max(ec, Math.hypot(curl.re[k] - lap.re[k], curl.im[k] - lap.im[k]));
    ed = Math.max(ed, Math.hypot(div.re[k], div.im[k]));
    scale = Math.max(scale, Math.hypot(lap.re[k], lap.im[k]));
  }
  check('spectral: curl of the wind of psi is the Laplacian of psi', ec / scale, 1e-10);
  check('spectral: the wind of psi has no divergence', ed / scale, 1e-10);

  // The Rossby-Haurwitz wave of wavenumber 4 (Williamson et al. 1992, test 6) is an
  // exact solution of the barotropic vorticity equation, travelling east at
  // nu = (R (3 + R) omega - 2 Omega) / ((1 + R)(2 + R)) (Haurwitz 1940).
  const R = 4, om = 7.848e-6, K = 7.848e-6;
  const nu = (R * (3 + R) * om - 2 * OMEGA) / ((1 + R) * (2 + R));
  const rh = (t) => {
    const g = new Float64Array(sh.nlat * sh.nlon);
    for (let j = 0; j < sh.nlat; j++) for (let k = 0; k < sh.nlon; k++) {
      const m = sh.mu[j], lam = (2 * Math.PI * k) / sh.nlon - nu * t;
      g[j * sh.nlon + k] = -RADIUS * RADIUS * om * m + RADIUS * RADIUS * K * (1 - m * m) ** 2 * m * Math.cos(R * lam);
    }
    return sh.analyse(g);
  };
  const model = new BarotropicModel({ truncation: 42, dampingHours: 0 });
  model.setStreamFunction(rh(0));
  for (let t = 0; t < 86400; t += 1800) model.step(1800);
  const exact = sh.scale(rh(86400), (n) => sh.laplacianEigen(n));
  let err = 0, mx = 0;
  for (let k = 0; k < sh.size; k++) {
    err = Math.max(err, Math.hypot(model.zeta.re[k] - exact.re[k], model.zeta.im[k] - exact.im[k]));
    mx = Math.max(mx, Math.hypot(exact.re[k], exact.im[k]));
  }
  check('barotropic: Rossby-Haurwitz wave after one day (relative error)', err / mx, 1e-6);
  const k54 = sh.index(4, 5);
  const phase = Math.atan2(model.zeta.im[k54], model.zeta.re[k54]) - Math.atan2(exact.im[k54], exact.re[k54]);
  check('barotropic: its phase speed, 12.2 degrees a day', Math.atan2(Math.sin(phase), Math.cos(phase)), 1e-6, 'rad');
  // A random flow: energy and enstrophy are kept (the spectral method keeps both
  // exactly; Runge-Kutta loses a little at the smallest scales).
  const m2 = new BarotropicModel({ truncation: 42, dampingHours: 0 });
  const z = sh.zeros();
  for (let k = 0; k < sh.size; k++) if (sh.n[k] > 0 && sh.n[k] <= 15) { z.re[k] = 1e-5 * (rnd() - 0.5); z.im[k] = sh.m[k] === 0 ? 0 : 1e-5 * (rnd() - 0.5); }
  m2.zeta = z;
  const i0 = m2.invariants();
  for (let t = 0; t < 86400; t += 1800) m2.step(1800);
  const i1 = m2.invariants();
  check('barotropic: energy kept over a day', (i1.energy - i0.energy) / i0.energy, 1e-4);
  check('barotropic: enstrophy kept over a day', (i1.enstrophy - i0.enstrophy) / i0.enstrophy, 2e-3);
  // Linear balance of a zonal flow u0 cos(phi): z = -(Omega a u0 / g) (mu^2 - 1/3).
  const zb = sh.synthesise(model.balancedHeight(psi));
  let eb = 0;
  for (let j = 0; j < sh.nlat; j++) eb = Math.max(eb, Math.abs(zb[j * sh.nlon] + ((OMEGA * RADIUS * u0) / GRAVITY) * (sh.mu[j] ** 2 - 1 / 3)));
  check('barotropic: balanced height of a zonal flow', eb, 1e-6, 'm');
}

// ---- joins on the real products ------------------------------------------------
{
  const geo = resolve(products, 'geography-ne50m.json');
  const ind = resolve(products, 'indicator-population.json');
  if (existsSync(geo) && existsSync(ind)) {
    const g = JSON.parse(readFileSync(geo, 'utf8'));
    const p = JSON.parse(readFileSync(ind, 'utf8'));
    const { codes, unmatched } = joinCountries(g.layers.countries, p.economies);
    const matched = codes.filter(Boolean).length;
    // Economies with no country of their own in Natural Earth at 1:50 million.
    const expected = ['GIB'];
    assert(`join: ${matched} of ${codes.length} countries carry a World Bank economy`, matched >= 200);
    assert('join: the unmatched economies are the expected ones', JSON.stringify(unmatched) === JSON.stringify(expected), unmatched.join(' '));
    const kosovo = codes[g.layers.countries.findIndex((c) => c.id === 'KOS')];
    assert('join: Kosovo by its ISO alpha-2 code (XK)', kosovo === 'XKX', String(kosovo));
    const france = codes[g.layers.countries.findIndex((c) => c.id === 'FRA')];
    assert('join: France, whose ISO_A3 is -99 in Natural Earth', france === 'FRA', String(france));
    const r = valueForYear(p.economies.JPN, p.years, 2030, 10);
    assert('join: value for a later year falls back to the latest', r && r.year === p.years[p.years.length - 1]);
    assert('join: a well-covered year exists', latestWellCoveredYear(p.economies, p.years) >= 2020);
  } else {
    console.log('skip  join: products missing (geography-ne50m, indicator-population)');
  }
}

// ---- the land fraction from the geography product -------------------------------
{
  const geo = resolve(products, 'geography-ne50m.json');
  if (existsSync(geo)) {
    const g = JSON.parse(readFileSync(geo, 'utf8'));
    const bin = readFileSync(resolve(products, g.coordinates.file));
    const coords = new Float32Array(bin.buffer, bin.byteOffset, bin.byteLength / 4);
    const W = 2048, H = 1024;
    const land = maskRaster(coords, g.layers.land, W, H);
    const rastered = areaFraction(land, W, H, (v) => v === 1);
    const exact = g.layers.land.reduce((s, f) => s + featureAreaFraction(coords, f.rings), 0);
    check('raster: land area from the raster equals that of the polygons', rastered - exact, 5e-4);
    // The land of the Earth is 29.2 % of its surface (148.9 of 510.1 million km^2, with
    // Antarctica's ice shelves, which Natural Earth's land omits: about 0.3 %).
    check('raster: Natural Earth land is 28.9 +/- 0.3 % of the sphere', exact - 0.289, 0.003);
  }
}

// ---- a field product -------------------------------------------------------------
{
  const meta = resolve(products, 'fields-gfs-forecast.json');
  if (existsSync(meta)) {
    const m = JSON.parse(readFileSync(meta, 'utf8'));
    const bin = readFileSync(resolve(products, m.file));
    const f = new FieldProduct(m, bin.buffer.slice(bin.byteOffset, bin.byteOffset + bin.byteLength));
    const t = f.slice('tmp2m', 0);
    const mean = areaMean(t, f.grid) - 273.15;
    // The global mean near-surface temperature is about 14-16 degrees C.
    assert('field: GFS global mean 2 m temperature is plausible', mean > 10 && mean < 19, `${mean.toFixed(2)} °C`);
    const p = f.slice('prmslmsl', 0);
    const pm = areaMean(p, f.grid) / 100;
    assert('field: GFS mean sea-level pressure is plausible', pm > 1005 && pm < 1015, `${pm.toFixed(1)} hPa`);
    const v = m.variables.find((x) => x.name === 'tmp2m');
    check('field: quantisation error at most half a step', v.scale / 2 - (v.max - v.min) / 130000, 1e-9);
  } else {
    console.log('skip  field: fields-gfs-forecast missing');
  }
}

// ---- the fast Fourier transform ------------------------------------------------------
{
  const n = 64, fft = new FFT(n), rnd = seeded(3);
  const x = Float64Array.from({ length: n }, () => rnd() - 0.5), y = Float64Array.from({ length: n }, () => rnd() - 0.5);
  const re = Float64Array.from(x), im = Float64Array.from(y);
  fft.transform(re, im, -1);
  let err = 0;
  for (let m = 0; m < n; m++) {
    let a = 0, b = 0;
    for (let k = 0; k < n; k++) { const t = (-2 * Math.PI * m * k) / n; a += x[k] * Math.cos(t) - y[k] * Math.sin(t); b += x[k] * Math.sin(t) + y[k] * Math.cos(t); }
    err = Math.max(err, Math.abs(a - re[m]), Math.abs(b - im[m]));
  }
  check('fft: 64 points agree with the direct sum', err, 1e-12);
  fft.transform(re, im, 1);
  check('fft: forward then backward returns n times the data', Math.max(...x.map((v, k) => Math.abs(re[k] / n - v))), 1e-14);
}

// ---- the shallow-water model ----------------------------------------------------------
{
  // Williamson et al. (1992), test case 2: steady, nonlinear, zonal geostrophic
  // flow, here also about an axis tilted by 45 degrees so that it crosses the
  // grid's poles. h and the wind are exact at any truncation (degrees <= 2).
  const case2 = (model, alpha) => {
    const sh = model.sh, n = sh.nlat * sh.nlon;
    const u0 = (2 * Math.PI * RADIUS) / (12 * 86400), gh0 = 2.94e4;
    const u = new Float64Array(n), v = new Float64Array(n), z = new Float64Array(n);
    for (let j = 0; j < sh.nlat; j++) {
      const phi = Math.asin(sh.mu[j]);
      for (let k = 0; k < sh.nlon; k++) {
        const lam = (2 * Math.PI * k) / sh.nlon, i = j * sh.nlon + k;
        u[i] = u0 * (Math.cos(phi) * Math.cos(alpha) + Math.cos(lam) * Math.sin(phi) * Math.sin(alpha));
        v[i] = -u0 * Math.sin(lam) * Math.sin(alpha);
        const s = -Math.cos(lam) * Math.cos(phi) * Math.sin(alpha) + Math.sin(phi) * Math.cos(alpha);
        z[i] = (gh0 - (RADIUS * OMEGA * u0 + (u0 * u0) / 2) * s * s) / GRAVITY;
      }
    }
    return { u, v, z };
  };
  const l2 = (sh, x, y) => {
    let s = 0, t = 0;
    for (let j = 0; j < sh.nlat; j++) for (let k = 0; k < sh.nlon; k++) { const i = j * sh.nlon + k; s += sh.w[j] * (x[i] - y[i]) ** 2; t += sh.w[j] * y[i] ** 2; }
    return Math.sqrt(s / t);
  };
  for (const [alpha, what] of [[0, 'about the pole'], [Math.PI / 4, 'across the poles']]) {
    const m = new ShallowWaterModel({ truncation: 21, dampingHours: 0, tilt: alpha });
    const { u, v, z } = case2(m, alpha);
    m.setState(u, v, z, { divergent: true });
    const h0 = m.sh.synthesise(m.state.phi), mass0 = m.invariants().mass;
    for (let s = 0; s < 96; s++) m.step(900);
    check(`shallow water: Williamson case 2, steady flow ${what}, height after a day`, l2(m.sh, m.sh.synthesise(m.state.phi), h0), 1e-12);
    if (alpha) check('shallow water: mass kept', (m.invariants().mass - mass0) / mass0, 1e-13);
  }
  // Energy over a day of a random, unbalanced flow (no damping).
  {
    const m = new ShallowWaterModel({ truncation: 21, dampingHours: 0 });
    const sh = m.sh, rnd = seeded(17), n = sh.nlat * sh.nlon;
    const u = new Float64Array(n), v = new Float64Array(n), z = new Float64Array(n);
    const waves = [...Array(6)].map(() => [1 + Math.floor(5 * rnd()), 2 * Math.PI * rnd(), 1 + Math.floor(4 * rnd()), 10 * (rnd() - 0.5)]);
    for (let j = 0; j < sh.nlat; j++) {
      const phi = Math.asin(sh.mu[j]);
      for (let k = 0; k < sh.nlon; k++) {
        const lam = (2 * Math.PI * k) / sh.nlon, i = j * sh.nlon + k;
        u[i] = 20 * Math.cos(phi); z[i] = 5500;
        for (const [mm, ph, l, amp] of waves) { u[i] += amp * Math.cos(phi) * Math.cos(mm * lam + ph) * Math.cos(l * phi); z[i] += 8 * amp * Math.cos(phi) ** 2 * Math.sin(mm * lam + ph); }
      }
    }
    m.setState(u, v, z);
    const e0 = m.invariants().energy;
    for (let s = 0; s < 96; s++) m.step(900);
    check('shallow water: energy kept over a day of unbalanced flow', (m.invariants().energy - e0) / e0, 1e-6);
  }
  // Gravity waves on a non-rotating sphere: a height of degree n oscillates with
  // frequency sqrt(g H n (n + 1)) / a (Lamb 1932, sec. 307).
  const wave = (H, n, truncation = 21) => {
    const m = new ShallowWaterModel({ truncation, omega: 0, dampingHours: 0 });
    m.state.phi.re[0] = (GRAVITY * H) / Math.sqrt(0.5);
    const k = m.sh.index(0, n);
    m.state.phi.re[k] = GRAVITY;
    return { m, k, period: (2 * Math.PI * RADIUS) / Math.sqrt(GRAVITY * H * n * (n + 1)) };
  };
  {
    const { m, k, period } = wave(1000, 5);
    for (let s = 0; s < 100; s++) m.step(period / 200);
    const half = m.state.phi.re[k] / GRAVITY;
    for (let s = 0; s < 100; s++) m.step(period / 200);
    check('shallow water: a gravity wave (n = 5, H = 1 km, 20.5 h) is reversed after half a period', half + 1, 1e-5);
    check('shallow water: and restored after a whole period', m.state.phi.re[k] / GRAVITY - 1, 1e-5);
  }
  // The digital filter (Lynch and Huang 1992) removes fast waves and keeps slow ones.
  {
    const fast = wave(10000, 20), slow = wave(1000, 2);
    fast.m.initialise({ dt: 300, spanHours: 6, cutoffHours: 6 });
    slow.m.initialise({ dt: 300, spanHours: 6, cutoffHours: 6 });
    const a = fast.m.state.phi.re[fast.k] / GRAVITY, b = slow.m.state.phi.re[slow.k] / GRAVITY;
    assert('digital filter: a gravity wave of 1.7 h is removed, one of 46 h kept', Math.abs(a) < 0.01 && b > 0.98, `amplitudes ${a.toFixed(4)} and ${b.toFixed(4)}`);
    check('digital filter: its weights sum to one', lanczosWeights(12, 0.26).reduce((x, y) => x + y, 0) - 1, 1e-14);
  }
}

// ---- tropical cyclones -----------------------------------------------------------------
{
  // A 1 degree global grid of sea-level pressure: 1013 hPa, with Gaussian lows
  // p = 1013 hPa - depth exp(-(r / 2.5 deg)^2) and, optionally, a long trough.
  const grid = { lon0: -180, dlon: 1, nlon: 360, lat0: 90, dlat: 1, nlat: 181, periodic: true };
  const field = (lows, trough = false) => {
    const v = new Float32Array(grid.nlon * grid.nlat);
    for (let r = 0; r < grid.nlat; r++) {
      const lat = 90 - r;
      for (let c = 0; c < grid.nlon; c++) {
        const lon = -180 + c;
        let p = 101300;
        for (const [lo, la, depth] of lows) p -= depth * Math.exp(-((arcDeg(lon, lat, lo, la) / 2.5) ** 2));
        // A trough along 40 S: low along a line, open at its ends.
        if (trough) p -= 900 * Math.exp(-(((lat + 40) / 3) ** 2));
        v[r * grid.nlon + c] = p;
      }
    }
    return v;
  };
  const found = detect(field([[130.4, 18.3, 8000]]), grid);
  const f0 = found[0] || { lon: NaN, lat: NaN, pa: NaN };
  // The paraboloid fit places the centre between grid points.
  check('cyclones: a low of 933 hPa is found where it is, within 0.15 degrees', arcDeg(f0.lon, f0.lat, 130.4, 18.3), 0.15, 'deg');
  check('cyclones: and its central pressure within 3 hPa', Math.abs(f0.pa - 93300) / 100, 3, 'hPa');
  assert('cyclones: one candidate, not several, for one low', found.length === 1, `${found.length}`);
  const trough = detect(field([], true), grid);
  assert('cyclones: a trough without a closed contour is not a cyclone', trough.length === 0, `${trough.length} candidates`);
  // A cyclone moving 3 degrees in 6 hours for two days, and a shallow, short-lived low.
  const times = [...Array(9).keys()].map((k) => Date.UTC(2026, 9, 20, 6 * k));
  const nodes = times.map((_, k) => detect(field([[150 - 3 * k, 12 + 0.8 * k, 4000 + 300 * k], ...(k < 3 ? [[-40, 20, 300]] : [])]), grid));
  const tracks = stitch(times, nodes);
  assert('cyclones: a moving cyclone becomes one track of 9 positions', tracks.length === 1 && tracks[0].points.length === 9,
    tracks.map((t) => t.points.length).join(', '));
  const born = stitch(times, times.map((_, k) => detect(field([[150, 35, 5000 + 100 * k]]), grid)));
  assert('cyclones: a low born at 35 degrees is not tropical', born.length === 0);
  const shallow = stitch(times, times.map((_, k) => detect(field([[150, 15, 800]]), grid)));
  assert('cyclones: a low that never deepens to 1000 hPa is left out', shallow.length === 0);
  const sea = stitch(times, nodes, { isOcean: (lon) => lon < 140 });
  assert('cyclones: and one born over land, when a land-sea test is given', sea.length === 0);
  // Positions between the stored times, across the date line.
  const q = positionAt([{ t: 0, lon: 179, lat: 10, pa: 100000 }, { t: 6, lon: -179, lat: 12, pa: 99000 }], 3);
  check('cyclones: positions interpolate across the date line', Math.abs(Math.abs(q.lon) - 180) + Math.abs(q.lat - 11) + Math.abs(q.pa - 99500) / 100, 1e-9);
}

// ---- the ocean beneath a place -------------------------------------------------------
{
  // At 30 degrees the inertial period is one sidereal day (86 164.1 s), since 2 Omega sin 30 = Omega.
  check('ocean: the inertial period at 30° is a sidereal day', Math.abs(inertialPeriodHours(30) * 3600 - 86164.1), 0.1, 's');
  check('ocean: and scales as 1 / sin φ (at 60° S)', Math.abs(inertialPeriodHours(-60) - inertialPeriodHours(30) * Math.sin(Math.PI / 6) / Math.sin(Math.PI / 3)), 1e-12, 'h');
  const toward = [flowToward(1, 0), flowToward(0, 1), flowToward(-1, 0), flowToward(0, -1), flowToward(1, 1)];
  assert('ocean: currents are named by where they flow: E 90, N 0, W 270, S 180, NE 45',
    [90, 0, 270, 180, 45].every((x, i) => Math.abs(toward[i] - x) < 1e-9), toward.join(', '));
  // Mixed to 50 m at 20 °C, then cooling by 0.01 °C per metre: the 0.2 °C criterion is met at 70 m,
  // between levels at 60 and 80 m, where linear interpolation is exact.
  const depths = [0, 2, 5, 10, 20, 30, 50, 60, 80, 100, 150, 200];
  const profile = depths.map((z) => (z <= 50 ? 20 : 20 - 0.01 * (z - 50)));
  check('ocean: the mixed layer of a constructed profile', Math.abs(mixedLayerDepth(profile, depths).depth - 70), 1e-9, 'm');
  const uniform = mixedLayerDepth(depths.map(() => 15).map((x, k) => (k > 9 ? NaN : x)), depths);
  assert('ocean: a column mixed to its deepest value reaches the bottom', uniform.reachesBottom && uniform.depth === 100, JSON.stringify(uniform));
  assert('ocean: no mixed layer without a temperature at 10 m', Number.isNaN(mixedLayerDepth(depths.map(() => NaN), depths).depth));
  check('ocean: a profile between levels', Math.abs(atDepth(profile, depths, 70) - 19.8), 1e-12, '°C');
  // Bilinear interpolation reproduces f = a + b t + c z + d t z exactly, on unequal steps of time and depth.
  const times = [0, 3, 9, 12, 24].map((h) => h * 3600000);
  const f = (t, z) => 1 + 2e-8 * t + 0.003 * z - 4e-11 * t * z;
  const sample = columnSampler(times, depths, times.map((t) => depths.map((z) => f(t, z))));
  let worst = 0;
  for (let i = 0; i < 200; i++) {
    const t = Math.random() * times[4], z = Math.random() * 200;
    worst = Math.max(worst, Math.abs(sample(t, z) - f(t, z)));
  }
  check('ocean: the column sampler is exact for bilinear fields', worst, 1e-12);
  // Beneath the deepest value (the sea floor) the sampler gives nothing.
  const floor = columnSampler(times, depths, times.map(() => depths.map((z) => (z <= 100 ? 1 : null))));
  assert('ocean: nothing beneath the sea floor', Number.isNaN(floor(5 * 3600000, 140)) && floor(5 * 3600000, 90) === 1);
  check('ocean: a symmetric range for currents', Math.abs(symmetricRange([-0.33, 0.1, 0.25], { q: 1 }) - 0.4), 1e-12);
}

// ---- the motion of the tectonic plates ------------------------------------------------
{
  // NNR-MORVEL56 velocities (mm/yr, east and north) from the EarthScope (formerly UNAVCO)
  // Plate Motion Calculator, model morvel56_2011, NNR frame, height 0; retrieved 2026-10-22.
  const reference = [['PA', -155.5, 19.5, -61.24, 32.23], ['NA', -100, 40, -15.49, -7.67], ['EU', 10, 50, 17.60, 14.59], ['PS', 135, 20, -45.29, 16.58]];
  let worst = 0;
  for (const [code, lon, lat, e, n] of reference) {
    const v = velocityAt(omega(poleOf(code)), lon, lat);
    worst = Math.max(worst, Math.abs(v.east - e), Math.abs(v.north - n));
  }
  check('plates: velocities of four plates against the EarthScope calculator', worst, 0.006, 'mm/yr');
  assert('plates: NNR-MORVEL56 has 56 plates, and PB2002 codes map to them', NNR_MORVEL56.length === 56
    && morvelCode('AF') === 'nb' && morvelCode('NB') === 'NB' && poleOf('NB').name === 'North Bismarck' && poleOf('AF').name === 'Nubia');
  // A pole recovered from its angular velocity.
  const pole = poleOfOmega(omega({ lat: -63.58, lon: 114.70, rate: 0.651 }));
  check('plates: a pole from its angular velocity', Math.abs(pole.lat + 63.58) + Math.abs(pole.lon - 114.70) + Math.abs(pole.rate - 0.651), 1e-12);
  // Relative motion is the difference of angular velocities: v_B - v_A at every point.
  const wA = omega(poleOf('PA')), wB = omega(poleOf('NA'));
  const vr = velocityAt(subtract(wB, wA), -130, 45), va = velocityAt(wA, -130, 45), vb = velocityAt(wB, -130, 45);
  check('plates: relative motion from the difference of angular velocities', Math.hypot(vr.east - (vb.east - va.east), vr.north - (vb.north - va.north)), 1e-12, 'mm/yr');
  // Across a step heading north, the right-hand plate is to the east: moving east it opens the
  // boundary; moving south it slides past right-laterally (as seen from either side).
  const open = acrossStep({ east: 3, north: 0 }, 0), slide = acrossStep({ east: 0, north: -2 }, 0);
  assert('plates: opening and right-lateral sliding by the conventions of Bird (2003)',
    Math.abs(open.divergence - 3) < 1e-12 && Math.abs(open.rightLateral) < 1e-12 && Math.abs(slide.divergence) < 1e-12 && Math.abs(slide.rightLateral - 2) < 1e-12,
    JSON.stringify([open, slide]));
  check('plates: the azimuth of a step due east along the equator', Math.abs(azimuth(10, 0, 11, 0) - 90), 1e-9, '°');

  const meta = resolve(products, 'plates-pb2002.json');
  if (existsSync(meta)) {
    const m = JSON.parse(readFileSync(meta, 'utf8'));
    const b = readFileSync(resolve(products, m.file));
    const coords = new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength));
    const W = 1440, H = 720, raster = new Uint8Array(W * H);
    m.plates.forEach((pl, i) => fillRings(coords, pl.rings, W, H, i + 1, raster));
    const total = m.plates.reduce((sum, pl) => sum + featureAreaFraction(coords, pl.rings), 0);
    assert('plates: the 52 plates of PB2002 tile the sphere: every cell a plate, areas summing to 4π',
      raster.every((x) => x > 0) && Math.abs(total - 1) < 1e-5, `${raster.filter((x) => !x).length} cells unassigned; areas ${total}`);
    // NNR-MORVEL56 on PB2002's boundaries, against Bird's classes and his own (PB2002) relative velocities.
    const codes = m.plates.map((pl) => pl.code), classes = m.classes.map((c) => c.code);
    const tally = { open: [0, 0], close: [0, 0], lateral: [0, 0] };
    for (const [lo1, la1, lo2raw, la2, k, L, R, , birdOpen, birdLateral] of m.steps) {
      let lo2 = lo2raw;
      if (lo2 - lo1 > 180) lo2 -= 360; else if (lo1 - lo2 > 180) lo2 += 360;
      const v = velocityAt(subtract(omega(poleOf(codes[R])), omega(poleOf(codes[L]))), (lo1 + lo2) / 2, (la1 + la2) / 2);
      const x = acrossStep(v, azimuth(lo1, la1, lo2, la2));
      const c = classes[k];
      if (c === 'OSR' || c === 'CRB') { tally.open[1]++; if (x.divergence > 0) tally.open[0]++; }
      if (c === 'SUB') { tally.close[1]++; if (x.divergence < 0) tally.close[0]++; }
      if ((c === 'OTF' || c === 'CTF') && Math.abs(birdLateral) > 2) { tally.lateral[1]++; if (Math.sign(birdLateral) === Math.sign(x.rightLateral)) tally.lateral[0]++; }
    }
    // No net rotation: the plates' motions, weighted by area, sum to no rotation of the whole
    // surface, w_net = (3 / 8 pi) sum_P int_P (I - r r^T) w_P dA (Argus and Gordon 1991).
    const net = [0, 0, 0], RAD = Math.PI / 180;
    for (let j = 0; j < H; j++) {
      const lat = 90 - ((j + 0.5) * 180) / H, dA = Math.cos(lat * RAD) * (Math.PI / H) * ((2 * Math.PI) / W);
      for (let i = 0; i < W; i++) {
        const lon = -180 + ((i + 0.5) * 360) / W;
        const r = [Math.cos(lat * RAD) * Math.cos(lon * RAD), Math.cos(lat * RAD) * Math.sin(lon * RAD), Math.sin(lat * RAD)];
        const w = omega(poleOf(m.plates[raster[j * W + i] - 1].code));
        const rw = r[0] * w[0] + r[1] * w[1] + r[2] * w[2];
        for (let k = 0; k < 3; k++) net[k] += (w[k] - r[k] * rw) * dA;
      }
    }
    check('plates: NNR-MORVEL56 on PB2002 has no net rotation (against the Pacific’s 0.651°/Myr)', (Math.hypot(...net) * 3) / (8 * Math.PI) / RAD, 0.002, '°/Myr');
    const share = ([a2, n]) => a2 / n;
    assert('plates: spreading ridges and rifts open (at least 99 %)', share(tally.open) >= 0.99, `${tally.open[0]} of ${tally.open[1]}`);
    assert('plates: subduction zones close (at least 99 %)', share(tally.close) >= 0.99, `${tally.close[0]} of ${tally.close[1]}`);
    assert('plates: transforms slip in the sense Bird (2003) gives (at least 98 %)', share(tally.lateral) >= 0.98, `${tally.lateral[0]} of ${tally.lateral[1]}`);
  } else {
    console.log('skip  plates: plates-pb2002 missing');
  }
}

// ---- GPS stations ---------------------------------------------------------------------
{
  // ITRF2020's plate motion model against the EarthScope calculator (model itrf2020, which
  // leaves out the origin rate bias): the Pacific plate at 19.5 N, 155.5 W, retrieved 2026-10-22.
  const v = velocityAt(itrfOmega('PA'), -155.5, 19.5);
  check('gnss: the ITRF2020 plate motion model against the EarthScope calculator', Math.max(Math.abs(v.east + 62.24), Math.abs(v.north - 33.90)), 0.006, 'mm/yr');
  const sm = resolve(products, 'stations-ngl-midas.json');
  if (existsSync(sm)) {
    const d = JSON.parse(readFileSync(sm, 'utf8'));
    const f = Object.fromEntries(d.fields.map((k, i) => [k, i]));
    // Stations deep in four stable plates move by tenths of a millimetre a year relative to them.
    const boxes = { Australia: ['AU', 120, 150, -35, -20], 'North America': ['NA', -95, -75, 30, 40], Europe: ['EU', 5, 20, 45, 52], 'South America': ['SA', -55, -40, -25, -10] };
    const medians = {};
    for (const [name, [code, w, e, so, n]] of Object.entries(boxes)) {
      const res = [];
      for (const r of d.stations) {
        const lon = r[f.lon], lat = r[f.lat];
        if (lon < w || lon > e || lat < so || lat > n) continue;
        const pv = itrfVelocityAt(code, lon, lat);
        res.push(Math.hypot(r[f.ve] - pv.east, r[f.vn] - pv.north));
      }
      res.sort((x, y) => x - y);
      medians[name] = +res[res.length >> 1].toFixed(2);
    }
    assert('gnss: stations in stable plate interiors move less than 1 mm/yr relative to their plate (median)',
      Object.values(medians).every((m) => m < 1), JSON.stringify(medians));
    const sp = resolve(products, 'series/series-ngl-p123.json');
    if (existsSync(sp)) {
      const sr = JSON.parse(readFileSync(sp, 'utf8'));
      const years = (Date.parse(sr.days[sr.days.length - 1]) - Date.parse(sr.days[0])) / (365.25 * 86400000);
      const st = d.stations.find((r) => r[f.id] === 'P123');
      const de = Math.abs(sr.components[0].values.at(-1) / years - st[f.ve]), dn = Math.abs(sr.components[1].values.at(-1) / years - st[f.vn]);
      check('gnss: P123’s 20-year record moves at its MIDAS velocity', Math.max(de, dn), 0.5, 'mm/yr');
    }
  } else {
    console.log('skip  gnss: stations-ngl-midas missing');
  }
}

// ---- earthquake mechanisms --------------------------------------------------------
{
  // A double couple from strike, dip and rake (Aki and Richards 2002, eqs. 4.88, 4.89), in the
  // catalogue's elements (r up, t south, p east): M = n d^T + d n^T in (north, east, down).
  const R = Math.PI / 180;
  const doubleCouple = (strike, dip, rake) => {
    const f = strike * R, dl = dip * R, l = rake * R;
    const n = [-Math.sin(dl) * Math.sin(f), Math.sin(dl) * Math.cos(f), -Math.cos(dl)];
    const d = [Math.cos(l) * Math.cos(f) + Math.sin(l) * Math.cos(dl) * Math.sin(f), Math.cos(l) * Math.sin(f) - Math.sin(l) * Math.cos(dl) * Math.cos(f), -Math.sin(l) * Math.sin(dl)];
    const M = (i, j) => n[i] * d[j] + d[i] * n[j];
    return { mtt: M(0, 0), mtp: -M(0, 1), mrt: M(0, 2), mpp: M(1, 1), mrp: -M(1, 2), mrr: M(2, 2) };
  };
  // (45, 70, 45): T plunges 45°, N 42°, P 13°, none steep enough for a pure class.
  const cases = [['thrust', 20, 30, 90], ['normal', 120, 55, -90], ['strike-slip', 300, 89, 5], ['oblique', 45, 70, 45]];
  let worst = 0;
  const kinds = [];
  for (const [kind, strike, dip, rake] of cases) {
    const m = mechanism(doubleCouple(strike, dip, rake));
    const target = { strike, dip, rake };
    worst = Math.max(worst, Math.min(planeDifference(m.planes[0], target), planeDifference(m.planes[1], target)));
    kinds.push(m.kind === kind);
  }
  // acos of a cosine within 1e-16 of 1 resolves angles to about 1e-6 degrees.
  check('mechanisms: a fault recovered from its double couple', worst, 1e-5, '°');
  assert('mechanisms: thrust, normal, strike-slip and oblique by the plunges of the axes (Frohlich 1992)', kinds.every(Boolean), JSON.stringify(kinds));
  // Eigenvalues of a double couple: 1, 0, -1; T in a quadrant of compression, P of dilatation.
  const dc = doubleCouple(20, 30, 90), e = eigenSymmetric(tensorNED(dc));
  check('mechanisms: the eigenvalues of a unit double couple are 1, 0, -1', Math.abs(e[0].value - 1) + Math.abs(e[1].value) + Math.abs(e[2].value + 1), 1e-12);
  assert('mechanisms: P waves push along T and pull along P', radiation(tensorNED(dc), e[0].vector) > 0 && radiation(tensorNED(dc), e[2].vector) < 0);
  let round = 0;
  for (let i = 0; i < 50; i++) {
    const x = random() * 1.4 - 0.7, y = random() * 1.4 - 0.7;
    const [x2, y2] = toEqualArea(fromEqualArea(x, y));
    round = Math.max(round, Math.hypot(x2 - x, y2 - y));
  }
  check('mechanisms: the equal-area projection, there and back', round, 1e-12);
  const gm = resolve(products, 'events-gcmt.json');
  if (existsSync(gm)) {
    const g = JSON.parse(readFileSync(gm, 'utf8'));
    const f = Object.fromEntries(g.fields.map((k, i) => [k, i]));
    const diffs = [];
    for (const r of g.events) {
      const m = mechanism({ mrr: r[f.mrr], mtt: r[f.mtt], mpp: r[f.mpp], mrt: r[f.mrt], mrp: r[f.mrp], mtp: r[f.mtp] });
      const theirs = [{ strike: r[f.strike1], dip: r[f.dip1] }, { strike: r[f.strike2], dip: r[f.dip2] }];
      diffs.push(Math.min(Math.max(planeDifference(m.planes[0], theirs[0]), planeDifference(m.planes[1], theirs[1])),
        Math.max(planeDifference(m.planes[0], theirs[1]), planeDifference(m.planes[1], theirs[0]))));
    }
    diffs.sort((x, y) => x - y);
    check(`mechanisms: nodal planes of all ${g.count.toLocaleString('en-GB')} Global CMT events against the catalogue's (99th percentile)`, diffs[Math.floor(0.99 * (diffs.length - 1))], 1.0, '°');
    const tohoku = g.events.find((r) => r[f.name] === 'M201103110546A');
    const mt = mechanism({ mrr: tohoku[f.mrr], mtt: tohoku[f.mtt], mpp: tohoku[f.mpp], mrt: tohoku[f.mrt], mrp: tohoku[f.mrp], mtp: tohoku[f.mtp] });
    assert('mechanisms: Tohoku 2011 (Mw 9.1) is a thrust on a plane dipping 10° west-north-west',
      mt.kind === 'thrust' && Math.abs(tohoku[f.mw] - 9.08) < 0.01 && planeDifference(mt.planes[0], { strike: 203, dip: 10 }) < 1, JSON.stringify(mt.planes[0]));
  } else {
    console.log('skip  mechanisms: events-gcmt missing');
  }
}

// ---- the crust, and sections through the Earth ----------------------------------------
{
  // A constructed oceanic column: 4 km of water, 0.5 km of sediment, 7 km of crust (no ice,
  // no middle or lower sediments, no middle crust).
  const tops = [0, -4, -4, -4.5, -4.5, -4.5, -6, -6, -11.5];
  const vp = [1.5, 0, 2, 0, 0, 5, 0, 7, 8.1];
  const c = crustColumn((q, k) => (q === 'top' ? tops[k] : q === 'vp' ? vp[k] : 0));
  assert('crust: a constructed oceanic column, layer by layer',
    c.layers.map((l) => l.name).join(',') === 'water,upper sediments,upper crystalline crust,lower crystalline crust'
    && Math.abs(c.moho - 11.5) < 1e-12 && Math.abs(c.crust - 7.5) < 1e-12 && Math.abs(c.sediments - 0.5) < 1e-12 && c.oceanic,
    JSON.stringify({ layers: c.layers.map((l) => [l.name, l.thickness]), moho: c.moho, crust: c.crust }));
  const g = { lon0: -179.5, dlon: 1, nlon: 360, lat0: 89.5, dlat: 1, nlat: 180 };
  assert('crust: the 1° cell containing a place (its corners and the antimeridian)',
    cellIndex(g, -180, 90) === 0 && cellIndex(g, 179.99, -89.99) === 360 * 180 - 1 && cellIndex(g, 0.5, 0.5) === 89 * 360 + 180 && cellIndex(g, 180, 0) === cellIndex(g, -180, 0));
  // Sections: a place 80 km to the right of a line, 400 km along it, is found there; one behind
  // the start has a negative distance along; a half-circumference closes the sphere.
  const lineAz = 100, [x1, y1] = travel(140, 38, lineAz, 400);
  // The line's direction there: the reverse of the azimuth back to its start.
  const local = azimuth(x1, y1, 140, 38) + 180;
  const [x2, y2] = travel(x1, y1, local + 90, 80);
  const r1 = relativeTo(140, 38, lineAz, x2, y2);
  const [x3, y3] = travel(140, 38, lineAz + 180, 250);
  const r2 = relativeTo(140, 38, lineAz, x3, y3);
  check('sections: along- and across-track distances of a place beside the line', Math.abs(r1.along - 400) + Math.abs(r1.across - 80) + Math.abs(r2.along + 250) + Math.abs(r2.across), 1e-6, 'km');
  const ends = sectionLine(0, 0, 90, Math.PI * R_KM / 2, Math.PI * R_KM / 2);
  assert('sections: a half circumference from 90° W to 90° E along the equator', ends.length === 3
    && Math.abs(ends[0].lon + 90) < 1e-9 && Math.abs(ends[2].lon - 90) < 1e-9 && Math.abs(ends[1].lon) < 1e-9);
  const cm = resolve(products, 'crust-crust1.json');
  if (existsSync(cm)) {
    const meta = JSON.parse(readFileSync(cm, 'utf8'));
    const buf = readFileSync(resolve(products, meta.file));
    const values = new Int16Array(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
    const q = Object.fromEntries(meta.quantities.map((x) => [x.name, x]));
    const at = (lon, lat) => { const cell = cellIndex(meta.grid, lon, lat); return crustColumn((name, k) => values[q[name].offsets[k] / 2 + cell] * q[name].scale); };
    const tibet = at(91.1, 29.65), pacific = at(-150, 10), shield = at(-100, 55);
    assert('crust: Tibet over 60 km thick; the central Pacific oceanic, its Moho 10 to 14 km below sea level; the Canadian Shield 35 to 45 km',
      tibet.crust > 60 && pacific.oceanic && pacific.moho > 10 && pacific.moho < 14 && shield.crust > 35 && shield.crust < 45 && !shield.oceanic,
      JSON.stringify({ tibet: tibet.crust, pacific: pacific.moho, shield: shield.crust }));
    assert('crust: P waves quicken to 7.6-8.4 km/s beneath the Moho everywhere checked (Pn)',
      [tibet, pacific, shield].every((x) => x.mantle.vp > 7.6 && x.mantle.vp < 8.4), [tibet, pacific, shield].map((x) => x.mantle.vp).join(', '));
  } else {
    console.log('skip  crust: crust-crust1 missing');
  }
}

// ---- water in shaded relief ---------------------------------------------------------
{
  // A 40 x 30 relief of random greys (60-200), with a lake of the water grey (88)
  // and single pixels of that grey scattered as shadows.
  const w = 40, h = 30, grey = new Uint8Array(w * h), rnd = seeded(5);
  for (let i = 0; i < w * h; i++) { grey[i] = 60 + Math.floor(140 * rnd()); if (grey[i] >= 87 && grey[i] <= 89) grey[i] = 120; }
  for (let y = 10; y < 20; y++) for (let x = 5; x < 18; x++) grey[y * w + x] = 88;
  const speckles = [[30, 4], [33, 25], [25, 12]];
  for (const [x, y] of speckles) grey[y * w + x] = 88;
  const m = waterMask(grey, w, h, 88);
  let lake = 0;
  for (let y = 10; y < 20; y++) for (let x = 5; x < 18; x++) lake += m[y * w + x];
  assert('water mask: the whole lake, to its edge, is water', lake === 130, `${lake} of 130 pixels`);
  assert('water mask: lone pixels of the water grey are not', speckles.every(([x, y]) => !m[y * w + x]));
}

// ---- climate: the Köppen–Geiger classification ------------------------------------
{
  // Climates constructed to meet one rule each (Beck et al. 2018, Table 2):
  // temperature a + b cos(2 pi (m - 6.5) / 12), warmest in July; precipitation by month.
  const season = (a, b) => [...Array(12).keys()].map((m) => a + b * Math.cos((2 * Math.PI * (m - 6.5)) / 12));
  const flat = (v) => Array(12).fill(v);
  // Precipitation in summer (April-September) and in winter (October-March), month by month.
  const halves = (summer, winter) => [...Array(12).keys()].map((m) => (m >= 3 && m <= 8 ? summer[m - 3] : winter[(m + 3) % 12]));
  const cases = [
    ['Af', season(27, 1), flat(250), 'rainforest: every month at least 60 mm'],
    ['Am', season(27, 1), [40, ...Array(11).fill(300)], 'monsoon: a dry month, made up by the wet ones'],
    ['Aw', season(26, 2), halves(flat(200), flat(10)), 'savannah: a long dry winter'],
    ['BWh', season(28, 6), flat(5), 'hot desert'],
    ['BSk', season(10, 12), flat(25), 'cold steppe'],
    ['BSk', season(1, 6), flat(8), 'B before E: a dry climate with no month above 10 C is arid, not tundra'],
    ['Csa', season(17, 8), halves(flat(5), flat(100)), 'Mediterranean: dry, hot summer'],
    ['Cfb', season(11, 7), flat(60), 'oceanic'],
    ['Cwa', season(20, 7), halves(flat(250), flat(10)), 'dry winter, hot summer'],
    ['Dfb', season(5, 14), flat(50), 'humid continental'],
    ['Dwd', season(-12, 30), halves(flat(60), flat(2)), 'dry, very cold winter'],
    ['ET', season(-2, 8), flat(30), 'tundra'],
    ['EF', season(-30, 10), flat(10), 'ice cap'],
    ['Csa', season(15, 8), halves([3, 200, 3, 3, 3, 3], [15, 60, 60, 60, 60, 60]), 'both s and w met, more rain in winter: s'],
    ['Cwa', season(15, 8), halves([3, 400, 3, 3, 3, 3], [15, 60, 60, 60, 60, 60]), 'both s and w met, more rain in summer: w'],
  ];
  const wrong = cases.filter(([code, T, P]) => koppen(T, P) !== code);
  assert('köppen: fifteen constructed climates take their classes', !wrong.length,
    wrong.map(([code, T, P, what]) => `${what}: ${koppen(T, P)} for ${code}`).join('; ') || cases.map((c) => c[0]).join(' '));
  assert('köppen: a missing month gives no class', koppen(season(10, 5), [NaN, ...flat(50).slice(1)]) === null);
  const q = koppenQuantities(season(-10, 20), flat(30));
  assert('köppen: summer is the warmer half-year', q.summer === 'April–September' && koppenQuantities(season(-10, -20), flat(30)).summer === 'October–March');

  // The classes are invariant under the half-year shift of the months (Z/2), for
  // climates whose half-years differ in temperature: summer and winter are
  // defined by temperature, not by the calendar.
  const shift = (a) => a.map((_, i) => a[(i + 6) % 12]);
  let broken = 0, n = 0;
  const rnd = seeded(7);
  for (let k = 0; k < 4000; k++) {
    const T = season(-25 + 55 * rnd(), (rnd() - 0.5) * 50).map((t) => t + 2 * (rnd() - 0.5));
    const P = [...Array(12)].map(() => 300 * rnd() ** 2);
    if (Math.abs(T.slice(3, 9).reduce((a, b) => a + b) - shift(T).slice(3, 9).reduce((a, b) => a + b)) < 1e-9) continue;
    n++;
    if (koppen(T, P) !== koppen(shift(T), shift(P))) broken++;
  }
  assert('köppen: classes are invariant under the half-year shift', broken === 0, `${n - broken} of ${n} random climates`);
  // More: within each half-year, the temperatures and the precipitations may be
  // permuted independently, and the halves exchanged: the class is a function on
  // the orbits of the group (S6 x S6) wr Z/2, of order 2 (6!)^4.
  const HALF = [[3, 4, 5, 6, 7, 8], [9, 10, 11, 0, 1, 2]];
  const shuffle = (r, xs) => { const a = [...xs]; for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(r() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; } return a; };
  const act = (r, x, swap) => {
    const out = new Array(12);
    HALF.forEach((half, h) => { const to = HALF[swap ? 1 - h : h]; shuffle(r, half).forEach((m, k) => { out[to[k]] = x[m]; }); });
    return out;
  };
  let moved = 0, kept = 0;
  const rnd3 = seeded(13);
  for (let k = 0; k < 4000; k++) {
    const T = season(-25 + 55 * rnd3(), (rnd3() - 0.5) * 50).map((t) => t + 3 * (rnd3() - 0.5));
    const P = [...Array(12)].map(() => 300 * rnd3() ** 2);
    const swap = rnd3() < 0.5;
    const T2 = act(rnd3, T, swap), P2 = act(rnd3, P, swap);
    moved++;
    if (koppen(T, P) === koppen(T2, P2)) kept++;
  }
  assert('köppen: classes are invariant under the wreath product (S6 x S6) wr Z/2', kept === moved, `${kept} of ${moved} random climates and group elements`);
  // Random climates with a wetter or drier half-year meet every class.
  const classes = new Map();
  const rnd2 = seeded(11);
  for (let k = 0; k < 60000; k++) {
    const T = season(-40 + 75 * rnd2(), 35 * rnd2());
    const wet = [300 * rnd2() ** 2, 300 * rnd2() ** 2];
    const c = koppen(T, [...Array(12).keys()].map((m) => wet[m >= 3 && m <= 8 ? 0 : 1] * rnd2()));
    classes.set(c, (classes.get(c) || 0) + 1);
  }
  const rarest = [...classes].sort((x, y) => x[1] - y[1])[0];
  assert('köppen: every one of the 30 classes is reachable', KOPPEN_CODES.every((c) => classes.has(c)), `${classes.size} classes met; the rarest, ${rarest[0]}, ${rarest[1]} times in 60 000`);

  // Dry months after Gaussen (P <= 2 T), and the Walter-Lieth precipitation scale.
  const dry = dryMonths(season(17, 8), halves(flat(5), flat(100)));
  assert('climate: the Mediterranean summer is dry, its winter wet', dry.slice(4, 9).every(Boolean) && !dry[0] && !dry[11]);
  check('climate: 100 mm stands at 50 degrees, 300 mm at 60', Math.abs(precipitationScale(100) - 50) + Math.abs(precipitationScale(300) - 60), 1e-12);
  const svg = climateDiagram({ T: season(-12, 30), P: halves(flat(60), flat(2)), Tmin: season(-18, 30), lat: 62 });
  const frost = (svg.match(/class="frost"/g) || []).length;
  assert('climate: the diagram marks the frost months', frost === 8, `${frost} months with mean daily minimum below 0 C`);
}

// ---- climate: the published maps and our classification -----------------------------
{
  const km = resolve(products, 'classes-koppen-geiger-1991-2020.json');
  const tm = resolve(products, 'fields-terraclimate-normals-medium.json');
  if (existsSync(km) && existsSync(tm)) {
    const meta = JSON.parse(readFileSync(km, 'utf8'));
    const beck = new ClassRaster(meta, readFileSync(resolve(products, meta.file)));
    assert('classes: the legend of the maps lists the classes in our order',
      meta.classes.map((c) => c.code).join(' ') === KOPPEN_CODES.join(' '));
    // Places of known climate on Beck et al.'s 1991-2020 map.
    const places = [['Paris', 2.35, 48.85, 'Cfb'], ['Cairo', 31.24, 30.04, 'BWh'], ['Singapore', 103.8, 1.35, 'Af'],
      ['Moscow', 37.6, 55.75, 'Dfb'], ['Rome', 12.5, 41.9, 'Csa'], ['Mumbai', 72.88, 19.08, 'Am'], ['Tokyo', 139.7, 35.68, 'Cfa']];
    const off = places.filter(([, lon, lat, code]) => beck.classAt(lon, lat)?.code !== code);
    assert('classes: seven cities in their classes on the 1991-2020 map', !off.length, off.map((p) => `${p[0]} ${beck.classAt(p[1], p[2])?.code}`).join(', '));
    const m = JSON.parse(readFileSync(tm, 'utf8'));
    const bin = readFileSync(resolve(products, m.file));
    const f = new FieldProduct(m, bin.buffer.slice(bin.byteOffset, bin.byteOffset + bin.byteLength));
    assert('normals: an instant is shown by its month', f.bracket(Date.UTC(2026, 6, 31, 23)).i0 === 6 && f.bracket(Date.UTC(1999, 0, 1)).i0 === 0);
    const months = [...Array(12).keys()];
    const ours = koppenGrid(months.map((i) => f.slice('tmax', i)), months.map((i) => f.slice('tmin', i)), months.map((i) => f.slice('ppt', i)));
    let both = 0, same = 0, group = 0;
    const g = f.grid;
    for (let r = 0; r < g.nlat; r++) {
      for (let c = 0; c < g.nlon; c++) {
        const a = ours[r * g.nlon + c], b = beck.valueAt(g.lon0 + c * g.dlon, g.lat0 - r * g.dlat);
        if (!a || !b) continue;
        both++;
        if (a === b) same++;
        if (KOPPEN_CODES[a - 1][0] === KOPPEN_CODES[b - 1][0]) group++;
      }
    }
    // The two rest on different climatologies (TerraClimate sampled every 1/3 degree;
    // Beck et al.'s ensemble at 1 km): agreement, not identity, is expected.
    assert('classes: ours agree with Beck et al. in at least 85 % of land cells', same / both >= 0.85, `${(100 * same / both).toFixed(1)} % of ${both}`);
    assert('classes: and in their main group in at least 93 %', group / both >= 0.93, `${(100 * group / both).toFixed(1)} %`);
    const at = (v) => months.map((i) => f.sample(v, Date.UTC(2005, i, 15), 2.35, 48.85));
    const T = meanTemperature(at('tmax'), at('tmin'));
    const mean = T.reduce((a, b) => a + b) / 12;
    const P = at('ppt').reduce((a, b) => a + b);
    // A plausibility check: the Paris basin has an oceanic climate, with a mean of
    // 11-13 C and 550-750 mm of precipitation a year (the cell is not the city).
    assert('normals: the year at Paris is mild and moderately wet', mean > 11 && mean < 13.5 && P > 550 && P < 750, `${mean.toFixed(2)} °C, ${P.toFixed(0)} mm`);
  } else {
    console.log('skip  classes: classes-koppen-geiger-1991-2020 or fields-terraclimate-normals-medium missing');
  }
}

console.log(`\n${checks - failures} of ${checks} checks passed.`);
process.exit(failures ? 1 : 0);
