// Checks of the shared visualization layer's provisional JavaScript: isolines and
// extrema, particles carried by a flow, the units in which quantities are shown,
// and the tiles of a map service.
// Run from the repository root:  node visualization/tests/run-tests.mjs

import { suite, seeded } from '../../development/tools/checks.mjs';
import {
  isolines, levels, extrema, ParticleSystem, uniformSpawn, displayUnit, toShown, toStored, format, tilesFor,
} from '../js/index.js';
import { DEG, EARTH_RADIUS_M } from '../../domains/geodesy/js/index.js';

const { check, assert, done } = suite();

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

// ---- units ------------------------------------------------------------------
{
  const k = displayUnit('K');
  check('units: 300 K in degrees Celsius', toShown(300, k) - 26.85, 1e-9);
  check('units: back to kelvin', toStored(26.85, k) - 300, 1e-9);
  check('units: 1e-3 kg m-2 s-1 is 3.6 mm/h', toShown(1e-3, displayUnit('kg m-2 s-1')) - 3.6, 1e-12);
  check('units: 10 m/s in knots', toShown(10, displayUnit('m s-1', 'kn')) - 19.4384, 1e-4);
  assert('units: numbers for reading', format(124516650) === '125 M' && format(3.14159) === '3.14' && format(1013.2) === '1013');
}

// ---- tiles ------------------------------------------------------------------
// Tiles of the Trek map services, for detail at close range: the tiles chosen
// cover the region asked for, at a level no coarser than needed (unless the
// count of tiles forbids it), and never more than the limit.
{
  const live = { tile: 256, levels: 8 };
  let uncovered = 0, tooMany = 0, coarse = 0;
  let seed = 9;
  const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
  for (let k = 0; k < 2000; k++) {
    const width = 1 + 200 * rnd() ** 2, height = Math.min(170, 1 + 120 * rnd() ** 2);
    const west = -180 + (360 - Math.min(width, 359)) * rnd(), south = -90 + (180 - height) * rnd();
    const box = [west, south, Math.min(width, 180 - west), height];
    const ppd = 2 + 600 * rnd() ** 3;
    const t = tilesFor(box, ppd, live);
    const [w, s, cw, ch] = t.box;
    if (w > box[0] + 1e-9 || s > box[1] + 1e-9 || w + cw < box[0] + box[2] - 1e-9 || s + ch < box[1] + box[3] - 1e-9) uncovered++;
    if ((t.c1 - t.c0 + 1) * (t.r1 - t.r0 + 1) > 48 && t.z > 0) tooMany++;
    const needed = Math.min(live.levels, Math.max(0, Math.ceil(Math.log2((ppd * 360) / 512))));
    if (t.z < needed && tilesFor(box, ppd, live, Infinity).z === t.z) coarse++;
    if (t.z === needed && needed < live.levels && (2 ** (t.z + 1) * 256) / 360 < ppd * 0.999) coarse++;
  }
  check('tiles cover the region asked for (2000 regions)', uncovered, 0, 'regions');
  check('tiles: never more than 48', tooMany, 0, 'regions');
  check('tiles: the level is the one needed, or limited by the count', coarse, 0, 'regions');
  const t = tilesFor([-10, -5, 20, 10], 100, live);
  check('tiles: 20 x 10 degrees at 100 pixels per degree is level 6, 8 x 4 tiles',
    Math.abs(t.z - 6) + Math.abs(t.c1 - t.c0 + 1 - 8) + Math.abs(t.r1 - t.r0 + 1 - 4), 0, '');
}

done();
