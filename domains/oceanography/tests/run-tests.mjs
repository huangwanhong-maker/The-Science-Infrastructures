// Checks of the oceanography domain: the ocean column beneath a place (inertial
// period, mixed layer, sampling in time and depth), the propagation of waves,
// and tsunami travel times over constructed oceans.
// Run from the repository root:  node domains/oceanography/tests/run-tests.mjs

import { suite } from '../../../development/tools/checks.mjs';
import {
  inertialPeriodHours, flowToward, mixedLayerDepth, columnSampler, symmetricRange, atDepth,
  deepWaterGroupSpeed, waveEnergyVelocity,
} from '../js/index.js';
import { TravelTimes, travelTimes, longWaveSpeed } from '../js/index.js';
import { distanceM } from '../../geodesy/js/index.js';

const { check, assert, done } = suite();

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

// ---- waves ------------------------------------------------------------------------------
{
  check('waves: deep-water group speed g T / 4 pi, T = 10 s', deepWaterGroupSpeed(10) - 7.80388, 1e-5, 'm/s');
  const [u, v] = waveEnergyVelocity(270, 10);           // waves from the west travel east
  check('waves: waves from the west travel east', Math.hypot(u - deepWaterGroupSpeed(10), v), 1e-12);
}

// ---- tsunami travel times ------------------------------------------------------------------
{
  check('tsunami: the long wave over 4000 m travels at √(g h) = 198.06 m/s', longWaveSpeed(4000) - 198.06, 0.005, 'm/s');
  // An ocean of uniform depth on a 1/2 degree grid: the first arrival is the great-circle
  // distance over √(g h). The bounds of the stencils (tsunami.js) are 2.1 % and 4.2 %.
  const step = 0.5, nlon = 720, nlat = 360;
  const grid = { lon0: -180 + step / 2, dlon: step, nlon, lat0: 90 - step / 2, dlat: step, nlat, periodic: true };
  const at = (i) => [grid.lon0 + (i % nlon) * step, grid.lat0 - Math.floor(i / nlon) * step];
  const c = longWaveSpeed(4000);
  const excess = (times, src, near = 1000e3, far = 9000e3) => {
    let worst = 0, sum = 0, n = 0;
    for (let i = 0; i < times.length; i++) {
      const d = distanceM(src[0], src[1], ...at(i));
      if (d < near || d > far) continue;
      const e = (times[i] * c) / d - 1;
      if (Math.abs(e) > Math.abs(worst)) worst = e;
      sum += e; n++;
    }
    return { worst, mean: sum / n };
  };
  const uniform = new Float32Array(nlon * nlat).fill(4000);
  for (const src of [[-150.3, 0.2], [140.1, 47.3], [10.4, -62.1]]) {
    const e32 = excess(travelTimes(uniform, grid, src, { neighbours: 32 }).times, src);
    const e16 = excess(travelTimes(uniform, grid, src, { neighbours: 16 }).times, src);
    check(`tsunami: uniform ocean, source at ${src[1]}°, 1000-9000 km: the 32-neighbour stencil within 3 % (mean ${(100 * e32.mean).toFixed(2)} %)`, e32.worst, 0.03);
    check(`tsunami: and the 16-neighbour stencil within 4.5 % (mean ${(100 * e16.mean).toFixed(2)} %)`, e16.worst, 0.045);
  }
  // Land is impassable. A wall of land along a meridian, from pole to pole but for a gap at
  // 30° N: behind the wall the wave arrives by the gap, as the great circles through it say.
  const walled = uniform.slice(), wallCol = Math.round((0 - grid.lon0) / step);
  for (let r = 0; r < nlat; r++) for (const dc of [0, 1]) walled[r * nlon + wallCol + dc] = Math.abs(grid.lat0 - r * step - 30) < 1 ? 4000 : NaN;
  // ... and on the far side of the globe, a second wall, so that the wave cannot go round.
  const backCol = Math.round((179.75 - grid.lon0) / step);
  for (let r = 0; r < nlat; r++) walled[r * nlon + backCol] = NaN;
  const src = [-20, 0], target = [20, 0], gap = [0.25, 30];
  const { times } = travelTimes(walled, grid, src, { neighbours: 32 });
  const ti = Math.round((grid.lat0 - target[1]) / step) * nlon + Math.round((target[0] - grid.lon0) / step);
  const detour = (distanceM(...src, ...gap) + distanceM(...gap, ...at(ti))) / c;
  check('tsunami: behind a wall, the wave comes round by its gap (relative to the path through the gap)', times[ti] / detour - 1, 0.03);
  assert('tsunami: and later than the straight path it cannot take', times[ti] > (1.2 * distanceM(...src, ...at(ti))) / c);
  assert('tsunami: land has no arrival time', Number.isNaN(times[Math.round((grid.lat0 - 0) / step) * nlon + wallCol]));
  // A closed sea whose shore is a diamond of single land cells joined only at their corners:
  // no step may cut a corner, so nothing outside is reached.
  const closed = uniform.slice(), r0 = Math.round((grid.lat0 - 10) / step), c0 = Math.round((100 - grid.lon0) / step), half = 12;
  for (let k = 0; k <= half; k++) for (const [dr, dc] of [[k - half, k], [k, half - k], [half - k, -k], [-k, k - half]]) closed[(r0 + dr) * nlon + c0 + dc] = NaN;
  const inside = new TravelTimes(closed, grid, [100, 10], { neighbours: 32 });
  const t = inside.solve();
  let outside = 0;
  for (let i = 0; i < t.length; i++) if (t[i] === t[i]) { const [lo, la] = at(i); if (Math.abs(lo - 100) + Math.abs(la - 10) > half * step) outside++; }
  assert('tsunami: a sea closed by land touching only at corners keeps the wave in', outside === 0 && inside.settled > 100, `${outside} points outside reached, ${inside.settled} inside`);
  // A source on an island of 5 x 5 cells, at its centre, moves to the nearest sea within
  // reach: three cells east or west, 1.5° of longitude at its latitude; or, beyond reach, nowhere.
  const island = uniform.slice(), [ilon, ilat] = at(r0 * nlon + c0);
  for (let r = r0 - 2; r <= r0 + 2; r++) for (let k = c0 - 2; k <= c0 + 2; k++) island[r * nlon + k] = NaN;
  const moved = new TravelTimes(island, grid, [ilon, ilat], { snapKm: 200 }).source;
  assert('tsunami: a source on an island starts from the nearest sea', moved && Math.abs(moved.movedKm - distanceM(ilon, ilat, ilon + 1.5, ilat) / 1000) < 1e-6 && moved.lat === ilat,
    moved ? `${moved.movedKm.toFixed(2)} km` : 'none');
  assert('tsunami: and none from an island wider than the reach', new TravelTimes(island, grid, [ilon, ilat], { snapKm: 100 }).source === null);
}

done();
