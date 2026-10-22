// Checks of the oceanography domain: the ocean column beneath a place (inertial
// period, mixed layer, sampling in time and depth) and the propagation of waves.
// Run from the repository root:  node domains/oceanography/tests/run-tests.mjs

import { suite } from '../../../development/tools/checks.mjs';
import {
  inertialPeriodHours, flowToward, mixedLayerDepth, columnSampler, symmetricRange, atDepth,
  deepWaterGroupSpeed, waveEnergyVelocity,
} from '../js/index.js';

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

done();
