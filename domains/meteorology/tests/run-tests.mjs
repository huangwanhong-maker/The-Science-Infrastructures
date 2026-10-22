// Checks of the meteorology domain: the wind; the barotropic vorticity model
// against the Rossby-Haurwitz wave; the shallow-water model against Williamson
// et al. (1992) and Lamb (1932); tropical cyclones found in constructed fields.
// Run from the repository root:  node domains/meteorology/tests/run-tests.mjs

import { suite, seeded } from '../../../development/tools/checks.mjs';
import {
  windFromDirection, BarotropicModel, GRAVITY, EARTH_RADIUS as RADIUS, OMEGA, ShallowWaterModel,
  lanczosWeights, detect, stitch, arcDeg, positionAt,
} from '../js/index.js';
import { SphericalHarmonics } from '../../../foundations/physis/js/index.js';

const { check, assert, done } = suite();

// ---- the wind ---------------------------------------------------------------------
{
  assert('wind from the west is 270', Math.abs(windFromDirection(5, 0) - 270) < 1e-12);
  assert('wind from the north is 0', Math.abs(windFromDirection(0, -5)) < 1e-12);
}

// ---- the barotropic model -------------------------------------------------------------
{
  const sh = new SphericalHarmonics(42);
  // The solid-body rotation u = u0 cos(phi), stream function psi = -a u0 mu.
  const u0 = 20, psi = sh.zeros();
  psi.re[sh.index(0, 1)] = (-RADIUS * u0) / Math.sqrt(1.5);
  const rnd = seeded(42);

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
  const exact = sh.scale(rh(86400), (n) => sh.laplacianEigen(n, RADIUS));
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

done();
