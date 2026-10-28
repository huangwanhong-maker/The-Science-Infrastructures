// Checks of the meteorology domain: the wind; the barotropic vorticity model
// against the Rossby-Haurwitz wave; the shallow-water model against Williamson
// et al. (1992) and Lamb (1932); tropical cyclones found in constructed fields;
// the Fire Weather Index System against its standard test data; the wet-bulb
// temperature and the heat index; the standard atmosphere and settling parcels.
// And its C++ (the U.S. Standard Atmosphere 1976 to 86 km, atmosphere/standard.hpp), built
// with CMake and run on calx's verification; the JavaScript standard atmosphere of
// settling.js is held to the same references (tests/fixtures/atmosphere-reference.txt).
// Run from the repository root:  node domains/meteorology/tests/run-tests.mjs

import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, seeded, cppChecks } from '../../../development/tools/checks.mjs';
import {
  windFromDirection, BarotropicModel, GRAVITY, EARTH_RADIUS as RADIUS, OMEGA, ShallowWaterModel,
  lanczosWeights, detect, stitch, arcDeg, positionAt,
} from '../js/index.js';
import { SphericalHarmonics } from '../../../foundations/physis/js/index.js';
import {
  STANDARD_START, fireWeatherDay, fineFuelMoistureCode, dcDayLength, dmcDayLength, localNoon, accumulated, forecastNoonDays,
  wetBulbStull, withinStullRange, heatIndexF, sampledMaximum, forecastLocalDays, localDay,
  standardAtmosphere, pressureAltitude, stokesSettling, windAtHeight, settlingPath,
} from '../js/index.js';

const { check, assert, skip, lines, done } = suite();
const HERE = dirname(fileURLToPath(import.meta.url));

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

// ---- fire weather ----------------------------------------------------------------------
{
  // The standard test data of Van Wagner and Pickett (1985), Forestry Technical Report 33
  // (also the first rows of the cffdrs package's test_fwi): from FFMC 85, DMC 6, DC 15 on
  // April 13, at 46° N; weather (T °C, RH %, wind km/h, rain mm) and the published
  // FFMC, DMC, DC, ISI, BUI, FWI, rounded to 0.1.
  const rows = [
    [17, 42, 25, 0, 87.7, 8.5, 19.0, 10.9, 8.5, 10.1],
    [20, 21, 25, 2.4, 86.2, 10.4, 23.6, 8.8, 10.4, 9.3],
    [8.5, 40, 17, 0, 87.0, 11.8, 26.1, 6.5, 11.7, 7.6],
    [6.5, 25, 6, 0, 88.8, 13.2, 28.2, 4.9, 13.1, 6.2],
    [13, 34, 24, 0, 89.1, 15.4, 31.5, 12.6, 15.3, 14.8],
  ];
  let codes = STANDARD_START, worst = 0, first = null;
  for (const [T, H, W, rain, ...published] of rows) {
    codes = fireWeatherDay(codes, { T, H, W, rain }, 4, 46);
    const ours = [codes.ffmc, codes.dmc, codes.dc, codes.isi, codes.bui, codes.fwi];
    first = first || ours;
    worst = Math.max(worst, ...ours.map((x, k) => Math.abs(x - published[k])));
  }
  check('fire weather: the first day of the standard test data (FFMC 87.7, DMC 8.5, DC 19.0, ISI 10.9, BUI 8.5, FWI 10.1)',
    Math.max(...first.map((x, k) => Math.abs(x - rows[0][4 + k]))), 0.05);
  check('fire weather: all six codes on five days, rain on the second, to the published rounding', worst, 0.05);
  // Equations (1) and (10), moisture content from FFMC and back: (10) inverts exactly to
  // m = 147.2 (101.05 - F) / (59.5 + F), since 59.5 × 250 / 147.2 = 101.05, which (1) rounds to 101.
  const m = (F) => (147.2 * (101 - F)) / (59.5 + F), F = (mc) => (59.5 * (250 - mc)) / (147.2 + mc);
  check('fire weather: the FFMC scale and moisture content are inverse to 0.05 (eq. 1 and 10)', Math.max(...[0, 30, 85, 101].map((x) => Math.abs(F(m(x)) - x))), 0.055);
  // Without rain or wind, the fine fuel approaches its equilibrium moisture: a day of 30 °C and 10 %
  // raises a moist FFMC toward its ceiling, a day of 5 °C and 95 % lowers a dry one.
  assert('fire weather: dry air raises the FFMC, humid air lowers it', fineFuelMoistureCode(70, 30, 10, 0, 0) > 85 && fineFuelMoistureCode(92, 5, 95, 0, 0) < 92);
  // The southern band's DC adjustment is the northern one six months on; the DMC's day length is 9 h at the equator.
  assert('fire weather: day-length factors of the southern hemisphere are the northern ones six months on',
    [...Array(12).keys()].every((k) => dcDayLength(k + 1, -45) === dcDayLength(((k + 6) % 12) + 1, 45)) && dmcDayLength(7, 0) === 9 && dmcDayLength(4, 46) === 12.8);
  // Local noon, and the rain of a day from a rate sampled every 6 hours.
  const H = 3600000, t0 = Date.UTC(2026, 9, 20, 18), times = [...Array(24).keys()].map((k) => t0 + 6 * k * H);
  check('daily weather: local mean noon at 90° E is 06:00 UTC', (localNoon(Date.UTC(2026, 9, 22), 90) - Date.UTC(2026, 9, 22, 6)) / H, 1e-12, 'h');
  const day = accumulated(times, () => 1e-3, Date.UTC(2026, 9, 21, 12), Date.UTC(2026, 9, 22, 12));
  check('daily weather: 1 g m-2 s-1 for 24 h is 86.4 mm', Math.abs(day.amount - 86.4) + Math.abs(day.covered - 86400), 1e-9, 'mm');
  const part = accumulated(times, (k) => (k === 2 ? 2e-3 : 0), Date.UTC(2026, 9, 21, 3), Date.UTC(2026, 9, 21, 9));
  check('daily weather: a window over half of a 6 h interval takes half of its rain', Math.abs(part.amount - 2e-3 * 3 * 3600), 1e-9, 'mm');
  const days = forecastNoonDays(times).map((d) => new Date(d).toISOString().slice(5, 10));
  assert('daily weather: a forecast from 20 Oct 18 UTC to 26 Oct 12 UTC has whole fire-weather days 22-25 Oct', days.join(' ') === '10-22 10-23 10-24 10-25', days.join(' '));
}

// ---- heat ----------------------------------------------------------------------------------
{
  check('heat: Stull (2011) wet-bulb temperature at 20 °C and 50 %, 13.7 °C', wetBulbStull(20, 50) - 13.7, 0.05, '°C');
  // Against the isobaric wet-bulb temperature at 1013.25 hPa, found here from the balance
  // c_p (T - Tw) = L (r_s(Tw) - r) with Bolton's (1980) saturation pressure: Stull states
  // errors of -1 to +0.65 °C against his own reference; against this one they reach 1.2 °C.
  const es = (t) => 611.2 * Math.exp((17.67 * t) / (t + 243.5));
  const isobaric = (t, rh) => {
    const p = 101325, e = (rh / 100) * es(t), r = (0.622 * e) / (p - e);
    let lo = -60, hi = t;
    for (let i = 0; i < 80; i++) {
      const w = (lo + hi) / 2, rs = (0.622 * es(w)) / (p - es(w));
      if ((1005.7 + 1875 * r) * (t - w) - (2.501e6 - 2370 * w) * (rs - r) > 0) lo = w; else hi = w;
    }
    return (lo + hi) / 2;
  };
  let worst = 0, sum = 0, n = 0, above = 0;
  for (let t = -20; t <= 50; t += 2) {
    for (let rh = 5; rh <= 99; rh += 2) {
      if (!withinStullRange(t, rh) || (t < 0 && rh < 40)) continue;              // not both cold and dry
      const d = wetBulbStull(t, rh) - isobaric(t, rh);
      worst = Math.max(worst, Math.abs(d)); sum += Math.abs(d); n++;
      if (wetBulbStull(t, rh) > t + 0.05) above++;
    }
  }
  check('heat: Stull against an isobaric psychrometric solution over his range (worst)', worst, 1.25, '°C');
  check('heat: and on average', sum / n, 0.4, '°C');
  assert('heat: the wet-bulb temperature does not exceed the temperature', above === 0, `${above} cases`);
  // The NWS heat index chart (weather.gov), whole °F, at cells away from the adjustments for
  // very dry and very humid air: [T °F, RH %, chart].
  const chart = [[80, 40, 80], [90, 40, 91], [100, 40, 109], [110, 40, 136], [90, 50, 95], [100, 50, 118], [96, 60, 116], [100, 65, 136], [88, 80, 106], [90, 90, 122]];
  check('heat: the NWS heat index against its chart at ten cells, to the chart’s rounding', Math.max(...chart.map(([T, RH, v]) => Math.abs(heatIndexF(T, RH) - v))), 0.5, '°F');
  check('heat: below 80 °F the simple formula: 70 °F at 50 % is 69.7 °F', heatIndexF(70, 50) - 0.5 * (70 + 61 + 1.2 * 2 + 4.7), 1e-12, '°F');
  // Daily maxima over the local solar day.
  const H = 3600000, t0 = Date.UTC(2026, 9, 20, 18), times = [...Array(24).keys()].map((k) => t0 + 6 * k * H);
  const [a, b] = localDay(Date.UTC(2026, 9, 23), 135);
  check('heat: the local day at 135° E runs from 15:00 UTC the day before', (Date.UTC(2026, 9, 23) - a) / H - 9 + (b - a) / H - 24, 1e-12, 'h');
  assert('heat: the greatest sample within a day', sampledMaximum(times, (k) => Math.sin(k), a, b) === Math.max(...times.map((t, k) => (t >= a && t < b ? Math.sin(k) : -9))));
  const days = forecastLocalDays(times).map((d) => new Date(d).toISOString().slice(5, 10));
  assert('heat: the same forecast has whole local days 22-25 Oct at every longitude', days.join(' ') === '10-22 10-23 10-24 10-25', days.join(' '));
}

// ---- settling particles ---------------------------------------------------------------------
{
  // ICAO (1993): 1.225 kg m-3 and 1.7894e-5 Pa s at sea level; 22 632 Pa and 216.65 K at 11 km.
  const s0 = standardAtmosphere(0), s11 = standardAtmosphere(11000);
  check('settling: the standard atmosphere at sea level (density, viscosity)', Math.abs(s0.rho - 1.225) / 1.225 + Math.abs(s0.mu - 1.7894e-5) / 1.7894e-5, 1e-4);
  check('settling: and at 11 km (pressure)', Math.abs(s11.p - 22632) / 22632 + Math.abs(s11.T - 216.65), 1e-4);
  check('settling: pressure altitudes of 850, 500 and 250 hPa (1457, 5574, 10 363 m)',
    Math.max(Math.abs(pressureAltitude(85000) - 1457.3), Math.abs(pressureAltitude(50000) - 5574.4), Math.abs(pressureAltitude(25000) - 10362.9)), 0.2, 'm');
  // By hand: (2500 - 1.225) kg m-3 × 9.80665 m s-2 × (30 µm)² / (18 × 1.7894e-5 Pa s) = 0.06847 m/s.
  const w0 = stokesSettling(30e-6, 2500, 0);
  check('settling: Stokes velocity of a 30 µm grain of 2500 kg m-3 at sea level, 0.06847 m/s', w0.w - 0.06847, 1e-5, 'm/s');
  assert('settling: its Reynolds number, 0.14, is at the edge of Stokes’ law', Math.abs(w0.reynolds - 0.1406) < 1e-3, w0.reynolds.toFixed(4));
  // A parcel in still air falls at w and stays put; in a uniform wind along the equator it
  // moves U t, and lands where its height meets the ground.
  const settle = (z) => stokesSettling(30e-6, 2500, z).w;
  const still = settlingPath({ lon: 10, lat: 45, z: 5000, t: 0 }, () => [0, 0], settle, { dt: 600, steps: 2000 });
  let expected = 0;
  for (let z = 5000; z > 0; z -= 1) expected += 1 / settle(z - 0.5);                // the time to fall 5 km, s
  check('settling: in still air a parcel lands where it started', Math.hypot(still.landed.lon - 10, still.landed.lat - 45), 1e-12, '°');
  check('settling: after the time ∫ dz / w(z) (relative)', (still.landed.t / 1000 - expected) / expected, 1e-4);
  const U = 10, windy = settlingPath({ lon: 0, lat: 0, z: 9000, t: 0 }, () => [U, 0], () => 0, { dt: 1800, steps: 48 });
  const n = windy.times.length - 1;
  check('settling: in a uniform wind along the equator, the displacement is U t', (windy.coords[2 * n] * Math.PI * 6371008.8) / 180 - U * 86400, 1e-6, 'm');
  check('settling: winds between levels are linear in height', Math.abs(windAtHeight(3000, [{ z: 10, u: 2, v: 0 }, { z: 5010, u: 12, v: -5 }])[0] - 7.98) + Math.abs(windAtHeight(9e4, [{ z: 10, u: 2, v: 0 }, { z: 5010, u: 12, v: -5 }])[1] + 5), 1e-12);
}

// ---- the JavaScript standard atmosphere against the C++'s references ------------------------------
{
  // fluids' U.S. Standard Atmosphere 1976 (Z T_M P rho C_s mu k g), against which the C++ is held to
  // 1e-12, from sea level to 20 km. settling.js takes geopotential heights and ICAO's R = 287.052 87
  // J kg-1 K-1, where the 1976 has R*/M0 = 287.053 07: its pressure and density differ by up to 3e-6.
  const R0 = 6356766;
  const rows = readFileSync(resolve(HERE, 'fixtures/atmosphere-reference.txt'), 'utf8').split('\n')
    .filter((l) => l.startsWith('fluids ')).map((l) => l.split(/\s+/).slice(1).map(Number)).filter((r) => r[0] >= 0 && r[0] <= 20000);
  let T = 0, p = 0, rho = 0, mu = 0;
  for (const [z, TM, P, density, , viscosity] of rows) {
    const s = standardAtmosphere((R0 * z) / (R0 + z));
    T = Math.max(T, Math.abs(s.T - TM));
    p = Math.max(p, Math.abs(s.p - P) / P);
    rho = Math.max(rho, Math.abs(s.rho - density) / density);
    mu = Math.max(mu, Math.abs(s.mu - viscosity) / viscosity);
  }
  assert(`standard atmosphere: ${rows.length} reference heights from 0 to 20 km`, rows.length > 40);
  check('standard atmosphere: settling.js temperature against the 1976 (fluids)', T, 1e-9, 'K');
  check('standard atmosphere: settling.js pressure, relative (ICAO\'s R)', p, 5e-6);
  check('standard atmosphere: settling.js density, relative (ICAO\'s R)', rho, 5e-6);
  check('standard atmosphere: settling.js viscosity, relative', mu, 1e-12);
}

// The C++.
cppChecks({ assert, skip, lines }, { project: resolve(HERE, '..'), name: 'meteorology' });

done();
