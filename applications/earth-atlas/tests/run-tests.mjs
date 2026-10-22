// Checks of what is the Earth Atlas's own: its reading of the data library's
// products, the checks that join two domains through those products (the plates'
// boundaries; our climate classes against the published maps), and its climate
// diagram. The science it shows is checked with its domains and the shared
// visualization layer.
// Run from the repository root:  node applications/earth-atlas/tests/run-tests.mjs
// The checks need the products of data/products/ (see the README); they are
// skipped, and say so, when a product is missing.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';

import { suite, PRODUCTS } from '../../../development/tools/checks.mjs';
import { FieldProduct, ClassRaster } from '../src/data/fields.js';
import { precipitationScale, climateDiagram } from '../src/ui/climate.js';
import { areaMean } from '../../../foundations/physis/js/index.js';
import { fillRings, featureAreaFraction } from '../../../domains/geography/js/index.js';
import { omega, poleOf, velocityAt, subtract, acrossStep, azimuth } from '../../../domains/geophysics/js/index.js';
import { koppenGrid, KOPPEN_CODES, meanTemperature } from '../../../domains/climate/js/index.js';
import { rmsDifference } from '../../../domains/meteorology/js/index.js';

const { check, assert, skip, done } = suite();

// ---- a field product -------------------------------------------------------------
{
  const meta = resolve(PRODUCTS, 'fields-gfs-forecast.json');
  if (existsSync(meta)) {
    const m = JSON.parse(readFileSync(meta, 'utf8'));
    const bin = readFileSync(resolve(PRODUCTS, m.file));
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
    skip('field: fields-gfs-forecast missing');
  }
}

// ---- the plates' boundaries (geophysics on geography's rasters) -----------------------
{
  const meta = resolve(PRODUCTS, 'plates-pb2002.json');
  if (existsSync(meta)) {
    const m = JSON.parse(readFileSync(meta, 'utf8'));
    const b = readFileSync(resolve(PRODUCTS, m.file));
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
    skip('plates: plates-pb2002 missing');
  }
}

// ---- the climate diagram ------------------------------------------------------------
{
  const season = (a, b) => [...Array(12).keys()].map((m) => a + b * Math.cos((2 * Math.PI * (m - 6.5)) / 12));
  const flat = (v) => Array(12).fill(v);
  const halves = (summer, winter) => [...Array(12).keys()].map((m) => (m >= 3 && m <= 8 ? summer[m - 3] : winter[(m + 3) % 12]));
  // The Walter-Lieth precipitation scale.
  check('climate: 100 mm stands at 50 degrees, 300 mm at 60', Math.abs(precipitationScale(100) - 50) + Math.abs(precipitationScale(300) - 60), 1e-12);
  const svg = climateDiagram({ T: season(-12, 30), P: halves(flat(60), flat(2)), Tmin: season(-18, 30), lat: 62 });
  const frost = (svg.match(/class="frost"/g) || []).length;
  assert('climate: the diagram marks the frost months', frost === 8, `${frost} months with mean daily minimum below 0 C`);
}

// ---- climate: the published maps and our classification -----------------------------
{
  const km = resolve(PRODUCTS, 'classes-koppen-geiger-1991-2020.json');
  const tm = resolve(PRODUCTS, 'fields-terraclimate-normals-medium.json');
  if (existsSync(km) && existsSync(tm)) {
    const meta = JSON.parse(readFileSync(km, 'utf8'));
    const beck = new ClassRaster(meta, readFileSync(resolve(PRODUCTS, meta.file)));
    assert('classes: the legend of the maps lists the classes in our order',
      meta.classes.map((c) => c.code).join(' ') === KOPPEN_CODES.join(' '));
    // Places of known climate on Beck et al.'s 1991-2020 map.
    const places = [['Paris', 2.35, 48.85, 'Cfb'], ['Cairo', 31.24, 30.04, 'BWh'], ['Singapore', 103.8, 1.35, 'Af'],
      ['Moscow', 37.6, 55.75, 'Dfb'], ['Rome', 12.5, 41.9, 'Csa'], ['Mumbai', 72.88, 19.08, 'Am'], ['Tokyo', 139.7, 35.68, 'Cfa']];
    const off = places.filter(([, lon, lat, code]) => beck.classAt(lon, lat)?.code !== code);
    assert('classes: seven cities in their classes on the 1991-2020 map', !off.length, off.map((p) => `${p[0]} ${beck.classAt(p[1], p[2])?.code}`).join(', '));
    const m = JSON.parse(readFileSync(tm, 'utf8'));
    const bin = readFileSync(resolve(PRODUCTS, m.file));
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
    skip('classes: classes-koppen-geiger-1991-2020 or fields-terraclimate-normals-medium missing');
  }
}

// ---- plug-ins (ADR 0027): the examples, run from synthetic fields -----------------------
{
  const api = await import('../src/plugins/index.js');
  const context = (time, pinned = null) => ({ time, pinned, progress() {}, fields: () => [], field: () => null, geography: () => null });
  const grid = { lon0: -180, dlon: 2.5, nlon: 144, lat0: 90, dlat: 2.5, nlat: 73, periodic: true };
  const HOUR = 3600000, t0 = Date.UTC(2026, 9, 20), times = [t0, t0 + 48 * HOUR];
  const U = 20, R = 6371008.8, DEG = Math.PI / 180;
  const fieldOf = (fn) => {
    const a = new Float32Array(grid.nlon * grid.nlat);
    for (let r = 0; r < grid.nlat; r++) for (let c = 0; c < grid.nlon; c++) a[r * grid.nlon + c] = fn(grid.lon0 + c * grid.dlon, grid.lat0 - r * grid.dlat);
    return a;
  };
  const input = (product, names, level = null) => ({ product, names, variables: names.map((n) => product.meta.variables.find((v) => v.name === n)), meta: product.meta, model: product.meta.title, level });
  // A solid-body rotation, u = U cos φ, v = 0: its vorticity is 2 U sin φ / a, its divergence 0,
  // and it is steady under the barotropic vorticity equation.
  const zonal = api.makeProduct({ name: 'test-zonal', title: 'Solid-body rotation', grid, times, variables: [
    { name: 'u', long_name: 'Eastward wind at 500 hPa', standard_name: 'eastward_wind', units: 'm s-1', level: '500 hPa', slices: times.map(() => fieldOf((lon, lat) => U * Math.cos(lat * DEG))) },
    { name: 'v', long_name: 'Northward wind at 500 hPa', standard_name: 'northward_wind', units: 'm s-1', level: '500 hPa', slices: times.map(() => fieldOf(() => 0)) },
  ] });
  const wind = input(zonal, ['u', 'v'], '500 hPa');

  const vorticity = (await import('../plugins/vorticity.js')).create(api);
  const [vort] = await vorticity.run({ wind }, {}, context(t0));
  const zeta = vort.product.slice('zeta', 0), delta = vort.product.slice('delta', 0), row45 = Math.round((grid.lat0 - 45) / grid.dlat) * grid.nlon;
  check('plug-in vorticity: a solid-body rotation has ζ = 2 U sin φ / a at 45° N', Math.abs(zeta[row45 + 10] - (2 * U * Math.sin(45 * DEG)) / R), 2e-8, '1/s');
  check('plug-in vorticity: and no divergence', Math.abs(delta[row45 + 10]), 1e-9, '1/s');
  assert('plug-in vorticity: the product carries the input level and both variables', vort.product.meta.variables.every((v) => v.level === '500 hPa') && vort.product.has('delta'));

  const parcels = (await import('../plugins/trajectories.js')).create(api);
  const [paths] = await parcels.run({ wind }, { hours: 24, step: 1, count: 1, radius: 0 }, context(t0, [10, 40]));
  const path = paths.polylines[0].coords, n = path.length / 2;
  check('plug-in trajectories: a parcel in a solid-body rotation keeps its latitude', Math.abs(path[2 * n - 1] - 40), 1e-6, '°');
  check('plug-in trajectories: and turns by U t / a in a day', Math.abs(path[2 * n - 2] - 10 - (U * 86400) / R / DEG), 0.01, '°');
  assert('plug-in trajectories: the moving point starts at the seed', paths.moving[0].times[0] === t0 && Math.abs(paths.moving[0].coords[0] - 10) < 1e-9);

  // Optimal interpolation with two observations at grid points, far apart, and σo → 0:
  // the analysis takes each observation's value where it is, and is the background far from both.
  const oi = (await import('../plugins/optimal-interpolation.js')).create(api);
  class OnGrid extends oi.constructor { async observations() { return [{ lon: 0, lat: 40, value: 283 }, { lon: 90, lat: -20, value: 278 }]; } }
  const background = api.makeProduct({ name: 'test-background', title: 'Background', grid, times: [t0], variables: [
    { name: 'tmp2m', long_name: 'air temperature at 2m', standard_name: 'air_temperature', units: 'K', slices: [fieldOf(() => 280)] }] });
  const bg = input(background, ['tmp2m'], '2 m');
  const [analysis] = await new OnGrid().run({ background: bg, truth: bg }, { count: 2, length: 500, sigmaB: 2, sigmaO: 0.01 }, context(t0));
  const an = analysis.product.slice('tmp2m', 0), inc = analysis.product.slice('tmp2m_increment', 0);
  const cell = (lon, lat) => Math.round((grid.lat0 - lat) / grid.dlat) * grid.nlon + Math.round((lon - grid.lon0) / grid.dlon);
  check('plug-in optimal interpolation: an exact observation is reproduced where it is', Math.abs(an[cell(0, 40)] - 283), 1e-3, 'K');
  check('plug-in optimal interpolation: and the second too', Math.abs(an[cell(90, -20)] - 278), 1e-3, 'K');
  check('plug-in optimal interpolation: the background is kept far from the observations', Math.abs(inc[cell(-150, -70)]), 1e-6, 'K');

  const barotropic = (await import('../plugins/barotropic.js')).create(api);
  const height = api.makeProduct({ name: 'test-height', title: 'Height', grid, times, variables: [
    { name: 'z', long_name: 'Geopotential height at 500 hPa', standard_name: 'geopotential_height', units: 'm', level: '500 hPa', slices: times.map(() => fieldOf(() => 5500)) }] });
  const [forecast] = await barotropic.run({ wind, height: input(height, ['z'], '500 hPa') }, { hours: 6, truncation: 21 }, context(t0));
  const z0 = forecast.product.slice('z@500hPa', 0), z6 = forecast.product.slice('z@500hPa', 1);
  check('plug-in barotropic: a solid-body rotation is steady (height change over 6 h)', rmsDifference(z0, z6, barotropic.grid()), 0.5, 'm');
  assert('plug-in barotropic: two times recorded, as a model of height at 500 hPa', forecast.product.times.length === 2 && forecast.product.meta.variables[0].standard_name === 'geopotential_height');

  const smallest = (await import('../plugins/wind-speed.js')).create(api);
  const [speed] = await smallest.run({ wind }, {}, context(t0));
  check('plug-in wind speed (the smallest): U cos φ at 45° N', Math.abs(speed.product.slice('speed', 0)[row45 + 10] - U * Math.cos(45 * DEG)), 1e-4, 'm/s');
}

done();
