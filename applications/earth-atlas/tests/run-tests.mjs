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

done();
