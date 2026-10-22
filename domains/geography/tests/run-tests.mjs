// Checks of the geography domain: map projections against Šavrič et al. (2018),
// polygons on rasters, joins of data to places, and water in shaded relief.
// Run from the repository root:  node domains/geography/tests/run-tests.mjs
// The checks on real data need the products of data/products/; they are
// skipped, and say so, when a product is missing.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';

import { suite, seeded, PRODUCTS } from '../../../development/tools/checks.mjs';
import {
  equalEarth, equirectangular, forward, inverse, extent, areaScale, fillRings, indexRaster, maskRaster,
  areaFraction, featureAreaFraction, waterMask, joinCountries, valueForYear, latestWellCoveredYear,
} from '../js/index.js';
import { DEG, wrapLongitude } from '../../geodesy/js/index.js';

const { check, assert, skip, done } = suite();
const random = seeded(20261022);

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

// ---- joins on the real products ------------------------------------------------
{
  const geo = resolve(PRODUCTS, 'geography-ne50m.json');
  const ind = resolve(PRODUCTS, 'indicator-population.json');
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
    skip('join: products missing (geography-ne50m, indicator-population)');
  }
}

// ---- the land fraction from the geography product -------------------------------
{
  const geo = resolve(PRODUCTS, 'geography-ne50m.json');
  if (existsSync(geo)) {
    const g = JSON.parse(readFileSync(geo, 'utf8'));
    const bin = readFileSync(resolve(PRODUCTS, g.coordinates.file));
    const coords = new Float32Array(bin.buffer, bin.byteOffset, bin.byteLength / 4);
    const W = 2048, H = 1024;
    const land = maskRaster(coords, g.layers.land, W, H);
    const rastered = areaFraction(land, W, H, (v) => v === 1);
    const exact = g.layers.land.reduce((s, f) => s + featureAreaFraction(coords, f.rings), 0);
    check('raster: land area from the raster equals that of the polygons', rastered - exact, 5e-4);
    // The land of the Earth is 29.2 % of its surface (148.9 of 510.1 million km^2, with
    // Antarctica's ice shelves, which Natural Earth's land omits: about 0.3 %).
    check('raster: Natural Earth land is 28.9 +/- 0.3 % of the sphere', exact - 0.289, 0.003);
  } else {
    skip('land fraction: geography-ne50m missing');
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

done();
