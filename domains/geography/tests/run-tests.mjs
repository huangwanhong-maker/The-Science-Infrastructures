// Checks of the geography domain: map projections against Šavrič et al. (2018),
// polygons on rasters, joins of data to places, water in shaded relief, zonal
// statistics on the sphere against closed forms, and drainage from elevation.
// Run from the repository root:  node domains/geography/tests/run-tests.mjs
// The checks on real data need the products of data/products/; they are
// skipped, and say so, when a product is missing.

import { existsSync, readFileSync } from 'node:fs';
import { resolve } from 'node:path';

import { suite, seeded, PRODUCTS } from '../../../development/tools/checks.mjs';
import {
  equalEarth, equirectangular, forward, inverse, extent, areaScale, fillRings, indexRaster, maskRaster,
  areaFraction, featureAreaFraction, waterMask, joinCountries, valueForYear, latestWellCoveredYear,
  cellAreas, overlapTable, zonalStatistics, regionAt, byCode,
  slope, outlets, fillDepressions, flowDirections, flowAccumulation, watershed, downstream, mainStem,
  cellOutline, pathOf, cellOf, neighbour, receiver,
} from '../js/index.js';
import { DEG, wrapLongitude, EARTH_RADIUS_M, distanceM as distance } from '../../geodesy/js/index.js';

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

// ---- zonal statistics on the sphere ------------------------------------------------
{
  const R = EARTH_RADIUS_M;
  // Regions on a 0.5° grid: the cells nearest each of 12 random seeds (a Voronoi
  // partition), with two latitude bands cut out of them and none beyond 80° S;
  // fields on grids not aligned with it, one coarser (0.75°, shifted) and one finer (0.2°).
  const rg = { lon0: -179.75, dlon: 0.5, nlon: 720, lat0: 89.75, dlat: 0.5, nlat: 360, periodic: true };
  const seeds = Array.from({ length: 12 }, () => [random() * 360 - 180, Math.asin(2 * random() - 1) / DEG]);
  const regions = new Uint16Array(rg.nlon * rg.nlat);
  for (let r = 0; r < rg.nlat; r++) {
    for (let c = 0; c < rg.nlon; c++) {
      const lon = rg.lon0 + c * rg.dlon, lat = rg.lat0 - r * rg.dlat;
      let best = 0, near = Infinity;
      seeds.forEach(([sl, sp], k) => { const d = distance(lon, lat, sl, sp); if (d < near) { near = d; best = k; } });
      let id = best + 3;
      if (lat > 10 && lat < 20.5) id = 1;          // the band 10°–20.5° N
      if (lat > -30 && lat < -5) id = 2;           // the band 30°–5° S
      if (lat < -80) id = 0;                       // no region beyond 80° S
      regions[r * rg.nlon + c] = id;
    }
  }
  const coarse = { lon0: -179.6, dlon: 0.75, nlon: 480, lat0: 89.625, dlat: 0.75, nlat: 240, periodic: true };
  const fine = { lon0: -179.9, dlon: 0.2, nlon: 1800, lat0: 89.9, dlat: 0.2, nlat: 900, periodic: true };
  const band = (p1, p2) => 2 * Math.PI * R * R * (Math.sin(p2 * DEG) - Math.sin(p1 * DEG));
  for (const [label, g] of [['a coarser field', coarse], ['a finer field', fine]]) {
    const table = overlapTable(g, rg);
    const constant = new Float32Array(g.nlon * g.nlat).fill(7.25);
    const s = zonalStatistics(constant, g, regions, rg, { count: 14, table });
    let worst = 0;
    for (let k = 1; k <= 14; k++) worst = Math.max(worst, Math.abs(s.mean[k] - 7.25), Math.abs(s.min[k] - 7.25), Math.abs(s.max[k] - 7.25));
    check(`zonal (${label}): a constant field has that mean, least and greatest value in every region`, worst, 1e-9);
    check(`zonal (${label}): the area of the band 10°–20.5° N is 2πR²(sin φ₂ − sin φ₁)`, s.area[1] / band(10, 20.5) - 1, 1e-9);
    check(`zonal (${label}): and that of the band 30°–5° S`, s.area[2] / band(-30, -5) - 1, 1e-9);
    check(`zonal (${label}): the regions and what lies in none cover the sphere, 4πR²`, s.area.reduce((a, x) => a + x, 0) / (4 * Math.PI * R * R) - 1, 1e-9);
    check(`zonal (${label}): what lies in no region is the cap beyond 80° S`, s.area[0] / band(-90, -80) - 1, 1e-9);
    // Counts, random, with cells of no value: every count is kept, however it is shared.
    const counts = new Float32Array(g.nlon * g.nlat);
    for (let i = 0; i < counts.length; i++) counts[i] = random() < 0.1 ? NaN : Math.floor(random() * 1000);
    for (const [share, reach] of [['regions', 0], ['regions', 4], ['area', 0]]) {
      const t = zonalStatistics(counts, g, regions, rg, { count: 14, share, reach, table });
      check(`zonal (${label}): counts shared by ${share}${reach ? `, reaching ${reach} cells,` : ''} sum to the total`, Math.abs(t.sum.reduce((a, x) => a + x, 0) / t.total - 1), 1e-12);
    }
    // The cap beyond 80° S is wider than 4 cells of 0.5°: only its rim reaches a region.
    const plain = zonalStatistics(counts, g, regions, rg, { count: 14, table }), reaching = zonalStatistics(counts, g, regions, rg, { count: 14, reach: 4, table });
    assert(`zonal (${label}): reaching 4 cells, the counts of the cap's rim go to the nearest region`,
      reaching.reached > 0 && reaching.unassigned > 0 && Math.abs(reaching.unassigned + reaching.reached - plain.unassigned) < 1e-6 * plain.unassigned);
  }
  const areas = cellAreas(coarse);
  check('zonal: the cells of a global grid cover 4πR²', (areas.reduce((a, x) => a + x, 0) * coarse.nlon) / (4 * Math.PI * R * R) - 1, 1e-12);
  assert('zonal: the region at a place is that of the cell containing it', regionAt(regions, rg, 0, 15) === 1 && regionAt(regions, rg, 179.99, -10) === 2 && regionAt(regions, rg, 0, -85) === 0);
}

// ---- terrain and water --------------------------------------------------------------------
{
  const R = EARTH_RADIUS_M;
  const make = (g, f) => {
    const z = new Float64Array(g.nlon * g.nlat);
    for (let r = 0; r < g.nlat; r++) for (let c = 0; c < g.nlon; c++) z[r * g.nlon + c] = f(r, c, g.lon0 + c * g.dlon, g.lat0 - r * g.dlat);
    return z;
  };
  const drain = (z, g, options = {}) => {
    const { filled, outlet, raised } = fillDepressions(z, g, options);
    const dirs = flowDirections(filled, g, outlet);
    return { filled, outlet, raised, dirs, acc: flowAccumulation(dirs, g) };
  };
  // A tilted plane drains in one direction: rising eastward it drains west; northward, south.
  const g = { lon0: 10, dlon: 0.1, nlon: 40, lat0: 45, dlat: 0.1, nlat: 30, periodic: false };
  for (const [label, f, dir] of [['east', (r, c) => 100 + 50 * c, 4], ['north', (r) => 100 + 50 * (g.nlat - r), 6]]) {
    const { dirs, outlet } = drain(make(g, f), g);
    let wrong = 0;
    dirs.forEach((d, i) => { if (!outlet[i] && d !== dir) wrong++; });
    assert(`hydrology: a plane rising to the ${label} drains the other way from every cell`, wrong === 0, `${wrong} cells otherwise`);
  }
  {
    const { dirs } = drain(make(g, (r, c) => 100 + 50 * c), g);
    const start = 15 * g.nlon + 25, path = downstream(dirs, g, start), { length } = pathOf(g, path);
    const lat = g.lat0 - 15 * g.dlat;
    assert('hydrology: downstream on that plane, straight west to the edge', path.length === 26 && path[25] === 15 * g.nlon);
    check('hydrology: and the path is 25 steps along the parallel on the sphere', length / (25 * distance(0, lat, g.dlon, lat)) - 1, 1e-12);
  }
  // Slopes of planes on the sphere: z = a R φ has slope a; z = b λ has slope b / (R cos φ).
  {
    const a = 0.02, b = 5e4;
    const s1 = slope(make(g, (r, c, lon, lat) => a * R * lat * DEG), g), s2 = slope(make(g, (r, c, lon) => b * lon * DEG), g);
    const i = 12 * g.nlon + 20, lat = g.lat0 - 12 * g.dlat;
    check('hydrology: the slope of z = a R φ is a', s1[i] / a - 1, 1e-6);
    check('hydrology: the slope of z = b λ is b / (R cos φ)', s2[i] / (b / (R * Math.cos(lat * DEG))) - 1, 1e-6);
  }
  // A pit is filled to its spill level: a bowl inside a rim of 200 m with one notch at 150 m.
  {
    const h = { lon0: 0, dlon: 0.1, nlon: 21, lat0: 1, dlat: 0.1, nlat: 21, periodic: false };
    const z = make(h, (r, c) => {
      const d = Math.max(Math.abs(r - 10), Math.abs(c - 10));
      if (d === 6) return r === 4 && c === 10 ? 150 : 200;
      return d < 6 ? 70 + 5 * d : 120 - 5 * (d - 6);
    });
    const { filled, raised } = drain(z, h);
    let worst = 0, inside = 0;
    for (let r = 5; r <= 15; r++) for (let c = 5; c <= 15; c++) { worst = Math.max(worst, Math.abs(filled[r * h.nlon + c] - 150)); inside++; }
    check('hydrology: a pit is filled to the level of its spill point (150 m)', worst, 1e-9, 'm');
    assert('hydrology: every cell of the pit was raised, and none outside', raised === inside, `${raised} of ${inside}`);
  }
  // A closed basin (its edges no outlet) draining through one ocean cell, with depressions
  // inside: the area draining through the outlet is the whole basin's, R² Δλ (sin φ₂ − sin φ₁).
  {
    const k = { lon0: 20.05, dlon: 0.1, nlon: 60, lat0: 59.95, dlat: 0.1, nlat: 50, periodic: false };
    const out = 25 * k.nlon + 30;
    const z = make(k, (r, c) => (r * k.nlon + c === out ? -5 : 10 + Math.hypot(r - 25, c - 30) + 30 * Math.sin(r / 3) ** 2 * Math.cos(c / 4) ** 2));
    const { acc, outlet, dirs } = drain(z, k, { edges: false, minSeaArea: 0 });
    const basin = (R * R * k.nlon * k.dlon * DEG * (Math.sin(60 * DEG) - Math.sin(55 * DEG))) / 1e6;
    assert('hydrology: the closed basin has one outlet', outlet.reduce((a, x) => a + x, 0) === 1);
    check('hydrology: the area draining through its outlet is the basin’s, R² Δλ (sin φ₂ − sin φ₁)', acc[out] / basin - 1, 1e-12);
    assert('hydrology: and its watershed is every cell', watershed(dirs, k, out).length === k.nlon * k.nlat);
    const below = receiver(k, dirs, 3 * k.nlon + 4), shed = watershed(dirs, k, below), rows = cellAreas(k);
    const area = Array.from(shed).reduce((a, i) => a + rows[Math.floor(i / k.nlon)] / 1e6, 0);
    check('hydrology: a watershed’s area is the upstream area at its outlet', area / acc[below] - 1, 1e-12);
    const stem = mainStem(dirs, acc, k, out), top = stem[stem.length - 1];
    assert('hydrology: the main stem climbs from the outlet to a cell into which none drains', stem[0] === out && stem.length > 10
      && [0, 1, 2, 3, 4, 5, 6, 7].every((d) => { const j = neighbour(k, top, d); return j < 0 || receiver(k, dirs, j) !== top; }));
  }
  // The watershed of a cone's apex is the apex alone, outlined by the four corners of its cell.
  {
    const c = { lon0: 0, dlon: 0.1, nlon: 31, lat0: 30, dlat: 0.1, nlat: 31, periodic: false };
    const { dirs } = drain(make(c, (r, q) => 1000 - 10 * Math.hypot(r - 15, q - 15)), c), apex = cellOf(c, 1.5, 28.5);
    const shed = watershed(dirs, c, apex), ring = cellOutline(c, shed)[0];
    assert('hydrology: the watershed of a cone’s apex is the apex alone', shed.length === 1 && shed[0] === 15 * c.nlon + 15);
    assert('hydrology: outlined by the four corners of its cell', ring.length === 10 && ring[0] === ring[8] && ring[1] === ring[9]);
    // An L of three cells: six corners, and the outline encloses the cells' area.
    const rings = cellOutline(c, [10 * c.nlon + 10, 11 * c.nlon + 10, 11 * c.nlon + 11]), cells = cellAreas(c);
    const area = featureAreaFraction(Float64Array.from(rings[0]), [[0, rings[0].length / 2]]) * 4 * Math.PI * R * R;
    assert('hydrology: the outline of an L of three cells has six corners', rings.length === 1 && rings[0].length === 14);
    check('hydrology: and encloses their area', area / (cells[10] + 2 * cells[11]) - 1, 1e-6);
  }
  // On a grid round the Earth flow crosses the antimeridian: a sea at one column, land rising away from it.
  {
    const w = { lon0: -179, dlon: 2, nlon: 180, lat0: 0, dlat: 2, nlat: 1, periodic: true };
    const { dirs } = drain(make(w, (r, c) => (c === 0 ? -1 : 100 + 10 * Math.min(c, 180 - c))), w, { minSeaArea: 0 });
    assert('hydrology: on a grid round the Earth, flow crosses the antimeridian', dirs[179] === 0 && dirs[1] === 4 && dirs[90] >= 0);
  }
}

// ---- real data: people per country, the regions' rasters, and the rivers of ETOPO --------------
{
  const read = (name) => {
    const m = JSON.parse(readFileSync(resolve(PRODUCTS, name), 'utf8'));
    return { m, b: readFileSync(resolve(PRODUCTS, m.file)) };
  };
  const decode = ({ m, b }, name) => {
    const v = m.variables.find((x) => x.name === name), q = new Int16Array(b.buffer, b.byteOffset, b.byteLength / 2);
    const n = m.grid.nlon * m.grid.nlat, out = new Float32Array(n);
    for (let i = 0; i < n; i++) { const x = q[v.offsets[0] / 2 + i]; out[i] = x === -32768 ? NaN : x * v.scale + v.offset; }
    return out;
  };
  const ids = ({ b }) => new Uint16Array(b.buffer, b.byteOffset, b.byteLength / 2);
  const have = (...names) => names.every((n) => existsSync(resolve(PRODUCTS, n)));
  if (have('fields-ghsl-population.json', 'regions-countries.json', 'regions-admin1.json')) {
    const p = read('fields-ghsl-population.json'), people = decode(p, 'population');
    for (const name of ['countries', 'admin1']) {
      const r = read(`regions-${name}.json`), t0 = Date.now();
      const s = zonalStatistics(people, p.m.grid, ids(r), r.m.grid, { count: r.m.regions.length, reach: 3 });
      const ms = Date.now() - t0;
      check(`zonal on GHS-POP by ${name} (${ms} ms): every person is counted once`, Math.abs(s.sum.reduce((a, x) => a + x, 0) / s.total - 1), 1e-9);
      check(`zonal on GHS-POP by ${name}: under 0.05 % of the people beyond 0.3° of every region`, s.unassigned / s.total, 5e-4);
      const by = byCode(s, r.m.regions);
      // Japan: 126 146 099 at the census of 1 October 2020 (Statistics Bureau of Japan,
      // 2020 Population Census, final count). GHS-POP's epoch is 2025, by when the UN's
      // estimates fall by about 2 %; and coasts and islands finer than the cells lose a little.
      const japan = name === 'countries' ? by.JPN : r.m.regions.filter((x) => x.country === 'JPN').reduce((a, x) => a + (by[x.code] || 0), 0);
      check(`zonal on GHS-POP by ${name}: Japan, ${(japan / 1e6).toFixed(2)} million, within 5 % of its census of 2020`, japan / 126146099 - 1, 0.05);
    }
  } else {
    skip('zonal on real data: fields-ghsl-population, regions-countries or regions-admin1 missing');
  }
  if (have('regions-countries.json', 'geography-ne50m.json')) {
    // The regions' raster (data/products/regions.py) and the domain's own fillRings, from the same polygons.
    const r = read('regions-countries.json'), stored = ids(r);
    const g = JSON.parse(readFileSync(resolve(PRODUCTS, 'geography-ne50m.json'), 'utf8'));
    const bin = readFileSync(resolve(PRODUCTS, g.coordinates.file)), coords = new Float32Array(bin.buffer, bin.byteOffset, bin.byteLength / 4);
    const ours = indexRaster(coords, g.layers.countries, r.m.grid.nlon, r.m.grid.nlat);
    let differ = 0;
    for (let i = 0; i < ours.length; i++) if (ours[i] !== stored[i]) differ++;
    check('regions: the countries’ raster is fillRings’s, but where float32 moves a vertex across a centre', differ / ours.length, 1e-5);
    assert('regions: the countries are those of the geography product, in its order', r.m.regions.every((x, k) => x.adm0_a3 === g.layers.countries[k].id));
  } else {
    skip('regions against fillRings: regions-countries or geography-ne50m missing');
  }
  if (have('fields-etopo-2022-medium.json')) {
    const e = read('fields-etopo-2022-medium.json'), z = decode(e, 'z'), grid = e.m.grid, t0 = Date.now();
    const { filled, outlet } = fillDepressions(z, grid);
    const t1 = Date.now(), dirs = flowDirections(filled, grid, outlet), acc = flowAccumulation(dirs, grid), t2 = Date.now();
    let stuck = 0;
    for (let i = 0; i < dirs.length; i++) if (!outlet[i] && dirs[i] < 0) stuck++;
    assert(`hydrology on ETOPO 10′ (${grid.nlon * grid.nlat} cells; filled in ${t1 - t0} ms, flow in ${t2 - t1} ms): every land cell drains`, stuck === 0, `${stuck} cells without a way down`);
    // Five great basins, at the cell of sea near each mouth that receives the most. Their
    // areas: Amazon 6.20 (without the Tocantins, which meets the sea apart), Murray–Darling
    // 1.05 million km² (World Resources Institute 2003, Watersheds of the World); Mackenzie
    // 1.81 and St Lawrence 1.34, with the Great Lakes, whose floors lie below sea level
    // and must not be outlets (Natural Resources Canada, The Atlas of Canada); Río de la
    // Plata 3.17 (Guerrero et al. 1997, Cont. Shelf Res. 17, 727). On a 10′ grid by D8, within 10 %.
    const mouth = (lon, lat, reach) => {
      const c0 = cellOf(grid, lon, lat), r0 = Math.floor(c0 / grid.nlon), q0 = c0 % grid.nlon;
      let best = -1;
      for (let dr = -reach; dr <= reach; dr++) for (let dc = -reach; dc <= reach; dc++) {
        const i = (r0 + dr) * grid.nlon + ((q0 + dc + grid.nlon) % grid.nlon);
        if (outlet[i] && (best < 0 || acc[i] > acc[best])) best = i;
      }
      return acc[best] / 1e6;
    };
    for (const [river, lon, lat, area] of [['Amazon', -50.5, 0, 6.2], ['Río de la Plata', -58.3, -34.5, 3.17], ['Mackenzie', -135.2, 69.2, 1.805],
      ['St Lawrence', -70.9, 46.9, 1.344], ['Murray–Darling', 138.9, -35.55, 1.05]]) {
      const a = mouth(lon, lat, 4);
      check(`hydrology on ETOPO 10′: the ${river}’s basin, ${a.toFixed(2)} million km², against ${area}`, a / area - 1, 0.1);
    }
  } else {
    skip('hydrology on real data: fields-etopo-2022-medium missing');
  }
}

done();
