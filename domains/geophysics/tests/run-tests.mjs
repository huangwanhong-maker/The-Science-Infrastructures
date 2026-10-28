// Checks of the geophysics domain: the motion of the plates against the
// EarthScope plate calculator; GPS velocities relative to their plates;
// earthquake mechanisms against double couples and the Global CMT catalogue;
// the crust (CRUST1.0); sections along great circles; and aftershock rates
// and probabilities (Reasenberg and Jones 1989).
// And its C++ (the geomagnetic field: the World Magnetic Model 2025), built with CMake and run on
// calx's verification, against NOAA's test values and the WMM2025 report's numerical example; its
// coefficients come through the data library (python data/acquisition/fetch.py noaa-wmm).
// Run from the repository root:  node domains/geophysics/tests/run-tests.mjs
// The checks on real data need the products of data/products/; they are
// skipped, and say so, when a product is missing. The checks of the plates'
// boundaries, which need geography's rasters as well, are the Earth Atlas's.

import { existsSync, readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

import { suite, seeded, cppChecks, PRODUCTS } from '../../../development/tools/checks.mjs';
import {
  NNR_MORVEL56, poleOf, omega, poleOfOmega, velocityAt, subtract, acrossStep, azimuth, morvelCode, itrfOmega,
  itrfVelocityAt, mechanism, tensorNED, eigenSymmetric, radiation, fromEqualArea, toEqualArea, planeDifference,
  cellIndex, column as crustColumn, travel, relativeTo, sectionLine, R_KM,
} from '../js/index.js';
import { GENERIC_CALIFORNIA, aftershockRate, expectedAftershocks, aftershockProbability, subsurfaceRuptureLengthKm } from '../js/index.js';

const { check, assert, skip, lines, done } = suite();
const random = seeded(20261022);

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
}

// ---- GPS stations ---------------------------------------------------------------------
{
  // ITRF2020's plate motion model against the EarthScope calculator (model itrf2020, which
  // leaves out the origin rate bias): the Pacific plate at 19.5 N, 155.5 W, retrieved 2026-10-22.
  const v = velocityAt(itrfOmega('PA'), -155.5, 19.5);
  check('gnss: the ITRF2020 plate motion model against the EarthScope calculator', Math.max(Math.abs(v.east + 62.24), Math.abs(v.north - 33.90)), 0.006, 'mm/yr');
  const sm = resolve(PRODUCTS, 'stations-ngl-midas.json');
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
    const sp = resolve(PRODUCTS, 'series/series-ngl-p123.json');
    if (existsSync(sp)) {
      const sr = JSON.parse(readFileSync(sp, 'utf8'));
      const years = (Date.parse(sr.days[sr.days.length - 1]) - Date.parse(sr.days[0])) / (365.25 * 86400000);
      const st = d.stations.find((r) => r[f.id] === 'P123');
      const de = Math.abs(sr.components[0].values.at(-1) / years - st[f.ve]), dn = Math.abs(sr.components[1].values.at(-1) / years - st[f.vn]);
      check('gnss: P123’s 20-year record moves at its MIDAS velocity', Math.max(de, dn), 0.5, 'mm/yr');
    } else {
      skip('gnss: the series of P123 is missing (requested on demand)');
    }
  } else {
    skip('gnss: stations-ngl-midas missing');
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
  const gm = resolve(PRODUCTS, 'events-gcmt.json');
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
    skip('mechanisms: events-gcmt missing');
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
  const cm = resolve(PRODUCTS, 'crust-crust1.json');
  if (existsSync(cm)) {
    const meta = JSON.parse(readFileSync(cm, 'utf8'));
    const buf = readFileSync(resolve(PRODUCTS, meta.file));
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
    skip('crust: crust-crust1 missing');
  }
}

// ---- aftershocks ---------------------------------------------------------------------------
{
  // The closed-form count against the rate integrated numerically (Simpson's rule on a
  // logarithmic variable, t + c = e^s, where the integrand is smooth).
  const simpson = (t1, t2, M, Mm, params) => {
    const { c } = params, a = Math.log(t1 + c), b = Math.log(t2 + c), n = 2000, h = (b - a) / n;
    let sum = 0;
    for (let k = 0; k <= n; k++) {
      const s = a + k * h, w = k === 0 || k === n ? 1 : k % 2 ? 4 : 2;
      sum += w * aftershockRate(Math.exp(s) - c, M, Mm, params) * Math.exp(s);
    }
    return (sum * h) / 3;
  };
  let worst = 0;
  for (const [t1, t2, M, Mm] of [[0, 1, 5, 6.5], [0.5, 7.5, 6, 7.7], [10, 40, 5, 7], [0, 365, 4, 6]]) {
    for (const params of [GENERIC_CALIFORNIA, { ...GENERIC_CALIFORNIA, p: 1 }, { a: -2.2, b: 1.0, p: 0.9, c: 0.02 }]) {
      worst = Math.max(worst, Math.abs(expectedAftershocks(t1, t2, M, Mm, params) / simpson(t1, t2, M, Mm, params) - 1));
    }
  }
  check('aftershocks: the closed-form count against numerical integration (p < 1, = 1, > 1)', worst, 1e-9);
  // Gutenberg-Richter scaling: a unit of magnitude more divides the count by 10^b.
  check('aftershocks: one unit of magnitude more, 10^-0.91 as many', expectedAftershocks(0, 7, 6, 7) / expectedAftershocks(0, 7, 5, 7) - 10 ** -0.91, 1e-12);
  // The probability of an M ≥ 5 aftershock of an M 7 over the next day falls with time since the mainshock.
  const later = [0, 1, 3, 7, 30, 100].map((t) => aftershockProbability(t, t + 1, 5, 7));
  assert('aftershocks: probabilities over the next day decrease with time since the mainshock', later.every((x, k) => !k || x < later[k - 1]),
    later.map((x) => x.toFixed(3)).join(', '));
  assert('aftershocks: and a week holds more than a day, a month more than a week', aftershockProbability(1, 2, 6, 7) < aftershockProbability(1, 8, 6, 7) && aftershockProbability(1, 8, 6, 7) < aftershockProbability(1, 31, 6, 7));
  // By hand: 10^(-1.67 + 0.91 × 2) [0.05^-0.08 - 1.05^-0.08] / 0.08 = 1.412538 × 0.274711 / 0.08 = 4.8505
  // aftershocks of M ≥ 5 in the first day after an M 7, so P = 1 - e^-4.8505 = 0.9922.
  check('aftershocks: the expected count of M ≥ 5 in the first day after an M 7 (by hand, 4.8505)', expectedAftershocks(0, 1, 5, 7) - 4.8505, 5e-5);
  // Wells and Coppersmith (1994), Table 2A, subsurface rupture length, all slip types: M 7 → 10^1.69 = 49.0 km.
  check('aftershocks: the subsurface rupture length of an M 7 is 49.0 km', subsurfaceRuptureLengthKm(7) - 48.98, 0.01, 'km');
}

// The C++.
cppChecks({ assert, skip, lines }, { project: resolve(dirname(fileURLToPath(import.meta.url)), '..'), name: 'geophysics' });

done();
