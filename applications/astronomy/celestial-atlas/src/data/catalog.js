// Loading of the star catalogue and the constellation figures prepared by data/.
//
// stars-v1 is built by data/products/star_catalog.py from the HYG database
// (CC BY-SA 4.0); constellations come from d3-celestial (BSD 3-Clause).

import { fromSpherical, DEG } from '../core/frames.js';
import {
  supergalactic, tidalIndices, tidalGroups, inferDisturberLuminosities, displayName,
} from '../core/galaxies.js';

async function json(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return r.json();
}

async function binary(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return r.arrayBuffer();
}

// The star catalogues built by data/products/star_catalog.py.
export async function listStarCatalogues(base) {
  try {
    return (await json(`${base}/products/index.json`)).stars;
  } catch {
    return [{ name: 'stars-v1', title: 'HYG v4.1 (Hipparcos, Yale, Gliese)', photometry: 'johnson' }];
  }
}

export async function loadStars(base, name = 'stars-v1') {
  const meta = await json(`${base}/products/${name}.json`);
  const [buffer, idBuffer] = await Promise.all([
    binary(`${base}/products/${name}.bin`),
    meta.ids ? binary(`${base}/products/${meta.ids}`) : null,
  ]);
  const ids = idBuffer ? new BigUint64Array(idBuffer) : null;
  const record = new Float32Array(buffer);
  const n = meta.count;
  const stride = meta.fields.length;
  const unit = new Float32Array(3 * n);   // ICRS unit vectors
  const xyz = new Float32Array(3 * n);    // equatorial position, pc (NaN if unknown)
  const mag = new Float32Array(n);
  const ci = new Float32Array(n);
  const dist = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const o = i * stride;
    const ra = record[o], dec = record[o + 1];
    const c = Math.cos(dec);
    unit[3 * i] = c * Math.cos(ra);
    unit[3 * i + 1] = c * Math.sin(ra);
    unit[3 * i + 2] = Math.sin(dec);
    mag[i] = record[o + 2];
    ci[i] = record[o + 3];
    dist[i] = record[o + 4];
    xyz[3 * i] = record[o + 5];
    xyz[3 * i + 1] = record[o + 6];
    xyz[3 * i + 2] = record[o + 7];
  }
  const names = new Map(meta.names);
  const designations = new Map(meta.designations);
  const spectra = new Map(meta.spectra);

  // Indices sorted from brightest to faintest, for magnitude-limited drawing.
  const byBrightness = Uint32Array.from({ length: n }, (_, i) => i)
    .sort((a, b) => mag[a] - mag[b]);

  const designation = (i) => designations.get(i) || (ids ? `Gaia DR3 ${ids[i]}` : '');
  return {
    name: meta.name,
    title: meta.title,
    photometry: meta.photometry || 'johnson',
    epoch: meta.epoch || 'J2000',
    count: n, unit, xyz, mag, ci, dist, byBrightness,
    starName: (i) => names.get(i) || '',
    designation,
    spectrum: (i) => spectra.get(i) || '',
    label: (i) => names.get(i) || designations.get(i) || '',
    namedIndices: [...names.keys()],
    source: meta.source,
  };
}

// ---- galaxies -----------------------------------------------------------------

// The galaxy catalogues built by data/products/galaxy_catalog.py.
export async function listGalaxyCatalogues(base) {
  try {
    return (await json(`${base}/products/index.json`)).galaxies || [];
  } catch {
    return [];
  }
}

const nanOf = (v) => (v === null || v === undefined ? NaN : v);

// A Local Volume catalogue with tidal relations: published (Karachentsev et
// al. 2013) and recomputed here from positions and luminosities.
async function loadTidal(base, name) {
  const meta = await json(`${base}/products/${name}.json`);
  const c = meta.columns;
  const n = meta.count;
  const positions = new Float32Array(3 * n);
  const exact = new Float64Array(3 * n);
  for (let i = 0; i < n; i++) exact.set(supergalactic(c.ra[i], c.dec[i], c.dist[i]), 3 * i);
  positions.set(exact);
  const logL = Float64Array.from(c.logL, nanOf);
  const byName = new Map(c.name.map((v, i) => [v, i]));
  const published = Int32Array.from(c.disturber, (d) => (byName.has(d) ? byName.get(d) : -1));
  const publishedTheta = Float64Array.from(c.theta1, nanOf);
  const filled = inferDisturberLuminosities(exact, logL, published, publishedTheta);
  const ours = tidalIndices(exact, filled.logL);
  const publishedGroups = tidalGroups(published, publishedTheta, filled.logL);
  const oursGroups = tidalGroups(ours.disturber, ours.theta1, filled.logL);
  let comparable = 0, agree = 0;
  for (let i = 0; i < n; i++) {
    if (published[i] < 0 || !Number.isFinite(publishedTheta[i])) continue;
    comparable++;
    if (ours.disturber[i] === published[i]) agree++;
  }
  return {
    name, kind: 'tidal', title: meta.title, units: 'Mpc', distances: 'measured', count: n,
    positions, exact,
    logL: Float64Array.from(filled.logL),
    catalogueLogL: logL,
    names: c.name.map(displayName),
    rawNames: c.name,
    dist: Float64Array.from(c.dist, nanOf),
    distMethod: c.distMethod,
    vLG: Float64Array.from(c.vLG, nanOf),
    type: Float64Array.from(c.type, nanOf),
    dwarfType: c.dwarfType,
    publishedDisturberName: c.disturber,
    relations: {
      published: { disturber: published, theta1: publishedTheta, groups: publishedGroups,
        label: 'Karachentsev et al. 2013' },
      ours: { disturber: ours.disturber, theta1: ours.theta1, theta5: ours.theta5, groups: oursGroups,
        label: 'recomputed here' },
    },
    inferred: filled.inferred,
    agreement: { comparable, agree },
    source: meta.source,
  };
}

// The 2MRS group catalogue: galaxies with their group membership.
async function loadGroups(base, name) {
  const [meta, galaxies, groupTable] = await Promise.all([
    json(`${base}/products/${name}.json`),
    binary(`${base}/products/${name}.bin`),
    binary(`${base}/products/${name}.groups.bin`),
  ]);
  const g = new Float32Array(galaxies);
  const t = new Float32Array(groupTable);
  const n = meta.count, gn = meta.groups;
  const fs = meta.fields.length, gs = meta.group_fields.length;
  const positions = new Float32Array(3 * n);
  const logL = new Float64Array(n), vCMB = new Float64Array(n), dist = new Float64Array(n);
  const type = new Float64Array(n), kmag = new Float64Array(n);
  const pgc = new Int32Array(n), groupOf = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    const o = i * fs;
    positions[3 * i] = g[o]; positions[3 * i + 1] = g[o + 1]; positions[3 * i + 2] = g[o + 2];
    logL[i] = g[o + 3]; vCMB[i] = g[o + 4]; dist[i] = g[o + 5]; type[i] = g[o + 6]; kmag[i] = g[o + 7];
    pgc[i] = g[o + 8]; groupOf[i] = g[o + 9];
  }
  const field = (k) => Float64Array.from({ length: gn }, (_, j) => t[j * gs + k]);
  const groups = {
    centre: Float32Array.from({ length: 3 * gn }, (_, j) => t[Math.floor(j / 3) * gs + (j % 3)]),
    dist: field(3), logLK: field(4), mLum: field(5), mVir: field(6), sigmaV: field(7), r2t: field(8),
    nHere: Int32Array.from(field(9)), nCatalog: Int32Array.from(field(10)),
    nest: Int32Array.from(field(11)), pgc1: Int32Array.from(field(12)), brightest: Int32Array.from(field(13)),
    names: new Map(meta.groupNames),
    count: gn,
  };
  const names = new Map(Object.entries(meta.names).map(([k, v]) => [Number(k), v]));
  return {
    name, kind: 'groups', title: meta.title, units: 'h⁻¹ Mpc', distances: 'kinematic', count: n,
    positions, logL, vCMB, dist, type, kmag, pgc, groupOf, groups,
    galaxyName: (i) => names.get(i) || `PGC ${pgc[i]}`,
    groupName: (j) => groups.names.get(j) || (groups.nHere[j] > 1 ? `Group ${groups.nest[j]}` : ''),
    source: meta.source,
  };
}

export async function loadGalaxies(base, entry) {
  return entry.kind === 'tidal' ? loadTidal(base, entry.name) : loadGroups(base, entry.name);
}

export async function loadConstellations(base) {
  const dir = `${base}/d3-celestial-constellations`;
  const [lines, names] = await Promise.all([
    json(`${dir}/constellations.lines.json`),
    json(`${dir}/constellations.json`),
  ]);
  const segments = [];
  for (const f of lines.features) {
    for (const line of f.geometry.coordinates) {
      for (let k = 1; k < line.length; k++) {
        const a = fromSpherical(line[k - 1][0] * DEG, line[k - 1][1] * DEG);
        const b = fromSpherical(line[k][0] * DEG, line[k][1] * DEG);
        segments.push(...a, ...b);
      }
    }
  }
  const labels = names.features.map((f) => ({
    id: f.id,
    rank: Number(f.properties.rank),
    names: { en: f.properties.en || f.properties.name, la: f.properties.la || f.properties.name,
      zh: f.properties.zh, ja: f.properties.ja },
    unit: fromSpherical(f.geometry.coordinates[0] * DEG, f.geometry.coordinates[1] * DEG),
  }));
  return { segments: new Float32Array(segments), labels };
}
