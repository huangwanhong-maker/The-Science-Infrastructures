// The laboratory's own forecasts, run in the page from the 500 hPa analysis of
// a forecast product and verified against that product's heights:
// - the barotropic vorticity model of 1950, from the wind;
// - the shallow-water model, from the wind and the height, started by a
//   digital filter.
// Each run becomes a product of its own, so the two can be compared with each
// other and with the forecast they start from.
//
// The model's output becomes a field product like any other (a ComputedField),
// so the atlas draws it, contours it, carries particles with it, and compares
// it with the other models, without code of its own.

import { bilinear, ComputedField, statistics } from '../core/grid.js';
import { gaussLegendre } from '../core/spectral.js';
import { rmsDifference } from '../core/barotropic.js';

export const LAB_MODELS = {
  barotropic: { label: 'Barotropic (laboratory)', title: 'Barotropic vorticity model (laboratory)', short: 'Barotropic', dt: 1800 },
  'shallow-water': { label: 'Shallow water (laboratory)', title: 'Shallow-water model (laboratory)', short: 'Shallow water', dt: 900 },
};
export const LAB_MODEL = LAB_MODELS.barotropic.label;
const OUTPUT = { lon0: -178.75, dlon: 2.5, nlon: 144, lat0: 88.75, dlat: 2.5, nlat: 72, periodic: true };
const TRUNCATION = 42;

// The product and variables that can start the model: geopotential height and
// wind at 500 hPa, preferring the GFS.
export function modelInputs(library) {
  const found = [];
  for (const [name, meta] of library.fieldMeta) {
    if (meta.computed) continue;
    const at500 = (std) => meta.variables.find((v) => v.standard_name === std && v.level === '500 hPa' && !/realization/.test(v.cell_methods || ''));
    const z = at500('geopotential_height'), u = at500('eastward_wind'), v = at500('northward_wind');
    if (z && u && v && meta.times && meta.times.length > 1) found.push({ product: name, model: meta.model || meta.title, z: z.name, u: u.name, v: v.name });
  }
  found.sort((a, b) => (a.model === 'GFS' ? -1 : 0) - (b.model === 'GFS' ? -1 : 0));
  return found[0] || null;
}

export class LabModel {
  constructor(library) {
    this.library = library;
    this.running = false;
    this.result = null;
  }

  // Runs a model ('barotropic' or 'shallow-water') from the product time at or
  // before `time`; returns { product, verification, start }.
  async run(time, { model: kind = 'barotropic', hours = 72, onProgress = () => {} } = {}) {
    const spec = LAB_MODELS[kind];
    const inputs = modelInputs(this.library);
    if (!inputs) throw new Error('no product with wind and height at 500 hPa');
    const product = await this.library.field(inputs.product);
    const t0 = [...product.times].reverse().find((t) => t <= time + 1000) ?? product.times[0];
    const i0 = product.times.indexOf(t0);
    const lastTime = product.times[product.times.length - 1];
    const span = Math.min(hours, Math.floor((lastTime - t0) / 3600000 / 6) * 6 || hours);

    // The 500 hPa wind at the start, on the model's Gaussian grid.
    const nlat = 2 * Math.ceil((3 * TRUNCATION + 1) / 4), nlon = 2 * Math.ceil((3 * TRUNCATION + 1) / 2);
    const { mu } = gaussLegendre(nlat);
    const u = new Float64Array(nlat * nlon), v = new Float64Array(nlat * nlon), z = new Float64Array(nlat * nlon);
    const su = product.slice(inputs.u, i0), sv = product.slice(inputs.v, i0), sz = product.slice(inputs.z, i0);
    for (let j = 0; j < nlat; j++) {
      const lat = (Math.asin(mu[j]) * 180) / Math.PI;
      for (let k = 0; k < nlon; k++) {
        const lon = ((360 * k) / nlon + 180) % 360 - 180;
        u[j * nlon + k] = bilinear(su, product.grid, lon, lat) || 0;
        v[j * nlon + k] = bilinear(sv, product.grid, lon, lat) || 0;
        z[j * nlon + k] = bilinear(sz, product.grid, lon, lat);
      }
    }

    this.running = true;
    const out = await new Promise((resolve, reject) => {
      const worker = new Worker(new URL('../workers/labmodel.worker.js', import.meta.url), { type: 'module' });
      worker.onmessage = (e) => {
        if (e.data.progress !== undefined) onProgress(e.data.progress);
        if (e.data.done) { worker.terminate(); resolve(e.data); }
      };
      worker.onerror = (e) => { worker.terminate(); reject(new Error(e.message || 'the model stopped')); };
      worker.postMessage({ model: kind, u, v, z, hours: span, dt: spec.dt, every: 6, grid: OUTPUT, truncation: TRUNCATION }, [u.buffer, v.buffer, z.buffer]);
    });
    this.running = false;

    // Heights: for the barotropic model, the analysis at the start plus the
    // model's balanced change; the shallow-water model carries its own.
    const zStart = regrid(product, inputs.z, i0);
    const heights = out.fields.z || out.fields.dz.map((dz) => dz.map((x, i) => x + zStart[i]));
    const times = out.times.map((s) => new Date(t0 + s * 1000).toISOString().replace('.000Z', 'Z'));
    const name = `lab-${kind}`;
    const variables = [
      { name: 'z@500hPa', long_name: 'Geopotential height at 500 hPa', standard_name: 'geopotential_height', units: 'm', level: '500 hPa', ...statistics(heights) },
      { name: 'u@500hPa', long_name: 'Eastward wind at 500 hPa', standard_name: 'eastward_wind', units: 'm s-1', level: '500 hPa', ...statistics(out.fields.u) },
      { name: 'v@500hPa', long_name: 'Northward wind at 500 hPa', standard_name: 'northward_wind', units: 'm s-1', level: '500 hPa', ...statistics(out.fields.v) },
      { name: 'zeta@500hPa', long_name: 'Relative vorticity at 500 hPa', standard_name: 'atmosphere_relative_vorticity', units: 's-1', level: '500 hPa', ...statistics(out.fields.zeta) },
    ];
    const meta = {
      schema: 'field-grid/1', name, computed: true,
      title: spec.title, model: spec.label, domain: 'meteorology',
      grid: OUTPUT, times, time_scale: 'UTC', variables,
      vectors: [{ name: 'wind@500hPa', u: 'u@500hPa', v: 'v@500hPa', units: 'm s-1', standard_name: 'wind', level: '500 hPa' }],
      source: { dataset: inputs.product, provider: `computed in the page from ${inputs.model} at ${times[0]}`, credit: `${spec.short} model T${TRUNCATION} (laboratory), from ${inputs.model}` },
    };
    const computed = new ComputedField(meta, { 'z@500hPa': heights, 'u@500hPa': out.fields.u, 'v@500hPa': out.fields.v, 'zeta@500hPa': out.fields.zeta });

    // Verification against the product's own heights, 20-80 N, with persistence for reference.
    const verification = [];
    const now = Date.now();
    out.times.forEach((s, k) => {
      if (k === 0) return;
      const t = t0 + s * 1000, idx = product.times.indexOf(t);
      if (idx < 0) return;
      const truth = regrid(product, inputs.z, idx);
      verification.push({
        hours: s / 3600,
        model: rmsDifference(heights[k], truth, OUTPUT, 20, 80),
        persistence: rmsDifference(zStart, truth, OUTPUT, 20, 80),
        against: t <= now ? 'analysis' : 'forecast',
      });
    });
    this.result = { kind, product: computed, verification, start: t0, against: inputs.model, hours: span };
    return this.result;
  }
}

function regrid(product, name, index) {
  const s = product.slice(name, index), out = new Float32Array(OUTPUT.nlat * OUTPUT.nlon);
  for (let r = 0; r < OUTPUT.nlat; r++) {
    const lat = OUTPUT.lat0 - r * OUTPUT.dlat;
    for (let c = 0; c < OUTPUT.nlon; c++) out[r * OUTPUT.nlon + c] = bilinear(s, product.grid, OUTPUT.lon0 + c * OUTPUT.dlon, lat);
  }
  return out;
}
