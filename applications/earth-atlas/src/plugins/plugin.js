// Plug-ins: algorithms written as scripts, loaded on request, run on the
// atlas's data and shown through its layers (ADR 0027). A script exports
//
//     export function create(api) { class Mine extends api.Diagnostic { ... } return new Mine(); }
//
// where `api` is ./index.js: these classes, the result makers and a few
// numerical helpers. The classes are named by what the algorithm is, and each
// writes run() from smaller hooks that a subclass overrides:
//   Plugin            identity, parameters, needs and run(): the general case;
//   Diagnostic        a field computed from fields, at every time of the input;
//   ForecastModel     a state advanced in time from an analysis, verified against later analyses;
//   DataAssimilation  a background and observations to an analysis and its increment;
//   Trajectories      parcels carried through a flow.
// Results are data, never drawing: a field product is drawn, contoured, compared
// and read like any product of the library; a sketch (polylines, markers,
// labels and moving points in longitude and latitude) is drawn as it is; tables
// and notes are shown in the panel. Nothing here touches the page or three.js,
// so a plug-in can be run and checked in Node (see tests/run-tests.mjs).

import { ComputedField, statistics } from '../data/fields.js';
import { bilinear } from '../../../../foundations/physis/js/index.js';
import { rmsDifference } from '../../../../domains/meteorology/js/index.js';
import { wrapLongitude, EARTH_RADIUS_M } from '../../../../domains/geodesy/js/index.js';

const HOUR = 3600000;
export const iso = (ms) => new Date(ms).toISOString().replace('.000Z', 'Z');

// ---- results -------------------------------------------------------------------------

// A field product to show. `show` names what to select once it is in the
// catalogue, by variable or vector name: { surface: 'zeta@500hPa', flow: 'wind@500hPa', lines: 'z@500hPa' }.
export const field = (product, show = {}) => ({ kind: 'field', product, show });
// A sketch over the globe: { title, polylines: [{ coords: [lon, lat, ...], colour, width, alpha, dashed }],
// markers: [{ lon, lat, colour, size }], labels: [{ lon, lat, text }],
// moving: [{ times: [ms], coords: [lon, lat, ...], colour, size }] }; a moving
// point is drawn where it is at the time shown.
export const sketch = (spec) => ({ kind: 'sketch', title: '', polylines: [], markers: [], labels: [], moving: [], ...spec });
// A table for the panel: headings, then rows of cells, the first cell of each a heading.
export const table = (title, columns, rows, note = '') => ({ kind: 'table', title, columns, rows, note });
export const note = (text) => ({ kind: 'note', text });

// ---- products and numerics ----------------------------------------------------------------

// A field product made in the page, with the interface of one read from the library.
// variables: [{ name, long_name, standard_name, units, level?, slices: [Float32Array per time] }];
// vectors: [{ name, u, v, units, standard_name, level? }]; times in ms or ISO.
export function makeProduct({ name, title, model = title, domain = 'meteorology', grid, times, variables, vectors = [], source = {}, verification = null }) {
  const slices = {};
  const described = variables.map(({ slices: list, ...v }) => { slices[v.name] = list; return { ...v, ...statistics(list) }; });
  const meta = {
    schema: 'field-grid/1', name, computed: true, title, model, domain, grid,
    times: times.map((t) => (typeof t === 'number' ? iso(t) : t)), time_scale: 'UTC',
    variables: described, vectors, source: { provider: 'computed in the page', credit: title, ...source },
  };
  if (verification) meta.verification = verification;
  return new ComputedField(meta, slices);
}

// A variable's values at one time, sampled bilinearly onto another grid.
export function sampleOnGrid(product, name, index, grid) {
  const s = product.slice(name, index), out = new Float32Array(grid.nlat * grid.nlon);
  for (let r = 0; r < grid.nlat; r++) {
    const lat = grid.lat0 - r * grid.dlat;
    for (let c = 0; c < grid.nlon; c++) out[r * grid.nlon + c] = bilinear(s, product.grid, grid.lon0 + c * grid.dlon, lat);
  }
  return out;
}

// The index of the last time at or before `ms`, or 0.
export const indexAtOrBefore = (times, ms) => Math.max(0, times.findLastIndex((t) => t <= ms + 1000));

// Points spread evenly over the sphere (a Fibonacci lattice): [{ lon, lat }].
export function lattice(n) {
  const golden = Math.PI * (3 - Math.sqrt(5)), out = [];
  for (let i = 0; i < n; i++) {
    const z = 1 - (2 * i + 1) / n, lon = wrapLongitude(((i * golden) * 180) / Math.PI);
    out.push({ lon, lat: (Math.asin(z) * 180) / Math.PI });
  }
  return out;
}

// One fourth-order Runge–Kutta step of a parcel on the sphere: from (lon, lat)
// in degrees at time ms, over dt seconds, through flow(lon, lat, ms) -> [u, v]
// in m/s, or null outside the data. Returns [lon, lat], or null.
export function parcelStep(flow, lon, lat, ms, dt) {
  const D = 180 / Math.PI;
  const rate = (lo, la, t) => {
    const uv = flow(wrapLongitude(lo), Math.max(-89.99, Math.min(89.99, la)), t);
    return uv ? [(uv[0] / (EARTH_RADIUS_M * Math.cos(la / D))) * D, (uv[1] / EARTH_RADIUS_M) * D] : null;
  };
  const k1 = rate(lon, lat, ms); if (!k1) return null;
  const k2 = rate(lon + 0.5 * dt * k1[0], lat + 0.5 * dt * k1[1], ms + 500 * dt); if (!k2) return null;
  const k3 = rate(lon + 0.5 * dt * k2[0], lat + 0.5 * dt * k2[1], ms + 500 * dt); if (!k3) return null;
  const k4 = rate(lon + dt * k3[0], lat + dt * k3[1], ms + 1000 * dt); if (!k4) return null;
  let nlon = lon + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
  let nlat = lat + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  if (nlat > 90) { nlat = 180 - nlat; nlon += 180; } else if (nlat < -90) { nlat = -180 - nlat; nlon += 180; }
  return [wrapLongitude(nlon), nlat];
}

// ---- the classes ---------------------------------------------------------------------------

export class Plugin {
  // Identity: a subclass sets these as class fields or in its constructor.
  id = 'plugin';
  title = 'A plug-in';
  description = '';
  domain = 'meteorology';
  // Parameters the panel offers and passes to run() as `params`:
  // { name: { label, default, min, max, step, unit } } or { name: { label, default, choices: [[value, label]] } }.
  parameters = {};

  // The fields the plug-in needs. Each { name, label, standard_name | standard_names: [u, v],
  // level?, times?, distinct? } is resolved by the host to the products that carry such
  // variables (one candidate per product and level); the user chooses the source, and run()
  // receives inputs[name] = { product, names, variables, meta, model, level }.
  // `distinct` asks for a source other than the needs before it; `times` for at least that many times.
  needs() { return []; }

  // The algorithm: returns a list of results made by field(), sketch(), table() and note().
  async run(inputs, params, context) { return []; } // eslint-disable-line no-unused-vars

  // A resting point in a long computation: reports progress and lets the page breathe.
  async breathe(context, fraction, text = '') {
    context.progress(fraction, text);
    await new Promise((resolve) => setTimeout(resolve, 0));
  }

  // The product name of this plug-in's output; a second run replaces the first.
  get productName() { return `plugin-${this.id}`; }
}

// A field computed from fields, at every time of the first need, on its grid.
// Override outputs(), compute(), and vectors() and show() if wanted.
export class Diagnostic extends Plugin {
  // [{ name, long_name, standard_name, units, level? }]: the variables compute() returns.
  outputs() { return []; }
  // [{ name, u, v, units, standard_name, level? }]: pairs among the outputs that form vectors.
  vectors() { return []; }
  // { output name: Float32Array on the first input's grid }, for time `index` of the first input.
  compute(inputs, params, index, context) { throw new Error(`${this.id}: compute() is not written`); } // eslint-disable-line no-unused-vars
  // What to select when the product is shown.
  show() { const first = this.outputs()[0]; return first ? { surface: first.name } : {}; }

  async run(inputs, params, context) {
    const first = inputs[this.needs()[0].name], product = first.product;
    const outputs = this.outputs(), lists = Object.fromEntries(outputs.map((o) => [o.name, []]));
    for (let k = 0; k < product.times.length; k++) {
      const values = this.compute(inputs, params, k, context);
      for (const o of outputs) lists[o.name].push(values[o.name]);
      await this.breathe(context, (k + 1) / product.times.length, `time ${k + 1} of ${product.times.length}`);
    }
    const made = makeProduct({
      name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain,
      grid: product.grid, times: product.times, vectors: this.vectors(),
      variables: outputs.map((o) => ({ level: first.level, ...o, slices: lists[o.name] })),
      source: { dataset: first.meta.name, provider: `computed in the page from ${first.model}`, credit: `${this.title} (plug-in), from ${first.model}` },
    });
    return [field(made, this.show())];
  }
}

// A state advanced in time from an analysis: the first need gives the start,
// at the input time at or before the atlas's instant. Override grid(), outputs(),
// initial(), advance() and record(); verify() names what to compare with the input.
export class ForecastModel extends Plugin {
  dt = 1800;      // the step, in seconds
  every = 6;      // hours between recorded times
  parameters = { hours: { label: 'Forecast length', default: 72, min: 6, max: 240, step: 6, unit: 'h' } };

  grid() { throw new Error(`${this.id}: grid() is not written`); }
  outputs() { return []; }
  vectors() { return []; }
  // { output: the output variable, need: the need holding the truth, label, unit }, or null for no verification.
  verify() { return null; }
  // The state at time index `index` of the inputs (any object the other hooks understand).
  async initial(inputs, index, params, context) { throw new Error(`${this.id}: initial() is not written`); } // eslint-disable-line no-unused-vars
  // One step of `seconds`.
  advance(state, seconds, params) { throw new Error(`${this.id}: advance() is not written`); } // eslint-disable-line no-unused-vars
  // { output name: Float32Array on grid() }: the state as fields.
  record(state, params) { throw new Error(`${this.id}: record() is not written`); } // eslint-disable-line no-unused-vars
  show() { const first = this.outputs()[0]; return first ? { surface: first.name } : {}; }

  async run(inputs, params, context) {
    const first = inputs[this.needs()[0].name], product = first.product;
    const i0 = indexAtOrBefore(product.times, context.time), t0 = product.times[i0];
    const hours = Number(params.hours ?? 72);
    const steps = Math.round((hours * 3600) / this.dt), recordEvery = Math.max(1, Math.round((this.every * 3600) / this.dt));
    const outputs = this.outputs(), lists = Object.fromEntries(outputs.map((o) => [o.name, []])), times = [t0];
    const state = await this.initial(inputs, i0, params, context);
    const take = () => { const r = this.record(state, params); for (const o of outputs) lists[o.name].push(r[o.name]); };
    take();
    for (let s = 1; s <= steps; s++) {
      this.advance(state, this.dt, params);
      if (s % recordEvery === 0) { take(); times.push(t0 + s * this.dt * 1000); }
      if (s % 4 === 0) await this.breathe(context, s / steps, `+${Math.round((s * this.dt) / 3600)} h`);
    }
    const grid = this.grid();
    // Verification against the input's own later times, 20–80° N, beside persistence.
    const v = this.verify();
    let verification = null, rows = [], truth = null;
    if (v && inputs[v.need] && inputs[v.need].product.times.includes(t0)) {
      truth = inputs[v.need];
      const name = truth.names[0], start = sampleOnGrid(truth.product, name, truth.product.times.indexOf(t0), grid), now = Date.now();
      times.forEach((t, k) => {
        const idx = truth.product.times.indexOf(t);
        if (k === 0 || idx < 0) return;
        const analysed = sampleOnGrid(truth.product, name, idx, grid);
        rows.push({ hours: (t - t0) / HOUR, model: rmsDifference(lists[v.output][k], analysed, grid, 20, 80),
          persistence: rmsDifference(start, analysed, grid, 20, 80), against: t <= now ? 'analysis' : 'forecast' });
      });
      verification = { quantity: v.label, against: truth.model, start: iso(t0), leads: rows };
    }
    const made = makeProduct({
      name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain,
      grid, times, vectors: this.vectors(), variables: outputs.map((o) => ({ ...o, slices: lists[o.name] })), verification,
      source: { dataset: first.meta.name, provider: `computed in the page from ${first.model} at ${iso(t0)}`, credit: `${this.title} (plug-in), from ${first.model}` },
    });
    const results = [field(made, this.show())];
    if (rows.length) {
      const unit = v.unit ? ` ${v.unit}` : '';
      const shown = rows.length > 8 ? rows.filter((r) => r.hours % 24 === 0) : rows;   // whole days, unless the run is short
      results.push(table(`Root-mean-square error of ${v.label}, 20–80° N, against ${truth.model}`, ['Lead', this.title, 'Persistence'],
        shown.map((r) => [`+${r.hours} h`, `${Math.round(r.model)}${unit}`, `${Math.round(r.persistence)}${unit}`]),
        rows.some((r) => r.against === 'forecast') ? 'Beyond the present, against the forecast rather than analyses.' : 'Against analyses.'));
    }
    return results;
  }
}

// A background field and observations to an analysis. The first need is the
// background, taken at the input time at or before the atlas's instant.
// Override observations() and analyse(); the result is the analysis and its
// increment as a product, with the observations sketched.
export class DataAssimilation extends Plugin {
  // [{ lon, lat, value, error? }]: what is observed, at time index `index` of the inputs.
  async observations(inputs, index, params, context) { return []; } // eslint-disable-line no-unused-vars
  // Float32Array: the analysis on the background's grid (NaN where the background is missing).
  async analyse(background, grid, observations, params, context) { return background; } // eslint-disable-line no-unused-vars
  // The markers drawn for the observations.
  markers(observations, innovations) {
    return observations.map((o, i) => ({ lon: o.lon, lat: o.lat, size: 2.6, colour: innovations[i] > 0 ? '#f6a07e' : '#8fb8e8' }));
  }

  async run(inputs, params, context) {
    const first = inputs[this.needs()[0].name], product = first.product, name = first.names[0], v = first.variables[0];
    const index = indexAtOrBefore(product.times, context.time), t = product.times[index];
    const background = product.slice(name, index);
    const observations = await this.observations(inputs, index, params, context);
    const innovations = observations.map((o) => o.value - product.sample(name, t, o.lon, o.lat));
    const analysis = await this.analyse(background, product.grid, observations, params, context);
    const increment = new Float32Array(analysis.length);
    for (let i = 0; i < analysis.length; i++) increment[i] = analysis[i] - background[i];
    const made = makeProduct({
      name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain,
      grid: product.grid, times: [t],
      variables: [
        { name: v.name, long_name: `${v.long_name || v.name}, analysed`, standard_name: v.standard_name, units: v.units, level: v.level, slices: [analysis] },
        { name: `${v.name}_increment`, long_name: `Analysis increment of ${(v.long_name || v.name).toLowerCase()}`, units: v.units, level: v.level, difference: true, slices: [increment] },
      ],
      source: { dataset: first.meta.name, provider: `analysed in the page from ${first.model} at ${iso(t)} and ${observations.length} observations`, credit: `${this.title} (plug-in), from ${first.model}` },
    });
    return [
      field(made, { surface: `${v.name}_increment` }),
      sketch({ title: `${observations.length} observations (warm where above the background)`, markers: this.markers(observations, innovations) }),
    ];
  }
}

// Parcels carried through a flow: the first need is a vector field. Override
// seeds(); velocity() samples the need by default. The result is a sketch of
// the paths, with the parcels moving along them as the time shown changes.
export class Trajectories extends Plugin {
  parameters = {
    hours: { label: 'Hours ahead', default: 48, min: 6, max: 240, step: 6, unit: 'h' },
    step: { label: 'Time step', default: 1, min: 0.25, max: 6, step: 0.25, unit: 'h' },
  };
  colour = '#e8c3dc';

  // [{ lon, lat, colour? }]: where the parcels start.
  async seeds(inputs, params, context) { return []; } // eslint-disable-line no-unused-vars
  // (lon, lat, ms) -> [u, v] in m/s, or null outside the data.
  velocity(inputs, params) { // eslint-disable-line no-unused-vars
    const w = inputs[this.needs()[0].name], [u, v] = w.names;
    return (lon, lat, ms) => {
      const a = w.product.sample(u, ms, lon, lat), b = w.product.sample(v, ms, lon, lat);
      return a === a && b === b ? [a, b] : null;
    };
  }

  async run(inputs, params, context) {
    const w = inputs[this.needs()[0].name], times = w.product.times, flow = this.velocity(inputs, params);
    const tEnd = times[times.length - 1], t0 = Math.max(times[0], Math.min(context.time, tEnd));
    const dt = Number(params.step ?? 1) * 3600;
    const steps = Math.max(0, Math.min(Math.round((Number(params.hours ?? 48) * 3600) / dt), Math.floor((tEnd - t0) / 1000 / dt)));
    const seeds = await this.seeds(inputs, params, context);
    const paths = [];
    for (const [s, seed] of seeds.entries()) {
      const coords = [seed.lon, seed.lat], when = [t0];
      let lon = seed.lon, lat = seed.lat;
      for (let k = 0; k < steps; k++) {
        const next = parcelStep(flow, lon, lat, t0 + k * dt * 1000, dt);
        if (!next) break;
        [lon, lat] = next;
        coords.push(lon, lat); when.push(t0 + (k + 1) * dt * 1000);
      }
      paths.push({ coords, times: when, colour: seed.colour || this.colour });
      if (s % 8 === 7) await this.breathe(context, (s + 1) / seeds.length, `parcel ${s + 1} of ${seeds.length}`);
    }
    return [sketch({
      title: `${paths.length} parcels through the wind of ${w.model}, from ${iso(t0)} for ${Math.round((steps * dt) / 3600)} h`,
      polylines: paths.map((p) => ({ coords: p.coords, colour: p.colour, width: 1.3, alpha: 0.85 })),
      moving: paths.map((p) => ({ times: p.times, coords: p.coords, colour: '#f6c47e', size: 3 })),
    })];
  }
}
