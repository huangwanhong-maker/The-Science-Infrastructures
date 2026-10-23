// The template's head, as the signals need it: the scalp cast from the T1 image, the
// electrodes of the 10-10 system placed on it by arcs (with the positions below its ring
// that recordings use), the scalp's surface, and the sphere that stands for the head in
// the forward models, with the radial normalisation that maps the head onto it. Built
// once, on first use, and shared by the Head and Signals modes.

import { loadVolume, FILES } from '../data/anatomy.js';
import { headMask, LANDMARKS, CENTRE } from '../modes/head.js';
import { maskScalp, tenTwenty, electrodeMap, radialNormaliser, sphericalSpline } from '../../../../domains/neuroscience/js/index.js';

let worker = null, requests = 0;
const pending = new Map();
export function surface(message) {
  if (!worker) {
    worker = new Worker(new URL('../workers/surface.worker.js', import.meta.url), { type: 'module' });
    worker.onmessage = (e) => { const p = pending.get(e.data.id); if (!p) return; pending.delete(e.data.id); if (e.data.error) p.reject(new Error(e.data.error)); else p.resolve(e.data); };
  }
  const id = ++requests;
  return new Promise((resolve, reject) => { pending.set(id, { resolve, reject }); worker.postMessage({ id, ...message }); });
}

let loading = null;

// {mask, scalp, landmarks, electrodes (Map name -> mm), mesh, centre, radius, normalise}
export function loadHeadModel() {
  if (!loading) loading = build();
  return loading;
}

async function build() {
  const mask = headMask(await loadVolume(FILES.head));
  const scalp = maskScalp(mask, { centre: CENTRE, axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], step: 2, maxRadius: 230 });
  const landmarks = {};
  for (const [k, p] of Object.entries(LANDMARKS)) landmarks[k] = scalp.onSurface(p);
  const placed = tenTwenty(scalp, landmarks, { system: '10-10', samples: 240 });
  const electrodes = electrodeMap(placed);
  const mesh = await surface({ kind: 'mask', dims: mask.dims, affine: mask.affine, data: mask.data, step: 2, smoothing: 4 });
  // The sphere: centred where the scalp is, of the median radius of the electrodes above the ring.
  const radii = placed.positions.filter((p) => !['Nz', 'Iz', 'LPA', 'RPA'].includes(p.name)).map((p) => Math.hypot(...p.position.map((x, a) => x - scalp.centre[a]))).sort((a, b) => a - b);
  const radius = radii[Math.floor(radii.length / 2)];
  return { mask, scalp, landmarks, placed, electrodes, mesh, centre: scalp.centre, radius, normalise: radialNormaliser(scalp, radius) };
}

// A spherical spline over the scalp (Perrin et al. 1989) for values at electrodes, quickly
// evaluated at many directions by a table of its kernel g. `directions` are the
// electrodes' directions from the centre.
export class ScalpSpline {
  constructor(directions) {
    this.directions = directions.map((d) => { const n = Math.hypot(...d); return d.map((x) => x / n); });
    this.table = null;
  }

  fit(values) {
    this.spline = sphericalSpline(this.directions, values, { regularisation: 1e-8 });
    if (!this.table) {
      const N = 2048;
      this.table = new Float64Array(N + 1);
      for (let i = 0; i <= N; i++) this.table[i] = this.spline.g(-1 + (2 * i) / N);
    }
    return this;
  }

  // The value in a direction (need not be a unit vector).
  at(d) {
    const n = Math.hypot(d[0], d[1], d[2]) || 1, u0 = d[0] / n, u1 = d[1] / n, u2 = d[2] / n;
    const c = this.spline.coefficients, e = this.directions, g = this.table, N = g.length - 1;
    let s = this.spline.constant;
    for (let i = 0; i < e.length; i++) {
      const t = u0 * e[i][0] + u1 * e[i][1] + u2 * e[i][2];
      const x = (Math.max(-1, Math.min(1, t)) + 1) * 0.5 * N, k = Math.min(N - 1, x | 0), f = x - k;
      s += c[i] * (g[k] + (g[k + 1] - g[k]) * f);
    }
    return s;
  }
}
