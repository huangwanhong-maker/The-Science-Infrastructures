// The source space: the places on the cortex where currents are estimated. The cortical
// surface of the anatomy product (its outer edge, two hemispheres) is thinned to one
// vertex a cell of a 6 mm grid; each source keeps the surface's normal there, the
// orientation of the pyramidal cells whose currents EEG and MEG see, and the gyrus of
// Harvard-Oxford that its vertex lies on. Every vertex of the drawn surface takes its
// value from the three nearest sources of its hemisphere, weighted by inverse distance,
// and every source knows its neighbours within 12 mm, for finding the maxima of a map.

import { decimate, structureOfKey } from '../../../../domains/neuroscience/js/index.js';

export function buildSourceSpace(instances, { spacing = 6 } = {}) {
  const positions = [], normals = [], where = [];        // where: [instance index, vertex]
  instances.forEach((inst, ii) => {
    const pos = inst.mesh.geometry.attributes.position.array, nrm = inst.mesh.geometry.attributes.normal.array;
    const n = pos.length / 3;
    const pts = [];
    for (let v = 0; v < n; v++) pts.push([pos[3 * v], pos[3 * v + 1], pos[3 * v + 2]]);
    for (const v of decimate(pts, spacing)) {
      positions.push(pts[v]);
      const nv = [nrm[3 * v], nrm[3 * v + 1], nrm[3 * v + 2]], l = Math.hypot(...nv) || 1;
      normals.push(nv.map((x) => x / l));
      where.push([ii, v]);
    }
  });
  const S = positions.length;
  const side = where.map(([ii]) => instances[ii].side);
  const structure = where.map(([ii, v]) => {
    const inst = instances[ii];
    const label = inst.labels ? inst.labels[v] : 0;
    return label ? structureOfKey(`ho:${label}`) : null;
  });
  // A hash of the sources by cells of 12 mm, for neighbours and for the vertices' weights.
  const cell = 12;
  const grid = new Map();
  const keyOf = (p) => `${Math.floor(p[0] / cell)},${Math.floor(p[1] / cell)},${Math.floor(p[2] / cell)}`;
  positions.forEach((p, s) => { const k = keyOf(p); if (!grid.has(k)) grid.set(k, []); grid.get(k).push(s); });
  const near = (p) => {
    const c = p.map((x) => Math.floor(x / cell)), out = [];
    for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) for (let d = -1; d <= 1; d++) {
      const list = grid.get(`${c[0] + a},${c[1] + b},${c[2] + d}`);
      if (list) out.push(...list);
    }
    return out;
  };
  const neighbours = positions.map((p, s) => near(p).filter((t) => t !== s && side[t] === side[s] && Math.hypot(...p.map((x, i) => x - positions[t][i])) <= 12));
  // The weights of every vertex of each instance.
  const weights = instances.map((inst) => {
    const pos = inst.mesh.geometry.attributes.position.array, n = pos.length / 3;
    const idx = new Int32Array(3 * n), w = new Float32Array(3 * n);
    for (let v = 0; v < n; v++) {
      const p = [pos[3 * v], pos[3 * v + 1], pos[3 * v + 2]];
      const cand = near(p).filter((t) => side[t] === inst.side);
      const best = [[-1, Infinity], [-1, Infinity], [-1, Infinity]];
      for (const t of cand) {
        const q = positions[t], d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2;
        if (d < best[2][1]) { best[2] = [t, d]; best.sort((x, y) => x[1] - y[1]); }
      }
      let total = 0;
      for (let k = 0; k < 3; k++) {
        const [t, d] = best[k];
        idx[3 * v + k] = t;
        const wk = t < 0 ? 0 : 1 / (Math.sqrt(d) + 1.5);
        w[3 * v + k] = wk; total += wk;
      }
      for (let k = 0; k < 3; k++) w[3 * v + k] /= total || 1;
    }
    return { idx, w };
  });
  return { positions, normals, side, structure, neighbours, weights, count: S, instances };
}

// The values of the sources spread over the vertices of one instance.
export function spread(space, i, values, out) {
  const { idx, w } = space.weights[i];
  const n = idx.length / 3;
  for (let v = 0; v < n; v++) {
    const a = idx[3 * v], b = idx[3 * v + 1], c = idx[3 * v + 2];
    out[v] = (a >= 0 ? w[3 * v] * values[a] : 0) + (b >= 0 ? w[3 * v + 1] * values[b] : 0) + (c >= 0 ? w[3 * v + 2] * values[c] : 0);
  }
  return out;
}

// The local maxima of |values| over the sources: those larger than every neighbour, the
// strongest first.
export function maxima(space, values, { count = 5, floor = 0 } = {}) {
  const out = [];
  for (let s = 0; s < space.count; s++) {
    const v = Math.abs(values[s]);
    if (v <= floor) continue;
    if (space.neighbours[s].every((t) => Math.abs(values[t]) <= v)) out.push(s);
  }
  out.sort((a, b) => Math.abs(values[b]) - Math.abs(values[a]));
  return out.slice(0, count);
}
