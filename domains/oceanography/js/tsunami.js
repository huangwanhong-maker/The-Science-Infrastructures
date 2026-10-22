// Tsunami travel times: the first arrival of the long wave from a source,
// over the bathymetry, by a shortest-path method on the grid.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).
//
// A tsunami is a long wave: its wavelength (100-500 km) is long beside the
// depth h, so it travels at c = √(g h), whatever its period (Lamb 1932, §170).
// The first-arrival time T from a source then satisfies the eikonal equation
// |∇T| = 1 / c on the sphere, and T(x) is the least, over paths from the source
// to x through the sea, of the integral of the slowness 1 / c along the path
// (Fermat's principle). Here the sea is a graph: each wet grid point is joined
// to the wet points of a stencil around it by great-circle steps, each costing
// its length times the mean of the slownesses at its ends, and T is the
// shortest-path distance from the source, found by Dijkstra's (1959) algorithm.
// This is the method of the travel-time charts of Shokin et al. (1987) and of
// Wessel's (2009) TTT, with the stencil adapted to latitude (below).
//
// The stencil. On cells of equal sides, the 16 directions (0, ±1), (±1, 0),
// (±1, ±1), (±1, ±2), (±2, ±1) leave gaps of at most 26.6° between them, so a
// shortest grid path is at most 1/cos(13.3°) - 1 = 2.75 % longer than the
// straight one. On a longitude-latitude grid the cells narrow toward the poles
// (east-west side cos φ times the north-south one), and the gaps open: 8 % at
// 60°. The column offsets of each row's stencil are therefore multiplied by
// 1/cos φ and rounded, with the unit steps (0, ±1) kept: at most 4.2 % at any
// latitude. Adding the moves (±1, ±3), (±3, ±1), (±2, ±3), (±3, ±2) (32
// neighbours) brings it to 2.1 %. These are bounds for the worst direction; the
// mean excess over all directions is about a third of them. Land is
// impassable: a step is not taken if it passes through a land cell, or
// through a corner between two land cells (a coast of single cells touching at
// their corners is closed).
//
// Range of validity: linear long waves, h ≪ wavelength; the source is a point
// (an earthquake's rupture is 100-1000 km long, and the near field is reached
// earlier where it faces the rupture). Observed first arrivals across an ocean
// come later than the long-wave times by up to a few per cent, since sea water
// is compressible, the Earth elastic and gravity perturbed by the wave (Watada
// et al. 2014); and the times depend on the bathymetry's resolution, most near
// the coasts, where the grid's depth is a mean over its cell.
//
// References:
//   Dijkstra, E. W. (1959), A note on two problems in connexion with graphs,
//   Numer. Math. 1, 269-271.
//   Shokin, Yu. I., Chubarov, L. B., Marchuk, An. G. and Simonov, K. V.
//   (1987), A computational experiment in the tsunami problem, Science of
//   Tsunami Hazards 5, 85-113 (travel-time charts of the Pacific).
//   Wessel, P. (2009), Analysis of observed and predicted tsunami travel
//   times for the Pacific and Indian oceans, Pure Appl. Geophys. 166, 301-324.
//   Watada, S., Kusumoto, S. and Satake, K. (2014), Traveltime delay and
//   initial phase reversal of distant tsunamis coupled with the
//   self-gravitating elastic Earth, J. Geophys. Res. Oceans 119, 4287-4310.
//   Lamb, H. (1932), Hydrodynamics, 6th ed., Cambridge.

import { STANDARD_GRAVITY } from './waves.js';
import { EARTH_RADIUS_M, DEG, centralAngle } from '../../geodesy/js/index.js';

// The speed (m/s) of a long wave over water of depth h (m).
export const longWaveSpeed = (h) => Math.sqrt(STANDARD_GRAVITY * h);

// The stencils: offsets (row, column) in the first quadrant, made symmetric below.
const QUADRANT = {
  16: [[0, 1], [1, 0], [1, 1], [1, 2], [2, 1]],
  32: [[0, 1], [1, 0], [1, 1], [1, 2], [2, 1], [1, 3], [3, 1], [2, 3], [3, 2]],
};

// The cells a straight step from the centre of cell (0, 0) to that of (dr, dc)
// passes through, in order, ends included (cell (r, c) spans r ± ½, c ± ½).
// Sampled finely between the points where the step meets a cell's side or
// corner, so that a step through a corner goes from one cell to its diagonal
// neighbour directly.
function stepCells(dr, dc) {
  const n = 64 * Math.max(Math.abs(dr), Math.abs(dc)), out = [[0, 0]];
  for (let j = 0; j < n; j++) {
    const f = (j + 0.5) / n, r = Math.round(dr * f), c = Math.round(dc * f);
    const [lr, lc] = out[out.length - 1];
    if (r !== lr || c !== lc) out.push([r, c]);
  }
  return out;
}

// A row's offsets, with column offsets scaled by 1/cos φ: [{ dr, dc, cells, corners }],
// `cells` the cells a step passes through between its ends, which must be wet;
// `corners` the pairs of cells beside each corner it passes through, of which
// one must be wet (for the diagonal step (±1, ±1), the two beside it).
const stencils = new Map();
function rowStencil(neighbours, scale) {
  const key = `${neighbours}:${QUADRANT[neighbours].map(([, c]) => (c ? Math.max(1, Math.round(c * scale)) : 0)).join(',')}`;
  if (!stencils.has(key)) stencils.set(key, makeStencil(neighbours, scale));
  return stencils.get(key);
}

function makeStencil(neighbours, scale) {
  const seen = new Set(), out = [];
  const add = (dr, dc) => {
    const key = `${dr},${dc}`;
    if ((!dr && !dc) || seen.has(key)) return;
    seen.add(key);
    const path = stepCells(dr, dc), corners = [];
    for (let k = 1; k < path.length; k++) {
      const [r0, c0] = path[k - 1], [r1, c1] = path[k];
      if (r0 !== r1 && c0 !== c1) corners.push([[r1, c0], [r0, c1]]);
    }
    out.push({ dr, dc, cells: path.slice(1, -1), corners });
  };
  for (const [r, c] of QUADRANT[neighbours]) {
    const sc = c ? Math.max(1, Math.round(c * scale)) : 0;
    for (const [sr, ss] of [[1, 1], [1, -1], [-1, 1], [-1, -1]]) add(sr * r, ss * sc);
  }
  add(0, 1); add(0, -1);
  return out;
}

// A binary min-heap of (key, index) pairs in typed arrays, growing as needed.
class Heap {
  constructor(capacity = 1 << 16) { this.keys = new Float64Array(capacity); this.items = new Int32Array(capacity); this.size = 0; this.top = NaN; }
  push(key, item) {
    if (this.size === this.keys.length) {
      const k = new Float64Array(this.size * 2), it = new Int32Array(this.size * 2);
      k.set(this.keys); it.set(this.items); this.keys = k; this.items = it;
    }
    const keys = this.keys, items = this.items;
    let i = this.size++;
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (keys[p] <= key) break;
      keys[i] = keys[p]; items[i] = items[p]; i = p;
    }
    keys[i] = key; items[i] = item;
  }
  // Removes the least item and returns it; its key is left in this.top.
  pop() {
    const keys = this.keys, items = this.items, item = items[0];
    this.top = keys[0];
    const key = keys[--this.size], last = items[this.size];
    let i = 0;
    for (;;) {
      let child = 2 * i + 1;
      if (child >= this.size) break;
      if (child + 1 < this.size && keys[child + 1] < keys[child]) child++;
      if (keys[child] >= key) break;
      keys[i] = keys[child]; items[i] = items[child]; i = child;
    }
    keys[i] = key; items[i] = last;
    return item;
  }
}

// First-arrival times from a source over a grid of depths.
//   depth: Float32Array (m, positive in the sea; land and missing as NaN or ≤ 0),
//   on a grid { lon0, dlon, nlon, lat0, dlat, nlat, periodic }, rows north to south;
//   source: [lon, lat] in degrees, moved to the nearest wet point within snapKm.
// The source's point and the points its stencil reaches start with their
// straight-line times from the source itself, so that the grid's offset from
// the source adds no time. advance(budget) settles up to `budget` points and
// says whether it is done, so that a page may rest between calls; result()
// gives seconds (NaN where the wave does not reach); solve() runs to the end.
export class TravelTimes {
  constructor(depth, grid, source, { neighbours = 16, minDepth = 1, snapKm = 100 } = {}) {
    if (!QUADRANT[neighbours]) throw new Error(`no stencil of ${neighbours} neighbours`);
    const { nlon, nlat, lat0, dlat, lon0, dlon } = grid, n = nlon * nlat;
    this.grid = grid;
    this.slowness = new Float32Array(n);
    for (let i = 0; i < n; i++) this.slowness[i] = depth[i] > 0 ? 1 / longWaveSpeed(Math.max(depth[i], minDepth)) : NaN;
    // Each row's stencil as typed arrays: the steps' offsets and great-circle lengths (m);
    // the cells each passes through (cstart[k] to cstart[k + 1] in cdr, cdc); and the
    // pairs beside each corner it passes through (kstart[k] to kstart[k + 1] in kdr, kdc, two by two).
    this.rows = [];
    let reach = 0;
    for (let r = 0; r < nlat; r++) {
      const lat = lat0 - r * dlat, scale = Math.min(nlon / 8, 1 / Math.max(Math.cos(lat * DEG), 1e-6));
      const steps = rowStencil(neighbours, scale).filter((s) => r + s.dr >= 0 && r + s.dr < nlat && (grid.periodic || Math.abs(s.dc) < nlon));
      const cells = steps.flatMap((s) => s.cells), corners = steps.flatMap((s) => s.corners.flat());
      const cstart = new Int32Array(steps.length + 1), kstart = new Int32Array(steps.length + 1);
      steps.forEach((s, k) => { cstart[k + 1] = cstart[k] + s.cells.length; kstart[k + 1] = kstart[k] + 2 * s.corners.length; });
      this.rows.push({
        dr: Int32Array.from(steps, (s) => s.dr), dc: Int32Array.from(steps, (s) => s.dc),
        length: Float64Array.from(steps, (s) => EARTH_RADIUS_M * centralAngle(0, lat, s.dc * dlon, lat0 - (r + s.dr) * dlat)),
        cstart, cdr: Int32Array.from(cells, (x) => x[0]), cdc: Int32Array.from(cells, (x) => x[1]),
        kstart, kdr: Int32Array.from(corners, (x) => x[0]), kdc: Int32Array.from(corners, (x) => x[1]),
        width: Math.max(...steps.map((s) => Math.abs(s.dc))),
      });
      reach = Math.max(reach, ...steps.map((s) => Math.abs(s.dr)));
    }
    // Open sea: every cell within reach of a point's stencil is wet, so that
    // none of its steps needs the test for land.
    this.open = new Uint8Array(n);
    const dryBefore = new Int32Array((nlon + 1) * nlat);         // dry cells of each row before each column
    for (let r = 0; r < nlat; r++) {
      for (let c = 0; c < nlon; c++) {
        const x = this.slowness[r * nlon + c];
        dryBefore[r * (nlon + 1) + c + 1] = dryBefore[r * (nlon + 1) + c] + (x === x ? 0 : 1);
      }
    }
    const dryIn = (q, a, b) => {                                  // dry cells of row q in columns a..b, across the wrap
      const at = (k) => dryBefore[q * (nlon + 1) + k];
      if (b - a + 1 >= nlon) return at(nlon);
      const a2 = ((a % nlon) + nlon) % nlon, b2 = ((b % nlon) + nlon) % nlon;
      return a2 <= b2 ? at(b2 + 1) - at(a2) : at(nlon) - at(a2) + at(b2 + 1);
    };
    for (let r = reach; r < nlat - reach; r++) {
      const w = this.rows[r].width;
      for (let c = 0; c < nlon; c++) {
        if (!grid.periodic && (c < w || c >= nlon - w)) continue;
        let dry = 0;
        for (let q = r - reach; q <= r + reach && !dry; q++) dry += dryIn(q, c - w, c + w);
        if (!dry) this.open[r * nlon + c] = 1;
      }
    }
    this.times = new Float64Array(n).fill(Infinity);
    this.done = new Uint8Array(n);
    this.heap = new Heap();
    this.settled = 0;
    // The source: the nearest wet point, searched over a box of rows and columns.
    const [slon, slat] = source;
    let best = -1, bestM = snapKm * 1000;
    const r0 = Math.round((lat0 - slat) / dlat), boxRows = Math.ceil((snapKm * 1000) / (EARTH_RADIUS_M * dlat * DEG)) + 1;
    for (let r = Math.max(0, r0 - boxRows); r <= Math.min(nlat - 1, r0 + boxRows); r++) {
      const lat = lat0 - r * dlat, cos = Math.max(Math.cos(lat * DEG), 1e-3);
      const boxCols = Math.min(nlon >> 1, Math.ceil((snapKm * 1000) / (EARTH_RADIUS_M * dlon * DEG * cos)) + 1);
      const c0 = Math.round((slon - lon0) / dlon);
      for (let c = c0 - boxCols; c <= c0 + boxCols; c++) {
        const cc = grid.periodic ? ((c % nlon) + nlon) % nlon : c;
        if (cc < 0 || cc >= nlon) continue;
        const i = r * nlon + cc;
        if (!(this.slowness[i] === this.slowness[i])) continue;
        const m = EARTH_RADIUS_M * centralAngle(slon, slat, lon0 + cc * dlon, lat);
        if (m < bestM) { bestM = m; best = i; }
      }
    }
    this.source = best < 0 ? null : {
      index: best, lon: lon0 + (best % nlon) * dlon, lat: lat0 - Math.floor(best / nlon) * dlat, depth: depth[best], movedKm: bestM / 1000,
    };
    if (best >= 0) {
      // The source's point and those its stencil reaches, at their straight-line times from the source.
      const start = (i) => {
        const lon = lon0 + (i % nlon) * dlon, lat = lat0 - Math.floor(i / nlon) * dlat;
        const t = 0.5 * EARTH_RADIUS_M * centralAngle(slon, slat, lon, lat) * (this.slowness[best] + this.slowness[i]);
        if (t < this.times[i]) { this.times[i] = t; this.heap.push(t, i); }
      };
      start(best);
      this.relax(best, start);
    }
  }

  // Calls f(j, k) for each wet point j reached from point i by a step k of its
  // row's stencil that does not pass over land.
  relax(i, f) {
    const { nlon, periodic } = this.grid, slowness = this.slowness;
    const r = (i / nlon) | 0, c = i - r * nlon;
    const { dr, dc, cstart, cdr, cdc, kstart, kdr, kdc } = this.rows[r], open = this.open[i];
    const wet = (rr, cc) => {
      if (cc < 0) cc += nlon; else if (cc >= nlon) cc -= nlon;
      const s = slowness[rr * nlon + cc];
      return s === s;
    };
    for (let k = 0; k < dr.length; k++) {
      let cc = c + dc[k];
      if (cc < 0 || cc >= nlon) { if (!periodic) continue; cc = cc < 0 ? cc + nlon : cc - nlon; }
      const j = (r + dr[k]) * nlon + cc;
      if (!open) {
        const sj = slowness[j];
        if (!(sj === sj)) continue;
        // Every cell the step passes through must be wet, and one of the two beside each corner.
        let pass = true;
        for (let q = cstart[k]; q < cstart[k + 1] && pass; q++) pass = wet(r + cdr[q], c + cdc[q]);
        for (let q = kstart[k]; q < kstart[k + 1] && pass; q += 2) pass = wet(r + kdr[q], c + kdc[q]) || wet(r + kdr[q + 1], c + kdc[q + 1]);
        if (!pass) continue;
      }
      f(j, k);
    }
  }

  advance(budget = Infinity) {
    const { nlon } = this.grid, { times, done, slowness, heap } = this;
    let count = 0, t = 0, si = 0, length = null;
    const step = (j, k) => {
      if (done[j]) return;
      const tj = t + 0.5 * length[k] * (si + slowness[j]);
      if (tj < times[j]) { times[j] = tj; heap.push(tj, j); }
    };
    while (heap.size && count < budget) {
      const i = heap.pop();
      if (done[i]) continue;
      done[i] = 1;
      this.settled++;
      count++;
      t = heap.top; si = slowness[i]; length = this.rows[(i / nlon) | 0].length;
      this.relax(i, step);
    }
    return heap.size === 0;
  }

  solve() { this.advance(Infinity); return this.result(); }

  // The times in seconds (Float32Array, NaN where not reached).
  result() {
    const out = new Float32Array(this.times.length);
    for (let i = 0; i < out.length; i++) out[i] = this.times[i] === Infinity ? NaN : this.times[i];
    return out;
  }
}

// First-arrival times (s) from a source, run to the end: { times, source }.
export function travelTimes(depth, grid, source, options = {}) {
  const solver = new TravelTimes(depth, grid, source, options);
  return { times: solver.solve(), source: solver.source };
}
