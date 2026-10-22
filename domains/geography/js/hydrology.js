// Terrain and water on the sphere: the slope of the surface; depressions filled
// by Priority-Flood; the direction of flow from each cell to the steepest of its
// eight neighbours (D8), with distances on the sphere; the area draining through
// each cell; the watershed of a place, with its outline; and the path of water
// downstream from a place, or up its main stem.
//
// Provisional JavaScript, standing in for the C++ planned for this domain
// (physical geography: rivers and drainage basins).
//
// Elevations are in metres on a grid as in zonal.js (rows north to south). The
// outlets, where water leaves the land, are the seas: the cells at or below sea
// level (elevation <= 0 in ETOPO) that form, with their neighbours at or below
// it, a body of water of at least 100 000 km² (the ocean, and on ETOPO's 10′ grid
// the Caspian and the Black Sea), or one that reaches the edge of a grid not
// round the Earth; with the cells without a value, and that grid's edges. The
// floors of deep lakes (Superior, Baikal) and dry depressions (Qattara, the
// Faiyum, the Dead Sea), below sea level but no sea, are filled instead, so that
// rivers pass through them to the sea. Every other cell drains, through a chain
// of cells each lower than the last, to one outlet.
//
// Sources:
//   Priority-Flood+ε: Barnes, R., Lehman, C., and Mulla, D. (2014), "Priority-flood:
//   an optimal depression-filling and watershed-labeling algorithm for digital
//   elevation models", Computers & Geosciences 62, 117-127, Algorithm 3.
//   D8: O'Callaghan, J. F., and Mark, D. M. (1984), "The extraction of drainage
//   networks from digital elevation data", Computer Vision, Graphics, and Image
//   Processing 28, 323-344.
//   Slope: Zevenbergen, L. W., and Thorne, C. R. (1987), "Quantitative analysis of
//   land surface topography", Earth Surface Processes and Landforms 12, 47-56: the
//   slope (G² + H²)^½ of their quadratic surface, G and H centred differences,
//   here with the metric of the sphere.
// See theory/drainage-from-elevation.md.

import { EARTH_RADIUS_M, centralAngle } from '../../geodesy/js/index.js';
import { cellAreas } from './zonal.js';

const RAD = Math.PI / 180;

// The eight neighbours, counter-clockwise from the east: [row step, column step].
// Rows run north to south, so north is a step of -1.
export const NEIGHBOURS = [[0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [1, 0], [1, 1]];

const DR = Int8Array.from(NEIGHBOURS, (n) => n[0]), DC = Int8Array.from(NEIGHBOURS, (n) => n[1]);

// The neighbours on one grid: (i, k) -> the index of the neighbour of cell i in
// direction k, or -1 beyond the grid. The grid's numbers are held in the
// closure, so that the loops over millions of cells below stay fast whatever
// the shape of the grid object passed.
export function stepper({ nlon, nlat, periodic }) {
  return (i, k) => {
    const r0 = (i / nlon) | 0, r = r0 + DR[k];
    if (r < 0 || r >= nlat) return -1;
    let c = i - r0 * nlon + DC[k];
    if (c < 0) { if (!periodic) return -1; c += nlon; } else if (c >= nlon) { if (!periodic) return -1; c -= nlon; }
    return r * nlon + c;
  };
}

// The index of the neighbour of cell i in direction k, or -1 beyond the grid.
export const neighbour = (grid, i, k) => stepper(grid)(i, k);

// The index of the cell containing a place, or -1 outside the grid.
export function cellOf(grid, lon, lat) {
  let c = Math.round((lon - grid.lon0) / grid.dlon);
  if (grid.periodic) c = ((c % grid.nlon) + grid.nlon) % grid.nlon;
  else if (c < 0 || c >= grid.nlon) return -1;
  const r = Math.round((grid.lat0 - lat) / grid.dlat);
  return r < 0 || r >= grid.nlat ? -1 : r * grid.nlon + c;
}

// The centre of cell i: [lon, lat].
export const centreOf = (grid, i) => [grid.lon0 + (i % grid.nlon) * grid.dlon, grid.lat0 - Math.floor(i / grid.nlon) * grid.dlat];

// Great-circle distances from a cell of each row to its eight neighbours:
// Float64Array(nlat × 8), metres.
export function neighbourDistances(grid, radius = EARTH_RADIUS_M) {
  const out = new Float64Array(grid.nlat * 8);
  for (let r = 0; r < grid.nlat; r++) {
    const lat = grid.lat0 - r * grid.dlat;
    NEIGHBOURS.forEach(([dr, dc], k) => { out[r * 8 + k] = radius * centralAngle(0, lat, dc * grid.dlon, lat - dr * grid.dlat); });
  }
  return out;
}

// ---- slope -------------------------------------------------------------------------

// The slope of the surface, |∇z| (rise over run, the tangent of the slope's
// angle), by centred differences: ∂z/∂x = (z_E − z_W) / (2 R cos φ Δλ),
// ∂z/∂y = (z_N − z_S) / (2 R Δφ); one-sided at the first and last rows and at the
// edges of a grid that does not go round the Earth. NaN where a value needed is
// missing, and at a pole.
// Accuracy: the slope of the surface smoothed over two cells (37 km on a 10′
// grid), which in mountains is several times smaller than the slopes of the
// ground itself; on a grid sampled from a finer one, aliased.
export function slope(z, grid, radius = EARTH_RADIUS_M) {
  const { nlon, nlat, dlon, dlat, lat0, periodic } = grid, out = new Float32Array(nlon * nlat).fill(NaN);
  const dy = radius * dlat * RAD;
  for (let r = 0; r < nlat; r++) {
    const cos = Math.cos((lat0 - r * dlat) * RAD);
    if (cos < 1e-9) continue;
    const dx = radius * cos * dlon * RAD;
    const rn = r > 0 ? r - 1 : r, rs = r < nlat - 1 ? r + 1 : r;
    for (let c = 0; c < nlon; c++) {
      let ce = c + 1, cw = c - 1;
      if (ce >= nlon) ce = periodic ? 0 : c;
      if (cw < 0) cw = periodic ? nlon - 1 : c;
      const gx = (z[r * nlon + ce] - z[r * nlon + cw]) / ((ce === c || cw === c ? 1 : 2) * dx);
      const gy = (z[rn * nlon + c] - z[rs * nlon + c]) / ((rs - rn) * dy);
      out[r * nlon + c] = Math.hypot(gx, gy);
    }
  }
  return out;
}

// ---- Priority-Flood ---------------------------------------------------------------

// The next float64 above x.
const f64 = new Float64Array(1), u64 = new BigUint64Array(f64.buffer);
export function nextUp(x) {
  if (!(x < Infinity)) return x;
  if (x === 0) return Number.MIN_VALUE;
  f64[0] = x;
  u64[0] += x > 0 ? 1n : -1n;
  return f64[0];
}

// A binary heap of cells by priority, least first.
class Heap {
  constructor(capacity = 1024) { this.key = new Float64Array(capacity); this.cell = new Int32Array(capacity); this.size = 0; }
  push(key, cell) {
    if (this.size === this.key.length) {
      const k = new Float64Array(this.size * 2), c = new Int32Array(this.size * 2);
      k.set(this.key); c.set(this.cell); this.key = k; this.cell = c;
    }
    let i = this.size++;
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (this.key[p] <= key) break;
      this.key[i] = this.key[p]; this.cell[i] = this.cell[p]; i = p;
    }
    this.key[i] = key; this.cell[i] = cell;
  }
  pop() {
    const top = this.cell[0], key = this.key[--this.size], cell = this.cell[this.size];
    let i = 0;
    for (;;) {
      let m = 2 * i + 1;
      if (m >= this.size) break;
      if (m + 1 < this.size && this.key[m + 1] < this.key[m]) m++;
      if (this.key[m] >= key) break;
      this.key[i] = this.key[m]; this.cell[i] = this.cell[m]; i = m;
    }
    this.key[i] = key; this.cell[i] = cell;
    return top;
  }
}

// Which cells are outlets: those without a value; those at or below `seaLevel`
// in a body of such cells (neighbours of eight) of at least `minSeaArea` km², or
// reaching the edge of a grid not round the Earth; and, when `edges`, the edges
// of the grid. Uint8Array, 1 for an outlet.
export function outlets(z, grid, { seaLevel = 0, minSeaArea = 1e5, edges = !grid.periodic, radius = EARTH_RADIUS_M } = {}) {
  const { nlon, nlat } = grid, n = nlon * nlat, out = new Uint8Array(n), seen = new Uint8Array(n), next = stepper(grid);
  const areas = cellAreas(grid, radius), stack = new Int32Array(n);
  for (let s = 0; s < n; s++) {
    if (z[s] !== z[s]) { out[s] = 1; continue; }
    if (seen[s] || z[s] > seaLevel) continue;
    // The body of water at or below sea level that holds cell s, by depth-first search.
    let top = 0, size = 0, area = 0, open = false;
    stack[top++] = s; seen[s] = 1;
    while (top) {
      const c = stack[--top], r = Math.floor(c / nlon), q = c % nlon;
      stack[n - 1 - size++] = c;                          // the body's cells, from the far end
      area += areas[r] / 1e6;
      if (!grid.periodic && (r === 0 || r === nlat - 1 || q === 0 || q === nlon - 1)) open = true;
      for (let k = 0; k < 8; k++) {
        const j = next(c, k);
        if (j >= 0 && !seen[j] && z[j] <= seaLevel) { seen[j] = 1; stack[top++] = j; }
      }
    }
    if (area >= minSeaArea || open) for (let k = 0; k < size; k++) out[stack[n - 1 - k]] = 1;
  }
  if (edges) {
    for (let c = 0; c < nlon; c++) { out[c] = 1; out[(nlat - 1) * nlon + c] = 1; }
    for (let r = 0; r < nlat; r++) { out[r * nlon] = 1; out[r * nlon + nlon - 1] = 1; }
  }
  return out;
}

// Depressions filled by Priority-Flood+ε (Barnes et al. 2014, Algorithm 3): the
// cells are flooded from the outlets inward in order of height, from a priority
// queue; a cell reached from a higher one lies in a depression, and is raised
// to the next float above that one and taken next, from a plain queue (the
// pit). Every cell that is not an outlet is then higher than the neighbour it
// was reached from, so that water leaves it, and a depression becomes a surface
// rising by the least possible steps from its spill point. When the next cell
// of the priority queue is no higher than the next of the pit, it goes first,
// so that cells are taken in order of their filled height.
// Returns { filled: Float64Array (outlets keep their height; a missing value is
// -Infinity), outlet: Uint8Array, raised: number of cells raised }.
export function fillDepressions(z, grid, options = {}) {
  const n = grid.nlon * grid.nlat, outlet = options.outlet || outlets(z, grid, options);
  const filled = new Float64Array(n), closed = new Uint8Array(n), open = new Heap(1 << 16), next = stepper(grid);
  const pit = new Int32Array(n);
  let head = 0, tail = 0, raised = 0;
  for (let i = 0; i < n; i++) filled[i] = z[i] === z[i] ? z[i] : -Infinity;
  // The outlets that touch land seed the flood; the rest are closed and never taken.
  for (let i = 0; i < n; i++) {
    if (!outlet[i]) continue;
    closed[i] = 1;
    for (let k = 0; k < 8; k++) {
      const j = next(i, k);
      if (j >= 0 && !outlet[j]) { open.push(filled[i], i); break; }
    }
  }
  while (open.size || head < tail) {
    let c;
    if (head < tail && !(open.size && open.key[0] <= filled[pit[head]])) c = pit[head++];
    else c = open.pop();
    const spill = nextUp(filled[c]);
    for (let k = 0; k < 8; k++) {
      const j = next(c, k);
      if (j < 0 || closed[j]) continue;
      closed[j] = 1;
      if (filled[j] <= spill) {
        if (filled[j] < filled[c]) raised++;
        filled[j] = spill;
        pit[tail++] = j;
      } else {
        open.push(filled[j], j);
      }
    }
  }
  return { filled, outlet, raised };
}

// ---- flow ---------------------------------------------------------------------------

// D8 flow directions (O'Callaghan and Mark 1984): from each cell that is not an
// outlet, to the neighbour of steepest descent, (h_c − h_n) / d_cn, with d the
// great-circle distance between the cells' centres. Int8Array: the index into
// NEIGHBOURS, or -1 for an outlet (and for a cell with no lower neighbour, which
// on a surface filled by fillDepressions does not occur).
export function flowDirections(filled, grid, outlet, radius = EARTH_RADIUS_M) {
  const n = grid.nlon * grid.nlat, dirs = new Int8Array(n).fill(-1), d = neighbourDistances(grid, radius), next = stepper(grid);
  for (let i = 0; i < n; i++) {
    if (outlet[i]) continue;
    const r = (i / grid.nlon) | 0, h = filled[i];
    let best = 0, dir = -1;
    for (let k = 0; k < 8; k++) {
      const j = next(i, k);
      if (j < 0) continue;
      const drop = (h - filled[j]) / d[r * 8 + k];
      if (drop > best) { best = drop; dir = k; }
    }
    dirs[i] = dir;
  }
  return dirs;
}

// The cell that cell i drains to, or -1.
export const receiver = (grid, dirs, i) => (dirs[i] < 0 ? -1 : neighbour(grid, i, dirs[i]));

// Every cell's receiver at once: Int32Array, -1 for an outlet.
export function receivers(dirs, grid) {
  const next = stepper(grid), down = new Int32Array(dirs.length);
  for (let i = 0; i < dirs.length; i++) down[i] = dirs[i] < 0 ? -1 : next(i, dirs[i]);
  return down;
}

// The area draining through each cell, its own included (km²): each cell's area on
// the sphere carried down its flow path. Cells are taken in topological order
// (Kahn 1962): a cell is passed on once every cell draining into it has been.
export function flowAccumulation(dirs, grid, radius = EARTH_RADIUS_M) {
  const n = grid.nlon * grid.nlat, areas = cellAreas(grid, radius), acc = new Float64Array(n);
  const inflow = new Int32Array(n), down = receivers(dirs, grid), queue = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    acc[i] = areas[(i / grid.nlon) | 0] / 1e6;
    if (down[i] >= 0) inflow[down[i]]++;
  }
  let head = 0, tail = 0;
  for (let i = 0; i < n; i++) if (!inflow[i]) queue[tail++] = i;
  while (head < tail) {
    const i = queue[head++], j = down[i];
    if (j < 0) continue;
    acc[j] += acc[i];
    if (--inflow[j] === 0) queue[tail++] = j;
  }
  return acc;
}

// ---- watersheds and paths -------------------------------------------------------------

// The watershed of cell `outletCell`: the cell and every cell draining to it.
// Int32Array of cell indices, the outlet first, by distance upstream.
export function watershed(dirs, grid, outletCell) {
  const seen = new Uint8Array(grid.nlon * grid.nlat), list = [outletCell], next = stepper(grid);
  seen[outletCell] = 1;
  for (let h = 0; h < list.length; h++) {
    const c = list[h];
    for (let k = 0; k < 8; k++) {
      const j = next(c, k);
      if (j >= 0 && !seen[j] && dirs[j] >= 0 && next(j, dirs[j]) === c) { seen[j] = 1; list.push(j); }
    }
  }
  return Int32Array.from(list);
}

// The cells from `cell` down to its outlet, both included.
export function downstream(dirs, grid, cell) {
  const path = [cell], limit = grid.nlon * grid.nlat, next = stepper(grid);
  for (let i = dirs[cell] < 0 ? -1 : next(cell, dirs[cell]); i >= 0 && path.length <= limit; i = dirs[i] < 0 ? -1 : next(i, dirs[i])) path.push(i);
  return Int32Array.from(path);
}

// The main stem above `cell`: from it, up to the neighbour draining into it
// with the greatest upstream area, until a cell into which none drains.
export function mainStem(dirs, accumulation, grid, cell) {
  const path = [cell], next = stepper(grid);
  for (let c = cell; ;) {
    let best = -1, area = 0;
    for (let k = 0; k < 8; k++) {
      const j = next(c, k);
      if (j >= 0 && dirs[j] >= 0 && next(j, dirs[j]) === c && accumulation[j] > area) { area = accumulation[j]; best = j; }
    }
    if (best < 0) break;
    path.push(best);
    c = best;
  }
  return Int32Array.from(path);
}

// The cell within `reach` cells of `cell` with the greatest upstream area, so
// that a place chosen near a river is taken on it (a pour point snapped to the
// stream, as in Lindsay, Rothwell and Davies 2008, Hydrol. Process. 22, 1701).
export function snapToStream(accumulation, grid, cell, reach = 1) {
  const r0 = Math.floor(cell / grid.nlon), c0 = cell % grid.nlon;
  let best = cell;
  for (let dr = -reach; dr <= reach; dr++) {
    const r = r0 + dr;
    if (r < 0 || r >= grid.nlat) continue;
    for (let dc = -reach; dc <= reach; dc++) {
      let c = c0 + dc;
      if (c < 0 || c >= grid.nlon) { if (!grid.periodic) continue; c = (c + grid.nlon) % grid.nlon; }
      const j = r * grid.nlon + c;
      if (accumulation[j] > accumulation[best]) best = j;
    }
  }
  return best;
}

// The centres of a list of cells as [lon, lat, ...], and the path's length on the sphere (m).
export function pathOf(grid, cells, radius = EARTH_RADIUS_M) {
  const coords = new Float64Array(cells.length * 2);
  let length = 0;
  cells.forEach((i, k) => {
    const [lon, lat] = centreOf(grid, i);
    coords[2 * k] = lon > 180 ? lon - 360 : lon; coords[2 * k + 1] = lat;
    if (k) length += radius * centralAngle(coords[2 * k - 2], coords[2 * k - 1], lon, lat);
  });
  return { coords, length };
}

// The outline of a set of cells: closed rings [lon, lat, ...] along the edges
// of the cells that separate the set from the rest, the set on the left of each
// ring (counter-clockwise round an outer edge), with the corners where the
// outline runs straight on left out.
export function cellOutline(grid, cells) {
  const { nlon, nlat } = grid, inside = new Uint8Array(nlon * nlat);
  for (const i of cells) inside[i] = 1;
  const at = (r, c) => {
    if (r < 0 || r >= nlat) return 0;
    if (c < 0 || c >= nlon) { if (!grid.periodic) return 0; c = (c + nlon) % nlon; }
    return inside[r * nlon + c];
  };
  // Corners (i, j): column edge i (0 to nlon, the same as 0 when periodic), row edge j (0 to nlat).
  const key = (i, j) => (grid.periodic ? i % nlon : i) * (nlat + 1) + j;
  const edges = new Map();
  const add = (i0, j0, i1, j1) => {
    const k = key(i0, j0);
    if (!edges.has(k)) edges.set(k, []);
    edges.get(k).push([i0, j0, i1, j1]);
  };
  for (const cell of cells) {
    const r = Math.floor(cell / nlon), c = cell % nlon;
    // The cell's corners counter-clockwise (east to the right, north up): SW, SE, NE, NW.
    if (!at(r + 1, c)) add(c, r + 1, c + 1, r + 1);       // south side, westward to eastward
    if (!at(r, c + 1)) add(c + 1, r + 1, c + 1, r);       // east side, northward
    if (!at(r - 1, c)) add(c + 1, r, c, r);               // north side, westward
    if (!at(r, c - 1)) add(c, r, c, r + 1);               // west side, southward
  }
  const west = grid.lon0 - grid.dlon / 2, north = grid.lat0 + grid.dlat / 2, rings = [];
  const point = (i, j) => [((west + i * grid.dlon + 540) % 360) - 180, Math.max(-90, Math.min(90, north - j * grid.dlat))];
  for (const list of edges.values()) {
    while (list.length) {
      const first = list.pop(), ring = [first];
      for (let e = first; ;) {
        const next = edges.get(key(e[2], e[3]));
        if (!next || !next.length) break;
        e = next.pop();
        ring.push(e);
        if (key(e[2], e[3]) === key(first[0], first[1])) break;
      }
      const coords = [];
      ring.forEach((e, k) => {
        const p = ring[(k + ring.length - 1) % ring.length];
        const straight = (e[2] - e[0]) === (p[2] - p[0]) && (e[3] - e[1]) === (p[3] - p[1]);
        if (!straight) coords.push(...point(e[0], e[1]));
      });
      if (coords.length) coords.push(coords[0], coords[1]);
      rings.push(coords);
    }
  }
  return rings;
}
