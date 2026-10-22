// Zonal statistics on the sphere: a field summarised over the regions of a
// regions raster (countries, provinces), as the areal mean, the sum of a count,
// the least and greatest value, and the area and number of cells.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (GIS).
//
// A grid is { lon0, dlon, nlon, lat0, dlat, nlat, periodic } (physis's grids):
// the centre of column c is at lon0 + c dlon, of row r at lat0 - r dlat, and the
// cell is the rectangle of longitude and latitude of those sizes about it, its
// edges clamped to the poles. The field and the regions need not share a grid:
// both are taken as constant on their cells, and the statistics are the exact
// integrals of the one over the other. The overlap of two cells is a rectangle
// of longitude and latitude, whose area on the sphere of radius R is
//
//   |A ∩ B| = R² Δλ Δ(sin φ),
//
// the product of an overlap in longitude and one in the sine of latitude, so the
// overlaps of every cell of one grid with every cell of the other are the
// products of two short lists computed once (overlapTable). This serves equally
// a field finer than the regions (ETOPO at 5′ on regions at 0.1°) and coarser
// (a forecast at 1°, whose cells each meet a hundred cells of regions), and it
// keeps sums: a count is divided among the regions a cell meets in proportion
// to the areas it shares with them. See theory/zonal-statistics.md.

import { EARTH_RADIUS_M } from '../../geodesy/js/index.js';

const RAD = Math.PI / 180;

// The areas of the cells of each row of a grid (m², or in units of R² for radius 1).
export function cellAreas(grid, radius = EARTH_RADIUS_M) {
  const out = new Float64Array(grid.nlat), dl = grid.dlon * RAD;
  for (let r = 0; r < grid.nlat; r++) {
    const [s, n] = rowEdges(grid, r);
    out[r] = radius * radius * dl * (Math.sin(n * RAD) - Math.sin(s * RAD));
  }
  return out;
}

// The southern and northern edges of row r, in degrees, clamped to the poles.
export function rowEdges(grid, r) {
  const c = grid.lat0 - r * grid.dlat;
  return [Math.max(-90, c - grid.dlat / 2), Math.min(90, c + grid.dlat / 2)];
}

// The overlaps of the cells of `grid` with those of `regionGrid`: for each row
// of the grid, the rows of the regions it meets and the overlap in sin φ; for
// each column, the columns of the regions it meets and the overlap in radians of
// longitude (wrapped when the regions go round the Earth). Lists are flat:
// rows.at[rows.start[r] .. rows.start[r + 1]) and rows.weight likewise.
export function overlapTable(grid, regionGrid) {
  const g = regionGrid;
  const north = g.lat0 + g.dlat / 2, west = g.lon0 - g.dlon / 2;
  const rows = axis(grid.nlat, (r) => {
    const [s, n] = rowEdges(grid, r), out = [];
    const k0 = Math.max(0, Math.floor((north - n) / g.dlat)), k1 = Math.min(g.nlat - 1, Math.ceil((north - s) / g.dlat) - 1);
    for (let k = k0; k <= k1; k++) {
      const [rs, rn] = rowEdges(g, k);
      const w = Math.sin(Math.min(n, rn) * RAD) - Math.sin(Math.max(s, rs) * RAD);
      if (w > 0) out.push([k, w]);
    }
    return out;
  });
  const cols = axis(grid.nlon, (c) => {
    const a = (grid.lon0 + (c - 0.5) * grid.dlon - west) / g.dlon, b = a + grid.dlon / g.dlon, out = [];
    for (let k = Math.floor(a); k < Math.ceil(b); k++) {
      const w = (Math.min(b, k + 1) - Math.max(a, k)) * g.dlon * RAD;
      if (w <= 0) continue;
      let kk = k;
      if (g.periodic) kk = ((k % g.nlon) + g.nlon) % g.nlon;
      else if (k < 0 || k >= g.nlon) continue;
      out.push([kk, w]);
    }
    return out;
  });
  return { rows, cols, regionGrid: g, grid };
}

function axis(n, overlaps) {
  const start = new Int32Array(n + 1), at = [], weight = [];
  for (let i = 0; i < n; i++) {
    for (const [k, w] of overlaps(i)) { at.push(k); weight.push(w); }
    start[i + 1] = at.length;
  }
  return { start, at: Int32Array.from(at), weight: Float64Array.from(weight) };
}

// Statistics of `values` (row-major on `grid`, NaN where there is no value) over
// the regions of `regions` (row-major ids on `regionGrid`, 0 for none). For each
// id k from 0 to `count` (0 collects what lies in no region):
//   area[k]   the area of region k where the field has a value (m²);
//   mean[k]   the area-weighted mean of the field over it;
//   min[k], max[k]   the least and greatest values of the cells that meet it;
//   cells[k]  the number of cells with a value that meet it;
//   sum[k]    for a count per cell (people, cells of a kind): the field's cells
//             divided among the regions in proportion to their shared areas.
// With share = 'regions' (the default) a cell's count is divided among the
// regions it meets only, as for people, who live on land when the cell is
// partly sea; the count of a cell meeting no region goes to sum[0], unless
// `reach` > 0: it then goes to the region of the cell of the regions' grid,
// within `reach` such cells of it, nearest to its centre (reached counts it), so
// that a coast or an island finer than the regions' cells keeps its people.
// With share = 'area' the part of the cell outside every region goes to sum[0].
// total is the sum of all values. Pass `table` from overlapTable to reuse it.
export function zonalStatistics(values, grid, regions, regionGrid, { count, share = 'regions', reach = 0, radius = EARTH_RADIUS_M, table } = {}) {
  const t = table || overlapTable(grid, regionGrid);
  const n = (count ?? maxId(regions)) + 1, rn = regionGrid.nlon;
  const area = new Float64Array(n), weighted = new Float64Array(n), sum = new Float64Array(n);
  const min = new Float64Array(n).fill(Infinity), max = new Float64Array(n).fill(-Infinity), cells = new Int32Array(n);
  const ids = new Int32Array(1024), ws = new Float64Array(1024);
  const byArea = share === 'area';
  let total = 0, reached = 0;
  for (let r = 0; r < grid.nlat; r++) {
    const r0 = t.rows.start[r], r1 = t.rows.start[r + 1];
    for (let c = 0; c < grid.nlon; c++) {
      const v = values[r * grid.nlon + c];
      if (v !== v) continue;
      total += v;
      const c0 = t.cols.start[c], c1 = t.cols.start[c + 1];
      // The ids this cell meets and the areas it shares with each (in units of R²).
      let m = 0, cellArea = 0, outside = 0;
      for (let a = r0; a < r1; a++) {
        const base = t.rows.at[a] * rn, wy = t.rows.weight[a];
        for (let b = c0; b < c1; b++) {
          const id = regions[base + t.cols.at[b]], w = wy * t.cols.weight[b];
          cellArea += w;
          if (id === 0) outside += w;
          let j = 0;
          while (j < m && ids[j] !== id) j++;
          if (j === m) { ids[m] = id; ws[m] = 0; m++; }
          ws[j] += w;
        }
      }
      if (!m) { sum[0] += v; continue; }                  // beyond the regions' grid
      const inRegions = cellArea - outside;
      for (let j = 0; j < m; j++) {
        const id = ids[j], w = ws[j];
        if (id >= n) throw new Error(`region id ${id} beyond the count ${n - 1}`);
        area[id] += w; weighted[id] += w * v; cells[id]++;
        if (v < min[id]) min[id] = v;
        if (v > max[id]) max[id] = v;
        if (byArea) sum[id] += (v * w) / cellArea;
        else if (id !== 0) sum[id] += (v * w) / inRegions;
      }
      if (!byArea && inRegions <= 0) {
        const id = reach > 0 && v !== 0
          ? nearestRegion(regions, regionGrid, t, r0, r1, c0, c1, reach, grid.lon0 + c * grid.dlon, grid.lat0 - r * grid.dlat) : 0;
        if (id >= n) throw new Error(`region id ${id} beyond the count ${n - 1}`);
        sum[id] += v;
        if (id) reached += v;
      }
    }
  }
  const mean = new Float64Array(n), R2 = radius * radius;
  for (let k = 0; k < n; k++) {
    mean[k] = area[k] > 0 ? weighted[k] / area[k] : NaN;
    if (!cells[k]) { min[k] = NaN; max[k] = NaN; }
    area[k] *= R2;
  }
  return { area, mean, sum, min, max, cells, total, unassigned: sum[0], reached };
}

// The region of the cell of `g` nearest to (lon, lat) among those within `reach`
// cells of the rows a0..a1 and columns b0..b1 of the table, or 0. Distances are
// those of the plane tangent at the place, enough to rank cells so near.
function nearestRegion(regions, g, t, a0, a1, b0, b1, reach, lon, lat) {
  const k0 = Math.max(0, t.rows.at[a0] - reach), k1 = Math.min(g.nlat - 1, t.rows.at[a1 - 1] + reach);
  const first = t.cols.at[b0], last = t.cols.at[b1 - 1];
  const width = (g.periodic ? (last - first + g.nlon) % g.nlon : last - first) + 2 * reach;
  const cos = Math.cos(lat * RAD);
  let best = 0, nearest = Infinity;
  for (let k = k0; k <= k1; k++) {
    const dlat = g.lat0 - k * g.dlat - lat;
    for (let s = 0; s <= width; s++) {
      let kc = first - reach + s;
      if (g.periodic) kc = ((kc % g.nlon) + g.nlon) % g.nlon;
      else if (kc < 0 || kc >= g.nlon) continue;
      const id = regions[k * g.nlon + kc];
      if (!id) continue;
      let dlon = g.lon0 + kc * g.dlon - lon;
      dlon -= 360 * Math.round(dlon / 360);
      const d = dlat * dlat + cos * cos * dlon * dlon;
      if (d < nearest) { nearest = d; best = id; }
    }
  }
  return best;
}

function maxId(regions) {
  let m = 0;
  for (let i = 0; i < regions.length; i++) if (regions[i] > m) m = regions[i];
  return m;
}

// The id of the region at a place (the cell containing it), or 0.
export function regionAt(regions, regionGrid, lon, lat) {
  const g = regionGrid;
  let c = Math.floor((lon - (g.lon0 - g.dlon / 2)) / g.dlon);
  if (g.periodic) c = ((c % g.nlon) + g.nlon) % g.nlon;
  else if (c < 0 || c >= g.nlon) return 0;
  const r = Math.floor((g.lat0 + g.dlat / 2 - lat) / g.dlat);
  if (r < 0 || r >= g.nlat) return 0;
  return regions[r * g.nlon + c];
}

// A statistic by region code: { code: value } for the regions listed
// ([{ id, code }], as in a regions product) that hold any of the field.
export function byCode(stats, list, key = 'sum') {
  const out = {};
  for (const region of list) {
    if (region.id < stats.cells.length && stats.cells[region.id] > 0) out[region.code] = stats[key][region.id];
  }
  return out;
}
