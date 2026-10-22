// Isolines of a gridded field by marching squares, and its local extrema
// (the highs and lows of a pressure chart).
//
// Provisional stand-in for core/domains/geography (cartography; laboratory ADR 0014).
//
// Each grid cell with corners on both sides of a level is crossed by one or two
// line segments, whose ends are found by linear interpolation along the cell's
// edges. The two ambiguous (saddle) cases are resolved by the mean of the four
// corners. Cells with a missing corner are skipped. On a periodic grid the cells
// between the last and the first column are included, so lines close round the
// Earth. Segments are returned in longitude and latitude (degrees); a segment
// from the last column to the first ends at longitude lon0 + 360, which the
// renderers treat as the same meridian.
//
// Sources: Lorensen and Cline (1987), "Marching cubes", Computer Graphics 21(4),
// 163-169, specialised to two dimensions; Maple (2003), "Geometric design and
// space planning using the marching squares and marching cube algorithms".

// For each of the 16 corner states (bit 0 top-left, 1 top-right, 2 bottom-right,
// 3 bottom-left set when the corner is at or above the level), the pairs of cell
// edges joined (edges: 0 top, 1 right, 2 bottom, 3 left).
const CASES = [
  [], [[3, 0]], [[0, 1]], [[3, 1]], [[1, 2]], null, [[0, 2]], [[3, 2]],
  [[2, 3]], [[0, 2]], null, [[1, 2]], [[1, 3]], [[0, 1]], [[3, 0]], [],
];

export function isolines(values, grid, level) {
  const { nlon, nlat, lon0, dlon, lat0, dlat } = grid;
  const columns = grid.periodic ? nlon : nlon - 1;
  const out = [];
  const point = (r, c, edge, v) => {
    // Corner values in order tl, tr, br, bl; the edge's two corners.
    const [a, b, ra, ca, rb, cb] = edge === 0 ? [v[0], v[1], r, c, r, c + 1]
      : edge === 1 ? [v[1], v[2], r, c + 1, r + 1, c + 1]
        : edge === 2 ? [v[3], v[2], r + 1, c, r + 1, c + 1]
          : [v[0], v[3], r, c, r + 1, c];
    const t = a === b ? 0.5 : (level - a) / (b - a);
    return [lon0 + (ca + (cb - ca) * t) * dlon, lat0 - (ra + (rb - ra) * t) * dlat];
  };
  const v = [0, 0, 0, 0];
  for (let r = 0; r < nlat - 1; r++) {
    for (let c = 0; c < columns; c++) {
      const c1 = c + 1 === nlon ? 0 : c + 1;
      v[0] = values[r * nlon + c]; v[1] = values[r * nlon + c1];
      v[2] = values[(r + 1) * nlon + c1]; v[3] = values[(r + 1) * nlon + c];
      if (!(v[0] === v[0] && v[1] === v[1] && v[2] === v[2] && v[3] === v[3])) continue;
      const state = (v[0] >= level ? 1 : 0) | (v[1] >= level ? 2 : 0) | (v[2] >= level ? 4 : 0) | (v[3] >= level ? 8 : 0);
      let pairs = CASES[state];
      if (pairs === null) {
        const centreHigh = (v[0] + v[1] + v[2] + v[3]) / 4 >= level;
        pairs = state === 5
          ? (centreHigh ? [[3, 2], [0, 1]] : [[3, 0], [1, 2]])
          : (centreHigh ? [[3, 0], [1, 2]] : [[0, 1], [2, 3]]);
      }
      for (const [e1, e2] of pairs) {
        const p = point(r, c, e1, v), q = point(r, c, e2, v);
        if (p[0] !== q[0] || p[1] !== q[1]) out.push(p[0], p[1], q[0], q[1]);
      }
    }
  }
  return new Float32Array(out);
}

// Levels at multiples of a step within the range of the values.
export function levels(values, step, base = 0) {
  let lo = Infinity, hi = -Infinity;
  for (const x of values) if (x === x) { if (x < lo) lo = x; if (x > hi) hi = x; }
  if (!(hi >= lo)) return [];
  const out = [];
  for (let k = Math.ceil((lo - base) / step); base + k * step <= hi; k++) out.push(base + k * step);
  return out;
}

// Grid points that are the highest (or lowest) within `radius` cells and differ
// from the mean of that neighbourhood by at least `prominence`. Rows within
// `polarRows` of either end are ignored, where the cells converge, and so are
// latitudes below `minAbsLat` (the tropics, where heights are flat).
export function extrema(values, grid, { radius = 6, prominence = 0, polarRows = 4, minAbsLat = 0 } = {}) {
  const { nlon, nlat } = grid;
  const found = [];
  for (let r = polarRows; r < nlat - polarRows; r++) {
    if (Math.abs(grid.lat0 - r * grid.dlat) < minAbsLat) continue;
    for (let c = 0; c < nlon; c++) {
      const x = values[r * nlon + c];
      if (x !== x) continue;
      let isMax = true, isMin = true, sum = 0, n = 0;
      for (let dr = -radius; dr <= radius && (isMax || isMin); dr++) {
        const rr = r + dr;
        if (rr < 0 || rr >= nlat) continue;
        for (let dc = -radius; dc <= radius; dc++) {
          let cc = c + dc;
          if (grid.periodic) cc = (cc + nlon) % nlon;
          else if (cc < 0 || cc >= nlon) continue;
          if (dr === 0 && dc === 0) continue;
          const y = values[rr * nlon + cc];
          if (y !== y) continue;
          if (y > x || (y === x && (dr < 0 || (dr === 0 && dc < 0)))) isMax = false;
          if (y < x || (y === x && (dr < 0 || (dr === 0 && dc < 0)))) isMin = false;
          sum += y; n++;
        }
      }
      if (!(isMax || isMin) || n === 0) continue;
      if (Math.abs(x - sum / n) < prominence) continue;
      found.push({ lon: grid.lon0 + c * grid.dlon, lat: grid.lat0 - r * grid.dlat, value: x, kind: isMax ? 'high' : 'low' });
    }
  }
  return found;
}
