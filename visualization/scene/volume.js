// The arithmetic of a field read in a volume (a VolumeGrid, grid.js): the cells a way of showing keeps,
// the linear fit of the values over the box (their gradient along its axes), the value at its centre,
// and the data of a level surface. Pure functions of the values and the grid, in the grid's axes and
// metres; a renderer places them.

import { volumeCell, volumeIndex, volumeOffset } from './grid.js';

// The cells a mode shows: every cell with a value ('cells'), the box's outer layer ('shell'), or three
// slices through the cell `slices` = [i, j, k] ('slices'); none for 'level', which draws a surface.
export function visibleCells(values, g, mode, slices = [0, 0, 0]) {
  const out = [];
  const [n0, n1, n2] = g.n, size = n0 * n1 * n2;
  for (let idx = 0; idx < size; idx++) {
    if (!Number.isFinite(values[idx])) continue;
    const [i, j, k] = volumeCell(g, idx);
    if (mode === 'shell' && !(i === 0 || j === 0 || k === 0 || i === n0 - 1 || j === n1 - 1 || k === n2 - 1)) continue;
    if (mode === 'slices' && !(i === slices[0] || j === slices[1] || k === slices[2])) continue;
    if (mode === 'level') continue;
    out.push(idx);
  }
  return out;
}

// The least-squares fit v ≈ a + g · x over the cells with a value, x the offset from the centre along the
// grid's axes (m): { mean: a, gradient: [g₀, g₁, g₂] }, or null with fewer than eight cells or a
// singular fit.
export function gradientOf(values, g) {
  const M = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], b = [0, 0, 0, 0];
  let count = 0;
  const size = g.n[0] * g.n[1] * g.n[2];
  for (let idx = 0; idx < size; idx++) {
    const v = values[idx];
    if (!Number.isFinite(v)) continue;
    const x = [1, ...volumeOffset(g, ...volumeCell(g, idx))];
    for (let p = 0; p < 4; p++) { b[p] += x[p] * v; for (let q = 0; q < 4; q++) M[p][q] += x[p] * x[q]; }
    count++;
  }
  if (count < 8) return null;
  // Gauss–Jordan elimination with partial pivoting on the 4 × 4 normal equations.
  const A = M.map((row, p) => [...row, b[p]]);
  for (let c = 0; c < 4; c++) {
    let piv = c;
    for (let r = c + 1; r < 4; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r;
    if (Math.abs(A[piv][c]) < 1e-30) return null;
    [A[c], A[piv]] = [A[piv], A[c]];
    for (let r = 0; r < 4; r++) {
      if (r === c) continue;
      const f = A[r][c] / A[c][c];
      for (let q = c; q < 5; q++) A[r][q] -= f * A[c][q];
    }
  }
  const sol = A.map((row, p) => row[4] / row[p]);
  return { mean: sol[0], gradient: [sol[1], sol[2], sol[3]] };
}

// The value at the box's centre: its middle cell's when every count is odd; else the mean of the cells
// about the centre (two along each even axis) that have a value.
export function centreValue(values, g) {
  const around = g.n.map((n) => (n % 2 ? [(n - 1) / 2] : [n / 2 - 1, n / 2]));
  let s = 0, c = 0;
  for (const i of around[0]) for (const j of around[1]) for (const k of around[2]) {
    const v = values[volumeIndex(g, i, j, k)];
    if (Number.isFinite(v)) { s += v; c++; }
  }
  return c ? s / c : NaN;
}

// The data of a level surface (visualization/js/isosurface.js): cells without a value counted below every
// level, as finite values, so that the surface's crossings interpolate and it closes on them.
export function levelData(values) {
  let lo = Infinity, hi = -Infinity;
  for (const v of values) if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  const below = lo - (hi - lo || 1);
  return Float32Array.from(values, (v) => (Number.isFinite(v) ? v : below));
}
