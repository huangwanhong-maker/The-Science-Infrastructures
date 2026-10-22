// Surfaces, away from the page's main thread: the brain's outline from its
// mask, a region's surface from the label index, and the scalp from the head.
//
// Message in:  { id, kind: 'mask' | 'region', dims, affine, data (Uint8Array mask,
//                or the Uint16Array label index with `indices`, the set of index
//                values that make the region), step (1 or 2: subsample the grid),
//                smoothing }
// Message out: { id, positions, normals, indices } or { id, error }

import { isosurface } from '../../../../visualization/js/index.js';

function subsample(dims, data, affine, step, test) {
  const nx = Math.ceil(dims[0] / step), ny = Math.ceil(dims[1] / step), nz = Math.ceil(dims[2] / step);
  // One voxel of margin so that the surface closes at the edges.
  const out = new Uint8Array((nx + 2) * (ny + 2) * (nz + 2));
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v = data[i * step + dims[0] * (j * step + dims[1] * k * step)];
    if (test(v)) out[(i + 1) + (nx + 2) * ((j + 1) + (ny + 2) * (k + 1))] = 1;
  }
  // The affine of the coarse grid: voxel (a, b, c) of it is voxel (step (a - 1), ...) of the fine one.
  const A = affine;
  const coarse = [
    A[0] * step, A[1] * step, A[2] * step, A[3] - (A[0] + A[1] + A[2]) * step,
    A[4] * step, A[5] * step, A[6] * step, A[7] - (A[4] + A[5] + A[6]) * step,
    A[8] * step, A[9] * step, A[10] * step, A[11] - (A[8] + A[9] + A[10]) * step,
    0, 0, 0, 1,
  ];
  return { dims: [nx + 2, ny + 2, nz + 2], data: out, affine: coarse };
}

self.onmessage = (event) => {
  const { id, kind, dims, affine, data, indices, step = 1, smoothing = 2 } = event.data;
  try {
    const test = kind === 'region' ? ((v) => indices.includes(v)) : ((v) => v > 0);
    const set = new Set(indices || []);
    const vol = subsample(dims, data, affine, step, kind === 'region' ? ((v) => set.has(v)) : test);
    const mesh = isosurface(vol, 0.5, { smoothing });
    self.postMessage({ id, positions: mesh.positions, normals: mesh.normals, indices: mesh.indices }, [mesh.positions.buffer, mesh.normals.buffer, mesh.indices.buffer]);
  } catch (e) {
    self.postMessage({ id, error: e.message });
  }
};
