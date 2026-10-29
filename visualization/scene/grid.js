// The samplings on which a field is read (physis fields/sampling.hpp), for the drawings of its samples:
// the same places and the same orders, so that a value and the place it is drawn at agree.
//
//   sphere   {centre, radius_m, nlat, nlon}   row i at latitude 90° − (i + ½)·180°/nlat, column j at
//                                             longitude −180° + (j + ½)·360°/nlon, index i·nlon + j
//                                             (an image's order: north up, west left);
//   plane    {origin, u, v, nu, nv}           cell (i, j) at origin + (i + ½)/nu u + (j + ½)/nv v,
//                                             index i + nu·j;
//   volume   {centre, axes, edges, n}         cell (i, j, k) offset ((i + ½)/n₀ − ½) e₀ along the first
//                                             axis (a row of `axes`) and so on, index i + n₀ (j + n₁ k);
//   points   {points}                         in their order.
//
// Points are on the axes of the field's frame, in metres. A sampling given by the engine's answer is
// its `sampling` member (fields.hpp, sampling_json); a volume without axes has the frame's.

const IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

export function sizeOf(g) {
  switch (g.kind) {
    case 'sphere': return g.nlat * g.nlon;
    case 'plane': return g.nu * g.nv;
    case 'volume': return g.n[0] * g.n[1] * g.n[2];
    case 'points': return g.points ? g.points.length : g.count;
    default: throw new Error(`no sampling of kind "${g.kind}"`);
  }
}

// ---- sphere ----

export const sphereLatitude = (g, i) => Math.PI / 2 - ((i + 0.5) * Math.PI) / g.nlat;
export const sphereLongitude = (g, j) => -Math.PI + ((j + 0.5) * 2 * Math.PI) / g.nlon;
export const sphereCell = (g, index) => [Math.floor(index / g.nlon), index % g.nlon];

// ---- volume ----

export const volumeIndex = (g, i, j, k) => i + g.n[0] * (j + g.n[1] * k);
export function volumeCell(g, index) {
  const i = index % g.n[0], rest = Math.floor(index / g.n[0]);
  return [i, rest % g.n[1], Math.floor(rest / g.n[1])];
}
// The offset of a cell's centre from the box's centre, along the box's axes (m).
export const volumeOffset = (g, i, j, k) => [((i + 0.5) / g.n[0] - 0.5) * g.edges[0], ((j + 0.5) / g.n[1] - 0.5) * g.edges[1], ((k + 0.5) / g.n[2] - 0.5) * g.edges[2]];
// A point given along the box's axes, on the frame's axes: centre + axesᵀ o.
export function fromAxes(g, o) {
  const A = g.axes || IDENTITY, c = g.centre || [0, 0, 0];
  return [0, 1, 2].map((x) => c[x] + A[0][x] * o[0] + A[1][x] * o[1] + A[2][x] * o[2]);
}

// ---- any sampling ----

// The place of a point of the sampling on the frame's axes (m).
export function pointOf(g, index) {
  switch (g.kind) {
    case 'sphere': {
      const [i, j] = sphereCell(g, index), p = sphereLatitude(g, i), l = sphereLongitude(g, j), c = g.centre || [0, 0, 0];
      return [c[0] + g.radius_m * Math.cos(p) * Math.cos(l), c[1] + g.radius_m * Math.cos(p) * Math.sin(l), c[2] + g.radius_m * Math.sin(p)];
    }
    case 'plane': {
      const i = index % g.nu, j = Math.floor(index / g.nu), s = (i + 0.5) / g.nu, r = (j + 0.5) / g.nv;
      return [0, 1, 2].map((x) => g.origin[x] + s * g.u[x] + r * g.v[x]);
    }
    case 'volume': return fromAxes(g, volumeOffset(g, ...volumeCell(g, index)));
    case 'points': return g.points[index];
    default: throw new Error(`no sampling of kind "${g.kind}"`);
  }
}
