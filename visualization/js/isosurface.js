// A surface from a volume: the isosurface at a level of a scalar field on a
// voxel grid, as a triangle mesh in world coordinates.
//
// The method is the surface net of Gibson (1998), in its naive form (Lysenko
// 2012): each cell of the grid whose eight corners are not all on one side of
// the level gets one vertex, the mean of the points where the level crosses
// the cell's edges; every grid edge that crosses the level joins the four cells
// around it by a quad, split into two triangles, wound so that the normals
// point toward the lower side (outward, for a mask whose inside is 1). It
// needs no case tables, closes the surface, and gives smoother meshes than
// marching cubes, which the optional relaxation (a few steps of Laplacian
// smoothing of the vertices toward their neighbours' mean) improves further.
//
// Input: a volume {dims, data, affine} as the domains represent it (data in
// Fortran order, affine voxel-to-world), the level, and options. Output:
// {positions: Float32Array (3 per vertex), indices: Uint32Array (3 per
// triangle), normals: Float32Array}. Any application may draw it.
//
// References: Gibson, S. F. F. (1998), Constrained elastic surface nets:
// generating smooth surfaces from binary segmented data, MICCAI 1998, LNCS
// 1496, 888-898; Lysenko, M. (2012), Smooth voxel terrain (0fps.net).

const EDGES = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];   // the 12 edges of a cell, corners numbered by bits (i, j, k)

export function isosurface(volume, level = 0.5, { smoothing = 2, smoothingWeight = 0.5 } = {}) {
  const [nx, ny, nz] = volume.dims;
  const d = volume.data, A = volume.affine;
  const cellIndex = new Int32Array((nx - 1) * (ny - 1) * (nz - 1)).fill(-1);
  const cell = (i, j, k) => i + (nx - 1) * (j + (ny - 1) * k);
  const vertices = [];                 // voxel coordinates, flat
  const at = (i, j, k) => d[i + nx * (j + ny * k)];
  // 1. a vertex for every cell the surface crosses
  for (let k = 0; k < nz - 1; k++) for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
    const c = [at(i, j, k), at(i + 1, j, k), at(i, j + 1, k), at(i + 1, j + 1, k), at(i, j, k + 1), at(i + 1, j, k + 1), at(i, j + 1, k + 1), at(i + 1, j + 1, k + 1)];
    let mask = 0;
    for (let n = 0; n < 8; n++) if (c[n] >= level) mask |= 1 << n;
    if (mask === 0 || mask === 255) continue;
    let sx = 0, sy = 0, sz = 0, count = 0;
    for (const [a, b] of EDGES) {
      const inA = (mask >> a) & 1, inB = (mask >> b) & 1;
      if (inA === inB) continue;
      const t = (level - c[a]) / (c[b] - c[a]);
      sx += (a & 1) + ((b & 1) - (a & 1)) * t;
      sy += ((a >> 1) & 1) + (((b >> 1) & 1) - ((a >> 1) & 1)) * t;
      sz += ((a >> 2) & 1) + (((b >> 2) & 1) - ((a >> 2) & 1)) * t;
      count++;
    }
    cellIndex[cell(i, j, k)] = vertices.length / 3;
    vertices.push(i + sx / count, j + sy / count, k + sz / count);
  }
  // 2. a quad for every edge the surface crosses, between the four cells around it
  const indices = [];
  const quad = (a, b, c, e, flip) => { if (flip) indices.push(a, b, c, a, c, e); else indices.push(a, c, b, a, e, c); };
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v0 = at(i, j, k) >= level;
    // the edge along x from (i, j, k): cells (i, j-1..j, k-1..k)
    if (i < nx - 1 && j > 0 && k > 0 && j < ny - 1 && k < nz - 1) {
      const v1 = at(i + 1, j, k) >= level;
      if (v0 !== v1) quad(cellIndex[cell(i, j - 1, k - 1)], cellIndex[cell(i, j, k - 1)], cellIndex[cell(i, j, k)], cellIndex[cell(i, j - 1, k)], v0);
    }
    if (j < ny - 1 && i > 0 && k > 0 && i < nx - 1 && k < nz - 1) {
      const v1 = at(i, j + 1, k) >= level;
      if (v0 !== v1) quad(cellIndex[cell(i - 1, j, k - 1)], cellIndex[cell(i - 1, j, k)], cellIndex[cell(i, j, k)], cellIndex[cell(i, j, k - 1)], v0);
    }
    if (k < nz - 1 && i > 0 && j > 0 && i < nx - 1 && j < ny - 1) {
      const v1 = at(i, j, k + 1) >= level;
      if (v0 !== v1) quad(cellIndex[cell(i - 1, j - 1, k)], cellIndex[cell(i, j - 1, k)], cellIndex[cell(i, j, k)], cellIndex[cell(i - 1, j, k)], v0);
    }
  }
  let pos = Float32Array.from(vertices);
  const nv = pos.length / 3;
  // 3. relaxation: each vertex toward the mean of its neighbours
  if (smoothing > 0 && nv) {
    const neighbours = Array.from({ length: nv }, () => new Set());
    for (let t = 0; t < indices.length; t += 3) {
      const [a, b, c] = [indices[t], indices[t + 1], indices[t + 2]];
      neighbours[a].add(b).add(c); neighbours[b].add(a).add(c); neighbours[c].add(a).add(b);
    }
    for (let step = 0; step < smoothing; step++) {
      const next = new Float32Array(pos.length);
      for (let v = 0; v < nv; v++) {
        const nb = neighbours[v];
        if (!nb.size) { next.set(pos.subarray(3 * v, 3 * v + 3), 3 * v); continue; }
        let mx = 0, my = 0, mz = 0;
        for (const u of nb) { mx += pos[3 * u]; my += pos[3 * u + 1]; mz += pos[3 * u + 2]; }
        const w = smoothingWeight;
        next[3 * v] = pos[3 * v] * (1 - w) + (mx / nb.size) * w;
        next[3 * v + 1] = pos[3 * v + 1] * (1 - w) + (my / nb.size) * w;
        next[3 * v + 2] = pos[3 * v + 2] * (1 - w) + (mz / nb.size) * w;
      }
      pos = next;
    }
  }
  // 4. to world coordinates, with normals from the faces
  const positions = new Float32Array(pos.length);
  for (let v = 0; v < nv; v++) {
    const x = pos[3 * v], y = pos[3 * v + 1], z = pos[3 * v + 2];
    positions[3 * v] = A[0] * x + A[1] * y + A[2] * z + A[3];
    positions[3 * v + 1] = A[4] * x + A[5] * y + A[6] * z + A[7];
    positions[3 * v + 2] = A[8] * x + A[9] * y + A[10] * z + A[11];
  }
  const index = Uint32Array.from(indices);
  return { positions, indices: index, normals: vertexNormals(positions, index), vertexCount: nv, triangleCount: index.length / 3 };
}

// Per-vertex normals: the sum of the face normals (weighted by area), normalised.
export function vertexNormals(positions, indices) {
  const normals = new Float32Array(positions.length);
  for (let t = 0; t < indices.length; t += 3) {
    const a = indices[t], b = indices[t + 1], c = indices[t + 2];
    const ax = positions[3 * a], ay = positions[3 * a + 1], az = positions[3 * a + 2];
    const ux = positions[3 * b] - ax, uy = positions[3 * b + 1] - ay, uz = positions[3 * b + 2] - az;
    const vx = positions[3 * c] - ax, vy = positions[3 * c + 1] - ay, vz = positions[3 * c + 2] - az;
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    for (const v of [a, b, c]) { normals[3 * v] += nx; normals[3 * v + 1] += ny; normals[3 * v + 2] += nz; }
  }
  for (let v = 0; v < normals.length; v += 3) {
    const n = Math.hypot(normals[v], normals[v + 1], normals[v + 2]) || 1;
    normals[v] /= n; normals[v + 1] /= n; normals[v + 2] /= n;
  }
  return normals;
}

// The area of a mesh and the volume it encloses (by the divergence theorem;
// positive when the faces are wound outward), and its Euler characteristic
// V - E + F (2 for a closed surface without holes or handles).
export function meshMeasures(positions, indices) {
  let area = 0, volume = 0;
  const edges = new Set();
  for (let t = 0; t < indices.length; t += 3) {
    const a = indices[t], b = indices[t + 1], c = indices[t + 2];
    const ax = positions[3 * a], ay = positions[3 * a + 1], az = positions[3 * a + 2];
    const bx = positions[3 * b], by = positions[3 * b + 1], bz = positions[3 * b + 2];
    const cx = positions[3 * c], cy = positions[3 * c + 1], cz = positions[3 * c + 2];
    const nx = (by - ay) * (cz - az) - (bz - az) * (cy - ay), ny = (bz - az) * (cx - ax) - (bx - ax) * (cz - az), nz = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    area += Math.hypot(nx, ny, nz) / 2;
    volume += (ax * nx + ay * ny + az * nz) / 6;
    for (const [p, q] of [[a, b], [b, c], [c, a]]) edges.add(p < q ? `${p},${q}` : `${q},${p}`);
  }
  const vertices = new Set(indices).size;
  return { area, volume, vertices, edges: edges.size, faces: indices.length / 3, euler: vertices - edges.size + indices.length / 3 };
}
