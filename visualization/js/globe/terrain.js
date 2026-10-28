// Terrain tiles: read from the data library's encoding, and built into meshes for drawing.
//
// A tile (data/tiles/terrain.py) holds (n + 2b)² heights and provenance codes: the tile's n x n
// samples at the corners of its cells, and a ring of b samples beyond its edges from which the
// normals at the edges are computed, so that neighbouring tiles shade alike across their seam.
// The mesh is built in the body-fixed frame, its positions relative to the tile's centre so that
// single precision keeps them to a millimetre (the centre itself is kept in double precision and
// the renderer subtracts the camera from it, the "relative to centre" method of Cozzi and Ring
// 2011, 3D Engine Design for Virtual Globes, §5.1). A skirt hangs from each edge, deep enough to
// hide the gap between a tile and a coarser neighbour (Ulrich 2002, Rendering Massive Terrains
// using Chunked Level of Detail Control, SIGGRAPH course notes), and from the two lines between its
// quarters, so that a quarter drawn beside a child that has taken the next quarter closes that
// seam too: the child's skirt hangs from the child's edge, and covers the gap only where the
// child stands higher.
//
// Provenance codes: 1 measured, 2 filled by the source's makers, 3 compiled (the source does not
// say), 4 finer than the data (interpolated here), 5 synthetic, 6 the sea's surface (the geoid).
//
// The sea. An Earth tile's heights are the surface one sees, the sea at its surface; the solid
// surface beneath, the sea floor, is a second layer of the tile ("floor"). A mesh is built of one or
// the other (buildMesh's `surface`); of the floor, the sea's samples are the source's (code 3).

import { cartesian, normal } from './figure.js';
import { cellDegrees, tileBounds } from './tiling.js';

const DEG = Math.PI / 180;

export const PROVENANCE = {
  1: 'measured',
  2: 'filled by the source',
  3: 'compiled source',
  4: 'finer than the data',
  5: 'synthetic',
  6: 'sea surface (geoid)',
};
export const SEA = 6;

export function decodeTile(buffer) {
  const bytes = new Uint8Array(buffer);
  if (bytes[0] !== 0x43 || bytes[1] !== 0x50 || bytes[2] !== 0x54 || bytes[3] !== 0x31) throw new Error('not a terrain tile (CPT1)');
  const view = new DataView(buffer);
  const length = view.getUint32(4, true);
  const header = JSON.parse(new TextDecoder().decode(bytes.subarray(8, 8 + length)));
  const m = header.n + 2 * (header.border || 0);
  const start = 8 + length;
  const heights = new Float32Array(buffer.slice(start, start + 4 * m * m));
  const end = start + 5 * m * m;
  const provenance = new Uint8Array(buffer.slice(start + 4 * m * m, end));
  let floor = null;
  if ((header.layers || []).includes('floor')) {
    const at = end + ((4 - (end % 4)) % 4);
    floor = new Float32Array(buffer.slice(at, at + 4 * m * m));
  }
  return { header, heights, provenance, floor, m };
}

// A flat tile at height zero, for drawing a body before its terrain arrives or where there is none.
export function flatTile(z, x, y, n = 65, border = 1) {
  const m = n + 2 * border;
  return { header: { z, x, y, n, border, min: 0, max: 0, sources: [], datum: 'none', flat: true }, heights: new Float32Array(m * m), provenance: new Uint8Array(m * m), floor: null, m };
}

// The mesh of a tile on a figure (figure.js), with `exaggeration` of the heights, of its visible
// surface or of its floor (`surface`: 'visible' or 'floor').
// Returns typed arrays for the GPU and the tile's centre (body-fixed, metres, double), its
// bounding radius, its geometric error (the cell's size, m) and points for the horizon test.
export function buildMesh(fig, tile, { exaggeration = 1, skirt = true, surface = 'visible' } = {}) {
  const { m } = tile;
  const ofFloor = surface === 'floor' && tile.floor;
  const heights = ofFloor ? tile.floor : tile.heights;
  const provenance = ofFloor ? tile.provenance.map((c) => (c === SEA ? 3 : c)) : tile.provenance;
  const header = ofFloor ? { ...tile.header, min: tile.header.floor_min, max: tile.header.floor_max } : tile.header;
  const n = header.n, b = header.border || 0;
  const [west, south, east, north] = tileBounds(header.z, header.x, header.y);
  const step = (east - west) / (n - 1);
  const at = (i, j) => heights[(j + b) * m + (i + b)] * exaggeration; // i east, j south; may reach into the border
  // The centre: the middle sample at the middle height.
  const midH = 0.5 * (header.min + header.max) * exaggeration;
  const centre = cartesian(fig, ((west + east) / 2) * DEG, ((south + north) / 2) * DEG, midH);
  const half = (n - 1) / 2, inner = Number.isInteger(half);
  const count = n * n, skirtCount = skirt ? 4 * n + (inner ? 2 * n : 0) : 0, total = count + skirtCount;
  const position = new Float32Array(total * 3);
  const normals = new Float32Array(total * 3);
  const uv = new Float32Array(total * 2);
  const prov = new Float32Array(total);
  const height = new Float32Array(total);
  const up = new Float32Array(total * 3); // the figure's normal: the local vertical, for day and night
  // Positions of the whole grid with its border, in double precision, for the normals.
  const full = new Float64Array(m * m * 3);
  for (let j = -b; j < n + b; j++) {
    let lat = north - step * j;
    if (lat > 90) lat = 180 - lat;
    if (lat < -90) lat = -180 - lat;
    for (let i = -b; i < n + b; i++) {
      cartesian(fig, (west + step * i) * DEG, lat * DEG, at(i, j), full, ((j + b) * m + (i + b)) * 3);
    }
  }
  const p = (i, j, k) => full[((j + b) * m + (i + b)) * 3 + k];
  let radius = 0;
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const v = j * n + i;
      for (let k = 0; k < 3; k++) position[3 * v + k] = p(i, j, k) - centre[k];
      radius = Math.max(radius, Math.hypot(position[3 * v], position[3 * v + 1], position[3 * v + 2]));
      // The normal from central differences across the cell, using the border at the edges; where
      // a tile has no border, one-sided.
      const i0 = Math.max(i - 1, -b), i1 = Math.min(i + 1, n - 1 + b), j0 = Math.max(j - 1, -b), j1 = Math.min(j + 1, n - 1 + b);
      const ex = [p(i1, j, 0) - p(i0, j, 0), p(i1, j, 1) - p(i0, j, 1), p(i1, j, 2) - p(i0, j, 2)];
      const ny = [p(i, j0, 0) - p(i, j1, 0), p(i, j0, 1) - p(i, j1, 1), p(i, j0, 2) - p(i, j1, 2)];
      let nx = ex[1] * ny[2] - ex[2] * ny[1], nyy = ex[2] * ny[0] - ex[0] * ny[2], nz = ex[0] * ny[1] - ex[1] * ny[0];
      let len = Math.hypot(nx, nyy, nz);
      if (!(len > 0)) {
        // At a pole the east difference vanishes: the figure's normal.
        const nn = normal((west + step * i) * DEG, (north - step * j) * DEG);
        [nx, nyy, nz] = nn;
        len = 1;
      }
      normals[3 * v] = nx / len, normals[3 * v + 1] = nyy / len, normals[3 * v + 2] = nz / len;
      uv[2 * v] = i / (n - 1);
      uv[2 * v + 1] = 1 - j / (n - 1);
      normal((west + step * i) * DEG, (north - step * j) * DEG, up, 3 * v);
      prov[v] = provenance[(j + b) * m + (i + b)];
      height[v] = heights[(j + b) * m + (i + b)];
    }
  }
  // Skirts: a copy of each edge sample lowered along the figure's normal.
  const cellM = cellDegrees(header.z, n) * DEG * fig.a;
  const depth = Math.max(cellM * 4, (header.max - header.min) * exaggeration * 0.5, 50);
  const edge = [];
  for (let i = 0; i < n; i++) edge.push([i, 0]);
  for (let j = 0; j < n; j++) edge.push([n - 1, j]);
  for (let i = n - 1; i >= 0; i--) edge.push([i, n - 1]);
  for (let j = n - 1; j >= 0; j--) edge.push([0, j]);
  // The lines between the quarters: the column i = half (north to south), then the row j = half
  // (west to east).
  const cross = [];
  if (inner) {
    for (let j = 0; j < n; j++) cross.push([half, j]);
    for (let i = 0; i < n; i++) cross.push([i, half]);
  }
  if (skirt) {
    edge.concat(cross).forEach(([i, j], s) => {
      const v = count + s, src = j * n + i;
      const nn = normal((west + step * i) * DEG, (north - step * j) * DEG);
      for (let k = 0; k < 3; k++) {
        position[3 * v + k] = position[3 * src + k] - depth * nn[k];
        normals[3 * v + k] = normals[3 * src + k];
        up[3 * v + k] = up[3 * src + k];
      }
      uv[2 * v] = uv[2 * src], uv[2 * v + 1] = uv[2 * src + 1];
      prov[v] = prov[src];
      height[v] = height[src];
    });
  }
  // Triangles: two per cell, split along alternating diagonals, and two per skirt segment, in the
  // order of the tile's quarters (north-west, north-east, south-west, south-east: its children's
  // order), each quarter's cells followed by its share of the skirts, so that a renderer can draw
  // any of the quarters alone (`quarters`: [start, count] of each in `index`).
  const perQuarter = [[], [], [], []];
  for (let j = 0; j < n - 1; j++) {
    for (let i = 0; i < n - 1; i++) {
      const a = j * n + i, bb = a + 1, c = a + n, d = c + 1;
      const q = (i >= half ? 1 : 0) + (j >= half ? 2 : 0);
      if ((i + j) % 2 === 0) perQuarter[q].push(a, c, bb, bb, c, d);
      else perQuarter[q].push(a, c, d, a, d, bb);
    }
  }
  if (skirt) {
    // The four sides, each n samples of the edge walk with their lowered copies: north (west to
    // east), east (north to south), south (east to west), west (south to north).
    const sideQuarters = [[0, 1], [1, 3], [3, 2], [2, 0]];
    for (let side = 0; side < 4; side++) {
      for (let k = 0; k < n - 1; k++) {
        const s0 = side * n + k, s1 = s0 + 1;
        const [i0, j0] = edge[s0], [i1, j1] = edge[s1];
        const top0 = j0 * n + i0, top1 = j1 * n + i1, low0 = count + s0, low1 = count + s1;
        perQuarter[sideQuarters[side][k < half ? 0 : 1]].push(top0, low0, top1, top1, low0, low1);
      }
    }
    if (inner) {
      // The walls under the lines between the quarters, each in both quarters beside it: drawn
      // whole, the tile hides them beneath its surface.
      for (let k = 0; k < n - 1; k++) {
        const c0 = k * n + half, c1 = c0 + n, cl0 = count + 4 * n + k, cl1 = cl0 + 1;
        for (const q of k < half ? [0, 1] : [2, 3]) perQuarter[q].push(c0, cl0, c1, c1, cl0, cl1);
        const r0 = half * n + k, r1 = r0 + 1, rl0 = count + 5 * n + k, rl1 = rl0 + 1;
        for (const q of k < half ? [0, 2] : [1, 3]) perQuarter[q].push(r0, rl0, r1, r1, rl0, rl1);
      }
    }
  }
  const index = new Uint16Array(perQuarter.reduce((m, q) => m + q.length, 0));
  const quarters = [];
  let t = 0;
  for (const q of perQuarter) {
    quarters.push([t, q.length]);
    index.set(q, t);
    t += q.length;
  }
  // Points for the horizon test: a 5 x 5 grid over the tile at its top (lod.js stands them by one).
  const top = header.max * exaggeration;
  const probes = [];
  for (const fx of [0, 0.25, 0.5, 0.75, 1]) for (const fy of [0, 0.25, 0.5, 0.75, 1]) probes.push(cartesian(fig, (west + (east - west) * fx) * DEG, (north - (north - south) * fy) * DEG, top));
  return {
    centre, radius, probes, geometricError: cellM, min: header.min * exaggeration, max: header.max * exaggeration, exaggeration, surface: ofFloor ? 'floor' : 'visible',
    position, normal: normals, up, uv, provenance: prov, height, index, quarters,
  };
}
