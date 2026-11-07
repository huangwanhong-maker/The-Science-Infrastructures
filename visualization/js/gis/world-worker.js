// A tile of the world made ready for the scene, off the page's thread: its text read, its footprints triangulated and its
// positions taken to the local frame, exactly (frame.js), with the ground under each from the tile's grid of heights.
// The result is typed arrays handed over without a copy.
//
// Everything stands on the real ground, except within and around a city's own model, whose ground the Earth's shader
// makes level and feathers back into the terrain (`level`), the same weight being used here so that tiles meet it. A
// building's foot goes to the lowest ground under its footprint and its roof is measured from the ground at its middle.

import { localFrameAt } from './frame.js';
import { edges, ringArea, triangulate } from './geometry.js';

const STEP_M = 110;           // a way follows the ground at least this often
const PEOPLE = new Set(['pedestrian', 'footway', 'path', 'steps', 'living_street', 'cycleway', 'corridor']);
const MAJOR = new Set(['motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'motorway_link', 'trunk_link', 'primary_link', 'secondary_link']);
const DEFAULT_HEIGHT = 7;
const LARGEST_RING = 4000;    // a polygon's points beyond this (a coast within a tile) are thinned

const bytes = (rgb) => [Math.round(Math.min(1, rgb[0]) * 255), Math.round(Math.min(1, rgb[1]) * 255), Math.round(Math.min(1, rgb[2]) * 255), 255];
const mix = (a, b, t) => a.map((v, i) => v + (b[i] - v) * t);
const smooth = (a, b, x) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };

// The ground of a tile: bilinear on its grid (rows north to south, on latitude), and levelled as the Earth's shader levels it.
function groundOf(box, heights, frame, level) {
  const [s, w, n, e] = box, N = heights.n, G = heights.ground;
  const raw = (lon, lat) => {
    const u = Math.min(N - 1.0001, Math.max(0, ((lon - w) / (e - w)) * (N - 1))), v = Math.min(N - 1.0001, Math.max(0, ((n - lat) / (n - s)) * (N - 1)));
    const i = Math.floor(u), j = Math.floor(v), fu = u - i, fv = v - j, a = j * N + i;
    return (1 - fv) * ((1 - fu) * G[a] + fu * G[a + 1]) + fv * ((1 - fu) * G[a + N] + fu * G[a + N + 1]);
  };
  if (!level) return raw;
  return (lon, lat) => {
    const h = raw(lon, lat), p = frame.toLocal(lat, lon, h);
    const beyond = Math.max(Math.abs(p[0] - level.cx) - level.halfEast, Math.abs(p[1] - level.cy) - level.halfNorth);
    const k = (1 - smooth(0, level.feather, beyond)) * (level.amount ?? 1);
    return h + (level.height - h) * k;
  };
}

// The buildings' average height on a tile's grid, from the same grid.
function builtOf(box, heights) {
  if (!heights.built) return null;
  const [s, w, n, e] = box, N = heights.n, G = heights.built;
  return (lon, lat) => {
    const u = Math.round(Math.min(N - 1, Math.max(0, ((lon - w) / (e - w)) * (N - 1)))), v = Math.round(Math.min(N - 1, Math.max(0, ((n - lat) / (n - s)) * (N - 1))));
    return G[v * N + u];
  };
}

function ringOf(flat, frame) {
  const out = [];
  const step = flat.length / 2 > LARGEST_RING ? Math.ceil(flat.length / 2 / LARGEST_RING) : 1;
  for (let i = 0; i < flat.length; i += 2 * step) out.push(frame.toLocal(flat[i + 1], flat[i], 0));
  return out;
}
const inBox = (b, x, y) => b && x > b.x0 && x < b.x1 && y > b.y0 && y < b.y1;

function solids(data, frame, ground, built, keep, C) {
  const low = bytes(C.building), high = bytes(C.tall), pale = bytes(C.guessed);
  const kept = [];
  let vertices = 0, indices = 0, points = 0;
  const counts = { measured: 0, floors: 0, estimated: 0 };
  for (const b of data.buildings) {
    const ring = b[3], holes = (b[4] || []).filter((q) => q.length >= 6);
    const n = ring.length / 2;
    if (n < 3) continue;
    let lon = 0, lat = 0;
    for (let i = 0; i < ring.length; i += 2) lon += ring[i], lat += ring[i + 1];
    lon /= n, lat /= n;
    const middle = frame.toLocal(lat, lon, 0);
    if (inBox(keep, middle[0], middle[1])) continue;
    let total = n;
    for (const q of holes) total += q.length / 2;
    kept.push({ b, holes, lon, lat, total });
    vertices += 2 * total;
    indices += 3 * (total - 2 + 2 * holes.length) + 6 * total;
    points += n;
  }
  if (!kept.length) return { solids: null, count: 0, ...counts };
  const position = new Float32Array(3 * vertices), colour = new Uint8Array(4 * vertices), rise = new Float32Array(vertices), index = new Uint32Array(indices);
  const boxes = new Float32Array(5 * kept.length), rings = new Float32Array(2 * points), starts = new Uint32Array(2 * kept.length);
  const features = [];
  let v = 0, t = 0, r = 0, min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  kept.forEach(({ b, holes, lon, lat }, k) => {
    let height = b[0], code = b[1];
    if (code === 2) {
      const average = built ? built(lon, lat) : 0;
      if (average > 0) { height = average; counts.estimated++; } else { height = DEFAULT_HEIGHT; code = 3; counts.estimated++; }
    } else if (code === 1) counts.floors++;
    else counts.measured++;
    const all = [b[3], ...holes].map((flat) => {
      const out = [];
      for (let i = 0; i < flat.length; i += 2) out.push(frame.toLocal(flat[i + 1], flat[i], ground(flat[i], flat[i + 1])));
      return out;
    });
    let foot = Infinity;
    for (const q of all) for (const p of q) foot = Math.min(foot, p[2]);
    const top = frame.toLocal(lat, lon, ground(lon, lat))[2] + height;
    const z0 = foot + (b[2] || 0) - (b[2] ? 0 : 0.5), z1 = Math.max(z0 + 1, top);
    const rgb = code >= 2 ? pale : mix(low, high, Math.min(1, height / 120)).map(Math.round);
    const first = v;
    // Each point of each ring, at the foot and at the roof.
    let total = 0, cx = 0, cy = 0, radius = 0;
    for (const q of all) total += q.length;
    const wallStart = [];
    for (const q of all) {
      wallStart.push(v - first);
      for (const p of q) {
        position[3 * v] = p[0], position[3 * v + 1] = z0, position[3 * v + 2] = -p[1];
        position[3 * (v + total)] = p[0], position[3 * (v + total) + 1] = z1, position[3 * (v + total) + 2] = -p[1];
        rise[v] = 0, rise[v + total] = z1 - z0;
        colour.set(rgb, 4 * v);
        colour.set(rgb, 4 * (v + total));
        v++;
      }
    }
    v += total;
    // The walls.
    all.forEach((q, a) => {
      const s = first + wallStart[a], m = q.length;
      for (let i = 0; i < m; i++) {
        const j = (i + 1) % m;
        index[t++] = s + i, index[t++] = s + j, index[t++] = s + j + total;
        index[t++] = s + i, index[t++] = s + j + total, index[t++] = s + i + total;
      }
    });
    // The roof.
    const flat2 = all.map((q) => q.map((p) => [p[0], p[1]]));
    const tri = triangulate(flat2[0], flat2.slice(1));
    for (let i = 0; i < tri.length; i++) index[t++] = first + total + tri[i];
    // For picking.
    for (const p of all[0]) cx += p[0], cy += p[1];
    cx /= all[0].length, cy /= all[0].length;
    starts[2 * k] = r, starts[2 * k + 1] = all[0].length;
    for (const p of all[0]) { rings[2 * r] = p[0], rings[2 * r + 1] = p[1]; r++; radius = Math.max(radius, Math.hypot(p[0] - cx, p[1] - cy)); }
    boxes.set([cx, cy, radius, z0, z1], 5 * k);
    features.push([Math.round(height * 10) / 10, code, b[5] || null, b[6] || null]);
    for (const p of all[0]) {
      min[0] = Math.min(min[0], p[0]), max[0] = Math.max(max[0], p[0]);
      min[1] = Math.min(min[1], p[1]), max[1] = Math.max(max[1], p[1]);
    }
    min[2] = Math.min(min[2], z0), max[2] = Math.max(max[2], z1);
  });
  return {
    solids: { position, colour, rise, index: index.subarray(0, t), boxes, rings, starts, features, sphere: sphereOf(min, max) },
    count: kept.length, ...counts,
  };
}

// The bounding sphere in the scene's axes (x east, y up, z south).
function sphereOf(min, max) {
  const c = [(min[0] + max[0]) / 2, (min[2] + max[2]) / 2, -(min[1] + max[1]) / 2];
  return { centre: c, radius: Math.hypot(max[0] - min[0], max[1] - min[1], max[2] - min[2]) / 2 + 1 };
}

function ways(data, frame, ground, keep, C) {
  const paint = bytes(C.pedestrian), road = bytes(C.road), minor = bytes(C.minor), rail = bytes(C.rail);
  const rank = (r) => (r[4] ? 0 : PEOPLE.has(r[0]) ? 3 : MAJOR.has(r[0]) ? 2 : 1);
  const order = data.ways.map((r, i) => [rank(r), i]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const lines = [], all = [];
  let vertices = 0, indices = 0;
  for (const [, i] of order) {
    const r = data.ways[i], q = r[3];
    const line = [];
    for (let a = 0; a + 3 < q.length; a += 2) {
      const A = frame.toLocal(q[a + 1], q[a], 0), B = frame.toLocal(q[a + 3], q[a + 2], 0);
      const m = Math.max(1, Math.ceil(Math.hypot(A[0] - B[0], A[1] - B[1]) / STEP_M));
      for (let k = a === 0 ? 0 : 1; k <= m; k++) {
        const lon = q[a] + ((q[a + 2] - q[a]) * k) / m, lat = q[a + 1] + ((q[a + 3] - q[a + 1]) * k) / m;
        line.push(frame.toLocal(lat, lon, ground(lon, lat) + 0.2));
      }
    }
    if (line.length < 2) continue;
    const xyz = new Float32Array(3 * line.length);
    line.forEach((p, k) => { xyz[3 * k] = p[0]; xyz[3 * k + 1] = p[1]; xyz[3 * k + 2] = p[2]; });
    all.push({ cls: r[0], width: r[1], name: r[2], rail: !!r[4], xyz });
    const mid = line[Math.floor(line.length / 2)];
    if (inBox(keep, mid[0], mid[1])) continue;
    lines.push({ r, line });
    vertices += 2 * line.length;
    indices += 6 * (line.length - 1);
  }
  if (!lines.length) return { paint: null, all };
  const position = new Float32Array(3 * vertices), colour = new Uint8Array(4 * vertices), index = new Uint32Array(indices);
  let v = 0, t = 0;
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (const { r, line } of lines) {
    const { left, right } = edges(line.map((p) => [p[0], p[1]]), Math.max(1, r[1]));
    const rgb = r[4] ? rail : PEOPLE.has(r[0]) ? paint : MAJOR.has(r[0]) ? road : minor;
    const start = v;
    for (let i = 0; i < line.length; i++) {
      const z = line[i][2];
      for (const p of [left[i], right[i]]) {
        position[3 * v] = p[0], position[3 * v + 1] = z, position[3 * v + 2] = -p[1];
        colour.set(rgb, 4 * v);
        min[0] = Math.min(min[0], p[0]), max[0] = Math.max(max[0], p[0]), min[1] = Math.min(min[1], p[1]), max[1] = Math.max(max[1], p[1]), min[2] = Math.min(min[2], z), max[2] = Math.max(max[2], z);
        v++;
      }
    }
    for (let i = 0; i < line.length - 1; i++) {
      const a = start + 2 * i;
      index[t++] = a, index[t++] = a + 1, index[t++] = a + 3;
      index[t++] = a, index[t++] = a + 3, index[t++] = a + 2;
    }
  }
  return { paint: { position, colour, index, sphere: sphereOf(min, max) }, all };
}

function groundMesh(heights, box, frame, ground, keep) {
  const N = heights.n, [s, w, n, e] = box;
  const position = new Float32Array(3 * N * N), index = [];
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    const lat = n + ((s - n) * j) / (N - 1), lon = w + ((e - w) * i) / (N - 1);
    const p = frame.toLocal(lat, lon, ground(lon, lat));
    const a = 3 * (j * N + i);
    position[a] = p[0], position[a + 1] = p[2], position[a + 2] = -p[1];
    for (let k = 0; k < 3; k++) { min[k] = Math.min(min[k], p[k]); max[k] = Math.max(max[k], p[k]); }
  }
  for (let j = 0; j < N - 1; j++) for (let i = 0; i < N - 1; i++) {
    const a = j * N + i;
    if (keep && inBox(keep, position[3 * a], -position[3 * a + 2]) && inBox(keep, position[3 * (a + N + 1)], -position[3 * (a + N + 1) + 2])) continue;
    index.push(a, a + N, a + 1, a + 1, a + N, a + N + 1);
  }
  return { position, index: Uint32Array.from(index), sphere: sphereOf(min, max) };
}

function spaces(data, frame, ground, keep, C) {
  const out = [];
  const [bs, bw, bn, be] = data.box;
  const corner = (lat, lon) => frame.toLocal(lat, lon, 0);
  const tileArea = Math.abs(ringArea([corner(bs, bw), corner(bs, be), corner(bn, be), corner(bn, bw)]));
  let vertices = 0, indices = 0;
  for (const sp of data.spaces) {
    const rgb = bytes(C[sp[0]] || C.plaza);
    const outer = ringOf(sp[1], frame), holes = sp[2].filter((q) => q.length >= 6).map((q) => ringOf(q, frame));
    if (outer.length < 3) continue;
    // A sea, filling much of its tile, is left to the imagery beneath: painted, it would end in a straight line where the tiles do.
    if (sp[0] === 'water' && Math.abs(ringArea(outer)) > 0.3 * tileArea) continue;
    const pts = [outer, ...holes].map((q) => q.map((p) => [p[0], p[1]]));
    if (keep && outer.every((p) => inBox(keep, p[0], p[1]))) continue;
    const tri = triangulate(pts[0], pts.slice(1));
    if (!tri.length) continue;
    const raw = [sp[1], ...sp[2].filter((q) => q.length >= 6)];
    out.push({ rgb, pts, tri, raw });
    for (const q of pts) vertices += q.length;
    indices += tri.length;
  }
  if (!out.length) return null;
  const position = new Float32Array(3 * vertices), colour = new Uint8Array(4 * vertices), index = new Uint32Array(indices);
  let v = 0, t = 0;
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (const { rgb, pts, tri } of out) {
    const start = v;
    for (const q of pts) for (const [x, y] of q) {
      const g = frame.toGeodetic(x, y, 0);
      const z = frame.toLocal(g.lat, g.lon, ground(g.lon, g.lat) + 0.08)[2];
      position[3 * v] = x, position[3 * v + 1] = z, position[3 * v + 2] = -y;
      colour.set(rgb, 4 * v);
      min[0] = Math.min(min[0], x), max[0] = Math.max(max[0], x), min[1] = Math.min(min[1], y), max[1] = Math.max(max[1], y), min[2] = Math.min(min[2], z), max[2] = Math.max(max[2], z);
      v++;
    }
    for (let i = 0; i < tri.length; i++) index[t++] = start + tri[i];
  }
  return { position, colour, index, sphere: sphereOf(min, max) };
}

const buffersOf = (o) => (o ? Object.values(o).filter((v) => v && v.buffer instanceof ArrayBuffer).map((v) => v.buffer) : []);

self.onmessage = async ({ data: job }) => {
  const { kind, z, x, y, base, origin, level, keep, colours, generation } = job;
  try {
    const frame = localFrameAt(origin.lat, origin.lon, origin.height);
    const url = (what) => `${base}/world/${what}/${z}/${x}/${y}`;
    const answers = await Promise.all((kind === 'city' ? [url('city'), url('heights')] : [url('base'), url('heights')]).map(async (u) => {
      // Low priority, as the globe's tiles: what a person asks of an engine goes first when a connection frees.
      const r = await fetch(u, { priority: 'low' });
      if (!r.ok) throw new Error(`${u} answered ${r.status}`);
      return r.json();
    }));
    const [data, heights] = answers;
    const ground = groundOf(data.box, heights, frame, level);
    if (kind === 'base') {
      const paint = spaces(data, frame, ground, keep, colours);
      self.postMessage({ generation, paint }, buffersOf(paint));
      return;
    }
    const { solids: s, count, measured, floors, estimated } = solids(data, frame, ground, builtOf(data.box, heights), keep, colours);
    const { paint: r, all: lines } = ways(data, frame, ground, keep, colours);
    const g = groundMesh(heights, data.box, frame, ground, keep);
    self.postMessage({ generation, solids: s, ways: r, lines, ground: g, count, measured, floors, estimated, pending: heights.pending, credit: data.credit }, [...buffersOf(s), ...buffersOf(r), ...buffersOf(g), ...lines.map((l) => l.xyz.buffer)]);
  } catch (error) {
    self.postMessage({ generation, error: String(error && error.message || error) });
  }
};
