// Any city in the world: as the camera comes near the ground, the tiles under it are asked for and drawn, buildings, ways, water
// and open land, each standing on the real ground, so that no place has to be prepared. The tiles are of the Web Mercator
// pyramid's level 14 (about 2.4 km at the equator); the server reduces them from a world's archives (data/acquisition/world_city.py),
// a pool of workers (world-worker.js) turns each into geometry in the scene's local frame, and each is added in a frame of its own,
// nearest first, rising from the ground as it comes. Tiles far behind are let go, and a limit on the buildings held
// keeps a dense region within a laptop's means.
//
//     const world = createWorld({ scene, frame, level, keep, colours, base, onStatus });
//     each few frames: world.want(lat, lon, range);
//
// The caller gives what changes with the place: `scene()` the city scene (renderer-three.js), `frame()` the frame's origin
// { lat, lon, height }, `level()` the ground a city's model stands on, levelled in the Earth's shader, so that tiles beyond it meet it,
// `keep()` the box of the local frame that a model of its own already fills, `colours()` the palette.

export const TILE_LEVEL = 14;
const BASE_LEVEL = 13;

// The tile holding a place, and what a tile covers [south, west, north, east] (degrees).
export function tileOf(lat, lon, z = TILE_LEVEL) {
  const n = 2 ** z, r = (lat * Math.PI) / 180;
  return [Math.floor(((lon + 180) / 360) * n), Math.floor(((1 - Math.asinh(Math.tan(r)) / Math.PI) / 2) * n)];
}
export function tileBox(x, y, z = TILE_LEVEL) {
  const n = 2 ** z, lat = (t) => (Math.atan(Math.sinh(Math.PI * (1 - (2 * t) / n))) * 180) / Math.PI;
  return [lat(y + 1), (x / n) * 360 - 180, lat(y), ((x + 1) / n) * 360 - 180];
}

// How many tiles out from the middle a view asks for: enough to fill what is seen, at most `most`.
export function reachOf(range, lat, most = 3) {
  const size = (40075016.7 / 2 ** TILE_LEVEL) * Math.cos((lat * Math.PI) / 180);
  return Math.min(most, Math.max(1, Math.ceil(Math.min(range * 1.25, 12000) / Math.max(400, size))));
}

export function createWorld({ scene, frame, level = () => null, keep = () => null, colours, form = () => 1, base = '', onStatus = null, workers: count = 3, maxBuildings = 420000, near = 26000 }) {
  const pool = Array.from({ length: count }, () => ({ worker: new Worker(new URL('./world-worker.js', import.meta.url), { type: 'module' }), busy: false }));
  const tiles = new Map();     // 'x/y' → { x, y, state, ids, buildings, measured, floors, estimated, tries }
  const bases = new Map();
  const arrivals = [];         // what is ready to be added, a tile a frame
  let generation = 0, centre = null, wantedReach = 1, idle = true;
  const status = { tiles: 0, buildings: 0, measured: 0, floors: 0, estimated: 0, pending: 0, failed: 0, credit: null };

  const tell = () => {
    status.pending = [...tiles.values(), ...bases.values()].filter((t) => t.state === 'queued' || t.state === 'loading').length + arrivals.length;
    onStatus && onStatus(status);
  };
  const distance = (t) => (centre ? Math.hypot(t.x / (t.scale || 1) - centre[0], t.y / (t.scale || 1) - centre[1]) : 0);

  function jobOf(t, kind) {
    return { kind, z: kind === 'city' ? TILE_LEVEL : BASE_LEVEL, x: t.x, y: t.y, base, origin: frame(), level: level(), keep: keep(), colours: colours(), generation };
  }

  function pump() {
    for (const w of pool) {
      if (w.busy) continue;
      // The nearest thing wanted; what has gone out of reach while it waited is not asked for.
      let next = null, best = Infinity, kind = null;
      for (const [map, k] of [[tiles, 'city'], [bases, 'base']]) for (const t of map.values()) {
        if (t.state !== 'queued') continue;
        const d = distance(t) + (k === 'base' ? 0.25 : 0);
        if (k === 'city' && d > 2.2 && status.buildings > maxBuildings) continue;
        if (d < best) { best = d; next = t; kind = k; }
      }
      if (!next) break;
      next.state = 'loading';
      w.busy = true;
      w.job = next;
      w.kind = kind;
      w.worker.postMessage(jobOf(next, kind));
    }
    tell();
  }

  for (const w of pool) {
    w.worker.onmessage = ({ data }) => {
      w.busy = false;
      const t = w.job, kind = w.kind;
      w.job = null;
      if (data.generation === generation && (kind === 'city' ? tiles : bases).get(`${t.x}/${t.y}`) === t) {
        if (data.error) {
          t.state = 'failed';
          t.error = data.error;
          status.failed++;
          console.warn(`world tile ${t.x}/${t.y}: ${data.error}`);
        } else arrivals.push({ t, kind, data });
      }
      pump();
    };
    w.worker.onerror = (e) => { console.warn(`world worker: ${e.message}`); w.busy = false; if (w.job) w.job.state = 'failed'; w.job = null; pump(); };
  }

  function drop(t) {
    const gis = scene();
    for (const id of t.ids || []) gis && gis.remove(id);
    t.ids = [];
    status.buildings -= t.buildings || 0, status.measured -= t.measured || 0, status.floors -= t.floors || 0, status.estimated -= t.estimated || 0;
    t.buildings = t.measured = t.floors = t.estimated = 0;
  }

  function add({ t, kind, data }) {
    const gis = scene();
    if (!gis) return;
    const key = `${t.x}/${t.y}`;
    if (kind === 'base') {
      if (data.paint) {
        gis.addBuiltPaint(`w:s:${key}`, data.paint, { order: -4 });
        t.ids = [`w:s:${key}`];
      }
      t.state = 'ready';
      return;
    }
    const again = !!t.waiting;
    drop(t);
    const ids = [];
    if (data.ground) { gis.addBuiltGround(`w:g:${key}`, data.ground); ids.push(`w:g:${key}`); }
    if (data.ways) { gis.addBuiltPaint(`w:r:${key}`, data.ways, { order: -2 }); ids.push(`w:r:${key}`); }
    if (data.solids) { gis.addBuiltSolids(`w:b:${key}`, data.solids, { facade: colours().facade, grow: !again, form: form(), pick: (f) => describe(f) }); ids.push(`w:b:${key}`); }
    t.ids = ids;
    t.state = 'ready';
    t.lines = data.lines;
    t.buildings = data.count, t.measured = data.measured, t.floors = data.floors, t.estimated = data.estimated;
    status.buildings += t.buildings, status.measured += t.measured, status.floors += t.floors, status.estimated += t.estimated;
    status.tiles = [...tiles.values()].filter((q) => q.state === 'ready' || q.state === 'again').length;
    status.credit = data.credit;
    onAdded && onAdded(ids);
    // The buildings' average heights of a region come once, some seconds after it is first seen: the tile is then built again.
    if (data.pending) {
      const gen = generation;
      setTimeout(() => {
        if (gen !== generation || tiles.get(key) !== t) return;
        t.state = 'queued';
        t.waiting = true;
        pump();
      }, 7000);
    }
    t.waiting = false;
  }
  let onAdded = null;

  function describe(f) {
    const how = { 0: 'measured', 1: 'from its floors', 2: 'estimated', 3: 'a default' }[f[1]];
    return `${f[2] || (f[3] && f[3] !== 'building' ? f[3].replace(/_/g, ' ') : 'A building')}\n${Math.round(f[0])} m (${how})`;
  }

  return {
    // Ask for what is under a view: the point of regard (degrees) and the camera's range (m).
    want(lat, lon, range) {
      if (range > near) { centre = null; tell(); return; }
      const [cx, cy] = tileOf(lat, lon);
      centre = [cx, cy];
      const reach = reachOf(range, lat);
      wantedReach = reach;
      for (let i = -reach; i <= reach; i++) for (let j = -reach; j <= reach; j++) {
        if (Math.hypot(i, j) > reach + 0.5) continue;
        const x = cx + i, y = cy + j, key = `${x}/${y}`;
        if (x < 0 || y < 0 || x >= 2 ** TILE_LEVEL || y >= 2 ** TILE_LEVEL) continue;
        if (!tiles.has(key)) tiles.set(key, { x, y, state: 'queued', ids: [] });
        const bx = Math.floor(x / 2), by = Math.floor(y / 2), bkey = `${bx}/${by}`;
        if (!bases.has(bkey)) bases.set(bkey, { x: bx, y: by, scale: 0.5, state: 'queued', ids: [] });
      }
      // Far behind: let go. What is held counts toward the limit, so a long walk across a city stays within it.
      for (const [key, t] of tiles) if (Math.hypot(t.x - cx, t.y - cy) > reach + 2.5) { drop(t); tiles.delete(key); }
      for (const [key, t] of bases) if (Math.hypot(t.x / t.scale - cx, t.y / t.scale - cy) > reach + 4) { const gis = scene(); for (const id of t.ids) gis && gis.remove(id); bases.delete(key); }
      for (const a of arrivals.splice(0)) if ((a.kind === 'city' ? tiles : bases).get(`${a.t.x}/${a.t.y}`) === a.t) arrivals.push(a);
      pump();
    },
    // Each frame: what is ready goes into the scene, one at a time, so that a frame is never long.
    tick() {
      const next = arrivals.shift();
      if (next) { add(next); tell(); }
      idle = !arrivals.length;
    },
    // Another place, another palette: what was asked for is forgotten, and asked again from the server's keep.
    reset() {
      generation++;
      arrivals.length = 0;
      tiles.clear();
      bases.clear();
      centre = null;
      Object.assign(status, { tiles: 0, buildings: 0, measured: 0, floors: 0, estimated: 0, pending: 0, failed: 0 });
      tell();
    },
    onAdded(f) { onAdded = f; },
    status: () => ({ ...status }),
    // Every way of the tiles held, whether or not a model covers it: { cls, width, name, rail, xyz (local, m) }.
    lines() {
      const out = [];
      for (const t of tiles.values()) if (t.lines) out.push(...t.lines);
      return out;
    },
    get reach() { return wantedReach; },
    get idle() { return idle; },
    tiles,
    dispose() { for (const w of pool) w.worker.terminate(); },
  };
}
