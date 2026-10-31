// What the plug-in's layers and analyses share: the district (its walking network cut from the city model in a
// box), the two stores, the lane, the day's run and the building, with listeners, so that the map and the drawer
// show one landscape. A module is loaded once by its URL, so every layer and tool holds this same object.
//
// The network: the city model's ways where people walk (streets, lanes, malls; not motorways, not the separate
// footway lines drawn beside streets, which would double every street), cut to the box, split at every point two
// ways share and at their ends: the junctions; each piece between two junctions a segment. A way that closes on
// itself is split once more so that no segment is a loop. Ways that cross without a shared point are left to
// the engine, which drops the shorter of each crossing pair and names it.
//
// The forcing of a day: the Sun's altitude and azimuth at each hour of the clock's date in the city's own time
// (the astronomy domain's subsolar point, as the clock lights the city), and the walkers entering the district
// each hour, half the hourly counts at the City of Melbourne's sensors inside the box on the clock's weekday
// (Street Oceanography's profiles), or a declared day's shape when they are not at hand.

import { instantOf } from '/src/clock.js';
import { subsolarPoint, sunAltAz } from '/domains/astronomy/js/index.js';

export const state = { box: null, network: null, district: null, stores: [], candidates: null, lane: null, run: null, building: null, active: new Set(),
  // The declared rates both analyses share (per hour): how much a store answers the other open (η) and the lane
  // (ξ); a permit granted, the lane built under it, recognised as it is used.
  rules: { eta: 0.5, xi: 0.8, permit: 0.3, build: 0.5, recognise: 0.4 } };
const listeners = new Set();
export function on(f) { listeners.add(f); return () => listeners.delete(f); }
export function emit(what) { for (const f of [...listeners]) try { f(what); } catch (e) { console.error(e); } }

const WALK = new Set(['primary', 'secondary', 'tertiary', 'residential', 'unclassified', 'living_street', 'pedestrian', 'service', 'primary_link', 'secondary_link', 'tertiary_link']);

// The Hoddle Grid's middle: Elizabeth Street to Russell Street, Flinders Street to Lonsdale Street.
export const DEFAULT_BOX = { south: -37.8173, west: 144.9617, north: -37.8117, east: 144.9697 };

export function boxOf(frame, b = DEFAULT_BOX) {
  const [x0, y0] = frame.toLocal(b.south, b.west), [x1, y1] = frame.toLocal(b.north, b.east);
  return { x0: Math.min(x0, x1), y0: Math.min(y0, y1), x1: Math.max(x0, x1), y1: Math.max(y0, y1), geo: b };
}
const inside = (p, B) => p[0] >= B.x0 && p[0] <= B.x1 && p[1] >= B.y0 && p[1] <= B.y1;

export function networkOf(model, B) {
  const runs = [];
  for (const r of model.roads || []) {
    if (!WALK.has(r.kind) || !r.line || r.line.length < 2) continue;
    let run = [];
    for (const p of r.line) {
      if (inside(p, B)) run.push(p);
      else { if (run.length >= 2) runs.push({ line: run, kind: r.kind, name: r.name || '' }); run = []; }
    }
    if (run.length >= 2) runs.push({ line: run, kind: r.kind, name: r.name || '' });
  }
  const key = (p) => `${Math.round(p[0] * 100)},${Math.round(p[1] * 100)}`;
  const uses = new Map();
  for (const r of runs) for (const p of new Set(r.line.map(key))) uses.set(p, (uses.get(p) || 0) + 1);
  const junctions = [], index = new Map();
  const junction = (p) => { const k = key(p); if (!index.has(k)) { index.set(k, junctions.length); junctions.push([+p[0].toFixed(2), +p[1].toFixed(2)]); } return index.get(k); };
  const segments = [], seen = new Set();
  const push = (line, kind, name) => {
    if (line.length < 2) return;
    let u = junction(line[0]), v = junction(line[line.length - 1]);
    if (u === v) {
      // A closed piece: split it at its middle point so that no segment is a loop.
      if (line.length < 3) return;
      const m = Math.floor(line.length / 2);
      push(line.slice(0, m + 1), kind, name);
      push(line.slice(m), kind, name);
      return;
    }
    const sig = `${Math.min(u, v)}-${Math.max(u, v)}-${line.length}-${key(line[Math.floor(line.length / 2)])}`;
    if (seen.has(sig)) return;
    seen.add(sig);
    const kindOf = kind === 'service' ? 'laneway' : kind === 'pedestrian' || kind === 'living_street' ? 'mall' : 'street';
    segments.push({ u, v, line: line.map((p) => [+p[0].toFixed(2), +p[1].toFixed(2)]), kind: kindOf, name });
  };
  for (const r of runs) {
    let start = 0;
    for (let k = 1; k < r.line.length; k++) {
      const last = k === r.line.length - 1;
      if (last || uses.get(key(r.line[k])) > 1) { push(r.line.slice(start, k + 1), r.kind, r.name); start = k; }
    }
  }
  // Drop junctions no segment uses (renumbering).
  const used = new Set(segments.flatMap((s) => [s.u, s.v]));
  const renumber = new Map();
  const J = [];
  junctions.forEach((p, i) => { if (used.has(i)) { renumber.set(i, J.length); J.push(p); } });
  for (const s of segments) s.u = renumber.get(s.u), s.v = renumber.get(s.v);
  return { junctions: J, segments };
}

// The district request: the network and the two stores, planarised.
export function request(extra = {}) {
  return { ...state.network, stores: state.stores.map((s) => ({ name: s.name, use: s.use, at: s.xy })), planarise: true, ...extra };
}

// Another plug-in's engine, when it is built (Street Oceanography's holdings and profiles); null otherwise.
export async function askOther(plugin, op, body = {}) {
  try {
    const r = await fetch(`/api/${plugin}/${op}`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
    return r.ok ? r.json() : null;
  } catch { return null; }
}

// The holdings inside the box, from Street Oceanography's streets.
export async function holdingsIn(ctx, B) {
  const streets = await askOther('street-oceanography', 'streets');
  if (!streets) return [];
  const mine = streets.filter((s) => s.area === ctx.area.id || (!s.area && s.city && ctx.area.name.includes(s.city)));
  const details = await Promise.all(mine.map((s) => askOther('street-oceanography', 'street', { street: s.id })));
  const out = [];
  for (const d of details) {
    if (!d) continue;
    for (const x of d.holdings)
      if (Number.isFinite(x.lat) && Number.isFinite(x.lon)) {
        const xy = ctx.frame.toLocal(x.lat, x.lon).slice(0, 2);
        if (inside(xy, B)) out.push({ id: x.id, name: x.name || '(unnamed)', use: (d.uses || [])[x.use] || '', street: d.name, xy });
      }
  }
  return out;
}

// The Sun at each civil hour of the clock's date: altitude and azimuth (radians, azimuth from north through east).
export function sunByHour(clock) {
  const alt = [], az = [];
  for (let h = 0; h < 24; h++) {
    const t = instantOf(clock.date, h + 0.5, clock.area.timezone);
    const s = sunAltAz(clock.area.lon, clock.area.lat, subsolarPoint(t));
    alt.push((s.altitude * Math.PI) / 180), az.push((s.azimuth * Math.PI) / 180);
  }
  return { altitude: alt, azimuth: az };
}

// Walkers entering the district each hour of the clock's weekday: half the counts at the sensors in the box.
export async function entriesByHour(ctx, B) {
  const list = await askOther('street-oceanography', 'profiles', { era: 'recent' });
  const wd = ctx.clock.weekday;
  if (list) {
    const inBox = list.filter((s) => inside(ctx.frame.toLocal(s.lat, s.lon).slice(0, 2), B));
    if (inBox.length) {
      const e = Array.from({ length: 24 }, (_, h) => 0.5 * inBox.reduce((a, s) => a + (Number.isFinite(s.profile[24 * wd + h]) ? s.profile[24 * wd + h] : 0), 0));
      return { entries: e, source: `half the hourly counts at the ${inBox.length} sensors in the district on ${['Mondays', 'Tuesdays', 'Wednesdays', 'Thursdays', 'Fridays', 'Saturdays', 'Sundays'][wd]} (City of Melbourne, 2024–2026)` };
    }
  }
  return { entries: Array.from({ length: 24 }, (_, h) => 150 + 1800 * Math.exp(-0.5 * ((h - 13) / 3.5) ** 2)), source: 'a declared day\'s shape (the sensors\' records are not at hand)' };
}

// ---- painters: the same drawing for a layer and for a tool ------------------------------------------------------

export const css = (rgb) => `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;

// The district's cells: blocks as pale flats (one hue each), segments as ribbons by kind, junctions as small posts,
// the stores and a planned lane. Returns the ids drawn.
export function paintCells(ctx, ids, cells, { highlight = new Set(), lane = null, stores = [] } = {}) {
  const { gis, G } = ctx;
  const cat = G.categorical();
  const night = ctx.theme === 'night';
  const pale = (c, t = night ? 0.25 : 0.55) => (night ? c.map((v) => v * (1 - t) + 0.1 * t) : c.map((v) => v + (1 - v) * t));
  const hot = [0.79, 0.55, 0.42];
  gis.addFlats(ids.blocks, cells.blocks, { ring: (f) => f.ring, color: (f, k) => (highlight.has(f.id) ? pale(hot, 0.2) : pale(cat(k % 10))), z: 0.12, pick: (f) => `a block, a 2-cell (#${f.id}): ${f.edges} streets round it, ${Math.round(f.area).toLocaleString('en')} m²` });
  const kindColour = { street: night ? [0.55, 0.66, 0.76] : [0.36, 0.56, 0.73], laneway: night ? [0.62, 0.6, 0.52] : [0.6, 0.56, 0.47], mall: [0.55, 0.7, 0.62], lane: [0.79, 0.55, 0.42] };
  gis.addRibbons(ids.edges, cells.edges, { line: (f) => f.line, width: (f) => (highlight.has(f.id) ? 4.5 : f.type === 'laneway' ? 2.2 : 3.2), color: (f) => (highlight.has(f.id) ? hot : kindColour[f.type] || kindColour.street), z: 0.3, order: 3,
    pick: (f) => `${f.name || 'a segment'}, a 1-cell (#${f.id}): ${Math.round(f.length)} m, ${f.blocks === 2 ? 'between two blocks' : f.blocks === 1 ? 'on the edge of one block' : 'bounding no block'}` });
  gis.addPoints(ids.vertices, cells.vertices, { position: (f) => f.at, radius: (f) => (highlight.has(f.id) ? 2.6 : 1.6), height: (f) => (highlight.has(f.id) ? 9 : 4), color: (f) => (highlight.has(f.id) ? hot : f.type === 'portal' ? [0.79, 0.55, 0.42] : night ? [0.8, 0.84, 0.88] : [0.3, 0.38, 0.45]),
    pick: (f) => `a ${f.type}, a 0-cell (#${f.id}) on ${f.degree} segments` });
  if (stores.length) gis.addPoints(ids.stores, stores, { position: (s) => s.at || s.xy, radius: () => 3.4, height: () => 14, color: () => [0.79, 0.55, 0.42], shadow: true, pick: (s) => `${s.name}${s.use ? ` — ${s.use}` : ''}` });
  if (lane && !lane.refused) gis.addArcs(ids.lane, [lane], { from: (l) => l.from, to: (l) => l.to, height: () => 6, radius: () => 0.9, base: 1, color: () => [0.79, 0.55, 0.42], opacity: 0.85 });
  return Object.values(ids);
}

// A frame of the run on its stratum: block warmth, walkers as columns at the junctions, flow arrows, the stores and
// the lane.
export function frameAt(run, hour) {
  if (!run || !run.frames.length) return null;
  const k = Math.max(0, Math.min(run.frames.length - 1, Math.floor(hour / 0.25 + 1e-9)));
  return run.frames[k];
}

// The middle of a polyline, by length; the unit direction from its first point to its last, reversed for a
// negative flow.
export const mid = (line) => {
  let L = 0;
  for (let i = 1; i < line.length; i++) L += Math.hypot(line[i][0] - line[i - 1][0], line[i][1] - line[i - 1][1]);
  let s = L / 2;
  for (let i = 1; i < line.length; i++) {
    const d = Math.hypot(line[i][0] - line[i - 1][0], line[i][1] - line[i - 1][1]);
    if (s <= d) return [line[i - 1][0] + ((line[i][0] - line[i - 1][0]) * s) / (d || 1), line[i - 1][1] + ((line[i][1] - line[i - 1][1]) * s) / (d || 1)];
    s -= d;
  }
  return line[line.length - 1];
};
export const dir = (line, f) => {
  const a = line[0], b = line[line.length - 1], d = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
  const s = f < 0 ? -1 : 1;
  return [(s * (b[0] - a[0])) / d, (s * (b[1] - a[1])) / d];
};

// A holding a passer-by walks into (not a vacant lot, an office, a hotel's rooms, a depot or a school).
export const visited = (c) => !/^(vacant|offices-and-professional|accommodation|transport|wholesale-and-industry|education)$/.test(c.use || '');

// The district, built once: the network in the box, the holdings in it, the nearest pair of them a lane could
// join (or, without holdings, two declared stores facing one block), and the engine's cells.
let building = null;
export function ensureDistrict(ctx) {
  if (state.district) return Promise.resolve(state.district);
  if (building) return building;
  building = (async () => {
    if (!ctx.model) throw new Error('the city model is not built yet');
    state.box = boxOf(ctx.frame);
    state.network = networkOf(ctx.model, state.box);
    state.candidates = await holdingsIn(ctx, state.box);
    // The two stores: of the holdings a passer-by walks into, the nearest pair a lane of 15 m or more could join
    // (shorter ones only cross a corner of the pavement).
    const shops = state.candidates.filter(visited);
    if (shops.length >= 2 && !state.stores.length) {
      const found = await ctx.ask('pairs', { ...state.network, planarise: true, stores: shops.map((c) => ({ name: c.name, use: c.use, at: c.xy })), within: 160, least: 15, limit: 8 });
      const best = found.pairs[0];
      if (best) state.stores = [shops[best.a], shops[best.b]];
    }
    if (state.stores.length < 2) {
      // Two declared stores on two opposite sides of the block nearest the box's middle.
      const d0 = await ctx.ask('district', { ...state.network, planarise: true });
      const mid = [(state.box.x0 + state.box.x1) / 2, (state.box.y0 + state.box.y1) / 2];
      const blocks = d0.cells.blocks.slice().sort((a, b) => dist(centre(a.ring), mid) - dist(centre(b.ring), mid));
      const ring = blocks[0].ring, n = ring.length;
      const at = (i) => { const a = ring[i], b = ring[(i + 1) % n]; return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]; };
      state.stores = [{ name: 'a declared store', use: '', xy: at(0) }, { name: 'another declared store', use: '', xy: at(Math.floor(n / 2)) }];
    }
    state.district = await ctx.ask('district', request());
    emit('district');
    return state.district;
  })();
  building.catch(() => { building = null; });
  return building;
}
const centre = (ring) => [ring.reduce((a, p) => a + p[0], 0) / ring.length, ring.reduce((a, p) => a + p[1], 0) / ring.length];
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);

// New stores (two holdings picked on the map): the district rebuilt, the lane and the run forgotten.
export async function setStores(ctx, stores) {
  state.stores = stores;
  state.lane = null, state.run = null;
  state.district = await ctx.ask('district', request());
  emit('district');
  emit('run');
  return state.district;
}

// A day of the district run as one process, forced by the clock's date: the Sun and the walkers entering.
export async function runDay(ctx, params = {}) {
  await ensureDistrict(ctx);
  const sun = sunByHour(ctx.clock);
  const { entries, source } = await entriesByHour(ctx, state.box);
  const run = await ctx.ask('run', request({ params: { seed: 3, ...params }, entries, altitude: sun.altitude }));
  state.run = { ...run, date: ctx.clock.date, entriesSource: source, entries, altitude: sun.altitude, params };
  emit('run');
  return state.run;
}

// ---- a building -------------------------------------------------------------------------------------------------

// The smallest rectangle about a footprint: of the directions of its convex hull's edges, the one whose box has the
// least area (Freeman and Shapira 1975: the minimal rectangle has a side on the hull). Returns its corner, heading
// (radians from east), width along the heading and depth across it.
export function rectangleOf(points) {
  const P = points.map((p) => [p[0], p[1]]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lower = [], upper = [];
  for (const p of P) { while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) lower.pop(); lower.push(p); }
  for (const p of P.slice().reverse()) { while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) upper.pop(); upper.push(p); }
  const hull = lower.slice(0, -1).concat(upper.slice(0, -1));
  let best = null;
  for (let i = 0; i < hull.length; i++) {
    const a = hull[i], b = hull[(i + 1) % hull.length], t = Math.atan2(b[1] - a[1], b[0] - a[0]);
    const c = Math.cos(t), s = Math.sin(t);
    let u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
    for (const p of hull) { const u = p[0] * c + p[1] * s, v = -p[0] * s + p[1] * c; u0 = Math.min(u0, u), u1 = Math.max(u1, u), v0 = Math.min(v0, v), v1 = Math.max(v1, v); }
    const area = (u1 - u0) * (v1 - v0);
    if (!best || area < best.area) best = { area, heading: t, width: u1 - u0, depth: v1 - v0, x: u0 * c - v0 * s, y: u0 * s + v0 * c };
  }
  return best;
}

// The building's box of rooms: about 7 m rooms, its floors from its levels or height, at most 160 rooms.
export function specOf(b) {
  const r = rectangleOf(b.footprint);
  const floors0 = Math.max(1, Math.min(8, Math.round(b.levels || (b.height_m || 9) / 3.5)));
  const nx = Math.max(1, Math.min(6, Math.round(r.width / 7))), ny = Math.max(1, Math.min(6, Math.round(r.depth / 7)));
  const floors = Math.max(1, Math.min(floors0, Math.floor(160 / (nx * ny))));
  return { x: r.x, y: r.y, heading: r.heading, width: r.width, depth: r.depth, storey: Math.max(2.8, Math.min(4.5, (b.height_m || 3.5 * floors0) / floors0)), nx, ny, floors, name: b.name || 'a building', id: b.id };
}

export async function runBuilding(ctx, spec, events) {
  const sun = sunByHour(ctx.clock);
  const result = await ctx.ask('building', { ...spec, altitude: sun.altitude, azimuth: sun.azimuth, events });
  state.building = { ...result, name: spec.name, id: spec.id, date: ctx.clock.date, requested: events };
  emit('building');
  return state.building;
}

// The world corners of a unit cube of the building's grid (x, y at its lower corner), as a footprint ring.
export function unitRing(spec, [i, j]) {
  const c = Math.cos(spec.heading), s = Math.sin(spec.heading), dx = spec.width / spec.nx, dy = spec.depth / spec.ny;
  const at = (u, v) => [spec.x + u * c - v * s, spec.y + u * s + v * c];
  return [at(i * dx, j * dy), at((i + 1) * dx, j * dy), at((i + 1) * dx, (j + 1) * dy), at(i * dx, (j + 1) * dy)];
}
