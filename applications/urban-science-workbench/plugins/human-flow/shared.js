// What the plug-in's layers and analyses share: the district (the shell's: its walking network and buildings cut from
// the city model in Landscape Dynamics' box), the City's sensors placed in the area's local frame with the way each
// stands on, the engine's preparation of the walkers, and the run the map draws, with listeners, so that the map and
// the drawer show one district's walkers. A module is loaded once by its URL, so every layer and tool holds this object.
//
// A run is named by its day and window of civil hours, a street's regime and an event's crowd (engine/runs.hpp). The
// map draws the run a tool chose while the clock is inside its window, and the day as built otherwise.

import { instantOf } from '/src/clock.js';
import { districtOf, boxOf, networkOf, footprintsIn } from '/src/district.js';

export const MODES = ['passing through', 'arriving at work', 'leaving work', 'between buildings', 'visiting shops', "an event's crowd"];
export const STATUS = ['about to leave', 'walking', 'waiting to enter a full street', 'waiting at a red signal', 'standing at the event', 'arrived'];
// Walkers by pattern: deep sea, sand, sea glass, lavender, gold and rose, every two apart in every vision (OKLab ΔE ≥ 0.09
// under the three deficiencies; checked in tests/run-tests.mjs).
export const MODE_COLOURS = ['#2E5C88', '#DA8C5E', '#96C8B8', '#8B6BC0', '#E9C24F', '#D45A86'];
// Walkers by what they are doing: walking sea blue, queued at a full street brick, held at red gold, standing at the event rose.
export const STATUS_COLOURS = ['#93A3AE', '#4C82B0', '#8E3B33', '#EDB45C', '#D45A86', '#93A3AE'];
export const rgb = (hex) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255);
export const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.round((x % 1) * 60) % 60).padStart(2, '0')}`;
export const DAY_NAMES = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'];

export const state = { gis: null, area: null, body: null, key: null, prep: null, scenario: null, colourBy: 'pattern', compare: null, record: null, picked: null };
const listeners = new Set();
export function on(f) { listeners.add(f); return () => listeners.delete(f); }
export function emit(what) { for (const f of [...listeners]) try { f(what); } catch (e) { console.error(e); } }

const WALK = new Set(['primary', 'secondary', 'tertiary', 'residential', 'unclassified', 'living_street', 'pedestrian', 'service', 'primary_link', 'secondary_link', 'tertiary_link']);
const wayKind = (k) => (k === 'service' ? 'laneway' : k === 'pedestrian' || k === 'living_street' ? 'mall' : 'street');

// The nearest way of the city model where people walk to a point: its name and kind (as the problem's prepare.py).
function nearestWay(roads, p) {
  let best = [Infinity, null, ''];
  for (const r of roads) for (let i = 0; i + 1 < r.line.length; i++) {
    const a = r.line[i], b = r.line[i + 1], ex = b[0] - a[0], ey = b[1] - a[1], L2 = ex * ex + ey * ey;
    const t = L2 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / L2)) : 0;
    const d = Math.hypot(p[0] - a[0] - t * ex, p[1] - a[1] - t * ey);
    if (d < best[0]) best = [d, r.kind, r.name || ''];
  }
  return { name: best[2], kind: wayKind(best[1]) };
}

// The district request: the box, the network and the buildings (the shell's), the sensors within 500 m of the box.
async function districtBody(ctx) {
  const B = boxOf(ctx.frame, districtOf(ctx.area));
  const net = networkOf(ctx.model, B);
  const buildings = footprintsIn(ctx.model, B).map((b) => ({ key: b.key, ring: b.ring, base: b.base, top: b.top }));
  const r = await fetch('/data/products/street-sensors.json');
  if (!r.ok) throw new Error('the City\'s sensors are not in the data library yet: the data lamp runs the recipe');
  const S = (await r.json()).sensors;
  const near = (p) => p[0] > B.x0 - 500 && p[0] < B.x1 + 500 && p[1] > B.y0 - 500 && p[1] < B.y1 + 500;
  const roads = (ctx.model.roads || []).filter((x) => WALK.has(x.kind) && x.line && x.line.length >= 2);
  const sensors = [];
  for (const s of S) {
    if (!Number.isFinite(s.lat) || !s.lat) continue;
    const p = ctx.frame.toLocal(s.lat, s.lon), at = [+p[0].toFixed(2), +p[1].toFixed(2)];
    if (near(at)) sensors.push({ id: s.id, at, way: nearestWay(roads, at) });
  }
  return { box: { x0: B.x0, y0: B.y0, x1: B.x1, y1: B.y1 }, ...net, buildings, sensors };
}

// The engine, asked of the prepared district by its key; sent whole when the engine does not know it (restarted).
export async function ask(ctx, op, body = {}) {
  await ensurePrepared(ctx);
  try {
    return await ctx.ask(op, { district: state.key, ...body });
  } catch (e) {
    if (!/unknown district/.test(String(e.message))) throw e;
    return ctx.ask(op, { ...state.body, ...body });
  }
}

let preparing = null;
export function ensurePrepared(ctx) {
  state.gis = ctx.gis;
  if (state.prep && state.area === ctx.area.id) return Promise.resolve(state.prep);
  if (preparing && state.area === ctx.area.id) return preparing;
  state.area = ctx.area.id, state.prep = null, state.key = null, state.scenario = null, state.compare = null, state.record = null;
  preparing = (async () => {
    if (!ctx.model) throw new Error('the city model is not built yet: the data lamp runs the recipe');
    state.body = await districtBody(ctx);
    const P = await ctx.ask('prepare', state.body);
    state.prep = P, state.key = P.key;
    emit('prepared');
    return P;
  })();
  preparing.catch(() => { preparing = null; });
  return preparing;
}

// The day as built on the clock's weekday: the whole day.
export const baseSpec = (clock) => ({ weekday: clock.weekday, from: 0, to: 24 });
// The run the map draws at the clock's hour: the chosen scenario inside its window, else the day as built.
export function specAt(clock) {
  const s = state.scenario;
  if (s && s.spec.weekday === clock.weekday && clock.hour >= s.spec.from && clock.hour < s.spec.to) return { spec: s.spec, label: s.label, scenario: true };
  return { spec: baseSpec(clock), label: `the ${DAY_NAMES[clock.weekday]} as built`, scenario: false };
}
export function setScenario(scenario) { state.scenario = scenario; emit('scenario'); }

// Starts a run (or finds it) and waits for it, reporting its progress; resolves with its result.
export async function runAndWait(ctx, spec, progress = () => {}, isAlive = () => true) {
  let st = await ask(ctx, 'run', spec);
  while (!st.done && isAlive()) {
    progress(st);
    await new Promise((r) => setTimeout(r, 1200));
    st = await ask(ctx, 'run', spec);
  }
  if (st.error) throw new Error(st.error);
  return ask(ctx, 'result', spec);
}

// The instant (MJD, UTC) a run's window begins on a civil date.
export const mjdOf = (date, hour, timezone) => instantOf(date, hour, timezone) / 86400000 + 40587;

export const streetOf = (name) => (state.prep ? state.prep.streets.find((s) => s.name === name) : null);
