// The district's spacetime (a model): the day of the landscape with time drawn as height above the district, from
// midnight at the ground, RISE metres an hour. Three kinds of emergence can be read in it.
//   The flow: at every hour the district's streets are a thin slice, each street coloured by the walkers it
//   carries; up the stack of slices the day's tide of walkers is a pattern of the whole that no street holds alone.
//   The structure: what has a history stands as a worldline or a sheet over the hours it exists. The two stores
//   rise as columns, thick while open and thin while shut; where the lane is built its sheet begins in mid-air, with
//   the worldlines of the two junctions it made, and the blocks are outlined at that instant, the new cells in sand.
//   The meaning: the lane's sheet deepens from the instant it is recognised. The cells are the same; what may happen
//   next is not, for a recognised lane can no longer be removed.
// Events are beads at their place and hour, each joined by a thread to the next events that read what it wrote:
// the record's causal order, its past all that can be reached down the threads. Faint frames mark every six hours,
// and a pale plane the clock's hour; beneath it the day as it went, above it, fainter, what the model goes on to do.
// The run is the engine's (urban-science, landscape/process.hpp), the same as the layer "the landscape running".

import { ensureDistrict, on, runDay, state } from '../shared.js';

export const RISE = 10;   // metres of height for an hour
const z = (t) => RISE * t;

// The events' kinds, in the order of the legend, with their colours (sRGB).
const KINDS = [
  { mechanism: 'a store opens', name: 'a store opens', day: [0.25, 0.48, 0.66], night: [0.55, 0.72, 0.86] },
  { mechanism: 'a store closes', name: 'a store shuts', day: [0.58, 0.64, 0.68], night: [0.66, 0.72, 0.76] },
  { mechanism: 'the permit is granted', name: 'a rule in force: the permit', day: [0.35, 0.44, 0.66], night: [0.62, 0.69, 0.89] },
  { mechanism: 'the lane is built', name: 'the rule carried out: the lane built', day: [0.79, 0.55, 0.42], night: [0.89, 0.67, 0.54] },
  { mechanism: 'the lane is removed', name: 'the lane removed', day: [0.65, 0.42, 0.31], night: [0.94, 0.76, 0.66] },
  { mechanism: 'the lane is recognised', name: 'a meaning adopted: the lane recognised', day: [0.79, 0.6, 0.29], night: [0.89, 0.74, 0.43] },
  { mechanism: 'the stores close for the night', name: 'the night (at the hour)', day: [0.5, 0.55, 0.6], night: [0.62, 0.68, 0.72] },
];

export async function layer(ctx) {
  const { dom, gis, G } = ctx;
  const { h, fmt } = dom;
  const night = ctx.theme === 'night';
  const ids = {
    slices: ctx.layerId(':slices'), lane: ctx.layerId(':lane'), lines: ctx.layerId(':worldlines'), outlines: ctx.layerId(':outlines'), frames: ctx.layerId(':frames'),
    events: ctx.layerId(':events'), threads: ctx.layerId(':threads'), now: ctx.layerId(':now'), nowEdge: ctx.layerId(':now-edge'),
  };
  state.active.add('spacetime');
  const status = h('div', { class: 'faint' }, 'running the day…');
  const nowLine = h('div', {}), rampBox = h('div', {});
  const kindOf = (m) => KINDS.find((k) => k.mechanism === m) || KINDS[6];
  const colourOf = (m) => kindOf(m)[night ? 'night' : 'day'];
  const css = (c) => `rgb(${c.map((v) => Math.round(v * 255)).join(',')})`;
  const ink = night ? [0.78, 0.84, 0.88] : [0.3, 0.38, 0.45], faint = night ? [0.42, 0.5, 0.56] : [0.68, 0.74, 0.78];
  const laneColour = night ? [0.89, 0.67, 0.54] : [0.79, 0.55, 0.42], laneDeep = night ? [0.94, 0.5, 0.38] : [0.6, 0.31, 0.2];
  const open = night ? [0.55, 0.72, 0.86] : [0.25, 0.48, 0.66];
  const legend = h('div', { class: 'card' }, h('b', {}, 'The district\'s spacetime (a model)'),
    h('div', { class: 'faint' }, `time rises ${RISE} m an hour from midnight at the ground: a slice of the streets each hour, faint frames every six`),
    rampBox,
    h('div', { class: 'legend', style: 'flex-direction:column;gap:2px;margin-top:4px' },
      h('span', {}, h('i', { style: `background:${css(open)};height:7px` }), 'a store\'s worldline, thick while open'),
      h('span', {}, h('i', { style: `background:${css(laneColour)};height:7px` }), 'the lane\'s sheet, deep once recognised'),
      KINDS.slice(0, 6).map((k) => h('span', {}, h('i', { class: 'dot', style: `background:${css(k[night ? 'night' : 'day'])}` }), k.name))),
    h('div', { class: 'faint' }, 'a thread joins an event to the next ones that read what it wrote'),
    nowLine, status);

  const paper = night ? [0.08, 0.11, 0.15] : [0.95, 0.965, 0.973];
  const fade = (c, past) => (past ? c : c.map((v, i) => v + (paper[i] - v) * 0.66));
  let hour = ctx.clock.hour, built = null;

  // Everything that does not move with the clock, built once for a run.
  function build() {
    const R = state.run;
    if (!R) return;
    const configs = R.strata;
    const configAt = (t) => configs.findIndex((c) => c.intervals.some(([a, b]) => t >= a - 1e-9 && t < b - 1e-9));
    const frameAt = (t) => R.frames[Math.max(0, Math.min(R.frames.length - 1, Math.round(t / 0.25)))];
    // The slices: the streets at each whole hour, coloured by the walkers they carry.
    const traffic = (F, id) => (F.forth[id] || 0) + (F.back[id] || 0);
    let top = 0;
    for (const f of R.frames) for (const id of Object.keys(f.flow)) top = Math.max(top, traffic(f, id));
    const scale = Math.max(1, top);
    const flowColour = G.sequential(0, Math.log1p(scale), 'sea', { from: night ? 0.14 : 0.16, to: 1 });
    rampBox.replaceChildren(ctx.chart.ramp('sea', 0, Math.round(scale), 'walkers an hour along a street, both ways (logarithmic)'));
    const slices = [];
    for (let t = 0; t < 24; t++) {
      const F = frameAt(t), C = configs[F.stratum];
      for (const e of C.cells.edges) slices.push({ e, t, flow: traffic(F, e.id), net: F.flow[e.id] || 0, lane: e.id === F.lane, recognised: e.id === F.lane && F.recognised });
    }
    gis.addSlabs(ids.slices, slices, {
      line: (s) => s.e.line, width: (s) => (s.lane ? 3.2 : s.e.type === 'laneway' ? 1.3 : 2.3), base: (s) => z(s.t), height: () => 0.45,
      color: (s) => fade(s.lane ? (s.recognised ? laneDeep : laneColour) : flowColour(Math.log1p(s.flow)), s.t <= hour), opacity: 0.88, order: 4,
      pick: (s) => `${s.e.name || (s.lane ? 'the lane' : 'a way')} at ${hhmm(s.t)}: ${fmt(s.flow, 3)} walkers an hour, ${fmt(Math.abs(s.net), 2)} more one way than the other${s.lane ? (s.recognised ? ' — the lane, recognised' : ' — the lane') : ''}`,
    });
    // The lane's sheet over the quarter hours it stands, cut where the cells change.
    const sheet = [];
    R.frames.forEach((F, k) => {
      if (!F.lane) return;
      const C = configs[F.stratum];
      const end = Math.min(F.time + 0.25, (C.intervals.find(([a, b]) => F.time >= a - 1e-9 && F.time < b - 1e-9) || [0, F.time + 0.25])[1]);
      const e = C.cells.edges.find((x) => x.id === F.lane);
      if (e) sheet.push({ e, a: F.time, b: end, recognised: F.recognised });
      // The quarter hour in which it is built, from the instant it is.
      const prev = R.frames[k - 1];
      if (prev && !prev.lane) {
        const from = C.intervals.find(([a]) => a <= F.time + 1e-9 && a > prev.time - 1e-9);
        if (from && e) sheet.push({ e, a: from[0], b: F.time, recognised: false });
      }
    });
    gis.addSlabs(ids.lane, sheet, { line: (s) => s.e.line, width: () => 2.2, base: (s) => z(s.a), height: (s) => Math.max(0.05, z(s.b) - z(s.a)),
      color: (s) => fade(s.recognised ? laneDeep : laneColour, s.a < hour), opacity: 0.82, order: 5,
      pick: (s) => `the lane, ${hhmm(s.a)}–${hhmm(s.b)}${s.recognised ? ': recognised' : ': built, not yet recognised'}` });
    // Worldlines: the stores, thick while open; the lane's two junctions, while they exist.
    const stores = configs[0].stores;
    const lines = [];
    R.frames.forEach((F) => {
      stores.forEach((s, i) => lines.push({ at: s.at, a: F.time, b: F.time + 0.25, open: i === 0 ? F.a : F.b, store: s }));
      if (F.lane) {
        const C = configs[F.stratum], e = C.cells.edges.find((x) => x.id === F.lane);
        if (e) for (const id of [e.tail, e.head]) { const v = C.cells.vertices.find((x) => x.id === id); if (v) lines.push({ at: v.at, a: F.time, b: F.time + 0.25, portal: v }); }
      }
    });
    gis.addPoints(ids.lines, lines, { position: (l) => l.at, radius: (l) => (l.portal ? 1.1 : l.open ? 3.2 : 0.9), base: (l) => z(l.a), height: (l) => z(l.b) - z(l.a),
      color: (l) => fade(l.portal ? laneColour : l.open ? open : faint, l.a < hour),
      pick: (l) => (l.portal ? `a junction the lane made (#${l.portal.id})` : `${l.store.name} at ${hhmm(l.a)}: ${l.open ? 'open' : 'shut'}`) });
    // The blocks outlined at each instant the cells change, the cells it made in sand.
    const outlines = [];
    for (const e of R.events) {
      if (!e.rewrites || !e.rewrites.length) continue;
      const c = configs[configAt(e.time + 1e-6)];
      if (!c) continue;
      const made = new Set(e.rewrites.flatMap((s) => s.created));
      for (const b of c.cells.blocks) if (made.has(b.id)) outlines.push({ ring: [...b.ring, b.ring[0]], t: e.time, made: true, event: e });
      for (const b of configs[configAt(e.time - 1e-6)]?.cells.blocks || []) if (e.rewrites.some((s) => s.deleted.includes(b.id))) outlines.push({ ring: [...b.ring, b.ring[0]], t: e.time, made: false, event: e });
    }
    gis.addTubes(ids.outlines, outlines, { path: (o) => o.ring.map(([x, y]) => [x, y, z(o.t) + 0.4]), radius: (o) => (o.made ? 0.55 : 0.3),
      color: (o) => (o.made ? laneColour : ink), opacity: 0.9, pick: (o) => `${hhmm(o.t)}: ${o.event.mechanism} — ${o.made ? 'a block it made' : 'the block it replaced'}` });
    // The frames of six hours and the box's four corner posts.
    const B = state.box;
    const rect = [[B.x0, B.y0], [B.x1, B.y0], [B.x1, B.y1], [B.x0, B.y1]];
    const frames = [6, 12, 18, 24].map((t) => ({ path: [...rect, rect[0]].map(([x, y]) => [x, y, z(t)]), t }));
    for (const [x, y] of rect) frames.push({ path: [[x, y, 0], [x, y, z(24)]], post: true });
    gis.addTubes(ids.frames, frames, { path: (f) => f.path, radius: (f) => (f.post ? 0.6 : 0.4), color: () => faint, opacity: 0.6, pick: (f) => (f.post ? 'a day, midnight to midnight' : `${hhmm(f.t % 24)}${f.t === 24 ? ' (the next midnight)' : ''}`) });
    // Events as beads, threads for the causal order.
    const places = eventPlaces(R);
    const beads = R.events.map((e, i) => ({ e, i, at: places[i] }));
    gis.addPoints(ids.events, beads, { position: (b) => b.at, radius: (b) => (b.e.rewrites && b.e.rewrites.length ? 7 : 5), base: (b) => z(b.e.time) - 1.2, height: () => 2.4,
      color: (b) => fade(colourOf(b.e.mechanism), b.e.time <= hour),
      pick: (b) => `${hhmm(b.e.time)} ${b.e.mechanism}${b.e.match && !/lane|permit|both/.test(b.e.match) ? `: ${b.e.match}` : ''}` });
    const threads = (R.order.covers || []).map(([a, b]) => ({ a: beads[a], b: beads[b] })).filter((t) => t.a && t.b);
    gis.addTubes(ids.threads, threads, { path: (t) => [[...t.a.at, z(t.a.e.time)], [...t.b.at, z(t.b.e.time)]], radius: () => 0.9,
      color: (t) => fade(colourOf(t.a.e.mechanism), t.b.e.time <= hour), opacity: 0.7,
      pick: (t) => `${hhmm(t.a.e.time)} ${t.a.e.mechanism} → ${hhmm(t.b.e.time)} ${t.b.e.mechanism}: the second read what the first wrote` });
    // The clock's plane, built at height zero and lifted to the hour.
    gis.addFlats(ids.now, [{ ring: rect }], { ring: (f) => f.ring, z: 0, color: () => (night ? [0.36, 0.5, 0.62] : [0.62, 0.76, 0.86]), opacity: 0.16, order: 6, pick: false });
    gis.addTubes(ids.nowEdge, [{ path: [...rect, rect[0]].map(([x, y]) => [x, y, 0]) }], { path: (f) => f.path, radius: () => 0.8, color: () => open, pick: false });
    built = { beads, threads, past: closure(R) };
    lift();
  }

  // What follows the clock: the plane, and the past drawn full and the future faint.
  function lift() {
    for (const id of [ids.now, ids.nowEdge]) { const l = gis.layers.get(id); if (l) l.object.position.y = z(hour); }
    const R = state.run;
    if (!R) return;
    const F = R.frames[Math.max(0, Math.min(R.frames.length - 1, Math.floor(hour / 0.25 + 1e-9)))];
    const done = R.events.filter((e) => e.time <= hour).length;
    const S = R.strata[F.stratum];
    nowLine.replaceChildren(
      h('div', { style: 'margin-top:4px' }, `${hhmm(hour)}, at the plane: ${done} of the day's ${R.events.length} events have happened`),
      h('div', { class: 'faint' }, `${S.counts[0]} junctions, ${S.counts[1]} ways, ${S.counts[2]} blocks${F.lane ? `; the lane ${F.recognised ? 'recognised' : 'standing'}` : F.permit ? '; a permit, no lane yet' : ''}; can happen now: ${(F.enabled || []).map((k) => kindOf(k).name.replace(/^.*: /, '')).join(', ') || 'nothing'}`));
    for (const id of [ids.slices, ids.lane, ids.lines, ids.events, ids.threads]) gis.restyle(id, {});
  }

  // The view: from the side, a little above, at the middle of the day, from the direction the map was seen.
  function look() {
    const B = state.box, T = gis.toWorld((B.x0 + B.x1) / 2, (B.y0 + B.y1) / 2, z(10));
    const o = gis.camera.position.clone().sub(gis.controls.target);
    const heading = Math.atan2(o.x, o.z), tilt = (64 * Math.PI) / 180, d = 1.45 * Math.max(B.x1 - B.x0, B.y1 - B.y0, z(24));
    gis.controls.target.copy(T);
    gis.camera.position.set(T.x + d * Math.sin(tilt) * Math.sin(heading), T.y + d * Math.cos(tilt), T.z + d * Math.sin(tilt) * Math.cos(heading));
    gis.controls.update();
  }
  const off = on((what) => { if (what === 'run') { build(); status.textContent = summary(); } });
  const summary = () => {
    const R = state.run;
    const changes = R.events.filter((e) => e.rewrites && e.rewrites.length).length;
    return `${R.events.length} events through ${R.date}; the cells changed ${changes} time${changes === 1 ? '' : 's'}, ${R.strata.length} configuration${R.strata.length === 1 ? '' : 's'}`;
  };
  try {
    await ensureDistrict(ctx);
    if (!state.run) await runDay(ctx);
    build();
    status.textContent = summary();
    look();
  } catch (e) {
    status.textContent = String(e.message || e);
  }
  return {
    legend,
    setTime(c) { hour = c.hour; if (built) lift(); },
    describe(hit) {
      if (hit.layer !== ids.events || !built) return null;
      const b = hit.feature, past = built.past;
      const before = past[b.i], after = past.map((p, j) => (p.has(b.i) ? j : -1)).filter((j) => j >= 0);
      const fp = b.e.footprint || { reads: [], writes: [] };
      const plain = (k) => ({ 'a:A': 'the first store\'s state', 'a:B': 'the second store\'s state', permit: 'the permit', lane: 'the lane', recognised: 'its recognition', day: 'the hour of the day' }[k] || `cell #${k}`);
      const named = (list) => [...new Set(list.map(plain).filter((x) => !x.startsWith('cell')))].join(', ') || '—';
      const cells = fp.writes.filter((k) => /^\d/.test(k)).length;
      return {
        kind: 'event', title: `${hhmm(b.e.time)} — ${b.e.mechanism}`,
        lines: [
          b.e.match || '',
          `read ${named(fp.reads)}; wrote ${named(fp.writes)}${cells ? `, and ${cells} cells and incidences` : ''}`,
          `its past: ${before.size} event${before.size === 1 ? '' : 's'} of the day; its future: ${after.length}`,
          b.e.rewrites && b.e.rewrites.length ? `${b.e.rewrites.length} rewrites of the cells: ${b.e.rewrites.map((s) => `${s.kind} (${s.refinement ? 'the same place, finer' : 'a new place'})`).join(', ')}` : '',
        ].filter(Boolean),
      };
    },
    dispose() { off(); state.active.delete('spacetime'); },
  };
}

// Where each event stands: a store's at the store, the permit's, the lane's and its recognition's at the lane's
// middle (between the stores when no lane can be built), the night's between the stores.
function eventPlaces(R) {
  const D = state.district;
  const stores = R.strata[0].stores;
  const mid = (a, b) => [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
  const lane = D && D.lane && !D.lane.refused ? mid(D.lane.from, D.lane.to) : mid(stores[0].at, stores[1].at);
  return R.events.map((e) => {
    if (/store/.test(e.mechanism) && !/night/.test(e.mechanism)) {
      const s = stores.find((x) => x.name === e.match);
      if (s) return s.at;
    }
    if (/night/.test(e.mechanism)) return mid(stores[0].at, stores[1].at);
    return lane;
  });
}

// The causal past of each event: everything below it along the covers.
function closure(R) {
  const n = R.events.length, below = Array.from({ length: n }, () => []);
  for (const [a, b] of R.order.covers || []) below[b].push(a);
  const past = [];
  for (let i = 0; i < n; i++) {
    const seen = new Set(), stack = [...below[i]];
    while (stack.length) { const j = stack.pop(); if (!seen.has(j)) { seen.add(j); stack.push(...below[j]); } }
    past.push(seen);
  }
  return past;
}

const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.floor((x % 1) * 60 + 1e-6)).padStart(2, '0')}`;
