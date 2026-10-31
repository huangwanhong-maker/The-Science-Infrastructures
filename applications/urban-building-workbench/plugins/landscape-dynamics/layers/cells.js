// The district as a complex: the walking network of the Hoddle Grid's middle as a regular cell complex, drawn as a
// complex. On the ground: each block (a 2-cell) a tile drawn a little inside its streets, so that the cells stand
// apart; each segment (a 1-cell) a ribbon, darker the more blocks it bounds (none, one, two), with an arrow from its
// tail to its head, its orientation; each junction (a 0-cell) a post; the holes (faces with a lane ending inside,
// which no single cell can be) outlined and faintly filled. By dimension: the junctions on the ground, the segments
// lifted above them and the blocks above those, each cell joined by a thread to each of its faces, sea where the
// incidence is +1 and sand where it is −1 — the complex's face poset with its signs, standing over the map. A cell
// picked shows its closure (it and its faces, sand) and its star (it and what it bounds, sea). The state set on the
// cells for a run stands on them: walkers as posts, a block's warmth as its tint. The edits made to the complex (the
// analysis "The district: complex, freedom, dynamics") are drawn as they are made.
// The buildings are cells too: each a 3-cell standing on the ground at its height, translucent, the walls two
// buildings share drawn in sand (2-cells with two cofaces); by dimension they float above the blocks as a fourth
// level, each joined to the neighbours it shares a wall with. After a run they are coloured by their warmth at the
// clock's hour.

import { airByHour, around, emit, ensureDistrict, ensureMassing, frameAt, on, state } from '../shared.js';

export const LIFT = { 0: 0, 1: 32, 2: 64, 3: 104 };   // metres above the ground of each dimension, drawn by dimension
let hourNow = 12;

export async function layer(ctx) {
  const { dom } = ctx;
  const { h } = dom;
  const ids = complexIds(ctx.layerId);
  const status = h('div', { class: 'faint' }, 'cutting the district from the city model…');
  const keys = h('div', { class: 'keys tight', style: 'margin:4px 0' });
  const legend = h('div', { class: 'card' }, h('b', {}, 'The district as a complex'), keys, status,
    h('div', { class: 'faint' }, 'buildings (3-cells): a soft hue for each terrace, the walls they share in sand; after a run, sea where cooler than the air, brick where warmer'));
  state.active.add('cells');
  const drawKeys = () => keys.replaceChildren(...[['ground', 'on the ground'], ['dimension', 'by dimension']].map(([v, name]) =>
    h('button', { class: `key${state.view === v ? ' on' : ''}`, onclick: () => { state.view = v; emit('view'); } }, name)));
  const draw = () => {
    if (!state.district) return;
    paintComplex(ctx, ids, state.district);
    status.replaceChildren(...summary(h, state.district));
    drawKeys();
  };
  const off = on((what) => { if (['district', 'picked', 'freedom', 'view', 'massing', 'run'].includes(what)) draw(); });
  try {
    await ensureDistrict(ctx);
    draw();
    await ensureMassing(ctx);
  } catch (e) {
    status.textContent = String(e.message || e);
  }
  return {
    legend,
    setTime(c) { recolourBuildings(ctx, ids, c.hour); },
    describe: (hit) => describeCell(hit, state.district),
    dispose() { off(); state.active.delete('cells'); },
  };
}

export function complexIds(layerId) {
  return { blocks: layerId(':blocks'), edges: layerId(':edges'), arrows: layerId(':arrows'), vertices: layerId(':vertices'), holes: layerId(':holes'),
    holeEdges: layerId(':hole-edges'), incidence: layerId(':incidence'), stores: layerId(':stores'), buildings: layerId(':buildings'), party: layerId(':party'),
    dual: layerId(':dual') };
}

// The counts, homology and holes of a complex, as lines.
export function summary(h, c) {
  const t = c.topology, [V, E, F] = t.counts, M = state.massing;
  return [
    h('div', {}, `${V} junctions (0-cells), ${E} segments (1-cells), ${F} blocks (2-cells)${M ? `; ${M.buildings.length} buildings (3-cells) sharing ${M.walls.length} walls, in ${M.terraces} terraces` : ''}`),
    h('div', {}, `b₀ = ${t.betti[0]}, b₁ = ${t.betti[1]}, χ = ${t.euler}: ${t.betti[0] === 1 ? 'one piece' : `${t.betti[0]} pieces`}, ${t.holes} hole${t.holes === 1 ? '' : 's'}`),
    h('div', { class: 'faint' }, `${c.edits && c.edits.length ? `${c.edits.length} edit${c.edits.length === 1 ? '' : 's'} made; ` : ''}regular: ${t.regular ? 'yes' : t.reason}`),
  ];
}

const inset = (ring, k) => {
  const n = ring.length, cx = ring.reduce((a, p) => a + p[0], 0) / n, cy = ring.reduce((a, p) => a + p[1], 0) / n;
  return ring.map(([x, y]) => [cx + k * (x - cx), cy + k * (y - cy)]);
};
const midOf = (line) => {
  let L = 0;
  for (let i = 1; i < line.length; i++) L += Math.hypot(line[i][0] - line[i - 1][0], line[i][1] - line[i - 1][1]);
  let s = L / 2;
  for (let i = 1; i < line.length; i++) {
    const d = Math.hypot(line[i][0] - line[i - 1][0], line[i][1] - line[i - 1][1]);
    if (s <= d) { const u = s / (d || 1); return [line[i - 1][0] + u * (line[i][0] - line[i - 1][0]), line[i - 1][1] + u * (line[i][1] - line[i - 1][1])]; }
    s -= d;
  }
  return line[line.length - 1];
};
export const centroid = (ring) => [ring.reduce((a, p) => a + p[0], 0) / ring.length, ring.reduce((a, p) => a + p[1], 0) / ring.length];

// Paints a complex on the map (both views), with the picked cell's closure and star and the state set on cells.
export function paintComplex(ctx, ids, c) {
  const { gis, G } = ctx;
  for (const id of Object.values(ids)) gis.remove(id);
  const night = ctx.theme === 'night';
  const lifted = state.view === 'dimension';
  const z = (d) => (lifted ? LIFT[d] : 0);
  const V = c.cells.vertices.map((x) => ({ ...x, dim: 0 })), E = c.cells.edges.map((x) => ({ ...x, dim: 1 })), B = c.cells.blocks.map((x) => ({ ...x, dim: 2 }));
  const byId = new Map([...V, ...E, ...B].map((x) => [x.id, x]));
  const picked = state.picked && (state.picked.dim === 3 ? state.picked : byId.get(state.picked.id));
  const near = picked ? around(c, picked) : null;
  const role = (id) => (!near ? 'none' : id === (picked.dim === 3 ? picked.key : picked.id) ? 'picked' : near.closure.has(id) ? 'closure' : near.star.has(id) ? 'star' : 'other');
  const sand = night ? [0.89, 0.67, 0.54] : [0.79, 0.55, 0.42], sea = night ? [0.55, 0.72, 0.86] : [0.25, 0.48, 0.66], deep = night ? [0.94, 0.5, 0.38] : [0.6, 0.31, 0.2];
  const paper = night ? [0.08, 0.11, 0.15] : [0.95, 0.965, 0.973];
  const mix = (a, b, t) => a.map((v, i) => v + (b[i] - v) * t);
  // By dimension the rest fades, so that the picked cell's neighbourhood in the face poset stands out; on the ground it stays.
  const tone = (base, id) => { const r = role(id); return r === 'picked' ? deep : r === 'closure' ? sand : r === 'star' ? sea : r === 'other' && lifted ? mix(base, paper, 0.6) : base; };
  const warmth = G.sequential(10, 34, 'gold', { from: 0.08, to: 1 });
  const blockBase = night ? [0.3, 0.44, 0.55] : [0.6, 0.74, 0.84];
  // 2-cells.
  gis.addFlats(ids.blocks, B, { ring: (f) => inset(f.ring, lifted ? 0.86 : 0.94), z: lifted ? LIFT[2] : 0.14, opacity: lifted ? 0.82 : 0.9, order: 1,
    color: (f) => tone(state.initial.temperature[f.id] !== undefined ? warmth(state.initial.temperature[f.id]) : blockBase, f.id),
    pick: (f) => `a block, a 2-cell (#${f.id}): bounded by ${f.boundary.length} segments, ${Math.round(f.area).toLocaleString('en')} m²${state.initial.temperature[f.id] !== undefined ? `, set to ${state.initial.temperature[f.id]} °C` : ''}` });
  // 1-cells, darker the more blocks they bound.
  const edgeBase = (e) => (e.cofaces.length === 2 ? (night ? [0.75, 0.83, 0.9] : [0.2, 0.32, 0.43]) : e.cofaces.length === 1 ? (night ? [0.55, 0.64, 0.72] : [0.42, 0.53, 0.62]) : (night ? [0.38, 0.45, 0.52] : [0.66, 0.72, 0.77]));
  gis.addRibbons(ids.edges, E, { line: (f) => f.line, z: lifted ? LIFT[1] : 0.32, order: 3, width: (f) => (role(f.id) === 'picked' ? 4.5 : f.type === 'lane' ? 3.4 : f.cofaces.length ? 2.6 : 1.6),
    color: (f) => (f.type === 'lane' && !near ? sand : tone(edgeBase(f), f.id)),
    pick: (f) => `${f.name || (f.type === 'lane' ? 'a passage' : 'a segment')}, a 1-cell (#${f.id}): ${Math.round(f.length)} m, bounding ${f.cofaces.length === 2 ? 'two blocks' : f.cofaces.length === 1 ? 'one block' : 'no block'}` });
  // Their orientations, tail to head.
  const longEnough = E.filter((e) => e.length > 14);
  gis.addArrows(ids.arrows, longEnough, {
    position: (e) => midOf(e.line), direction: (e) => { const a = e.line[0], b = e.line[e.line.length - 1], d = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1; return [(b[0] - a[0]) / d, (b[1] - a[1]) / d]; },
    length: (e) => Math.min(9, 3 + e.length / 20), width: () => 3, z: lifted ? LIFT[1] + 0.4 : 0.7, opacity: 0.9, color: (e) => tone(night ? [0.62, 0.7, 0.78] : [0.36, 0.46, 0.55], e.id) });
  // 0-cells, with the walkers set on them.
  gis.addPoints(ids.vertices, V, { position: (f) => f.at, radius: (f) => (role(f.id) === 'picked' ? 2.6 : f.type === 'portal' ? 1.9 : 1.4), base: () => z(0),
    height: (f) => 3 + (state.initial.walkers[f.id] || 0) / 8, color: (f) => (state.initial.walkers[f.id] ? (role(f.id) === 'other' ? mix(sea, paper, 0.5) : sea) : tone(night ? [0.8, 0.84, 0.88] : [0.3, 0.38, 0.45], f.id)),
    pick: (f) => `a junction, a 0-cell (#${f.id}) on ${f.degree} segment${f.degree === 1 ? '' : 's'}${state.initial.walkers[f.id] ? `, ${state.initial.walkers[f.id]} walkers set` : ''}` });
  // The holes: faces no single cell can be.
  const holes = (c.holes || []).map((ring, k) => ({ ring, k }));
  gis.addFlats(ids.holes, holes, { ring: (f) => f.ring, z: lifted ? LIFT[2] - 0.3 : 0.08, opacity: lifted ? 0.18 : 0.34, order: 0, color: () => (night ? [0.62, 0.46, 0.36] : [0.86, 0.66, 0.52]),
    pick: () => 'a hole: a face with a lane ending inside it, which no single cell can be' });
  gis.addTubes(ids.holeEdges, holes, { path: (f) => [...f.ring, f.ring[0]].map(([x, y]) => [x, y, (lifted ? LIFT[2] : 0) + 0.5]), radius: () => 0.45, color: () => (night ? [0.7, 0.58, 0.48] : [0.66, 0.52, 0.42]), opacity: 0.7, pick: false });
  // By dimension: each cell joined to its faces, + sea and − sand.
  if (lifted) {
    const threads = [];
    const vAt = new Map(V.map((v) => [v.id, v.at]));
    const eMid = new Map(E.map((e) => [e.id, midOf(e.line)]));
    for (const e of E) {
      const m = eMid.get(e.id);
      threads.push({ from: [...m, LIFT[1]], to: [...vAt.get(e.tail), 3], sign: -1, a: e.id, b: e.tail });
      threads.push({ from: [...m, LIFT[1]], to: [...vAt.get(e.head), 3], sign: 1, a: e.id, b: e.head });
    }
    for (const f of B) {
      const cc = centroid(f.ring);
      for (const [e, sign] of f.boundary) if (eMid.has(e)) threads.push({ from: [...cc, LIFT[2]], to: [...eMid.get(e), LIFT[1]], sign, a: f.id, b: e });
    }
    gis.addTubes(ids.incidence, threads, { path: (t) => [t.from, t.to], radius: (t) => (near && (role(t.a) !== 'other' && role(t.b) !== 'other') ? 0.55 : 0.22),
      color: (t) => { const base = t.sign > 0 ? sea : sand; return near && (role(t.a) === 'other' || role(t.b) === 'other') ? mix(base, paper, 0.75) : base; }, opacity: 0.75,
      pick: (t) => `[${t.a} : ${t.b}] = ${t.sign > 0 ? '+1' : '−1'}: ${byId.get(t.b) ? ['the junction', 'the segment'][byId.get(t.b).dim] : 'a face'} #${t.b} is a face of #${t.a}, ${t.sign > 0 ? 'with' : 'against'} its orientation` });
  }
  // 3-cells: the buildings, the walls they share, and by dimension the joins between neighbours.
  const M = state.massing;
  if (M) {
    const bodies = M.buildings.map((b) => ({ ...b, dim: 3, id: b.key }));
    // Each terrace (buildings joined by shared walls) in its own soft hue; a building's warmth, when it has one,
    // instead, against the air of the hour: sea cooler, brick warmer, within 6 K.
    const againstAir = G.diverging(6, 'sea-brick');
    const cat = G.categorical();
    const terraceHue = (b) => { const c = cat(b.terrace % 10); return night ? c.map((v) => v * 0.75) : c.map((v) => v + (1 - v) * 0.35); };
    gis.addSolids(ids.buildings, bodies, {
      ring: (b) => (lifted ? inset(b.ring, 0.6) : b.ring), base: (b) => (lifted ? LIFT[3] : b.base), height: (b) => (lifted ? LIFT[3] + 3 : b.top), opacity: lifted ? 0.8 : 0.62,
      color: (b) => { const T = warmthOf(b.key); return tone(T !== undefined ? againstAir(T - airByHour()[Math.min(23, Math.floor(hourNow))]) : terraceHue(b), b.key); },
      pick: (b) => `${b.name || 'a building'}, a 3-cell: ${Math.round(b.top - b.base)} m high, ${Math.round(b.volume).toLocaleString('en')} m³; it shares walls with ${b.neighbours.length} building${b.neighbours.length === 1 ? '' : 's'}${warmthOf(b.key) !== undefined ? `; ${warmthOf(b.key).toFixed(1)} °C` : ''}` });
    if (!lifted) {
      gis.addSlabs(ids.party, M.walls, { line: (w) => w.line, width: () => 1.4, base: (w) => w.z0, height: (w) => Math.max(0.1, w.z1 - w.z0), color: (w) => (near && (near.star.has(w.a) || near.star.has(w.b)) && (role(w.a) === 'picked' || role(w.b) === 'picked') ? deep : sand), opacity: 0.92, order: 5,
        pick: (w) => `a wall two buildings share, a 2-cell with two cofaces: ${Math.round(w.area)} m²` });
    } else {
      const at = new Map(bodies.map((b) => [b.key, centroid(b.ring)]));
      const pairs = [];
      for (const b of bodies) for (const n of b.neighbours) if (b.key < n && at.has(n)) pairs.push({ a: b.key, b: n });
      gis.addTubes(ids.dual, pairs, { path: (q) => [[...at.get(q.a), LIFT[3] + 5], [...at.get(q.b), LIFT[3] + 5]], radius: (q) => (near && (role(q.a) === 'picked' || role(q.b) === 'picked') ? 1.1 : 0.6),
        color: () => sand, opacity: 0.85, pick: () => 'two buildings that share a wall' });
    }
  }

  // The stores.
  if (c.stores && c.stores.length) gis.addPoints(ids.stores, c.stores, { position: (s) => s.at, radius: () => 3.2, height: (s, k) => ((k === 0 ? state.initial.openA : state.initial.openB) ? 16 : 6), base: () => z(0),
    color: (s, k) => ((k === 0 ? state.initial.openA : state.initial.openB) ? sand : (night ? [0.6, 0.62, 0.64] : [0.66, 0.66, 0.66])), shadow: true,
    pick: (s, k) => `${s.name}${s.use ? ` — ${s.use.replace(/-/g, ' ')}` : ''}: ${(k === 0 ? state.initial.openA : state.initial.openB) ? 'open' : 'shut'} at the start` });
}

// A building's warmth: at the clock's hour in a run, else as set by hand.
function warmthOf(key) {
  if (state.heat) {
    const F = frameAt(state.heat, hourNow);
    if (F && F.temperature[key] !== undefined) return F.temperature[key];
  }
  return state.initial.buildings[key];
}
// The buildings recoloured at an hour, without being rebuilt.
export function recolourBuildings(ctx, ids, hour) {
  hourNow = hour;
  if (state.heat && ctx.gis.layers.get(ids.buildings)) ctx.gis.restyle(ids.buildings, {});
}

// What a picked cell is, for the info card and the analyses: its dimension, its faces and cofaces, what it carries.
export function describeCell(hit, c) {
  const f = hit.feature;
  if (!f || !c || f.dim === undefined) return null;
  const point = hit.point;
  if (f.dim === 3) return { kind: 'cell', cell: f, point, title: f.name || 'a building, a 3-cell',
    lines: [`${Math.round(f.top - f.base)} m high, ${Math.round(f.volume).toLocaleString('en')} m³ on ${Math.round(f.area).toLocaleString('en')} m²`,
      f.neighbours.length ? `it shares walls with ${f.neighbours.length} building${f.neighbours.length === 1 ? '' : 's'}, in a terrace of the district` : 'it stands alone: no wall shared',
      warmthOf(f.key) !== undefined ? `${warmthOf(f.key).toFixed(1)} °C` : 'at the air\'s temperature'] };
  if (f.dim === 0) return { kind: 'cell', cell: f, point, title: `a junction, a 0-cell (#${f.id})`, lines: [`on ${f.degree} segment${f.degree === 1 ? '' : 's'}`, state.initial.walkers[f.id] ? `${state.initial.walkers[f.id]} walkers set on it` : 'no walkers set on it'] };
  if (f.dim === 1) return { kind: 'cell', cell: f, point, title: f.name || (f.type === 'lane' ? 'a passage' : 'a segment'),
    lines: [`a 1-cell (#${f.id}), ${Math.round(f.length)} m, from junction #${f.tail} to #${f.head}`, f.cofaces.length ? `it bounds block${f.cofaces.length === 1 ? '' : 's'} #${f.cofaces.join(' and #')}` : 'it bounds no block'] };
  return { kind: 'cell', cell: f, point, title: `a block, a 2-cell (#${f.id})`,
    lines: [`bounded by ${f.boundary.length} segments, ${Math.round(f.area).toLocaleString('en')} m²`, state.initial.temperature[f.id] !== undefined ? `set to ${state.initial.temperature[f.id]} °C` : `at ${f.temperature} °C`] };
}
