// The landscape and its events: the district as cells; two stores picked on the map (the holdings of Street
// Oceanography) or the nearest pair a lane could join; the lane between them built as three rewrites of the cells,
// in the order they need, and closed again; and a day of the district run as one process — the walkers, the
// warmth of the blocks, the stores, a permit, the lane and its recognition — with its record, the order of its
// events and the balances it keeps. The map shows the same through the layers "junctions, streets and blocks"
// and "the landscape running", drawn here when they are off.

import { emit, ensureDistrict, on, request, runDay, setStores, state } from '../shared.js';

const big = (x) => Math.round(x).toLocaleString('en');
const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.floor((x % 1) * 60)).padStart(2, '0')}`;

export async function render(root, ctx) {
  const { h, svg, busy, failure, fmt, css } = ctx.dom;
  const { lines, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The landscape and its events'),
    h('span', { class: 'lede' }, 'A district of Melbourne as cells — junctions, streets, blocks — two stores in it, and what happens there in a day, each event a change of the cells or of what they carry.')));
  const head = h('div', {}), laneBox = h('div', {}), dayBox = h('div', {});
  root.append(head, laneBox, dayBox);
  busy(head, 'cutting the district from the city model and finding two stores a lane could join…');
  try { await ensureDistrict(ctx); } catch (e) { return failure(head, e); }

  // The map: the layers' drawing when they are off.
  const shown = [];
  const showOnMap = async (which) => {
    if (state.active.has(which) || shown.includes(which)) return;
    const mod = await import(`../layers/${which === 'cells' ? 'cells' : 'landscape'}.js`);
    const handle = await mod.layer({ ...ctx, layerId: (s = '') => `landscape-dynamics/tool-${which}${s}` });
    shown.push(which);
    let alive = true;
    ctx.clock.on((c) => { if (alive && handle.setTime) handle.setTime(c); });
    ctx.onDispose(() => {
      alive = false;
      if (handle.dispose) handle.dispose();
      for (const id of [...ctx.gis.layers.keys()]) if (id.startsWith(`landscape-dynamics/tool-${which}`)) ctx.gis.remove(id);
    });
  };
  await showOnMap('cells');

  // ---- the district and the stores ----
  const pickNote = h('p', { class: 'note' });
  const drawHead = () => {
    const D = state.district, t = D.topology;
    pickNote.textContent = `The two stores: ${D.stores[0].name}${D.stores[0].street ? ` on ${D.stores[0].street}` : ''} and ${D.stores[1].name}${D.stores[1].street ? ` on ${D.stores[1].street}` : ''}. Pick two others on the map (turn on Street Oceanography's "shops and services"), facing one block.`;
    head.replaceChildren(h('h2', {}, 'The district as cells'),
      ctx.finding({ title: `${t.counts[0]} junctions, ${t.counts[1]} segments and ${t.counts[2]} blocks make one regular cell complex`, kind: 'proved',
        statement: `Every block is a disc bounded by a circuit of streets, every segment joins two distinct junctions: the complex is regular, so its cells and their faces determine it. ${t.betti[0] === 1 ? 'It is connected' : `It falls into ${t.betti[0]} pieces within the box (b₀ = ${t.betti[0]}: ways the box cuts off, or joined to the streets only by footways left out)`}, and it has ${t.betti[1]} hole${t.betti[1] === 1 ? '' : 's'} (b₁ = ${t.betti[1]}): the block${t.holes === 1 ? '' : 's'} with a dead end or a loop of lanes inside, which no single cell can be. Euler characteristic ${t.euler} = b₀ − b₁.`,
        method: 'OpenStreetMap\'s ways where people walk, cut to the box Elizabeth–Russell, Flinders–Lonsdale and split at every shared point; the faces traced from the directions at each junction; each face bounded by a circuit a 2-cell; regularity tested cell by cell; homology over the integers (morphe).',
        limits: `${D.dropped && D.dropped.length ? `Left out, crossing another way without a junction (a ramp, a way at another level): ${[...new Set(D.dropped)].join(', ')}. ` : ''}The footways drawn beside streets are left out so that no street counts twice.` }),
      pickNote);
  };
  drawHead();
  const picked = [];
  ctx.onSelect(async (sel) => {
    if (!sel || sel.kind !== 'holding' || !state.candidates) return;
    const c = state.candidates.find((x) => x.id === sel.id);
    if (!c) { pickNote.textContent = `${sel.name} lies outside the district.`; return; }
    if (picked.length === 2) picked.length = 0;
    picked.push(c);
    if (picked.length === 1) { pickNote.textContent = `First: ${c.name}. Now the second store.`; return; }
    busy(laneBox, 'the district with the two stores…');
    try { await setStores(ctx, picked.slice()); } catch (e) { return failure(laneBox, e); }
    drawHead();
    drawLane();
  });

  // ---- the lane ----
  const drawLane = () => {
    const D = state.district;
    laneBox.replaceChildren(h('h2', {}, 'A lane between them'));
    if (!D.lane || D.lane.refused) { laneBox.append(h('p', { class: 'note' }, `No lane: ${D.lane ? D.lane.refused : 'two stores are needed'}.`)); return; }
    const build = h('button', { class: 'primary' }, state.lane ? 'built — show it again' : 'build the lane');
    const out = h('div', {});
    laneBox.append(h('p', {}, `A passage ${Math.round(D.lane.length)} m long through the block between them, from the street each faces.`), h('p', {}, build), out);
    build.addEventListener('click', async () => {
      busy(out, 'rewriting the cells…');
      let L;
      try { L = await ctx.ask('lane', request()); } catch (e) { return failure(out, e); }
      state.lane = L;
      emit('lane');
      const names = { subdivide: 'a street subdivided where the lane meets it', split: 'the block split by the lane' };
      out.replaceChildren(
        h('ol', {}, L.steps.map((s) => h('li', {}, `${names[s.kind] || s.kind}${s.kind === 'subdivide' ? ` (${s.match})` : ''}: ${s.deleted.length} cell${s.deleted.length === 1 ? '' : 's'} removed, ${s.created.length} made `,
          h('span', { class: `chip ${s.refinement ? 'proved' : 'model'}` }, s.refinement ? 'a refinement: the same street, described finer' : 'an intervention: a new place')))),
        orderDiagram(h, svg, css, ['subdivide one street', 'subdivide the other', 'split the block'], L.order),
        ctx.finding({ title: 'Three rewrites; the first two in either order, the split after both', kind: 'proved',
          statement: `Each rewrite replaces a few cells by others and keeps the rest: a street becomes two with a junction between them, its blocks bounded by both halves; the block becomes two with the lane between them. The two streets can be cut in either order — neither reads what the other changes — and the block can be split only once both are cut, so the ${L.order.records} orders of the three are one history. A subdivision is a refinement — the street is the same street, its description finer, and nothing that held of it changes — while the split is an intervention: a passage that was not there. The district keeps its Betti numbers (b₀ = ${L.after.betti[0]}, b₁ = ${L.after.betti[1]}): a lane through a block adds two junctions, three segments and a block, and no hole. Closed again — the two blocks merged across it, the streets merged at its ends — the district is the same up to isomorphism: ${L.undo.isomorphic ? 'yes' : 'no'}.`,
          method: 'Double pushouts on the complex\'s incidence graph (morphe\'s graph rewriting), each result tested for regularity; the order from what each rewrite read and wrote.',
          against: `the blocks' heat capacity ${L.capacityKept === 0 ? 'kept exactly' : `kept to a relative ${fmt(Math.abs(L.capacityKept), 2)}`}; the walkers kept (a new junction starts with none)` }));
    });
    if (state.lane) build.click();
  };
  drawLane();

  // ---- a day ----
  const P = { ...state.rules, seed: 3 };
  const slider = (key, name, min, max, step, unit = '') => {
    const input = h('input', { type: 'range', min, max, step, value: P[key] });
    const o = h('output', {}, `${fmt(P[key])}${unit}`);
    input.addEventListener('input', () => { o.textContent = `${fmt(+input.value)}${unit}`; });
    input.addEventListener('change', () => { P[key] = +input.value; if (key in state.rules) state.rules[key] = P[key]; });
    return h('label', { class: 'slider' }, name, input, o);
  };
  const go = h('button', { class: 'primary' }, 'run the day');
  const result = h('div', {});
  dayBox.append(h('h2', {}, 'A day of the district'),
    h('p', { class: 'caption' }, 'The stores open and shut by the hour, each more readily when the other is open (η) and when the lane is there (ξ); a permit may be granted, a lane built under it, removed while nobody has recognised it, and recognised as it is used. At 22:00 both stores shut.'),
    h('div', { class: 'panel' }, slider('eta', 'the other store open', 0, 2, 0.1), slider('xi', 'the lane there', 0, 2, 0.1), slider('permit', 'a permit, per hour', 0, 1, 0.05), slider('build', 'built, per hour', 0, 2, 0.1),
      slider('recognise', 'recognised, per hour', 0, 2, 0.1), slider('seed', 'which day', 1, 40, 1)), h('p', {}, go), result);
  const showRun = () => {
    const R = state.run;
    if (!R) return;
    const rows = R.events.slice(0, 60).map((e) => h('tr', {}, h('td', { class: 'mono' }, hhmm(e.time)), h('td', { class: 'left' }, e.mechanism, e.forced ? h('span', { class: 'chip model', style: 'margin-left:6px' }, 'at the hour') : null, e.rewrites && e.rewrites.length ? h('span', { class: 'chip proved', style: 'margin-left:6px' }, `${e.rewrites.length} rewrites`) : null), h('td', { class: 'left faint' }, e.match)));
    const total = R.frames.map((f) => Object.values(f.walkers).reduce((a, x) => a + x, 0));
    const warm = R.frames.map((f) => { const v = Object.values(f.temperature); return v.reduce((a, x) => a + x, 0) / Math.max(1, v.length); });
    const times = R.frames.map((f) => f.time);
    const changes = R.events.filter((e) => e.rewrites && e.rewrites.length).length;
    const lane = R.events.find((e) => e.mechanism === 'the lane is built');
    result.replaceChildren(
      lines([{ name: 'walkers in the district', x: times, y: total, color: css('--accent') }], { xlabel: 'hour', ylabel: 'walkers', height: 170 }),
      lines([{ name: 'mean warmth of the blocks', x: times, y: warm, color: css('--minor') }], { xlabel: 'hour', ylabel: '°C', height: 150 }),
      lines([{ name: 'first store open', x: times, y: R.frames.map((f) => (f.a ? 1 : 0)), color: css('--series-1') }, { name: 'second store open', x: times, y: R.frames.map((f) => (f.b ? 0.9 : -0.1)), color: css('--series-4'), dash: '4 3' },
        { name: 'the lane there', x: times, y: R.frames.map((f) => (f.lane ? 0.5 : -0.2)), color: css('--minor') }], { xlabel: 'hour', ylabel: 'open (1) or not', height: 140 }),
      legend([{ name: R.strata[0].stores[0].name, color: css('--series-1') }, { name: R.strata[0].stores[1].name, color: css('--series-4') }, { name: 'the lane', color: css('--minor') }]),
      ctx.finding({ title: `${R.events.length} events; ${lane ? `the lane built at ${hhmm(lane.time)}` : 'no lane this day'}; ${R.strata.every((s) => s.regular) ? 'every configuration regular' : 'a configuration not regular'}`, kind: 'model',
        statement: `The cells changed ${changes} time${changes === 1 ? '' : 's'} as events came, through ${R.strata.length} configuration${R.strata.length === 1 ? '' : 's'} — ${R.strata.filter((s) => s.lane).length ? 'with the lane the district has two junctions, three segments and a block more' : 'the lane never built'} — and every configuration was a regular complex with b₁ = ${R.strata[0].betti[1]}; the other events changed only what the cells carry. Walkers are kept exactly: ${big(R.balance.end)} at midnight = ${big(R.balance.start)} + ${big(R.balance.entered)} entered − ${big(R.balance.left)} left (residual ${fmt(R.balance.residual, 2)}), the rewrites carrying them across.`,
        method: `A piecewise-deterministic process on the strata of the complex: between events the walkers hop along the streets toward open stores and leave, and the blocks warm toward the air and the Sun of ${R.date}; events at rates that read the state, the closing at 22:00 a boundary of the domain. Walkers entering: ${R.entriesSource}.`,
        limits: 'A representation: every rate is declared, not measured, to explore how structure, use and recognition can change one another.' }),
      h('h3', {}, 'The record'), h('table', { class: 'data', style: 'width:100%' }, h('tbody', {}, rows)), R.events.length > 60 ? h('p', { class: 'faint' }, `…and ${R.events.length - 60} more.`) : null,
      h('h3', {}, 'What could happen, and what did'), runOrder(h, svg, css, R),
      ctx.finding({ title: `${R.order.steps.length} steps of events that could happen together`, kind: 'proved',
        statement: `An event comes after another only if it reads what the other wrote: a store opening after the lane it answers to was built, the lane after the permit, the recognition after the lane. Events that touch nothing in common — the permit and the two stores' first openings — are in no order; ${R.order.records ? `the record is one of ${Number(R.order.records).toLocaleString('en')} orderings of the same history` : 'many orderings of the record are the same history'}. The steps group what could have happened at once.`,
        method: 'Each event\'s reads and writes (the stores\' activity, the permit, the lane, its recognition, and every cell and incidence a rewrite touched); dependent events ordered, the order closed transitively; its Foata steps (morphe\'s traces).' }));
  };
  go.addEventListener('click', async () => {
    busy(result, 'running the day: walkers, warmth, stores, permit, lane…');
    try { await runDay(ctx, { eta: P.eta, xi: P.xi, permit: P.permit, build: P.build, recognise: P.recognise, seed: P.seed }); } catch (e) { return failure(result, e); }
    await showOnMap('landscape');
    showRun();
  });
  const off = on((what) => { if (what === 'run') showRun(); });
  ctx.onDispose(off);
  if (state.run) { await showOnMap('landscape'); showRun(); }
}

// Three steps and their order: the two subdivisions side by side, the split beneath both.
function orderDiagram(h, svg, css, names, order) {
  const W = 420, H = 120, pos = [[90, 30], [330, 30], [210, 95]];
  const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  s.append(svg('defs', {}, svg('marker', { id: 'ld-arrow', viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: 'auto-start-reverse' }, svg('path', { d: 'M0,0 L10,5 L0,10 z', fill: css('--ink-2') }))));
  for (const [a, b] of order.covers) s.append(svg('line', { x1: pos[a][0], y1: pos[a][1] + 10, x2: pos[b][0], y2: pos[b][1] - 12, stroke: css('--ink-2'), 'stroke-width': 1.2, 'marker-end': 'url(#ld-arrow)' }));
  names.forEach((n, k) => s.append(svg('text', { x: pos[k][0], y: pos[k][1], 'text-anchor': 'middle', 'font-size': 13, fill: css('--ink') }, n)));
  return h('div', { class: 'panel', style: 'margin:8px 0' }, s, h('div', { class: 'faint' }, `${order.records} orders of the three steps are one history`));
}

// The record's events on the day's hours, by kind, with the arrows of its order (its Hasse diagram).
function runOrder(h, svg, css, R) {
  const kinds = ['a store opens', 'a store closes', 'the permit is granted', 'the lane is built', 'the lane is removed', 'the lane is recognised', 'the stores close for the night'];
  const label = ['opens', 'shuts', 'permit', 'lane built', 'lane removed', 'recognised', 'night'];
  const W = 600, H = 34 + 26 * kinds.length, x = (t) => 70 + (t / 24) * (W - 90), y = (e) => 22 + 26 * Math.max(0, kinds.indexOf(e.mechanism));
  const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%' });
  kinds.forEach((k, i) => s.append(svg('text', { x: 4, y: 26 + 26 * i, 'font-size': 11, fill: css('--ink-2') }, label[i])));
  // The registry: the hours each kind of event could happen, as the state allowed it.
  for (const F of R.frames)
    for (const k of F.enabled || []) {
      const i = kinds.indexOf(k);
      if (i >= 0) s.append(svg('rect', { x: x(F.time), y: 22 + 26 * i - 7, width: Math.max(0.5, x(F.time + 0.25) - x(F.time)), height: 14, fill: css('--accent-soft') }));
    }
  for (let t = 0; t <= 24; t += 6) s.append(svg('text', { x: x(t), y: H - 4, 'font-size': 10, 'text-anchor': 'middle', fill: css('--ink-faint') }, `${t}:00`));
  for (const [a, b] of R.order.covers) {
    const A = R.events[a], B = R.events[b];
    s.append(svg('line', { x1: x(A.time), y1: y(A), x2: x(B.time), y2: y(B), stroke: css('--accent'), 'stroke-opacity': 0.35, 'stroke-width': 1 }));
  }
  R.events.forEach((e) => s.append(svg('circle', { cx: x(e.time), cy: y(e), r: 3.6, fill: e.forced ? css('--minor') : css('--accent-deep') })));
  return h('div', { class: 'panel', style: 'margin:8px 0' }, s, h('div', { class: 'faint' }, 'each kind of event on its row: shaded while the state allowed it, a dot when it happened; a line joins an event to the next ones that read what it wrote'));
}
