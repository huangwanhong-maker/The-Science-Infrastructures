// The district: its complex, its freedom, its dynamics. One place to work on the landscape as the formalism has it,
// in three tabs. The complex: the district's cells — junctions, streets, blocks and its buildings, which share their
// walls — drawn on the map on the ground or by dimension, a picked cell with its faces and what it bounds. The
// freedom: what the cells carry at the start, set by hand, and the cells themselves rewritten (a street cut, a passage
// opened through a block, two blocks joined), each edit checked and kept as a history. The dynamics: the landscape
// and the buildings' warmth run on that complex from that state from the clock's hour, drawn on the map as the
// clock plays. Each tab is drawn from the shared state, so an edit, a pick or a run redraws only what is shown.

import { clearFreedom, edit, emit, ensureDistrict, ensureMassing, on, resetEdits, runFrom, setFreedom, state, undoEdit } from '../shared.js';
import { complexIds, paintComplex, recolourBuildings } from '../layers/cells.js';
import { brief, figures, fold, hint, seg, slider, swatches, tabs, toggle } from '../ui.js';

const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.round((x % 1) * 60) % 60).padStart(2, '0')}`;
const MODES = [
  ['look', 'Look', 'Pick any cell to see it.'],
  ['cut', 'Cut a street', 'Pick a street where it should be cut.'],
  ['passage', 'Open a passage', 'Pick a block, then two points on its edge.'],
  ['join', 'Join two blocks', 'Pick the street between them.'],
];

export async function render(root, ctx) {
  const { h, svg, fmt, css } = ctx.dom;
  const { lines } = ctx.chart;
  const { gis } = ctx;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The district'), h('span', { class: 'lede' }, 'Its cells, what they carry, and how they change.')));
  const body = h('div', {});
  root.append(body);
  body.append(hint(h, 'cutting the district from the city model…'));
  try { await ensureDistrict(ctx); } catch (e) { return ctx.dom.failure(body, e); }

  // The complex on the map, drawn here when its layer is off.
  const ids = complexIds((s) => `landscape-dynamics/tool-complex${s}`);
  let alive = true;
  const paint = () => { if (alive && !state.active.has('cells') && state.district) paintComplex(ctx, ids, state.district); };
  ctx.onDispose(() => { alive = false; for (const id of Object.values(ids)) gis.remove(id); });

  let mode = 'look', passage = null, notice = '', running = false, P = { span: 6, seed: 3 };

  // ---- the complex ----
  const describe = (p) => {
    if (p.dim === 3) return { title: p.name || 'A building', kind: '3-cell', lines: [`${Math.round(p.top - p.base)} m high · ${Math.round(p.volume).toLocaleString('en')} m³`, p.neighbours.length ? `shares walls with ${p.neighbours.length}` : 'stands alone'] };
    if (p.dim === 2) return { title: 'A block', kind: `2-cell #${p.id}`, lines: [`${p.boundary.length} streets round it · ${Math.round(p.area).toLocaleString('en')} m²`] };
    if (p.dim === 1) return { title: p.name || 'A street', kind: `1-cell #${p.id}`, lines: [`${Math.round(p.length)} m · bounds ${p.cofaces.length === 2 ? 'two blocks' : p.cofaces.length === 1 ? 'one block' : 'no block'}`] };
    return { title: 'A junction', kind: `0-cell #${p.id}`, lines: [`on ${p.degree} street${p.degree === 1 ? '' : 's'}`] };
  };
  const pickedCard = () => {
    const p = state.picked;
    if (!p) return hint(h, 'Pick a cell on the map to see its faces and what it bounds.');
    const d = describe(p);
    const detail = p.dim === 2 ? `its boundary: ${p.boundary.map(([e, sg]) => `${sg > 0 ? '+' : '−'}${e}`).join(' ')}` : p.dim === 1 ? `∂ = #${p.head} − #${p.tail}` : null;
    return h('div', { class: 'ld-card' }, h('div', { class: 'ld-row' }, h('span', { class: 'ld-title' }, d.title), h('span', { class: 'ld-sub' }, d.kind)),
      ...d.lines.map((l) => h('div', { class: 'ld-sub' }, l)), detail ? fold(h, 'orientation', h('p', {}, detail)) : null);
  };
  const complexPanel = () => {
    const c = state.district, t = c.topology, [V, E, F] = t.counts, M = state.massing;
    const night = ctx.theme === 'night', rgb = (a) => `rgb(${a.map((v) => Math.round(v * 255)).join(',')})`;
    return h('div', {},
      figures(h, [{ value: V, label: 'junctions' }, { value: E, label: 'streets' }, { value: F, label: 'blocks' }, M && { value: M.buildings.length, label: 'buildings' }, { value: t.holes, label: 'holes' }]),
      hint(h, `Betti numbers ${t.betti[0]}, ${t.betti[1]} · χ = ${String(t.euler).replace('-', '−')}${M ? ` · ${M.walls.length} shared walls in ${M.terraces} terraces` : ''}${state.edits.length ? ` · ${state.edits.length} edit${state.edits.length === 1 ? '' : 's'}` : ''}`),
      h('div', { class: 'ld-row', style: 'margin:10px 0 12px' }, seg(h, [['ground', 'On the ground'], ['dimension', 'By dimension']], state.view, (v) => { state.view = v; emit('view'); })),
      pickedCard(),
      fold(h, 'how to read the map',
        swatches(h, [{ color: rgb(night ? [0.3, 0.44, 0.55] : [0.6, 0.74, 0.84]), label: 'block' }, { color: rgb(night ? [0.75, 0.83, 0.9] : [0.2, 0.32, 0.43]), label: 'street', kind: 'line' },
          { color: rgb(night ? [0.8, 0.84, 0.88] : [0.3, 0.38, 0.45]), label: 'junction', kind: 'dot' }, { color: css('--series-5'), label: 'building (a hue per terrace)' },
          { color: css('--minor'), label: 'shared wall', kind: 'line' }, { color: rgb([0.86, 0.66, 0.52]), label: 'hole' }]),
        h('p', {}, 'A hole is a face with a lane ending inside it: no single cell can be such a face. By dimension, each cell is joined to its faces — sea +1, sand −1 — and each building to those it shares a wall with.')));
  };

  // ---- the freedom ----
  const numberInput = (value, onSet, step = 1) => {
    const input = h('input', { type: 'number', value: value ?? '', step, class: 'ld-number', style: 'width:6.5em' });
    input.addEventListener('change', () => onSet(input.value === '' ? null : +input.value));
    return input;
  };
  const carryControl = () => {
    const p = state.picked;
    if (!p) return hint(h, 'Pick a junction, a block or a building to set what it carries.');
    const d = describe(p);
    if (p.dim === 0) return h('div', { class: 'ld-row' }, h('span', {}, d.title), numberInput(state.initial.walkers[p.id], (v) => setFreedom('walkers', p.id, v)), h('span', { class: 'ld-sub' }, 'walkers'));
    if (p.dim === 2) return h('div', { class: 'ld-row' }, h('span', {}, d.title), numberInput(state.initial.temperature[p.id] ?? p.temperature, (v) => setFreedom('temperature', p.id, v), 0.5), h('span', { class: 'ld-sub' }, '°C'));
    if (p.dim === 3) return h('div', { class: 'ld-row' }, h('span', {}, d.title), numberInput(state.initial.buildings[p.key], (v) => setFreedom('buildings', p.key, v), 0.5), h('span', { class: 'ld-sub' }, '°C (the air if empty)'));
    return hint(h, 'A street carries the walkers its junctions send along it: set walkers at a junction.');
  };
  const freedomPanel = () => {
    const c = state.district, names = (c.stores || []).map((s) => s.name);
    const set = [];
    const nW = Object.keys(state.initial.walkers).length, nT = Object.keys(state.initial.temperature).length, nB = Object.keys(state.initial.buildings).length;
    if (nW) set.push(`walkers at ${nW} junction${nW === 1 ? '' : 's'}`);
    if (nT) set.push(`${nT} block${nT === 1 ? '' : 's'} warmed`);
    if (nB) set.push(`${nB} building${nB === 1 ? '' : 's'} warmed`);
    const items = state.edits.map((e, k) => {
      const done = c.edits[k];
      const what = e.kind === 'cut' ? 'A street cut' : e.kind === 'passage' ? 'A passage opened' : 'Two blocks joined';
      return h('li', {}, h('span', {}, what), done ? h('span', { class: `chip ${done.refinement ? 'proved' : 'model'}` }, done.refinement ? 'refinement' : 'intervention') : null,
        h('span', { class: 'faint' }, done ? `${done.steps.length} rewrite${done.steps.length === 1 ? '' : 's'}` : ''));
    });
    const modeHint = passage ? (passage.from ? 'Now a second point, on another of its streets.' : 'Now a point on one of its streets.') : MODES.find((m) => m[0] === mode)[2];
    return h('div', {},
      h('div', { class: 'ld-label' }, 'What the cells carry at the start'),
      carryControl(),
      h('div', { style: 'margin-top:10px' }, names[0] ? toggle(h, `${names[0]} open`, state.initial.openA, (v) => setFreedom('openA', null, v)) : null,
        names[1] ? toggle(h, `${names[1]} open`, state.initial.openB, (v) => setFreedom('openB', null, v)) : null, toggle(h, 'A permit for the lane', state.initial.permit, (v) => setFreedom('permit', null, v))),
      set.length ? h('div', { class: 'ld-row ld-hint' }, `Set: ${set.join(', ')}.`, h('button', { class: 'ld-link', onclick: () => clearFreedom() }, 'clear')) : null,
      h('div', { class: 'ld-label', style: 'margin-top:22px' }, 'Change the cells'),
      seg(h, MODES.map(([k, n]) => [k, n]), mode, (v) => { mode = v; passage = null; notice = ''; redraw(); }),
      hint(h, modeHint),
      notice ? h('p', { class: 'ld-hint', style: 'color:var(--bad)' }, notice) : null,
      items.length ? h('ul', { class: 'ld-list' }, items) : null,
      items.length ? h('div', { class: 'ld-row' }, h('button', { class: 'ld-link', onclick: async () => { notice = ''; await undoEdit(ctx); } }, 'undo the last'),
        h('button', { class: 'ld-link', onclick: async () => { notice = ''; await resetEdits(ctx); } }, 'back to the map as it is')) : null,
      items.length > 1 ? fold(h, 'their order', h('p', {}, `These ${items.length} edits need ${(c.order && c.order.steps ? c.order.steps.length : 1)} step${c.order && c.order.steps && c.order.steps.length === 1 ? '' : 's'}: an edit must come after another only if it rewrites cells the other made.`)) : null);
  };

  // ---- the dynamics ----
  const startHour = () => Math.min(23.5, Math.floor(ctx.clock.hour * 4) / 4);
  const run = async () => {
    const s = startHour(), e = Math.min(24, s + P.span);
    running = true; redraw();
    try { await runFrom(ctx, { start: s, end: e, seed: P.seed }); } catch (x) { notice = String(x.message || x); }
    running = false;
    await showOnMap();
    ctx.clock.set({ hour: s });
    const sl = document.querySelector('.clock input[type=range]');
    if (sl) sl.value = s;
    redraw();
  };
  const dynamicsPanel = () => {
    const s = startHour(), e = Math.min(24, s + P.span), R = state.run && state.run.start !== undefined ? state.run : null, Hh = state.heat;
    const out = [
      h('div', { class: 'ld-row' }, h('button', { class: 'primary', onclick: run, disabled: running || null }, running ? 'running…' : `Run ${hhmm(s)} – ${hhmm(e)}`),
        h('span', { class: 'ld-hint' }, `${ctx.clock.date}, from what is set`)),
      fold(h, 'conditions',
        slider(h, fmt, { label: 'for', min: 1, max: 12, step: 0.5, value: P.span, unit: ' h', onChange: (v) => { P.span = v; redraw(); } }),
        slider(h, fmt, { label: 'which run', min: 1, max: 40, step: 1, value: P.seed, onChange: (v) => { P.seed = v; } }),
        slider(h, fmt, { label: 'the air, mean', min: -5, max: 35, step: 0.5, value: state.air.mean, unit: ' °C', onChange: (v) => { state.air.mean = v; } }),
        slider(h, fmt, { label: 'warmest above it', min: 0, max: 12, step: 0.5, value: state.air.swing, unit: ' °C', onChange: (v) => { state.air.swing = v; } })),
    ];
    if (R) {
      const times = R.frames.map((f) => f.time);
      const total = R.frames.map((f) => Object.values(f.walkers).reduce((a, x) => a + x, 0));
      const changes = R.events.filter((x) => x.rewrites && x.rewrites.length).length;
      const ht = Hh ? Hh.frames.map((f) => Object.values(f.temperature)) : null;
      out.push(
        figures(h, [{ value: R.events.length, label: 'events' }, { value: changes, label: 'rewrites' }, { value: Math.round(Math.max(...total)).toLocaleString('en'), label: 'walkers at most' },
          ht && { value: `${fmt(Math.min(...ht.flat()), 3)}–${fmt(Math.max(...ht.flat()), 3)} °C`, label: 'the buildings' }]),
        h('div', {}, h('div', { class: 'ld-label' }, 'Walkers in the district'), lines([{ name: 'walkers in the district', x: times, y: total, color: css('--accent-deep'), width: 2 }], { xlabel: 'hour', height: 120 })),
        Hh ? buildingsChart(Hh) : null,
        hint(h, 'Press ▶ on the clock to watch it on the map; the layer "the district\'s spacetime" shows the whole run at once.'),
        fold(h, 'what could happen, and what did', registry(h, svg, css, R)),
        brief(ctx, { title: 'Walkers kept exactly', kind: 'model',
          statement: `${Math.round(R.balance.end).toLocaleString('en')} at the end = ${Math.round(R.balance.start).toLocaleString('en')} at the start + ${Math.round(R.balance.entered).toLocaleString('en')} entered − ${Math.round(R.balance.left).toLocaleString('en')} left.`,
          method: `A piecewise-deterministic process on the strata of the complex, from the state set; walkers entering: ${R.entriesSource}. The buildings: Δ₃ on their massing, through shared walls and to the air and the ground, under the Sun.`,
          limits: 'Every rate is declared, not measured.' }));
    } else if (!running) out.push(hint(h, 'Nothing run yet.'));
    return h('div', {}, ...out.filter(Boolean));
  };
  const buildingsChart = (Hh) => {
    const times = Hh.frames.map((f) => f.time);
    const v = Hh.frames.map((f) => Object.values(f.temperature));
    const series = [{ name: 'coolest to warmest', x: times, y: v.map((a) => Math.min(...a)), y2: v.map((a) => Math.max(...a)), color: css('--minor'), opacity: 0.18 },
      { name: 'their mean', x: times, y: v.map((a) => a.reduce((s, x) => s + x, 0) / a.length), color: css('--minor-deep'), width: 2 }];
    const p = state.picked;
    if (p && p.dim === 3) series.push({ name: p.name || 'the building picked', x: times, y: Hh.frames.map((f) => f.temperature[p.key]), color: css('--accent-deep'), width: 2 });
    return h('div', {}, h('div', { class: 'ld-label' }, `The buildings, °C${p && p.dim === 3 ? ` · ${p.name || 'the one picked'}` : ''}`), lines(series, { xlabel: 'hour', height: 120 }));
  };
  let shown = false;
  const showOnMap = async () => {
    if (shown || state.active.has('landscape')) return;
    shown = true;
    const mod = await import('../layers/landscape.js');
    const handle = await mod.layer({ ...ctx, layerId: (s = '') => `landscape-dynamics/tool-landscape${s}` });
    ctx.clock.on((c) => { if (alive && handle.setTime) handle.setTime(c); });
    ctx.onDispose(() => { if (handle.dispose) handle.dispose(); for (const id of [...gis.layers.keys()]) if (id.startsWith('landscape-dynamics/tool-landscape')) gis.remove(id); });
  };

  // ---- the tabs ----
  let current = 'Complex';
  const T = tabs(h, [{ name: 'Complex', panel: complexPanel }, { name: 'Freedom', panel: freedomPanel }, { name: 'Dynamics', panel: dynamicsPanel }],
    { remember: 'landscape-workspace', onShow: (n) => { current = n; } });
  const redraw = () => { if (alive) T.show(current); };
  body.replaceChildren(T.el);
  paint();
  // The buildings as glass so that the cells show through (the rail's own key, so that its state shows); the district
  // framed in the part of the map the drawer leaves free.
  const glass = [...document.querySelectorAll('[data-role="building-form"] button')].find((b) => b.textContent === 'glass');
  if (glass && !glass.classList.contains('on')) glass.click();
  {
    const B = state.box, cx = (B.x0 + B.x1) / 2, cy = (B.y0 + B.y1) / 2, size = Math.max(B.x1 - B.x0, B.y1 - B.y0);
    const o = gis.camera.position.clone().sub(gis.controls.target), n = Math.hypot(o.x, o.z) || 1;
    gis.flyTo([cx + (o.z / n) * 0.32 * size, cy + (o.x / n) * 0.32 * size], 1.5 * size);
  }
  const off = on((what) => {
    if (!alive) return;
    if (['district', 'freedom', 'picked', 'view', 'massing'].includes(what)) { paint(); redraw(); }
  });
  ctx.onDispose(off);
  ctx.clock.on((c) => { if (alive) { recolourBuildings(ctx, ids, c.hour); if (current === 'Dynamics' && !running) redraw(); } });
  ensureMassing(ctx).catch((x) => console.error(x));

  // Picks on the complex, wherever it is drawn.
  gis.onPick(async (hit) => {
    if (!alive || !hit || !hit.layer || !/^landscape-dynamics\/(cells|tool-complex):(blocks|edges|vertices|buildings)$/.test(hit.layer)) return;
    const cell = hit.feature, at = hit.point ? [hit.point[0], hit.point[1]] : null;
    const attempt = async (e) => {
      notice = '';
      let r;
      try { r = await edit(ctx, e); } catch (x) { notice = String(x.message || x); }
      if (r && r.refused) notice = `Not done: ${r.refused}.`;
      passage = null;
      redraw();
    };
    if (mode === 'look') { state.picked = cell; emit('picked'); return; }
    if (mode === 'cut') { if (cell.dim !== 1) { notice = 'A street is cut: pick a street.'; return redraw(); } return attempt({ kind: 'cut', segment: cell.id, at }); }
    if (mode === 'join') { if (cell.dim !== 1) { notice = 'Pick the street between the two blocks.'; return redraw(); } return attempt({ kind: 'join', segment: cell.id }); }
    if (mode === 'passage') {
      if (!passage) {
        if (cell.dim !== 2) { notice = 'First the block the passage goes through.'; return redraw(); }
        passage = { block: cell.id };
        state.picked = cell; notice = '';
        emit('picked');
        return;
      }
      if (!passage.from) { passage.from = at; return redraw(); }
      return attempt({ kind: 'passage', block: passage.block, from: passage.from, to: at });
    }
  });
}

// The run's events on its hours, by kind: shaded while the state allowed the kind, a dot when it happened, a line to
// the next events that read what it wrote.
function registry(h, svg, css, R) {
  const kinds = ['a store opens', 'a store closes', 'the permit is granted', 'the lane is built', 'the lane is removed', 'the lane is recognised', 'the stores close for the night'];
  const label = ['opens', 'shuts', 'permit', 'lane built', 'lane removed', 'recognised', 'night'];
  const t0 = R.start ?? 0, t1 = R.end ?? 24;
  const W = 600, H = 34 + 24 * kinds.length, x = (t) => 76 + ((t - t0) / (t1 - t0 || 1)) * (W - 92), y = (e) => 20 + 24 * Math.max(0, kinds.indexOf(e.mechanism));
  const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%' });
  kinds.forEach((k, i) => s.append(svg('text', { x: 4, y: 24 + 24 * i, 'font-size': 11, fill: css('--ink-faint') }, label[i])));
  for (const F of R.frames)
    for (const k of F.enabled || []) {
      const i = kinds.indexOf(k);
      if (i >= 0) s.append(svg('rect', { x: x(F.time), y: 20 + 24 * i - 6, width: Math.max(0.5, x(F.time + 0.25) - x(F.time)), height: 12, fill: css('--accent-soft') }));
    }
  for (let t = Math.ceil(t0); t <= t1; t += Math.max(1, Math.round((t1 - t0) / 6))) s.append(svg('text', { x: x(t), y: H - 4, 'font-size': 10, 'text-anchor': 'middle', fill: css('--ink-faint') }, `${t}:00`));
  for (const [a, b] of R.order.covers || []) {
    const A = R.events[a], B = R.events[b];
    s.append(svg('line', { x1: x(A.time), y1: y(A), x2: x(B.time), y2: y(B), stroke: css('--accent'), 'stroke-opacity': 0.3, 'stroke-width': 1 }));
  }
  R.events.forEach((e) => s.append(svg('circle', { cx: x(e.time), cy: y(e), r: 3.2, fill: e.forced ? css('--minor') : css('--accent-deep') })));
  return h('div', {}, s, h('p', { class: 'ld-hint' }, 'Shaded while the state allowed each kind of event; a dot when it happened; a line to the events that read what it wrote.'));
}
