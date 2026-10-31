// The district: its complex, its freedom, its dynamics. One place to work on the landscape as the formalism has it:
// the district as a cell complex — its streets and blocks, and its buildings as 3-cells sharing their walls — shown
// as a complex (on the ground or by dimension, a picked cell with its faces and what it bounds); the freedom on top of it — the state each cell carries, set by hand (walkers at a junction,
// a block's warmth, the stores open or shut, a rule in force), and the complex itself rewritten (a street cut, a
// passage opened through a block, two blocks joined), each edit checked and kept as a history whose order the
// engine gives; and the dynamics — the landscape run on that complex from that state, from the clock's hour, drawn
// on the map as the clock plays, with its record, what could happen at each hour and what did.

import { clearFreedom, edit, emit, ensureDistrict, ensureMassing, frameAt, on, resetEdits, runFrom, setFreedom, state, undoEdit } from '../shared.js';
import { complexIds, describeCell, paintComplex, recolourBuildings, summary } from '../layers/cells.js';

const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.round((x % 1) * 60) % 60).padStart(2, '0')}`;
const MODES = {
  look: { name: 'look', help: 'Pick any cell on the map: its faces and what it bounds light up, and what it carries can be set below.' },
  cut: { name: 'cut a street', help: 'Pick a segment where it should be cut: a junction is made there. The same street, described finer — a refinement.' },
  passage: { name: 'open a passage', help: 'Pick a block, then two points on its edge (two of its streets): a passage is opened between them, and the block becomes two — an intervention.' },
  join: { name: 'join two blocks', help: 'Pick a segment between two blocks: it is taken away and the two blocks become one — an intervention.' },
};

export async function render(root, ctx) {
  const { h, svg, busy, failure, fmt, css } = ctx.dom;
  const { lines, legend } = ctx.chart;
  const { gis } = ctx;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The district: its complex, its freedom, its dynamics'),
    h('span', { class: 'lede' }, 'A district of Melbourne as a cell complex; what can be changed on it — the state its cells carry, and the cells themselves; and how it then runs.')));
  const head = h('div', {});
  root.append(head);
  busy(head, 'cutting the district from the city model…');
  try { await ensureDistrict(ctx); } catch (e) { return failure(head, e); }

  // The complex on the map, drawn here when its layer is off; picks routed to the workspace either way.
  const ids = complexIds((s) => `landscape-dynamics/tool-complex${s}`);
  let alive = true;
  const paint = () => { if (alive && !state.active.has('cells') && state.district) paintComplex(ctx, ids, state.district); };
  ctx.onDispose(() => { alive = false; for (const id of Object.values(ids)) gis.remove(id); });

  // ---- 1. the complex ----
  const stats = h('div', {});
  const viewKeys = h('div', { class: 'keys' });
  const pickedBox = h('div', {});
  const one = h('div', { class: 'section' }, h('h2', {}, 'The complex'), stats, viewKeys,
    h('p', { class: 'caption' }, 'On the ground, each block is drawn a little inside its streets so that the cells stand apart; a segment is darker the more blocks it bounds, and its arrow runs from its tail to its head. The buildings are 3-cells, each bounded by its ground, its roof and its walls; where two stand against each other the wall between them is one 2-cell they share, drawn in sand, and buildings so joined make a terrace. By dimension, the junctions stay on the ground, the segments rise above them, the blocks above those and the buildings highest, each cell joined to its faces — sea where the incidence is +1, sand where it is −1 — and each building to those it shares a wall with. Holes are the faces with a lane ending inside: no single cell can be such a face.'),
    pickedBox);
  const drawStats = () => {
    const c = state.district, t = c.topology, [V, E, F] = t.counts;
    stats.replaceChildren(h('div', { class: 'stat-row' },
      h('div', { class: 'count' }, h('b', {}, String(V)), h('span', {}, 'junctions, 0-cells')),
      h('div', { class: 'count' }, h('b', {}, String(E)), h('span', {}, 'segments, 1-cells')),
      h('div', { class: 'count' }, h('b', {}, String(F)), h('span', {}, 'blocks, 2-cells')),
      h('div', { class: 'count' }, h('b', {}, String(t.holes)), h('span', {}, 'holes')),
      state.massing ? h('div', { class: 'count' }, h('b', {}, String(state.massing.buildings.length)), h('span', {}, 'buildings, 3-cells')) : null,
      state.massing ? h('div', { class: 'count' }, h('b', {}, String(state.massing.walls.length)), h('span', {}, 'walls they share')) : null),
    h('p', { class: 'faint' }, summary(h, c).slice(1).map((x) => x.textContent).join(' · ')));
    viewKeys.replaceChildren(...[['ground', 'on the ground'], ['dimension', 'by dimension']].map(([v, name]) =>
      h('button', { class: `key${state.view === v ? ' on' : ''}`, onclick: () => { state.view = v; paint(); drawStats(); } }, name)));
  };

  // ---- 2. freedom ----
  let mode = 'look', passage = null;
  const modeKeys = h('div', { class: 'keys' });
  const modeHelp = h('p', { class: 'note' });
  const history = h('div', {});
  const freedomBox = h('div', {});
  const two = h('div', { class: 'section' }, h('h2', {}, 'The freedom on top of it'),
    h('h3', {}, 'What the cells carry at the start'), freedomBox,
    h('h3', { style: 'margin-top:12px' }, 'The cells themselves'), modeKeys, modeHelp, history);
  const drawModes = () => {
    modeKeys.replaceChildren(...Object.entries(MODES).map(([k, m]) => h('button', { class: `key${k === mode ? ' on' : ''}`, onclick: () => { mode = k; passage = null; drawModes(); } }, m.name)));
    modeHelp.textContent = passage ? (passage.from ? 'Now the second point, on another of its streets.' : `Block #${passage.block}: now a point on one of its streets.`) : MODES[mode].help;
  };
  const drawHistory = (message = '') => {
    const c = state.district;
    const items = state.edits.map((e, k) => {
      const done = c.edits[k];
      const what = e.kind === 'cut' ? `segment #${e.segment} cut` : e.kind === 'passage' ? `a passage through block #${e.block}` : `the blocks either side of segment #${e.segment} joined`;
      return h('li', {}, what, ' ', done ? h('span', { class: `chip ${done.refinement ? 'proved' : 'model'}` }, done.refinement ? 'a refinement' : 'an intervention') : null,
        done ? h('span', { class: 'faint' }, ` · ${done.steps.length} rewrite${done.steps.length === 1 ? '' : 's'}: ${done.steps.map((s) => s.kind).join(', ')}`) : null);
    });
    const steps = c.order && c.order.steps ? c.order.steps.length : 0;
    history.replaceChildren(...[
      message ? h('p', { class: 'error-box' }, message) : null,
      items.length ? h('ol', {}, items) : h('p', { class: 'faint' }, 'The district as mapped: no edit yet.'),
      items.length > 1 ? h('p', { class: 'faint' }, `These ${items.length} edits need ${steps} step${steps === 1 ? '' : 's'}: an edit comes after another only if it rewrites cells the other made or changed; the rest could have been made in any order.`) : null,
      items.length ? h('div', { class: 'keys' }, h('button', { class: 'plain', onclick: async () => { busy(history, 'undoing…'); await undoEdit(ctx); } }, 'undo the last'),
        h('button', { class: 'plain', onclick: async () => { busy(history, 'back to the district as mapped…'); await resetEdits(ctx); } }, 'back to the district as mapped')) : null].filter(Boolean));
  };
  const number = (value, onSet, unit, step = 1) => {
    const input = h('input', { type: 'number', value: value ?? '', step, style: 'width:7em' });
    input.addEventListener('change', () => onSet(input.value === '' ? null : +input.value));
    return h('label', { class: 'pick' }, input, unit);
  };
  const drawFreedom = () => {
    const c = state.district, p = state.picked;
    const set = [];
    for (const [id, n] of Object.entries(state.initial.walkers)) set.push(`${n} walkers at junction #${id}`);
    for (const [id, t] of Object.entries(state.initial.temperature)) set.push(`block #${id} at ${t} °C`);
    for (const [k, t] of Object.entries(state.initial.buildings)) set.push(`${(state.massing && (state.massing.buildings.find((b) => b.key === k) || {}).name) || 'a building'} at ${t} °C`);
    const names = (c.stores || []).map((s) => s.name);
    const toggle = (key, label) => {
      const box = h('input', { type: 'checkbox', checked: state.initial[key] || null });
      box.addEventListener('change', () => setFreedom(key, box.checked));
      return h('label', { class: 'pick', style: 'margin-right:14px' }, box, label);
    };
    let cellControl = h('p', { class: 'faint' }, 'Pick a junction to set walkers on it, a block or a building to set its warmth.');
    if (p && p.dim === 3) cellControl = h('div', {}, `${p.name || 'The building'}: `, number(state.initial.buildings[p.key], (v) => setFreedom('buildings', p.key, v), '°C at the start (the air\'s if left empty)', 0.5));
    if (p && p.dim === 0) cellControl = h('div', {}, `Junction #${p.id}: `, number(state.initial.walkers[p.id], (v) => setFreedom('walkers', p.id, v), 'walkers'));
    if (p && p.dim === 2) cellControl = h('div', {}, `Block #${p.id}: `, number(state.initial.temperature[p.id] ?? p.temperature, (v) => setFreedom('temperature', p.id, v), '°C', 0.5));
    if (p && p.dim === 1) cellControl = h('p', { class: 'faint' }, `Segment #${p.id} carries the flow of walkers, which follows from the junctions at its ends: set walkers there.`);
    freedomBox.replaceChildren(h('div', { class: 'panel' }, cellControl,
      h('div', { style: 'margin-top:8px' }, names[0] ? toggle('openA', `${names[0]} open`) : null, names[1] ? toggle('openB', `${names[1]} open`) : null, toggle('permit', 'a rule in force (a permit for the lane)')),
      set.length ? h('div', { class: 'faint', style: 'margin-top:6px' }, `Set: ${set.join('; ')}. `, h('button', { class: 'plain', onclick: () => clearFreedom() }, 'clear')) : null));
  };
  const drawPicked = () => {
    const p = state.picked;
    if (!p) { pickedBox.replaceChildren(h('p', { class: 'faint' }, 'Pick a cell on the map.')); return; }
    const d = describeCell({ feature: p }, state.district);
    const orient = p.dim === 3 ? h('div', { class: 'faint' }, `its faces: its ground, its roof and its walls; ${p.neighbours.length ? `${p.neighbours.length} of its walls shared, each a 2-cell with two cofaces` : 'none shared'}`) : p.dim === 2 ? h('div', { class: 'faint' }, `its boundary, counter-clockwise: ${p.boundary.map(([e, sg]) => `${sg > 0 ? '+' : '−'}#${e}`).join(' ')}`)
      : p.dim === 1 ? h('div', { class: 'faint' }, `∂ = #${p.head} − #${p.tail}`) : null;
    pickedBox.replaceChildren(h('div', { class: 'panel' }, h('b', {}, d.title), ...d.lines.map((l) => h('div', {}, l)), orient));
  };

  // ---- 3. dynamics ----
  const P = { span: 6, seed: 3 };
  const result = h('div', {});
  const startLabel = h('span', {});
  const runButton = h('button', { class: 'primary' }, 'run');
  const slider = (key, name, min, max, step, unit = '') => {
    const input = h('input', { type: 'range', min, max, step, value: P[key] });
    const o = h('output', {}, `${fmt(P[key])}${unit}`);
    input.addEventListener('input', () => { o.textContent = `${fmt(+input.value)}${unit}`; });
    input.addEventListener('change', () => { P[key] = +input.value; drawStart(); });
    return h('label', { class: 'slider' }, name, input, o);
  };
  const airSlider = (key, name, min, max, step) => {
    const input = h('input', { type: 'range', min, max, step, value: state.air[key] });
    const o = h('output', {}, `${fmt(state.air[key])} °C`);
    input.addEventListener('input', () => { o.textContent = `${fmt(+input.value)} °C`; });
    input.addEventListener('change', () => { state.air[key] = +input.value; });
    return h('label', { class: 'slider' }, name, input, o);
  };
  const startHour = () => Math.min(23.5, Math.floor(ctx.clock.hour * 4) / 4);
  const drawStart = () => {
    const s = startHour(), e = Math.min(24, s + P.span);
    startLabel.textContent = `from ${hhmm(s)} to ${hhmm(e)} of ${ctx.clock.date}`;
    runButton.textContent = `run from ${hhmm(s)}`;
  };
  const three = h('div', { class: 'section' }, h('h2', {}, 'The dynamics'),
    h('p', { class: 'caption' }, 'The landscape run on this complex from the state set above, from the clock\'s hour: walkers entering at the district\'s edge as Melbourne\'s counts say, walking each way along the streets, the blocks warming with the Sun; the stores opening and shutting; a permit, a lane built under it, and its recognition. The map draws it as the clock plays (the layers "the landscape running" and "the district\'s spacetime").'),
    h('div', { class: 'panel' }, slider('span', 'for', 1, 12, 0.5, ' h'), slider('seed', 'which run', 1, 40, 1), airSlider('mean', 'the day\'s air', -5, 35, 0.5), airSlider('swing', 'warmest above it', 0, 12, 0.5),
      h('div', { class: 'faint' }, startLabel)),
    h('p', {}, runButton), result);
  let shown = false;
  const showOnMap = async () => {
    if (shown || state.active.has('landscape')) return;
    shown = true;
    const mod = await import('../layers/landscape.js');
    const handle = await mod.layer({ ...ctx, layerId: (s = '') => `landscape-dynamics/tool-landscape${s}` });
    ctx.clock.on((c) => { if (alive && handle.setTime) handle.setTime(c); });
    ctx.onDispose(() => { if (handle.dispose) handle.dispose(); for (const id of [...gis.layers.keys()]) if (id.startsWith('landscape-dynamics/tool-landscape')) gis.remove(id); });
  };
  runButton.addEventListener('click', async () => {
    const s = startHour(), e = Math.min(24, s + P.span);
    busy(result, `running ${hhmm(s)}–${hhmm(e)} on the complex as it is, from the state set…`);
    let R;
    try { R = await runFrom(ctx, { start: s, end: e, seed: P.seed }); } catch (x) { return failure(result, x); }
    await showOnMap();
    drawRun(R);
    ctx.clock.set({ hour: s });
    const sl = document.querySelector('.clock input[type=range]');
    if (sl) sl.value = s;
  });
  // The buildings' warmth through the run: their mean, the warmest and coolest, and the building picked.
  const heatChart = () => {
    const Hh = state.heat;
    if (!Hh) return null;
    const times = Hh.frames.map((f) => f.time);
    const stats = Hh.frames.map((f) => { const v = Object.values(f.temperature); return { lo: Math.min(...v), hi: Math.max(...v), mean: v.reduce((a, x) => a + x, 0) / v.length }; });
    const series = [{ name: 'coolest to warmest building', x: times, y: stats.map((x) => x.lo), y2: stats.map((x) => x.hi), color: css('--minor'), opacity: 0.2 },
      { name: 'their mean', x: times, y: stats.map((x) => x.mean), color: css('--minor-deep'), width: 2 }];
    const p = state.picked;
    if (p && p.dim === 3) series.push({ name: p.name || 'the building picked', x: times, y: Hh.frames.map((f) => f.temperature[p.key]), color: css('--accent-deep'), width: 2 });
    return h('div', {}, lines(series, { xlabel: 'hour', ylabel: '°C, the buildings', height: 150 }),
      legend([{ name: 'the buildings, coolest to warmest', color: css('--minor') }, { name: 'their mean', color: css('--minor-deep') }, ...(p && p.dim === 3 ? [{ name: p.name || 'the building picked', color: css('--accent-deep') }] : [])]),
      h('p', { class: 'caption' }, `The ${Hh.buildings.length} buildings' warmth through the run: each a 3-cell exchanging heat with its neighbours through the ${Hh.walls.length} walls they share, with the air through its facades and roof, and with the ground; the Sun on what it sees of the sky; the air declared above.`));
  };
  const drawRun = (R) => {
    const times = R.frames.map((f) => f.time);
    const total = R.frames.map((f) => Object.values(f.walkers).reduce((a, x) => a + x, 0));
    const traffic = R.frames.map((f) => Object.keys(f.flow).reduce((a, id) => a + (f.forth[id] || 0) + (f.back[id] || 0), 0) / Math.max(1, Object.keys(f.flow).length));
    const rows = R.events.slice(0, 40).map((e) => h('tr', {}, h('td', { class: 'mono' }, hhmm(e.time)),
      h('td', { class: 'left' }, e.mechanism, e.rewrites && e.rewrites.length ? h('span', { class: 'chip proved', style: 'margin-left:6px' }, `${e.rewrites.length} rewrites`) : null), h('td', { class: 'left faint' }, e.match)));
    const changes = R.events.filter((e) => e.rewrites && e.rewrites.length).length;
    result.replaceChildren(...[
      lines([{ name: 'walkers in the district', x: times, y: total, color: css('--accent-deep'), width: 2 }], { xlabel: 'hour', ylabel: 'walkers', height: 160 }),
      lines([{ name: 'a street\'s traffic, on average', x: times, y: traffic, color: css('--series-2'), width: 2 }], { xlabel: 'hour', ylabel: 'walkers an hour, both ways', height: 130 }),
      heatChart(),
      ctx.finding({ title: `${R.events.length} events from ${hhmm(R.start)} to ${hhmm(R.end)}; the cells changed ${changes} time${changes === 1 ? '' : 's'}`, kind: 'model',
        statement: `Walkers are kept exactly: ${Math.round(R.balance.end).toLocaleString('en')} at the end = ${Math.round(R.balance.start).toLocaleString('en')} at the start + ${Math.round(R.balance.entered).toLocaleString('en')} entered − ${Math.round(R.balance.left).toLocaleString('en')} left (residual ${fmt(R.balance.residual, 2)}), every rewrite carrying them across. Every configuration the run passed through was a regular complex.`,
        method: `A piecewise-deterministic process on the strata of the complex, from the state set; walkers entering: ${R.entriesSource}.`,
        limits: 'Every rate is declared, not measured: a model to explore how the cells, what they carry and how they are read change one another.' }),
      h('h3', {}, 'What could happen, and what did'), registry(h, svg, css, R),
      h('table', { class: 'data', style: 'width:100%' }, h('tbody', {}, rows)), R.events.length > 40 ? h('p', { class: 'faint' }, `…and ${R.events.length - 40} more.`) : null,
      h('p', { class: 'note' }, 'Press the clock\'s ▶ to watch it on the map; turn on "the district\'s spacetime" to see the whole run at once, time rising above the district.')].filter(Boolean));
  };

  head.replaceChildren(one, two, three);
  // The buildings stand on the blocks; as glass (the rail's own key, so that its state shows) the cells show through.
  const glass = [...document.querySelectorAll('[data-role="building-form"] button')].find((b) => b.textContent === 'glass');
  if (glass && !glass.classList.contains('on')) glass.click();
  // The district framed in the part of the map the drawer leaves free: the view's target moved to its right.
  {
    const B = state.box, cx = (B.x0 + B.x1) / 2, cy = (B.y0 + B.y1) / 2, size = Math.max(B.x1 - B.x0, B.y1 - B.y0);
    const o = gis.camera.position.clone().sub(gis.controls.target), n = Math.hypot(o.x, o.z) || 1;
    gis.flyTo([cx + (o.z / n) * 0.32 * size, cy + (o.x / n) * 0.32 * size], 1.5 * size);
  }
  const refresh = () => { drawStats(); drawPicked(); drawFreedom(); drawHistory(); drawModes(); drawStart(); paint(); };
  refresh();
  if (state.run && state.run.start !== undefined) drawRun(state.run);
  const off = on((what) => {
    if (!alive) return;
    if (what === 'district') { drawStats(); drawPicked(); drawFreedom(); drawHistory(); paint(); }
    if (what === 'freedom' || what === 'picked') { drawFreedom(); drawPicked(); paint(); }
    if (what === 'view') { drawStats(); paint(); }
  });
  ctx.onDispose(off);
  ctx.clock.on((c) => { if (alive) { drawStart(); recolourBuildings(ctx, ids, c.hour); } });
  ensureMassing(ctx).then(() => { if (alive) { drawStats(); paint(); } }).catch((x) => console.error(x));

  // Picks on the complex, wherever it is drawn.
  gis.onPick(async (hit) => {
    if (!alive || !hit || !hit.layer || !/^landscape-dynamics\/(cells|tool-complex):(blocks|edges|vertices|buildings)$/.test(hit.layer)) return;
    const cell = hit.feature, at = hit.point ? [hit.point[0], hit.point[1]] : null;
    if (mode === 'look') { state.picked = cell; emit('picked'); if (state.run && state.run.start !== undefined) drawRun(state.run); return; }
    const attempt = async (e) => {
      busy(history, 'rewriting the cells…');
      let r;
      try { r = await edit(ctx, e); } catch (x) { return drawHistory(String(x.message || x)); }
      if (r.refused) drawHistory(`Not done: ${r.refused}.`);
      passage = null;
      drawModes();
    };
    if (mode === 'cut') { if (cell.dim !== 1) return drawHistory('A street is cut: pick a segment.'); return attempt({ kind: 'cut', segment: cell.id, at }); }
    if (mode === 'join') { if (cell.dim !== 1) return drawHistory('Blocks are joined across a segment: pick the segment between them.'); return attempt({ kind: 'join', segment: cell.id }); }
    if (mode === 'passage') {
      if (!passage) {
        if (cell.dim !== 2) return drawHistory('First pick the block the passage goes through.');
        passage = { block: cell.id };
        state.picked = cell; emit('picked');
        return drawModes();
      }
      if (!passage.from) { passage.from = at; return drawModes(); }
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
  const W = 600, H = 34 + 26 * kinds.length, x = (t) => 70 + ((t - t0) / (t1 - t0 || 1)) * (W - 90), y = (e) => 22 + 26 * Math.max(0, kinds.indexOf(e.mechanism));
  const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%' });
  kinds.forEach((k, i) => s.append(svg('text', { x: 4, y: 26 + 26 * i, 'font-size': 11, fill: css('--ink-2') }, label[i])));
  for (const F of R.frames)
    for (const k of F.enabled || []) {
      const i = kinds.indexOf(k);
      if (i >= 0) s.append(svg('rect', { x: x(F.time), y: 22 + 26 * i - 7, width: Math.max(0.5, x(F.time + 0.25) - x(F.time)), height: 14, fill: css('--accent-soft') }));
    }
  for (let t = Math.ceil(t0); t <= t1; t += Math.max(1, Math.round((t1 - t0) / 6))) s.append(svg('text', { x: x(t), y: H - 4, 'font-size': 10, 'text-anchor': 'middle', fill: css('--ink-faint') }, `${t}:00`));
  for (const [a, b] of R.order.covers || []) {
    const A = R.events[a], B = R.events[b];
    s.append(svg('line', { x1: x(A.time), y1: y(A), x2: x(B.time), y2: y(B), stroke: css('--accent'), 'stroke-opacity': 0.35, 'stroke-width': 1 }));
  }
  R.events.forEach((e) => s.append(svg('circle', { cx: x(e.time), cy: y(e), r: 3.6, fill: e.forced ? css('--minor') : css('--accent-deep') })));
  return h('div', { class: 'panel', style: 'margin:8px 0' }, s, h('div', { class: 'faint' }, 'each kind of event on its row: shaded while the state allowed it, a dot when it happened; a line joins an event to the next ones that read what it wrote'));
}
