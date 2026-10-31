// A building's rooms: a building picked on the map filled with a box of rooms fitted to its footprint and floors —
// each room a 3-cell, each wall, floor, facade and roof a 2-cell of one regular complex — and its heat through the
// clock's day under the real Sun and a declared day's air; how long it takes to forget a disturbance; and the
// events that change it: doors opened in its walls and stairs in its floors (the same cells, carrying another
// meaning), two rooms merged into one and parted again (new cells), each with the heat it keeps and the rooms one
// can reach from the street. The map shows the rooms at the clock's hour (the layer "a building's rooms").

import { forArea, runBuilding, specOf, state } from '../shared.js';
import { brief, figures, fold, hint, slider, swatches, tabs, toggle } from '../ui.js';

const power = (w, fmt) => (w >= 1e6 ? `${fmt(w / 1e6, 3)} MW` : `${fmt(w / 1e3, 3)} kW`);

const KINDS = [
  { kind: 'doors', name: 'Doors in the inside walls', short: 'doors', what: 'the same walls, read as doors' },
  { kind: 'stairs', name: 'A stair between floors', short: 'stairs', what: 'a corner of each slab, read as a stair' },
  { kind: 'merge', name: 'Two rooms merged', short: 'merged', what: 'a wall removed, two rooms made one' },
  { kind: 'partition', name: 'And parted again', short: 'parted', what: 'a wall put back where it stood' },
];
const hh = (x) => `${String(x).padStart(2, '0')}:00`;

export async function render(root, ctx) {
  forArea(ctx);
  const { h, failure, fmt, css } = ctx.dom;
  const { lines, bars } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'A building\'s rooms'),
    h('span', { class: 'lede' }, 'Its rooms as one complex, their heat under the Sun, and what changes them.')));
  const plan = { doors: { on: true, hour: 8 }, stairs: { on: true, hour: 8 }, merge: { on: true, hour: 12 }, partition: { on: true, hour: 17 } };
  const air = state.air;
  const name = h('div', { class: 'ld-title' }, 'Pick a building on the map.');
  const sub = h('div', { class: 'ld-sub' }, 'Its rooms are fitted to its footprint and floors.');
  const out = h('div', {});
  root.append(h('div', { class: 'ld-card' }, name, sub), out);

  let picked = null, B = null, airByHour = null;
  const choose = (sel) => {
    if (!sel || sel.kind !== 'building' || !ctx.model) return false;
    const id = (sel.lines || [])[1];
    const b = (ctx.model.buildings || []).find((x) => x.id === id);
    if (!b || !b.footprint || b.footprint.length < 3) return false;
    picked = b;
    return true;
  };
  async function go() {
    if (!picked) return;
    const spec = specOf(picked);
    name.textContent = spec.name.charAt(0).toUpperCase() + spec.name.slice(1);
    sub.textContent = `about ${Math.round(spec.width)} × ${Math.round(spec.depth)} m · ${spec.floors} floor${spec.floors === 1 ? '' : 's'} · ${spec.nx * spec.ny * spec.floors} rooms`;
    const events = KINDS.filter((k) => plan[k.kind].on && (k.kind !== 'partition' || plan.merge.on)).map((k) => ({ hour: plan[k.kind].hour, kind: k.kind }));
    airByHour = Array.from({ length: 24 }, (_, hr) => air.mean + air.swing * Math.cos((2 * Math.PI * (hr + 0.5 - 15)) / 24));
    if (!B) out.replaceChildren(hint(h, 'filling the building with rooms and following its heat through the day…'));
    try { B = await runBuilding(ctx, { ...spec, air: airByHour[0], airByHour }, events); } catch (e) { return failure(out, e); }
    await showOnMap();
    const on = out.querySelector('.ld-tab.on');
    out.replaceChildren(T.el);
    T.show(on ? on.dataset.name : undefined);
  }

  // The map: the layer's drawing when it is off.
  let shown = false;
  const showOnMap = async () => {
    if (state.active.has('rooms') || shown) return;
    shown = true;
    const mod = await import('../layers/rooms.js');
    const handle = await mod.layer({ ...ctx, layerId: (s = '') => `landscape-dynamics/tool-rooms${s}` });
    let alive = true;
    ctx.clock.on((c) => { if (alive && handle.setTime) handle.setTime(c); });
    ctx.onDispose(() => {
      alive = false;
      if (handle.dispose) handle.dispose();
      for (const id of [...ctx.gis.layers.keys()]) if (id.startsWith('landscape-dynamics/tool-rooms')) ctx.gis.remove(id);
    });
  };

  // ---- the tabs ----
  const roomsPanel = () => {
    const [V, E, F, C] = B.counts;
    const ball = B.betti.join(',') === '1,0,0,0';
    return h('div', {},
      figures(h, [{ value: C, label: 'rooms' }, { value: F, label: 'walls, floors, facades' }, { value: E, label: 'edges' }, { value: V, label: 'corners' }]),
      brief(ctx, { title: ball ? 'One complex, a ball' : `One complex, Betti numbers ${B.betti.join(', ')}`, kind: 'proved',
        statement: `Each room is a 3-cell bounded by its walls, floor and ceiling; each wall is shared by at most two rooms.${ball ? ' No hole through it, no enclosed void.' : ''}`,
        method: 'A cubical complex of rooms fitted to the footprint\'s smallest rectangle; its regularity and homology checked (morphe).',
        limits: 'The rooms are a regular grid fitted to the building, not its plan.' }),
      fold(h, 'what the map shows',
        swatches(h, [{ color: css('--minor'), label: 'cooler rooms' }, { color: css('--accent-deep'), label: 'warmer rooms' }]),
        hint(h, 'The rooms at the clock\'s hour; move the clock to follow them through the day.')));
  };

  const heatPanel = () => {
    const hoursOf = B.hours.map((x) => x.hour + 0.5);
    const stats = B.hours.map((x) => {
      const t = x.rooms.map((r) => r.temperature);
      return { lo: Math.min(...t), hi: Math.max(...t), mean: t.reduce((a, v) => a + v, 0) / t.length };
    });
    const slowest = B.timeConstants[0];
    const marks = B.events.map((e, k) => ({ x: e.hour, label: KINDS.find((x) => x.kind === e.kind).short, row: k % 3 }));
    const warmest = stats.reduce((m, s, k) => (s.mean > stats[m].mean ? k : m), 0);
    return h('div', {},
      figures(h, [{ value: `${fmt(slowest, 3)} h`, label: 'to forget a disturbance' }, { value: `${fmt(stats[warmest].mean, 3)} °C`, label: `warmest, near ${hh(B.hours[warmest].hour)}` },
        { value: power(Math.max(...B.hours.map((x) => x.sun)), fmt), label: 'from the Sun at most' }]),
      lines([{ name: 'coolest to warmest room', x: hoursOf, y: stats.map((s) => s.lo), y2: stats.map((s) => s.hi), color: css('--minor'), opacity: 0.22 },
        { name: 'mean of the rooms', x: hoursOf, y: stats.map((s) => s.mean), color: css('--minor-deep'), width: 2 },
        { name: 'the air outside (declared)', x: hoursOf, y: airByHour, color: css('--accent'), dash: '5 4' }],
      { xlabel: 'hour', ylabel: '°C', height: 180, xlim: [0, 24], marks }),
      swatches(h, [{ color: css('--minor'), label: 'rooms, coolest to warmest' }, { color: css('--minor-deep'), label: 'their mean', kind: 'line' }, { color: css('--accent'), label: 'the air', kind: 'line' }]),
      brief(ctx, { title: `It forgets in ${fmt(slowest, 2)} hours`, kind: 'proved',
        statement: 'Heat lives in the rooms and passes through their shared faces. A disturbance is a sum of the building\'s modes: the slowest is the whole building with the air, the faster ones its rooms evening out.',
        method: 'The Hodge Laplacian of the rooms, faces\' conductances by transmittance and area, rooms\' capacities by volume; advanced exactly over each hour (physis). The Sun of the date at this place on each facade and the roof (the astronomy domain).',
        limits: 'Transmittances are typical values, not this building\'s; the air is declared.' }),
      fold(h, 'the day\'s air',
        slider(h, fmt, { label: 'its mean', min: -5, max: 35, step: 0.5, value: air.mean, unit: ' °C', onChange: (v) => { air.mean = v; go(); } }),
        slider(h, fmt, { label: 'warmest above it', min: 0, max: 12, step: 0.5, value: air.swing, unit: ' °C', onChange: (v) => { air.swing = v; go(); } })),
      fold(h, 'the Sun, and every mode',
        lines([{ name: 'the Sun\'s heat reaching its rooms', x: hoursOf, y: B.hours.map((x) => x.sun / 1000), color: css('--series-4'), width: 2 }], { xlabel: 'hour', ylabel: 'kW reaching the rooms', height: 120, xlim: [0, 24] }),
        bars(B.timeConstants.map((t, k) => ({ x: k, y: t, color: k ? css('--series-3') : css('--accent-deep'), label: `mode ${k + 1}: ${fmt(t, 3)} h` })), { labels: B.timeConstants.map((_, k) => `${k + 1}`), ylabel: 'hours to 1/e', height: 130 })));
  };

  const changesPanel = () => {
    const rows = KINDS.map((k) => {
      const e = B.events.find((x) => x.kind === k.kind);
      const disabled = k.kind === 'partition' && !plan.merge.on;
      const head = h('div', { class: 'ld-row' },
        toggle(h, k.name, plan[k.kind].on && !disabled, (v) => { plan[k.kind].on = v; go(); }),
        plan[k.kind].on && !disabled ? slider(h, fmt, { label: '', min: 0, max: 23, step: 1, value: plan[k.kind].hour, unit: ':00', onChange: (v) => { plan[k.kind].hour = v; go(); } }) : null);
      const kept = e && Math.abs(e.heatChange) <= 1e-12 * Math.abs(B.hours[0].heat);
      const refines = e && e.steps.every((s) => s.kind === 'relabel');
      const what = e ? h('div', { class: 'ld-quiet' }, `${k.what} · ${refines ? 'the same cells' : 'new cells'} · ${e.rooms} rooms, ${e.reachable} reached · heat ${kept ? 'kept' : `${fmt(e.heatChange, 3)} J`}`)
        : h('div', { class: 'ld-quiet' }, disabled ? 'only after a merge' : k.what);
      return h('li', {}, head, what);
    });
    const doors = B.events.find((e) => e.kind === 'doors');
    return h('div', {},
      h('ul', { class: 'ld-list stack' }, rows),
      B.events.length ? brief(ctx, { title: 'Openings change how heat passes, not how much', kind: 'proved',
        statement: `A door or a stair is a wall or a floor read otherwise: the same cell, another conductance, a way through${doors ? ` — with doors, ${doors.reachable} rooms are reached from the entrance` : ''}. Merging makes two rooms' heat one; parting copies its temperature to both. Heat is kept to rounding, the complex stays regular.`,
        method: 'Relabelling of faces; merging and parting by rewriting the complex (morphe); heat extensive, temperature intensive; rooms reached by a search through doors and stairs.' }) : hint(h, 'Nothing happens to it today: switch a change on.'));
  };

  const T = tabs(h, [{ name: 'Rooms', panel: () => roomsPanel() }, { name: 'Heat', panel: () => heatPanel() }, { name: 'Changes', panel: () => changesPanel() }], { remember: 'landscape-building', lazy: true });

  ctx.onSelect((sel) => { if (choose(sel)) go(); });
  if (choose(ctx.selection)) go();
  else if (state.building && ctx.model) {
    const b = (ctx.model.buildings || []).find((x) => x.id === state.building.id);
    if (b) { picked = b; go(); }
  }
}
