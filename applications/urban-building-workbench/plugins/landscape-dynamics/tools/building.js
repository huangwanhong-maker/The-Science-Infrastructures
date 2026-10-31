// A building's rooms: a building picked on the map filled with a box of rooms fitted to its footprint and floors —
// each room a 3-cell, each wall, floor, facade and roof a 2-cell of one regular complex — and its heat through the
// clock's day under the real Sun and a declared day's air; how long it takes to forget a disturbance; and the
// events that change it: doors opened in its walls and stairs in its floors (the same cells, carrying another
// meaning), two rooms merged into one and parted again (new cells), each with the heat it keeps and the rooms one
// can reach from the street. The map shows the rooms at the clock's hour (the layer "a building's rooms").

import { runBuilding, specOf, state } from '../shared.js';

const power = (w, fmt) => (w >= 1e6 ? `${fmt(w / 1e6, 3)} MW` : `${fmt(w / 1e3, 3)} kW`);

const KINDS = [
  { kind: 'doors', name: 'doors in every inside wall', what: 'every wall between two rooms relabelled a door: the same faces, conducting as doors do' },
  { kind: 'stairs', name: 'a stair between every two floors', what: 'a floor slab in one corner of each storey relabelled a stair' },
  { kind: 'merge', name: 'two ground-floor rooms merged', what: 'the wall between two rooms removed and the two rooms made one cell' },
  { kind: 'partition', name: 'and parted again', what: 'the merged room parted by a wall where the old one stood' },
];

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'A building\'s rooms'),
    h('span', { class: 'lede' }, 'A building as a complex of rooms: its heat through the day under the Sun, how long it remembers a disturbance, and the openings and walls that change how heat and people pass through it.')));
  const pickNote = h('p', { class: 'note' }, 'Pick a building on the map.');
  const plan = { doors: { on: true, hour: 8 }, stairs: { on: true, hour: 8 }, merge: { on: true, hour: 12 }, partition: { on: true, hour: 17 } };
  const air = { mean: 18, swing: 5 };
  const panel = h('div', { class: 'panel' });
  const slider = (o, key, name, min, max, step, unit) => {
    const input = h('input', { type: 'range', min, max, step, value: o[key] });
    const out = h('output', {}, `${fmt(o[key])}${unit}`);
    input.addEventListener('input', () => { out.textContent = `${fmt(+input.value)}${unit}`; });
    input.addEventListener('change', () => { o[key] = +input.value; go(); });
    return h('label', { class: 'slider' }, name, input, out);
  };
  panel.append(h('h3', {}, 'What happens to it, and when'),
    ...KINDS.map((k) => {
      const box = h('input', { type: 'checkbox', checked: plan[k.kind].on || null });
      box.addEventListener('change', () => { plan[k.kind].on = box.checked; go(); });
      return h('div', { style: 'display:grid;grid-template-columns:1.6em 1fr;align-items:center' }, box, slider(plan[k.kind], 'hour', k.name, 0, 23, 1, ':00'));
    }),
    h('h3', { style: 'margin-top:8px' }, 'The day\'s air (declared)'), slider(air, 'mean', 'its mean', -5, 35, 0.5, ' °C'), slider(air, 'swing', 'warmest above it', 0, 12, 0.5, ' °C'));
  const out = h('div', {});
  root.append(pickNote, panel, out);

  let picked = null;
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
    pickNote.textContent = `${spec.name}: about ${Math.round(spec.width)} × ${Math.round(spec.depth)} m on ${spec.floors} floor${spec.floors === 1 ? '' : 's'}, filled with ${spec.nx * spec.ny * spec.floors} rooms of about ${fmt(spec.width / spec.nx, 2)} × ${fmt(spec.depth / spec.ny, 2)} m. Pick another to compare.`;
    const events = KINDS.filter((k) => plan[k.kind].on && (k.kind !== 'partition' || plan.merge.on)).map((k) => ({ hour: plan[k.kind].hour, kind: k.kind }));
    const airByHour = Array.from({ length: 24 }, (_, hr) => air.mean + air.swing * Math.cos((2 * Math.PI * (hr + 0.5 - 15)) / 24));
    busy(out, 'filling the building with rooms and following its heat through the day…');
    let B;
    try { B = await runBuilding(ctx, { ...spec, air: airByHour[0], airByHour }, events); } catch (e) { return failure(out, e); }
    await showOnMap();
    draw(B, airByHour);
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

  function draw(B, airByHour) {
    const hoursOf = B.hours.map((x) => x.hour + 0.5);
    const stats = B.hours.map((x) => {
      const t = x.rooms.map((r) => r.temperature);
      return { lo: Math.min(...t), hi: Math.max(...t), mean: t.reduce((a, v) => a + v, 0) / t.length };
    });
    const [V, E, F, C] = B.counts;
    const slowest = B.timeConstants[0];
    const marks = B.events.map((e, k) => ({ x: e.hour, label: KINDS.find((x) => x.kind === e.kind).name.split(' ').slice(0, 2).join(' '), row: k % 3 }));
    const evRows = B.events.map((e) => {
      const k = KINDS.find((x) => x.kind === e.kind);
      const refines = e.steps.every((s) => s.kind === 'relabel');
      const kept = Math.abs(e.heatChange) <= 1e-12 * Math.abs(B.hours[0].heat);
      return h('tr', {}, h('td', { class: 'mono' }, `${String(e.hour).padStart(2, '0')}:00`),
        h('td', { class: 'left' }, k.name, ' ', h('span', { class: `chip ${refines ? 'proved' : 'model'}` }, refines ? 'the same cells, another meaning' : 'new cells'), h('div', { class: 'faint' }, k.what)),
        h('td', {}, String(e.rooms)), h('td', {}, String(e.reachable)), h('td', {}, kept ? 'kept' : `${fmt(e.heatChange, 3)} J`));
    });
    out.replaceChildren(
      h('div', { class: 'stat-row' },
        h('div', { class: 'count' }, h('b', {}, String(C)), h('span', {}, 'rooms (3-cells)')),
        h('div', { class: 'count' }, h('b', {}, String(F)), h('span', {}, 'walls, floors, facades')),
        h('div', { class: 'count' }, h('b', {}, `${fmt(slowest, 3)} h`), h('span', {}, 'to forget a disturbance')),
        h('div', { class: 'count' }, h('b', {}, power(Math.max(...B.hours.map((x) => x.sun)), fmt)), h('span', {}, 'from the Sun at most'))),
      ctx.finding({ title: `${C} rooms, ${F} faces, ${E} edges and ${V} corners make one regular complex`, kind: 'proved',
        statement: `Each room is a box bounded by its walls, floor and ceiling; each wall is shared by at most two rooms. ${B.betti.join(',') === '1,0,0,0' ? 'The complex is a ball (Betti numbers 1, 0, 0, 0): no hole through it and no enclosed void.' : `Its Betti numbers are ${B.betti.join(', ')}.`} Heat lives in the rooms and passes through their shared faces, so the building's conduction is the Hodge Laplacian of its rooms, its faces' conductances and its rooms' capacities the stars; the slowest of its modes decays in ${fmt(slowest, 3)} hours, the next in ${fmt(B.timeConstants[1] ?? slowest, 3)}.`,
        method: 'A cubical complex of the rooms fitted to the footprint\'s smallest rectangle (morphe); conductances of walls, slabs, facades, roof and ground by their transmittance and area, capacities by the rooms\' volume; the heat advanced exactly over each hour (physis\'s Hodge Laplacian, the matrix exponential).',
        limits: 'The rooms are a regular grid fitted to the building, not its plan; transmittances are typical values, not this building\'s.' }),
      h('h3', {}, 'Its rooms through the day'),
      lines([{ name: 'coolest to warmest room', x: hoursOf, y: stats.map((s) => s.lo), y2: stats.map((s) => s.hi), color: css('--minor'), opacity: 0.22 },
        { name: 'mean of the rooms', x: hoursOf, y: stats.map((s) => s.mean), color: css('--minor-deep'), width: 2 },
        { name: 'the air outside (declared)', x: hoursOf, y: airByHour, color: css('--accent'), dash: '5 4' }],
      { xlabel: 'hour', ylabel: '°C', height: 200, xlim: [0, 24], marks }),
      legend([{ name: 'the rooms, coolest to warmest', color: css('--minor') }, { name: 'their mean', color: css('--minor-deep') }, { name: 'the air (declared)', color: css('--accent') }]),
      lines([{ name: 'the Sun\'s heat reaching its rooms', x: hoursOf, y: B.hours.map((x) => x.sun / 1000), color: css('--series-4'), width: 2 }], { xlabel: 'hour', ylabel: 'kW reaching the rooms', height: 140, xlim: [0, 24] }),
      h('p', { class: 'caption' }, `The Sun of ${B.date} at this place (the astronomy domain), on each facade by its orientation and on the roof: what the glass lets in, and the small share of what walls and roof absorb that reaches the rooms rather than the air outside. The clock's slider moves the hour drawn on the map.`),
      h('h3', {}, 'How long it remembers'),
      bars(B.timeConstants.map((t, k) => ({ x: k, y: t, color: k ? css('--series-3') : css('--accent-deep'), label: `mode ${k + 1}: ${fmt(t, 3)} h` })), { labels: B.timeConstants.map((_, k) => `${k + 1}`), ylabel: 'hours to fall to 1/e, the slowest modes first', height: 160 }),
      h('p', { class: 'caption' }, 'A disturbance — a sunny hour, a door left open — is a sum of the building\'s modes, each fading at its own rate; the slowest is the building as a whole exchanging heat with the air, the faster ones its rooms evening out among themselves.'),
      B.events.length ? h('h3', {}, 'What happened to it') : null,
      B.events.length ? h('table', { class: 'data', style: 'width:100%' }, h('thead', {}, h('tr', {}, h('th', { class: 'left' }, 'at'), h('th', { class: 'left' }, 'the event'),
        h('th', {}, 'rooms'), h('th', {}, 'reached from the street'), h('th', {}, 'heat changed'))), h('tbody', {}, evRows)) : null,
      B.events.length ? ctx.finding({ title: 'Openings change how heat and people pass, never how much heat there is', kind: 'proved',
        statement: `A door or a stair is a wall or a floor read otherwise: the same cell, carrying another conductance, and a way through for people — ${B.events.find((e) => e.kind === 'doors') ? `with doors, ${B.events.find((e) => e.kind === 'doors').reachable} rooms can be reached from the entrance` : 'the rooms reached from the entrance change with them'}. Merging two rooms removes a wall and makes their heat one room's, at the mean temperature their capacities weigh; parting them again copies that temperature to both. Each event keeps the building's heat to rounding and its complex regular.`,
        method: 'Relabelling of faces; merging and parting of rooms by rewriting the complex (morphe); heat carried as an extensive quantity, temperatures as intensive ones; rooms reached by a search through doors and stairs.' }) : null);
  }

  ctx.onSelect((sel) => { if (choose(sel)) go(); });
  if (choose(ctx.selection)) go();
  else if (state.building && ctx.model) {
    const b = (ctx.model.buildings || []).find((x) => x.id === state.building.id);
    if (b) { picked = b; go(); }
  }
}
