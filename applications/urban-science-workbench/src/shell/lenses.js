// What can be looked at, by the question it answers or by the work that holds it: two orders over one list, chosen by a key at their head. Each opens in place and shows
// its layers (each with a small copy of the altitude tape, marking the scales at which it means something, and dimmed when the camera is elsewhere) and its analyses. One is open at
// a time. What is of another place is not hidden: it says where it is, and asking for it goes there.

import { h, svg } from '../dom.js';
import { along } from './navigator.js';

const GLYPHS = {
  // People: cells of a grid rising.
  people: (g) => { g.append(svg('path', { d: 'M3.5 14l4-2.3 4 2.3v4.6l-4 2.3-4-2.3zM12.5 14l4-2.3 4 2.3v4.6l-4 2.3-4-2.3zM8 5.6l4-2.3 4 2.3v4.6l-4 2.3-4-2.3z' })); },
  // Movement: currents two ways.
  movement: (g) => { g.append(svg('path', { d: 'M3 8.5c3-2.4 5-2.4 8 0s5 2.4 8 0M16.5 5.6l3.2 2.5-2.6 3' }), svg('path', { d: 'M21 15.5c-3 2.4-5 2.4-8 0s-5-2.4-8 0M7.5 18.4l-3.2-2.5 2.6-3' })); },
  // Fabric: a building in axonometric.
  fabric: (g) => { g.append(svg('path', { d: 'M12 3.2l7.5 4.3v9L12 20.8l-7.5-4.3v-9zM12 11.8l7.5-4.3M12 11.8l-7.5-4.3M12 11.8v9' })); },
  // Risk: water rising against a mark.
  risk: (g) => { g.append(svg('path', { d: 'M12 3.4l8.4 14.4H3.6z' }), svg('path', { d: 'M2.5 20.6c2.4-1.6 4.8-1.6 7.2 0s4.8 1.6 7.2 0 3-1.5 4.6-.6M12 8.8v4.6M12 15.2v.6' })); },
  // Change: rings of years.
  change: (g) => { g.append(svg('path', { d: 'M3.4 12a8.6 8.6 0 1 1 3.2 6.7' }), svg('path', { d: 'M7.4 12a4.6 4.6 0 1 1 1.9 3.7' }), svg('path', { d: 'M3 15.8l3.6 2.9.9-4.1M12 12h.2' })); },
  // Yours: a stone dropped in still water, its rings widening.
  yours: (g) => { g.append(svg('circle', { cx: 12, cy: 12, r: 1.6 }), svg('circle', { cx: 12, cy: 12, r: 5.4 }), svg('path', { d: 'M12 2.8a9.2 9.2 0 0 1 9.2 9.2M12 21.2A9.2 9.2 0 0 1 2.8 12' })); },
};

const remembered = (key, fallback) => { try { return localStorage.getItem(key) ?? fallback; } catch { return fallback; } };
const remember = (key, value) => { try { localStorage.setItem(key, value); } catch { /* a private window */ } };

export function createLenses(el, { lensList, onToggle, onTool, plannedPanel }) {
  let open = null, mode = remembered('urban-lens-order', 'questions'), data = [], planned = [], tool = null, works = [];
  const rows = new Map();   // key → { row, input, marker }
  const counts = new Map(); // id → the count element

  function miniScale(near, far) {
    const a = along(far) * 100, b = along(near) * 100;
    const marker = h('b', {});
    return { el: h('span', { class: 'scale-mini', title: 'the scales at which this means something' }, h('i', { style: `left:${a}%;width:${Math.max(4, b - a)}%` }), marker), marker };
  }

  function bodyOf(items, { byGroup = true } = {}) {
    const layers = items.filter((d) => d.kind === 'layer'), tools = items.filter((d) => d.kind === 'tool');
    const body = h('div', { class: 'lens-body' });
    let group = null;
    for (const d of layers) {
      if (byGroup && d.group && d.group !== group) { group = d.group; body.append(h('div', { class: 'lens-group' }, group)); }
      const input = h('input', { type: 'checkbox', checked: d.isOn() ? true : null, onchange: () => { onToggle(d, input.checked); refresh(); } });
      const mini = miniScale(d.scale[0], d.scale[1]);
      const row = h('label', { class: `layer-row${d.here ? '' : ' elsewhere'}`, title: d.about || '' }, input, h('span', { class: 'box' }), h('span', { class: 'name' }, d.title, d.here ? null : h('em', { class: 'where' }, ` in ${d.places}`)), mini.el);
      rows.set(d.key, { row, input, marker: mini.marker, d });
      body.append(row);
    }
    if (tools.length) body.append(h('div', { class: 'lens-group' }, 'Analyses'));
    for (const t of tools) {
      body.append(h('button', { class: `tool-row${tool === t.key ? ' active' : ''}${t.here ? '' : ' elsewhere'}`, title: t.about || '', onclick: () => onTool(t) }, h('span', {}, t.title, t.here ? null : h('em', { class: 'where' }, ` in ${t.places}`))));
    }
    if (!layers.length && !tools.length) body.append(h('div', { class: 'lens-foot' }, 'Nothing here yet.'));
    return body;
  }

  function section(id, title, question, glyphId, items) {
    const glyph = svg('svg', { viewBox: '0 0 24 24' });
    GLYPHS[glyphId]?.(glyph);
    const count = h('span', { class: 'on-count' });
    counts.set(id, count);
    const head = h('button', { class: 'lens-head', onclick: () => { open = open === id ? null : id; draw(); } }, glyph, h('span', {}, h('b', {}, title), h('span', { class: 'q' }, question)), count);
    const wrap = h('section', { class: `lens${open === id ? ' open' : ''}` }, head);
    if (open === id) wrap.append(bodyOf(items, { byGroup: true }));
    return wrap;
  }

  function draw() {
    el.replaceChildren();
    rows.clear();
    counts.clear();
    el.append(h('div', { class: 'lens-order', role: 'group', 'aria-label': 'the order of what can be looked at' },
      ...[['questions', 'by question'], ['works', 'by work']].map(([m, name]) => h('button', { class: m === mode ? 'on' : '', onclick: () => { mode = m; remember('urban-lens-order', m); open = null; draw(); } }, name))));
    if (mode === 'questions') {
      for (const lens of lensList) el.append(section(lens.id, lens.title, lens.question, lens.id, data.filter((d) => d.lens === lens.id)));
    } else {
      for (const w of works) el.append(section(w.name, w.title, w.tagline, w.lens, data.filter((d) => d.plugin.name === w.name)));
    }
    if (planned.length) {
      el.append(h('div', { class: 'about' }, h('button', { title: planned.map((p) => p.title).join(', '), onclick: () => plannedPanel(planned[0]) }, `Planned: ${planned.map((p) => p.title).join(' · ')}`)));
    }
    refresh();
  }

  let lastRange = 0;
  function refresh() {
    if (mode === 'questions') {
      for (const lens of lensList) {
        const n = data.filter((d) => d.lens === lens.id && d.kind === 'layer' && d.isOn()).length;
        const c = counts.get(lens.id);
        if (c) c.textContent = n ? String(n) : '';
      }
    } else {
      for (const w of works) {
        const n = data.filter((d) => d.plugin.name === w.name && d.kind === 'layer' && d.isOn()).length;
        const c = counts.get(w.name);
        if (c) c.textContent = n ? String(n) : '';
      }
    }
    for (const [, r] of rows) {
      const on = r.d.isOn();
      r.input.checked = on;
      r.row.classList.toggle('on', on);
      r.row.classList.toggle('loading', !!r.d.isLoading?.());
    }
    update(lastRange, true);
  }
  function update(range, force = false) {
    if (!force && Math.abs(range - lastRange) / (lastRange || 1) < 0.01) return;
    lastRange = range;
    const u = along(range) * 100;
    for (const [, r] of rows) {
      r.marker.style.left = `${u}%`;
      r.row.classList.toggle('away', !(range >= r.d.scale[0] && range <= r.d.scale[1]) || (r.d.isAway?.() ?? false));
    }
  }

  return {
    set(items, plan, workList) { data = items; planned = plan; works = workList; draw(); },
    open(id) { open = id; draw(); },
    // Open a work by name, in the order of works.
    showWork(name) { mode = 'works'; remember('urban-lens-order', 'works'); open = name; draw(); },
    setTool(key) { tool = key; draw(); },
    refresh,
    update,
  };
}
