// What can be looked at, ordered by the question it answers rather than by the plug-in that holds it. Each lens opens in
// place and shows its layers (each with a small copy of the altitude tape, marking the scales at which it means
// something, and dimmed when the camera is elsewhere) and its analyses. One lens is open at a time.

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
};

export function createLenses(el, { lensList, onToggle, onTool, plannedPanel, inArea }) {
  let open = null, data = [], planned = [], tool = null;
  const rows = new Map();   // key → { row, input, marker }
  const counts = new Map(); // lens → the count element

  function miniScale(near, far) {
    const a = along(far) * 100, b = along(near) * 100;
    const marker = h('b', {});
    return { el: h('span', { class: 'scale-mini', title: 'the scales at which this means something' }, h('i', { style: `left:${a}%;width:${Math.max(4, b - a)}%` }), marker), marker };
  }

  function draw() {
    el.replaceChildren();
    rows.clear();
    counts.clear();
    for (const lens of lensList) {
      const items = data.filter((d) => d.lens === lens.id);
      const layers = items.filter((d) => d.kind === 'layer'), tools = items.filter((d) => d.kind === 'tool');
      const glyph = svg('svg', { viewBox: '0 0 24 24' });
      GLYPHS[lens.id]?.(glyph);
      const count = h('span', { class: 'on-count' });
      counts.set(lens.id, count);
      const head = h('button', { class: 'lens-head', onclick: () => { open = open === lens.id ? null : lens.id; draw(); } }, glyph, h('span', {}, h('b', {}, lens.title), h('span', { class: 'q' }, lens.question)), count);
      const wrap = h('section', { class: `lens${open === lens.id ? ' open' : ''}` }, head);
      if (open === lens.id) {
        const body = h('div', { class: 'lens-body' });
        let group = null;
        for (const d of layers) {
          if (d.group && d.group !== group) { group = d.group; body.append(h('div', { class: 'lens-group' }, group)); }
          const input = h('input', { type: 'checkbox', checked: d.isOn() ? true : null, onchange: () => { onToggle(d, input.checked); refresh(); } });
          const mini = miniScale(d.scale[0], d.scale[1]);
          const row = h('label', { class: 'layer-row', title: d.about || '' }, input, h('span', { class: 'box' }), h('span', { class: 'name' }, d.title), mini.el);
          rows.set(d.key, { row, input, marker: mini.marker, d });
          body.append(row);
        }
        if (tools.length) body.append(h('div', { class: 'lens-group' }, 'Analyses'));
        for (const t of tools) {
          body.append(h('button', { class: `tool-row${tool === t.key ? ' active' : ''}`, title: t.about || '', onclick: () => onTool(t) }, t.title));
        }
        if (!layers.length && !tools.length) body.append(h('div', { class: 'lens-foot' }, 'Nothing here yet for this place.'));
        wrap.append(body);
      }
      el.append(wrap);
    }
    if (planned.length) {
      el.append(h('div', { class: 'about' }, h('button', { title: planned.map((p) => p.title).join(', '), onclick: () => plannedPanel(planned[0]) }, `Planned: ${planned.map((p) => p.title).join(' · ')}`)));
    }
    refresh();
  }

  let lastRange = 0;
  function refresh() {
    for (const lens of lensList) {
      const n = data.filter((d) => d.lens === lens.id && d.kind === 'layer' && d.isOn()).length;
      const c = counts.get(lens.id);
      if (c) c.textContent = n ? String(n) : '';
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
    set(items, plan) { data = items; planned = plan; draw(); },
    open(id) { open = id; draw(); },
    setTool(key) { tool = key; draw(); },
    refresh,
    update,
  };
}
