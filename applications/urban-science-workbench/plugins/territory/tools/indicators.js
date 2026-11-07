// Any indicator on the globe: the World Bank's Data360 searched by what an indicator means, one chosen and raised over the
// world's economies, in the mosaic of the cells their people live in. It is the same search that data360-mcp offers to language
// models, asked here directly.

import { compact, nearest, tone, world } from '../shared.js';

const SETS = [{ res: 2, fade: [5.5e6, 9e7], top: 1.5e6 }, { res: 3, fade: [9e5, 5.5e6], top: 7.5e5 }];

export async function render(root, ctx) {
  const { earth, dom: { h } } = ctx;
  const input = h('input', { type: 'search', class: 'wide', placeholder: 'what it means: access to electricity, maternal mortality, internet users…', 'aria-label': 'what the indicator means' });
  const go = h('button', { class: 'primary' }, 'Search');
  const list = h('div', { class: 'choices' });
  const shown = h('div', {});
  root.append(h('h1', {}, 'Any indicator on the globe'),
    h('p', { class: 'muted' }, 'Say what you are after in your own words: Data360 searches its indicators by meaning. Choose one, and the economies rise over the globe in its values, the latest year each one reports.'),
    h('form', { class: 'search-row', onsubmit: (e) => { e.preventDefault(); search(); } }, input, go), list, shown);
  const keys = [];
  const clear = () => { for (const k of keys.splice(0)) earth.removeColumns(k); };
  ctx.onDispose(clear);

  async function search() {
    const q = input.value.trim();
    if (!q) return;
    list.replaceChildren(h('p', { class: 'busy' }, 'searching…'));
    const r = await fetch(`/world/data360?q=${encodeURIComponent(q)}`);
    const a = await r.json();
    if (!r.ok) { list.replaceChildren(h('p', { class: 'error-box' }, a.message)); return; }
    if (!a.found.length) { list.replaceChildren(h('p', { class: 'muted' }, 'Nothing by that meaning.')); return; }
    list.replaceChildren(...a.found.map((d) => h('button', { class: 'choice', onclick: () => raise(d) }, h('b', {}, d.name), h('span', {}, `${d.database}${d.unit ? ` · ${d.unit}` : ''}`), d.definition ? h('em', {}, d.definition) : null)));
  }

  async function raise(d) {
    shown.replaceChildren(h('p', { class: 'busy' }, `reading every economy's ${d.name}…`));
    const r = await fetch(`/world/data360?db=${encodeURIComponent(d.database)}&id=${encodeURIComponent(d.id)}`);
    const a = await r.json();
    if (!r.ok) { shown.replaceChildren(h('p', { class: 'error-box' }, a.message)); return; }
    const latest = {};
    for (const [code, series] of Object.entries(a.values)) {
      const ys = Object.keys(series).sort();
      if (ys.length) latest[code] = { year: ys[ys.length - 1], value: series[ys[ys.length - 1]] };
    }
    const values = Object.values(latest).map((x) => x.value).filter(Number.isFinite).sort((p, q) => p - q);
    if (!values.length) { shown.replaceChildren(h('p', { class: 'muted' }, 'This indicator gives no value by economy.')); return; }
    const low = values[Math.floor(values.length * 0.02)], high = values[Math.floor(values.length * 0.98)];
    const u = (v) => (high > low ? Math.min(1, Math.max(0, (v - low) / (high - low))) : 0.5);
    const colour = tone('plum', ctx.theme, { from: 0.14 });
    clear();
    let named = 0;
    const names = {};
    for (const S of SETS) {
      const w = await world(ctx, S.res);
      w.code.forEach((c, k) => { if (c && w.name[k]) names[c] = w.name[k]; });
      const key = `indicator:${S.res}`;
      keys.push(key);
      const set = earth.addColumns(key, w.rings.map((ring) => ({ ring })), { fade: S.fade, shrink: 0.9 });
      const got = w.code.map((c) => latest[c]?.value);
      named = new Set(w.code.filter((c) => latest[c])).size;
      set.prisms.set(got.map((v) => (Number.isFinite(v) ? S.top * Math.max(0.02, u(v)) : S.top * 0.012)), got.map((v) => (Number.isFinite(v) ? colour(u(v)) : (ctx.theme === 'night' ? [0.13, 0.09, 0.16] : [0.86, 0.82, 0.86]))));
    }
    if (earth.view.range < 9e5) earth.flyTo({ range: 1.6e7, tilt: 0 }, 2.4);
    // Economies only, by the names of the regions product: a database's aggregates are left out of what is ranked.
    const sorted = Object.entries(latest).filter(([c, x]) => names[c] && Number.isFinite(x.value)).sort((a, b) => b[1].value - a[1].value);
    const name = (code) => names[code];
    shown.replaceChildren(ctx.finding({
      title: d.name, kind: 'measured',
      statement: `${named} economies on the globe, each at the latest year it reports. Highest: ${sorted.slice(0, 3).map(([c, x]) => `${name(c)} ${compact(x.value)}`).join(', ')}. Lowest: ${sorted.slice(-3).map(([c, x]) => `${name(c)} ${compact(x.value)}`).join(', ')}.`,
      method: `${d.database}, indicator ${d.id}; column height and colour from the 2nd to the 98th percentile${d.unit ? `, in ${d.unit}` : ''}.`,
      limits: 'Years differ from economy to economy, and an economy\'s value says nothing of the cells within it: the hexagons only show where its people live. A cell holding any part of an area in dispute takes no economy\'s value.',
      against: 'World Bank Data360, retrieved through its public API when asked.',
    }));
  }
  input.focus();
}
