// The mix of uses: along a street (OpenStreetMap's holdings) or about a Melbourne sensor (the census of land
// use): how far the street departs from an evenly mixed one, and at what scale; the distance over which its uses
// change; how the mix has moved over the census years; and, across Melbourne, why no single score ranks streets.

export async function render(root, ctx) {
  const { h, svg, busy, failure, fmt, css, hover } = ctx.dom;
  const { lines, scatter, legend, category } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Mix of uses'),
    h('span', { class: 'lede' }, 'Pick a sensor or a shop on the map, or choose a street here.')));
  let streets = [], sensors = [];
  [streets, sensors] = await Promise.all([ctx.ask('streets').catch(() => []), ctx.ask('sensors').catch(() => [])]);
  streets = streets.filter((s) => s.area === ctx.area.id);
  const sel = h('select', {});
  for (const s of streets) sel.append(h('option', { value: `street:${s.id}` }, `${s.name} (${s.holdings} holdings)`));
  if (ctx.area.id === 'melbourne') for (const s of sensors.filter((x) => x.terrain)) sel.append(h('option', { value: `sensor:${s.id}` }, `about ${s.name}`));
  const walk = h('input', { type: 'range', min: 100, max: 1200, step: 50, value: 400 });
  const walkOut = h('output', {}, '400 m');
  root.append(h('div', { class: 'keys' }, h('label', { class: 'pick' }, sel), h('label', { class: 'slider', style: 'min-width: 24em' }, 'a visit walks', walk, walkOut)));
  const out = h('div', {});
  const dom = h('div', { class: 'section' });
  root.append(out, dom);
  if (!sel.options.length) return failure(out, new Error('no streets or sensors here yet: the records are not built'));
  const choose = (s) => {
    if (!s) return;
    const v = s.kind === 'sensor' ? `sensor:${s.id}` : s.kind === 'holding' ? `street:${s.street}` : null;
    if (v && [...sel.options].some((o) => o.value === v)) { sel.value = v; return true; }
    return false;
  };
  choose(ctx.selection);
  ctx.onSelect((s) => { if (choose(s)) run(); });
  const useLegend = (uses, present) => legend(uses.map((u, i) => ({ name: u, color: category(i), dot: true })).filter((_, i) => present[i]));

  const run = async () => {
    busy(out);
    const [kind, id] = sel.value.split(':');
    const w = +walk.value;
    walkOut.textContent = `${w} m`;
    let d;
    try { d = kind === 'street' ? await ctx.ask('street', { street: id, walk: w }) : await ctx.ask('terrain', { sensor: +id, walk: w }); }
    catch (e) { return failure(out, e); }
    const blocks = [];
    if (kind === 'street') {
      const W = 600, H = 130, L = 14, R = 14, mid = H / 2;
      const X = (s) => L + (s / d.length) * (W - L - R);
      const maxOff = Math.max(20, ...d.holdings.map((x) => Math.abs(x.offset)));
      const strip = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%' });
      strip.append(svg('line', { x1: L, x2: W - R, y1: mid, y2: mid, stroke: css('--ink-faint'), 'stroke-width': 1.5 }));
      for (let s = 0; s <= d.length; s += 200) strip.append(svg('text', { x: X(s), y: H - 3, 'text-anchor': 'middle', class: 'mark-label' }, `${s} m`));
      for (const x of d.holdings) {
        const c = svg('circle', { cx: X(x.s), cy: mid - (x.offset / maxOff) * (mid - 14), r: 3.6, fill: category(x.use), stroke: css('--surface'), 'stroke-width': 0.8 });
        hover(c, `${x.name || '(unnamed)'}\n${d.uses[x.use]}\n${Math.round(x.s)} m along`);
        strip.append(c);
      }
      blocks.push(h('h3', {}, `${d.name}: ${d.holdings.length} holdings along ${fmt(d.length, 4)} m`), strip, useLegend(d.uses, d.uses.map((_, i) => d.holdings.some((x) => x.use === i))));
    } else {
      const pts = d.points.map(([x, y, u]) => ({ x, y, color: category(u), r: 3, label: d.uses[u] }));
      blocks.push(h('h3', {}, `About ${d.name}: ${d.points.length} businesses within 200 m (${d.year})`),
        scatter(pts, { xlabel: 'east (m)', ylabel: 'north (m)', xlim: [-210, 210], ylim: [-210, 210], height: 380, width: 420 }), useLegend(d.uses, d.uses.map((_, i) => d.points.some((p) => p[2] === i))));
    }
    if (d.mixing) {
      const m = d.mixing, c = d.correlation;
      const pts = c.r.map((r, i) => ({ r, rho: c.rho[i] })).filter((p) => Number.isFinite(p.rho));
      blocks.push(h('div', { class: 'section' },
        ctx.finding({ title: `${fmt(100 * m.normalised, 3)}% away from an even mix`, kind: 'measured',
          statement: `How much knowing where you are tells you about which use you meet (0: every place holds the uses in the street's proportions; 100%: every place holds one use). Of it, ${fmt(100 * m.between / Math.max(1e-12, m.departure), 3)}% lies between ${m.stretch}-metre ${kind === 'street' ? 'stretches' : 'blocks'} and the rest within them, among places of ${m.place} m. Exchange that ignores what things are can only lower it, so whatever keeps a street unmixed depends on the uses themselves.`,
          method: 'The mutual information of place and use, split exactly between the two scales.', limits: 'Uses only; price, language and standing are not registered.' }),
        lines([{ name: 'correlation of uses', x: pts.map((p) => p.r), y: pts.map((p) => p.rho), color: css('--series-1') }, { x: [0, pts[pts.length - 1].r], y: [0, 0], color: css('--hair-2'), dash: '3 3' }],
          { xlabel: 'distance between two holdings (m)', ylabel: 'how much more often they share a use', height: 190, marks: [{ x: c.firstZero, label: `${fmt(c.firstZero, 3)} m` }] }),
        ctx.finding({ title: c.singleRegime ? 'One use only' : `Uses change every ${fmt(c.length, 3)} m; a visit of ${w} m crosses ${fmt(c.compression, 3)} of them`, kind: 'measured',
          statement: c.singleRegime ? 'No change to measure.' : 'The distance over which two holdings stop sharing their use more than chance would, and how many such changes a visit of the chosen length passes.',
          method: 'The correlation of uses at each separation, integrated to its first zero.', limits: 'The length of a visit is chosen, not measured.' })));
    }
    if (d.record && d.record.length > 1) {
      blocks.push(h('div', { class: 'section' }, h('h3', {}, 'Over the census years'),
        lines([{ name: 'away from an even mix', x: d.record.map((x) => x.year), y: d.record.map((x) => 100 * x.normalised), color: css('--series-2') }],
          { xlabel: 'census year', ylabel: '% away from an even mix', height: 180, points: d.record.map((x) => ({ x: x.year, y: 100 * x.normalised, label: `${x.year}: ${x.establishments} businesses`, color: css('--series-2') })) })));
    }
    out.replaceChildren(...blocks);
  };

  const dominance = async () => {
    if (ctx.area.id !== 'melbourne') return;
    busy(dom);
    let d;
    try { d = await ctx.ask('dominance', {}); } catch { return dom.replaceChildren(); }
    const n = d.items.length, cols = d.rankings.length;
    const W = 600, H = Math.max(240, 11 * n + 40), L = 150, R = 150, T = 22;
    const X = (c) => L + (c / Math.max(1, cols - 1)) * (W - L - R), Y = (k) => T + (k / Math.max(1, n - 1)) * (H - T - 16);
    const bump = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%' });
    const rankOf = d.rankings.map((rk) => { const p = new Array(n); rk.rank.forEach((item, k) => { p[item] = k; }); return p; });
    const maximal = new Set(d.maximal);
    for (let i = 0; i < n; i++) {
      let path = '';
      for (let c = 0; c < cols; c++) path += `${c ? 'L' : 'M'}${X(c)},${Y(rankOf[c][i])}`;
      const p = svg('path', { d: path, fill: 'none', stroke: maximal.has(i) ? css('--accent') : css('--ink-faint'), 'stroke-width': maximal.has(i) ? 1.8 : 0.8, opacity: maximal.has(i) ? 0.95 : 0.5 });
      hover(p, d.items[i].name);
      bump.append(p);
    }
    d.rankings.forEach((rk, c) => bump.append(svg('text', { x: X(c), y: 12, 'text-anchor': 'middle', class: 'mark-label' }, `reader ${c + 1}`)));
    rankOf[0].forEach((k, i) => bump.append(svg('text', { x: L - 5, y: Y(k) + 3, 'text-anchor': 'end', class: 'mark-label' }, d.items[i].name.slice(0, 24))));
    dom.replaceChildren(h('h2', {}, 'No single score ranks streets'), bump,
      h('ol', { class: 'caption' }, d.rankings.map((rk) => h('li', {}, rk.name))),
      ctx.finding({ title: 'Readers disagree except where one street has more of everything', kind: 'proved',
        statement: `Each street's surroundings counted in unlike units (${d.components.join('; ')}). Four readers weigh them differently, and their rankings cross; the last two are the same reader measuring length in metres and in kilometres. Only where one street has at least as much of every component as another do all readers agree (rose lines: the streets nothing outdoes in every respect).`,
        method: 'Rankings by weighted sums; the order every non-negative weighting agrees on.' }));
  };
  sel.addEventListener('change', run);
  walk.addEventListener('change', run);
  walk.addEventListener('input', () => { walkOut.textContent = `${walk.value} m`; });
  run();
  dominance();
}
