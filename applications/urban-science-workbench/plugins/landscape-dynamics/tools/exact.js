// The landscape's events, exactly: the discrete part of the district's state — the two stores open or shut, a rule
// in force (the permit for the lane), the rule carried out (the lane built), a meaning adopted (the lane
// recognised) — composed from its five parts at their interfaces and solved exactly: how often each joint state
// holds and whether the whole runs round cycles; which coarser descriptions forget something exactly and which do
// not, with the two states each confuses and how long it stays close; how three interventions interact on an
// outcome chosen; and what a record of the events says of how much the stores answer each other. The rates are
// the ones the analysis "The landscape and its events" runs the day with.

import { state } from '../shared.js';
import { brief, figures, fold, hint, seg, slider, tabs } from '../ui.js';

const PART = { A: 'the first store', B: 'the second store', P: 'the permit', W: 'the lane', R: 'its recognition' };
const MECHANISM = {
  'A opens': 'the first store opens', 'B opens': 'the second store opens', 'A/shuts': 'the first store shuts', 'B/shuts': 'the second store shuts',
  'P/granted': 'the permit is granted', 'P/withdrawn': 'the permit is withdrawn', 'W built': 'the lane is built', 'W/removed': 'the lane is removed',
  'R gained': 'the lane is recognised', 'R/lapses': 'its recognition lapses',
};
const STORES = [[0, 0], [1, 0], [0, 1], [1, 1]];
const INSTITUTIONS = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]];
const ROWS = ['no permit', 'a rule in force', 'the rule carried out', 'a meaning adopted'];
const ROW_NOTE = ['', 'the permit', 'the lane built', 'the lane recognised'];

export async function render(root, ctx) {
  const { h, svg, failure, fmt, css } = ctx.dom;
  const { lines, bars } = ctx.chart;
  const names = state.district ? state.district.stores.map((s) => s.name) : ['the first store', 'the second store'];
  const COLS = [['both', 'shut'], [tiny(names[0]), 'open'], [tiny(names[1]), 'open'], ['both', 'open']];
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The events, exactly'),
    h('span', { class: 'lede' }, 'Two stores, a rule, the lane it brings and the meaning it takes on — solved without simulation.')));

  const R = state.rules;
  const P = { alphaA: 1.6, alphaB: 1.6, beta: 0.7, etaA: R.eta, etaB: R.eta, xi: R.xi, draw: 0.5, grant: R.permit, revoke: 0.05, build: R.build, remove: 0.06, recognise: R.recognise, lapse: 0.02, gain: 0.5, hours: 500, seed: 5 };
  const share = () => Object.assign(state.rules, { eta: (P.etaA + P.etaB) / 2, xi: P.xi, permit: P.grant, build: P.build, recognise: P.recognise });
  const rate = (key, label, min, max, step, after) => slider(h, fmt, { label, min, max, step, value: P[key], unit: '/h', onChange: (v) => { P[key] = v; if (after) after(); run(); } });
  root.append(fold(h, 'the rates',
    h('div', { class: 'ld-label' }, 'The stores'),
    rate('alphaA', `${short(names[0])} opens`, 0.1, 4, 0.1), rate('alphaB', `${short(names[1])} opens`, 0.1, 4, 0.1), rate('beta', 'each shuts', 0.1, 3, 0.1),
    rate('etaA', `${tiny(names[0])} answers the other`, 0, 3, 0.05, share), rate('etaB', `${tiny(names[1])} answers the other`, 0, 3, 0.05, share), rate('xi', 'each answers the lane', 0, 3, 0.05, share),
    h('div', { class: 'ld-label', style: 'margin-top:10px' }, 'The rule, the lane, its meaning'),
    rate('grant', 'a permit granted', 0.01, 2, 0.01, share), rate('revoke', 'withdrawn', 0.01, 1, 0.01), rate('build', 'the lane built', 0.05, 2, 0.05, share), rate('remove', 'removed', 0.01, 1, 0.01),
    rate('recognise', 'recognised', 0.05, 2, 0.05, share), rate('draw', 'use per open store', 0, 3, 0.05), rate('lapse', 'meaning lapses', 0.005, 0.5, 0.005),
    hint(h, 'A store opens faster when the other is open and when the lane stands. The lane stands only under a permit, and once recognised it cannot be removed.')));
  const out = h('div', {});
  root.append(out);
  let r = null, outcome = 'both';

  const partsPanel = () => h('div', {},
    compositionDiagram(h, svg, css, names),
    brief(ctx, { title: r.coherent ? 'Composed either way, one process' : 'The two orders of composition disagree', kind: 'proved',
      statement: `Five parts, joined where their states must agree: no lane without a permit, no meaning without the lane. ${r.states.length} joint states.`,
      method: 'Open systems with ports, joined by fibre products of their states at the relations of their ports; each change that reads several parts added once, at the first composition holding them all; the two orders compared rate by rate (physis composition).' }));

  const lawPanel = () => {
    const cellOf = (s) => [STORES.findIndex(([a, b]) => a === s.A && b === s.B), INSTITUTIONS.findIndex(([p, w, q]) => p === s.P && w === s.W && q === s.R)];
    const top = [...r.cycles].sort((a, b) => Math.abs(b.current * b.affinity) - Math.abs(a.current * a.affinity))[0];
    return h('div', {},
      lawGrid(h, svg, css, fmt, r, cellOf, COLS, null),
      brief(ctx, { title: r.detailedBalance ? 'In balance' : `Driven: entropy produced at ${fmt(r.entropyProduction, 3)} an hour`, kind: 'proved',
        statement: r.detailedBalance ? 'Every cycle of states is run as often one way as the other.'
          : `It runs round its cycles${top ? `, the busiest ${fmt(Math.abs(top.current), 3)} times an hour net` : ''}: the stores answer the lane, and the lane's meaning answers the stores.`,
        method: 'The stationary law by the GTH algorithm; the currents and affinities of a spanning tree\'s chords; Schnakenberg\'s entropy production (physis).',
        limits: 'It would be in balance only if neither store answered the lane, the lane\'s use did not depend on them, and each answered the other equally.' }));
  };

  const descriptionsPanel = () => {
    const cellOf = (s) => [STORES.findIndex(([a, b]) => a === s.A && b === s.B), INSTITUTIONS.findIndex(([p, w, q]) => p === s.P && w === s.W && q === s.R)];
    const lawBox = h('div', {}, lawGrid(h, svg, css, fmt, r, cellOf, COLS, null));
    const rows = r.views.map((v) => {
      const verdict = v.eventwise ? 'exact' : v.whole ? 'exact as a whole' : 'not exact';
      const kind = v.eventwise ? 'pass' : v.whole ? 'amended' : 'model';
      const a = r.states[v.state], b = r.states[v.other];
      const differ = a && b ? Object.keys(PART).filter((k) => a[k] !== b[k]) : [];
      const confused = differ.length === 2 && differ.includes('A') && differ.includes('B') ? 'which store is open' : differ.map((k) => PART[k]).join(' and ');
      const tr = h('tr', { title: v.eventwise ? '' : `"${MECHANISM[v.mechanism] || v.mechanism}" tells apart two states it cannot: they differ in ${confused}` },
        h('td', { class: 'left' }, v.name), h('td', { class: 'left' }, h('span', { class: `chip ${kind}` }, verdict)), h('td', {}, v.horizon == null ? '—' : hours(v.horizon, fmt)));
      if (!v.eventwise && a && b) {
        tr.addEventListener('mouseenter', () => lawBox.replaceChildren(lawGrid(h, svg, css, fmt, r, cellOf, COLS, [cellOf(a), cellOf(b)])));
        tr.addEventListener('mouseleave', () => lawBox.replaceChildren(lawGrid(h, svg, css, fmt, r, cellOf, COLS, null)));
      }
      return tr;
    });
    const alone = r.views[0];
    return h('div', {},
      h('table', { class: 'data', style: 'width:100%' }, h('thead', {}, h('tr', {}, h('th', { class: 'left' }, 'a description that keeps'), h('th', { class: 'left' }, ''), h('th', {}, 'within 10% for'))), h('tbody', {}, rows)),
      hint(h, 'Hover a row: the two states it confuses light up below.'),
      lawBox,
      brief(ctx, { title: alone.eventwise ? 'The stores alone forget the institutions exactly' : 'The stores alone remember the lane', kind: 'proved',
        statement: 'A description is a process of its own exactly when no change it keeps reads what it forgets. The permit and the meaning can never be forgotten: building reads the one, removal waits on the other.',
        method: 'Strong lumpability for each mechanism and for the whole; the aggregation by the stationary law, its outward discrepancy δ, and the coupling bound 1 − exp(−δT) (physis lumping).' }),
      fold(h, 'how long both stores stay shut',
        lines([{ name: 'the true law', x: r.dwell.tau, y: r.dwell.density, color: css('--accent-deep'), width: 2 }, { name: 'an exponential of the same mean', x: r.dwell.tau, y: r.dwell.exponential, color: css('--minor'), dash: '5 4' }],
          { xlabel: 'hours', ylabel: 'density (log)', height: 150, logy: true }),
        h('p', {}, `Seen by the stores alone, ${fmt(60 * r.dwell.mean, 3)} minutes on average; its spread (squared coefficient of variation ${fmt(r.dwell.cv2, 4)}, an exponential's 1) is what the description remembers.`)));
  };

  const interventionsPanel = () => {
    const g = r.interventions[outcome];
    const main = [[1, 'the first'], [2, 'the second'], [4, 'permits']];
    const joint = [[3, 'both stores'], [5, 'first, permits'], [6, 'second, permits'], [7, 'all three']];
    const strongest = joint.reduce((m, j) => (Math.abs(g[j[0]]) > Math.abs(g[m[0]]) ? j : m), joint[0]);
    return h('div', {},
      h('div', { class: 'ld-row' }, seg(h, [['both', 'Both open'], ['open', 'Stores open'], ['recognised', 'Lane recognised']], outcome, (v) => { outcome = v; T.show('Interventions'); })),
      slider(h, fmt, { label: 'each one adds', min: 0.05, max: 2, step: 0.05, value: P.gain, onChange: (v) => { P.gain = v; run(); } }),
      h('div', { class: 'grid cols-2' },
        h('div', {}, h('div', { class: 'ld-label' }, 'each alone'), bars(main.map(([s], i) => ({ x: i, y: g[s], color: css('--series-1'), label: `${main[i][1]}: ${fmt(g[s], 3)}` })), { labels: main.map((m) => m[1]), height: 170, width: 300 })),
        h('div', {}, h('div', { class: 'ld-label' }, 'beyond each alone'), bars(joint.map(([s], i) => ({ x: i, y: g[s], color: g[s] >= 0 ? css('--series-4') : css('--series-5'), label: `${joint[i][1]}: ${fmt(g[s], 3)}` })), { labels: joint.map((j) => j[1]), height: 170, width: 300 }))),
      brief(ctx, { title: `Strongest together: ${strongest[1]} (${fmt(g[strongest[0]], 3)})`, kind: 'proved',
        statement: 'What doing several together adds beyond each alone. It belongs to the outcome: stores that ignore each other show none on how many are open, yet one on whether both are.',
        method: 'The anchored Möbius inversion over the subsets of the three interventions (promote the first store, the second, grant permits faster), each outcome at stationarity (morphe).' }));
  };

  const evidencePanel = () => {
    const ev = r.evidence;
    return h('div', {},
      bars(ev.worlds.map((w, k) => ({ x: k, y: ev.posterior[k], color: w === 0 ? css('--series-3') : css('--accent'), label: `${fmt(w)}: ${fmt(100 * ev.posterior[k], 3)}%` })),
        { labels: ev.worlds.map((w) => fmt(w)), ylabel: 'how much they answer each other, after the record', height: 180 }),
      figures(h, [{ value: fmt(ev.coupled, 3), label: 'P(they answer each other)' }, { value: fmt(ev.strong, 3), label: 'P(by 1 or more)' }, { value: ev.events, label: `changes in ${fmt(ev.hours)} h` }]),
      fold(h, 'the record',
        slider(h, fmt, { label: 'hours of record', min: 50, max: 3000, step: 50, value: P.hours, onChange: (v) => { P.hours = v; run(); } }),
        slider(h, fmt, { label: 'which record', min: 1, max: 40, step: 1, value: P.seed, onChange: (v) => { P.seed = v; run(); } }),
        h('p', {}, 'A record made by the model itself, every change seen; each answer weighed by the record\'s exact path density under it, from a uniform prior.')));
  };

  const T = tabs(h, [{ name: 'Parts', panel: () => partsPanel() }, { name: 'Law', panel: () => lawPanel() }, { name: 'Descriptions', panel: () => descriptionsPanel() },
    { name: 'Interventions', panel: () => interventionsPanel() }, { name: 'Evidence', panel: () => evidencePanel() }], { remember: 'landscape-exact', lazy: true });
  let current = null;
  async function run() {
    const shownTab = out.querySelector('.ld-tab.on');
    current = shownTab ? shownTab.dataset.name : current;
    if (!r) out.replaceChildren(hint(h, 'solving the sixteen joint states…'));
    try {
      r = await ctx.ask('skeleton', { alphaA: P.alphaA, alphaB: P.alphaB, betaA: P.beta, betaB: P.beta, etaA: P.etaA, etaB: P.etaB, xiA: P.xi, xiB: P.xi, draw: P.draw,
        grant: P.grant, revoke: P.revoke, build: P.build, remove: P.remove, recognise: P.recognise, lapse: P.lapse, gain: P.gain, hours: P.hours, seed: P.seed });
    } catch (e) { return failure(out, e); }
    out.replaceChildren(T.el);
    T.show(current);
  }
  run();
}

const short = (name) => (name.length > 22 ? `${name.slice(0, 20)}…` : name);
const tiny = (name) => (name === 'the first store' ? 'the first' : name === 'the second store' ? 'the second' : name.length > 15 ? `${name.slice(0, 14)}…` : name);
const hours = (x, fmt) => (x == null ? 'always' : x < 1 ? `${fmt(60 * x, 2)} minutes` : `${fmt(x, 3)} hours`);

// The sixteen joint states as a table: the stores across, the institutions down, shaded by the share of time.
function lawGrid(h, svg, css, fmt, r, cellOf, cols, mark) {
  const W = 600, H = 250, L = 150, T = 34, cw = (W - L - 10) / 4, ch = (H - T - 10) / 4;
  const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  cols.forEach(([a, b], i) => { s.append(svg('text', { x: L + cw * (i + 0.5), y: 13, 'text-anchor': 'middle', 'font-size': 11.5, fill: css('--ink-2') }, a)); s.append(svg('text', { x: L + cw * (i + 0.5), y: 27, 'text-anchor': 'middle', 'font-size': 11.5, fill: css('--ink-faint') }, b)); });
  ROWS.forEach((name, j) => {
    s.append(svg('text', { x: L - 10, y: T + ch * (j + 0.5) - 2, 'text-anchor': 'end', 'font-size': 12.5, fill: css('--ink') }, name));
    if (ROW_NOTE[j]) s.append(svg('text', { x: L - 10, y: T + ch * (j + 0.5) + 13, 'text-anchor': 'end', 'font-size': 10.5, fill: css('--ink-faint') }, ROW_NOTE[j]));
  });
  const top = Math.max(...r.states.map((x) => x.p));
  const accent = css('--accent-deep');
  for (const st of r.states) {
    const [i, j] = cellOf(st);
    if (i < 0 || j < 0) continue;
    const t = Math.sqrt(st.p / top);
    const x = L + cw * i, y = T + ch * j;
    s.append(svg('rect', { x: x + 2, y: y + 2, width: cw - 4, height: ch - 4, rx: 3, fill: accent, 'fill-opacity': (0.06 + 0.8 * t).toFixed(3), stroke: css('--hair-2'), 'stroke-width': 0.6 }));
    s.append(svg('text', { x: x + cw / 2, y: y + ch / 2 + 4, 'text-anchor': 'middle', 'font-size': 12, fill: t > 0.55 ? css('--surface') : css('--ink') }, `${fmt(100 * st.p, 2)}%`));
  }
  if (mark) for (const [i, j] of mark) s.append(svg('rect', { x: L + cw * i + 1, y: T + ch * j + 1, width: cw - 2, height: ch - 2, rx: 4, fill: 'none', stroke: css('--minor'), 'stroke-width': 2.4 }));
  return h('div', { class: 'ld-card', style: 'margin:6px 0' }, s);
}

// The five parts, where they meet, and the joint changes on the parts they read.
function compositionDiagram(h, svg, css, names) {
  const W = 600, H = 236;
  const at = { A: [80, 55], B: [80, 165], P: [280, 110], W: [420, 110], R: [560, 110] };
  const label = { A: names[0], B: names[1], P: 'the permit', W: 'the lane', R: 'its meaning' };
  const s = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  // Joint changes: a small diamond joined to the parts each reads.
  const joint = [
    { name: 'opening', at: [190, 40], parts: ['A', 'B', 'W'] },
    { name: 'opening', at: [190, 180], parts: ['A', 'B', 'W'] },
    { name: 'building', at: [350, 36], parts: ['P', 'W'] },
    { name: 'recognition', at: [300, 205], parts: ['A', 'B', 'W', 'R'] },
  ];
  for (const j of joint) {
    for (const p of j.parts) s.append(svg('line', { x1: j.at[0], y1: j.at[1], x2: at[p][0], y2: at[p][1], stroke: css('--minor'), 'stroke-opacity': 0.55, 'stroke-width': 1.1, 'stroke-dasharray': '3 3' }));
  }
  // Interfaces: the lane under the permit, the meaning under the lane.
  for (const [a, b, rel] of [['P', 'W', 'lane ≤ permit'], ['W', 'R', 'meaning ≤ lane']]) {
    s.append(svg('line', { x1: at[a][0] + 26, y1: at[a][1], x2: at[b][0] - 26, y2: at[b][1], stroke: css('--accent-deep'), 'stroke-width': 2 }));
    s.append(svg('text', { x: (at[a][0] + at[b][0]) / 2, y: at[a][1] - 8, 'text-anchor': 'middle', 'font-size': 10, fill: css('--accent-deep') }, rel));
  }
  for (const j of joint) {
    const [x, y] = j.at;
    s.append(svg('path', { d: `M${x},${y - 7} L${x + 7},${y} L${x},${y + 7} L${x - 7},${y} Z`, fill: css('--minor'), stroke: css('--surface'), 'stroke-width': 1 }));
    s.append(svg('text', { x, y: y + 20, 'text-anchor': 'middle', 'font-size': 10.5, fill: css('--ink-2') }, j.name));
  }
  for (const [k, [x, y]] of Object.entries(at)) {
    s.append(svg('circle', { cx: x, cy: y, r: 24, fill: css('--surface'), stroke: css('--accent'), 'stroke-width': 1.6 }));
    s.append(svg('text', { x, y: y + 4.5, 'text-anchor': 'middle', 'font-size': 13, 'font-weight': 600, fill: css('--ink') }, k));
    const nm = label[k].length > 18 ? `${label[k].slice(0, 17)}…` : label[k];
    s.append(svg('text', { x, y: k === 'A' ? y - 32 : y + 40, 'text-anchor': 'middle', 'font-size': 11, fill: css('--ink-2') }, nm));
  }
  return h('div', { class: 'ld-card' }, s,
    h('p', { class: 'ld-hint' }, 'circles: the parts · solid: where they meet · diamonds: changes that read several parts'));
}
