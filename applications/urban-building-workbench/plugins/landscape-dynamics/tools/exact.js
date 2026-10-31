// The landscape's events, exactly: the discrete part of the district's state — the two stores open or shut, a rule
// in force (the permit for the lane), the rule carried out (the lane built), a meaning adopted (the lane
// recognised) — composed from its five parts at their interfaces and solved exactly: how often each joint state
// holds and whether the whole runs round cycles; which coarser descriptions forget something exactly and which do
// not, with the two states each confuses and how long it stays close; how three interventions interact on an
// outcome chosen; and what a record of the events says of how much the stores answer each other. The rates are
// the ones the analysis "The landscape and its events" runs the day with.

import { state } from '../shared.js';

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
  const { h, svg, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, legend } = ctx.chart;
  const names = state.district ? state.district.stores.map((s) => s.name) : ['the first store', 'the second store'];
  const COLS = [['both', 'shut'], [tiny(names[0]), 'open'], [tiny(names[1]), 'open'], ['both', 'open']];
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The landscape\'s events, exactly'),
    h('span', { class: 'lede' }, `${names[0]} and ${names[1]} open and shut; a rule may come into force, be carried out as a lane between them, and the lane take on a meaning as it is used. Five parts, composed where they meet, and solved without simulation.`)));

  const R = state.rules;
  const P = { alphaA: 1.6, alphaB: 1.6, beta: 0.7, etaA: R.eta, etaB: R.eta, xi: R.xi, draw: 0.5, grant: R.permit, revoke: 0.05, build: R.build, remove: 0.06, recognise: R.recognise, lapse: 0.02, gain: 0.5, hours: 500, seed: 5 };
  const slider = (key, name, min, max, step, after = () => {}) => {
    const input = h('input', { type: 'range', min, max, step, value: P[key] });
    const o = h('output', {}, fmt(P[key]));
    input.addEventListener('input', () => { o.textContent = fmt(+input.value); });
    input.addEventListener('change', () => { P[key] = +input.value; after(); run(); });
    return h('label', { class: 'slider' }, name, input, o);
  };
  const share = () => Object.assign(state.rules, { eta: (P.etaA + P.etaB) / 2, xi: P.xi, permit: P.grant, build: P.build, recognise: P.recognise });
  root.append(h('div', { class: 'grid', style: 'margin-top:10px;gap:10px' },
    h('div', { class: 'panel' }, h('h3', {}, 'The stores'),
      slider('alphaA', `${short(names[0])} opens`, 0.1, 4, 0.1), slider('alphaB', `${short(names[1])} opens`, 0.1, 4, 0.1), slider('beta', 'each shuts', 0.1, 3, 0.1),
      slider('etaA', `${short(names[0])} answers the other`, 0, 3, 0.05, share), slider('etaB', `${short(names[1])} answers the other`, 0, 3, 0.05, share),
      slider('xi', 'each answers the lane', 0, 3, 0.05, share)),
    h('div', { class: 'panel' }, h('h3', {}, 'The rule, the lane, its meaning'),
      slider('grant', 'a permit granted', 0.01, 2, 0.01, share), slider('revoke', 'withdrawn', 0.01, 1, 0.01), slider('build', 'the lane built', 0.05, 2, 0.05, share),
      slider('remove', 'removed', 0.01, 1, 0.01), slider('recognise', 'recognised', 0.05, 2, 0.05, share), slider('draw', 'each open store adds to its use', 0, 3, 0.05),
      slider('lapse', 'its meaning lapses', 0.005, 0.5, 0.005))),
    h('p', { class: 'caption' }, 'Rates per hour. A store opens at its rate times 1 + how much it answers the other (when that one is open) + how much it answers the lane (when the lane stands); the lane is recognised at its rate times its use, which grows with each store open. The lane can stand only under a permit and be recognised only while it stands, so a permit cannot be withdrawn under a standing lane, nor a recognised lane removed.'));
  const out = h('div', {});
  root.append(out);

  async function run() {
    busy(out, 'composing the five parts and solving the sixteen joint states…');
    let r;
    try {
      r = await ctx.ask('skeleton', { alphaA: P.alphaA, alphaB: P.alphaB, betaA: P.beta, betaB: P.beta, etaA: P.etaA, etaB: P.etaB, xiA: P.xi, xiB: P.xi, draw: P.draw,
        grant: P.grant, revoke: P.revoke, build: P.build, remove: P.remove, recognise: P.recognise, lapse: P.lapse, gain: P.gain, hours: P.hours, seed: P.seed });
    } catch (e) { return failure(out, e); }
    const cellOf = (s) => [STORES.findIndex(([a, b]) => a === s.A && b === s.B), INSTITUTIONS.findIndex(([p, w, q]) => p === s.P && w === s.W && q === s.R)];
    const lawBox = h('div', {});
    const drawLaw = (mark = null) => lawBox.replaceChildren(lawGrid(h, svg, css, fmt, r, cellOf, COLS, mark));
    drawLaw();

    // ---- the parts ----
    const parts = h('div', { class: 'section' }, h('h2', {}, 'Five parts, composed where they meet'),
      compositionDiagram(h, svg, css, names),
      ctx.finding({ title: r.coherent ? 'Composed in either order, the same process' : 'The two orders of composition disagree', kind: 'proved',
        statement: `Each part keeps its own changes; each joint change is declared once, on the parts it reads — a store's opening reads both stores and the lane, the lane's building reads the permit, its recognition reads the lane and both stores. Composing ((the stores, the permit), the lane), its meaning and the stores with (the permit with (the lane with its meaning)) give the same rate for every change between every two of the ${r.states.length} joint states. Where two parts meet, only states that agree are joined: no lane without a permit, no meaning without the lane.`,
        method: 'Open systems with ports, joined by fibre products of their states at the relations of their ports; each joint mechanism added at the first composition that holds all its parts (physis composition).' }));

    // ---- the law ----
    const top = [...r.cycles].sort((a, b) => Math.abs(b.current * b.affinity) - Math.abs(a.current * a.affinity))[0];
    const law = h('div', { class: 'section' }, h('h2', {}, 'How often each joint state holds, and whether it runs round'),
      lawBox,
      h('p', { class: 'caption' }, 'Each square a joint state, shaded by the share of time it holds in the long run. Hover a row of the table below to see the two states a description confuses.'),
      ctx.finding({ title: r.detailedBalance ? 'In balance: every cycle of states is run as often one way as the other' : `Driven: it runs round its cycles, producing entropy at ${fmt(r.entropyProduction, 3)} an hour`, kind: 'proved',
        statement: r.detailedBalance
          ? 'The shares of time are those of an energy, the stores\' and the institutions\' separately: the stores open more readily together, and nothing runs one way round.'
          : `The stores answer the lane, and the lane's meaning answers the stores (it is recognised faster when they are open): influences that do not return in kind. Round its busiest cycle (${top ? top.states.length : 0} states) the process runs ${fmt(Math.abs(top ? top.current : 0), 3)} times an hour net, against an affinity of ${fmt(Math.abs(top ? top.affinity : 0), 3)}. It would be in balance only if neither store answered the lane, the lane's use did not depend on them, and each answered the other equally.`,
        method: 'The stationary law by the GTH algorithm; the cycles of a spanning tree\'s chords with their currents and affinities; Schnakenberg\'s entropy production (physis).' }));

    // ---- coarser descriptions ----
    const rows = r.views.map((v) => {
      const verdict = v.eventwise ? 'exact for every change' : v.whole ? 'exact as a whole' : 'not exact';
      const kind = v.eventwise ? 'pass' : v.whole ? 'amended' : 'model';
      const a = r.states[v.state], b = r.states[v.other];
      const differ = a && b ? Object.keys(PART).filter((k) => a[k] !== b[k]) : [];
      const confused = differ.length === 2 && differ.includes('A') && differ.includes('B') ? 'which store is open' : differ.map((k) => PART[k]).join(' and ');
      const mech = MECHANISM[v.mechanism] || v.mechanism;
      const why = v.eventwise ? 'no change it keeps reads what it forgets'
        : v.whole ? `the whole is a process of its own, but "${mech}" alone is not: it leaves two states it cannot tell apart, differing in ${confused}, at rates ${fmt(v.by, 3)} an hour apart`
          : `"${mech}" goes at rates ${fmt(v.by, 3)} an hour apart from two states it cannot tell apart, which differ in ${confused}`;
      const tr = h('tr', {}, h('td', { class: 'left' }, v.name), h('td', {}, String(v.blocks)), h('td', { class: 'left' }, h('span', { class: `chip ${kind}` }, verdict)),
        h('td', {}, v.horizon == null ? 'always' : hours(v.horizon, fmt)), h('td', { class: 'left faint' }, why));
      if (!v.eventwise && a && b) {
        tr.addEventListener('mouseenter', () => drawLaw([cellOf(a), cellOf(b)]));
        tr.addEventListener('mouseleave', () => drawLaw());
      }
      return tr;
    });
    const alone = r.views[0];
    const coarse = h('div', { class: 'section' }, h('h2', {}, 'What a coarser description forgets'),
      h('table', { class: 'data', style: 'width:100%' }, h('thead', {}, h('tr', {}, h('th', { class: 'left' }, 'description'), h('th', {}, 'states'), h('th', { class: 'left' }, 'verdict'), h('th', {}, 'within 10%'), h('th', { class: 'left' }, 'why'))), h('tbody', {}, rows)),
      ctx.finding({ title: alone.eventwise ? 'The stores alone can forget the institutions exactly' : 'The stores alone remember the lane', kind: 'proved',
        statement: `A description that keeps some parts and forgets others is a process of its own exactly when no change it keeps reads what it forgets. ${alone.eventwise ? 'Neither store answers the lane, so the two stores are a process of their own.' : `A store opens faster when the lane stands, so the stores alone are not a process of their own; the best they can do, averaging what they forget by how often it holds, stays within 10% of the whole for ${hours(alone.horizon, fmt)}.`} The permit cannot be forgotten, for building the lane reads it; nor the recognition, for while it holds the lane cannot be removed: a meaning adopted changes what may happen next. Counting the open stores without saying which is exact as a whole only when the two stores are alike, and never change by change.`,
        method: 'Strong lumpability tested for each mechanism and for the whole; the aggregation by the stationary law, its outward discrepancy δ, and the coupling bound 1 − exp(−δT) on its paths (physis lumping).' }),
      lines([{ name: 'both shut, the true law', x: r.dwell.tau, y: r.dwell.density, color: css('--accent-deep'), width: 2 }, { name: 'an exponential of the same mean', x: r.dwell.tau, y: r.dwell.exponential, color: css('--minor'), dash: '5 4' }],
        { xlabel: 'hours both stores stay shut', ylabel: 'density (logarithmic)', height: 170, logy: true }),
      legend([{ name: 'how long both stay shut', color: css('--accent-deep') }, { name: 'exponential, the same mean', color: css('--minor') }]),
      h('p', { class: 'caption' }, `Seen by the stores alone, a stretch with both shut lasts ${fmt(60 * r.dwell.mean, 3)} minutes on average; its spread (squared coefficient of variation ${fmt(r.dwell.cv2, 4)}, 1 for an exponential) is the memory of what the description forgot: ${Math.abs(r.dwell.cv2 - 1) < 1e-6 ? 'none' : 'whether the lane stood when the stretch began'}.`));

    // ---- interventions ----
    const outcomeKeys = [['both', 'both open'], ['open', 'stores open'], ['recognised', 'the lane recognised']];
    let outcome = 'both';
    const interBox = h('div', {});
    const drawInter = () => {
      const g = r.interventions[outcome];
      const main = [[1, 'the first'], [2, 'the second'], [4, 'permits']];
      const joint = [[3, 'both stores'], [5, 'first, permits'], [6, 'second, permits'], [7, 'all three']];
      const strongest = joint.reduce((m, j) => (Math.abs(g[j[0]]) > Math.abs(g[m[0]]) ? j : m), joint[0]);
      const jointName = { 3: 'the two promotions', 5: `promoting ${short(names[0])} with faster permits`, 6: `promoting ${short(names[1])} with faster permits`, 7: 'all three together' };
      interBox.replaceChildren(
        h('div', { class: 'grid cols-2' },
          h('div', {}, h('h3', {}, 'each alone'), h('div', { class: 'faint' }, `promoting ${short(names[0])}, ${short(names[1])}, or granting permits faster`), bars(main.map(([s], i) => ({ x: i, y: g[s], color: css('--series-1'), label: `${main[i][1]}: ${fmt(g[s], 3)}` })), { labels: main.map((m) => m[1]), height: 190, width: 300 })),
          h('div', {}, h('h3', {}, 'beyond each alone'), bars(joint.map(([s], i) => ({ x: i, y: g[s], color: g[s] >= 0 ? css('--series-4') : css('--series-5'), label: `${joint[i][1]}: ${fmt(g[s], 3)}` })), { labels: joint.map((j) => j[1]), height: 190, width: 300 }))),
        ctx.finding({ title: `On "${outcomeKeys.find((k) => k[0] === outcome)[1]}", the strongest interaction is ${jointName[strongest[0]]} (${fmt(g[strongest[0]], 3)})`, kind: 'proved',
          statement: `Each promotion multiplies its store's opening by ${fmt(1 + P.gain)}, and permits granted faster multiply the grant by as much. An interaction is what doing several together adds beyond doing each: for two, Y(both) − Y(one) − Y(other) + Y(neither); for three, the alternating sum over all eight. It belongs to the outcome, not to the stores: stores that answer nothing have none on how many are open, yet one on the chance both are open, which needs both. Here, with the outcome ${fmt(g[0], 4)} without any intervention, ${outcome === 'recognised' ? 'the promotions reach the lane\'s meaning only through its use, which the open stores make' : 'permits reach the stores only through the lane they answer'}.`,
          method: 'The anchored Möbius inversion over the subsets of the three interventions, each outcome at stationarity (morphe mobius).' }));
    };
    const keys = h('div', { class: 'keys' }, outcomeKeys.map(([k, name]) => h('button', { class: `key${k === outcome ? ' on' : ''}`, onclick: (e) => { outcome = k; for (const b of keys.children) b.classList.toggle('on', b === e.target); drawInter(); } }, name)));
    drawInter();
    const inter = h('div', { class: 'section' }, h('h2', {}, 'Three interventions and how they interact'),
      h('div', { class: 'panel' }, slider('gain', 'each intervention adds', 0.05, 2, 0.05)), keys, interBox);

    // ---- evidence ----
    const ev = r.evidence;
    const evidence = h('div', { class: 'section' }, h('h2', {}, 'What a record says of the stores\' answer to each other'),
      bars(ev.worlds.map((w, k) => ({ x: k, y: ev.posterior[k], color: w === 0 ? css('--series-3') : css('--accent'), label: `${fmt(w)}: ${fmt(100 * ev.posterior[k], 3)}%` })),
        { labels: ev.worlds.map((w) => `answers by ${fmt(w)}`), ylabel: 'probability after the record', height: 200 }),
      h('div', { class: 'panel' }, slider('hours', 'hours of record', 50, 3000, 50), slider('seed', 'which record', 1, 40, 1)),
      ctx.finding({ title: `After ${fmt(ev.hours)} hours: they answer each other with probability ${fmt(ev.coupled, 3)}, strongly (by 1 or more) ${fmt(ev.strong, 3)}`, kind: 'model',
        statement: `A record of ${ev.events} changes, made with the stores answering each other by ${fmt(P.etaA)} and ${fmt(P.etaB)}, weighs seven possible answers, alike for both stores, by how likely each makes the record. Overlapping propositions are weighed through the answers they hold for: "strongly" holds in fewer of them than "at all", and never has the larger probability.`,
        method: 'Each answer\'s exact path density of the record under the composed process; a uniform prior over the seven.',
        limits: 'A record made by the model itself, every change seen; a real one would leave the other rates to be learnt too, and some changes unseen.' }));

    out.replaceChildren(parts, law, coarse, inter, evidence);
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
  return h('div', { class: 'panel', style: 'margin:8px 0' }, s);
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
  return h('div', { class: 'panel', style: 'margin:8px 0' }, s,
    h('div', { class: 'faint' }, 'the parts (circles), where they meet (solid, with the relation their states must keep) and the changes that read several parts (diamonds, joined to the parts each reads)'));
}
