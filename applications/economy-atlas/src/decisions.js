// Decisions: perception and decision on the economies' growth (ADR 0074). The World Bank's growth of real GDP,
// 1961–2024, read as a hidden regime (contraction, normal growth, boom) by a hidden Markov model fitted to every
// economy at once, each economy against its own normal; a policymaker's POMDP on that belief; planning online; a fund
// under a budget; a budget on support; two policymakers each seeing only its own growth; and agents that herd and
// sort. Every figure is computed by the engine (tyche's decision and processes, nomos's agents, the economics domain's
// regimes and policy); what policy does and is worth is declared, and shown with the figures.

import * as api from './api.js';
import { el, svg, put, fill } from './dom.js';
import { onTheme, token } from './palette.js';
import { lineChart, formatNumber as fmt } from './charts.js';
import { shared, emit } from './state.js';
import { namingPolicy } from '/visualization/js/map/index.js';

const SECTIONS = [
  ['perception', 'Perception', 'The regime of an economy, read from its growth'],
  ['decision', 'Decision', 'A policymaker acting on that belief'],
  ['search', 'Search', 'Planning online by Monte Carlo tree search'],
  ['fund', 'A fund', 'Support among economies under a year’s budget'],
  ['budget', 'A budget', 'A budget on support itself'],
  ['coordination', 'Two policymakers', 'Each seeing only its own growth'],
  ['herding', 'Herding', 'Kirman’s agents: two views, no preference'],
  ['sorting', 'Sorting', 'Schelling’s agents: mild preferences'],
];
const REGIME = ['contraction', 'normal growth', 'boom'];
const REGIME_COLOUR = () => [token('local'), '#A9B0B4', token('driven')];
const ACTION = ['wait', 'support', 'restrain'];
const ACTION_COLOUR = () => ['#DCDDD6', token('sends'), token('magnitude')];
// The interface's own colours (the canvas, the chosen), which the data's tokens do not hold.
const ui = (name, fallback) => (typeof document !== 'undefined'
  && getComputedStyle(document.documentElement).getPropertyValue(`--${name}`).trim()) || fallback;
const DEFAULTS = { welfare_recession: -4, welfare_normal: 1, welfare_boom: -0.5, cost_support: 0.15, cost_restrain: 0.2, support_effect: 0.5, restraint_effect: 0.35, discount: 0.9 };

const state = {
  section: 'perception', economy: 'USA', year: 2009, fundYear: 2020, budget: 20, supportBudget: 1, threshold: 1 / 3,
  assumptions: { ...DEFAULTS }, answers: new Map(), pending: new Map(), policy: namingPolicy(), names: {}, wired: false,
};

const nameOf = (code) => state.policy.name(code, state.names[code] || code);
const pct = (x, d = 0) => (Number.isFinite(x) ? `${(100 * x).toFixed(d)} %` : '—');

// ---- the engine's answers, kept ---------------------------------------------------------------------------------

function ask(kind, params = {}) {
  const key = JSON.stringify([kind, params]);
  if (state.answers.has(key)) return Promise.resolve(state.answers.get(key));
  if (state.pending.has(key)) return state.pending.get(key);
  const p = api.decide(kind, params).then((a) => { state.answers.set(key, a); state.pending.delete(key); return a; },
    (e) => { state.pending.delete(key); throw e; });
  state.pending.set(key, p);
  return p;
}
const policyParams = (extra = {}) => ({ ...state.assumptions, ...extra });

// ---- figures ------------------------------------------------------------------------------------------------------

function figure(title, W, H) {
  const host = svg('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': title });
  return { box: el('figure', { class: 'chart' }, el('figcaption', { text: title }), host), host };
}

const block = (title, ...children) => el('div', { class: 'tb-block' }, title ? el('div', { class: 'heading' }, el('span', { text: title })) : null, ...children);
const wide = (title, ...children) => el('div', { class: 'tb-block dc-wide' }, title ? el('div', { class: 'heading' }, el('span', { text: title })) : null, ...children);
const legend = (items) => el('div', { class: 'chart-legend' }, items.map(([label, colour]) => el('span', {}, el('i', { style: { background: colour } }), label)));
const figures = (pairs) => el('div', { class: 'figures' }, pairs.map(([v, label]) => el('div', { class: 'figure-cell' }, el('div', { class: 'value', text: v }), el('div', { class: 'label', text: label }))));

// An economy's years: growth as bars coloured by its most probable regime, and the filtered and smoothed chance of
// contraction as lines on the right axis.
function yearsFigure(r, highlight) {
  const W = 1120, H = 330, m = { l: 48, r: 60, t: 26, b: 34 };
  const { box, host } = figure(`${nameOf(r.economy)}: growth of real GDP, and the chance of contraction`, W, H);
  const n = r.growth.length, years = r.growth.map((_, t) => r.first + t);
  const g = r.growth.filter(Number.isFinite);
  const lo = Math.min(-4, ...g), hi = Math.max(6, ...g);
  const X = (t) => m.l + ((W - m.l - m.r) * (t + 0.5)) / n, bw = Math.max(1.5, (W - m.l - m.r) / n - 1.5);
  const Y = (v) => m.t + ((H - m.t - m.b) * (hi - v)) / (hi - lo), P = (p) => m.t + (H - m.t - m.b) * (1 - p);
  const colours = REGIME_COLOUR();
  for (let v = Math.ceil(lo / 2) * 2; v <= hi; v += 2) put(host, svg('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: v === 0 ? 'zero' : 'grid' }),
    svg('text', { x: m.l - 6, y: Y(v) + 4, 'text-anchor': 'end', class: 'tick', text: `${v}` }));
  for (const p of [0, 0.5, 1]) put(host, svg('text', { x: W - m.r + 6, y: P(p) + 4, class: 'tick', text: pct(p) }));
  r.growth.forEach((v, t) => {
    if (!Number.isFinite(v)) return;
    const y0 = Y(Math.max(0, v)), y1 = Y(Math.min(0, v));
    put(host, svg('rect', { x: X(t) - bw / 2, y: y0, width: bw, height: Math.max(0.5, y1 - y0), fill: colours[r.path[t]], 'fill-opacity': 0.85 },
      svg('title', { text: `${years[t]}: ${v.toFixed(2)} %; chance of contraction ${pct(r.filtered[t][0])} (filtered), ${pct(r.smoothed[t][0])} (smoothed)` })));
  });
  const line = (rows, dashed, colour) => {
    let d = '';
    rows.forEach((row, t) => { if (row) d += `${d ? 'L' : 'M'}${X(t).toFixed(1)},${P(row[0]).toFixed(1)}`; });
    put(host, svg('path', { d, fill: 'none', stroke: colour, 'stroke-width': dashed ? 1.4 : 2, 'stroke-dasharray': dashed ? '5 3' : null }));
  };
  line(r.smoothed, true, 'currentColor');
  line(r.filtered, false, colours[0]);
  for (let t = 0; t < n; t++) if (years[t] % 10 === 0) put(host, svg('text', { x: X(t), y: H - m.b + 16, 'text-anchor': 'middle', class: 'tick', text: `${years[t]}` }));
  if (highlight && highlight >= r.first && highlight < r.first + n) {
    const t = highlight - r.first;
    put(host, svg('line', { x1: X(t), x2: X(t), y1: m.t, y2: H - m.b, stroke: ui('accent', '#B08A3E'), 'stroke-width': 1.2 }));
  }
  put(host, svg('text', { x: 12, y: m.t - 10, class: 'axis-label', text: 'growth, %' }), svg('text', { x: W - m.r - 6, y: m.t - 10, 'text-anchor': 'end', class: 'axis-label', text: 'chance of contraction' }));
  put(box, legend([[REGIME[0], colours[0]], [REGIME[1], colours[1]], [REGIME[2], colours[2]], ['filtered (the years so far)', colours[0]], ['smoothed (the whole record), dashed', 'currentColor']]));
  return box;
}

// The share of the world's economies more likely than not in contraction, year by year.
function worldFigure(readings) {
  const first = Math.min(...readings.map((r) => r.first)), last = Math.max(...readings.map((r) => r.first + r.growth.length - 1));
  const x = [], share = [], mean = [];
  for (let y = first; y <= last; y++) {
    let n = 0, k = 0, s = 0;
    for (const r of readings) {
      const t = y - r.first;
      if (t < 0 || t >= r.growth.length || !Number.isFinite(r.growth[t])) continue;
      n++, s += r.filtered[t][0], k += r.filtered[t][0] > 0.5;
    }
    if (n >= 20) x.push(y), share.push(k / n), mean.push(s / n);
  }
  return lineChart({ title: 'The world’s contractions: economies more likely than not in contraction, and the mean chance', x,
    series: [{ name: 'share above ½', y: share }, { name: 'mean chance', y: mean, dashed: true }], xlabel: 'year', ylabel: 'share of economies' });
}

// The simplex of beliefs: the policy's action at each belief, and an economy's beliefs year by year.
function simplexFigure(grid, years) {
  const W = 520, H = 470, pad = 46;
  const { box, host } = figure('The policy over beliefs, and an economy’s beliefs year by year', W, H);
  const A = [pad, H - pad], B = [W - pad, H - pad], C = [W / 2, pad + 6];  // contraction, normal, boom
  const at = (b) => [b[0] * A[0] + b[1] * B[0] + b[2] * C[0], b[0] * A[1] + b[1] * B[1] + b[2] * C[1]];
  const n = grid.n, colours = ACTION_COLOUR();
  const point = new Map(grid.points.map((p) => [`${Math.round(p[0] * n)},${Math.round(p[1] * n)}`, p]));
  // Small triangles of the grid, each coloured by its first corner's action.
  for (let i = 0; i < n; i++)
    for (let j = 0; i + j < n; j++) {
      const tris = [[[i, j], [i + 1, j], [i, j + 1]]];
      if (i + j + 1 < n) tris.push([[i + 1, j], [i + 1, j + 1], [i, j + 1]]);
      for (const t of tris) {
        // The triangle's action: the commonest of its corners' (the first corner's on a tie).
        const acts = t.map(([a, b]) => point.get(`${a},${b}`)?.[3]).filter((v) => v !== undefined);
        const count = (v) => acts.filter((q) => q === v).length;
        const act = acts.reduce((best, v) => (count(v) > count(best) ? v : best), acts[0]);
        const pts = t.map(([a, b]) => at([a / n, b / n, (n - a - b) / n]));
        put(host, svg('polygon', { points: pts.map((q) => q.map((v) => v.toFixed(1)).join(',')).join(' '), fill: colours[act ?? 0], stroke: colours[act ?? 0], 'stroke-width': 0.4 }));
      }
    }
  put(host, svg('polygon', { points: [A, B, C].map((q) => q.join(',')).join(' '), fill: 'none', stroke: 'currentColor', 'stroke-width': 1 }));
  put(host, svg('text', { x: A[0] - 4, y: A[1] + 32, 'text-anchor': 'start', class: 'tick', text: 'certain of contraction' }),
    svg('text', { x: B[0] + 4, y: B[1] + 32, 'text-anchor': 'end', class: 'tick', text: 'certain of normal growth' }),
    svg('text', { x: C[0], y: C[1] - 10, 'text-anchor': 'middle', class: 'tick', text: 'certain of a boom' }));
  if (years) {
    const pts = years.years.map((y) => at(y[2]));
    put(host, svg('polyline', { points: pts.map((q) => q.map((v) => v.toFixed(1)).join(',')).join(' '), fill: 'none', stroke: 'currentColor', 'stroke-width': 0.8, 'stroke-opacity': 0.55 }));
    years.years.forEach((y, k) => {
      const [px, py] = pts[k];
      put(host, svg('circle', { cx: px, cy: py, r: 3.2, fill: ui('paper', '#FAFAF7'), stroke: 'currentColor', 'stroke-width': 1.2 },
        svg('title', { text: `${y[0]}: growth ${Number.isFinite(y[1]) ? y[1].toFixed(2) + ' %' : '—'}; belief ${y[2].map((v) => pct(v)).join(', ')}; ${ACTION[y[3]]}` })));
      // The years of contraction and recovery named: beneath the lower edge when the belief lies on it.
      if ([1982, 2009, 2020, 2021].includes(y[0])) {
        const low = y[2][2] < 0.05;
        put(host, svg('text', { x: low ? px : px + 7, y: low ? A[1] + 15 : py - 4, 'text-anchor': low ? 'middle' : 'start', class: 'tick', text: `${y[0]}` }));
      }
    });
  }
  put(box, legend(ACTION.map((a, k) => [a, colours[k]])));
  return box;
}

// An economy's years as a stripe: the action the policy takes on each year's belief.
function stripeFigure(years) {
  const W = 1120, H = 64, m = { l: 48, r: 60 };
  const { box, host } = figure(`${nameOf(years.economy)}: the action on each year’s belief`, W, H);
  const n = years.years.length, w = (W - m.l - m.r) / n, colours = ACTION_COLOUR();
  years.years.forEach((y, t) => {
    put(host, svg('rect', { x: m.l + t * w, y: 8, width: Math.max(1, w - 1), height: 26, fill: colours[y[3]] }, svg('title', { text: `${y[0]}: ${ACTION[y[3]]} (chance of contraction ${pct(y[2][0])})` })));
    if (y[0] % 10 === 0) put(host, svg('text', { x: m.l + (t + 0.5) * w, y: 52, 'text-anchor': 'middle', class: 'tick', text: `${y[0]}` }));
  });
  return box;
}

// Two policy trees of depth 3, each an action at the root and one subtree a signal (below or above normal).
function treesFigure(trees, title) {
  const W = 620, H = 210;
  const { box, host } = figure(title, W, H);
  const colours = ACTION_COLOUR();
  trees.forEach((tree, k) => {
    const x0 = 20 + k * (W / 2), w = W / 2 - 40;
    const pos = (node) => {
      let depth = 0, first = 0, width = 1;
      while (node >= first + width) first += width, width *= 2, depth++;
      const idx = node - first;
      return [x0 + (w * (idx + 0.5)) / width, 34 + depth * 64];
    };
    put(host, svg('text', { x: x0 + w / 2, y: 14, 'text-anchor': 'middle', class: 'axis-label', text: `policymaker ${k + 1}` }));
    tree.forEach((a, node) => {
      const [px, py] = pos(node);
      for (const o of [0, 1]) {
        const child = node * 2 + 1 + o;
        if (child < tree.length) {
          const [cx, cy] = pos(child);
          put(host, svg('line', { x1: px, y1: py, x2: cx, y2: cy, stroke: 'currentColor', 'stroke-opacity': 0.35 }),
            svg('text', { x: (px + cx) / 2 + (o ? 6 : -6), y: (py + cy) / 2, 'text-anchor': o ? 'start' : 'end', class: 'tick', text: o ? 'above' : 'below' }));
        }
      }
      put(host, svg('circle', { cx: px, cy: py, r: 11, fill: colours[a], stroke: 'currentColor', 'stroke-width': 0.8 }),
        svg('text', { x: px, y: py + 4, 'text-anchor': 'middle', class: 'tick', fill: a ? '#FFFFFF' : null, text: a ? 'S' : 'W' }));
    });
  });
  put(box, legend([['W: wait', colours[0]], ['S: support', colours[1]]]));
  return box;
}

// Schelling's frames: small multiples of the lattice.
function framesFigure(answer) {
  const side = answer.side, cell = 4, gap = 14, per = side * cell;
  const W = answer.frames.length * (per + gap) + gap, H = per + 40;
  const { box, host } = figure(`Schelling’s lattice, ${side} × ${side}, ${answer.agents} agents, threshold ${fmt(answer.threshold)}, tick by tick`, W, H);
  const colours = [token('local'), token('driven')];
  answer.frames.forEach((f, k) => {
    const x0 = gap + k * (per + gap);
    put(host, svg('rect', { x: x0, y: 8, width: per, height: per, fill: ui('paper-2', '#F3F3EF') }));
    f.cells.forEach((c, i) => {
      if (c < 0) return;
      put(host, svg('rect', { x: x0 + (i % side) * cell, y: 8 + Math.floor(i / side) * cell, width: cell, height: cell, fill: colours[c] }));
    });
    put(host, svg('text', { x: x0 + per / 2, y: per + 26, 'text-anchor': 'middle', class: 'tick', text: `tick ${f.tick}` }));
  });
  return box;
}

// Kirman's law: the share of time with k agents of view 1, simulated (bars) against the exact law (line).
function lawFigure(answer) {
  const W = 640, H = 280, m = { l: 52, r: 20, t: 18, b: 36 };
  const N = answer.agents, { box, host } = figure(`The law of the number holding view 1: over a run of ${fmt(answer.long_run)} (bars) and exact (line), ${N} agents`, W, H);
  const top = Math.max(...answer.seen, ...answer.exact);
  const X = (k) => m.l + ((W - m.l - m.r) * (k + 0.5)) / (N + 1), bw = (W - m.l - m.r) / (N + 1) - 1;
  const Y = (v) => m.t + ((H - m.t - m.b) * (top - v)) / top;
  for (const v of [0, top / 2, top]) put(host, svg('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: v === 0 ? 'zero' : 'grid' }),
    svg('text', { x: m.l - 6, y: Y(v) + 4, 'text-anchor': 'end', class: 'tick', text: fmt(v) }));
  answer.seen.forEach((v, k) => put(host, svg('rect', { x: X(k) - bw / 2, y: Y(v), width: bw, height: Y(0) - Y(v), fill: token('sends'), 'fill-opacity': 0.55 })));
  put(host, svg('path', { d: answer.exact.map((v, k) => `${k ? 'L' : 'M'}${X(k).toFixed(1)},${Y(v).toFixed(1)}`).join(''), fill: 'none', stroke: token('local'), 'stroke-width': 2 }));
  for (const k of [0, Math.round(N / 2), N]) put(host, svg('text', { x: X(k), y: H - m.b + 16, 'text-anchor': 'middle', class: 'tick', text: `${k}` }));
  put(host, svg('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', class: 'axis-label', text: 'agents holding view 1' }));
  return box;
}

// ---- the sections -------------------------------------------------------------------------------------------------

async function perception() {
  const R = await ask('regimes');
  state.names = R.names || {};
  const m = R.model, readings = R.readings;
  const r = readings.find((x) => x.economy === state.economy) || readings[0];
  const options = [...readings].sort((a, b) => nameOf(a.economy).localeCompare(nameOf(b.economy)));
  const choose = el('select', { onchange: (e) => { state.economy = e.target.value; render(); } },
    options.map((x) => el('option', { value: x.economy, selected: x.economy === r.economy ? true : null, text: nameOf(x.economy) })));
  const colours = REGIME_COLOUR();
  const table = el('table', { class: 'results' }, el('thead', {}, el('tr', {}, ['regime', 'mean growth', 'spread', 'mean duration', 'share of years'].map((h) => el('th', { text: h })))),
    el('tbody', {}, REGIME.map((name, k) => el('tr', {}, el('td', {}, el('i', { class: 'swatch', style: { background: colours[k] } }), name),
      el('td', { text: `${fmt(m.mean[k])} σ` }), el('td', { text: `${fmt(m.sd[k])} σ` }), el('td', { text: `${fmt(m.duration[k])} years` }), el('td', { text: pct(m.stationary[k], 1) })))));
  const trans = el('table', { class: 'results' }, el('thead', {}, el('tr', {}, el('th', { text: 'from \\ to' }), REGIME.map((h) => el('th', { text: h })))),
    el('tbody', {}, m.transition.map((row, i) => el('tr', {}, el('th', { text: REGIME[i] }), row.map((v) => el('td', { text: pct(v, 1) }))))));
  const year = state.year, layer = () => {
    const values = {};
    for (const x of readings) {
      const t = year - x.first;
      if (t >= 0 && t < x.growth.length && Number.isFinite(x.growth[t])) values[x.economy] = x.filtered[t][0];
    }
    putLayer(`decide-contraction-${year}`, `The chance of contraction in ${year}, filtered`, values, 'sequential', 'probability');
  };
  return [
    wide('', el('div', { class: 'dc-controls' }, el('label', {}, 'economy ', choose),
      el('label', {}, 'year ', el('input', { type: 'number', min: 1962, max: 2024, value: year, onchange: (e) => { state.year = +e.target.value; render(); } })),
      el('button', { class: 'quiet', text: `the chance of contraction in ${year} on the map`, onclick: layer })), yearsFigure(r, year)),
    block('The regimes, fitted to every economy at once',
      figures([[`${m.economies}`, 'economies'], [`${m.observations.toLocaleString('en-GB')}`, 'years of growth'], [`${m.iterations}`, 'iterations of EM']]),
      table, el('p', { class: 'note', text: 'In each economy’s own deviations σ from its median growth (a robust deviation, 1.4826 × the median absolute deviation); one σ for the three regimes, as in Hamilton’s model, so that they differ by their means alone.' })),
    block('How the regimes follow one another', trans, el('p', { class: 'note', text: 'A contraction rarely outlasts its year: the chance that it continues is the first row’s first entry.' })),
    block('', worldFigure(readings)),
  ];
}

function assumptionsCard(onChange) {
  const A = state.assumptions;
  const row = (key, label, step = 0.05) => el('label', { class: 'dc-param' }, el('span', { text: label }),
    el('input', { type: 'number', step, value: A[key], onchange: (e) => { A[key] = +e.target.value; onChange(); } }));
  return block('What policy does and is worth (declared)',
    el('div', { class: 'dc-params' }, row('welfare_recession', 'a year of contraction'), row('welfare_normal', 'a year of normal growth'), row('welfare_boom', 'a year of boom'),
      row('cost_support', 'the cost of support'), row('cost_restrain', 'the cost of restraint'), row('support_effect', 'support: share of a contraction’s persistence ended'),
      row('restraint_effect', 'restraint: share of a boom’s persistence ended'), row('discount', 'the discount a year', 0.01)),
    el('div', { class: 'row-actions' }, el('button', { class: 'quiet', text: 'the defaults again', onclick: () => { state.assumptions = { ...DEFAULTS }; onChange(); } })),
    el('p', { class: 'note', text: 'The data give the regimes, how they follow one another and how each shows in growth. What support and restraint do, and what a year is worth, are assumptions, and every figure here carries them.' }));
}

async function decision() {
  const P = await ask('policy', policyParams({ economy: state.economy }));
  const h = P.hsvi, colours = ACTION_COLOUR();
  const edge = P.edges[0];
  const x = edge.points.map((p) => p[0]);
  const slice = lineChart({ title: 'The value along the edge from contraction to normal growth: HSVI’s bounds and three bounds from the states', x,
    series: [{ name: 'HSVI, upper', y: edge.points.map((p) => p[2]), dashed: true }, { name: 'HSVI, lower', y: edge.points.map((p) => p[1]) },
      { name: 'fast informed bound', y: edge.points.map((p) => p[3]), dashed: true }, { name: 'QMDP', y: edge.points.map((p) => p[4]), dashed: true },
      { name: 'blind policies', y: edge.points.map((p) => p[5]) }], xlabel: 'the belief in normal growth (the rest in contraction)', ylabel: 'value' });
  const trace = lineChart({ title: 'HSVI closing on the optimal value at the starting belief, trial by trial', x: h.trace.map((t) => t[0]),
    series: [{ name: 'upper bound', y: h.trace.map((t) => t[2]) }, { name: 'lower bound', y: h.trace.map((t) => t[1]) }], xlabel: 'trial', ylabel: 'value' });
  return [
    wide('', el('div', { class: 'dc-controls' }, el('span', { class: 'note', text: `${nameOf(state.economy)}: its beliefs after each year are drawn over the policy; choose the economy under Perception.` })),
      el('div', { class: 'dc-pair' }, simplexFigure(P.grid, P.economy), el('div', {},
        figures([[`${fmt(h.lower)} … ${fmt(h.upper)}`, 'the optimal value at an even belief'], [`${h.trials}`, 'trials of HSVI'], [`${h.vectors.length}`, 'plans in its lower bound']]),
        el('p', { class: 'lede', text: `Each point of the triangle is a belief: how sure the policymaker is of contraction, normal growth or a boom. Its colour is the action the solved POMDP takes there. Near certainty of contraction it supports; near certainty of a boom it restrains; between, it waits. The dots are the economy’s beliefs after each year’s growth, read by the regimes’ filter.` }),
        trace))),
    P.economy ? wide('', stripeFigure(P.economy)) : null,
    block('', slice),
    assumptionsCard(() => render()),
  ];
}

async function search() {
  const S = await ask('online', policyParams({ economy: state.economy, year: state.year }));
  const x = S.ladder.map((l) => Math.log10(l.simulations));
  const series = ACTION.map((a, k) => ({ name: a, y: S.ladder.map((l) => l.q[k]) }));
  const chart = lineChart({ title: `${nameOf(S.economy)} after ${S.year}: POMCP’s action values against the number of simulations`, x, series, xlabel: 'simulations (log₁₀)', ylabel: 'value' });
  const rows = ACTION.map((a, k) => el('tr', {}, el('td', { text: a }), el('td', { text: `${fmt(S.q_lower[k])} … ${fmt(S.q_upper[k])}` }), el('td', { text: fmt(S.ladder.at(-1).q[k]) })));
  return [
    block('', chart),
    block(`${nameOf(S.economy)} after ${S.year}: the belief planned from`, figures([[pct(S.belief[0]), 'contraction'], [pct(S.belief[1]), 'normal growth'], [pct(S.belief[2]), 'boom']]),
      el('table', { class: 'results' }, el('thead', {}, el('tr', {}, ['action', 'HSVI’s interval', 'POMCP, 30 000 simulations'].map((t) => el('th', { text: t })))), el('tbody', {}, rows)),
      el('p', { class: 'note', text: `HSVI takes ${ACTION[S.hsvi_action]}; POMCP, by the end, ${ACTION[S.ladder.at(-1).action]}. POMCP plans from particles by sampling, without the tables: its values approach the interval from below, its leaf the value of the best blind policy.` })),
  ];
}

function putLayer(key, title, values, scale, unit) {
  shared.layers = shared.layers.filter((l) => l.key !== key);
  shared.layers.push({ key, file: `decisions:${key}`, title, plugin: 'Decisions', language: 'tool', year: undefined, dataset: 'library:trade',
    result: { kind: 'values', title, values, scale, unit } });
  shared.colour = key;
  emit('layers');
  document.getElementById('status').textContent = `${title}: on the map`;
}

async function fund() {
  const F = await ask('fund', policyParams({ year: state.fundYear, budget: state.budget }));
  const supported = F.economies.filter((e) => e[3] === 1);  // in the order of choice
  const curve = lineChart({ title: 'The Lagrangian bound on what any allocation within the budget can reach, against the multiplier λ', x: F.curve.map((c) => c[0]),
    series: [{ name: 'the bound', y: F.curve.map((c) => c[1]) }], xlabel: 'λ, the price of a unit of the budget', ylabel: 'discounted welfare',
    marks: [{ x: F.lambda, label: `the least bound, λ = ${fmt(F.lambda)}` }] });
  const tie = F.sure > F.supported
    ? ` ${F.sure} economies were as sure of contraction after ${F.year}, and the model cannot tell them apart: among equal gains the budget goes to the deepest falls, in each economy’s own deviations σ, a declared rule.`
    : '';
  const layer = () => {
    const values = {};
    for (const e of F.economies) values[e[0]] = e[3];
    putLayer(`decide-fund-${F.year}`, `The fund’s support after ${F.year}`, values, 'sequential', '1 supported, 0 not');
  };
  return [
    wide('', el('div', { class: 'dc-controls' },
      el('label', {}, 'after the year ', el('input', { type: 'number', min: 1962, max: 2024, value: state.fundYear, onchange: (e) => { state.fundYear = +e.target.value; render(); } })),
      el('label', {}, 'economies the budget supports ', el('input', { type: 'number', min: 1, max: 100, value: state.budget, onchange: (e) => { state.budget = +e.target.value; render(); } })),
      el('button', { class: 'quiet', text: 'the allocation on the map', onclick: layer })),
      figures([[`${F.supported} of ${F.economies.length}`, 'economies supported'], [fmt(F.lambda), 'λ at the least bound'], [fmt(F.upper), 'the bound'], [fmt(F.waiting), 'waiting everywhere (exact)']])),
    block(`Supported after ${F.year}, in the order chosen`, el('ol', { class: 'dc-list' }, supported.map((e) => el('li', {}, el('span', { text: nameOf(e[0]) }),
      el('span', { class: 'note', text: `${pct(e[1])} in contraction; growth ${Number.isFinite(e[4]) ? e[4].toFixed(1) + ' %' : '—'}, ${Number.isFinite(e[5]) ? e[5].toFixed(1) + ' σ' : '—'}` })))),
      el('p', { class: 'note', text: `Each economy is one POMDP of the same model, from its own belief; the budget joins them. The bound relaxes the budget into a price λ, solving each economy alone; at the least bound, the economies whose support gains most for its cost are supported until the budget is spent.${tie}` })),
    block('', curve),
  ];
}

async function budget() {
  const B = await ask('budget', policyParams({ economy: state.economy, year: state.year, budget: state.supportBudget }));
  const W = 560, H = 330, m = { l: 58, r: 20, t: 34, b: 40 };
  const { box, host } = figure(`What support buys: the plans for each price λ of support, ${nameOf(B.economy)} after ${B.year}`, W, H);
  // The plans of the price path, with the two the budget mixes, by their support; repeated plans drawn once.
  const pts = [...B.frontier, [B.lambda, ...B.low], [B.lambda, ...B.high]].sort((a, b) => a[2] - b[2])
    .filter((p, i, all) => i === 0 || Math.abs(p[2] - all[i - 1][2]) > 1e-6);
  const cs = [...pts.map((p) => p[2]), B.budget, B.cost], rs = [...pts.map((p) => p[1]), B.reward, B.dual];
  const c0 = 0, c1 = Math.max(...cs) * 1.08 || 1, r0 = Math.min(...rs), r1 = Math.max(...rs);
  const X = (c) => m.l + ((W - m.l - m.r) * (c - c0)) / (c1 - c0), Y = (r) => m.t + ((H - m.t - m.b) * (r1 - r)) / (r1 - r0 || 1);
  put(host, svg('line', { x1: X(B.budget), x2: X(B.budget), y1: m.t, y2: H - m.b, stroke: ui('accent', '#B08A3E'), 'stroke-dasharray': '4 3' }),
    svg('text', { x: X(B.budget) + 4, y: H - m.b - 8, class: 'tick', text: `the budget: ${fmt(B.budget)}` }));
  put(host, svg('polyline', { points: pts.map((p) => `${X(p[2]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(' '), fill: 'none', stroke: token('sends'), 'stroke-width': 1.5 }));
  pts.forEach((p) => put(host, svg('circle', { cx: X(p[2]), cy: Y(p[1]), r: 3, fill: token('sends') }, svg('title', { text: `λ ${fmt(p[0])}: welfare ${fmt(p[1])}, support ${fmt(p[2])}` }))));
  put(host, svg('circle', { cx: X(B.cost), cy: Y(B.reward), r: 6, fill: 'none', stroke: token('local'), 'stroke-width': 2 }),
    svg('text', { x: X(B.cost) - 9, y: Y(B.reward) - 9, 'text-anchor': 'end', class: 'tick', text: 'the mixture' }));
  put(host, svg('line', { x1: m.l, x2: W - m.r, y1: Y(B.dual), y2: Y(B.dual), stroke: 'currentColor', 'stroke-dasharray': '2 3', 'stroke-opacity': 0.6 }),
    svg('text', { x: m.l + 6, y: Y(B.dual) - 5, class: 'tick', text: 'the dual bound' }));
  for (const c of [0, c1 / 2, c1]) put(host, svg('text', { x: X(c), y: H - m.b + 16, 'text-anchor': 'middle', class: 'tick', text: fmt(c) }));
  for (const r of [r0, (r0 + r1) / 2, r1]) put(host, svg('text', { x: m.l - 6, y: Y(r) + 4, 'text-anchor': 'end', class: 'tick', text: fmt(r) }));
  put(host, svg('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', class: 'axis-label', text: 'expected discounted years of support' }),
    svg('text', { x: 12, y: m.t - 18, class: 'axis-label', text: 'discounted welfare' }));
  return [
    wide('', el('div', { class: 'dc-controls' }, el('label', {}, 'the budget, in discounted years of support ', el('input', { type: 'number', step: 0.25, min: 0, value: state.supportBudget, onchange: (e) => { state.supportBudget = +e.target.value; render(); } })),
      el('span', { class: 'note', text: `${nameOf(B.economy)} after ${B.year}: the belief the plans start from` })),
      el('div', { class: 'dc-pair' }, box, el('div', {}, figures([[fmt(B.reward), 'the mixture’s welfare'], [fmt(B.cost), 'its support'], [fmt(B.dual), 'the dual bound'], [`${fmt(B.free[0])} … ${fmt(B.free[1])}`, 'without the budget']]),
        el('p', { class: 'lede', text: `Each point is the best plan when a year of support is priced at λ: dearer support, less of it. Where the budget binds, at λ = ${fmt(B.lambda)}, the policy draws at the start between the two plans on either side of the budget, the more generous with probability ${pct(B.weight, 1)}, and spends the budget exactly in expectation. Each plan is a finite-state controller whose welfare and support are computed exactly.` })))),
  ];
}

async function coordination() {
  const C = await ask('coordination');
  return [
    wide('', el('div', { class: 'dc-pair' }, treesFigure(C.optimum.trees, 'The best joint policy over three years: each policymaker’s tree'), el('div', {},
      figures([[fmt(C.optimum.value), 'the best joint policy'], [fmt(C.jesp.value), 'JESP’s equilibrium'], [fmt(C.alone), 'the best of one, the other waiting'], [fmt(C.never), 'neither ever supporting']]),
      el('p', { class: 'lede', text: `Two neighbouring economies, each in contraction with chance ${pct(C.shock)} after a shared shock. Each policymaker sees only whether its own growth fell below its normal, and both share one welfare. Support by both shortens both contractions by a further ${pct(C.spillover)}, a declared spillover. The best joint policy is found by searching every tree of the one against the other’s exact best response; JESP, alternating best responses, reaches an equilibrium.${Math.abs(C.alone - C.never) < 1e-9 && C.optimum.value > C.never ? ' Alone, neither finds support worth its cost; together, at once after the shock, they do.' : ''}` }),
      el('p', { class: 'note', text: `From the regimes: a contraction continues with chance ${pct(C.recession_stay)}, begins with chance ${pct(C.recession_from_rest)}; growth falls below normal with chance ${pct(C.low_in_recession)} in contraction and ${pct(C.low_otherwise)} otherwise.` })))),
  ];
}

async function herding() {
  const H = await ask('herding');
  const x = H.series.map((_, t) => t);
  return [
    block('', lineChart({ title: `Kirman’s market: the share of ${H.agents} agents holding view 1, a unit of time apart`, x, series: [{ name: 'share holding view 1', y: H.series }], xlabel: 'time', ylabel: 'share' })),
    block('', lawFigure(H)),
    block('Herding without a preference', figures([[fmt(H.epsilon), 'ε, changing one’s mind alone'], [fmt(H.delta), 'δ, recruitment'], [fmt(H.threshold), 'δ / (N − 1)']]),
      el('p', { class: 'lede', text: 'No agent prefers either view, yet the market spends most of its time near one or the other and switches suddenly: ε lies below δ / (N − 1), where the law of the number holding a view, exact for this birth–death chain, is bimodal. Optimism and pessimism in a market, or two suppliers, herd so (Kirman 1993).' })),
  ];
}

async function sorting() {
  const S = await ask('sorting', { threshold: state.threshold });
  const x = S.segregation.map((_, t) => t);
  return [
    wide('', el('div', { class: 'dc-controls' }, el('label', {}, 'the threshold of like neighbours ', el('input', { type: 'number', step: 0.05, min: 0, max: 1, value: state.threshold.toFixed(2), onchange: (e) => { state.threshold = +e.target.value; render(); } }))),
      framesFigure(S)),
    block('', lineChart({ title: 'Like neighbours, and the share of agents content, tick by tick', x, series: [{ name: 'like neighbours', y: S.segregation }, { name: 'content', y: S.content }], xlabel: 'tick', ylabel: 'share' })),
    block('Sorting from mild preferences', el('p', { class: 'lede', text: `Each agent moves when fewer than ${pct(S.threshold)} of its neighbours are of its kind. From a random start, about half of each agent’s neighbours are like it; within a few ticks most are. Residential sorting by income, origin or tenure may follow so from preferences far milder than the pattern suggests (Schelling 1971).` })),
  ];
}

const RENDER = { perception, decision, search, fund, budget, coordination, herding, sorting };

function renderRail() {
  fill(document.getElementById('dc-rail'), el('div', { class: 'tb-group' }, el('div', { class: 'heading' }, el('span', { text: 'perception and decision' })),
    SECTIONS.map(([id, title, sub]) => el('button', { class: `tb-tool${state.section === id ? ' on' : ''}`, onclick: () => { state.section = id; history.replaceState(null, '', `#decisions&section=${id}`); renderRail(); render(); } },
      el('span', {}, el('span', { text: title }), el('br'), el('span', { class: 'note', text: sub }))))),
  el('p', { class: 'note dc-rail-note', text: 'The World Bank’s growth of real GDP, 1961–2024. The regimes are fitted; what policy does and is worth is declared. Engine: tyche (decision, processes), nomos (agents), the economics domain (regimes, policy).' }));
}

let rendering = 0;
async function render() {
  const main = document.getElementById('dc-main');
  const [, title, sub] = SECTIONS.find((s) => s[0] === state.section);
  const ticket = ++rendering;
  fill(main, el('div', { class: 'tb-head' }, el('h2', { text: title }), el('p', { class: 'lede', text: sub })), el('p', { class: 'note', text: 'computing…' }));
  try {
    if (state.section !== 'perception' && !Object.keys(state.names).length) { const R = await ask('regimes'); state.names = R.names || {}; }
    const blocks = await RENDER[state.section]();
    if (ticket !== rendering) return;
    fill(main, el('div', { class: 'tb-head' }, el('h2', { text: title }), el('p', { class: 'lede', text: sub })), el('div', { class: 'tb-answer dc-answer' }, blocks));
  } catch (error) {
    if (ticket !== rendering) return;
    fill(main, el('div', { class: 'tb-head' }, el('h2', { text: title })), el('p', { class: 'note', text: error.message }));
  }
}

// The address may name a section: #decisions&section=fund, and an economy: &economy=GBR.
export async function showDecisions(data, params = {}) {
  if (!state.wired) {
    state.wired = true;
    onTheme(() => render());
    try { state.policy = namingPolicy(await api.names()); } catch { /* the names as given */ }
  }
  if (params.section && RENDER[params.section]) state.section = params.section;
  if (params.economy) state.economy = params.economy;
  if (params.year) state.year = +params.year;
  renderRail();
  await render();
}
