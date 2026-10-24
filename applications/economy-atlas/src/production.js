// Production: an input–output table as a circular flow (theory note §4). The table is edited here;
// value added balances each column (output by rows less the inputs). The engine gives the Leontief
// model (technical coefficients A, inverse L, multipliers, the price dual, linkages) and the
// decomposition of the net deliveries between sectors into what is driven and what circulates.

import * as api from './api.js';
import { el, svg, clear, fixed, percent, shareBar } from './dom.js';
import { PART, seriesColour } from './palette.js';
import { TABLES } from './presets.js';

const state = { table: structuredClone(TABLES[0].table), preset: TABLES[0].id, answer: null, shock: { sector: 0, dc: 0, dy: 0 } };

function balanced(t) {
  const n = t.sectors.length;
  const x = t.sectors.map((_, i) => t.z[i].reduce((s, v) => s + v, 0) + t.final_demand[i]);
  const value_added = t.sectors.map((_, j) => x[j] - t.z.reduce((s, row) => s + row[j], 0));
  return { sectors: t.sectors, z: t.z.map((r) => r.slice(0, n)), final_demand: t.final_demand, value_added };
}

function editor() {
  const t = state.table, n = t.sectors.length;
  const cell = (value, set) => el('input', { type: 'number', value, step: 'any', onchange: (e) => { set(Number(e.target.value)); run(); } });
  const rows = t.sectors.map((name, i) => el('tr', {},
    el('th', { class: 'row' }, el('input', { class: 'name', value: name, onchange: (e) => { t.sectors[i] = e.target.value; run(); } })),
    t.sectors.map((_, j) => el('td', {}, cell(t.z[i][j], (v) => { t.z[i][j] = v; }))),
    el('td', {}, cell(t.final_demand[i], (v) => { t.final_demand[i] = v; }))));
  return el('div', {},
    el('table', { class: 'grid' }, el('thead', {}, el('tr', {}, el('th', { class: 'row', text: 'from \\ to' }), t.sectors.map((s) => el('th', { text: s })), el('th', { text: 'final demand' }))), el('tbody', {}, rows)),
    el('div', { class: 'row-actions' },
      el('button', { onclick: () => { t.sectors.push(`sector ${n + 1}`); t.z.forEach((r) => r.push(0)); t.z.push(new Array(n + 1).fill(0)); t.final_demand.push(100); render(); run(); } }, '+ a sector'),
      n > 2 ? el('button', { onclick: () => { t.sectors.pop(); t.z.pop(); t.z.forEach((r) => r.pop()); t.final_demand.pop(); render(); run(); } }, '− the last sector') : null));
}

function matrix(M, digits = 4) {
  return el('table', { class: 'matrix' }, el('tbody', {}, M.map((row) => el('tr', {}, row.map((v) => el('td', { text: fixed(v, digits) }))))));
}

function circle(A) {
  const n = A.sectors.length, R = 120, cx = 200, cy = 160;
  const pos = A.sectors.map((_, i) => [cx + R * Math.cos((2 * Math.PI * i) / n - Math.PI / 2), cy + R * Math.sin((2 * Math.PI * i) / n - Math.PI / 2)]);
  const legs = A.deliveries.legs;
  const max = Math.max(1, ...legs.map((l) => l[2]));
  const g = svg('svg', { class: 'figure', width: 400, height: 320, viewBox: '0 0 400 320' },
    svg('defs', {}, svg('marker', { id: 'pr-head', viewBox: '0 0 10 10', refX: 8, refY: 5, markerWidth: 5, markerHeight: 5, orient: 'auto-start-reverse' }, svg('path', { d: 'M0,1L9,5L0,9z', fill: 'context-stroke' }))));
  for (const l of legs) {
    const [a, b] = [pos[l[0]], pos[l[1]]];
    const mx = (a[0] + b[0]) / 2 + (b[1] - a[1]) * 0.18, my = (a[1] + b[1]) / 2 - (b[0] - a[0]) * 0.18;
    const circ = Math.abs(l[4] + l[5]), drv = Math.abs(l[3]);
    g.append(svg('path', { d: `M${a[0]},${a[1]}Q${mx},${my} ${b[0]},${b[1]}`, fill: 'none', stroke: circ >= drv ? PART.local : PART.driven,
      'stroke-width': (1 + 6 * Math.sqrt(l[2] / max)).toFixed(1), 'stroke-opacity': 0.75, 'marker-end': 'url(#pr-head)' }, svg('title', { text: `${A.sectors[l[0]]} → ${A.sectors[l[1]]}: ${fixed(l[2], 1)} net, driven ${fixed(l[3], 1)}, circulating ${fixed(l[4] + l[5], 1)}` })));
  }
  A.sectors.forEach((s, i) => {
    g.append(svg('circle', { cx: pos[i][0], cy: pos[i][1], r: 7, fill: seriesColour(i) }));
    g.append(svg('text', { x: pos[i][0], y: pos[i][1] + (pos[i][1] < cy ? -13 : 22), 'text-anchor': 'middle', text: s }));
  });
  return g;
}

function results() {
  const box = document.getElementById('pr-panel');
  clear(box);
  const A = state.answer;
  if (!A) return;
  const d = A.deliveries;
  box.append(
    el('h2', { text: 'The Leontief model' }),
    el('p', { class: 'lede', text: A.productive ? 'The economy is productive: its Leontief inverse is non-negative, so any final demand can be met.' : 'The economy is not productive: its Leontief inverse has negative entries, and some final demand cannot be met.' }),
    el('h3', { text: 'Technical coefficients A' }), matrix(A.A),
    el('h3', { text: 'Leontief inverse L = (I − A)⁻¹' }), matrix(A.L),
    el('table', {}, el('thead', {}, el('tr', {}, el('th', { text: 'sector' }), el('th', { class: 'num', text: 'output' }), el('th', { class: 'num', text: 'multiplier' }), el('th', { class: 'num', text: 'price' }))),
      el('tbody', {}, A.sectors.map((s, i) => el('tr', {}, el('td', { text: s }), el('td', { class: 'num', text: fixed(A.output[i], 1) }), el('td', { class: 'num', text: fixed(A.multipliers[i], 3) }),
        el('td', { class: 'num', text: fixed(A.prices[i], 3) }))))),
    el('table', {}, el('thead', {}, el('tr', {}, el('th', { text: 'linkages' }), el('th', { class: 'num', text: 'backward' }), el('th', { class: 'num', text: 'forward' }))),
      el('tbody', {}, A.sectors.map((s, i) => el('tr', {}, el('td', { text: s }), el('td', { class: 'num', text: fixed(A.backward[i], 2) }), el('td', { class: 'num', text: fixed(A.forward[i], 2) }))))),
    el('p', { class: 'note', text: 'Multipliers: the output of all sectors one unit of final demand for a sector calls for. Prices: p = c L, which a table in money makes 1. Linkages over 1 are above the average (Rasmussen).' }),
    el('div', { class: 'section' }, el('h3', { text: 'The circular flow of deliveries' }),
      shareBar([{ label: 'reciprocal', value: d.reciprocal, colour: PART.reciprocal }, { label: 'net', value: d.net, colour: PART.driven }]),
      shareBar([{ label: 'driven', value: d.driven, colour: PART.driven }, { label: 'circulating', value: d.local + d.global, colour: PART.local }]),
      el('p', {}, 'Deliveries both ways between two sectors are exchange; of the net deliveries, the driven part runs from sectors that supply more than they use to those that use more than they are supplied, and the rest runs round the cycles of production.')),
  );
  if (A.shocked_prices) box.append(el('div', { class: 'section' }, el('h3', { text: 'Cost-push' }), el('p', {}, `With ${A.sectors[state.shock.sector]}'s value added per unit raised by ${state.shock.dc}: prices ${A.shocked_prices.map((p) => fixed(p, 3)).join(', ')}.`)));
  if (A.shocked_output) box.append(el('p', {}, `With ${A.sectors[state.shock.sector]}'s final demand raised by ${state.shock.dy}: outputs ${A.shocked_output.map((x) => fixed(x, 1)).join(', ')}.`));
}

function render() {
  const main = document.getElementById('pr-main');
  clear(main);
  const preset = TABLES.find((t) => t.id === state.preset);
  main.append(
    el('h2', { text: 'Production as a circular flow' }),
    el('p', { class: 'lede', text: 'An input–output table records what each sector delivers to each other in a period, and to final demand; value added balances each column. The Leontief model solves it, and the deliveries are read as a flow between sectors whose circulation the decomposition separates.' }),
    el('div', { class: 'presets' }, TABLES.map((t) => el('button', { class: t.id === state.preset ? 'active' : '', onclick: () => { state.preset = t.id; state.table = structuredClone(t.table); render(); run(); } }, t.title))),
    preset ? el('p', { class: 'note', text: preset.about }) : null,
    editor(),
    el('div', { class: 'row-actions' },
      el('label', {}, 'Shock sector ', el('select', { onchange: (e) => { state.shock.sector = Number(e.target.value); run(); } }, state.table.sectors.map((s, i) => el('option', { value: i, selected: i === state.shock.sector, text: s })))),
      el('label', {}, ' value added per unit +', el('input', { type: 'number', step: '0.01', value: state.shock.dc, style: { width: '70px' }, onchange: (e) => { state.shock.dc = Number(e.target.value); run(); } })),
      el('label', {}, ' final demand +', el('input', { type: 'number', step: '1', value: state.shock.dy, style: { width: '70px' }, onchange: (e) => { state.shock.dy = Number(e.target.value); run(); } }))),
    el('div', { id: 'pr-figure' }),
  );
}

async function run() {
  const t = balanced(state.table);
  const extra = {};
  try {
    if (state.answer && (state.shock.dc || state.shock.dy)) {
      if (state.shock.dc) extra.value_added_coefficients = state.answer.value_added_coefficients.map((c, i) => c + (i === state.shock.sector ? state.shock.dc : 0));
      if (state.shock.dy) extra.final_demand = t.final_demand.map((y, i) => y + (i === state.shock.sector ? state.shock.dy : 0));
    }
    document.getElementById('status').textContent = 'solving the table…';
    state.answer = await api.production(t, extra);
    document.getElementById('status').textContent = '';
  } catch (error) {
    document.getElementById('status').textContent = error.message;
    return;
  }
  const fig = document.getElementById('pr-figure');
  if (fig) clear(fig).append(circle(state.answer));
  results();
}

// The address may ask for an example table: #production&table=three.
export async function showProduction(data, params = {}) {
  const wanted = TABLES.find((t) => t.id === params.table);
  if (wanted && wanted.id !== state.preset) {
    state.preset = wanted.id;
    state.table = structuredClone(wanted.table);
    state.shock = { sector: 0, dc: 0, dy: 0 };
    render();
  }
  if (!document.getElementById('pr-main').firstChild) render();
  await run();
}
