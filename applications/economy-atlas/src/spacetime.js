// Spacetime: currencies held through time at their policy rates and exchanged with a base each
// month (the euro, for the library's record; one's own record names its own), as the economics
// domain's spacetime complex (theory note §5).
//
// Each lens of the complex, the base and one currency between two consecutive months, has the
// curvature κ = (r_B − r_X) Δt + log S_q − log S_p, S the units of X per unit of the base: the excess
// return of holding the base and exchanging later over exchanging at once and holding X, the departure from
// uncovered interest parity in that month. A sheet shows every lens; lines show the curvature
// summed in time, a currency's cumulative excess return against the base; the panel measures the
// curvature against the exchange rates' alone (κ₀ = log S_q − log S_p, nothing earned by holding):
// the share of it the rates of holding remove. Curvature is invariant under every re-description
// of value, unlike the share of the transports in their exact part, which a change of a
// currency's unit alters, and is therefore no measure.

import * as api from './api.js';
import { shared, datasetsOf, dataset, on } from './state.js';
import { openDataDialog, describe } from './data.js';
import { el, svg, clear, put, fill, percent, fixed, signed } from './dom.js';
import { diverging, seriesColour, onTheme } from './palette.js';

const state = { answer: null, selected: null, months: [], wired: false };

// How the base is named: the euro for the library's record, its code otherwise.
const base = (A) => (A.base || 'EUR') === 'EUR' ? { name: 'the euro', its: "the euro's", holding: 'euros', per: 'per euro' }
  : { name: A.base, its: `${A.base}'s`, holding: A.base, per: `per ${A.base}` };

function quantile(values, q) {
  const v = [...values].sort((a, b) => a - b);
  return v.length ? v[Math.min(v.length - 1, Math.floor(q * v.length))] : 0;
}

function sheet() {
  const A = state.answer;
  const host = document.getElementById('st-sheet');
  clear(host);
  const rows = A.currencies.length - 1, cols = A.months.length - 1;
  const box = host.getBoundingClientRect();
  const W = Math.max(600, box.width), left = 92, right = 74, top = 26, rowH = Math.max(14, Math.min(24, (box.height - top - 30) / rows));
  const cw = (W - left - right) / cols;
  host.setAttribute('viewBox', `0 0 ${W} ${top + rows * rowH + 30}`);
  host.style.height = rows * 24 + 56 < 520 ? `${top + rows * rowH + 30}px` : '';
  const scale = Math.max(1e-6, quantile(A.lenses.map((l) => Math.abs(l[3])), 0.97));
  const byCurrency = new Map();
  for (const l of A.lenses) {
    if (!byCurrency.has(l[0])) byCurrency.set(l[0], []);
    byCurrency.get(l[0]).push(l);
  }
  // Years along the top.
  for (let k = 0; k < A.months.length; k++)
    if (A.months[k].endsWith('-01')) put(host, svg('text', { x: left + k * cw, y: top - 8, text: A.months[k].slice(0, 4) }));
  const stats = new Map(A.stats.map((s) => [s.currency, s]));
  for (let r = 0; r < rows; r++) {
    const x = r + 1, y = top + r * rowH;
    const s = stats.get(x);
    put(host, svg('text', { class: 'row-label', x: left - 8, y: y + rowH * 0.7, 'text-anchor': 'end', text: A.currencies[x],
      style: state.selected === x ? 'font-weight:600' : '', onclick: () => select(x) }));
    for (const l of byCurrency.get(x) || [])
      put(host, svg('rect', { class: 'cell', x: (left + l[1] * cw).toFixed(1), y: (y + 1).toFixed(1), width: Math.max(0.6, cw * (l[2] - l[1]) - 0.4).toFixed(2), height: (rowH - 2).toFixed(1),
        fill: diverging(l[3] / scale), onclick: () => select(x) },
        svg('title', { text: `${A.currencies[x]}, ${A.months[l[1]]} to ${A.months[l[2]]}: κ = ${signed(100 * l[3], 2)}% (${base(A).its} excess return)` })));
    if (s) put(host, svg('text', { x: W - right + 8, y: y + rowH * 0.7, text: `${signed(100 * s.mean, 1)}%/yr` }));
  }
  put(host, svg('text', { x: left, y: top + rows * rowH + 20, text: `curvature of each lens, ${base(A).its} excess return over the month (colour scale ±${(100 * scale).toFixed(1)}%); right: its mean a year` }));
}

function lines() {
  const A = state.answer;
  const host = document.getElementById('st-lines');
  clear(host);
  const box = host.getBoundingClientRect();
  const W = Math.max(600, box.width), H = Math.max(220, box.height), left = 92, right = 74, top = 14, bottom = 26;
  host.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const all = A.stats.flatMap((s) => s.cumulative);
  const lo = Math.min(0, ...all), hi = Math.max(0, ...all);
  const x = (k) => left + ((W - left - right) * k) / (A.months.length - 1);
  const y = (v) => top + ((H - top - bottom) * (hi - v)) / Math.max(1e-9, hi - lo);
  put(host, svg('path', { d: `M${left},${y(0).toFixed(1)}H${W - right}`, stroke: 'var(--hair-2)' }));
  put(host, svg('text', { x: left - 8, y: y(0) + 4, 'text-anchor': 'end', text: '0' }));
  put(host, svg('text', { x: left - 8, y: y(hi) + 10, 'text-anchor': 'end', text: signed(hi, 2) }));
  put(host, svg('text', { x: left - 8, y: y(lo), 'text-anchor': 'end', text: signed(lo, 2) }));
  for (const s of A.stats) {
    const c = s.currency;
    let d = `M${x(0).toFixed(1)},${y(0).toFixed(1)}`;
    s.cumulative.forEach((v, k) => { d += `L${x(k + 1).toFixed(1)},${y(v).toFixed(1)}`; });
    put(host, svg('path', { class: `line${state.selected === c ? ' selected' : ''}`, d, stroke: seriesColour(c - 1),
      'stroke-opacity': state.selected === null || state.selected === c ? 0.9 : 0.25, onclick: () => select(c) }, svg('title', { text: A.currencies[c] })));
    const end = s.cumulative.at(-1);
    put(host, svg('text', { x: W - right + 6, y: y(end) + 4, text: A.currencies[c], fill: seriesColour(c - 1), 'fill-opacity': state.selected === null || state.selected === c ? 1 : 0.4 }));
  }
  put(host, svg('text', { x: left, y: H - 6, text: `the curvature summed in time: the cumulative log excess return of ${base(A).name} over each currency (Stokes: the deficit round the loop of all its lenses)` }));
}

function select(c) {
  state.selected = state.selected === c ? null : c;
  sheet();
  lines();
  panel();
}

function panel() {
  const A = state.answer;
  const box = document.getElementById('st-panel');
  clear(box);
  const K = A.curvature;
  put(box, 
    el('h2', { text: 'Currencies in spacetime' }),
    el('p', { class: 'lede', text: `${A.currencies.length} currencies, ${A.months[0]} to ${A.months.at(-1)}: each held from month to month at its policy rate and exchanged with ${base(A).name} at the month's rate. A month's lens is a loop: ${base(A).holding} held and exchanged, against exchanged and held.` }),
    el('dl', {},
      el('dt', { text: 'presences (currency, month)' }), el('dd', { text: A.cells.presences }),
      el('dt', { text: 'holdings and exchanges' }), el('dd', { text: A.cells.edges }),
      el('dt', { text: 'lenses' }), el('dd', { text: A.cells.lenses })),
    el('div', { class: 'heading' }, el('span', { text: 'how far from a price system' })),
    el('p', {}, `A single system of intertemporal prices, one log price for each currency in each month, would flatten every lens. The lenses' curvature is ${(100 * K.rms).toFixed(2)}% a month (root mean square). With nothing earned by holding, the exchange rates alone would curve them by ${(100 * K.rms_exchange).toFixed(2)}%; the policy rates remove ${percent(K.removed)} of that (1 − Σκ²/Σκ₀²). ${K.removed < 0.05 ? 'Month by month the rates offset almost none of the movement of exchange rates, which parity says they should offset in expectation; what they leave on average is the carry in the table below.' : 'The rates offset part of the movement of exchange rates, as parity says they should in expectation; what they leave on average is the carry in the table below.'}`),
    el('p', { class: 'note', text: 'Curvature does not change when value is re-described, by a change of a currency\'s unit for instance; the share of the transports that a price system would carry does, and is therefore not shown.' }),
  );
  const rows = [...A.stats].sort((a, b) => a.mean - b.mean);
  put(box, el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: `${base(A).name} against each currency` })),
    el('table', {}, el('thead', {}, el('tr', {}, el('th', { text: 'currency' }), el('th', { class: 'num', text: 'a year' }), el('th', { class: 'num', text: 'volatility' }), el('th', { class: 'num', text: 'ratio' }))),
      el('tbody', {}, rows.map((s) => el('tr', { class: 'pick', style: state.selected === s.currency ? 'background: var(--gold-wash)' : '', onclick: () => select(s.currency) },
        el('td', { text: `${A.currencies[s.currency]}` }), el('td', { class: 'num', text: `${signed(100 * s.mean, 1)}%` }),
        el('td', { class: 'num', text: `${(100 * s.volatility).toFixed(1)}%` }), el('td', { class: 'num', text: fixed(s.volatility > 0 ? s.mean / s.volatility : 0, 2) })))))),
    el('p', { class: 'note', text: `A negative mean: holding the currency, at its higher rate, beat holding ${base(A).holding} even after its exchange rate moved: the carry that interest parity says should vanish.` }));
  if (state.selected !== null) {
    const c = state.selected;
    const policy = A.policy[c], ofBase = A.policy[0];
    const mean = (v) => v.reduce((s, x) => s + x, 0) / v.length;
    put(box, el('div', { class: 'section' }, el('h3', { text: `${A.names[c]} (${A.currencies[c]})` }),
      el('dl', {},
        el('dt', { text: 'policy rate, mean' }), el('dd', { text: `${mean(policy).toFixed(2)}%` }),
        el('dt', { text: `${base(A).its}, mean` }), el('dd', { text: `${mean(ofBase).toFixed(2)}%` }),
        el('dt', { text: `rate ${base(A).per}, ${A.months[0]}` }), el('dd', { text: fixed(A.rates[c][0], 4) }),
        el('dt', { text: `rate ${base(A).per}, ${A.months.at(-1)}` }), el('dd', { text: fixed(A.rates[c].at(-1), 4) }),
        el('dt', { text: 'its exchange-rate curvature the rates remove' }), el('dd', { text: percent(A.stats.find((s) => s.currency === c).removed) }))));
  }
  put(box, el('p', { class: 'note', text: `The nearest price system lies ${K.distance.toFixed(3)} from the record (the norm of the transports' co-exact and harmonic parts, found ${A.method === 'sparse' ? `by conjugate gradients in ${A.iterations} iterations` : 'exactly'}). ${shared.currencies === 'library:currencies' ? 'Rates: ECB; policy rates: BIS, the euro\u2019s the ECB\u2019s deposit facility rate.' : `From ${dataset(shared.currencies)?.title || 'your record'}.`}` }));
}

async function load() {
  const first = document.getElementById('st-first').value, last = document.getElementById('st-last').value;
  const status = document.getElementById('status');
  // The currency shown is kept by its code: the currencies with complete records change with the span.
  const keep = state.wanted ?? (state.answer && state.selected !== null ? state.answer.currencies[state.selected] : null);
  state.wanted = undefined;
  status.textContent = 'building the spacetime of currencies…';
  try {
    state.answer = await api.spacetime(first, last, shared.currencies);
    status.textContent = '';
  } catch (error) {
    status.textContent = error.message;
    return;
  }
  const i = keep ? state.answer.currencies.indexOf(keep) : -1;
  state.selected = i > 0 ? i : null;
  document.getElementById('st-note').textContent = `${state.answer.currencies.length} currencies with complete records over the span`;
  sheet();
  lines();
  panel();
}

// The months of the chosen record, in the two selects; a span kept where the record allows it.
function months(preferFirst) {
  const d = dataset(shared.currencies);
  const [m0, m1] = d?.months || [];
  const all = [];
  if (m0 && m1) {
    let [y, m] = [Number(m0.slice(0, 4)), Number(m0.slice(5, 7))];
    while (`${y}-${String(m).padStart(2, '0')}` <= m1) {
      all.push(`${y}-${String(m).padStart(2, '0')}`);
      [y, m] = m === 12 ? [y + 1, 1] : [y, m + 1];
    }
  }
  state.months = all;
  const first = document.getElementById('st-first'), last = document.getElementById('st-last');
  for (const s of [first, last]) fill(s, ...all.map((x) => el('option', { value: x, text: x })));
  first.value = all.includes(preferFirst) ? preferFirst : all[0];
  last.value = all.at(-1);
}

function datasetSelect() {
  const select = document.getElementById('st-dataset');
  fill(select, ...datasetsOf('currencies').map((d) => el('option', { value: d.id, text: d.library ? d.title : `${d.title} (yours)`, selected: d.id === shared.currencies,
    title: describe(d) })), el('option', { value: '+', text: '+ your data…' }));
}

async function changeRecord() {
  state.answer = null;
  state.selected = null;
  months(shared.currencies === 'library:currencies' ? '2019-01' : null);
  datasetSelect();
  await load();
}

// The address may ask for a record, a span and a currency: #spacetime&from=2015-01&to=2026-06&currency=MXN.
export async function showSpacetime(data, params = {}) {
  if (!state.wired) {
    state.wired = true;
    const first = document.getElementById('st-first'), last = document.getElementById('st-last');
    first.addEventListener('change', load);
    last.addEventListener('change', load);
    document.getElementById('st-dataset').addEventListener('change', (e) => {
      if (e.target.value === '+') { datasetSelect(); openDataDialog('currencies'); return; }
      shared.currencies = e.target.value;
      changeRecord();
    });
    on('dataset-chosen', (kind) => { if (kind === 'currencies') changeRecord(); });
    on('datasets', datasetSelect);
    onTheme(() => { if (state.answer) { sheet(); lines(); panel(); } });
    window.addEventListener('resize', () => { if (state.answer && !document.querySelector('.view[data-view="spacetime"]').hidden) { sheet(); lines(); } });
    if (params.dataset && dataset(params.dataset)) shared.currencies = params.dataset;
    months(shared.currencies === 'library:currencies' ? '2019-01' : null);
    datasetSelect();
  } else if (params.dataset && params.dataset !== shared.currencies && dataset(params.dataset)) {
    shared.currencies = params.dataset;
    state.answer = null;
    months(null);
    datasetSelect();
  }
  const first = document.getElementById('st-first'), last = document.getElementById('st-last');
  let changed = false;
  if (state.months.includes(params.from) && params.from !== first.value) { first.value = params.from; changed = true; }
  if (state.months.includes(params.to) && params.to !== last.value) { last.value = params.to; changed = true; }
  if ('currency' in params) state.wanted = params.currency.toUpperCase();
  if (!state.answer || changed) await load();
  else {
    if (state.wanted !== undefined) {
      const i = state.answer.currencies.indexOf(state.wanted);
      state.selected = i > 0 ? i : null;
      state.wanted = undefined;
    }
    sheet(); lines(); panel();
  }
}
