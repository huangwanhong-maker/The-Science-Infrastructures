// Small charts for the toolbox's answers, in SVG and in the atlas's colours: lines (with bands),
// bars (grouped or stacked, with a band), a scatter (with a line), and a matrix of values. Each
// answers the pointer with a card, has its legend, and keeps its marks thin.

import { el, svg, put } from './dom.js';
import { token, diverging } from './palette.js';

// The series' colours in a fixed order, validated as neighbours for colour vision deficiencies:
// indigo, rose, gold, sky, green, plum; more series than six share a legend and a table instead.
export const LINE_COLOURS = () => [token('driven'), token('rose'), token('global'), token('sends'), '#5E9A4C', '#A05C9C'];

const fmt = (x) => {
  if (!Number.isFinite(x)) return '—';
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return x.toExponential(2);
  return x.toLocaleString('en-GB', { maximumFractionDigits: a >= 100 ? 1 : a >= 1 ? 3 : 4 });
};
export { fmt as formatNumber };

// Ticks at 1, 2 or 5 times a power of ten.
function ticks(lo, hi, n = 5) {
  if (!(hi > lo)) { const d = Math.abs(lo) || 1; lo -= d / 2; hi += d / 2; }
  const raw = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 5, 10].map((m) => m * p).find((s) => s >= raw) || 10 * p;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) out.push(+t.toPrecision(12));
  return out;
}

function frame(title, width, height) {
  const box = el('figure', { class: 'chart' }, el('figcaption', { text: title }));
  const host = svg('svg', { viewBox: `0 0 ${width} ${height}`, role: 'img', 'aria-label': title });
  put(box, host);
  return { box, host };
}

function hoverCard(box) {
  const card = el('div', { class: 'chart-tip' });
  put(box, card);
  return {
    show(x, y, rows) {
      card.replaceChildren(...rows.map(([k, v, c]) => el('div', { class: 'row' }, el('span', {},
        c ? el('i', { style: { background: c } }) : null, k), el('span', { text: v }))));
      card.style.left = `${x + 14}px`;
      card.style.top = `${y + 8}px`;
      card.classList.add('on');
    },
    hide() { card.classList.remove('on'); },
  };
}

function axes(host, X, Y, xt, yt, W, H, m, xlabel, ylabel, xfmt = fmt) {
  for (const t of yt) {
    put(host, svg('line', { x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), class: t === 0 ? 'zero' : 'grid' }),
      svg('text', { x: m.l - 6, y: Y(t) + 4, 'text-anchor': 'end', class: 'tick', text: fmt(t) }));
  }
  for (const t of xt) put(host, svg('text', { x: X(t), y: H - m.b + 16, 'text-anchor': 'middle', class: 'tick', text: xfmt(t) }));
  if (xlabel) put(host, svg('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', class: 'axis-label', text: xlabel }));
  if (ylabel) put(host, svg('text', { x: 12, y: m.t - 6, class: 'axis-label', text: ylabel }));
}

// `compact` draws it for a narrow panel: fewer pixels, so that its text stays readable there.
export function lineChart({ title, x, series, xlabel = '', ylabel = '', compact = false }) {
  const W = compact ? 360 : 640, H = compact ? 200 : 300, m = { l: compact ? 50 : 58, r: compact ? 64 : 110, t: 22, b: compact ? 34 : 40 };
  const { box, host } = frame(title, W, H);
  const colours = LINE_COLOURS();
  const all = series.flatMap((s) => [...s.y, ...(s.lower || []), ...(s.upper || [])]).filter(Number.isFinite);
  const xs = x.filter(Number.isFinite);
  if (!all.length || !xs.length) return el('p', { class: 'note', text: `${title}: nothing to draw` });
  const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
  const yt = ticks(Math.min(...all), Math.max(...all));
  const [y0, y1] = [Math.min(yt[0], ...all), Math.max(yt.at(-1), ...all)];
  const X = (v) => m.l + ((W - m.l - m.r) * (v - x0)) / (x1 - x0 || 1);
  const Y = (v) => m.t + ((H - m.t - m.b) * (y1 - v)) / (y1 - y0 || 1);
  const integer = xs.every((v) => Number.isInteger(v));
  // Years, whole or decimal, are written without a thousands separator.
  const years = x0 >= 1000 && x1 <= 3000;
  axes(host, X, Y, ticks(x0, x1, 6).filter((t) => t >= x0 && t <= x1), yt, W, H, m, xlabel, ylabel,
    integer || years ? (v) => String(+v.toFixed(2)) : fmt);
  series.forEach((s, k) => {
    const c = colours[k % colours.length];
    if (s.lower && s.upper) {
      let d = '';
      x.forEach((v, i) => { if (Number.isFinite(s.upper[i])) d += `${d ? 'L' : 'M'}${X(v).toFixed(1)},${Y(s.upper[i]).toFixed(1)}`; });
      for (let i = x.length - 1; i >= 0; i--) if (Number.isFinite(s.lower[i])) d += `L${X(x[i]).toFixed(1)},${Y(s.lower[i]).toFixed(1)}`;
      put(host, svg('path', { d: `${d}Z`, fill: c, 'fill-opacity': 0.14, stroke: 'none' }));
    }
    let d = '', pen = false;
    x.forEach((v, i) => {
      if (!Number.isFinite(s.y[i])) { pen = false; return; }
      d += `${pen ? 'L' : 'M'}${X(v).toFixed(1)},${Y(s.y[i]).toFixed(1)}`;
      pen = true;
    });
    put(host, svg('path', { d, fill: 'none', stroke: c, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-dasharray': s.dashed ? '6 4' : null }));
    // A direct label at the series' last point, for two to four series (the title names a single one);
    // a series that stops before the right edge (the history before a forecast) is labelled above its end.
    if (series.length > 1 && series.length <= 4) {
      const last = [...s.y.keys()].reverse().find((i) => Number.isFinite(s.y[i]));
      const edge = x[last] >= x1;
      if (last !== undefined) {
        put(host, svg('text', { x: X(x[last]) + (edge ? 6 : -4), y: Y(s.y[last]) + (edge ? 4 : -10), 'text-anchor': edge ? 'start' : 'end',
          class: 'series-label', fill: 'currentColor', text: s.name }));
      }
    }
  });
  if (series.length > 4) {
    put(box, el('div', { class: 'chart-legend' }, series.map((s, k) => el('span', {}, el('i', { style: { background: colours[k % colours.length] } }), s.name))));
  }
  // The pointer: the nearest x, every series' value there.
  const tip = hoverCard(box);
  const cursor = svg('line', { y1: m.t, y2: H - m.b, class: 'cursor', visibility: 'hidden' });
  put(host, cursor);
  host.addEventListener('pointermove', (e) => {
    const r = host.getBoundingClientRect();
    const vx = x0 + ((e.clientX - r.left) * (W / r.width) - m.l) / (W - m.l - m.r) * (x1 - x0);
    let i = 0;
    x.forEach((v, j) => { if (Math.abs(v - vx) < Math.abs(x[i] - vx)) i = j; });
    cursor.setAttribute('x1', X(x[i]));
    cursor.setAttribute('x2', X(x[i]));
    cursor.setAttribute('visibility', 'visible');
    tip.show(e.clientX - box.getBoundingClientRect().left, e.clientY - box.getBoundingClientRect().top,
      [[xlabel || 'x', integer ? String(x[i]) : fmt(x[i])], ...series.map((s, k) => [s.name, fmt(s.y[i]), colours[k % colours.length]])]);
  });
  host.addEventListener('pointerleave', () => { tip.hide(); cursor.setAttribute('visibility', 'hidden'); });
  return box;
}

export function barChart({ title, labels, series, stacked = false, band = null }) {
  const W = 640, H = 280, m = { l: 58, r: 16, t: 22, b: 52 };
  const { box, host } = frame(title, W, H);
  const colours = LINE_COLOURS();
  const n = labels.length;
  let lo = 0, hi = 0;
  for (let i = 0; i < n; i++) {
    if (stacked) {
      let pos = 0, neg = 0;
      for (const s of series) { const v = s.values[i] || 0; if (v >= 0) pos += v; else neg += v; }
      hi = Math.max(hi, pos);
      lo = Math.min(lo, neg);
    } else for (const s of series) { const v = s.values[i]; if (Number.isFinite(v)) { hi = Math.max(hi, v); lo = Math.min(lo, v); } }
  }
  if (band) { lo = Math.min(lo, band[0]); hi = Math.max(hi, band[1]); }
  const yt = ticks(lo, hi);
  const [y0, y1] = [Math.min(yt[0], lo), Math.max(yt.at(-1), hi)];
  const Y = (v) => m.t + ((H - m.t - m.b) * (y1 - v)) / (y1 - y0 || 1);
  const slot = (W - m.l - m.r) / Math.max(1, n);
  axes(host, () => 0, Y, [], yt, W, H, m, '', '');
  const every = Math.ceil(n / 24);
  labels.forEach((lab, i) => {
    if (i % every === 0) put(host, svg('text', { x: m.l + slot * (i + 0.5), y: H - m.b + 16, 'text-anchor': 'middle', class: 'tick', text: lab.length > 8 ? lab.slice(0, 8) : lab }));
  });
  if (band) for (const b of band) put(host, svg('line', { x1: m.l, x2: W - m.r, y1: Y(b), y2: Y(b), class: 'band' }));
  const width = stacked ? Math.max(2, slot - 2) : Math.max(1.5, (slot - 2) / series.length - 2);
  for (let i = 0; i < n; i++) {
    let pos = 0, neg = 0;
    series.forEach((s, k) => {
      const v = s.values[i];
      if (!Number.isFinite(v)) return;
      let top, bottom;
      if (stacked) { if (v >= 0) { bottom = pos; pos += v; top = pos; } else { top = neg; neg += v; bottom = neg; } }
      else { top = Math.max(0, v); bottom = Math.min(0, v); }
      const xLeft = m.l + slot * i + 1 + (stacked ? 0 : k * (width + 2));
      put(host, svg('rect', { x: xLeft.toFixed(1), y: Y(top).toFixed(1), width: width.toFixed(1), height: Math.max(0.5, Y(bottom) - Y(top)).toFixed(1),
        rx: 2, fill: colours[k % colours.length], 'data-i': i }));
    });
  }
  if (series.length > 1) put(box, el('div', { class: 'chart-legend' }, series.map((s, k) => el('span', {}, el('i', { style: { background: colours[k % colours.length] } }), s.name))));
  const tip = hoverCard(box);
  host.addEventListener('pointermove', (e) => {
    const r = host.getBoundingClientRect();
    const i = Math.floor(((e.clientX - r.left) * (W / r.width) - m.l) / slot);
    if (i < 0 || i >= n) return tip.hide();
    tip.show(e.clientX - box.getBoundingClientRect().left, e.clientY - box.getBoundingClientRect().top,
      [[labels[i], ''], ...series.map((s, k) => [s.name, fmt(s.values[i]), colours[k % colours.length]])]);
  });
  host.addEventListener('pointerleave', () => tip.hide());
  return box;
}

export function scatterChart({ title, x, y, xlabel = '', ylabel = '', labels = null, line = null }) {
  const W = 640, H = 320, m = { l: 58, r: 16, t: 22, b: 40 };
  const { box, host } = frame(title, W, H);
  const pts = x.map((v, i) => [v, y[i], labels ? labels[i] : '']).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
  if (!pts.length) return el('p', { class: 'note', text: `${title}: nothing to draw` });
  const xt = ticks(Math.min(...pts.map((p) => p[0])), Math.max(...pts.map((p) => p[0])));
  const yt = ticks(Math.min(...pts.map((p) => p[1])), Math.max(...pts.map((p) => p[1])));
  const X = (v) => m.l + ((W - m.l - m.r) * (v - xt[0])) / (xt.at(-1) - xt[0] || 1);
  const Y = (v) => m.t + ((H - m.t - m.b) * (yt.at(-1) - v)) / (yt.at(-1) - yt[0] || 1);
  axes(host, X, Y, xt, yt, W, H, m, xlabel, ylabel);
  const dots = svg('g', {});
  pts.forEach(([a, b], i) => put(dots, svg('circle', { cx: X(a).toFixed(1), cy: Y(b).toFixed(1), r: 3.2, fill: token('sends'), 'fill-opacity': 0.55, stroke: 'var(--card)', 'stroke-width': 0.8, 'data-i': i })));
  put(host, dots);
  if (line) {
    let d = '';
    line.x.forEach((v, i) => { d += `${i ? 'L' : 'M'}${X(v).toFixed(1)},${Y(line.y[i]).toFixed(1)}`; });
    put(host, svg('path', { d, fill: 'none', stroke: token('rose'), 'stroke-width': 2 }));
  }
  const tip = hoverCard(box);
  dots.addEventListener('pointermove', (e) => {
    const i = Number(e.target.dataset?.i);
    if (!Number.isInteger(i)) return;
    const [a, b, lab] = pts[i];
    tip.show(e.clientX - box.getBoundingClientRect().left, e.clientY - box.getBoundingClientRect().top,
      [[lab || 'point', ''], [xlabel || 'x', fmt(a)], [ylabel || 'y', fmt(b)]]);
  });
  dots.addEventListener('pointerleave', () => tip.hide());
  return box;
}

export function matrixChart({ title, rows, columns, values }) {
  const flat = values.flat().filter(Number.isFinite);
  const big = Math.max(1e-12, ...flat.map(Math.abs));
  return el('figure', { class: 'chart' }, el('figcaption', { text: title }),
    el('div', { class: 'matrix-wrap' }, el('table', { class: 'matrix-chart' },
      el('thead', {}, el('tr', {}, el('th', {}), columns.map((c) => el('th', { text: c })))),
      el('tbody', {}, rows.map((r, i) => el('tr', {}, el('th', { text: r }),
        columns.map((_, j) => {
          const v = values[i][j];
          return el('td', { style: { background: Number.isFinite(v) ? diverging(v / big) : 'transparent' }, text: fmt(v) });
        })))))));
}
