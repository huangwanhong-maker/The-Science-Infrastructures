// Line charts and bar charts in SVG, on hairline axes, coloured by the page's tokens.

import { css, fmt, hover, svg } from '../dom.js';

function ticks(lo, hi, n = 5) {
  if (hi === lo) return [lo];
  const step0 = (hi - lo) / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= step0) || step0;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-12; t += step) out.push(+t.toFixed(12));
  return out;
}

// series: [{name, x: [], y: [], color, width?, dash?, opacity?}]; options: {width, height, xlabel, ylabel, ylim, xlim}
export function lines(series, o = {}) {
  const W = o.width || 560, H = o.height || 260, L = 52, R = 12, T = 26, B = 36;
  const xs = series.flatMap((s) => s.x), ys = series.flatMap((s) => s.y).filter(Number.isFinite);
  const [x0, x1] = o.xlim || [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = o.ylim || [Math.min(...ys), Math.max(...ys)];
  if (y0 === y1) { y0 -= 1; y1 += 1; }
  const X = (x) => L + ((x - x0) / (x1 - x0 || 1)) * (W - L - R);
  const Y = (y) => H - B - ((y - y0) / (y1 - y0)) * (H - T - B);
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  const axis = svg('g', { class: 'axis' });
  for (const t of ticks(y0, y1)) axis.append(svg('line', { x1: L, x2: W - R, y1: Y(t), y2: Y(t), 'stroke-width': 0.6 }), svg('text', { x: L - 6, y: Y(t) + 4, 'text-anchor': 'end' }, fmt(t, 3)));
  for (const t of ticks(x0, x1)) axis.append(svg('text', { x: X(t), y: H - B + 16, 'text-anchor': 'middle' }, fmt(t, 3)));
  if (o.xlabel) axis.append(svg('text', { x: (L + W - R) / 2, y: H - 4, 'text-anchor': 'middle' }, o.xlabel));
  if (o.ylabel) axis.append(svg('text', { x: L, y: 12, 'text-anchor': 'start' }, o.ylabel));
  root.append(axis);
  for (const s of series) {
    let d = '';
    s.x.forEach((x, i) => { if (Number.isFinite(s.y[i])) d += `${d ? 'L' : 'M'}${X(x).toFixed(1)},${Y(s.y[i]).toFixed(1)}`; });
    const p = svg('path', { d, fill: 'none', stroke: s.color || css('--series-1'), 'stroke-width': s.width || 1.8, 'stroke-dasharray': s.dash || null, opacity: s.opacity ?? 1, 'stroke-linejoin': 'round' });
    if (s.name) hover(p, s.name);
    root.append(p);
  }
  for (const m of o.marks || []) root.append(svg('circle', { cx: X(m.x), cy: Y(m.y), r: m.r || 4, fill: m.fill || 'none', stroke: m.color || css('--ink'), 'stroke-width': 1.5 }));
  return root;
}

// bars: [{x, y, color?, label?}] on a numeric x axis.
export function bars(data, o = {}) {
  const W = o.width || 560, H = o.height || 240, L = 52, R = 12, T = 26, B = 36;
  const xs = data.map((d) => d.x), ys = data.map((d) => d.y);
  const x0 = Math.min(...xs) - 0.5, x1 = Math.max(...xs) + 0.5, y1 = Math.max(...ys, 1e-12) * 1.08;
  const X = (x) => L + ((x - x0) / (x1 - x0)) * (W - L - R);
  const Y = (y) => H - B - (y / y1) * (H - T - B);
  const bw = Math.max(1, (W - L - R) / (x1 - x0) - 2);
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  const axis = svg('g', { class: 'axis' });
  for (const t of ticks(0, y1)) axis.append(svg('line', { x1: L, x2: W - R, y1: Y(t), y2: Y(t), 'stroke-width': 0.6 }), svg('text', { x: L - 6, y: Y(t) + 4, 'text-anchor': 'end' }, fmt(t, 3)));
  for (const t of ticks(x0 + 0.5, x1 - 0.5)) axis.append(svg('text', { x: X(t), y: H - B + 16, 'text-anchor': 'middle' }, fmt(t, 3)));
  if (o.xlabel) axis.append(svg('text', { x: (L + W - R) / 2, y: H - 4, 'text-anchor': 'middle' }, o.xlabel));
  if (o.ylabel) axis.append(svg('text', { x: L, y: 12 }, o.ylabel));
  root.append(axis);
  for (const d of data) {
    const r = svg('rect', { x: X(d.x) - bw / 2, y: Y(d.y), width: bw, height: Y(0) - Y(d.y), fill: d.color || css('--series-1'), rx: 1.5 });
    hover(r, d.label || `${fmt(d.x)}: ${fmt(d.y)}`);
    root.append(r);
  }
  return root;
}

export function legend(items) {
  const el = document.createElement('div');
  el.className = 'legend';
  for (const it of items) {
    const s = document.createElement('span');
    const i = document.createElement('i');
    i.style.background = it.color;
    s.append(i, it.name);
    el.append(s);
  }
  return el;
}
