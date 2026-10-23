// Charts of curves, in SVG: several series against a common axis, linear or logarithmic,
// with reference lines (a critical point, an expected dimension, the classical value) and
// marked points (a bounce); and histograms of small integers (the spins of a spin foam).

import { svg, el } from './dom.js';
import { number } from './format.js';

const COLOURS = { rose: '#C08081', gold: '#f2c27c', lilac: '#b9a2e0', blush: '#e8c3dc', sage: '#9fc2b8', faint: 'rgba(240,226,238,0.35)' };

// series: [{x, y, colour, dash, label, width}]; options: {xlabel, ylabel, logx, logy, hlines:
// [{y, label, colour}], vlines: [{x, label, colour}], points: [{x, y, label, colour}], height,
// ylim: [lo, hi], xlim}.
export function xyChart(series, options = {}) {
  const W = 360, H = options.height ?? 190, L = 46, R = 10, T = 12, B = 32;
  const fx = options.logx ? Math.log10 : (v) => v, fy = options.logy ? Math.log10 : (v) => v;
  const xs = [], ys = [];
  for (const s of series) s.x.forEach((x, i) => { const y = s.y[i]; if (Number.isFinite(x) && Number.isFinite(y) && (!options.logx || x > 0) && (!options.logy || y > 0)) { xs.push(fx(x)); ys.push(fy(y)); } });
  for (const h of options.hlines || []) if (Number.isFinite(h.y)) ys.push(fy(h.y));
  for (const p of options.points || []) if (Number.isFinite(p.y)) { xs.push(fx(p.x)); ys.push(fy(p.y)); }
  let [x0, x1] = options.xlim ? options.xlim.map(fx) : [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = options.ylim ? options.ylim.map(fy) : [Math.min(...ys), Math.max(...ys)];
  if (!Number.isFinite(x0)) { x0 = 0; x1 = 1; }
  if (!Number.isFinite(y0)) { y0 = 0; y1 = 1; }
  if (x1 - x0 < 1e-12) { x0 -= 0.5; x1 += 0.5; }
  if (y1 - y0 < 1e-12) { y0 -= 0.5; y1 += 0.5; }
  if (!options.ylim) { const pad = (y1 - y0) * 0.07; y0 -= pad; y1 += pad; }
  const X = (v) => L + ((fx(v) - x0) / (x1 - x0)) * (W - L - R), Y = (v) => T + (1 - (fy(v) - y0) / (y1 - y0)) * (H - T - B);
  const out = svg('svg', { class: 'chart xy', viewBox: `0 0 ${W} ${H}`, role: 'img' });
  out.append(svg('rect', { x: L, y: T, width: W - L - R, height: H - T - B, fill: 'rgba(240,226,238,0.025)', stroke: 'rgba(240,226,238,0.12)' }));
  // Ticks: ends and middle of each axis, in the units of the data.
  const back = (v, log) => (log ? 10 ** v : v);
  for (const t of [0, 0.5, 1]) {
    const vx = back(x0 + t * (x1 - x0), options.logx), vy = back(y0 + t * (y1 - y0), options.logy);
    out.append(svg('text', { x: X(vx), y: H - B + 13, 'text-anchor': t === 0 ? 'start' : t === 1 ? 'end' : 'middle' }, number(vx, 3)));
    out.append(svg('text', { x: L - 5, y: Y(vy) + 4, 'text-anchor': 'end' }, number(vy, 3)));
  }
  if (options.xlabel) out.append(svg('text', { x: (L + W - R) / 2, y: H - 4, 'text-anchor': 'middle', class: 'axis' }, options.xlabel));
  if (options.ylabel) out.append(svg('text', { x: 4, y: T + 2, 'text-anchor': 'start', class: 'axis', transform: `rotate(90 4 ${T + 2})` }, options.ylabel));
  const clip = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  for (const h of options.hlines || []) {
    if (!Number.isFinite(h.y) || fy(h.y) < y0 || fy(h.y) > y1) continue;
    out.append(svg('line', { x1: L, x2: W - R, y1: Y(h.y), y2: Y(h.y), stroke: COLOURS[h.colour] || h.colour || COLOURS.gold, 'stroke-dasharray': '4 3', 'stroke-width': 1.1 }));
    if (h.label) out.append(svg('text', { x: W - R - 3, y: Y(h.y) - 4, 'text-anchor': 'end', class: 'mark' }, h.label));
  }
  for (const v of options.vlines || []) {
    if (!Number.isFinite(v.x) || fx(v.x) < x0 || fx(v.x) > x1) continue;
    out.append(svg('line', { x1: X(v.x), x2: X(v.x), y1: T, y2: H - B, stroke: COLOURS[v.colour] || v.colour || COLOURS.gold, 'stroke-dasharray': '4 3', 'stroke-width': 1.1 }));
    if (v.label) { const right = X(v.x) > W - R - 70; out.append(svg('text', { x: X(v.x) + (right ? -4 : 4), y: T + 12, 'text-anchor': right ? 'end' : 'start', class: 'mark' }, v.label)); }
  }
  for (const s of series) {
    const pts = [];
    s.x.forEach((x, i) => {
      const y = s.y[i];
      if (!Number.isFinite(x) || !Number.isFinite(y) || (options.logx && x <= 0) || (options.logy && y <= 0)) return;
      pts.push(`${X(x).toFixed(1)},${clip(Y(y), T, H - B).toFixed(1)}`);
    });
    if (pts.length > 1) out.append(svg('polyline', { points: pts.join(' '), fill: 'none', stroke: COLOURS[s.colour] || s.colour || COLOURS.rose, 'stroke-width': s.width ?? 1.5, 'stroke-dasharray': s.dash || null, 'stroke-linejoin': 'round' }));
  }
  for (const p of options.points || []) {
    if (!Number.isFinite(p.x) || !Number.isFinite(p.y)) continue;
    out.append(svg('circle', { cx: X(p.x), cy: Y(p.y), r: 3.4, fill: COLOURS[p.colour] || p.colour || COLOURS.gold }));
    if (p.label) out.append(svg('text', { x: X(p.x) + 6, y: Y(p.y) - 6, class: 'mark' }, p.label));
  }
  return out;
}

// A legend line under a chart: [{label, colour, dash}].
export function chartLegend(items) {
  return el('div', { class: 'chart-legend' }, ...items.map((i) => el('span', {},
    svg('svg', { width: 22, height: 8 }, svg('line', { x1: 1, x2: 21, y1: 4, y2: 4, stroke: COLOURS[i.colour] || i.colour, 'stroke-width': 2, 'stroke-dasharray': i.dash || null })), ` ${i.label}`)));
}

// A histogram of values (counts by value), bars in Old Rose with their value beneath.
export function histogram(values, { label = (v) => String(v), height = 110 } = {}) {
  const counts = new Map();
  for (const v of values) counts.set(v, (counts.get(v) || 0) + 1);
  const keys = [...counts.keys()].sort((a, b) => a - b);
  const W = 360, H = height, B = 22, T = 14, gap = 4;
  const w = Math.min(46, (W - gap * (keys.length + 1)) / Math.max(1, keys.length)), max = Math.max(1, ...counts.values());
  const out = svg('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}` });
  keys.forEach((k, i) => {
    const h = ((H - B - T) * counts.get(k)) / max, x = gap + i * (w + gap);
    out.append(svg('rect', { x, y: H - B - h, width: w, height: h, rx: 2, fill: '#C08081', 'fill-opacity': 0.75 }));
    out.append(svg('text', { x: x + w / 2, y: H - B - h - 3, 'text-anchor': 'middle', class: 'mark' }, String(counts.get(k))));
    out.append(svg('text', { x: x + w / 2, y: H - 7, 'text-anchor': 'middle' }, label(k)));
  });
  return out;
}
