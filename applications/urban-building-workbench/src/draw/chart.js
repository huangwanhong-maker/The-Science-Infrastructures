// Charts in SVG on hairline axes, coloured by the page's tokens, and rasters on canvas coloured by the
// visualization layer's validated maps: lines (linear or logarithmic in y, with bands and marked
// frequencies), bars, scatter plots with labels, and a raster of a field over two axes.

import { css, fmt, h, hover, svg } from '../dom.js';
import { colourMap } from '/visualization/styles/index.js';

function ticks(lo, hi, n = 5) {
  if (hi === lo) return [lo];
  const step0 = (hi - lo) / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= step0) || step0;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-12; t += step) out.push(+t.toFixed(12));
  return out;
}

function decades(lo, hi) {
  const out = [];
  for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) out.push(10 ** e);
  return out.filter((t) => t >= lo * 0.999 && t <= hi * 1.001);
}

// series: [{name, x, y, color, width?, dash?, opacity?, fill? (y2 for a band)}]
// o: {width, height, xlabel, ylabel, ylim, xlim, logy, marks: [{x, label, color?, dash?}], points: [{x, y, label, color}]}
export function lines(series, o = {}) {
  const W = o.width || 620, H = o.height || 260, L = 56, R = 14, T = 26, B = 36;
  const xs = series.flatMap((s) => s.x).filter(Number.isFinite);
  const ys = series.flatMap((s) => [...s.y, ...(s.y2 || [])]).filter((v) => Number.isFinite(v) && (!o.logy || v > 0));
  const [x0, x1] = o.xlim || [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = o.ylim || [Math.min(...ys), Math.max(...ys)];
  if (y0 === y1) { y0 = o.logy ? y0 / 10 : y0 - 1; y1 = o.logy ? y1 * 10 : y1 + 1; }
  const ly = (v) => (o.logy ? Math.log10(Math.max(v, y0)) : v);
  const X = (x) => L + ((x - x0) / (x1 - x0 || 1)) * (W - L - R);
  const Y = (y) => H - B - ((ly(y) - ly(y0)) / (ly(y1) - ly(y0) || 1)) * (H - T - B);
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  const axis = svg('g', { class: 'axis' });
  for (const t of o.logy ? decades(y0, y1) : ticks(y0, y1))
    axis.append(svg('line', { x1: L, x2: W - R, y1: Y(t), y2: Y(t), 'stroke-width': 0.6 }), svg('text', { x: L - 6, y: Y(t) + 4, 'text-anchor': 'end' }, o.logy ? t.toExponential(0) : fmt(t, 3)));
  // Years are written plainly (2010, not 2,010).
  const xfmt = o.xfmt || (/year/.test(o.xlabel || '') ? (t) => String(Math.round(t)) : (t) => fmt(t, 3));
  for (const t of ticks(x0, x1, o.xticks || 6)) axis.append(svg('text', { x: X(t), y: H - B + 16, 'text-anchor': 'middle' }, xfmt(t)));
  if (o.xlabel) axis.append(svg('text', { x: (L + W - R) / 2, y: H - 4, 'text-anchor': 'middle' }, o.xlabel));
  if (o.ylabel) axis.append(svg('text', { x: L, y: 12, 'text-anchor': 'start' }, o.ylabel));
  root.append(axis);
  for (const m of o.marks || []) {
    if (m.x < x0 || m.x > x1) continue;
    root.append(svg('line', { x1: X(m.x), x2: X(m.x), y1: T, y2: H - B, stroke: m.color || css('--hair-2'), 'stroke-dasharray': m.dash || '3 3', 'stroke-width': 1 }));
    if (m.label) root.append(svg('text', { x: X(m.x) + 3, y: T + 10 + (m.row || 0) * 11, class: 'mark-label' }, m.label));
  }
  for (const s of series) {
    if (s.y2) {
      let d = '';
      s.x.forEach((x, i) => { if (Number.isFinite(s.y[i])) d += `${d ? 'L' : 'M'}${X(x).toFixed(1)},${Y(s.y[i]).toFixed(1)}`; });
      for (let i = s.x.length - 1; i >= 0; i--) if (Number.isFinite(s.y2[i])) d += `L${X(s.x[i]).toFixed(1)},${Y(s.y2[i]).toFixed(1)}`;
      root.append(svg('path', { d: d + 'Z', fill: s.color, opacity: s.opacity ?? 0.18, stroke: 'none' }));
      continue;
    }
    let d = '';
    s.x.forEach((x, i) => { if (Number.isFinite(s.y[i]) && (!o.logy || s.y[i] > 0)) d += `${d ? 'L' : 'M'}${X(x).toFixed(1)},${Y(s.y[i]).toFixed(1)}`; else if (s.gaps) d += ''; });
    const p = svg('path', { d, fill: 'none', stroke: s.color || css('--series-1'), 'stroke-width': s.width || 1.6, 'stroke-dasharray': s.dash || null, opacity: s.opacity ?? 1, 'stroke-linejoin': 'round' });
    if (s.name) hover(p, s.name);
    root.append(p);
  }
  for (const pt of o.points || []) {
    const c = svg('circle', { cx: X(pt.x), cy: Y(pt.y), r: pt.r || 3.5, fill: pt.color || css('--accent'), stroke: css('--surface'), 'stroke-width': 1 });
    if (pt.label) hover(c, pt.label);
    root.append(c);
  }
  return root;
}

// bars: [{x, y, color?, label?}] on a numeric or categorical x (x an index; o.labels the category names).
export function bars(data, o = {}) {
  const W = o.width || 620, H = o.height || 240, L = 56, R = 14, T = 26, B = o.labels ? 64 : 36;
  const xs = data.map((d) => d.x), ys = data.map((d) => d.y);
  const x0 = Math.min(...xs) - 0.5, x1 = Math.max(...xs) + 0.5;
  const ylo = Math.min(0, ...ys), yhi = Math.max(...ys, 1e-12) * 1.08;
  const X = (x) => L + ((x - x0) / (x1 - x0)) * (W - L - R);
  const Y = (y) => H - B - ((y - ylo) / (yhi - ylo)) * (H - T - B);
  const bw = Math.max(1, (W - L - R) / (x1 - x0) - 2);
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  const axis = svg('g', { class: 'axis' });
  for (const t of ticks(ylo, yhi)) axis.append(svg('line', { x1: L, x2: W - R, y1: Y(t), y2: Y(t), 'stroke-width': 0.6 }), svg('text', { x: L - 6, y: Y(t) + 4, 'text-anchor': 'end' }, fmt(t, 3)));
  if (o.labels) o.labels.forEach((name, i) => axis.append(svg('text', { x: X(i), y: H - B + 12, 'text-anchor': 'end', transform: `rotate(-35 ${X(i)} ${H - B + 12})` }, name)));
  else for (const t of ticks(x0 + 0.5, x1 - 0.5)) axis.append(svg('text', { x: X(t), y: H - B + 16, 'text-anchor': 'middle' }, fmt(t, 3)));
  if (o.xlabel) axis.append(svg('text', { x: (L + W - R) / 2, y: H - 4, 'text-anchor': 'middle' }, o.xlabel));
  if (o.ylabel) axis.append(svg('text', { x: L, y: 12 }, o.ylabel));
  root.append(axis);
  for (const d of data) {
    const top = Math.min(Y(d.y), Y(0)), height = Math.abs(Y(0) - Y(d.y));
    const r = svg('rect', { x: X(d.x) - bw / 2, y: top, width: bw, height, fill: d.color || css('--series-1'), rx: 1.5 });
    hover(r, d.label || `${fmt(d.x)}: ${fmt(d.y)}`);
    root.append(r);
  }
  return root;
}

// points: [{x, y, label, color?, r?, name?}]; o: {xlabel, ylabel, xlim, ylim, names: true to print names}
export function scatter(points, o = {}) {
  const W = o.width || 620, H = o.height || 380, L = 56, R = 18, T = 22, B = 40;
  const xs = points.map((p) => p.x), ys = points.map((p) => p.y);
  const pad = (lo, hi) => { const d = (hi - lo) || 1; return [lo - 0.06 * d, hi + 0.06 * d]; };
  const [x0, x1] = o.xlim || pad(Math.min(...xs), Math.max(...xs));
  const [y0, y1] = o.ylim || pad(Math.min(...ys), Math.max(...ys));
  const X = (x) => L + ((x - x0) / (x1 - x0)) * (W - L - R);
  const Y = (y) => H - B - ((y - y0) / (y1 - y0)) * (H - T - B);
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
  const axis = svg('g', { class: 'axis' });
  for (const t of ticks(y0, y1)) axis.append(svg('line', { x1: L, x2: W - R, y1: Y(t), y2: Y(t), 'stroke-width': 0.6 }), svg('text', { x: L - 6, y: Y(t) + 4, 'text-anchor': 'end' }, fmt(t, 3)));
  for (const t of ticks(x0, x1)) axis.append(svg('line', { x1: X(t), x2: X(t), y1: T, y2: H - B, 'stroke-width': 0.4 }), svg('text', { x: X(t), y: H - B + 16, 'text-anchor': 'middle' }, fmt(t, 3)));
  if (o.xlabel) axis.append(svg('text', { x: (L + W - R) / 2, y: H - 4, 'text-anchor': 'middle' }, o.xlabel));
  if (o.ylabel) axis.append(svg('text', { x: L, y: 12 }, o.ylabel));
  root.append(axis);
  for (const p of points) {
    const c = svg('circle', { cx: X(p.x), cy: Y(p.y), r: p.r || 4.5, fill: p.color || css('--accent'), 'fill-opacity': 0.85, stroke: css('--surface'), 'stroke-width': 1 });
    if (p.label) hover(c, p.label);
    root.append(c);
    if (o.names && p.name) root.append(svg('text', { x: X(p.x) + 6, y: Y(p.y) - 5, class: 'mark-label' }, p.name));
  }
  return root;
}

// A raster of values[row][col] on canvas, NaN transparent; o: {map: name of a colour map, lo, hi, height, rows: labels, cols: labels}
export function raster(values, o = {}) {
  const rows = values.length, cols = rows ? values[0].length : 0;
  const flat = values.flat().filter(Number.isFinite);
  const lo = o.lo ?? Math.min(...flat), hi = o.hi ?? Math.max(...flat);
  const colour = colourMap(o.map || 'plum').colour;
  const canvas = h('canvas', { class: 'raster', width: cols, height: rows, style: `height:${o.height || 220}px` });
  const ctx = canvas.getContext('2d');
  const img = ctx.createImageData(cols, rows);
  for (let r = 0; r < rows; r++)
    for (let c = 0; c < cols; c++) {
      const v = values[r][c], k = 4 * (r * cols + c);
      if (!Number.isFinite(v)) { img.data[k + 3] = 0; continue; }
      const rgb = colour(Math.min(1, Math.max(0, (v - lo) / (hi - lo || 1))));
      img.data[k] = Math.round(rgb[0] * 255), img.data[k + 1] = Math.round(rgb[1] * 255), img.data[k + 2] = Math.round(rgb[2] * 255), img.data[k + 3] = 255;
    }
  ctx.putImageData(img, 0, 0);
  if (o.hover) hover(canvas, () => o.hover.text || '');
  return canvas;
}

// A legend for a colour map: a gradient between two labels.
export function ramp(map, lo, hi, label = '') {
  const colour = colourMap(map).colour;
  const stops = [];
  for (let i = 0; i <= 8; i++) {
    const rgb = colour(i / 8);
    stops.push(`rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`);
  }
  return h('span', { class: 'legend' }, h('span', {}, `${label} ${fmt(lo, 3)}`), h('span', { class: 'ramp', style: `background: linear-gradient(to right, ${stops.join(',')})` }), h('span', {}, fmt(hi, 3)));
}

export function legend(items) {
  return h('div', { class: 'legend' }, items.map((it) => h('span', {}, h('i', { class: it.dot ? 'dot' : '', style: `background:${it.color}` }), it.name)));
}

// The hue of a category by index, from the validated categorical map (its order keeps neighbours apart).
export function category(i) {
  const rgb = colourMap('categorical-apart').colour(i);
  return `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
}
