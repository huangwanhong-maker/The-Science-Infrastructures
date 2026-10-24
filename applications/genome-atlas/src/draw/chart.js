// Small charts as SVG strings, in the atlas's hand: hairline axes, serif labels, few ticks.
//
//   chart({ width, height, x: [lo, hi], y: [lo, hi], xlabel, ylabel, series: [...], bars: [...] })
//
// A series is { points: [[x, y], ...], colour, dash, width, label, area, dots, line }; bars are
// { x0, x1, y, colour, title }. Scales are linear, or logarithmic with log: 'x' or 'y'.

import { esc } from '../dom.js';
import { token } from '../palette.js';
import { tickStep, number } from '../format.js';

export function chart({ width = 520, height = 260, x, y, xlabel = '', ylabel = '', series = [], bars = [], marks = [], log = '', xticks, yticks, title = '' }) {
  const m = { l: 56, r: 14, t: title ? 26 : 12, b: 40 };
  const W = width - m.l - m.r, H = height - m.t - m.b;
  const lg = (v) => Math.log10(Math.max(v, 1e-300));
  const sx = log.includes('x') ? (v) => m.l + ((lg(v) - lg(x[0])) / (lg(x[1]) - lg(x[0]))) * W : (v) => m.l + ((v - x[0]) / (x[1] - x[0] || 1)) * W;
  const sy = log.includes('y') ? (v) => m.t + H - ((lg(v) - lg(y[0])) / (lg(y[1]) - lg(y[0]))) * H : (v) => m.t + H - ((v - y[0]) / (y[1] - y[0] || 1)) * H;
  const faint = token('--ink-faint'), hair = token('--hair'), hair2 = token('--hair-2');
  const out = [];
  if (title) out.push(`<text x="${m.l}" y="14" style="font-size:13px" fill="${token('--ink-2')}">${esc(title)}</text>`);
  // Grid and ticks.
  const ticks = (lo, hi, isLog, given) => {
    if (given) return given;
    if (isLog) {
      const t = [];
      for (let e = Math.ceil(lg(lo)); e <= Math.floor(lg(hi)); e++) t.push(10 ** e);
      return t;
    }
    const step = tickStep(hi - lo, 5), t = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) t.push(+v.toPrecision(12));
    return t;
  };
  for (const v of ticks(y[0], y[1], log.includes('y'), yticks)) {
    const yy = sy(v);
    out.push(`<line x1="${m.l}" y1="${yy.toFixed(1)}" x2="${m.l + W}" y2="${yy.toFixed(1)}" stroke="${hair}" />`);
    out.push(`<text class="tick-label" x="${m.l - 6}" y="${(yy + 3.5).toFixed(1)}" text-anchor="end">${esc(number(v, 3))}</text>`);
  }
  for (const v of ticks(x[0], x[1], log.includes('x'), xticks)) {
    const xx = sx(v);
    out.push(`<line x1="${xx.toFixed(1)}" y1="${m.t + H}" x2="${xx.toFixed(1)}" y2="${m.t + H + 4}" stroke="${faint}" />`);
    out.push(`<text class="tick-label" x="${xx.toFixed(1)}" y="${m.t + H + 16}" text-anchor="middle">${esc(number(v, 3))}</text>`);
  }
  out.push(`<line x1="${m.l}" y1="${m.t + H}" x2="${m.l + W}" y2="${m.t + H}" stroke="${hair2}" />`);
  out.push(`<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${m.t + H}" stroke="${hair2}" />`);
  if (xlabel) out.push(`<text class="tick-label" x="${m.l + W / 2}" y="${height - 6}" text-anchor="middle" style="font-style:italic">${esc(xlabel)}</text>`);
  if (ylabel) out.push(`<text class="tick-label" transform="translate(13 ${m.t + H / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">${esc(ylabel)}</text>`);
  // Bars.
  for (const b of bars) {
    const x0 = sx(b.x0), x1 = sx(b.x1), y0 = sy(b.y0 ?? (log.includes('y') ? y[0] : Math.max(0, y[0]))), y1 = sy(b.y);
    out.push(`<rect x="${Math.min(x0, x1).toFixed(1)}" y="${Math.min(y0, y1).toFixed(1)}" width="${Math.max(0.5, Math.abs(x1 - x0) - (b.gap ?? 1)).toFixed(1)}" height="${Math.max(0.5, Math.abs(y1 - y0)).toFixed(1)}" fill="${b.colour || token('--accent')}" fill-opacity="${b.opacity ?? 0.8}">${b.title ? `<title>${esc(b.title)}</title>` : ''}</rect>`);
  }
  // Series.
  for (const s of series) {
    const pts = s.points.filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b) && (!log.includes('y') || b > 0) && (!log.includes('x') || a > 0));
    if (!pts.length) continue;
    const d = pts.map(([a, b], k) => `${k ? 'L' : 'M'}${sx(a).toFixed(1)},${sy(b).toFixed(1)}`).join('');
    if (s.area) out.push(`<path d="${d}V${(m.t + H).toFixed(1)}H${sx(pts[0][0]).toFixed(1)}Z" fill="${s.colour}" fill-opacity="0.12" />`);
    if (s.line !== false) out.push(`<path d="${d}" fill="none" stroke="${s.colour || token('--accent')}" stroke-width="${s.width || 1.6}" ${s.dash ? `stroke-dasharray="${s.dash}"` : ''} stroke-linejoin="round" />`);
    if (s.dots) for (const [a, b] of pts) out.push(`<circle cx="${sx(a).toFixed(1)}" cy="${sy(b).toFixed(1)}" r="3" fill="${s.colour}" />`);
  }
  for (const k of marks) out.push(`<text class="tick-label" x="${sx(k.x).toFixed(1)}" y="${(sy(k.y) - 6).toFixed(1)}" text-anchor="middle">${esc(k.text)}</text>`);
  // A legend of the labelled series, top right.
  let ly = m.t + 4;
  for (const s of series.filter((q) => q.label)) {
    out.push(`<line x1="${m.l + W - 120}" y1="${ly + 4}" x2="${m.l + W - 104}" y2="${ly + 4}" stroke="${s.colour}" stroke-width="2" ${s.dash ? `stroke-dasharray="${s.dash}"` : ''} />`);
    out.push(`<text class="tick-label" x="${m.l + W - 98}" y="${ly + 8}">${esc(s.label)}</text>`);
    ly += 15;
  }
  return `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}

// A timeline: rows of bars [start, end] with colours, for workers and threads.
export function timeline({ width = 760, rows, bars, seconds, rowHeight = 22, label = (r) => r }) {
  const m = { l: 118, r: 14, t: 8, b: 30 };
  const W = width - m.l - m.r;
  const height = m.t + rows.length * rowHeight + m.b;
  const sx = (t) => m.l + (t / Math.max(1e-9, seconds)) * W;
  const out = [];
  rows.forEach((r, k) => {
    const y = m.t + k * rowHeight;
    out.push(`<rect x="${m.l}" y="${y + 2}" width="${W}" height="${rowHeight - 4}" fill="${token('--band')}" />`);
    out.push(`<text class="tick-label" x="${m.l - 8}" y="${y + rowHeight / 2 + 4}" text-anchor="end">${esc(label(r))}</text>`);
  });
  for (const b of bars) {
    const y = m.t + b.row * rowHeight;
    out.push(`<rect x="${sx(b.start).toFixed(1)}" y="${y + 4}" width="${Math.max(1, sx(b.end) - sx(b.start) - 1).toFixed(1)}" height="${rowHeight - 8}" rx="2" fill="${b.colour}" fill-opacity="${b.opacity ?? 0.85}">${b.title ? `<title>${esc(b.title)}</title>` : ''}</rect>`);
  }
  const step = tickStep(seconds, 6);
  for (let t = 0; t <= seconds + 1e-9; t += step) {
    const x = sx(t);
    out.push(`<line x1="${x.toFixed(1)}" y1="${m.t + rows.length * rowHeight}" x2="${x.toFixed(1)}" y2="${m.t + rows.length * rowHeight + 4}" stroke="${token('--ink-faint')}" />`);
    out.push(`<text class="tick-label" x="${x.toFixed(1)}" y="${m.t + rows.length * rowHeight + 16}" text-anchor="middle">${esc(number(t, 3))} s</text>`);
  }
  return `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg">${out.join('')}</svg>`;
}
