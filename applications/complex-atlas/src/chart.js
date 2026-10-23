// A small line chart in SVG: a series against its index, with an exact value and the mean drawn
// across it.

import { svg } from './dom.js';
import { number } from './format.js';

export function lineChart(values, { exact = null, mean = null, label = '' } = {}) {
  const W = 340, H = 150, L = 40, R = 8, T = 10, B = 22;
  const all = [...values, ...(exact !== null ? [exact] : []), ...(mean !== null ? [mean] : [])];
  let lo = Math.min(...all), hi = Math.max(...all);
  if (hi - lo < 1e-9) { lo -= 0.5; hi += 0.5; }
  const padY = (hi - lo) * 0.08;
  lo -= padY; hi += padY;
  const x = (i) => L + (i / Math.max(1, values.length - 1)) * (W - L - R), y = (v) => T + (1 - (v - lo) / (hi - lo)) * (H - T - B);
  const out = svg('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}`, preserveAspectRatio: 'none' });
  out.append(svg('line', { x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: 'rgba(240,226,238,0.2)' }));
  for (const v of [lo + padY, hi - padY]) out.append(svg('text', { x: L - 4, y: y(v) + 4, 'text-anchor': 'end' }, number(v, 3)));
  out.append(svg('text', { x: W - R, y: H - 6, 'text-anchor': 'end' }, `${values.length} ${label}`));
  if (exact !== null) out.append(svg('line', { x1: L, x2: W - R, y1: y(exact), y2: y(exact), stroke: '#f2c27c', 'stroke-dasharray': '4 3', 'stroke-width': 1.2 }));
  if (mean !== null) out.append(svg('line', { x1: L, x2: W - R, y1: y(mean), y2: y(mean), stroke: '#b9a2e0', 'stroke-width': 1 }));
  const step = Math.max(1, Math.floor(values.length / 600));
  const pts = [];
  for (let i = 0; i < values.length; i += step) pts.push(`${x(i).toFixed(1)},${y(values[i]).toFixed(1)}`);
  out.append(svg('polyline', { points: pts.join(' '), fill: 'none', stroke: '#C08081', 'stroke-width': 1.3, 'vector-effect': 'non-scaling-stroke' }));
  return out;
}
