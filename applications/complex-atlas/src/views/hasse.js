// The Hasse diagram: the cells in layers by dimension, vertices at the bottom, each joined to the
// cells of its boundary (hasse() in visualization/js/complex.js). An incidence of multiplicity
// other than ±1 is written on its link (a 2-cell attached by a a meets a twice); a loop meets its
// vertex with incidence 0, drawn dashed. It pans and zooms as the net does.

import { svg, clear } from '../dom.js';
import { attachPanZoom } from './net.js';
import { selectionSets } from '../selection.js';
import { cellName } from '../format.js';

const MARK = '#f6cf7e';
const LAYER = ['vertices', 'edges', 'faces'];
export const HASSE_LIMIT = 3000;

export function drawHasse(host, complex, result, { names, style, selection, onHover, fit }) {
  clear(host);
  const counts = result.layers.map((l) => l.length);
  const total = counts.reduce((s, n) => s + n, 0);
  const W = Math.max(900, Math.min(4000, Math.max(...counts) * 26)), gap = Math.max(170, W * 0.55 / Math.max(1, counts.length - 1)), H = gap * (counts.length - 1) + 120;
  if (fit || !host.__view) host.__view = { x: -150, y: -30, w: W + 200, h: H + 60 };
  host.setAttribute('viewBox', `${host.__view.x} ${host.__view.y} ${host.__view.w} ${host.__view.h}`);
  attachPanZoom(host);
  const marks = selectionSets(selection, complex);
  const markSet = [marks.vertices, marks.edges, marks.faces];
  const X = (d, c) => 40 + result.nodes[d][c].x * (W - 80), Y = (d) => H - 60 - d * gap;
  const linkLayer = svg('g', { stroke: '#cfb9ca', 'stroke-opacity': total > 600 ? 0.14 : 0.3, 'stroke-width': 1 });
  const labelLayer = svg('g', {});
  for (const l of result.links) {
    const lower = l.point ? 0 : l.dim - 1;
    if (l.dim >= counts.length || (l.point ? l.face >= counts[0] : l.face >= counts[lower])) continue;
    const a = [X(l.dim, l.cell), Y(l.dim)], b = [X(lower, l.face), Y(lower)];
    const marked = markSet[l.dim]?.has(l.cell) && (l.point || markSet[lower]?.has(l.face) || l.dim === 1);
    linkLayer.append(svg('line', { x1: a[0], y1: a[1], x2: b[0], y2: b[1], 'stroke-dasharray': l.coefficient === 0 ? '4 4' : null, stroke: marked ? MARK : null, 'stroke-opacity': marked ? 0.9 : null }));
    if (Math.abs(l.coefficient) > 1 && total <= 600) labelLayer.append(svg('text', { x: (a[0] + b[0]) / 2 + 4, y: (a[1] + b[1]) / 2, 'font-size': 12, fill: '#e9b872', 'font-style': 'italic' }, `${l.coefficient > 0 ? '' : '−'}${Math.abs(l.coefficient)}`));
  }
  host.append(linkLayer);
  const nodeLayer = svg('g', {});
  const r = Math.max(2.5, Math.min(7, 7 * style.vertexSize * (total > 600 ? 0.5 : 1)));
  result.layers.forEach((order, d) => {
    labelLayer.append(svg('text', { x: -20, y: Y(d) + 5, 'font-size': 16, fill: '#9a8496', 'text-anchor': 'end', 'font-style': 'italic' }, `${LAYER[d] || `${d}-cells`}`));
    for (const c of order) {
      const marked = markSet[d]?.has(c);
      nodeLayer.append(svg('circle', { cx: X(d, c), cy: Y(d), r: marked ? r * 1.6 : r, fill: marked ? MARK : ['#f2c27c', '#ecd4e2', '#c08081'][d] ?? '#b9a2e0', stroke: '#17111a', 'stroke-width': 1, onpointerenter: () => onHover?.({ dim: d, index: c }) }));
      if (style.labels && order.length <= 60) labelLayer.append(svg('text', { x: X(d, c), y: Y(d) - r - 6, 'font-size': 15, fill: '#cfb9ca', 'text-anchor': 'middle', 'font-style': d === 1 ? 'italic' : null }, cellName(d, c, names)));
    }
  });
  host.append(nodeLayer, labelLayer);
  host.onpointerleave = () => onHover?.(null);
}
