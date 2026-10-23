// Small drawings, in SVG, of the graphs of the generation of spacetime: the bubbles of an action
// (its interactions), the Feynman graphs of its expansion, the vertex of a simplicial model, and
// the boundary spin network of a spin-foam vertex. Black vertices are the fields φ, white ones
// φ̄; each index (colour) is drawn in its warm hue, the propagators (colour 0) as pale hairlines.

import { svg } from './dom.js';
import { colouredGraphLayout, indexColour, cssColour } from '../../../visualization/js/index.js';

const BLACK = '#2a1f2f', WHITE = '#f6eef5', RING = '#e8c3dc';

// A coloured graph, as a bubble (propagator = false: its colours are the indices 1..r) or as a
// Feynman graph (colour 0 the propagators).
export function graphSvg(cg, { size = 120, propagator = true, seed = 1, title = null } = {}) {
  const out = svg('svg', { class: 'graph', viewBox: `0 0 ${size} ${size}`, width: size, height: size, role: 'img' });
  if (title) out.append(svg('title', {}, title));
  const pad = size * 0.14, span = size - 2 * pad;
  if (!cg.colours) {
    // Rank 0: the vertices alone, pairs side by side.
    for (let b = 0; b < cg.p; b++) {
      const y = pad + (span * (b + 0.5)) / cg.p;
      out.append(svg('line', { x1: pad + span * 0.25, x2: pad + span * 0.75, y1: y, y2: y, stroke: 'rgba(240,226,238,0.25)', 'stroke-dasharray': '2 3' }));
      out.append(vertex(pad + span * 0.25, y, true, size), vertex(pad + span * 0.75, y, false, size));
    }
    return out;
  }
  const L = colouredGraphLayout(cg, { seed });
  const P = L.nodes.map((n) => [pad + n.x * span, pad + n.y * span]);
  for (const e of L.edges) {
    const [x1, y1] = P[e.u], [x2, y2] = P[e.v];
    const colour = cssColour(indexColour(e.colour, { propagator }));
    const mx = (x1 + x2) / 2, my = (y1 + y2) / 2, dx = x2 - x1, dy = y2 - y1, len = Math.hypot(dx, dy) || 1;
    const bend = e.bend * Math.min(len * 0.9, size * 0.32);
    const cx = mx - (dy / len) * bend, cy = my + (dx / len) * bend;
    out.append(svg('path', { d: `M${x1},${y1} Q${cx},${cy} ${x2},${y2}`, fill: 'none', stroke: colour,
      'stroke-width': propagator && e.colour === 0 ? 1.1 : 1.9, 'stroke-opacity': propagator && e.colour === 0 ? 0.75 : 0.95 }));
  }
  L.nodes.forEach((n, i) => out.append(vertex(P[i][0], P[i][1], n.black, size)));
  return out;
}

function vertex(x, y, black, size) {
  const r = Math.max(2.6, size * 0.035);
  return svg('circle', { cx: x, cy: y, r, fill: black ? BLACK : WHITE, stroke: RING, 'stroke-width': 1.1 });
}

// The vertex of a coloured simplicial model: a black (λ ∏ φ_c) and a white (λ̄ ∏ φ̄_c) vertex
// with one half-edge of each colour 0..D.
export function simplicialVertexSvg(colours, { size = 120 } = {}) {
  const out = svg('svg', { class: 'graph', viewBox: `0 0 ${size} ${size}`, width: size, height: size, role: 'img' });
  [[size * 0.3, true], [size * 0.7, false]].forEach(([cx, black]) => {
    const cy = size / 2, R = size * 0.2;
    for (let c = 0; c < colours; c++) {
      const a = -Math.PI / 2 + (2 * Math.PI * c) / colours;
      out.append(svg('line', { x1: cx, y1: cy, x2: cx + R * Math.cos(a), y2: cy + R * Math.sin(a), stroke: cssColour(indexColour(c, { propagator: false })), 'stroke-width': 2 }));
      out.append(svg('text', { x: cx + (R + 7) * Math.cos(a), y: cy + (R + 7) * Math.sin(a) + 3, 'text-anchor': 'middle', class: 'mark' }, String(c)));
    }
    out.append(vertex(cx, cy, black, size));
  });
  return out;
}

// The boundary spin network of a spin-foam vertex in three dimensions: the tetrahedral graph
// whose nodes are the four triangles T_k (opposite vertex k of the tetrahedron) and whose six
// links are the edges ab of the tetrahedron, joining the two triangles that contain ab, each
// labelled with its spin. Edges in the order 01, 02, 03, 12, 13, 23.
export const TETRAHEDRON_EDGES = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

export function spinNetworkSvg(spins, { size = 200, spinText = (s) => String(s), colourOf = null } = {}) {
  const out = svg('svg', { class: 'graph spin-network', viewBox: `0 0 ${size} ${size}`, width: size, height: size, role: 'img' });
  // T0 at the centre, T1, T2, T3 at the corners of a triangle. The spins of the three spokes
  // are written two thirds of the way out, those of the outer sides at their middles.
  const c = size / 2, R = size * 0.4, cy = c + size * 0.05;
  const node = [[c, cy], ...[0, 1, 2].map((k) => [c + R * Math.cos(-Math.PI / 2 + (2 * Math.PI * k) / 3), cy + R * Math.sin(-Math.PI / 2 + (2 * Math.PI * k) / 3)])];
  TETRAHEDRON_EDGES.forEach(([a, b], i) => {
    const [u, v] = [0, 1, 2, 3].filter((k) => k !== a && k !== b);  // the triangles containing edge ab
    const [x1, y1] = node[u], [x2, y2] = node[v];
    const colour = colourOf ? colourOf(spins[i]) : '#C08081';
    out.append(svg('line', { x1, y1, x2, y2, stroke: colour, 'stroke-width': 2.2, 'stroke-linecap': 'round' }));
    const t = u === 0 ? 0.64 : 0.5;
    const mx = x1 + (x2 - x1) * t, my = y1 + (y2 - y1) * t;
    out.append(svg('rect', { x: mx - 15, y: my - 8, width: 30, height: 15, rx: 7, fill: '#1d1520', stroke: 'rgba(240,226,238,0.2)' }));
    out.append(svg('text', { x: mx, y: my + 3, 'text-anchor': 'middle', class: 'spin' }, spinText(spins[i])));
  });
  node.forEach(([x, y], k) => {
    out.append(svg('circle', { cx: x, cy: y, r: 6.5, fill: '#2a1f2f', stroke: '#f2c27c', 'stroke-width': 1.4 }));
    // Outer labels outward from the centre; the centre's below-left of it.
    const dx = k ? (x - c) / R : -0.75, dy = k ? (y - cy) / R : 0.95;
    out.append(svg('text', { x: x + dx * 15, y: y + dy * 15 + 3.5, 'text-anchor': 'middle', class: 'mark' }, `T${k}`));
  });
  return out;
}
