// Labelled graphs in SVG, laid out by a seeded spring embedding: atoms as discs coloured by element,
// nonterminals (labels starting with a capital and holding a lower-case or digit suffix that is no
// element, e.g. A0, Tail1, A_tet) as dashed boxes, as the manuscript draws them; double bonds doubled,
// ports dashed, directed edges with arrows.

import { css, hover, svg } from '../dom.js';

const ELEMENTS = ['H', 'B', 'C', 'N', 'O', 'F', 'S', 'Cl', 'Br'];

export function parseLabel(label) {
  let s = label, charge = '', tau = '';
  const at = s.indexOf('@');
  if (at >= 0) { tau = s.slice(at + 1); s = s.slice(0, at); }
  if (s.endsWith('+') || s.endsWith('-')) { charge = s.slice(-1); s = s.slice(0, -1); }
  const element = ELEMENTS.includes(s) ? s : null;
  const free = s.endsWith('•');
  return { text: s, element, charge, tau, free, nonterminal: !element && /^(A|Tail|A_)/.test(s) };
}

// The connected components of a graph, each a list of its vertices.
function components(n, edges) {
  const p = Array.from({ length: n }, (_, i) => i);
  const root = (x) => (p[x] === x ? x : (p[x] = root(p[x])));
  for (const [u, v] of edges) p[root(u)] = root(v);
  const groups = new Map();
  for (let i = 0; i < n; i++) {
    const r = root(i);
    if (!groups.has(r)) groups.set(r, []);
    groups.get(r).push(i);
  }
  return [...groups.values()];
}

// A spring embedding of one connected component (its vertices `vs`, edges in global numbers), with a
// weak pull toward its centre so that nothing drifts; bonds settle near unit length.
function springs(vs, edges, rand, iterations = 500) {
  const local = new Map(vs.map((v, i) => [v, i]));
  const n = vs.length;
  const P = vs.map((_, i) => [Math.cos((2 * Math.PI * i) / n) * 0.8 + 0.1 * rand(), Math.sin((2 * Math.PI * i) / n) * 0.8 + 0.1 * rand()]);
  const E = edges.filter(([u, v]) => local.has(u) && u !== v).map(([u, v]) => [local.get(u), local.get(v)]);
  for (let it = 0; it < iterations && n > 1; it++) {
    const F = P.map(() => [0, 0]);
    const t = 0.08 * (1 - it / iterations) + 0.004;
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) {
        const dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1];
        const d2 = dx * dx + dy * dy + 1e-4, f = 0.35 / d2;
        F[i][0] += f * dx, F[i][1] += f * dy, F[j][0] -= f * dx, F[j][1] -= f * dy;
      }
    for (const [a, b] of E) {
      const dx = P[a][0] - P[b][0], dy = P[a][1] - P[b][1], d = Math.hypot(dx, dy) + 1e-6, f = (d - 1) * 0.9;
      F[a][0] -= (f * dx) / d, F[a][1] -= (f * dy) / d, F[b][0] += (f * dx) / d, F[b][1] += (f * dy) / d;
    }
    const cx = P.reduce((s, q) => s + q[0], 0) / n, cy = P.reduce((s, q) => s + q[1], 0) / n;
    for (let i = 0; i < n; i++) {
      F[i][0] -= 0.02 * (P[i][0] - cx), F[i][1] -= 0.02 * (P[i][1] - cy);
      const m = Math.hypot(F[i][0], F[i][1]) || 1;
      P[i][0] += (F[i][0] / m) * Math.min(m, t), P[i][1] += (F[i][1] / m) * Math.min(m, t);
    }
  }
  return P;
}

// Turns a component so that its long axis (the principal axis of its points) lies horizontal.
function level(P) {
  const n = P.length;
  const cx = P.reduce((s, q) => s + q[0], 0) / n, cy = P.reduce((s, q) => s + q[1], 0) / n;
  let sxx = 0, syy = 0, sxy = 0;
  for (const [x, y] of P) sxx += (x - cx) ** 2, syy += (y - cy) ** 2, sxy += (x - cx) * (y - cy);
  const t = -0.5 * Math.atan2(2 * sxy, sxx - syy), c = Math.cos(t), s = Math.sin(t);
  return P.map(([x, y]) => [c * (x - cx) - s * (y - cy), s * (x - cx) + c * (y - cy)]);
}

// Each component laid out on its own and levelled, then the components packed side by side, larger
// first, so that a molecule never collapses because another one has drifted far away; single free
// vertices (a released water, an unused unit) are stacked in a column at the end.
function layout(G, seed = 7) {
  let r = seed;
  const rand = () => ((r = (r * 1103515245 + 12345) % 2147483648) / 2147483648);
  const n = G.vertices.length;
  const out = Array.from({ length: n }, () => [0, 0]);
  const parts = components(n, G.edges).sort((a, b) => b.length - a.length);
  const singles = parts.filter((vs) => vs.length === 1).map((vs) => vs[0]);
  let x = 0;
  for (const vs of parts.filter((p) => p.length > 1)) {
    const P = level(springs(vs, G.edges, rand));
    const xs = P.map((q) => q[0]), ys = P.map((q) => q[1]);
    const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
    vs.forEach((v, i) => { out[v] = [x + (P[i][0] - x0), P[i][1] - (y0 + y1) / 2]; });
    x += (x1 - x0) + 1.3;
  }
  singles.forEach((v, i) => { out[v] = [x, (i - (singles.length - 1) / 2) * 0.8]; });
  return out;
}

// G: {vertices: [label], edges: [[u, v, label]], directed}; o: {width, height, highlight: Set of vertex
// ids, ids: show vertex numbers, positions}.
export function draw(G, o = {}) {
  const W = o.width || 260, H = o.height || 180, pad = 26;
  const P = o.positions || layout(G, o.seed);
  const xs = P.map((p) => p[0]), ys = P.map((p) => p[1]);
  const [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
  const s = Math.min((W - 2 * pad) / (x1 - x0 || 1), (H - 2 * pad) / (y1 - y0 || 1), 58);
  const X = (x) => W / 2 + (x - (x0 + x1) / 2) * s, Y = (y) => H / 2 + (y - (y0 + y1) / 2) * s;
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H });
  const ink = css('--ink-2'), accent = css('--accent');
  const marker = svg('marker', { id: `arrow-${Math.random().toString(36).slice(2, 8)}`, viewBox: '0 0 10 10', refX: 15, refY: 5, markerWidth: 6, markerHeight: 6, orient: 'auto-start-reverse' },
    svg('path', { d: 'M0,0 L10,5 L0,10 z', fill: ink }));
  root.append(svg('defs', {}, marker));
  const counts = new Map();
  for (const [u, v, label] of G.edges) {
    const key = u < v ? `${u}-${v}` : `${v}-${u}`;
    const k = counts.get(key) || 0;
    counts.set(key, k + 1);
    const [ax, ay, bx, by] = [X(P[u][0]), Y(P[u][1]), X(P[v][0]), Y(P[v][1])];
    const dx = bx - ax, dy = by - ay, d = Math.hypot(dx, dy) || 1, nx = -dy / d, ny = dx / d;
    const off = k * 6 - (label === 'double' ? 2.2 : 0);
    const hl = o.highlightEdges && o.highlightEdges.has(G.edges.indexOf(G.edges.find((e) => e[0] === u && e[1] === v && e[2] === label)));
    const attrs = { stroke: hl ? accent : ink, 'stroke-width': label === 'amide' ? 3 : 1.5, 'stroke-dasharray': label === 'port' || label.startsWith('port') ? '4 3' : null, 'stroke-linecap': 'round' };
    const line = (o2) => svg('line', { x1: ax + nx * o2, y1: ay + ny * o2, x2: bx + nx * o2, y2: by + ny * o2, ...attrs, 'marker-end': G.directed ? `url(#${marker.id})` : null });
    root.append(line(off));
    if (label === 'double') root.append(line(off + 4.4));
    const tag = svg('title', {}, label);
    root.lastChild.append(tag);
  }
  G.vertices.forEach((label, i) => {
    const L = parseLabel(label);
    const x = X(P[i][0]), y = Y(P[i][1]);
    const g = svg('g', {});
    const hl = o.highlight && o.highlight.has(i);
    if (L.nonterminal) {
      const w = Math.max(30, 7.5 * L.text.length + 10);
      g.append(svg('rect', { x: x - w / 2, y: y - 11, width: w, height: 22, rx: 3, fill: css('--sunken'), stroke: hl ? accent : ink, 'stroke-dasharray': '4 3', 'stroke-width': 1.2 }));
      g.append(svg('text', { x, y: y + 4.5, 'text-anchor': 'middle', 'font-size': 12, 'font-style': 'italic' }, L.text));
    } else if (L.element) {
      const r = L.element === 'H' ? 8 : 11;
      g.append(svg('circle', { cx: x, cy: y, r, fill: css(`--atom-${L.element}`) || css('--surface'), stroke: hl ? accent : css('--hair-2'), 'stroke-width': hl ? 2.4 : 1 }));
      const light = ['H', 'F', 'S', 'B'].includes(L.element);
      g.append(svg('text', { x, y: y + 4, 'text-anchor': 'middle', 'font-size': L.element.length > 1 ? 10 : 11.5, style: `fill:${light ? '#2E2638' : '#FFFFFF'}` }, L.element));
      if (L.charge) g.append(svg('text', { x: x + r + 1, y: y - r + 4, 'font-size': 12, style: `fill:${css('--minor')}` }, L.charge === '+' ? '+' : '−'));
      if (L.tau) {
        const tag = `τ = ${L.tau === '+' ? '+1' : '−1'}`;
        g.append(svg('rect', { x: x + r + 3, y: y + 4, width: 7 * tag.length, height: 16, rx: 3, fill: css('--accent-soft') }));
        g.append(svg('text', { x: x + r + 6, y: y + 16, 'font-size': 12, style: `fill:${css('--accent-deep')}` }, tag));
      }
    } else {
      const w = Math.max(26, 7.2 * L.text.length + 10);
      g.append(svg('rect', { x: x - w / 2, y: y - 10, width: w, height: 20, rx: 4, fill: css('--surface'), stroke: hl ? accent : css('--hair-2'), 'stroke-width': hl ? 2 : 1 }));
      g.append(svg('text', { x, y: y + 4.5, 'text-anchor': 'middle', 'font-size': 11.5 }, L.text));
    }
    if (o.ids) g.append(svg('text', { x, y: y + 23, 'text-anchor': 'middle', 'font-size': 9, style: `fill:${css('--ink-faint')}` }, String(i)));
    hover(g, `${i}: ${label}`);
    root.append(g);
  });
  return root;
}
