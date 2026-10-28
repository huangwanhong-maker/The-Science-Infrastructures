// Small 2-complexes in SVG: vertices on a circle (or a line for two), edges as arcs (parallel edges
// fanned out, loops as small circles), each 2-cell listed by its attaching word and highlighted on its
// edges when hovered; a filled polygon when a single 2-cell runs once round a cycle of all the
// vertices.

import { css, h, svg } from '../dom.js';

export function draw(X, o = {}) {
  const W = o.width || 300, H = o.height || 220, n = X.vertices;
  const R = Math.min(W, H) / 2 - 30;
  const pos = n === 2 ? [[W / 2 - R, H / 2], [W / 2 + R, H / 2]] :
    Array.from({ length: n }, (_, i) => [W / 2 + R * Math.cos((2 * Math.PI * i) / n - Math.PI / 2), H / 2 + R * Math.sin((2 * Math.PI * i) / n - Math.PI / 2)]);
  const root = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H });
  const ink = css('--ink-2'), minor = css('--minor');
  // A polygon filled when one face runs round all the vertices once.
  if (X.faces.length === 1 && X.faces[0].word.length === n && n > 2)
    root.append(svg('polygon', { points: pos.map((p) => p.join(',')).join(' '), fill: css('--accent-soft'), stroke: 'none' }));
  const seen = new Map();
  const paths = X.edges.map(([t, hd], e) => {
    const key = t < hd ? `${t}-${hd}` : `${hd}-${t}`;
    const k = seen.get(key) || 0;
    seen.set(key, k + 1);
    const total = X.edges.filter(([a, b]) => (a < b ? `${a}-${b}` : `${b}-${a}`) === key).length;
    const [ax, ay] = pos[t], [bx, by] = pos[hd];
    let d;
    if (t === hd) {
      // A loop: a teardrop pointing away from the centre, the k-th of a vertex turned and enlarged.
      let ox = ax - W / 2, oy = ay - H / 2;
      const len = Math.hypot(ox, oy);
      if (len < 1) { ox = 0; oy = -1; } else { ox /= len; oy /= len; }
      const turn = (k - (total - 1) / 2) * 0.9, c = Math.cos(turn), s = Math.sin(turn);
      const ux = ox * c - oy * s, uy = ox * s + oy * c, px = -uy, py = ux, reach = 46 + 6 * k, wide = 26 + 4 * k;
      d = `M${ax},${ay} C${ax + ux * reach + px * wide},${ay + uy * reach + py * wide} ${ax + ux * reach - px * wide},${ay + uy * reach - py * wide} ${ax},${ay}`;
    } else {
      const bend = total > 1 ? (k - (total - 1) / 2) * (n === 2 ? 46 : 22) : 0;
      const mx = (ax + bx) / 2, my = (ay + by) / 2, dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy) || 1;
      const sign = t < hd ? 1 : -1;
      d = `M${ax},${ay} Q${mx - (dy / L) * bend * 2 * sign},${my + (dx / L) * bend * 2 * sign} ${bx},${by}`;
    }
    const p = svg('path', { d, fill: 'none', stroke: ink, 'stroke-width': 1.6 });
    root.append(p);
    return p;
  });
  // Edge labels at the middle of each path, and a dot toward its head (after insertion, so that the
  // path can be measured).
  requestAnimationFrame(() => {
    paths.forEach((p, e) => {
      try {
        const L = p.getTotalLength(), m = p.getPointAtLength(L * 0.5), a = p.getPointAtLength(L * 0.62);
        root.append(svg('circle', { cx: a.x, cy: a.y, r: 2.4, fill: ink }));
        root.append(svg('text', { x: m.x + 4, y: m.y - 5, 'font-size': 11, style: `fill:${css('--ink-faint')}` }, `e${e + 1}`));
      } catch { /* not rendered */ }
    });
  });
  pos.forEach(([x, y], i) => {
    root.append(svg('circle', { cx: x, cy: y, r: 9, fill: css('--surface'), stroke: css('--ink'), 'stroke-width': 1.3 }));
    root.append(svg('text', { x, y: y + 4, 'text-anchor': 'middle', 'font-size': 11 }, o.vertexNames ? o.vertexNames[i] : `v${i + 1}`));
  });
  const words = h('div', { class: 'keys' });
  X.faces.forEach((f, i) => {
    const text = f.word.map(([e, s]) => `e${e + 1}${s < 0 ? '⁻¹' : ''}`).join(' ') || '(empty)';
    const chip = h('span', { class: 'chip', style: 'cursor:default' }, `f${i + 1}: ${text}`);
    chip.addEventListener('mouseenter', () => { for (const [e] of f.word) { paths[e].setAttribute('stroke', minor); paths[e].setAttribute('stroke-width', 3); } });
    chip.addEventListener('mouseleave', () => { for (const p of paths) { p.setAttribute('stroke', ink); p.setAttribute('stroke-width', 1.6); } });
    words.append(chip);
  });
  return h('div', {}, root, X.faces.length ? words : h('p', { class: 'caption' }, 'no 2-cells'));
}

// An integer matrix as a small table, signs coloured.
export function matrix(rows, rowNames, colNames) {
  const t = h('table', { class: 'matrix' });
  if (colNames) t.append(h('tr', {}, h('td', {}), ...colNames.map((c) => h('td', { class: 'faint' }, c))));
  rows.forEach((r, i) => t.append(h('tr', {}, h('td', { class: 'faint' }, rowNames ? rowNames[i] : ''), ...r.map((v) => h('td', { class: v < 0 ? 'neg' : v > 0 ? 'pos' : 'faint' }, v === 0 ? '·' : String(v))))));
  return t;
}
