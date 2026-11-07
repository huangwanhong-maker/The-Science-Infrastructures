// The plate's neatline: a hairline round the scene, with the graticule's ticks where it crosses it, as on an engraved map.
// Each edge is walked in the page's pixels, the ground under each step found on the globe, and a tick set where the
// longitude or latitude passes a round value. It is drawn when the view has moved, at a step chosen so that about half a
// dozen ticks stand on the longest edge; from far off, where the edges miss the globe, only the line stands.

import { svg } from '../dom.js';

const STEPS = [30, 15, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005];

const format = (v, step, pos, neg) => {
  const a = Math.abs(v);
  const digits = step >= 1 ? 0 : step >= 0.1 ? 1 : step >= 0.01 ? 2 : step >= 0.001 ? 3 : 4;
  return `${a.toFixed(digits)}°${v >= 0 ? pos : neg}`;
};

export function createNeatline(el, earth) {
  let key = '';
  const INSET = 8;
  function draw(W, H, reserve = {}) {
    const v = earth.view;
    const k = [v.lon.toFixed(6), v.lat.toFixed(6), v.range.toFixed(1), v.heading.toFixed(4), v.tilt.toFixed(4), W, H].join();
    if (k === key) return;
    key = k;
    const out = [svg('rect', { class: 'frame', x: INSET + 0.5, y: INSET + 0.5, width: W - 2 * INSET - 1, height: H - 2 * INSET - 1 })];
    const rect = earth.canvas.getBoundingClientRect();
    const N = 80;
    const walk = (edge, x0, y0, x1, y1) => {
      const pts = [];
      for (let i = 0; i <= N; i++) {
        const x = x0 + ((x1 - x0) * i) / N, y = y0 + ((y1 - y0) * i) / N;
        const g = earth.groundAt(rect.left + x, rect.top + y);
        pts.push(g ? { x, y, lon: g.lon, lat: g.lat } : null);
      }
      return { edge, pts };
    };
    const edges = [walk('top', INSET, INSET, W - INSET, INSET), walk('bottom', INSET, H - INSET, W - INSET, H - INSET), walk('left', INSET, INSET, INSET, H - INSET), walk('right', W - INSET, INSET, W - INSET, H - INSET)];
    // The step: from the span of latitude and longitude seen along the edges.
    let lonMin = 400, lonMax = -400, latMin = 100, latMax = -100, hits = 0;
    for (const e of edges) for (const p of e.pts) if (p) { hits++; lonMin = Math.min(lonMin, p.lon), lonMax = Math.max(lonMax, p.lon), latMin = Math.min(latMin, p.lat), latMax = Math.max(latMax, p.lat); }
    if (hits > 40 && lonMax - lonMin < 200) {
      const span = Math.max(lonMax - lonMin, latMax - latMin);
      const step = STEPS.find((s) => span / s <= 9) ?? STEPS[STEPS.length - 1];
      const usedAt = [];
      for (const e of edges) {
        const horizontal = e.edge === 'top' || e.edge === 'bottom';
        for (let i = 1; i < e.pts.length; i++) {
          const a = e.pts[i - 1], b = e.pts[i];
          if (!a || !b) continue;
          for (const [q, pos, neg] of [['lon', 'E', 'W'], ['lat', 'N', 'S']]) {
            if (Math.abs(a[q] - b[q]) > step) continue;
            const ka = Math.floor(a[q] / step), kb = Math.floor(b[q] / step);
            if (ka === kb) continue;
            const v0 = Math.max(ka, kb) * step, t = (v0 - a[q]) / (b[q] - a[q]);
            const x = a.x + (b.x - a.x) * t, y = a.y + (b.y - a.y) * t;
            if (e.edge === 'left' && y < (reserve.leftTop ?? 0)) continue;
            if (usedAt.some(([ux, uy]) => Math.hypot(ux - x, uy - y) < 74)) continue;
            usedAt.push([x, y]);
            const dx = e.edge === 'left' ? 1 : e.edge === 'right' ? -1 : 0, dy = e.edge === 'top' ? 1 : e.edge === 'bottom' ? -1 : 0;
            out.push(svg('line', { class: 'tick', x1: x, y1: y, x2: x + dx * 6, y2: y + dy * 6 }));
            const label = format(v0, step, pos, neg);
            out.push(svg('text', { x: horizontal ? x : x + dx * 9, y: horizontal ? y + dy * 16 + (dy < 0 ? 0 : 0) : y + 3.5, 'text-anchor': horizontal ? 'middle' : dx > 0 ? 'start' : 'end' }, label));
          }
        }
      }
    }
    el.replaceChildren(...out);
    el.setAttribute('viewBox', `0 0 ${W} ${H}`);
  }
  return { draw, invalidate: () => { key = ''; } };
}
