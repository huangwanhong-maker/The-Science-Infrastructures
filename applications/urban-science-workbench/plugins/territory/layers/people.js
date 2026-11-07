// Where people live, on the whole Earth: hexagons of the H3 grid with the world's people binned on them, each column as
// tall as the square root of its density, rising when they arrive. Resolution 2 from orbit, resolution 3 from nearer: a
// coarser grid holds its shape where a finer one would only shimmer.

import { card, compact, nearest, tone, world } from '../shared.js';

const SETS = [{ res: 2, fade: [5.5e6, 9e7], top: 1.5e6 }, { res: 3, fade: [9e5, 5.5e6], top: 7.5e5 }];

export async function layer(ctx) {
  const { earth, dom: { h } } = ctx;
  const id = ctx.layerId();
  const colour = tone('plum', ctx.theme, { from: 0.16 });
  const made = [];
  const pick = [];
  for (const S of SETS) {
    const w = await world(ctx, S.res);
    const density = w.people.map((p, k) => p / Math.max(1, w.area_km2[k]));
    const sorted = [...density].sort((a, b) => a - b);
    const top = Math.sqrt(sorted[Math.floor(sorted.length * 0.995)] || 1);
    const key = `${id}:${S.res}`;
    const set = earth.addColumns(key, w.rings.map((ring) => ({ ring })), { fade: S.fade, shrink: 0.86 });
    const u = density.map((d) => Math.min(1, Math.sqrt(d) / top));
    set.prisms.set(u.map((x) => S.top * x), u.map((x) => colour(x)));
    made.push(key);
    pick.push({ w, S, near: nearest(w.centres, S.res === 2 ? 8 : 4) });
  }
  const find = (lat, lon) => {
    const range = earth.view.range;
    if (range < 6e5 || range > 9e7) return null;
    const P = pick.find((q) => range >= q.S.fade[0]) || pick[pick.length - 1];
    const k = P.near(lon, lat, P.S.res === 2 ? 3.2 : 1.3);
    return k >= 0 ? { w: P.w, k, P } : null;
  };
  return {
    legend: card(h, 'Where people live', {
      colour, low: 'few', high: 'most crowded', note: 'Columns by the square root of people per km², 2025, on H3 cells.', source: 'GHS-POP R2023A, European Commission JRC',
    }),
    hoverAt(lat, lon) {
      const f = find(lat, lon);
      if (!f) return null;
      const { w, k } = f;
      return `${compact(w.people[k])} people on this cell${w.name[k] ? `, most in ${w.name[k]}` : ''}\n${Math.round(w.people[k] / w.area_km2[k]).toLocaleString('en')} per km²`;
    },
    dispose() { for (const key of made) earth.removeColumns(key); },
  };
}
