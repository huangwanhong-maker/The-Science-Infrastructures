// Poverty and inequality over the world's economies: the mosaic of the H3 cells each economy's people live in, the height of a
// column its people, its colour the economy's measure in its latest survey. An economy without a survey stands uncoloured:
// no survey is not a low value.

import { card, compact, nearest, tone, world } from '../shared.js';

const SETS = [{ res: 2, fade: [5.5e6, 9e7], top: 1.5e6 }, { res: 3, fade: [9e5, 5.5e6], top: 7.5e5 }];
const MEASURES = {
  headcount: { title: 'Poverty at $3 a day', map: 'plum', low: '0', high: 'over 75 %', top: 0.75, label: (v) => `${(v * 100).toFixed(1)} % of its people live below $3 a day`, note: 'Each economy, in its latest survey: the share below the international line, in 2021 PPP dollars. Columns by people.' },
  gini: { title: 'Inequality', map: 'orchid', low: 'Gini 25', high: '60', bottom: 0.25, top: 0.6, label: (v) => `Gini index ${(v * 100).toFixed(1)}`, note: 'Each economy, in its latest survey, of income or of consumption, which are not strictly alike. Columns by people.' },
};

export async function layer(ctx) {
  const { earth, dom: { h } } = ctx;
  const M = MEASURES[ctx.spec.measure];
  const id = ctx.layerId();
  const r = await fetch('/data/products/poverty-latest.json');
  if (!r.ok) throw new Error('the survey estimates are not built: the records lamp runs the recipe that fetches them');
  const { economies, source } = await r.json();
  const colour = tone(M.map, ctx.theme, { from: 0.14 });
  const made = [], pick = [];
  for (const S of SETS) {
    const w = await world(ctx, S.res);
    const density = w.people.map((p, k) => p / Math.max(1, w.area_km2[k]));
    const sorted = [...density].sort((a, b) => a - b);
    const top = Math.sqrt(sorted[Math.floor(sorted.length * 0.995)] || 1);
    const key = `${id}:${S.res}`;
    const set = earth.addColumns(key, w.rings.map((ring) => ({ ring })), { fade: S.fade, shrink: 0.9 });
    const values = w.code.map((c) => economies[c]?.[ctx.spec.measure] ?? null);
    set.prisms.set(density.map((d, k) => (values[k] === null ? S.top * 0.02 : S.top * Math.min(1, Math.sqrt(d) / top))),
      values.map((v) => (v === null ? (ctx.theme === 'night' ? [0.13, 0.09, 0.16] : [0.86, 0.82, 0.86]) : colour((v - (M.bottom ?? 0)) / (M.top - (M.bottom ?? 0))))));
    made.push(key);
    pick.push({ w, S, near: nearest(w.centres, S.res === 2 ? 8 : 4) });
  }
  const find = (lat, lon) => {
    const range = earth.view.range;
    if (range < 6e5 || range > 9e7) return null;
    const P = pick.find((q) => range >= q.S.fade[0]) || pick[pick.length - 1];
    const k = P.near(lon, lat, P.S.res === 2 ? 3.2 : 1.3);
    return k >= 0 ? { w: P.w, k } : null;
  };
  return {
    legend: card(h, M.title, { colour, low: M.low, high: M.high, note: M.note, source: source.credit }),
    hoverAt(lat, lon) {
      const f = find(lat, lon);
      if (!f) return null;
      const name = f.w.name[f.k], e = economies[f.w.code[f.k]];
      if (!e) return name ? `${name}: no survey in the Platform` : null;
      return `${name || e.name}, ${e.year}\n${M.label(e[ctx.spec.measure])}\n${compact(f.w.people[f.k])} people on this cell`;
    },
    dispose() { for (const key of made) earth.removeColumns(key); },
  };
}
