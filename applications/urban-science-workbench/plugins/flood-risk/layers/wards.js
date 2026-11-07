// What each ward stands to lose: the assessment's average annual loss by ward, raised on its own ground, the height by the loss per square
// kilometre so that a large ward does not stand tall only for being large.

import { assessment, money } from '../shared.js';
import { colourMap } from '/visualization/styles/index.js';
import { ringArea } from '/visualization/js/gis/geometry.js';

export async function layer(ctx) {
  const { gis, frame, dom: { h } } = ctx;
  const P = await assessment();
  const id = ctx.layerId();
  const map = colourMap('orchid').colour;
  const colour = (u) => map(ctx.theme === 'night' ? 0.25 + 0.75 * u : 1 - (0.1 + 0.8 * u));
  const features = [];
  for (const w of P.wards) for (const p of w.polygons) {
    const outer = p.outer.map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2));
    features.push({ ward: w, outer, holes: p.holes.map((q) => q.map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2))), area: Math.abs(ringArea(outer)) });
  }
  const area = new Map();
  for (const f of features) area.set(f.ward, (area.get(f.ward) || 0) + f.area);
  const density = (f) => (f.ward.aal || 0) / Math.max(1e4, area.get(f.ward)) * 1e6;
  const sorted = features.map(density).sort((a, b) => a - b);
  const top = sorted[Math.floor(sorted.length * 0.98)] || 1;
  const u = (f) => Math.sqrt(Math.min(1, density(f) / top));
  gis.addSolids(id, features, {
    ring: (f) => f.outer, holes: (f) => f.holes, base: () => 0.5, height: (f) => 2 + u(f) * 420, color: (f) => colour(u(f)), opacity: 0.86,
    pick: (f) => `${f.ward.name}\n${money(f.ward.aal || 0)} a year on average\n${money(density(f))} a year per km²`,
  });
  const ranked = [...P.wards].sort((a, b) => (b.aal || 0) - (a.aal || 0));
  const total = P.wards.reduce((a, w) => a + (w.aal || 0), 0);
  return {
    legend: h('div', { class: 'card' }, h('b', {}, 'What each ward stands to lose'),
      h('div', { class: 'ld-line' }, `${money(total)} a year over ${P.wards.length} wards; most in ${ranked.slice(0, 3).map((w) => w.name).join(', ')}`),
      h('div', { class: 'ld-quiet' }, 'Height: the average annual loss per km², by the square root. From the assessment\'s model, of its modelled buildings.')),
    describe: (hit) => ({ kind: 'ward', title: hit.feature.ward.name, lines: [`${money(hit.feature.ward.aal || 0)} a year on average`, `${money(density(hit.feature))} a year per km²`, hit.feature.ward.code] }),
    dispose() { gis.remove(id); },
  };
}
