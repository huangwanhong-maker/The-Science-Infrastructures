// Where the City Council records floods: its own places and areas of flash flooding, drawn over whatever else is on, as a check of the model
// by what the city has seen.

import { assessment } from '../shared.js';

export async function layer(ctx) {
  const { gis, frame, dom: { h } } = ctx;
  const P = await assessment();
  const points = ctx.layerId(':points'), areas = ctx.layerId(':areas');
  const rose = ctx.theme === 'night' ? [1, 0.62, 0.78] : [0.76, 0.2, 0.46];
  const polys = [];
  for (const a of P.reported_areas) for (const p of a.polygons) polys.push({ ring: p.outer.map(([lon, lat]) => frame.toLocal(lat, lon).slice(0, 2)) });
  gis.addFlats(areas, polys, { ring: (f) => f.ring, color: () => rose, z: 0.7, opacity: 0.32, order: 4, pick: () => 'An area the Council records as flooding' });
  gis.addPoints(points, P.reported_points, { position: (f) => frame.toLocal(f.at[0], f.at[1]).slice(0, 2), radius: () => 14, base: () => 1, height: () => 90, color: () => rose, opacity: 0.95, pick: (f) => f.name || 'A place the Council records as flooding' });
  return {
    legend: h('div', { class: 'card' }, h('b', {}, 'Where the Council records floods'),
      h('div', { class: 'ld-line' }, `${P.reported_points.length} places and ${polys.length} areas of flash flooding`),
      h('div', { class: 'ld-quiet' }, 'Freetown City Council, in the assessment\'s data.')),
    describe: (hit) => ({ kind: 'report', title: hit.feature.name || 'A recorded flood', lines: ['Recorded by the Freetown City Council'] }),
    dispose() { gis.remove(points); gis.remove(areas); },
  };
}
