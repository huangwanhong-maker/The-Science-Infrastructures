// The base of the GIS: a city's buildings, streets and open spaces from the data library's city model (in the
// area's local east–north–up frame; Melbourne's from OpenStreetMap, New York's from Valdrada's compiled city of
// the city's records), drawn as extruded solids, ribbons of their widths, flat areas and, where the model has
// them, trees. A building whose height was not mapped stands at a default height, drawn paler so that the guess
// shows as one. The legend credits each area's own sources.

const PALETTES = {
  day: { building: [0.95, 0.91, 0.95], tall: [0.84, 0.75, 0.89], guessed: [0.98, 0.96, 0.98], pedestrian: [0.97, 0.80, 0.87], road: [0.88, 0.84, 0.90], minor: [0.93, 0.90, 0.94],
    park: [0.78, 0.86, 0.76], water: [0.70, 0.70, 0.92], plaza: [0.97, 0.92, 0.94], rail: [0.84, 0.80, 0.86], tree: [0.58, 0.74, 0.60], background: '#faf4fa', ground: '#f3e9f3' },
  night: { building: [0.30, 0.18, 0.38], tall: [0.62, 0.34, 0.66], guessed: [0.215, 0.14, 0.27], pedestrian: [0.62, 0.20, 0.42], road: [0.075, 0.045, 0.105], minor: [0.05, 0.032, 0.072],
    park: [0.055, 0.125, 0.11], water: [0.03, 0.026, 0.115], plaza: [0.11, 0.065, 0.14], rail: [0.08, 0.055, 0.105], tree: [0.09, 0.2, 0.16], background: '#150d1b', ground: '#1c1129', css: { building: 'rgb(100,74,120)', guessed: 'rgb(78,58,96)', pedestrian: 'rgb(214,120,170)', park: 'rgb(62,104,96)', water: 'rgb(48,52,120)' } },
};

// Walls by theme: at night their windows come on, in the rose of the page and a little lamplight.
export const facadeOf = (theme) => (theme === 'night' ? { foot: 0.55, glow: 0.62 } : { foot: 0.78, glow: 0.5, pane: [0.76, 0.28, 0.5], warm: [0.9, 0.55, 0.3] });
export const palette = (theme) => ({ ...PALETTES[theme === 'night' ? 'night' : 'day'], facade: facadeOf(theme) });

const PEDESTRIAN = new Set(['pedestrian', 'footway', 'path', 'steps', 'living_street', 'cycleway', 'corridor']);
const MAJOR = new Set(['motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'motorway_link', 'trunk_link', 'primary_link', 'secondary_link']);

export function addCity(gis, model, theme) {
  const P = palette(theme);
  const mix = (a, b, t) => a.map((v, i) => v + (b[i] - v) * t);
  gis.addFlats('spaces', model.spaces || [], {
    ring: (f) => f.ring,
    holes: (f) => f.holes || [],
    color: (f) => P[f.kind] || P.plaza,
    z: 0.03,
    pick: false,
  });
  gis.addRibbons('streets', (model.roads || []).filter((r) => r.line && r.line.length >= 2), {
    line: (f) => f.line,
    width: (f) => Math.max(1.2, f.width_m || 6),
    color: (f) => (PEDESTRIAN.has(f.kind) ? P.pedestrian : MAJOR.has(f.kind) ? P.road : P.minor),
    z: 0.08,
    pick: (f) => f.name || f.kind,
  });
  gis.addSolids('buildings', model.buildings || [], {
    ring: (f) => f.footprint,
    holes: (f) => f.holes || [],
    base: (f) => f.min_height_m || 0,
    height: (f) => f.height_m || 9,
    color: (f) => (f.height_source === 'default' ? P.guessed : mix(P.building, P.tall, Math.min(1, (f.height_m || 0) / 120))),
    facade: P.facade,
  });
  // Trees, where the model has them (New York's): a crown on its trunk, its breadth by the trunk's diameter.
  if ((model.trees || []).length)
    gis.addPoints('trees', model.trees, {
      position: (f) => f.at,
      radius: (f) => Math.min(4.5, 1.2 + 0.11 * (f.dbh_in || 8)),
      base: () => 2.6,
      height: (f) => Math.min(9, 3 + 0.18 * (f.dbh_in || 8)),
      color: (f, k) => P.tree.map((v) => v * (0.94 + 0.12 * ((k * 2654435761) % 97) / 97)),
      opacity: 0.92,
      pick: (f) => `a tree: ${f.species}${f.dbh_in ? `, ${f.dbh_in} in across the trunk` : ''} (${f.kind})`,
    });
  const all = [];
  for (const b of model.buildings || []) for (const p of b.footprint || []) all.push(p);
  for (const r of model.roads || []) for (const p of r.line || []) all.push(p);
  return all.length ? gis.extentOf(all) : null;
}

export function legendOfCity(dom, model, theme) {
  const { h } = dom;
  const P = palette(theme);
  const css = (rgb, key) => (P.css && P.css[key]) || `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
  const n = (model.buildings || []).length, guessed = (model.buildings || []).filter((b) => b.height_source === 'default').length;
  return h('div', { class: 'card' }, h('b', {}, 'The city'),
    h('div', { class: 'legend' }, ...[['buildings', P.building, 'building'], guessed && ['height not mapped', P.guessed, 'guessed'], ['pedestrian ways', P.pedestrian, 'pedestrian'], ['parks', P.park, 'park'], ['water', P.water, 'water']]
      .filter(Boolean).map(([name, c, k]) => h('span', {}, h('i', { class: 'dot', style: `background:${css(c, k)}` }), name))),
    h('div', { class: 'faint', style: 'margin-top:3px' }, `${n.toLocaleString('en')} buildings${guessed ? ` (${Math.round((100 * guessed) / Math.max(1, n))}% at a default height)` : ''}, ${(model.roads || []).length.toLocaleString('en')} ways`),
    h('div', { class: 'faint' }, model.credit || '© OpenStreetMap contributors'));
}

// The card of a place drawn from the world's archives: what is standing in view, how many of those heights are known, and whose they are.
export function worldCard(dom, theme) {
  const { h } = dom;
  const P = palette(theme);
  const css = (rgb, key) => (P.css && P.css[key]) || `rgb(${rgb.map((c) => Math.round(Math.min(1, c) * 255)).join(',')})`;
  const count = h('div', { class: 'world-count' }, 'asking for the buildings here…');
  const share = h('div', { class: 'world-share' });
  const card = h('div', { class: 'card city-card' }, h('b', {}, 'The city'),
    h('div', { class: 'legend' }, ...[['height known', P.tall, 'building'], ['height estimated', P.guessed, 'guessed'], ['ways on foot', P.pedestrian, 'pedestrian'], ['parks', P.park, 'park'], ['water', P.water, 'water']]
      .map(([name, c, k]) => h('span', {}, h('i', { class: 'dot', style: `background:${css(c, k)}` }), name))),
    count, share,
    h('div', { class: 'faint' }, '© OpenStreetMap contributors, Overture Maps Foundation; footprints also from Google Open Buildings and Microsoft. Estimated heights: GHS-BUILT-H R2023A, the average of the buildings on a 100 m cell in 2018.'));
  return {
    card,
    update(s) {
      if (!s.buildings) { count.textContent = s.pending ? 'asking for the buildings here…' : 'no buildings are recorded here'; share.replaceChildren(); return; }
      count.textContent = `${s.buildings.toLocaleString('en')} buildings in view${s.pending ? ', more coming' : ''}`;
      const known = s.measured + s.floors, all = Math.max(1, s.buildings);
      share.replaceChildren(h('i', { style: `width:${(100 * s.measured) / all}%` }), h('i', { style: `width:${(100 * s.floors) / all}%` }), h('em', {}, `${Math.round((100 * known) / all)}% with a recorded height or floors`));
    },
  };
}
