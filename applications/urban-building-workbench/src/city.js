// The base of the GIS: a city's buildings, streets and open spaces from the data library's city model (in the
// area's local east–north–up frame; Melbourne's from OpenStreetMap, New York's from Valdrada's compiled city of
// the city's records), drawn as extruded solids, ribbons of their widths, flat areas and, where the model has
// them, trees. A building whose height was not mapped stands at a default height, drawn paler so that the guess
// shows as one. The legend credits each area's own sources.

const PALETTES = {
  day: { building: [0.91, 0.905, 0.89], tall: [0.75, 0.78, 0.8], guessed: [0.955, 0.955, 0.945], pedestrian: [0.9, 0.84, 0.76], road: [0.84, 0.86, 0.87], minor: [0.89, 0.9, 0.9],
    park: [0.76, 0.84, 0.7], water: [0.6, 0.75, 0.86], plaza: [0.9, 0.9, 0.87], rail: [0.81, 0.81, 0.81], tree: [0.56, 0.7, 0.52], background: '#eef3f6', ground: '#e4eaee' },
  night: { building: [0.3, 0.33, 0.37], tall: [0.39, 0.43, 0.48], guessed: [0.24, 0.27, 0.3], pedestrian: [0.45, 0.4, 0.36], road: [0.26, 0.3, 0.34], minor: [0.22, 0.25, 0.29],
    park: [0.22, 0.3, 0.26], water: [0.17, 0.27, 0.35], plaza: [0.27, 0.29, 0.31], rail: [0.25, 0.27, 0.3], tree: [0.25, 0.36, 0.28], background: '#141d26', ground: '#1b2630' },
};

export const palette = (theme) => PALETTES[theme === 'night' ? 'night' : 'day'];

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
  const css = (rgb) => `rgb(${rgb.map((c) => Math.round(c * 255)).join(',')})`;
  const n = (model.buildings || []).length, guessed = (model.buildings || []).filter((b) => b.height_source === 'default').length;
  return h('div', { class: 'card' }, h('b', {}, 'The city'),
    h('div', { class: 'legend' }, ...[['buildings', P.building], guessed && ['height not mapped', P.guessed], ['pedestrian ways', P.pedestrian], ['parks', P.park], ['water', P.water]]
      .filter(Boolean).map(([name, c]) => h('span', {}, h('i', { class: 'dot', style: `background:${css(c)}` }), name))),
    h('div', { class: 'faint', style: 'margin-top:3px' }, `${n.toLocaleString('en')} buildings${guessed ? ` (${Math.round((100 * guessed) / Math.max(1, n))}% at a default height)` : ''}, ${(model.roads || []).length.toLocaleString('en')} ways`),
    h('div', { class: 'faint' }, model.credit || '© OpenStreetMap contributors'));
}
