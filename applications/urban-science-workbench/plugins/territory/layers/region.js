// A region, cell by cell: Space2Stats' H3 level 6 cells (about 36 km²) wherever the camera looks, asked a tile at a time
// and kept, raised as columns. What a cell holds over years has a scale of years above the clock. One variable is shown at a
// time (the manifest's rows are exclusive): the columns stand in the same cells.

import { card, compact, nearest, tilesInView, tileAt, tone } from '../shared.js';

const year = (prefix, y) => `${prefix}${y}`;
const range = (a, b, step = 1) => Array.from({ length: Math.floor((b - a) / step) + 1 }, (_, i) => a + i * step);

export const VARIABLES = {
  people: {
    title: 'People, cell by cell', years: range(2015, 2030), start: 2025, field: (y) => year('sum_pop_', y), map: 'plum', scale: 'sqrt', top: 95000, source: 'WorldPop and the UN’s projections, through Space2Stats (World Bank)',
    label: 'people', note: 'Columns by the square root of people on a cell of about 36 km².',
  },
  lights: {
    title: 'Night lights, cell by cell', years: range(2012, 2024), start: 2024, field: (y) => year('sum_viirs_ntl_', y), map: 'gold', scale: 'sqrt', top: 95000, source: 'VIIRS annual composites, through Space2Stats (World Bank)',
    label: 'the sum of radiance', note: 'Columns by the square root of the summed radiance on a cell.',
  },
  built: {
    title: 'Built-up land, cell by cell', years: range(1975, 2030, 5), start: 2020, field: (y) => year('sum_built_area_m_', y), map: 'orchid', scale: 'linear', top: 70000, source: 'GHS-BUILT-S R2023A, through Space2Stats (World Bank); 2025 and 2030 are projections',
    label: 'm² built', note: 'Columns by the surface built on a cell.',
  },
  urban: {
    title: 'Urban centres and suburbs', field: () => null, value: (t, k) => {
      const f = t.fields, all = f.ghs_total_pop[k];
      if (!(all > 0)) return null;
      return (f.ghs_30_pop[k] + f.ghs_23_pop[k] + f.ghs_22_pop[k] + 0.6 * f.ghs_21_pop[k]) / all;
    }, height: (t, k) => t.fields.ghs_total_pop[k], map: 'orchid', scale: 'share', top: 90000, source: 'GHS-SMOD R2023A, through Space2Stats (World Bank)',
    label: 'of its people in urban land', note: 'Colour: the share of a cell’s people in urban centres, dense clusters and, at 0.6, suburbs. Height: its people.', percent: true,
  },
  flood: {
    title: 'People in flood zones', field: () => 'pop_flood', map: 'sea', scale: 'sqrt', top: 95000, source: 'Fathom flood maps with WorldPop, through Space2Stats (World Bank)',
    label: 'people', note: 'People on a cell living where a 1-in-100-year flood reaches over 15 cm.',
  },
  cyclone: {
    title: 'Tropical cyclones', field: () => 'cy_frequency_mean', map: 'lavender', scale: 'linear', top: 60000, source: 'Cyclone frequency, through Space2Stats (World Bank)',
    label: 'mean frequency', note: 'How often tropical cyclones pass, averaged over a cell.',
  },
  landslide: {
    title: 'Landslide susceptibility', field: () => 'landslide_susceptibility_mean_2023', map: 'gold', scale: 'linear', top: 50000, source: 'Global landslide susceptibility, through Space2Stats (World Bank)',
    label: 'susceptibility', note: 'A cell’s mean susceptibility to landslides.',
  },
  drought: {
    title: 'Drought', field: () => 'drought_spei_1_5_rp100_mean', map: 'gold', scale: 'magnitude', top: 50000, source: 'SPEI at 1.5 °C, a 100-year event, through Space2Stats (World Bank)',
    label: 'SPEI', note: 'The depth of a 100-year drought at 1.5 °C of warming, averaged over a cell.',
  },
  fires: {
    title: 'Wildfires', field: () => 'fires_density_mean', map: 'gold', scale: 'sqrt', top: 50000, source: 'Fire density, through Space2Stats (World Bank)',
    label: 'density', note: 'A cell’s mean density of fires.',
  },
};

export async function layer(ctx) {
  const { earth, dom: { h } } = ctx;
  const V = VARIABLES[ctx.spec.variable];
  const id = ctx.layerId();
  const fade = ctx.spec.scale || [18000, 2800000];
  let current = V.years ? V.start : null;
  let loaded = [];
  let set = null;
  let cells = [];   // [{ tile, k }]
  let near = null;
  let dead = false;
  const colour = tone(V.map === 'lavender' ? 'orchid' : V.map, ctx.theme, { from: 0.12 });
  const status = h('div', { class: 'ld-quiet' }, 'asking for the cells in view…');
  const ends = h('div', { class: 'ramp-ends' }, h('span', {}, ''), h('span', {}, ''));

  const valueOf = (c) => {
    if (V.value) return V.value(c.tile, c.k);
    const f = c.tile.fields[V.field(current)];
    const v = f ? f[c.k] : null;
    return v === null || v === undefined ? null : Math.abs(v);
  };
  const heightOf = (c) => (V.height ? c.tile.fields.ghs_total_pop[c.k] : valueOf(c));

  function paint() {
    if (!set) return;
    const values = cells.map(valueOf), heights = cells.map(heightOf);
    const known = values.filter((v) => Number.isFinite(v) && v > 0).sort((a, b) => a - b);
    const hk = heights.filter((v) => Number.isFinite(v) && v > 0).sort((a, b) => a - b);
    const top = known.length ? known[Math.floor(known.length * 0.995)] : 1;
    const htop = hk.length ? hk[Math.floor(hk.length * 0.995)] : 1;
    const f = (v, t) => (V.scale === 'sqrt' ? Math.sqrt(Math.max(0, v) / t) : Math.max(0, v) / t);
    set.prisms.set(heights.map((v) => (Number.isFinite(v) && v > 0 ? V.top * Math.min(1, V.height ? Math.sqrt(v / htop) : f(v, top)) : -1)),
      values.map((v) => (Number.isFinite(v) && v > 0 ? colour(Math.min(1, V.percent ? v : f(v, top))) : null)));
    ends.children[0].textContent = V.percent ? '0' : 'little';
    ends.children[1].textContent = V.percent ? '100 %' : `${compact(top)} ${V.label}`;
  }

  async function refresh() {
    const v = earth.view;
    if (dead || v.range < fade[0] * 0.5 || v.range > fade[1] * 1.8) return;
    const want = tilesInView(v);
    const keys = want.map(([la, lo]) => `${la},${lo}`);
    if (keys.every((k) => loaded.includes(k))) return;
    status.textContent = 'asking for the cells in view…';
    const got = await Promise.allSettled(want.map(([la, lo]) => tileAt(la + 2, lo + 2)));
    if (dead) return;
    const fresh = [];
    got.forEach((g, i) => { if (g.status === 'fulfilled') fresh.push([keys[i], g.value]); });
    const before = new Set(loaded);
    if (fresh.every(([k]) => before.has(k)) && fresh.length === loaded.length) { status.textContent = `${cells.length.toLocaleString('en')} cells of H3 level 6`; return; }
    loaded = fresh.map(([k]) => k);
    cells = [];
    const rings = [], centres = [];
    for (const [, t] of fresh) t.cells.forEach((_, k) => { cells.push({ tile: t, k }); rings.push({ ring: t.rings[k] }); centres.push(t.centres[k]); });
    if (!cells.length) { status.textContent = 'no cells here (the service answers for land)'; return; }
    set = earth.addColumns(id, rings, { fade, shrink: 0.9, grow: !set });
    near = nearest(centres, 1);
    paint();
    status.textContent = `${cells.length.toLocaleString('en')} cells of H3 level 6`;
  }
  const timer = setInterval(() => refresh().catch((e) => { status.textContent = String(e.message || e); }), 900);
  refresh().catch((e) => { status.textContent = String(e.message || e); });
  const years = V.years ? ctx.years({ label: V.title.split(',')[0], years: V.years, value: current, onChange: (y) => { current = y; paint(); } }) : null;

  return {
    legend: card(h, V.title, { colour, note: V.note, extra: [], source: V.source, low: undefined }),
    afterLegend: null,
    onLegend(c) { c.insertBefore(ends, c.querySelector('.ld-quiet')); c.append(status); },
    hoverAt(lat, lon) {
      if (!near) return null;
      const range = earth.view.range;
      if (range < fade[0] * 0.5 || range > fade[1] * 1.7) return null;
      const k = near(lon, lat, 0.06);
      if (k < 0) return null;
      const c = cells[k], v = valueOf(c);
      if (!Number.isFinite(v)) return null;
      return `${V.percent ? `${Math.round(v * 100)} %` : compact(v)} ${V.label}${current ? ` in ${current}` : ''}`;
    },
    dispose() {
      dead = true;
      clearInterval(timer);
      earth.removeColumns(id);
      if (years) years.close();
    },
  };
}
