// A composite of place: your own index over a region, cell by cell. Each factor, from Space2Stats' grid, is scaled to 0–5
// between its 2nd and 98th percentiles in the cells in view, and they are averaged with the weights you set, a factor missing on a
// cell left out of its mean (nomos's composite). The idea is that of GeoE3; the arithmetic is written here, and it is a
// description of a region, by one's own weights, not a measure of anything by itself.

import { card, compact, nearest, tilesInView, tileAt, tone } from '../shared.js';

const FACTORS = [
  { key: 'pop', name: 'People', field: (f, k) => f.sum_pop_2025[k], weight: 1 },
  { key: 'lights', name: 'Night light', field: (f, k) => f.sum_viirs_ntl_2024[k], weight: 1 },
  { key: 'built', name: 'Built-up land', field: (f, k) => f.sum_built_area_m_2020[k], weight: 1 },
  { key: 'urban', name: 'Urbanisation', field: (f, k) => (f.ghs_total_pop[k] > 0 ? (f.ghs_30_pop[k] + f.ghs_23_pop[k] + f.ghs_22_pop[k]) / f.ghs_total_pop[k] : null), weight: 0 },
  { key: 'flood', name: 'Flood exposure', field: (f, k) => f.pop_flood_pct[k], weight: 0, invert: false },
  { key: 'landslide', name: 'Landslide susceptibility', field: (f, k) => f.landslide_susceptibility_mean_2023[k], weight: 0 },
  { key: 'cyclone', name: 'Cyclones', field: (f, k) => f.cy_frequency_mean[k], weight: 0 },
  { key: 'fires', name: 'Wildfires', field: (f, k) => f.fires_density_mean[k], weight: 0 },
];

export async function render(root, ctx) {
  const { earth, dom: { h } } = ctx;
  const id = 'composite';
  let cells = [], near = null, set = null, token = 0;
  const status = h('p', { class: 'muted' }, 'asking for the cells in view…');
  const legendBox = h('div', {});
  const list = h('div', { class: 'factors' });
  const colour = tone('plum', ctx.theme, { from: 0.14 });
  const out = h('div', {});
  root.append(h('h1', {}, 'A composite of place'),
    h('p', { class: 'muted' }, 'Weigh the factors; the columns over the region follow as you move them. Each factor is scaled to 0–5 within what is in view, so that the index compares the cells here with one another, and nothing else.'),
    list, status, legendBox, out);
  ctx.onDispose(() => earth.removeColumns(id));

  const state = FACTORS.map((F) => ({ ...F }));
  for (const F of state) {
    const value = h('span', { class: 'mono' }, F.weight.toFixed(1));
    const input = h('input', { type: 'range', min: 0, max: 3, step: 0.1, value: F.weight, 'aria-label': `the weight of ${F.name}`, oninput: () => { F.weight = +input.value; value.textContent = F.weight.toFixed(1); schedule(); } });
    list.append(h('label', { class: 'factor' }, h('span', {}, F.name), input, value));
  }
  let timer = null;
  const schedule = () => { clearTimeout(timer); timer = setTimeout(paint, 120); };

  async function load() {
    const want = tilesInView(earth.view, 9);
    const got = await Promise.allSettled(want.map(([la, lo]) => tileAt(la + 2, lo + 2)));
    const tiles = got.filter((g) => g.status === 'fulfilled').map((g) => g.value);
    cells = [];
    const rings = [], centres = [];
    for (const t of tiles) t.cells.forEach((_, k) => { cells.push({ t, k }); rings.push({ ring: t.rings[k] }); centres.push(t.centres[k]); });
    if (!cells.length) { status.textContent = 'No cells here: the service answers for land. Move the view over a country and open this again.'; return; }
    set = earth.addColumns(id, rings, { fade: [15000, 3.2e6], shrink: 0.9 });
    near = nearest(centres, 1);
    paint();
  }

  async function paint() {
    if (!set) return;
    const me = ++token;
    const factors = state.filter((F) => F.weight > 0).map((F) => ({ values: cells.map((c) => { const v = F.field(c.t.fields, c.k); return Number.isFinite(v) ? v : null; }), weight: F.weight, invert: !!F.invert }));
    if (!factors.length) { set.prisms.set(cells.map(() => -1), cells.map(() => null)); status.textContent = 'Give some factor a weight.'; return; }
    const { index } = await ctx.ask('composite', { factors });
    if (me !== token) return;
    set.prisms.set(index.map((v) => (v === null ? -1 : 1.5e3 + (v / 5) * 4.2e4)), index.map((v) => (v === null ? null : colour(v / 5))));
    const known = index.filter((v) => v !== null);
    const mean = known.reduce((a, b) => a + b, 0) / Math.max(1, known.length);
    status.textContent = `${known.length.toLocaleString('en')} cells of about 36 km² in view; the mean index is ${mean.toFixed(2)} of 5.`;
    if (!legendBox.firstChild) legendBox.append(card(h, 'The index', { colour, low: '0', high: '5' }));
    set.index = index;
  }
  // Move the view and it asks again, once the camera has settled somewhere new.
  let at = null;
  const follow = setInterval(() => {
    const v = earth.view, here = [Math.floor((v.lat * 57.3) / 4), Math.floor((v.lon * 57.3) / 4), Math.round(Math.log2(v.range))];
    if (at && here.join() === at.join()) return;
    at = here;
    load().catch((e) => { status.textContent = String(e.message || e); });
  }, 1200);
  ctx.onDispose(() => clearInterval(follow));
  ctx.onDispose(() => clearTimeout(timer));
  out.append(ctx.finding({
    title: 'What the index is', kind: 'model',
    statement: 'A description of the region by your own weights, after the composite indicators of GeoE3 and of the OECD and JRC handbook. It is not a measure: other weights give another order.',
    method: 'Each factor to 0–5 between its 2nd and 98th percentiles in view; the weighted mean of those present on a cell (nomos measurement::composite).',
    limits: 'The percentiles move with the view, so that an index read here cannot be compared with one read elsewhere.',
  }));
  void compact;
}
