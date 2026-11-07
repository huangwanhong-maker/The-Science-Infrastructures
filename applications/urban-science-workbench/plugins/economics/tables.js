// The tables a tool can run on from where one stands: the place's own, the cells about the point looked at, the world's economies,
// what is shown in the studio, or a table pasted. A table is { id, title, source, columns: [{name, kind, label}], data: {name: [...]},
// rows, places }, its rows being places where `places` says how:
//
//   { kind: 'cells', column, time?, time_kind?, rings: {id: ring}, names: {id: name} }   a polygon for each id
//   { kind: 'points', column, time?, time_kind?, at: {id: [lon, lat]}, names }           a point for each id
//   { kind: 'economies', column, time, time_kind, at, names }                            an economy: its cells of the world's mosaic
//   { kind: 'row', lon, lat, time?, time_kind? }                                         each row its own point
//
// A table with places gains the columns place_lon and place_lat, for the tools that ask where its rows are.

import { spacetime as ST } from '/visualization/js/gis/index.js';
import { parsePaste } from '/visualization/js/toolbox/form.js';
import { tileAt } from '../territory/shared.js';
import { VARIABLES } from '../territory/layers/region.js';
import { centreOf } from './answers.js';

const DEG = 180 / Math.PI;
const files = new Map();

function product(id) {
  if (!files.has(id)) {
    files.set(id, fetch(`/data/products/tables/${id}.json`).then(async (r) => {
      if (!r.ok) throw new Error(`The data library has no table ${id} yet: the records lamp runs the recipe that makes it.`);
      return r.json();
    }).catch((e) => { files.delete(id); throw e; }));
  }
  return files.get(id);
}


function positioned(t) {
  const P = t.places;
  if (!P || P.kind === 'row' || !t.data[P.column]) return t;
  const where = (id) => (P.kind === 'cells' ? (P.rings[id] ? centreOf(P.rings[id]) : null) : P.at?.[id] ?? null);
  const seen = new Map();
  const lon = [], lat = [];
  for (const raw of t.data[P.column]) {
    const id = String(raw);
    if (!seen.has(id)) seen.set(id, where(id));
    const c = seen.get(id);
    lon.push(c ? c[0] : null);
    lat.push(c ? c[1] : null);
  }
  return {
    ...t,
    data: { ...t.data, place_lon: lon, place_lat: lat },
    columns: [...t.columns, { name: 'place_lon', kind: 'number', label: 'where this row is: longitude' }, { name: 'place_lat', kind: 'number', label: 'where this row is: latitude' }],
  };
}

// How tall an answer stands over each table's places, in metres: under the towers for a city's blocks, above them for its counters.
const TOPS = { 'melbourne-blocks': { top: 46 }, 'melbourne-counters-day': { top: 240, base: 0.3 }, 'melbourne-counters-month': { top: 240, base: 0.3 }, 'freetown-wards': { top: 320 } };

async function library(id) {
  const t = positioned(await product(id));
  return t.places ? { ...t, places: { ...t.places, ...(TOPS[id] || {}) } } : t;
}

// The economies, whose places are the world's own: each economy's point, its name by the data library's policy.
async function economies() {
  const t = await product('economies');
  const at = {}, names = {};
  const { iso3, name, longitude, latitude } = t.data;
  iso3.forEach((id, k) => {
    names[id] = name[k];
    if (Number.isFinite(longitude[k]) && Number.isFinite(latitude[k])) at[id] = [longitude[k], latitude[k]];
  });
  return { ...t, places: { kind: 'economies', column: 'iso3', time: 'year', time_kind: 'year', at, names, top: 3.2e5 } };
}

// The cells of the world's mosaic, kept: every economy's cells, those in dispute having none (the territory's engine).
let mosaic = null;
export function worldMosaic() {
  if (!mosaic) {
    mosaic = fetch('/api/territory/earth_people', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ res: 2 }) }).then(async (r) => {
      const w = await r.json();
      if (!r.ok) throw new Error(w.message || 'the territory\'s engine did not answer');
      const by = new Map();
      w.code.forEach((c, k) => { if (c) { if (!by.has(c)) by.set(c, []); by.get(c).push(w.rings[k]); } });
      return by;
    }).catch((e) => { mosaic = null; throw e; });
  }
  return mosaic;
}

// Space2Stats' cells of the tile under the point looked at.
async function cellsHere(ctx, through) {
  const v = ctx.earth.view, tile = await tileAt(v.lat * DEG, v.lon * DEG), f = tile.fields;
  const rings = {}, names = {};
  const rows = [];
  const deg = (v, pos, neg) => `${Math.abs(v).toFixed(2)}°${v >= 0 ? pos : neg}`;
  tile.cells.forEach((cell, k) => {
    rings[cell] = tile.rings[k];
    const [x, y] = tile.centres[k];
    names[cell] = `The cell at ${deg(y, 'N', 'S')} ${deg(x, 'E', 'W')}`;
    if (through) {
      for (let y = 2015; y <= 2024; y++) rows.push({ cell, year: y, people: f[`sum_pop_${y}`][k], lights: f[`sum_viirs_ntl_${y}`][k] });
    } else {
      rows.push({
        cell, people: f.sum_pop_2024[k], lights: f.sum_viirs_ntl_2024[k], built_m2: f.sum_built_area_m_2020[k], urban_share: VARIABLES.urban.value(tile, k),
        people_in_flood_zones: f.pop_flood[k], fires: f.fires_density_mean[k], cyclones: f.cy_frequency_mean[k], landslides: f.landslide_susceptibility_mean_2023[k], drought: f.drought_spei_1_5_rp100_mean[k],
      });
    }
  });
  const spec = through ? [['cell', 'text', 'an H3 cell of about 36 km²'], ['year', 'number', ''], ['people', 'number', 'people on the cell, that year (WorldPop)'], ['lights', 'number', 'the sum of radiance, that year (VIIRS)']]
    : [['cell', 'text', 'an H3 cell of about 36 km²'], ['people', 'number', 'people on the cell, 2024 (WorldPop)'], ['lights', 'number', 'the sum of radiance, 2024 (VIIRS)'], ['built_m2', 'number', 'built-up surface, m², 2020 (GHS-BUILT-S)'],
      ['urban_share', 'number', 'the share of its people in urban land (GHS-SMOD)'], ['people_in_flood_zones', 'number', 'people where a 1-in-100-year flood reaches over 15 cm (Fathom)'], ['fires', 'number', 'wildfires, density'],
      ['cyclones', 'number', 'tropical cyclones, frequency'], ['landslides', 'number', 'landslide susceptibility'], ['drought', 'number', 'drought, at a 1-in-100-year return']];
  const t = {
    id: through ? 'cells-years' : 'cells-here',
    title: through ? 'The cells about here, year by year' : 'The cells about here',
    source: 'Space2Stats (World Bank), from WorldPop, VIIRS, GHSL, Fathom and others',
    rows: rows.length,
    columns: spec.map(([name, kind, label]) => ({ name, kind, label })),
    data: Object.fromEntries(spec.map(([n]) => [n, rows.map((r) => r[n] ?? null)])),
    places: { kind: 'cells', column: 'cell', ...(through ? { time: 'year', time_kind: 'year' } : {}), rings, names },
  };
  return positioned(t);
}

// A record shown in the studio, as a table: each event a row; each cell and step of time a row; each sample of a track a row.
function tableOfRecord(e) {
  const r = e.rec, variables = r.variables || Object.keys(r.values || {});
  const base = { id: `shown:${e.id}`, title: r.title, source: 'Shown in the studio' };
  const column = (name, kind = 'number', label = '') => ({ name, kind, label });
  if (r.kind === 'events') {
    return {
      ...base, rows: r.count,
      columns: [column('time_ms', 'number', 'milliseconds since 1970, UTC'), column('lat'), column('lon'), ...variables.map((v) => column(v, 'number', r.units?.[v] || ''))],
      data: { time_ms: [...r.time], lat: [...r.lat], lon: [...r.lon], ...Object.fromEntries(variables.map((v) => [v, Array.from(r.values[v], (x) => (Number.isFinite(x) ? x : null))])) },
      places: { kind: 'row', lon: 'lon', lat: 'lat', time: 'time_ms', time_kind: 'ms' },
    };
  }
  if (r.kind === 'cells') {
    const n = r.count, steps = r.time.length;
    const rings = {}, names = {};
    r.cells.forEach((c) => { rings[c.id] = c.ring; names[c.id] = c.name || c.id; });
    const cell = [], time = [];
    for (let s = 0; s < steps; s++) for (let k = 0; k < n; k++) { cell.push(String(r.cells[k].id)); time.push(r.time[s]); }
    return positioned({
      ...base, rows: cell.length,
      columns: [column('cell', 'text'), column('time_ms', 'number', 'milliseconds since 1970, UTC'), ...variables.map((v) => column(v, 'number', r.units?.[v] || ''))],
      data: { cell, time_ms: time, ...Object.fromEntries(variables.map((v) => [v, Array.from(r.values[v], (x) => (Number.isFinite(x) ? x : null))])) },
      places: { kind: 'cells', column: 'cell', time: 'time_ms', time_kind: 'ms', rings, names },
    });
  }
  return null;
}

// A table pasted: its rows are places if it has columns that say where.
export function pasted(text) {
  const { columns, data } = parsePaste(text);
  const rows = columns.length ? data[columns[0].name].length : 0;
  const t = { id: '__paste', title: 'A table pasted', source: 'Pasted here', columns, data, rows };
  const g = ST.guessColumns({ columns: columns.map((c) => c.name), rows: [] });
  if (g.lat && g.lon) t.places = { kind: 'row', lon: g.lon, lat: g.lat };
  return t;
}

// What can be chosen from where one stands: grouped, each with how to read it.
export function tablesHere(ctx) {
  const area = ctx.area?.id;
  const out = [];
  if (area === 'melbourne') {
    out.push({ id: 'melbourne-blocks', group: 'Melbourne', title: 'Blocks by census year: floor space, businesses, dwellings', load: () => library('melbourne-blocks') });
    out.push({ id: 'melbourne-counters-day', group: 'Melbourne', title: 'Counters by day: walkers and the weather', load: () => library('melbourne-counters-day') });
    out.push({ id: 'melbourne-counters-month', group: 'Melbourne', title: 'Counters by month since 2009: walkers', load: () => library('melbourne-counters-month') });
  }
  if (area === 'freetown') out.push({ id: 'freetown-wards', group: 'Freetown', title: 'Wards: the annual loss from floods', load: () => library('freetown-wards') });
  out.push({ id: 'cells-here', group: 'About the point looked at', title: 'The cells about here: people, lights, built land, hazards', load: () => cellsHere(ctx, false) });
  out.push({ id: 'cells-years', group: 'About the point looked at', title: 'The cells about here, year by year: people and lights', load: () => cellsHere(ctx, true) });
  out.push({ id: 'economies', group: 'The world', title: 'Economies by year: trade and the World Bank\'s indicators', load: economies });
  out.push({ id: 'currencies', group: 'The world', title: 'Currencies by month, long', load: () => product('currencies') });
  out.push({ id: 'currency-rates', group: 'The world', title: 'Currencies by month, wide', load: () => product('currency-rates') });
  for (const e of ctx.st.entries()) {
    if (e.style?.toolbox) continue;
    const t = tableOfRecord(e);
    if (t) out.push({ id: t.id, group: 'Shown in the studio', title: `${t.title} (${t.rows.toLocaleString('en')} rows)`, load: async () => t });
  }
  return out;
}

// A table the questions may name that this place does not offer: the world's and the cells' are everywhere, the rest by their city.
export function findTable(ctx, id) {
  return tablesHere(ctx).find((t) => t.id === id) || null;
}
