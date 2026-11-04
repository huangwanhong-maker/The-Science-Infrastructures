// The map: flows between places in a year on an Equal Earth map (the geography domain's projection),
// from the library's world trade or from one's own data.
//
// The engine splits the year's flows into reciprocal exchange (what two places send each other both
// ways, a round trip) and net flows, one for each pair; and the net flows, on a complex whose faces are
// the spherical Delaunay triangles of the places, into the driven part, local circulation round the
// triangles, and global circulation round the holes the map leaves. Each link conducts in proportion
// to its gross flow (Ohm's law), so the driven net flow between two places is their gross flow times
// their difference in attraction α: α is a share of gross flow, without unit, the potential from which
// value is drawn to where it is drawn to.
//
// Places are coloured by α, by the flows' other measures, by the World Bank's indicators (from the
// toolbox's table of economies), or by the values a plug-in returned; places with a country of their
// own colour it, the others are circles. Arcs are the largest net flows, their driven or circulating
// parts, or a plug-in's arcs; small circular arrows on the triangles are the eddies of local
// circulation. The strip beneath gives every year's shares and plays through them; the map zooms and
// pans, finds a place, and gives its numbers away as files and itself as a figure composed for
// publication. Its instruments name the projection and the place under the pointer; its side panels
// fold away; its address keeps what it shows, so that a view can be linked and drawn again.

import * as api from './api.js';
import { el, svg, clear, put, fill, money, quantity, percent, signed, shareBar, download, csv } from './dom.js';
import { lineChart } from './charts.js';
import { diverging, sequential, rampCss, PART, onTheme } from './palette.js';
import { shared, on, dataset } from './state.js';
import { datasetChooser } from './data.js';
import { pluginsSection } from './plugins.js';
import { equalEarth, extent, outline, inverse } from '/domains/geography/js/projections.js';

const COUNTRIES = '/data/natural-earth-50m/ne_50m_admin_0_countries.geojson';
const SCALE = 300;   // plane units to SVG units
const [XM, YM] = extent(equalEarth);
const W = XM * SCALE, H = YM * SCALE;
const DEG = Math.PI / 180;
// Trading enough to rank: over 20 bn $ of exports and imports, for the library's trade.
const RANKED = 2e4;

const state = {
  geo: null, answer: null, selected: null, wanted: undefined, cache: new Map(), year: null, years: [],
  timeline: new Map(), view: { k: 1, x: 0, y: 0 }, opts: { faces: true, eddies: true, labels: true, flowing: true, shown: 150 },
  drawn: [], built: false, covered: new Set(), indicators: null, playing: null,
};

// The raw projection, not `forward`, which wraps +180° to −180°: Natural Earth splits its polygons at
// the antimeridian, and a ring's vertices on +180° must stay on the map's right edge.
const project = (lon, lat) => {
  const [x, y] = equalEarth.raw(Math.max(-180, Math.min(180, lon)) * DEG, lat * DEG);
  return [x * SCALE, -y * SCALE];
};

function economyCode(p) {
  for (const k of ['ISO_A3', 'ISO_A3_EH', 'ADM0_A3']) if (p[k] && p[k] !== '-99' && p[k].length === 3) return p[k];
  return null;
}

function ringPath(ring) {
  let d = '';
  ring.forEach(([lon, lat], i) => {
    const [x, y] = project(lon, lat);
    d += `${i ? 'L' : 'M'}${x.toFixed(1)},${y.toFixed(1)}`;
  });
  return `${d}Z`;
}

function featurePath(geometry) {
  const polygons = geometry.type === 'Polygon' ? [geometry.coordinates] : geometry.type === 'MultiPolygon' ? geometry.coordinates : [];
  return polygons.map((poly) => poly.map(ringPath).join('')).join('');
}

// The great-circle arc between two places, projected, split where it crosses the map's edge, lifted a
// little so that arcs between the same places part.
function arcPath(a, b, bulge = 0.12) {
  const toVec = (lat, lon) => [Math.cos(lat * DEG) * Math.cos(lon * DEG), Math.cos(lat * DEG) * Math.sin(lon * DEG), Math.sin(lat * DEG)];
  const u = toVec(a.latitude, a.longitude), v = toVec(b.latitude, b.longitude);
  const omega = Math.acos(Math.max(-1, Math.min(1, u[0] * v[0] + u[1] * v[1] + u[2] * v[2])));
  if (omega < 1e-6) return '';
  const n = Math.max(8, Math.ceil(omega * 40));
  let d = '', prev = null;
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    const s0 = Math.sin((1 - t) * omega) / Math.sin(omega), s1 = Math.sin(t * omega) / Math.sin(omega);
    const p = [s0 * u[0] + s1 * v[0], s0 * u[1] + s1 * v[1], s0 * u[2] + s1 * v[2]];
    const lat = Math.asin(Math.max(-1, Math.min(1, p[2]))) / DEG, lon = Math.atan2(p[1], p[0]) / DEG;
    let [x, y] = project(lon, lat);
    y -= bulge * SCALE * Math.sin(Math.PI * t) * omega * 0.35;
    const jump = prev && Math.abs(x - prev[0]) > W;
    d += `${!prev || jump ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`;
    prev = [x, y];
  }
  return d;
}

// ---- layers: what colours the places, which arcs ----------------------------------------------------------

const layer = (key) => shared.layers.find((l) => l.key === key) || null;
// Amounts of the flows: money for the library's trade (millions of US dollars), plain numbers in one's own unit.
const amount = (x) => (shared.flows === 'library:trade' ? money(x) : quantity(x));
const fmt = (x) => (!Number.isFinite(x) ? '—' : Math.abs(x) >= 1000 ? quantity(x) : Math.abs(x) >= 1 ? x.toFixed(2) : x.toPrecision(2));

// The q-quantile of sorted values, linear between order statistics.
const quantile = (s, q) => { const h = (s.length - 1) * q, i = Math.floor(h); return i + 1 < s.length ? s[i] + (h - i) * (s[i + 1] - s[i]) : s[i]; };

// Values by place id coloured on a scale: diverging about zero, or sequential from the least to the
// most, on a log scale when the values are positive and span more than two decades. With twenty
// places or more the scale runs between the 2nd and 98th percentiles (of |v| for a diverging one),
// so that one outlying value does not wash the rest into one end; the values beyond are clamped,
// and the legend's ends say so.
function scaled(title, label, values, sequentialScale, show = fmt) {
  const xs = Object.values(values).filter(Number.isFinite).sort((a, b) => a - b);
  if (!xs.length) return null;
  const robust = xs.length >= 20;
  const lo = robust ? quantile(xs, 0.02) : xs[0], hi = robust ? quantile(xs, 0.98) : xs.at(-1);
  const mags = xs.map(Math.abs).sort((a, b) => a - b);
  const big = Math.max(1e-12, robust ? quantile(mags, 0.98) : mags.at(-1));
  const log = sequentialScale && lo > 0 && hi / lo > 100;
  const clamp = (u) => Math.max(0, Math.min(1, u));
  const t = (v) => clamp(log ? Math.log(v / lo) / Math.log(hi / lo) : hi > lo ? (v - lo) / (hi - lo) : 0.5);
  const below = xs[0] < lo, above = xs.at(-1) > hi, beyond = mags.at(-1) > big;
  return {
    title: `${title}${log ? ', on a log scale' : ''}`, kind: sequentialScale ? 'sequential' : 'diverging',
    ends: sequentialScale ? [`${below ? '≤ ' : ''}${show(lo)}`, `${above ? '≥ ' : ''}${show(hi)}`] : [`${beyond ? '≤ ' : ''}−${show(big)}`, `${beyond ? '≥ ' : ''}+${show(big)}`],
    of: (e) => { const v = values[e.iso3]; return Number.isFinite(v) ? (sequentialScale ? sequential(t(v)) : diverging(Math.max(-1, Math.min(1, v / big)))) : null; },
    value: (e) => values[e.iso3], label, show,
  };
}

// The flows' own measures of each place, besides α.
const MEASURES = {
  balance: ['Net inflow: imports − exports', 'net inflow', (e) => e.balance, false],
  exports: ['Exports', 'exports', (e) => e.exports, true],
  imports: ['Imports', 'imports', (e) => e.imports, true],
  circulating: ['Circulating, of the parts on its links', 'circulating', (e) => e.circulating, true],
};

function colouring(A) {
  const key = shared.colour;
  const byPlace = (f) => Object.fromEntries(A.economies.map((e) => [e.iso3, f(e)]));
  if (MEASURES[key]) {
    const [title, label, f, seq] = MEASURES[key];
    const C = scaled(title, label, byPlace(f), seq, key === 'circulating' ? (v) => percent(v, 0) : amount);
    if (C) return C;
  }
  if (key.startsWith('wb:') && state.indicators) {
    const name = key.slice(3), column = state.indicators.columns.find((c) => c.name === name);
    const values = {}, years = new Set();
    for (const e of A.economies) {
      const found = indicator(name, e.iso3, state.year);
      if (found) { values[e.iso3] = found.value; years.add(found.year); }
    }
    const span = years.size ? [Math.min(...years), Math.max(...years)] : [];
    const when = span.length ? (span[0] === span[1] ? `${span[0]}` : `${span[0]}–${span[1]}, the latest by ${state.year}`) : '';
    const C = scaled(`${column?.label || name} (World Bank${when ? `, ${when}` : ''})`, column?.label || name, values, true);
    if (C) return C;
  }
  const L = key !== 'attraction' ? layer(key) : null;
  if (L && L.result.kind === 'values') {
    const C = scaled(`${L.title}${L.result.unit ? ` (${L.result.unit})` : ''}`, L.title, L.result.values || {}, L.result.scale === 'sequential');
    if (C) return C;
  }
  const maxA = Math.max(1e-9, ...A.economies.map((e) => Math.abs(e.attraction)));
  return { title: 'Attraction α, a share of gross flow', kind: 'diverging', ends: [`sends, −${maxA.toFixed(2)}`, `draws, +${maxA.toFixed(2)}`],
    of: (e) => diverging(e.attraction / maxA), value: (e) => e.attraction, label: 'α', show: fmt };
}

// ---- the World Bank's indicators, from the toolbox's table of economies by year ---------------------------------

// The table indexed by economy and year; the columns that are not the flows' own.
const FLOW_COLUMNS = new Set(['year', 'exports', 'imports', 'balance', 'trade']);
async function loadIndicators() {
  if (state.indicators !== null) return state.indicators;
  state.indicators = false;   // asked for once
  try {
    const t = await api.table('economies');
    const index = new Map();
    t.data.iso3.forEach((code, i) => index.set(`${code}|${t.data.year[i]}`, i));
    state.indicators = { ...t, index, columns: t.columns.filter((c) => c.kind === 'number' && !FLOW_COLUMNS.has(c.name)) };
  } catch {
    // Without the table the map colours by the flows alone.
  }
  return state.indicators;
}

// An indicator's value for an economy in a year, or in the latest of the five years before it.
function indicator(name, code, year) {
  const T = state.indicators;
  if (!T || !T.data[name]) return null;
  for (let y = year; y >= year - 5; y--) {
    const i = T.index.get(`${code}|${y}`);
    const v = i === undefined ? null : T.data[name][i];
    if (Number.isFinite(v)) return { value: v, year: y };
  }
  return null;
}

// An economy's years in the table: its trade and, when the map is coloured by one, that indicator.
function placeYears(code) {
  const T = state.indicators;
  if (!T) return null;
  const rows = [];
  T.data.iso3.forEach((c, i) => { if (c === code) rows.push(i); });
  rows.sort((a, b) => T.data.year[a] - T.data.year[b]);
  return { years: rows.map((i) => T.data.year[i]), column: (name) => rows.map((i) => T.data[name]?.[i] ?? null) };
}

// The arcs drawn: [{from, to, value, colour, parts?, title?}], from and to indices into the economies.
function arcsOf(A) {
  if (shared.arcs === 'none') return [];
  const L = ['all', 'driven', 'circulating'].includes(shared.arcs) ? null : layer(shared.arcs);
  const keep = (a) => state.selected === null || a.from === state.selected || a.to === state.selected;
  if (L && L.result.kind === 'arcs') {
    const index = new Map(A.economies.map((e, i) => [e.iso3, i]));
    return (L.result.arcs || []).map((a) => ({ from: index.get(a.from), to: index.get(a.to), value: Number(a.value), colour: PART.plugin, title: L.title }))
      .filter((a) => a.from !== undefined && a.to !== undefined && Number.isFinite(a.value) && keep(a));
  }
  const part = shared.arcs;
  const weight = (l) => (part === 'driven' ? Math.abs(l[3]) : part === 'circulating' ? Math.abs(l[4] + l[5]) : l[2]);
  return A.legs.filter((l) => state.selected === null || l[0] === state.selected || l[1] === state.selected)
    .sort((a, b) => weight(b) - weight(a)).slice(0, state.opts.shown).map((l) => {
      const value = part === 'driven' ? l[3] : part === 'circulating' ? l[4] + l[5] : l[2];
      const [from, to] = value >= 0 ? [l[0], l[1]] : [l[1], l[0]];
      // An arc wears the colour of its largest part: driven, round the triangles, or round the holes.
      const [d, lo, g] = [Math.abs(l[3]), Math.abs(l[4]), Math.abs(l[5])];
      const colour = part === 'driven' ? PART.driven : part === 'circulating' ? (lo >= g ? PART.local : PART.global)
        : d >= lo && d >= g ? PART.driven : lo >= g ? PART.local : PART.global;
      return { from, to, value: Math.abs(value), colour, parts: l };
    });
}

// ---- drawing ------------------------------------------------------------------------------------------------

function legend(C) {
  const box = document.getElementById('map-legend');
  const key = (colour, text, dash) => el('div', { class: 'key' }, svg('svg', { width: 28, height: 10 },
    svg('path', { d: 'M1 5H27', stroke: colour, 'stroke-width': 2.4, 'stroke-linecap': 'round', 'stroke-dasharray': dash || null })), text);
  const L = layer(shared.arcs);
  const src = sourceLine();
  fill(box,
    el('div', { class: 'title', text: C.title }),
    el('div', { class: 'ramp', style: { background: rampCss(C.kind) } }),
    el('div', { class: 'ends' }, el('span', { text: C.ends[0] }), el('span', { text: C.ends[1] })),
    L ? key(PART.plugin, `arcs: ${L.title}`) : shared.arcs === 'none' ? null : [
      key(PART.driven, 'mostly driven: from senders to drawers'),
      key(PART.local, 'mostly circulating round the triangles of neighbours'),
      key(PART.global, 'mostly circulating round holes, through distant links'),
      state.opts.eddies ? key(PART.local, 'eddies round the triangles', '2 3') : null],
    el('div', { class: 'source', text: src }),
  );
}

// What the map shows, in a line: the flows, their year and unit, and the projection.
function sourceLine() {
  const A = state.answer;
  const d = dataset(shared.flows);
  const what = shared.flows === 'library:trade' ? 'Trade in goods: OECD BIMTS' : d?.title || 'Flows of one’s own';
  const unit = shared.flows === 'library:trade' ? 'millions of US$' : A?.unit || '';
  return `${what}${A?.year ? `, ${A.year}` : ''}${unit ? ` · ${unit}` : ''} · Equal Earth`;
}

function buildBase() {
  const host = document.getElementById('map');
  host.setAttribute('viewBox', `${-W - 8} ${-H - 8} ${2 * W + 16} ${2 * H + 16}`);
  clear(host);
  const world = svg('g', { id: 'world' });
  const sphere = outline(equalEarth, 160).map(([x, y], i) => `${i ? 'L' : 'M'}${(x * SCALE).toFixed(1)},${(-y * SCALE).toFixed(1)}`).join('') + 'Z';
  put(world, svg('defs', {}, svg('clipPath', { id: 'sphere-clip' }, svg('path', { d: sphere }))));
  put(world, svg('path', { class: 'sphere', d: sphere }));
  // A graticule every 30°.
  let g = '';
  for (let lon = -180; lon <= 180; lon += 30) {
    for (let lat = -90, i = 0; lat <= 90; lat += 3, i++) { const [x, y] = project(lon, lat); g += `${i ? 'L' : 'M'}${x.toFixed(1)},${y.toFixed(1)}`; }
  }
  for (let lat = -60; lat <= 60; lat += 30) {
    for (let lon = -180, i = 0; lon <= 180; lon += 3, i++) { const [x, y] = project(lon, lat); g += `${i ? 'L' : 'M'}${x.toFixed(1)},${y.toFixed(1)}`; }
  }
  put(world, svg('path', { class: 'graticule', d: g, 'clip-path': 'url(#sphere-clip)' }));
  const countries = svg('g', { id: 'countries' });
  for (const f of state.geo.features) {
    put(countries, svg('path', { class: 'country', d: featurePath(f.geometry), 'data-code': economyCode(f.properties) || '', 'data-name': f.properties.NAME || '' }));
  }
  put(world, countries, svg('g', { id: 'over' }));
  put(host, world);
  applyView();
  state.built = true;
}

function draw() {
  const A = state.answer;
  if (!state.geo || !A) return;
  if (!state.built) buildBase();
  const C = colouring(A);
  legend(C);
  const byCode = new Map(A.economies.map((e, i) => [e.iso3, i]));
  state.covered = new Set();
  for (const path of document.querySelectorAll('#countries .country')) {
    const i = byCode.get(path.dataset.code);
    const e = i === undefined ? null : A.economies[i];
    const colour = e ? C.of(e) : null;
    if (e) state.covered.add(i);
    path.setAttribute('fill', colour || 'var(--land)');
    path.dataset.i = e ? String(i) : '';
    path.classList.toggle('selected', e !== null && state.selected === i);
  }
  drawOver(C);
  keepAddress();
}

// The address keeps what the map shows (ADR 0068): a link draws the same view again.
function keepAddress() {
  if (document.querySelector('.view[data-view="map"]').hidden) return;
  const A = state.answer;
  const q = [];
  if (shared.flows !== 'library:trade') q.push(['dataset', shared.flows]);
  if (state.year !== null && state.year !== undefined) q.push(['year', state.year]);
  if (shared.colour !== 'attraction' && (MEASURES[shared.colour] || shared.colour.startsWith('wb:'))) q.push(['colour', shared.colour]);
  if (['driven', 'circulating', 'none'].includes(shared.arcs)) q.push(['arcs', shared.arcs]);
  if (state.opts.shown !== 150) q.push(['shown', state.opts.shown]);
  for (const k of ['faces', 'labels', 'flowing']) if (!state.opts[k]) q.push([k, 0]);
  if (A && state.selected !== null) q.push(['economy', A.economies[state.selected].iso3]);
  // Only what would break the address is encoded, so that a link reads as it is (colour=wb:gini).
  const hash = `#map${q.map(([k, v]) => `&${k}=${String(v).replace(/[&#=%\s]/g, encodeURIComponent)}`).join('')}`;
  if (location.hash !== hash) window.history.replaceState(null, '', hash);
}

// What lies over the countries, drawn again when the view zooms (its sizes are kept on the screen).
function drawOver(C = state.answer ? colouring(state.answer) : null) {
  const A = state.answer;
  const over = document.getElementById('over');
  if (!over || !A) return;
  clear(over);
  const k = state.view.k, sk = Math.sqrt(k);
  // Sizes are given in pixels of the screen: the map's units per pixel, at this zoom.
  const px = (2 * W + 16) / Math.max(1, document.getElementById('map').clientWidth) / k;
  // Places without a country of their own: circles, by their flows.
  const through = A.economies.map((e) => e.exports + e.imports);
  const big = Math.max(1e-12, ...through);
  const places = svg('g', {});
  A.economies.forEach((e, i) => {
    if (state.covered.has(i)) return;
    const [x, y] = project(e.longitude, e.latitude);
    const r = (2.5 + 6.5 * Math.sqrt(through[i] / big)) * px * Math.sqrt(sk);
    put(places, svg('circle', { class: `place${state.selected === i ? ' selected' : ''}`, cx: x.toFixed(2), cy: y.toFixed(2), r: r.toFixed(2),
      fill: C.of(e) || 'var(--land)', 'data-i': i }));
  });
  // The map's triangles and their eddies.
  const faces = svg('g', {});
  if (state.opts.faces || state.opts.eddies) {
    const maxCurl = Math.max(1e-12, ...A.triangles.map((t) => Math.abs(t.curl)));
    for (const t of A.triangles) {
      const pts = t.places.map((i) => project(A.economies[i].longitude, A.economies[i].latitude));
      if (Math.max(...pts.map((p) => p[0])) - Math.min(...pts.map((p) => p[0])) > W) continue;  // across the map's edge
      if (state.opts.faces) put(faces, svg('path', { class: 'face', d: `M${pts.map((p) => p.map((x) => x.toFixed(1)).join(',')).join('L')}Z` }));
      const s = Math.abs(t.curl) / maxCurl;
      if (!state.opts.eddies || s < 0.05) continue;
      const cx = (pts[0][0] + pts[1][0] + pts[2][0]) / 3, cy = (pts[0][1] + pts[1][1] + pts[2][1]) / 3, r = (3 + 9 * Math.sqrt(s)) * px * Math.sqrt(sk);
      const sweep = t.curl > 0 ? 0 : 1;  // anticlockwise on the sphere seen from outside is anticlockwise on the map
      const a0 = 0.3, a1 = 5.3, flip = sweep ? -1 : 1;
      const p0 = [cx + r * Math.cos(a0), cy - r * Math.sin(a0) * flip], p1 = [cx + r * Math.cos(a1), cy - r * Math.sin(a1) * flip];
      put(faces, svg('path', { class: 'eddy', stroke: PART.local, 'stroke-width': (1 + 1.2 * s).toFixed(2), 'stroke-opacity': (0.45 + 0.5 * s).toFixed(2),
        'marker-end': 'url(#eddy-head)', d: `M${p0[0].toFixed(2)},${p0[1].toFixed(2)}A${r.toFixed(2)},${r.toFixed(2)} 0 1 ${sweep} ${p1[0].toFixed(2)},${p1[1].toFixed(2)}` }));
    }
  }
  // Arcs.
  const arcs = arcsOf(A);
  state.drawn = arcs;
  const maxW = Math.max(1e-12, ...arcs.map((a) => a.value));
  const arcLayer = svg('g', {});
  arcs.forEach((a, n) => {
    const w = Math.sqrt(a.value / maxW);
    put(arcLayer, svg('path', { class: `arc${state.opts.flowing ? ' flowing' : ''}`, d: arcPath(A.economies[a.from], A.economies[a.to]), stroke: a.colour,
      'stroke-width': (0.5 + 4 * w).toFixed(2), 'stroke-opacity': (0.35 + 0.5 * w).toFixed(2), 'marker-end': state.opts.flowing ? null : 'url(#head)', 'data-arc': n }));
  });
  // Sketches of plug-ins: polylines, markers and labels in longitude and latitude.
  const sketches = svg('g', {});
  for (const L of shared.layers.filter((x) => x.result.kind === 'sketch' && x.dataset === shared.flows)) {
    for (const p of L.result.polylines || []) {
      const c = p.coords || [];
      let d = '';
      for (let i = 0; i + 1 < c.length; i += 2) { const [x, y] = project(c[i], c[i + 1]); d += `${i ? 'L' : 'M'}${x.toFixed(1)},${y.toFixed(1)}`; }
      put(sketches, svg('path', { class: 'sketch-line', d, stroke: p.colour || PART.plugin, 'stroke-width': p.width || 1.5, 'stroke-opacity': p.alpha ?? 0.9,
        'stroke-dasharray': p.dashed ? '4 3' : null }));
    }
    for (const m of L.result.markers || []) {
      const [x, y] = project(m.lon, m.lat);
      put(sketches, svg('circle', { cx: x.toFixed(1), cy: y.toFixed(1), r: ((m.size || 4) * px).toFixed(2), fill: m.colour || PART.plugin, 'fill-opacity': 0.85 }));
    }
    for (const t of L.result.labels || []) {
      const [x, y] = project(t.lon, t.lat);
      put(sketches, svg('text', { class: 'label', x: x.toFixed(1), y: (y - 5 * px).toFixed(1), 'font-size': (12 * px).toFixed(2), 'stroke-width': (3 * px).toFixed(2), text: t.text }));
    }
  }
  // The largest flows' places, named.
  const labels = svg('g', {});
  if (state.opts.labels) {
    // From the largest down, each name where it touches none placed before it.
    const order = A.economies.map((e, i) => [through[i], i]).sort((a, b) => b[0] - a[0]).slice(0, Math.round(16 * sk));
    const placed = [];
    for (const [, i] of order) {
      const e = A.economies[i];
      const [x, y] = project(e.longitude, e.latitude);
      const text = e.name.length > 18 ? e.iso3 : e.name;
      const w = 0.5 * 12 * px * text.length, h = 13 * px, box = [x - w / 2, y - 6 * px - h, x + w / 2, y - 6 * px + 2 * px];
      if (placed.some((b) => box[0] < b[2] && b[0] < box[2] && box[1] < b[3] && b[1] < box[3])) continue;
      placed.push(box);
      put(labels, svg('text', { class: 'label', x: x.toFixed(1), y: (y - 6 * px).toFixed(1), 'text-anchor': 'middle', 'font-size': (12 * px).toFixed(2),
        'stroke-width': (3 * px).toFixed(2), text }));
    }
  }
  const m = 4 * px;
  put(over, svg('defs', {},
    svg('marker', { id: 'head', viewBox: '0 0 10 10', refX: 8, refY: 5, markerUnits: 'userSpaceOnUse', markerWidth: (m * 2).toFixed(2), markerHeight: (m * 2).toFixed(2), orient: 'auto-start-reverse' },
      svg('path', { d: 'M0,1L9,5L0,9z', fill: 'context-stroke', 'fill-opacity': 0.85 })),
    svg('marker', { id: 'eddy-head', viewBox: '0 0 10 10', refX: 6, refY: 5, markerUnits: 'userSpaceOnUse', markerWidth: (m * 1.6).toFixed(2), markerHeight: (m * 1.6).toFixed(2), orient: 'auto' },
      svg('path', { d: 'M0,1L9,5L0,9z', fill: PART.local }))),
    faces, places, arcLayer, sketches, labels);
}

// ---- zoom and pan ---------------------------------------------------------------------------------------------

function applyView() {
  const world = document.getElementById('world');
  if (world) world.setAttribute('transform', `translate(${state.view.x.toFixed(2)},${state.view.y.toFixed(2)}) scale(${state.view.k.toFixed(4)})`);
}

let settle = null;
function zoomTo(k, px, py) {
  const v = state.view;
  const k2 = Math.max(1, Math.min(18, k));
  v.x = px - (px - v.x) * (k2 / v.k);
  v.y = py - (py - v.y) * (k2 / v.k);
  v.k = k2;
  if (k2 === 1) { v.x = 0; v.y = 0; }
  applyView();
  readout(null);
  clearTimeout(settle);
  settle = setTimeout(() => drawOver(), 140);
}

// A dataset of a region: the view fitted to its places; the world's, the whole map.
function fitToPlaces() {
  const A = state.answer;
  if (!A) return;
  const pts = A.economies.map((e) => project(e.longitude, e.latitude));
  const xs = pts.map((p) => p[0]), ys = pts.map((p) => p[1]);
  const bw = Math.max(...xs) - Math.min(...xs), bh = Math.max(...ys) - Math.min(...ys);
  if (bw > 1.2 * W || bh > 1.2 * H) { state.view = { k: 1, x: 0, y: 0 }; applyView(); return; }
  const k = Math.max(1, Math.min(18, 0.8 * (2 * W) / Math.max(bw, 40), 0.8 * (2 * H) / Math.max(bh, 30)));
  const cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
  state.view = { k, x: -cx * k, y: -cy * k };
  applyView();
}

function resetView() {
  state.view = { k: 1, x: 0, y: 0 };
  applyView();
  drawOver();
}

function svgPoint(host, event) {
  const p = host.createSVGPoint();
  p.x = event.clientX;
  p.y = event.clientY;
  return p.matrixTransform(host.getScreenCTM().inverse());
}

// The panels fold away, so that the map takes the room; the choice is kept in this browser.
function panes() {
  const layout = document.querySelector('.map-layout');
  const stage = document.querySelector('.map-stage');
  const remembered = (k) => { try { return localStorage.getItem(`economy-atlas-${k}`) === 'closed'; } catch { return false; } };
  const set = (k, closed) => {
    layout.classList.toggle(`${k}-closed`, closed);
    try { localStorage.setItem(`economy-atlas-${k}`, closed ? 'closed' : 'open'); } catch { /* for this visit only */ }
    render();
    drawOver();
  };
  const button = (k, open, shut) => el('button', { 'data-pane': k, onclick: () => set(k, !layout.classList.contains(`${k}-closed`)) });
  const dockB = button('dock'), panelB = button('panel');
  const render = () => {
    dockB.textContent = layout.classList.contains('dock-closed') ? '› controls' : '‹ controls';
    dockB.title = layout.classList.contains('dock-closed') ? 'Show the controls' : 'Hide the controls';
    panelB.textContent = layout.classList.contains('panel-closed') ? 'inspector ‹' : 'inspector ›';
    panelB.title = layout.classList.contains('panel-closed') ? 'Show the inspector' : 'Hide the inspector';
  };
  layout.classList.toggle('dock-closed', remembered('dock'));
  layout.classList.toggle('panel-closed', remembered('panel'));
  render();
  put(stage, el('div', { class: 'panes' }, dockB, panelB), el('div', { class: 'readout', id: 'readout', 'aria-live': 'off' }));
}

// The projection, the place under the pointer and the zoom.
function readout(e) {
  const box = document.getElementById('readout');
  if (!box) return;
  const host = document.getElementById('map');
  let where = '';
  if (e) {
    const p = svgPoint(host, e);
    const x = (p.x - state.view.x) / state.view.k / SCALE, y = -(p.y - state.view.y) / state.view.k / SCALE;
    const ll = inverse(equalEarth, x, y);
    if (ll) {
      const [lon, lat] = ll;
      where = ` · ${Math.abs(lat).toFixed(1)}° ${lat >= 0 ? 'N' : 'S'}, ${Math.abs(lon).toFixed(1)}° ${lon >= 0 ? 'E' : 'W'}`;
    }
  }
  fill(box, el('b', { text: 'Equal Earth' }), `${where} · ×${state.view.k.toFixed(1)}`);
}

function wire() {
  const host = document.getElementById('map');
  host.addEventListener('wheel', (e) => {
    e.preventDefault();
    const p = svgPoint(host, e);
    zoomTo(state.view.k * Math.exp(-e.deltaY * 0.0016), p.x, p.y);
  }, { passive: false });
  let drag = null, moved = false;
  host.addEventListener('pointerdown', (e) => {
    if (e.button !== 0) return;
    drag = { p: svgPoint(host, e), x: state.view.x, y: state.view.y };
    moved = false;
  });
  host.addEventListener('pointermove', (e) => {
    if (drag) {
      const p = svgPoint(host, e);
      const dx = p.x - drag.p.x, dy = p.y - drag.p.y;
      if (!moved && Math.hypot(dx, dy) > 4) { moved = true; host.setPointerCapture(e.pointerId); host.classList.add('dragging'); hideTip(); }
      if (moved) {
        state.view.x = drag.x + dx;
        state.view.y = drag.y + dy;
        applyView();
        return;
      }
    }
    tip(e);
    readout(e);
  });
  const end = () => { drag = null; host.classList.remove('dragging'); };
  host.addEventListener('pointerup', (e) => {
    if (drag && !moved) choose(e);
    end();
  });
  host.addEventListener('pointerleave', () => { hideTip(); readout(null); });
  host.addEventListener('dblclick', resetView);
  document.getElementById('zoom-in').onclick = () => zoomTo(state.view.k * 1.6, 0, 0);
  document.getElementById('zoom-out').onclick = () => zoomTo(state.view.k / 1.6, 0, 0);
  document.getElementById('zoom-reset').onclick = resetView;
  window.addEventListener('keydown', (e) => {
    if (document.querySelector('.view[data-view="map"]').hidden || e.target.closest?.('input, select, textarea, dialog')) return;
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
      const i = state.years.indexOf(state.year) + (e.key === 'ArrowLeft' ? -1 : 1);
      if (i >= 0 && i < state.years.length) { e.preventDefault(); setYear(state.years[i]); }
    }
    if (e.key === 'Escape' && state.selected !== null) { state.selected = null; draw(); panel(); }
  });
  onTheme(() => { if (state.answer) { draw(); strip(); panel(); } });
}

function choose(e) {
  const t = e.target.closest('[data-i]');
  const i = t && t.dataset.i !== '' ? Number(t.dataset.i) : null;
  if (i === null) return;
  state.selected = state.selected === i ? null : i;
  draw();
  panel();
}

// ---- the tooltip --------------------------------------------------------------------------------------------

function tip(e) {
  const A = state.answer;
  const box = document.getElementById('tip');
  if (!A) return;
  const row = (a, b) => el('div', { class: 'row' }, el('span', { text: a }), el('span', { text: b }));
  const place = e.target.closest('[data-i]');
  const arc = e.target.closest('[data-arc]');
  let content = null;
  if (place && place.dataset.i !== '') {
    const x = A.economies[Number(place.dataset.i)];
    const C = colouring(A);
    const v = C.value(x);
    content = [el('b', { text: x.name }),
      C.label !== 'α' ? row(C.label, Number.isFinite(v) ? C.show(v) : '—') : null,
      row('attraction α', signed(x.attraction, 3)), row('exports', amount(x.exports)), row('imports', amount(x.imports)),
      row('net inflow', amount(x.balance)), row('circulating', percent(x.circulating))];
  } else if (arc) {
    const a = state.drawn[Number(arc.dataset.arc)];
    if (a) {
      const from = A.economies[a.from], to = A.economies[a.to];
      content = [el('b', { text: `${from.name} → ${to.name}` }), row(a.title || 'net flow', amount(a.value))];
      if (a.parts) content.push(row('driven', amount(Math.abs(a.parts[3]))), row('round triangles', amount(Math.abs(a.parts[4]))), row('round holes', amount(Math.abs(a.parts[5]))));
    }
  } else if (e.target.classList?.contains('country')) {
    content = [el('b', { text: e.target.dataset.name }), el('span', { class: 'note', text: 'not in these data' })];
  }
  if (!content) return hideTip();
  fill(box, ...content.filter(Boolean));
  box.classList.add('on');
  const x = Math.min(window.innerWidth - box.offsetWidth - 12, e.clientX + 16), y = Math.min(window.innerHeight - box.offsetHeight - 12, e.clientY + 14);
  box.style.left = `${x}px`;
  box.style.top = `${y}px`;
}
const hideTip = () => document.getElementById('tip').classList.remove('on');

// ---- the strip of years ---------------------------------------------------------------------------------------

function strip() {
  const box = document.getElementById('years');
  if (state.years.length < 2) { box.hidden = true; return; }
  box.hidden = false;
  const byYear = new Map((state.timeline.get(shared.flows) || []).map((t) => [t.year, t]));
  const short = state.years.length > 11;
  fill(box,
    el('div', { class: 'caption' },
      el('button', { class: 'play', 'aria-pressed': String(!!state.playing), title: state.playing ? 'Stop' : 'Play through the years', onclick: play },
        state.playing ? '❚❚ stop' : '▶ play'),
      el('div', {}, 'The net flows’ energy, year by year: ', el('span', { style: { color: PART.driven }, text: 'driven' }), ', ',
        el('span', { style: { color: PART.local }, text: 'local' }), ', ', el('span', { style: { color: PART.global }, text: 'global' }), '.')),
    el('div', { class: 'strip', role: 'listbox', 'aria-label': 'Years' }, state.years.map((y, n) => {
      const t = byYear.get(y);
      const stack = el('div', { class: 'stack' }, t ? [['driven', t.driven], ['local', t.local], ['global', t.global]]
        .map(([k, v]) => el('span', { style: { height: `${(100 * v).toFixed(1)}%`, background: PART[k] } })) : []);
      const label = !y ? '—' : short && n > 0 && n < state.years.length - 1 ? `’${String(y).slice(-2)}` : String(y);
      return el('div', { class: `year${y === state.year ? ' on' : ''}${t ? '' : ' pending'}`, role: 'option', 'aria-selected': String(y === state.year),
        title: t ? `${y}: driven ${percent(t.driven)}, round triangles ${percent(t.local)}, round holes ${percent(t.global)}; reciprocal ${percent(t.reciprocal / t.total)} of all` : String(y),
        onclick: () => setYear(y) }, stack, el('div', { class: 'label', text: label }));
    })));
}

// The years one after another, from the one shown (or from the first, at the last), each held a moment;
// pressed again, it stops.
async function play() {
  if (state.playing) { state.playing = null; strip(); return; }
  const run = (state.playing = {});
  let i = state.years.indexOf(state.year);
  if (i < 0 || i >= state.years.length - 1) i = -1;
  strip();
  while (state.playing === run && ++i < state.years.length) {
    const t0 = performance.now();
    await setYear(state.years[i]);
    await new Promise((r) => setTimeout(r, Math.max(0, 1400 - (performance.now() - t0))));
  }
  if (state.playing === run) { state.playing = null; strip(); }
}

// The next year's decomposition asked for while this one is looked at.
function prefetch() {
  const next = state.years[state.years.indexOf(state.year) + 1];
  const key = `${shared.flows}|${next}`;
  if (next === undefined || state.cache.has(key)) return;
  const flows = shared.flows;
  api.map(next, flows).then((A) => state.cache.set(`${flows}|${next}`, A)).catch(() => {});
}

async function loadTimeline() {
  const id = shared.flows;
  if (state.timeline.has(id) || state.years.length < 2) return strip();
  strip();
  try {
    const T = await api.timeline(id);
    state.timeline.set(id, T.years);
    if (id === shared.flows) strip();
  } catch {
    // The strip stays pending; the map does not need it.
  }
}

// ---- the panel ------------------------------------------------------------------------------------------------

function list(title, rows, onPick) {
  return el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: title })),
    el('table', {}, el('tbody', {}, rows.map(([a, b, i]) => el('tr', { class: onPick && i !== undefined ? 'pick' : '', onclick: onPick && i !== undefined ? () => onPick(i) : null },
      el('td', {}, a), el('td', { class: 'num', text: b }))))));
}

const pick = (i) => { state.selected = i; draw(); panel(); };

function partsBar(l) {
  const d = Math.abs(l[3]), c = Math.abs(l[4] + l[5]), s = d + c || 1;
  return el('span', { class: 'mini', title: `driven ${percent(d / s, 0)}, circulating ${percent(c / s, 0)}` },
    el('span', { style: { width: `${(100 * d) / s}%`, background: PART.driven } }), el('span', { style: { width: `${(100 * c) / s}%`, background: PART.global } }));
}

function pluginResults() {
  const mine = shared.layers.filter((l) => l.dataset === shared.flows);
  if (!mine.length) return null;
  const byFile = new Map();
  for (const l of mine) byFile.set(l.file, [...(byFile.get(l.file) || []), l]);
  const show = (fn) => () => { fn(); refreshDock(); draw(); panel(); };
  return el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: 'from the plug-ins and the toolbox' }),
    el('button', { text: 'clear', onclick: show(() => { shared.layers = shared.layers.filter((l) => l.dataset !== shared.flows); shared.colour = 'attraction'; shared.arcs = 'all'; }) })),
  [...byFile.values()].map((ls) => el('div', { style: { marginBottom: '12px' } },
    el('h3', {}, ls[0].plugin, ' ', el('span', { class: `badge ${ls[0].language === 'cpp' ? 'cpp' : ls[0].language === 'tool' ? 'economics' : 'py'}`,
      text: ls[0].language === 'cpp' ? 'C++' : ls[0].language === 'tool' ? 'toolbox' : 'Py' }),
      ls[0].year ? el('span', { class: 'note', text: ` · ${ls[0].year}` }) : null),
    ls.map((l) => {
      const r = l.result;
      if (r.kind === 'note') return el('p', { class: 'note', text: r.text });
      if (r.kind === 'table') {
        return el('div', {}, el('div', { class: 'note', text: r.title }),
          el('table', {}, el('thead', {}, el('tr', {}, (r.columns || []).map((c, i) => el('th', { class: i ? 'num' : '', text: String(c) })))),
            el('tbody', {}, (r.rows || []).map((row) => el('tr', {}, row.map((v, i) => el('td', { class: i ? 'num' : '', text: typeof v === 'number' ? fmt(v) : String(v ?? '') })))))),
          r.note ? el('p', { class: 'note', text: r.note }) : null);
      }
      const shown = r.kind === 'values' ? shared.colour === l.key : r.kind === 'arcs' ? shared.arcs === l.key : true;
      return el('p', {}, el('button', { class: shown ? 'primary' : 'quiet', text: `${r.kind === 'values' ? 'colour by' : r.kind === 'arcs' ? 'draw' : 'sketch'}: ${l.title}`,
        onclick: show(() => { if (r.kind === 'values') shared.colour = shown ? 'attraction' : l.key; if (r.kind === 'arcs') shared.arcs = shown ? 'all' : l.key; }) }),
      r.note ? el('span', { class: 'note', text: ` ${r.note}` }) : null);
    }))));
}

function panel() {
  const box = document.getElementById('map-panel');
  const A = state.answer;
  clear(box);
  if (!A) return;
  const T = A.totals;
  const library = shared.flows === 'library:trade';
  const d = dataset(shared.flows);
  if (state.selected !== null) {
    const i = state.selected, e = A.economies[i];
    const mine = A.legs.filter((l) => l[0] === i || l[1] === i).sort((a, b) => b[2] - a[2]);
    put(box, el('h2', { text: e.name }), el('p', { class: 'lede' }, `${e.continent ? `${e.continent} · ` : ''}${A.year || ''} `,
      el('button', { class: 'quiet', onclick: () => { state.selected = null; draw(); panel(); } }, 'the whole map')),
    el('dl', {},
      el('dt', { text: 'exports' }), el('dd', { text: amount(e.exports) }),
      el('dt', { text: 'imports' }), el('dd', { text: amount(e.imports) }),
      el('dt', { text: 'net inflow (imports − exports)' }), el('dd', { text: amount(e.balance) }),
      el('dt', { text: 'attraction α' }), el('dd', { text: signed(e.attraction, 3) }),
      el('dt', { text: 'circulating, of the parts on its links' }), el('dd', { text: percent(e.circulating) })),
    el('p', {}, e.attraction > 0 ? 'Value is drawn to it: the driven part of net flows runs up to it from the places of lower attraction, along each link in proportion to the link’s gross flow.'
      : 'Value is sent from it: the driven part of net flows runs down from it to the places of higher attraction, along each link in proportion to the link’s gross flow.'),
    list('its largest net flows: driven, circulating', mine.slice(0, 14).map((l) => {
      const j = l[0] === i ? l[1] : l[0];
      return [[`${l[0] === i ? 'to' : 'from'} ${A.economies[j].name}`, partsBar(l)], amount(l[2]), j];
    }), pick),
    placeHistory(e),
    pluginResults());
    return;
  }
  put(box,
    el('h2', { text: `${library ? 'World trade in goods' : d?.title || 'Flows'}${A.year ? `, ${A.year}` : ''}` }),
    el('p', { class: 'lede', text: `${A.economies.length} places, ${amount(T.total)} flowing${!library && A.unit ? ` (${A.unit})` : ''}. What two places send each other both ways is exchange; what is left, one net flow for each pair, is driven or circulates.` }),
    el('div', { class: 'heading' }, el('span', { text: 'exchange and net transfer' })),
    shareBar([{ label: 'reciprocal exchange', value: T.reciprocal, colour: PART.reciprocal }, { label: 'net flows', value: T.net, colour: PART.driven }]),
    el('div', { class: 'heading', style: { marginTop: '14px' } }, el('span', { text: 'the net flows, by energy' })),
    shareBar([{ label: 'driven', value: T.driven, colour: PART.driven }, { label: 'round triangles', value: T.local, colour: PART.local }, { label: 'round holes', value: T.global, colour: PART.global }]),
    el('p', {}, `Each link conducts in proportion to its gross flow. The driven part runs up a single potential, the attraction α: it is the flow the balances alone would drive through those links, and between two places it is their gross flow times their difference in α. The rest circulates: round the map's ${A.faces.triangles} triangles of neighbours (local eddies) and round the loops through distant partners that no triangle fills (global circulation). Value round a loop returns to where it started; it is moved by the structure of the flows, not by any one place's balance.`),
    pluginResults(),
    list(`drawing value most (α)${library ? ', of those trading over 20 bn $' : ''}`, ranked(A, (a, b) => b.attraction - a.attraction, library).map((x) => [x.name, signed(x.attraction, 3), x.i]), pick),
    list('sending value most (α)', ranked(A, (a, b) => a.attraction - b.attraction, library).map((x) => [x.name, signed(x.attraction, 3), x.i]), pick),
    list('circulating most, of the parts on their links', ranked(A, (a, b) => b.circulating - a.circulating, library).map((x) => [x.name, percent(x.circulating), x.i]), pick),
    el('p', { class: 'note', text: `Decomposed ${A.method === 'sparse' ? 'by conjugate gradients' : 'exactly'} on the complex of ${A.economies.length} places, ${A.legs.length} net flows and ${A.faces.triangles} triangles. Source: ${A.source?.citation || ''}` }),
    downloads(A),
  );
}

// A place over the years, from the toolbox's table of economies: its trade, and the indicator the
// map is coloured by, each on its own axes.
function placeHistory(e) {
  if (shared.flows !== 'library:trade') return null;
  const H = placeYears(e.iso3);
  if (!H || H.years.length < 2) return null;
  const toBn = (v) => (Number.isFinite(v) ? v / 1e3 : null);
  const charts = [lineChart({ title: 'Its trade in goods, bn US$', x: H.years, compact: true,
    series: [{ name: 'exports', y: H.column('exports').map(toBn) }, { name: 'imports', y: H.column('imports').map(toBn) }] })];
  if (shared.colour.startsWith('wb:')) {
    const name = shared.colour.slice(3), label = state.indicators.columns.find((c) => c.name === name)?.label || name;
    charts.push(lineChart({ title: label, x: H.years, compact: true, series: [{ name: 'World Bank', y: H.column(name) }] }));
  }
  return el('div', { class: 'section history' }, el('div', { class: 'heading' }, el('span', { text: 'over the years' })), charts);
}

// The year's numbers and the map itself, to take away.
function downloads(A) {
  const stem = `${shared.flows === 'library:trade' ? 'world-trade' : 'flows'}-${A.year || 'all'}`;
  const C = colouring(A);
  const places = () => download(`${stem}-places.csv`, csv(
    ['id', 'name', 'latitude', 'longitude', 'exports', 'imports', 'net_inflow', 'attraction', 'circulating', ...(C.label !== 'α' ? [C.label] : [])],
    A.economies.map((e) => [e.iso3, e.name, e.latitude, e.longitude, e.exports, e.imports, e.balance, e.attraction, e.circulating, ...(C.label !== 'α' ? [C.value(e)] : [])])));
  const flows = () => download(`${stem}-net-flows.csv`, csv(['from', 'to', 'net', 'driven', 'round_triangles', 'round_holes'],
    A.legs.map((l) => [A.economies[l[0]].iso3, A.economies[l[1]].iso3, l[2], l[3], l[4], l[5]])));
  return el('div', { class: 'section downloads' }, el('div', { class: 'heading' }, el('span', { text: 'take away' })),
    el('div', { class: 'chips' },
      el('button', { class: 'quiet', onclick: places, text: 'the places, CSV' }),
      el('button', { class: 'quiet', onclick: flows, text: 'the net flows, CSV' }),
      el('button', { class: 'quiet', onclick: async () => download(`${stem}-figure.svg`, await figureSvg(), 'image/svg+xml'), text: 'the figure, SVG' }),
      el('button', { class: 'quiet', onclick: figurePng(stem), text: 'the figure, PNG' })),
    el('p', { class: 'note', text: 'A figure is the map as shown, with its title, legend, source and link, at 1800 pixels wide; the PNG at twice that.' }));
}

// The map's drawing with the page's styles written into it, its arcs standing still.
function mapClone() {
  const host = document.getElementById('map');
  const copy = host.cloneNode(true);
  const from = host.querySelectorAll('*'), to = copy.querySelectorAll('*');
  const keep = ['fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-opacity', 'stroke-dasharray', 'stroke-linecap', 'stroke-linejoin', 'paint-order', 'opacity',
    'font-family', 'font-style', 'font-weight', 'font-variant', 'letter-spacing', 'vector-effect'];
  from.forEach((node, n) => {
    if (!(node instanceof SVGGraphicsElement) || node.closest('defs')) return;
    const s = getComputedStyle(node);
    const flowing = node.classList.contains('flowing');
    to[n].setAttribute('style', keep.map((p) => `${p}:${p === 'stroke-dasharray' && flowing ? 'none' : s.getPropertyValue(p)}`).join(';'));
    // On paper an arc does not move: it is drawn whole, with its head.
    if (flowing) to[n].setAttribute('marker-end', 'url(#head)');
    to[n].removeAttribute('class');
  });
  return copy;
}

// Source Serif 4 written into the figure, so that it is set in the page's type wherever it is
// opened; without the network, the figure names the font and falls back to Georgia.
let fontFaces = null;
function embeddedFonts() {
  fontFaces ??= (async () => {
    const link = document.querySelector('link[href*="fonts.googleapis.com/css2"]');
    if (!link) return '';
    const css = await (await fetch(link.href)).text();
    const faces = [...css.matchAll(/@font-face\s*{[^}]*}/g)].map((m) => m[0]).filter((f) => /latin \*\/|unicode-range: U\+0000-00FF/.test(f) || !/unicode-range/.test(f));
    const out = [];
    for (const face of (faces.length ? faces : [...css.matchAll(/@font-face\s*{[^}]*}/g)].map((m) => m[0])).slice(0, 4)) {
      const url = face.match(/url\(([^)]+)\)/)?.[1];
      if (!url) continue;
      const bytes = new Uint8Array(await (await fetch(url)).arrayBuffer());
      let bin = '';
      for (let i = 0; i < bytes.length; i += 0x8000) bin += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
      out.push(face.replace(url, `data:font/woff2;base64,${btoa(bin)}`));
    }
    return out.join('\n');
  })().catch(() => '');
  return fontFaces;
}

// The figure: a title, the map as shown, its legend, and a line of source, projection and link.
async function figureSvg() {
  const A = state.answer;
  const C = colouring(A);
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  const map = document.getElementById('map');
  const vb = map.viewBox.baseVal;
  const Wf = 1800, pad = 48, top = 104, mapW = Wf - 2 * pad, mapH = Math.round((mapW * vb.height) / vb.width), keyH = 150;
  const Hf = top + mapH + keyH;
  const font = "'Source Serif 4', Georgia, serif";
  const text = (x, y, t, size, extra = {}) => svg('text', { x, y, 'font-family': font, 'font-size': size, fill: css('--ink'), ...extra, text: t });
  const d = dataset(shared.flows);
  const title = `${shared.flows === 'library:trade' ? 'World trade in goods' : d?.title || 'Flows'}${A.year ? `, ${A.year}` : ''}`;
  const sub = `Places coloured by ${C.title.charAt(0).toLowerCase()}${C.title.slice(1)}; arcs: ${shared.arcs === 'none' ? 'none' : layer(shared.arcs) ? layer(shared.arcs).title : `the ${state.opts.shown} largest ${shared.arcs === 'all' ? 'net flows' : `${shared.arcs} parts of the net flows`}`}.`;
  const out = svg('svg', { width: Wf, height: Hf, viewBox: `0 0 ${Wf} ${Hf}` });
  const fonts = await embeddedFonts();
  if (fonts) put(out, svg('style', {}, fonts));
  // The ramp of the legend, as a gradient of the colours the map uses.
  const stops = C.kind === 'sequential' ? [0, 0.25, 0.5, 0.75, 1].map((t) => [t, sequential(t)]) : [-1, -0.5, -0.2, 0, 0.2, 0.5, 1].map((t) => [(t + 1) / 2, diverging(t)]);
  put(out, svg('defs', {}, svg('linearGradient', { id: 'figure-ramp' }, stops.map(([o, c]) => svg('stop', { offset: o, 'stop-color': c })))));
  put(out, svg('rect', { width: Wf, height: Hf, fill: css('--paper') }),
    text(pad, 56, title, 30, { 'font-weight': 600 }),
    text(pad, 86, sub, 17, { fill: css('--ink-2') }));
  const inner = mapClone();
  inner.setAttribute('x', pad);
  inner.setAttribute('y', top);
  inner.setAttribute('width', mapW);
  inner.setAttribute('height', mapH);
  inner.removeAttribute('id');
  put(out, inner);
  // The legend beneath: the ramp and its ends, then the arcs' keys.
  const ly = top + mapH + 30;
  put(out, svg('rect', { x: pad, y: ly, width: 360, height: 12, fill: 'url(#figure-ramp)' }),
    text(pad, ly - 10, C.title, 15, { 'font-weight': 500 }),
    text(pad, ly + 32, C.ends[0], 14, { fill: css('--ink-2') }),
    text(pad + 360, ly + 32, C.ends[1], 14, { fill: css('--ink-2'), 'text-anchor': 'end' }));
  const keys = layer(shared.arcs) ? [[PART.plugin, layer(shared.arcs).title]] : shared.arcs === 'none' ? []
    : [[PART.driven, 'mostly driven: from senders to drawers'], [PART.local, 'mostly circulating round the triangles of neighbours'], [PART.global, 'mostly circulating round holes, through distant links']];
  keys.forEach(([c, t], k) => put(out, svg('path', { d: `M${pad + 440} ${ly + 6 + 24 * k}h40`, stroke: c, 'stroke-width': 3, 'stroke-linecap': 'round' }),
    text(pad + 492, ly + 11 + 24 * k, t, 14, { fill: css('--ink-2') })));
  put(out, svg('path', { d: `M${pad} ${Hf - 40}H${Wf - pad}`, stroke: css('--border'), 'stroke-width': 1 }),
    text(pad, Hf - 18, `${sourceLine()} projection (Šavrič, Patterson and Jenny 2018). Decomposed by the economics domain of The Science Infrastructures. ${new Date().toISOString().slice(0, 10)}.`, 13, { fill: css('--ink-faint') }),
    text(Wf - pad, Hf - 18, location.href, 13, { fill: css('--ink-faint'), 'text-anchor': 'end' }));
  return `<?xml version="1.0" encoding="utf-8"?>\n${new XMLSerializer().serializeToString(out)}\n`;
}

// The figure drawn at twice its size into an image.
function figurePng(stem) {
  return async () => {
    const text = await figureSvg();
    const url = URL.createObjectURL(new Blob([text], { type: 'image/svg+xml' }));
    const img = new Image();
    img.onload = () => {
      const c = document.createElement('canvas');
      c.width = img.width * 2;
      c.height = img.height * 2;
      const g = c.getContext('2d');
      g.scale(2, 2);
      g.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      c.toBlob((b) => {
        const a = el('a', { href: URL.createObjectURL(b), download: `${stem}-figure.png` });
        document.body.append(a);
        a.click();
        a.remove();
      }, 'image/png');
    };
    img.src = url;
  };
}

const ranked = (A, order, library) => A.economies.map((e, i) => ({ ...e, i }))
  .filter((e) => !library || e.exports + e.imports > RANKED).sort(order).slice(0, 8);

// ---- the dock -------------------------------------------------------------------------------------------------

let pluginsBox = null;

// A place found by its name or code: chosen, and brought into view when the map is drawn closer.
function findPlace(text) {
  const A = state.answer;
  const q = text.trim().toLowerCase();
  if (!A || !q) return;
  const i = A.economies.findIndex((e) => e.name.toLowerCase() === q || e.iso3.toLowerCase() === q);
  const j = i >= 0 ? i : A.economies.findIndex((e) => e.name.toLowerCase().startsWith(q));
  if (j < 0) return;
  if (state.view.k > 1) {
    const [x, y] = project(A.economies[j].longitude, A.economies[j].latitude);
    state.view.x = -x * state.view.k;
    state.view.y = -y * state.view.k;
    applyView();
  }
  pick(j);
}

async function refreshDock() {
  const dock = document.getElementById('map-dock');
  const plugins = shared.layers.filter((l) => l.result.kind === 'values' && l.dataset === shared.flows).map((l) => [l.key, `${l.title} (${l.plugin})`]);
  const indicators = state.indicators ? state.indicators.columns.map((c) => [`wb:${c.name}`, c.label || c.name]) : [];
  const colourGroups = [
    ['the flows', [['attraction', 'Attraction α'], ...Object.entries(MEASURES).map(([k, [title]]) => [k, title])]],
    ...(indicators.length ? [['the World Bank’s indicators', indicators]] : []),
    ...(plugins.length ? [['plug-ins', plugins]] : []),
  ];
  const arcOptions = [['all', 'all net flows'], ['driven', 'their driven part'], ['circulating', 'their circulating part'], ['none', 'none'],
    ...shared.layers.filter((l) => l.result.kind === 'arcs' && l.dataset === shared.flows).map((l) => [l.key, `${l.title} (${l.plugin})`])];
  if (!colourGroups.some(([, os]) => os.some(([k]) => k === shared.colour)) && !(shared.colour.startsWith('wb:') && !state.indicators)) shared.colour = 'attraction';
  if (!arcOptions.some(([k]) => k === shared.arcs)) shared.arcs = 'all';
  const option = ([k, label], value) => el('option', { value: k, text: label, selected: k === value });
  const select = (options, value, onChange) => el('select', { onchange: (e) => onChange(e.target.value) }, options.map((o) => option(o, value)));
  const grouped = (groups, value, onChange) => el('select', { onchange: (e) => onChange(e.target.value) },
    groups.map(([label, os]) => el('optgroup', { label }, os.map((o) => option(o, value)))));
  const toggle = (key, label) => el('label', {}, el('input', { type: 'checkbox', checked: state.opts[key], onchange: (e) => { state.opts[key] = e.target.checked; draw(); } }), label);
  const names = el('datalist', { id: 'map-places' }, (state.answer?.economies || []).map((e) => el('option', { value: e.name })));
  if (!pluginsBox) pluginsBox = await pluginsSection(() => ({ dataset: shared.flows, year: state.year }));
  fill(dock,
    el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: 'data' })),
      datasetChooser('flows', () => changeDataset())),
    el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: 'layers' })),
      el('label', { class: 'field' }, el('span', { text: 'Find a place' }), el('input', { type: 'search', list: 'map-places', placeholder: 'a name or a code',
        onchange: (e) => findPlace(e.target.value), onkeydown: (e) => { if (e.key === 'Enter') findPlace(e.target.value); } }), names),
      el('label', { class: 'field' }, el('span', { text: 'Colour the places by' }), grouped(colourGroups, shared.colour, (v) => { shared.colour = v; draw(); panel(); })),
      el('label', { class: 'field' }, el('span', { text: 'Arcs' }), select(arcOptions, shared.arcs, (v) => { shared.arcs = v; draw(); })),
      el('label', { class: 'field' }, el('span', { text: 'How many' }), select([['60', '60'], ['150', '150'], ['400', '400'], ['1000', '1000']], String(state.opts.shown),
        (v) => { state.opts.shown = Number(v); draw(); })),
      el('div', { class: 'toggles' }, toggle('faces', 'the map’s triangles'), toggle('eddies', 'eddies round them'), toggle('labels', 'names of the largest'),
        toggle('flowing', 'arcs flowing'))),
    el('div', { class: 'section' }, pluginsBox),
  );
}

// ---- loading --------------------------------------------------------------------------------------------------

function yearsOf(id) {
  const d = dataset(id);
  return d?.years ? [...d.years] : [];
}

async function setYear(year) {
  state.year = year;
  strip();
  await load();
}

async function load() {
  const status = document.getElementById('status');
  const key = `${shared.flows}|${state.year}`;
  if (state.cache.has(key)) {
    state.answer = state.cache.get(key);
  } else {
    status.textContent = `decomposing the flows of ${state.year || 'the data'}…`;
    try {
      state.answer = await api.map(state.year, shared.flows);
      state.cache.set(key, state.answer);
      status.textContent = '';
    } catch (error) {
      status.textContent = error.message;
      return;
    }
  }
  if (state.wanted !== undefined) {
    const i = state.answer.economies.findIndex((e) => e.iso3 === state.wanted);
    state.selected = i >= 0 ? i : null;
    state.wanted = undefined;
  }
  if (state.fit) {
    state.fit = false;
    if (!state.built && state.geo) buildBase();
    fitToPlaces();
  }
  const names = document.getElementById('map-places');
  if (names) fill(names, state.answer.economies.map((e) => el('option', { value: e.name })));
  draw();
  panel();
  if (state.playing) prefetch();
}

async function changeDataset() {
  state.selected = null;
  state.years = yearsOf(shared.flows);
  state.year = state.years.at(-1) ?? null;
  state.fit = true;
  await refreshDock();
  loadTimeline();
  await load();
}

// The address may ask for a dataset, a year, the arcs and how many, the triangles, the names, the
// flowing, a place (its id; the ISO alpha-3 code for the library's economies) and the colouring:
// #map&year=2010&arcs=driven&shown=400&faces=0&economy=CHN&colour=wb:gdp_per_capita.
export async function showMap(data, params = {}) {
  if (!state.geo) {
    state.geo = await (await fetch(COUNTRIES)).json();
    panes();
    wire();
    readout(null);
    on('layers', () => { refreshDock(); draw(); panel(); });
    on('dataset-chosen', (kind) => { if (kind === 'flows') changeDataset(); });
    on('datasets', () => { if (!document.querySelector('.view[data-view="map"]').hidden) refreshDock(); });
    if (params.dataset && dataset(params.dataset)) shared.flows = params.dataset;
    state.fit = shared.flows !== 'library:trade';
    state.years = yearsOf(shared.flows);
    state.year = state.years.at(-1) ?? null;
    await refreshDock();
    loadTimeline();
    loadIndicators().then((T) => { if (T) { refreshDock(); if (shared.colour.startsWith('wb:')) { draw(); panel(); } } });
  } else if (params.dataset && params.dataset !== shared.flows && dataset(params.dataset)) {
    shared.flows = params.dataset;
    state.fit = true;
    state.years = yearsOf(shared.flows);
    state.year = state.years.at(-1) ?? null;
    loadTimeline();
  }
  const year = Number(params.year);
  if (state.years.includes(year)) state.year = year;
  if (['all', 'driven', 'circulating', 'none'].includes(params.arcs)) shared.arcs = params.arcs;
  if (['60', '150', '400', '1000'].includes(params.shown)) state.opts.shown = Number(params.shown);
  if (params.faces) state.opts.faces = params.faces !== '0';
  if (params.labels) state.opts.labels = params.labels !== '0';
  if (params.flowing) state.opts.flowing = params.flowing !== '0';
  if ('economy' in params) state.wanted = params.economy.toUpperCase();
  if (params.colour && (params.colour === 'attraction' || MEASURES[params.colour] || params.colour.startsWith('wb:'))) shared.colour = params.colour;
  if (Object.keys(params).length) await refreshDock();
  strip();
  await load();
}
