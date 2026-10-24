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
// Places are coloured by α, or by the values a plug-in returned; places with a country of their own
// colour it, the others are circles. Arcs are the largest net flows, their driven or circulating
// parts, or a plug-in's arcs; small circular arrows on the triangles are the eddies of local
// circulation. The strip beneath gives every year's shares; the map zooms and pans.

import * as api from './api.js';
import { el, svg, clear, put, fill, money, quantity, percent, signed, shareBar } from './dom.js';
import { diverging, sequential, rampCss, PART, onTheme } from './palette.js';
import { shared, on, dataset } from './state.js';
import { datasetChooser } from './data.js';
import { pluginsSection } from './plugins.js';
import { equalEarth, extent, outline } from '/domains/geography/js/projections.js';

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
  drawn: [], built: false, covered: new Set(),
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

function colouring(A) {
  const L = shared.colour !== 'attraction' ? layer(shared.colour) : null;
  if (L && L.result.kind === 'values') {
    const values = L.result.values || {};
    const xs = Object.values(values).filter(Number.isFinite);
    const seq = L.result.scale === 'sequential';
    const lo = Math.min(...xs), hi = Math.max(...xs), big = Math.max(1e-12, ...xs.map(Math.abs));
    return {
      title: `${L.title}${L.result.unit ? ` (${L.result.unit})` : ''}`, kind: seq ? 'sequential' : 'diverging',
      ends: seq ? [fmt(lo), fmt(hi)] : [`−${fmt(big)}`, `+${fmt(big)}`],
      of: (e) => { const v = values[e.iso3]; return Number.isFinite(v) ? (seq ? sequential(hi > lo ? (v - lo) / (hi - lo) : 0.5) : diverging(v / big)) : null; },
      value: (e) => values[e.iso3], label: L.title,
    };
  }
  const maxA = Math.max(1e-9, ...A.economies.map((e) => Math.abs(e.attraction)));
  return { title: 'Attraction α, a share of gross flow', kind: 'diverging', ends: [`sends, −${maxA.toFixed(2)}`, `draws, +${maxA.toFixed(2)}`],
    of: (e) => diverging(e.attraction / maxA), value: (e) => e.attraction, label: 'α' };
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
  fill(box, 
    el('div', { class: 'title', text: C.title }),
    el('div', { class: 'ramp', style: { background: rampCss(C.kind) } }),
    el('div', { class: 'ends' }, el('span', { text: C.ends[0] }), el('span', { text: C.ends[1] })),
    L ? key(PART.plugin, `arcs: ${L.title}`) : shared.arcs === 'none' ? null : [
      key(PART.driven, 'mostly driven: from senders to drawers'),
      key(PART.local, 'mostly circulating round the triangles of neighbours'),
      key(PART.global, 'mostly circulating round holes, through distant links'),
      state.opts.eddies ? key(PART.local, 'eddies round the triangles', '2 3') : null],
  );
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
  });
  const end = () => { drag = null; host.classList.remove('dragging'); };
  host.addEventListener('pointerup', (e) => {
    if (drag && !moved) choose(e);
    end();
  });
  host.addEventListener('pointerleave', () => { hideTip(); });
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
      C.label !== 'α' ? row(C.label, fmt(v)) : null,
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
    el('div', { class: 'caption' }, 'The net flows’ energy, year by year: ', el('span', { style: { color: PART.driven }, text: 'driven' }), ', ',
      el('span', { style: { color: PART.local }, text: 'local' }), ', ', el('span', { style: { color: PART.global }, text: 'global' }), '.'),
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
  return el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: 'from the plug-ins' }),
    el('button', { text: 'clear', onclick: show(() => { shared.layers = shared.layers.filter((l) => l.dataset !== shared.flows); shared.colour = 'attraction'; shared.arcs = 'all'; }) })),
  [...byFile.values()].map((ls) => el('div', { style: { marginBottom: '12px' } },
    el('h3', {}, ls[0].plugin, ' ', el('span', { class: `badge ${ls[0].language === 'cpp' ? 'cpp' : 'py'}`, text: ls[0].language === 'cpp' ? 'C++' : 'Py' }),
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
  );
}

const ranked = (A, order, library) => A.economies.map((e, i) => ({ ...e, i }))
  .filter((e) => !library || e.exports + e.imports > RANKED).sort(order).slice(0, 8);

// ---- the dock -------------------------------------------------------------------------------------------------

let pluginsBox = null;

async function refreshDock() {
  const dock = document.getElementById('map-dock');
  const colourOptions = [['attraction', 'Attraction α'], ...shared.layers.filter((l) => l.result.kind === 'values' && l.dataset === shared.flows).map((l) => [l.key, `${l.title} (${l.plugin})`])];
  const arcOptions = [['all', 'all net flows'], ['driven', 'their driven part'], ['circulating', 'their circulating part'], ['none', 'none'],
    ...shared.layers.filter((l) => l.result.kind === 'arcs' && l.dataset === shared.flows).map((l) => [l.key, `${l.title} (${l.plugin})`])];
  if (!colourOptions.some(([k]) => k === shared.colour)) shared.colour = 'attraction';
  if (!arcOptions.some(([k]) => k === shared.arcs)) shared.arcs = 'all';
  const select = (options, value, onChange) => el('select', { onchange: (e) => onChange(e.target.value) },
    options.map(([k, label]) => el('option', { value: k, text: label, selected: k === value })));
  const toggle = (key, label) => el('label', {}, el('input', { type: 'checkbox', checked: state.opts[key], onchange: (e) => { state.opts[key] = e.target.checked; draw(); } }), label);
  if (!pluginsBox) pluginsBox = await pluginsSection(() => ({ dataset: shared.flows, year: state.year }));
  fill(dock, 
    el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: 'data' })),
      datasetChooser('flows', () => changeDataset())),
    el('div', { class: 'section' }, el('div', { class: 'heading' }, el('span', { text: 'layers' })),
      el('label', { class: 'field' }, el('span', { text: 'Colour the places by' }), select(colourOptions, shared.colour, (v) => { shared.colour = v; draw(); panel(); })),
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
  draw();
  panel();
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
// flowing, and a place (its id; the ISO alpha-3 code for the library's economies):
// #map&year=2010&arcs=driven&shown=400&faces=0&economy=CHN.
export async function showMap(data, params = {}) {
  if (!state.geo) {
    state.geo = await (await fetch(COUNTRIES)).json();
    wire();
    on('layers', () => { refreshDock(); draw(); panel(); });
    on('dataset-chosen', (kind) => { if (kind === 'flows') changeDataset(); });
    on('datasets', () => { if (!document.querySelector('.view[data-view="map"]').hidden) refreshDock(); });
    if (params.dataset && dataset(params.dataset)) shared.flows = params.dataset;
    state.fit = shared.flows !== 'library:trade';
    state.years = yearsOf(shared.flows);
    state.year = state.years.at(-1) ?? null;
    await refreshDock();
    loadTimeline();
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
  if (Object.keys(params).length) await refreshDock();
  strip();
  await load();
}
