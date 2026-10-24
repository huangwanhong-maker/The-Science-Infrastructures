// The map: the world's trade in a year on an Equal Earth map (the geography domain's projection).
//
// The engine splits the year's exports into reciprocal exchange (what two economies send each other
// both ways, a round trip) and net flows, one for each pair; and the net flows, on a complex whose
// faces are the spherical Delaunay triangles of the economies' places, into the driven part, local
// circulation round the triangles, and global circulation round the holes the map leaves. Each
// link conducts in proportion to its gross trade (Ohm's law), so the driven net flow between two
// economies is their gross trade times their difference in attraction α: α is a share of gross
// trade, without unit, the potential from which value is drawn to where it is drawn to. Countries
// are coloured by α; arcs are the largest net flows, or their driven or circulating parts; small
// circular arrows on the triangles are the eddies of local circulation.

import * as api from './api.js';
import { el, svg, clear, money, percent, fixed, signed, shareBar } from './dom.js';
import { diverging, rampCss, PART } from './palette.js';
import { equalEarth, extent, outline } from '/domains/geography/js/projections.js';

const COUNTRIES = '/data/natural-earth-50m/ne_50m_admin_0_countries.geojson';
const SCALE = 300;   // plane units to SVG units
const [XM, YM] = extent(equalEarth);

const state = { geo: null, answer: null, selected: null, cache: new Map() };

// The raw projection, not `forward`, which wraps +180° to −180°: Natural Earth splits its polygons at
// the antimeridian, and a ring's vertices on +180° must stay on the map's right edge.
const DEG = Math.PI / 180;
const project = (lon, lat) => {
  const [x, y] = equalEarth.raw(Math.max(-180, Math.min(180, lon)) * DEG, lat * DEG);
  return [x * SCALE, -y * SCALE];
};

function economyCode(p) {
  for (const k of ['ISO_A3', 'ISO_A3_EH']) if (p[k] && p[k] !== '-99' && p[k].length === 3) return p[k];
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

// The great-circle arc between two places, projected, split where it crosses the map's edge.
function arcPath(a, b, bulge = 0.12) {
  const toVec = (lat, lon) => {
    const f = (lat * Math.PI) / 180, l = (lon * Math.PI) / 180;
    return [Math.cos(f) * Math.cos(l), Math.cos(f) * Math.sin(l), Math.sin(f)];
  };
  const u = toVec(a.latitude, a.longitude), v = toVec(b.latitude, b.longitude);
  const omega = Math.acos(Math.max(-1, Math.min(1, u[0] * v[0] + u[1] * v[1] + u[2] * v[2])));
  if (omega < 1e-6) return '';
  const n = Math.max(8, Math.ceil(omega * 40));
  let d = '', prev = null;
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    const s0 = Math.sin((1 - t) * omega) / Math.sin(omega), s1 = Math.sin(t * omega) / Math.sin(omega);
    const p = [s0 * u[0] + s1 * v[0], s0 * u[1] + s1 * v[1], s0 * u[2] + s1 * v[2]];
    const lat = (Math.asin(Math.max(-1, Math.min(1, p[2]))) * 180) / Math.PI, lon = (Math.atan2(p[1], p[0]) * 180) / Math.PI;
    let [x, y] = project(lon, lat);
    y -= bulge * SCALE * Math.sin(Math.PI * t) * omega * 0.35;  // a gentle lift, so that arcs separate
    const jump = prev && Math.abs(x - prev[0]) > XM * SCALE;
    d += `${!prev || jump ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`;
    prev = [x, y];
  }
  return d;
}

// Trading enough to rank: over 20 bn $ of exports and imports.
const RANKED = 2e4;

function legend(maxAttraction) {
  const box = document.getElementById('map-legend');
  clear(box).append(
    el('div', {}, 'Attraction α, a share of gross trade'),
    el('div', { class: 'ramp', style: { background: rampCss() } }),
    el('div', { class: 'ends' }, el('span', {}, `sends, −${maxAttraction.toFixed(2)}`), el('span', {}, `draws, +${maxAttraction.toFixed(2)}`)),
    el('div', { class: 'key' }, svg('svg', { width: 28, height: 10 }, svg('path', { d: 'M1 5H27', stroke: PART.driven, 'stroke-width': 2.4 })), 'driven: down from senders, up to drawers'),
    el('div', { class: 'key' }, svg('svg', { width: 28, height: 10 }, svg('path', { d: 'M1 5H27', stroke: PART.local, 'stroke-width': 2.4 })), 'circulating: round triangles and holes'),
  );
}

function draw() {
  const host = document.getElementById('map');
  const A = state.answer;
  if (!state.geo || !A) return;
  const W = XM * SCALE, H = YM * SCALE;
  host.setAttribute('viewBox', `${-W - 6} ${-H - 6} ${2 * W + 12} ${2 * H + 12}`);
  clear(host);
  const sphere = outline(equalEarth, 120).map(([x, y], i) => `${i ? 'L' : 'M'}${(x * SCALE).toFixed(1)},${(-y * SCALE).toFixed(1)}`).join('') + 'Z';
  host.append(svg('path', { class: 'sphere', d: sphere }));
  const byCode = new Map(A.economies.map((e, i) => [e.iso3, i]));
  const maxA = Math.max(1e-9, ...A.economies.map((e) => Math.abs(e.attraction)));
  legend(maxA);
  const countries = svg('g', {});
  for (const f of state.geo.features) {
    const code = economyCode(f.properties);
    const i = code ? byCode.get(code) : undefined;
    const e = i === undefined ? null : A.economies[i];
    countries.append(svg('path', {
      class: `country${state.selected === i && i !== undefined ? ' selected' : ''}`,
      d: featurePath(f.geometry),
      fill: e ? diverging(e.attraction / maxA) : 'var(--land)',
      onclick: () => { if (i !== undefined) { state.selected = i; draw(); panel(); } },
    }, svg('title', { text: e ? `${e.name}: α ${signed(e.attraction, 3)}` : f.properties.NAME })));
  }
  host.append(countries);
  // The map's triangles and their eddies.
  if (document.getElementById('show-faces').checked) {
    const faces = svg('g', {});
    const maxCurl = Math.max(1, ...A.triangles.map((t) => Math.abs(t.curl)));
    for (const t of A.triangles) {
      const pts = t.places.map((i) => project(A.economies[i].longitude, A.economies[i].latitude));
      if (Math.max(...pts.map((p) => p[0])) - Math.min(...pts.map((p) => p[0])) > W) continue;  // across the map's edge
      faces.append(svg('path', { class: 'face', d: `M${pts.map((p) => p.map((x) => x.toFixed(1)).join(',')).join('L')}Z` }));
      const s = Math.abs(t.curl) / maxCurl;
      if (s < 0.06) continue;
      const cx = (pts[0][0] + pts[1][0] + pts[2][0]) / 3, cy = (pts[0][1] + pts[1][1] + pts[2][1]) / 3, r = 3 + 9 * Math.sqrt(s);
      const sweep = t.curl > 0 ? 0 : 1;  // anticlockwise on the sphere seen from outside is anticlockwise on the map
      const a0 = 0.3, a1 = 5.3;
      const p0 = [cx + r * Math.cos(a0), cy - r * Math.sin(a0) * (sweep ? -1 : 1)], p1 = [cx + r * Math.cos(a1), cy - r * Math.sin(a1) * (sweep ? -1 : 1)];
      faces.append(svg('path', { class: 'eddy', stroke: PART.local, 'stroke-opacity': 0.35 + 0.6 * s, 'marker-end': 'url(#eddy-head)',
        d: `M${p0[0].toFixed(1)},${p0[1].toFixed(1)}A${r.toFixed(1)},${r.toFixed(1)} 0 1 ${sweep} ${p1[0].toFixed(1)},${p1[1].toFixed(1)}` },
        svg('title', { text: `circulation ${t.curl > 0 ? 'anticlockwise' : 'clockwise'} round ${t.places.map((i) => A.economies[i].iso3).join(', ')}: curl ${money(Math.abs(t.curl))}` })));
    }
    host.append(faces);
  }
  // Arcs.
  const part = document.getElementById('arc-part').value, count = Number(document.getElementById('arc-count').value);
  const weight = (l) => (part === 'driven' ? Math.abs(l[3]) : part === 'circulating' ? Math.abs(l[4] + l[5]) : l[2]);
  const legs = A.legs.filter((l) => state.selected === null || l[0] === state.selected || l[1] === state.selected)
    .sort((a, b) => weight(b) - weight(a)).slice(0, count);
  const maxW = Math.max(1, ...legs.map(weight));
  const arcs = svg('g', {});
  for (const l of legs) {
    const w = weight(l);
    if (w <= 0) continue;
    // The direction drawn: the net flow's, or the part's (a part may run against the net flow).
    const value = part === 'driven' ? l[3] : part === 'circulating' ? l[4] + l[5] : l[2];
    const [from, to] = value >= 0 ? [l[0], l[1]] : [l[1], l[0]];
    const colour = part === 'driven' ? PART.driven : part === 'circulating' ? PART.local : (Math.abs(l[3]) >= Math.abs(l[4] + l[5]) ? PART.driven : PART.local);
    arcs.append(svg('path', { class: 'arc', d: arcPath(A.economies[from], A.economies[to]), stroke: colour,
      'stroke-width': (0.4 + 4.2 * Math.sqrt(w / maxW)).toFixed(2), 'stroke-opacity': (0.35 + 0.5 * Math.sqrt(w / maxW)).toFixed(2), 'marker-end': 'url(#head)' }));
  }
  host.append(svg('defs', {},
    svg('marker', { id: 'head', viewBox: '0 0 10 10', refX: 8, refY: 5, markerWidth: 4, markerHeight: 4, orient: 'auto-start-reverse' },
      svg('path', { d: 'M0,1L9,5L0,9z', fill: 'context-stroke', 'fill-opacity': 0.8 })),
    svg('marker', { id: 'eddy-head', viewBox: '0 0 10 10', refX: 6, refY: 5, markerWidth: 5, markerHeight: 5, orient: 'auto' },
      svg('path', { d: 'M0,1L9,5L0,9z', fill: PART.local }))));
  host.append(arcs);
  const places = svg('g', {});
  for (const e of A.economies) {
    const [x, y] = project(e.longitude, e.latitude);
    places.append(svg('circle', { class: 'place', cx: x.toFixed(1), cy: y.toFixed(1), r: 0.9 }));
  }
  host.append(places);
}

function list(title, rows) {
  return el('div', { class: 'section' }, el('h3', { text: title }),
    el('table', {}, el('tbody', {}, rows.map(([a, b]) => el('tr', {}, el('td', { text: a }), el('td', { class: 'num', text: b }))))));
}

const ranked = (A, order) => A.economies.filter((e) => e.exports + e.imports > RANKED).sort(order).slice(0, 8);

function panel() {
  const box = document.getElementById('map-panel');
  const A = state.answer;
  clear(box);
  if (!A) return;
  const T = A.totals;
  if (state.selected !== null) {
    const i = state.selected, e = A.economies[i];
    const mine = A.legs.filter((l) => l[0] === i || l[1] === i).sort((a, b) => b[2] - a[2]);
    box.append(el('h2', { text: e.name }), el('p', { class: 'lede' }, `${e.continent} · ${A.year}. `,
      el('button', { class: 'primary', onclick: () => { state.selected = null; draw(); panel(); } }, 'the whole world')),
    el('dl', {},
      el('dt', { text: 'exports' }), el('dd', { text: money(e.exports) }),
      el('dt', { text: 'imports' }), el('dd', { text: money(e.imports) }),
      el('dt', { text: 'net inflow (imports − exports)' }), el('dd', { text: money(e.balance) }),
      el('dt', { text: 'attraction α' }), el('dd', { text: signed(e.attraction, 3) }),
      el('dt', { text: 'circulating, of the parts on its links' }), el('dd', { text: percent(e.circulating) })),
    el('p', {}, e.attraction > 0 ? 'Value is drawn to it: the driven part of net trade runs up to it from the economies of lower attraction, along each link in proportion to the link\'s gross trade.' : 'Value is sent from it: the driven part of net trade runs down from it to the economies of higher attraction, along each link in proportion to the link\'s gross trade.'),
    list('Its largest net flows', mine.slice(0, 14).map((l) => {
      const other = A.economies[l[0] === i ? l[1] : l[0]];
      const circ = Math.abs(l[4] + l[5]), parts = Math.abs(l[3]) + circ;
      return [`${l[0] === i ? 'to' : 'from'} ${other.name}`, `${money(l[2])}, ${percent(parts > 0 ? circ / parts : 0, 0)} circ.`];
    })));
    return;
  }
  box.append(
    el('h2', { text: `World trade in goods, ${A.year}` }),
    el('p', { class: 'lede', text: `${A.economies.length} economies, ${money(T.total)} of exports. What two economies send each other both ways is exchange; what is left, one net flow for each pair, is driven or circulates.` }),
    el('h3', { text: 'Exchange and net transfer' }),
    shareBar([{ label: 'reciprocal exchange', value: T.reciprocal, colour: PART.reciprocal }, { label: 'net flows', value: T.net, colour: PART.driven }]),
    el('h3', { text: 'The net flows, by energy', title: 'shares of Σ F²/c, the energy of the net flows with each link conducting its gross trade', style: { marginTop: '14px' } }),
    shareBar([{ label: 'driven', value: T.driven, colour: PART.driven }, { label: 'round triangles', value: T.local, colour: PART.local }, { label: 'round holes', value: T.global, colour: PART.global }]),
    el('p', {}, `Each link conducts in proportion to its gross trade. The driven part runs up a single potential, the attraction α: it is the flow the balances alone would drive through those links, and between two economies it is their gross trade times their difference in α. The rest circulates: round the map's ${A.faces.triangles} triangles of neighbouring economies (local eddies) and round the loops through distant partners that no triangle fills (global circulation). Value round a loop returns to where it started; it is moved by the structure of trade, not by any one economy's balance.`),
    list('Drawing value most (α), of those trading over 20 bn $', ranked(A, (a, b) => b.attraction - a.attraction).map((e) => [e.name, signed(e.attraction, 3)])),
    list('Sending value most (α)', ranked(A, (a, b) => a.attraction - b.attraction).map((e) => [e.name, signed(e.attraction, 3)])),
    list('Circulating most, of the parts on their links', ranked(A, (a, b) => b.circulating - a.circulating).map((e) => [e.name, percent(e.circulating)])),
    el('p', { class: 'note', text: `Decomposed ${A.method === 'sparse' ? 'by conjugate gradients' : 'exactly'} on the complex of ${A.economies.length} places, ${A.legs.length} net flows and ${A.faces.triangles} triangles. Source: ${A.source.citation}` }),
  );
}

async function load(year) {
  const status = document.getElementById('status');
  if (state.cache.has(year)) {
    state.answer = state.cache.get(year);
  } else {
    status.textContent = `decomposing the trade of ${year}…`;
    try {
      state.answer = await api.map(year);
      state.cache.set(year, state.answer);
      status.textContent = '';
    } catch (error) {
      status.textContent = error.message;
      return;
    }
  }
  draw();
  panel();
}

export async function showMap(data) {
  if (!state.geo) {
    const response = await fetch(COUNTRIES);
    state.geo = await response.json();
    const slider = document.getElementById('year');
    const years = data?.trade?.years || [2010, 2024];
    slider.min = years[0];
    slider.max = years[1];
    slider.value = years[1];
    document.getElementById('year-out').textContent = years[1];
    slider.addEventListener('input', () => { document.getElementById('year-out').textContent = slider.value; });
    slider.addEventListener('change', () => load(Number(slider.value)));
    for (const id of ['show-faces', 'arc-part', 'arc-count']) document.getElementById(id).addEventListener('change', draw);
  }
  await load(Number(document.getElementById('year').value));
}
