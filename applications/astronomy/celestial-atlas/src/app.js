// Celestial Atlas: data loading, the clock, the four views and the interface.

import * as THREE from 'three';
import { SPK } from './core/spk.js';
import { scalesFromUnixMs, unixMsFromJdTDB, formatUtc, J2000, DAY } from './core/time.js';
import { DEG, AU_KM, ICRS_TO_ECLIPTIC, ICRS_TO_GALACTIC, apply, sub, norm, toSpherical } from './core/frames.js';
import { BODIES, BY_KEY, PLANETS, orbitalMu } from './core/bodies.js';
import { snapshot } from './core/sky.js';
import { osculating } from './core/kepler.js';
import { loadStars, loadConstellations, listStarCatalogues } from './data/catalog.js';
import { loadMaps } from './render/maps.js';
import { GlobePicker } from './render/globe-picker.js';
import { nearestPlace } from './core/geodesy.js';
import { SolarSystemView } from './views/solar-system.js';
import { NeighbourhoodView, LY_PER_PC, absoluteMagnitude } from './views/neighbourhood.js';
import { SkyChartView } from './views/sky-chart.js';
import { OrbitDiagramView } from './views/orbit-diagram.js';
import { GalaxiesView } from './views/galaxies.js';
import { HubbleDiagram, MassLightDiagram, HUBBLE_CONSTANT } from './views/galaxy-diagrams.js';
import { listGalaxyCatalogues, loadGalaxies } from './data/catalog.js';

const DATA = 'data';
const STORAGE_KEY = 'celestial-atlas/v1';

const SPEEDS = [
  [1, 'real time'], [60, '1 minute per second'], [3600, '1 hour per second'],
  [86400, '1 day per second'], [604800, '1 week per second'],
  [2629800, '1 month per second'], [31557600, '1 year per second'],
];
const RATES = [...SPEEDS.slice().reverse().map(([s, l]) => [-s, `${l}, backwards`]), ...SPEEDS];
const REAL_TIME = SPEEDS.length;
const SIZE_STEPS = [1, 10, 50, 100, 300, 600, 1500, 3000];
const RADIUS_STEPS = [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
const PLACES = [
  { name: 'Tokyo', lat: 35.6812, lon: 139.7671, height: 0.04 },
  { name: 'Beijing', lat: 39.9042, lon: 116.4074, height: 0.05 },
  { name: 'Greenwich', lat: 51.4769, lon: -0.0005, height: 0.05 },
  { name: 'Mauna Kea', lat: 19.8207, lon: -155.4681, height: 4.2 },
  { name: 'Cerro Paranal', lat: -24.6272, lon: -70.4042, height: 2.64 },
  { name: 'Sydney', lat: -33.8688, lon: 151.2093, height: 0.03 },
];
const VIEW_NAMES = { solar: 'Solar System', near: 'Neighbourhood', galaxies: 'Galaxies', sky: 'Sky Chart', orbits: 'Orbits' };
const MLY_PER_MPC = 3.261563777;

const $ = (selector) => document.querySelector(selector);
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => `&#${c.charCodeAt(0)};`);
const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));

// ---- persistence (a per-viewer convenience; the page works without it) ----

function readSaved() {
  try { return JSON.parse(localStorage.getItem(STORAGE_KEY)) || {}; } catch { return {}; }
}

function save() {
  try {
    localStorage.setItem(STORAGE_KEY, JSON.stringify({
      view: state.view, place: state.place, options: state.options,
    }));
  } catch { /* storage unavailable */ }
}

// An address may name a view, a place and an instant, as another application
// links here: #view=sky&lat=35.68&lon=139.77&name=Tokyo&t=2026-10-22T12:00Z.
function readLink() {
  const q = new URLSearchParams(location.hash.slice(1));
  const out = {};
  if (VIEW_NAMES[q.get('view')]) out.view = q.get('view');
  const lat = Number(q.get('lat')), lon = Number(q.get('lon'));
  if (q.has('lat') && q.has('lon') && Math.abs(lat) <= 90 && Math.abs(lon) <= 180) {
    out.place = { name: q.get('name') || '', lat, lon, height: 0, timezone: q.get('tz') || undefined };
  }
  const t = Date.parse(q.get('t') || '');
  if (Number.isFinite(t)) out.ms = t;
  return out;
}

const saved = readSaved();
const linked = readLink();
const state = {
  view: linked.view || (VIEW_NAMES[saved.view] ? saved.view : 'solar'),
  ms: linked.ms ?? Date.now(),
  live: linked.ms === undefined,
  paused: linked.ms !== undefined,
  rate: REAL_TIME,
  place: linked.place || saved.place || PLACES[0],
  selection: {
    solar: { kind: 'body', key: 'earth' },
    near: { kind: 'star', index: null },
    sky: { kind: 'body', key: 'moon' },
    orbits: { kind: 'body', key: 'mars' },
    galaxies: { kind: 'galaxy', index: null },
  },
  options: {
    solar: { size: 5, orbits: true, labels: true, stars: true, milkyWay: true, constellations: false, grid: true },
    near: { catalogue: 'stars-v1', radius: 2, labels: true, grid: true, milkyWay: true },
    sky: { magLimit: 5.5, lines: true, names: true, nameLang: 'la', ecliptic: true, grid: true,
      refraction: true, starNames: true },
    orbits: { range: 'inner', sightLines: true },
    galaxies: { catalogue: 'galaxies-ungc', mode: 'map', relations: 'published', field: false,
      colourBy: 'relation', minMembers: 3, labels: true, grid: true },
  },
};
for (const view of Object.keys(state.options)) Object.assign(state.options[view], saved.options?.[view]);

const observer = () => ({
  name: state.place.name, lat: state.place.lat * DEG, lon: state.place.lon * DEG, height: state.place.height || 0,
  timezone: state.place.timezone,
});

// Coordinates as words: 35.68° N, 139.77° E.
const placeCoords = (lat, lon) => `${Math.abs(lat).toFixed(2)}° ${lat >= 0 ? 'N' : 'S'}, ${Math.abs(lon).toFixed(2)}° ${lon >= 0 ? 'E' : 'W'}`;

// The Earth Atlas at the same place and instant: on the local servers it has a
// port of its own; on a site holding the repository, a path of its own.
function earthAtlasUrl(place, ms) {
  const base = location.port === '8770' ? `${location.protocol}//${location.hostname}:8790/` : '../../geography/earth-atlas/';
  const { lat, lon } = place;
  return `${base}#p=globe&c=${lon.toFixed(1)},${lat.toFixed(1)},2.4&at=${lon.toFixed(3)},${lat.toFixed(3)}`
    + `&t=${new Date(ms).toISOString().slice(0, 16)}Z`;
}

// ---- formatting --------------------------------------------------------------

const pad = (n) => String(n).padStart(2, '0');
const signed = (x, digits) => `${x < 0 ? '−' : ''}${Math.abs(x).toFixed(digits)}`;

function formatRA(rad) {
  const total = Math.round(((rad / DEG / 15 + 24) % 24) * 36000);   // tenths of a second
  const h = Math.floor(total / 36000), m = Math.floor((total % 36000) / 600), s = (total % 600) / 10;
  return `${h}h ${pad(m)}m ${s.toFixed(1).padStart(4, '0')}s`;
}

function formatDec(rad) {
  const total = Math.round(Math.abs(rad / DEG) * 3600);
  const d = Math.floor(total / 3600), m = Math.floor((total % 3600) / 60), s = total % 60;
  return `${rad < 0 ? '−' : '+'}${d}° ${pad(m)}′ ${pad(s)}″`;
}

function formatDistance(km) {
  const au = km / AU_KM;
  if (au < 0.01) return `${Math.round(km).toLocaleString('en')} km`;
  return `${au.toFixed(au < 10 ? 4 : 3)} AU`;
}

function formatDuration(seconds) {
  if (seconds < 60) return `${seconds.toFixed(1)} s`;
  if (seconds < 3600) return `${Math.floor(seconds / 60)} min ${Math.round(seconds % 60)} s`;
  if (seconds < 86400 * 2) return `${Math.floor(seconds / 3600)} h ${Math.round((seconds % 3600) / 60)} min`;
  const days = seconds / 86400;
  return days < 800 ? `${days.toFixed(2)} days` : `${(days / 365.25).toFixed(2)} years`;
}

function compass(azDeg) {
  const points = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'];
  return points[Math.round(azDeg / 22.5) % 16];
}

function eclipticLongitude(v) {
  return toSpherical(apply(ICRS_TO_ECLIPTIC, v)).lon;
}

// ---- information shown in the panel -----------------------------------------

function bodyInfo(snap, key) {
  const b = BY_KEY[key];
  const e = snap.bodies[key];
  const rows = [];
  if (key !== 'earth') {
    rows.push(['Distance from Earth', formatDistance(e.distKm)]);
    rows.push(['Light-time', formatDuration(e.lightTime)]);
  }
  if (key !== 'sun') rows.push(['Distance from Sun', formatDistance(norm(e.helio.p))]);
  if (key !== 'earth') {
    rows.push(['Right ascension', formatRA(e.ra)]);
    rows.push(['Declination', formatDec(e.dec)]);
    const alt = e.altApparent / DEG, az = e.az / DEG;
    rows.push(['Altitude', `${signed(alt, 1)}°${alt < 0 ? ' · below horizon' : ''}`]);
    rows.push(['Azimuth', `${az.toFixed(1)}° ${compass(az)}`]);
  }
  if (key !== 'sun' && key !== 'earth') {
    let d = (eclipticLongitude(e.geo) - eclipticLongitude(snap.bodies.sun.geo)) / DEG;
    d = ((d + 540) % 360) - 180;
    rows.push(['Elongation', `${(e.elongation / DEG).toFixed(1)}° ${d >= 0 ? 'E · evening' : 'W · morning'}`]);
    rows.push(['Illuminated', `${Math.round(e.illuminated * 100)} %`]);
  }
  if (key !== 'sun') {
    let state6;
    if (key === 'earth') state6 = snap.emb;
    else if (b.parent) {
      const p = snap.bodies[b.parent].ssb;
      state6 = { p: sub(e.ssb.p, p.p), v: sub(e.ssb.v, p.v) };
    } else state6 = e.helio;
    const o = osculating(state6.p, state6.v, orbitalMu(b));
    rows.push(['Semi-major axis', formatDistance(o.a)]);
    rows.push(['Eccentricity', o.e.toFixed(4)]);
    rows.push(['Orbital period', formatDuration(o.period)]);
  }
  const kind = key === 'sun' ? 'Star' : key === 'moon' ? 'Natural satellite of the Earth'
    : key === 'pluto' ? 'Dwarf planet' : 'Planet';
  if (key === 'earth') {
    const p = state.place;
    rows.push(['Observer', p.name || placeCoords(p.lat, p.lon)]);
    if (p.name) rows.push(['', placeCoords(p.lat, p.lon)]);
    return { title: b.name, subtitle: kind, rows, actions: [['sky', 'The sky from here']] };
  }
  return { title: b.name, subtitle: kind, rows };
}

function starInfo(stars, i) {
  const gaia = stars.photometry === 'gaia';
  const band = gaia ? 'G' : 'V';
  const rows = [];
  const dist = stars.dist[i];
  if (Number.isFinite(dist)) {
    const digits = dist < 100 ? 2 : 0;
    rows.push(['Distance', `${dist.toFixed(digits)} pc · ${(dist * LY_PER_PC).toFixed(dist < 100 ? 1 : 0)} ly`]);
  }
  rows.push([`Apparent magnitude (${band})`, signed(stars.mag[i], 2)]);
  if (Number.isFinite(dist)) rows.push([`Absolute magnitude (${band})`, signed(absoluteMagnitude(stars.mag[i], dist), 2)]);
  if (stars.spectrum(i)) rows.push(['Spectral type', stars.spectrum(i)]);
  if (Number.isFinite(stars.ci[i])) rows.push([gaia ? 'BP − RP' : 'B − V', signed(stars.ci[i], 2)]);
  const u = [stars.unit[3 * i], stars.unit[3 * i + 1], stars.unit[3 * i + 2]];
  const eq = toSpherical(u);
  rows.push([`Right ascension (${stars.epoch})`, formatRA(eq.lon)]);
  rows.push([`Declination (${stars.epoch})`, formatDec(eq.lat)]);
  const g = toSpherical(apply(ICRS_TO_GALACTIC, u));
  rows.push(['Galactic l, b', `${(g.lon / DEG).toFixed(1)}°, ${(g.lat / DEG).toFixed(1)}°`]);
  const designation = stars.designation(i);
  return {
    title: stars.label(i) || designation || 'Unnamed star',
    subtitle: stars.label(i) && designation && designation !== stars.label(i) ? designation : stars.title,
    rows,
  };
}

function morphology(t) {
  if (!Number.isFinite(t)) return 'Galaxy';
  if (t <= -4) return 'Elliptical galaxy';
  if (t <= -1) return 'Lenticular galaxy';
  if (t <= 8) return 'Spiral galaxy';
  if (t <= 9) return 'Magellanic spiral';
  if (t <= 10) return 'Irregular galaxy';
  return `Galaxy, type ${t}`;
}

function galaxyInfo(cat, i, options) {
  const rows = [];
  if (cat.kind === 'tidal') {
    const rel = cat.relations[options.relations];
    const pub = cat.relations.published, ours = cat.relations.ours;
    const d = cat.dist[i];
    rows.push(['Distance', `${d.toFixed(2)} Mpc`]);
    rows.push(['Light-travel', `${(d * MLY_PER_MPC).toFixed(1)} million years`]);
    if (cat.distMethod[i]) rows.push(['Distance method', cat.distMethod[i]]);
    if (Number.isFinite(cat.catalogueLogL[i])) rows.push(['log L_K (solar)', cat.catalogueLogL[i].toFixed(2)]);
    else if (Number.isFinite(cat.logL[i])) rows.push(['log L_K, inferred', cat.logL[i].toFixed(2)]);
    if (Number.isFinite(cat.vLG[i])) {
      rows.push(['Velocity (Local Group)', `${signed(cat.vLG[i], 0)} km/s`]);
      rows.push(['Peculiar velocity', `${signed(cat.vLG[i] - HUBBLE_CONSTANT * d, 0)} km/s`]);
    }
    if (Number.isFinite(pub.theta1[i])) rows.push(['Θ₁ published', signed(pub.theta1[i], 1)]);
    if (Number.isFinite(ours.theta1[i])) rows.push(['Θ₁ recomputed', signed(ours.theta1[i], 2)]);
    const md = pub.disturber[i];
    const sep = (k) => Math.hypot(cat.exact[3 * k] - cat.exact[3 * i], cat.exact[3 * k + 1] - cat.exact[3 * i + 1],
      cat.exact[3 * k + 2] - cat.exact[3 * i + 2]);
    if (md >= 0) rows.push(['Main disturber', `${cat.names[md]}, ${sep(md).toFixed(2)} Mpc`]);
    else if (cat.publishedDisturberName[i]) {
      rows.push(['Main disturber', `${cat.publishedDisturberName[i]}, outside the catalogue`]);
    }
    const od = ours.disturber[i];
    if (od >= 0 && od !== md) rows.push(['Recomputed disturber', `${cat.names[od]}, ${sep(od).toFixed(2)} Mpc`]);
    const g = rel.groups.groups[rel.groups.groupOf[i]];
    rows.push(['Tidal group', g.size > 1 ? `${cat.names[g.leader]} · ${g.size}` : 'none (field)']);
    const dwarf = cat.dwarfType[i] ? ` · dwarf ${cat.dwarfType[i]}` : '';
    return { title: cat.names[i], subtitle: `${morphology(cat.type[i])}${dwarf}`, rows };
  }
  const gr = cat.groups, j = cat.groupOf[i];
  rows.push(['PGC', String(cat.pgc[i])]);
  rows.push(['Velocity (CMB)', `${cat.vCMB[i].toFixed(0)} km/s`]);
  rows.push(['Distance', `${cat.dist[i].toFixed(1)} h⁻¹ Mpc`]);
  if (Number.isFinite(cat.logL[i])) rows.push(['log L_K (solar)', cat.logL[i].toFixed(2)]);
  if (Number.isFinite(cat.kmag[i])) rows.push(['Ks magnitude', cat.kmag[i].toFixed(2)]);
  if (gr.nHere[j] > 1) {
    rows.push(['Group members', `${gr.nHere[j]}`]);
    if (Number.isFinite(gr.sigmaV[j])) rows.push(['Dispersion σ', `${gr.sigmaV[j].toFixed(0)} km/s`]);
    if (Number.isFinite(gr.r2t[j])) rows.push(['Radius R₂ₜ', `${gr.r2t[j].toFixed(2)} Mpc`]);
    if (gr.mVir[j] > 0) rows.push(['Virial mass', `${(gr.mVir[j]).toPrecision(3)} × 10¹² M☉`]);
    if (gr.mLum[j] > 0) rows.push(['Mass from light', `${(gr.mLum[j]).toPrecision(3)} × 10¹² M☉`]);
  }
  return {
    title: cat.galaxyName(i),
    subtitle: gr.nHere[j] > 1 ? (cat.groupName(j) || `Group ${gr.nest[j]}`) : `${morphology(cat.type[i])} · field`,
    rows,
  };
}

function renderInfo(info) {
  const html = info
    ? `<h2${info.title.length > 18 ? ' class="long"' : ''}>${esc(info.title)}</h2><p class="sub">${esc(info.subtitle)}</p>`
      + `<dl>${info.rows.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('')}</dl>`
      + (info.actions ? `<div class="actions">${info.actions.map(([id, label]) => `<button type="button" class="action" data-action="${id}">${esc(label)}</button>`).join('')}</div>` : '')
    : `<p class="hint">Select a ${state.view === 'galaxies' ? 'galaxy' : 'star'} to see its details.</p>`;
  const target = $('#info');
  if (target.dataset.html !== html) {
    target.innerHTML = html;
    target.dataset.html = html;
  }
}

// ---- application --------------------------------------------------------------

async function main() {
  const progress = $('#loading-bar');
  const [spk, stars, constellations, catalogues, galaxyCatalogues] = await Promise.all([
    SPK.fromUrl(`${DATA}/jpl-de440s/de440s.bsp`, (p) => { progress.style.width = `${(p * 100).toFixed(1)}%`; }),
    loadStars(DATA),
    loadConstellations(DATA),
    listStarCatalogues(DATA),
    listGalaxyCatalogues(DATA),
  ]);
  if (!catalogues.some((c) => c.name === state.options.near.catalogue)) state.options.near.catalogue = 'stars-v1';
  if (!galaxyCatalogues.some((c) => c.name === state.options.galaxies.catalogue)) {
    state.options.galaxies.catalogue = galaxyCatalogues[0]?.name || '';
  }
  const loadedCatalogues = new Map([[stars.name, stars]]);
  const loadedGalaxies = new Map();
  let galaxies = null;
  const minMs = unixMsFromJdTDB(J2000 + spk.coverage.startEt / DAY) + 2 * 86400000;
  const maxMs = unixMsFromJdTDB(J2000 + spk.coverage.endEt / DAY) - 2 * 86400000;

  const stage = $('#stage');
  const renderer = new THREE.WebGLRenderer({ canvas: $('#gl'), antialias: true, logarithmicDepthBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  renderer.setClearColor(0x120d14, 1);

  const views = {
    solar: new SolarSystemView({ renderer, container: stage, stars, constellations,
      onSelect: (key, how) => selectBody(key, how) }),
    near: new NeighbourhoodView({ renderer, container: stage, stars, onSelect: (i) => selectStar(i) }),
    sky: new SkyChartView({ canvas: $('#chart'), stars, constellations }),
    orbits: new OrbitDiagramView({ canvas: $('#chart') }),
    galaxies: new GalaxiesView({ renderer, container: stage, onSelect: (i) => selectGalaxy(i, { focus: true }) }),
    hubble: new HubbleDiagram($('#chart')),
    massLight: new MassLightDiagram($('#chart')),
  };
  // Maps of the bodies and the sky arrive after the first frame; they are optional.
  loadMaps(DATA, renderer).then((maps) => {
    if (!maps) return;
    views.solar.setMaps(maps);
    globe.setImages(maps.earthDay && maps.earthDay.image, maps.earthNight && maps.earthNight.image);
    if (maps.sky) views.near.setSky(maps.sky);
    dirty = true;
  });
  const globe = new GlobePicker((lat, lon) => choosePlace(lat, lon));
  const galaxyMap = () => state.options.galaxies.mode === 'map';
  const is3d = (v) => v === 'solar' || v === 'near' || (v === 'galaxies' && galaxyMap());
  // The view object that draws the current tab.
  const active = () => {
    if (state.view !== 'galaxies') return views[state.view];
    if (galaxyMap()) return views.galaxies;
    return galaxies?.kind === 'groups' ? views.massLight : views.hubble;
  };
  let dirty = true;
  let lastDrawnMs = 0;
  let snap = null;

  function applyOptions() {
    const o = state.options;
    views.solar.setOptions({ sizeScale: SIZE_STEPS[o.solar.size], orbits: o.solar.orbits, labels: o.solar.labels,
      stars: o.solar.stars, milkyWay: o.solar.milkyWay, constellations: o.solar.constellations, grid: o.solar.grid });
    if (views.near.radius !== RADIUS_STEPS[o.near.radius]) views.near.setRadius(RADIUS_STEPS[o.near.radius]);
    views.near.setOptions({ labels: o.near.labels, grid: o.near.grid, milkyWay: o.near.milkyWay });
    if (views.near.stars.name !== o.near.catalogue) useCatalogue(o.near.catalogue);
    views.sky.setOptions(o.sky);
    views.orbits.setOptions(o.orbits);
    const g = o.galaxies;
    views.galaxies.setOptions({ relations: g.relations, field: g.field, colourBy: g.colourBy,
      minMembers: g.minMembers, labels: g.labels, grid: g.grid });
    if (state.view === 'galaxies' && g.catalogue && galaxies?.name !== g.catalogue) useGalaxies(g.catalogue);
    dirty = true;
  }

  // Galaxy catalogues are loaded when the Galaxies view first needs them.
  async function useGalaxies(name) {
    let cat = loadedGalaxies.get(name);
    if (!cat) {
      const entry = galaxyCatalogues.find((c) => c.name === name);
      if (!entry) return;
      toast(`Loading ${entry.title}…`);
      try {
        cat = await loadGalaxies(DATA, entry);
      } catch (error) {
        toast(`Could not load ${name}: ${error.message}`);
        return;
      }
      loadedGalaxies.set(name, cat);
    }
    if (state.options.galaxies.catalogue !== name) return;
    galaxies = cat;
    views.galaxies.setCatalogue(cat);
    state.selection.galaxies = { kind: 'galaxy', index: null };
    views.hubble.selected = views.massLight.selected = null;
    if (state.view === 'galaxies') showView('galaxies');
    dirty = true;
  }

  function selectGalaxy(index, { focus = false } = {}) {
    if (!galaxies || index == null) return;
    state.selection.galaxies = { kind: 'galaxy', index };
    views.galaxies.select(index);
    views.hubble.selected = views.massLight.selected = index;
    if (focus && galaxyMap()) views.galaxies.focusGalaxy(index);
    dirty = true;
    updatePanel(true);
  }

  // Star catalogues for the neighbourhood are loaded when first chosen.
  async function useCatalogue(name) {
    let catalogue = loadedCatalogues.get(name);
    if (!catalogue) {
      toast(`Loading ${catalogues.find((c) => c.name === name)?.title || name}…`);
      try {
        catalogue = await loadStars(DATA, name);
      } catch (error) {
        toast(`Could not load ${name}: ${error.message}`);
        return;
      }
      loadedCatalogues.set(name, catalogue);
    }
    if (state.options.near.catalogue !== name) return;   // another choice was made meanwhile
    views.near.setStars(catalogue);
    state.selection.near = { kind: 'star', index: null };
    updateLegend();
    updatePanel(true);
  }

  function selectBody(key, { focus = false } = {}) {
    state.selection[state.view] = { kind: 'body', key };
    if (state.view === 'solar') {
      views.solar.select(key);
      if (focus) views.solar.focus(key);
    }
    if (state.view === 'sky') views.sky.selectedKey = key;
    if (state.view === 'orbits') views.orbits.selectedKey = key;
    dirty = true;
    updatePanel(true);
  }

  // ---- where the observer stands ----

  // Names of places from the geography product (Natural Earth), read when first needed.
  let gazetteer = null;
  function loadGazetteer() {
    if (!gazetteer) {
      gazetteer = fetch(`${DATA}/products/index.json`).then((r) => r.json()).then((index) => {
        const entry = (index.geography || [])[0];
        return entry ? fetch(`${DATA}/products/${entry.file}`).then((r) => r.json()) : null;
      }).then((g) => {
        if (!g) return null;
        // Natural Earth's short names in English: China, not People's Republic of China.
        const countries = new Map(g.layers.countries.map((c) => [c.id, c.names.default || c.names.en]));
        return g.layers.places.map((p) => ({ name: p.names.default || p.names.en, lat: p.lat, lon: p.lon,
          country: countries.get(p.country) || '', timezone: p.timezone }));
      }).catch(() => null);
    }
    return gazetteer;
  }

  // Nautical time at a longitude (radians): the zone of whole hours from UTC,
  // as the IANA zone Etc/GMT-k (whose sign is reversed by convention).
  const nauticalZone = (lon) => {
    const k = Math.round(lon / DEG / 15);
    return k === 0 ? 'Etc/GMT' : `Etc/GMT${k > 0 ? '-' : '+'}${Math.abs(k)}`;
  };

  // The observer at a place (radians), named after the nearest town within
  // 120 km, and keeping the time zone of the nearest town within 500 km, or
  // else nautical time.
  async function choosePlace(lat, lon) {
    const latDeg = +(lat / DEG).toFixed(4), lonDeg = +(((lon / DEG + 540) % 360) - 180).toFixed(4);
    state.place = { name: '', lat: latDeg, lon: lonDeg, height: 0 };
    placeChanged(true);
    renderOptions();
    dirty = true;
    const places = await loadGazetteer();
    const near = places && nearestPlace(places, lat, lon);
    if (state.place.lat === latDeg && state.place.lon === lonDeg) {
      state.place.timezone = near && near.km < 500 && near.place.timezone ? near.place.timezone : nauticalZone(lon);
      if (near && near.km < 120) state.place.name = `near ${near.place.name}${near.place.country ? `, ${near.place.country}` : ''}`;
      renderOptions();
      updatePanel(true);
    }
    save();
  }

  // Keeps the marks of the observer's place in step with it.
  function placeChanged(centre = false) {
    const place = { lat: state.place.lat * DEG, lon: state.place.lon * DEG };
    views.solar.setObserver(place);
    globe.setPlace(place, centre);
    updatePanel(true);
  }

  function selectStar(index) {
    state.selection[state.view] = { kind: 'star', index };
    if (state.view === 'near') views.near.select(index);
    updatePanel(true);
  }

  // ---- views and options ----

  function showView(view) {
    state.view = view;
    for (const button of document.querySelectorAll('[data-view]')) {
      button.setAttribute('aria-selected', String(button.dataset.view === view));
    }
    $('#gl').hidden = !is3d(view);
    $('#chart').hidden = is3d(view);
    views.solar.setActive(view === 'solar');
    views.near.setActive(view === 'near');
    views.galaxies.setActive(view === 'galaxies' && galaxyMap() && !!galaxies);
    if (view === 'galaxies' && !galaxies && state.options.galaxies.catalogue) useGalaxies(state.options.galaxies.catalogue);
    const sel = state.selection[view];
    if (view === 'solar') { views.solar.select(sel.key); }
    if (view === 'sky') views.sky.selectedKey = sel.key;
    if (view === 'orbits') views.orbits.selectedKey = sel.key;
    if (view === 'near') views.near.select(sel.index);
    renderOptions();
    resize();
    dirty = true;
    updatePanel(true);
    save();
  }

  function check(view, key, label) {
    return `<label class="check"><input type="checkbox" data-view-option="${view}" data-key="${key}"`
      + `${state.options[view][key] ? ' checked' : ''}><span>${label}</span></label>`;
  }

  function renderOptions() {
    const v = state.view, o = state.options[v];
    let html = '';
    if (v === 'solar') {
      html = `<h3>Display</h3>
        <label class="range"><span>Body size</span>
          <input type="range" min="0" max="${SIZE_STEPS.length - 1}" step="1" value="${o.size}" data-view-option="solar" data-key="size"></label>
        <p class="legend" id="legend"></p>
        ${check('solar', 'orbits', 'Orbits')}${check('solar', 'labels', 'Names')}
        ${check('solar', 'grid', 'Distance rings')}${check('solar', 'stars', 'Stars')}${check('solar', 'milkyWay', 'Milky Way')}
        ${check('solar', 'constellations', 'Constellation figures')}
        <h3>Go to</h3>
        <div class="chips">${BODIES.map((b) => `<button type="button" data-focus="${b.key}">`
          + `<span class="dot" style="background:${b.colour}"></span>${b.name}</button>`).join('')}</div>`;
    } else if (v === 'near') {
      html = `<h3>Catalogue</h3>
        <select data-view-option="near" data-key="catalogue" class="wide">
          ${catalogues.map((c) => `<option value="${esc(c.name)}"${c.name === o.catalogue ? ' selected' : ''}>`
            + `${esc(c.title)} · ${Number(c.count).toLocaleString('en')} stars</option>`).join('')}</select>
        <p class="legend">More catalogues: <code>python data/acquisition/fetch.py --list</code></p>
        <h3>Display</h3>
        <label class="range"><span>Radius</span>
          <input type="range" min="0" max="${RADIUS_STEPS.length - 1}" step="1" value="${o.radius}" data-view-option="near" data-key="radius"></label>
        <p class="legend" id="legend"></p>
        ${check('near', 'labels', 'Names')}${check('near', 'grid', 'Rings and axes')}${check('near', 'milkyWay', 'Milky Way')}`;
    } else if (v === 'sky') {
      const custom = !PLACES.some((p) => p.name === state.place.name);
      html = `<h3>Observer</h3>
        <canvas id="globe" class="globe-picker" aria-label="Globe: click a place to stand there"></canvas>
        <select id="place">${PLACES.map((p) => `<option${p.name === state.place.name ? ' selected' : ''}>${p.name}</option>`).join('')}
          <option value="custom"${custom ? ' selected' : ''}>${custom ? esc(state.place.name) : 'Custom'}</option></select>
        <div class="coords">
          <label><span>Latitude</span><input id="lat" type="number" step="0.0001" min="-90" max="90" value="${state.place.lat}"></label>
          <label><span>Longitude</span><input id="lon" type="number" step="0.0001" min="-180" max="180" value="${state.place.lon}"></label>
        </div>
        <div class="observer-links"><button type="button" class="link" id="locate">Use my location</button>
          <a class="link" id="weather-here" href="${esc(earthAtlasUrl(state.place, state.ms))}" target="_blank" rel="noopener">The weather here ↗</a></div>
        <h3>Display</h3>
        <label class="range"><span>Faintest stars</span>
          <input type="range" min="3" max="7" step="0.5" value="${o.magLimit}" data-view-option="sky" data-key="magLimit"></label>
        <p class="legend" id="legend"></p>
        ${check('sky', 'lines', 'Constellation figures')}${check('sky', 'names', 'Constellation names')}
        <label class="select"><span>Names in</span><select data-view-option="sky" data-key="nameLang">
          ${[['la', 'Latin'], ['en', 'English'], ['zh', '中文'], ['ja', '日本語']].map(([k, l]) =>
            `<option value="${k}"${o.nameLang === k ? ' selected' : ''}>${l}</option>`).join('')}</select></label>
        ${check('sky', 'starNames', 'Bright star names')}${check('sky', 'ecliptic', 'Ecliptic')}
        ${check('sky', 'grid', 'Altitude grid')}${check('sky', 'refraction', 'Refraction')}`;
    } else if (v === 'galaxies') {
      html = galaxyOptions(o);
    } else {
      html = `<h3>Range</h3>
        <div class="segmented">
          <button type="button" data-range="inner" aria-pressed="${o.range === 'inner'}">Inner planets</button>
          <button type="button" data-range="all" aria-pressed="${o.range === 'all'}">All planets</button>
        </div>
        ${check('orbits', 'sightLines', 'Lines of sight from the Earth')}
        <h3>Configuration</h3><table class="config" id="config"></table>`;
    }
    $('#options').innerHTML = html;
    if (v === 'sky') globe.attach($('#globe'));
    updateLegend();
  }

  function galaxyOptions(o) {
    if (!galaxyCatalogues.length) {
      return `<h3>Catalogue</h3><p class="legend">No galaxy catalogue yet. From the repository root:<br>
        <code>python data/acquisition/fetch.py ungc</code><br>
        <code>python data/products/galaxy_catalog.py --source ungc</code></p>`;
    }
    const tidal = galaxies?.kind === 'tidal';
    let html = `<h3>Catalogue</h3>
      <select data-view-option="galaxies" data-key="catalogue" class="wide">
        ${galaxyCatalogues.map((c) => `<option value="${esc(c.name)}"${c.name === o.catalogue ? ' selected' : ''}>`
          + `${esc(c.title)} · ${Number(c.count).toLocaleString('en')}</option>`).join('')}</select>
      <div class="segmented spaced">
        <button type="button" data-mode="map" aria-pressed="${o.mode === 'map'}">3D map</button>
        <button type="button" data-mode="diagram" aria-pressed="${o.mode === 'diagram'}">${tidal ? 'Hubble flow' : 'Mass and light'}</button>
      </div>
      <p class="legend" id="legend"></p>
      <h3>Relations</h3>`;
    if (!galaxies) return `${html}<p class="legend">Loading…</p>`;
    if (tidal) {
      const a = galaxies.agreement;
      html += `<label class="select"><span>Main disturbers</span><select data-view-option="galaxies" data-key="relations">
          <option value="published"${o.relations === 'published' ? ' selected' : ''}>Published</option>
          <option value="ours"${o.relations === 'ours' ? ' selected' : ''}>Recomputed here</option></select></label>
        ${check('galaxies', 'field', 'Relations of field galaxies (Θ₁ ≤ 0)')}
        <p class="legend">Recomputed main disturbers agree with the published ones for
          ${(100 * a.agree / a.comparable).toFixed(1)} % of ${a.comparable} galaxies.</p>`;
    } else {
      html += `<label class="range"><span>Link members of groups with at least ${o.minMembers} galaxies</span>
        <input type="range" min="2" max="30" step="1" value="${o.minMembers}" data-view-option="galaxies" data-key="minMembers"></label>`;
    }
    html += `<h3>Display</h3>
      <label class="select"><span>Colour</span><select data-view-option="galaxies" data-key="colourBy">
        <option value="relation"${o.colourBy === 'relation' ? ' selected' : ''}>${tidal ? 'Tidal index' : 'Group size'}</option>
        <option value="type"${o.colourBy === 'type' ? ' selected' : ''}>Morphological type</option></select></label>
      ${check('galaxies', 'labels', 'Names')}${check('galaxies', 'grid', 'Supergalactic plane')}
      <h3>${tidal ? 'Largest tidal groups' : 'Richest groups'}</h3><div class="chips">`;
    if (tidal) {
      const rel = galaxies.relations[o.relations];
      html += rel.groups.groups.filter((g) => g.size >= 3).slice(0, 10).map((g) =>
        `<button type="button" data-galaxy="${g.leader}">${esc(galaxies.names[g.leader])} · ${g.size}</button>`).join('');
    } else {
      const gr = galaxies.groups;
      const ranked = [...gr.names.keys()].sort((x, y) => gr.nHere[y] - gr.nHere[x]).slice(0, 12);
      html += ranked.map((j) => `<button type="button" data-group="${j}">${esc(gr.names.get(j))} · ${gr.nHere[j]}</button>`).join('');
    }
    return `${html}</div>`;
  }

  function updateLegend() {
    const el = $('#legend');
    if (!el) return;
    if (state.view === 'solar') el.textContent = views.solar.legend();
    if (state.view === 'near') el.textContent = views.near.legend();
    if (state.view === 'sky') el.textContent = `Stars to magnitude ${state.options.sky.magLimit}`;
    if (state.view === 'galaxies') el.textContent = views.galaxies.legend();
  }

  $('#options').addEventListener('input', (event) => {
    const t = event.target;
    if (!t.dataset.viewOption) return;
    const value = t.type === 'checkbox' ? t.checked : t.type === 'range' ? Number(t.value) : t.value;
    state.options[t.dataset.viewOption][t.dataset.key] = value;
    applyOptions();
    if (t.dataset.viewOption === 'galaxies' && t.dataset.key === 'relations') renderOptions();
    if (t.dataset.viewOption === 'galaxies' && t.dataset.key === 'minMembers') {
      t.closest('label').querySelector('span').textContent = `Link members of groups with at least ${value} galaxies`;
    }
    updateLegend();
    updatePanel(true);
    save();
  });

  $('#info').addEventListener('click', (event) => {
    const action = event.target.closest('[data-action]');
    if (action && action.dataset.action === 'sky') showView('sky');
  });

  $('#options').addEventListener('click', (event) => {
    // The link to the Earth Atlas carries the instant shown when it is followed.
    if (event.target.id === 'weather-here') event.target.href = earthAtlasUrl(state.place, state.ms);
    const focus = event.target.closest('[data-focus]');
    if (focus) selectBody(focus.dataset.focus, { focus: true });
    const mode = event.target.closest('[data-mode]');
    if (mode) {
      state.options.galaxies.mode = mode.dataset.mode;
      showView('galaxies');
    }
    const galaxy = event.target.closest('[data-galaxy]');
    if (galaxy) selectGalaxy(Number(galaxy.dataset.galaxy), { focus: true });
    const group = event.target.closest('[data-group]');
    if (group && galaxies?.kind === 'groups') {
      const j = Number(group.dataset.group);
      if (galaxies.groups.brightest[j] >= 0) selectGalaxy(galaxies.groups.brightest[j]);
      if (galaxyMap()) views.galaxies.focusGroupCentre(j);
    }
    const range = event.target.closest('[data-range]');
    if (range) {
      state.options.orbits.range = range.dataset.range;
      applyOptions();
      renderOptions();
      save();
    }
    if (event.target.id === 'locate' && navigator.geolocation) {
      navigator.geolocation.getCurrentPosition((pos) => {
        state.place = { name: 'My location', lat: +pos.coords.latitude.toFixed(4),
          lon: +pos.coords.longitude.toFixed(4), height: (pos.coords.altitude || 0) / 1000,
          timezone: Intl.DateTimeFormat().resolvedOptions().timeZone };
        placeChanged(true);
        renderOptions(); dirty = true; save();
      }, () => toast('Location is not available.'));
    }
  });

  $('#options').addEventListener('change', (event) => {
    const t = event.target;
    if (t.id === 'place' && t.value !== 'custom') {
      state.place = { ...PLACES.find((p) => p.name === t.value) };
      placeChanged();
      renderOptions();
    } else if (t.id === 'lat' || t.id === 'lon') {
      const lat = clamp(Number($('#lat').value), -90, 90), lon = clamp(Number($('#lon').value), -180, 180);
      if (Number.isFinite(lat) && Number.isFinite(lon)) {
        choosePlace(lat * DEG, lon * DEG);
      }
    } else return;
    dirty = true;
    save();
  });

  for (const button of document.querySelectorAll('[data-view]')) {
    button.addEventListener('click', () => showView(button.dataset.view));
  }

  // ---- pointer interaction ----

  let down = null;
  $('#gl').addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
  $('#gl').addEventListener('pointerup', (e) => {
    if (!down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 5) return;
    const rect = $('#gl').getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1,
      -((e.clientY - rect.top) / rect.height) * 2 + 1);
    if (state.view === 'solar') {
      const key = views.solar.pick(ndc);
      // On the focused Earth, a click chooses where to stand.
      const place = key === 'earth' && views.solar.focusKey === 'earth' ? views.solar.pickPlace(ndc) : null;
      if (place) choosePlace(place.lat, place.lon);
      if (key) selectBody(key);
    } else if (state.view === 'galaxies') {
      const index = views.galaxies.pick(ndc);
      if (index != null) selectGalaxy(index);
    } else {
      const index = views.near.pick(ndc);
      if (index != null) selectStar(index);
    }
  });

  const chart = $('#chart');
  const chartPoint = (e) => {
    const rect = chart.getBoundingClientRect();
    return [e.clientX - rect.left, e.clientY - rect.top];
  };
  const hitId = (hit) => (hit ? `${hit.kind}:${hit.key ?? hit.group ?? hit.index}` : '');
  chart.addEventListener('pointermove', (e) => {
    const view = active();
    const hit = view.hitTest(...chartPoint(e));
    if (hitId(hit) !== hitId(view.hover)) {
      view.hover = hit;
      chart.style.cursor = hit ? 'pointer' : 'default';
      dirty = true;
    }
  });
  chart.addEventListener('pointerleave', () => { active().hover = null; dirty = true; });
  chart.addEventListener('click', (e) => {
    const hit = active().hitTest(...chartPoint(e));
    if (!hit) return;
    if (hit.kind === 'body') selectBody(hit.key);
    else if (hit.kind === 'galaxy' || hit.kind === 'group') selectGalaxy(hit.index);
    else selectStar(hit.index);
  });

  // ---- clock ----

  function running() { return state.live || !state.paused; }

  function setRate(index) {
    state.rate = clamp(index, 0, RATES.length - 1);
    state.live = false;
    state.paused = false;
    updateTimebar(true);
  }

  $('#now').addEventListener('click', () => {
    state.live = true; state.paused = false; state.rate = REAL_TIME;
    updateTimebar(true);
  });
  $('#toggle').addEventListener('click', () => {
    if (running()) { state.live = false; state.paused = true; } else state.paused = false;
    updateTimebar(true);
  });
  $('#faster').addEventListener('click', () => setRate(state.rate + 1));
  $('#slower').addEventListener('click', () => setRate(state.rate - 1));
  $('#date').addEventListener('change', (e) => {
    const m = /^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2})(?::(\d{2}))?/.exec(e.target.value);
    if (!m) return;
    const ms = Date.UTC(+m[1], +m[2] - 1, +m[3], +m[4], +m[5], +(m[6] || 0));
    state.ms = clamp(ms, minMs, maxMs);
    state.live = false;
    state.paused = true;
    dirty = true;
    updateTimebar(true);
  });
  window.addEventListener('keydown', (e) => {
    if (e.target.closest('input, select, textarea')) return;
    if (e.key === ' ') { e.preventDefault(); $('#toggle').click(); }
    if (e.key === 'n') $('#now').click();
    if (e.key === ']') $('#faster').click();
    if (e.key === '[') $('#slower').click();
  });

  let lastClock = '';
  function updateTimebar(force = false) {
    const utc = `${formatUtc(state.ms)} UTC`;
    if (!force && utc === lastClock) return;
    lastClock = utc;
    $('#time-utc').textContent = utc;
    $('#time-local').textContent = new Date(state.ms).toLocaleString('en-GB', {
      weekday: 'short', day: 'numeric', month: 'short', year: 'numeric',
      hour: '2-digit', minute: '2-digit', second: '2-digit', timeZoneName: 'short',
    });
    $('#rate').textContent = state.live ? 'live' : state.paused ? 'paused' : RATES[state.rate][1];
    $('#rate').classList.toggle('live', state.live);
    $('#toggle').innerHTML = running()
      ? '<svg viewBox="0 0 24 24"><path d="M6.5 5h3.6v14H6.5zM13.9 5h3.6v14h-3.6z"/></svg>'
      : '<svg viewBox="0 0 24 24"><path d="M7 4.8v14.4L19 12z"/></svg>';
    $('#toggle').setAttribute('aria-label', running() ? 'Pause' : 'Play');
    const input = $('#date');
    if (document.activeElement !== input) input.value = new Date(state.ms).toISOString().slice(0, 19);
  }

  let lastPanel = 0;
  function updatePanel(force = false) {
    const now = performance.now();
    if (!snap || (!force && now - lastPanel < 250)) return;
    lastPanel = now;
    const sel = state.selection[state.view];
    if (state.view === 'galaxies') {
      renderInfo(galaxies && sel.index != null ? galaxyInfo(galaxies, sel.index, state.options.galaxies) : null);
      return;
    }
    const starSource = state.view === 'near' ? views.near.stars : stars;
    renderInfo(sel.kind === 'body' ? bodyInfo(snap, sel.key) : sel.index != null ? starInfo(starSource, sel.index) : null);
    if (state.view === 'orbits') renderConfiguration();
  }

  function renderConfiguration() {
    const table = $('#config');
    if (!table) return;
    const sunLon = eclipticLongitude(snap.bodies.sun.geo);
    const rows = PLANETS.filter((b) => b.key !== 'earth').map((b) => {
      const e = snap.bodies[b.key];
      let d = (eclipticLongitude(e.geo) - sunLon) / DEG;
      d = ((d + 540) % 360) - 180;
      return `<tr><th><span class="dot" style="background:${b.colour}"></span>${b.name}</th>`
        + `<td>${(norm(e.helio.p) / AU_KM).toFixed(2)}</td><td>${(e.distKm / AU_KM).toFixed(2)}</td>`
        + `<td>${Math.abs(d).toFixed(0)}° ${d >= 0 ? 'E' : 'W'}</td></tr>`;
    }).join('');
    const html = '<thead><tr><th></th><th>Sun, AU</th><th>Earth, AU</th><th>Elong.</th></tr></thead>'
      + `<tbody>${rows}</tbody>`;
    if (table.dataset.html !== html) { table.innerHTML = html; table.dataset.html = html; }
  }

  let toastTimer = 0;
  function toast(text) {
    const el = $('#toast');
    el.textContent = text;
    el.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => { el.hidden = true; }, 4000);
  }

  // ---- layout ----

  function resize() {
    const { width, height } = stage.getBoundingClientRect();
    if (!width || !height) return;
    renderer.setSize(width, height, false);
    views.solar.resize(width, height);
    views.near.resize(width, height);
    views.sky.resize(width, height);
    views.orbits.resize(width, height);
    views.galaxies.resize(width, height);
    views.hubble.resize(width, height);
    views.massLight.resize(width, height);
    dirty = true;
  }
  new ResizeObserver(resize).observe(stage);

  // ---- frame loop ----

  // The point under the Sun, for the globe: the Sun's right ascension and
  // declination of date, less the sidereal time at Greenwich.
  let lastGlobeSun = -Infinity, lastGlobeMs = NaN;
  function updateGlobeSun(now) {
    if (now - lastGlobeSun < 1500 && Math.abs(state.ms - lastGlobeMs) < 600000) return;
    lastGlobeSun = now; lastGlobeMs = state.ms;
    const d = toSpherical(apply(snap.precession, snap.bodies.sun.geo));
    globe.setSun({ lat: d.lat, lon: ((d.lon - snap.gmst) % (2 * Math.PI) + 3 * Math.PI) % (2 * Math.PI) - Math.PI });
  }

  let last = performance.now();
  function frame(now) {
    const dt = Math.min(0.25, (now - last) / 1000);
    last = now;
    if (state.live) state.ms = Date.now();
    else if (!state.paused) state.ms += RATES[state.rate][0] * dt * 1000;
    if (state.ms < minMs || state.ms > maxMs) {
      state.ms = clamp(state.ms, minMs, maxMs);
      state.live = false;
      state.paused = true;
      toast('The DE440 ephemeris covers 1849-12-26 to 2150-01-22.');
      updateTimebar(true);
    }
    snap = snapshot(spk, scalesFromUnixMs(state.ms), observer());
    const view = active();
    if (state.view === 'galaxies') {
      if (!galaxies) {
        if (galaxyMap()) renderer.clear();
      } else if (galaxyMap()) {
        view.update(snap, dt);
        view.render();
      } else if (dirty) {
        view.draw(galaxies, galaxies.relations?.[state.options.galaxies.relations]);
        dirty = false;
      }
    } else if (is3d(state.view)) {
      view.update(snap, dt);
      view.render();
    } else {
      if (state.view === 'sky') updateGlobeSun(now);
      const animating = !state.live && !state.paused;
      if (dirty || animating || Math.abs(state.ms - lastDrawnMs) >= 1000) {
        view.draw(snap);
        lastDrawnMs = state.ms;
        dirty = false;
      }
    }
    updateTimebar();
    updatePanel();
    requestAnimationFrame(frame);
  }

  applyOptions();
  placeChanged(true);
  // A place given by a link, without its time zone: the zone of the nearest town, or nautical time.
  if (!state.place.timezone) {
    const { lat, lon } = state.place;
    loadGazetteer().then((places) => {
      if (state.place.lat !== lat || state.place.lon !== lon) return;
      const near = places && nearestPlace(places, lat * DEG, lon * DEG);
      state.place.timezone = near && near.km < 500 && near.place.timezone ? near.place.timezone : nauticalZone(lon * DEG);
      dirty = true;
    });
  }
  showView(state.view);
  $('#loading').hidden = true;
  // For inspection from the browser's console and for automated checks.
  window.celestial = { state, views, selectBody, showView, applyOptions };
  requestAnimationFrame(frame);
}

main().catch((error) => {
  console.error(error);
  $('#loading-text').innerHTML = `<strong>The atlas could not start.</strong><br>${esc(error.message)}<br><br>`
    + 'Fetch the data from the repository root with<br>'
    + '<code>python data/acquisition/fetch.py</code><br><code>python data/products/star_catalog.py</code>';
  $('#loading').classList.add('failed');
});
