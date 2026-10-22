// The Earth Atlas: the atmosphere, the ocean, the land and its people on one
// globe, each from its own dataset, at one instant.
//
// The application composes layers (src/layers/) over a base map, all drawn by
// the stage (src/render/stage.js), and keeps them at the time shown. Every
// number it shows is computed by the provisional core (src/core/) from data
// products of the laboratory's data library (data/products/).

import { Stage } from './render/stage.js';
import { Labels } from './render/labels.js';
import { Library } from './data/products.js';
import { buildCatalogue, DOMAINS, INDICATORS, modelRank } from './layers/catalogue.js';
import { BaseLayer } from './layers/base.js';
import { SurfaceLayer } from './layers/surface.js';
import { FlowLayer } from './layers/flow.js';
import { LinesLayer } from './layers/lines.js';
import { EventsLayer } from './layers/events.js';
import { TracksLayer } from './layers/tracks.js';
import { PlatesLayer, MOTION_PRODUCT } from './layers/plates.js';
import { StationsLayer } from './layers/stations.js';
import { MechanismsLayer, KIND_COLOURS } from './layers/mechanisms.js';
import { mechanism } from './core/moment.js';
import { stereonet } from './ui/stereonet.js';
import { arcDeg } from './core/cyclones.js';
import { cellIndex, column as crustColumn } from './core/crust.js';
import { crustDiagram } from './ui/crust.js';
import { SectionFigure } from './ui/section.js';
import { sectionLine } from './core/section.js';
import { SegmentLayer, polylineSegments } from './render/segments.js';
import { SeriesFigure } from './ui/series.js';
import { Requests } from './data/requests.js';
import { ImageryLayer, DetailLayer } from './layers/imagery.js';
import { LabModel, modelInputs, LAB_MODELS } from './layers/labmodel.js';
import { setClouds } from './render/surface.js';
import { Timeline, formatUTC, formatShortUTC, formatDay, formatOffset } from './ui/timeline.js';
import { renderLegend } from './ui/legend.js';
import { renderPlace } from './ui/place.js';
import { buildPanel } from './ui/panel.js';
import { toVector } from './core/sphere.js';
import { subsolarPoint, sunAltAz, daylight, localMeanSolarTime } from './core/sun.js';
import { windFromDirection, compassPoint } from './core/derived.js';
import { flowToward } from './core/ocean.js';
import { displayUnit, toShown, format } from './core/units.js';
import { areaMean, ClassRaster } from './core/grid.js';
import { koppen, koppenGrid, meanTemperature } from './core/climate.js';
import { climateDiagram } from './ui/climate.js';
import { joinCountries, valueForYear, latestWellCoveredYear } from './core/join.js';
import { OceanBeneath } from './ui/ocean.js';

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const HOUR = 3600000, DAY = 86400000;

const state = {
  projection: 'globe',
  surface: null, compare: null, flow: null, lines: 'none', events: null, indicator: null, tracks: [], plates: false,
  stations: false, stationFrame: 'plate', station: null,
  mechanisms: false, mechanismMode: 'year', quake: null,
  style: 'relief', night: true, lights: true, clouds: false,
  coast: true, borders: true, rivers: true, graticule: false, labels: true, reliefOnData: true,
  lang: 'en',
  time: Date.now(), live: true, playing: false,
  pinned: null,
};

let library, catalogue, stage, base, surface, flow, lines, events, labels, timeline, geography, elevation, tracks;
let imagery, detail, labModel, detailSource = null, ocean, plates, stations, series, mechanisms, crust = null, section, sectionLayer;
const asked = new Requests();       // data requested on demand for the place card, other than the ocean's
const lab = { input: null, running: null, progress: 0, runs: {}, start: null, against: '', error: '' };
const MODEL_COLOURS = ['#CA8DBE', '#f6c47e', '#9fb6e8', '#8fd1b4', '#e8a0a0'];
const loaded = new Map();          // field products in memory, by name
const classMaps = new Map();       // class products in memory, by name
const MONTH_NAMES = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'];
const indicators = new Map();      // indicator products, by name, with their joins
let hover = null;
let dirty = { labels: true, legend: true, place: true };
let lastPlace = 0;

// ---- loading ------------------------------------------------------------------

function progress(fraction, text) {
  $('loading-bar').style.width = `${Math.round(fraction * 100)}%`;
  if (text) $('loading-text').textContent = text;
}

function fail(html) {
  const el = $('loading');
  el.classList.add('failed');
  $('loading-text').innerHTML = html;
}

const MISSING_DATA = `The data products are not built yet. From the repository root, run:<br><br>
<code>python data/acquisition/fetch.py natural-earth-50m etopo-2022 gfs oisst ww3</code><br>
<code>python data/products/geography.py</code><br>
<code>python data/products/field_grid.py --source etopo-2022</code><br>
<code>python data/products/field_grid.py --source gfs</code><br><br>
and see the README for the other datasets.`;

async function start() {
  try {
    progress(0.05, 'Reading the data library…');
    library = await new Library('data/products').open();
    if (!library.list('geography').length) throw new Error('no geography');
    geography = await library.geography();
  } catch (e) {
    console.error(e);
    fail(MISSING_DATA);
    return;
  }
  progress(0.2, 'Drawing the map…');
  stage = new Stage($('gl'));
  base = new BaseLayer(stage, geography);
  // The plates and their motion, a field computed here, before the catalogue is built.
  plates = new PlatesLayer(stage, library);
  try { await plates.load(); } catch (e) { console.error('plates', e); }
  catalogue = buildCatalogue(library);
  restore();
  chooseDefaults();

  progress(0.3, 'Reading the relief…');
  // The finest relief in the library.
  const relief = [...library.fieldMeta].filter(([, m]) => m.variables.some((v) => v.standard_name === 'altitude'))
    .sort((a, b) => b[1].grid.nlon - a[1].grid.nlon)[0];
  if (relief) {
    elevation = await library.field(relief[0]);
    loaded.set(relief[0], elevation);
  }
  base.painted = state.style === 'plain' ? 'plain' : 'relief';
  const painting = base.paint(elevation, base.painted);
  imagery = new ImageryLayer(stage, library);
  detail = new DetailLayer(stage);
  labModel = new LabModel(library);
  const inputs = modelInputs(library);
  if (inputs) lab.input = `${inputs.model} at the start`;

  surface = new SurfaceLayer(stage, library, geography.meta.layers.countries);
  flow = new FlowLayer(stage, library);
  lines = new LinesLayer(stage, library);
  events = new EventsLayer(stage, library);
  tracks = new TracksLayer(stage, library);
  labels = new Labels($('labels'));
  ocean = new OceanBeneath($('figure'), { onChange: () => { dirty.place = true; } });
  stations = new StationsLayer(stage, library, plates);
  stations.setMode(state.stationFrame);
  mechanisms = new MechanismsLayer(stage, library);
  mechanisms.setMode(state.mechanismMode);
  // The section figure, and its line on the map while it is open.
  sectionLayer = new SegmentLayer({ color: '#f6c47e', width: 2, opacity: 0.95, order: 6 });
  stage.overlay.add(sectionLayer.mesh);
  section = new SectionFigure($('section'), {
    quakes: () => mechanisms.all(),
    elevation: (lon, lat) => (elevation ? elevation.sample(elevation.meta.variables[0].name, 0, lon, lat) : NaN),
    crust: (lon, lat) => crustAt(lon, lat),
    crossings: (points) => plates.crossings(points),
  }, { onChange: drawSectionLine });
  series = new SeriesFigure($('series'));

  progress(0.5, 'Reading the forecast…');
  await Promise.all([
    selectSurface(state.surface, false), selectFlow(state.flow, false),
    selectLines(state.lines, false), selectEvents(state.events, false), selectTracks(state.tracks, false),
  ]);
  await painting;
  progress(0.95);

  timeline = new Timeline($('timeline'), (t) => { state.live = false; state.playing = false; setTime(t); updatePlay(); remember(); });
  bind();
  updateTimeline();
  applyBase('style');
  applyBase('lights');
  if (state.clouds) applyBase('clouds');
  plates.show(state.plates);
  if (state.stations) await stations.show(true);
  if (state.mechanisms) await mechanisms.show(true);
  setProjection(state.projection, true);
  setTime(state.live ? Date.now() : state.time);
  renderPanel();
  $('loading').hidden = true;
  // For inspection from the browser's console and for automated checks.
  window.atlas = { state, stage, base, surface, flow, lines, events, library, catalogue: () => catalogue, labels, imagery, detail, lab, runLab, ocean, plates,
    select: { surface: (id) => selectSurface(id), compare: (id) => selectCompare(id), flow: (id) => selectFlow(id), lines: (id) => selectLines(id) } };
  requestAnimationFrame(loop);

  // The rest of the data, for the place card.
  library.crust().then((c) => { crust = c; dirty.place = true; }).catch((e) => console.error('crust', e));
  for (const [name] of library.fieldMeta) track(name);
  for (const e of library.list('indicators')) loadIndicator(e.name);
  const present = presentClasses();
  if (present) library.classes(present.product).then((r) => { classMaps.set(present.product, r); dirty.place = true; });
  // Our own classification of the climate normals, offered beside the published maps.
  const normals = catalogue.normals[0];
  if (normals && present) track(normals.product).then((p) => addOurClasses(normals, p, present));
}

// The published map of climate classes for the present period (1991-2020).
function presentClasses() {
  return catalogue.surfaces.find((d) => d.kind === 'classes' && !d.computed && !d.scenario && d.period[1] === 2020) || null;
}

function addOurClasses(normals, p, present) {
  const months = [...Array(12).keys()], v = normals.variables;
  const values = koppenGrid(months.map((i) => p.slice(v.tmax, i)), months.map((i) => p.slice(v.tmin, i)), months.map((i) => p.slice(v.ppt, i)));
  const ref = library.classMeta.get(present.product);
  const period = normals.period.split('/').map((d) => Number(d.slice(0, 4)));
  library.addClasses(new ClassRaster({
    name: `${ref.series}-laboratory`, title: `${ref.title.split(',')[0]}, classified here from ${normals.model} normals`,
    series: ref.series, period, scenario: null, computed: true, model: `ours, from ${normals.model}`, domain: ref.domain,
    grid: p.grid, classes: ref.classes, no_value: 0,
    source: { credit: `classified here from ${p.meta.source.credit}` },
  }, values));
  catalogue = buildCatalogue(library);
  renderPanel();
}

function track(name) {
  if (loaded.has(name)) return Promise.resolve(loaded.get(name));
  return library.field(name).then((p) => { loaded.set(name, p); dirty.place = true; return p; });
}

function loadIndicator(name) {
  if (indicators.has(name)) return indicators.get(name);
  const promise = library.indicator(name).then((p) => {
    const join = joinCountries(geography.meta.layers.countries, p.economies);
    const entry = { product: p, join, latest: latestWellCoveredYear(p.economies, p.years) };
    indicators.set(name, entry);
    dirty.place = true;
    return entry;
  });
  indicators.set(name, promise);
  return promise;
}

function chooseDefaults() {
  const ids = (list) => list.map((x) => x.id);
  if (!state.surface || (state.surface !== 'none' && !ids(catalogue.surfaces).includes(state.surface))) {
    state.surface = (catalogue.surfaces.find((s) => s.key === 'speed:wind|10 m') || catalogue.surfaces[0] || { id: 'none' }).id;
  }
  if (state.compare && !ids(catalogue.surfaces).includes(state.compare)) state.compare = null;
  if (!state.flow || (state.flow !== 'none' && !ids(catalogue.flows).includes(state.flow))) {
    state.flow = (catalogue.flows.find((f) => f.key === 'wind|10 m') || { id: 'none' }).id;
  }
  if (state.lines !== 'none' && !ids(catalogue.lines).includes(state.lines)) state.lines = 'none';
  if (state.events && !ids(catalogue.events).includes(state.events)) state.events = null;
  state.tracks = state.tracks.filter((id) => ids(catalogue.tracks).includes(id));
  const first = catalogue.surfaces.find((s) => s.kind === 'indicator');
  if (!state.indicator && first) state.indicator = first.id;
}

// ---- selections ---------------------------------------------------------------

const find = (list, id) => list.find((x) => x.id === id) || null;

async function selectSurface(id, render = true) {
  state.surface = id;
  const def = find(catalogue.surfaces, id);
  if (def && def.kind === 'indicator') state.indicator = id;
  if (def && def.kind === 'field') track(def.product);
  let compare = find(catalogue.surfaces, state.compare);
  if (!compare || !def || compare.key !== def.key || compare.id === def.id) { state.compare = null; compare = null; }
  await surface.select(def, compare, state.time);
  if (def && def.kind === 'indicator') loadIndicator(def.product);
  surface.update(state.time);
  afterChange(render);
}

async function selectCompare(id) {
  state.compare = id;
  await selectSurface(state.surface);
}

async function selectFlow(id, render = true) {
  state.flow = id;
  const def = find(catalogue.flows, id);
  if (def) track(def.product);
  await flow.select(def);
  flow.update(state.time);
  afterChange(render);
}

async function selectLines(id, render = true) {
  state.lines = id;
  await lines.select(find(catalogue.lines, id));
  lines.update(state.time, true);
  afterChange(render);
}

// The colour of a model, the same in the place card, the panel and the tracks.
const modelColour = (model) => MODEL_COLOURS[modelRank(model) % MODEL_COLOURS.length];

// Cyclone tracks of the chosen sources; genesis over the sea is judged by the relief.
async function selectTracks(list, render = true) {
  state.tracks = list;
  const name = elevation && elevation.meta.variables[0].name;
  const isOcean = elevation ? (lon, lat) => elevation.sample(name, 0, lon, lat) < 0 : null;
  const defs = list.map((id) => find(catalogue.tracks, id)).filter(Boolean)
    .map((d) => ({ ...d, colour: modelColour(d.model) }));
  await tracks.select(defs, isOcean);
  tracks.update(state.time, true);
  afterChange(render);
}

async function selectEvents(id, render = true) {
  state.events = id;
  await events.select(find(catalogue.events, id));
  events.update(state.time, true);
  afterChange(render);
}

function afterChange(render) {
  dirty = { labels: true, legend: true, place: true };
  if (render) { renderPanel(); updateTimeline(); }
  remember();
}

function renderPanel() {
  const def = find(catalogue.surfaces, state.surface);
  const years = def && def.kind === 'indicator' && surface.year ? { year: surface.year } : null;
  const extra = {
    imagery: imagery ? imagery.dayImagery() : [],
    lights: Boolean(imagery && imagery.entry('night-lights')),
    clouds: catalogue.clouds.length > 0,
    lab: lab.input ? lab : null,
    modelColour,
    plates: Boolean(plates && plates.plates),
    stations: library.list('stations').length > 0,
    mechanisms: library.list('events').some((e) => e.kind === 'mechanism'),
    ai: [...library.fieldMeta.values()].filter((m) => m.verification && m.verification.start).map((m) => ({
      label: m.model, short: m.model.replace(/ \(laboratory\)$/, '').replace(/ hindcast$/, ''),
      start: Date.parse(m.verification.start), leads: m.verification.leads,
    })),
  };
  buildPanel($('layers'), catalogue, state, years, extra, {
    surface: (id) => selectSurface(id),
    compare: (id) => selectCompare(id),
    runLab: (kind) => runLab(kind),
    flow: (id) => selectFlow(id),
    lines: (id) => selectLines(id),
    events: (id) => selectEvents(id),
    track: (id, on) => selectTracks(on ? [...state.tracks, id] : state.tracks.filter((x) => x !== id)),
    plates: (on) => { state.plates = on; plates.show(on); dirty.labels = dirty.legend = true; remember(); },
    stations: async (on) => {
      state.stations = on;
      if (!on) state.station = null;
      await stations.show(on);
      dirty.legend = dirty.place = true; renderPanel(); remember();
    },
    mechanisms: async (on) => {
      state.mechanisms = on;
      if (!on) state.quake = null;
      await mechanisms.show(on);
      mechanisms.layout(state.time);
      dirty.legend = dirty.place = true; renderPanel(); remember();
    },
    mechanismMode: (mode) => { state.mechanismMode = mode; mechanisms.setMode(mode); mechanisms.layout(state.time); dirty.legend = true; renderPanel(); remember(); },
    stationFrame: (mode) => { state.stationFrame = mode; stations.setMode(mode); dirty.legend = dirty.place = true; renderPanel(); remember(); },
    year: (step) => {
      const p = surface.product;
      if (!p || !p.years) return;
      const y = Math.max(p.years[0], Math.min(p.years[p.years.length - 1], surface.year + step));
      surface.setYear(y);
      dirty.legend = dirty.place = true;
      renderPanel();
    },
    base: (key, value) => { state[key] = value; applyBase(key); renderPanel(); remember(); },
  });
}

function applyBase(key) {
  const u = stage.material.uniforms;
  u.uNight.value = state.night ? 1 : 0;
  u.uReliefOnData.value = state.reliefOnData ? 0.55 : 0;
  if (key === 'style') {
    const imageryStyle = state.style !== 'relief' && state.style !== 'plain';
    if (!imageryStyle || base.painted !== 'relief') {
      const want = state.style === 'plain' ? 'plain' : 'relief';
      if (base.painted !== want) { base.paint(elevation, want); base.painted = want; }
    }
    imagery.setStyle(state.style).then(() => imagery.detailSource()).then((source) => {
      detailSource = source;
      if (source) detail.request(source); else detail.clear();
      dirty.legend = true;
    });
  }
  if (key === 'lights') imagery.setLights(state.lights);
  if (key === 'clouds') {
    const c = catalogue.clouds[0];
    if (state.clouds && c) track(c.product).then(() => setTime(state.time));
    else u.uCloudOn.value = 0;
    dirty.legend = true;
  }
  base.show(state);
  $('labels').hidden = !state.labels;
  dirty.labels = true;
}

function setProjection(mode, initial = false) {
  state.projection = mode;
  document.querySelectorAll('[data-projection]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.projection === mode)));
  if (!initial || mode !== stage.mode) stage.setProjection(mode);
  base.rebuild();
  base.show(state);
  if (flow.product) flow.reset();
  lines.update(state.time, true);
  events.update(state.time, true);
  tracks.rebuild();
  tracks.update(state.time, true);
  if (plates.built) plates.rebuild();
  drawSectionLine();
  stations.stale = true;
  mechanisms.stale = true;
  dirty.labels = true;
  remember();
}

// ---- time ---------------------------------------------------------------------

function timeDomain() {
  const domains = [], steps = [];
  const add = (product) => {
    if (product && product.times && product.times.length > 1) {
      domains.push([product.times[0], product.times[product.times.length - 1]]);
      if (!steps.length) steps.push(...product.times);
    }
  };
  if (surface.def && surface.def.kind === 'field' && surface.product && !surface.product.climatology) add(surface.product);
  add(flow.product);
  add(lines.product);
  if (events.product) domains.push(events.window());
  if (tracks.sources.length && tracks.window()) domains.push(tracks.window());
  const now = Date.now();
  if (!domains.length) return { span: [now - DAY, now + DAY], steps: [], coverage: [] };
  const t0 = Math.min(...domains.map((d) => d[0])), t1 = Math.max(...domains.map((d) => d[1]));
  const pad = (t1 - t0) * 0.02;
  return { span: [Math.min(t0, now) - pad, Math.max(t1, now) + pad], steps, coverage: domains };
}

function updateTimeline() {
  const d = timeDomain();
  timeline.set({ ...d, time: state.time, now: Date.now() });
  stepTimes = d.steps;
}
let stepTimes = [];

function setTime(t) {
  state.time = t;
  surface.update(t);
  if (state.clouds && catalogue.clouds.length) {
    const c = catalogue.clouds[0], p = loaded.get(c.product);
    if (p) setClouds(stage.material, surface.textures, p, c.variable, p.bracket(t, c.variable), true);
  }
  flow.update(t);
  lines.update(t);
  events.update(t);
  tracks.update(t);
  ocean.figure.setNow(t);
  const sub = subsolarPoint(t);
  stage.material.uniforms.uSun.value.fromArray(toVector(sub.lon, sub.lat));
  timeline && timeline.set({ time: t, now: Date.now() });
  $('time-utc').textContent = formatUTC(t);
  const note = $('time-note');
  note.classList.toggle('live', state.live);
  const ahead = t - Date.now();
  note.textContent = state.live ? 'now' : Math.abs(ahead) < 60000 ? 'the present' : `${formatOffset(-ahead).replace(' the time shown', ' now').replace('before now', 'ago').replace('after now', 'from now')}`;
  dirty.legend = dirty.place = dirty.labels = true;
}

function step(direction) {
  state.live = false; state.playing = false;
  const t = state.time;
  const list = stepTimes.length ? stepTimes : [];
  let next = direction > 0 ? list.find((s) => s > t + 1000) : [...list].reverse().find((s) => s < t - 1000);
  if (next === undefined) next = t + direction * HOUR;
  setTime(next);
  updatePlay();
  remember();
}

function updatePlay() {
  const b = $('play');
  b.innerHTML = state.playing
    ? '<svg viewBox="0 0 24 24"><path d="M6.5 5h4v14h-4zM13.5 5h4v14h-4z"/></svg>'
    : '<svg viewBox="0 0 24 24"><path d="M7 4.5v15l12.5-7.5z"/></svg>';
  b.setAttribute('aria-label', state.playing ? 'Pause' : 'Play');
}

// ---- interaction ----------------------------------------------------------------

function bind() {
  document.querySelectorAll('[data-projection]').forEach((b) => b.addEventListener('click', () => setProjection(b.dataset.projection)));
  $('now').addEventListener('click', () => { state.live = true; state.playing = false; setTime(Date.now()); updatePlay(); remember(); });
  $('back').addEventListener('click', () => step(-1));
  $('forward').addEventListener('click', () => step(1));
  $('play').addEventListener('click', () => {
    state.playing = !state.playing;
    state.live = false;
    const [a, b] = timeline.span;
    if (state.playing && state.time >= b - HOUR) setTime(a);
    updatePlay();
  });
  updatePlay();
  $('data-link').addEventListener('click', () => {
    renderAbout();
    $('about').showModal();
    $('about').scrollTop = 0;
  });
  $('place').addEventListener('click', (e) => {
    if (e.target.closest('[data-action="unpin"]')) { state.pinned = null; dirty.place = dirty.labels = true; remember(); }
    if (e.target.closest('[data-action="ocean-column"]') && state.pinned) ocean.ask(state.pinned[0], state.pinned[1]);
    if (e.target.closest('[data-action="ocean-show"]') && state.pinned) {
      series.hide();
      ocean.show(state.pinned[0], state.pinned[1], state.time, seaFloorAt(state.pinned[0], state.pinned[1]));
    }
    if (e.target.closest('[data-action="station-series"]') && state.station) showStationSeries(state.station);
    if (e.target.closest('[data-action="section"]') && state.pinned) openSection(state.pinned[0], state.pinned[1]);
  });
  stage.on('hover', (place, e) => {
    if (!place || !e) { hover = null; $('tooltip').hidden = true; return; }
    const r = $('gl').getBoundingClientRect();
    hover = { lon: place[0], lat: place[1], x: e.clientX - r.left, y: e.clientY - r.top };
    showTooltip();
  });
  stage.on('click', (place) => {
    if (!place) return;
    // An earthquake's mechanism, or else a GPS station, under the pointer is chosen with its place.
    const s = stage.screen(place[0], place[1]);
    const quake = s && mechanisms.shown ? mechanisms.pick(s[0], s[1]) : null;
    const st = s && !quake && stations.shown ? stations.pick(s[0], s[1]) : null;
    state.quake = quake ? quake.row[mechanisms.f.name] : null;
    state.station = st ? st.id : null;
    if (quake) place = [quake.lon, quake.lat];
    if (st) place = [st.lon, st.lat];
    state.pinned = [place[0], place[1]];
    dirty.place = dirty.labels = true;
    remember();
  });
  stage.on('change', () => { dirty.labels = true; });
  window.addEventListener('keydown', (e) => {
    if (e.target.closest('input, select, textarea')) return;
    if (e.key === 'ArrowRight') step(1);
    else if (e.key === 'ArrowLeft') step(-1);
    else if (e.key === ' ') { e.preventDefault(); $('play').click(); }
  });
}

function placeName(lon, lat) {
  const i = base.countryAt(lon, lat);
  const lang = state.lang;
  const pick = (names) => names[lang] || names.en || names.default;
  if (i >= 0) return { kind: 'country', index: i, country: geography.meta.layers.countries[i], name: pick(geography.meta.layers.countries[i].names) };
  const sea = base.seaAt(lon, lat);
  if (sea) return { kind: 'sea', name: pick(sea.names), sea };
  return { kind: 'none', name: '' };
}

function showTooltip() {
  const tip = $('tooltip');
  if (!hover) { tip.hidden = true; return; }
  const quake = events.product ? events.pick(hover.x, hover.y) : null;
  let html;
  if (quake && events.product.kind === 'disaster') {
    html = `<span class="v">${esc(quake.name || '')}</span> · ${esc(quake.alert || '')}<br><span class="n">${esc(quake.severity || events.product.types[quake.type] || '')}</span><br>${formatShortUTC(quake.time_ms)}`;
  } else if (quake) {
    html = `<span class="v">M ${quake.mag.toFixed(1)}</span> · ${format(quake.depth_km)} km deep<br><span class="n">${esc(quake.place || '')}</span><br>${formatShortUTC(quake.time_ms)}`;
  } else {
    const name = placeName(hover.lon, hover.lat).name;
    const value = surfaceValueText(hover.lon, hover.lat, hover);
    if (!value && !name) { tip.hidden = true; return; }
    html = `${value ? `<span class="v">${value}</span>` : ''}${value && name ? ' · ' : ''}${name ? `<span class="n">${esc(name)}</span>` : ''}`;
  }
  tip.innerHTML = html;
  tip.style.left = `${hover.x}px`;
  tip.style.top = `${hover.y}px`;
  tip.hidden = false;
}

function surfaceValueText(lon, lat) {
  const def = surface.def;
  if (!def || !surface.product) return '';
  if (def.kind === 'classes') {
    const c = surface.classesAt(lon, lat);
    if (!c || !c.now) return '';
    const then = c.then ? ` <span class="n">${c.then.code === c.now.code ? 'unchanged' : `was ${esc(c.then.code)}`}</span>` : '';
    return `${esc(c.now.code)} <span class="n">${esc(c.now.name)}</span>${then}`;
  }
  if (def.kind === 'indicator') {
    const i = base.countryAt(lon, lat);
    const v = i >= 0 && surface.values && surface.values[i];
    return v ? `${format(v.value)} ${esc(def.units || '')}${v.year !== surface.year ? ` (${v.year})` : ''}` : '';
  }
  const x = surface.sample(lon, lat, state.time);
  if (!Number.isFinite(x)) return '';
  return `${format(toShown(x, surface.unit))} ${esc(surface.unit.unit)}`;
}

// ---- the place card ---------------------------------------------------------------

function coords(lon, lat) {
  return `${Math.abs(lat).toFixed(2)}° ${lat >= 0 ? 'N' : 'S'}, ${Math.abs(lon).toFixed(2)}° ${lon >= 0 ? 'E' : 'W'}`;
}

// The Celestial Atlas's Sky Chart at a place and instant: on the local servers
// it has a port of its own; on a site holding the repository, a path of its own.
function skyUrl(lon, lat, t, name) {
  const base = location.port === '8790' ? `${location.protocol}//${location.hostname}:8770/` : '../../astronomy/celestial-atlas/';
  const q = new URLSearchParams({ view: 'sky', lat: lat.toFixed(4), lon: lon.toFixed(4), t: new Date(t).toISOString().slice(0, 16) + 'Z' });
  if (name) q.set('name', name);
  return `${base}#${q}`;
}

function hhmm(hours) {
  const m = Math.round(hours * 60) % 1440;
  return `${String(Math.floor(m / 60)).padStart(2, '0')}:${String(m % 60).padStart(2, '0')}`;
}

function valueAt(probe, product, lon, lat, t) {
  const [a, b] = probe.variables;
  if (probe.kind === 'vector') {
    const u = product.sample(a, t, lon, lat), v = product.sample(b, t, lon, lat);
    return { u, v, speed: Math.hypot(u, v) };
  }
  if (probe.kind === 'direction') return { value: product.sampleDirection(a, t, lon, lat) };
  return { value: product.sample(a, t, lon, lat) };
}

function probeUnit(probe) {
  const unit = displayUnit(probe.units);
  return probe.spread ? { ...unit, b: 0 } : unit;
}

function formatProbe(probe, r) {
  const unit = probeUnit(probe);
  if (probe.kind === 'vector') {
    if (!Number.isFinite(r.speed)) return null;
    // Wind is named by where it comes from; currents by where they flow.
    const angle = probe.toward ? flowToward(r.u, r.v) : windFromDirection(r.u, r.v);
    const dir = Number.isFinite(angle) ? ` <span class="dir">${probe.toward ? 'toward' : 'from'} ${compassPoint(angle)}</span>` : '';
    return `${format(toShown(r.speed, unit))} ${esc(unit.unit)}${dir}`;
  }
  if (!Number.isFinite(r.value)) return null;
  if (probe.kind === 'direction') return `from ${compassPoint(r.value)} <span class="dir">${Math.round(r.value)}°</span>`;
  if (probe.percent) return `${format(r.value * 100)} %`;
  const shown = toShown(r.value, unit);
  // Rates of precipitation below 0.05 mm/h are below what the model resolves.
  if (unit.unit === 'mm/h' && shown < 0.05) return `0 ${esc(unit.unit)}`;
  return `${format(shown)} ${esc(unit.unit)}`;
}

function timeNote(product, t) {
  if (!product.times.length) return '';
  const b = product.bracket(t);
  if (b.offset && Math.abs(b.offset) > HOUR) return ` <span class="dir">${formatDay(b.shown).replace(/ \d{4}$/, '')}</span>`;
  return '';
}

function placeCard() {
  const t = state.time;
  const sub = subsolarPoint(t);
  if (!state.pinned) {
    const rows = [['Sun overhead at', coords(sub.lon, sub.lat)]];
    const def = surface.def;
    if (def && def.kind === 'field' && surface.product) {
      const mean = (p, name) => {
        const b = p.bracket(t, name);
        const a = areaMean(p.slice(name, Math.max(0, b.i0)), p.grid);
        return b.w > 0 ? a + (areaMean(p.slice(name, b.i1), p.grid) - a) * b.w : a;
      };
      const label = `${def.label.charAt(0).toLowerCase()}${def.label.slice(1)}`;
      if (def.variables.length === 1 && surface.compare) {
        // A difference of models: both means are over the whole globe, so their
        // difference is the mean of the difference (in units without offset).
        const c = surface.compare;
        const m = mean(surface.product, def.variables[0]) - mean(c.product, c.def.variables[0]);
        if (Number.isFinite(m)) rows.push([`Mean difference in ${label}`, `${format(toShown(m, surface.unit))} ${esc(surface.unit.unit)}`]);
      } else if (def.variables.length === 1) {
        const m = mean(surface.product, def.variables[0]);
        if (Number.isFinite(m)) rows.push([`Mean ${label}`, `${format(toShown(m, surface.unit))} ${esc(surface.unit.unit)}`]);
      }
    } else if (def && def.kind === 'indicator' && surface.product) {
      const w = surface.product.aggregates.WLD;
      const r = w && valueForYear(w, surface.product.years, surface.year, 0);
      if (r) rows.push([`World, ${r.year}`, `${format(r.value)} ${esc(def.units || '')}`]);
    }
    if (events.product && events.product.kind === 'disaster') {
      const f = events.index;
      rows.push(['Disaster alerts shown', String(events.shown.length)]);
      for (const e of events.shown.filter((x) => x[f.alert] !== 'Green')) rows.push([e[f.alert], `${esc(e[f.name] || '')}`]);
    } else if (events.product) {
      const f = events.index;
      const shown = events.shown;
      if (shown.length) {
        const big = shown.reduce((m, e) => (e[f.mag] > m[f.mag] ? e : m), shown[0]);
        rows.push(['Earthquakes shown', String(shown.length)]);
        rows.push(['Largest', `M ${big[f.mag].toFixed(1)}`]);
        // The place, which can be long, on a line of its own across the card.
        if (big[f.place]) rows.push(['', esc(big[f.place]), 'gloss']);
      }
    }
    return { title: 'The Earth', sub: formatUTC(t), rows, hint: 'Choose a place on the map.' };
  }

  const [lon, lat] = state.pinned;
  const where = placeName(lon, lat);
  const sun = sunAltAz(lon, lat, sub);
  const rows = [
    ['Sun’s altitude', `${sun.altitude.toFixed(0).replace('-', '−')}° <span class="dir">${daylight(sun.altitude)}</span>`],
    ['Local mean solar time', hhmm(localMeanSolarTime(t, lon))],
    ['Sky', `<a class="link" href="${esc(skyUrl(lon, lat, t, where.name))}" target="_blank" rel="noopener">The sky here ↗</a>`],
  ];
  const sections = new Map();
  const add = (domain, term, html, kind = '') => {
    const key = DOMAINS[domain] || domain;
    if (!sections.has(key)) sections.set(key, []);
    sections.get(key).push([term, html, kind]);
  };
  const series = [];
  // Probes of the same quantity from several models are shown together.
  const groups = new Map();
  for (const probe of catalogue.probes) {
    const k = probe.key || `${probe.product}/${probe.variables.join('/')}`;
    if (!groups.has(k)) groups.set(k, []);
    groups.get(k).push(probe);
  }
  const models = new Set();
  const cells = (probe) => { const g = library.fieldMeta.get(probe.product).grid; return g.nlon * g.nlat; };
  for (const list of groups.values()) {
    // Models in their order; for relief, only the finest product.
    if (list[0].id === 'elevation') list.splice(0, list.length, list.reduce((a, b) => (cells(b) > cells(a) ? b : a)));
    list.sort((a, b) => modelRank(a.model) - modelRank(b.model));
    const rows = [];
    for (const probe of list) {
      const p = loaded.get(probe.product);
      if (!p) continue;
      if (probe.id === 'elevation') {
        const z = p.sample(probe.variables[0], t, lon, lat);
        if (Number.isFinite(z)) add('geography', z >= 0 ? 'Elevation' : 'Depth', `${format(Math.abs(z))} m`);
        continue;
      }
      const r = valueAt(probe, p, lon, lat, t);
      const text = formatProbe(probe, r);
      if (text) rows.push({ probe, p, r, text });
    }
    if (!rows.length) continue;
    const [main, ...others] = rows;
    add(main.probe.domain, main.probe.label, main.text + timeNote(main.p, t) + (others.length ? ` <span class="dir">${esc(main.probe.model)}</span>` : ''));
    // Other models of the same quantity, beneath it.
    for (const o of others) add(main.probe.domain, o.probe.model, o.text + timeNote(o.p, t), 'model');
    const lines = [];
    for (const row of rows) {
      const { probe, p, r } = row;
      if (p.times.length < 2 || probe.kind === 'direction') continue;
      const unit = probeUnit(probe);
      let values = p.series(probe.variables[0], lon, lat);
      if (probe.kind === 'vector') {
        const v = p.series(probe.variables[1], lon, lat);
        values = values.map((u, i) => Math.hypot(u, v[i]));
      }
      const shown = values.map((v) => (probe.percent ? v * 100 : toShown(v, unit)));
      const cursor = { time: t, value: probe.kind === 'vector' ? toShown(r.speed, unit) : probe.percent ? r.value * 100 : toShown(r.value, unit) };
      models.add(probe.model);
      lines.push({ times: p.times, values: shown, cursor, model: probe.model });
    }
    if (lines.length) {
      series.push({ label: main.probe.kind === 'vector' ? `${main.probe.label} speed` : main.probe.label, now: main.text.replace(/ <span.*$/, ''), lines });
    }
  }
  const climate = climateAt(lon, lat);
  for (const [term, html, kind] of climate.rows) add('climate', term, html, kind);
  // The solid Earth: an earthquake chosen, a GPS station chosen, the plate, its motion, and the nearest plate boundary.
  const quake = state.quake && mechanisms.find(state.quake);
  let quakeFigure = null;
  if (quake) {
    const f = mechanisms.f, r = quake.row, mech = mechanism(quake.elements);
    const when = new Date(quake.t).toISOString().slice(0, 16).replace('T', ' ');
    add('geophysics', 'Earthquake', `M ${quake.mw.toFixed(1)} <span class="dir">${when} UTC</span>`);
    add('geophysics', '', `${esc(r[f.region])}, ${r[f.depth_km].toFixed(0)} km deep (centroid)${r[f.quick] ? '; a quick solution' : ''}`, 'gloss');
    const words = { thrust: 'thrust: one side pushed up over the other, where the crust is squeezed', normal: 'normal: one side dropped, where the crust is stretched', 'strike-slip': 'strike-slip: the two sides slid past each other', oblique: 'oblique: thrust or normal and sliding at once' };
    add('geophysics', 'Faulting', `<span style="color:${KIND_COLOURS[mech.kind]}">${mech.kind}</span>`);
    add('geophysics', '', words[mech.kind], 'gloss');
    const plane = (p) => `${p.strike.toFixed(0)}° / ${p.dip.toFixed(0)}° / ${p.rake.toFixed(0)}°`;
    add('geophysics', 'Nodal planes', `${plane(mech.planes[0])} <span class="dir">and</span> ${plane(mech.planes[1])}`);
    add('geophysics', '', 'strike / dip / rake: one plane is the fault, the other auxiliary; the waves alone cannot tell which', 'gloss');
    const axis = (a) => `${a.plunge.toFixed(0)}° toward ${a.azimuth.toFixed(0)}°`;
    add('geophysics', 'Pressure axis', axis(mech.axes.P));
    add('geophysics', 'Tension axis', axis(mech.axes.T));
    quakeFigure = { svg: stereonet(mech, quake.m, KIND_COLOURS[mech.kind]), caption: 'Lower hemisphere, equal area, north up: coloured where the first motion of P waves is a push; the two nodal planes; P and T axes. Global CMT.' };
  }
  // The crust beneath, from CRUST1.0: its layers to the Moho.
  let crustFigure = null;
  if (crust) {
    const c = crustAt(lon, lat);
    if (c) {
      add('geophysics', 'The crust beneath', `${c.crust.toFixed(1)} km thick <span class="dir">${c.oceanic ? 'oceanic' : 'continental'}</span>`);
      add('geophysics', 'Moho', `${c.moho.toFixed(1)} km below sea level`);
      if (c.sediments > 0.05) add('geophysics', 'Sediments', `${c.sediments.toFixed(1)} km`);
      crustFigure = { svg: crustDiagram(c), caption: 'CRUST1.0 (Laske et al. 2013), the 1° cell containing this place: layers to scale with their P-wave speeds; the Moho, where P waves quicken to about 8 km/s, divides crust from mantle.' };
    }
  }
  // A GPS station chosen.
  const st = state.station && stations.find(state.station);
  if (st) {
    for (const [term, html, kind] of stations.rows(st)) add('geophysics', term, html, kind);
    const entry = asked.get(`station:${st.id}`);
    if (entry && !entry.product) {
      add('geophysics', '', entry.error ? esc(entry.error) : `asking the Nevada Geodetic Laboratory… ${Math.round((Date.now() - entry.asked) / 1000)} s`, 'gloss');
    }
  }
  for (const [term, html, kind] of plates.rows(lon, lat)) add('geophysics', term, html, kind);
  if (plates.plates && library.list('events').some((e) => e.kind === 'mechanism')) {
    add('geophysics', 'Through the Earth', '<button type="button" class="action" data-action="section">A section through here</button>');
  }
  // The ocean beneath: currents and temperature by depth, asked for on demand.
  for (const [term, html, kind] of ocean.rows(lon, lat, t, seaFloorAt(lon, lat))) add('oceanography', term, html, kind);
  const palette = [...models].sort((a, b) => modelRank(a) - modelRank(b));
  for (const s2 of series) for (const l of s2.lines) l.colour = MODEL_COLOURS[palette.indexOf(l.model) % MODEL_COLOURS.length];
  const order = ['meteorology', 'oceanography', 'climate', 'geography', 'geophysics'].map((d) => DOMAINS[d]);
  const ordered = [...sections].sort((a, b) => (order.indexOf(a[0]) + 1 || 99) - (order.indexOf(b[0]) + 1 || 99));
  if (where.kind === 'country') {
    const rowsC = [];
    for (const e of library.list('indicators')) {
      const entry = indicators.get(e.name);
      if (!entry || entry instanceof Promise) continue;
      const code = entry.join.codes[where.index];
      if (!code) continue;
      const def = find(catalogue.surfaces, `indicator:${e.name}`);
      const year = def && def.id === state.surface ? surface.year : entry.latest;
      const r = valueForYear(entry.product.economies[code], entry.product.years, year, 4);
      const [label, unit] = INDICATORS[e.code] || [e.quantity || e.title, e.unit || ''];
      if (r) rowsC.push([`${label}, ${r.year}`, `${format(r.value)}${unit ? ` <span class="dir">${esc(unit)}</span>` : ''}`]);
    }
    if (rowsC.length) ordered.push([`${where.name}, World Bank`, rowsC]);
  }
  const sub2 = where.kind === 'country' ? [coords(lon, lat), where.country.subregion].filter(Boolean).join(' · ') : coords(lon, lat);
  return {
    title: where.name || (lat < -60 ? 'Antarctica' : 'Open water'),
    sub: sub2,
    pinned: true,
    rows,
    sections: ordered.map(([heading, r]) => ({ heading, rows: r, figure: heading === DOMAINS.climate ? climate.figure : heading === DOMAINS.geophysics ? (quakeFigure || crustFigure) : null })),
    series,
    models: palette.length > 1 ? palette.map((m, i) => ({ model: m, colour: MODEL_COLOURS[i % MODEL_COLOURS.length] })) : null,
  };
}

// The crust at a place from CRUST1.0, or null.
function crustAt(lon, lat) {
  if (!crust) return null;
  const g = crust.meta.grid, cell = cellIndex(g, lon, lat);
  const q = Object.fromEntries(crust.meta.quantities.map((x) => [x.name, x]));
  const c = crustColumn((name, k) => crust.values[q[name].offsets[k] / 2 + cell] * q[name].scale);
  return Number.isFinite(c.moho) ? c : null;
}

// A section through the Earth at a place: across the nearest plate boundary.
async function openSection(lon, lat) {
  await mechanisms.load();
  ocean.figure.hide(); series.hide();
  section.show(plates.sectionAt(lon, lat));
}

// The section's line on the map, with the edges of its band, while it is open.
function drawSectionLine() {
  if (!section || !section.open) { sectionLayer.clear(); return; }
  const proj = stage.projector();
  const lift = (lon, lat, out) => proj.project(lon, lat, out, 1.0022);
  const pts = sectionLine(section.lon, section.lat, section.az, section.half, 20);
  const coords = new Float32Array(pts.flatMap((p) => [p.lon, p.lat]));
  const seg = polylineSegments(coords, [[0, pts.length]], lift, proj.breaks, { subdivide: 0.5 });
  sectionLayer.set(seg.a, seg.b, seg.info);
}

// The daily positions of a GPS station, asked for on demand, then drawn.
function showStationSeries(id) {
  const st = stations.find(id);
  if (!st) return;
  const key = `station:${id}`;
  const open = (entry) => {
    dirty.place = true;
    if (entry.product && state.station === id) {
      ocean.figure.hide();
      // Earthquakes near enough to have moved the station measurably, if the catalogue is at hand.
      const quakes = mechanisms.near(st.lon, st.lat, (a, b, c, d) => arcDeg(a, b, c, d) * 111.195);
      Object.assign(quakes, { minMw: 6, radius: 'near enough to move the ground here measurably (within 40 km at magnitude 6, 400 km at 8)' });
      series.show(entry.product, st, quakes);
    }
  };
  const known = asked.get(key);
  if (known && known.product) { open(known); return; }
  asked.ask(key, 'ngl-station', { station: id }, open);
}

// The depth of the sea floor beneath a place (m, positive down), from the relief.
function seaFloorAt(lon, lat) {
  if (!elevation) return NaN;
  const z = elevation.sample(elevation.meta.variables[0].name, 0, lon, lat);
  return Number.isFinite(z) && z < 0 ? -z : NaN;
}

// The climate of a place: its class in the present and in the period shown,
// and its normal year, as figures and as a climate diagram.
function climateAt(lon, lat) {
  const rows = [];
  let figure = null;
  // The code in the row; the class's name beneath, across the card.
  const classRow = (term, c) => rows.push([term, `<span class="code">${esc(c.code)}</span>`], ['', esc(c.name), 'gloss']);
  const present = presentClasses();
  const now = present && classMaps.get(present.product);
  const c0 = now && now.classAt(lon, lat);
  if (c0) classRow(`Climate class, ${present.model}`, c0);
  if (surface.def && surface.def.kind === 'classes' && surface.product && surface.product !== now) {
    const c = surface.product.classAt(lon, lat);
    if (c) classRow(`Climate class, ${surface.def.model}`, c);
  }
  const n = catalogue.normals[0], p = n && loaded.get(n.product);
  if (p) {
    const at = (v) => [...Array(12).keys()].map((i) => p.sample(v, Date.UTC(p.climatology.year, i, 15), lon, lat));
    const tmax = at(n.variables.tmax), tmin = at(n.variables.tmin), P = at(n.variables.ppt);
    if ([...tmax, ...tmin, ...P].every(Number.isFinite)) {
      const T = meanTemperature(tmax, tmin);
      const code = koppen(T, P);
      rows.push(['Mean temperature, year', `${format(T.reduce((a, b) => a + b, 0) / 12)} °C`]);
      rows.push(['Precipitation, year', `${format(P.reduce((a, b) => a + b, 0))} mm`]);
      if (code && (!c0 || c0.code !== code)) rows.push(['Class of these normals', `<span class="code">${esc(code)}</span>`]);
      const [a, b] = n.period.split('/').map((d) => d.slice(0, 4));
      figure = { svg: climateDiagram({ T, P, Tmin: tmin, lat }), caption: `${n.model}, ${a}–${b}: temperature and precipitation by month, after Walter and Lieth` };
    }
  }
  return { rows, figure };
}

// ---- legend, labels, credits, about ---------------------------------------------

function legendData() {
  const out = {};
  const lg = surface.legend();
  if (lg) {
    const def = surface.def;
    if (def.kind === 'classes') {
      const credit = surface.product.meta.source && surface.product.meta.source.credit;
      lg.note = [credit, surface.compare ? 'faint where the class is the same' : '0.1° cells'].filter(Boolean).join(' · ');
    } else if (def.kind === 'field' && surface.product.climatology) {
      const p = surface.product, b = p.bracket(state.time);
      lg.title = `${def.label.replace(/, normal$/, '').replace(/ in the month$/, '')} in ${MONTH_NAMES[b.month]}, normal`;
      lg.note = [p.meta.source && p.meta.source.credit, `the normal year of ${p.climatology.period.slice(0, 4)}–${p.climatology.period.slice(11, 15)}`].filter(Boolean).join(' · ');
    } else if (def.kind === 'field') {
      const p = surface.product, b = surface.shown || p.bracket(state.time);
      const credit = [p.meta.source && p.meta.source.credit, surface.compare && surface.compare.product.meta.source.credit].filter(Boolean).join(' · ');
      const when = p.times.length ? (b.offset ? formatShortUTC(b.shown) : 'interpolated in time') : '';
      lg.note = [credit, p.times.length === 1 ? `daily mean, ${formatDay(p.times[0])}` : when].filter(Boolean).join(' · ');
      if (b.offset && Math.abs(b.offset) > HOUR) lg.offset = `Data from ${formatOffset(b.offset)}`;
    } else {
      lg.title = `${def.label}, ${surface.year}`;
      lg.note = `${(surface.product.source && surface.product.source.credit) || 'World Bank'} · by country`;
    }
    out.surface = lg;
  }
  if (flow.def && flow.product) {
    // Screen speed (m/s of the Earth's surface) for one unit of the field, against its real speed in m/s.
    const realPerFrame = (0.16 * (flow.def.pace || 1)) / stage.pixelsPerMetre();
    const factor = (realPerFrame * 60) / (flow.def.siPerUnit || 1);
    const b = flow.product.bracket(state.time);
    // A large factor as a power of ten: plates move a few centimetres a year.
    const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
    const e = Math.floor(Math.log10(factor));
    const times = factor < 1e7 ? (Math.round(factor / 1000) * 1000).toLocaleString('en-GB')
      : `${(factor / 10 ** e).toFixed(1)} × 10${String(e).split('').map((c) => SUP[c]).join('')}`;
    out.flow = {
      title: flow.def.label,
      note: `${flow.def.kind === 'waves' ? 'At the group speed g T / 4π, ' : ''}shown ${times} times faster than real`,
      offset: b.offset && Math.abs(b.offset) > HOUR ? `Data from ${formatOffset(b.offset)}` : '',
    };
  }
  if (lines.def && lines.product) {
    const step = lines.def.step * displayUnit(lines.def.units, lines.def.labelUnits).a;
    out.lines = { title: lines.def.label, unit: `every ${+step.toPrecision(3)} ${lines.def.labelUnits}` };
  }
  if (state.plates && plates.plates) out.plates = plates.legend();
  if (state.stations && stations.product) out.stations = stations.legend();
  if (state.mechanisms && mechanisms.product) out.mechanisms = mechanisms.legend();
  if (tracks.sources.length) {
    const keys = tracks.sources.flatMap((s) => (s.def.kind === 'official'
      ? [{ label: 'observed', colour: s.colour }, { label: 'official forecast', colour: s.colour, dashed: true }]
      : [{ label: s.def.model, colour: s.colour }]));
    const count = new Set(tracks.marks.map((m) => m.name || `${m.lon.toFixed(0)},${m.lat.toFixed(0)}`)).size;
    const models = tracks.sources.some((s) => s.def.kind !== 'official');
    const official = tracks.sources.find((s) => s.def.kind === 'official');
    out.tracks = {
      title: 'Tropical cyclones', count: count ? `${count} now` : '', keys,
      note: [models ? 'models: found here in sea-level pressure (a closed 2 hPa contour; Ullrich and Zarzycki 2017)' : '',
        official ? 'official: GDACS, automatic; not a warning' : ''].filter(Boolean).join(' · '),
    };
  }
  if (events.def && events.product && events.product.kind === 'disaster') {
    out.events = { kind: 'disaster', title: 'Disaster alerts', count: `${events.shown.length}`,
      note: `GDACS alert levels, until a week after the last update · automatic; not an official warning` };
  } else if (events.def && events.product) {
    const w = events.product.window;
    out.events = { title: 'Earthquakes', count: `${events.shown.length}`, note: `M ≥ ${events.product.selection.min_magnitude} · ${formatDay(Date.parse(w.start))} to the time shown · ${(events.product.source && events.product.source.credit) || 'USGS'}` };
  }
  return out;
}

function updateCredits() {
  const credits = new Set(['Natural Earth']);
  if (elevation && elevation.meta.source.credit) credits.add(elevation.meta.source.credit);
  for (const p of [surface.def && surface.def.kind === 'field' && surface.product, surface.compare && surface.compare.product, flow.product, lines.product]) {
    if (p && p.meta && p.meta.source.credit) credits.add(p.meta.source.credit);
  }
  if (imagery && imagery.active) credits.add(`${imagery.active.source.credit}${imagery.active.date ? ` (${imagery.active.date})` : ''}`);
  if (detail && detail.source && detail.source.role === 'relief-shading' && detail.target) credits.add(detail.source.source.credit);
  if (state.lights && imagery && imagery.entry('night-lights')) credits.add('NASA Black Marble');
  if (state.clouds && catalogue.clouds.length) credits.add(`clouds: ${catalogue.clouds[0].model}`);
  if (surface.def && surface.def.kind === 'indicator' && surface.product) credits.add(surface.product.source.credit);
  if (surface.def && surface.def.kind === 'classes' && surface.product) {
    for (const r of [surface.product, surface.compare && surface.compare.product]) if (r && r.meta.source) credits.add(r.meta.source.credit);
  }
  if (events.product) credits.add(events.product.source.credit);
  for (const s of tracks.sources) if (s.def.kind === 'official' && s.def.source) credits.add(s.def.source.credit);
  if (state.plates && plates.meta) credits.add(plates.meta.source.credit);
  if (state.stations && stations.product) credits.add(stations.product.source.credit);
  if (state.mechanisms && mechanisms.product) credits.add(mechanisms.product.source.credit);
  $('credits').textContent = [...credits].filter(Boolean).join(' · ');
}

function updateLabels() {
  const items = [];
  if (state.pinned) {
    items.push({ key: 'pin', text: '', html: '', lon: state.pinned[0], lat: state.pinned[1], className: 'pin', font: '12px serif', size: 18, always: true });
  }
  if (lines.def && lines.marks.length) {
    for (const m of lines.marks) {
      items.push({ key: `x${m.lon},${m.lat}`, text: `${m.text} ${m.value}`, html: `${m.text}<small>${m.value}</small>`, lon: m.lon, lat: m.lat, className: `pressure ${m.kind}`, font: '600 15px serif', size: 26, always: true });
    }
  }
  // One label per cyclone: its name, at the official position where there is one,
  // and each model's central pressure (hPa) in that model's colour.
  const cyclones = new Map();
  tracks.marks.forEach((m, i) => {
    const key = m.name || `unnamed${i}`;
    if (!cyclones.has(key)) cyclones.set(key, { name: m.name, at: m, pressures: [] });
    const c = cyclones.get(key);
    if (m.source.kind === 'official') c.at = m;
    if (m.hPa !== undefined) c.pressures.push(`<span style="color:${m.colour}">${Math.round(m.hPa)}</span>`);
  });
  for (const [key, c] of cyclones) {
    const name = c.name ? esc(c.name.replace(/-\d+$/, '')) : '';
    const html = `${name}${c.pressures.length ? `<small>${c.pressures.join(' ')}</small>` : ''}`;
    if (html) items.push({ key: `tc${key}`, text: html, html, lon: c.at.lon, lat: c.at.lat, className: 'cyclone', font: '600 12px serif', size: 18, always: true });
  }
  items.push(...plates.labelItems(stage.zoomLevel()));
  if (state.labels) items.push(...base.labelItems(stage.zoomLevel(), state.lang));
  labels.update(items, (lon, lat) => stage.screen(lon, lat), stage.width, stage.height);
}

function renderAbout() {
  const rows = [];
  const add = (title, s, extra = '') => {
    if (!s) return;
    rows.push(`<dt>${esc(title)}</dt><dd>${esc(s.provider || '')}${s.licence ? `. Licence: ${esc(s.licence)}` : ''}${s.citation ? `<br><i>${esc(s.citation)}</i>` : ''}${s.retrieved_at ? `<br>Retrieved ${esc(s.retrieved_at)}${s.profile ? `, profile ${esc(s.profile)}` : ''}.` : ''}${extra}</dd>`);
  };
  add('Base geography', geography.meta.source);
  const readable = (url) => { try { return decodeURIComponent(url); } catch { return url; } };
  for (const [, m] of library.fieldMeta) add(m.title, m.source, m.source.request ? `<br><code>${esc(readable(m.source.request))}</code>` : '');
  const ind = library.list('indicators');
  if (ind.length) add(`Indicators by country (${ind.length})`, { ...ind[0].source, licence: 'CC BY 4.0' });
  for (const e of library.list('events')) add(e.title, e.source);
  for (const e of library.list('imagery')) add(e.title, { ...e.source, provider: `${e.layer}${e.date ? `, ${e.date}` : ''}` });
  const series = new Set();
  for (const [, m] of library.classMeta) {
    if (m.computed || series.has(m.series)) continue;
    series.add(m.series);
    add(`${m.title.split(',')[0]}, ${[...library.classMeta.values()].filter((x) => x.series === m.series && !x.computed).length} maps`, m.source);
  }
  if (lab.input) rows.push(`<dt>Laboratory models</dt><dd>Two models of the flow at 500 hPa, solved by the spectral transform method at triangular truncation T42 (a Gaussian grid of 128 × 64) with fourth-order Runge–Kutta steps. The barotropic vorticity equation, as in the first numerical forecast (Charney, Fjørtoft and von Neumann 1950), starts from the wind of the GFS, in steps of 30 min; heights follow by linear balance and are added, as changes, to the analysed height. The shallow-water equations start from its wind and height, in steps of 15 min, after a digital filter over six hours removes the gravity waves of an unbalanced start (Lynch and Huang 1992); their finite radius of deformation slows the longest waves, which the barotropic model moves westward too fast. Tested against the Rossby–Haurwitz wave, Williamson et al.'s steady flow across the poles, and gravity waves of known period.</dd>`);
  if (catalogue.normals.length) rows.push(`<dt>Climate</dt><dd>Climate classes follow Köppen (1936) as revised by Peel et al. (2007) and Beck et al. (2018): the arid group first, then the tropical, polar, temperate and cold groups, with summer the warmer half-year. Besides the published maps, the atlas classifies the ${esc(catalogue.normals[0].model)} normals itself, with the mean temperature of a month taken as the mean of its daily maximum and minimum; the two agree in about nine land cells in ten. Climate diagrams follow Walter and Lieth (1960): 10 °C to 20 mm, dry where precipitation falls below twice the temperature (Gaussen).</dd>`);
  if (plates.meta) rows.push(`<dt>Tectonic plates</dt><dd>The outlines of 52 plates and the classes of their boundaries are Bird's model PB2002 (as converted by H. Ahlenius, Nordpil; ODC Attribution licence). Each plate turns rigidly about its pole, at the rate of NNR-MORVEL56 (Argus, Gordon and DeMets 2011), relative to the no-net-rotation frame, the mean motion of the whole surface: a point moves at ω × r, computed here on the WGS 84 ellipsoid. Bird's Africa takes MORVEL's Nubia; MORVEL's Capricorn, Lwandle, Macquarie and Sur, which PB2002 does not outline, move with the plates containing them. Across a boundary the relative motion is resolved into opening (or closing) and sliding, right- or left-lateral. <i>Bird, P. (2003), G-cubed 4(3), 1027; Argus, D. F., Gordon, R. G. and DeMets, C. (2011), G-cubed 12, Q11001.</i></dd>`);
  if (library.list('crust').length) rows.push(`<dt>The crust</dt><dd>CRUST1.0 (Laske, Masters, Ma and Pasyanos 2013): on 1° cells, water, ice, three layers of sediment and three of crystalline crust down to the Moho, each with its P- and S-wave speed and density, and the speed of P waves just below the Moho. A place takes the values of the cell containing it, as the model's authors advise. No licence is stated; cite Laske et al. (2013), EGU2013-2658.</dd>`);
  if (library.list('stations').length) rows.push(`<dt>GPS stations</dt><dd>The velocities of some 20,000 continuous GNSS stations, estimated by MIDAS (Blewitt et al. 2016) at the Nevada Geodetic Laboratory, in IGS20, aligned with ITRF2020; kept where a station has run for two years or more with horizontal uncertainties of 2 mm/yr or less. Relative to its plate, a station's velocity is less the motion of the plate it stands on (NNR-MORVEL56; ITRF2020's own rotation against the no-net-rotation frame, a fraction of a millimetre a year, remains). A station's daily positions are asked for on demand. No licence is stated; data are supplied as is. <i>Blewitt, G., Hammond, W. C. and Kreemer, C. (2018), Eos 99, doi:10.1029/2018EO104623.</i></dd>`);
  if (library.list('events').some((e) => e.kind === 'mechanism')) rows.push(`<dt>Earthquake mechanisms</dt><dd>The centroid moment tensors of the Global CMT Project, about 70,000 earthquakes since 1976 (final solutions, and the quick ones of the present year). Each is drawn as a beach ball: the lower hemisphere of the source in equal-area projection, north up, coloured where the first motion of P waves is a push. Its axes and nodal planes are computed here from the tensor's eigenvectors (they agree with the catalogue's own to within a degree), and it is called thrust, normal or strike-slip by the plunges of its axes (Frohlich 1992). <i>Dziewonski, Chou and Woodhouse (1981), J. Geophys. Res. 86, 2825; Ekström, Nettles and Dziewonski (2012), Phys. Earth Planet. Inter. 200-201, 1.</i></dd>`);
  rows.push(`<dt>The ocean beneath a place</dt><dd>Asked for on demand, through this atlas's server and the laboratory's data library, from the US Navy's global ocean model (ESPC-D-V02, its ocean HYCOM on a 1/12° grid with 40 levels to 5000 m) on HYCOM.org's OPeNDAP server: the column at the model's grid point nearest the place, every 6 hours from two days back to four ahead. Approved for public release, distribution unlimited; HYCOM was developed with funding from the National Ocean Partnership Program and the Office of Naval Research, and its data-assimilative products are funded by the U.S. Navy. <i>Barton, N., et al. (2021), The Navy's Earth System Prediction Capability, Earth and Space Science 8, e2020EA001199.</i> The mixed layer is where temperature first differs from that at 10 m by 0.2 °C (de Boyer Montégut et al. 2004).</dd>`);
  rows.push(`<dt>Method</dt><dd>Fields are interpolated bilinearly in space and linearly in time between the stored steps; a place near missing values takes the mean of the valid neighbours carrying at least half the weight. Wind direction is the direction it comes from. Wave energy travels at the deep-water group speed g T / 4π along the peak direction. The Sun's position follows the Astronomical Almanac's low-precision formulas (0.01°); day, civil, nautical and astronomical twilight are drawn as bands. Countries are joined to indicators by their World Bank and ISO codes; boundaries are Natural Earth's de facto depiction. The Equal Earth projection keeps areas true.</dd>`);
  $('about-body').innerHTML = `<dl>${rows.join('')}</dl>`;
}

// ---- the laboratory's model ---------------------------------------------------------

let labPanelTimer = 0;
async function runLab(kind = 'barotropic') {
  if (lab.running) return;
  lab.running = kind; lab.progress = 0; lab.error = '';
  renderPanel();
  try {
    const result = await labModel.run(state.time, {
      model: kind,
      onProgress: (f) => {
        lab.progress = f;
        clearTimeout(labPanelTimer);
        labPanelTimer = setTimeout(renderPanel, 60);
      },
    });
    const meta = result.product.meta;
    library.fieldMeta.set(meta.name, meta);
    library.fields.set(meta.name, Promise.resolve(result.product));
    loaded.set(meta.name, result.product);
    surface.textures.forget(meta.name);
    catalogue = buildCatalogue(library);
    // Runs from another start are not comparable: they are forgotten.
    if (lab.start !== result.start) lab.runs = {};
    lab.runs[kind] = result.verification;
    lab.start = result.start;
    lab.against = result.against;
    state.live = false;
    state.time = result.start;
    const label = LAB_MODELS[kind].label;
    const height = catalogue.surfaces.find((d) => d.key === 'height|500 hPa|' && d.model === label);
    const contours = catalogue.lines.find((d) => d.key === 'heights|500 hPa' && d.model === label);
    const wind = catalogue.flows.find((d) => d.key === 'wind|500 hPa' && d.model === label);
    lab.running = null;
    await Promise.all([
      height ? selectSurface(height.id, false) : null,
      contours ? selectLines(contours.id, false) : null,
      wind ? selectFlow(wind.id, false) : null,
    ]);
    setTime(state.time);
    updatePlay();
  } catch (e) {
    console.error(e);
    lab.error = `The model could not run: ${e.message}`;
  }
  lab.running = null;
  renderPanel();
  updateTimeline();
  remember();
}

// ---- the URL remembers the view ------------------------------------------------------

let rememberTimer = 0;
function remember() {
  clearTimeout(rememberTimer);
  rememberTimer = setTimeout(() => {
    const v = stage.globeView, q = new URLSearchParams();
    q.set('p', state.projection);
    q.set('c', `${v.lon.toFixed(1)},${v.lat.toFixed(1)},${v.dist.toFixed(2)}`);
    q.set('s', state.surface || 'none'); q.set('f', state.flow || 'none'); q.set('l', state.lines || 'none');
    if (state.compare) q.set('m', state.compare);
    if (!state.lights) q.set('nl', '0');
    if (state.clouds) q.set('cl', '1');
    if (state.events) q.set('e', state.events);
    if (state.tracks.length) q.set('k', state.tracks.join(' '));
    if (state.plates) q.set('pb', '1');
    if (state.stations) q.set('gp', state.stationFrame);
    if (state.station) q.set('st', state.station);
    if (state.mechanisms) q.set('cmt', state.mechanismMode);
    if (state.quake) q.set('eq', state.quake);
    if (!state.live) q.set('t', new Date(state.time).toISOString().slice(0, 16) + 'Z');
    if (state.pinned) q.set('at', `${state.pinned[0].toFixed(3)},${state.pinned[1].toFixed(3)}`);
    if (state.lang !== 'en') q.set('lang', state.lang);
    if (state.style !== 'relief') q.set('style', state.style);
    history.replaceState(null, '', `#${q.toString()}`);
  }, 400);
}

function restore() {
  const q = new URLSearchParams(location.hash.slice(1));
  if (['globe', 'equal-earth', 'equirectangular'].includes(q.get('p'))) state.projection = q.get('p');
  if (q.get('c')) {
    const [lon, lat, dist] = q.get('c').split(',').map(Number);
    if ([lon, lat, dist].every(Number.isFinite)) {
      Object.assign(stage.globeView, { lon, lat, dist });
      stage.updateCamera();
    }
  }
  for (const [k, key] of [['s', 'surface'], ['f', 'flow'], ['l', 'lines'], ['e', 'events'], ['m', 'compare']]) if (q.get(k)) state[key] = q.get(k);
  if (q.get('nl') === '0') state.lights = false;
  if (q.get('k')) state.tracks = q.get('k').split(' ').filter(Boolean);
  if (q.get('pb') === '1') state.plates = true;
  if (['plate', 'igs'].includes(q.get('gp'))) { state.stations = true; state.stationFrame = q.get('gp'); }
  if (/^[A-Z0-9]{4}$/.test(q.get('st') || '')) state.station = q.get('st');
  if (['year', 'large'].includes(q.get('cmt'))) { state.mechanisms = true; state.mechanismMode = q.get('cmt'); }
  if (/^[A-Z0-9]{6,16}$/.test(q.get('eq') || '')) state.quake = q.get('eq');
  if (q.get('cl') === '1') state.clouds = true;
  if (q.get('t')) { const t = Date.parse(q.get('t')); if (Number.isFinite(t)) { state.time = t; state.live = false; } }
  if (q.get('at')) { const a = q.get('at').split(',').map(Number); if (a.every(Number.isFinite)) state.pinned = a; }
  if (q.get('lang')) state.lang = q.get('lang');
  if (q.get('style')) state.style = q.get('style');
}

// ---- the frame loop ------------------------------------------------------------------

let lastFrame = performance.now(), lastLive = 0, lastCoverage = 0, lastRevision = -1, lastDetailTarget = 0, stationsMoved = 0;
function loop(now) {
  requestAnimationFrame(loop);
  const dt = Math.min(100, now - lastFrame);
  lastFrame = now;
  if (state.live && now - lastLive > 20000) { lastLive = now; setTime(Date.now()); }
  if (state.playing) {
    const [a, b] = timeline.span;
    let t = state.time + ((b - a) / 45000) * dt;
    if (t > b) t = a;
    setTime(t);
  }
  stage.coast();
  if (flow.product) {
    const cov = stage.coverage();
    if (Math.abs(cov - lastCoverage) > 0.25 * Math.max(lastCoverage, 0.05)) { lastCoverage = cov; flow.reset(); }
    flow.step();
  }
  // GPS arrows and beach balls are laid out again once the view has settled, and
  // beach balls when the time shown has moved.
  const settled = stationsMoved && now - stationsMoved > 140;
  if (stations.shown && (stations.stale || settled)) { stations.layout(); dirty.legend = true; }
  if (mechanisms.shown && (mechanisms.stale || settled || Math.abs(state.time - mechanisms.time) > (state.playing ? 6 : 0.5) * HOUR)) {
    mechanisms.layout(state.time); dirty.legend = true;
  }
  if (settled) stationsMoved = 0;
  // Names follow the camera, whatever moved it.
  if (stage.revision !== lastRevision) {
    lastRevision = stage.revision;
    stationsMoved = now;
    dirty.labels = true;
    if (detailSource) detail.request(detailSource);
  }
  detail.frame();
  if ((detail.target || 0) !== lastDetailTarget) { lastDetailTarget = detail.target || 0; dirty.legend = true; }
  if (dirty.labels) { dirty.labels = false; flow.viewChanged(); updateLabels(); }
  const res = stage.resolution;
  base.frame(res); flow.frame(res); lines.frame(res); events.frame(res); tracks.frame(res); plates.frame(res); stations.frame(res); mechanisms.frame(res);
  sectionLayer.frame(res, 0);
  stage.render();
  if (dirty.legend) { dirty.legend = false; renderLegend($('legend'), legendData()); updateCredits(); }
  if (dirty.place && now - lastPlace > (state.playing ? 400 : 120)) {
    dirty.place = false; lastPlace = now;
    renderPlace($('place'), placeCard());
    if (hover) showTooltip();
  }
}

start();
