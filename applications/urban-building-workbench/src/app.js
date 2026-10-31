// The Urban Building Workbench's shell: a three-dimensional GIS of a city (the visualization layer's gis: the
// local frame on WGS 84, extruded buildings, streets and spaces from the data library's city model, the Sun of
// the chosen instant), a clock in the city's own time zone, and on top of it the plug-ins: each declares
// layers (modules that draw on the map and answer for what is picked there) and tools (analyses that open in
// a drawer and can act on the map). A plug-in joins by being a directory with a manifest; one of status
// "planned" names what is to come. Nothing here names a plug-in.

import * as THREE from 'three';
import { MapControls } from 'three/addons/controls/MapControls.js';
import * as G from '/visualization/js/gis/index.js';

import { asker, plugins, runSetup, setup } from './api.js';
import * as dom from './dom.js';
import * as chart from './draw/chart.js';
import { finding } from './finding.js';
import { createClock, DAYS } from './clock.js';
import { addCity, legendOfCity, palette } from './city.js';

const { h } = dom;
const state = { plugins: [], setup: null, areas: [], area: null, model: null, gis: null, frame: null, clock: null, layers: new Map(), tool: null, selection: null,
  selectListeners: [], disposers: [], buildingForm: (() => { try { return localStorage.getItem('urban-buildings') || 'low'; } catch { return 'low'; } })() };

function theme(set) {
  let t = set;
  if (!t) { try { t = localStorage.getItem('urban-theme'); } catch { t = null; } }
  t = t === 'night' ? 'night' : 'day';
  document.documentElement.dataset.theme = t;
  document.getElementById('theme').textContent = t === 'night' ? 'day' : 'night';
  try { localStorage.setItem('urban-theme', t); } catch { /* private window */ }
  return t;
}

function lamps() {
  const el = document.getElementById('lamp');
  const engines = state.plugins.filter((p) => p.engine);
  const current = engines.every((p) => p.engineStatus && p.engineStatus.current);
  el.className = `lamp ${current ? 'on' : 'off'}`;
  el.textContent = current ? 'engines current' : 'engines to build';
  const d = document.getElementById('data-lamp');
  const steps = (state.setup && state.setup.steps || []).filter((s) => !/engine/.test(s.id));
  const done = steps.filter((s) => s.done).length;
  d.className = `lamp ${steps.length && done === steps.length ? 'on' : done ? 'part' : 'off'}`;
  d.textContent = steps.length ? `data ${done}/${steps.length}` : 'data';
  d.title = steps.map((s) => `${s.done ? '✓' : '·'} ${s.title}`).join('\n');
  d.style.cursor = 'pointer';
  d.onclick = () => buildPanel();
}

// ---- the drawer ---------------------------------------------------------------------------------------------

function openPanel(title) {
  for (const f of state.disposers.splice(0)) try { f(); } catch { /* gone */ }
  const panel = document.getElementById('panel');
  panel.hidden = false;
  const body = h('div', {});
  panel.replaceChildren(h('button', { class: 'plain close', onclick: closePanel, title: 'close' }, '×'), body);
  document.getElementById('where').textContent = title || '';
  return body;
}

function closePanel() {
  for (const f of state.disposers.splice(0)) try { f(); } catch { /* gone */ }
  document.getElementById('panel').hidden = true;
  state.tool = null;
  for (const b of document.querySelectorAll('.rail button.tool')) b.classList.remove('active');
  document.getElementById('where').textContent = '';
}

async function buildPanel() {
  const body = openPanel('the engine and the records');
  const log = h('pre', { class: 'mono', style: 'white-space: pre-wrap; max-height: 320px; overflow: auto; font-size: 0.8rem' });
  const button = h('button', { class: 'primary' }, 'Build the engine and fetch the records');
  const steps = h('ul', {}, (state.setup && state.setup.steps || []).map((s) => h('li', {}, `${s.done ? '✓' : '·'} ${s.title}`)));
  body.append(h('h2', {}, 'The engine and the records'),
    h('p', { class: 'muted' }, 'The analyses compute in C++ on the urban-science domain and the foundations; the map and the analyses read the data library\'s records: the city model (OpenStreetMap\'s buildings, streets and spaces), holdings and their declared hours, Melbourne\'s pedestrian counts, land use, weather and holidays. Compiling takes a minute or two (g++, natively or in WSL); fetching the records several minutes.'),
    steps, h('p', {}, button), log);
  const poll = async () => {
    const s = await setup();
    state.setup = s;
    log.textContent = s.log || '';
    if (s.running) return setTimeout(poll, 1500);
    state.plugins = await plugins();
    lamps();
    button.disabled = false;
  };
  button.addEventListener('click', async () => { button.disabled = true; const s = await runSetup(); if (s.message) log.textContent = s.message; poll(); });
  poll();
}

function plannedPanel(plugin) {
  const body = openPanel(`${plugin.title} · planned`);
  body.append(h('div', { class: 'view-head' }, h('h1', {}, plugin.title), h('span', { class: 'chip model' }, 'planned')),
    h('p', { class: 'muted' }, plugin.summary || ''),
    ...(plugin.planned || []).map((item) => h('div', { class: 'finding' }, h('div', { class: 'finding-head' }, h('b', {}, item.title)), h('p', {}, item.about))),
    plugin.begins ? h('p', { class: 'note' }, plugin.begins) : null);
}

// ---- the context given to layers and tools ---------------------------------------------------------------------

function select(sel) {
  state.selection = sel;
  showInfo(sel);
  for (const f of state.selectListeners) try { f(sel); } catch (e) { console.error(e); }
}

function context(plugin, extra = {}) {
  return {
    plugin, ask: asker(plugin.name), dom, chart, finding, gis: state.gis, G, THREE, frame: state.frame, area: state.area, clock: state.clock,
    model: state.model, selection: state.selection, select,
    onSelect: (f) => { state.selectListeners.push(f); state.disposers.push(() => { state.selectListeners = state.selectListeners.filter((x) => x !== f); }); },
    openTool: (id, sel) => { if (sel) select(sel); openTool(plugin, plugin.tools.find((t) => t.id === id)); },
    theme: document.documentElement.dataset.theme,
    base: `/plugins/${plugin.name}/`,
    onDispose: (f) => state.disposers.push(f),
    ...extra,
  };
}

async function openTool(plugin, tool) {
  if (!tool) return;
  if (plugin.engine && !(plugin.engineStatus && plugin.engineStatus.current)) return buildPanel();
  state.tool = { plugin, tool };
  for (const b of document.querySelectorAll('.rail button.tool')) b.classList.toggle('active', b.dataset.key === `${plugin.name}/${tool.id}`);
  const body = openPanel(`${plugin.title} · ${tool.title}`);
  body.append(h('p', { class: 'busy' }, 'loading…'));
  try {
    const mod = await import(`/plugins/${plugin.name}/${tool.module}`);
    body.replaceChildren();
    await mod.render(body, context(plugin));
  } catch (error) {
    dom.failure(body, error);
    console.error(error);
  }
}

// ---- the info card --------------------------------------------------------------------------------------------

function showInfo(sel) {
  const info = document.getElementById('info');
  if (!sel) { info.hidden = true; return; }
  info.hidden = false;
  const actions = h('div', { class: 'actions' });
  for (const a of sel.actions || []) actions.append(h('button', { class: 'plain', onclick: a.run }, a.title));
  info.replaceChildren(h('h3', {}, sel.title || ''), ...(sel.lines || []).map((l) => h('div', { class: 'muted' }, l)), sel.extra || null, sel.actions && sel.actions.length ? actions : null);
}

// ---- the rail ---------------------------------------------------------------------------------------------------

function rail() {
  const r = document.getElementById('rail');
  r.replaceChildren();
  r.append(h('h4', {}, 'The city'));
  for (const [id, title] of [['buildings', 'buildings'], ['streets', 'streets and ways'], ['spaces', 'parks, water, squares']]) {
    const box = h('input', { type: 'checkbox', checked: true, onchange: () => state.gis && state.gis.setVisible(id, box.checked) });
    r.append(h('label', { class: 'layer' }, box, title));
  }
  // The buildings as built, lowered to a fifth of their height, or of glass: the streets' water and walkers
  // show between towers.
  const forms = [['full', 'as built'], ['low', 'low'], ['glass', 'glass']];
  const formKeys = h('div', { class: 'keys tight' }, forms.map(([k, name]) => h('button', { class: `key${k === state.buildingForm ? ' on' : ''}`, onclick: () => setBuildingForm(k) }, name)));
  formKeys.dataset.role = 'building-form';
  r.append(h('div', { class: 'faint form-label' }, 'buildings'), formKeys);
  for (const p of state.plugins.filter((x) => x.status !== 'planned')) {
    const layers = (p.layers || []).filter((l) => !l.areas || l.areas.includes(state.area.id));
    if (layers.length) {
      r.append(h('h4', {}, p.title));
      let group = null;
      for (const l of layers) {
        if (l.group && l.group !== group) { group = l.group; r.append(h('h5', {}, group)); }
        const key = `${p.name}/${l.id}`;
        const on = state.layers.get(key)?.on ?? !!l.default;
        const box = h('input', { type: 'checkbox', checked: on });
        box.addEventListener('change', () => toggleLayer(p, l, box.checked));
        r.append(h('label', { class: 'layer', title: l.about || '' }, box, l.title));
      }
    }
    if ((p.tools || []).length) {
      r.append(h('h4', {}, 'Analyses'));
      for (const t of p.tools) r.append(h('button', { class: 'tool', 'data-key': `${p.name}/${t.id}`, title: t.about || '', onclick: () => openTool(p, t) }, t.title));
    }
  }
  const planned = state.plugins.filter((x) => x.status === 'planned');
  if (planned.length) {
    r.append(h('h4', {}, 'Planned'));
    for (const p of planned) r.append(h('button', { class: 'tool', title: p.summary || '', onclick: () => plannedPanel(p) }, p.title));
  }
  r.append(h('div', { class: 'about' }, 'Plug-ins live in plugins/<name>/: a manifest of their layers and analyses, and an engine in C++ on the foundations and the urban-science domain.'));
}

function setBuildingForm(form) {
  state.buildingForm = form;
  try { localStorage.setItem('urban-buildings', form); } catch { /* private window */ }
  for (const b of document.querySelectorAll('[data-role="building-form"] button')) b.classList.toggle('on', b.textContent === { full: 'as built', low: 'low', glass: 'glass' }[form]);
  const l = state.gis && state.gis.layers.get('buildings');
  if (!l) return;
  l.object.scale.y = form === 'low' ? 0.2 : 1;
  l.object.material.transparent = form === 'glass';
  l.object.material.opacity = form === 'glass' ? 0.28 : 1;
  l.object.material.depthWrite = form !== 'glass';
  l.object.castShadow = form !== 'glass';
  l.object.material.needsUpdate = true;
}

async function toggleLayer(plugin, layer, on) {
  const key = `${plugin.name}/${layer.id}`;
  const entry = state.layers.get(key) || { on: false, handle: null, legend: null };
  entry.on = on;
  state.layers.set(key, entry);
  if (!on) {
    if (entry.handle && entry.handle.dispose) entry.handle.dispose();
    for (const id of entry.ids || []) state.gis.remove(id);
    entry.handle = null;
    if (entry.legend) entry.legend.remove();
    return;
  }
  if (plugin.engine && !(plugin.engineStatus && plugin.engineStatus.current)) return buildPanel();
  try {
    const mod = await import(`/plugins/${plugin.name}/${layer.module}`);
    const ids = [];
    const ctx = context(plugin, { layerId: (suffix = '') => { const id = `${key}${suffix}`; ids.push(id); return id; } });
    entry.handle = await mod.layer(ctx);
    entry.ids = ids;
    if (entry.handle && entry.handle.legend) {
      entry.legend = entry.handle.legend;
      document.getElementById('legend').append(entry.legend);
    }
    if (entry.handle && entry.handle.setTime) entry.handle.setTime(state.clock);
  } catch (error) {
    console.error(error);
    const card = h('div', { class: 'card' }, h('b', {}, layer.title), h('span', { class: 'faint' }, String(error.message || error)));
    entry.legend = card;
    document.getElementById('legend').append(card);
  }
}

// ---- the clock ---------------------------------------------------------------------------------------------------

function clockBar() {
  const c = state.clock;
  const bar = document.getElementById('clock');
  const date = h('input', { type: 'date', value: c.date });
  const slider = h('input', { type: 'range', min: 0, max: 23.75, step: 0.25, value: c.hour });
  const label = h('span', { class: 'time' }, '');
  const sun = h('span', { class: 'sun' }, '');
  let playing = null;
  const play = h('button', { class: 'plain', title: 'let the day run' }, '▶');
  const show = () => {
    label.textContent = c.label();
    const s = c.sun;
    sun.textContent = s.altitude > 0 ? `Sun ${Math.round(s.altitude)}° high, ${Math.round(s.azimuth)}°` : `Sun ${Math.round(-s.altitude)}° below`;
  };
  date.addEventListener('change', () => c.set({ date: date.value }));
  slider.addEventListener('input', () => c.set({ hour: +slider.value }));
  play.addEventListener('click', () => {
    if (playing) { clearInterval(playing); playing = null; play.textContent = '▶'; return; }
    play.textContent = '❚❚';
    playing = setInterval(() => { const next = (c.hour + 0.25) % 24; slider.value = next; c.set({ hour: next }); }, 400);
  });
  if (state.stopClock) state.stopClock();
  state.stopClock = () => { if (playing) clearInterval(playing); playing = null; };
  bar.replaceChildren(play, date, slider, label, sun);
  c.on(show);
  show();
}

// ---- an area -------------------------------------------------------------------------------------------------------

async function cityModel() {
  if (state.modelAll !== undefined) return state.modelAll;
  try {
    const r = await fetch('/data/products/city-model.json');
    state.modelAll = r.ok ? await r.json() : null;
  } catch { state.modelAll = null; }
  return state.modelAll;
}

async function openArea(area) {
  // Leave the previous city: its layers, its scene.
  for (const [, e] of state.layers) { if (e.handle && e.handle.dispose) e.handle.dispose(); if (e.legend) e.legend.remove(); e.handle = null; }
  if (state.gis) state.gis.dispose();
  document.getElementById('legend').replaceChildren();
  select(null);
  const all = await cityModel();
  const model = all && (all.areas || []).find((a) => a.id === area.id);
  state.model = model || null;
  const origin = model && model.origin ? model.origin : { lat: area.lat, lon: area.lon };
  state.area = { ...area, lat: origin.lat, lon: origin.lon };
  state.frame = G.localFrameAt(origin.lat, origin.lon);
  const t = document.documentElement.dataset.theme;
  const P = palette(t);
  const mapEl = document.getElementById('map');
  mapEl.replaceChildren();
  state.gis = G.cityScene(THREE, { MapControls }, mapEl, { background: P.background, ground: P.ground });
  let ext = null;
  if (model) {
    ext = addCity(state.gis, model, t);
    setBuildingForm(state.buildingForm);
    document.getElementById('legend').append(legendOfCity(dom, model, t));
  } else {
    mapEl.append(h('div', { class: 'empty-map' }, 'The city model is not built yet: the data lamp, above, runs the recipe that fetches it.'));
  }
  // The first view frames the area's box if it names one (where its records are), else the whole model.
  const box = area.view && area.view.box;
  if (box) {
    const corners = [[box[0], box[1]], [box[0], box[3]], [box[2], box[1]], [box[2], box[3]]].map(([la, lo]) => state.frame.toLocal(la, lo).slice(0, 2));
    ext = state.gis.extentOf(corners);
  }
  state.gis.fit(ext || { x0: -1500, y0: -1500, x1: 1500, y1: 1500 }, area.view || {});
  state.gis.onPick((hit) => {
    if (!hit || !hit.layer) return select(null);
    if (hit.layer === 'buildings') {
      const b = hit.feature;
      return select({ kind: 'building', title: b.name || `a building (${b.kind || 'yes'})`, lines: [`${Math.round(b.height_m)} m${b.height_source === 'default' ? ' (height not mapped: a default)' : b.height_source === 'levels' ? ` (${b.levels} levels)` : ''}`, b.id] });
    }
    if (hit.layer === 'streets') return select({ kind: 'way', title: hit.feature.name || hit.feature.kind, lines: [`${hit.feature.kind}, ${Math.round(hit.feature.width_m)} m wide${hit.feature.width_source === 'default' ? ' (a default width)' : ''}`] });
    for (const [, e] of state.layers) if (e.handle && e.ids && e.ids.includes(hit.layer) && e.handle.describe) return select(e.handle.describe(hit));
  });
  state.gis.onHover((hit, ev) => {
    const box = document.getElementById('hover');
    if (!hit || !hit.layer || !ev) { box.hidden = true; return; }
    const l = state.gis.layers.get(hit.layer);
    const text = l && typeof l.style.pick === 'function' ? l.style.pick(hit.feature, hit.index) : null;
    if (!text) { box.hidden = true; return; }
    box.textContent = text;
    box.hidden = false;
    box.style.left = `${Math.min(ev.clientX + 14, window.innerWidth - 340)}px`;
    box.style.top = `${ev.clientY + 14}px`;
  });
  state.clock = createClock(state.area, state.clock ? { date: state.clock.date, hour: state.clock.hour } : {});
  const applySun = (c) => { const s = c.sun; state.gis.setSun({ azimuth: (s.azimuth * Math.PI) / 180, altitude: (s.altitude * Math.PI) / 180 }); };
  state.clock.on(applySun);
  state.clock.on((c) => { for (const [, e] of state.layers) if (e.handle && e.handle.setTime) e.handle.setTime(c); });
  applySun(state.clock);
  clockBar();
  rail();
  for (const p of state.plugins.filter((x) => x.status !== 'planned'))
    for (const l of (p.layers || []).filter((x) => !x.areas || x.areas.includes(area.id))) {
      const key = `${p.name}/${l.id}`;
      const on = state.layers.get(key)?.on ?? !!l.default;
      if (on) toggleLayer(p, l, true);
    }
  if (state.tool) openTool(state.tool.plugin, state.tool.tool);
}

async function start() {
  theme();
  document.getElementById('theme').addEventListener('click', () => {
    const t = theme(document.documentElement.dataset.theme === 'night' ? 'day' : 'night');
    const P = palette(t);
    if (state.gis) state.gis.setPalette({ background: P.background, ground: P.ground });
    if (state.area) openArea(state.areas.find((a) => a.id === state.area.id));
  });
  try {
    const [ps, s, areas] = await Promise.all([plugins(), setup().catch(() => null), fetch('areas.json').then((r) => r.json())]);
    state.plugins = ps, state.setup = s, state.areas = areas.areas;
  } catch (error) {
    return dom.failure(document.getElementById('map'), error);
  }
  state.plugins.sort((a, b) => (a.status === 'planned') - (b.status === 'planned') || (a.order ?? 99) - (b.order ?? 99));
  lamps();
  const sel = document.getElementById('area');
  sel.replaceChildren(...state.areas.map((a) => h('option', { value: a.id }, a.name)));
  const wanted = (location.hash.slice(1) || '').split('/')[0];
  if (state.areas.some((a) => a.id === wanted)) sel.value = wanted;
  sel.addEventListener('change', () => { location.hash = sel.value; openArea(state.areas.find((a) => a.id === sel.value)); });
  await openArea(state.areas.find((a) => a.id === sel.value));
  const engines = state.plugins.filter((p) => p.engine);
  if (!engines.every((p) => p.engineStatus && p.engineStatus.current)) buildPanel();
}

start();
