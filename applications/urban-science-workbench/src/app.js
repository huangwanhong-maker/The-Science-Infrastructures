// The Urban Science Workbench's shell: one three-dimensional GIS from the whole Earth to a street. The globe of the
// shared visualization layer (terrain and imagery by level of detail, graded into the workbench's palette, lit by the
// real Sun) carries what is known of the world's people and hazards; a city (its buildings, streets and spaces in the
// local frame on WGS 84) stands on it, drawn on a second canvas that follows the same view; a clock in the city's own
// time zone lights both. On top are the plug-ins: each declares layers (modules that draw and answer for what is picked)
// and tools (analyses that open in a sheet and can act on the map), and may add readings to the probe. Each names the
// lens it belongs to, the question it answers. A plug-in joins by being a directory with a manifest; one of status
// "planned" names what is to come. Nothing here names a plug-in.

import * as THREE from 'three';
import { MapControls } from 'three/addons/controls/MapControls.js';
import * as G from '/visualization/js/gis/index.js';
import { subsolarPoint } from '/domains/astronomy/js/index.js';

import { asker, plugins, runSetup, setup } from './api.js';
import * as dom from './dom.js';
import * as chart from './draw/chart.js';
import { finding } from './finding.js';
import { createClock } from './clock.js';
import { addCity, legendOfCity, palette, worldCard } from './city.js';
import { look, } from './shell/look.js';
import { createNavigator, distance, scaleName } from './shell/navigator.js';
import { createRibbon, createYears } from './shell/ribbon.js';
import { createNeatline } from './shell/neatline.js';
import { createLenses } from './shell/lenses.js';
import { createSpaceTimeHost } from './shell/spacetime-host.js';

const { h } = dom;
const DEG = Math.PI / 180;
const remembered = (key, fallback) => { try { return localStorage.getItem(key) ?? fallback; } catch { return fallback; } };
const remember = (key, value) => { try { localStorage.setItem(key, String(value)); } catch { /* private window */ } };

const state = {
  plugins: [], setup: null, areas: [], lenses: [], area: null, model: null, gis: null, earth: null, frame: null, clock: null, layers: new Map(), tool: null, selection: null,
  selectListeners: [], disposers: [], readers: [], places: null, modelIndex: undefined, modelAreas: new Map(),
  buildingForm: remembered('urban-buildings', 'full'),
  // How much of the city is seen (1 solid, toward 0 a ghost), over its buildings, ways, spaces and trees, so that a
  // model drawn above it (a foam in spacetime) can be seen through it; the buildings' form multiplies with it.
  cityOpacity: (() => { const v = parseFloat(remembered('urban-city-solid', '1')); return Number.isFinite(v) ? Math.min(1, Math.max(0.04, v)) : 1; })(),
  separation: 0, groundPlain: remembered('urban-ground', 'image') === 'plain',
  ui: {},
};
const CITY_LAYERS = ['buildings', 'streets', 'spaces', 'trees'];
const themeNow = () => document.documentElement.dataset.theme;

// ---- the theme: the page's, and the scene's look --------------------------------------------------------------------

function applyLook() {
  const L = look(themeNow());
  state.earth.setLook({ ...L, haze: L.sky(0).horizon, tintAmount: state.groundPlain ? 0.96 : L.tintAmount });
}
function theme(set) {
  const t = set || (remembered('urban-theme', 'night') === 'day' ? 'day' : 'night');
  document.documentElement.dataset.theme = t;
  remember('urban-theme', t);
  return t;
}

// ---- lamps, status ---------------------------------------------------------------------------------------------------

function lamps() {
  const el = document.getElementById('status');
  const engines = state.plugins.filter((p) => p.engine);
  const current = engines.every((p) => p.engineStatus && p.engineStatus.current);
  const steps = (state.setup && state.setup.steps || []).filter((s) => !/engine/.test(s.id));
  const done = steps.filter((s) => s.done).length;
  const lamp = (cls, text, title, onclick) => h('button', { class: `lamp button ${cls}`, title, onclick }, h('i', {}), text);
  const key = h('button', { class: 'lamp button theme-key', title: 'a plum night, or lilac paper', onclick: () => {
    theme(themeNow() === 'night' ? 'day' : 'night');
    applyLook();
    state.ribbon?.redraw();
    if (state.area) openArea(state.area, { view: 'keep', keepLayers: true });
    lamps();
  } }, h('i', { style: themeNow() === 'night' ? 'background:var(--paper-2)' : 'background:var(--accent)' }), themeNow() === 'night' ? 'night' : 'day');
  el.replaceChildren(
    lamp(steps.length && done === steps.length ? 'on' : done ? 'part' : 'off', steps.length ? `records ${done}/${steps.length}` : 'records', steps.map((s) => `${s.done ? '✓' : '·'} ${s.title}`).join('\n'), buildPanel),
    lamp(current ? 'on' : 'off', current ? 'engines current' : 'engines to build', 'the analyses compute in C++ on the domains and the foundations', buildPanel),
    key);
}

// ---- the sheet -----------------------------------------------------------------------------------------------------

function openPanel(title) {
  for (const f of state.disposers.splice(0)) try { f(); } catch { /* gone */ }
  const panel = document.getElementById('panel');
  panel.hidden = false;
  document.body.classList.add('has-folio');
  const body = h('div', {});
  panel.replaceChildren(h('button', { class: 'close', onclick: closePanel, title: 'close' }, '×'), title ? h('div', { class: 'folio-kicker' }, title) : null, body);
  panel.scrollTop = 0;
  return body;
}

function closePanel() {
  for (const f of state.disposers.splice(0)) try { f(); } catch { /* gone */ }
  document.getElementById('panel').hidden = true;
  document.body.classList.remove('has-folio');
  state.tool = null;
  state.ui.lenses?.setTool(null);
}

async function buildPanel() {
  const body = openPanel('The engines and the records');
  const log = h('pre', { class: 'mono', style: 'white-space: pre-wrap; max-height: 320px; overflow: auto; font-size: 0.8rem' });
  const button = h('button', { class: 'primary' }, 'Build the engines and fetch the records');
  const steps = h('ul', {}, (state.setup && state.setup.steps || []).map((s) => h('li', {}, `${s.done ? '✓' : '·'} ${s.title}`)));
  body.append(h('h1', {}, 'The engines and the records'),
    h('p', { class: 'muted' }, 'The analyses compute in C++ on the domains and the foundations; the map and the analyses read the data library\'s records. Compiling takes a few minutes (g++, natively or in WSL); fetching the records several more. Terrain, imagery and the data asked for on demand come as they are needed.'),
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

function plannedPanel() {
  const body = openPanel('Planned');
  for (const plugin of state.plugins.filter((p) => p.status === 'planned')) {
    body.append(h('h1', {}, plugin.title), h('p', { class: 'muted' }, plugin.summary || ''),
      ...(plugin.planned || []).map((item) => h('div', { class: 'finding' }, h('div', { class: 'finding-head' }, h('b', {}, item.title)), h('p', {}, item.about))),
      plugin.begins ? h('p', { class: 'note' }, plugin.begins) : null);
  }
}

// ---- the context given to layers and tools ---------------------------------------------------------------------

function select(sel) {
  state.selection = sel;
  showInfo(sel);
  for (const f of state.selectListeners) try { f(sel); } catch (e) { console.error(e); }
}

function context(plugin, extra = {}) {
  return {
    plugin, ask: asker(plugin.name), dom, chart, finding, gis: state.gis, G, THREE, frame: state.frame, area: state.area, clock: state.clock, world: state.world, st: state.spacetime,
    model: state.model, selection: state.selection, select, earth: state.earth,
    onSelect: (f) => { state.selectListeners.push(f); state.disposers.push(() => { state.selectListeners = state.selectListeners.filter((x) => x !== f); }); },
    openTool: (id, sel) => { if (sel) select(sel); openTool(plugin, plugin.tools.find((t) => t.id === id)); },
    cityOpacity: () => state.cityOpacity, setCityOpacity: (v, o) => setCityOpacity(v, o),
    theme: themeNow(),
    base: `/plugins/${plugin.name}/`,
    onDispose: (f) => state.disposers.push(f),
    // A scale of years above the clock, for a layer that holds a record over years.
    years: (spec) => state.years.open(spec),
    ...extra,
  };
}

async function openTool(plugin, tool) {
  if (!tool) return;
  if (plugin.engine && !(plugin.engineStatus && plugin.engineStatus.current)) return buildPanel();
  state.tool = { plugin, tool };
  state.ui.lenses?.setTool(`${plugin.name}/${tool.id}`);
  const lensTitle = state.lenses.find((l) => l.id === (tool.lens || plugin.lens))?.title;
  const body = openPanel(lensTitle ? `${lensTitle} · ${tool.title}` : tool.title);
  body.append(h('p', { class: 'busy' }, 'loading…'));
  try {
    const mod = await import(`/plugins/${plugin.name}/${tool.module}`);
    body.replaceChildren();
    await mod.render(body, context(plugin, { spec: tool }));
  } catch (error) {
    dom.failure(body, error);
    console.error(error);
  }
}

// ---- the probe: what is known at the point picked ---------------------------------------------------------------------

const dms = (v, pos, neg) => {
  const a = Math.abs(v), d = Math.floor(a), mf = (a - d) * 60, m = Math.floor(mf), s = Math.round((mf - m) * 60);
  const [mm, ss] = s === 60 ? [m + 1, 0] : [m, s];
  return `${d}°${String(mm).padStart(2, '0')}′${String(ss).padStart(2, '0')}″${v >= 0 ? pos : neg}`;
};
const coords = (lat, lon) => `${dms(lat, 'N', 'S')} ${dms(lon, 'E', 'W')}`;

function showInfo(sel) {
  const info = document.getElementById('info');
  state.probe = sel && sel.at ? { ...sel.at } : null;
  if (!sel) { info.hidden = true; document.getElementById('leader').replaceChildren(); state.ui.probeMark && (state.ui.probeMark.style.display = 'none'); return; }
  info.hidden = false;
  const actions = h('div', { class: 'actions' });
  for (const a of sel.actions || []) actions.append(h('button', { onclick: a.run }, a.title));
  const readings = h('div', {});
  info.replaceChildren(...[
    h('button', { class: 'close', title: 'close', onclick: () => select(null) }, '×'),
    h('h3', {}, sel.title || 'Here'),
    sel.at ? h('div', { class: 'where' }, coords(sel.at.lat, sel.at.lon)) : null,
    ...(sel.lines || []).map((l) => h('div', { class: 'muted' }, l)),
    sel.extra, sel.actions && sel.actions.length ? actions : null, readings,
  ].filter(Boolean));
  if (sel.at) readAt(sel.at, readings);
}

// Each plug-in that offers a reader adds what it knows of the point, as it arrives.
let probeToken = 0;
function readAt(at, holder) {
  const token = ++probeToken;
  for (const { plugin, read } of state.readers) {
    const slot = h('div', {});
    holder.append(slot);
    Promise.resolve().then(() => read(context(plugin), at)).then((answer) => {
      if (token !== probeToken || !answer) return;
      const nodes = [];
      for (const part of [].concat(answer)) {
        if (part.title) nodes.push(h('div', { class: 'reading-head' }, part.title));
        for (const r of part.rows || []) nodes.push(h('div', { class: 'reading' }, h('span', { class: 'label' }, r.label), h('span', { class: 'value' }, r.value), r.source ? h('span', { class: 'source' }, r.source) : null));
        if (part.node) nodes.push(part.node);
      }
      slot.replaceChildren(...nodes);
    }).catch((e) => console.warn(e));
  }
}

// ---- the city's look ---------------------------------------------------------------------------------------------------

// A layer or a tool for the city in view. A place opened from the world's archives has only what holds anywhere.
const inArea = (x) => (state.area?.generic ? !!(x.global || x.anywhere) : !x.areas || (state.area && x.areas.includes(state.area.id)));

function setBuildingForm(form) {
  state.buildingForm = form;
  remember('urban-buildings', form);
  for (const b of document.querySelectorAll('[data-role="building-form"] button')) b.classList.toggle('on', b.textContent === { full: 'as built', low: 'low', glass: 'glass' }[form]);
  if (!state.gis) return;
  const l = state.gis.layers.get('buildings');
  if (l) {
    l.object.scale.y = form === 'low' ? 0.2 : 1;
    l.object.material.userData.formOpacity = form === 'glass' ? 0.28 : 1;
    l.object.castShadow = form !== 'glass';
  }
  for (const [id, b] of state.gis.layers) {
    if (b.kind !== 'built-solids') continue;
    state.gis.setBuiltForm(id, form === 'low' ? 0.2 : 1);
    b.object.castShadow = form !== 'glass';
  }
  applyOpacity();
}

// Remembered only when a person chose it: a layer that dims the city for itself gives it back.
function setCityOpacity(v, { keep = true } = {}) {
  state.cityOpacity = Math.min(1, Math.max(0.04, Number.isFinite(v) ? v : 1));
  if (keep) remember('urban-city-solid', state.cityOpacity);
  for (const el of document.querySelectorAll('[data-role="city-opacity"]')) if (parseFloat(el.value) !== state.cityOpacity) el.value = String(state.cityOpacity);
  applyOpacity();
}

// Each of the city's layers at its own opacity (the buildings' form, or the style it was drawn with) times the
// city's; a layer drawn transparent stops writing depth, so that what lies within or above it shows.
function applyOpacity() {
  if (!state.gis) return;
  for (const [id, l] of state.gis.layers) {
    const painted = l.kind === 'built-paint', solid = l.kind === 'built-solids';
    if (!CITY_LAYERS.includes(id) && !painted && !solid) continue;
    if (!l.object || !l.object.material) continue;
    const m = l.object.material;
    if (m.userData.styleOpacity === undefined) m.userData.styleOpacity = m.opacity;
    const glass = state.buildingForm === 'glass' ? 0.28 : 1;
    const own = solid ? glass : id === 'buildings' && m.userData.formOpacity !== undefined ? m.userData.formOpacity : m.userData.styleOpacity;
    const a = own * state.cityOpacity;
    if (m.opacity === a) continue;
    m.transparent = a < 1;
    m.opacity = a;
    if (!painted) m.depthWrite = a >= 0.999;
    m.needsUpdate = true;
  }
}

// The layers lifted from the ground as sheets, each above the last: the exploded view of a GIS's stack.
function applySeparation() {
  if (!state.gis) return;
  const step = state.separation * state.earth.view.range * 0.16;
  let k = 0;
  for (const [, e] of state.layers) {
    if (!e.on || !e.ids || !e.ids.length) continue;
    k++;
    for (const id of e.ids) state.gis.setLift(id, step * k);
  }
}

function viewKeys() {
  const el = document.getElementById('view-keys');
  const forms = [['full', 'as built'], ['low', 'low'], ['glass', 'glass']];
  const formKeys = h('div', { class: 'keys-tight', title: 'the buildings: as built, lowered to a fifth of their height, or of glass, so that the streets show between towers' },
    forms.map(([k, name]) => h('button', { class: k === state.buildingForm ? 'on' : '', onclick: () => setBuildingForm(k) }, name)));
  formKeys.dataset.role = 'building-form';
  const seen = h('input', { type: 'range', min: '0.04', max: '1', step: '0.02', value: String(state.cityOpacity), 'aria-label': 'how solid the city is drawn', oninput: () => setCityOpacity(parseFloat(seen.value)) });
  seen.dataset.role = 'city-opacity';
  const apart = h('input', { type: 'range', min: '0', max: '1', step: '0.01', value: '0', 'aria-label': 'how far the layers stand apart', oninput: () => { state.separation = parseFloat(apart.value); applySeparation(); } });
  const ground = h('div', { class: 'keys-tight', title: 'the ground beneath the city' },
    [['image', 'imagery'], ['plain', 'plain']].map(([k, name]) => h('button', { class: (k === 'plain') === state.groundPlain ? 'on' : '', onclick: (e) => {
      state.groundPlain = k === 'plain';
      remember('urban-ground', k);
      for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget);
      applyLook();
    } }, name)));
  el.replaceChildren(formKeys, ground,
    h('label', { class: 'view-slider', title: 'Lower it to see what a layer draws above or within the city' }, h('span', {}, 'city solid'), seen),
    h('label', { class: 'view-slider', title: 'Lift each layer from the ground as a sheet, to see how they stand over one another' }, h('span', {}, 'layers apart'), apart));
}

// ---- layers --------------------------------------------------------------------------------------------------------

async function toggleLayer(plugin, layer, on, { quiet = false } = {}) {
  const key = `${plugin.name}/${layer.id}`;
  // Rows of one `exclusive` group stand in the same place: one at a time.
  if (on && layer.exclusive) for (const p of state.plugins) for (const other of p.layers || []) {
    if (other.exclusive === layer.exclusive && other !== layer && state.layers.get(`${p.name}/${other.id}`)?.on) toggleLayer(p, other, false);
  }
  const entry = state.layers.get(key) || { on: false, handle: null, legend: null, global: !!layer.global };
  entry.spec = layer;
  entry.on = on;
  // Each load has a token: a layer switched off (or on again) while it loads drops what that load drew.
  const token = (entry.token = {});
  state.layers.set(key, entry);
  if (!on) {
    if (entry.handle && entry.handle.dispose) entry.handle.dispose();
    for (const id of entry.ids || []) state.gis.remove(id);
    entry.handle = null;
    if (entry.legend) entry.legend.remove();
    applySeparation();
    state.ui.lenses?.refresh();
    return;
  }
  if (plugin.engine && !(plugin.engineStatus && plugin.engineStatus.current) && !layer.withoutEngine) { entry.on = false; return quiet ? undefined : buildPanel(); }
  entry.loading = true;
  state.ui.lenses?.refresh();
  try {
    const mod = await import(`/plugins/${plugin.name}/${layer.module}`);
    const ids = [];
    entry.ids = ids;
    const ctx = context(plugin, { spec: layer, layerId: (suffix = '') => { const id = `${key}${suffix}`; ids.push(id); return id; } });
    const handle = await mod.layer(ctx);
    if (entry.token !== token) {
      if (handle && handle.dispose) handle.dispose();
      if (!entry.on) for (const id of ids) state.gis.remove(id);  // on again: the newer load draws under the same ids
      return;
    }
    entry.handle = handle;
    if (entry.handle && entry.handle.legend) {
      entry.legend = tidyLegend(entry.handle.legend);
      entry.legend.dataset.global = layer.global ? '1' : '';
      if (entry.handle.onLegend) entry.handle.onLegend(entry.legend);
      document.getElementById('legend').append(entry.legend);
      // Opened, unless it was turned on by itself somewhere it means nothing.
      const range = state.earth.view.range;
      if (!quiet || !layer.scale || (range >= layer.scale[0] * 0.5 && range <= layer.scale[1] * 1.7)) openLegend(entry.legend);
      else entry.legend.classList.add('folded');
    }
    if (entry.handle && entry.handle.setTime) entry.handle.setTime(state.clock);
    applySeparation();
  } catch (error) {
    console.error(error);
    const card = h('div', { class: 'card' }, h('b', {}, layer.title), h('span', { class: 'faint' }, String(error.message || error)));
    entry.legend = card;
    document.getElementById('legend').append(card);
  } finally {
    entry.loading = false;
    state.ui.lenses?.refresh();
  }
}

// Legends fold to their titles, one open at a time, so that they do not stand over the map: the latest turned on, or the one asked for.
function openLegend(card) {
  for (const c of document.getElementById('legend').children) c.classList.toggle('folded', c !== card);
}
// A legend's prose folds away: its title, swatches and controls stay, as the eye wants them.
function tidyLegend(card) {
  const head = card.querySelector(':scope > b');
  if (head && !head.dataset.folds) {
    head.dataset.folds = '1';
    head.addEventListener('click', () => openLegend(card));
  }
  const long = [...card.children].filter((c) => !c.matches('b, details, .ld-seg, .ramp, .ramp-ends, .world-count, .world-share') && !c.matches('button, input, select, svg, canvas') && !c.querySelector('button, input, select, svg, canvas, i.dot, i[style]') && c.textContent.length > 64);
  if (long.length) card.append(h('details', { class: 'ld-fold' }, h('summary', {}, 'about this layer'), h('div', { class: 'ld-fold-body' }, long)));
  return card;
}

// The default scales: a city's layers from the doorstep to the metropolis.
const CITY_SCALE = [40, 60000];

// Where a layer or a tool can be had: the names of the places, or of the workbench's cities, when it is not of anywhere.
function placesOf(x) {
  if (x.global || x.anywhere) return 'anywhere';
  const ids = x.areas || state.areas.filter((a) => a.id !== 'freetown').map((a) => a.id);
  return ids.map((id) => state.areas.find((a) => a.id === id)?.short || id).join(' and ');
}
function lensItems() {
  const items = [], plan = [], works = [];
  for (const p of state.plugins) {
    if (p.status === 'planned') { plan.push(p); continue; }
    for (const l of p.layers || []) {
      const key = `${p.name}/${l.id}`;
      items.push({
        kind: 'layer', key, plugin: p, spec: l, lens: l.lens || p.lens || 'fabric', group: l.group, title: l.title, about: l.about, scale: l.scale || CITY_SCALE,
        here: inArea(l), places: placesOf(l), areas: l.areas,
        isOn: () => state.layers.get(key)?.on ?? false, isLoading: () => state.layers.get(key)?.loading ?? false,
        isAway: () => !l.global && state.ui.cityFar === true,
      });
    }
    for (const t of p.tools || []) {
      items.push({ kind: 'tool', key: `${p.name}/${t.id}`, plugin: p, spec: t, lens: t.lens || p.lens || 'fabric', title: t.title, about: t.about, here: inArea(t), places: placesOf(t), areas: t.areas });
    }
    works.push({ name: p.name, title: p.title, tagline: p.tagline || '', lens: p.lens || (p.layers?.[0]?.lens) || 'fabric' });
  }
  return { items, plan, works };
}

// The works of the place open, under its name: what was written for it, and what holds anywhere is not named again.
function placeWorks() {
  const el = document.getElementById('place-works');
  if (!el) return;
  const near = state.area && !state.ui.cityFar;
  if (!near) { if (el.childElementCount) { el.replaceChildren(); document.body.classList.remove('has-works'); } return; }
  const key = `${state.area.id}`;
  if (el.dataset.area === key && el.childElementCount) return;
  el.dataset.area = key;
  const here = state.plugins.filter((p) => p.status !== 'planned' && [...(p.layers || []), ...(p.tools || [])].some((x) => x.areas ? x.areas.includes(state.area.id) : !state.area.generic && !x.global && !x.anywhere));
  const list = here.length ? here : state.plugins.filter((p) => ['territory', 'access', 'studio'].includes(p.name));
  el.replaceChildren(...list.map((p) => h('button', { title: p.tagline || p.summary || '', onclick: () => state.ui.lenses.showWork(p.name) }, p.title)));
  document.body.classList.toggle('has-works', list.length > 0);
}

// A layer or a tool of another place: go there, and then turn it on.
async function goThenDo(d, act) {
  const id = d.areas?.[0] || state.areas[0].id;
  const area = state.areas.find((a) => a.id === id);
  if (!area) return;
  await openArea(area, { view: 'fly' });
  act();
  state.ui.lenses.refresh();
}

function onToolRow(t) {
  if (!inArea(t.spec)) return goThenDo(t, () => openTool(t.plugin, t.spec));
  return openTool(t.plugin, t.spec);
}

function onLayerRow(d, on) {
  if (on && !inArea(d.spec)) { goThenDo(d, () => toggleLayer(d.plugin, d.spec, true)); return; }
  toggleLayer(d.plugin, d.spec, on);
  if (!on) return;
  // A layer asked for where it cannot be seen: go to where it can.
  const v = state.earth.view;
  if (!d.spec.global && state.ui.cityFar) return goToArea(state.area);
  if (v.range < d.scale[0] || v.range > d.scale[1]) {
    const target = Math.sqrt(d.scale[0] * d.scale[1]);
    state.earth.flyTo({ range: Math.min(d.scale[1] * 0.7, Math.max(d.scale[0] * 1.4, target)), tilt: target > 2e5 ? 0 : v.tilt }, 1.8);
  }
}

// ---- places: the pins, the gazetteer --------------------------------------------------------------------------------

function drawPins() {
  const pins = document.getElementById('pins');
  state.ui.pins = state.areas.map((a) => {
    const el = h('button', { class: 'pin', onclick: () => goToArea(a), title: a.about || '' },
      h('span', { class: 'flag' }, h('b', {}, a.short || a.name), a.tag ? h('span', {}, a.tag) : null), h('span', { class: 'stem' }), h('span', { class: 'foot' }));
    pins.append(el);
    return { a, el };
  });
  const mark = h('div', { class: 'probe-mark', style: 'display:none' });
  pins.append(mark);
  state.ui.probeMark = mark;
}

function updatePins(v) {
  const near = v.range < 1.1e6;
  for (const { a, el } of state.ui.pins) {
    const p = state.earth.project(a.lon, a.lat, 0);
    const show = p && !(near && v.range < 2.5e5 && state.area && a.id === state.area.id);
    el.style.display = show ? 'flex' : 'none';
    if (show) el.style.transform = `translate(${p.x.toFixed(1)}px, ${(p.y - el.offsetHeight + 4).toFixed(1)}px)`;
  }
  const probe = state.probe;
  const leader = document.getElementById('leader');
  if (probe) {
    const p = state.earth.project(probe.lon, probe.lat, state.earth.city ? state.earth.city.height : 0);
    state.ui.probeMark.style.display = p ? 'block' : 'none';
    if (p) {
      state.ui.probeMark.style.transform = `translate(${p.x.toFixed(1)}px, ${p.y.toFixed(1)}px) rotate(45deg)`;
      const r = document.getElementById('info').getBoundingClientRect();
      const x2 = r.right, y2 = Math.min(r.bottom, r.top + 26);
      leader.replaceChildren(dom.svg('line', { x1: p.x, y1: p.y, x2, y2 }));
    } else leader.replaceChildren();
  }
}

function goToArea(a, { seconds = 3.2 } = {}) {
  if (!state.area || state.area.id !== a.id) return openArea(a, { view: 'fly' });
  const ext = state.ui.cityExtent;
  state.gis.fit(ext, { ...a.view, seconds });
}

function placeTitle(v) {
  const lat = v.lat / DEG, lon = v.lon / DEG;
  const near = state.areas.find((a) => Math.hypot((a.lat - lat) * 111e3, (a.lon - lon) * 111e3 * Math.cos(a.lat * DEG)) < Math.max(40e3, v.range));
  const here = state.area?.generic && v.range < 2.5e5 ? state.area.title : null;
  const title = here || (near && v.range < 1.2e6 ? (v.range > 2.5e4 ? near.short || near.name : near.title || near.name) : nameOfView(lat, lon, v.range));
  const el = document.getElementById('place-title');
  if (el.textContent !== title) el.textContent = title;
  const line = document.getElementById('place-line');
  line.replaceChildren(`${coords(lat, lon)}  `, h('b', {}, distance(v.range)), `  ${scaleName(v.range).toLowerCase()}`);
}

async function gazetteer() {
  const box = document.getElementById('gazetteer');
  if (!box.hidden) { box.hidden = true; return; }
  box.hidden = false;
  if (!state.places) {
    try { state.places = (await (await fetch('/data/products/places/earth.json')).json()).places; } catch { state.places = []; }
  }
  const input = h('input', { type: 'search', placeholder: 'A place on Earth, or a way of seeing it', 'aria-label': 'a place on Earth, or a way of seeing it' });
  const list = h('ul', {});
  const fly = (lon, lat, range) => { box.hidden = true; state.earth.flyTo({ lon: lon * DEG, lat: lat * DEG, height: state.earth.tree.heightAt(lon * DEG, lat * DEG) || 0, range, tilt: range < 3e5 ? 56 * DEG : 0 }, 4.2); };
  const draw = () => {
    const q = input.value.trim().toLowerCase();
    const rows = [];
    for (const a of state.areas) if (!q || (a.title || a.name).toLowerCase().includes(q) || a.name.toLowerCase().includes(q)) rows.push(h('li', { onclick: () => { box.hidden = true; goToArea(a); } }, h('i', {}, a.title || a.name), h('span', {}, a.tag || 'in three dimensions')));
    if (q) for (const d of lensItems().items) {
      if (rows.length > 10) break;
      if (!`${d.title} ${d.plugin.title} ${d.about || ''}`.toLowerCase().includes(q)) continue;
      rows.push(h('li', { class: 'way', onclick: () => { box.hidden = true; if (d.kind === 'tool') onToolRow(d); else onLayerRow(d, true); state.ui.lenses.refresh(); } }, d.title, h('span', {}, `${d.plugin.title}${d.here ? '' : `, in ${d.places}`}`)));
    }
    if (q) for (const p of state.places) {
      if (rows.length > 24) break;
      if (p[0].toLowerCase().includes(q)) rows.push(h('li', { onclick: () => fly(p[1], p[2], p[3] === 'Country' ? 3.5e6 : 3200) }, p[0], h('span', {}, p[3] === 'Country' ? (p[4] || '') : (p[4] || p[3]))));
    }
    if (!q) rows.push(h('li', { onclick: () => { box.hidden = true; state.earth.flyTo({ range: 2e7, tilt: 0, heading: 0 }, 3); } }, h('i', {}, 'The whole Earth'), h('span', {}, '')));
    list.replaceChildren(...rows);
  };
  input.addEventListener('input', draw);
  box.replaceChildren(input, list);
  draw();
  input.focus();
}

// ---- any place ---------------------------------------------------------------------------------------------------------

// The gazetteer, read once: names by the data library's policy (ADR 0069).
async function placesList() {
  if (!state.places) {
    try { state.places = (await (await fetch('/data/products/places/earth.json')).json()).places; } catch { state.places = []; }
  }
  return state.places;
}
const metres = (la0, lo0, la1, lo1) => Math.hypot((la0 - la1) * 111e3, (((lo0 - lo1 + 540) % 360) - 180) * 111e3 * Math.cos(la0 * DEG));
// A nearby place's name, if the gazetteer has one within `within` metres: the largest kinds first have no claim, the nearest does.
function placeNear(lat, lon, within) {
  let best = null, d = within;
  for (const p of state.places || []) {
    if (p[3] === 'Country') continue;
    const q = metres(lat, lon, p[2], p[1]);
    if (q < d) { d = q; best = p; }
  }
  return best;
}
function nameOfView(lat, lon, range) {
  if (range > 1.2e6) return 'The Earth';
  const p = placeNear(lat, lon, Math.max(60e3, Math.min(range * 0.5, 2e5)));
  return p ? p[0] : 'The Earth';
}
// A place with no model of its own. Its clock keeps the hour of its longitude, since the workbench holds no map of time zones, and says so.
function placeAt(lat, lon) {
  const n = Math.round(lon / 15), name = placeNear(lat, lon, 50e3);
  const zone = n === 0 ? 'Etc/UTC' : `Etc/GMT${n > 0 ? '-' : '+'}${Math.abs(n)}`;
  const title = name ? name[0] : coords(lat, lon);
  return { id: `@${lat.toFixed(4)},${lon.toFixed(4)}`, generic: true, name: title, title, short: title, lat, lon, timezone: zone, zoneNote: `UTC${n >= 0 ? '+' : '−'}${Math.abs(n)}, by longitude`, ground_m: state.earth.tree.heightAt(lon * DEG, lat * DEG) || 0 };
}

// The world's tiles about the view, in whichever city is open: those beyond its model, or all of them where it has none.
function worldLevel() {
  const c = state.earth.city;
  return c && c.level ? { cx: c.cx, cy: c.cy, halfEast: c.halfEast, halfNorth: c.halfNorth, feather: c.feather, height: c.height, amount: 1 } : null;
}
function createWorldTiles() {
  state.world = G.createWorld({
    base: location.origin,
    scene: () => state.gis,
    frame: () => ({ lat: state.area.lat, lon: state.area.lon, height: state.earth.city.height }),
    level: worldLevel,
    keep: () => state.ui.modelBox || null,
    colours: () => palette(themeNow()),
    form: () => (state.buildingForm === 'low' ? 0.2 : 1),
    onStatus: (s) => { state.ui.worldCard?.update(s); },
  });
  state.world.onAdded(() => applyOpacity());
}

// ---- an area -------------------------------------------------------------------------------------------------------

// Whose this is: said where it is drawn, as their terms ask. The imagery and terrain at any height, a city's ways and buildings near the ground.
function credits(v) {
  const el = document.getElementById('credits');
  const near = v.range < 6e4;
  const text = `Imagery: Sentinel-2 cloudless by EOX; Black Marble, NASA. Terrain: Mapzen, ETOPO 2022.${near ? ' Buildings and ways: © OpenStreetMap contributors, Overture Maps Foundation; Google Open Buildings, Microsoft. Heights estimated: GHS-BUILT-H, JRC.' : ''}`;
  if (el.textContent !== text) el.textContent = text;
}

// The city model: an index of its areas (city-model.json), each area in a file of its own, read when it is opened.
async function cityModel(id) {
  if (state.modelIndex === undefined) {
    try {
      const r = await fetch('/data/products/city-model.json');
      state.modelIndex = r.ok ? await r.json() : null;
    } catch { state.modelIndex = null; }
  }
  const entry = state.modelIndex && (state.modelIndex.areas || []).find((a) => a.id === id);
  if (!entry) return null;
  if (!state.modelAreas.has(id)) {
    try {
      const r = await fetch(`/data/products/${entry.file}`);
      state.modelAreas.set(id, r.ok ? await r.json() : null);
    } catch { state.modelAreas.set(id, null); }
  }
  return state.modelAreas.get(id);
}

async function openArea(area, { view = 'fly', keepLayers = false } = {}) {
  const token = (state.ui.openToken = {});
  state.ui.opening = true;
  const model = area.generic ? null : await cityModel(area.id);
  if (token !== state.ui.openToken) return;
  state.world.reset();
  // Leave the previous city: its layers, its scene. Global layers stay where they are.
  for (const [, e] of state.layers) {
    if (e.global && !keepLayers) continue;
    e.token = {};
    if (e.handle && e.handle.dispose) e.handle.dispose();
    if (e.legend) e.legend.remove();
    e.handle = null;
    e.loading = false;
  }
  if (state.gis) state.gis.dispose();
  document.getElementById('legend').querySelector('.city-card')?.remove();
  select(null);
  state.model = model || null;
  const origin = model && model.origin ? model.origin : { lat: area.lat, lon: area.lon };
  state.area = { ...area, lat: origin.lat, lon: origin.lon };
  state.frame = G.localFrameAt(origin.lat, origin.lon);
  const t = themeNow(), P = palette(t), L = look(t);
  const mapEl = document.getElementById('map');
  mapEl.querySelector('.empty-map')?.remove();
  // The city's extent is known before its scene: the Earth's ground is made level under it.
  const pts = [];
  for (const b of (model && model.buildings) || []) { const p = b.footprint && b.footprint[0]; if (p) pts.push(p); }
  for (const r of (model && model.roads) || []) { const p = r.line && r.line[0]; if (p) pts.push(p); }
  const e0 = pts.length ? G.extent(pts) : null;
  state.earth.bindCity({ lat: origin.lat, lon: origin.lon, height: area.ground_m ?? 0, box: e0 ? [e0.x0, e0.y0, e0.x1, e0.y1] : null, level: !!e0 });
  state.ui.modelBox = e0 || null;
  state.gis = G.cityScene(THREE, { MapControls }, mapEl, { driver: state.earth.driver(), eventElement: mapEl, background: P.background, ground: P.ground, shadowColour: L.shadow, shadowOpacity: L.shadowOpacity, lightScale: L.lightScale });
  state.gis.renderer.domElement.classList.add('city-canvas');
  // The flat ground beneath a model takes its shadows; beyond it, and everywhere else, each of the world's tiles brings its own.
  state.gis.setGroundBox(e0 ? { x0: e0.x0 - 40, y0: e0.y0 - 40, x1: e0.x1 + 40, y1: e0.y1 + 40 } : null);
  state.ui.shadowKey = null;
  let ext = null;
  if (model) {
    ext = addCity(state.gis, model, t);
    setBuildingForm(state.buildingForm);
    applyOpacity();
    const card = legendOfCity(dom, model, t);
    card.classList.add('city-card');
    state.ui.cityCard = card;
    state.ui.worldCard = null;
  } else {
    // No model here: the world's archives stand in, wherever this is.
    state.ui.worldCard = worldCard(dom, t);
    state.ui.cityCard = state.ui.worldCard.card;
  }
  // The first view frames the area's box if it names one (where its records are), else the whole model.
  const box = area.view && area.view.box;
  if (box) {
    const corners = [[box[0], box[1]], [box[0], box[3]], [box[2], box[1]], [box[2], box[3]]].map(([la, lo]) => state.frame.toLocal(la, lo).slice(0, 2));
    ext = state.gis.extentOf(corners);
  }
  ext = ext || { x0: -1500, y0: -1500, x1: 1500, y1: 1500 };
  state.ui.cityExtent = ext;
  state.gis.fit(ext, { ...(area.view || {}), skipView: view === 'keep', instant: view === 'instant', seconds: 4 });
  state.gis.onPick(onPick);
  state.gis.onHover(onHover);
  state.gis.onFrame(onFrame);
  state.spacetime.sceneChanged();
  state.clock = createClock(state.area, state.clock ? { date: state.clock.date, hour: state.clock.hour } : {});
  const applySun = (c) => {
    const s = c.sun;
    state.gis.setSun({ azimuth: (s.azimuth * Math.PI) / 180, altitude: (s.altitude * Math.PI) / 180 });
    const sub = subsolarPoint(c.instant);
    state.earth.setSubsolar(sub.lat, sub.lon);
  };
  state.clock.on(applySun);
  state.clock.on((c) => { for (const [, e] of state.layers) if (e.handle && e.handle.setTime) e.handle.setTime(c); });
  applySun(state.clock);
  state.ribbon.bind(state.clock);
  const { items, plan, works } = lensItems();
  state.ui.lenses.set(items, plan, works);
  document.getElementById('place-works')?.replaceChildren();
  viewKeys();
  history.replaceState(null, '', `#${area.id}`);
  for (const p of state.plugins.filter((x) => x.status !== 'planned'))
    for (const l of (p.layers || []).filter((x) => !x.areas || x.areas.includes(area.id))) {
      const key = `${p.name}/${l.id}`;
      const existing = state.layers.get(key);
      if (existing?.global && existing.handle && !keepLayers) continue;
      if (!inArea(l)) continue;
      if (existing?.on ?? !!l.default) toggleLayer(p, l, true, { quiet: true });
    }
  if (state.tool && inArea(state.tool.tool)) openTool(state.tool.plugin, state.tool.tool);
  else if (state.tool) closePanel();
  state.ui.neatline.invalidate();
  state.ui.citySeen = undefined;
  state.ui.opening = false;
}

// ---- pointer ------------------------------------------------------------------------------------------------------------

function onPick(hit, ev) {
  const at = ev ? state.earth.groundAt(ev.clientX, ev.clientY) : null;
  // A layer that takes the place itself (reach, say) comes before what stands there.
  if (at) for (const [, e] of state.layers) {
    if (e.handle?.pickAt && e.spec?.pickFirst) {
      const sel = e.handle.pickAt(at.lat, at.lon);
      if (sel) { if (!sel.at) sel.at = at; return select(sel); }
    }
  }
  if (hit && hit.layer) {
    let sel = null;
    if (hit.layer.startsWith('w:b:')) {
      const [m, code, name, kind] = hit.feature;
      const how = ['measured', 'from its floors', 'estimated: the average of the buildings on its 100 m cell, 2018', 'not known: a default'][code];
      sel = { kind: 'building', title: name || (kind && kind !== 'building' ? kind.replace(/_/g, ' ').replace(/^./, (c) => c.toUpperCase()) : 'A building'), lines: [`${Math.round(m)} m, ${how}`] };
    } else if (hit.layer === 'buildings') {
      const b = hit.feature;
      sel = { kind: 'building', title: b.name || `A building (${b.kind || 'yes'})`, lines: [`${Math.round(b.height_m)} m${b.height_source === 'default' ? ' (height not mapped: a default)' : b.height_source === 'levels' ? ` (${b.levels} levels)` : ''}`, b.id] };
    } else if (hit.layer === 'streets') {
      sel = { kind: 'way', title: hit.feature.name || hit.feature.kind, lines: [`${hit.feature.kind}, ${Math.round(hit.feature.width_m)} m wide${hit.feature.width_source === 'default' ? ' (a default width)' : ''}`] };
    } else {
      for (const [, e] of state.layers) if (e.handle && e.ids && e.ids.includes(hit.layer) && e.handle.describe) { sel = e.handle.describe(hit); break; }
    }
    if (sel) { if (at && !sel.at) sel.at = at; return select(sel); }
  }
  if (!at) return select(null);
  const inTime = state.spacetime.pickAt(at.lat, at.lon);
  if (inTime) return select(inTime);
  for (const [, e] of state.layers) {
    if (!e.handle || !e.handle.pickAt) continue;
    const sel = e.handle.pickAt(at.lat, at.lon);
    if (sel) { if (!sel.at) sel.at = at; return select(sel); }
  }
  select({ kind: 'place', title: 'Here', at });
}

function onHover(hit, ev) {
  const box = document.getElementById('hover');
  let text = null;
  if (hit && hit.layer) {
    const l = state.gis.layers.get(hit.layer);
    text = l && typeof l.style.pick === 'function' ? l.style.pick(hit.feature, hit.index) : null;
  }
  if (!text && ev) {
    const at = state.earth.groundAt(ev.clientX, ev.clientY);
    if (at) text = state.spacetime.hoverAt(at.lat, at.lon);
    if (at && !text) for (const [, e] of state.layers) if (e.handle && e.handle.hoverAt) { text = e.handle.hoverAt(at.lat, at.lon); if (text) break; }
  }
  if (!text || !ev) { box.hidden = true; return; }
  box.textContent = text;
  box.hidden = false;
  box.style.left = `${Math.min(ev.clientX + 14, window.innerWidth - 340)}px`;
  box.style.top = `${ev.clientY + 14}px`;
}

// ---- each frame ---------------------------------------------------------------------------------------------------------

const smooth = (a, b, x) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };
let titleAt = 0;
function onFrame(now) {
  const E = state.earth;
  E.frame(now);
  const v = E.view;
  // The city seen from near, fading out as the camera rises or leaves it.
  const a = state.area;
  const off = a ? Math.hypot((a.lat - v.lat / DEG) * 111e3, (a.lon - v.lon / DEG) * 111e3 * Math.cos(a.lat * DEG)) : Infinity;
  const seen = 1 - Math.max(smooth(70e3, 170e3, v.range), smooth(80e3, 220e3, off));
  const city = state.gis.renderer.domElement;
  if (state.ui.citySeen !== seen) {
    state.ui.citySeen = seen;
    city.style.opacity = seen.toFixed(3);
    state.gis.setPaused(seen <= 0.002);
    const far = seen <= 0.002;
    const flipped = state.ui.cityFar !== far;
    state.ui.cityFar = far;
    if (flipped && state.spacetime.needsRebuild()) state.spacetime.sceneChanged();
    E.city && (E.city.amount = seen);
    const legend = document.getElementById('legend');
    if (state.ui.cityCard) {
      if (seen > 0.5 && !state.ui.cityCard.isConnected) { legend.prepend(tidyLegend(state.ui.cityCard)); if (legend.children.length > 1) state.ui.cityCard.classList.add('folded'); }
      else if (seen <= 0.5 && state.ui.cityCard.isConnected) state.ui.cityCard.remove();
    }
    legend.classList.toggle('from-afar', seen <= 0.5);
    state.ui.lenses.update(v.range, true);
  }
  state.world.tick();
  if (now - titleAt > 120) {
    titleAt = now;
    credits(v);
    placeWorks();
    // A legend stands while its layer means something at this height.
    for (const [, e] of state.layers) {
      if (!e.legend || !e.spec?.scale) continue;
      const within = v.range >= e.spec.scale[0] * 0.5 && v.range <= e.spec.scale[1] * 1.7;
      e.legend.style.display = within ? '' : 'none';
    }
    placeTitle(v);
    state.ui.lenses.update(v.range);
    const lat = v.lat / DEG, lon = v.lon / DEG;
    if (!E.controls.flight && !state.ui.opening) {
      // Coming near another city of the workbench, it opens; coming down anywhere else, that place does.
      const near = v.range < 4e5 && state.areas.find((q) => q.id !== a?.id && metres(q.lat, q.lon, lat, lon) < Math.max(60e3, v.range));
      if (near) openArea(near, { view: 'keep' });
      else if (v.range < 3e4 && (!a || metres(a.lat, a.lon, lat, lon) > (a.generic ? 16e3 : 60e3))) openArea(placeAt(lat, lon), { view: 'keep', keepLayers: false });
    }
    if (!state.ui.opening && E.city) {
      state.world.want(lat, lon, v.range);
      const [ox, oy] = state.frame.toLocal(lat, lon, 0);
      state.gis.setDissolve(ox, oy, state.world.reach * (40075016.7 / 2 ** G.TILE_LEVEL) * Math.cos(v.lat) * 1.02);
      // The Sun's shadows follow what is looked at, as fine as it asks: stepping, so that they do not crawl.
      const s = Math.min(3600, Math.max(380, v.range * 1.2)), step = s / 5;
      const [x, y] = state.frame.toLocal(lat, lon, 0);
      const cx = Math.round(x / step) * step, cy = Math.round(y / step) * step, z = v.height - E.city.height, key = `${cx},${cy},${s.toFixed(0)}`;
      if (key !== state.ui.shadowKey && (!state.ui.shadowS || Math.abs(Math.log(s / state.ui.shadowS)) > 0.2 || state.ui.shadowCx !== cx || state.ui.shadowCy !== cy)) {
        state.ui.shadowKey = key, state.ui.shadowS = s, state.ui.shadowCx = cx, state.ui.shadowCy = cy;
        state.gis.setShadowBox({ x0: cx - s, y0: cy - s, x1: cx + s, y1: cy + s, z });
      }
    }
  }
  state.ui.navigator.update(v);
  const stage = document.getElementById('stage');
  state.ui.neatline.draw(stage.clientWidth, stage.clientHeight, { leftTop: 0 });
  updatePins(v);
}

// ---- the beginning ------------------------------------------------------------------------------------------------------

async function start() {
  theme();
  const mapEl = document.getElementById('map');
  try {
    const [ps, s, areas, lenses] = await Promise.all([plugins(), setup().catch(() => null), fetch('areas.json').then((r) => r.json()), fetch('lenses.json').then((r) => r.json())]);
    state.plugins = ps, state.setup = s, state.areas = areas.areas, state.lenses = lenses.lenses;
  } catch (error) {
    return dom.failure(mapEl, error);
  }
  state.plugins.sort((a, b) => (a.status === 'planned') - (b.status === 'planned') || (a.order ?? 99) - (b.order ?? 99));
  state.earth = G.createEarth(THREE, mapEl, { eventElement: mapEl, terrainUrl: (z, x, y) => `${location.origin}/tiles/terrain/earth/${z}/${x}/${y}`, preserveDrawingBuffer: /[?&]capture/.test(location.search) });
  applyLook();
  createWorldTiles();
  placesList();
  state.spacetime = createSpaceTimeHost({
    THREE, earth: state.earth, gis: () => state.gis, theme: themeNow, cityFar: () => state.ui.cityFar,
    area: { current: () => state.area, clock: () => state.clock },
    legendEl: () => document.getElementById('legend'), openLegend,
  });
  // A file dropped anywhere goes to the studio, which opens to take it.
  window.addEventListener('dragover', (e) => { if (e.dataTransfer?.types?.includes('Files')) { e.preventDefault(); document.body.classList.add('st-dragging'); } });
  window.addEventListener('dragleave', (e) => { if (!e.relatedTarget) document.body.classList.remove('st-dragging'); });
  window.addEventListener('drop', (e) => {
    const file = e.dataTransfer?.files?.[0];
    if (!file) return;
    e.preventDefault();
    document.body.classList.remove('st-dragging');
    if (state.spacetime.dropped(file)) return;
    const studio = state.plugins.find((p) => p.name === 'studio');
    if (studio) openTool(studio, studio.tools.find((t) => t.id === 'studio'));
  });
  state.years = createYears();
  state.ribbon = createRibbon(document.getElementById('clock'), { theme: themeNow });
  state.ui.navigator = createNavigator({ ladderEl: document.getElementById('ladder'), compassEl: document.getElementById('compass'), scaleEl: document.getElementById('scalebar'), earth: state.earth });
  state.ui.neatline = createNeatline(document.getElementById('neatline'), state.earth);
  state.ui.lenses = createLenses(document.getElementById('lenses'), {
    lensList: state.lenses, onToggle: onLayerRow, plannedPanel,
    onTool: onToolRow,
  });
  drawPins();
  lamps();
  document.getElementById('place-title').addEventListener('click', gazetteer);
  // Readings the plug-ins add to the probe, whatever layers are on.
  for (const p of state.plugins) if (p.reader) {
    import(`/plugins/${p.name}/${p.reader}`).then((m) => state.readers.push({ plugin: p, read: m.read })).catch((e) => console.warn(e));
  }
  const wanted = (location.hash.slice(1) || '').split('/')[0];
  const area = state.areas.find((a) => a.id === wanted);
  const there = /^@(-?[\d.]+),(-?[\d.]+)$/.exec(decodeURIComponent(wanted));
  if (area) await openArea(area, { view: 'instant' });
  else if (there) {
    // A place named by its coordinates: the camera is put above it, and it opens beneath.
    await placesList();
    const lat = +there[1], lon = +there[2];
    await openArea(state.areas[0], { view: 'keep' });
    state.earth.setView({ lon: lon * DEG, lat: lat * DEG, height: 0, range: 3600, tilt: 56 * DEG, heading: 0.35 });
  } else await openArea(state.areas[0], { view: 'keep' });
  window.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') { if (!document.getElementById('gazetteer').hidden) document.getElementById('gazetteer').hidden = true; else if (!document.getElementById('panel').hidden) closePanel(); else select(null); }
    if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'k') { e.preventDefault(); gazetteer(); }
  });
  // The page's state for its own checks and for its pictures.
  window.workbench = { state, goToArea, openArea, select, toggleLayer, openTool, closePanel };
}

start();
