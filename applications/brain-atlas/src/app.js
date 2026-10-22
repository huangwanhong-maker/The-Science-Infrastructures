// The Brain Atlas: a template brain in three slices and in the round, its
// parcellations as a tree, and data mapped onto it by region, by coordinate
// and from electrodes on the scalp.
//
// The application composes the neuroscience domain (domains/neuroscience/js)
// with the shared visualization layer (surfaces) and draws with the stage
// (src/render/stage.js). The data are the files the data library fetched into
// its cache (src/data/atlases.js).

import { Stage } from './render/stage.js';
import { loadTemplate, loadParcellation, regionIndices, MISSING_DATA, FILES, fetchBytes } from './data/atlases.js';
import { RegionTree } from './ui/tree.js';
import { valuesByRegion, parsePoints, channelValues, valueRange, colourOfValue, tidy } from './ui/data.js';
import { colourAt, categoryColour, hexToRgb } from './render/colormaps.js';
import { pointCard, parcellationPanel, dataPanel, headPanel, legend, aboutHtml, esc } from './ui/panel.js';
import {
  ATLASES, mniToTalairach, talairachToMni, worldToVoxel, maskScalp, tenTwenty, scaleToHead, sphericalSpline, apply, invert,
  readNifti, threshold, robustRange, volume as makeVolume, rasToCcf, ccfToRas, headFrame,
} from '../../../domains/neuroscience/js/index.js';

const $ = (id) => document.getElementById(id);

const state = {
  atlas: 'human', parcellation: 'julich-brain-3.1', layout: 'four',
  overlay: 0.55, selected: null, showBrain: true, showRegion: true,
  colormap: 'orchid', valuesText: '', valuesReport: '', values: null, valueRange: null, categories: null,
  pointsText: '', pointsReport: '', points: [], pointSpace: 'mni',
  head: { available: false, on: false, system: '10-20', measurements: {}, template: {}, scales: [1, 1, 1], channelsText: '', channelsReport: '', spline: null, positions: [], landmarks: null },
};

let stage, tree, template, entry, worker, requestId = 0;
const pending = new Map();
const entries = new Map();     // parcellation id -> loaded entry

// ---- loading ------------------------------------------------------------------------------

function progress(fraction, text) {
  $('loading-bar').style.width = `${Math.round(fraction * 100)}%`;
  if (text) $('loading-text').textContent = text;
}
function fail(html) { $('loading').classList.add('failed'); $('loading').classList.remove('done'); $('loading-text').innerHTML = html; }
function toast(text, ms = 2600) { const el = $('toast'); el.textContent = text; el.hidden = false; clearTimeout(toast.timer); toast.timer = setTimeout(() => { el.hidden = true; }, ms); }

function surface(kind, dims, affine, data, indices, step, smoothing = 2) {
  const id = ++requestId;
  return new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, kind, dims, affine, data, indices, step, smoothing });
  });
}

async function start() {
  worker = new Worker(new URL('./workers/surface.worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = (e) => { const p = pending.get(e.data.id); if (!p) return; pending.delete(e.data.id); if (e.data.error) p.reject(new Error(e.data.error)); else p.resolve(e.data); };
  const views = Object.fromEntries([...document.querySelectorAll('#views .view')].map((el) => [el.dataset.view, el]));
  stage = new Stage($('gl'), views);
  tree = new RegionTree($('tree'), { onSelect: selectRegion, colourOf: regionColour });
  restore();
  bindHeader();
  bindPanel();
  stage.on('point', onPoint);
  stage.on('hover', onHover);
  stage.on('render', drawOverlays);
  await loadAtlas(state.atlas);
}

async function loadAtlas(atlasId) {
  state.atlas = atlasId;
  state.selected = null; state.values = null; state.points = []; state.head.on = false; state.head.spline = null; state.head.positions = [];
  entries.clear();
  stage.setSurface('region', null); stage.setSurface('scalp', null); stage.setMarkers([]); stage.setScale([1, 1, 1]);
  $('loading').classList.remove('done', 'failed');
  progress(0.05, 'Reading the data library…');
  try {
    template = await loadTemplate(atlasId, progress);
  } catch (e) {
    console.error(e);
    fail(MISSING_DATA);
    return;
  }
  stage.setVolume(template.display, template.volume.dims, template.volume.affine);
  stage.setWindow(0.02, 0.98);
  $('credits').textContent = template.credit;
  document.querySelectorAll('.views [data-atlas]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.atlas === atlasId)));
  const available = ATLASES[atlasId].parcellations.filter((p) => FILES[atlasId].parcellations[p.id]);
  if (!available.some((p) => p.id === state.parcellation)) state.parcellation = available[0].id;
  // The head, for the human: measurements from the template's landmarks.
  state.head.available = Boolean(template.landmarks);
  if (template.landmarks) {
    const t = scaleToHead(template.landmarks, {}, template.headCentre).template;
    state.head.template = t;
    state.head.measurements = { ...t };
    state.head.scales = [1, 1, 1];
  }
  renderPanels();
  await loadParcellationById(state.parcellation);
  progress(0.95, 'The brain’s surface…');
  try {
    const mesh = await surface('mask', template.mask.dims, template.mask.affine, template.mask.data, null, atlasId === 'human' ? 2 : 1, 3);
    stage.setSurface('brain', mesh, { opacity: 0.28 });
    stage.surfaceVisible('brain', state.showBrain);
  } catch (e) { console.error(e); }
  $('loading').classList.add('done');
  if (state.restorePoint) { stage.setPoint(state.restorePoint); state.restorePoint = null; } else onPoint(stage.point);
  if (state.restoreSelected) { const r = entry.parcellation.get(state.restoreSelected); if (r) selectRegion(r); state.restoreSelected = null; }
}

async function loadParcellationById(id) {
  state.parcellation = id;
  if (!entries.has(id)) {
    $('loading').classList.remove('done');
    progress(0.3, 'Reading the parcellation…');
    try {
      entries.set(id, await loadParcellation(state.atlas, id, template, progress));
    } catch (e) {
      console.error(e);
      fail(`The parcellation could not be read: ${esc(e.message)}.<br><br>${MISSING_DATA}`);
      return;
    }
    $('loading').classList.add('done');
  }
  entry = entries.get(id);
  if (entry.unknownLabels) console.warn(`${entry.unknownLabels} labels of the map are not in the parcellation`);
  stage.setLabels(entry.index);
  state.selected = null;
  stage.setSurface('region', null);
  tree.setParcellation(entry.parcellation);
  state.values = null; state.valuesReport = '';
  updatePalette();
  renderPanels();
  onPoint(stage.point);
  save();
}

// ---- colours ---------------------------------------------------------------------------------

// A region's own colour, or the nearest ancestor's; a parcellation without
// colours (Julich-Brain's configuration has none) gets one hue per area, the
// two hemispheric parts alike, spaced by the golden angle and kept soft.
function regionColour(r) {
  for (let x = r; x; x = x.parent) if (x.colour) return x.colour;
  if (!r.isLeaf && r.depth < 2) return null;
  const area = r.hemisphere && r.parent ? r.parent : r;
  const hue = (area.index * 137.508) % 360;
  return hslToRgb(hue, 0.42, 0.64);
}

function hslToRgb(h, s, l) {
  const k = (n) => (n + h / 30) % 12;
  const a = s * Math.min(l, 1 - l);
  const f = (n) => l - a * Math.max(-1, Math.min(k(n) - 3, 9 - k(n), 1));
  return [f(0), f(8), f(4)].map((x) => Math.round(255 * x));
}

// The value a region takes from the mapped data: its own, or its nearest ancestor's.
function valueOf(r) { if (!state.values) return null; for (let x = r; x; x = x.parent) if (state.values.has(x)) return state.values.get(x); return null; }

function updatePalette() {
  const rgba = new Uint8Array(256 * 256 * 4);
  const selectedIndices = state.selected ? regionIndices(entry, state.selected) : null;
  for (let k = 1; k < entry.regions.length; k++) {
    const r = entry.regions[k];
    let colour = regionColour(r) || [160, 150, 170], alpha = 255;
    if (state.values) {
      const v = valueOf(r);
      if (v === null) { colour = [110, 104, 112]; alpha = 70; } else if (state.categories) { colour = hexToRgb(categoryColour(state.categories.get(v) ?? 0)); } else { colour = colourOfValue(v, state.valueRange, state.colormap); }
    }
    if (selectedIndices && selectedIndices.has(k)) { colour = colour.map((c) => Math.min(255, Math.round(c * 0.6 + 110))); alpha = 255; }
    else if (selectedIndices) { alpha = Math.round(alpha * 0.7); }
    rgba.set([Math.round(colour[0]), Math.round(colour[1]), Math.round(colour[2]), alpha], 4 * k);
  }
  stage.setPalette(rgba);
  renderLegend();
}

function renderLegend() {
  let spec = null;
  if (state.values) {
    spec = state.categories
      ? { title: state.valuesTitle || 'Mapped data', keys: [...state.categories.entries()].slice(0, 24).map(([name, i]) => ({ name: String(name), colour: categoryColour(i) })), note: `${state.values.size} regions` }
      : { title: state.valuesTitle || 'Mapped data', unit: state.valuesUnit || '', colormap: state.colormap, range: state.valueRange, note: `${state.values.size} regions coloured; the rest dimmed` };
  } else if (state.head.spline) {
    spec = { title: 'Scalp map', unit: state.head.unit || '', colormap: 'balance', range: state.head.range, note: `${state.head.channels.size} channels, spherical spline` };
  }
  $('legend').innerHTML = spec ? legend(spec) : '';
}

// ---- the point and the region ------------------------------------------------------------

function pointInfo(point) {
  const at = entry ? entry.map.at(point) : null;
  const region = at ? at.region : null;
  const voxel = worldToVoxel(template.volume, point);
  const info = { atlas: state.atlas, point, voxel, region, path: region ? region.path() : [], selected: state.selected, stats: null, talairach: null, ccf: null, value: region ? valueOf(region) : null, scaled: null, pointsHere: [], hemisphere: region ? region.side() : null };
  if (state.atlas === 'human') info.talairach = { brett: mniToTalairach(point), spm: mniToTalairach(point, { method: 'lancaster', software: 'spm' }) };
  if (state.atlas === 'mouse') info.ccf = rasToCcf(point);
  if (state.head.scales.some((s) => s !== 1)) info.scaled = apply(state.head.affine, point);
  if (state.selected) info.stats = entry.map.statistics(state.selected);
  info.pointsHere = state.points.filter((p) => Math.hypot(p.position[0] - point[0], p.position[1] - point[1], p.position[2] - point[2]) < 6);
  return info;
}

function onPoint(point, from) {
  if (!template) return;
  $('point').innerHTML = pointCard(pointInfo(point));
  const at = entry ? entry.map.at(point) : null;
  stage.setSelected(at && at.region ? entry.indexOf.get(at.region) ?? -1 : -1);
  if (from && from.marker && from.marker.region) selectRegion(from.marker.region);
  save();
}

function onHover(view, x, y) {
  const el = stage.views[view];
  let label = el.querySelector('.label.hover');
  if (x === null || view === 'three' || !entry) { if (label) label.remove(); return; }
  const p = stage.unproject(view, x, y);
  const at = p ? entry.map.at(p) : null;
  if (!at || !at.region) { if (label) label.remove(); return; }
  if (!label) { label = document.createElement('div'); label.className = 'label hover'; el.querySelector('.labels').appendChild(label); }
  label.textContent = at.region.name;
  label.style.left = `${x}px`; label.style.top = `${y - 14}px`;
}

async function selectRegion(region) {
  state.selected = region;
  tree.select(region);
  updatePalette();
  onPoint(stage.point);
  stage.setSurface('region', null);
  if (!region) return;
  const indices = [...regionIndices(entry, region)];
  if (!indices.length) return;
  try {
    const mesh = await surface('region', template.volume.dims, template.volume.affine, entry.index, indices, indices.length > 60 || state.atlas === 'mouse' ? 1 : 1, 2);
    if (state.selected !== region) return;
    const c = regionColour(region) || [246, 196, 126];
    stage.setSurface('region', mesh, { colour: `rgb(${c[0]}, ${c[1]}, ${c[2]})`, opacity: 1 });
    stage.surfaceVisible('region', state.showRegion);
  } catch (e) { console.error(e); }
}

// ---- data mapping ---------------------------------------------------------------------------

function mapValues(values, { title, unit, categories = null } = {}) {
  state.values = values; state.valuesTitle = title; state.valuesUnit = unit; state.categories = categories;
  state.valueRange = categories ? null : valueRange(values.values(), { diverging: state.colormap === 'balance' });
  updatePalette();
  onPoint(stage.point);
}

function builtin(kind) {
  const p = entry.parcellation;
  const values = new Map();
  if (kind === 'volume') {
    for (const r of p.leaves()) { const s = entry.map.statistics(r); if (s) values.set(r, s.volume / 1000); }
    mapValues(values, { title: 'Volume of each region', unit: 'cm³' });
    state.valuesReport = `<span class="good">${values.size} regions, from the voxels of the map.</span>`;
  } else if (kind === 'hemisphere') {
    const cats = new Map([['left', 0], ['right', 1], ['midline', 2]]);
    for (const r of p.leaves()) values.set(r, r.side() || 'midline');
    mapValues(values, { title: 'Hemisphere', categories: cats });
    state.valuesReport = '<span class="good">Each leaf by the side its name states.</span>';
  } else if (kind === 'depth') {
    for (const r of p.leaves()) values.set(r, r.depth);
    mapValues(values, { title: 'Depth in the tree', unit: 'levels' });
    state.valuesReport = '<span class="good">The depth of each leaf below the root.</span>';
  } else if (kind === 'group') {
    // The ancestor two levels below the root, as a category: a network, a lobe, a cerebral division.
    const cats = new Map();
    for (const r of p.leaves()) { const path = r.path(); const g = path[Math.min(path.length - 1, state.atlas === 'mouse' ? 2 : 0)]; if (!cats.has(g.name)) cats.set(g.name, cats.size); values.set(r, g.name); }
    mapValues(values, { title: 'Group', categories: cats });
    state.valuesReport = `<span class="good">${cats.size} groups.</span>`;
  }
  renderPanels();
}

function showPoints() {
  const { points, bad } = parsePoints(state.pointsText);
  const convert = (q) => {
    if (state.atlas === 'mouse') return state.pointSpace === 'ccf' ? ccfToRas(q) : q;
    if (state.pointSpace === 'mni') return q;
    if (state.pointSpace === 'brett') return talairachToMni(q);
    return talairachToMni(q, { method: 'lancaster', software: state.pointSpace });
  };
  state.points = points.map((p) => ({ ...p, given: p.position, position: convert(p.position), region: entry.map.regionAt(convert(p.position)) }));
  const range = valueRange(state.points.map((p) => p.value).filter((v) => v !== null));
  const markers = state.points.map((p, i) => ({ kind: 'datum', position: p.position, radius: state.atlas === 'mouse' ? 0.25 : 2.6, colour: p.value === null ? '#e8a0a0' : `rgb(${colourOfValue(p.value, range, state.colormap).join(',')})`, name: p.name || `${i + 1}`, region: p.region, value: p.value }));
  stage.setMarkers([...markers, ...electrodeMarkers()]);
  state.pointsReport = `<span class="good">${points.length} points shown</span>${bad.length ? `; <span class="bad">${bad.length} lines not read</span>` : ''}.` + (state.points.length ? `<br>${state.points.slice(0, 12).map((p) => `${esc(p.name || p.given.join(', '))}: ${p.region ? esc(p.region.name) : 'outside the parcellation'}`).join('<br>')}${state.points.length > 12 ? '<br>…' : ''}` : '');
  renderPanels();
  stage.invalidate();
}

// ---- the head ---------------------------------------------------------------------------------

let headMesh = null, scalpModel = null;

async function buildHead() {
  const h = state.head;
  if (!headMesh) {
    progress(0.2, 'Reading the template’s head…');
    $('loading').classList.remove('done');
    const spec = FILES.human.template;
    const head = readNifti(await fetchBytes(`${spec.dataset}/${spec.head}`, (f) => progress(0.2 + 0.6 * f, 'Reading the template’s head…'))).volume;
    const [lo, hi] = robustRange(head.data, 0.02, 0.98);
    const mask = makeVolume({ dims: head.dims, data: threshold(head, lo + 0.12 * (hi - lo)), affine: head.affine });
    headMesh = { mask };
    $('loading').classList.add('done');
  }
  // The scalp, scaled to the subject's head: cast through the mask with the scaled affine.
  const affine = h.affine || [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  const scaledMask = { ...headMesh.mask, affine: composeAffine(affine, headMesh.mask.affine) };
  const centre = apply(affine, template.headCentre);
  scalpModel = maskScalp(scaledMask, { centre, axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], step: 2, maxRadius: 220 });
  const landmarks = {};
  for (const [k, p] of Object.entries(template.landmarks)) landmarks[k] = scalpModel.onSurface(apply(affine, p));
  h.landmarks = landmarks;
  const system = tenTwenty(scalpModel, landmarks, { system: h.system });
  h.positions = system.positions;        // in the scaled space
  h.czDiscrepancy = system.czDiscrepancy;
  // The scalp surface for drawing, in atlas coordinates (the stage scales the group).
  try {
    const mesh = await surface('mask', headMesh.mask.dims, headMesh.mask.affine, headMesh.mask.data, null, 2, 3);
    colourScalp(mesh);
  } catch (e) { console.error(e); }
  stage.setMarkers([...pointMarkers(), ...electrodeMarkers()]);
  renderPanels();
}

const composeAffine = (A, B) => {
  const m = new Array(16).fill(0);
  for (let r = 0; r < 4; r++) for (let c = 0; c < 4; c++) { let s = 0; for (let k = 0; k < 4; k++) s += A[4 * r + k] * B[4 * k + c]; m[4 * r + c] = s; }
  return m;
};

// Markers for the electrodes, in atlas coordinates (the inverse of the head's scaling).
function electrodeMarkers() {
  const h = state.head;
  if (!h.on || !h.positions.length) return [];
  const inv = h.affine ? invert(h.affine) : null;
  const toAtlas = (p) => (inv ? apply(inv, p) : p);
  const out = h.positions.map((e) => ({ kind: 'electrode', name: e.name, position: toAtlas(e.position), radius: 3.2, colour: h.spline && h.channels.has(e.name) ? `rgb(${colourOfValue(h.channels.get(e.name), h.range, 'balance').join(',')})` : '#f6c47e' }));
  for (const [k, p] of Object.entries(h.landmarks || {})) out.push({ kind: 'landmark', name: { nasion: 'nasion', inion: 'inion', lpa: 'LPA', rpa: 'RPA', vertex: 'vertex' }[k], position: toAtlas(p), radius: 2.2, colour: '#e8c3dc' });
  return out;
}
function pointMarkers() { return (stage.markerList || []).filter((m) => m.kind === 'datum'); }

// The scalp mesh coloured by the spline (or cream), and the brain's surface likewise.
function colourScalp(mesh) {
  const h = state.head;
  const colours = h.spline ? vertexColours(mesh.positions) : null;
  stage.setSurface('scalp', mesh, { opacity: h.spline ? 0.9 : 0.16, colours, colour: '#f4dbc8' });
  stage.surfaceVisible('scalp', h.on);
  if (stage.surfaces.brain) {
    const brainColours = h.spline ? vertexColours(stage.surfaces.brain.geometry.attributes.position.array) : null;
    stage.setSurfaceColours('brain', brainColours, { opacity: brainColours ? 0.95 : 0.28 });
  }
  headMesh.mesh = mesh;
  stage.invalidate();
}

function vertexColours(positions) {
  const h = state.head;
  const inv = h.affine ? invert(h.affine) : null;
  const centre = inv ? apply(inv, scalpModel.centre) : scalpModel.centre;
  const out = new Float32Array(positions.length);
  for (let v = 0; v < positions.length; v += 3) {
    const d = [positions[v] - centre[0], positions[v + 1] - centre[1], positions[v + 2] - centre[2]];
    const value = h.spline.at(d);
    const [r, g, b] = colourOfValue(value, h.range, 'balance');
    out[v] = r / 255; out[v + 1] = g / 255; out[v + 2] = b / 255;
  }
  return out;
}

async function mapChannels() {
  const h = state.head;
  if (!h.on) { h.on = true; await buildHead(); }
  const { values, unknown } = channelValues(h.channelsText, h.positions.map((p) => p.name));
  if (values.size < 3) { h.channelsReport = '<span class="bad">At least three known channels are needed.</span>'; renderPanels(); return; }
  const dirs = [], vals = [];
  for (const [name, value] of values) { const e = h.positions.find((p) => p.name === name); dirs.push([e.position[0] - scalpModel.centre[0], e.position[1] - scalpModel.centre[1], e.position[2] - scalpModel.centre[2]]); vals.push(value); }
  h.channels = values;
  h.spline = sphericalSpline(dirs, vals, { regularisation: 1e-9 });
  h.range = valueRange(vals, { diverging: true });
  h.channelsReport = `<span class="good">${values.size} channels mapped</span>${unknown.length ? `; <span class="bad">unknown: ${esc(unknown.join(', '))}</span>` : ''}.`;
  if (headMesh && headMesh.mesh) colourScalp(headMesh.mesh);
  stage.setMarkers([...pointMarkers(), ...electrodeMarkers()]);
  renderLegend();
  renderPanels();
}

function applyHead(reset = false) {
  const h = state.head;
  if (reset) { h.measurements = { ...h.template }; }
  else for (const k of ['length', 'breadth', 'height']) { const v = Number($(`head-${k}`).value); if (v > 0) h.measurements[k] = v; }
  const s = scaleToHead(template.landmarks, h.measurements, template.headCentre);
  h.scales = s.scales; h.affine = s.affine;
  const earMid = [(template.landmarks.lpa[0] + template.landmarks.rpa[0]) / 2, (template.landmarks.lpa[1] + template.landmarks.rpa[1]) / 2, (template.landmarks.lpa[2] + template.landmarks.rpa[2]) / 2];
  stage.setScale(h.scales, template.headCentre || earMid);
  if (h.on) buildHead(); else renderPanels();
  onPoint(stage.point);
}

// ---- overlays -------------------------------------------------------------------------------

function drawOverlays() {
  for (const view of ['axial', 'coronal', 'sagittal', 'three']) {
    const el = stage.views[view];
    if (!el) continue;
    const r = stage.rect(view);
    const svg = el.querySelector('.overlay'), labels = el.querySelector('.labels');
    if (!r) continue;
    const parts = [], tags = [];
    if (view !== 'three') {
      const p = stage.project(view, stage.point);
      if (p) parts.push(`<line x1="${p.x}" y1="0" x2="${p.x}" y2="${r.h}"/><line x1="0" y1="${p.y}" x2="${r.w}" y2="${p.y}"/>`);
      // Markers within 6 mm of the slice plane, as rings.
      const axis = { axial: 2, coronal: 1, sagittal: 0 }[view];
      for (const m of stage.markerList || []) {
        const dz = Math.abs(m.position[axis] - stage.point[axis]);
        const reach = m.kind === 'datum' ? 6 : 5;
        if (dz > reach) continue;
        const q = stage.project(view, m.position);
        if (!q) continue;
        const radius = Math.max(2, (m.radius || 2.5) * 1.2 * Math.sqrt(1 - (dz / reach) ** 2) / (stage.mmPerPx ? stage.mmPerPx[view] : 1));
        parts.push(`<circle class="${m.kind === 'datum' ? 'datum' : 'point'}" cx="${q.x}" cy="${q.y}" r="${radius}" style="stroke:${m.colour}"/>`);
        tags.push(`<div class="label ${m.kind}" style="left:${q.x}px;top:${q.y}px">${esc(m.name)}</div>`);
      }
    } else {
      for (const m of stage.markerList || []) {
        if (m.kind === 'datum' && !m.name) continue;
        const q = stage.project(view, m.position);
        if (!q || !q.inside) continue;
        tags.push(`<div class="label ${m.kind}" style="left:${q.x}px;top:${q.y}px">${esc(m.name)}</div>`);
      }
    }
    svg.innerHTML = parts.join('');
    const hover = labels.querySelector('.label.hover');
    labels.innerHTML = tags.join('');
    if (hover) labels.appendChild(hover);
  }
}

// ---- the panel --------------------------------------------------------------------------------

function renderPanels() {
  const available = ATLASES[state.atlas].parcellations.filter((p) => FILES[state.atlas].parcellations[p.id] || p.reader === 'built-in');
  $('parcellations').innerHTML = parcellationPanel(state.atlas, available, state.parcellation, state.overlay);
  $('parcellations').querySelector('[data-action="toggle-brain"]').setAttribute('aria-pressed', String(state.showBrain));
  $('parcellations').querySelector('[data-action="toggle-region"]').setAttribute('aria-pressed', String(state.showRegion));
  $('data').innerHTML = dataPanel(state, entry ? entry.parcellation.name : '');
  $('head').innerHTML = headPanel(state.head);
}

function bindHeader() {
  document.querySelectorAll('.views [data-atlas]').forEach((b) => b.addEventListener('click', () => { if (b.dataset.atlas !== state.atlas) { state.parcellation = null; loadAtlas(b.dataset.atlas); } }));
  document.querySelectorAll('.layouts [data-layout]').forEach((b) => b.addEventListener('click', () => setLayout(b.dataset.layout)));
  $('data-link').addEventListener('click', () => { $('about-body').innerHTML = aboutHtml(state.atlas, ATLASES[state.atlas].parcellations); $('about').showModal(); });
  $('region-search').addEventListener('input', (e) => tree.search(e.target.value));
  window.addEventListener('keydown', (e) => {
    if (e.target.matches('input, textarea, select')) return;
    const step = e.shiftKey ? 5 : 1;
    const p = stage.point.slice();
    const moves = { ArrowLeft: [0, -step], ArrowRight: [0, step], ArrowUp: [1, step], ArrowDown: [1, -step], PageUp: [2, step], PageDown: [2, -step] };
    if (moves[e.key]) { e.preventDefault(); p[moves[e.key][0]] += moves[e.key][1]; stage.setPoint(p); }
    if (e.key === 'Escape') selectRegion(null);
  });
}

function setLayout(layout) {
  state.layout = layout;
  $('views').className = layout === 'four' ? '' : layout;
  document.querySelectorAll('.layouts [data-layout]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.layout === layout)));
  stage.setLayout(layout);
  save();
}

function bindPanel() {
  const panel = $('panel');
  panel.addEventListener('change', (e) => {
    const t = e.target;
    if (t.name === 'parcellation') loadParcellationById(t.value);
    if (t.id === 'overlay') { state.overlay = Number(t.value); $('overlay-value').textContent = `${Math.round(state.overlay * 100)}%`; stage.setOverlay(state.overlay); }
    if (t.id === 'colormap') { state.colormap = t.value; if (state.values && !state.categories) { state.valueRange = valueRange(state.values.values(), { diverging: state.colormap === 'balance' }); } updatePalette(); }
    if (t.id === 'points-space') state.pointSpace = t.value;
  });
  panel.addEventListener('input', (e) => {
    const t = e.target;
    if (t.id === 'overlay') { state.overlay = Number(t.value); $('overlay-value').textContent = `${Math.round(state.overlay * 100)}%`; stage.setOverlay(state.overlay); }
    if (t.id === 'values-text') state.valuesText = t.value;
    if (t.id === 'points-text') state.pointsText = t.value;
    if (t.id === 'channels-text') state.head.channelsText = t.value;
  });
  panel.addEventListener('click', async (e) => {
    const b = e.target.closest('button');
    if (!b) return;
    const action = b.dataset.action;
    if (b.dataset.builtin) builtin(b.dataset.builtin);
    else if (b.dataset.system) { state.head.system = b.dataset.system; if (state.head.on) await buildHead(); else renderPanels(); }
    else if (action === 'map-values') {
      const { values, unmatched } = valuesByRegion(entry.parcellation, state.valuesText);
      if (!values.size) { state.valuesReport = '<span class="bad">No region of this parcellation matched a line.</span>'; renderPanels(); return; }
      mapValues(values, { title: 'Pasted values' });
      state.valuesReport = `<span class="good">${values.size} regions matched</span>${unmatched.length ? `; <span class="bad">not matched: ${esc(unmatched.slice(0, 8).join(', '))}${unmatched.length > 8 ? ', …' : ''}</span>` : ''}.`;
      renderPanels();
    } else if (action === 'clear-values') { state.values = null; state.categories = null; state.valuesReport = ''; updatePalette(); onPoint(stage.point); renderPanels(); }
    else if (action === 'show-points') showPoints();
    else if (action === 'clear-points') { state.points = []; state.pointsReport = ''; stage.setMarkers(electrodeMarkers()); renderPanels(); onPoint(stage.point); }
    else if (action === 'centre' && state.selected) { const s = entry.map.statistics(state.selected); if (s) stage.setPoint(s.centroid); }
    else if (action === 'deselect') selectRegion(null);
    else if (action === 'toggle-brain') { state.showBrain = !state.showBrain; stage.surfaceVisible('brain', state.showBrain); renderPanels(); }
    else if (action === 'toggle-region') { state.showRegion = !state.showRegion; stage.surfaceVisible('region', state.showRegion); renderPanels(); }
    else if (action === 'toggle-head') { state.head.on = !state.head.on; if (state.head.on) await buildHead(); else { stage.surfaceVisible('scalp', false); stage.setMarkers(pointMarkers()); renderPanels(); } }
    else if (action === 'apply-head') applyHead(false);
    else if (action === 'reset-head') applyHead(true);
    else if (action === 'map-channels') mapChannels();
    else if (action === 'clear-channels') { state.head.spline = null; state.head.channelsReport = ''; if (headMesh && headMesh.mesh) colourScalp(headMesh.mesh); stage.setMarkers([...pointMarkers(), ...electrodeMarkers()]); renderLegend(); renderPanels(); }
  });
}

// ---- the address bar ----------------------------------------------------------------------

function save() {
  if (!template) return;
  const q = new URLSearchParams();
  q.set('atlas', state.atlas); q.set('parcellation', state.parcellation); q.set('layout', state.layout);
  q.set('point', stage.point.map((x) => Math.round(x * 10) / 10).join(','));
  if (state.selected) q.set('region', state.selected.name);
  history.replaceState(null, '', `#${q.toString()}`);
}

function restore() {
  const q = new URLSearchParams(location.hash.slice(1));
  if (q.get('atlas') && ATLASES[q.get('atlas')]) state.atlas = q.get('atlas');
  if (q.get('parcellation')) state.parcellation = q.get('parcellation');
  if (q.get('layout')) { state.layout = q.get('layout'); $('views').className = state.layout === 'four' ? '' : state.layout; document.querySelectorAll('.layouts [data-layout]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.layout === state.layout))); }
  const p = (q.get('point') || '').split(',').map(Number);
  if (p.length === 3 && p.every(Number.isFinite)) state.restorePoint = p;
  if (q.get('region')) state.restoreSelected = q.get('region');
}

start().catch((e) => { console.error(e); fail(`The atlas could not start: ${esc(e.message)}`); });
