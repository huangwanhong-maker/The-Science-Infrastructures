// The Brain Atlas: the human brain in the round, every part of it clickable for
// its anatomy; sections through it; atlases and data mapped onto it; the head and
// its electrodes around it.
//
// The anatomy is the data library's product (data/products/brain_anatomy.py); the
// names, the tree of structures and what is known of each are the neuroscience
// domain's (brain.js, knowledge.js); the drawing is the viewer's (scene/viewer.js)
// and the sections' (sections/sections.js).

import { Viewer } from './scene/viewer.js';
import { Sections, SectionView, PLANES } from './sections/sections.js';
import { loadProduct, loadPack, loadTemplate, loadLabels, MISSING } from './data/anatomy.js';
import { StructureTree } from './ui/tree.js';
import { cardHtml, esc, volumeText } from './ui/card.js';
import { layoutPlate, clearPlate, samplesOf } from './ui/plate.js';
import { drawCompass, toolbarHtml, searchResults, resultsHtml, VIEW_DIRECTIONS, VIEW_UPS } from './ui/controls.js';
import { TOURS } from './ui/tours.js';
import { AtlasesMode } from './modes/atlases.js';
import { DataMode } from './modes/data.js';
import { HeadMode } from './modes/head.js';
import {
  BRAIN_STRUCTURES, structure, childrenOf, ancestorsOf, descendantsOf, keysOf, structureOfKey, LOBE_COLOURS, lobeOf,
  SYSTEMS, system, knowledgeOf, referencesOf, mniToTalairach,
} from '../../../domains/neuroscience/js/index.js';

const $ = (id) => document.getElementById(id);

// Structures drawn by default are the gross ones; their finer divisions wait to be asked for.
const HIDDEN_AT_FIRST = new Set(['ca1', 'ca23', 'dentate-gyrus', 'subiculum', 'presubiculum', 'uncus',
  'thalamus-ventral-anterior', 'thalamus-ventrolateral', 'thalamus-ventral-posterior-ventrolateral', 'thalamus-lateral-posterior-ventral-posterior',
  'thalamus-mediodorsal', 'thalamus-pulvinar', 'thalamus-pulvinar-medial-centrolateral', 'globus-pallidus-externus', 'globus-pallidus-internus',
  'internal-capsule', 'cerebral-white-matter']);
const SHELLS = new Set(['cerebral-cortex', 'cerebral-white-matter', 'cerebellar-cortex', 'brainstem', 'cerebellar-white-matter']);
const NATURAL = '#e9cbc4';

const state = {
  mode: 'anatomy', colouring: 'gyri', cortexStyle: 'solid', cut: null, explode: 0, labels: true,
  showWhite: false, showBundles: false, selected: null, hover: null, isolate: false, system: null,
  point: [0, -18, 10], mainPlane: 'coronal', overlay: 0.55,
};
const app = {
  state, viewer: null, product: null, sections: null, views: {}, tree: null,
  instances: new Map(),        // structure id -> [instance]
  hidden: new Set(HIDDEN_AT_FIRST),
  palettes: {},
  viewHooks: [],               // functions called when the view moves (kind: 'change' | 'end')
  labelValues: new Map(),      // joined label value -> key
  cortexColour: null,          // (inst, vertex) -> [r, g, b] when a mode colours the cortex
  cortexRegion: null,          // (inst, vertex) -> region descriptor, for atlases
};
export default app;
window.brainAtlas = app;       // for inspection from the console
const INITIAL_HASH = location.hash;   // the address as the page was opened, before any view is saved

// ---- loading ------------------------------------------------------------------------------------

function progress(f, text) { $('loading-bar').style.width = `${Math.round(f * 100)}%`; if (text) $('loading-text').textContent = text; }
function fail(html) { $('loading').classList.add('failed'); $('loading').classList.remove('done'); $('loading-text').innerHTML = html; }
export function toast(text, ms = 2600) { const el = $('toast'); el.textContent = text; el.hidden = false; clearTimeout(toast.t); toast.t = setTimeout(() => { el.hidden = true; }, ms); }

async function start() {
  app.viewer = new Viewer($('gl'), { onHover, onPick, onView });
  app.viewer.screenBox = screenBox;
  app.viewer.afterRender = () => drawCompass($('compass'), app.viewer.basis());
  bindChrome();
  try {
    progress(0.04, 'Reading the catalogue of the anatomy…');
    app.product = await loadProduct();
  } catch (e) {
    console.error(e);
    fail(MISSING);
    return;
  }
  for (const [value, key] of app.product.labelTable) app.labelValues.set(value, key);
  progress(0.1, 'Reading the surfaces of the brain…');
  const core = await loadPack(app.product, 'core', (f) => progress(0.1 + 0.7 * f));
  progress(0.82, 'Building the surfaces…');
  addPack('core', core);
  restore();
  app.tree = new StructureTree(document.createElement('div'), { nodes: treeNodes(), isVisible: treeVisible, onSelect: (id) => select({ id, side: null }), onToggle: toggleStructure });
  renderLeft();
  paintAll();
  applyStyles();
  computeExplodeOffsets();
  setInitialView();
  renderToolbar();
  renderCard();
  $('loading').classList.add('done');
  app.ready = true;
  // The sections follow, in the background.
  try {
    const [t1, labels] = await Promise.all([loadTemplate(), loadLabels(app.product)]);
    app.sections = new Sections({ t1, labels, palette: labelPalette() });
    app.sections.overlay = state.overlay;
    setupSectionViews();
    if (state.selected) renderCard();
    if (state.cut) updateCut();
  } catch (e) { console.error('sections', e); }
  if (app.pendingMode) setMode(app.pendingMode);
}

function addPack(pack, bytes) {
  const scale = app.product.positionScale;
  for (const entry of app.product.meshes) {
    if (entry.pack !== pack) continue;
    const at = structureOfKey(entry.key, 'mesh');
    if (!at) continue;
    const s = structure(at.id);
    const kind = entry.key.startsWith('cortex:') ? 'cortex' : s.kind === 'ventricle' ? 'ventricle' : s.kind === 'tract' ? 'tract' : 'tissue';
    const inst = app.viewer.addMesh(entry, bytes, scale, { kind, id: s.id, side: at.side, colour: s.colour || '#d9c3cf', clipped: SHELLS.has(s.id) });
    inst.region = app.product.regions[entry.key] || null;
    if (!app.instances.has(s.id)) app.instances.set(s.id, []);
    app.instances.get(s.id).push(inst);
  }
}

async function ensurePack(pack) {
  if (app.loadedPacks?.has(pack)) return;
  app.loadedPacks = app.loadedPacks || new Set(['core']);
  app.loadedPacks.add(pack);
  toast(pack === 'white' ? 'Reading the white matter…' : 'Reading the fibre bundles…');
  const bytes = await loadPack(app.product, pack);
  addPack(pack, bytes);
  computeExplodeOffsets();
  paintAll();
  applyStyles();
  app.tree.render();
}

// ---- the tree of structures and its visibility ---------------------------------------------------------

function treeNodes() {
  const k = (id) => knowledgeOf(id).entry;
  return BRAIN_STRUCTURES.filter((s) => s.id !== 'brain' && s.id !== 'other' && !['vessel', 'white-matter-hypointensities', 'csf'].includes(s.id))
    .map((s) => ({ id: s.id, name: s.name, parent: s.parent, colour: colourOf(s.id), also: (k(s.id) && k(s.id).also) || [] }));
}

// Structures with a surface somewhere beneath them can be shown or hidden.
function drawable(id) { return [id, ...descendantsOf(id).map((d) => d.id)].some((x) => app.instances.has(x) || isGyrus(x)); }
const isGyrus = (id) => { const s = structure(id); return Boolean(s && s.label && s.label.left && s.label.left.startsWith('ho:')); };

function treeVisible(id) {
  if (!drawable(id)) return null;
  if (isGyrus(id) || ancestorsOf(id).some((a) => a.id === 'cerebral-cortex') || id === 'cerebral-cortex') return state.cortexStyle !== 'hidden' && !app.hidden.has('cerebral-cortex');
  const own = app.instances.has(id) ? !app.hidden.has(id) : false;
  const kids = descendantsOf(id).filter((d) => app.instances.has(d.id));
  if (!app.instances.has(id)) return kids.some((d) => !app.hidden.has(d.id));
  return own || kids.some((d) => !app.hidden.has(d.id));
}

async function toggleStructure(id) {
  if (id === 'cerebral-white-matter' || descendantsOf('fibre-bundles').some((d) => d.id === id) || id === 'fibre-bundles') {
    await ensurePack(id === 'cerebral-white-matter' ? 'white' : 'bundles');
    if (id === 'cerebral-white-matter') state.showWhite = !treeVisible(id) || app.hidden.has('cerebral-white-matter');
  }
  if (isGyrus(id) || id === 'cerebral-cortex' || ancestorsOf(id).some((a) => a.id === 'cerebral-cortex')) {
    state.cortexStyle = state.cortexStyle === 'hidden' ? 'solid' : 'hidden';
    app.hidden.delete('cerebral-cortex');
  } else {
    const ids = [id, ...descendantsOf(id).map((d) => d.id)].filter((x) => app.instances.has(x));
    const shown = treeVisible(id);
    if (shown) ids.forEach((x) => app.hidden.add(x));
    else {
      // Showing a structure shows it, not its finer divisions, unless it has no surface of its own.
      if (app.instances.has(id)) app.hidden.delete(id);
      else ids.forEach((x) => { if (!HIDDEN_AT_FIRST.has(x) || structure(x).parent === id) app.hidden.delete(x); });
      if (app.instances.has(id)) for (const a of ancestorsOf(id)) if (app.instances.has(a.id) && !SHELLS.has(a.id)) { /* the parent is drawn as glass */ }
    }
    if (id === 'cerebral-white-matter') state.showWhite = !app.hidden.has(id);
    if (id === 'fibre-bundles' || descendantsOf('fibre-bundles').some((d) => d.id === id)) state.showBundles = descendantsOf('fibre-bundles').some((d) => app.instances.has(d.id) && !app.hidden.has(d.id));
  }
  applyStyles();
  app.tree.render();
  renderToolbar();
  schedulePlate();
}

// ---- colours ----------------------------------------------------------------------------------

function colourOf(id) {
  const s = structure(id);
  if (!s) return null;
  if (s.colour) return s.colour;
  const kids = childrenOf(id).map((c) => colourOf(c.id)).filter(Boolean);
  return kids[0] || null;
}

const rgbOf = (hex) => [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255);

// The cortex's colour at each Harvard-Oxford label, by the colouring chosen.
function cortexPalette(colouring) {
  const pal = new Map();
  for (let v = 1; v <= 96; v++) {
    const at = structureOfKey(`ho:${v}`);
    if (!at) continue;
    const hex = colouring === 'lobes' ? LOBE_COLOURS[lobeOf(at.id)] : colouring === 'natural' ? NATURAL : structure(at.id).colour;
    pal.set(v, rgbOf(hex || NATURAL));
  }
  pal.set(0, rgbOf(colouring === 'natural' ? NATURAL : '#e2c9c4'));
  return pal;
}

// The palette of the joined labels for the sections: colour and overlay strength per value.
function labelPalette() {
  const own = app.modes[state.mode]?.sectionPalette?.();
  if (own) return own;
  const n = app.product.labelTable.length;
  const pal = new Uint8Array(4 * (n + 1));
  for (const [value, key] of app.product.labelTable) {
    const at = structureOfKey(key);
    const s = at && structure(at.id);
    let hex = (s && colourOf(s.id)) || '#cccccc';
    let alpha = 255;
    if (s && (s.id === 'cerebral-white-matter' || s.id === 'cerebellar-white-matter')) alpha = 0;
    if (s && s.id === 'cerebral-cortex') alpha = 110;
    if (s && (s.id === 'csf' || s.id === 'vessel' || s.id === 'white-matter-hypointensities')) alpha = 0;
    if (state.mode === 'anatomy' && state.colouring === 'lobes' && s && isGyrus(s.id)) hex = LOBE_COLOURS[lobeOf(s.id)];
    const [r, g, b] = rgbOf(hex).map((x) => Math.round(255 * x));
    pal.set([r, g, b, alpha], 4 * value);
  }
  return pal;
}

// ---- painting and styles ----------------------------------------------------------------------------

function selectionKeys() {
  const sel = state.selected;
  if (!sel) return null;
  if (sel.system) {
    const out = new Set();
    for (const id of system(sel.system).members) for (const side of ['left', 'right', 'mid']) {
      for (const which of ['label', 'mesh']) for (const k of keysOf(id, which)[side]) out.add(k);
    }
    return out;
  }
  const out = new Set();
  for (const which of ['label', 'mesh']) {
    const k = keysOf(sel.id, which);
    for (const side of ['left', 'right', 'mid']) if (!sel.side || side === sel.side || side === 'mid') for (const x of k[side]) out.add(x);
  }
  return out;
}

function paintCortex() {
  const sel = selectionKeys();
  const hover = state.hover && state.hover.key;
  const pal = app.cortexColour ? null : cortexPalette(state.colouring);
  for (const inst of app.instances.get('cerebral-cortex') || []) {
    const labels = inst.labels;
    const side = inst.side;
    const lit = new Set();
    let any = false;
    if (sel) for (let v = 1; v <= 96; v++) if (sel.has(`ho:${v}`)) { lit.add(v); any = true; }
    const whole = sel && sel.has(inst.key);
    const hoverLabel = hover && hover.startsWith('ho:') ? Number(hover.slice(3)) : -1;
    const regionLit = app.regionLit;            // a set of vertex indices lit by an atlas region
    const emphasis = (i) => {
      const l = labels[i];
      if (regionLit) return regionLit.has(inst.key) ? (regionLit.get(inst.key)[i] ? 1.22 : 0.62) : 0.62;
      if (any || (sel && !whole)) return lit.has(l) ? 1.2 : 0.62;
      return l === hoverLabel ? 1.14 : 1;
    };
    const colour = app.cortexColour ? (i) => app.cortexColour(inst, i) : (i) => pal.get(labels[i]) || pal.get(0);
    // A fine engraved line along the borders of the gyri (or of the regions an atlas draws).
    const engrave = state.colouring !== 'natural' || app.cortexColour;
    if (engrave && !inst.borderFlags) inst.borderFlags = app.viewer.borders(inst, app.cortexKey ? (i) => app.cortexKey(inst, i) : undefined);
    const flags = engrave ? inst.borderFlags : null;
    const withBorders = flags ? (i) => emphasis(i) * (flags[i] ? 0.8 : 1) : emphasis;
    app.viewer.paint(inst, colour, { emphasis: withBorders });
    void side;
  }
}

function paintAll() {
  for (const [id, list] of app.instances) {
    if (id === 'cerebral-cortex') continue;
    for (const inst of list) {
      inst.colour = (app.meshColour && app.meshColour(inst)) || colourOf(id) || inst.colour;
      app.viewer.paint(inst);
    }
  }
  paintCortex();
}

// The style of every instance: drawn, glass, ghost or hidden, and its emphasis.
function applyStyles() {
  const sel = selectionKeys();
  const focus = (state.isolate || state.system) && sel;
  for (const [id, list] of app.instances) {
    for (const inst of list) {
      let visible = !app.hidden.has(id) && !app.otherSpecies;
      let style = 'solid';
      if (id === 'cerebral-cortex') {
        visible = state.cortexStyle !== 'hidden' && !app.hidden.has(id) && !app.otherSpecies;
        style = state.cortexStyle === 'glass' ? 'glass' : 'solid';
        if (focus && !sel.has(inst.key) && ![...sel].some((k) => k.startsWith('ho:'))) style = 'glass';
      } else {
        // A structure whose finer divisions are shown is seen through.
        const kids = descendantsOf(id).filter((d) => app.instances.has(d.id) && !app.hidden.has(d.id));
        if (kids.length) style = 'ghost';
        if (focus && !sel.has(inst.key)) style = 'ghost';
        if (focus && sel.has(inst.key) && !app.otherSpecies) visible = true;
      }
      if (inst.mesh.visible !== visible || inst.visible !== visible) app.viewer.setVisible(inst, visible);
      if (inst.style !== style) app.viewer.setStyle(inst, style);
      const selected = sel && sel.has(inst.key) && id !== 'cerebral-cortex';
      const hovered = state.hover && state.hover.key === inst.key && id !== 'cerebral-cortex';
      app.viewer.setEmphasis(inst, selected ? 0.22 : hovered ? 0.12 : 0, selected ? '#f6d39a' : null);
    }
  }
  app.viewer.dirty = true;
}

// ---- hovering and choosing ------------------------------------------------------------------

// The structure, side and key of a picked surface (a gyrus, for the cortex).
function whatIs(hit) {
  if (!hit) return null;
  const inst = hit.inst;
  if (app.cortexRegion && inst.id === 'cerebral-cortex') {
    const region = app.cortexRegion(inst, hit);
    if (region) return { ...region, inst, point: hit.point };
  }
  if (inst.labels && hit.label) {
    const key = `ho:${hit.label}`;
    const at = structureOfKey(key);
    if (at) return { id: at.id, side: at.side, key, inst, point: hit.point };
  }
  return { id: inst.id, side: inst.side, key: inst.key, inst, point: hit.point };
}

function onHover(hit, pointer) {
  const tip = $('tip');
  // A mode's markers (electrodes, points) come before the surfaces.
  const info = pointer && app.modes[state.mode]?.hoverInfo?.(pointer);
  if (info) {
    const r = $('viewport').getBoundingClientRect();
    tip.style.left = `${pointer.x - r.left}px`; tip.style.top = `${pointer.y - r.top}px`;
    tip.innerHTML = `<div class="t"><i style="background:${esc(info.colour || '#ccc')}"></i>${esc(info.name)}</div><div class="s">${esc(info.sub || '')}</div>`;
    tip.hidden = false;
    $('gl').classList.add('pointing');
    return;
  }
  const w = whatIs(hit);
  const key = w ? w.key : null;
  $('gl').classList.toggle('pointing', Boolean(w));
  if (!w || !pointer) { tip.hidden = true; if (state.hover) { state.hover = null; applyStyles(); paintCortex(); } return; }
  const r = $('viewport').getBoundingClientRect();
  tip.style.left = `${pointer.x - r.left}px`; tip.style.top = `${pointer.y - r.top}px`;
  const name = w.name || structure(w.id).name;
  const sub = w.sub || [w.side && w.side !== 'mid' ? `${w.side} hemisphere` : '', ancestorsOf(w.id).filter((a) => a.id !== 'brain')[0]?.name].filter(Boolean).join(' · ');
  tip.innerHTML = `<div class="t"><i style="background:${esc(w.colour || colourOf(w.id) || '#ccc')}"></i>${esc(name)}</div><div class="s">${esc(sub)}</div>`;
  tip.hidden = false;
  if (!state.hover || state.hover.key !== key) {
    const was = state.hover;
    state.hover = { key, id: w.id };
    applyStyles();
    if ((was && was.key && was.key.startsWith('ho:')) || key.startsWith('ho:')) paintCortex();
  }
}

function onPick(hit, event, double) {
  if (app.modes[state.mode]?.pick?.(event)) return;
  const w = whatIs(hit);
  if (!w) { if (!double) select(null); return; }
  if (w.region) { app.modes[state.mode]?.selectRegion?.(w.region, w); return; }
  select({ id: w.id, side: w.side === 'mid' ? null : w.side }, { point: w.point });
  if (double) focus();
}

// Whether a structure lies within the brain, hidden by a solid cortex.
function isDeep(id) {
  if (!id || SHELLS.has(id) || id === 'cerebral-cortex' || isGyrus(id)) return false;
  return !ancestorsOf(id).some((a) => a.id === 'cerebral-cortex');
}

export function select(sel, { point = null } = {}) {
  state.selected = sel;
  // A deep structure chosen while the cortex is solid is shown through glass.
  if (sel && sel.id && isDeep(sel.id) && state.cortexStyle === 'solid' && state.mode === 'anatomy' && app.instances.size && !app.tour) {
    state.cortexStyle = 'glass';
    renderToolbar();
  }
  state.system = sel && sel.system ? sel.system : null;
  if (!sel || !sel.system) state.isolate = state.isolate && Boolean(sel);
  app.regionLit = null;
  if (point) setPoint(point, { quiet: true });
  else if (sel && sel.id && app.instances.size) {
    const c = centroidOf(sel);
    if (c) setPoint(c, { quiet: true });
  }
  if (app.tree) app.tree.select(sel && sel.id ? sel.id : null);
  renderSystems();
  applyStyles();
  paintCortex();
  renderCard();
  updateSectionSelection();
  schedulePlate();
  save();
}

function centroidOf(sel) {
  const keys = keysOf(sel.id, 'label');
  const sides = sel.side ? [sel.side, 'mid'] : ['left', 'right', 'mid'];
  let n = 0;
  const c = [0, 0, 0];
  for (const side of sides) for (const k of keys[side]) {
    const r = app.product.regions[k];
    if (!r) continue;
    const v = r.joined ? r.joined : r;
    for (let a = 0; a < 3; a++) c[a] += v.centroid[a] * v.voxels;
    n += v.voxels;
  }
  if (!n) return null;
  return c.map((x) => x / n);
}

function boxOf() {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  const keys = selectionKeys();
  let found = false;
  for (const k of keys || []) {
    const r = app.product.regions[k];
    if (!r) continue;
    const b = r.joined ? r.joined.bounds : r.bounds;
    for (let a = 0; a < 3; a++) { min[a] = Math.min(min[a], b.min[a]); max[a] = Math.max(max[a], b.max[a]); }
    found = true;
  }
  return found ? { min, max } : null;
}

async function focusBox(min, max, margin = 1.6) {
  const { Box3, Vector3 } = await import('three');
  app.viewer.frame(new Box3(new Vector3(...min), new Vector3(...max)), { margin });
}

async function focus() {
  const b = boxOf();
  if (!b) return;
  const { Box3, Vector3 } = await import('three');
  const box = new Box3(new Vector3(...b.min), new Vector3(...b.max));
  app.viewer.frame(box, { margin: 1.6 });
}

// ---- the point and the sections -----------------------------------------------------------------

function setPoint(p, { quiet = false } = {}) {
  state.point = p.map((x) => Math.round(x * 10) / 10);
  for (const v of Object.values(app.views)) { v.point = state.point; v.redraw(); }
  if (!quiet) save();
}

function setupSectionViews() {
  const els = [...document.querySelectorAll('#sections .section')];
  const flip = (plane) => plane === 'sagittal';
  for (const el of els) {
    const plane = el.dataset.plane;
    const v = new SectionView(el, plane, app.sections, {
      flipU: flip(plane),
      onPoint: (w, how) => {
        setPoint(w);
        const mode = app.modes[state.mode];
        if (how.click && mode?.sectionAt) { const at = mode.sectionAt(w); if (at) at.choose(); return; }
        if (how.click) {
          const L = app.sections.labelAt(w);
          const key = app.labelValues.get(L);
          const at = key && structureOfKey(key);
          if (at && !['cerebral-white-matter', 'csf', 'vessel', 'white-matter-hypointensities'].includes(at.id)) select({ id: at.id, side: at.side === 'mid' ? null : at.side }, { point: w });
        }
      },
      onHover: (w, e) => sectionHover(w, e),
      onSwap: (plane2) => { state.mainPlane = plane2; arrangeSections(); renderToolbar(); },
    });
    v.point = state.point;
    app.views[el.dataset.slot = el.classList.contains('main') ? 'main' : el.dataset.plane] = v;
  }
  arrangeSections();
}

function arrangeSections() {
  const planes = ['axial', 'coronal', 'sagittal'];
  const main = state.mainPlane;
  const others = planes.filter((p) => p !== main);
  const els = [...document.querySelectorAll('#sections .section')];
  const order = [main, ...others];
  const views = Object.values(app.views);
  els.forEach((el, i) => {
    const v = views.find((x) => x.el === el);
    v.setPlane(order[i], order[i] === 'sagittal');
    v.point = state.point;
    v.redraw();
  });
}

function sectionHover(w, e) {
  const tip = $('tip');
  if (!w || !app.sections) { tip.hidden = true; return; }
  const mode = app.modes[state.mode];
  if (mode?.sectionAt) {
    const at = mode.sectionAt(w);
    if (!at) { tip.hidden = true; return; }
    const r = $('stage').getBoundingClientRect();
    tip.style.left = `${e.clientX - r.left}px`; tip.style.top = `${e.clientY - r.top}px`;
    tip.innerHTML = `<div class="t"><i style="background:${esc(at.colour)}"></i>${esc(at.name)}</div><div class="s">${esc(at.sub || '')} · MNI ${w.map((x) => Math.round(x)).join(', ')} mm</div>`;
    tip.hidden = false;
    return;
  }
  const L = app.sections.labelAt(w);
  const key = app.labelValues.get(L);
  const at = key && structureOfKey(key);
  if (!at) { tip.hidden = true; return; }
  const s = structure(at.id);
  const r = $('stage').getBoundingClientRect();
  tip.style.left = `${e.clientX - r.left}px`; tip.style.top = `${e.clientY - r.top}px`;
  tip.innerHTML = `<div class="t"><i style="background:${esc(colourOf(s.id) || '#ccc')}"></i>${esc(s.name)}</div><div class="s">${esc(at.side !== 'mid' ? `${at.side} hemisphere · ` : '')}MNI ${w.map((x) => Math.round(x)).join(', ')} mm</div>`;
  tip.hidden = false;
}

function updateSectionSelection() {
  if (!app.sections) return;
  const mode = app.modes[state.mode];
  if (mode?.sectionSelection) {
    app.sections.selected = mode.sectionSelection();
    for (const v of Object.values(app.views)) v.redraw();
    if (state.cut && state.cut.cap) updateCap();
    return;
  }
  const keys = selectionKeys();
  const lit = new Set();
  if (keys) for (const [value, key] of app.labelValues) if (keys.has(key)) lit.add(value);
  app.sections.selected = lit.size ? lit : null;
  app.sections.dimOthers = Boolean(state.isolate || state.system);
  for (const v of Object.values(app.views)) v.redraw();
  if (state.cut && state.cut.cap) updateCap();
}

// ---- the cut ----------------------------------------------------------------------------------

const CUT_RANGES = { sagittal: [-75, 75], coronal: [-110, 75], axial: [-50, 85] };

function setCut(axis) {
  if (!axis || axis === 'off') { state.cut = null; app.viewer.setCut(null); renderToolbar(); schedulePlate(); return; }
  const [min, max] = CUT_RANGES[axis];
  const position = state.cut && state.cut.axis === axis ? state.cut.position : axis === 'sagittal' ? -4 : axis === 'coronal' ? -18 : 12;
  state.cut = { axis, position, keep: state.cut?.keep ?? (axis === 'sagittal' ? 1 : axis === 'coronal' ? -1 : -1), cap: state.cut?.cap ?? true, min, max };
  updateCut();
  renderToolbar();
}

function updateCut() {
  const c = state.cut;
  if (!c) return;
  app.viewer.setCut({ axis: c.axis, position: c.position, keep: c.keep, cap: c.cap && Boolean(app.sections), all: c.cap });
  if (c.cap) updateCap();
  schedulePlate();
}

function updateCap() {
  const c = state.cut;
  if (!c || !app.sections) return;
  if (!app.capCanvas) app.capCanvas = document.createElement('canvas');
  app.sections.draw(app.capCanvas, c.axis, c.position, { scale: 3, boundaries: true });
  app.viewer.setCapImage(app.capCanvas, app.sections.extent(c.axis));
}

// ---- explosion --------------------------------------------------------------------------------

function computeExplodeOffsets() {
  for (const [id, list] of app.instances) {
    const anc = new Set(ancestorsOf(id).map((a) => a.id));
    for (const inst of list) {
      const s = inst.side === 'left' ? -1 : inst.side === 'right' ? 1 : 0;
      const c = inst.centroid;
      let o;
      if (id === 'cerebral-cortex' || id === 'cerebral-white-matter') o = [s * 58, 0, 0];
      else if (anc.has('cerebral-cortex')) o = [s * 58, 0, 0];
      else if (id === 'cerebellar-cortex' || anc.has('cerebellum') || id === 'cerebellar-white-matter') o = [s * 26, -46, -40];
      else if (id === 'brainstem' || anc.has('brainstem')) o = [c.x * 0.6, -14, -52];
      else if (anc.has('fibre-bundles')) o = [s * 34 + c.x * 0.4, c.y * 0.15, c.z * 0.15];
      else o = [s * 18 + (c.x - s * 14) * 1.1, (c.y + 10) * 0.75, (c.z - 4) * 0.7 + 10];
      inst.offset.set(o[0], o[1], o[2]);
    }
  }
}

// ---- the plate of labels --------------------------------------------------------------------------

let plateTimer = null;
function schedulePlate(delay = 120) {
  clearTimeout(plateTimer);
  plateTimer = setTimeout(drawPlate, delay);
}

function drawPlate() {
  const svg = $('plate'), el = $('plate-labels');
  if (!state.labels || app.plateOff || inSections() || !app.instances.size) { clearPlate(svg, el); return; }
  const sel = selectionKeys();
  const candidates = [];
  const cortexShown = state.cortexStyle === 'solid' && !app.hidden.has('cerebral-cortex');
  if (cortexShown && !app.cortexRegion) {
    for (const inst of app.instances.get('cerebral-cortex') || []) {
      if (!inst.visible) continue;
      const index = app.viewer.labelIndex(inst);
      // Divisions of one gyrus are labelled once, by the gyrus.
      const groups = new Map();
      for (const [label, verts] of index) {
        if (!label) continue;
        const at = structureOfKey(`ho:${label}`);
        if (!at) continue;
        const s = structure(at.id);
        const parent = structure(s.parent);
        const whole = parent && parent.kind === 'group' && !['inferior-frontal-gyrus'].includes(parent.id) ? parent : s;
        if (!groups.has(whole.id)) groups.set(whole.id, { id: whole.id, side: at.side, labels: [], verts: [] });
        const g = groups.get(whole.id);
        g.labels.push(label);
        g.verts.push(verts);
      }
      inst.samples = inst.samples || new Map();
      for (const g of groups.values()) {
        if (!inst.samples.has(g.id)) {
          const all = new Uint32Array(g.verts.reduce((n, v) => n + v.length, 0));
          let o = 0;
          for (const v of g.verts) { all.set(v, o); o += v.length; }
          inst.samples.set(g.id, samplesOf(inst, all, 200));
        }
        candidates.push({ id: g.id, side: g.side, name: structure(g.id).short, colour: colourOf(g.id), inst, samples: inst.samples.get(g.id), selected: Boolean(sel && g.labels.some((l) => sel.has(`ho:${l}`))), weight: 1 });
      }
    }
  }
  for (const [id, list] of app.instances) {
    if (id === 'cerebral-cortex') continue;
    for (const inst of list) {
      if (!inst.visible || inst.style !== 'solid') continue;
      inst.allSamples = inst.allSamples || samplesOf(inst, null, 220);
      const selected = Boolean(sel && sel.has(inst.key));
      // The segments of the corpus callosum are named as one, unless one is chosen.
      const s = !selected && structure(id).parent === 'corpus-callosum' ? structure('corpus-callosum') : structure(id);
      candidates.push({ id: s.id, side: inst.side, name: s.short, colour: colourOf(s.id), inst, samples: inst.allSamples, selected, weight: SHELLS.has(id) ? 0.6 : 2.2 });
    }
  }
  layoutPlate({ viewer: app.viewer, candidates, max: 14, svg, labelsEl: el });
}

function screenBox() {
  const v = app.viewer;
  let x0 = Infinity, x1 = -Infinity;
  for (const list of app.instances.values()) for (const inst of list) {
    if (!inst.visible || inst.style === 'ghost') continue;
    inst.allSamples = inst.allSamples || samplesOf(inst, null, 220);
    const s = inst.allSamples;
    for (let i = 0; i < s.length; i += 6) {
      const q = v.project([s[i], s[i + 1], s[i + 2]], inst);
      if (q) { x0 = Math.min(x0, q.x); x1 = Math.max(x1, q.x); }
    }
  }
  if (!Number.isFinite(x0)) return { x0: v.width * 0.3, x1: v.width * 0.7 };
  return { x0: Math.max(0, x0), x1: Math.min(v.width, x1) };
}

// The labels fade while the view moves and return when it rests.
let restTimer = null;
function onView(kind) {
  for (const fn of app.viewHooks) fn(kind);
  $('plate-labels').classList.add('moving');
  $('plate').style.opacity = 0;
  clearTimeout(restTimer);
  restTimer = setTimeout(() => { $('plate-labels').classList.remove('moving'); $('plate').style.opacity = 1; drawPlate(); save(); }, 260);
}

// ---- the card -----------------------------------------------------------------------------------

function sideText(sel) {
  if (sel.side === 'left') return 'Left hemisphere';
  if (sel.side === 'right') return 'Right hemisphere';
  const s = structure(sel.id);
  const k = s.mesh || s.label;
  if (k && k.mid) return 'On the midline';
  return 'Both hemispheres';
}

// The volume and centre of a selection, side by side: {volume, centroid, sides: {left, right, mid}}.
// The centroid is that of the chosen side, or of the left when both are chosen, so that a
// section through it cuts the structure rather than the midline between its halves.
function measures(sel) {
  const keys = keysOf(sel.id, 'label');
  const sides = sel.side ? [sel.side, 'mid'] : ['left', 'right', 'mid'];
  const by = {};
  for (const side of sides) {
    let vol = 0, n = 0;
    const c = [0, 0, 0];
    for (const k of keys[side]) {
      const r = app.product.regions[k];
      if (!r) continue;
      const v = r.joined ? r.joined : r;
      vol += v.volume; n += v.voxels;
      for (let a = 0; a < 3; a++) c[a] += v.centroid[a] * v.voxels;
    }
    if (n) by[side] = { volume: vol, voxels: n, centroid: c.map((x) => x / n) };
  }
  const parts = Object.values(by);
  if (!parts.length) return null;
  const volume = parts.reduce((t, p) => t + p.volume, 0);
  const main = by[sel.side] || (by.left && by.right ? by.left : null) || by.mid || parts[0];
  return { volume, centroid: main.centroid, sides: by };
}

const signed = (x) => (x < 0 ? `−${-x}` : `${x}`);
const coords = (c) => c.map((x) => signed(Math.round(x))).join(', ');

function contactsOf(sel) {
  const keys = keysOf(sel.id, 'label');
  const own = new Set([...keys.left, ...keys.right, ...keys.mid]);
  const sides = sel.side ? [sel.side, 'mid'] : ['left', 'right', 'mid'];
  const tally = new Map();
  for (const side of sides) for (const k of keys[side]) {
    const r = app.product.regions[k];
    if (!r || !r.contacts) continue;
    for (const [other, faces] of r.contacts) {
      if (own.has(other)) continue;
      const at = structureOfKey(other);
      if (!at || ['cerebral-white-matter', 'csf', 'vessel', 'white-matter-hypointensities'].includes(at.id)) continue;
      const anc = ancestorsOf(sel.id).map((a) => a.id);
      if (anc.includes(at.id)) continue;
      const t = tally.get(at.id) || 0;
      tally.set(at.id, t + faces);
    }
  }
  return [...tally].sort((a, b) => b[1] - a[1]).slice(0, 12).map(([id]) => ({ id, name: structure(id).name, colour: colourOf(id) }));
}

function nearOf(sel) {
  const keys = keysOf(sel.id, 'mesh');
  const sides = sel.side ? [sel.side] : ['left', 'right'];
  const tally = new Map();
  for (const side of sides) for (const k of keys[side]) {
    const r = app.product.regions[k];
    if (!r || !r.near) continue;
    for (const [other, count] of r.near) {
      const at = structureOfKey(other);
      if (at) tally.set(at.id, (tally.get(at.id) || 0) + count);
    }
  }
  return [...tally].sort((a, b) => b[1] - a[1]).slice(0, 10).map(([id]) => ({ id, name: structure(id).name, colour: colourOf(id) }));
}

export function renderCard() {
  const el = $('card');
  const sel = state.selected;
  if (app.cardOverride) { el.innerHTML = app.cardOverride(); bindCard(); return; }
  if (!sel) { el.innerHTML = cardHtml({ kind: 'welcome', count: BRAIN_STRUCTURES.filter((s) => s.mesh || s.label).length }); return; }
  if (sel.system) {
    const sys = system(sel.system);
    el.innerHTML = cardHtml({ kind: 'system', name: sys.name, summary: sys.summary, members: sys.members.map((id) => ({ id, name: structure(id).name, colour: colourOf(id) })) });
    bindCard();
    return;
  }
  const s = structure(sel.id);
  const { entry, inherited } = knowledgeOf(s.id);
  const inh = {};
  for (const [f, v] of Object.entries(inherited)) inh[f] = { name: structure(v.id)?.name || v.id, value: v.value };
  const m = measures(sel);
  const facts = [];
  if (m) {
    const { left, right, mid } = m.sides;
    if (left && right) {
      facts.push(['Volume, left', volumeText(left.volume)]);
      facts.push(['Volume, right', volumeText(right.volume)]);
      // The asymmetry index 2(L − R)/(L + R), positive when the left is larger.
      const ai = (2 * (left.volume - right.volume)) / (left.volume + right.volume);
      facts.push(['Asymmetry', Math.abs(ai) < 0.02 ? 'symmetric <small>within 2 %</small>' : `${ai > 0 ? 'left' : 'right'} larger by ${Math.round(Math.abs(ai) * 100)} <small>%</small>`]);
      facts.push(['Centre, MNI', `<span class="pair">L ${coords(left.centroid)}</span><span class="pair">R ${coords(right.centroid)}</span>`]);
    } else {
      facts.push(['Volume', volumeText(m.volume)]);
      facts.push(['Centre, MNI', `${coords((left || right || mid).centroid)} <small>mm</small>`]);
      facts.push(['Talairach', coords(mniToTalairach((left || right || mid).centroid))]);
    }
  }
  if (app.cardFacts) facts.push(...app.cardFacts(sel));
  const path = ancestorsOf(s.id).filter((a) => a.id !== 'brain').reverse();
  const kids = childrenOf(s.id).map((c) => ({ id: c.id, name: c.name, colour: colourOf(c.id) }));
  const systems = SYSTEMS.filter((x) => x.members.includes(s.id) || x.members.some((mm) => ancestorsOf(mm).some((a) => a.id === s.id)));
  const sources = sourceOf(s);
  el.innerHTML = cardHtml({
    kind: 'structure', name: s.name, colour: colourOf(s.id), knowledge: entry, inherited: inh, path, sideText: sideText(sel),
    facts, contacts: contactsOf(sel), near: nearOf(sel), children: kids, systems, isolated: state.isolate,
    hasSections: Boolean(app.sections && m), canHide: drawable(s.id) && treeVisible(s.id) !== null,
    fallback: entry ? null : fallbackText(s), sources: [...sources, ...referencesOf(entry)],
  });
  bindCard();
  drawMiniSections(m);
}

function sourceOf(s) {
  const k = (s.mesh || s.label);
  const key = (k && (k.left || k.mid || k.right)) || '';
  const src = key.split(':')[0];
  return {
    ho: ['Delineation: Harvard–Oxford cortical atlas (FSL), at 25 % probability, on the MNI 152 2009c template.'],
    aseg: ['Delineation: FreeSurfer’s segmentation of the MNI 152 2009c template (TemplateFlow).'],
    massp: ['Delineation: MASSP 2.0, 7 T quantitative MRI of the AHEAD database (Bazin et al. 2020, 2025).'],
    mial: ['Delineation: thalamic nuclei from the diffusion MRI of 67 adults (Najdenovska et al. 2018).'],
    dwm: ['Delineation: atlas of deep white matter bundles, tractography of 78 adults (Guevara et al. 2012), EBRAINS.'],
    jlh: ['Delineation: Julich-Brain 3.1 cytoarchitectonic maps (Amunts et al. 2020), EBRAINS.'],
    cortex: ['Surface: from the tissue probability maps of the MNI 152 2009c template.'],
    white: ['Surface: from the tissue probability maps of the MNI 152 2009c template.'],
    cerebellum: ['Surface: from the tissue probability maps of the MNI 152 2009c template.'],
  }[src] || [];
}

function fallbackText(s) {
  const path = ancestorsOf(s.id).filter((a) => a.id !== 'brain').map((a) => a.name.toLowerCase());
  return `Part of the ${path[0] || 'brain'}${path.length > 1 ? `, in the ${path.slice(1).join(', in the ')}` : ''}.`;
}

function drawMiniSections(m) {
  if (!app.sections || !m) return;
  for (const c of document.querySelectorAll('#card canvas[data-mini]')) {
    const plane = c.dataset.mini;
    const ax = PLANES[plane].axis;
    app.sections.draw(c, plane, m.centroid[ax], { scale: 2, flipU: plane === 'sagittal' });
    c.addEventListener('click', () => { setPoint(m.centroid); state.mainPlane = plane; setMode('sections'); });
  }
}

function bindCard() {
  const el = $('card');
  el.onclick = (e) => {
    const b = e.target.closest('button');
    if (!b) return;
    if (b.dataset.select) { select({ id: b.dataset.select, side: b.dataset.side || state.selected?.side || null }); return; }
    if (b.dataset.system) { select({ system: b.dataset.system }); return; }
    const a = b.dataset.action;
    if (a === 'deselect') select(null);
    if (a === 'focus') focus();
    if (a === 'isolate') { state.isolate = !state.isolate; applyStyles(); paintCortex(); renderCard(); updateSectionSelection(); schedulePlate(); }
    if (a === 'sections') { const m = measures(state.selected); if (m) setPoint(m.centroid); setMode('sections'); }
    if (a === 'hide') { toggleStructure(state.selected.id); }
    app.modes[state.mode]?.cardAction?.(a, b);
  };
}

// ---- the left panel ------------------------------------------------------------------------------

function renderLeft() {
  const left = $('left');
  const mode = app.modes[state.mode];
  if (mode && mode.panel) { left.innerHTML = mode.panel(); mode.bind && mode.bind(left); return; }
  const seg = (attr, value, label, on) => `<button type="button" class="opt" data-${attr}="${value}" aria-pressed="${on}">${label}</button>`;
  left.innerHTML = `<h3>The cortex, coloured by</h3><div class="pad seg">${seg('colouring', 'gyri', 'Gyri', state.colouring === 'gyri')}${seg('colouring', 'lobes', 'Lobes', state.colouring === 'lobes')}${seg('colouring', 'natural', 'Specimen', state.colouring === 'natural')}</div>
<h3>Layers</h3><div class="pad seg">${seg('layer', 'white', 'White matter', state.showWhite)}${seg('layer', 'bundles', 'Fibre bundles', state.showBundles)}${seg('layer', 'fine', 'Finer nuclei', Boolean(state.fine))}</div>
<h3>Systems</h3><div class="systems" id="systems"></div><h3 style="margin-top:20px">Structures</h3><div id="tree-host"></div>`;
  left.querySelectorAll('.seg').forEach((el) => el.addEventListener('click', (e) => onLeftOption(e)));
  $('tree-host').appendChild(app.tree.el);
  app.tree.el.className = 'tree';
  app.tree.render();
  renderSystems();
}

async function onLeftOption(e) {
  const b = e.target.closest('button');
  if (!b) return;
  if (b.dataset.colouring) { state.colouring = b.dataset.colouring; for (const inst of app.instances.get('cerebral-cortex') || []) inst.borderFlags = null; paintCortex(); if (app.sections) { app.sections.palette = labelPalette(); updateSectionSelection(); } }
  if (b.dataset.layer === 'white') {
    await ensurePack('white');
    state.showWhite = !state.showWhite;
    if (state.showWhite) { app.hidden.delete('cerebral-white-matter'); if (state.cortexStyle === 'solid') state.cortexStyle = 'glass'; } else app.hidden.add('cerebral-white-matter');
  }
  if (b.dataset.layer === 'bundles') {
    await ensurePack('bundles');
    state.showBundles = !state.showBundles;
    for (const d of descendantsOf('fibre-bundles')) if (app.instances.has(d.id)) { if (state.showBundles) app.hidden.delete(d.id); else app.hidden.add(d.id); }
    if (state.showBundles && state.cortexStyle === 'solid') state.cortexStyle = 'glass';
  }
  if (b.dataset.layer === 'fine') {
    state.fine = !state.fine;
    for (const id of HIDDEN_AT_FIRST) if (id !== 'cerebral-white-matter' && id !== 'internal-capsule') { if (state.fine) app.hidden.delete(id); else app.hidden.add(id); }
  }
  applyStyles(); renderToolbar(); renderLeft(); schedulePlate(); save();
}

function renderSystems() {
  const el = $('systems');
  if (!el) return;
  el.innerHTML = SYSTEMS.map((s) => `<button type="button" class="chip" data-system="${esc(s.id)}" aria-pressed="${state.system === s.id}">${esc(s.name.replace(/:.*$/, ''))}</button>`).join('');
  el.onclick = (e) => {
    const b = e.target.closest('[data-system]');
    if (!b) return;
    if (state.system === b.dataset.system) select(null);
    else {
      select({ system: b.dataset.system });
      if (state.cortexStyle === 'solid') { state.cortexStyle = 'glass'; applyStyles(); renderToolbar(); }
      // The parts of a system that are hidden at first are shown.
      for (const id of system(b.dataset.system).members) if (app.instances.has(id)) app.hidden.delete(id);
      if (system(b.dataset.system).members.some((id) => descendantsOf('fibre-bundles').some((d) => d.id === id))) ensurePack('bundles').then(() => { state.showBundles = true; for (const id of system(b.dataset.system).members) app.hidden.delete(id); applyStyles(); renderToolbar(); });
      applyStyles();
      app.tree.render();
    }
  };
}

// ---- the toolbar ----------------------------------------------------------------------------------

function renderToolbar() {
  const el = $('toolbar');
  const mode = app.modes[state.mode];
  el.innerHTML = mode?.toolbar ? mode.toolbar() : toolbarHtml(state);
}

function bindChrome() {
  const tb = $('toolbar');
  tb.addEventListener('click', async (e) => {
    const b = e.target.closest('button');
    if (!b) return;
    if (b.dataset.modeAction) { app.modes[state.mode]?.toolbarAction?.(b.dataset.modeAction, b); return; }
    if (b.dataset.show) { state.show = b.dataset.show; showStage(); }
    if (b.dataset.cortex) { state.cortexStyle = b.dataset.cortex; app.hidden.delete('cerebral-cortex'); applyStyles(); paintCortex(); app.tree?.render(); }
    if (b.dataset.colouring) { state.colouring = b.dataset.colouring; paintCortex(); if (app.sections) { app.sections.palette = labelPalette(); updateSectionSelection(); } }
    if (b.dataset.cut) setCut(b.dataset.cut);
    if (b.dataset.flip) { state.cut.keep *= -1; updateCut(); }
    if (b.dataset.cap) { state.cut.cap = !state.cut.cap; updateCut(); }
    if (b.dataset.main) { state.mainPlane = b.dataset.main; arrangeSections(); }
    if (b.dataset.toggle === 'labels') { state.labels = !state.labels; schedulePlate(0); }
    if (b.dataset.toggle === 'white') {
      await ensurePack('white');
      state.showWhite = !state.showWhite;
      if (state.showWhite) { app.hidden.delete('cerebral-white-matter'); if (state.cortexStyle === 'solid') state.cortexStyle = 'glass'; } else app.hidden.add('cerebral-white-matter');
      applyStyles(); app.tree.render();
    }
    if (b.dataset.toggle === 'bundles') {
      await ensurePack('bundles');
      state.showBundles = !state.showBundles;
      for (const d of descendantsOf('fibre-bundles')) if (app.instances.has(d.id)) { if (state.showBundles) app.hidden.delete(d.id); else app.hidden.add(d.id); }
      if (state.showBundles && state.cortexStyle === 'solid') state.cortexStyle = 'glass';
      applyStyles(); app.tree.render();
    }
    renderToolbar();
    schedulePlate();
    save();
  });
  tb.addEventListener('input', (e) => {
    const r = e.target.closest('[data-range]');
    if (!r) return;
    const v = Number(r.value);
    if (r.dataset.range === 'explode') { state.explode = v; app.viewer.setExplode(v); }
    if (r.dataset.range === 'cut') { state.cut.position = v; r.nextElementSibling.textContent = `${v} mm`; updateCut(); }
    if (r.dataset.range === 'overlay') { state.overlay = v; app.sections.overlay = v; r.nextElementSibling.textContent = `${Math.round(v * 100)}%`; for (const s of Object.values(app.views)) s.redraw(); }
  });
  document.querySelectorAll('.modes [data-mode]').forEach((b) => b.addEventListener('click', () => setMode(b.dataset.mode)));
  $('compass').addEventListener('click', (e) => {
    const t = e.target.closest('[data-view]');
    if (!t) return;
    const l = t.dataset.view;
    app.viewer.view(VIEW_DIRECTIONS[l], VIEW_UPS[l] || [0, 0, 1]);
  });
  $('plate-labels').addEventListener('click', (e) => {
    const l = e.target.closest('[data-select]');
    if (l) select({ id: l.dataset.select, side: l.dataset.side && l.dataset.side !== 'mid' ? l.dataset.side : null });
  });
  // The search.
  const find = $('find'), results = $('find-results');
  const entries = () => BRAIN_STRUCTURES.filter((s) => s.id !== 'brain' && s.id !== 'other').map((s) => ({ id: s.id, name: s.name, colour: colourOf(s.id), also: knowledgeOf(s.id).entry?.also || [], path: ancestorsOf(s.id).filter((a) => a.id !== 'brain')[0]?.name || '' }));
  find.addEventListener('input', () => {
    const list = searchResults(find.value, entries());
    results.innerHTML = resultsHtml(list);
    results.hidden = !find.value.trim();
  });
  find.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') { const b = results.querySelector('button'); if (b) b.click(); }
    if (e.key === 'Escape') { find.value = ''; results.hidden = true; find.blur(); }
  });
  find.addEventListener('blur', () => setTimeout(() => { results.hidden = true; }, 180));
  results.addEventListener('mousedown', (e) => {
    const b = e.target.closest('[data-select]');
    if (!b) return;
    e.preventDefault();
    find.value = ''; results.hidden = true;
    revealAndSelect(b.dataset.select);
  });
  $('about-link').addEventListener('click', () => { $('about-body').innerHTML = aboutHtml(); $('about').showModal(); });
  $('tours-link').addEventListener('click', () => { $('about-body').innerHTML = toursHtml(); $('about').showModal(); });
  $('about-body').addEventListener('click', (e) => {
    const b = e.target.closest('[data-tour]');
    if (b) { $('about').close(); startTour(b.dataset.tour); }
  });
  window.addEventListener('keydown', (e) => {
    if (e.target.matches('input, textarea, select')) return;
    if (e.key === 'Escape') { if (app.tour) endTour(); else select(null); }
    if (e.key === 'f' && state.selected) focus();
  });
}

// Choose a structure from the search: show it if it is hidden, and fly to it.
function revealAndSelect(id) {
  if (app.instances.has(id) && app.hidden.has(id)) app.hidden.delete(id);
  for (const d of descendantsOf(id)) if (app.instances.has(d.id) && !HIDDEN_AT_FIRST.has(d.id)) app.hidden.delete(d.id);
  const inBundles = id === 'fibre-bundles' || descendantsOf('fibre-bundles').some((d) => d.id === id);
  const go = () => {
    const deep = !isGyrus(id) && !ancestorsOf(id).some((a) => a.id === 'cerebral-cortex') && id !== 'cerebral-cortex';
    if (deep && state.cortexStyle === 'solid' && !SHELLS.has(id)) state.cortexStyle = 'glass';
    applyStyles(); renderToolbar();
    select({ id, side: null });
    focus();
  };
  if (inBundles) ensurePack('bundles').then(() => { state.showBundles = true; app.hidden.delete(id); for (const d of descendantsOf(id)) app.hidden.delete(d.id); go(); });
  else if (id === 'cerebral-white-matter') ensurePack('white').then(() => { state.showWhite = true; app.hidden.delete(id); go(); });
  else go();
}
export { revealAndSelect };

// ---- modes --------------------------------------------------------------------------------------

// The sections are shown in their own mode, and in the Atlases and Data modes when chosen there.
function inSections() { return state.mode === 'sections' || (['atlases', 'data'].includes(state.mode) && state.show === 'sections'); }
function showStage() {
  const sections = inSections();
  $('sections').hidden = !sections;
  $('viewport').style.visibility = sections ? 'hidden' : 'visible';
  if (sections && app.sections) { app.sections.palette = labelPalette(); arrangeSections(); }
  schedulePlate();
}

app.modes = {};
function setMode(mode) {
  if (!app.product) { app.pendingMode = mode; return; }
  const was = app.modes[state.mode];
  if (was && was.leave && state.mode !== mode) was.leave();
  state.mode = mode;
  document.querySelectorAll('.modes [data-mode]').forEach((b) => b.setAttribute('aria-selected', String(b.dataset.mode === mode)));
  const sections = inSections();
  $('sections').hidden = !sections;
  $('viewport').style.visibility = sections ? 'hidden' : 'visible';
  const m = app.modes[mode];
  if (m && m.enter) m.enter();
  if (sections && app.sections) { app.sections.palette = labelPalette(); arrangeSections(); }
  if (!m || !m.enter) { if (app.sections && !sections) app.sections.palette = labelPalette(); }
  renderLeft();
  renderToolbar();
  renderCard();
  schedulePlate();
  save();
}

// ---- tours ----------------------------------------------------------------------------------------

function toursHtml() {
  return `<h2>Guided tours</h2><p>Short walks through the brain, each a few steps, with the structures shown, cut and coloured as the text needs.</p><div class="tour-list">${TOURS.map((t) => `<button type="button" data-tour="${esc(t.id)}"><b>${esc(t.title)}</b><span>${esc(t.summary)}</span></button>`).join('')}</div>`;
}

function startTour(id) {
  const tour = TOURS.find((t) => t.id === id);
  if (!tour) return;
  app.tour = { tour, step: 0 };
  showStep();
}

async function showStep() {
  const { tour, step } = app.tour;
  const s = tour.steps[step];
  if (state.mode !== 'anatomy') setMode('anatomy');
  if (s.bundles) { await ensurePack('bundles'); }
  // A later step asked for while the bundles were read wins.
  if (!app.tour || app.tour.tour !== tour || app.tour.step !== step) return;
  // What the step shows: the gross structures, and those it names.
  app.hidden = new Set(HIDDEN_AT_FIRST);
  for (const d of descendantsOf('fibre-bundles')) if (app.instances.has(d.id)) app.hidden.add(d.id);
  for (const id of s.show || []) { app.hidden.delete(id); for (const d of descendantsOf(id)) if (app.instances.has(d.id)) app.hidden.delete(d.id); }
  for (const id of s.hide || []) { app.hidden.add(id); for (const d of descendantsOf(id)) app.hidden.add(d.id); }
  state.showBundles = Boolean(s.bundles);
  state.cortexStyle = s.cortex || 'solid';
  state.colouring = s.colouring || 'gyri';
  state.isolate = Boolean(s.isolate);
  state.cut = null; app.viewer.setCut(null);
  if (s.cut) { setCut(s.cut.axis); state.cut.position = s.cut.position; state.cut.keep = s.cut.keep ?? state.cut.keep; state.cut.cap = s.cut.cap ?? true; updateCut(); }
  state.explode = s.explode || 0; app.viewer.setExplode(state.explode);
  paintCortex();
  if (s.select) select(typeof s.select === 'string' ? { id: s.select, side: s.side || null } : s.select);
  else select(null);
  applyStyles();
  if (s.view) app.viewer.flyTo({ target: s.view.target, direction: s.view.direction, distance: s.view.distance, up: s.view.up || [0, 0, 1], duration: 1100 });
  renderToolbar();
  app.tree.render();
  const el = $('tour');
  el.hidden = false;
  el.innerHTML = `<button type="button" class="close" data-tour-action="end" aria-label="End the tour"><svg><use href="#i-close"/></svg></button>
<div class="over">${esc(tour.title)} · ${step + 1} of ${tour.steps.length}</div><h4>${esc(s.title)}</h4><p>${s.text}</p>
<div class="nav"><button type="button" data-tour-action="back"${step ? '' : ' disabled style="opacity:.3"'}>‹ Back</button><span class="dots">${tour.steps.map((_, i) => `<i class="${i === step ? 'on' : ''}"></i>`).join('')}</span><button type="button" data-tour-action="next">${step === tour.steps.length - 1 ? 'Finish' : 'Next ›'}</button></div>`;
  el.onclick = (e) => {
    const b = e.target.closest('[data-tour-action]');
    if (!b) return;
    const a = b.dataset.tourAction;
    if (a === 'end') endTour();
    if (a === 'back' && app.tour.step > 0) { app.tour.step--; showStep(); }
    if (a === 'next') { if (app.tour.step < tour.steps.length - 1) { app.tour.step++; showStep(); } else endTour(); }
  };
}

function endTour() { app.tour = null; $('tour').hidden = true; }

// ---- the initial view and the address --------------------------------------------------------------

function setInitialView() {
  const v = app.viewer;
  const r = state.restoreView;
  if (r) { v.camera.position.set(...r.position); v.controls.target.set(...r.target); v.camera.lookAt(...r.target); v.dirty = true; return; }
  v.controls.target.set(0, -16, 12);
  const dir = [-0.78, 0.38, 0.5];
  const len = Math.hypot(...dir);
  const d = v.fitDistance(98, 1.5);
  v.camera.position.set(0 + (dir[0] / len) * d, -16 + (dir[1] / len) * d, 12 + (dir[2] / len) * d);
  v.camera.lookAt(0, -16, 12);
  v.dirty = true;
}

function save() {
  if (!app.ready) return;
  const q = new URLSearchParams();
  q.set('mode', state.mode);
  if (state.selected?.id) q.set('s', `${state.selected.id}${state.selected.side ? `@${state.selected.side}` : ''}`);
  if (state.selected?.system) q.set('system', state.selected.system);
  q.set('cortex', state.cortexStyle);
  q.set('colour', state.colouring);
  if (state.cut) q.set('cut', `${state.cut.axis},${state.cut.position},${state.cut.keep},${state.cut.cap ? 1 : 0}`);
  if (state.explode) q.set('explode', state.explode.toFixed(2));
  const c = app.viewer.camera.position, t = app.viewer.controls.target;
  q.set('view', [c.x, c.y, c.z, t.x, t.y, t.z].map((x) => Math.round(x)).join(','));
  q.set('p', state.point.map((x) => Math.round(x)).join(','));
  history.replaceState(null, '', `#${q.toString()}`);
}

function restore() {
  const q = new URLSearchParams(INITIAL_HASH.slice(1));
  if (q.get('cortex')) state.cortexStyle = q.get('cortex');
  if (q.get('colour')) state.colouring = q.get('colour');
  const v = (q.get('view') || '').split(',').map(Number);
  if (v.length === 6 && v.every(Number.isFinite)) state.restoreView = { position: v.slice(0, 3), target: v.slice(3) };
  const p = (q.get('p') || '').split(',').map(Number);
  if (p.length === 3 && p.every(Number.isFinite)) state.point = p;
  const cut = (q.get('cut') || '').split(',');
  if (cut.length >= 2 && CUT_RANGES[cut[0]]) {
    const [min, max] = CUT_RANGES[cut[0]];
    state.cut = { axis: cut[0], position: Number(cut[1]), keep: Number(cut[2]) || 1, cap: cut[3] !== '0', min, max };
    setTimeout(updateCut, 0);
  }
  if (q.get('explode')) { state.explode = Number(q.get('explode')) || 0; app.viewer.setExplode(state.explode); }
  const s = q.get('s');
  if (s) { const [id, side] = s.split('@'); if (structure(id)) state.selected = { id, side: side || null }; }
  if (q.get('system') && system(q.get('system'))) state.selected = { system: q.get('system') }, state.system = q.get('system');
  const mode = q.get('mode');
  if (mode && mode !== 'anatomy') app.pendingMode = mode;
  setTimeout(() => {
    if (state.selected) select(state.selected);
    if (app.pendingMode) { const m = app.pendingMode; app.pendingMode = null; setMode(m); } else setMode('anatomy');
  }, 0);
}

// ---- about --------------------------------------------------------------------------------------

function aboutHtml() {
  return `<h2>About the Brain Atlas</h2>
<p>One brain, the unbiased average of 152 adults (MNI 152 ICBM 2009c nonlinear asymmetric, Fonov et al. 2011), with its parts delineated by several atlases on the same 1 mm grid and drawn as smooth surfaces. Everything is computed on this computer from the data library’s files; nothing is sent anywhere.</p>
<h4>The parts</h4>
<ul>
<li><b>Cortex and white matter:</b> surfaces from the template’s tissue probability maps, the cortex coloured by the gyri of the Harvard–Oxford atlas (lateralised, 25 %).</li>
<li><b>Deep structures, ventricles, brainstem, cerebellum:</b> FreeSurfer’s segmentation of the template (aseg).</li>
<li><b>Finer nuclei and tracts:</b> MASSP 2.0 (Bazin et al. 2020, 2025), from 7 T quantitative MRI; the thalamic nuclei of Najdenovska et al. (2018); the cerebellar nuclei of Julich-Brain 3.1 (Amunts et al. 2020).</li>
<li><b>Fibre bundles:</b> the atlas of deep white matter bundles (Guevara et al. 2012), EBRAINS.</li>
</ul>
<h4>How it is drawn</h4>
<p>Surfaces are extracted by surface nets from blurred and upsampled masks, relaxed by Taubin’s smoothing, and shaded by a baked volumetric obscurance, so that sulci and creases fall into shadow (data/products/brain_anatomy.py). The sections show the template’s T1 image with every structure washed in its colour.</p>
<h4>What is said of each part</h4>
<p>The descriptions follow the standard references: Kandel et al., <i>Principles of Neural Science</i> (2021); Standring, <i>Gray’s Anatomy</i> (2020); Haines and Mihailoff, <i>Fundamental Neuroscience</i> (2018); Purves et al., <i>Neuroscience</i> (2018); Vanderah and Gould, <i>Nolte’s The Human Brain</i> (2021); Duvernoy, <i>The Human Brain</i> (1999); Tatu et al. on arterial territories (1996, 1998); Catani and Thiebaut de Schotten on tracts (2008). Latin names follow <i>Terminologia Neuroanatomica</i> (FIPAT 2017). They are for learning, not for clinical decisions.</p>
<h4>Licences</h4>
<p>The software is MIT-licensed. The data keep their own licences: the MNI templates (permissive), MASSP 2.0 and the thalamic nuclei (CC BY 4.0), the Harvard–Oxford atlas (FSL licence, non-commercial), Julich-Brain and the white matter bundles (CC BY-NC-SA 4.0). They are fetched into the cache and not redistributed.</p>`;
}

// ---- the modes that colour and map --------------------------------------------------------------------

const ctx = {
  app, state, structure, colourOf, select, renderCard, paintCortex, paintAll, applyStyles, schedulePlate, toast, setPoint, renderLeft, ensurePack, labelPalette, updateSectionSelection,
  setMode, renderToolbar, focusBox, onViewHook: (fn) => { app.viewHooks.push(fn); },
  sectionsReady: () => Boolean(app.sections),
};
app.modes.atlases = new AtlasesMode(ctx);
app.modes.data = new DataMode(ctx);
app.modes.head = new HeadMode(ctx);

start().catch((e) => { console.error(e); fail(`The atlas could not start: ${esc(e.message)}`); });
