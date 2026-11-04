// The Body Atlas: the whole body in the round, by system, every part measured by the anatomy domain; sections
// in the planes of imaging; tours. This file joins the parts: loading, choosing, the cut, the sections, the
// tours, the search and the address, which keeps the mode, the part chosen, the cut and the view.

import * as THREE from 'three';

import { renderCard, renderProfile } from './card.js';
import { ask, decodeElement, loadCatalogue, loadKnowledge, loadPacks, loadProduct, loadStatus, setCentre } from './data.js';
import { SYSTEM_STYLE } from './palette.js';
import { attachPlateControls, describeSection, renderPlate } from './plate.js';
import { TOURS } from './tours.js';
import { PartsTree, capital, renderFound, renderSystems, search, searchIndex, updateSystems } from './tree.js';
import { Viewer } from './viewer.js';

const $ = (id) => document.getElementById(id);
const AXIS = { sagittal: 0, coronal: 1, transverse: 2 };
// The half each cut keeps (along its LPS axis), so that its face turns to the view it is seen from: a
// transverse cut keeps the lower half, seen from above; a coronal cut the back, seen from the front; a
// sagittal cut the right, seen from the left.
const KEEP = { transverse: -1, coronal: 1, sagittal: -1 };
const FACE = { transverse: 'above', coronal: 'front', sagittal: 'left' };
const VIEWS = {
  front: [0, 0.05, 1], back: [0, 0.05, -1], left: [1, 0.05, 0], right: [-1, 0.05, 0], above: [0, 1, 0.12],
  'front-left': [0.6, 0.1, 0.8], 'front-right': [-0.6, 0.1, 0.8], 'back-left': [0.6, 0.1, -0.8], 'back-right': [-0.6, 0.1, -0.8],
};

const state = {
  mode: 'body',
  styles: {},
  selection: null,      // { id, relation, elements, part }
  cut: { plane: '', level: 0, keep: -1 },
  plate: { plane: 'transverse', level: 1137, convention: 'radiological' },
  plateView: { zoom: 1, du: 0, dv: 0 },
  plateCtx: {},
  section: null,
  tour: null,
};
let catalogue, knowledge, atlas, viewer, tree, index, systemsByKey, extents;

// ---- loading ------------------------------------------------------------------------------------------------

function progress(text, done = false) {
  const li = document.createElement('li');
  li.textContent = text;
  if (done) li.className = 'done';
  $('progress').appendChild(li);
  return li;
}

async function boot() {
  const theme = (() => { try { return localStorage.getItem('body-atlas-theme'); } catch (_) { return null; } })();
  if (theme === 'night') setTheme('night');
  let status;
  try { status = await loadStatus(); } catch (e) { progress(`The server did not answer: ${e.message}`); return; }
  if (!status.product || !status.engine.current) { showSetup(status); return; }
  const p1 = progress('Reading the body and its measures…');
  const [product, cat, kn] = await Promise.all([loadProduct(), loadCatalogue(), loadKnowledge()]);
  catalogue = cat;
  knowledge = kn;
  p1.className = 'done';
  p1.textContent = `${product.counts.elements.toLocaleString('en-GB')} parts, ${(product.counts.triangles / 1e6).toFixed(1)} million triangles`;
  buildAtlas(product);
  setCentre((catalogue.height.top + catalogue.height.bottom) / 2);
  viewer = new Viewer($('viewport'));
  viewer.onCompass = (c) => {
    const k = $('compass');
    k.querySelector('.c-up').textContent = c.up;
    k.querySelector('.c-down').textContent = c.down;
    k.querySelector('.c-left').textContent = c.left;
    k.querySelector('.c-right').textContent = c.right;
  };
  for (const s of catalogue.systems) state.styles[s.key] = (SYSTEM_STYLE[s.key] || {}).show || 'solid';
  for (const [k, v] of Object.entries(state.styles)) viewer.styles[k] = v;
  for (const e of catalogue.elements) if (e.copyOf !== null) viewer.hidden.add(e.i);
  const p2 = progress('Fetching the meshes…');
  const byPack = product.elements.reduce((m, e, i) => { (m[e.pack] ||= []).push(i); return m; }, {});
  await loadPacks(product, (k, bytes) => {
    for (const i of byPack[k]) {
      const e = product.elements[i], c = catalogue.elements[i];
      viewer.addElement(i, decodeElement(e, bytes), c.system, { inward: c.measures.inward, solid: c.measures.topology.consistent });
    }
  }, (got, total) => { p2.textContent = `Fetching the meshes… ${total ? Math.round((100 * got) / total) : 0} %`; });
  p2.className = 'done';
  p2.textContent = 'The meshes';
  wire();
  readAddress();
  $('veil').classList.add('gone');
}

function showSetup(status) {
  const el = $('setup');
  el.hidden = false;
  el.innerHTML = `<p>The atlas needs its data and its engine, prepared once on this computer:</p><ol>`
    + `<li>BodyParts3D's body, about 145 MB${status.product ? ' (ready)' : ''}</li>`
    + `<li>the body as a product, a minute${status.product ? ' (ready)' : ''}</li>`
    + `<li>the engine on the anatomy domain, compiled with g++ (in WSL on Windows), a minute${status.engine.current ? ' (ready)' : ''}</li></ol>`
    + '<button class="key" id="prepare">Prepare</button><p class="hint" id="setup-state"></p>';
  $('prepare').addEventListener('click', async () => {
    await fetch('/setup/run', { method: 'POST' });
    const poll = async () => {
      const s = await (await fetch('/setup')).json();
      const steps = s.steps || [];
      $('setup-state').textContent = steps.map((x) => `${x.title}: ${x.state || (x.done ? 'done' : 'waiting')}`).join(' · ');
      if (steps.length && steps.every((x) => x.done)) location.reload();
      else setTimeout(poll, 2000);
    };
    poll();
  });
}

function buildAtlas(product) {
  const concepts = new Map(catalogue.concepts.map((c) => [c.id, c]));
  const childrenIsa = new Map(), childrenPartof = new Map();
  for (const c of concepts.values()) {
    if (c.inIsa && c.isa && concepts.has(c.isa) && c.elementsIsa.length) (childrenIsa.get(c.isa) || childrenIsa.set(c.isa, []).get(c.isa)).push(c.id);
    if (c.inPartof && c.partof && concepts.has(c.partof) && c.elementsPartof.length) (childrenPartof.get(c.partof) || childrenPartof.set(c.partof, []).get(c.partof)).push(c.id);
  }
  const byName = (a, b) => concepts.get(a).name.localeCompare(concepts.get(b).name);
  for (const m of [childrenIsa, childrenPartof]) for (const l of m.values()) l.sort(byName);
  const rootsIsa = [...concepts.values()].filter((c) => c.inIsa && !c.isa && c.elementsIsa.length).map((c) => c.id);
  const rootsPartof = [...concepts.values()].filter((c) => c.inPartof && !c.partof && c.elementsPartof.length).map((c) => c.id);
  const touching = new Map();
  for (const [a, b, d] of catalogue.touching) {
    (touching.get(a) || touching.set(a, []).get(a)).push({ j: b, d });
    (touching.get(b) || touching.set(b, []).get(b)).push({ j: a, d });
  }
  atlas = { concepts, childrenIsa, childrenPartof, rootsIsa, rootsPartof, touching, product };
  systemsByKey = Object.fromEntries(catalogue.systems.map((s) => [s.key, s]));
  index = searchIndex(atlas);
  extents = [[Infinity, -Infinity], [Infinity, -Infinity], [Infinity, -Infinity]];
  for (const e of catalogue.elements) {
    if (e.copyOf !== null) continue;
    const b = e.measures.bounds;
    for (let k = 0; k < 3; k++) { extents[k][0] = Math.min(extents[k][0], b[k]); extents[k][1] = Math.max(extents[k][1], b[k + 3]); }
  }
}

// ---- the parts chosen -----------------------------------------------------------------------------------------

function elementsOf(id, relation) {
  const c = atlas.concepts.get(id);
  if (!c) return [];
  const list = relation === 'isa' ? c.elementsIsa : c.elementsPartof;
  return list.filter((i) => catalogue.elements[i].copyOf === null);
}

function choose(id, relation = null, { fly = false, quiet = false } = {}) {
  const c = atlas.concepts.get(id);
  if (!c) return;
  const rel = relation || (c.elementsIsa.length ? 'isa' : 'partof');
  const elements = elementsOf(id, rel);
  state.selection = { id, relation: rel, elements, part: null };
  viewer.select(elements);
  tree.reveal(id);
  drawCard();
  if (elements.length > 1) {
    ask('part', { concept: id, relation: rel }).then((part) => {
      if (state.selection?.id !== id) return;
      state.selection.part = part;
      drawCard();
    }).catch(() => {});
  }
  if (fly) viewer.frame(viewer.boxOf(elements));
  if (state.section) { viewer.setCaps(state.section); if (state.mode === 'sections') drawPlate(); }
  if (!quiet) writeAddress();
}

function chooseElement(i, opts) {
  const e = catalogue.elements[i];
  choose(e.concept, null, opts);
}

function clearChoice() {
  state.selection = null;
  viewer.select([]);
  $('card-body').hidden = true;
  $('card-empty').hidden = false;
  writeAddress();
}

function drawCard() {
  const s = state.selection;
  if (!s) return;
  $('card-empty').hidden = true;
  const el = $('card-body');
  el.hidden = false;
  renderCard(el, { atlas, knowledge, catalogue, selection: s, systemsByKey, touching: atlas.touching, isolated: !!viewer.isolated });
  el.querySelectorAll('[data-concept]').forEach((a) => a.addEventListener('click', () => choose(a.dataset.concept, a.dataset.relation || null, { fly: true })));
  el.querySelectorAll('[data-act]').forEach((b) => b.addEventListener('click', () => act(b.dataset.act)));
}

function act(what) {
  const s = state.selection;
  if (!s) return;
  if (what === 'fly') viewer.frame(viewer.boxOf(s.elements));
  if (what === 'isolate') {
    viewer.isolate(viewer.isolated ? null : s.elements);
    if (viewer.isolated) viewer.frame(viewer.boxOf(s.elements));
    drawCard();
    refreshSection();
  }
  if (what === 'hide') {
    for (const i of s.elements) viewer.hidden.add(i);
    viewer.restyleAll();
    clearChoice();
    refreshSection();
  }
  if (what === 'cut') {
    const box = viewer.boxOf(s.elements);
    const m = s.elements.length === 1 ? catalogue.elements[s.elements[0]].measures : s.part?.measures;
    const centre = m ? m.centre : [0, 0, 0];
    if (state.mode === 'sections') setPlate({ level: centre[AXIS[state.plate.plane]] });
    else setCut(state.cut.plane || 'transverse', centre[AXIS[state.cut.plane || 'transverse']]);
    if (!box.isEmpty()) viewer.frame(box);
  }
  if (what === 'calibre') {
    const i = s.elements.reduce((a, b) => (catalogue.elements[b].measures.extents[0] > catalogue.elements[a].measures.extents[0] ? b : a), s.elements[0]);
    const el = document.getElementById('profile');
    el.innerHTML = '<p class="seed">Cutting it along its length…</p>';
    ask('profile', { element: i, axis: 'long', samples: 120 }).then((p) => renderProfile(el, p)).catch((e) => { el.innerHTML = `<p class="seed">${e.message}</p>`; });
  }
}

// ---- systems and visibility -------------------------------------------------------------------------------------

function setStyle(key, style) {
  state.styles[key] = style;
  viewer.setStyle(key, style);
  updateSystems($('systems'), state.styles);
  refreshSection();
  writeAddress();
}

function applySystems(spec) {
  if (!spec) return;
  for (const s of catalogue.systems) {
    const style = spec.only ? spec[s.key] || 'off' : spec[s.key] || state.styles[s.key];
    state.styles[s.key] = style;
    viewer.styles[s.key] = style;
  }
  viewer.restyleAll();
  updateSystems($('systems'), state.styles);
}

// The parts the engine should not cut: those not shown and, in the round, those drawn as glass (whose faces
// would hide the cut). The plate cuts glass parts too: the skin's section is the body's outline.
function hiddenElements() {
  const glassToo = state.mode === 'sections';
  const out = [];
  for (let i = 0; i < catalogue.elements.length; i++) {
    const glass = state.styles[catalogue.elements[i].system] === 'glass';
    const shown = viewer.isolated ? viewer.isolated.has(i) && !viewer.hidden.has(i) : viewer.isVisible(i) && (glassToo || !glass);
    if (!shown) out.push(i);
  }
  return out;
}

// ---- the cut and the sections -----------------------------------------------------------------------------------

let sectionTimer = null, sectionSerial = 0;
function refreshSection() {
  clearTimeout(sectionTimer);
  sectionTimer = setTimeout(requestSection, 90);
}

async function requestSection() {
  const plane = state.mode === 'sections' ? state.plate.plane : state.cut.plane;
  if (!plane) { state.section = null; viewer.setCaps(null); return; }
  const level = state.mode === 'sections' ? state.plate.level : state.cut.level;
  const convention = state.mode === 'sections' ? state.plate.convention : 'radiological';
  const serial = ++sectionSerial;
  try {
    const section = await ask('section', { plane, level, convention, hidden: hiddenElements() });
    if (serial !== sectionSerial) return;
    state.section = section;
    viewer.setCaps(section);
    if (state.mode === 'sections') drawPlate();
    readout();
  } catch (e) {
    $('plate-note').textContent = e.message;
  }
}

function setCut(plane, level = null, keep = null) {
  state.cut.plane = plane;
  if (plane) {
    const [lo, hi] = extents[AXIS[plane]];
    state.cut.level = level ?? (state.cut.level && lo <= state.cut.level && state.cut.level <= hi ? state.cut.level : (lo + hi) / 2);
    if (keep !== null) state.cut.keep = keep;
    const slider = $('cut-level');
    slider.disabled = false;
    slider.min = lo; slider.max = hi; slider.step = 0.5; slider.value = state.cut.level;
  } else {
    $('cut-level').disabled = true;
  }
  $('cut-flip').disabled = !plane;
  for (const b of $('cut-plane').children) b.classList.toggle('on', b.dataset.plane === plane);
  viewer.setCut(plane || null, state.cut.level, state.cut.keep);
  refreshSection();
  readout();
  writeAddress();
}

// A level in the reader's words: a height above the soles (and its vertebra); a depth from the body's middle,
// halfway between its front and its back; a distance from the midline.
function levelWords(plane, level, vertebra) {
  if (plane === 'transverse') return `${Math.round(level - catalogue.height.bottom)} mm above the soles${vertebra ? `, at ${vertebra}` : ''}`;
  if (plane === 'coronal') {
    const d = level - (extents[1][0] + extents[1][1]) / 2;
    return Math.abs(d) < 2 ? 'through the body\'s middle' : `${Math.abs(Math.round(d))} mm ${d < 0 ? 'in front of' : 'behind'} the body\'s middle`;
  }
  return Math.abs(level) < 2 ? 'on the midline' : `${Math.abs(Math.round(level))} mm ${level > 0 ? 'left' : 'right'} of the midline`;
}

function readout() {
  const sec = state.section;
  if (state.cut.plane) $('cut-read').textContent = levelWords(state.cut.plane, state.cut.level, sec && sec.plane === 'transverse' && state.mode === 'body' ? sec.vertebra : '');
  else $('cut-read').textContent = '';
  if (state.mode === 'sections') {
    $('plate-where').textContent = levelWords(state.plate.plane, state.plate.level, sec?.vertebra);
    $('plate-note').textContent = sec ? `${describeSection(sec)} · ${sec.pieces.length} parts cut` : '';
  }
}

function setPlate(change) {
  Object.assign(state.plate, change);
  const [lo, hi] = extents[AXIS[state.plate.plane]];
  if (!(state.plate.level >= lo && state.plate.level <= hi)) state.plate.level = (lo + hi) / 2;
  const slider = $('plate-level');
  slider.min = lo; slider.max = hi; slider.step = 0.5; slider.value = state.plate.level;
  for (const b of $('plate-plane').children) b.classList.toggle('on', b.dataset.plane === state.plate.plane);
  for (const b of $('plate-convention').children) b.classList.toggle('on', b.dataset.convention === state.plate.convention);
  // In the round, keep the half that shows the cut's face to the usual viewpoint.
  viewer.setCut(state.plate.plane, state.plate.level, KEEP[state.plate.plane]);
  refreshSection();
  readout();
  writeAddress();
}

function drawPlate() {
  state.plateCtx = {
    catalogue, atlas, selected: new Set(state.selection?.elements || []), view: state.plateView,
    onChoose: (i) => chooseElement(i),
    onHover: (i, ev) => hoverAt(i, ev),
  };
  renderPlate($('plate-svg'), state.section, state.plateCtx);
}

// The round view in the Sections mode: the trunk about the cut, seen from where its face shows.
function frameCut() {
  const plane = state.plate.plane;
  const level = state.plate.level;
  const centre = [0, 0, 0];
  centre[AXIS[plane]] = level;
  const c = new THREE.Vector3(centre[0], centre[2] - (catalogue.height.top + catalogue.height.bottom) / 2, -centre[1]);
  const half = new THREE.Vector3(260, 300, 200);
  const box = new THREE.Box3(c.clone().sub(half), c.clone().add(half));
  const dir = { transverse: [0, 0.75, 0.66], coronal: [0, 0.2, 1], sagittal: [-1, 0.2, 0.15] }[plane];
  viewer.frame(box, { direction: new THREE.Vector3(...dir), duration: 600 });
}

function vertebraLevel(label) {
  return catalogue.levels.find((l) => l.label === label)?.z;
}

// ---- modes ----------------------------------------------------------------------------------------------------------

function setMode(mode, { quiet = false } = {}) {
  state.mode = mode;
  for (const b of document.querySelectorAll('.modes .key')) b.classList.toggle('on', b.dataset.mode === mode);
  const main = $('main');
  main.classList.toggle('mode-sections', mode === 'sections');
  main.classList.toggle('mode-body', mode === 'body');
  main.classList.toggle('mode-tours', mode === 'tours');
  $('plate').hidden = mode !== 'sections';
  if (mode !== 'tours') closeTour();
  // The plate shows every system but the fascia, until changed on the left; the round view's choice returns after.
  if (mode === 'sections' && !state.stash) {
    state.stash = { ...state.styles };
    const all = { only: true };
    for (const s of catalogue.systems) all[s.key] = s.key === 'fascia' ? 'off' : s.key === 'skin' ? 'glass' : 'solid';
    applySystems(all);
  } else if (mode !== 'sections' && state.stash) {
    applySystems({ only: true, ...state.stash });
    state.stash = null;
  }
  if (mode === 'sections') {
    setPlate({});
    frameCut();
  } else {
    viewer.setCut(state.cut.plane || null, state.cut.level, state.cut.keep);
    refreshSection();
  }
  if (mode === 'tours') openTours();
  setTimeout(() => { viewer.resize(); if (mode === 'sections') drawPlate(); }, 30);
  if (!quiet) writeAddress();
}

function lookFrom(name, what = 'body') {
  const d = new THREE.Vector3(...(VIEWS[name] || VIEWS.front));
  let box;
  if (what === 'select' && state.selection) box = viewer.boxOf(state.selection.elements);
  else if (what === 'isolate' && viewer.isolated) box = viewer.boxOf([...viewer.isolated]);
  else if (what === 'trunk') box = new THREE.Box3(new THREE.Vector3(-220, 0, -180), new THREE.Vector3(220, 700, 180));
  else box = new THREE.Box3(new THREE.Vector3(-350, -870, -200), new THREE.Vector3(350, 870, 200));
  if (box.isEmpty()) return;
  viewer.frame(box, { direction: d });
}

// ---- tours ------------------------------------------------------------------------------------------------------------

function openTours() {
  const box = $('tour');
  box.hidden = false;
  box.innerHTML = `<h4>Tours</h4><div class="tour-choices">${TOURS.map((t) => `<button class="key small" data-tour="${t.key}" title="${t.about}">${t.title}</button>`).join(' ')}</div><p class="text" style="margin-top:10px">${TOURS.map((t) => `<b>${t.title}</b>: ${t.about}`).join('<br>')}</p>`;
  box.querySelectorAll('[data-tour]').forEach((b) => b.addEventListener('click', () => startTour(b.dataset.tour)));
}

function closeTour() {
  state.tour = null;
  $('tour').hidden = true;
}

function startTour(key, step = 0) {
  const tour = TOURS.find((t) => t.key === key);
  if (!tour) return;
  state.tour = { key, step };
  showStep();
}

function showStep() {
  const { key, step } = state.tour;
  const tour = TOURS.find((t) => t.key === key);
  const s = tour.steps[step];
  const box = $('tour');
  box.hidden = false;
  box.innerHTML = `<h4>${tour.title}</h4><p class="title">${s.title}</p><p class="text">${s.text}</p>`
    + `<div class="nav"><span class="count">${step + 1} of ${tour.steps.length}</span>`
    + `<button class="key small" data-go="-1"${step ? '' : ' disabled'}>Back</button>`
    + `<button class="key small" data-go="1">${step + 1 < tour.steps.length ? 'Next' : 'End'}</button>`
    + '<button class="key small" data-go="0">Tours</button></div>';
  box.querySelectorAll('[data-go]').forEach((b) => b.addEventListener('click', () => {
    const g = Number(b.dataset.go);
    if (g === 0) { state.tour = null; openTours(); return; }
    const next = step + g;
    if (next >= tour.steps.length) { state.tour = null; openTours(); return; }
    state.tour.step = next;
    showStep();
  }));
  // The step's view.
  const sectioned = s.mode === 'sections';
  if (sectioned !== (state.mode === 'sections')) {
    state.mode = sectioned ? 'sections' : 'tours';
    $('main').classList.toggle('mode-sections', sectioned);
    $('main').classList.toggle('mode-tours', !sectioned);
    $('plate').hidden = !sectioned;
  }
  applySystems(s.systems || (sectioned ? { only: true, skin: 'solid', skeleton: 'solid', muscles: 'solid', heart: 'solid', arteries: 'solid', veins: 'solid', respiratory: 'solid', digestive: 'solid', urinary: 'solid', nerves: 'solid', brain: 'solid', endocrine: 'solid', lymphoid: 'solid', genital: 'solid' } : null));
  viewer.isolate(null);
  if (s.isolate) {
    const els = s.isolate.flatMap((id) => elementsOf(id, atlas.concepts.get(id)?.elementsIsa.length ? 'isa' : 'partof'));
    viewer.isolate(els);
  }
  if (s.select) choose(s.select, null, { quiet: true });
  else if (!sectioned) clearChoice();
  if (sectioned) {
    const level = s.plate.vertebra ? vertebraLevel(s.plate.vertebra) : s.plate.level;
    Object.assign(state.plateView, { zoom: 1, du: 0, dv: 0 });
    setPlate({ plane: s.plate.plane, level, convention: 'radiological' });
    frameCut();
  } else {
    viewer.setCut(null);
    state.section = null;
    viewer.setCaps(null);
    lookFrom(s.view || 'front', s.frame || 'body');
  }
  setTimeout(() => { viewer.resize(); if (sectioned) drawPlate(); }, 30);
}

// ---- hover and the address --------------------------------------------------------------------------------------

function hoverAt(i, ev) {
  const h = $('hover');
  if (i === null || i === undefined || !ev) { h.hidden = true; return; }
  const e = catalogue.elements[i];
  const c = atlas.concepts.get(e.concept);
  h.innerHTML = `${capital(c?.name || e.name)}<span class="k">${systemsByKey[e.system]?.name || ''}</span>`;
  const stage = $('stage').getBoundingClientRect();
  const inStage = ev.clientX >= stage.left && ev.clientX <= stage.right;
  h.style.left = `${(inStage ? ev.clientX - stage.left : stage.width - 260) + 14}px`;
  h.style.top = `${(inStage ? ev.clientY - stage.top : 12) + 14}px`;
  h.hidden = false;
}

let addressTimer = null, reading = false;
function writeAddress() {
  if (reading) return;
  clearTimeout(addressTimer);
  addressTimer = setTimeout(() => {
    const q = new URLSearchParams();
    q.set('m', state.mode === 'tours' ? 'body' : state.mode);
    if (state.selection) q.set('p', state.selection.id);
    if (state.cut.plane && state.mode === 'body') q.set('cut', `${state.cut.plane}:${Math.round(state.cut.level)}:${state.cut.keep}`);
    if (state.mode === 'sections') q.set('sec', `${state.plate.plane}:${Math.round(state.plate.level)}:${state.plate.convention[0]}`);
    q.set('sys', catalogue.systems.map((s) => `${s.key}.${state.styles[s.key][0]}`).join(','));
    const v = viewer.view();
    q.set('v', [...v.p, ...v.t].join(','));
    history.replaceState(null, '', `#${q.toString()}`);
  }, 250);
}

function readAddress() {
  reading = true;
  const q = new URLSearchParams(location.hash.slice(1));
  if (q.get('sys')) {
    const spec = {};
    for (const part of q.get('sys').split(',')) {
      const [k, s] = part.split('.');
      spec[k] = { s: 'solid', g: 'glass', o: 'off' }[s] || 'solid';
    }
    applySystems(spec);
  }
  if (q.get('v')) {
    const n = q.get('v').split(',').map(Number);
    if (n.length === 6 && n.every(isFinite)) viewer.setView({ p: n.slice(0, 3), t: n.slice(3) });
  } else lookFrom('front', 'body');
  if (q.get('cut')) {
    const [plane, level, keep] = q.get('cut').split(':');
    if (AXIS[plane] !== undefined) setCut(plane, Number(level), Number(keep) || -1);
  }
  if (q.get('sec')) {
    const [plane, level, conv] = q.get('sec').split(':');
    if (AXIS[plane] !== undefined) Object.assign(state.plate, { plane, level: Number(level), convention: conv === 'n' ? 'neurological' : 'radiological' });
  }
  const mode = q.get('m');
  if (mode === 'sections') setMode('sections', { quiet: true });
  if (q.get('p')) choose(q.get('p'), null, { quiet: true, fly: !q.get('v') && state.mode === 'body' });
  reading = false;
}

// ---- wiring ------------------------------------------------------------------------------------------------------

function setTheme(t) {
  document.documentElement.dataset.theme = t === 'night' ? 'night' : '';
  $('theme').innerHTML = `<svg class="i"><use href="#i-${t === 'night' ? 'sun' : 'moon'}"/></svg>`;
  try { localStorage.setItem('body-atlas-theme', t); } catch (_) { /* no storage */ }
}

function wire() {
  renderSystems($('systems'), catalogue.systems, state.styles, {
    onToggle: (k) => setStyle(k, state.styles[k] === 'off' ? ((SYSTEM_STYLE[k] || {}).show === 'glass' ? 'glass' : 'solid') : 'off'),
    onGlass: (k) => setStyle(k, state.styles[k] === 'glass' ? 'solid' : 'glass'),
  });
  tree = new PartsTree($('tree'), atlas, { onChoose: (id, rel) => choose(id, rel, { fly: true }) });
  for (const r of [...atlas.rootsPartof, ...atlas.rootsIsa]) tree.open.add(r);
  tree.render();
  for (const b of $('tree-relation').children) b.addEventListener('click', () => {
    for (const x of $('tree-relation').children) x.classList.toggle('on', x === b);
    tree.setRelation(b.dataset.relation);
  });
  for (const b of document.querySelectorAll('.modes .key')) b.addEventListener('click', () => setMode(b.dataset.mode));
  $('theme').addEventListener('click', () => setTheme(document.documentElement.dataset.theme === 'night' ? 'day' : 'night'));

  // The cut in the round.
  for (const b of $('cut-plane').children) b.addEventListener('click', () => {
    const plane = b.dataset.plane;
    if (!plane) { setCut(''); return; }
    setCut(plane, null, KEEP[plane]);
    lookFrom(FACE[plane], 'body');
  });
  $('cut-level').addEventListener('input', (ev) => { state.cut.level = Number(ev.target.value); viewer.setCut(state.cut.plane, state.cut.level, state.cut.keep); refreshSection(); readout(); writeAddress(); });
  $('cut-flip').addEventListener('click', () => setCut(state.cut.plane, state.cut.level, -state.cut.keep));

  // The plate.
  for (const b of $('plate-plane').children) b.addEventListener('click', () => {
    Object.assign(state.plateView, { zoom: 1, du: 0, dv: 0 });
    setPlate({ plane: b.dataset.plane, level: NaN });
    frameCut();
  });
  for (const b of $('plate-convention').children) b.addEventListener('click', () => setPlate({ convention: b.dataset.convention }));
  $('plate-level').addEventListener('input', (ev) => setPlate({ level: Number(ev.target.value) }));
  attachPlateControls($('plate-svg'), () => state.plateCtx, drawPlate, (dir, big) => setPlate({ level: state.plate.level + dir * (big ? 25 : 5) }));
  new ResizeObserver(() => { if (state.mode === 'sections') drawPlate(); }).observe($('plate'));

  // Picking in the round.
  const canvas = viewer.renderer.domElement;
  let pending = null, down = null;
  canvas.addEventListener('pointermove', (ev) => {
    if (ev.buttons) { $('hover').hidden = true; return; }
    if (pending) return;
    pending = requestAnimationFrame(() => {
      pending = null;
      const hit = viewer.pick(ev.clientX, ev.clientY);
      hoverAt(hit ? hit.i : null, ev);
    });
  });
  canvas.addEventListener('pointerleave', () => { $('hover').hidden = true; });
  canvas.addEventListener('pointerdown', (ev) => { down = [ev.clientX, ev.clientY]; });
  canvas.addEventListener('click', (ev) => {
    if (down && Math.hypot(ev.clientX - down[0], ev.clientY - down[1]) > 4) return;
    const hit = viewer.pick(ev.clientX, ev.clientY);
    if (hit) chooseElement(hit.i); else clearChoice();
  });
  canvas.addEventListener('dblclick', (ev) => {
    const hit = viewer.pick(ev.clientX, ev.clientY);
    if (hit) chooseElement(hit.i, { fly: true });
  });
  viewer.controls.addEventListener('end', () => writeAddress());

  // Finding.
  const input = $('find'), found = $('found');
  let results = [];
  input.addEventListener('input', () => { results = search(index, input.value); renderFound(found, results, systemsByKey, pick); });
  input.addEventListener('keydown', (ev) => {
    if (ev.key === 'Enter' && results.length) pick(results[0].id);
    if (ev.key === 'Escape') { input.value = ''; found.hidden = true; }
  });
  input.addEventListener('blur', () => setTimeout(() => { found.hidden = true; }, 150));
  function pick(id) { found.hidden = true; input.value = ''; input.blur(); choose(id, null, { fly: true }); }
  document.addEventListener('keydown', (ev) => {
    if (ev.key === '/' && document.activeElement !== input) { ev.preventDefault(); input.focus(); }
    if (ev.key === 'Escape' && document.activeElement !== input) { viewer.isolate(null); clearChoice(); }
  });
}

boot().catch((e) => { progress(`Something failed: ${e.message}`); console.error(e); });
