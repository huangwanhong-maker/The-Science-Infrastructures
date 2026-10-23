// The Complex Atlas: cell complexes, their invariants and the fields on them.
//
// The page holds one complex (Complex JSON: README.md, the protocol). It comes from the engine's
// catalogue, from the examples, or from the user (words, facets or JSON, read in input.js). It is
// drawn in three views, all computed here in the page (visualization/js/complex.js): realised in
// space, unfolded into a net, and as a Hasse diagram. Its invariants, a gauge field on it and a
// state sum come from the engine, through the server; the views do not need the engine.

import * as api from './api.js';
import { el, clear } from './dom.js';
import { realise, net, hasse, edgeNames, cyclicColour, divergingColour, cssColour } from '../../../visualization/js/index.js';
import { Realisation } from './views/realisation.js';
import { drawNet } from './views/net.js';
import { drawHasse, HASSE_LIMIT } from './views/hasse.js';
import { renderInvariants } from './invariants.js';
import { renderFields } from './fields.js';
import { renderSums } from './statesums.js';
import { renderCatalogue } from './catalogue.js';
import { readComplex } from './input.js';
import { EXAMPLES } from './examples.js';
import { cellName, wordText, countsText, number } from './format.js';

const $ = (id) => document.getElementById(id);
const STORE = 'complex-atlas:v1';

const state = {
  engine: null, families: null, family: null,
  complex: null, names: [], about: '', source: null,
  geometry: null, net: null, hasse: null,
  analysis: null, analysisState: null,
  fields: null, fieldsShown: true, fieldsRunning: false, fieldsError: null,
  sums: null, sumsRunning: false, sumsError: null,
  selection: null,
  view: 'realisation', panel: 'invariants',
  style: { vertexSize: 1, edgeWidth: 2.5, arrowSize: 1, faceOpacity: 0.38, rings: 3, labels: true, arrows: true, positions: true, seed: 1 },
  input: { text: 'abAB', kind: 'auto' },
};
let realisation = null;
let token = 0;

function save() {
  try {
    localStorage.setItem(STORE, JSON.stringify({ style: state.style, view: state.view, panel: state.panel, input: state.input, source: state.source }));
  } catch { /* storage unavailable */ }
}
function restore() {
  try { return JSON.parse(localStorage.getItem(STORE) ?? 'null'); } catch { return null; }
}
const status = (text) => { $('status').textContent = text || ''; };

// ---- start -------------------------------------------------------------------------------------
async function start() {
  const saved = restore();
  if (saved) {
    Object.assign(state.style, saved.style || {});
    state.view = saved.view || state.view;
    state.panel = saved.panel || state.panel;
    Object.assign(state.input, saved.input || {});
  }
  for (const b of document.querySelectorAll('.views button')) b.onclick = () => { state.view = b.dataset.view; save(); tabs(); draw(); };
  for (const b of document.querySelectorAll('.tabs button')) b.onclick = () => { state.panel = b.dataset.panel; save(); tabs(); };
  tabs();
  renderInput();
  renderDisplay();
  await connect();
  // The first complex: from the address (#example=, #family=, #words=, #facets=), or as left.
  const hash = new URLSearchParams(location.hash.slice(1));
  const source = hash.get('example') ? { example: hash.get('example') }
    : hash.get('family') ? { family: hash.get('family'), params: Object.fromEntries([...hash].filter(([k]) => k !== 'family').map(([k, v]) => [k, Number.isNaN(Number(v)) ? v : Number(v)])) }
    : hash.get('words') ? { text: hash.get('words'), kind: 'words' }
    : hash.get('facets') ? { text: hash.get('facets'), kind: 'facets' }
    : saved?.source || { example: 'torus' };
  await load(source);
}

async function connect() {
  const pill = $('engine');
  try {
    state.engine = await api.engineState();
  } catch (error) {
    state.engine = null;
    pill.textContent = 'no server';
    pill.className = 'engine stale';
    status(`The server does not answer: start it with python applications/complex-atlas/serve.py`);
    renderCatalogue($('catalogue'), { families: null });
    return;
  }
  const e = state.engine;
  pill.textContent = e.fixtures ? 'fixtures' : e.current ? `engine ready${e.toolchain === 'wsl' ? ' · WSL' : ''}` : e.built ? 'engine out of date' : 'no engine yet';
  pill.className = 'engine ' + (e.current ? 'ok' : 'stale');
  pill.title = e.current ? '' : 'Build the engine';
  pill.style.cursor = e.current ? '' : 'pointer';
  pill.onclick = e.current ? null : offerBuild;
  if (e.current) {
    try { state.families = (await api.catalogue()).families; } catch { state.families = null; }
  }
  drawCatalogue();
}

const ready = () => Boolean(state.engine?.current);

function drawCatalogue() {
  renderCatalogue($('catalogue'), {
    families: state.families,
    selected: state.family,
    onSelect: (f) => { state.family = state.family?.id === f.id ? null : { id: f.id, params: {} }; drawCatalogue(); },
    onBuild: (f, params) => load({ family: f.id, params }),
  });
}

// A complex from a source: {example}, {family, params} or {text, kind}.
async function load(source) {
  try {
    if (source.example) {
      const ex = EXAMPLES.find((e) => e.id === source.example) || EXAMPLES[0];
      const { complex } = readComplex(ex.text, ex.kind);
      complex.name = ex.id;
      complex.title = ex.title;
      state.input = { text: ex.text, kind: ex.kind };
      renderInput();
      setComplex(complex, { ...source, about: ex.about });
    } else if (source.family) {
      if (!ready()) throw new Error('the catalogue needs the engine');
      status('Building…');
      const answer = await api.build(source.family, source.params || {});
      const f = state.families?.find((x) => x.id === source.family);
      state.family = { id: source.family, params: source.params || {} };
      drawCatalogue();
      setComplex(answer.complex, { ...source, about: f?.description || '' });
      status('');
    } else {
      const { complex } = readComplex(source.text, source.kind || 'auto');
      setComplex(complex, { ...source, about: 'Your complex, read from what you wrote.' });
    }
  } catch (error) {
    status('');
    if (source.example || source.family) { if (!state.complex) return load({ example: 'torus' }); }
    showInputError(error.message);
  }
}

function setComplex(complex, source) {
  state.complex = complex;
  state.names = edgeNames(complex);
  state.about = source.about || '';
  state.source = { ...source, about: undefined };
  state.geometry = state.net = state.hasse = null;
  state.analysis = null;
  state.fields = state.sums = null;
  state.fieldsError = state.sumsError = null;
  state.selection = null;
  save();
  renderAbout();
  draw(true);
  renderPanels();
  analyse();
}

async function analyse() {
  const mine = ++token;
  if (!ready()) { state.analysisState = null; renderPanels(); return; }
  state.analysisState = 'loading';
  renderPanels();
  try {
    const answer = await api.analyse(state.complex);
    if (mine !== token) return;
    state.analysis = answer.analysis;
    state.analysisState = null;
  } catch (error) {
    if (mine !== token) return;
    state.analysisState = { error: error.message };
  }
  renderPanels();
}

// ---- the views ----------------------------------------------------------------------------------
function tabs() {
  for (const b of document.querySelectorAll('.views button')) b.setAttribute('aria-selected', String(b.dataset.view === state.view));
  for (const b of document.querySelectorAll('.tabs button')) b.setAttribute('aria-selected', String(b.dataset.panel === state.panel));
  // toggleAttribute, since an <svg> has no hidden property.
  for (const v of ['realisation', 'net', 'hasse']) $(`view-${v}`).toggleAttribute('hidden', v !== state.view);
  for (const p of ['invariants', 'fields', 'sums']) $(`panel-${p}`).toggleAttribute('hidden', p !== state.panel);
}

function renderAbout() {
  const c = state.complex;
  const counts = [c.vertices, c.edges.length, (c.faces || []).length, ...(c.cells || []).map((l) => l.length)];
  while (counts.length > 2 && !counts[counts.length - 1]) counts.pop();
  clear($('about')).append(el('h2', { text: c.title || c.name || 'A complex' }), el('p', {}, state.about ? `${state.about} ` : '', el('span', { class: 'note', text: countsText(counts) })));
}

const shownFields = () => (state.fieldsShown && state.fields ? state.fields : null);

function draw(fit = false) {
  const c = state.complex;
  if (!c) return;
  const notice = $('notice');
  notice.hidden = true;
  const onHover = hover;
  if (state.view === 'realisation') {
    if (!state.geometry) {
      const t0 = performance.now();
      const complex = state.style.positions ? c : { ...c, positions: undefined };
      state.geometry = realise(complex, { seed: state.style.seed, rings: state.style.rings });
      const ms = performance.now() - t0;
      if (ms > 400) status(`laid out in ${(ms / 1000).toFixed(1)} s`);
    }
    if (!realisation) realisation = new Realisation($('view-realisation'), { onHover });
    realisation.setNames(state.names);
    realisation.fields = shownFields();
    realisation.selection = state.selection;
    realisation.show(c, state.geometry, state.style, fit);
  } else if (state.view === 'net') {
    if (!(c.faces || []).length) return showNotice('This complex has no 2-cells, so it has no net. Its realisation and Hasse diagram show it.');
    const fresh = !state.net;
    if (fresh) state.net = net(c);
    drawNet($('view-net'), c, state.net, { names: state.names, style: state.style, fields: shownFields(), selection: state.selection, onHover, fit: fit || fresh });
  } else {
    const total = c.vertices + c.edges.length + (c.faces || []).length + (c.cells || []).reduce((s, l) => s + l.length, 0);
    if (total > HASSE_LIMIT) return showNotice(`${total.toLocaleString()} cells are too many to read as a Hasse diagram (at most ${HASSE_LIMIT.toLocaleString()} are drawn).`);
    const fresh = !state.hasse;
    if (fresh) state.hasse = hasse(c);
    drawHasse($('view-hasse'), c, state.hasse, { names: state.names, style: state.style, selection: state.selection, onHover, fit: fit || fresh });
  }
  legend();
}

function showNotice(text, ...more) {
  const n = clear($('notice'));
  n.hidden = false;
  n.append(el('p', { class: 'lead', text }), ...more);
}

function hover(what) {
  const h = $('hover');
  if (!what || !state.complex) { h.hidden = true; return; }
  const c = state.complex, f = shownFields();
  let text = '';
  if (what.dim === 0) text = `vertex v${what.index}`;
  else if (what.dim === 1) {
    const [t, hd] = c.edges[what.index];
    text = `edge ${state.names[what.index]}: v${t} → v${hd}${t === hd ? ', a loop' : ''}`;
    if (f?.link_angle) text += ` · link angle ${number(f.link_angle[what.index], 3)}`;
  } else if (what.dim === 2) {
    const face = c.faces[what.index];
    text = `face ${cellName(2, what.index)}: ∂ = ${wordText(face.boundary, state.names, { powers: false })}`;
    if (f?.plaquettes) text += ` · Re tr U / N = ${number(f.plaquettes[what.index], 3)}`;
  }
  h.textContent = text;
  h.hidden = false;
}

function legend() {
  const box = $('legend');
  const f = shownFields();
  if (!f) { box.hidden = true; return; }
  const ramp = (fn, label, from, to) => {
    const canvas = el('canvas', { width: 120, height: 10 });
    const ctx = canvas.getContext('2d');
    for (let x = 0; x < 120; x++) { ctx.fillStyle = cssColour(fn(x / 119)); ctx.fillRect(x, 0, 1, 10); }
    return el('div', { class: 'ramp' }, el('span', { text: label }), el('span', { text: from }), canvas, el('span', { text: to }));
  };
  clear(box).append(ramp((t) => cyclicColour(t * 2 * Math.PI), 'link angle', '0', '2π'), ramp((t) => divergingColour(2 * t - 1), 'Re tr U / N', '−1', '1'));
  box.hidden = false;
}

function select(selection) {
  state.selection = selection;
  if (realisation && state.view === 'realisation') realisation.highlight(selection);
  else draw();
  renderPanels();
}

// ---- the panels ----------------------------------------------------------------------------------
function renderPanels() {
  const c = state.complex;
  renderInvariants($('panel-invariants'), { complex: c, analysis: state.analysis, names: state.names, selection: state.selection, onSelect: select, state: state.analysisState });
  renderFields($('panel-fields'), {
    complex: c, analysis: state.analysis, names: state.names, result: state.fields, running: state.fieldsRunning, error: state.fieldsError,
    showing: state.fieldsShown, selection: state.selection, onSelect: select,
    onShow: (on) => { state.fieldsShown = on; redrawColours(); },
    onRun: runGauge,
  });
  renderSums($('panel-sums'), { complex: c, result: state.sums, running: state.sumsRunning, error: state.sumsError, onRun: runSums });
  if (!ready() && c) {
    const note = () => el('div', {}, el('p', { class: 'note', text: state.engine ? 'The engine is not built on this computer yet.' : 'The server does not answer.' }),
      state.engine ? el('button', { class: 'button', onclick: offerBuild }, 'Build the engine') : null);
    for (const id of ['panel-invariants', 'panel-fields', 'panel-sums']) $(id).append(note());
  }
}

function redrawColours() {
  if (realisation && state.view === 'realisation') { realisation.setFields(shownFields()); legend(); } else draw();
  renderPanels();
}

async function runGauge(params) {
  if (!ready()) return offerBuild();
  state.fieldsRunning = true;
  state.fieldsError = null;
  renderPanels();
  // Wilson loops along the generators of π1 after simplification (on a lattice torus, the two
  // Polyakov loops), or along the first raw generators when there is no simplified form.
  const pi1 = state.analysis?.fundamental_group;
  const loops = (pi1?.simplified?.loops ?? pi1?.loops ?? []).slice(0, 8);
  try {
    state.fields = await api.gauge({ complex: state.complex, ...params, loops });
    state.fieldsShown = true;
  } catch (error) {
    state.fieldsError = error.message;
  }
  state.fieldsRunning = false;
  redrawColours();
}

async function runSums(r) {
  if (!ready()) return offerBuild();
  state.sumsRunning = true;
  state.sumsError = null;
  renderPanels();
  try { state.sums = await api.turaevViro(state.complex, r); } catch (error) { state.sumsError = error.message; }
  state.sumsRunning = false;
  renderPanels();
}

// ---- your complex --------------------------------------------------------------------------------
function renderInput() {
  const host = clear($('input'));
  const text = el('textarea', { rows: 4, spellcheck: 'false', 'aria-label': 'Your complex: words, facets or JSON' });
  text.value = state.input.text;
  const kind = el('select', { 'aria-label': 'What the text is' },
    ...[['auto', 'read as it looks'], ['words', 'surface words'], ['facets', 'facets'], ['json', 'Complex JSON']].map(([v, t]) => el('option', { value: v, selected: state.input.kind === v ? true : null }, t)));
  const drawIt = () => { state.input = { text: text.value, kind: kind.value }; save(); load({ text: text.value, kind: kind.value }); };
  text.addEventListener('keydown', (e) => { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) drawIt(); });
  const file = el('input', { type: 'file', accept: '.json,application/json,.txt', hidden: true, onchange: async () => {
    const f = file.files[0];
    if (!f) return;
    text.value = await f.text();
    kind.value = 'auto';
    drawIt();
  } });
  const download = () => {
    if (!state.complex) return;
    const blob = new Blob([JSON.stringify(state.complex, null, 1)], { type: 'application/json' });
    const a = el('a', { href: URL.createObjectURL(blob), download: `${state.complex.name || 'complex'}.json` });
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  host.append(
    el('div', { class: 'examples' }, ...EXAMPLES.map((ex) => el('button', { class: 'chip', title: ex.about, onclick: () => load({ example: ex.id }) }, ex.title))),
    text,
    el('div', { class: 'row' }, kind, el('button', { class: 'button primary', title: 'Draw (Ctrl+Enter)', onclick: drawIt }, 'Draw')),
    el('div', { class: 'row' }, el('label', { class: 'button' }, 'Open a file', file), el('button', { class: 'button', onclick: download }, 'Save as JSON')),
    el('p', { class: 'note', text: 'Words: polygons glued along their sides, lowercase forward and uppercase backward (abAB, aa, a1b1A1B1a2b2A2B2). Facets: the simplices of a triangulation (0 1 2; 0 2 3). Or the Complex JSON of the README, which may carry positions for its vertices.' }),
    el('div', { id: 'input-error', class: 'error' }));
}
function showInputError(message) { const e = $('input-error'); if (e) e.textContent = message; }

// ---- display -------------------------------------------------------------------------------------
function renderDisplay() {
  const host = clear($('display'));
  const st = state.style;
  const slider = (label, key, min, max, step, relayout = false) => {
    const out = el('output', { text: String(st[key]) });
    const input = el('input', { type: 'range', min, max, step, value: st[key], oninput: () => {
      st[key] = Number(input.value);
      out.textContent = String(st[key]);
      save();
      if (relayout) { state.geometry = null; draw(); } else restyle();
    } });
    return el('label', { class: 'slider' }, el('span', { text: label }), input, out);
  };
  const check = (label, key, relayout = false) => el('label', { class: 'check' }, el('input', { type: 'checkbox', checked: st[key] ? true : null, onchange: (e) => {
    st[key] = e.target.checked;
    save();
    if (relayout) { state.geometry = null; draw(true); } else restyle();
  } }), label);
  const seed = el('input', { type: 'number', value: st.seed, min: 1, style: 'width: 70px' });
  host.append(
    slider('vertex size', 'vertexSize', 0.2, 3, 0.1),
    slider('edge width', 'edgeWidth', 0.5, 7, 0.5),
    slider('arrow size', 'arrowSize', 0, 2.5, 0.1),
    slider('face opacity', 'faceOpacity', 0, 1, 0.02),
    slider('face smoothing', 'rings', 1, 6, 1, true),
    check('names of the cells', 'labels'),
    check('arrows on the edges', 'arrows'),
    check('the complex\'s own positions, when it has them', 'positions', true),
    el('div', { class: 'row' }, el('span', { class: 'note', text: 'layout seed' }), seed,
      el('button', { class: 'button', onclick: () => { st.seed = Number(seed.value) || 1; save(); state.geometry = null; draw(true); } }, 'Lay out again'),
      el('button', { class: 'button', title: 'Save the realisation as a picture', onclick: snapshot }, 'Picture')));
}

function restyle() {
  if (state.view === 'realisation' && realisation) realisation.setStyle(state.style);
  else draw();
}

function snapshot() {
  if (state.view !== 'realisation' || !realisation) return;
  const a = el('a', { href: realisation.snapshot(), download: `${state.complex?.name || 'complex'}.png` });
  a.click();
}

// ---- building the engine ---------------------------------------------------------------------------
function offerBuild() {
  const log = el('pre', { class: 'log' });
  const button = el('button', { class: 'button primary', onclick: () => build(button, log) }, state.engine?.built ? 'Build it again' : 'Build the engine');
  showNotice(state.engine?.built ? 'The engine was built from older sources.' : 'The engine has not been built on this computer.',
    el('p', { class: 'note', text: 'It is the foundations\' C++ for complexes, compiled once with g++ (natively, or in WSL on Windows); it takes a minute or two. The views need no engine: close this to keep drawing.' }),
    el('div', { class: 'row' }, button, el('button', { class: 'button', onclick: () => { $('notice').hidden = true; } }, 'Close')), log);
}

async function build(button, log) {
  button.disabled = true;
  try {
    await api.buildEngine();
  } catch (error) {
    log.textContent = error.message;
    button.disabled = false;
    return;
  }
  const poll = async () => {
    const s = await api.setup();
    log.textContent = s.log || '';
    if (s.running) return setTimeout(poll, 1000);
    if (s.state === 'done') { $('notice').hidden = true; await connect(); analyse(); return; }
    button.disabled = false;
  };
  poll();
}

start();
