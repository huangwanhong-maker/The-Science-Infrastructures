// The Spiking Simulator: compose a network of spiking neurons, run it in the neuroscience
// domain's C++ engine, and read what it did.
//
// The page holds one specification (spec.js), edited on the left and shown as JSON; a run sends
// it to the server, which runs the engine and returns the results. The engine is built by the
// application's recipe, from the page when it is missing.

import * as api from './api.js';
import { Catalogue, renderCatalogue } from './catalogue.js';
import { el, clear } from './dom.js';
import { renderEditor } from './editor.js';
import { PRESETS, preset } from './presets.js';
import { renderResults } from './results.js';
import * as views from './views.js';
import { clone, empty, forEngine, problems, addPopulation, neurons, normalise } from './spec.js';
import { seconds } from './units.js';

const $ = (id) => document.getElementById(id);
const STORE = 'spiking-simulator:spec';

const state = { spec: null, ran: null, catalogue: null, results: null, running: false, selected: 'lif_psc_exp', presetId: null, figures: [] };

function save() {
  try { localStorage.setItem(STORE, JSON.stringify({ spec: state.spec, preset: state.presetId })); } catch { /* storage unavailable */ }
}
function restore() {
  try { return JSON.parse(localStorage.getItem(STORE) ?? 'null'); } catch { return null; }
}

// ---- the engine and its catalogue ----
async function start() {
  let engine;
  try {
    engine = await api.engineState();
  } catch (error) {
    return fail(`The server does not answer (${error.message}). Start it with: python applications/spiking-simulator/serve.py`);
  }
  $('engine').textContent = engine.current ? `engine ready${engine.toolchain === 'wsl' ? ' · WSL' : ''}` : engine.built ? 'engine out of date' : 'no engine yet';
  $('engine').className = 'engine ' + (engine.current ? 'ok' : 'stale');
  state.limits = engine.limits;
  const threads = clear($('threads'));
  threads.append(el('option', { value: '', text: `all ${engine.limits.threads} threads` }));
  for (let n = 1; n < engine.limits.threads; n *= 2) threads.append(el('option', { value: n, text: `${n} ${n === 1 ? 'thread' : 'threads'}` }));
  if (!engine.current) return offerBuild(engine);
  state.catalogue = new Catalogue(await api.catalogue());
  ready();
}

function offerBuild(engine) {
  const main = clear($('activity'));
  const log = el('pre', { class: 'log' });
  const button = el('button', { class: 'button primary', onclick: () => build(button, log) }, engine.built ? 'Build it again' : 'Build the engine');
  main.append(el('div', { class: 'empty' },
    el('p', { class: 'lead', text: engine.built ? 'The engine was built from older sources.' : 'The engine has not been built on this computer.' }),
    el('p', { class: 'note', text: 'It is the neuroscience domain\'s C++ simulator, compiled once with g++ (natively, or in WSL on Windows). It takes a minute or two.' }),
    button, log));
  showTab('activity');
}

async function build(button, log) {
  button.disabled = true;
  try {
    await api.build();
  } catch (error) {
    log.textContent = error.message;
    button.disabled = false;
    return;
  }
  const poll = async () => {
    const s = await api.setup();
    log.textContent = s.log || '';
    if (s.running) return setTimeout(poll, 1000);
    if (s.state === 'done') return start();
    button.disabled = false;
  };
  poll();
}

function fail(message) {
  clear($('results')).append(el('div', { class: 'empty' }, el('p', { class: 'lead', text: message })));
}

// ---- the page ----
function ready() {
  const params = new URLSearchParams(location.hash.slice(1));
  // A link may open the activity on a span of time (in ms) and on a tab of its own.
  if (params.has('from') && params.has('to')) state.domain = [Number(params.get('from')) / 1e3, Number(params.get('to')) / 1e3];
  state.tab = params.get('tab');
  if (params.get('threads')) $('threads').value = params.get('threads');
  const saved = restore();
  if (params.get('preset') && preset(params.get('preset'))) load(params.get('preset'));
  else if (saved?.spec) { state.spec = normalise(saved.spec); state.presetId = saved.preset; }
  else load('brunel-ai');
  fillPresets();
  draw();
  renderCatalogue($('catalogue'), state.catalogue, {
    selected: state.selected,
    onSelect: (name) => { state.selected = name; },
    onAdd: (name) => { addPopulation(state.spec, name); changed(true); },
  });
  if (params.get('run') === '1') run();
}

function load(id) {
  const p = preset(id);
  state.spec = normalise(clone(p.spec()));
  state.presetId = id;
  state.results = null;
}

function fillPresets() {
  const select = $('preset');
  clear(select).append(el('option', { value: '', text: 'Examples…' }),
    ...PRESETS.map((p) => el('option', { value: p.id, text: p.title, selected: p.id === state.presetId ? true : null })),
    el('option', { value: '__new', text: 'A new network' }));
  select.onchange = () => {
    if (select.value === '__new') { state.spec = empty(); state.presetId = null; addPopulation(state.spec, 'lif_psc_exp'); }
    else if (select.value) load(select.value);
    state.results = null;
    draw();
    save();
  };
}

function draw() {
  renderEditor($('editor'), state.spec, state.catalogue, {
    onChange: () => changed(false),
    onStructure: () => changed(false),
    onFocusModel: (name) => { state.selected = name; },
  });
  drawAbout();
  drawSpec();
  drawProblems();
  drawViews();
}

// The views of the last run, drawn from the specification as it was run; the network's diagram
// from the specification as it is now.
function drawViews() {
  for (const f of state.figures) f.destroy();
  state.figures = [];
  renderResults($('results'), state.results, { spec: state.ran, onDownload: download });
  const ids = ['activity', 'statistics', 'states', 'weights', 'profile'];
  if (!state.results) {
    for (const id of ids) clear($(id)).append(el('div', { class: 'empty' }, el('p', { class: 'lead', text: id === 'activity' ? 'Run the network to see its activity.' : 'Run the network first.' }),
      id === 'activity' ? el('p', { class: 'note', text: 'Choose an example above, or compose a network on the left from the models on the right; Run, or Ctrl+Enter, sends it to the engine.' }) : null));
  } else {
    views.activity($('activity'), state.results, state.ran, state.figures, { domain: state.domain });
    views.statistics($('statistics'), state.results, state.ran, state.figures, state.catalogue);
    views.states($('states'), state.results, state.ran, state.figures);
    views.weights($('weights'), state.results, state.ran, state.figures, state.catalogue);
    views.profile($('profile'), state.results);
  }
  views.network($('network'), state.spec, state.catalogue);
}

function changed(redraw) {
  state.presetId = null;
  save();
  if (redraw) draw();
  else { drawSpec(); drawProblems(); views.network($('network'), state.spec, state.catalogue); }
}

function drawAbout() {
  const p = preset(state.presetId);
  const box = clear($('about'));
  if (p) box.append(el('h2', { text: p.title }), el('p', { text: p.about }));
  else box.append(el('h2', { text: 'Your network' }), el('p', { text: `${neurons(state.spec).toLocaleString('en')} neurons in ${state.spec.populations.length} ${state.spec.populations.length === 1 ? 'population' : 'populations'}.` }));
}

function drawSpec() {
  const text = $('spec-text');
  if (document.activeElement !== text) text.value = JSON.stringify(state.spec, null, 2);
}

function drawProblems() {
  const list = problems(state.spec, state.catalogue);
  const box = clear($('problems'));
  const T = seconds(state.spec.duration);
  if (state.limits) {
    if (neurons(state.spec) > state.limits.neurons) list.push(`this server runs at most ${state.limits.neurons.toLocaleString('en')} neurons`);
    if (T > state.limits.duration) list.push(`this server simulates at most ${state.limits.duration} s in a run`);
  }
  box.hidden = !list.length;
  if (list.length) box.append(el('ul', {}, list.map((p) => el('li', { text: p }))));
  $('run').disabled = state.running || list.length > 0;
}

// ---- running ----
async function run() {
  if (state.running || problems(state.spec, state.catalogue).length) return;
  state.running = true;
  $('run').disabled = true;
  $('stop').hidden = false;
  const t0 = performance.now();
  const status = $('status');
  const tick = setInterval(() => { status.textContent = `running… ${((performance.now() - t0) / 1000).toFixed(1)} s`; }, 100);
  document.body.classList.add('running');
  try {
    const threads = Number($('threads').value) || undefined;
    const ran = clone(state.spec);
    state.results = await api.run(forEngine(ran), threads);
    state.ran = ran;
    status.textContent = `done in ${((performance.now() - t0) / 1000).toFixed(1)} s`;
    drawViews();
    showTab(state.tab ?? 'activity');
    state.tab = null;
  } catch (error) {
    status.textContent = error.status === 499 ? 'stopped' : '';
    if (error.status !== 499) clear($('results')).append(el('div', { class: 'empty' }, el('p', { class: 'lead error', text: `The engine refused: ${error.message}` })));
    showTab('results');
  } finally {
    clearInterval(tick);
    state.running = false;
    $('stop').hidden = true;
    document.body.classList.remove('running');
    drawProblems();
  }
}

function download(what) {
  const data = what === 'results' ? state.results : forEngine(state.spec);
  if (!data) return;
  const name = (state.presetId ?? 'network') + (what === 'results' ? '-results.json' : '.json');
  const url = URL.createObjectURL(new Blob([JSON.stringify(data, null, what === 'results' ? 0 : 2)], { type: 'application/json' }));
  el('a', { href: url, download: name }).click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// ---- tabs, the specification as text ----
function showTab(name) {
  for (const b of document.querySelectorAll('.tabs button')) b.setAttribute('aria-selected', b.dataset.tab === name ? 'true' : 'false');
  for (const p of document.querySelectorAll('.tab')) p.hidden = p.id !== `tab-${name}`;
}

function wire() {
  for (const b of document.querySelectorAll('.tabs button')) b.addEventListener('click', () => showTab(b.dataset.tab));
  $('run').addEventListener('click', run);
  $('stop').addEventListener('click', () => api.cancel());
  document.addEventListener('keydown', (e) => { if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') { e.preventDefault(); run(); } });
  $('spec-apply').addEventListener('click', () => {
    try {
      const spec = JSON.parse($('spec-text').value);
      if (!spec || !Array.isArray(spec.populations)) throw new Error('a specification has a list of populations');
      state.spec = normalise({ ...empty(), ...spec });
      state.presetId = null;
      $('spec-error').textContent = '';
      draw();
      save();
    } catch (error) {
      $('spec-error').textContent = error.message;
    }
  });
  $('spec-download').addEventListener('click', () => download('spec'));
  $('spec-open').addEventListener('change', async (e) => {
    const file = e.target.files?.[0];
    if (!file) return;
    $('spec-text').value = await file.text();
    $('spec-apply').click();
  });
}

wire();
showTab('activity');
start();
