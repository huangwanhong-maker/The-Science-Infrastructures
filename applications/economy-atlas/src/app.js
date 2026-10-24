// The Economy Atlas: four views on the economics domain, computed by the engine through the server.

import * as api from './api.js';
import { el, clear } from './dom.js';
import { showMap } from './map.js';
import { showSpacetime } from './spacetime.js';
import { showProduction } from './production.js';
import { showCircuits } from './circuits.js';

const VIEWS = { map: showMap, spacetime: showSpacetime, production: showProduction, circuits: showCircuits };
const state = { view: null, data: null };

function viewFromHash() {
  const v = location.hash.slice(1);
  return VIEWS[v] ? v : 'map';
}

async function show(view) {
  state.view = view;
  for (const b of document.querySelectorAll('.views button')) b.setAttribute('aria-selected', String(b.dataset.view === view));
  for (const s of document.querySelectorAll('.view')) s.hidden = s.dataset.view !== view;
  if (location.hash.slice(1) !== view) history.replaceState(null, '', `#${view}`);
  try {
    await VIEWS[view](state.data);
  } catch (error) {
    document.getElementById('status').textContent = error.message;
  }
}

async function engineBadge() {
  const badge = document.getElementById('engine');
  try {
    const e = await api.engineState();
    badge.className = `engine ${e.current ? 'ok' : 'stale'}`;
    badge.textContent = e.current ? 'engine ready' : e.built ? 'engine stale' : 'no engine';
    return e;
  } catch {
    badge.className = 'engine stale';
    badge.textContent = 'no server';
    return null;
  }
}

// The recipe, when the engine or the data are missing.
async function setup() {
  const box = document.getElementById('setup');
  const view = await api.setup();
  const steps = document.getElementById('setup-steps');
  clear(steps).append(...view.steps.map((s) => el('li', { class: s.done ? 'done' : s.state === 'running' ? 'running' : '', text: s.title })));
  const missing = view.steps.some((s) => !s.done);
  box.hidden = !missing;
  document.getElementById('setup-note').textContent = view.running ? 'making…' : view.writable === false ? view.refusal : '';
  if (view.running) setTimeout(setup, 3000);
  else if (!missing) await start();
  return missing;
}

async function start() {
  state.data = await api.dataState();
  await engineBadge();
  await show(viewFromHash());
}

document.getElementById('setup-run').addEventListener('click', async () => {
  try {
    await api.runSetup();
  } catch (error) {
    document.getElementById('setup-note').textContent = error.message;
  }
  setTimeout(setup, 1500);
});
for (const b of document.querySelectorAll('.views button')) b.addEventListener('click', () => show(b.dataset.view));
window.addEventListener('hashchange', () => show(viewFromHash()));

(async () => {
  const e = await engineBadge();
  let missing = false;
  try {
    missing = await setup();
  } catch {
    missing = false;
  }
  if (!missing && e) await start();
})();
