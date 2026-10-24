// The Economy Atlas: four views on the economics domain, computed by the engine through the server,
// on the library's data or one's own, with plug-ins in Python and C++ (ADR 0044, ADR 0046).

import * as api from './api.js';
import { el, clear } from './dom.js';
import { themeChanged } from './palette.js';
import { setDatasets } from './state.js';
import { showMap } from './map.js';
import { showSpacetime } from './spacetime.js';
import { showProduction } from './production.js';
import { showCircuits } from './circuits.js';

const VIEWS = { map: showMap, spacetime: showSpacetime, production: showProduction, circuits: showCircuits };
const state = { view: null, data: null };

// The address: the view, then what it shows, so that a picture can be linked and taken again:
// #map&year=2010&arcs=driven&economy=CHN, #spacetime&from=2015-01&to=2026-06&currency=MXN,
// #production&table=three, #circuits&network=five; and on the map or in spacetime, dataset=<id>.
function fromHash() {
  const [view, ...rest] = location.hash.slice(1).split('&');
  const params = Object.fromEntries(rest.filter((p) => p.includes('=')).map((p) => p.split('=').map(decodeURIComponent)));
  return { view: VIEWS[view] ? view : 'map', params };
}

async function show(view, params = {}) {
  state.view = view;
  for (const b of document.querySelectorAll('.views button')) b.setAttribute('aria-selected', String(b.dataset.view === view));
  for (const s of document.querySelectorAll('.view')) s.hidden = s.dataset.view !== view;
  if (location.hash.slice(1).split('&')[0] !== view) history.replaceState(null, '', `#${view}`);
  try {
    await VIEWS[view](state.data, params);
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

// Day and night: the system's unless chosen; the choice is kept in this browser.
const MOON = 'M11.8 10.4A5.2 5.2 0 0 1 5.6 4.2a5.4 5.4 0 1 0 6.2 6.2z';
const SUN = 'M8 4.6a3.4 3.4 0 1 0 0 6.8 3.4 3.4 0 0 0 0-6.8zM8 1v1.6M8 13.4V15M1 8h1.6M13.4 8H15M3 3l1.1 1.1M11.9 11.9 13 13M3 13l1.1-1.1M11.9 4.1 13 3';
const dark = window.matchMedia('(prefers-color-scheme: dark)');
const theme = () => document.documentElement.dataset.theme || (dark.matches ? 'night' : 'day');

function themeButton() {
  const b = document.getElementById('theme');
  const night = theme() === 'night';
  b.querySelector('path').setAttribute('d', night ? SUN : MOON);
  b.title = night ? 'By day' : 'By night';
  b.setAttribute('aria-label', b.title);
}

document.getElementById('theme').addEventListener('click', () => {
  const next = theme() === 'night' ? 'day' : 'night';
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem('economy-atlas-theme', next); } catch { /* kept for this visit only */ }
  themeButton();
  themeChanged();
});
dark.addEventListener('change', () => { if (!document.documentElement.dataset.theme) { themeButton(); themeChanged(); } });
themeButton();

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
  setDatasets(state.data.datasets || []);
  await engineBadge();
  const { view, params } = fromHash();
  await show(view, params);
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
window.addEventListener('hashchange', () => { const { view, params } = fromHash(); show(view, params); });

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
