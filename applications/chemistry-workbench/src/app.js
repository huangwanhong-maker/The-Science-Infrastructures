// The Chemistry Workbench's shell. Every plug-in (plugins/<name>/plugin.json) declares its views; the
// rail lists them, the address (#<plugin>/<view>[/<args>…]) chooses one, and its module is loaded and given a
// context: the plug-in's engine (ask), the drawing kit and the page's helpers. A plug-in joins the
// workbench by being a directory with a manifest; nothing here names one.

import { asker, plugins, runSetup, setup } from './api.js';
import * as dom from './dom.js';
import * as chart from './draw/chart.js';
import * as graph from './draw/graph.js';
import * as complex from './draw/complex.js';

const { h } = dom;
const state = { plugins: [], current: null, disposers: [] };

function theme(set) {
  let t = set;
  if (!t) {
    try { t = localStorage.getItem('chemistry-theme'); } catch { t = null; }
  }
  t = t === 'night' ? 'night' : 'day';
  document.documentElement.dataset.theme = t;
  document.getElementById('theme').textContent = t === 'night' ? 'day' : 'night';
  try { localStorage.setItem('chemistry-theme', t); } catch { /* private window */ }
  return t;
}

function rail() {
  const r = document.getElementById('rail');
  r.replaceChildren();
  for (const p of state.plugins) {
    const block = h('div', { class: 'plugin' }, h('div', { class: 'plugin-title' }, p.title), h('div', { class: 'plugin-summary' }, p.summary || ''));
    for (const v of p.views || []) block.append(h('a', { class: 'view', href: `#${p.name}/${v.id}`, 'data-key': `${p.name}/${v.id}`, title: v.about || '' }, v.title));
    r.append(block);
  }
  r.append(h('div', { class: 'about' }, 'Plug-ins live in plugins/<name>/, each with a manifest of its views and, if it computes, an engine in C++ on the foundations.'));
}

function lamp() {
  const el = document.getElementById('lamp');
  const engines = state.plugins.filter((p) => p.engine);
  const current = engines.every((p) => p.engineStatus && p.engineStatus.current);
  el.className = `lamp ${current ? 'on' : 'off'}`;
  el.textContent = current ? 'engines current' : 'engines to build';
}

async function buildPanel(main, plugin) {
  const log = h('pre', { class: 'mono', style: 'white-space: pre-wrap; max-height: 320px; overflow: auto; font-size: 0.8rem' });
  const button = h('button', { class: 'primary' }, 'Build the engines');
  const panel = h('div', { class: 'panel', style: 'max-width: 760px' },
    h('h2', {}, 'The engine is not built from the present sources'),
    h('p', { class: 'muted' }, `The plug-in ${plugin.title} computes in C++ on the foundations. Building compiles it with g++ (natively, or in WSL), which takes a minute or two.`),
    h('p', {}, button), log);
  main.replaceChildren(panel);
  const poll = async () => {
    const s = await setup();
    log.textContent = s.log || '';
    if (s.running) return setTimeout(poll, 1500);
    state.plugins = await plugins();
    lamp();
    const p = state.plugins.find((q) => q.name === plugin.name);
    if (p && p.engineStatus && p.engineStatus.current) route();
    else button.disabled = false;
  };
  button.addEventListener('click', async () => {
    button.disabled = true;
    const s = await runSetup();
    if (s.message) log.textContent = s.message;
    poll();
  });
  poll();
}

async function route() {
  for (const d of state.disposers.splice(0)) try { d(); } catch { /* already gone */ }
  const main = document.getElementById('main');
  const [name, view, ...args] = (location.hash.slice(1) || '').split('/');
  const plugin = state.plugins.find((p) => p.name === name) || state.plugins[0];
  if (!plugin) return main.replaceChildren(h('p', {}, 'No plug-ins are installed.'));
  const v = (plugin.views || []).find((x) => x.id === view) || plugin.views[0];
  for (const a of document.querySelectorAll('.rail a.view')) a.classList.toggle('active', a.dataset.key === `${plugin.name}/${v.id}`);
  document.getElementById('where').textContent = `${plugin.title} · ${v.title}`;
  if (plugin.engine && !(plugin.engineStatus && plugin.engineStatus.current)) return buildPanel(main, plugin);
  main.replaceChildren(h('p', { class: 'busy' }, 'loading…'));
  try {
    const mod = await import(`/plugins/${plugin.name}/${v.module}`);
    const root = h('div', {});
    main.replaceChildren(root);
    await mod.render(root, {
      plugin, view: v, args, ask: asker(plugin.name), dom, chart, graph, complex,
      base: `/plugins/${plugin.name}/`,
      onDispose: (f) => state.disposers.push(f),
      go: (id) => { location.hash = `${plugin.name}/${id}`; },
    });
  } catch (error) {
    dom.failure(main, error);
    console.error(error);
  }
}

async function start() {
  theme();
  document.getElementById('theme').addEventListener('click', () => {
    theme(document.documentElement.dataset.theme === 'night' ? 'day' : 'night');
    route();
  });
  try {
    state.plugins = await plugins();
  } catch (error) {
    return dom.failure(document.getElementById('main'), error);
  }
  rail();
  lamp();
  window.addEventListener('hashchange', route);
  route();
}

start();
