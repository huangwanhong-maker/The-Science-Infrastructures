// The Genome Atlas: the shell. It reads the address, asks the server for the engine and the
// genomes, fills the header, and shows one view at a time; each view is a module with
// mount(section, ctx), called once, and show(ctx), called whenever it becomes visible or the
// genome changes. The context carries the genome, the view's parameters from the address, and
// setParams(), which writes them back.

import { $, el, clear, notice } from './dom.js';
import * as api from './api.js';
import { readHash, writeHash, VIEWS } from './state.js';
import { setupPanel } from './setup.js';
import * as genome from './views/genome.js';
import * as run from './views/run.js';
import * as align from './views/align.js';
import * as indexView from './views/index.js';
import * as assembly from './views/assembly.js';
import * as populations from './views/populations.js';
import * as trees from './views/trees.js';
import * as readsView from './views/reads.js';
import * as expression from './views/expression.js';

const MODULES = {
  genome,
  align,
  index: indexView,
  reads: readsView,
  expression,
  assembly,
  trees,
  populations,
  run,
};

const ctx = {
  view: 'genome', genome: null, params: {}, catalogue: [], engine: null,
  setParams(params, replace = true) {
    ctx.params = { ...ctx.params, ...params };
    for (const k of Object.keys(ctx.params)) if (ctx.params[k] === null || ctx.params[k] === undefined) delete ctx.params[k];
    const h = writeHash(ctx);
    if (h !== location.hash) (replace ? history.replaceState : history.pushState).call(history, null, '', h);
  },
  status(text) { $('status').textContent = text || ''; },
  genomeInfo() { return ctx.catalogue.find((g) => g.id === ctx.genome) || null; },
};
const mounted = new Set();

function fillGenomes() {
  const select = clear($('genome'));
  const groups = [['bacteria', 'Bacteria'], ['viruses', 'Viruses'], ['mitochondria', 'Mitochondria'], ['simulated', 'Simulated, with their truth']];
  for (const [key, label] of groups) {
    const items = ctx.catalogue.filter((g) => g.group === key);
    if (!items.length) continue;
    const og = el('optgroup', { label });
    for (const g of items) og.append(el('option', { value: g.id, disabled: !g.fetched }, `${g.title.replace(/, complete genome.*$/, '').replace(/ \(.*\)$/, '')}${g.fetched ? '' : ' (not fetched)'}`));
    select.append(og);
  }
  select.value = ctx.genome || '';
  select.onchange = () => {
    ctx.genome = select.value;
    ctx.params = {};
    history.pushState(null, '', writeHash(ctx));
    show();
  };
}

function engineBadge() {
  const b = $('engine'), e = ctx.engine;
  const workers = e?.cluster?.workers?.filter((w) => w.alive) || [];
  const target = e?.target?.built_for || workers[0]?.built_for;
  b.className = `engine ${!e?.built ? 'off' : e.current ? 'ok' : 'stale'}`;
  b.textContent = !e?.built ? 'engine not built' : `${e.current ? 'engine' : 'engine (stale)'}${target ? ` · ${target.replace('avx512bw', 'AVX-512').replace('avx2', 'AVX2')}` : ''}${workers.length ? ` · ${workers.length} worker${workers.length > 1 ? 's' : ''}` : ''}`;
  b.onclick = () => { ctx.view = 'run'; ctx.params = {}; history.pushState(null, '', writeHash(ctx)); show(); };
}

async function refreshEngine() {
  ctx.engine = await api.engineState().catch(() => null);
  engineBadge();
}

function show() {
  for (const v of VIEWS) {
    $(`view-${v}`).hidden = v !== ctx.view;
    document.querySelector(`.views button[data-view="${v}"]`)?.setAttribute('aria-selected', String(v === ctx.view));
  }
  $('genome').value = ctx.genome || '';
  const section = $(`view-${ctx.view}`);
  const needsEngine = ctx.view !== 'run' && !ctx.engine?.built;
  if (needsEngine) {
    section.hidden = false;
    clear(section).append(el('div', { class: 'page' }, setupPanel(async () => { await refreshEngine(); await loadCatalogue(); mounted.clear(); show(); })));
    mounted.delete(ctx.view);
    return;
  }
  const mod = MODULES[ctx.view];
  if (!mounted.has(ctx.view)) {
    clear(section);
    mod.mount(section, ctx);
    mounted.add(ctx.view);
  }
  mod.show(ctx);
}

async function loadCatalogue() {
  const g = await api.genomes().catch(() => ({ genomes: [] }));
  ctx.catalogue = g.genomes;
  if (!ctx.genome || !ctx.catalogue.some((x) => x.id === ctx.genome && x.fetched)) {
    ctx.genome = (['ncbi:NC_000913.3', 'ncbi:NC_001422.1'].map((id) => ctx.catalogue.find((x) => x.id === id && x.fetched)).find(Boolean) || ctx.catalogue.find((x) => x.fetched))?.id || null;
  }
  fillGenomes();
}

async function init() {
  const h = readHash();
  ctx.view = h.view;
  ctx.genome = h.genome;
  ctx.params = h.params;
  for (const b of document.querySelectorAll('.views button')) {
    b.onclick = () => {
      if (ctx.view === b.dataset.view) return;
      ctx.view = b.dataset.view;
      ctx.params = {};
      history.pushState(null, '', writeHash(ctx));
      show();
    };
  }
  await Promise.all([refreshEngine(), loadCatalogue()]);
  window.addEventListener('popstate', () => {
    const s = readHash();
    const changedGenome = s.genome !== ctx.genome;
    ctx.view = s.view;
    ctx.genome = s.genome || ctx.genome;
    ctx.params = s.params;
    if (changedGenome) $('genome').value = ctx.genome;
    show();
  });
  show();
  setInterval(refreshEngine, 15000);
}

init().catch((e) => notice(`The atlas could not start: ${e.message}`, true, 20000));
