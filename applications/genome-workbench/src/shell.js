// The Genome Workbench's shell: the rail of views (from the registry), the top bar (where you
// are, the command palette, the genome in hand), the workspace, and the status bar (the engine,
// the activity, what just happened). It reads and writes the address, mounts each view once and
// shows it on demand, takes files dropped anywhere into the library, and keeps the theme and the
// rail's width for this browser.
//
// A view is a module with mount(section, ctx), called once, and show(ctx), called whenever it is
// shown or the genome changes. The context it receives: the genome and the catalogue of genomes,
// the view's parameters from the address and setParams() to write them back, go() to another view,
// openTool(), setGenome(), status() for the status bar, and the engine's state.

import { $, el, clear } from './dom.js';
import * as api from './api.js';
import * as activity from './core/activity.js';
import * as registry from './core/registry.js';
import * as library from './core/library.js';
import * as palette from './core/palette.js';
import { readHash, writeHash } from './state.js';
import { setupPanel } from './setup.js';
import { icon } from './ui/icons.js';
import { iconBtn, popover, toast, fmtBytes } from './ui/kit.js';
import * as fmt from './format.js';

const store = {
  get(k, d) { try { return localStorage.getItem(`genome-workbench:${k}`) ?? d; } catch { return d; } },
  set(k, v) { try { localStorage.setItem(`genome-workbench:${k}`, v); } catch { /* private mode */ } },
};

const ctx = {
  view: 'genome', genome: null, params: {}, catalogue: [], engine: null, pendingTool: null,
  setParams(params, replace = true) {
    ctx.params = { ...ctx.params, ...params };
    for (const k of Object.keys(ctx.params)) if (ctx.params[k] === null || ctx.params[k] === undefined) delete ctx.params[k];
    const h = writeHash(ctx);
    if (h !== location.hash) (replace ? history.replaceState : history.pushState).call(history, null, '', h);
  },
  status(text) { const m = $('status-msg'); if (m) m.textContent = text || ''; },
  genomeInfo() { return ctx.catalogue.find((g) => g.id === ctx.genome) || null; },
  go(view, params = {}) {
    ctx.view = view;
    ctx.params = params;
    history.pushState(null, '', writeHash(ctx));
    show();
  },
  openTool(id, values = null) { ctx.pendingTool = { id, values }; ctx.go('toolbox', { tool: id }); },
  setGenome(id, toGenomeView = false) {
    ctx.genome = id;
    if (toGenomeView) { ctx.view = 'genome'; ctx.params = {}; }
    history.pushState(null, '', writeHash(ctx));
    drawDataset();
    show();
  },
  async refreshCatalogue() { await loadCatalogue(); drawDataset(); },
};
const mounted = new Set();
const sections = new Map();

// ---- the frame ---------------------------------------------------------------------------------------
function frame() {
  const shell = el('div', { class: `shell${store.get('rail', '') === 'collapsed' ? ' rail-collapsed' : ''}`, id: 'shell' });
  // the rail
  const nav = el('nav', { 'aria-label': 'Views' });
  for (const [g, label] of registry.RAIL_GROUPS) {
    const vs = registry.views.filter((v) => v.group === g);
    if (!vs.length) continue;
    nav.append(el('div', { class: 'rail-group' }, label));
    for (const v of vs) nav.append(el('button', { class: 'rail-item', 'data-view': v.id, title: v.title, onclick: () => { if (ctx.view !== v.id) ctx.go(v.id); } }, icon(v.icon), el('span', {}, v.title), v.count ? el('span', { class: 'count', id: `count-${v.id}` }) : null));
  }
  const theme = iconBtn(document.documentElement.dataset.theme === 'night' ? 'sun' : 'moon', { title: 'day or night', size: 'sm', onClick: () => setTheme(document.documentElement.dataset.theme === 'night' ? 'day' : 'night') });
  theme.id = 'theme-btn';
  const rail = el('aside', { class: 'rail' },
    el('div', { class: 'rail-brand' }, mark(), el('div', { class: 'name' }, el('b', {}, 'Genome Workbench'), el('i', {}, 'genomes, computed and drawn'))),
    nav,
    el('div', { class: 'rail-foot' },
      iconBtn('sidebar', { title: 'narrow or widen the rail', size: 'sm', onClick: () => { shell.classList.toggle('rail-collapsed'); store.set('rail', shell.classList.contains('rail-collapsed') ? 'collapsed' : ''); } }),
      theme,
      iconBtn('keyboard', { title: 'keys', size: 'sm', onClick: (e) => keysHelp(e.currentTarget) })));
  // the top bar
  const crumb = el('div', { class: 'crumb' }, el('h1', { id: 'crumb-title' }), el('span', { class: 'sub', id: 'crumb-sub' }));
  const top = el('header', { class: 'topbar' }, crumb, el('div', { class: 'spacer' }),
    el('button', { class: 'palette-trigger', onclick: () => palette.open(ctx) }, icon('search', 'sm'), el('span', {}, 'Find anything'), el('kbd', {}, 'Ctrl K')),
    el('button', { class: 'dataset', id: 'dataset', title: 'the genome in hand', onclick: (e) => datasetMenu(e.currentTarget) }));
  // the workspace
  const work = el('main', { class: 'workspace', id: 'workspace' });
  for (const v of registry.views) {
    const s = el('section', { id: `view-${v.id}`, class: `view ${v.full ? '' : 'page'} ${v.cls || ''}`.trim(), hidden: true });
    sections.set(v.id, s);
    work.append(s);
  }
  // the status bar
  const status = el('footer', { class: 'statusbar' },
    el('button', { id: 'engine-status', onclick: () => ctx.go('run') }, el('span', { class: 'dot' }), el('span', {}, 'engine')),
    el('button', { id: 'activity-btn', onclick: (e) => trayToggle(e.currentTarget) }, icon('activity', 'sm'), el('span', { id: 'activity-text' }, 'idle')),
    el('span', { class: 'msg', id: 'status-msg' }),
    el('span', { style: { marginLeft: 'auto' }, id: 'status-right' }));
  shell.append(rail, top, work, status);
  document.body.prepend(shell);
  activity.subscribe(drawActivity);
}

function mark() {
  const s = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  s.setAttribute('viewBox', '0 0 32 32');
  s.setAttribute('class', 'mark');
  s.innerHTML = '<circle cx="16" cy="16" r="10.5" fill="none" stroke="currentColor" stroke-opacity="0.75" stroke-width="1.3"/><path d="M16 3.6A12.4 12.4 0 0 1 28.1 13.4" fill="none" stroke="#a8577e" stroke-width="2.6" stroke-linecap="round"/><path d="M26.6 21.6A12.4 12.4 0 0 1 9.8 26.8" fill="none" stroke="#c9a043" stroke-width="2.6" stroke-linecap="round"/><path d="M8.2 16.6A7.8 7.8 0 0 1 13.4 8.6" fill="none" stroke="#6e9c7a" stroke-width="2" stroke-linecap="round"/>';
  return s;
}

function setTheme(t) {
  if (t === 'night') document.documentElement.dataset.theme = 'night';
  else delete document.documentElement.dataset.theme;
  store.set('theme', t);
  const b = $('theme-btn');
  if (b) b.replaceChildren(icon(t === 'night' ? 'sun' : 'moon', 'sm'));
  // Drawings read colours from the tokens: redraw the view in hand.
  mounted.delete(ctx.view);
  show();
}

function keysHelp(anchor) {
  const rows = [['Ctrl K or /', 'find anything'], ['Esc', 'close'], ['drop files', 'add them to your library'], ['wheel on the strip', 'zoom the genome'], ['drag on the circle', 'choose a stretch']];
  popover(anchor, el('div', { style: { padding: '6px 8px' } }, el('div', { class: 'menu-group' }, 'keys and gestures'),
    ...rows.map(([k, v]) => el('div', { class: 'row', style: { justifyContent: 'space-between', padding: '4px 4px', gap: '24px' } }, el('kbd', {}, k), el('span', { class: 'small faint' }, v)))), { place: 'above' });
}

// ---- the genome in hand ------------------------------------------------------------------------------
const GROUPS = [['bacteria', 'Bacteria'], ['viruses', 'Viruses'], ['mitochondria', 'Mitochondria'], ['files', 'Your files'], ['simulated', 'Simulated, with their truth']];
const shortTitle = (g) => g.title.replace(/, complete genome.*$/, '').replace(/ \(.*\)$/, '');

function drawDataset() {
  const b = $('dataset');
  if (!b) return;
  const g = ctx.genomeInfo();
  b.replaceChildren(el('span', { class: 'glyph' }, icon(g?.group === 'files' ? 'file' : 'genome', 'sm')),
    el('span', { class: 'what' }, el('b', {}, g ? shortTitle(g) : 'no genome'), el('i', {}, g ? `${g.group === 'files' ? 'your file' : g.group}${g.size ? ` · ${fmtBytes(g.size)}` : ''}` : 'choose one')), icon('down', 'sm'));
}

function datasetMenu(anchor) {
  const list = el('div', { style: { maxHeight: '64vh', overflowY: 'auto' } });
  const p = popover(anchor, list, { width: 420 });
  for (const [key, label] of GROUPS) {
    const items = ctx.catalogue.filter((g) => g.group === key);
    if (!items.length) continue;
    list.append(el('div', { class: 'menu-group' }, label));
    for (const g of items) {
      list.append(el('button', { class: `menu-item${g.id === ctx.genome ? ' active' : ''}`, disabled: !g.fetched || null, onclick: () => { p.close(); ctx.setGenome(g.id); } },
        icon(key === 'files' ? 'file' : key === 'simulated' ? 'spark' : 'genome', 'sm'),
        el('span', { style: { minWidth: 0, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' } }, shortTitle(g)),
        el('span', { class: 'meta' }, g.fetched ? (g.accession || g.kind || '') : 'not fetched')));
    }
  }
  list.append(el('div', { class: 'menu-sep' }),
    el('button', { class: 'menu-item', onclick: () => { p.close(); ctx.go('data'); } }, icon('data', 'sm'), el('span', {}, 'Your data: add a genome from a file…')));
}

async function loadCatalogue() {
  const g = await api.genomes().catch(() => ({ genomes: [] }));
  ctx.catalogue = g.genomes;
  if (!ctx.genome || !ctx.catalogue.some((x) => x.id === ctx.genome && x.fetched)) {
    ctx.genome = (['ncbi:NC_000913.3', 'ncbi:NC_001422.1'].map((id) => ctx.catalogue.find((x) => x.id === id && x.fetched)).find(Boolean) || ctx.catalogue.find((x) => x.fetched))?.id || null;
  }
}

// ---- the engine and the activity ---------------------------------------------------------------------
function drawEngine() {
  const b = $('engine-status'), e = ctx.engine;
  if (!b) return;
  const workers = e?.cluster?.workers?.filter((w) => w.alive) || [];
  const target = e?.target?.built_for || workers[0]?.built_for;
  const state = !e?.built ? 'off' : e.current ? 'ok' : 'warn';
  b.replaceChildren(el('span', { class: `dot ${state}` }),
    el('span', {}, !e?.built ? 'engine not built' : `${e.current ? 'engine' : 'engine (sources changed)'}${target ? ` · ${target.replace('avx512bw', 'AVX-512').replace('avx2', 'AVX2')}` : ''}${workers.length ? ` · ${workers.length} worker${workers.length > 1 ? 's' : ''}` : ''}`));
  const right = $('status-right');
  if (right && e?.target) right.textContent = `${e.target.threads || ''} threads · ${e.target.brand?.replace(/\(R\)|\(TM\)|CPU|@.*$/g, '').replace(/\s+/g, ' ').trim() || ''}`;
}

async function refreshEngine() {
  ctx.engine = await api.engineState().catch(() => null);
  drawEngine();
}

function drawActivity() {
  const t = $('activity-text'), dot = $('activity-btn');
  if (!t) return;
  const run = activity.running();
  const last = activity.list().find((e) => e.ended !== null);
  t.textContent = run.length ? `${run.length} running: ${run[0].op}` : last ? `${last.op} ${last.ok ? `in ${fmt.seconds((last.ended - last.started) / 1000)}` : 'failed'}` : 'idle';
  dot.classList.toggle('busy', run.length > 0);
  if (tray) fillTray();
}

let tray = null;
function trayToggle() {
  if (tray) { tray.remove(); tray = null; return; }
  tray = el('div', { class: 'tray' });
  document.body.append(tray);
  fillTray();
  setTimeout(() => document.addEventListener('pointerdown', function away(e) {
    if (tray && !tray.contains(e.target) && !$('activity-btn').contains(e.target)) { tray.remove(); tray = null; document.removeEventListener('pointerdown', away, true); }
  }, true));
}
function fillTray() {
  const items = activity.list();
  tray.replaceChildren(el('div', { class: 'tray-head' }, icon('activity', 'sm'), el('h3', {}, 'Activity'), el('span', { class: 'small faint' }, `${items.length} requests this session`), el('span', { style: { flex: 1 } }),
    iconBtn('x', { title: 'close', size: 'sm', onClick: () => { tray.remove(); tray = null; } })));
  const list = el('div', { class: 'tray-list' });
  for (const e of items.slice(0, 80)) {
    const dur = e.ended === null ? performance.now() - e.started : e.ended - e.started;
    list.append(el('div', { class: 'tray-item', title: e.error || '' },
      e.ended === null ? el('span', { class: 'spin' }) : el('span', { class: `dot ${e.ok ? 'ok' : 'off'}` }),
      el('span', { class: 'op' }, `${e.kind === 'job' ? 'job ' : ''}${e.op}`),
      el('span', { class: 'when' }, `${fmt.seconds(dur / 1000)}${e.worker ? ` · ${e.worker}` : ''}${e.shards ? ` · ${e.shards} shards` : ''}`),
      el('span', { class: 'about' }, e.error ? `failed: ${e.error}` : e.about || '—')));
  }
  if (!items.length) list.append(el('p', { class: 'note', style: { padding: '10px' } }, 'Nothing has been asked of the engine yet.'));
  tray.append(list);
}

// ---- showing a view ----------------------------------------------------------------------------------
function show() {
  const def = registry.viewById(ctx.view) || registry.views[0];
  ctx.view = def.id;
  for (const [id, s] of sections) s.hidden = id !== def.id;
  for (const b of document.querySelectorAll('.rail-item')) {
    if (b.dataset.view === def.id) b.setAttribute('aria-current', 'page');
    else b.removeAttribute('aria-current');
  }
  $('crumb-title').textContent = def.title;
  $('crumb-sub').textContent = def.tagline || '';
  document.title = `${def.title} · Genome Workbench`;
  const section = sections.get(def.id);
  if (def.needsEngine !== false && !ctx.engine?.built) {
    clear(section).append(el('div', { class: 'page-narrow' }, setupPanel(async () => { await refreshEngine(); await loadCatalogue(); drawDataset(); mounted.clear(); show(); })));
    mounted.delete(def.id);
    return;
  }
  if (!mounted.has(def.id)) {
    clear(section);
    def.module.mount(section, ctx);
    adoptHead(section, def);
    mounted.add(def.id);
  }
  def.module.show(ctx);
}

// Views written before the shell open with their own title and long text: the title is the top
// bar's now, and the long text is kept behind the view's "about".
function adoptHead(section, def) {
  const h = section.querySelector(':scope > h2.view-title');
  const lede = section.querySelector(':scope > p.lede');
  if (h) h.remove();
  if (lede) {
    lede.classList.remove('lede');
    lede.classList.add('about');
    lede.hidden = true;
    const b = iconBtn('info', { title: 'about this view', size: 'sm', onClick: () => { lede.hidden = !lede.hidden; b.setAttribute('aria-pressed', String(!lede.hidden)); } });
    b.style.marginLeft = '6px';
    b.dataset.forView = def.id;
    section.prepend(lede);
    section.aboutButton = b;
  }
  const crumb = document.querySelector('.crumb');
  crumb.querySelectorAll('[data-for-view]').forEach((x) => x.remove());
  if (section.aboutButton) crumb.append(section.aboutButton);
}

// ---- files dropped anywhere --------------------------------------------------------------------------
function dropAnywhere() {
  let veil = null, depth = 0;
  const hasFiles = (e) => [...(e.dataTransfer?.types || [])].includes('Files');
  window.addEventListener('dragenter', (e) => {
    if (!hasFiles(e)) return;
    depth++;
    if (!veil) {
      veil = el('div', { class: 'drop-veil' }, el('div', { class: 'card' }, icon('upload', 'lg'), el('h2', {}, 'Drop to add to your data'),
        el('p', { class: 'small faint', style: { margin: 0 } }, 'FASTA, FASTQ, GenBank, VCF, Newick, SAM, plain or gzipped')));
      document.body.append(veil);
    }
  });
  window.addEventListener('dragleave', () => { depth = Math.max(0, depth - 1); if (!depth && veil) { veil.remove(); veil = null; } });
  window.addEventListener('dragover', (e) => { if (hasFiles(e)) e.preventDefault(); });
  window.addEventListener('drop', async (e) => {
    depth = 0;
    if (veil) { veil.remove(); veil = null; }
    if (!hasFiles(e)) return;
    e.preventDefault();
    const added = await library.add([...e.dataTransfer.files]);
    if (added.some((r) => r.genome)) await ctx.refreshCatalogue();
    if (added.length === 1 && added[0].genome) {
      const t = toast(el('span', {}, `${added[0].name} can be read as a genome. `, el('button', { class: 'link', onclick: () => { t.remove(); ctx.setGenome(added[0].id, true); } }, 'Open it')), 'good', 9000);
    }
    if (ctx.view === 'data') { mounted.delete('data'); show(); }
  });
}

// ---- keys ----------------------------------------------------------------------------------------------
function keys() {
  window.addEventListener('keydown', (e) => {
    const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement?.tagName || '') || document.activeElement?.isContentEditable;
    if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'k') { e.preventDefault(); palette.open(ctx); }
    else if (e.key === '/' && !typing && !document.querySelector('.dialog-veil')) { e.preventDefault(); palette.open(ctx); }
  });
}

// ---- start ----------------------------------------------------------------------------------------------
export async function start() {
  if (store.get('theme', 'day') === 'night') document.documentElement.dataset.theme = 'night';
  frame();
  keys();
  dropAnywhere();
  const h = readHash();
  ctx.view = h.view;
  ctx.genome = h.genome;
  ctx.params = h.params;
  await Promise.all([refreshEngine(), loadCatalogue(), library.ready()]);
  drawDataset();
  window.addEventListener('popstate', () => {
    const s = readHash();
    const changedGenome = s.genome !== ctx.genome;
    ctx.view = s.view;
    ctx.genome = s.genome || ctx.genome;
    ctx.params = s.params;
    if (changedGenome) drawDataset();
    show();
  });
  show();
  setInterval(refreshEngine, 15000);
}

export { ctx };
