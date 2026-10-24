// The Genome view: a genome as a circle and as a strip. The circle shows the whole: features by
// strand, GC content and skew, and the cumulative skew whose extremes mark the origin and the
// terminus of replication. The strip shows a window: features in rows, the six frames' ORFs,
// composition in bins, the letters at base resolution, the hits of a search. A column of cards
// shows the genome, or the feature chosen, with its protein checked against the record's own and
// its codons.
//
// The toolbar above switches the strip's tracks (each chip carries the colours its track is drawn
// in), computes what the sequence alone says (genes, CpG islands, repeats), and sets the shortest
// ORF. The strip's bar holds the window, its navigation and the search, which takes a gene, a
// sequence of four or more bases, a position or a range.
//
// Drag on the circle to choose a window, or click it to move there; on the strip, the wheel
// zooms about the pointer and dragging pans; + and − zoom, the arrows pan, 0 shows the whole, f
// finds, Esc returns to the genome's card. The window and the chosen feature live in the address
// (begin, end, contig, feature, orf); a search asked from elsewhere arrives as find.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { hashNumber } from '../state.js';
import { circleSvg, positionAt, angleOf } from '../draw/circle.js';
import { linearSvg, LEFT } from '../draw/linear.js';
import * as fmt from '../format.js';
import { featureColour, warm, token, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, iconBtn, select as selectBox, card, metrics, kv, callout, empty, popover, menu, copy, toast } from '../ui/kit.js';

const ask = api.cached(64);
let view = null;   // the view's elements and state
let keysBound = false;

const LEDE = 'A genome as a circle and as a strip. The circle shows the whole: its features by strand, GC content and skew, and the cumulative skew, whose extremes mark the origin and the terminus of replication in most bacteria. The strip shows a window: the features in rows, the open reading frames of the six frames, composition in bins, and the letters themselves once there is room. Drag on the circle to choose a stretch, or click it to move there; on the strip the wheel zooms about the pointer and a drag pans. A feature clicked opens its card, with its protein checked against the record’s own and the codons it uses.';

// The kit's icons under the names this view uses for them.
const GLYPHS = { left: 'left', right: 'chevron', up: 'up', down: 'down', zoomIn: 'zoom-in', zoomOut: 'zoom-out', fit: 'fit', target: 'target' };
const glyph = (name) => icon(GLYPHS[name], 'sm');
const glyphBtn = (name, title, onClick) => iconBtn(GLYPHS[name], { title, size: 'sm', onClick });
const withGlyph = (b, name) => { b.prepend(glyph(name)); return b; };

// The strip's tracks, each shown with the colours draw/linear.js draws it in.
const split = (a, b) => `linear-gradient(135deg, ${a} 50%, ${b} 50%)`;
const TRACKS = [
  ['orfs', 'ORFs', 'open reading frames in the six frames', 'color-mix(in srgb, var(--accent) 45%, transparent)'],
  ['gc', 'GC', 'GC content about the genome’s mean: ochre above it, blue below', split('var(--high)', 'var(--low)')],
  ['skew', 'skew', 'GC skew, (G − C)/(G + C): mulberry where G leads, sage where C leads', split('var(--pos)', 'var(--neg)')],
  ['cpg', 'CpG o/e', 'CpG dinucleotides observed over expected', 'var(--ncrna)'],
  ['entropy', 'entropy', 'the entropy of 3-mers in each bin', 'var(--C)'],
];
const ANNOTATIONS = [
  ['genes', 'predicted genes', 'a self-trained gene finder: periodic Markov chains of coding', 'color-mix(in srgb, var(--cds) 40%, transparent)'],
  ['islands', 'CpG islands', 'the HMM of Durbin et al.', 'var(--G)'],
  ['repeats', 'repeats ≥ 100 bp', 'maximal repeats on both strands, by the suffix array (drawn as chords on the circle)', `conic-gradient(${CATEGORIES.slice(0, 4).map((c, k) => `${c} ${k * 25}% ${(k + 1) * 25}%`).join(', ')})`],
];
const ORF_MIN = [30, 50, 100, 200, 300];

export function mount(section, ctx) {
  const tools = el('div', { class: 'gv-toolbar', role: 'toolbar', 'aria-label': 'tracks of the strip' });
  const body = el('div', { class: 'genome-body' });
  const circlePane = el('div', { class: 'circle-pane' });
  const linearPane = el('div', { class: 'linear-pane' });
  const panel = el('aside', { class: 'card-panel', 'aria-live': 'polite' });
  const bar = el('div', { class: 'linear-bar' });
  const stage = el('div', { class: 'linear-stage' });
  linearPane.append(bar, stage);
  body.append(circlePane, linearPane, panel);
  section.append(el('p', { class: 'lede' }, LEDE), tools, body);
  view = { ctx, section, tools, circlePane, linearPane, bar, stage, panel, genome: null, summary: null, whole: null, tracks: null, contig: 0,
           window: null, selected: null, hits: [], hitAt: -1, hitTotal: 0, query: '', keepHits: false,
           show: { orfs: true, gc: true, skew: true, cpg: true, entropy: false }, request: 0, layout: null,
           annotate: { genes: false, islands: false, repeats: false }, counts: {}, busy: {}, extra: {} };
  new ResizeObserver(() => drawLinear()).observe(stage);
  stage.addEventListener('wheel', onWheel, { passive: false });
  stage.addEventListener('pointerdown', onStripDown);
  stage.addEventListener('pointermove', onStripHover);
  stage.addEventListener('pointerleave', () => hover(null));
  circlePane.addEventListener('pointerdown', onCircleDown);
  circlePane.addEventListener('pointermove', onCircleHover);
  circlePane.addEventListener('pointerleave', () => hover(null));
  if (!keysBound) { window.addEventListener('keydown', onKey); keysBound = true; }
  buildTools();
  buildBar();
}

export async function show(ctx) {
  view.ctx = ctx;
  if (!ctx.genome) { noGenome(); return; }
  const contig = hashNumber(ctx.params, 'contig', 0, 0);
  if (view.genome !== ctx.genome || view.contig !== contig) {
    view.genome = ctx.genome;
    view.contig = contig;
    // A search stepping onto another contig keeps its hits; anything else starts afresh.
    if (!view.keepHits) { view.hits = []; view.hitTotal = 0; view.hitAt = -1; }
    view.keepHits = false;
    view.selected = null;
    view.tracks = null;
    view.extra = {};
    view.counts = {};
    for (const k of Object.keys(view.annotate)) view.annotate[k] = false;
    await load();
  }
  if (!view.whole || view.genome !== ctx.genome) return;
  const L = view.length;
  const begin = hashNumber(ctx.params, 'begin', 0, 0, L - 1), end = hashNumber(ctx.params, 'end', L, begin + 1, L);
  view.window = { begin: Math.floor(begin), end: Math.ceil(end) };
  const feature = ctx.params.feature !== undefined ? Number(ctx.params.feature) : null;
  if (feature !== view.selected) select(Number.isFinite(feature) ? feature : null, false);
  drawCircle();
  updateBar();
  fetchWindow();
  // A search asked from elsewhere (the command palette), run here once.
  if (ctx.params.find) {
    const q = String(ctx.params.find);
    ctx.setParams({ find: null });
    view.searchInput.value = q;
    syncClear();
    search(q);
  }
}

function noGenome() {
  clear(view.circlePane).append(empty({ icon: 'genome', title: 'No genome in hand', text: 'Fetch the catalogue’s genomes, or add your own FASTA or GenBank file: it opens here.',
    action: btn('Your data', { icon: 'data', kind: 'primary', onClick: () => view.ctx.go('data') }) }));
  clear(view.stage);
  clear(view.panel);
}

// ---- loading ----------------------------------------------------------------------------------------
async function load() {
  const id = view.genome;
  view.ctx.status('reading the genome…');
  view.summary = null;
  view.whole = null;
  clear(view.circlePane).append(el('div', { class: 'gv-ring' }, el('div', { class: 'wait' }, el('span', { class: 'spin' }), 'Reading the genome…')));
  clear(view.stage).append(stripSkeleton());
  clear(view.panel).append(cardSkeleton());
  buildTools();
  try {
    const summary = await ask('genome.summary', { genome: id });
    if (view.genome !== id) return;
    view.summary = summary;
    const c = summary.contigs[view.contig] || summary.contigs[0];
    view.length = c.length;
    const whole = await ask('genome.tracks', { genome: id, contig: view.contig, bins: 720, orfs: false, max_features: 20000 });
    if (view.genome !== id) return;
    view.whole = whole;
    view.ctx.status(`read in ${fmt.seconds(summary._seconds)}`);
    buildTools();
    if (view.selected === null) showGenomeCard();
  } catch (e) {
    if (view.genome !== id) return;
    view.ctx.status('');
    view.genome = null;   // asked again on the next showing
    clear(view.circlePane).append(el('div', { class: 'gv-pad' }, callout(`The genome could not be read: ${e.message}`, 'bad')));
    clear(view.stage);
    clear(view.panel);
  }
}

let windowTimer = null;
function fetchWindow() {
  clearTimeout(windowTimer);
  view.linearPane.classList.add('loading');
  windowTimer = setTimeout(async () => {
    const id = view.genome, { begin, end } = view.window;
    const width = Math.max(200, view.stage.clientWidth - LEFT - 14);
    const bins = Math.max(20, Math.min(4000, Math.floor(width / 2)));
    const n = ++view.request;
    try {
      const t = await ask('genome.tracks', { genome: id, contig: view.contig, begin, end, bins, orf_min: Number(view.ctx.params.orf || 100), max_features: 12000 });
      if (n !== view.request) return;
      view.tracks = t;
      drawLinear();
    } catch (e) {
      if (n === view.request) toast(e.message, 'error');
    } finally {
      if (n === view.request) view.linearPane.classList.remove('loading');
    }
  }, 90);
}

const stripSkeleton = () => el('div', { class: 'gv-strip-skel', 'aria-hidden': 'true' },
  ...[12, 34, 54, 46, 38, 30].map((h, k) => el('div', { class: 'skeleton', style: { height: `${h}px`, opacity: String(1 - k * 0.1) } })));
const cardSkeleton = () => el('div', { class: 'card gv-card-skel', 'aria-hidden': 'true' },
  ...[['30%', 10], ['78%', 20], ['56%', 12], ['100%', 44], ['90%', 12], ['70%', 12], ['84%', 12]].map(([w, h]) => el('div', { class: 'skeleton', style: { width: w, height: `${h}px` } })));

// ---- the toolbar ------------------------------------------------------------------------------------
// A chip that switches a track: its swatch is the track's colours, hollow while it is off.
function chip(label, title, swatch, on, onClick) {
  const sw = el('span', { class: 'gv-sw', style: { background: swatch } });
  const b = el('button', { type: 'button', class: `chip gv-chip${on ? ' on' : ''}`, title, 'aria-pressed': String(on), onclick: () => onClick(b) }, sw, el('span', {}, label), el('span', { class: 'gv-count' }));
  b.swatch = sw;
  b.set = (v) => { b.classList.toggle('on', v); b.setAttribute('aria-pressed', String(v)); };
  return b;
}

const group = (label, ...content) => el('div', { class: 'gv-group' }, el('span', { class: 'gv-label' }, label), ...content);
const sep = () => el('span', { class: 'gv-sep', 'aria-hidden': 'true' });

function buildTools() {
  const t = clear(view.tools);
  const s = view.summary;
  if (s && s.contigs.length > 1) {
    t.append(group('contig', selectBox(s.contigs.map((c, k) => [k, `${c.name} (${fmt.length(c.length)})`]), view.contig,
      (v) => { view.ctx.setParams({ contig: v, begin: null, end: null, feature: null }); show(view.ctx); }, { title: 'the contig shown' })), sep());
  }
  // The ORFs' chip carries their shortest length beside it, chosen from a menu.
  const orfMin = () => Number(view.ctx.params.orf || 100);
  const orfOpt = el('button', { type: 'button', class: 'chip gv-chip-opt', title: 'the shortest open reading frame drawn', onclick: () => menu(orfOpt,
    [{ group: 'ORFs of at least' }, ...ORF_MIN.map((v) => ({ label: `${v} codons`, icon: v === orfMin() ? 'check' : null, meta: `${3 * v + 3} bp`, onClick: () => {
      view.ctx.setParams({ orf: v });
      orfOpt.firstChild.textContent = `≥ ${v}`;
      fetchWindow();
    } }))], { width: 250 }) }, el('span', {}, `≥ ${orfMin()}`), icon('down', 'sm'));
  const trackChips = TRACKS.map(([k, label, title, sw]) => chip(label, title, sw, view.show[k], (b) => {
    view.show[k] = !view.show[k];
    b.set(view.show[k]);
    if (k === 'orfs') orfOpt.classList.toggle('dim', !view.show.orfs);
    drawLinear();
  }));
  if (!view.show.orfs) orfOpt.classList.add('dim');
  t.append(group('tracks', el('div', { class: 'chips' }, el('span', { class: 'gv-split' }, trackChips[0], orfOpt), ...trackChips.slice(1))), sep());
  t.append(group('from the sequence', el('div', { class: 'chips' }, ...ANNOTATIONS.map(([k, label, title, sw]) => {
    const b = chip(label, title, sw, view.annotate[k], (b) => { view.annotate[k] = !view.annotate[k]; b.set(view.annotate[k]); annotate(k, b); });
    if (view.annotate[k] && view.counts[k] !== undefined) b.querySelector('.gv-count').textContent = fmt.count(view.counts[k]);
    if (view.busy[k]) busyChip(b, true);
    b.dataset.kind = k;
    return b;
  }))));
  t.append(el('span', { class: 'gv-spacer' }),
    iconBtn('keyboard', { title: 'keys and gestures', size: 'sm', onClick: (e) => keysHelp(e.currentTarget) }));
}

function busyChip(b, on) {
  if (!b) return;
  b.setAttribute('aria-busy', String(on));
  if (on && !b.spinner) { b.spinner = el('span', { class: 'spin' }); b.swatch.replaceWith(b.spinner); }
  else if (!on && b.spinner) { b.spinner.replaceWith(b.swatch); b.spinner = null; }
}

function keysHelp(anchor) {
  const rows = [['drag', 'on the circle: choose a stretch'], ['click', 'on the circle: move there'], ['wheel', 'on the strip: zoom about the pointer'], ['drag', 'on the strip: pan'],
    ['+  −', 'zoom in, out'], ['←  →', 'pan a quarter (Shift: a whole window)'], ['0', 'the whole contig'], ['o', 'the minimum of the cumulative skew'], ['f', 'find in this genome'],
    ['Enter', 'in the search: the next hit (Shift: the previous)'], ['Esc', 'back to the genome’s card']];
  popover(anchor, el('div', { class: 'gv-keys' }, el('div', { class: 'menu-group' }, 'keys and gestures'),
    ...rows.map(([k, v]) => el('div', { class: 'gv-key' }, el('kbd', {}, k), el('span', {}, v)))), { width: 340 });
}

// Predicted genes, CpG islands or repeats of the contig, fetched once and drawn on the strip (and
// the repeats as chords on the circle).
async function annotate(kind, chipEl = null) {
  const find = () => chipEl?.isConnected ? chipEl : view.tools.querySelector(`.gv-chip[data-kind="${kind}"]`);
  if (!view.annotate[kind]) {
    delete view.extra[kind];
    const b = find();
    if (b) b.querySelector('.gv-count').textContent = '';
    drawCircle(); drawLinear();
    return;
  }
  const genome = view.genome, contig = view.contig;
  view.busy[kind] = true;
  busyChip(find(), true);
  view.ctx.status(kind === 'genes' ? 'predicting genes…' : kind === 'islands' ? 'finding CpG islands…' : 'finding repeats…');
  try {
    const r = await ask('genome.annotate', { genome, contig, what: kind, min_length: 100 });
    if (view.genome !== genome || view.contig !== contig || !view.annotate[kind]) return;
    view.extra[kind] = kind === 'genes' ? r.genes : kind === 'islands' ? r.islands : r.repeats;
    view.counts[kind] = kind === 'genes' ? r.genes.length : kind === 'islands' ? r.islands.length : r.count;
    const b = find();
    if (b) b.querySelector('.gv-count').textContent = fmt.count(view.counts[kind]);
    // What the numbers need said beside them: a model trained on too little, or the wrong background.
    const L = view.summary?.contigs?.[view.contig]?.length ?? Infinity;
    const mito = /mitochondri/i.test(`${view.summary?.title || ''} ${view.summary?.definition || ''}`);
    const small = L < 200000 ? ` (only ${fmt.length(L)} to train the coding chains on${mito ? '; mitochondrial genes also overlap and end in stops completed by polyadenylation' : ''})` : '';
    const bg = bacterial() ? ' (a bacterial genome is not depleted of CpG, and the model\'s background is the vertebrate one, so it flags much of the genome)' : '';
    view.ctx.status(kind === 'genes' ? `${fmt.count(r.genes.length)} genes predicted${r.sensitivity !== null ? `: ${fmt.percent(r.sensitivity)} of the ${fmt.count(r.annotated)} annotated CDSs found, ${fmt.percent(r.specificity)} of the predictions annotated` : ''}${small}`
      : kind === 'islands' ? `${fmt.count(r.islands.length)} CpG islands${bg}` : `${fmt.count(r.count)} repeats of at least 100 bp`);
    drawCircle();
    drawLinear();
  } catch (e) {
    toast(e.message, 'error');
    view.ctx.status('');
    view.annotate[kind] = false;
    find()?.set(false);
  } finally {
    view.busy[kind] = false;
    busyChip(find(), false);
  }
}

// ---- the circle ----------------------------------------------------------------------------------------
function drawCircle() {
  const s = view.summary, w = view.whole;
  if (!s || !w) return;
  const c = s.contigs[view.contig];
  const genes = (c.features.CDS || 0) + (c.features.rRNA || 0) + (c.features.tRNA || 0);
  const title = (s.organism || s.title || '').replace(/ str\. .*$/, '');
  view.circlePane.innerHTML = circleSvg({
    L: view.length, circular: c.circular, tracks: w, features: w.features, window: view.window, selected: view.selected,
    title: title.length > 34 ? `${title.slice(0, 32)}…` : title,
    subtitle: `${fmt.length(c.length)} · GC ${fmt.percent(c.gc)}${genes ? ` · ${fmt.count(genes)} genes` : ''}`,
    markers: bacterial() ? ['origin', 'terminus'] : ['skew min', 'skew max'],
    links: (view.extra.repeats || []).slice(0, 300).flatMap((r, k) => r.at.slice(1, 6).map(([q]) => [r.at[0][0], q, CATEGORIES[k % CATEGORIES.length], `a repeat of ${r.length} bp, ${r.at.length} copies`])),
  });
  view.circlePane.append(el('div', { class: 'circle-caption' }, `${c.circular ? 'circular, read clockwise from the top' : 'linear: its two ends meet at the top'} · drag to choose a stretch`));
  const key = (colour, text) => el('span', { class: 'key-item' }, el('span', { class: 'swatch', style: { background: colour } }), text);
  view.circlePane.append(el('div', { class: 'legend' },
    key(featureColour('CDS'), 'coding'), key(featureColour('rRNA'), 'rRNA'), key(featureColour('tRNA'), 'tRNA'), key(featureColour('ncRNA'), 'other RNA'), key(featureColour('mobile_element'), 'mobile'),
    el('span', { class: 'legend-gap' }),
    key(token('--high'), 'GC above its mean'), key(token('--low'), 'below'), key(token('--pos'), 'G leads (skew)'), key(token('--neg'), 'C leads'),
    key(token('--accent'), 'cumulative skew')));
}

function circlePoint(e) {
  const svg = view.circlePane.querySelector('svg');
  if (!svg) return null;
  const r = svg.getBoundingClientRect();
  const sx = 1230 / r.width, sy = 1096 / r.height;
  const x = (e.clientX - r.left) * sx - 615, y = (e.clientY - r.top) * sy - 548;
  return { x, y, r: Math.hypot(x, y) };
}

function onCircleDown(e) {
  if (!view.whole) return;
  const p = circlePoint(e);
  if (!p || p.r > 520) return;
  const target = e.target.closest?.('.feature');
  if (target) { selectAndShow(Number(target.dataset.i)); return; }
  const L = view.length;
  const start = positionAt(p.x, p.y, L);
  let last = start, moved = false;
  const wedge = () => {
    const [b, z] = last >= start ? [start, last] : [last, start];
    view.window = { begin: b, end: Math.max(b + 50, z) };
    drawCircle();
    updateBar();
  };
  const move = (ev) => { const q = circlePoint(ev); if (q) { last = positionAt(q.x, q.y, L); moved = true; wedge(); } };
  const up = () => {
    window.removeEventListener('pointermove', move);
    window.removeEventListener('pointerup', up);
    if (!moved || Math.abs(last - start) < L / 720) {
      // A click: the same width of window, centred here.
      const width = Math.min(L, view.window.end - view.window.begin === L ? Math.max(2000, Math.round(L / 40)) : view.window.end - view.window.begin);
      setWindow(start - width / 2, start + width / 2);
    } else {
      setWindow(view.window.begin, view.window.end);
    }
  };
  window.addEventListener('pointermove', move);
  window.addEventListener('pointerup', up);
}

function onCircleHover(e) {
  const p = circlePoint(e);
  if (!p || !view.whole) return;
  const f = e.target.closest?.('.feature');
  const pos = positionAt(p.x, p.y, view.length);
  if (f) {
    const feat = view.whole.features.find((x) => x.i === Number(f.dataset.i));
    if (feat) { hover(featureTip(feat), e); return; }
  }
  if (p.r < 480 && p.r > 140) hover(`<b>${fmt.position(pos)}</b>${p.r < 290 ? '<br>the cumulative GC skew' : p.r < 350 ? '<br>GC skew' : p.r < 425 ? '<br>GC content against its mean' : ''}`, e);
  else hover(null);
}

// ---- the strip ----------------------------------------------------------------------------------------
function drawLinear() {
  const t = view.tracks;
  if (!t || !view.summary) return;
  const width = view.stage.clientWidth;
  const layout = linearSvg({ width, tracks: t, selected: view.selected, hits: view.hits, contig: view.contig, show: view.show, extra: view.extra });
  view.layout = layout;
  view.stage.innerHTML = layout.svg;
  updateBar();
}

// The strip's bar, built once: where the window is, its navigation, and the search.
function buildBar() {
  const zoom = (f) => { if (!view.window) return; const { begin, end } = view.window; const m = (begin + end) / 2, h = ((end - begin) * f) / 2; setWindow(m - h, m + h); };
  const pan = (f) => { if (!view.window) return; const { begin, end } = view.window; const d = (end - begin) * f; setWindow(begin + d, end + d); };
  view.zoom = zoom;
  view.pan = pan;
  view.whereText = el('div', { class: 'gv-where-text' });
  view.locator = el('span');
  const track = el('div', { class: 'gv-locator', title: 'the window within the contig: click to move there', onpointerdown: (e) => {
    if (!view.window) return;
    const r = track.getBoundingClientRect(), at = ((e.clientX - r.left) / r.width) * view.length, w = view.window.end - view.window.begin;
    setWindow(at - w / 2, at + w / 2);
  } }, view.locator);
  const where = el('div', { class: 'gv-where' }, view.whereText, track);
  view.originBtn = glyphBtn('target', 'the origin: the minimum of the cumulative GC skew (o)', () => goOrigin());
  const nav = el('div', { class: 'gv-cluster', role: 'group', 'aria-label': 'move along the genome' },
    glyphBtn('left', 'pan left (←)', () => pan(-0.25)),
    glyphBtn('zoomOut', 'zoom out (−)', () => zoom(2)),
    glyphBtn('zoomIn', 'zoom in (+)', () => zoom(0.5)),
    glyphBtn('right', 'pan right (→)', () => pan(0.25)),
    el('span', { class: 'gv-cluster-sep' }),
    glyphBtn('fit', 'the whole contig (0)', () => setWindow(0, view.length)),
    view.originBtn);
  // The search: a gene, a sequence, a position or a range; Enter again steps through the hits.
  const input = el('input', { type: 'search', placeholder: 'Find a gene, ACGT…, or 1200-4800', spellcheck: false, autocomplete: 'off', 'aria-label': 'find in this genome' });
  input.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      const q = input.value.trim();
      if (q && q === view.query && view.hits.length) stepHit(e.shiftKey ? -1 : 1);
      else search(q);
    } else if (e.key === 'Escape') {
      e.stopPropagation();
      if (input.value) { input.value = ''; syncClear(); search(''); } else input.blur();
    }
  });
  input.addEventListener('input', () => { syncClear(); if (!input.value) search(''); });
  view.searchInput = input;
  view.searchIcon = el('span', { class: 'gv-search-icon' }, icon('search', 'sm'));
  view.clearBtn = iconBtn('x', { title: 'clear the search (Esc)', size: 'sm', onClick: () => { input.value = ''; syncClear(); search(''); input.focus(); } });
  view.clearBtn.hidden = true;
  view.hitText = el('span', { class: 'gv-hits-text' });
  view.hitBox = el('span', { class: 'gv-hits', hidden: true }, view.hitText,
    glyphBtn('up', 'the previous hit (Shift Enter)', () => stepHit(-1)), glyphBtn('down', 'the next hit (Enter)', () => stepHit(1)));
  const field = el('label', { class: 'gv-search' }, view.searchIcon, input, view.clearBtn, view.hitBox);
  clear(view.bar).append(where, nav, field);
}

function syncClear() { if (view.clearBtn) view.clearBtn.hidden = !view.searchInput.value; }

function updateBar() {
  if (!view.window || !view.whereText) return;
  const { begin, end } = view.window, L = view.length || 1;
  const whole = end - begin >= L;
  view.whereText.replaceChildren(el('b', {}, fmt.interval(begin, end)),
    el('span', {}, `${fmt.length(end - begin)} · ${whole ? 'the whole' : `${fmt.number((100 * (end - begin)) / L, 2)}% of it`}`));
  Object.assign(view.locator.style, { left: `${(100 * begin) / L}%`, width: `${(100 * (end - begin)) / L}%` });
  view.originBtn.hidden = !(view.whole && bacterial());
  view.hitBox.hidden = !view.hits.length;
  if (view.hits.length) view.hitText.textContent = `${view.hitAt >= 0 ? fmt.count(view.hitAt + 1) : '–'} of ${fmt.count(view.hitTotal)}${view.hitTotal > view.hits.length ? ` (first ${fmt.count(view.hits.length)})` : ''}`;
}

function setWindow(b, e, push = false) {
  const L = view.length;
  if (!L) return;
  let w = Math.max(40, Math.min(L, e - b));
  b = Math.round(Math.max(0, Math.min(L - w, b)));
  w = Math.round(w);
  view.window = { begin: b, end: b + w };
  view.ctx.setParams({ begin: b, end: b + w }, !push);
  drawCircle();
  updateBar();
  fetchWindow();
}

const goOrigin = () => { if (view.whole) setWindow(view.whole.origin - 5000, view.whole.origin + 5000); };

function onWheel(e) {
  if (!view.layout) return;
  e.preventDefault();
  const r = view.stage.getBoundingClientRect();
  const at = view.layout.position(e.clientX - r.left);
  const { begin, end } = view.window;
  const f = e.deltaY > 0 ? 1.3 : 1 / 1.3;
  setWindow(at - (at - begin) * f, at + (end - at) * f);
}

function onStripDown(e) {
  if (!view.layout) return;
  const target = e.target.closest?.('.feature');
  if (target) { selectAndShow(Number(target.dataset.i)); return; }
  const svg = view.stage.querySelector('svg');
  const x0 = e.clientX, { begin, end } = view.window;
  const perPx = (end - begin) / Math.max(1, view.stage.clientWidth - LEFT - 14);
  svg?.classList.add('dragging');
  const move = (ev) => {
    const d = (ev.clientX - x0) * perPx;
    view.window = { begin: Math.round(begin - d), end: Math.round(end - d) };
    // Slide the drawing until the new tracks come.
    if (svg) svg.style.transform = `translateX(${ev.clientX - x0}px)`;
  };
  const up = (ev) => {
    window.removeEventListener('pointermove', move);
    window.removeEventListener('pointerup', up);
    svg?.classList.remove('dragging');
    const d = (ev.clientX - x0) * perPx;
    if (Math.abs(ev.clientX - x0) > 2) setWindow(begin - d, end - d);
  };
  window.addEventListener('pointermove', move);
  window.addEventListener('pointerup', up);
}

function onStripHover(e) {
  if (!view.layout || !view.tracks) return;
  const f = e.target.closest?.('.feature');
  if (f) {
    const feat = view.tracks.features.find((x) => x.i === Number(f.dataset.i));
    if (feat) { hover(featureTip(feat), e); return; }
  }
  const r = view.stage.getBoundingClientRect();
  const px = e.clientX - r.left;
  if (px < LEFT) { hover(null); return; }
  const p = Math.floor(view.layout.position(px));
  const t = view.tracks;
  const k = Math.floor(((p - t.begin) / (t.end - t.begin)) * t.gc.length);
  const base = t.letters ? t.letters[p - t.begin] : null;
  hover(`<b>${fmt.position(p)}</b>${base ? ` · ${base}` : ''}<br>GC ${fmt.percent(t.gc[k])} · skew ${fmt.number(t.skew[k], 2)} · CpG o/e ${fmt.number(t.cpg[k], 2)}`, e);
}

// The view's keys, while it is shown and nothing is being typed.
function onKey(e) {
  if (!view || view.section.hidden || !view.window || e.ctrlKey || e.metaKey || e.altKey) return;
  const a = document.activeElement;
  if (/^(INPUT|TEXTAREA|SELECT)$/.test(a?.tagName || '') || a?.isContentEditable) return;
  if (document.querySelector('.dialog-veil, .popover')) return;
  const k = e.key;
  if (k === '+' || k === '=') view.zoom(0.5);
  else if (k === '-' || k === '_') view.zoom(2);
  else if (k === 'ArrowLeft') view.pan(e.shiftKey ? -1 : -0.25);
  else if (k === 'ArrowRight') view.pan(e.shiftKey ? 1 : 0.25);
  else if (k === '0' || k === 'Home') setWindow(0, view.length);
  else if (k === 'o') goOrigin();
  else if (k === 'f') view.searchInput.focus();
  else if (k === 'Escape' && view.selected !== null) { view.ctx.setParams({ feature: null }); select(null, true); }
  else return;
  e.preventDefault();
}

// Whether the genome is a bacterial or archaeal chromosome, whose skew marks its origin and terminus.
const bacterial = () => (view.summary?.lineage || []).some((x) => x === 'Bacteria' || x === 'Archaea');

const featureTip = (f) => `<b>${esc(f.gene || f.name)}</b> <span class="faint">${esc(f.key)}</span>${f.product ? `<br>${esc(f.product)}` : ''}<br>${f.parts.map(([b, e]) => fmt.interval(b, e)).join(', ')} (${f.strand})`;

// ---- searching ----------------------------------------------------------------------------------------
function stepHit(d) {
  const n = view.hits.length;
  if (!n) return;
  view.hitAt = view.hitAt < 0 ? (d > 0 ? 0 : n - 1) : (view.hitAt + d + n) % n;
  const [c, p] = view.hits[view.hitAt];
  const b = p - 60, e = p + 60 + (view.query || '').length;
  if (c !== view.contig) {
    view.keepHits = true;
    view.ctx.setParams({ contig: c, begin: Math.max(0, Math.round(b)), end: Math.round(e), feature: null });
    show(view.ctx);
    return;
  }
  setWindow(b, e);
}

function searching(on) {
  view.searchIcon.replaceChildren(on ? el('span', { class: 'spin' }) : icon('search', 'sm'));
}

async function search(q) {
  view.query = q;
  if (!q) { view.hits = []; view.hitTotal = 0; view.hitAt = -1; updateBar(); drawLinear(); return; }
  if (!view.whole) return;
  // A range, 1-based and inclusive as written (1200-4800, 1,200..4,800); or one position.
  const range = q.replace(/,/g, '').match(/^(\d+)\s*(?:-|–|\.\.|:)\s*(\d+)$/);
  if (range) {
    const a = Number(range[1]), z = Number(range[2]);
    setWindow(Math.min(a, z) - 1, Math.max(a, z));
    return;
  }
  if (/^\d[\d,]*$/.test(q)) {
    const p = Number(q.replace(/,/g, '')) - 1, w = view.window.end - view.window.begin;
    const width = w >= view.length ? Math.max(2000, Math.round(view.length / 40)) : w;
    setWindow(p - width / 2, p + width / 2);
    return;
  }
  if (/^[ACGTacgt]+$/.test(q) && q.length >= 4) {
    searching(true);
    try {
      const r = await api.run('genome.search', { genome: view.genome, pattern: q });
      if (view.query !== q) return;
      view.hits = r.hits;
      view.hitTotal = r.total;
      view.hitAt = -1;
      toast(`${fmt.count(r.total)} occurrence${r.total === 1 ? '' : 's'} of ${q.length} bases on both strands, found in ${fmt.seconds(r._seconds)}${r.hits.length ? '; Enter steps through them' : ''}`);
      updateBar();
      drawLinear();
    } catch (e) { toast(e.message, 'error'); } finally { searching(false); }
    return;
  }
  // A gene or a locus tag: the features of the whole genome.
  view.hits = []; view.hitTotal = 0; view.hitAt = -1;
  const lower = q.toLowerCase();
  const hit = view.whole.features.find((f) => [f.gene, f.locus_tag, f.name].some((v) => v && v.toLowerCase() === lower))
    || view.whole.features.find((f) => [f.gene, f.locus_tag, f.name, f.product].some((v) => v && v.toLowerCase().includes(lower)));
  if (!hit) { updateBar(); toast(`No feature named ${q}`, 'error'); return; }
  const lo = Math.min(...hit.parts.map((p) => p[0])), hi = Math.max(...hit.parts.map((p) => p[1]));
  setWindow(lo - (hi - lo) * 0.6 - 200, hi + (hi - lo) * 0.6 + 200);
  selectAndShow(hit.i);
}

// ---- the cards ---------------------------------------------------------------------------------------
function selectAndShow(i) {
  view.ctx.setParams({ feature: i });
  select(i, true);
}

function select(i, redraw) {
  view.selected = i;
  if (redraw) { drawCircle(); drawLinear(); }
  if (i === null) { showGenomeCard(); return; }
  showFeatureCard(i);
}

const body = (...children) => el('div', { class: 'card-body' }, ...children);

function showGenomeCard() {
  const s = view.summary, panel = clear(view.panel);
  if (!s) return;
  const c = s.contigs[view.contig] || s.contigs[0];
  const info = view.ctx.genomeInfo?.();
  const genes = (c.features.CDS || 0) + (c.features.rRNA || 0) + (c.features.tRNA || 0);
  panel.append(card({ cls: 'gv-card', body: body(
    el('p', { class: 'eyebrow' }, icon(s.source === 'simulated' ? 'spark' : info?.group === 'files' ? 'file' : 'genome', 'sm'),
      s.source === 'simulated' ? 'a simulated genome' : info?.group === 'files' ? 'your genome' : 'genome'),
    el('h3', { class: 'gv-title' }, s.organism || s.title),
    el('div', { class: 'gv-sub' }, c.definition || s.title),
    metrics([[fmt.length(c.length), c.circular ? 'circular' : 'linear'], [fmt.percent(c.gc), 'GC'], genes ? [fmt.count(genes), 'genes'] : null]),
    kv([['accession', c.accession || null], ['length', `${fmt.count(c.length)} bp`], ['GC', fmt.percent(c.gc, 2)], ['genetic code', `${s.table}, ${s.table_name}`],
      ['repeats', s.repeats?.length ? `${s.repeats.length} copies (the simulation's truth)` : null]])) }));
  // The features of the record, by their key.
  const counts = Object.entries(c.features).sort((a, b) => b[1] - a[1]).slice(0, 7);
  if (counts.length) {
    const max = counts[0][1];
    panel.append(card({ title: 'Features', sub: 'by their key in the record', pad: true, body: el('div', { class: 'gv-bars' },
      ...counts.map(([k, n]) => el('div', { class: 'gv-bar-row' },
        el('span', { class: 'k' }, el('span', { class: 'swatch', style: { background: featureColour(k) } }), k),
        el('span', { class: 'bar' }, el('span', { style: { width: `${Math.max(1.5, (100 * n) / max)}%`, background: featureColour(k) } })),
        el('span', { class: 'n' }, fmt.count(n))))) }));
  }
  if (view.whole) {
    const lineage = s.lineage || [];
    const extremes = `The cumulative GC skew falls to its minimum at ${fmt.position(view.whole.origin)} and rises to its maximum at ${fmt.position(view.whole.terminus)}.`;
    const reading = lineage.includes('Bacteria') || lineage.includes('Archaea')
      ? ' In most bacteria these mark the origin and the terminus of replication: the leading strand is richer in G than in C, so the skew changes sign at both (Lobry 1996; Grigoriev 1998).'
      : /mitochondri/i.test(`${s.title} ${c.definition}`)
        ? ' A mitochondrion is copied by strand displacement, which leaves one strand single for long and wears its cytosines into thymines: the skew keeps one sign along nearly the whole genome, and its extremes lie near the control region rather than at a bacterial origin and terminus (Tanaka and Ozawa 1994; Reyes et al. 1998).'
        : s.source === 'simulated'
          ? ' This genome was drawn without strand asymmetry, so its skew is noise and its extremes mean nothing: a control for reading the others.'
          : ' Outside bacterial chromosomes the extremes need not mark replication; read them against what is known of how this genome is copied.';
    const to = (p) => () => setWindow(p - 5000, p + 5000);
    panel.append(card({ title: 'Replication', sub: 'read from the cumulative GC skew', body: body(
      el('p', { class: 'note gv-prose' }, extremes + reading),
      el('div', { class: 'gv-acts' },
        withGlyph(btn(bacterial() ? 'The origin' : 'The minimum', { size: 'sm', onClick: to(view.whole.origin) }), 'target'),
        withGlyph(btn(bacterial() ? 'The terminus' : 'The maximum', { size: 'sm', onClick: to(view.whole.terminus) }), 'target'))) }));
  }
  if (s.lineage?.length) {
    const path = el('div', { class: 'gv-lineage' });
    s.lineage.forEach((x, k) => { if (k) path.append(glyph('right')); path.append(el('span', {}, x)); });
    panel.append(card({ title: 'Lineage', pad: true, body: path }));
  }
  if (s.problems?.length) panel.append(card({ title: 'Problems in the file', pad: true, body: el('div', { class: 'stack' }, ...s.problems.slice(0, 8).map((p) => callout(`line ${p.line}: ${p.message}`, 'bad'))) }));
  panel.append(el('p', { class: 'small faint gv-own' }, 'Your own genome: drop a FASTA or GenBank file anywhere on the page, or ',
    el('button', { class: 'link', onclick: () => view.ctx.go('data') }, 'add it in Your data'), '; it opens here.'));
}

async function showFeatureCard(i) {
  const panel = clear(view.panel);
  panel.append(cardSkeleton());
  let f;
  try { f = await ask('genome.feature', { genome: view.genome, contig: view.contig, index: i }); } catch (e) {
    if (view.selected === i) clear(panel).append(callout(e.message, 'bad'));
    return;
  }
  if (view.selected !== i) return;
  clear(panel);
  const back = withGlyph(btn('The genome', { size: 'sm', kind: 'ghost', title: 'back to the genome’s card (Esc)', onClick: () => { view.ctx.setParams({ feature: null }); select(null, true); } }), 'left');
  panel.append(el('div', { class: 'gv-back' }, back, el('span', { class: 'small faint' }, `feature ${fmt.count(i + 1)}`)));
  const lo = Math.min(...f.parts.map((p) => p[0])), hi = Math.max(...f.parts.map((p) => p[1]));
  const aa = f.key === 'CDS' ? Math.floor(f.length / 3) - 1 : null;
  panel.append(card({ cls: 'gv-card', body: body(
    el('p', { class: 'eyebrow' }, el('span', { class: 'swatch', style: { background: featureColour(f.key) } }), f.key,
      el('span', { class: 'badge' }, f.strand === '-' ? 'reverse (−)' : 'forward (+)')),
    el('h3', { class: 'gv-title' }, f.gene || f.name),
    f.product ? el('div', { class: 'gv-sub' }, f.product) : null,
    metrics([[fmt.count(f.length), 'bp'], aa !== null ? [fmt.count(aa), 'amino acids'] : null, [fmt.percent(f.gc), 'GC']]),
    kv([['locus tag', f.locus_tag || null], ['location', f.location.replace(/([(,])/g, '$1\u200b')], ['strand', f.strand === '-' ? 'reverse (−)' : 'forward (+)'],
      ['GC at third positions', f.gc3 !== undefined ? fmt.percent(f.gc3) : null], ['effective codons (Wright)', f.nc !== undefined ? fmt.number(f.nc, 3) : null]]),
    el('div', { class: 'gv-acts' },
      withGlyph(btn('Show it', { size: 'sm', title: 'the feature on the strip, with some of its neighbourhood', onClick: () => setWindow(lo - (hi - lo) * 0.4 - 100, hi + (hi - lo) * 0.4 + 100) }), 'target'),
      f.translation ? btn('Align its protein', { size: 'sm', icon: 'align', title: 'with another sequence, in the Toolbox', onClick: () => view.ctx.openTool('align-pair', { a: f.translation, alphabet: 'protein' }) }) : null)) }));
  if (f.translation) {
    panel.append(card({ title: ['Protein', f.agrees === true ? el('span', { class: 'badge good' }, 'as RefSeq translates it') : f.agrees === false ? el('span', { class: 'badge bad' }, 'differs from RefSeq’s') : null],
      actions: [iconBtn('copy', { title: 'copy the protein', size: 'sm', onClick: () => copy(f.translation) })],
      body: body(el('div', { class: 'protein' }, f.translation)) }));
  }
  if (f.codons) {
    const max = Math.max(1, ...f.codons.map((c) => c[2]));
    const grid = el('div', { class: 'codon-grid' });
    // NCBI's order: the first base by quarter, the second by column, the third by row.
    for (let first = 0; first < 4; first++)
      for (let third = 0; third < 4; third++)
        for (let second = 0; second < 4; second++) {
          const [codon, aa_, n] = f.codons[16 * first + 4 * second + third];
          grid.append(el('div', { style: { background: warm(0.12 + 0.8 * (n / max)), color: n / max > 0.6 ? '#fff' : '' }, title: `${codon} ${aa_}: ${n}` }, el('span', {}, `${codon} ${aa_}`), el('span', {}, String(n))));
        }
    panel.append(card({ title: 'Codons used', sub: 'darker, used more', body: body(grid) }));
  }
  if (f.qualifiers?.length) {
    const list = el('dl', { class: 'kv small' });
    for (const [k, v] of f.qualifiers) if (k !== 'translation') list.append(el('dt', {}, k), el('dd', { style: { wordBreak: 'break-word' } }, v || '✓'));
    panel.append(card({ title: 'Qualifiers', sub: 'as the record gives them', body: body(list) }));
  }
}

export { angleOf };
