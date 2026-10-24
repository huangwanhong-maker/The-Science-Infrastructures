// The Genome view: a genome as a circle and as a strip. The circle shows the whole: features by
// strand, GC content and skew, and the cumulative skew whose extremes mark the origin and the
// terminus of replication. The strip shows a window: features in rows, the six frames' ORFs,
// composition in bins, the letters at base resolution, the hits of a search. A card shows the
// genome, or the feature chosen, with its protein checked against the record's own and its
// codons.
//
// Drag on the circle to choose a window, or click it to move there; on the strip, the wheel
// zooms about the pointer and dragging pans. The window and the chosen feature live in the
// address (begin, end, contig, feature).

import { el, clear, esc, hover, notice, toggles, field } from '../dom.js';
import * as api from '../api.js';
import { hashNumber } from '../state.js';
import { circleSvg, positionAt, angleOf } from '../draw/circle.js';
import { linearSvg, LEFT } from '../draw/linear.js';
import * as fmt from '../format.js';
import { featureColour, warm, token, CATEGORIES } from '../palette.js';

const ask = api.cached(64);
let view = null;   // the view's elements and state

export function mount(section, ctx) {
  const tools = el('div', { class: 'tools' });
  const body = el('div', { class: 'genome-body' });
  const circlePane = el('div', { class: 'circle-pane' });
  const linearPane = el('div', { class: 'linear-pane' });
  const card = el('aside', { class: 'card-panel', 'aria-live': 'polite' });
  const bar = el('div', { class: 'linear-bar' });
  const stage = el('div', { class: 'linear-stage' });
  linearPane.append(bar, stage);
  body.append(circlePane, linearPane, card);
  section.append(tools, body);
  view = { ctx, tools, circlePane, linearPane, bar, stage, card, genome: null, summary: null, whole: null, tracks: null, contig: 0,
           window: null, selected: null, hits: [], hitAt: -1, show: { orfs: true, gc: true, skew: true, cpg: true, entropy: false }, request: 0, layout: null,
           annotate: { genes: false, islands: false, repeats: false }, extra: {} };
  new ResizeObserver(() => drawLinear()).observe(stage);
  stage.addEventListener('wheel', onWheel, { passive: false });
  stage.addEventListener('pointerdown', onStripDown);
  stage.addEventListener('pointermove', onStripHover);
  stage.addEventListener('pointerleave', () => hover(null));
  circlePane.addEventListener('pointerdown', onCircleDown);
  circlePane.addEventListener('pointermove', onCircleHover);
  circlePane.addEventListener('pointerleave', () => hover(null));
  buildTools();
}

export async function show(ctx) {
  view.ctx = ctx;
  if (!ctx.genome) { clear(view.card).append(el('p', { class: 'wait' }, 'No genome is fetched yet: run Setup.')); return; }
  const contig = hashNumber(ctx.params, 'contig', 0, 0);
  if (view.genome !== ctx.genome || view.contig !== contig) {
    view.genome = ctx.genome;
    view.contig = contig;
    view.hits = [];
    view.selected = null;
    view.extra = {};
    for (const k of Object.keys(view.annotate)) view.annotate[k] = false;
    await load();
  }
  const L = view.length;
  const begin = hashNumber(ctx.params, 'begin', 0, 0, L - 1), end = hashNumber(ctx.params, 'end', L, begin + 1, L);
  view.window = { begin: Math.floor(begin), end: Math.ceil(end) };
  const feature = ctx.params.feature !== undefined ? Number(ctx.params.feature) : null;
  if (feature !== view.selected) select(Number.isFinite(feature) ? feature : null, false);
  drawCircle();
  fetchWindow();
}

// ---- loading ----------------------------------------------------------------------------------------
async function load() {
  const id = view.genome;
  view.ctx.status('reading the genome…');
  clear(view.circlePane).append(el('p', { class: 'wait' }, 'Reading the genome…'));
  try {
    const summary = await ask('genome.summary', { genome: id });
    if (view.genome !== id) return;
    view.summary = summary;
    const c = summary.contigs[view.contig] || summary.contigs[0];
    view.length = c.length;
    view.whole = await ask('genome.tracks', { genome: id, contig: view.contig, bins: 720, orfs: false, max_features: 20000 });
    view.ctx.status(`read in ${fmt.seconds(summary._seconds)}`);
    buildTools();
    if (view.selected === null) showGenomeCard();
  } catch (e) {
    view.ctx.status('');
    clear(view.circlePane).append(el('p', { class: 'error' }, e.message));
  }
}

let windowTimer = null;
function fetchWindow() {
  clearTimeout(windowTimer);
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
      notice(e.message, true);
    }
  }, 90);
}

// ---- tools ----------------------------------------------------------------------------------------
function buildTools() {
  const t = clear(view.tools);
  const s = view.summary;
  if (s && s.contigs.length > 1) {
    const sel = el('select', { onchange: () => view.ctx.setParams({ contig: sel.value, begin: null, end: null, feature: null }) || show(view.ctx) });
    s.contigs.forEach((c, k) => sel.append(el('option', { value: k, selected: k === view.contig }, `${c.name} (${fmt.length(c.length)})`)));
    t.append(field('contig', sel));
  }
  const tracks = [['orfs', 'ORFs'], ['gc', 'GC'], ['skew', 'skew'], ['cpg', 'CpG'], ['entropy', 'entropy']];
  const box = el('div', { class: 'row' }, el('span', { class: 'faint small' }, 'tracks'));
  for (const [k, label] of tracks) {
    box.append(el('label', { class: 'field' }, el('input', { type: 'checkbox', checked: view.show[k], onchange: (e) => { view.show[k] = e.target.checked; drawLinear(); } }), label));
  }
  t.append(box);
  const ann = el('div', { class: 'row' }, el('span', { class: 'faint small' }, 'from the sequence alone'));
  for (const [k, label, title] of [['genes', 'predicted genes', 'a self-trained gene finder: periodic Markov chains of coding'], ['islands', 'CpG islands', 'the HMM of Durbin et al.'], ['repeats', 'repeats ≥ 100 bp', 'maximal repeats on both strands, by the suffix array']]) {
    ann.append(el('label', { class: 'field', title }, el('input', { type: 'checkbox', checked: view.annotate[k], onchange: (e) => { view.annotate[k] = e.target.checked; annotate(k); } }), label));
  }
  t.append(ann);
  const orf = el('select', { onchange: () => { view.ctx.setParams({ orf: orf.value }); fetchWindow(); } });
  for (const v of [30, 50, 100, 200, 300]) orf.append(el('option', { value: v, selected: Number(view.ctx.params.orf || 100) === v }, `${v} codons`));
  t.append(field('ORFs from', orf));
  t.append(el('span', { class: 'hint' }, 'drag on the circle to choose a stretch · the wheel zooms the strip · drag to pan'));
}

// Predicted genes, CpG islands or repeats of the contig, fetched once and drawn on the strip (and
// the repeats as chords on the circle).
async function annotate(kind) {
  if (!view.annotate[kind]) { delete view.extra[kind]; drawCircle(); drawLinear(); return; }
  view.ctx.status(kind === 'genes' ? 'predicting genes…' : kind === 'islands' ? 'finding CpG islands…' : 'finding repeats…');
  try {
    const r = await ask('genome.annotate', { genome: view.genome, contig: view.contig, what: kind, min_length: 100 });
    view.extra[kind] = kind === 'genes' ? r.genes : kind === 'islands' ? r.islands : r.repeats;
    // What the numbers need said beside them: a model trained on too little, or the wrong background.
    const L = view.summary?.contigs?.[view.contig]?.length ?? Infinity;
    const small = L < 200000 ? ` (only ${fmt.length(L)} to train the coding chains on; mitochondrial genes also overlap and end in stops completed by polyadenylation)` : '';
    const bg = bacterial() ? ' (a bacterial genome is not depleted of CpG, and the model\'s background is the vertebrate one, so it flags much of the genome)' : '';
    view.ctx.status(kind === 'genes' ? `${fmt.count(r.genes.length)} genes predicted${r.sensitivity !== null ? `: ${fmt.percent(r.sensitivity)} of the ${fmt.count(r.annotated)} annotated CDSs found, ${fmt.percent(r.specificity)} of the predictions annotated` : ''}${small}`
      : kind === 'islands' ? `${fmt.count(r.islands.length)} CpG islands${bg}` : `${fmt.count(r.count)} repeats of at least 100 bp`);
    drawCircle();
    drawLinear();
  } catch (e) { notice(e.message, true); view.ctx.status(''); }
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
  view.circlePane.append(el('div', { class: 'circle-caption' }, c.circular ? 'circular, read clockwise from the top' : 'linear: its two ends meet at the top'));
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
  drawBar();
}

function drawBar() {
  const b = clear(view.bar), { begin, end } = view.window;
  const L = view.length;
  b.append(el('span', { class: 'where' }, `${fmt.interval(begin, end)} · ${fmt.length(end - begin)}`));
  const zoom = (f) => { const m = (begin + end) / 2, h = ((end - begin) * f) / 2; setWindow(m - h, m + h); };
  b.append(el('button', { class: 'button small', onclick: () => zoom(0.5), title: 'zoom in' }, '+'));
  b.append(el('button', { class: 'button small', onclick: () => zoom(2), title: 'zoom out' }, '−'));
  b.append(el('button', { class: 'button small', onclick: () => setWindow(0, L) }, 'whole'));
  if (view.whole && bacterial()) b.append(el('button', { class: 'button small', title: 'the minimum of the cumulative GC skew', onclick: () => setWindow(view.whole.origin - 5000, view.whole.origin + 5000) }, 'origin'));
  const input = el('input', { class: 'search', placeholder: 'search: ACGT… or a gene', value: view.query || '', onkeydown: (e) => { if (e.key === 'Enter') search(input.value.trim()); } });
  b.append(input);
  if (view.hits.length) {
    const step = (d) => {
      view.hitAt = (view.hitAt + d + view.hits.length) % view.hits.length;
      const [c, p] = view.hits[view.hitAt];
      if (c !== view.contig) view.ctx.setParams({ contig: c });
      setWindow(p - 60, p + 60 + (view.query || '').length);
    };
    b.append(el('span', { class: 'faint small' }, `${view.hitAt + 1 || '–'} of ${fmt.count(view.hitTotal)}${view.hitTotal > view.hits.length ? ` (first ${view.hits.length})` : ''}`));
    b.append(el('button', { class: 'button small', onclick: () => step(-1) }, '‹'), el('button', { class: 'button small', onclick: () => step(1) }, '›'));
  }
}

function setWindow(b, e, push = false) {
  const L = view.length;
  let w = Math.max(40, Math.min(L, e - b));
  b = Math.round(Math.max(0, Math.min(L - w, b)));
  w = Math.round(w);
  view.window = { begin: b, end: b + w };
  view.ctx.setParams({ begin: b, end: b + w }, !push);
  drawCircle();
  fetchWindow();
}

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

// Whether the genome is a bacterial or archaeal chromosome, whose skew marks its origin and terminus.
const bacterial = () => (view.summary?.lineage || []).some((x) => x === 'Bacteria' || x === 'Archaea');

const featureTip = (f) => `<b>${esc(f.gene || f.name)}</b> <span class="faint">${esc(f.key)}</span>${f.product ? `<br>${esc(f.product)}` : ''}<br>${f.parts.map(([b, e]) => fmt.interval(b, e)).join(', ')} (${f.strand})`;

// ---- searching ----------------------------------------------------------------------------------------
async function search(q) {
  view.query = q;
  if (!q) { view.hits = []; drawLinear(); return; }
  if (/^[ACGTacgt]+$/.test(q) && q.length >= 4) {
    try {
      const r = await api.run('genome.search', { genome: view.genome, pattern: q });
      view.hits = r.hits;
      view.hitTotal = r.total;
      view.hitAt = -1;
      notice(`${fmt.count(r.total)} occurrence${r.total === 1 ? '' : 's'} of ${q.length} bases on both strands, found in ${fmt.seconds(r._seconds)}`);
      drawLinear();
    } catch (e) { notice(e.message, true); }
    return;
  }
  // A gene or a locus tag: the features of the whole genome.
  const lower = q.toLowerCase();
  const hit = view.whole.features.find((f) => [f.gene, f.locus_tag, f.name].some((v) => v && v.toLowerCase() === lower))
    || view.whole.features.find((f) => [f.gene, f.locus_tag, f.name, f.product].some((v) => v && v.toLowerCase().includes(lower)));
  if (!hit) { notice(`No feature named ${q}`, true); return; }
  const lo = Math.min(...hit.parts.map((p) => p[0])), hi = Math.max(...hit.parts.map((p) => p[1]));
  setWindow(lo - (hi - lo) * 0.6 - 200, hi + (hi - lo) * 0.6 + 200);
  selectAndShow(hit.i);
}

// ---- the card ----------------------------------------------------------------------------------------
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

function showGenomeCard() {
  const s = view.summary, card = clear(view.card);
  if (!s) return;
  const c = s.contigs[view.contig];
  card.append(el('p', { class: 'eyebrow' }, s.source === 'simulated' ? 'a simulated genome' : 'genome'));
  card.append(el('h3', {}, s.organism || s.title));
  card.append(el('div', { class: 'sub' }, c.definition || s.title));
  const kv = el('dl', { class: 'kv' });
  const row = (k, v) => kv.append(el('dt', {}, k), el('dd', {}, v));
  if (c.accession) row('accession', c.accession);
  row('length', `${fmt.count(c.length)} bp, ${c.circular ? 'circular' : 'linear'}`);
  row('GC', fmt.percent(c.gc, 2));
  row('genetic code', `${s.table}, ${s.table_name}`);
  for (const [k, n] of Object.entries(c.features).sort((a, b) => b[1] - a[1]).slice(0, 7)) row(k, fmt.count(n));
  if (s.repeats?.length) row('repeats', `${s.repeats.length} copies (the simulation's truth)`);
  card.append(kv);
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
    card.append(el('div', { class: 'panel' },
      el('p', { class: 'eyebrow' }, 'replication'),
      el('p', { class: 'note' }, extremes + reading),
      el('div', { class: 'row' }, el('button', { class: 'button small', onclick: () => setWindow(view.whole.origin - 5000, view.whole.origin + 5000) }, bacterial() ? 'the origin' : 'the minimum'),
        el('button', { class: 'button small', onclick: () => setWindow(view.whole.terminus - 5000, view.whole.terminus + 5000) }, bacterial() ? 'the terminus' : 'the maximum'))));
  }
  if (s.lineage?.length) card.append(el('div', { class: 'panel' }, el('p', { class: 'eyebrow' }, 'lineage'), el('p', { class: 'small faint' }, s.lineage.join(' › '))));
  if (s.problems?.length) card.append(el('div', { class: 'panel' }, el('p', { class: 'eyebrow' }, 'problems in the file'), ...s.problems.slice(0, 8).map((p) => el('p', { class: 'small error' }, `line ${p.line}: ${p.message}`))));
}

async function showFeatureCard(i) {
  const card = clear(view.card);
  card.append(el('p', { class: 'wait' }, '…'));
  let f;
  try { f = await ask('genome.feature', { genome: view.genome, contig: view.contig, index: i }); } catch (e) { clear(card).append(el('p', { class: 'error' }, e.message)); return; }
  if (view.selected !== i) return;
  clear(card);
  card.append(el('div', { class: 'row', style: { justifyContent: 'space-between' } },
    el('p', { class: 'eyebrow' }, el('span', { class: 'swatch', style: { background: featureColour(f.key) } }), f.key),
    el('button', { class: 'link small', onclick: () => { view.ctx.setParams({ feature: null }); select(null, true); } }, 'the genome')));
  card.append(el('h3', {}, f.gene || f.name));
  if (f.product) card.append(el('div', { class: 'sub' }, f.product));
  const kv = el('dl', { class: 'kv' });
  const row = (k, v) => kv.append(el('dt', {}, k), el('dd', {}, v));
  if (f.locus_tag) row('locus tag', f.locus_tag);
  row('location', f.location);
  row('strand', f.strand === '-' ? 'reverse (−)' : 'forward (+)');
  row('length', f.key === 'CDS' ? `${fmt.count(f.length)} bp, ${fmt.count(Math.floor(f.length / 3) - 1)} aa` : `${fmt.count(f.length)} bp`);
  row('GC', fmt.percent(f.gc));
  if (f.gc3 !== undefined) row('GC at third positions', fmt.percent(f.gc3));
  if (f.nc !== undefined) row('effective codons (Wright)', fmt.number(f.nc, 3));
  card.append(kv);
  const lo = Math.min(...f.parts.map((p) => p[0])), hi = Math.max(...f.parts.map((p) => p[1]));
  card.append(el('div', { class: 'row' }, el('button', { class: 'button small', onclick: () => setWindow(lo - (hi - lo) * 0.4 - 100, hi + (hi - lo) * 0.4 + 100) }, 'show it')));
  if (f.translation) {
    card.append(el('div', { class: 'panel' },
      el('p', { class: 'eyebrow' }, 'protein', f.agrees === true ? el('span', { class: 'badge good' }, 'as RefSeq translates it') : f.agrees === false ? el('span', { class: 'badge bad' }, 'differs from RefSeq’s') : null),
      el('div', { class: 'protein' }, f.translation)));
  }
  if (f.codons) {
    const max = Math.max(1, ...f.codons.map((c) => c[2]));
    const grid = el('div', { class: 'codon-grid' });
    // NCBI's order: the first base by quarter, the second by column, the third by row.
    for (let first = 0; first < 4; first++)
      for (let third = 0; third < 4; third++)
        for (let second = 0; second < 4; second++) {
          const [codon, aa, n] = f.codons[16 * first + 4 * second + third];
          grid.append(el('div', { style: { background: warm(0.12 + 0.8 * (n / max)), color: n / max > 0.6 ? '#fff' : '' }, title: `${codon} ${aa}: ${n}` }, el('span', {}, `${codon} ${aa}`), el('span', {}, String(n))));
        }
    card.append(el('div', { class: 'panel' }, el('p', { class: 'eyebrow' }, 'codons used'), grid));
  }
  if (f.qualifiers?.length) {
    const list = el('dl', { class: 'kv small' });
    for (const [k, v] of f.qualifiers) if (k !== 'translation') list.append(el('dt', {}, k), el('dd', { style: { wordBreak: 'break-word' } }, v || '✓'));
    card.append(el('div', { class: 'panel' }, el('p', { class: 'eyebrow' }, 'qualifiers'), list));
  }
}

export { angleOf };
