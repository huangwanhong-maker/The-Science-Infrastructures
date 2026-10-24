// The Assembly view: a genome rebuilt from its reads. Reads are simulated from a genome (with
// errors), cut into k-mers and counted; the spectrum of counts separates the genome's k-mers
// (around the coverage) from the errors' (at one); the solid k-mers make a de Bruijn graph, whose
// non-branching paths are the unitigs; errors leave tips and bubbles, which cleaning removes step
// by step; what remains are the contigs, drawn on the genome they came from. And, for learning,
// the graph of a short text's k-mers, its Eulerian circuits, and how many texts share its k-mer
// composition: the BEST theorem, exactly. The instruments stand in the bench's column and the
// results on its stage; the values are kept in the address.

import { el, clear, esc } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, warm } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, iconBtn, seg, select, fld, card, metrics, callout, empty, table, kv, working, busy, onEnter } from '../ui/kit.js';

const GENOMES = [
  ['sim:asm-30k', 'simulated, 30 kb with repeats'],
  ['ncbi:NC_001422.1', 'phiX174 (5.4 kb)'],
  ['ncbi:NC_001416.1', 'phage lambda (48.5 kb)'],
  ['sim:repeats-200k', 'simulated, 200 kb with repeats'],
  ['ncbi:NC_000913.3', 'Escherichia coli (4.6 Mb, some seconds)'],
];
const COVERAGES = [5, 10, 20, 30, 60];
const ERRORS = [0, 0.002, 0.005, 0.01, 0.02];
// k is odd, so that no k-mer is its own reverse complement, and at most 31, two bits a base in a word.
const KS = [11, 15, 19, 21, 25, 27, 31];
const EULER_KS = [3, 4, 5, 6];
const EXAMPLE = 'TAATGCCATGGGATGTT';

// The values kept in the address: [field, key, read].
const oneOf = (list) => (v) => (list.includes(v) ? v : undefined);
const numberOf = (list) => (v) => (list.includes(Number(v)) ? Number(v) : undefined);
const HASH = [
  ['genome', 'from', oneOf(GENOMES.map(([g]) => g))], ['coverage', 'cov', numberOf(COVERAGES)], ['error', 'err', numberOf(ERRORS)],
  ['k', 'k', numberOf(KS)], ['seed', 'seed', (v) => (/^\d{1,6}$/.test(v) ? Number(v) : undefined)],
];

let view = null;
let keysBound = false;

export function mount(section, ctx) {
  view = { ctx, section, p: { genome: 'sim:asm-30k', coverage: 30, error: 0.005, k: 31, seed: 1 }, e: { text: EXAMPLE, k: 3 },
           result: null, ran: null, step: null, wantStep: null, euler: null, running: false, eulering: false };
  adopt(ctx.params);
  view.controls = el('div', { class: 'controls' });
  view.state = el('div', { class: 'results-head' });
  view.box = el('div', { class: 'stack', style: { gap: '16px' } });
  view.eulerCard = buildEuler();
  section.append(
    el('p', { class: 'lede' }, 'A sequencer reads a genome as millions of short, overlapping, slightly wrong pieces. Cut them into k-mers and count: the genome\'s k-mers appear about as often as the coverage, the errors\' once or twice. The solid k-mers, joined where k − 1 letters overlap, make a de Bruijn graph; walking it without choice gives unitigs; a repeat longer than k is a junction the reads cannot resolve, which is why assemblies break where genomes repeat. A genome is an Eulerian path of its own graph, and the BEST theorem counts how many texts share its k-mers.'),
    el('div', { class: 'bench' }, view.controls, el('div', { class: 'stage' }, view.state, view.box, view.eulerCard)));
  buildControls();
  if (!keysBound) { keysBound = true; window.addEventListener('keydown', onKey); }
}

export function show(ctx) {
  view.ctx = ctx;
  const { changed, euler } = adopt(ctx.params);
  if (changed || euler) buildControls();
  if (euler && !view.eulering) eulerian();
  if (view.running) return;
  if (!view.result || changed) { assemble(); return; }
  if (view.wantStep !== null && view.result.graphs.length) { view.step = view.wantStep; view.wantStep = null; drawGraph(); }
  remember();
}

// The address's values, taken where they are valid: whether the reads differ, and whether a
// text's graph is asked for that is not the one drawn.
function adopt(params) {
  let changed = false;
  for (const [f, k, read] of HASH) {
    if (params[k] === undefined) continue;
    const v = read(params[k]);
    if (v !== undefined && v !== view.p[f]) { view.p[f] = v; changed = true; }
  }
  if (/^\d{1,3}$/.test(params.step || '') && Number(params.step) !== view.step) view.wantStep = Number(params.step);
  let euler = false;
  const text = /^[ACGT]{2,400}$/i.test(params.text || '') ? params.text.toUpperCase() : null;
  const ek = EULER_KS.includes(Number(params.ek)) ? Number(params.ek) : view.e.k;
  if (text && (text !== view.euler?.text || ek !== view.euler?.k)) { view.e.text = text; view.e.k = ek; euler = true; }
  return { changed, euler };
}

function remember() {
  if (view.ctx.view !== 'assembly') return;
  const p = view.p;
  view.ctx.setParams({ from: p.genome, cov: p.coverage, err: p.error, k: p.k, seed: p.seed, step: view.result?.graphs.length ? view.step : null,
    text: view.euler ? view.euler.text : null, ek: view.euler ? view.euler.k : null });
}

const swatch = (markup, w = 18) => { const s = el('span', { class: 'sw' }); s.innerHTML = `<svg width="${w}" height="10" viewBox="0 0 ${w} 10">${markup}</svg>`; return s; };
const key = (markup, label, w) => el('span', { class: 'k' }, swatch(markup, w), label);
const genomeLabel = (g) => GENOMES.find(([x]) => x === g)?.[1].replace(/ \(.*\)$/, '') || g;

// ---- the instruments -------------------------------------------------------------------------------------
function buildControls() {
  const p = view.p, e = view.e;
  view.runBtn = btn('Assemble', { icon: 'play', kind: 'primary', onClick: () => assemble() });
  view.againBtn = btn('Other reads', { icon: 'refresh', title: 'other reads of the same genome: the next seed', onClick: () => { p.seed++; assemble(); } });
  if (view.running) { busy(view.runBtn, true, 'Assembling…'); view.againBtn.disabled = true; }
  const reads = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'The reads'),
      fld('from the genome', select(GENOMES, p.genome, (v) => { p.genome = v; })),
      fld('coverage', seg(COVERAGES.map((v) => [v, `${v}×`]), p.coverage, (v) => { p.coverage = v; }, { block: true }), { help: '150 bases a read' }),
      fld('sequencing errors', seg(ERRORS.map((v) => [v, v ? `${(v * 100).toFixed(1)}%` : 'none']), p.error, (v) => { p.error = v; }, { block: true }))),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'The graph'),
      fld('k', seg(KS.map((v) => [v, String(v)]), p.k, (v) => { p.k = v; }, { block: true }), { hint: 'Odd, so that no k-mer is its own reverse complement; at most 31, two bits a base in a word.' })),
    el('div', { class: 'run-bar' }, view.runBtn, view.againBtn));
  onEnter(reads, assemble);
  view.eulerBtn = btn('Build its graph', { icon: 'play', kind: 'primary', onClick: () => eulerian() });
  if (view.eulering) busy(view.eulerBtn, true, 'Counting…');
  const text = el('input', { type: 'text', class: 'asm-text', value: e.text, spellcheck: 'false', placeholder: 'a short text of A, C, G and T', oninput: () => { e.text = text.value.trim().toUpperCase(); } });
  const circuits = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Eulerian circuits'),
      fld('a text', text, { hint: el('span', {}, 'Compeau and Pevzner\'s example has two texts with the same 3-mers. ', el('button', { class: 'link', onclick: () => { text.value = EXAMPLE; e.text = EXAMPLE; } }, 'Bring it back')) }),
      fld('k', seg(EULER_KS.map((v) => [v, String(v)]), e.k, (v) => { e.k = v; }, { block: true }))),
    el('div', { class: 'run-bar' }, view.eulerBtn));
  onEnter(circuits, eulerian);
  const own = el('div', { class: 'toolbox-hint' }, icon('toolbox', 'sm'), el('span', {},
    'These reads are simulated from a genome, so the contigs can be placed on it. With your own reads, ',
    el('button', { class: 'link', onclick: () => view.ctx.openTool('assembly-reads', { k: p.k }) }, 'assemble them'), ' in the Toolbox.'));
  clear(view.controls).append(reads, circuits, own);
}

// ---- assembling ---------------------------------------------------------------------------------------------
async function assemble() {
  if (view.running) return;
  view.running = true;
  const p = { ...view.p };
  busy(view.runBtn, true, 'Assembling…');
  view.againBtn.disabled = true;
  view.state.replaceChildren(working(`Simulating ${p.coverage}× of reads, counting their ${p.k}-mers, building and cleaning the graph…`));
  if (view.result) view.box.classList.add('stale');
  else view.box.replaceChildren(el('div', { class: 'skeleton', style: { height: '74px' } }), el('div', { class: 'skeleton', style: { height: '310px' } }), el('div', { class: 'skeleton', style: { height: '480px' } }));
  remember();
  try {
    view.result = await api.run('assembly.run', { genome: p.genome, coverage: p.coverage, error: p.error, k: p.k, seed: p.seed, length: 150 });
    view.ran = { p, at: new Date() };
    view.step = view.wantStep ?? Math.max(0, view.result.graphs.length - 1);
    view.wantStep = null;
    drawResults();
    remember();
  } catch (e) {
    view.state.replaceChildren(callout(`The reads could not be assembled: ${e.message}`, 'bad'));
    if (!view.result) clear(view.box);
  } finally {
    view.running = false;
    view.box.classList.remove('stale');
    busy(view.runBtn, false);
    view.againBtn.disabled = false;
  }
}

function drawResults() {
  const r = view.result, q = view.ran.p, st = r.stats;
  view.state.replaceChildren(el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Result'),
    el('span', { class: 'when' }, `${view.ran.at.toLocaleTimeString()} · ${genomeLabel(q.genome)}, ${q.coverage}× at ${q.error ? `${(q.error * 100).toFixed(1)}% errors` : 'no errors'}, k = ${q.k}, seed ${q.seed}`));
  view.graphCard = card({ title: 'The de Bruijn graph', sub: 'cleaned step by step: choose a step, or press ← and →', body: el('div') });
  view.graphSub = view.graphCard.querySelector('.card-head .sub');
  clear(view.box).append(
    metrics([
      [fmt.count(st.contigs), 'contigs'],
      [fmt.length(st.n50), `N50${st.ng50 ? `, NG50 ${fmt.length(st.ng50)}` : ''}`],
      [fmt.length(st.largest), 'the largest'],
      [fmt.length(st.total), 'in all'],
      [fmt.length(r.size_estimate), r.length ? `the size, from the spectrum (it is ${fmt.length(r.length)})` : 'the size, from the spectrum'],
      r.placements?.length ? [fmt.percent(r.covered), 'of the genome under a contig'] : null,
    ]),
    spectrumCard(), view.graphCard, contigsCard());
  drawGraph();
}

function spectrumCard() {
  const r = view.result, h = r.spectrum;
  const pts = h.map(([c, n]) => [c, n]).filter(([c, n]) => c >= 1 && n > 0);
  const maxC = Math.min(Math.max(...pts.map((x) => x[0])), Math.max(8, Math.ceil(r.peak * 2.5)));
  const shown = pts.filter(([c]) => c <= maxC);
  const c = el('div', { class: 'asm-chart' });
  c.innerHTML = chart({ width: 440, height: 250, x: [1, maxC], y: [1, Math.max(...shown.map((x) => x[1])) * 1.5], log: 'y', xlabel: 'times a k-mer is seen', ylabel: 'k-mers',
    bars: shown.map(([cc, n]) => ({ x0: cc - 0.5, x1: cc + 0.5, y: n, colour: cc < r.threshold ? token('--ink-faint') : token('--accent'), gap: 0.5, title: `${fmt.count(n)} k-mers seen ${cc} time${cc === 1 ? '' : 's'}` })),
    marks: [{ x: r.peak, y: Math.max(...shown.map((x) => x[1])), text: `peak ${r.peak}` }] });
  return card({ title: 'The k-mer spectrum', sub: `the errors cut at ${r.threshold}, the peak at ${r.peak}`, pad: true, body: el('div', { class: 'asm-split' }, c,
    el('div', { class: 'stack' },
      kv([['reads', `${fmt.count(r.reads)} of 150 bases`], ['distinct k-mers of the errors', fmt.count(r.error_kmers)], ['solid k-mers', fmt.count(r.solid_kmers)], ['coverage in k-mers', fmt.number(r.kmer_coverage, 3)]]),
      el('p', { class: 'note' }, `The errors' k-mers (grey) are cut at ${r.threshold}, where the fitted mixture of errors and genome changes hands. The peak at ${r.peak} is the coverage in k-mers, c(L − k + 1)/L. The genome's size, read from the spectrum as the solid k-mers' occurrences over that coverage, is ${fmt.length(r.size_estimate)}${r.length ? ` (it is ${fmt.length(r.length)})` : ''}.`))) });
}

// ---- the graph, by a force layout ----------------------------------------------------------------------
function layout(nodes, links, iterations = 260) {
  const n = nodes.length;
  const pos = nodes.map((_, i) => ({ x: Math.cos((2 * Math.PI * i) / n) * 200 + (i % 7), y: Math.sin((2 * Math.PI * i) / n) * 200 + (i % 5), vx: 0, vy: 0 }));
  const index = new Map(nodes.map((d, i) => [d.id, i]));
  const edges = links.map((l) => [index.get(l.from), index.get(l.to)]).filter(([a, b]) => a !== undefined && b !== undefined && a !== b);
  for (let it = 0; it < iterations; it++) {
    const t = 1 - it / iterations;
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) {
        const dx = pos[i].x - pos[j].x, dy = pos[i].y - pos[j].y;
        const d2 = Math.max(4, dx * dx + dy * dy), f = 900 / d2;
        pos[i].vx += dx * f; pos[i].vy += dy * f; pos[j].vx -= dx * f; pos[j].vy -= dy * f;
      }
    for (const [a, b] of edges) {
      const dx = pos[b].x - pos[a].x, dy = pos[b].y - pos[a].y, d = Math.sqrt(dx * dx + dy * dy) + 1e-9;
      const want = 14 + 0.6 * (Math.log2(1 + nodes[a].length) + Math.log2(1 + nodes[b].length));
      const f = 0.06 * (d - want) / d;
      pos[a].vx += dx * f; pos[a].vy += dy * f; pos[b].vx -= dx * f; pos[b].vy -= dy * f;
    }
    for (const q of pos) { q.vx -= q.x * 0.004; q.vy -= q.y * 0.004; q.x += Math.max(-12, Math.min(12, q.vx * t)); q.y += Math.max(-12, Math.min(12, q.vy * t)); q.vx *= 0.5; q.vy *= 0.5; }
  }
  return pos;
}

function drawGraph() {
  const r = view.result, box = clear(view.graphCard.body);
  const st = r.stats;
  const summary = el('p', { class: 'note' }, `After cleaning: ${fmt.count(st.contigs)} contigs, ${fmt.length(st.total)} in all, the largest ${fmt.length(st.largest)}, N50 ${fmt.length(st.n50)}${st.ng50 ? `, NG50 ${fmt.length(st.ng50)}` : ''}; counted, built and cleaned in ${fmt.seconds(r.seconds.count + r.seconds.spectrum + r.seconds.graph + r.seconds.clean)} on ${r.threads} threads.`);
  if (!r.graphs.length) {
    // Too large to draw: the steps alone.
    view.graphSub.textContent = `${fmt.count(r.raw_unitigs)} unitigs before cleaning: too many to draw one by one`;
    box.append(el('div', { class: 'card-body' }, el('p', { class: 'eyebrow' }, 'the cleaning, step by step'),
      table([{ key: 'name', label: 'step' }, { key: 'removed', label: 'removed', num: true, format: fmt.count }, { key: 'before', label: 'unitigs before', num: true, format: fmt.count }, { key: 'after', label: 'after', num: true, format: fmt.count }], r.steps, { wide: false }), summary));
    return;
  }
  view.step = Math.max(0, Math.min(r.graphs.length - 1, view.step));
  const step = view.step;   // 0: the raw graph; i: after the i-th step of cleaning
  const frame = r.graphs[step];
  const nodes = frame.unitigs.map(([length, coverage], id) => ({ id, length, coverage }));
  const links = frame.links.map(([from, to]) => ({ from, to }));
  const doomed = new Set(frame.removed);
  const label = step === 0 ? 'the raw graph' : `after ${r.steps[step - 1].name}`;
  const next = step < r.steps.length ? `; the next step, ${r.steps[step].name}, removes ${fmt.count(r.steps[step].removed)} (red)` : r.steps.length ? '; cleaned' : '; nothing to clean';
  view.graphSub.textContent = `${label}: ${fmt.count(nodes.length)} unitigs, ${fmt.count(links.length)} links${next}`;
  const fig = el('div', { class: 'asm-graph' });
  const steps = el('div', { class: 'asm-steps' }, el('div', { class: 'eyebrow' }, 'the cleaning, step by step'),
    ...r.graphs.map((g, i) => el('button', { class: `asm-step${i === step ? ' on' : ''}`, onclick: () => { view.step = i; drawGraph(); remember(); } },
      el('span', { class: 'n' }, String(i)),
      el('span', { class: 't' }, el('span', {}, i === 0 ? 'the raw graph' : r.steps[i - 1].name),
        el('i', {}, i === 0 ? `${fmt.count(g.unitigs.length)} unitigs` : `${fmt.count(r.steps[i - 1].removed)} removed, ${fmt.count(g.unitigs.length)} unitigs left`)))),
    el('div', { class: 'row', style: { gap: '4px', marginTop: '6px' } },
      flipped(iconBtn('chevron', { title: 'a step back (←)', size: 'sm', onClick: () => goStep(-1) })), iconBtn('chevron', { title: 'a step on (→)', size: 'sm', onClick: () => goStep(1) }),
      el('span', { class: 'small faint' }, `${step} of ${r.graphs.length - 1}`)));
  box.append(el('div', { class: 'asm-graph-body' }, fig, steps));
  if (nodes.length > 700) {
    fig.append(el('p', { class: 'note' }, `${fmt.count(nodes.length)} unitigs: too many to draw one by one; choose fewer errors, a larger k or a smaller genome.`));
  } else {
    const pos = layout(nodes, links);
    const xs = pos.map((q) => q.x), ys = pos.map((q) => q.y);
    const x0 = Math.min(...xs) - 20, x1 = Math.max(...xs) + 20, y0 = Math.min(...ys) - 20, y1 = Math.max(...ys) + 20;
    const W = Math.max(380, Math.min(720, fig.clientWidth > 300 ? fig.clientWidth - 4 : 640)), H = 420;
    // The largest disc's radius kept inside the frame, and the drawing centred in it.
    const pad = Math.min(48, 2 + Math.sqrt(Math.max(...nodes.map((u) => u.length))) * 0.35 + 2);
    const sc = Math.min((W - 2 * pad) / (x1 - x0), (H - 2 * pad) / (y1 - y0));
    const ox = pad + (W - 2 * pad - (x1 - x0) * sc) / 2, oy = pad + (H - 2 * pad - (y1 - y0) * sc) / 2;
    const X = (x) => ox + (x - x0) * sc, Y = (y) => oy + (y - y0) * sc;
    const maxCov = Math.max(...nodes.map((u) => u.coverage));
    let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" class="dbg">`;
    for (const l of links) {
      const a = pos[l.from], b = pos[l.to];
      if (l.from === l.to) { svg += `<circle cx="${X(a.x).toFixed(1)}" cy="${(Y(a.y) - 9).toFixed(1)}" r="7" fill="none" stroke="${token('--hair-2')}" />`; continue; }
      svg += `<line x1="${X(a.x).toFixed(1)}" y1="${Y(a.y).toFixed(1)}" x2="${X(b.x).toFixed(1)}" y2="${Y(b.y).toFixed(1)}" stroke="${token('--hair-2')}" stroke-width="1" />`;
    }
    nodes.forEach((u, i) => {
      const rr = 2 + Math.sqrt(u.length) * 0.35;
      const colour = doomed.has(u.id) ? token('--bad') : warm(0.25 + 0.7 * Math.min(1, u.coverage / maxCov));
      svg += `<circle cx="${X(pos[i].x).toFixed(1)}" cy="${Y(pos[i].y).toFixed(1)}" r="${rr.toFixed(1)}" fill="${colour}" stroke="${token('--paper-2')}" stroke-width="0.8"><title>unitig ${u.id + 1}: ${fmt.length(u.length)}, coverage ${fmt.number(u.coverage, 3)}${doomed.has(u.id) ? `; the next step removes it (${esc(r.steps[step].kind)})` : ''}</title></circle>`;
    });
    svg += '</svg>';
    fig.innerHTML = svg;
  }
  box.append(el('div', { class: 'key-row' },
    key(`<circle cx="5" cy="5" r="4" fill="${warm(0.35)}"/><circle cx="14" cy="5" r="4" fill="${warm(0.9)}"/>`, 'a unitig, sized by its length; deeper where its coverage is higher', 20),
    key(`<circle cx="5" cy="5" r="4" fill="${token('--bad')}"/>`, 'removed by the next step', 10),
    key(`<line x1="1" y1="5" x2="17" y2="5" stroke="${token('--hair-3')}"/>`, 'an overlap of k − 1 bases')),
  el('div', { class: 'card-body', style: { paddingTop: '4px' } }, summary,
    el('p', { class: 'note' }, 'A tip is a dead end an error makes near a read\'s end; a bubble, two short arms between the same two unitigs, an error inside a read (or a heterozygous site); a low-coverage unitig, an error the reads repeat.')));
}

function flipped(b) { b.classList.add('flip'); return b; }

function goStep(d) {
  const r = view.result;
  if (!r?.graphs.length) return;
  const s = Math.max(0, Math.min(r.graphs.length - 1, view.step + d));
  if (s === view.step) return;
  view.step = s;
  drawGraph();
  remember();
}

// The contigs placed on the genome they came from (by their exact k-mers).
function contigsCard() {
  const r = view.result;
  if (!r.placements?.length) return card({ title: 'The contigs, on the genome', pad: true, body: el('p', { class: 'note', style: { margin: 0 } }, 'The contigs could not be placed on the genome.') });
  const fig = el('div', { class: 'asm-contigs' });
  const c = card({ title: 'The contigs, on the genome', sub: `${fmt.percent(r.covered)} of the genome under a contig placed by a k-mer it holds once`, body: [fig,
    el('div', { class: 'key-row' },
      key(`<rect y="3" width="18" height="4" fill="${token('--ink')}"/>`, 'the genome'),
      key(`<rect y="1" width="18" height="8" fill="${token('--G')}" fill-opacity="0.7"/>`, 'its repeats'),
      key(`<rect y="2" width="18" height="6" fill="${token('--accent')}" fill-opacity="0.8"/>`, 'a contig, forward'),
      key(`<rect y="2" width="18" height="6" fill="${token('--C')}" fill-opacity="0.8"/>`, 'reverse')),
    el('div', { class: 'card-body', style: { paddingTop: '4px' } }, el('p', { class: 'note', style: { margin: 0 } }, `${fmt.count(r.placements.length)} of ${fmt.count(r.stats.contigs)} contigs are placed; the others lie within repeats. Contigs end where a repeat longer than k makes the graph branch.`))] });
  // Drawn once the card is in the page, at the width it has there.
  requestAnimationFrame(() => drawContigs(fig));
  return c;
}

function drawContigs(fig) {
  const r = view.result;
  const L = r.length, W = Math.min(1180, fig.clientWidth > 300 ? fig.clientWidth - 8 : 980), rowH = 7;
  const rows = [];
  const placed = [...r.placements].sort((a, b) => a.begin - b.begin).map((p) => {
    let row = rows.findIndex((end) => end < p.begin);
    if (row < 0) { row = rows.length; rows.push(-1); }
    rows[row] = p.end;
    return { ...p, row };
  });
  const H = 40 + Math.min(40, rows.length) * rowH;
  const x = (p) => 10 + (p / L) * (W - 20);
  let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`;
  svg += `<line x1="10" y1="16" x2="${W - 10}" y2="16" stroke="${token('--ink')}" stroke-width="2" />`;
  for (const [b, e, tandem] of r.repeats || []) svg += `<rect x="${x(b).toFixed(1)}" y="10" width="${Math.max(1, x(e) - x(b)).toFixed(1)}" height="12" fill="${token('--G')}" fill-opacity="0.7"><title>${tandem ? 'a tandem array' : 'a repeat copy'} (${e - b} bp)</title></rect>`;
  for (const p of placed) {
    if (p.row >= 40) continue;
    svg += `<rect x="${x(p.begin).toFixed(1)}" y="${30 + p.row * rowH}" width="${Math.max(1, x(p.end) - x(p.begin)).toFixed(1)}" height="${rowH - 2}" fill="${p.reverse ? token('--C') : token('--accent')}" fill-opacity="0.8"><title>a contig of ${fmt.length(p.length)} at ${fmt.interval(p.begin, p.end)}${p.reverse ? ' (reverse strand)' : ''}</title></rect>`;
  }
  svg += '</svg>';
  fig.innerHTML = svg;
}

// ---- Eulerian circuits ----------------------------------------------------------------------------------------
function buildEuler() {
  view.eulerState = el('div');
  view.eulerOut = el('div', {}, empty({ icon: 'assembly', title: 'No text\'s graph yet', text: `The graph of ${EXAMPLE}'s 3-mers is Compeau and Pevzner's example: two texts share them.`,
    action: btn('Build its graph', { icon: 'play', onClick: () => eulerian() }) }));
  return card({ title: 'Eulerian circuits and the BEST theorem', sub: 'how many texts share a text\'s k-mers', pad: true, body: [
    el('p', { class: 'note', style: { marginTop: 0 } }, 'The de Bruijn graph of a text\'s k-mers has the (k − 1)-mers as nodes and one edge for each k-mer; the text is an Eulerian trail through it, using every edge once. When the graph has several such trails, the k-mers cannot tell the texts apart. The circuits number ec(G) = t_w(G) ∏_v (deg⁺ v − 1)!, the BEST theorem (van Aardenne-Ehrenfest and de Bruijn 1951), with t_w, the arborescences towards a vertex, counted exactly by the matrix-tree theorem as a determinant over the integers (morphe).'),
    view.eulerState, view.eulerOut] });
}

async function eulerian() {
  if (view.eulering) return;
  const e = view.e;
  if (!/^[ACGT]+$/.test(e.text) || e.text.length <= e.k) { view.eulerState.replaceChildren(callout(`A text of A, C, G and T, longer than k = ${e.k}.`, 'bad')); return; }
  view.eulering = true;
  busy(view.eulerBtn, true, 'Counting…');
  view.eulerState.replaceChildren(working(`Building the graph of ${e.text.length < 24 ? e.text : 'the text'}'s ${e.k}-mers and counting its circuits…`));
  if (view.euler) view.eulerOut.classList.add('stale');
  try {
    const r = await api.run('assembly.eulerian', { text: e.text, k: e.k });
    view.euler = { r, text: e.text, k: e.k };
    clear(view.eulerState);
    drawEuler(view.eulerOut, r);
    remember();
  } catch (err) { view.eulerState.replaceChildren(callout(err.message, 'bad')); }
  finally { view.eulering = false; view.eulerOut.classList.remove('stale'); busy(view.eulerBtn, false); }
}

function drawEuler(box, r) {
  clear(box);
  const n = r.nodes.length, R = 150, W = 420, H = 380;
  const pos = r.nodes.map((_, i) => ({ x: W / 2 + R * Math.cos((2 * Math.PI * i) / n - Math.PI / 2), y: H / 2 + R * Math.sin((2 * Math.PI * i) / n - Math.PI / 2) }));
  let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}"><defs><marker id="arrow" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0L8,4L0,8z" fill="${token('--ink-2')}"/></marker></defs>`;
  const order = new Map((r.trail || []).map((e, k) => [e, k]));
  r.edges.forEach(([a, b, label], e) => {
    const p = pos[a], q = pos[b];
    const k = order.get(e);
    const colour = k === undefined ? token('--hair-2') : warm(0.25 + 0.7 * (k / Math.max(1, r.edges.length - 1)));
    if (a === b) {
      svg += `<circle cx="${p.x + (p.x - W / 2) * 0.12}" cy="${p.y + (p.y - H / 2) * 0.12}" r="12" fill="none" stroke="${colour}" stroke-width="1.6"><title>${esc(label)}</title></circle>`;
    } else {
      const mx = (p.x + q.x) / 2 + (q.y - p.y) * 0.12, my = (p.y + q.y) / 2 - (q.x - p.x) * 0.12;
      svg += `<path d="M${p.x},${p.y}Q${mx},${my} ${q.x},${q.y}" fill="none" stroke="${colour}" stroke-width="1.6" marker-end="url(#arrow)"><title>${esc(label)}${k !== undefined ? `, step ${k + 1}` : ''}</title></path>`;
    }
  });
  r.nodes.forEach((label, i) => {
    svg += `<circle cx="${pos[i].x}" cy="${pos[i].y}" r="13" fill="${token('--paper-2')}" stroke="${token('--ink-2')}" />`;
    svg += `<text x="${pos[i].x}" y="${pos[i].y + 4}" text-anchor="middle" style="font: 10.5px var(--mono)">${esc(label)}</text>`;
  });
  svg += '</svg>';
  const fig = el('div');
  fig.innerHTML = svg;
  const facts = kv([
    ['nodes, edges', `${n} (k − 1)-mers, ${r.edges.length} k-mers`],
    ['an Eulerian path spells', el('span', { class: 'mono' }, r.spelled || '–')],
    ['the text', r.spelled === r.text ? el('span', { class: 'badge good' }, 'recovered') : el('span', { class: 'badge warn' }, 'another text with the same k-mers')],
  ]);
  box.append(el('div', { class: 'asm-split' }, fig, el('div', { style: { maxWidth: '460px' } },
    metrics([[String(r.trails), 'Eulerian trails (BEST)'], [String(r.texts), 'texts with these k-mers', Number(r.texts) === 1 ? 'good' : '']]),
    facts,
    el('p', { class: 'note' }, 'The edges are coloured in the order the trail takes them, light to dark. The BEST theorem counts Eulerian circuits; a trail is a circuit of the graph closed by one more edge, from the trail\'s end to its start, opened there. Trails that differ only in the order they take the copies of a repeated k-mer spell the same text, so the texts are the trails divided by ∏ m(x)! over the k-mers seen m(x) times.'))));
}

// ---- keys ------------------------------------------------------------------------------------------------
function onKey(e) {
  if (!view || view.section.hidden || !view.result?.graphs.length) return;
  const a = document.activeElement;
  if (/^(INPUT|TEXTAREA|SELECT)$/.test(a?.tagName || '') || a?.isContentEditable || e.ctrlKey || e.metaKey || e.altKey || document.querySelector('.dialog-veil, .popover')) return;
  if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') { e.preventDefault(); goStep(e.key === 'ArrowLeft' ? -1 : 1); }
}
