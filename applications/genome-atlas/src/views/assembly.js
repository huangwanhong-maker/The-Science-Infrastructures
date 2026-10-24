// The Assembly view: a genome rebuilt from its reads. Reads are simulated from a genome (with
// errors), cut into k-mers and counted; the spectrum of counts separates the genome's k-mers
// (around the coverage) from the errors' (at one); the solid k-mers make a de Bruijn graph, whose
// non-branching paths are the unitigs; errors leave tips and bubbles, which cleaning removes step
// by step; what remains are the contigs, drawn on the genome they came from. And, for learning,
// the graph of a short text's k-mers, its Eulerian circuits, and how many texts share its k-mer
// composition: the BEST theorem, exactly.

import { el, clear, esc, notice, toggles, field, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, warm, CATEGORIES } from '../palette.js';

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, controls: el('div'), spectrum: el('div'), graph: el('div'), contigs: el('div'), euler: el('div'),
           p: { genome: 'sim:asm-30k', coverage: 30, error: 0.005, k: 31, seed: 1 }, result: null, step: null };
  section.append(
    el('h2', { class: 'view-title' }, 'Assembly'),
    el('p', { class: 'lede' }, 'A sequencer reads a genome as millions of short, overlapping, slightly wrong pieces. Cut them into k-mers and count: the genome\'s k-mers appear about as often as the coverage, the errors\' once or twice. The solid k-mers, joined where k − 1 letters overlap, make a de Bruijn graph; walking it without choice gives unitigs; a repeat longer than k is a junction the reads cannot resolve, which is why assemblies break where genomes repeat. A genome is an Eulerian path of its own graph, and the BEST theorem counts how many texts share its k-mers.'),
    view.controls,
    el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '28px' } }, el('div', {}, el('h3', { class: 'section-title' }, 'The k-mer spectrum'), view.spectrum), el('div', { style: { flex: 1, minWidth: '420px' } }, el('h3', { class: 'section-title' }, 'The de Bruijn graph'), view.graph)),
    el('h3', { class: 'section-title' }, 'The contigs, on the genome'), view.contigs,
    el('h3', { class: 'section-title' }, 'Eulerian circuits and the BEST theorem'), view.euler);
  buildControls();
  buildEuler();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.result) assemble();
}

function buildControls() {
  const p = view.p, box = clear(view.controls);
  const genome = el('select', { onchange: () => { p.genome = genome.value; } },
    el('option', { value: 'sim:asm-30k', selected: p.genome === 'sim:asm-30k' }, 'simulated, 30 kb with repeats'),
    el('option', { value: 'ncbi:NC_001422.1', selected: p.genome === 'ncbi:NC_001422.1' }, 'phiX174 (5.4 kb)'),
    el('option', { value: 'ncbi:NC_001416.1', selected: p.genome === 'ncbi:NC_001416.1' }, 'phage lambda (48.5 kb)'),
    el('option', { value: 'sim:repeats-200k', selected: p.genome === 'sim:repeats-200k' }, 'simulated, 200 kb with repeats'),
    el('option', { value: 'ncbi:NC_000913.3', selected: p.genome === 'ncbi:NC_000913.3' }, 'Escherichia coli (4.6 Mb, some seconds)'));
  // k is odd, so that no k-mer is its own reverse complement, and at most 31, two bits a base in a word.
  const k = el('select', { onchange: () => { p.k = Number(k.value); } });
  for (const v of [11, 15, 19, 21, 25, 27, 31]) k.append(el('option', { value: v, selected: v === p.k }, `k = ${v}`));
  const cov = el('select', { onchange: () => { p.coverage = Number(cov.value); } });
  for (const v of [5, 10, 20, 30, 60]) cov.append(el('option', { value: v, selected: v === p.coverage }, `${v}×`));
  const err = el('select', { onchange: () => { p.error = Number(err.value); } });
  for (const v of [0, 0.002, 0.005, 0.01, 0.02]) err.append(el('option', { value: v, selected: v === p.error }, v ? `${(v * 100).toFixed(1)}% errors` : 'no errors'));
  box.append(el('div', { class: 'row' }, field('genome', genome), field('reads', cov), err, k, el('button', { class: 'button primary', onclick: assemble }, 'Assemble')));
}

async function assemble() {
  const p = view.p;
  clear(view.graph).append(el('p', { class: 'wait' }, 'Counting k-mers and building the graph…'));
  try {
    view.result = await api.run('assembly.run', { genome: p.genome, coverage: p.coverage, error: p.error, k: p.k, seed: p.seed, length: 150 });
    view.step = Math.max(0, view.result.graphs.length - 1);
    drawSpectrum();
    drawGraph();
    drawContigs();
  } catch (e) { clear(view.graph).append(el('p', { class: 'error' }, e.message)); }
}

function drawSpectrum() {
  const r = view.result, h = r.spectrum;
  const pts = h.map(([c, n]) => [c, n]).filter(([c, n]) => c >= 1 && n > 0);
  const maxC = Math.min(Math.max(...pts.map((x) => x[0])), Math.max(8, Math.ceil(r.peak * 2.5)));
  const shown = pts.filter(([c]) => c <= maxC);
  const box = clear(view.spectrum);
  const c = el('div');
  c.innerHTML = chart({ width: 440, height: 250, x: [1, maxC], y: [1, Math.max(...shown.map((x) => x[1])) * 1.5], log: 'y', xlabel: 'times a k-mer is seen', ylabel: 'k-mers',
    bars: shown.map(([cc, n]) => ({ x0: cc - 0.5, x1: cc + 0.5, y: n, colour: cc < r.threshold ? token('--ink-faint') : token('--accent'), gap: 0.5, title: `${fmt.count(n)} k-mers seen ${cc} time${cc === 1 ? '' : 's'}` })),
    marks: [{ x: r.peak, y: Math.max(...shown.map((x) => x[1])), text: `peak ${r.peak}` }] });
  box.append(c, el('p', { class: 'note', style: { maxWidth: '440px' } },
    `${fmt.count(r.reads)} reads of 150 bases. The errors' k-mers (grey, ${fmt.count(r.error_kmers)} distinct) are cut at ${r.threshold}, where the fitted mixture of errors and genome changes hands; ${fmt.count(r.solid_kmers)} solid k-mers remain. The peak at ${r.peak} is the coverage in k-mers, c(L − k + 1)/L. The genome's size, read from the spectrum as the solid k-mers' occurrences over that coverage, is ${fmt.length(r.size_estimate)}${r.length ? ` (it is ${fmt.length(r.length)})` : ''}.`));
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
  const r = view.result, box = clear(view.graph);
  const st = r.stats;
  const summary = el('p', { class: 'note' }, `After cleaning: ${fmt.count(st.contigs)} contigs, ${fmt.length(st.total)} in all, the largest ${fmt.length(st.largest)}, N50 ${fmt.length(st.n50)}${st.ng50 ? `, NG50 ${fmt.length(st.ng50)}` : ''}; counted, built and cleaned in ${fmt.seconds(r.seconds.count + r.seconds.spectrum + r.seconds.graph + r.seconds.clean)} on ${r.threads} threads.`);
  if (!r.graphs.length) {
    // Too large to draw: the steps alone.
    const table = el('table', { class: 'data' }, el('tr', {}, ...['step', 'removed', 'unitigs before', 'after'].map((h) => el('th', { class: h === 'step' ? '' : 'num' }, h))),
      ...r.steps.map((x) => el('tr', {}, el('td', {}, x.name), el('td', { class: 'num' }, fmt.count(x.removed)), el('td', { class: 'num' }, fmt.count(x.before)), el('td', { class: 'num' }, fmt.count(x.after)))));
    box.append(el('p', { class: 'note' }, `${fmt.count(r.raw_unitigs)} unitigs before cleaning: too many to draw one by one. The cleaning, step by step:`), table, summary);
    return;
  }
  view.step = Math.max(0, Math.min(r.graphs.length - 1, view.step));
  const step = view.step;   // 0: the raw graph; i: after the i-th step of cleaning
  const frame = r.graphs[step];
  const nodes = frame.unitigs.map(([length, coverage], id) => ({ id, length, coverage }));
  const links = frame.links.map(([from, to]) => ({ from, to }));
  const doomed = new Set(frame.removed);
  if (nodes.length > 700) {
    box.append(el('p', { class: 'note' }, `${fmt.count(nodes.length)} unitigs: too many to draw one by one; choose fewer errors, a larger k or a smaller genome.`));
  } else {
    const pos = layout(nodes, links);
    const xs = pos.map((q) => q.x), ys = pos.map((q) => q.y);
    const x0 = Math.min(...xs) - 20, x1 = Math.max(...xs) + 20, y0 = Math.min(...ys) - 20, y1 = Math.max(...ys) + 20;
    const W = Math.max(380, Math.min(720, view.section.clientWidth - 520)), H = 420;
    const sc = Math.min(W / (x1 - x0), H / (y1 - y0));
    const X = (x) => (x - x0) * sc, Y = (y) => (y - y0) * sc;
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
    const fig = el('div');
    fig.innerHTML = svg;
    box.append(fig);
  }
  // The steps of cleaning, one at a time.
  const label = step === 0 ? 'the raw graph' : `after ${r.steps[step - 1].name}`;
  const next = step < r.steps.length ? `; the next step, ${r.steps[step].name}, removes ${fmt.count(r.steps[step].removed)} (red)` : r.steps.length ? '; cleaned' : '; nothing to clean';
  box.prepend(el('div', { class: 'row' }, el('span', { class: 'faint small' }, 'cleaning:'),
    el('button', { class: 'button small', disabled: step <= 0, onclick: () => { view.step--; drawGraph(); } }, '‹'),
    el('button', { class: 'button small', disabled: step >= r.graphs.length - 1, onclick: () => { view.step++; drawGraph(); } }, '›'),
    el('span', { class: 'small' }, `${label}: ${fmt.count(nodes.length)} unitigs, ${fmt.count(links.length)} links${next}`)));
  box.append(summary, el('p', { class: 'note' }, 'Each disc is a unitig, sized by its length and coloured by its coverage; lines join unitigs that overlap by k − 1 bases. A tip is a dead end an error makes near a read\'s end; a bubble, two short arms between the same two unitigs, an error inside a read (or a heterozygous site); a low-coverage unitig, an error the reads repeat.'));
}

// The contigs placed on the genome they came from (by their exact k-mers).
function drawContigs() {
  const r = view.result, box = clear(view.contigs);
  if (!r.placements?.length) { box.append(el('p', { class: 'note' }, 'The contigs could not be placed on the genome.')); return; }
  const L = r.length, W = Math.min(980, view.section.clientWidth - 70), rowH = 7;
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
  const fig = el('div');
  fig.innerHTML = svg;
  box.append(el('p', { class: 'note' }, `The genome (the dark line), its repeats (ochre) and the contigs placed on it (mulberry forward, blue reverse): ${fmt.percent(r.covered)} of the genome is under a contig placed by a k-mer it holds once (${fmt.count(r.placements.length)} of ${fmt.count(r.stats.contigs)} contigs; the others lie within repeats). Contigs end where a repeat longer than k makes the graph branch.`), fig);
}

// ---- Eulerian circuits ----------------------------------------------------------------------------------------
function buildEuler() {
  const box = clear(view.euler);
  const text = el('input', { value: 'TAATGCCATGGGATGTT', style: { width: '320px', font: '14px var(--mono)' } });
  const k = el('select', {});
  for (const v of [3, 4, 5, 6]) k.append(el('option', { value: v, selected: v === 3 }, `k = ${v}`));
  const out = el('div');
  const go = el('button', { class: 'button', onclick: async () => {
    try {
      const r = await api.run('assembly.eulerian', { text: text.value.trim(), k: Number(k.value) });
      drawEuler(out, r);
    } catch (e) { notice(e.message, true); }
  } }, 'Build its graph');
  box.append(el('p', { class: 'note' }, 'The de Bruijn graph of a text\'s k-mers has the (k − 1)-mers as nodes and one edge for each k-mer; the text is an Eulerian trail through it, using every edge once. When the graph has several such trails, the k-mers cannot tell the texts apart. The circuits number ec(G) = t_w(G) ∏_v (deg⁺ v − 1)!, the BEST theorem (van Aardenne-Ehrenfest and de Bruijn 1951), with t_w, the arborescences towards a vertex, counted exactly by the matrix-tree theorem as a determinant over the integers (morphe). Compeau and Pevzner\'s example below has two texts.'),
    el('div', { class: 'row' }, field('text', text), k, go), out);
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
  const facts = el('dl', { class: 'kv' });
  const row = (k, v) => facts.append(el('dt', {}, k), el('dd', {}, v));
  row('nodes, edges', `${n} (k − 1)-mers, ${r.edges.length} k-mers`);
  row('an Eulerian path spells', el('span', { class: 'mono' }, r.spelled || '–'));
  row('the text', r.spelled === r.text ? el('span', { class: 'badge good' }, 'recovered') : el('span', { class: 'badge warn' }, 'another text with the same k-mers'));
  row('Eulerian trails (BEST)', el('span', { class: 'mono' }, r.trails));
  row('texts with these k-mers', el('span', { class: 'mono' }, r.texts));
  box.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '28px' } }, fig, el('div', { style: { maxWidth: '460px' } }, facts,
    el('p', { class: 'note' }, 'The edges are coloured in the order the trail takes them, light to dark. The BEST theorem counts Eulerian circuits; a trail is a circuit of the graph closed by one more edge, from the trail\'s end to its start, opened there. Trails that differ only in the order they take the copies of a repeated k-mer spell the same text, so the texts are the trails divided by ∏ m(x)! over the k-mers seen m(x) times.'))));
}
