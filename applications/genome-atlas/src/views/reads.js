// The Reads view: a diploid made from a genome (SNVs and short indels, with their truth), read
// at a coverage, mapped back and called. Its parts: the calls against the truth, type by type; a
// stretch of the genome read by read (the pileup), with the calls and the truth above it; how
// well the mapper places reads at each mapping quality, measured on a job spread over the
// workers (Li 2018, Fig. 1); and the reads and calls as SAM and VCF.

import { el, clear, esc, notice, toggles, field, hover } from '../dom.js';
import * as api from '../api.js';
import { chart, timeline } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, BASE_COLOURS, CATEGORIES } from '../palette.js';

const GENOMES = [
  ['ncbi:NC_001416.1', 'phage lambda (48.5 kb)'],
  ['sim:repeats-200k', 'simulated, 200 kb with repeats'],
  ['ncbi:NC_000913.3', 'Escherichia coli (4.6 Mb, half a minute)'],
];

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, controls: el('div'), summary: el('div'), browser: el('div'), mapping: el('div'), files: el('div'),
           p: { genome: GENOMES[0][0], kind: 'short', coverage: 30, snv: 0.001, indel: 0.0001, seed: 5 }, result: null, at: null, width: 240 };
  section.append(
    el('h2', { class: 'view-title' }, 'Reads'),
    el('p', { class: 'lede' }, 'A read is placed on its genome in three steps, each paid for only when the one before has narrowed the choice: its minimizers, looked up in an index, give exact anchors; a dynamic programme over the anchors finds the chains that run in the same order on read and genome; only the best chains are aligned base by base. Whether the place is right is a posterior, the mapping quality: a read in a repeat has two places equally good and a quality of zero. At each position of a diploid the reads then vote, each with its base\'s error probability, among three genotypes; an indel must first be written in one canonical way, for in a repeat it has many.'),
    view.controls,
    el('h3', { class: 'section-title' }, 'The calls against the truth'), view.summary,
    el('h3', { class: 'section-title' }, 'A stretch, read by read'), view.browser,
    el('h3', { class: 'section-title' }, 'Mapping quality, measured across the workers'), view.mapping,
    el('h3', { class: 'section-title' }, 'As files'), view.files);
  buildControls();
  buildMapping();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.result) call();
}

function params() {
  const p = view.p;
  return { genome: p.genome, kind: p.kind, coverage: p.coverage, snv: p.snv, indel: p.indel, seed: p.seed };
}

function buildControls() {
  const p = view.p, box = clear(view.controls);
  const genome = el('select', { onchange: () => { p.genome = genome.value; } }, ...GENOMES.map(([v, l]) => el('option', { value: v, selected: v === p.genome }, l)));
  const cov = el('select', { onchange: () => { p.coverage = Number(cov.value); } }, ...[10, 20, 30, 60].map((v) => el('option', { value: v, selected: v === p.coverage }, `${v}×`)));
  const snv = el('select', { onchange: () => { p.snv = Number(snv.value); } }, ...[0.0005, 0.001, 0.002, 0.005].map((v) => el('option', { value: v, selected: v === p.snv }, `SNVs ${fmt.number(v * 1000, 2)} per kb`)));
  const indel = el('select', { onchange: () => { p.indel = Number(indel.value); } }, ...[0.00005, 0.0001, 0.0003].map((v) => el('option', { value: v, selected: v === p.indel }, `indels ${fmt.number(v * 1000, 2)} per kb`)));
  box.append(el('div', { class: 'row' }, field('genome', genome),
    toggles([['short', 'short reads, 150 b'], ['long', 'long reads, 8 kb at 8%']], p.kind, (v) => { p.kind = v; }),
    cov, snv, indel,
    el('button', { class: 'button', onclick: () => { p.seed++; call(); } }, 'another diploid'),
    el('button', { class: 'button primary', onclick: call }, 'Read, map and call')));
}

async function call() {
  clear(view.summary).append(el('p', { class: 'wait' }, 'Making the diploid, reading it, mapping the reads and calling the variants…'));
  clear(view.browser);
  try {
    const r = await api.run('reads.call', params());
    view.result = r;
    // Open the browser on the first call that is an indel, or the first call.
    const first = r.calls.find((c) => c.indel) || r.calls[0];
    view.at = first ? { contig: first.contig, position: first.position } : { contig: 0, position: 1000 };
    drawSummary();
    drawFiles();
    browse();
  } catch (e) { clear(view.summary).append(el('p', { class: 'error' }, e.message)); }
}

// ---- the calls against the truth ------------------------------------------------------------------------
function drawSummary() {
  const r = view.result, box = clear(view.summary);
  const s = r.seconds;
  const rows = [['snv', 'SNVs'], ['insertion', 'insertions'], ['deletion', 'deletions'], ['all', 'all']];
  const table = el('table', { class: 'data' }, el('tr', {}, ...['', 'true', 'false', 'missed', 'precision', 'recall', 'F₁', 'genotypes'].map((h, k) => el('th', { class: k ? 'num' : '' }, h))),
    ...rows.map(([k, label]) => {
      const c = r.types[k];
      return el('tr', {}, el('td', {}, label), el('td', { class: 'num' }, fmt.count(c.tp)), el('td', { class: 'num' }, fmt.count(c.fp)), el('td', { class: 'num' }, fmt.count(c.fn)),
        el('td', { class: 'num' }, fmt.number(c.precision, 4)), el('td', { class: 'num' }, fmt.number(c.recall, 4)), el('td', { class: 'num' }, fmt.number(c.f1, 4)), el('td', { class: 'num' }, fmt.number(c.genotype, 4)));
    }));
  // Where the time went.
  const parts = [['simulate', 'reads simulated'], ['index', 'the index'], ['map', 'mapping'], ['call', 'calling'], ['compare', 'comparing']];
  const total = parts.reduce((a, [k]) => a + (s[k] || 0), 0);
  const bar = el('div', { class: 'stackbar' }, ...parts.map(([k, label], i) => el('span', { style: { width: `${(100 * (s[k] || 0)) / total}%`, background: CATEGORIES[i] }, title: `${label}: ${fmt.seconds(s[k] || 0)}` })));
  const legend = el('div', { class: 'row small faint' }, ...parts.map(([k, label], i) => el('span', {}, el('span', { class: 'swatch', style: { background: CATEGORIES[i] } }), `${label} ${fmt.seconds(s[k] || 0)}`)));
  box.append(el('p', { class: 'note' }, `${fmt.count(r.reads)} reads (${fmt.length(r.mean_length)} on average) of a diploid made from ${fmt.length(r.length)}, at a mean depth of ${fmt.number(r.mean_depth, 3)}; ${fmt.percent(r.correct / Math.max(1, r.reads), 2)} placed where they came from. ${fmt.count(r.truth)} variants were made and ${fmt.count(r.called)} called. A call is true when, written in canonical form (left-aligned, parsimonious: Tan, Abecasis and Kang 2015), it is one of the variants made; a genotype agrees when it has their number of copies.`),
    table, el('p', { class: 'eyebrow', style: { marginTop: '14px' } }, `the time, ${fmt.seconds(total)} on ${r.threads} threads`), bar, legend);
}

// ---- the pileup ------------------------------------------------------------------------------------------
async function browse() {
  const r = view.result, box = clear(view.browser);
  const at = view.at, half = Math.round(view.width / 2);
  const begin = Math.max(0, at.position - half), end = begin + view.width;
  // Navigation: the calls, in order, and the window's width.
  const calls = r.calls;
  const index = calls.findIndex((c) => c.contig === at.contig && c.position >= at.position);
  const go = (k) => { const c = calls[(k + calls.length) % calls.length]; if (c) { view.at = { contig: c.contig, position: c.position }; browse(); } };
  const pick = el('select', { onchange: () => go(Number(pick.value)) },
    ...calls.slice(0, 2000).map((c, k) => el('option', { value: k, selected: k === index }, `${fmt.position(c.position)} ${c.ref}→${c.alts.join(',')} ${c.genotype}${c.true === false ? ' (false)' : ''}`)));
  const missed = el('select', { onchange: () => { const v = r.missed[Number(missed.value)]; if (v) { view.at = { contig: v.contig, position: v.position }; browse(); } } },
    el('option', { value: '', disabled: true, selected: true }, `${r.missed.length} missed…`), ...r.missed.map((v, k) => el('option', { value: k }, `${fmt.position(v.position)} ${v.ref}→${v.alt}`)));
  box.append(el('div', { class: 'row' },
    el('button', { class: 'button small', onclick: () => go(index - 1) }, '‹ call'), el('button', { class: 'button small', onclick: () => go(index + 1) }, 'call ›'),
    field('calls', pick), missed,
    toggles([[120, '120 b'], [240, '240 b'], [1000, '1 kb'], [4000, '4 kb']], view.width, (v) => { view.width = v; browse(); })));
  const fig = el('div', { class: 'pileup' }, el('p', { class: 'wait' }, 'Reading the stretch…'));
  box.append(fig);
  let w;
  try { w = await api.run('reads.window', { ...params(), contig: at.contig, begin, end, max_reads: 1500 }); } catch (e) { clear(fig).append(el('p', { class: 'error' }, e.message)); return; }
  drawPileup(fig, w);
}

function drawPileup(fig, w) {
  clear(fig);
  const colours = BASE_COLOURS();
  const W = Math.min(1180, view.section.clientWidth - 40), left = 8, right = 8, span = w.end - w.begin;
  const x = (p) => left + ((p - w.begin) / span) * (W - left - right);
  const px = (W - left - right) / span;
  // Rows: each read in the first row whose last read ended before it begins.
  const reads = [...w.reads].sort((a, b) => a[0] - b[0]);
  const ends = [], row = [];
  for (const rd of reads) {
    let k = ends.findIndex((e) => e + Math.max(2, 3 / px) <= rd[0]);
    if (k < 0) { k = ends.length; ends.push(0); }
    ends[k] = rd[1];
    row.push(k);
  }
  const rowH = span <= 300 ? 9 : span <= 1200 ? 6 : 4, top = 70;
  const H = top + Math.min(ends.length, 120) * rowH + 12;
  const out = [`<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`];
  // The scale and the reference.
  const step = fmt.tickStep(span, 8);
  for (let p = Math.ceil(w.begin / step) * step; p < w.end; p += step) out.push(`<line x1="${x(p)}" y1="12" x2="${x(p)}" y2="17" stroke="${token('--ink-faint')}" /><text class="tick-label" x="${x(p)}" y="10" text-anchor="middle">${fmt.position(p)}</text>`);
  if (px >= 6) {
    for (let k = 0; k < span; k++) out.push(`<text x="${(x(w.begin + k) + px / 2).toFixed(1)}" y="30" text-anchor="middle" style="font:${Math.min(12, px * 1.3).toFixed(1)}px var(--mono);fill:${colours[w.reference[k]] || token('--N')}">${w.reference[k]}</text>`);
  } else {
    for (let k = 0; k < span; k++) out.push(`<rect x="${x(w.begin + k).toFixed(2)}" y="22" width="${Math.max(0.6, px).toFixed(2)}" height="8" fill="${colours[w.reference[k]] || token('--N')}" />`);
  }
  // The truth (ochre) and the calls (green when true, red when not), above the reads.
  const mark = (p, y, colour, title, shape) => out.push(shape === 'tri'
    ? `<path d="M${x(p) + px / 2 - 5},${y + 8}L${x(p) + px / 2 + 5},${y + 8}L${x(p) + px / 2},${y}Z" fill="${colour}"><title>${esc(title)}</title></path>`
    : `<path d="M${x(p) + px / 2},${y}L${x(p) + px / 2 + 5},${y + 5}L${x(p) + px / 2},${y + 10}L${x(p) + px / 2 - 5},${y + 5}Z" fill="${colour}"><title>${esc(title)}</title></path>`);
  out.push(`<text class="tick-label" x="${W - right}" y="44" text-anchor="end">made</text><text class="tick-label" x="${W - right}" y="58" text-anchor="end">called</text>`);
  for (const v of w.truth) mark(v.position, 36, token('--G'), `made: ${v.ref}→${v.alt} at ${fmt.position(v.position)}, ${v.dosage === 2 ? 'homozygous' : 'heterozygous'}`, 'tri');
  for (const c of w.calls) mark(c.position, 50, c.true === false ? token('--bad') : c.passed ? token('--good') : token('--ink-faint'),
    `called: ${c.ref}→${c.alts.join(',')} ${c.genotype}, QUAL ${c.qual}, GQ ${c.gq}, depth ${c.depth}, AD ${c.ad.join(',')}${c.filters.length ? `, filtered: ${c.filters.join(', ')}` : ''}${c.true === false ? ' — not among the variants made' : ''}`, 'diamond');
  // The reads.
  reads.forEach((rd, k) => {
    const r = row[k];
    if (r >= 120) return;
    const [b, e, reverse, mapq, mism, ins, del] = rd;
    const y = top + r * rowH;
    const a = Math.max(b, w.begin), z = Math.min(e, w.end);
    out.push(`<rect x="${x(a).toFixed(1)}" y="${y}" width="${Math.max(1, x(z) - x(a)).toFixed(1)}" height="${rowH - 1.5}" fill="${reverse ? token('--C') : token('--accent')}" fill-opacity="${mapq >= 20 ? 0.22 : 0.08}" stroke="${mapq === 0 ? token('--bad') : 'none'}" stroke-width="0.6"><title>${esc(rd[8])}: ${fmt.interval(b, e)}, ${reverse ? 'reverse' : 'forward'} strand, MAPQ ${mapq}</title></rect>`);
    for (const [t, base] of mism) out.push(`<rect x="${x(t).toFixed(1)}" y="${y}" width="${Math.max(1, px).toFixed(1)}" height="${rowH - 1.5}" fill="${colours[base] || token('--N')}" />`);
    for (const [t, len] of del) {
      const d0 = Math.max(t, w.begin), d1 = Math.min(t + len, w.end);
      out.push(`<rect x="${x(d0).toFixed(1)}" y="${y}" width="${Math.max(1, x(d1) - x(d0)).toFixed(1)}" height="${rowH - 1.5}" fill="${token('--paper')}" /><line x1="${x(d0)}" y1="${y + (rowH - 1.5) / 2}" x2="${x(d1)}" y2="${y + (rowH - 1.5) / 2}" stroke="${token('--ink')}" stroke-width="1" />`);
    }
    for (const [t, len] of ins) out.push(`<rect x="${(x(t) - 1).toFixed(1)}" y="${y - 0.5}" width="2" height="${rowH}" fill="${CATEGORIES[4]}"><title>${len} inserted</title></rect>`);
  });
  out.push('</svg>');
  const svg = el('div');
  svg.innerHTML = out.join('');
  fig.append(svg, el('p', { class: 'note' }, `${fmt.count(w.total)} reads over ${fmt.interval(w.begin, w.end)}${w.reads.length < w.total ? ` (the first ${fmt.count(w.reads.length)} drawn)` : ''}: mulberry forward, blue reverse, pale below MAPQ 20 and outlined in red at MAPQ 0; a coloured base differs from the reference, a line is a deletion, a violet tick an insertion. Above, the variants made (ochre) and those called (green when made, red when not).`));
}

// ---- mapping quality, as a job ---------------------------------------------------------------------------
function buildMapping() {
  const box = clear(view.mapping);
  const reads = el('select', {}, ...[100000, 400000, 1000000].map((v) => el('option', { value: v, selected: v === 400000 }, `${fmt.count(v)} reads`)));
  const shards = el('select', {}, ...[4, 8, 16, 32].map((v) => el('option', { value: v, selected: v === 16 }, `${v} shards`)));
  const out = el('div');
  const go = el('button', { class: 'button', onclick: async () => {
    go.disabled = true;
    clear(out).append(el('p', { class: 'wait' }, 'Cutting the reads into shards, simulating and mapping each on a worker…'));
    try {
      const r = await api.job('reads.map', { ...params(), diploid: true, reads: Number(reads.value), shards: Number(shards.value) });
      drawMapping(out, r);
    } catch (e) { clear(out).append(el('p', { class: 'error' }, e.message)); }
    go.disabled = false;
  } }, 'Map them');
  box.append(el('p', { class: 'note' }, 'The reads of the diploid, simulated and mapped in shards: each shard is two numbers, its first read and their count, so any worker, here or on another machine, makes the same reads; each worker builds its index of the genome once. Every read knows where it came from, so each mapping is right or wrong; sorting them by mapping quality gives Li\'s curve: how many reads are kept, and how many of those are wrong, at each threshold.'),
    el('div', { class: 'row' }, reads, shards, go), out);
}

function drawMapping(box, r) {
  clear(box);
  // Li's curve: thresholds from 60 down; x the error rate among the reads kept, y the fraction kept.
  const pts = r.curve.filter(([, mapped, wrong]) => mapped > 0).map(([q, mapped, wrong, frac, err]) => ({ q, frac, err: Math.max(err, 0.5 / Math.max(1, mapped)) }));
  const xs = pts.map((p) => p.err);
  const c = el('div');
  c.innerHTML = chart({ width: 460, height: 300, x: [Math.min(...xs) / 2, Math.max(...xs) * 2], y: [Math.min(...pts.map((p) => p.frac)) - 0.01, 1.002], log: 'x',
    xlabel: 'error rate among the reads kept', ylabel: 'fraction of reads kept',
    series: [{ points: pts.map((p) => [p.err, p.frac]), colour: token('--accent'), dots: true }],
    marks: pts.filter((p) => [60, 20, 10, 1, 0].includes(p.q)).map((p) => ({ x: p.err, y: p.frac, text: `≥${p.q}` })) });
  const hist = el('div');
  const by = r.by_mapq, maxN = Math.max(...by);
  hist.innerHTML = chart({ width: 420, height: 300, x: [-0.5, 60.5], y: [0.5, maxN * 2], log: 'y', xlabel: 'mapping quality', ylabel: 'reads',
    bars: by.map((n, q) => ({ x0: q - 0.5, x1: q + 0.5, y: Math.max(n, 0.5), colour: r.wrong_by_mapq[q] ? token('--bad') : token('--accent'), gap: 0.3, title: `MAPQ ${q}: ${fmt.count(n)} reads, ${fmt.count(r.wrong_by_mapq[q])} wrong` })) });
  const at = (q) => r.curve.find((x) => x[0] === q);
  const q20 = at(20), q0 = at(0);
  const job = r._job;
  box.append(el('p', { class: 'note' }, `${fmt.count(r.reads)} reads, ${fmt.count(r.unmapped)} unmapped. At MAPQ ≥ 20, ${fmt.percent(q20[3], 2)} of the reads are kept and ${fmt.count(q20[2])} of them are wrong (${fmt.number(q20[4], 3)}); with every read kept, ${fmt.percent(q0[4], 2)} are wrong, nearly all at MAPQ 0, in the genome's repeats, where the mapper said it could not tell. The mapping took ${fmt.seconds(r.map_seconds)} of engine time in all (${fmt.count(Math.round(r.reads / Math.max(1e-9, r.map_seconds)))} reads a second for a worker's threads), the job ${fmt.seconds(job.seconds)}.`),
    el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '24px' } }, c, hist));
  const workers = [...new Set(job.trace.map((t) => t.worker))].sort();
  const colour = (wk) => CATEGORIES[workers.indexOf(wk) % CATEGORIES.length];
  const bars = job.trace.map((t) => ({ row: workers.indexOf(t.worker), start: t.start, end: t.end, colour: colour(t.worker), opacity: 0.55 + 0.4 * (t.shard % 2 === 0), title: `shard ${t.shard}: ${fmt.count(t.size)} reads, ${fmt.seconds(t.engine_seconds)} in the engine` }));
  const tl = el('div');
  tl.innerHTML = timeline({ width: Math.min(980, view.section.clientWidth - 70), rows: workers, bars, seconds: job.seconds });
  box.append(el('p', { class: 'eyebrow' }, 'the shards on the workers, in seconds from the plan (the first shard on each worker builds its index)'), tl);
}

// ---- as files --------------------------------------------------------------------------------------------
function drawFiles() {
  const r = view.result, box = clear(view.files);
  const save = (name, text, type) => el('button', { class: 'button small', onclick: () => {
    const a = el('a', { href: URL.createObjectURL(new Blob([text], { type })), download: name });
    document.body.append(a); a.click(); a.remove();
  } }, `save ${name}`);
  box.append(el('p', { class: 'note' }, 'The first records of each, as the domain writes them: SAM 1.6 (Li et al. 2009) with the tags NM, AS, tp and s1; VCF 4.2 with the genotype, its quality, the depth and the depth of each allele.'),
    el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '18px' } },
      el('div', { style: { flex: 1, minWidth: 0 } }, el('p', { class: 'eyebrow' }, 'SAM'), el('pre', { class: 'file' }, r.sam), save('reads.sam', r.sam, 'text/plain')),
      el('div', { style: { flex: 1, minWidth: 0 } }, el('p', { class: 'eyebrow' }, 'VCF'), el('pre', { class: 'file' }, r.vcf), save('calls.vcf', r.vcf, 'text/plain'))));
}
