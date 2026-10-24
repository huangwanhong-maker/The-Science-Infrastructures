// The Reads view: a diploid made from a genome (SNVs and short indels, with their truth), read
// at a coverage, mapped back and called. Its parts: the calls against the truth, type by type; a
// stretch of the genome read by read (the pileup), with the calls and the truth above it; how
// well the mapper places reads at each mapping quality, measured on a job spread over the
// workers (Li 2018, Fig. 1); and the reads and calls as SAM and VCF. The instruments stand in the
// bench's column and the results on its stage; the values are kept in the address.

import { el, clear, esc } from '../dom.js';
import * as api from '../api.js';
import { chart, timeline } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, BASE_COLOURS, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, iconBtn, seg, select, fld, card, metrics, callout, empty, table, download, copy, working, busy, onEnter } from '../ui/kit.js';

const GENOMES = [
  ['ncbi:NC_001416.1', 'phage lambda (48.5 kb)'],
  ['sim:repeats-200k', 'simulated, 200 kb with repeats'],
  ['ncbi:NC_000913.3', 'Escherichia coli (4.6 Mb, half a minute)'],
];
const KINDS = [['short', 'short, 150 b'], ['long', 'long, 8 kb at 8%']];
const COVERAGES = [10, 20, 30, 60];
const SNVS = [0.0005, 0.001, 0.002, 0.005];
const INDELS = [0.00005, 0.0001, 0.0003];
const WIDTHS = [[120, '120 b'], [240, '240 b'], [1000, '1 kb'], [4000, '4 kb']];
const MAP_READS = [[100000, '100 k'], [400000, '400 k'], [1000000, '1 M']];
const SHARDS = [4, 8, 16, 32];

// The values kept in the address: [field, key, read].
const oneOf = (list) => (v) => (list.includes(v) ? v : undefined);
const numberOf = (list) => (v) => (list.includes(Number(v)) ? Number(v) : undefined);
const HASH = [
  ['genome', 'from', oneOf(GENOMES.map(([g]) => g))], ['kind', 'kind', oneOf(KINDS.map(([k]) => k))], ['coverage', 'cov', numberOf(COVERAGES)],
  ['snv', 'snv', numberOf(SNVS)], ['indel', 'indel', numberOf(INDELS)], ['seed', 'seed', (v) => (/^\d{1,6}$/.test(v) ? Number(v) : undefined)],
];

let view = null;
let keysBound = false;

export function mount(section, ctx) {
  view = { ctx, section, p: { genome: GENOMES[0][0], kind: 'short', coverage: 30, snv: 0.001, indel: 0.0001, seed: 5 }, m: { reads: 400000, shards: 16 },
           result: null, ran: null, at: null, wantAt: null, width: 240, file: 'sam', running: false, mapping: false, mapRun: null, browseToken: 0 };
  adopt(ctx.params);
  view.controls = el('div', { class: 'controls' });
  view.state = el('div', { class: 'results-head' });
  view.callBox = el('div', { class: 'stack', style: { gap: '16px' } });
  view.filesBox = el('div');
  view.mapCard = buildMapping();
  section.append(
    el('p', { class: 'lede' }, 'A read is placed on its genome in three steps, each paid for only when the one before has narrowed the choice: its minimizers, looked up in an index, give exact anchors; a dynamic programme over the anchors finds the chains that run in the same order on read and genome; only the best chains are aligned base by base. Whether the place is right is a posterior, the mapping quality: a read in a repeat has two places equally good and a quality of zero. At each position of a diploid the reads then vote, each with its base\'s error probability, among three genotypes; an indel must first be written in one canonical way, for in a repeat it has many.'),
    el('div', { class: 'bench' }, view.controls, el('div', { class: 'stage' }, view.state, view.callBox, view.mapCard, view.filesBox)));
  buildControls();
  if (!keysBound) { keysBound = true; window.addEventListener('keydown', onKey); }
}

export function show(ctx) {
  view.ctx = ctx;
  const changed = adopt(ctx.params);
  if (changed) buildControls();
  if (view.running) return;
  if (!view.result || changed) { call(); return; }
  if (view.wantAt) { view.at = view.wantAt; view.wantAt = null; browse(); return; }
  remember();
}

// The address's values, taken where they are valid; true when the diploid differs.
function adopt(params) {
  let changed = false;
  for (const [f, k, read] of HASH) {
    if (params[k] === undefined) continue;
    const v = read(params[k]);
    if (v !== undefined && v !== view.p[f]) { view.p[f] = v; changed = true; }
  }
  const w = Number(params.w);
  if (WIDTHS.some(([x]) => x === w)) view.width = w;
  const at = /^(\d+):(\d+)$/.exec(params.at || '');
  if (at && (!view.at || view.at.contig !== Number(at[1]) || view.at.position !== Number(at[2]))) view.wantAt = { contig: Number(at[1]), position: Number(at[2]) };
  return changed;
}

function remember() {
  if (view.ctx.view !== 'reads') return;
  const p = view.p;
  view.ctx.setParams({ from: p.genome, kind: p.kind, cov: p.coverage, snv: p.snv, indel: p.indel, seed: p.seed, w: view.width, at: view.at ? `${view.at.contig}:${view.at.position}` : null });
}

function params() {
  const p = view.p;
  return { genome: p.genome, kind: p.kind, coverage: p.coverage, snv: p.snv, indel: p.indel, seed: p.seed };
}

const swatch = (markup, w = 18) => { const s = el('span', { class: 'sw' }); s.innerHTML = `<svg width="${w}" height="10" viewBox="0 0 ${w} 10">${markup}</svg>`; return s; };
const key = (markup, label, w) => el('span', { class: 'k' }, swatch(markup, w), label);
const genomeLabel = (g) => GENOMES.find(([x]) => x === g)?.[1].replace(/ \(.*\)$/, '') || g;

// ---- the instruments -------------------------------------------------------------------------------------
function buildControls() {
  const p = view.p, m = view.m;
  view.runBtn = btn('Read, map and call', { icon: 'play', kind: 'primary', onClick: () => call() });
  view.againBtn = btn('Another', { icon: 'refresh', title: 'another diploid of the same genome: the next seed', onClick: () => { p.seed++; call(); } });
  if (view.running) { busy(view.runBtn, true, 'Reading…'); view.againBtn.disabled = true; }
  const diploid = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'The diploid'),
      fld('made from', select(GENOMES, p.genome, (v) => { p.genome = v; })),
      el('div', { class: 'ctl-row' },
        fld('SNVs', select(SNVS.map((v) => [v, `${fmt.number(v * 1000, 2)} per kb`]), p.snv, (v) => { p.snv = Number(v); })),
        fld('indels', select(INDELS.map((v) => [v, `${fmt.number(v * 1000, 2)} per kb`]), p.indel, (v) => { p.indel = Number(v); })))),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Its reads'),
      fld('kind', seg(KINDS, p.kind, (v) => { p.kind = v; }, { block: true })),
      fld('coverage', seg(COVERAGES.map((v) => [v, `${v}×`]), p.coverage, (v) => { p.coverage = v; }, { block: true }))),
    el('div', { class: 'run-bar' }, view.runBtn, view.againBtn));
  onEnter(diploid, call);
  view.mapBtn = btn('Map them', { icon: 'play', kind: 'primary', onClick: mapJob });
  if (view.mapping) busy(view.mapBtn, true, 'Mapping…');
  const job = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Mapping quality, as a job'),
      fld('reads', seg(MAP_READS, m.reads, (v) => { m.reads = v; }, { block: true })),
      fld('shards', seg(SHARDS.map((v) => [v, String(v)]), m.shards, (v) => { m.shards = v; }, { block: true }), { hint: 'The diploid above, read many times over; each shard is simulated and mapped on a worker.' })),
    el('div', { class: 'run-bar' }, view.mapBtn));
  onEnter(job, mapJob);
  const tool = (label, id, values) => el('button', { class: 'link', onclick: () => view.ctx.openTool(id, values()) }, label);
  const genome = () => (view.ctx.catalogue.some((g) => g.id === p.genome && g.fetched) ? { genome: p.genome } : {});
  const own = el('div', { class: 'toolbox-hint' }, icon('toolbox', 'sm'), el('span', {},
    'These reads are simulated, so each carries its truth. With your own reads, ', tool('map them', 'reads-map', () => ({ ...genome(), kind: p.kind })),
    ' or ', tool('call their variants', 'reads-call', genome), ' in the Toolbox.'));
  clear(view.controls).append(diploid, job, own);
}

// ---- reading, mapping and calling ---------------------------------------------------------------------------
async function call() {
  if (view.running) return;
  view.running = true;
  const p = params();
  busy(view.runBtn, true, 'Reading…');
  view.againBtn.disabled = true;
  view.state.replaceChildren(working('Making the diploid, reading it, mapping the reads and calling the variants…'));
  if (view.result) view.callBox.classList.add('stale');
  else view.callBox.replaceChildren(el('div', { class: 'skeleton', style: { height: '74px' } }), el('div', { class: 'skeleton', style: { height: '300px' } }), el('div', { class: 'skeleton', style: { height: '420px' } }));
  remember();
  try {
    const r = await api.run('reads.call', p);
    view.result = r;
    view.ran = { p, at: new Date() };
    // Open the browser where the address says, or on the first call that is an indel, or the first call.
    const first = r.calls.find((c) => c.indel) || r.calls[0];
    view.at = view.wantAt || (first ? { contig: first.contig, position: first.position } : { contig: 0, position: 1000 });
    view.wantAt = null;
    drawResults();
    browse();
  } catch (e) {
    view.state.replaceChildren(callout(`The diploid could not be read: ${e.message}`, 'bad'));
    if (!view.result) clear(view.callBox);
  } finally {
    view.running = false;
    view.callBox.classList.remove('stale');
    busy(view.runBtn, false);
    view.againBtn.disabled = false;
  }
}

function drawResults() {
  const r = view.result, q = view.ran.p, all = r.types.all;
  view.state.replaceChildren(el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Result'),
    el('span', { class: 'when' }, `${view.ran.at.toLocaleTimeString()} · ${genomeLabel(q.genome)}, ${q.kind} reads at ${q.coverage}×, seed ${q.seed}`));
  clear(view.callBox).append(...[
    metrics([
      [fmt.count(r.reads), `reads, ${fmt.length(r.mean_length)} on average`],
      [fmt.percent(r.correct / Math.max(1, r.reads), 2), 'placed where they came from'],
      [`${fmt.number(r.mean_depth, 3)}×`, 'mean depth'],
      [`${fmt.count(r.called)} / ${fmt.count(r.truth)}`, 'variants called / made'],
      [fmt.number(all.f1, 4), 'F₁ of all the calls', all.f1 >= 0.99 ? 'good' : ''],
      [fmt.percent(all.genotype, 1), 'of the true calls with the right genotype'],
    ]),
    q.kind === 'long' ? callout('Long reads are mapped well, but the caller is a caller of short reads: its genotype likelihoods (Li 2011) take each base\'s error as a substitution of known probability, and errors of 8%, mostly insertions and deletions, break that model, so most of its calls from long reads are false. Calling from long reads needs a model of their errors (a pair HMM over the read, or a consensus of the reads first), which is not yet written.', 'warn') : null,
    callsCard(r),
    stretchCard(),
  ].filter(Boolean));
  clear(view.filesBox).append(filesCard(r));
}

// ---- the calls against the truth ------------------------------------------------------------------------
function callsCard(r) {
  const s = r.seconds;
  const rows = [['snv', 'SNVs'], ['insertion', 'insertions'], ['deletion', 'deletions'], ['all', 'all']].map(([k, label]) => ({ label, ...r.types[k] }));
  const f = (d) => (x) => fmt.number(x, d);
  const t = table([{ key: 'label', label: 'variants' }, { key: 'tp', label: 'true', num: true, format: fmt.count }, { key: 'fp', label: 'false', num: true, format: fmt.count }, { key: 'fn', label: 'missed', num: true, format: fmt.count },
    { key: 'precision', label: 'precision', num: true, format: f(4) }, { key: 'recall', label: 'recall', num: true, format: f(4) }, { key: 'f1', label: 'F₁', num: true, format: f(4) }, { key: 'genotype', label: 'genotypes', num: true, format: f(4) }], rows);
  // Where the time went.
  const parts = [['simulate', 'reads simulated'], ['index', 'the index'], ['map', 'mapping'], ['call', 'calling'], ['compare', 'comparing']];
  const total = parts.reduce((a, [k]) => a + (s[k] || 0), 0);
  const bar = el('div', { class: 'stackbar' }, ...parts.map(([k, label], i) => el('span', { style: { width: `${(100 * (s[k] || 0)) / total}%`, background: CATEGORIES[i] }, title: `${label}: ${fmt.seconds(s[k] || 0)}` })));
  const legend = el('div', { class: 'key-row flush' }, ...parts.map(([k, label], i) => el('span', { class: 'k' }, el('span', { class: 'swatch', style: { background: CATEGORIES[i] } }), `${label} ${fmt.seconds(s[k] || 0)}`)));
  return card({ title: 'The calls against the truth', sub: `${fmt.count(r.truth)} variants made in ${fmt.length(r.length)}, ${fmt.count(r.called)} called`, pad: true, body: [
    t,
    el('p', { class: 'note' }, 'A call is true when, written in canonical form (left-aligned, parsimonious: Tan, Abecasis and Kang 2015), it is one of the variants made; a genotype agrees when it has their number of copies.'),
    el('p', { class: 'eyebrow', style: { marginTop: '16px' } }, `where the time went: ${fmt.seconds(total)} on ${r.threads} threads`), bar, legend] });
}

// ---- the pileup ------------------------------------------------------------------------------------------
function stretchCard() {
  view.bar = el('div', { class: 'rd-bar' });
  view.fig = el('div', { class: 'pileup' });
  view.count = el('span', { class: 'k count' });
  const accent = token('--accent');
  const legend = el('div', { class: 'key-row' }, view.count,
    key(`<rect y="1" width="18" height="8" fill="${accent}" fill-opacity="0.3"/>`, 'forward'),
    key(`<rect y="1" width="18" height="8" fill="${token('--C')}" fill-opacity="0.3"/>`, 'reverse'),
    key(`<rect y="1" width="18" height="8" fill="${accent}" fill-opacity="0.1"/>`, 'below MAPQ 20'),
    key(`<rect x="0.5" y="1.5" width="17" height="7" fill="${accent}" fill-opacity="0.08" stroke="${token('--bad')}" stroke-width="0.8"/>`, 'MAPQ 0'),
    key(['A', 'C', 'G', 'T'].map((b, i) => `<rect x="${i * 4.5}" y="1" width="3.5" height="8" fill="${token(`--${b}`)}"/>`).join(''), 'a base unlike the reference'),
    key(`<line x1="1" y1="5" x2="17" y2="5" stroke="${token('--ink')}"/>`, 'deletion'),
    key(`<rect x="8" y="0" width="2" height="10" fill="${CATEGORIES[4]}"/>`, 'insertion', 10),
    key(`<path d="M2,9L12,9L7,1Z" fill="${token('--G')}"/>`, 'made', 14),
    key(`<path d="M7,0L12,5L7,10L2,5Z" fill="${token('--good')}"/>`, 'called, and made', 14),
    key(`<path d="M7,0L12,5L7,10L2,5Z" fill="${token('--bad')}"/>`, 'called, not made', 14),
    key(`<path d="M7,0L12,5L7,10L2,5Z" fill="${token('--ink-faint')}"/>`, 'filtered', 14));
  return card({ title: 'A stretch, read by read', sub: 'the variants made and called above the reads; ← and → step through the calls', body: [view.bar, view.fig, legend] });
}

// The call at the stretch's centre, or the first after it.
function callIndex() {
  const at = view.at;
  return view.result.calls.findIndex((c) => c.contig > at.contig || (c.contig === at.contig && c.position >= at.position));
}

function goCall(k) {
  const calls = view.result.calls;
  const c = calls[(k + calls.length) % calls.length];
  if (c) { view.at = { contig: c.contig, position: c.position }; browse(); }
}

function step(d) {
  const calls = view.result?.calls;
  if (!calls?.length) return;
  let k = callIndex();
  if (k < 0) k = calls.length;
  const here = calls[k] && calls[k].contig === view.at.contig && calls[k].position === view.at.position;
  goCall(d > 0 ? (here ? k + 1 : k) : k - 1);
}

function drawBar() {
  const r = view.result, calls = r.calls, index = callIndex();
  const prev = iconBtn('chevron', { title: 'the call before (←)', size: 'sm', onClick: () => step(-1) });
  prev.classList.add('flip');
  const next = iconBtn('chevron', { title: 'the next call (→)', size: 'sm', onClick: () => step(1) });
  const pick = select(calls.slice(0, 2000).map((c, k) => [k, `${fmt.position(c.position)} ${c.ref}→${c.alts.join(',')} ${c.genotype}${c.true === false ? ' (false)' : ''}`]), index, (v) => goCall(Number(v)), { title: 'the calls, in order' });
  if (!calls.length) pick.disabled = true;
  const missed = select([{ value: '', label: r.missed.length ? `${r.missed.length} missed…` : 'none missed', disabled: true }, ...r.missed.map((v, k) => [k, `${fmt.position(v.position)} ${v.ref}→${v.alt}`])], '',
    (v) => { const x = r.missed[Number(v)]; if (x) { view.at = { contig: x.contig, position: x.position }; browse(); } }, { title: 'the variants made and not called' });
  if (!r.missed.length) missed.disabled = true;
  // What is known of the call at the centre.
  const c = calls[index], here = c && c.contig === view.at.contig && c.position === view.at.position;
  const facts = here ? [c.true === false ? el('span', { class: 'badge bad' }, 'not among those made') : c.passed ? el('span', { class: 'badge good' }, 'made') : el('span', { class: 'badge' }, 'filtered'),
    el('span', { class: 'facts' }, `QUAL ${c.qual} · GQ ${c.gq} · depth ${c.depth} · AD ${c.ad.join(',')}${c.filters.length ? ` · ${c.filters.join(', ')}` : ''}`)]
    : r.missed.some((v) => v.contig === view.at.contig && v.position === view.at.position) ? [el('span', { class: 'badge warn' }, 'made, not called')] : [];
  view.bar.replaceChildren(el('div', { class: 'nav' }, prev, pick, next), missed, ...facts,
    seg(WIDTHS, view.width, (v) => { view.width = v; browse(); }));
}

async function browse() {
  if (!view.result || !view.fig) return;
  const at = view.at, half = Math.round(view.width / 2);
  const begin = Math.max(0, at.position - half), end = begin + view.width;
  drawBar();
  remember();
  if (view.fig.querySelector('svg')) view.fig.classList.add('stale');
  else view.fig.replaceChildren(working('Reading the stretch…'));
  const t = ++view.browseToken;
  let w;
  try { w = await api.run('reads.window', { ...view.ran.p, contig: at.contig, begin, end, max_reads: 1500 }); } catch (e) {
    if (t === view.browseToken) { view.fig.classList.remove('stale'); clear(view.fig).append(callout(e.message, 'bad')); }
    return;
  }
  if (t !== view.browseToken) return;   // a later stretch was asked for
  view.fig.classList.remove('stale');
  drawPileup(view.fig, w);
  view.count.textContent = `${fmt.count(w.total)} reads over ${fmt.interval(w.begin, w.end)}${w.reads.length < w.total ? ` (the first ${fmt.count(w.reads.length)} drawn)` : ''}`;
}

function drawPileup(fig, w) {
  clear(fig);
  const colours = BASE_COLOURS();
  const W = Math.min(1180, fig.clientWidth > 300 ? fig.clientWidth - 16 : 1000), left = 8, right = 8, span = w.end - w.begin;
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
  const step_ = fmt.tickStep(span, 8);
  for (let p = Math.ceil(w.begin / step_) * step_; p < w.end; p += step_) out.push(`<line x1="${x(p)}" y1="12" x2="${x(p)}" y2="17" stroke="${token('--ink-faint')}" /><text class="tick-label" x="${x(p)}" y="10" text-anchor="middle">${fmt.position(p)}</text>`);
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
  fig.append(svg);
}

// ---- mapping quality, as a job ---------------------------------------------------------------------------
function buildMapping() {
  view.mapState = el('div');
  view.mapRes = el('div', { class: 'stack', style: { gap: '14px', marginTop: '12px' } },
    empty({ icon: 'activity', title: 'Not measured yet', text: 'Run the job from the column at the left: the diploid\'s reads, many times over, mapped in shards on every worker there is.', action: btn('Run the job', { icon: 'play', onClick: mapJob }) }));
  view.mapSub = el('span');
  return card({ title: 'Mapping quality, measured across the workers', sub: 'how many reads are kept, and how many of those are wrong, at each threshold', pad: true, body: [
    el('p', { class: 'note', style: { marginTop: 0 } }, 'The reads of the diploid, simulated and mapped in shards: each shard is two numbers, its first read and their count, so any worker, here or on another machine, makes the same reads; each worker builds its index of the genome once. Every read knows where it came from, so each mapping is right or wrong; sorting them by mapping quality gives Li\'s curve: how many reads are kept, and how many of those are wrong, at each threshold.'),
    view.mapState, view.mapRes] });
}

async function mapJob() {
  if (view.mapping) return;
  view.mapping = true;
  const p = { ...params(), diploid: true, reads: view.m.reads, shards: view.m.shards };
  busy(view.mapBtn, true, 'Mapping…');
  view.mapState.replaceChildren(working(`Cutting ${fmt.count(p.reads)} reads into ${p.shards} shards, simulating and mapping each on a worker…`));
  if (view.mapRun) view.mapRes.classList.add('stale');
  else view.mapRes.replaceChildren(el('div', { class: 'skeleton', style: { height: '74px' } }), el('div', { class: 'skeleton', style: { height: '300px' } }));
  const r0 = view.mapCard.getBoundingClientRect();
  if (r0.top > window.innerHeight - 120 || r0.bottom < 60) view.mapCard.scrollIntoView({ behavior: 'smooth', block: 'start' });
  try {
    const r = await api.job('reads.map', p);
    view.mapRun = { r, p, at: new Date() };
    clear(view.mapState);
    drawMapping();
  } catch (e) { view.mapState.replaceChildren(callout(`The job failed: ${e.message}`, 'bad')); }
  finally { view.mapping = false; view.mapRes.classList.remove('stale'); busy(view.mapBtn, false); }
}

function drawMapping() {
  const { r, p, at: when } = view.mapRun, box = clear(view.mapRes);
  // Li's curve: thresholds from 60 down; x the error rate among the reads kept, y the fraction kept.
  const pts = r.curve.filter(([, mapped, wrong]) => mapped > 0).map(([q, mapped, wrong, frac, err]) => ({ q, frac, err: Math.max(err, 0.5 / Math.max(1, mapped)) }));
  const xs = pts.map((x) => x.err);
  const c = el('div');
  c.innerHTML = chart({ width: 460, height: 300, x: [Math.min(...xs) / 2, Math.max(...xs) * 2], y: [Math.min(...pts.map((x) => x.frac)) - 0.01, 1.002], log: 'x',
    xlabel: 'error rate among the reads kept', ylabel: 'fraction of reads kept',
    series: [{ points: pts.map((x) => [x.err, x.frac]), colour: token('--accent'), dots: true }],
    marks: pts.filter((x) => [60, 20, 10, 1, 0].includes(x.q)).map((x) => ({ x: x.err, y: x.frac, text: `≥${x.q}` })) });
  const hist = el('div');
  const by = r.by_mapq, maxN = Math.max(...by);
  hist.innerHTML = chart({ width: 420, height: 300, x: [-0.5, 60.5], y: [0.5, maxN * 2], log: 'y', xlabel: 'mapping quality', ylabel: 'reads',
    bars: by.map((n, q) => ({ x0: q - 0.5, x1: q + 0.5, y: Math.max(n, 0.5), colour: r.wrong_by_mapq[q] ? token('--bad') : token('--accent'), gap: 0.3, title: `MAPQ ${q}: ${fmt.count(n)} reads, ${fmt.count(r.wrong_by_mapq[q])} wrong` })) });
  const atQ = (q) => r.curve.find((x) => x[0] === q);
  const q20 = atQ(20), q0 = atQ(0);
  const job = r._job;
  box.append(
    el('div', { class: 'results-head', style: { margin: 0 } }, el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Measured'),
      el('span', { class: 'when' }, `${when.toLocaleTimeString()} · ${genomeLabel(p.genome)}, ${p.kind} reads at ${p.coverage}×, seed ${p.seed} · ${p.shards} shards`)),
    metrics([
      [fmt.count(r.reads), 'reads mapped'],
      [fmt.percent(q20[3], 2), 'kept at MAPQ ≥ 20'],
      [fmt.count(q20[2]), 'of those, wrong'],
      [fmt.percent(q0[4], 2), 'wrong, with every read kept'],
      [fmt.count(Math.round(r.reads / Math.max(1e-9, r.map_seconds))), 'reads a second, a worker'],
      [fmt.seconds(job.seconds), 'the job, start to end'],
    ]),
    el('p', { class: 'note', style: { margin: 0 } }, `${fmt.count(r.reads)} reads, ${fmt.count(r.unmapped)} unmapped. At MAPQ ≥ 20, ${fmt.percent(q20[3], 2)} of the reads are kept and ${fmt.count(q20[2])} of them are wrong (${fmt.number(q20[4], 3)}); with every read kept, ${fmt.percent(q0[4], 2)} are wrong, nearly all at MAPQ 0, in the genome's repeats, where the mapper said it could not tell. The mapping took ${fmt.seconds(r.map_seconds)} of engine time in all (${fmt.count(Math.round(r.reads / Math.max(1e-9, r.map_seconds)))} reads a second for a worker's threads), the job ${fmt.seconds(job.seconds)}.`),
    el('div', { class: 'rd-charts' }, el('div', {}, el('p', { class: 'eyebrow' }, 'Li\'s curve, thresholds from 60 down'), c), el('div', {}, el('p', { class: 'eyebrow' }, 'reads by mapping quality, wrong ones in red'), hist)));
  const workers = [...new Set(job.trace.map((t) => t.worker))].sort();
  const colour = (wk) => CATEGORIES[workers.indexOf(wk) % CATEGORIES.length];
  const bars = job.trace.map((t) => ({ row: workers.indexOf(t.worker), start: t.start, end: t.end, colour: colour(t.worker), opacity: 0.55 + 0.4 * (t.shard % 2 === 0), title: `shard ${t.shard}: ${fmt.count(t.size)} reads, ${fmt.seconds(t.engine_seconds)} in the engine` }));
  const tl = el('div', { style: { overflowX: 'auto' } });
  tl.innerHTML = timeline({ width: Math.min(980, Math.max(480, box.clientWidth - 8)), rows: workers, bars, seconds: job.seconds });
  box.append(el('div', {}, el('p', { class: 'eyebrow' }, 'the shards on the workers, in seconds from the plan (the first shard on each worker builds its index)'), tl));
}

// ---- as files --------------------------------------------------------------------------------------------
function filesCard(r) {
  const pre = el('pre', { class: 'file', style: { border: 0, borderRadius: 0, borderTop: '1px solid var(--hair)', maxHeight: '340px' } });
  const text = () => (view.file === 'sam' ? r.sam : r.vcf);
  const name = () => (view.file === 'sam' ? 'reads.sam' : 'calls.vcf');
  const pick = (k) => { view.file = k; pre.textContent = text(); };
  pick(view.file);
  return card({ title: 'As files', sub: 'the first records of each, as the domain writes them',
    actions: [seg([['sam', 'SAM'], ['vcf', 'VCF']], view.file, pick), btn('Save', { icon: 'download', size: 'sm', onClick: () => download(name(), text()) }), iconBtn('copy', { title: 'copy', size: 'sm', onClick: () => copy(text()) })],
    body: [el('p', { class: 'note', style: { margin: 0, padding: '10px 16px' } }, 'SAM 1.6 (Li et al. 2009) with the tags NM, AS, tp and s1; VCF 4.2 with the genotype, its quality, the depth and the depth of each allele.'), pre] });
}

// ---- keys ------------------------------------------------------------------------------------------------
function onKey(e) {
  if (!view || view.section.hidden || !view.result || !view.fig) return;
  const a = document.activeElement;
  if (/^(INPUT|TEXTAREA|SELECT)$/.test(a?.tagName || '') || a?.isContentEditable || e.ctrlKey || e.metaKey || e.altKey || document.querySelector('.dialog-veil, .popover')) return;
  if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') { e.preventDefault(); step(e.key === 'ArrowLeft' ? -1 : 1); }
}
