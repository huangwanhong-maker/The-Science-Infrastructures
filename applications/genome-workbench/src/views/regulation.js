// The Regulation view: how much of a gene's expression its sequence predicts (genomics::codon and
// genomics::regulation, tyche's learning; ADR 0064). Three readings, each against measurement:
//
//   codons      a gene's codons, against those of the genome's most expressed genes (CAI) and its
//               tRNA genes (tAI), painted round the genome and set against measured protein
//               abundance (PaxDb);
//   promoters   σ70's two boxes and their spacing, learnt by EM from Escherichia coli's measured
//               transcription start sites, and the promoters they find along the genome;
//   learnt      expression learnt from measured promoters (the DREAM Challenge's yeast data): words of
//               k letters with an elastic net, and a convolutional network whose filters are read back
//               as motifs and matched to JASPAR's.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import { sector, angleOf } from '../draw/circle.js';
import { logoSvg } from '../draw/logo.js';
import * as fmt from '../format.js';
import { token, warm, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, seg, select, fld, card, metrics, callout, empty, table, working, busy } from '../ui/kit.js';

const GENOMES = [
  ['ncbi:NC_000913.3', 'Escherichia coli K-12', { abundance: 'data:paxdb-ecoli/511145-WHOLE_ORGANISM-integrated.txt' }],
  ['ncbi:NC_001416.1', 'phage lambda, on its host', { host: 'ncbi:NC_000913.3' }],
  ['ncbi:NC_001422.1', 'phage phiX174, on its host', { host: 'ncbi:NC_000913.3' }],
  ['ncbi:NC_012920.1', 'the human mitochondrion', { table: 2 }],
];
const TSS = { tss: 'data:ecoli-tss-cappable/ettwiller-2016-tss.gtf', genome_fasta: 'data:ecoli-tss-cappable/U00096.2.fasta' };
const DREAM = { data: 'data:dream-promoters/filtered_test_data_with_MAUDE_expression.txt', motifs: 'data:jaspar-fungi/JASPAR2024_CORE_fungi_non-redundant_pfms_meme.txt' };

const kept = {
  tab: 'codons', genome: GENOMES[0][0], index: 'cai',
  codons: null, codonsError: null, codonsBusy: false,
  sigma: null, sigmaError: null, sigmaBusy: false, region: 0, promoters: null,
  learn: null, learnError: null, learnBusy: false, epochs: 8, filters: 32,
};
let view = null;

export function mount(section, ctx) {
  view = { ctx, section };
  if (['codons', 'promoters', 'learnt'].includes(ctx.params.tab)) kept.tab = ctx.params.tab;
  view.tabs = el('div', { class: 'tabs', role: 'tablist' });
  for (const [v, l] of [['codons', 'Codons and abundance'], ['promoters', 'Promoters'], ['learnt', 'Expression learnt from sequence']]) {
    view.tabs.append(el('button', { role: 'tab', 'data-v': v, class: v === kept.tab ? 'on' : '', onclick: () => { kept.tab = v; ctx.setParams({ tab: v }); for (const b of view.tabs.children) b.classList.toggle('on', b.dataset.v === v); draw(); } }, l));
  }
  view.body = el('div');
  section.append(
    el('p', { class: 'lede' }, 'How much of a gene\'s expression is written in its sequence? In bacteria a highly expressed gene prefers the codons read by its genome\'s most abundant tRNAs, so its codons alone, compared with those of the ribosomal proteins (CAI) or weighed by the genome\'s tRNA genes (tAI), predict its protein\'s abundance. A promoter\'s strength is written in two short boxes and their spacing, read by σ70; their model can be learnt from measured start sites alone. And expression can be learnt from many measured promoters at once, by words of a few letters or by a convolutional network, whose filters are then read back as motifs. Every reading here is set against measurement.'),
    view.tabs, view.body);
  draw();
}

export function show(ctx) { view.ctx = ctx; }

function draw() {
  clear(view.body);
  if (kept.tab === 'codons') drawCodons();
  else if (kept.tab === 'promoters') drawPromoters();
  else drawLearnt();
}

// ---- codons --------------------------------------------------------------------------------------------
function drawCodons() {
  const runBtn = btn('Read the codons', { icon: 'play', kind: 'primary', onClick: () => codons() });
  if (kept.codonsBusy) busy(runBtn, true, 'Reading…');
  const controls = el('div', { class: 'ctl-row', style: { alignItems: 'flex-end', gap: '12px', flexWrap: 'wrap' } },
    fld('genome', select(GENOMES.map(([v, l]) => [v, l]), kept.genome, (v) => { kept.genome = v; codons(); })),
    fld('painted by', seg([['cai', 'CAI'], ['tai', 'tAI'], ['nc_wright', 'Nc'], ['gc3s', 'GC3s']], kept.index, (v) => { kept.index = v; drawCodonsResult(); })),
    runBtn);
  view.codonStage = el('div', { class: 'stack', style: { gap: '16px', marginTop: '14px' } });
  view.body.append(card({ pad: true, body: controls }), view.codonStage);
  if (kept.codons && kept.codons.request === kept.genome) drawCodonsResult();
  else if (!kept.codonsBusy) codons();
  else view.codonStage.append(working('Reading every gene\'s codons…'));
}

async function codons() {
  if (kept.codonsBusy) return;
  kept.codonsBusy = true; kept.codonsError = null;
  const [id, , extra] = GENOMES.find(([v]) => v === kept.genome);
  if (view.codonStage) clear(view.codonStage).append(working('Reading every gene\'s codons, the genome\'s tRNA genes, and the S value of translational selection…'));
  try {
    const r = await api.run('codon.indices', { genome: id, reference: 'ribosomal', s_permutations: 199, table: 11, ...extra });
    r.request = id;
    kept.codons = r;
  } catch (e) {
    kept.codonsError = e.message;
  } finally {
    kept.codonsBusy = false;
    if (kept.tab === 'codons') drawCodonsResult();
  }
}

const col = (r, name) => r.columns.indexOf(name);

function drawCodonsResult() {
  const box = view.codonStage && clear(view.codonStage);
  if (!box) return;
  if (kept.codonsError) { box.append(callout(kept.codonsError, 'bad')); return; }
  const r = kept.codons;
  if (!r) return;
  const genes = r.genes;
  const ag = r.abundance?.agreement || {};
  box.append(metrics([
    [fmt.count(genes.length), 'genes read'],
    [String(r.reference?.genes ?? '–'), `reference genes (${r.reference?.kind || ''})`],
    [r.trna ? String(r.trna.genes) : '–', 'tRNA genes'],
    [r.s_value ? fmt.number(r.s_value.s, 3) : '–', r.s_value ? `S, translational selection (p ${fmt.number(r.s_value.p, 2)})` : 'S value'],
    ag.cai ? [fmt.number(ag.cai.spearman, 3), 'ρ of CAI with abundance'] : null,
    ag.tai ? [fmt.number(ag.tai.spearman, 3), 'ρ of tAI with abundance'] : null,
  ]));
  const row = el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(440px, 1fr))', gap: '16px' } });
  row.append(card({ title: 'Round the genome', sub: `each gene painted by its ${indexName(kept.index)}; the outer ring forward, the inner reverse; dots mark the genes predicted highly expressed (PHX)`, pad: true, body: ring(r) }));
  if (r.abundance) row.append(card({ title: 'Against measured abundance', sub: `${r.abundance.source || 'PaxDb'}; ${r.abundance.points.length} genes matched`, pad: true, body: abundancePlot(r) }));
  box.append(row);
  const row2 = el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(440px, 1fr))', gap: '16px' } });
  row2.append(card({ title: 'The effective number of codons', sub: 'each gene\'s Nc against its GC3s; the curve is Wright\'s expectation without selection (Wright 1990); genes far below it are biased by more than mutation', pad: true, body: ncPlot(r) }));
  row2.append(card({ title: 'Which codons the expressed genes prefer', sub: 'the genetic code as a table, each codon shaded by its relative adaptiveness w in the reference set (CAI\'s weights)', pad: true, body: codeTable(r) }));
  box.append(row2);
  if (Object.keys(ag).length) {
    const rows = Object.entries(ag).map(([k, v]) => ({ index: indexName(k), spearman: fmt.number(v.spearman, 3), pearson: fmt.number(v.pearson_log, 3), matched: v.matched })).sort((a, b) => Math.abs(b.spearman) - Math.abs(a.spearman));
    box.append(card({ title: 'Each index against abundance', sub: 'Spearman\'s ρ, and Pearson\'s r against the logarithm of abundance', pad: true,
      body: table([{ key: 'index', label: 'index' }, { key: 'spearman', label: 'ρ', num: true }, { key: 'pearson', label: 'r (log)', num: true }, { key: 'matched', label: 'genes', num: true }], rows) }));
  }
  if (r.trna?.note) box.append(callout(r.trna.note, 'info'));
}

const indexName = (k) => ({ cai: 'CAI', tai: 'tAI', fop: 'Fop', nc_wright: 'Nc (Wright)', nc_sun: 'Nc (Sun)', gc3s: 'GC3s', milc: 'MILC', milc_negated: '−MILC', e: 'E (Karlin)' }[k] || k);

function ring(r) {
  const genes = r.genes, L = r.length;
  const k = col(r, kept.index), s0 = col(r, 'start'), s1 = col(r, 'stop'), rev = col(r, 'reverse'), phx = col(r, 'phx'), name = col(r, 'gene'), prod = col(r, 'product');
  const vals = genes.map((g) => g[k]).filter(Number.isFinite).sort((a, b) => a - b);
  const lo = vals[Math.floor(vals.length * 0.02)] ?? 0, hi = vals[Math.floor(vals.length * 0.98)] ?? 1;
  const flip = kept.index === 'nc_wright';  // a low Nc means strong bias: drawn dark like a high CAI
  const shade = (v) => { const t = (v - lo) / (hi - lo || 1); return warm(Math.min(1, Math.max(0, flip ? 1 - t : t))); };
  const out = [`<circle r="300" fill="none" stroke="${token('--hair')}" />`];
  genes.forEach((g, i) => {
    const v = g[k];
    if (!Number.isFinite(v)) return;
    let a0 = angleOf(g[s0], L), a1 = angleOf(g[s1], L);
    if (a1 < a0) a1 += 2 * Math.PI;
    a1 = Math.max(a1, a0 + 0.004);
    const [r0, r1] = g[rev] ? [262, 296] : [304, 338];
    out.push(`<path data-i="${i}" d="${sector(r0, r1, a0, a1)}" fill="${shade(v)}" />`);
    if (g[phx]) { const a = (a0 + a1) / 2, rr = g[rev] ? 254 : 346; out.push(`<circle cx="${(rr * Math.cos(a)).toFixed(1)}" cy="${(rr * Math.sin(a)).toFixed(1)}" r="2" fill="${token('--accent-deep')}" />`); }
  });
  const ag = r.abundance?.agreement?.[kept.index];
  out.push(`<text text-anchor="middle" y="-14" style="font: 500 22px var(--serif); fill:${token('--ink')}">${esc(indexName(kept.index))}</text>`);
  out.push(`<text text-anchor="middle" y="12" style="font: italic 13px var(--serif); fill:${token('--ink-2')}">${esc(r.accession || '')}</text>`);
  if (ag) out.push(`<text text-anchor="middle" y="36" style="font: 13px var(--serif); fill:${token('--ink-2')}">ρ ${fmt.number(ag.spearman, 3)} with abundance</text>`);
  const fig = el('div', { style: { display: 'flex', justifyContent: 'center' } });
  fig.innerHTML = `<svg width="500" height="500" viewBox="-360 -360 720 720">${out.join('')}</svg>`;
  fig.addEventListener('mousemove', (e) => {
    const p = e.target.closest('path[data-i]');
    if (!p) { hover(null); return; }
    const g = genes[Number(p.dataset.i)];
    hover(`<b>${esc(g[name] || g[0])}</b> · ${esc(g[prod] || '')}<br>CAI ${fmt.number(g[col(r, 'cai')], 3)} · tAI ${fmt.number(g[col(r, 'tai')], 3)} · Nc ${fmt.number(g[col(r, 'nc_wright')], 3)}${g[phx] ? ' · PHX' : ''}`, e);
  });
  fig.addEventListener('mouseleave', () => hover(null));
  return fig;
}

function abundancePlot(r) {
  const pts = r.abundance.points;  // [gene, cai, tai, ppm]
  const which = kept.index === 'tai' ? 2 : 1;
  const xs = pts.map((p) => p[which]).filter(Number.isFinite);
  const series = [{ points: pts.filter((p) => p[3] > 0 && Number.isFinite(p[which])).map((p) => [p[which], p[3]]), colour: token('--accent'), line: false, dots: true }];
  const fig = el('div');
  fig.innerHTML = chart({ width: 440, height: 320, x: [Math.min(...xs), Math.max(...xs)], y: [0.01, Math.max(...pts.map((p) => p[3])) * 1.5], log: 'y', xlabel: which === 2 ? 'tAI' : 'CAI', ylabel: 'protein abundance (ppm)', series }).replace(/r="3"/g, 'r="1.8" fill-opacity="0.45"');
  return fig;
}

function ncPlot(r) {
  const g3 = col(r, 'gc3s'), nc = col(r, 'nc_wright');
  const pts = r.genes.map((g) => [g[g3], g[nc]]).filter(([a, b]) => Number.isFinite(a) && Number.isFinite(b));
  const curve = [];
  for (let s = 0; s <= 1.0001; s += 0.02) curve.push([s, 2 + s + 29 / (s * s + (1 - s) * (1 - s))]);
  const fig = el('div');
  fig.innerHTML = chart({ width: 440, height: 300, x: [0, 1], y: [20, 62], xlabel: 'GC3s', ylabel: 'Nc', series: [
    { points: pts, colour: token('--C'), line: false, dots: true },
    { points: curve, colour: token('--accent-deep'), width: 1.6, label: 'Wright\'s expectation' }] }).replace(/r="3"/g, 'r="1.6" fill-opacity="0.4"');
  return fig;
}

// The standard table of the code: first base by rows, second by columns, third within each cell.
function codeTable(r) {
  const w = new Map(r.weights.map((x) => [x[0], x]));
  const B = 'TCAG';
  const t = el('table', { class: 'data', style: { fontSize: '12px', borderCollapse: 'separate', borderSpacing: '2px' } });
  for (const a of B) for (const c of B) {
    const tr = el('tr');
    for (const b of B) {
      const codon = a + b + c, x = w.get(codon);
      const v = x ? x[2] : null;
      tr.append(el('td', { style: { background: Number.isFinite(v) ? warm(0.08 + 0.75 * v) : 'transparent', color: Number.isFinite(v) && v > 0.62 ? 'var(--on-accent)' : 'var(--ink)', padding: '3px 7px', borderRadius: '4px', whiteSpace: 'nowrap' },
        title: x ? `${codon} (${x[1]}): w ${fmt.number(v, 3)}${x[3] != null ? `, tAI's w ${fmt.number(x[3], 3)}` : ''}; ${x[4]} in all genes, ${x[5]} in the reference set` : codon },
        el('span', { class: 'mono', style: { color: 'inherit' } }, codon), ' ', el('b', { style: { fontWeight: 500 } }, x ? x[1] : '')));
    }
    t.append(tr);
  }
  return el('div', { style: { overflowX: 'auto' } }, t);
}

// ---- promoters ---------------------------------------------------------------------------------------------
function drawPromoters() {
  const runBtn = btn('Learn σ70 from start sites', { icon: 'play', kind: 'primary', onClick: () => sigma() });
  if (kept.sigmaBusy) busy(runBtn, true, 'Learning…');
  view.promStage = el('div', { class: 'stack', style: { gap: '16px', marginTop: '14px' } });
  view.body.append(card({ pad: true, body: [el('p', { class: 'hint', style: { margin: '0 0 10px' } }, 'From 2,000 of Escherichia coli\'s transcription start sites measured by Cappable-seq (Ettwiller et al. 2016): the windows upstream of each, searched by EM for two boxes at a variable distance, each box present or not.'), runBtn] }), view.promStage);
  if (kept.sigma) drawSigma();
  else if (!kept.sigmaBusy) sigma();
}

async function sigma() {
  if (kept.sigmaBusy) return;
  kept.sigmaBusy = true; kept.sigmaError = null;
  if (view.promStage) clear(view.promStage).append(working('Learning the two boxes and their spacing by EM from the start sites…'));
  try {
    kept.sigma = await api.run('regulation.sigma70', { ...TSS, sites: 2000 });
    await promoters();
  } catch (e) {
    kept.sigmaError = e.message;
  } finally {
    kept.sigmaBusy = false;
    if (kept.tab === 'promoters') drawSigma();
  }
}

async function promoters() {
  const begin = kept.region, end = begin + 20000;
  kept.promoters = await api.run('regulation.promoters', { ...TSS, sites: 2000, genome: 'ncbi:NC_000913.3', contig: 0, begin, end, threshold: 8 });
}

const rowsOf = (box) => (Array.isArray(box.probability?.[0]) ? box.probability : (box.logo || []).map((x) => { const s = x.heights.reduce((a, b) => a + b, 0) || 1; return x.heights.map((h) => h / s); }));

function drawSigma() {
  const box = view.promStage && clear(view.promStage);
  if (!box) return;
  if (kept.sigmaError) { box.append(callout(kept.sigmaError, 'bad')); return; }
  const s = kept.sigma;
  if (!s) return;
  const m35 = logoSvg({ rows: rowsOf(s.minus35), title: `−35 box · ${s.minus35.consensus}` });
  const m10 = logoSvg({ rows: rowsOf(s.minus10), title: `−10 box · ${s.minus10.consensus}` });
  const spacerBars = (s.spacers.p || []).map((p, i) => ({ x0: s.spacers.from + i - 0.4, x1: s.spacers.from + i + 0.4, y: p, colour: token('--accent-2') }));
  const spacer = chart({ width: 250, height: 150, x: [s.spacers.from - 1, s.spacers.from + s.spacers.p.length], y: [0, Math.max(...s.spacers.p) * 1.15], xlabel: 'spacer (bases)', bars: spacerBars });
  const fig = el('div', { style: { display: 'flex', alignItems: 'flex-end', gap: '6px', flexWrap: 'wrap' } });
  fig.innerHTML = `<div>${m35}</div><div>${spacer}</div><div>${m10}</div>`;
  const trace = (s.trace || []).map((v, i) => [i + 1, v]);
  const tr = el('div');
  tr.innerHTML = chart({ width: 300, height: 150, x: [1, Math.max(2, trace.length)], y: [Math.min(...trace.map((p) => p[1])), Math.max(...trace.map((p) => p[1]))], xlabel: 'EM iteration', ylabel: 'log likelihood', series: [{ points: trace, colour: token('--accent') }] });
  box.append(card({ title: 'σ70\'s promoter, learnt', sub: `two boxes with the spacing between them; ${fmt.percent(s.gamma, 0)} of the windows hold a promoter by the model (γ); ${s.iterations} iterations of EM on ${s.learnt?.windows ?? '?'} windows`, pad: true,
    body: el('div', { style: { display: 'flex', gap: '20px', flexWrap: 'wrap', alignItems: 'flex-end' } }, fig, tr) }));
  if (kept.promoters) box.append(promoterStrip(kept.promoters));
}

function promoterStrip(p) {
  const c = p.columns, hits = p.hits || [];
  const L = p.end - p.begin, W = 980, H = 120, m = 20;
  const X = (x) => m + ((x - p.begin) / L) * (W - 2 * m);
  const out = [`<line x1="${m}" y1="60" x2="${W - m}" y2="60" stroke="${token('--hair-3')}" />`];
  const bi = c.indexOf('bits'), st = c.indexOf('start'), rv = c.indexOf('reverse');
  const maxBits = Math.max(1, ...hits.map((h) => h[bi]));
  for (const h of hits) {
    const x = X(h[st]), up = !h[rv], len = 10 + 26 * (h[bi] / maxBits);
    const y = up ? 58 : 62;
    out.push(`<path d="M${x},${y} v${up ? -len : len} h${up ? 8 : -8}" fill="none" stroke="${token('--accent')}" stroke-width="1.5" stroke-opacity="${0.35 + 0.65 * h[bi] / maxBits}"><title>${esc(`${fmt.number(h[bi], 3)} bits · −35 ${h[c.indexOf('minus35')]} · −10 ${h[c.indexOf('minus10')]} · spacer ${h[c.indexOf('spacer')]} · +1 at ${fmt.position(h[st])}`)}</title></path>`);
  }
  for (let k = 0; k <= 4; k++) { const x = p.begin + (k * L) / 4; out.push(`<text class="tick-label" x="${X(x)}" y="${H - 6}" text-anchor="middle">${esc(fmt.number(x / 1000, 4))} kb</text>`); }
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${out.join('')}</svg>`;
  const prev = btn('', { icon: 'left', size: 'sm', title: 'the 20 kb before', onClick: async () => { kept.region = Math.max(0, kept.region - 20000); await promoters(); drawSigma(); } });
  const next = btn('', { icon: 'chevron', size: 'sm', title: 'the next 20 kb', onClick: async () => { kept.region += 20000; await promoters(); drawSigma(); } });
  return card({ title: 'Promoters along the genome', sub: `${hits.length} above ${p.threshold} bits in ${fmt.length(L)} of Escherichia coli K-12 from ${fmt.position(p.begin)}; arrows point the way they transcribe, taller where stronger`, actions: [prev, next], pad: true, body: fig });
}

// ---- learnt ----------------------------------------------------------------------------------------------------
function drawLearnt() {
  const runBtn = btn('Learn expression', { icon: 'play', kind: 'primary', onClick: () => learn() });
  if (kept.learnBusy) busy(runBtn, true, 'Learning…');
  const controls = el('div', { class: 'ctl-row', style: { alignItems: 'flex-end', gap: '12px', flexWrap: 'wrap' } },
    fld('epochs', seg([[4, '4'], [8, '8'], [16, '16']], kept.epochs, (v) => { kept.epochs = Number(v); })),
    fld('filters', seg([[16, '16'], [32, '32'], [64, '64']], kept.filters, (v) => { kept.filters = Number(v); })),
    runBtn);
  view.learnStage = el('div', { class: 'stack', style: { gap: '16px', marginTop: '14px' } });
  view.body.append(card({ pad: true, body: [el('p', { class: 'hint', style: { margin: '0 0 10px' } }, 'Yeast promoters of 80 random bases, each measured for expression (the DREAM Challenge, Rafi et al.): 8,500 to learn from, the rest held out by family, so that near-copies never straddle the split. About a minute.'), controls] }), view.learnStage);
  if (kept.learn) drawLearn();
  else view.learnStage.append(card({ body: empty({ icon: 'regulation', title: 'Expression from sequence', text: 'Two models of the same measurements: an elastic net on words of a few letters, and a convolutional network. Learn them to compare.' }) }));
}

async function learn() {
  if (kept.learnBusy) return;
  kept.learnBusy = true; kept.learnError = null;
  draw();
  if (view.learnStage) clear(view.learnStage).append(working('Fitting the elastic net on k-mers, then training the network…'));
  try {
    kept.learn = await api.run('regulation.learn', { ...DREAM, train: 8500, epochs: kept.epochs, filters: kept.filters, width: 15, seed: 7 });
  } catch (e) {
    kept.learnError = e.message;
  } finally {
    kept.learnBusy = false;
    if (kept.tab === 'learnt') draw();
  }
}

function drawLearn() {
  const box = view.learnStage && clear(view.learnStage);
  if (!box) return;
  if (kept.learnError) { box.append(callout(kept.learnError, 'bad')); return; }
  const r = kept.learn;
  if (!r) return;
  box.append(metrics([
    [fmt.number(r.baseline.pearson, 3), `r, words (${r.baseline.nonzero} kept)`],
    [fmt.number(r.cnn.pearson, 3), `r, network (${fmt.count(r.cnn.parameters)} parameters)`],
    [fmt.number(r.baseline.spearman, 3), 'ρ, words'], [fmt.number(r.cnn.spearman, 3), 'ρ, network'],
    [fmt.count(r.sizes?.held_out ?? r.held_out.length), 'held out'],
  ]));
  const ho = r.held_out;
  const xs = ho.map((x) => x[0]);
  const lim = [Math.min(...xs), Math.max(...xs)];
  const scatter = (k, label, colour) => chart({ width: 360, height: 300, x: lim, y: lim, xlabel: 'measured', ylabel: `predicted by ${label}`,
    series: [{ points: ho.map((x) => [x[0], x[k]]), colour, line: false, dots: true }, { points: [lim, lim].map((_, i) => [lim[i], lim[i]]), colour: token('--hair-3'), dash: '3 3', width: 1 }] }).replace(/r="3"/g, 'r="1.5" fill-opacity="0.35"');
  const curve = r.cnn.curve || [];
  const lc = chart({ width: 360, height: 300, x: [1, Math.max(2, curve.length)], y: [0, Math.max(...curve.flatMap((c) => [c[1], c[2]])) * 1.05], xlabel: 'epoch', ylabel: 'loss',
    series: [{ points: curve.map((c) => [c[0], c[1]]), colour: token('--accent'), dots: true, label: 'learning' }, { points: curve.map((c) => [c[0], c[2]]), colour: token('--C'), dots: true, label: 'held out' }] });
  const fig = el('div', { style: { display: 'flex', gap: '10px', flexWrap: 'wrap' } });
  fig.innerHTML = `<div>${scatter(1, 'words', token('--accent'))}</div><div>${scatter(2, 'the network', token('--C'))}</div><div>${lc}</div>`;
  box.append(card({ title: 'Held out: predicted against measured', sub: 'the dashed line is perfect prediction; at this size of data the words do better than the network, as the measurements say', pad: true, body: fig }));
  // The words the elastic net weighs most.
  const words = (r.baseline.words || []).slice(0, 16);
  if (words.length) box.append(card({ title: 'The words that move expression', sub: 'the k-mers with the largest coefficients in the elastic net, up in mulberry, down in blue', pad: true,
    body: el('div', { class: 'chips' }, ...words.map(([w, c]) => el('span', { class: 'chip', style: { borderColor: c > 0 ? 'var(--accent-ring)' : 'var(--C)', fontFamily: 'var(--mono)' }, title: `coefficient ${fmt.number(c, 3)}` }, w, el('i', { style: { color: c > 0 ? 'var(--accent-deep)' : 'var(--C)' } }, ` ${c > 0 ? '+' : '−'}${fmt.number(Math.abs(c), 2)}`)))) }));
  // The network's filters, read back as motifs.
  const grid = el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(290px, 1fr))', gap: '10px' } });
  for (const f of (r.filters || []).filter((x) => x.sites > 0).sort((a, b) => (b.match ? 1 : 0) - (a.match ? 1 : 0) || b.bits - a.bits).slice(0, 12)) {
    const rows = f.logo.map((x) => { const s = x.heights.reduce((a, b) => a + b, 0) || 1; return x.heights.map((h) => h / s); });
    const c = el('div', { class: 'form' });
    c.innerHTML = logoSvg({ rows, height: 80 });
    c.append(el('div', { class: 'note' }, f.match ? `like ${f.match.name} (${f.match.id}), similarity ${fmt.number(f.match.similarity, 2)}, p ${fmt.number(f.match.p, 2)}` : `filter ${f.filter}: ${f.consensus}, no JASPAR match`));
    grid.append(c);
  }
  box.append(card({ title: 'What the network learnt', sub: 'each filter read back as the motif of the windows it fires on, and matched to JASPAR 2024\'s fungal profiles', pad: true, body: grid }));
}
