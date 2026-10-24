// The Expression view: RNA-seq, from reads to what changed. An RNA-seq sample is fragments of a
// tissue's transcripts in proportion to their abundance; pseudoalignment finds each read's
// equivalence class (the transcripts that hold all its k-mers) by hashing alone, and the EM shares
// each class's reads among its transcripts. Between two conditions, counts of a gene are
// overdispersed (negative binomial); dispersions are shrunk toward a trend in the mean, fold
// changes tested by Wald's statistic and the p-values adjusted for the false discovery rate. The
// two instruments stand in the bench's column and their results on its stage; the values are
// kept in the address.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, CATEGORIES, warm } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, seg, fld, card, metrics, callout, empty, working, busy, onEnter } from '../ui/kit.js';
import { tableCard } from '../ui/results.js';

const GENES = [100, 400, 1500];
const FRAGMENTS = [[100000, '100 k'], [400000, '400 k'], [2000000, '2 M']];
const PARALOGS = [0, 0.15, 0.4];
const SAMPLES = [2, 3, 4, 6, 10];
const FOLDS = [0.5, 1, 2];
const CHANGED = [0, 0.05, 0.1, 0.3];

// The values kept in the address: [part, field, key, read].
const numberOf = (list) => (v) => (list.includes(Number(v)) ? Number(v) : undefined);
const seedOf = (v) => (/^\d{1,6}$/.test(v) ? Number(v) : undefined);
const HASH = [
  ['q', 'genes', 'genes', numberOf(GENES)], ['q', 'fragments', 'frags', numberOf(FRAGMENTS.map(([v]) => v))], ['q', 'paralogs', 'paralogs', numberOf(PARALOGS)],
  ['q', 'paired', 'pairs', (v) => (v === '1' ? true : v === '0' ? false : undefined)], ['q', 'seed', 'qseed', seedOf],
  ['d', 'samples', 'samples', numberOf(SAMPLES)], ['d', 'fold', 'fold', numberOf(FOLDS)], ['d', 'changed', 'changed', numberOf(CHANGED)], ['d', 'seed', 'dseed', seedOf],
];

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, q: { genes: 400, fragments: 400000, paralogs: 0.15, seed: 1, paired: false, estimate: 'em' },
           d: { genes: 4000, samples: 4, changed: 0.1, fold: 1, seed: 1 }, qr: null, dr: null, qRan: null, dRan: null, quanting: false, testing: false };
  adopt(ctx.params);
  view.controls = el('div', { class: 'controls' });
  view.qState = el('div', { class: 'results-head' });
  view.qBox = el('div', { class: 'stack', style: { gap: '16px' } });
  view.dState = el('div', { class: 'results-head', style: { marginTop: '18px' } });
  view.dBox = el('div', { class: 'stack', style: { gap: '16px' } });
  section.append(
    el('p', { class: 'lede' }, 'An RNA-seq sample is a few million fragments of a tissue\'s transcripts, each transcript copied in proportion to its abundance and its length. The likelihood of the abundances needs only which transcripts each read could have come from: the intersection of the sets holding its k-mers, its equivalence class, found by hashing alone (Bray et al. 2016). The EM then shares each class\'s reads among its transcripts in proportion to their current abundances. Between conditions, a gene\'s counts vary more than Poisson allows; the negative binomial\'s dispersion, shrunk toward its trend in the mean, lets a few samples be tested (Love, Huber and Anders 2014).'),
    el('div', { class: 'bench' }, view.controls, el('div', { class: 'stage' }, view.qState, view.qBox, view.dState, view.dBox)));
  buildControls();
  drawDiffEmpty();
}

export function show(ctx) {
  view.ctx = ctx;
  const asked = adopt(ctx.params);
  if (asked.q || asked.d) buildControls();
  if (!view.quanting && (!view.qr || asked.q)) quantify();
  else if (asked.est && view.qr) drawQuant();
  if (!view.testing && asked.d) test();
  if (!asked.q && !asked.d) remember();
}

// The address's values, taken where they are valid: which instruments' values differ from those
// in hand (a part counts as asked for only when one of its values is in the address).
function adopt(params) {
  const asked = { q: false, d: false, est: false };
  for (const [part, f, k, read] of HASH) {
    if (params[k] === undefined) continue;
    const v = read(params[k]);
    if (v === undefined) continue;
    if (part === 'd' && !view.dr) asked.d = true;
    if (v !== view[part][f]) { view[part][f] = v; asked[part] = true; }
  }
  if (['em', 'split'].includes(params.est) && params.est !== view.q.estimate) { view.q.estimate = params.est; asked.est = true; }
  return asked;
}

function remember() {
  if (view.ctx.view !== 'expression') return;
  const q = view.q, d = view.d;
  view.ctx.setParams({ genes: q.genes, frags: q.fragments, paralogs: q.paralogs, pairs: q.paired ? 1 : 0, qseed: q.seed, est: q.estimate,
    ...(view.dr || view.testing ? { samples: d.samples, fold: d.fold, changed: d.changed, dseed: d.seed } : {}) });
}

const dot = (colour, label) => el('span', { class: 'k' }, el('span', { class: 'swatch', style: { background: colour, borderRadius: '50%' } }), label);

// ---- the instruments -------------------------------------------------------------------------------------
function buildControls() {
  const q = view.q, d = view.d;
  view.qBtn = btn('Quantify', { icon: 'play', kind: 'primary', onClick: () => quantify() });
  view.qAgain = btn('Another sample', { icon: 'refresh', title: 'another sample of the same design: the next seed', onClick: () => { q.seed++; quantify(); } });
  if (view.quanting) { busy(view.qBtn, true, 'Quantifying…'); view.qAgain.disabled = true; }
  const sample = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'One sample'),
      fld('genes', seg(GENES.map((v) => [v, fmt.count(v)]), q.genes, (v) => { q.genes = v; }, { block: true }), { help: 'with their isoforms' }),
      fld('fragments', seg(FRAGMENTS, q.fragments, (v) => { q.fragments = v; }, { block: true })),
      fld('paralogs', seg(PARALOGS.map((v) => [v, v ? fmt.percent(v, 0) : 'none']), q.paralogs, (v) => { q.paralogs = v; }, { block: true }), { help: 'genes that share sequence' }),
      fld('reads', seg([[false, 'single'], [true, 'pairs']], q.paired, (v) => { q.paired = v; }, { block: true }))),
    el('div', { class: 'run-bar' }, view.qBtn, view.qAgain));
  onEnter(sample, quantify);
  view.dBtn = btn('Test', { icon: 'play', kind: 'primary', onClick: () => test() });
  view.dAgain = btn('Another study', { icon: 'refresh', title: 'another study of the same design: the next seed', onClick: () => { d.seed++; test(); } });
  if (view.testing) { busy(view.dBtn, true, 'Testing…'); view.dAgain.disabled = true; }
  const study = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Two conditions'),
      fld('samples a condition', seg(SAMPLES.map((v) => [v, String(v)]), d.samples, (v) => { d.samples = v; }, { block: true }), { help: `${fmt.count(d.genes)} genes` }),
      fld('a changed gene changes', seg(FOLDS.map((v) => [v, `${fmt.number(2 ** v, 3)}-fold`]), d.fold, (v) => { d.fold = v; }, { block: true })),
      fld('genes changed', seg(CHANGED.map((v) => [v, v ? fmt.percent(v, 0) : 'none']), d.changed, (v) => { d.changed = v; }, { block: true }))),
    el('div', { class: 'run-bar' }, view.dBtn, view.dAgain));
  onEnter(study, test);
  const own = el('div', { class: 'toolbox-hint' }, icon('toolbox', 'sm'), el('span', {},
    'These samples are simulated, so every transcript\'s truth is known. With your own transcripts and reads, ',
    el('button', { class: 'link', onclick: () => view.ctx.openTool('expr-quant') }, 'quantify them'), ' in the Toolbox.'));
  clear(view.controls).append(sample, study, own);
}

// A scatter on a canvas, for thousands of points: log axes optional, a diagonal optional.
function scatter({ width = 420, height = 340, points, x, y, logx = false, logy = false, xlabel, ylabel, diagonal = false, lines = [], hoverText }) {
  const m = { l: 54, r: 12, t: 10, b: 40 };
  const W = width - m.l - m.r, H = height - m.t - m.b;
  const fx = logx ? Math.log10 : (v) => v, fy = logy ? Math.log10 : (v) => v;
  const sx = (v) => m.l + ((fx(v) - fx(x[0])) / (fx(x[1]) - fx(x[0]))) * W, sy = (v) => m.t + H - ((fy(v) - fy(y[0])) / (fy(y[1]) - fy(y[0]))) * H;
  const wrap = el('div', { style: { position: 'relative', width: `${width}px`, height: `${height}px` } });
  const axes = el('div');
  axes.innerHTML = chart({ width, height, x, y, log: (logx ? 'x' : '') + (logy ? 'y' : ''), xlabel, ylabel,
    series: [...(diagonal ? [{ points: [[x[0], x[0]], [x[1], x[1]]], colour: token('--ink-faint'), dash: '4 3', width: 1 }] : []), ...lines] });
  const canvas = el('canvas', { width: W * 2, height: H * 2, style: { position: 'absolute', left: `${m.l}px`, top: `${m.t}px`, width: `${W}px`, height: `${H}px` } });
  const g = canvas.getContext('2d');
  g.scale(2, 2);
  for (const p of points) {
    if (!(p.x > (logx ? 0 : -Infinity)) || !(p.y > (logy ? 0 : -Infinity))) continue;
    g.fillStyle = p.colour;
    g.globalAlpha = p.alpha ?? 0.6;
    g.beginPath();
    g.arc(sx(p.x) - m.l, sy(p.y) - m.t, p.r ?? 1.6, 0, 2 * Math.PI);
    g.fill();
  }
  wrap.append(axes, canvas);
  if (hoverText) {
    canvas.addEventListener('pointermove', (e) => {
      const rect = canvas.getBoundingClientRect();
      const cx = e.clientX - rect.left + m.l, cy = e.clientY - rect.top + m.t;
      let best = null, bd = 64;
      for (const p of points) {
        const dx = sx(p.x) - cx, dy = sy(p.y) - cy, d2 = dx * dx + dy * dy;
        if (d2 < bd) { bd = d2; best = p; }
      }
      hover(best ? hoverText(best) : null, e);
    });
    canvas.addEventListener('pointerleave', () => hover(null));
  }
  return wrap;
}

// ---- quantification --------------------------------------------------------------------------------------
async function quantify() {
  if (view.quanting) return;
  view.quanting = true;
  const q = { ...view.q };
  busy(view.qBtn, true, 'Quantifying…');
  view.qAgain.disabled = true;
  view.qState.replaceChildren(working(`Simulating ${fmt.count(q.fragments)} fragments of ${fmt.count(q.genes)} genes' transcripts, pseudoaligning, and running the EM…`));
  if (view.qr) view.qBox.classList.add('stale');
  else view.qBox.replaceChildren(el('div', { class: 'skeleton', style: { height: '74px' } }), el('div', { class: 'skeleton', style: { height: '470px' } }));
  remember();
  try {
    view.qr = await api.run('expression.quantify', { genes: q.genes, fragments: q.fragments, paralogs: q.paralogs, paired: q.paired, seed: q.seed, bootstrap: 30 });
    view.qRan = { q, at: new Date() };
    drawQuant();
  } catch (e) {
    view.qState.replaceChildren(callout(`The sample could not be quantified: ${e.message}`, 'bad'));
    if (!view.qr) clear(view.qBox);
  } finally {
    view.quanting = false;
    view.qBox.classList.remove('stale');
    busy(view.qBtn, false);
    view.qAgain.disabled = false;
  }
}

function drawQuant() {
  const r = view.qr, q = view.qRan.q, out = clear(view.qBox);
  view.qState.replaceChildren(el('span', { class: 'eyebrow', style: { margin: 0 } }, 'One sample'),
    el('span', { class: 'when' }, `${view.qRan.at.toLocaleTimeString()} · ${fmt.count(q.genes)} genes, ${fmt.percent(q.paralogs, 0)} paralogs, ${q.paired ? 'pairs' : 'single reads'}, seed ${q.seed}`));
  // [name, gene, length, true, em, split, true tpm, tpm, sd, own, exons]
  const tx = r.transcripts;
  const est = view.q.estimate;
  const pick = (t) => (est === 'em' ? t[4] : t[5]);
  const colour = (own) => warm(0.25 + 0.7 * own);
  const pts = tx.filter((t) => t[3] > 0).map((t) => ({ x: t[3] + 1, y: pick(t) + 1, colour: colour(t[9]), t }));
  const max = Math.max(...pts.map((p) => Math.max(p.x, p.y))) * 1.5;
  const s = r.seconds, quality = r[est === 'em' ? 'em' : 'split'];
  out.append(metrics([
    [fmt.count(r.fragments), `fragments of ${fmt.count(tx.length)} transcripts`],
    [fmt.percent(r.aligned / r.fragments, 1), 'pseudoaligned'],
    [fmt.count(r.classes), 'equivalence classes'],
    [fmt.count(r.rounds), 'rounds of the EM'],
    [fmt.number(quality.correlation, 4), `correlation of log counts, ${est === 'em' ? 'EM' : 'equal split'}`],
    [fmt.number(quality.median_error, 3), `median relative error, ${est === 'em' ? 'EM' : 'equal split'}`],
  ]));
  const sc = scatter({ width: 460, height: 380, points: pts, x: [1, max], y: [1, max], logx: true, logy: true, diagonal: true, xlabel: 'true fragments + 1', ylabel: `estimated fragments + 1 (${est === 'em' ? 'EM' : 'equal split'})`,
    hoverText: (p) => `${esc(p.t[0])} (gene ${p.t[1] + 1}, ${fmt.length(p.t[2])}, exons ${p.t[10].map((e) => e + 1).join(' ')})<br>true ${fmt.number(p.t[3], 4)}, EM ${fmt.number(p.t[4], 4)} ± ${fmt.number(2 * p.t[8], 3)}, split ${fmt.number(p.t[5], 4)}<br>${fmt.percent(p.t[9], 0)} of its fragments in classes of its own` });
  // Class sizes.
  const sizes = r.class_sizes, frs = r.class_fragments;
  const bars = sizes.map((n, k) => ({ x0: k - 0.4, x1: k + 0.4, y: n, colour: token('--accent'), title: `${fmt.count(n)} classes of ${k}${k === 8 ? ' or more' : ''} transcripts, holding ${fmt.count(frs[k])} fragments` })).slice(1);
  const cs = el('div');
  cs.innerHTML = chart({ width: 360, height: 220, x: [0.4, 8.6], y: [0, Math.max(...sizes.slice(1)) * 1.15], xlabel: 'transcripts in the class', ylabel: 'classes', bars, xticks: [1, 2, 3, 4, 5, 6, 7, 8] });
  // Where the time went.
  const parts = [['simulate', 'the sample'], ['index', 'the index'], ['pseudoalign', 'pseudoalignment'], ['em', 'the EM'], ['bootstrap', `${r.bootstrap} bootstraps`]];
  const total = parts.reduce((a, [k]) => a + (s[k] || 0), 0);
  const time = el('div', {}, el('p', { class: 'eyebrow', style: { marginTop: '14px' } }, `where the time went: ${fmt.seconds(total)} on ${r.threads} threads, ${fmt.count(Math.round(r.fragments / s.pseudoalign))} fragments a second pseudoaligned`),
    el('div', { class: 'stackbar' }, ...parts.map(([k, label], i) => el('span', { style: { width: `${(100 * (s[k] || 0)) / total}%`, background: CATEGORIES[i] }, title: `${label}: ${fmt.seconds(s[k] || 0)}` }))),
    el('div', { class: 'key-row flush' }, ...parts.map(([k, label], i) => el('span', { class: 'k' }, el('span', { class: 'swatch', style: { background: CATEGORIES[i] } }), `${label} ${fmt.seconds(s[k] || 0)}`))));
  const ramp = el('span', { class: 'ex-ramp', style: { background: `linear-gradient(90deg, ${colour(0)}, ${colour(0.5)}, ${colour(1)})` } });
  out.append(card({ title: 'How much of each transcript', sub: 'each dot a transcript: its true fragments against the estimate',
    actions: [seg([['em', 'the EM'], ['split', 'an equal split of each class']], est, (v) => { view.q.estimate = v; drawQuant(); remember(); })], body: [
      el('div', { class: 'ex-split card-body' }, sc,
        el('div', { class: 'ex-side' }, el('p', { class: 'eyebrow' }, 'equivalence classes by size'), cs,
          el('p', { class: 'note' }, `${fmt.count(r.kmers)} distinct 31-mers index the transcriptome in ${fmt.bytes(r.index_bytes)}; the EM converged in ${r.rounds} rounds.`))),
      el('div', { class: 'key-row' }, el('span', { class: 'k' }, el('span', {}, 'shared'), ramp, el('span', {}, 'its own')), el('span', { class: 'k faint' }, 'the share of a transcript\'s fragments in classes of its own'),
        el('span', { class: 'k' }, el('span', { class: 'ex-dash' }), 'estimate = truth')),
      el('div', { class: 'card-body', style: { paddingTop: '6px' } },
        el('p', { class: 'note', style: { marginTop: 0 } }, 'Pale transcripts scatter: from classes alone, an isoform e1e2e3 beside e1e2, e2e3 and e1e3 cannot be told apart (the theory note, §5.3), and the EM can only share. The equal split, which gives each transcript of a class the same share, shows what the EM\'s likelihood adds.'),
        time)] }));
  // Genes with the most isoforms: their transcripts, true and estimated, with the bootstrap's spread.
  const byGene = new Map();
  for (const t of tx) { if (!byGene.has(t[1])) byGene.set(t[1], []); byGene.get(t[1]).push(t); }
  const genes = [...byGene.values()].filter((ts) => ts.length >= 3 && ts.some((t) => t[3] > 50)).sort((a, b) => b.length - a.length).slice(0, 6);
  const rows = genes.flat().map((t) => ({ name: t[0], exons: t[10].map((e) => e + 1).join('·'), truth: Math.round(t[3]), em: Math.round(t[4]), sd2: Math.round(2 * t[8]), split: Math.round(t[5]), own: t[9] }));
  if (rows.length) {
    out.append(tableCard('Genes with several isoforms', [
      { key: 'name', label: 'transcript' }, { key: 'exons', label: 'exons', format: (x) => el('span', { class: 'mono' }, x) },
      { key: 'truth', label: 'true', num: true, format: fmt.count }, { key: 'em', label: 'EM ± 2 sd', num: true, format: (x, row) => `${fmt.count(x)} ± ${fmt.count(row.sd2)}` },
      { key: 'split', label: 'split', num: true, format: fmt.count }, { key: 'own', label: 'own', num: true, format: (x) => fmt.percent(x, 0) },
    ], rows, { sub: 'each isoform\'s exons, its fragments, and the estimates', file: 'isoforms.csv' }));
  }
}

// ---- differential expression -----------------------------------------------------------------------------
function drawDiffEmpty() {
  view.dState.replaceChildren(el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Two conditions'), el('span', { class: 'when' }, 'not tested yet'));
  view.dBox.replaceChildren(card({ body: empty({ icon: 'expression', title: 'What changed between two conditions',
    text: `A study of ${fmt.count(view.d.genes)} genes in two conditions, some of them changed: each gene's counts are fitted by a negative binomial, its dispersion shrunk toward the trend, and its change tested.`,
    action: btn('Test', { icon: 'play', onClick: () => test() }) }) }));
}

async function test() {
  if (view.testing) return;
  view.testing = true;
  const d = { ...view.d };
  busy(view.dBtn, true, 'Testing…');
  view.dAgain.disabled = true;
  view.dState.replaceChildren(working(`Drawing the counts of ${fmt.count(d.genes)} genes in ${d.samples} samples a condition and fitting a model to every gene…`));
  if (view.dr) view.dBox.classList.add('stale');
  else view.dBox.replaceChildren(el('div', { class: 'skeleton', style: { height: '74px' } }), el('div', { class: 'skeleton', style: { height: '760px' } }));
  remember();
  const r0 = view.dState.getBoundingClientRect();
  if (r0.top > window.innerHeight - 120) view.dState.scrollIntoView({ behavior: 'smooth', block: 'start' });
  try {
    view.dr = await api.run('expression.differential', d);
    view.dRan = { d, at: new Date() };
    drawDiff();
  } catch (e) {
    view.dState.replaceChildren(callout(`The study could not be tested: ${e.message}`, 'bad'));
    if (!view.dr) drawDiffEmpty();
  } finally {
    view.testing = false;
    view.dBox.classList.remove('stale');
    busy(view.dBtn, false);
    view.dAgain.disabled = false;
    remember();
  }
}

function drawDiff() {
  const r = view.dr, d = view.dRan.d, out = clear(view.dBox);
  view.dState.replaceChildren(el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Two conditions'),
    el('span', { class: 'when' }, `${view.dRan.at.toLocaleTimeString()} · ${d.samples} samples a condition, ${fmt.percent(d.changed, 0)} of genes changed ${fmt.number(2 ** d.fold, 3)}-fold, seed ${d.seed}`));
  // [base mean, lfc, p, padj, disp gene, disp trend, disp, outlier, true lfc, true disp]
  const G = r.genes;
  const col = (g) => (g[3] !== null && g[3] < 0.1 ? (g[8] !== 0 ? token('--good') : token('--bad')) : g[8] !== 0 ? token('--C') : token('--ink-faint'));
  const tip = (p) => `base mean ${fmt.number(p.g[0], 4)}, log₂ fold ${fmt.number(p.g[1], 3)} (true ${p.g[8]})<br>p ${p.g[2] === null ? '–' : fmt.number(p.g[2], 3)}, adjusted ${p.g[3] === null ? '–' : fmt.number(p.g[3], 3)}<br>dispersion: the gene's own ${fmt.number(p.g[4], 3)}, the trend ${fmt.number(p.g[5], 3)}, shrunk ${fmt.number(p.g[6], 3)} (true ${fmt.number(p.g[9], 3)})`;
  const ok = G.filter((g) => g[0] > 0);
  const means = ok.map((g) => g[0]);
  const maxL = Math.max(1, ...ok.map((g) => Math.abs(g[1]))) * 1.1;
  const ma = scatter({ width: 420, height: 320, points: ok.map((g) => ({ x: g[0], y: g[1], colour: col(g), g, r: 1.4, alpha: 0.65 })), x: [Math.min(...means) * 0.8, Math.max(...means) * 1.2], y: [-maxL, maxL], logx: true,
    xlabel: 'mean of normalised counts', ylabel: 'log₂ fold change', hoverText: tip });
  const ps = ok.filter((g) => g[2] !== null);
  const maxP = Math.max(...ps.map((g) => -Math.log10(Math.max(g[2], 1e-300))));
  const vol = scatter({ width: 420, height: 320, points: ps.map((g) => ({ x: g[1], y: -Math.log10(Math.max(g[2], 1e-300)), colour: col(g), g, r: 1.4, alpha: 0.65 })), x: [-maxL, maxL], y: [0, maxP * 1.05],
    xlabel: 'log₂ fold change', ylabel: '−log₁₀ p', hoverText: tip });
  // Dispersions: each gene's own estimate, and the shrunk ones, against the mean, with the trend.
  const trend = (mu) => r.trend.a0 + r.trend.a1 / mu, truth = (mu) => r.trend.true_a0 + r.trend.true_a1 / mu;
  const lo = Math.min(...means), hi = Math.max(...means);
  const curve = (f) => Array.from({ length: 60 }, (_, k) => { const mu = lo * (hi / lo) ** (k / 59); return [mu, f(mu)]; });
  const disps = ok.flatMap((g) => [{ x: g[0], y: g[4], colour: token('--ink-faint'), g, r: 1.2, alpha: 0.45 }, { x: g[0], y: g[6], colour: g[7] ? token('--bad') : token('--C'), g, r: 1.2, alpha: 0.6 }]);
  const dvals = disps.map((p) => p.y).filter((v) => v > 0);
  const disp = scatter({ width: 420, height: 320, points: disps, x: [lo * 0.8, hi * 1.2], y: [Math.max(1e-4, Math.min(...dvals) * 0.8), Math.max(...dvals) * 1.3], logx: true, logy: true,
    xlabel: 'mean of normalised counts', ylabel: 'dispersion α', hoverText: tip,
    lines: [{ points: curve(trend), colour: token('--accent'), width: 2, label: 'fitted trend' }, { points: curve(truth), colour: token('--ink'), dash: '4 3', width: 1, label: 'the true trend' }] });
  // The p-values: uniform under the null, a peak at 0 from the genes that changed.
  const bins = Array(20).fill(0), binsDe = Array(20).fill(0);
  for (const g of ps) { const k = Math.min(19, Math.floor(g[2] * 20)); bins[k]++; if (g[8] !== 0) binsDe[k]++; }
  const ph = el('div');
  ph.innerHTML = chart({ width: 420, height: 320, x: [0, 1], y: [0, Math.max(...bins) * 1.1], xlabel: 'p', ylabel: 'genes',
    bars: [...bins.map((n, k) => ({ x0: k / 20, x1: (k + 1) / 20, y: n, colour: token('--ink-faint'), title: `${fmt.count(n)} genes with p in [${k / 20}, ${(k + 1) / 20})` })),
      ...binsDe.map((n, k) => ({ x0: k / 20, x1: (k + 1) / 20, y: n, colour: token('--good'), title: `${fmt.count(n)} of them changed` }))] });
  out.append(metrics([
    [fmt.count(G.length), `genes in ${r.samples} samples`],
    [fmt.count(r.de), 'changed'],
    [fmt.count(r.discoveries), 'called at an adjusted p below 0.1'],
    [fmt.count(r.true_discoveries), 'of those, truly changed'],
    [fmt.percent(r.power, 1), 'power'],
    [fmt.percent(r.fdr, 1), 'false discovery rate (the target, 10%)'],
  ]));
  const panel = (title, figure, keys = null) => el('div', { class: 'ex-panel' }, el('p', { class: 'eyebrow' }, title), figure, keys);
  out.append(card({ title: 'What changed between two conditions', sub: 'each dot a gene; point at one for its numbers', body: [
    el('div', { class: 'key-row', style: { borderTop: 0, borderBottom: '1px solid var(--hair)' } },
      dot(token('--good'), 'called, and changed'), dot(token('--bad'), 'called, but unchanged'), dot(token('--C'), 'changed, but missed'), dot(token('--ink-faint'), 'unchanged, not called')),
    el('div', { class: 'ex-grid card-body' },
      panel('the MA plot', ma), panel('the volcano', vol),
      panel('dispersions against the mean', disp, el('div', { class: 'key-row flush' }, dot(token('--ink-faint'), 'each gene\'s own'), dot(token('--C'), 'shrunk toward the trend'), dot(token('--bad'), 'an outlier, not shrunk'))),
      panel('the p-values: the null\'s are uniform', ph, el('div', { class: 'key-row flush' }, dot(token('--ink-faint'), 'all genes'), dot(token('--good'), 'those that changed')))),
    el('div', { class: 'card-body', style: { paddingTop: 0 } }, el('p', { class: 'note', style: { marginTop: 0 } },
      `${fmt.count(G.length)} genes in ${r.samples} samples, ${fmt.count(r.de)} of them changed. At an adjusted p below 0.1 (Benjamini and Hochberg), ${fmt.count(r.discoveries)} genes are called, ${fmt.count(r.true_discoveries)} of them truly changed: power ${fmt.percent(r.power, 1)}, false discovery rate ${fmt.percent(r.fdr, 1)} (the target, 10%). The dispersion's trend is fitted as ${fmt.number(r.trend.a0, 3)} + ${fmt.number(r.trend.a1, 3)}/μ (the truth ${r.trend.true_a0} + ${r.trend.true_a1}/μ). Fitted in ${fmt.seconds(r.seconds)} on ${r.threads} threads.`))] }));
}
