// The Expression view: RNA-seq, from reads to what changed. An RNA-seq sample is fragments of a
// tissue's transcripts in proportion to their abundance; pseudoalignment finds each read's
// equivalence class (the transcripts that hold all its k-mers) by hashing alone, and the EM shares
// each class's reads among its transcripts. Between two conditions, counts of a gene are
// overdispersed (negative binomial); dispersions are shrunk toward a trend in the mean, fold
// changes tested by Wald's statistic and the p-values adjusted for the false discovery rate.

import { el, clear, esc, toggles, field, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, CATEGORIES, warm } from '../palette.js';

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, quant: el('div'), diff: el('div'), q: { genes: 400, fragments: 400000, paralogs: 0.15, seed: 1, paired: false, estimate: 'em' },
           d: { genes: 4000, samples: 4, changed: 0.1, fold: 1, seed: 1 }, qr: null, dr: null };
  section.append(
    el('h2', { class: 'view-title' }, 'Expression'),
    el('p', { class: 'lede' }, 'An RNA-seq sample is a few million fragments of a tissue\'s transcripts, each transcript copied in proportion to its abundance and its length. The likelihood of the abundances needs only which transcripts each read could have come from: the intersection of the sets holding its k-mers, its equivalence class, found by hashing alone (Bray et al. 2016). The EM then shares each class\'s reads among its transcripts in proportion to their current abundances. Between conditions, a gene\'s counts vary more than Poisson allows; the negative binomial\'s dispersion, shrunk toward its trend in the mean, lets a few samples be tested (Love, Huber and Anders 2014).'),
    el('h3', { class: 'section-title' }, 'How much of each transcript'), view.quant,
    el('h3', { class: 'section-title' }, 'What changed between two conditions'), view.diff);
  buildQuant();
  buildDiff();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.qr) quantify();
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
function buildQuant() {
  const q = view.q, box = clear(view.quant);
  const genes = el('select', { onchange: () => { q.genes = Number(genes.value); } }, ...[100, 400, 1500].map((v) => el('option', { value: v, selected: v === q.genes }, `${v} genes`)));
  const frags = el('select', { onchange: () => { q.fragments = Number(frags.value); } }, ...[100000, 400000, 2000000].map((v) => el('option', { value: v, selected: v === q.fragments }, `${fmt.count(v)} fragments`)));
  const par = el('select', { onchange: () => { q.paralogs = Number(par.value); } }, ...[0, 0.15, 0.4].map((v) => el('option', { value: v, selected: v === q.paralogs }, `${fmt.percent(v, 0)} paralogs`)));
  box.append(el('div', { class: 'row' }, genes, frags, par, toggles([[false, 'single reads'], [true, 'pairs']], q.paired, (v) => { q.paired = v; }),
    el('button', { class: 'button', onclick: () => { q.seed++; quantify(); } }, 'another sample'), el('button', { class: 'button primary', onclick: quantify }, 'Quantify')),
    el('div', { class: 'quant-out' }));
}

async function quantify() {
  const out = view.quant.querySelector('.quant-out');
  clear(out).append(el('p', { class: 'wait' }, 'Simulating the transcriptome and its sample, pseudoaligning, and running the EM…'));
  try {
    const q = view.q;
    view.qr = await api.run('expression.quantify', { genes: q.genes, fragments: q.fragments, paralogs: q.paralogs, paired: q.paired, seed: q.seed, bootstrap: 30 });
    drawQuant();
  } catch (e) { clear(out).append(el('p', { class: 'error' }, e.message)); }
}

function drawQuant() {
  const r = view.qr, out = clear(view.quant.querySelector('.quant-out'));
  // [name, gene, length, true, em, split, true tpm, tpm, sd, own, exons]
  const tx = r.transcripts;
  const est = view.q.estimate;
  const pick = (t) => (est === 'em' ? t[4] : t[5]);
  const colour = (own) => warm(0.25 + 0.7 * own);
  const pts = tx.filter((t) => t[3] > 0).map((t) => ({ x: t[3] + 1, y: pick(t) + 1, colour: colour(t[9]), t }));
  const max = Math.max(...pts.map((p) => Math.max(p.x, p.y))) * 1.5;
  const s = r.seconds;
  out.append(el('p', { class: 'note' }, `${fmt.count(r.fragments)} fragments of ${fmt.count(tx.length)} transcripts of ${fmt.count(r.genes)} genes; ${fmt.percent(r.aligned / r.fragments, 1)} pseudoaligned to ${fmt.count(r.classes)} equivalence classes by ${fmt.count(r.kmers)} distinct 31-mers (${fmt.bytes(r.index_bytes)}). The EM converged in ${r.rounds} rounds. Times on ${r.threads} threads: the sample ${fmt.seconds(s.simulate)}, the index ${fmt.seconds(s.index)}, pseudoalignment ${fmt.seconds(s.pseudoalign)} (${fmt.count(Math.round(r.fragments / s.pseudoalign))} fragments a second), the EM ${fmt.seconds(s.em)}, ${r.bootstrap} bootstrap replicates ${fmt.seconds(s.bootstrap)}.`),
    el('div', { class: 'row' }, toggles([['em', 'the EM'], ['split', 'an equal split of each class']], est, (v) => { view.q.estimate = v; drawQuant(); }),
      el('span', { class: 'small' }, `correlation of log counts ${fmt.number(r[est === 'em' ? 'em' : 'split'].correlation, 4)}, median relative error ${fmt.number(r[est === 'em' ? 'em' : 'split'].median_error, 3)}`)));
  const sc = scatter({ width: 460, height: 380, points: pts, x: [1, max], y: [1, max], logx: true, logy: true, diagonal: true, xlabel: 'true fragments + 1', ylabel: `estimated fragments + 1 (${est === 'em' ? 'EM' : 'equal split'})`,
    hoverText: (p) => `${esc(p.t[0])} (gene ${p.t[1] + 1}, ${fmt.length(p.t[2])}, exons ${p.t[10].map((e) => e + 1).join(' ')})<br>true ${fmt.number(p.t[3], 4)}, EM ${fmt.number(p.t[4], 4)} ± ${fmt.number(2 * p.t[8], 3)}, split ${fmt.number(p.t[5], 4)}<br>${fmt.percent(p.t[9], 0)} of its fragments in classes of its own` });
  // Class sizes.
  const sizes = r.class_sizes, frs = r.class_fragments;
  const bars = sizes.map((n, k) => ({ x0: k - 0.4, x1: k + 0.4, y: n, colour: token('--accent'), title: `${fmt.count(n)} classes of ${k}${k === 8 ? ' or more' : ''} transcripts, holding ${fmt.count(frs[k])} fragments` })).slice(1);
  const cs = el('div');
  cs.innerHTML = chart({ width: 360, height: 220, x: [0.4, 8.6], y: [0, Math.max(...sizes.slice(1)) * 1.15], xlabel: 'transcripts in the class', ylabel: 'classes', bars, xticks: [1, 2, 3, 4, 5, 6, 7, 8] });
  out.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '28px' } }, sc,
    el('div', { style: { maxWidth: '400px' } }, el('p', { class: 'eyebrow' }, 'equivalence classes by size'), cs,
      el('p', { class: 'note' }, 'Each dot a transcript, coloured by the share of its fragments that fall in classes of its own: deep where many do, pale where almost every k-mer it has is shared with an isoform or a paralog. Pale transcripts scatter: from classes alone, an isoform e1e2e3 beside e1e2, e2e3 and e1e3 cannot be told apart (the theory note, §5.3), and the EM can only share. The equal split, which gives each transcript of a class the same share, shows what the EM\'s likelihood adds.'))));
  // A gene with the most isoforms: its transcripts, true and estimated, with the bootstrap's spread.
  const byGene = new Map();
  for (const t of tx) { if (!byGene.has(t[1])) byGene.set(t[1], []); byGene.get(t[1]).push(t); }
  const genes = [...byGene.values()].filter((ts) => ts.length >= 3 && ts.some((t) => t[3] > 50)).sort((a, b) => b.length - a.length).slice(0, 6);
  const tab = el('table', { class: 'data' }, el('tr', {}, ...['transcript', 'exons', 'true', 'EM ± 2 sd', 'split', 'own'].map((h, k) => el('th', { class: k > 1 ? 'num' : '' }, h))));
  for (const ts of genes) for (const t of ts) tab.append(el('tr', {}, el('td', { class: 'small' }, t[0]), el('td', { class: 'mono small' }, t[10].map((e) => e + 1).join('·')),
    el('td', { class: 'num' }, fmt.count(Math.round(t[3]))), el('td', { class: 'num' }, `${fmt.count(Math.round(t[4]))} ± ${fmt.count(Math.round(2 * t[8]))}`), el('td', { class: 'num' }, fmt.count(Math.round(t[5]))), el('td', { class: 'num' }, fmt.percent(t[9], 0))));
  out.append(el('p', { class: 'eyebrow', style: { marginTop: '14px' } }, 'genes with several isoforms: each isoform\'s exons, its fragments, and the estimates'), el('div', { style: { maxHeight: '320px', overflowY: 'auto' } }, tab));
}

// ---- differential expression -----------------------------------------------------------------------------
function buildDiff() {
  const d = view.d, box = clear(view.diff);
  const samples = el('select', { onchange: () => { d.samples = Number(samples.value); } }, ...[2, 3, 4, 6, 10].map((v) => el('option', { value: v, selected: v === d.samples }, `${v} samples a condition`)));
  const fold = el('select', { onchange: () => { d.fold = Number(fold.value); } }, ...[0.5, 1, 2].map((v) => el('option', { value: v, selected: v === d.fold }, `changed ${fmt.number(2 ** v, 3)}-fold`)));
  const changed = el('select', { onchange: () => { d.changed = Number(changed.value); } }, ...[0, 0.05, 0.1, 0.3].map((v) => el('option', { value: v, selected: v === d.changed }, `${fmt.percent(v, 0)} of genes changed`)));
  box.append(el('div', { class: 'row' }, samples, fold, changed, el('button', { class: 'button', onclick: () => { d.seed++; test(); } }, 'another study'), el('button', { class: 'button primary', onclick: test }, 'Test')),
    el('div', { class: 'diff-out' }));
}

async function test() {
  const out = view.diff.querySelector('.diff-out');
  clear(out).append(el('p', { class: 'wait' }, 'Drawing the counts and fitting a model to every gene…'));
  try {
    view.dr = await api.run('expression.differential', view.d);
    drawDiff();
  } catch (e) { clear(out).append(el('p', { class: 'error' }, e.message)); }
}

function drawDiff() {
  const r = view.dr, out = clear(view.diff.querySelector('.diff-out'));
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
  out.append(el('p', { class: 'note' }, `${fmt.count(G.length)} genes in ${r.samples} samples, ${fmt.count(r.de)} of them changed. At an adjusted p below 0.1 (Benjamini and Hochberg), ${fmt.count(r.discoveries)} genes are called, ${fmt.count(r.true_discoveries)} of them truly changed: power ${fmt.percent(r.power, 1)}, false discovery rate ${fmt.percent(r.fdr, 1)} (the target, 10%). The dispersion's trend is fitted as ${fmt.number(r.trend.a0, 3)} + ${fmt.number(r.trend.a1, 3)}/μ (the truth ${r.trend.true_a0} + ${r.trend.true_a1}/μ). Fitted in ${fmt.seconds(r.seconds)} on ${r.threads} threads. Green: called and changed; red: called but unchanged; blue: changed but missed.`),
    el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '20px', flexWrap: 'wrap' } },
      el('div', {}, el('p', { class: 'eyebrow' }, 'the MA plot'), ma), el('div', {}, el('p', { class: 'eyebrow' }, 'the volcano'), vol),
      el('div', {}, el('p', { class: 'eyebrow' }, 'dispersions: each gene\'s own (grey), shrunk toward the trend (blue)'), disp), el('div', {}, el('p', { class: 'eyebrow' }, 'the p-values: the null\'s are uniform'), ph)));
}
