// The Populations view: variation as the trace of a history. A sample's genealogy, drawn by the
// coalescent backwards in time, with the mutations that fell on its branches; along a
// recombining chromosome the genealogy changes from place to place. The site frequency
// spectrum and the estimators of θ read the genealogy's shape: Tajima's D compares two of them.
// Forwards in time, Wright–Fisher drift and selection move frequencies to loss or fixation as
// Kimura's diffusion predicts. And a sample from several populations, mixed, shows its history
// in its principal components and in the ancestry an EM recovers.

import { el, clear, esc, notice, field, toggles, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, warm, CATEGORIES } from '../palette.js';

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, coal: el('div'), wf: el('div'), structure: el('div'),
           c: { samples: 16, theta: 12, rho: 0, demography: 'constant', seed: 1, at: 0.5 },
           w: { individuals: 100, s: 0.01, h: 0.5, p0: 0.1, demes: 1, m: 0, generations: 1500, replicates: 1000, seed: 1 },
           st: { populations: 3, each: 25, admixed: 30, fst: 0.08, sites: 4000, seed: 1 } };
  section.append(
    el('h2', { class: 'view-title' }, 'Populations'),
    el('p', { class: 'lede' }, 'Look backwards from a sample of genes and their lineages meet, two at a time, at a rate that grows with the number of pairs: Kingman\'s coalescent, the genealogy that a Wright–Fisher population draws. Mutations fall on its branches as a Poisson process, so the variation of the sample is a picture of its tree. The number of segregating sites and the mean number of differences both estimate θ = 4Nμ; when the tree is stretched at its root or squeezed at its tips (a bottleneck, an expansion, selection) they disagree, and Tajima\'s D measures by how much.'),
    el('h3', { class: 'section-title' }, 'The coalescent'), view.coal,
    el('h3', { class: 'section-title' }, 'Wright–Fisher: drift and selection'), view.wf,
    el('h3', { class: 'section-title' }, 'Structure: populations and their mixture'), view.structure);
  buildCoalescent();
  buildWF();
  buildStructure();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.coalDone) { view.coalDone = true; runCoalescent(); runWF(); runStructure(); }
}

const num = (label, value, onchange, attrs = {}) => field(label, el('input', { type: 'number', value, step: 'any', onchange: (e) => onchange(Number(e.target.value)), ...attrs }));

// ---- the coalescent ----------------------------------------------------------------------------------------
function buildCoalescent() {
  const c = view.c, box = clear(view.coal);
  const dem = toggles([['constant', 'constant size'], ['growth', 'expansion'], ['bottleneck', 'a bottleneck']], c.demography, (v) => { c.demography = v; runCoalescent(); });
  box.append(el('div', { class: 'row' },
    num('samples n', c.samples, (v) => { c.samples = v; runCoalescent(); }, { min: 2, max: 60 }),
    num('θ', c.theta, (v) => { c.theta = v; runCoalescent(); }, { min: 0, max: 200 }),
    num('ρ', c.rho, (v) => { c.rho = v; runCoalescent(); }, { min: 0, max: 200 }),
    dem,
    el('button', { class: 'button', onclick: () => { c.seed++; runCoalescent(); } }, 'another sample')));
  view.coalBody = el('div');
  box.append(view.coalBody);
}

async function runCoalescent() {
  const c = view.c;
  const demography = c.demography === 'growth' ? { kind: 'growth', alpha: 8 } : c.demography === 'bottleneck' ? { kind: 'bottleneck', start: 0.05, duration: 0.05, size: 0.04 } : { kind: 'constant' };
  try {
    const r = await api.run('pop.coalescent', { samples: c.samples, theta: c.theta, rho: c.rho, demography, seed: c.seed, at: c.at });
    drawCoalescent(r);
  } catch (e) { clear(view.coalBody).append(el('p', { class: 'error' }, e.message)); }
}

// A tree with its root at the top and time downwards to the samples; returns the leaf order.
function treeSvg(t, W, H) {
  const n = t.nodes.length;
  const kids = Array.from({ length: n }, () => []);
  t.nodes.forEach(([, parent], v) => { if (parent >= 0) kids[parent].push(v); });
  const x = new Array(n).fill(0), order = [];
  const place = (v) => {
    if (!kids[v].length) { x[v] = order.length; order.push(v); return; }
    kids[v].sort((a, b) => t.nodes[b][0] - t.nodes[a][0]);
    for (const c of kids[v]) place(c);
    x[v] = kids[v].reduce((s, c) => s + x[c], 0) / kids[v].length;
  };
  place(t.root);
  const leaves = order.length;
  const pad = 26, left = 44;
  const X = (v) => left + (x[v] / Math.max(1, leaves - 1)) * (W - left - pad);
  const Y = (time) => H - 24 - (time / Math.max(1e-9, t.tmrca)) * (H - 48);
  let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`;
  // A time axis.
  for (const f of [0, 0.25, 0.5, 0.75, 1]) {
    const yy = Y(f * t.tmrca);
    s += `<line x1="${left - 6}" y1="${yy}" x2="${W - pad}" y2="${yy}" stroke="${token('--hair')}" />`;
    s += `<text class="tick-label" x="${left - 8}" y="${yy + 3}" text-anchor="end">${fmt.number(f * t.tmrca, 2)}</text>`;
  }
  t.nodes.forEach(([time, parent], v) => {
    if (parent < 0) return;
    const yp = Y(t.nodes[parent][0]);
    s += `<path d="M${X(v).toFixed(1)},${Y(time).toFixed(1)}V${yp.toFixed(1)}H${X(parent).toFixed(1)}" fill="none" stroke="${token('--ink-2')}" stroke-width="1.3" />`;
  });
  for (const [node, time, pos] of t.mutations) {
    const xx = X(node), yy = Y(time);
    s += `<line x1="${xx - 5}" y1="${yy}" x2="${xx + 5}" y2="${yy}" stroke="${token('--accent')}" stroke-width="2.4"><title>a mutation at ${fmt.number(pos, 3)} of the chromosome, ${fmt.number(time, 3)} back</title></line>`;
  }
  order.forEach((v) => { s += `<circle cx="${X(v).toFixed(1)}" cy="${Y(0)}" r="3.2" fill="${token('--C')}" />`; });
  s += `<text class="tick-label" x="${W - pad}" y="${H - 6}" text-anchor="end">time back, in units of 2N generations</text></svg>`;
  return { svg: s, order };
}

function drawCoalescent(r) {
  const box = clear(view.coalBody), c = view.c;
  const W = Math.min(600, view.section.clientWidth / 2 - 30), H = 340;
  const { svg, order } = treeSvg(r.tree, W, H);
  const left = el('div');
  left.innerHTML = svg;
  const facts = el('dl', { class: 'kv' });
  const row = (k, v) => facts.append(el('dt', {}, k), el('dd', {}, v));
  row('segregating sites S', fmt.count(r.segregating));
  row('θ from S (Watterson)', fmt.number(r.theta_w, 4));
  row('θ from π (Tajima)', fmt.number(r.pi, 4));
  row('θ from θ_H (Fay and Wu)', fmt.number(r.theta_h, 4));
  row("Tajima's D", r.tajima_d === null ? '–' : fmt.number(r.tajima_d, 3));
  row('T_MRCA here', `${fmt.number(r.tree.tmrca, 3)} (expected ${fmt.number(r.expected_tmrca, 3)} at constant size)`);
  if (c.rho > 0) row('along the chromosome', `${r.trees} marginal trees after ${r.recombinations} recombinations`);
  // The spectrum against θ/i.
  const n = r.samples, xi = r.spectrum;
  const bars = [];
  for (let i = 1; i < n; i++) bars.push({ x0: i - 0.42, x1: i + 0.42, y: xi[i], colour: token('--accent'), gap: 0, title: `${fmt.number(xi[i], 3)} sites with ${i} derived` });
  const expect = [];
  for (let i = 1; i < n; i++) expect.push([i, c.theta / i]);
  const sfs = el('div');
  sfs.innerHTML = chart({ width: 380, height: 210, x: [0.4, n - 0.4], y: [0, Math.max(c.theta, ...xi.slice(1, n)) * 1.1 || 1], xlabel: 'derived copies i', ylabel: 'sites', bars,
    series: [{ points: expect, colour: token('--ink-2'), dash: '4 3', label: 'θ / i' }], xticks: [...Array(Math.min(n - 1, 10)).keys()].map((k) => Math.round(1 + (k * (n - 2)) / Math.max(1, Math.min(n - 2, 9)))) });
  const right = el('div', {}, facts, el('p', { class: 'eyebrow' }, 'the site frequency spectrum'), sfs);
  box.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '30px' } }, left, right));
  // The marginal trees' T_MRCA along the chromosome; a click chooses where the tree is drawn.
  if (c.rho > 0 && r.tmrca_along?.length) {
    const pts = [];
    for (const [a, b, t] of r.tmrca_along) { pts.push([a, t], [b, t]); }
    const along = el('div', { style: { cursor: 'pointer' } });
    along.innerHTML = chart({ width: Math.min(900, view.section.clientWidth - 80), height: 150, x: [0, 1], y: [0, Math.max(...pts.map((p) => p[1])) * 1.1], xlabel: 'position on the chromosome (click to draw its tree)', ylabel: 'T_MRCA',
      series: [{ points: pts, colour: token('--accent'), area: true }], marks: [{ x: c.at, y: 0, text: '▲' }] });
    along.addEventListener('click', (e) => {
      const rect = along.querySelector('svg').getBoundingClientRect();
      const fx = (e.clientX - rect.left - 56) / (rect.width - 70);
      if (fx >= 0 && fx <= 1) { c.at = fx; runCoalescent(); }
    });
    box.append(el('p', { class: 'eyebrow' }, 'the genealogy changes along a recombining chromosome'), along);
  }
  // The haplotypes, in the order of the tree's leaves.
  if (r.haplotypes && r.haplotypes[0]?.length) {
    const S = r.haplotypes[0].length, rows = r.haplotypes.length;
    const cw = Math.max(2, Math.min(10, Math.floor(Math.min(900, view.section.clientWidth - 80) / S))), ch = 9;
    const canvas = el('canvas', { width: S * cw, height: rows * ch, class: 'dp' });
    const g = canvas.getContext('2d');
    order.forEach((leaf, k) => {
      const h = r.haplotypes[leaf];
      for (let j = 0; j < S; j++) { g.fillStyle = h[j] === '1' ? token('--accent') : token('--paper-3'); g.fillRect(j * cw, k * ch, cw - (cw > 3 ? 1 : 0), ch - 1); }
    });
    box.append(el('p', { class: 'eyebrow' }, `the sample: ${rows} sequences × ${S} segregating sites, derived alleles dark, in the order of the tree's leaves`), canvas);
  }
}

// ---- Wright–Fisher ----------------------------------------------------------------------------------------
function buildWF() {
  const w = view.w, box = clear(view.wf);
  box.append(el('div', { class: 'row' },
    num('N', w.individuals, (v) => { w.individuals = v; }, { min: 2 }),
    num('s', w.s, (v) => { w.s = v; }),
    num('h', w.h, (v) => { w.h = v; }),
    num('p₀', w.p0, (v) => { w.p0 = v; }, { min: 0, max: 1 }),
    num('demes', w.demes, (v) => { w.demes = v; }, { min: 1, max: 20 }),
    num('m', w.m, (v) => { w.m = v; }, { min: 0, max: 1 }),
    num('generations', w.generations, (v) => { w.generations = v; }),
    el('button', { class: 'button primary', onclick: runWF }, 'Run'), el('button', { class: 'button', onclick: () => { w.seed++; runWF(); } }, 'other replicates')));
  view.wfBody = el('div');
  box.append(view.wfBody);
}

async function runWF() {
  try {
    const r = await api.run('pop.wright_fisher', { ...view.w, replicates: view.w.replicates, draw: 60 });
    const box = clear(view.wfBody);
    const T = r.generations, every = r.every;
    const series = r.paths.map((p, k) => ({ points: p.map((v, i) => [i * every, v]), colour: CATEGORIES[k % CATEGORIES.length], width: 0.9 }));
    series.push({ points: r.mean_p.map((v, i) => [i * every, v]), colour: token('--ink'), width: 2.2, label: `mean of ${fmt.count(r.replicates)}` });
    const c = el('div');
    c.innerHTML = chart({ width: Math.min(640, view.section.clientWidth - 420), height: 280, x: [0, T], y: [0, 1], xlabel: 'generations', ylabel: 'frequency of A₁', series });
    const facts = el('dl', { class: 'kv' });
    const row = (k, v) => facts.append(el('dt', {}, k), el('dd', {}, v));
    row('fixed', `${fmt.count(r.fixed)} of ${fmt.count(r.replicates)}: ${fmt.percent(r.fixation, 1)}`);
    row('Kimura predicts', fmt.percent(r.kimura, 1));
    row('lost', fmt.count(r.lost));
    row('settled after', `${fmt.number(r.mean_settling, 4)} generations on average`);
    if (r.fst !== null) row('F_ST among demes', `${fmt.number(r.fst, 3)}${r.island_fst !== null ? ` (the island model: ${fmt.number(r.island_fst, 3)})` : ''}`);
    box.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '30px' } }, c, el('div', { style: { maxWidth: '380px' } }, facts,
      el('p', { class: 'note' }, 'Kimura\'s probability of fixation, u(p) = (1 − e^{−4Nsp}) / (1 − e^{−4Ns}) for genic selection (Kimura 1962), solves the diffusion\'s backward equation; with dominance it is a ratio of integrals, computed by quadrature. A neutral allele fixes with probability p₀.'))));
  } catch (e) { clear(view.wfBody).append(el('p', { class: 'error' }, e.message)); }
}

// ---- structure ----------------------------------------------------------------------------------------
function buildStructure() {
  const st = view.st, box = clear(view.structure);
  box.append(el('div', { class: 'row' },
    num('populations K', st.populations, (v) => { st.populations = v; }, { min: 2, max: 6 }),
    num('of each', st.each, (v) => { st.each = v; }, { min: 2, max: 100 }),
    num('admixed', st.admixed, (v) => { st.admixed = v; }, { min: 0, max: 200 }),
    num('F_ST', st.fst, (v) => { st.fst = v; }, { min: 0.001, max: 0.9 }),
    num('sites', st.sites, (v) => { st.sites = v; }, { min: 200, max: 20000 }),
    el('button', { class: 'button primary', onclick: runStructure }, 'Simulate and infer'), el('button', { class: 'button', onclick: () => { st.seed++; runStructure(); } }, 'another sample')));
  view.stBody = el('div');
  box.append(view.stBody);
}

const mixColour = (weights) => {
  const hex = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16));
  const rgb = [0, 0, 0];
  weights.forEach((w, k) => { const c = hex(CATEGORIES[k % CATEGORIES.length]); for (let i = 0; i < 3; i++) rgb[i] += w * c[i]; });
  return `rgb(${rgb.map((v) => Math.round(v)).join(',')})`;
};

async function runStructure() {
  clear(view.stBody).append(el('p', { class: 'wait' }, 'Simulating the sample, then its principal components and the EM of ancestry…'));
  try {
    const r = await api.run('pop.structure', view.st);
    const box = clear(view.stBody);
    const coords = r.pca.coordinates;
    const xs = coords.map((p) => p[0]), ys = coords.map((p) => p[1]);
    const W = 380, H = 320, pad = 34;
    const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
    const X = (v) => pad + ((v - x0) / Math.max(1e-9, x1 - x0)) * (W - 2 * pad), Y = (v) => H - pad - ((v - y0) / Math.max(1e-9, y1 - y0)) * (H - 2 * pad);
    let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}"><line x1="${pad}" y1="${H - pad}" x2="${W - pad}" y2="${H - pad}" stroke="${token('--hair-2')}"/><line x1="${pad}" y1="${pad}" x2="${pad}" y2="${H - pad}" stroke="${token('--hair-2')}"/>`;
    coords.forEach((p, i) => { s += `<circle cx="${X(p[0]).toFixed(1)}" cy="${Y(p[1]).toFixed(1)}" r="${r.admixed[i] ? 3.4 : 4.2}" fill="${mixColour(r.truth[i])}" fill-opacity="0.85" stroke="${r.admixed[i] ? token('--ink') : 'none'}" stroke-width="0.6"><title>individual ${i + 1}${r.admixed[i] ? ', admixed' : ''}</title></circle>`; });
    s += `<text class="tick-label" x="${W / 2}" y="${H - 8}" text-anchor="middle" style="font-style:italic">PC 1 (eigenvalue ${fmt.number(r.pca.values[0], 3)})</text><text class="tick-label" transform="translate(12 ${H / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">PC 2 (${fmt.number(r.pca.values[1], 3)})</text></svg>`;
    const pca = el('div');
    pca.innerHTML = s;
    // Ancestry: the EM's bars above, the truth's below; components matched to the truth's by best agreement.
    const K = r.k, n = r.samples;
    const perm = matchComponents(r.admixture.q, r.truth, K, r.populations);
    const idx = [...Array(n).keys()].sort((a, b) => r.population[a] - r.population[b] || r.truth[b][r.population[b]] - r.truth[a][r.population[a]]);
    const BW = Math.min(820, view.section.clientWidth - 470), BH = 90;
    const bw = BW / n;
    let b = `<svg width="${BW}" height="${2 * BH + 30}" viewBox="0 0 ${BW} ${2 * BH + 30}">`;
    let mae = 0;
    idx.forEach((i, k) => {
      let y = 0;
      for (let c = 0; c < K; c++) {
        const v = r.admixture.q[i][c], tc = perm[c];
        b += `<rect x="${(k * bw).toFixed(2)}" y="${(y * BH).toFixed(2)}" width="${(bw + 0.3).toFixed(2)}" height="${(v * BH).toFixed(2)}" fill="${CATEGORIES[tc % CATEGORIES.length]}" />`;
        y += v;
      }
      y = 0;
      for (let c = 0; c < r.populations; c++) {
        const v = r.truth[i][c];
        b += `<rect x="${(k * bw).toFixed(2)}" y="${(BH + 20 + y * BH).toFixed(2)}" width="${(bw + 0.3).toFixed(2)}" height="${(v * BH).toFixed(2)}" fill="${CATEGORIES[c % CATEGORIES.length]}" fill-opacity="0.85" />`;
        y += v;
      }
      for (let c = 0; c < Math.min(K, r.populations); c++) mae += Math.abs(r.admixture.q[i][c] - r.truth[i][perm[c]]);
    });
    b += `<text class="tick-label" x="0" y="${BH + 13}">the EM's ancestry, above; the truth, below</text></svg>`;
    const bars = el('div');
    bars.innerHTML = b;
    box.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '30px' } }, pca, el('div', {}, bars,
      el('p', { class: 'note' }, `${n} individuals at ${fmt.count(r.sites)} sites, from ${r.populations} populations diverged to F_ST ≈ ${view.st.fst} (Balding–Nichols), ${view.st.admixed} of them admixed. The principal components (Patterson, Price and Reich 2006) place each individual between the populations it descends from; the EM of ancestry (Tang et al. 2005, accelerated) recovers the proportions to a mean absolute error of ${fmt.number(mae / (n * Math.min(K, r.populations)), 3)} per component, in ${r.admixture.iterations} steps.`))));
  } catch (e) { clear(view.stBody).append(el('p', { class: 'error' }, e.message)); }
}

// The permutation of the EM's components that best matches the truth's (greedy on agreement).
function matchComponents(q, truth, K, T) {
  const score = Array.from({ length: K }, () => new Array(T).fill(0));
  q.forEach((row, i) => { for (let a = 0; a < K; a++) for (let b = 0; b < T; b++) score[a][b] += row[a] * truth[i][b]; });
  const perm = new Array(K).fill(0), used = new Set();
  const pairs = [];
  for (let a = 0; a < K; a++) for (let b = 0; b < T; b++) pairs.push([score[a][b], a, b]);
  pairs.sort((x, y) => y[0] - x[0]);
  const done = new Set();
  for (const [, a, b] of pairs) if (!done.has(a) && !used.has(b)) { perm[a] = b; done.add(a); used.add(b); }
  for (let a = 0; a < K; a++) if (!done.has(a)) perm[a] = (T + a) % CATEGORIES.length;
  return perm;
}
