// The Populations view: variation as the trace of a history. A sample's genealogy, drawn by the
// coalescent backwards in time, with the mutations that fell on its branches; along a
// recombining chromosome the genealogy changes from place to place. The site frequency
// spectrum and the estimators of θ read the genealogy's shape: Tajima's D compares two of them.
// Forwards in time, Wright–Fisher drift and selection move frequencies to loss or fixation as
// Kimura's diffusion predicts. And a sample from several populations, mixed, shows its history
// in its principal components and in the ancestry an EM recovers.
//
// Three simulators, each a card: its parameters in a row along its top, its figures beneath. The
// coalescent answers every change at once; the other two run when asked (or on Enter).

import { el, clear } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, seg, numberField, fld, card, metrics, callout, kv } from '../ui/kit.js';

let view = null;
// The parameters and the last results, kept when the view is mounted again (a change of theme).
const kept = {
  c: { samples: 16, theta: 12, rho: 0, demography: 'constant', seed: 1, at: 0.5 },
  w: { individuals: 100, s: 0.01, h: 0.5, p0: 0.1, demes: 1, m: 0, generations: 1500, replicates: 1000, seed: 1 },
  st: { populations: 3, each: 25, admixed: 30, fst: 0.08, sites: 4000, seed: 1 },
  coal: null, wf: null, structure: null, started: false,
};

export function mount(section, ctx) {
  view = { ctx, section, tickets: { coal: 0, wf: 0, st: 0 } };
  section.append(
    el('p', { class: 'lede' }, 'Look backwards from a sample of genes and their lineages meet, two at a time, at a rate that grows with the number of pairs: Kingman\'s coalescent, the genealogy that a Wright–Fisher population draws. Mutations fall on its branches as a Poisson process, so the variation of the sample is a picture of its tree. The number of segregating sites and the mean number of differences both estimate θ = 4Nμ; when the tree is stretched at its root or squeezed at its tips (a bottleneck, an expansion, selection) they disagree, and Tajima\'s D measures by how much.'),
    el('p', { class: 'toolbox-hint top' }, icon('toolbox', 'sm'), el('span', {}, 'Each simulation carries its truth, against which its estimates are read. A population of your own, from its variants: ',
      el('button', { type: 'button', class: 'link', onclick: () => view.ctx.openTool('pop-vcf') }, 'A population from its VCF'), ', in the Toolbox.')),
    el('div', { class: 'stack sims' }, coalescentCard(), wrightFisherCard(), structureCard()));
  if (kept.coal) drawCoalescent();
  if (kept.wf) drawWF();
  if (kept.structure) drawStructure();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!kept.started) { kept.started = true; runCoalescent(); runWF(); runStructure(); }
}

// ---- the pieces of a simulator ------------------------------------------------------------------------
// A number in the row of parameters; steps that are not whole are rounded, so that 0.1 + 0.05 shows
// as 0.15. `live` runs the simulation on every change.
function num(label, state, key, { min = null, max = null, step = 1, title = null, unit = '' }, onChange = () => {}) {
  const f = numberField({ value: state[key], min, max, step, unit, onChange: (v) => {
    const r = +v.toPrecision(10);
    if (r !== v) f.input.value = r;
    if (state[key] === r) return;
    state[key] = r;
    onChange(r);
  } });
  const node = fld(label, f);
  if (title) node.title = title;
  return node;
}

// Enter in a field of the row commits it and runs.
function onEnter(strip, run) {
  strip.addEventListener('keydown', (e) => {
    if (e.key !== 'Enter' || e.target.tagName !== 'INPUT') return;
    e.preventDefault();
    e.target.dispatchEvent(new Event('change'));
    run();
  });
}

// A run in progress: the card's figures dimmed and kept until the new ones come, a line with the
// seconds since it began in the card's head, its button disabled. Returns the function that ends it.
function working(sim, text, button = null) {
  const elapsed = el('span', { class: 'elapsed' });
  const line = el('span', { class: 'wait busy-line' }, el('span', { class: 'spin' }), el('span', {}, text), elapsed);
  sim.status.replaceChildren(line);
  sim.figs.classList.add('is-busy');
  const keep = button ? [...button.childNodes] : null;
  if (button) { button.disabled = true; button.replaceChildren(el('span', { class: 'spin' }), el('span', {}, 'Running…')); }
  const t0 = performance.now();
  const tick = () => { elapsed.textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`; };
  tick();
  const timer = setInterval(tick, 100);
  sim.line = line;
  return () => {
    clearInterval(timer);
    if (button) { button.replaceChildren(...keep); button.disabled = false; }
    if (sim.line !== line) return;   // a later run of the same card is still going
    sim.line = null;
    sim.status.replaceChildren();
    sim.figs.classList.remove('is-busy');
  };
}

function simCard(title, sub, strip) {
  const status = el('div', { class: 'status' });
  const figs = el('div', { class: 'figs' }, el('div', { class: 'wait', style: { padding: '30px 4px' } }, el('span', { class: 'spin' }), 'Simulating…'));
  const c = card({ cls: 'sim', title, sub, actions: [status], body: [strip, figs] });
  return { card: c, status, figs };
}

function failed(sim, e) { clear(sim.figs).append(el('div', { style: { padding: '2px 0' } }, callout(e.message, 'bad'))); }

// ---- the coalescent ----------------------------------------------------------------------------------------
function coalescentCard() {
  const c = kept.c, again = () => runCoalescent();
  const strip = el('div', { class: 'ctl-strip' },
    num('samples n', c, 'samples', { min: 2, max: 60, title: 'genes sampled' }, again),
    num('θ = 4Nμ', c, 'theta', { min: 0, max: 200, title: 'the scaled rate of mutation of the whole stretch' }, again),
    num('ρ = 4Nr', c, 'rho', { min: 0, max: 200, title: 'the scaled rate of recombination; 0 for a single tree' }, again),
    fld('the population\'s size through time', seg([['constant', 'constant'], ['growth', 'an expansion'], ['bottleneck', 'a bottleneck']], c.demography, (v) => { c.demography = v; runCoalescent(); })),
    el('span', { class: 'grow' }),
    btn('Another sample', { icon: 'refresh', title: 'the same parameters, another seed', onClick: () => { c.seed++; runCoalescent(); } }));
  onEnter(strip, () => {});
  view.coal = simCard('The coalescent', 'a sample\'s genealogy, backwards in time, with the mutations that fell on it', strip);
  return view.coal.card;
}

async function runCoalescent() {
  const c = kept.c, ticket = ++view.tickets.coal;
  const demography = c.demography === 'growth' ? { kind: 'growth', alpha: 8 } : c.demography === 'bottleneck' ? { kind: 'bottleneck', start: 0.05, duration: 0.05, size: 0.04 } : { kind: 'constant' };
  const stop = working(view.coal, 'drawing a genealogy…');
  try {
    const r = await api.run('pop.coalescent', { samples: c.samples, theta: c.theta, rho: c.rho, demography, seed: c.seed, at: c.at });
    if (ticket !== view.tickets.coal) return;
    kept.coal = r;
    drawCoalescent();
  } catch (e) { if (ticket === view.tickets.coal) failed(view.coal, e); } finally { stop(); }
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

function drawCoalescent() {
  const r = kept.coal, c = kept.c, box = clear(view.coal.figs);
  const width = view.coal.card.clientWidth || 1200;
  const W = Math.max(380, Math.min(600, width / 2 - 40)), H = 340;
  const { svg, order } = treeSvg(r.tree, W, H);
  const tree = el('div', { class: 'fig' });
  tree.innerHTML = svg;
  // The spectrum against θ/i.
  const n = r.samples, xi = r.spectrum;
  const bars = [];
  for (let i = 1; i < n; i++) bars.push({ x0: i - 0.42, x1: i + 0.42, y: xi[i], colour: token('--accent'), gap: 0, title: `${fmt.number(xi[i], 3)} sites with ${i} derived` });
  const expect = [];
  for (let i = 1; i < n; i++) expect.push([i, c.theta / i]);
  const sfs = el('div', { class: 'fig' });
  sfs.innerHTML = chart({ width: 400, height: 210, x: [0.4, n - 0.4], y: [0, Math.max(c.theta, ...xi.slice(1, n)) * 1.1 || 1], xlabel: 'derived copies i', ylabel: 'sites', bars,
    series: [{ points: expect, colour: token('--ink-2'), dash: '4 3', label: 'θ / i' }], xticks: [...Array(Math.min(n - 1, 10)).keys()].map((k) => Math.round(1 + (k * (n - 2)) / Math.max(1, Math.min(n - 2, 9)))) });
  const D = r.tajima_d;
  const facts = [['θ from θ_H (Fay and Wu)', fmt.number(r.theta_h, 4)], ['T_MRCA here', `${fmt.number(r.tree.tmrca, 3)} (expected ${fmt.number(r.expected_tmrca, 3)} at constant size)`]];
  if (c.rho > 0) facts.push(['along the chromosome', `${r.trees} marginal trees after ${r.recombinations} recombinations`]);
  box.append(el('div', { class: 'fig-row' }, tree, el('div', { class: 'side' },
    metrics([[fmt.count(r.segregating), 'segregating sites S'], [fmt.number(r.theta_w, 4), 'θ from S (Watterson)'], [fmt.number(r.pi, 4), 'θ from π (Tajima)'],
      [D === null ? '–' : fmt.number(D, 3), 'Tajima\'s D', D !== null && Math.abs(D) > 2 ? 'strong' : '']]),
    kv(facts),
    el('div', { class: 'eyebrow', style: { marginTop: '8px' } }, 'the site frequency spectrum, against θ / i'), sfs)));
  // The marginal trees' T_MRCA along the chromosome; a click chooses where the tree is drawn.
  if (c.rho > 0 && r.tmrca_along?.length) {
    const pts = [];
    for (const [a, b, t] of r.tmrca_along) { pts.push([a, t], [b, t]); }
    const along = el('div', { class: 'fig along', title: 'click to draw the tree at that place' });
    along.innerHTML = chart({ width: Math.min(900, width - 60), height: 150, x: [0, 1], y: [0, Math.max(...pts.map((p) => p[1])) * 1.1], xlabel: 'position on the chromosome (click to draw its tree)', ylabel: 'T_MRCA',
      series: [{ points: pts, colour: token('--accent'), area: true }], marks: [{ x: c.at, y: 0, text: '▲' }] });
    along.addEventListener('click', (e) => {
      const rect = along.querySelector('svg').getBoundingClientRect();
      const fx = (e.clientX - rect.left - 56) / (rect.width - 70);
      if (fx >= 0 && fx <= 1) { c.at = fx; runCoalescent(); }
    });
    box.append(el('div', { class: 'block' }, el('div', { class: 'eyebrow' }, 'the genealogy changes along a recombining chromosome'), along));
  }
  // The haplotypes, in the order of the tree's leaves.
  if (r.haplotypes && r.haplotypes[0]?.length) {
    const S = r.haplotypes[0].length, rows = r.haplotypes.length;
    const cw = Math.max(2, Math.min(10, Math.floor(Math.min(900, width - 60) / S))), ch = 9;
    const canvas = el('canvas', { width: S * cw, height: rows * ch, class: 'dp' });
    const g = canvas.getContext('2d');
    order.forEach((leaf, k) => {
      const h = r.haplotypes[leaf];
      for (let j = 0; j < S; j++) { g.fillStyle = h[j] === '1' ? token('--accent') : token('--paper-3'); g.fillRect(j * cw, k * ch, cw - (cw > 3 ? 1 : 0), ch - 1); }
    });
    box.append(el('div', { class: 'block' }, el('div', { class: 'eyebrow' }, `the sample: ${rows} sequences × ${S} segregating sites, derived alleles dark, in the order of the tree's leaves`), el('div', { class: 'fig' }, canvas)));
  }
}

// ---- Wright–Fisher ----------------------------------------------------------------------------------------
function wrightFisherCard() {
  const w = kept.w;
  view.wfBtn = btn('Run', { icon: 'play', kind: 'primary', onClick: runWF });
  const strip = el('div', { class: 'ctl-strip' },
    num('N individuals', w, 'individuals', { min: 2, max: 100000, step: 10, title: 'the individuals of each deme (2N copies of the gene)' }),
    num('s selection', w, 's', { min: -1, max: 1, step: 0.005, title: 'the advantage of A₁A₁ over A₂A₂' }),
    num('h dominance', w, 'h', { min: 0, max: 1, step: 0.1, title: 'the heterozygote\'s share of the advantage' }),
    num('p₀ at the start', w, 'p0', { min: 0, max: 1, step: 0.05, title: 'the frequency of A₁ at the start' }),
    num('demes', w, 'demes', { min: 1, max: 20, title: 'populations, exchanging migrants as islands' }),
    num('m migration', w, 'm', { min: 0, max: 1, step: 0.01, title: 'the share of each deme replaced by migrants each generation' }),
    num('generations', w, 'generations', { min: 10, max: 100000, step: 100 }),
    el('span', { class: 'grow' }),
    btn('Other replicates', { icon: 'refresh', title: 'the same parameters, another seed', onClick: () => { if (!view.wfBtn.disabled) { w.seed++; runWF(); } } }),
    view.wfBtn);
  onEnter(strip, runWF);
  view.wf = simCard('Wright–Fisher: drift and selection', `${fmt.count(w.replicates)} populations followed forwards in time, against Kimura's diffusion`, strip);
  return view.wf.card;
}

async function runWF() {
  if (view.wfBtn.disabled) return;
  const ticket = ++view.tickets.wf;
  const stop = working(view.wf, `following ${fmt.count(kept.w.replicates)} populations…`, view.wfBtn);
  try {
    const r = await api.run('pop.wright_fisher', { ...kept.w, replicates: kept.w.replicates, draw: 60 });
    if (ticket !== view.tickets.wf) return;
    kept.wf = r;
    drawWF();
  } catch (e) { if (ticket === view.tickets.wf) failed(view.wf, e); } finally { stop(); }
}

function drawWF() {
  const r = kept.wf, box = clear(view.wf.figs);
  const width = view.wf.card.clientWidth || 1200;
  const T = r.generations, every = r.every;
  const series = r.paths.map((p, k) => ({ points: p.map((v, i) => [i * every, v]), colour: CATEGORIES[k % CATEGORIES.length], width: 0.9 }));
  series.push({ points: r.mean_p.map((v, i) => [i * every, v]), colour: token('--ink'), width: 2.2, label: `mean of ${fmt.count(r.replicates)}` });
  const c = el('div', { class: 'fig' });
  c.innerHTML = chart({ width: Math.max(420, Math.min(720, width - 480)), height: 290, x: [0, T], y: [0, 1], xlabel: 'generations', ylabel: 'frequency of A₁', series });
  const facts = [['lost', fmt.count(r.lost)], ['settled after', `${fmt.number(r.mean_settling, 4)} generations on average`]];
  if (r.fst !== null) facts.push(['F_ST among demes', `${fmt.number(r.fst, 3)}${r.island_fst !== null ? ` (the island model: ${fmt.number(r.island_fst, 3)})` : ''}`]);
  box.append(el('div', { class: 'fig-row' }, c, el('div', { class: 'side' },
    metrics([[fmt.percent(r.fixation, 1), `fixed: ${fmt.count(r.fixed)} of ${fmt.count(r.replicates)}`], [fmt.percent(r.kimura, 1), 'Kimura predicts']]),
    kv(facts),
    el('p', { class: 'note' }, 'Kimura\'s probability of fixation, u(p) = (1 − e^{−4Nsp}) / (1 − e^{−4Ns}) for genic selection (Kimura 1962), solves the diffusion\'s backward equation; with dominance it is a ratio of integrals, computed by quadrature. A neutral allele fixes with probability p₀.'))));
}

// ---- structure ----------------------------------------------------------------------------------------
function structureCard() {
  const st = kept.st;
  view.stBtn = btn('Simulate and infer', { icon: 'play', kind: 'primary', onClick: runStructure });
  const strip = el('div', { class: 'ctl-strip' },
    num('populations K', st, 'populations', { min: 2, max: 6 }),
    num('of each', st, 'each', { min: 2, max: 100, step: 5, title: 'unmixed individuals from each population' }),
    num('admixed', st, 'admixed', { min: 0, max: 200, step: 5, title: 'individuals with ancestors in several populations' }),
    num('F_ST', st, 'fst', { min: 0.001, max: 0.9, step: 0.01, title: 'how far the populations have diverged (Balding–Nichols)' }),
    num('sites', st, 'sites', { min: 200, max: 20000, step: 500 }),
    el('span', { class: 'grow' }),
    btn('Another sample', { icon: 'refresh', title: 'the same parameters, another seed', onClick: () => { if (!view.stBtn.disabled) { st.seed++; runStructure(); } } }),
    view.stBtn);
  onEnter(strip, runStructure);
  view.st = simCard('Structure: populations and their mixture', 'principal components and the EM of ancestry, each against the truth', strip);
  return view.st.card;
}

const mixColour = (weights) => {
  const hex = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16));
  const rgb = [0, 0, 0];
  weights.forEach((w, k) => { const c = hex(CATEGORIES[k % CATEGORIES.length]); for (let i = 0; i < 3; i++) rgb[i] += w * c[i]; });
  return `rgb(${rgb.map((v) => Math.round(v)).join(',')})`;
};

async function runStructure() {
  if (view.stBtn.disabled) return;
  const ticket = ++view.tickets.st;
  const stop = working(view.st, 'simulating the sample, then its principal components and the EM of ancestry…', view.stBtn);
  try {
    const r = await api.run('pop.structure', kept.st);
    if (ticket !== view.tickets.st) return;
    kept.structure = { r, st: { ...kept.st } };
    drawStructure();
  } catch (e) { if (ticket === view.tickets.st) failed(view.st, e); } finally { stop(); }
}

function drawStructure() {
  const { r, st } = kept.structure, box = clear(view.st.figs);
  const width = view.st.card.clientWidth || 1200;
  const coords = r.pca.coordinates;
  const xs = coords.map((p) => p[0]), ys = coords.map((p) => p[1]);
  const W = 380, H = 320, pad = 34;
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  const X = (v) => pad + ((v - x0) / Math.max(1e-9, x1 - x0)) * (W - 2 * pad), Y = (v) => H - pad - ((v - y0) / Math.max(1e-9, y1 - y0)) * (H - 2 * pad);
  let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}"><line x1="${pad}" y1="${H - pad}" x2="${W - pad}" y2="${H - pad}" stroke="${token('--hair-2')}"/><line x1="${pad}" y1="${pad}" x2="${pad}" y2="${H - pad}" stroke="${token('--hair-2')}"/>`;
  coords.forEach((p, i) => { s += `<circle cx="${X(p[0]).toFixed(1)}" cy="${Y(p[1]).toFixed(1)}" r="${r.admixed[i] ? 3.4 : 4.2}" fill="${mixColour(r.truth[i])}" fill-opacity="0.85" stroke="${r.admixed[i] ? token('--ink') : 'none'}" stroke-width="0.6"><title>individual ${i + 1}${r.admixed[i] ? ', admixed' : ''}</title></circle>`; });
  s += `<text class="tick-label" x="${W / 2}" y="${H - 8}" text-anchor="middle" style="font-style:italic">PC 1 (eigenvalue ${fmt.number(r.pca.values[0], 3)})</text><text class="tick-label" transform="translate(12 ${H / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">PC 2 (${fmt.number(r.pca.values[1], 3)})</text></svg>`;
  const pca = el('div', { class: 'fig' });
  pca.innerHTML = s;
  // Ancestry: the EM's bars above, the truth's below; components matched to the truth's by best agreement.
  const K = r.k, n = r.samples;
  const perm = matchComponents(r.admixture.q, r.truth, K, r.populations);
  const idx = [...Array(n).keys()].sort((a, b) => r.population[a] - r.population[b] || r.truth[b][r.population[b]] - r.truth[a][r.population[a]]);
  const BW = Math.max(360, Math.min(820, width - 480)), BH = 90;
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
  const bars = el('div', { class: 'fig' });
  bars.innerHTML = b;
  const err = mae / (n * Math.min(K, r.populations));
  box.append(el('div', { class: 'fig-row' }, pca, el('div', { class: 'side' },
    metrics([[fmt.count(n), 'individuals'], [fmt.count(r.sites), 'sites'], [fmt.number(err, 3), 'mean error of ancestry'], [fmt.count(r.admixture.iterations), 'steps of the EM']]),
    bars,
    el('p', { class: 'note' }, `${n} individuals at ${fmt.count(r.sites)} sites, from ${r.populations} populations diverged to F_ST ≈ ${st.fst} (Balding–Nichols), ${st.admixed} of them admixed. The principal components (Patterson, Price and Reich 2006) place each individual between the populations it descends from; the EM of ancestry (Tang et al. 2005, accelerated) recovers the proportions to a mean absolute error of ${fmt.number(err, 3)} per component, in ${r.admixture.iterations} steps.`))));
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
