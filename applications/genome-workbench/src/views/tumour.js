// The Tumour view: a cancer genome against its normal (genomics::somatic, ADR 0064). A tumour is a
// mixture: normal cells and the clones of a tree, each with its mutations and its copy number. From
// the reads of tumour and normal, the view recovers what the mixture hides:
//
//   the somatic mutations, by MuTect's two likelihood tests (Cibulskis et al. 2013);
//   purity and ploidy, and every segment's integer copy numbers, by ASCAT's fit (Van Loo et al. 2010),
//     drawn as its "sunrise", the distance to an integer solution over ploidy and purity;
//   the clones, from each mutation's cancer cell fraction, clustered, and the trees they allow;
//   the mutational processes, from the 96 trinucleotide classes, refitted to their signatures;
//   microsatellite instability, by MSIsensor's test site by site.
//
// The tumour is simulated on a real genome and labelled so; its truth is drawn beside every estimate.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import { genomeTracksSvg } from '../draw/tracks.js';
import * as fmt from '../format.js';
import { token, warm, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, seg, select, toggle, fld, card, metrics, callout, empty, table, working, busy, onEnter } from '../ui/kit.js';

const REGIONS = [[300000, '300 kb'], [600000, '600 kb'], [1200000, '1.2 Mb']];
const PURITIES = [0.3, 0.5, 0.7, 0.9];
const BURDENS = [[500, '500 per Mb'], [2000, '2,000 per Mb'], [5000, '5,000 per Mb']];
const DEPTHS = [[40, '40×'], [60, '60×'], [100, '100×']];
const SEGMENTS = [3, 7, 11];

// SBS-96's six substitution classes in softened conventional colours (COSMIC's blue, black, red,
// grey, green and pink).
const SBS6 = [['C>A', '--C'], ['C>G', '--ink-2'], ['C>T', '--T'], ['T>A', '--other'], ['T>C', '--A'], ['T>G', '--accent-2']];

const kept = {
  p: { length: 600000, purity: 0.7, mutations: 2000, subclone: 0.55, subsub: 0.3, wgd: false, segments: 7, msi: 0.3, tumour: 60, normal: 35, seed: 3 },
  c: { tumours: 40, max_rank: 6, mean: 400, seed: 1 },
  r: null, rc: null, error: null, errorC: null, running: false, runningC: false, site: null,
};
let view = null;

export function mount(section, ctx) {
  view = { ctx, section };
  const p = ctx.params;
  if (p.purity && PURITIES.includes(Number(p.purity))) kept.p.purity = Number(p.purity);
  if (/^\d{1,6}$/.test(p.seed || '')) kept.p.seed = Number(p.seed);
  view.controls = el('div', { class: 'controls' });
  view.stage = el('div', { class: 'stage' });
  section.append(
    el('p', { class: 'lede' }, 'A tumour sample is a mixture: normal cells, and the clones of a tree, each carrying the mutations of its ancestors and its own copy number. Its reads therefore show every change diluted. A mutation present in every cancer cell, on one of two copies, is read at half the purity; one in a subclone, at less. ASCAT recovers purity and ploidy by asking which pair turns every segment\'s measured depth and allele balance into whole numbers of copies. With them, each mutation\'s fraction of cancer cells follows, the clones are found as clusters of those fractions, and the tree is constrained by the rule that children cannot outnumber their parent. The mutations\' trinucleotide contexts record the processes that made them.'),
    el('div', { class: 'bench' }, view.controls, view.stage));
  buildControls();
  draw();
  if (!kept.r && !kept.running) grow();
}

export function show(ctx) { view.ctx = ctx; }

function buildControls() {
  const p = kept.p, c = kept.c;
  view.runBtn = btn('Grow, read and call', { icon: 'play', kind: 'primary', onClick: grow });
  if (kept.running) busy(view.runBtn, true, 'Growing…');
  const again = btn('Another', { icon: 'refresh', title: 'another tumour of the same kind: the next seed', onClick: () => { p.seed++; grow(); } });
  const tumour = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'The tumour'),
      fld('purity', seg(PURITIES.map((v) => [v, fmt.percent(v, 0)]), p.purity, (v) => { p.purity = Number(v); }, { block: true }), { hint: 'the fraction of cells in the sample that are cancer cells' }),
      el('div', { class: 'ctl-row' },
        fld('mutations', select(BURDENS, p.mutations, (v) => { p.mutations = Number(v); })),
        fld('copy-number segments', select(SEGMENTS.map((v) => [v, String(v)]), p.segments, (v) => { p.segments = Number(v); }))),
      el('div', { class: 'ctl-row' },
        fld('subclone', select([[0, 'none'], [0.3, '30%'], [0.55, '55%'], [0.75, '75%']], p.subclone, (v) => { p.subclone = Number(v); }), { help: 'of cancer cells' }),
        fld('its subclone', select([[0, 'none'], [0.2, '20%'], [0.3, '30%']], p.subsub, (v) => { p.subsub = Number(v); }))),
      toggle('whole-genome doubling', p.wgd, (v) => { p.wgd = v; }),
      fld('microsatellite slippage', seg([[0, 'stable'], [0.3, 'unstable']], p.msi, (v) => { p.msi = Number(v); }, { block: true }))),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Its reads'),
      el('div', { class: 'ctl-row' },
        fld('tumour', select(DEPTHS, p.tumour, (v) => { p.tumour = Number(v); })),
        fld('normal', select([[30, '30×'], [35, '35×'], [50, '50×']], p.normal, (v) => { p.normal = Number(v); }))),
      fld('reference', select(REGIONS.map(([v, l]) => [v, `Escherichia coli K-12, ${l}`]), p.length, (v) => { p.length = Number(v); }), { hint: 'a real genome, so that its repeats and homopolymers are real; the tumour on it is simulated' })),
    el('div', { class: 'run-bar' }, view.runBtn, again));
  onEnter(tumour, grow);
  view.cohortBtn = btn('Extract', { icon: 'play', kind: 'primary', onClick: cohort });
  if (kept.runningC) busy(view.cohortBtn, true, 'Extracting…');
  const sigs = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Signatures across a cohort'),
      fld('tumours', seg([[20, '20'], [40, '40'], [80, '80']], c.tumours, (v) => { c.tumours = Number(v); }, { block: true })),
      fld('mutations each', select([[100, '100'], [400, '400'], [1500, '1,500']], c.mean, (v) => { c.mean = Number(v); }),
        { hint: 'Catalogues drawn from known signatures; non-negative matrix factorisation must find how many there are, and which, from the catalogues alone.' })),
    el('div', { class: 'run-bar' }, view.cohortBtn));
  clear(view.controls).append(tumour, sigs);
}

async function grow() {
  if (kept.running) return;
  kept.running = true; kept.error = null;
  busy(view.runBtn, true, 'Growing…');
  view.ctx.setParams({ purity: kept.p.purity, seed: kept.p.seed });
  draw();
  const p = kept.p;
  try {
    kept.r = await api.run('somatic.tumour', {
      genome: 'ncbi:NC_000913.3', contig: 0, begin: 0, length: p.length, seed: p.seed, purity: p.purity, mutations_per_mb: p.mutations,
      subclone_ccf: p.subclone, subsubclone_ccf: p.subclone ? p.subsub : 0, segments: p.segments, whole_genome_doubling: p.wgd, msi_rate: p.msi,
      tumour_depth: p.tumour, normal_depth: p.normal, bootstrap: 100 });
    kept.site = null;
  } catch (e) {
    kept.error = e.message;
  } finally {
    kept.running = false;
    busy(view.runBtn, false);
    draw();
  }
}

async function cohort() {
  if (kept.runningC) return;
  kept.runningC = true; kept.errorC = null;
  busy(view.cohortBtn, true, 'Extracting…');
  draw();
  const c = kept.c;
  try {
    kept.rc = await api.run('somatic.signatures', { seed: c.seed, tumours: c.tumours, max_rank: c.max_rank, mean_mutations: c.mean, replicates: 20 });
  } catch (e) {
    kept.errorC = e.message;
  } finally {
    kept.runningC = false;
    busy(view.cohortBtn, false);
    draw();
  }
}

// ---- drawing ------------------------------------------------------------------------------------------
function draw() {
  const box = clear(view.stage);
  if (kept.running) box.append(el('div', { class: 'results-head' }, working('Growing the tumour, reading tumour and normal, mapping, calling, fitting…')));
  if (kept.error) box.append(callout(kept.error, 'bad'));
  const r = kept.r;
  if (r) box.append(...tumourCards(r));
  else if (!kept.running && !kept.error) box.append(card({ body: empty({ icon: 'tumour', title: 'A tumour and its normal', text: 'Choose its purity and its clones, then grow it.' }) }));
  if (kept.runningC) box.append(el('div', { class: 'results-head', style: { marginTop: '18px' } }, working('Drawing a cohort\'s catalogues and factorising them at every rank…')));
  if (kept.errorC) box.append(callout(kept.errorC, 'bad'));
  if (kept.rc) box.append(cohortCard(kept.rc));
}

const W = () => Math.max(700, (view.stage?.clientWidth || 1000) - 34);

function tumourCards(r) {
  const out = [];
  const ev = r.evaluation || {}, fit = r.copy_number?.fit || {}, t = r.truth || {};
  out.push(el('div', { class: 'results-head' }, el('span', { class: 'badge warn' }, 'simulated'),
    el('span', { class: 'when' }, `${r.reference?.title || ''}, ${fmt.length(r.reference?.length || 0)} · tumour ${r.reads?.tumour_depth ?? ''}×, normal ${r.reads?.normal_depth ?? ''}× · ${fmt.seconds(Object.values(r.timing || {}).reduce((a, b) => a + (b || 0), 0))}`)));
  out.push(metrics([
    [fmt.percent(fit.purity, 0), `purity (truth ${fmt.percent(t.purity, 0)})`],
    [fmt.number(fit.ploidy, 3), `ploidy (truth ${fmt.number(t.ploidy, 3)})`],
    [fmt.percent(ev.snv?.precision, 1), `precision of substitutions`],
    [fmt.percent(ev.snv?.recall, 1), `recall of substitutions`],
    [String(r.clonality?.clones?.length ?? '–'), `clones found (truth ${(t.clones || []).length})`],
    [fmt.number(r.burden?.per_megabase, 3), 'mutations per megabase'],
    [r.msi ? (r.msi.high ? 'MSI-high' : 'stable') : '–', r.msi ? `${r.msi.somatic} of ${r.msi.tested} sites unstable` : 'microsatellites'],
  ]));
  out.push(profileCard(r));
  out.push(el('div', { class: 'two', style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', gap: '16px' } }, sunriseCard(r), clonesCard(r)));
  out.push(signatureCard(r));
  if (r.msi) out.push(msiCard(r));
  out.push(evaluationCard(r));
  return out;
}

// The copy-number profile: logR and BAF with their segments, the integer copy numbers, the truth behind.
function profileCard(r) {
  const cn = r.copy_number, L = r.reference.length, contig = r.reference.contig || 'contig';
  const bins = cn.bins, sites = cn.sites;
  const lr = bins.logr.filter(Number.isFinite);
  const lo = Math.min(-1.2, ...lr.map((v) => Math.max(-3, v))), hi = Math.max(1.2, ...lr.map((v) => Math.min(3, v)));
  const faint = token('--ink-faint'), acc = token('--accent'), blue = token('--C'), ghost = token('--hair-3');
  const top = Math.max(5, ...cn.segments.map((s) => s.major + 1), ...(r.truth?.segments || []).map((s) => s.major + 1));
  const tracks = [
    { kind: 'points', label: 'log R', height: 110, y: [lo, hi], lines: [0], ticks: [lo, 0, hi], opacity: 0.45,
      points: bins.begin.map((b, i) => ({ contig, pos: b + cn.bin / 2, value: bins.logr[i], colour: faint })),
      segments: cn.segments.map((s) => ({ contig, start: s.begin, end: s.end, value: s.logr, colour: acc, width: 2.2, title: `log R ${fmt.number(s.logr, 3)} over ${s.sites} sites` })) },
    { kind: 'points', label: 'B-allele fraction', height: 100, y: [0, 1], lines: [0.5], ticks: [0, 0.5, 1], opacity: 0.5,
      points: sites.position.map((x, i) => ({ contig, pos: x, value: sites.baf[i], colour: blue })),
      segments: cn.segments.flatMap((s) => [s.baf, 1 - s.baf].map((v) => ({ contig, start: s.begin, end: s.end, value: v, colour: acc, width: 1.8, title: `mirrored BAF ${fmt.number(s.baf, 3)}` }))) },
    { kind: 'segments', label: 'copies', height: 96, y: [0, top], lines: [2], ticks: [0, 2, 4],
      segments: [
        ...(r.truth?.segments || []).flatMap((s) => [
          { contig, start: s.begin, end: s.end, value: s.major + 0.1, colour: ghost, width: 6, title: `truth: ${s.major} + ${s.minor}` },
          { contig, start: s.begin, end: s.end, value: s.minor - 0.1, colour: ghost, width: 6, title: `truth: ${s.major} + ${s.minor}` }]),
        ...cn.segments.flatMap((s) => [
          { contig, start: s.begin, end: s.end, value: s.major + 0.1, colour: acc, width: 2.4, title: `major ${s.major} (${fmt.number(s.n_a_raw, 3)} before rounding)` },
          { contig, start: s.begin, end: s.end, value: s.minor - 0.1, colour: blue, width: 2.4, title: `minor ${s.minor} (${fmt.number(s.n_b_raw, 3)} before rounding)` }]),
      ] },
  ];
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = genomeTracksSvg({ contigs: [{ name: contig, length: L }], width: W(), window: { contig, start: 0, end: L, offset: r.reference.begin || 0 }, tracks });
  return card({ title: 'Copy number, allele by allele', sub: 'depth against the normal (log R) and the B-allele fraction at the normal\'s heterozygous sites, each with its segments; then ASCAT\'s integer copies, major in mulberry and minor in blue, the truth in grey behind', pad: true, body: fig });
}

// ASCAT's sunrise: the distance to integer copy numbers over ploidy ψ and purity ρ.
function sunriseCard(r) {
  const s = r.copy_number.sunrise, fit = r.copy_number.fit, t = r.truth;
  const nx = s.psi.length, ny = s.rho.length;
  const w = 430, h = 300, m = { l: 46, r: 10, t: 10, b: 36 };
  const cw = (w - m.l - m.r) / nx, ch = (h - m.t - m.b) / ny;
  const ds = s.distance.filter((v) => Number.isFinite(v) && v > 0);
  const dlo = Math.log(Math.min(...ds)), dhi = Math.log(Math.max(...ds));
  const out = [];
  for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
    const d = s.distance[i * ny + j];
    const v = Number.isFinite(d) && d > 0 ? 1 - (Math.log(d) - dlo) / (dhi - dlo || 1) : 0;
    out.push(`<rect x="${(m.l + i * cw).toFixed(2)}" y="${(m.t + (ny - 1 - j) * ch).toFixed(2)}" width="${(cw + 0.3).toFixed(2)}" height="${(ch + 0.3).toFixed(2)}" fill="${warm(v ** 1.6)}" />`);
  }
  const X = (psi) => m.l + ((psi - s.psi[0]) / (s.psi[nx - 1] - s.psi[0])) * (w - m.l - m.r);
  const Y = (rho) => m.t + (1 - (rho - s.rho[0]) / (s.rho[ny - 1] - s.rho[0])) * (h - m.t - m.b);
  for (const o of s.optima || []) out.push(`<circle cx="${X(o.psi).toFixed(1)}" cy="${Y(o.rho).toFixed(1)}" r="3.2" fill="none" stroke="${token('--ink-2')}" stroke-width="1"><title>${esc(`a local optimum: purity ${fmt.percent(o.rho, 0)}, ploidy ${fmt.number(o.ploidy, 3)}, goodness ${fmt.percent(o.goodness / 100, 1)}`)}</title></circle>`);
  if (fit.found) out.push(`<circle cx="${X(fit.psi).toFixed(1)}" cy="${Y(fit.purity).toFixed(1)}" r="5.5" fill="none" stroke="${token('--surface')}" stroke-width="2.4" /><circle cx="${X(fit.psi).toFixed(1)}" cy="${Y(fit.purity).toFixed(1)}" r="5.5" fill="none" stroke="${token('--accent-deep')}" stroke-width="1.4" />`);
  if (t) out.push(`<path d="M${X(t.ploidy) - 5},${Y(t.purity)}h10M${X(t.ploidy)},${Y(t.purity) - 5}v10" stroke="${token('--C')}" stroke-width="1.6" />`);
  for (const v of [1.5, 2, 3, 4, 5]) if (v >= s.psi[0] && v <= s.psi[nx - 1]) out.push(`<text class="tick-label" x="${X(v).toFixed(1)}" y="${h - m.b + 14}" text-anchor="middle">${v}</text>`);
  for (const v of [0.2, 0.4, 0.6, 0.8, 1]) if (v >= s.rho[0] && v <= s.rho[ny - 1]) out.push(`<text class="tick-label" x="${m.l - 6}" y="${(Y(v) + 3.5).toFixed(1)}" text-anchor="end">${v}</text>`);
  out.push(`<text class="tick-label" x="${(m.l + w - m.r) / 2}" y="${h - 4}" text-anchor="middle" style="font-style:italic">ploidy ψ</text>`);
  out.push(`<text class="tick-label" transform="translate(12 ${(m.t + h - m.b) / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">purity ρ</text>`);
  const fig = el('div');
  fig.innerHTML = `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">${out.join('')}</svg>`;
  return card({ title: 'Purity and ploidy, by ASCAT', sub: `the distance from integer copy numbers, darker where nearer; the fit circled (goodness ${fmt.percent((fit.goodness || 0) / 100, 1)}), the truth a blue cross`, pad: true, body: fig });
}

// The clones: each call's cancer cell fraction against its allele fraction, the clusters, the tree.
function clonesCard(r) {
  const calls = (r.calls || []).filter((c) => Number.isFinite(c.ccf) && !c.filters?.length);
  const clusters = r.clonality?.clones || [];
  const series = clusters.map((cl, k) => ({ points: calls.filter((c) => c.cluster === k).map((c) => [c.ccf, c.vaf]), colour: CATEGORIES[k % CATEGORIES.length], line: false, dots: true, label: `clone at ${fmt.percent(cl.ccf, 0)} · ${cl.mutations}` }));
  const scatter = chart({ width: 430, height: 230, x: [0, 1.3], y: [0, Math.max(0.6, ...calls.map((c) => c.vaf))], xlabel: 'cancer cell fraction', ylabel: 'allele fraction', series });
  // A histogram of the fractions.
  const nb = 26, hist = new Array(nb).fill(0);
  for (const c of calls) hist[Math.min(nb - 1, Math.floor((Math.min(1.3, c.ccf) / 1.3) * nb))]++;
  const bars = hist.map((v, i) => ({ x0: (i / nb) * 1.3, x1: ((i + 1) / nb) * 1.3, y: v, colour: token('--accent-2') }));
  const truthMarks = (r.truth?.clones || []).filter((c) => c.ccf < 1.01).map((c) => ({ points: [[c.ccf, 0], [c.ccf, Math.max(...hist) * 1.05]], colour: token('--C'), dash: '3 3', width: 1.2 }));
  const histogram = chart({ width: 430, height: 160, x: [0, 1.3], y: [0, Math.max(...hist, 1) * 1.1], xlabel: 'cancer cell fraction (true clones dashed)', bars, series: truthMarks });
  const tree = treeSvg(r.clonality);
  const fig = el('div');
  fig.innerHTML = scatter + histogram + tree;
  const n = r.clonality?.trees_consistent;
  return card({ title: 'Its clones', sub: `fractions of cancer cells, clustered; ${n === 1 ? 'one tree' : `${n} trees`} of ${r.clonality?.trees_examined ?? '?'} keep the rule that a parent outnumbers its children`, pad: true, body: fig });
}

function treeSvg(cl) {
  if (!cl || !cl.best_tree) return '';
  const parent = cl.best_tree.parent || [], ccf = cl.tree_nodes_ccf || [], clone = cl.tree_nodes_clone || [];
  const n = parent.length;
  if (!n) return '';
  const kids = Array.from({ length: n }, () => []);
  let root = 0;
  parent.forEach((p, i) => { if (p < 0) root = i; else kids[p].push(i); });
  const depth = new Array(n).fill(0), order = [];
  const walk = (i, d) => { depth[i] = d; order.push(i); for (const k of kids[i]) walk(k, d + 1); };
  walk(root, 0);
  const leafX = new Map();
  let next = 0;
  const place = (i) => { if (!kids[i].length) { leafX.set(i, next++); return leafX.get(i); } const xs = kids[i].map(place); const x = (Math.min(...xs) + Math.max(...xs)) / 2; leafX.set(i, x); return x; };
  place(root);
  const w = 430, h = 120, dx = (w - 60) / Math.max(1, next - 1 || 1), dy = 34;
  const X = (i) => 30 + leafX.get(i) * (next > 1 ? dx : 0) + (next > 1 ? 0 : (w - 60) / 2);
  const Y = (i) => 22 + depth[i] * dy;
  const out = [];
  for (const i of order) for (const k of kids[i]) out.push(`<path d="M${X(i)},${Y(i)} L${X(k)},${Y(k)}" stroke="${token('--hair-3')}" stroke-width="1.4" />`);
  for (const i of order) {
    const c = clone[i];
    const col = c >= 0 ? CATEGORIES[c % CATEGORIES.length] : token('--ink-faint');
    const rr = 5 + 9 * Math.sqrt(Math.max(0, ccf[i] || 0));
    out.push(`<circle cx="${X(i)}" cy="${Y(i)}" r="${rr.toFixed(1)}" fill="${col}" fill-opacity="0.75" stroke="${token('--surface')}" stroke-width="1.5"><title>${esc(c >= 0 ? `clone at ${fmt.percent(ccf[i], 0)} of cancer cells` : 'the trunk')}</title></circle>`);
    out.push(`<text class="tick-label" x="${X(i) + rr + 4}" y="${Y(i) + 4}">${esc(fmt.percent(ccf[i], 0))}</text>`);
  }
  return `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">${out.join('')}</svg>`;
}

// SBS-96 with its reconstruction, and the exposures with their bootstrap intervals.
function sbsSvg(values, { recon = null, width = 760, height = 170, title = '' } = {}) {
  const m = { l: 40, r: 8, t: title ? 22 : 10, b: 36 };
  const n = values.length, bw = (width - m.l - m.r) / n;
  const max = Math.max(...values, ...(recon || []), 1e-9);
  const H = height - m.t - m.b;
  const out = [];
  if (title) out.push(`<text x="${m.l}" y="14" style="font: italic 12.5px var(--serif); fill:${token('--ink-2')}">${esc(title)}</text>`);
  for (let i = 0; i < n; i++) {
    const g = Math.floor(i / 16), col = token(SBS6[g][1]);
    const hh = (values[i] / max) * H;
    out.push(`<rect x="${(m.l + i * bw + 0.6).toFixed(1)}" y="${(m.t + H - hh).toFixed(1)}" width="${(bw - 1.2).toFixed(1)}" height="${Math.max(0, hh).toFixed(1)}" fill="${col}" fill-opacity="0.85" />`);
  }
  if (recon) out.push(`<path d="${recon.map((v, i) => `${i ? 'L' : 'M'}${(m.l + (i + 0.5) * bw).toFixed(1)},${(m.t + H - (v / max) * H).toFixed(1)}`).join('')}" fill="none" stroke="${token('--ink')}" stroke-width="1" stroke-opacity="0.65" />`);
  for (let g = 0; g < 6; g++) {
    const x0 = m.l + g * 16 * bw, x1 = x0 + 16 * bw;
    out.push(`<rect x="${(x0 + 1).toFixed(1)}" y="${m.t + H + 4}" width="${(x1 - x0 - 2).toFixed(1)}" height="5" rx="1.5" fill="${token(SBS6[g][1])}" />`);
    out.push(`<text class="tick-label" x="${((x0 + x1) / 2).toFixed(1)}" y="${m.t + H + 22}" text-anchor="middle">${SBS6[g][0]}</text>`);
  }
  out.push(`<line x1="${m.l}" y1="${m.t + H}" x2="${width - m.r}" y2="${m.t + H}" stroke="${token('--hair-2')}" />`);
  return `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">${out.join('')}</svg>`;
}

function signatureCard(r) {
  const cat = r.catalogue, sg = r.signatures, t = r.truth;
  const names = (t.signatures || []).map((s) => s.name);
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = sbsSvg(cat.sbs96, { recon: sg?.refit?.reconstruction, width: Math.min(900, W() + 160), title: `${cat.substitutions} substitutions in their 96 trinucleotide classes; the line is their refit (cosine ${fmt.number(sg?.refit?.cosine, 4)})` });
  // Exposures: estimated with intervals beside the truth.
  const ex = sg?.refit?.exposures || [], bs = sg?.bootstrap || [], te = sg?.truth_exposures || [];
  const total = ex.reduce((a, b) => a + b, 0) || 1, ttotal = te.reduce((a, b) => a + b, 0) || 1;
  const rows = names.map((nm, i) => el('div', { class: 'change', style: { gridTemplateColumns: '210px 1fr' } },
    el('span', { class: 'what' }, nm),
    el('span', {}, el('b', {}, `${fmt.count(Math.round(ex[i] || 0))} mutations`), ` · ${fmt.percent((ex[i] || 0) / total, 1)}`,
      bs[i] ? el('span', { class: 'hint' }, ` (bootstrap ${fmt.count(Math.round(bs[i].low))}–${fmt.count(Math.round(bs[i].high))})`) : '',
      el('span', { class: 'hint' }, `  ·  truth ${fmt.count(Math.round(te[i] || 0))}, ${fmt.percent((te[i] || 0) / ttotal, 1)}`))));
  return card({ title: 'The processes that made its mutations', sub: 'the catalogue refitted to known signatures by non-negative least squares, with bootstrap intervals; the signatures here are synthetic, built from their mechanisms (COSMIC\'s may not be redistributed)', pad: true,
    body: [fig, el('div', { class: 'change-list', style: { marginTop: '10px' } }, ...rows)] });
}

function msiCard(r) {
  const m = r.msi;
  const sites = m.sites || [];
  const L = r.reference.length;
  const pts = sites.filter((s) => s.tested).map((s) => [s.position, Math.min(12, -Math.log10(Math.max(1e-12, s.p)))]);
  const series = [
    { points: pts.filter((_, i) => !sites.filter((s) => s.tested)[i].somatic), colour: token('--ink-faint'), line: false, dots: true },
    { points: pts.filter((_, i) => sites.filter((s) => s.tested)[i].somatic), colour: token('--T'), line: false, dots: true, label: 'unstable' },
  ];
  const strip = chart({ width: Math.min(900, W() + 160), height: 190, x: [0, L], y: [0, 12.5], xlabel: 'position', ylabel: '−log₁₀ p', series });
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = strip;
  // The most unstable site's lengths, normal against tumour.
  const top = sites.filter((s) => s.tested).sort((a, b) => a.p - b.p)[0];
  if (top) {
    const lens = [...new Set([...Object.keys(top.normal || {}), ...Object.keys(top.tumour || {})].map(Number))].sort((a, b) => a - b);
    const nt = Object.values(top.normal).reduce((a, b) => a + b, 0) || 1, tt = Object.values(top.tumour).reduce((a, b) => a + b, 0) || 1;
    const bars = [];
    lens.forEach((len) => {
      bars.push({ x0: len - 0.42, x1: len - 0.02, y: (top.normal[len] || 0) / nt, colour: token('--C') });
      bars.push({ x0: len + 0.02, x1: len + 0.42, y: (top.tumour[len] || 0) / tt, colour: token('--T') });
    });
    const h = el('div');
    h.innerHTML = chart({ width: 430, height: 180, x: [lens[0] - 1, lens[lens.length - 1] + 1], y: [0, 1], xlabel: `repeat length at ${top.motif}×${top.copies} (normal blue, tumour rose)`, ylabel: 'fraction of reads', bars });
    fig.append(h);
  }
  return card({ title: 'Microsatellites', sub: `${m.tested} sites tested by MSIsensor's χ² of repeat lengths, tumour against normal; ${m.somatic} unstable, ${fmt.percent(m.tested ? m.somatic / m.tested : 0, 1)}: ${m.high ? 'MSI-high' : 'stable'}`, pad: true, body: fig });
}

function evaluationCard(r) {
  const ev = r.evaluation || {};
  const bands = ev.recall_by_vaf || [];
  const bars = bands.map((b) => ({ x0: b.low, x1: b.high, y: b.truth ? b.found / b.truth : 0, colour: token('--accent-2'), title: `${b.found} of ${b.truth}` }));
  const fig = el('div');
  fig.innerHTML = chart({ width: 430, height: 190, x: [0, Math.max(0.6, ...bands.map((b) => b.high))], y: [0, 1.05], xlabel: 'true allele fraction', ylabel: 'found', bars });
  const rows = [
    ['substitutions', ev.snv], ['insertions and deletions', ev.indel],
  ].filter(([, x]) => x).map(([k, x]) => ({ kind: k, tp: x.tp, fp: x.fp, fn: x.fn, precision: fmt.percent(x.precision, 1), recall: fmt.percent(x.recall, 1) }));
  return card({ title: 'Against the truth', sub: `copy numbers right on ${fmt.percent(ev.copy_number_bases_right, 1)} of bases; purity off by ${fmt.number(ev.purity_error, 2)}, ploidy by ${fmt.number(ev.ploidy_error, 2)}`, pad: true,
    body: el('div', { style: { display: 'grid', gridTemplateColumns: 'minmax(0, 1fr) 440px', gap: '16px', alignItems: 'start' } },
      table([{ key: 'kind', label: 'calls' }, { key: 'tp', label: 'right', num: true }, { key: 'fp', label: 'false', num: true }, { key: 'fn', label: 'missed', num: true }, { key: 'precision', label: 'precision' }, { key: 'recall', label: 'recall' }], rows), fig) });
}

function cohortCard(rc) {
  const ranks = rc.ranks || [];
  const stab = chart({ width: 430, height: 210, x: [1, Math.max(2, ...ranks.map((x) => x.rank))], y: [0, 1.05], xlabel: 'rank', ylabel: 'stability · error',
    series: [
      { points: ranks.map((x) => [x.rank, x.stability]), colour: token('--accent'), dots: true, label: 'stability' },
      { points: ranks.map((x) => [x.rank, x.relative_error]), colour: token('--C'), dots: true, label: 'relative error' },
      { points: [[1, 0.8], [ranks.length + 1, 0.8]], colour: token('--hair-3'), dash: '3 3', width: 1 },
    ] });
  const fig = el('div');
  fig.innerHTML = stab;
  const pairs = el('div', { class: 'stack', style: { gap: '4px' } });
  (rc.match || []).forEach((mt, k) => {
    const truth = (rc.truth || []).find((t) => t.name === mt.name) || rc.truth?.[mt.truth];
    const box = el('div', { style: { overflowX: 'auto' } });
    box.innerHTML = sbsSvg(rc.extracted[k], { width: 430, height: 110, title: `extracted ${k + 1} ↔ ${mt.name}, cosine ${fmt.number(mt.cosine, 4)}` }) + (truth ? sbsSvg(truth.profile, { width: 430, height: 100, title: 'its truth' }) : '');
    pairs.append(box);
  });
  return card({ title: 'Signatures found from catalogues alone', sub: `${rc.truth?.length ?? '?'} signatures behind ${rc.exposures?.length ?? '?'} tumours; the rank chosen: ${rc.chosen_rank}, the largest whose solutions are stable across resamplings (SigProfiler's criterion)`, pad: true,
    body: el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(440px, 1fr))', gap: '14px' } }, fig, pairs) });
}
