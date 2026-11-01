// The Association view: from genotype to trait and to expression (evolution::association and tyche;
// ADR 0064). A simulated cohort, clearly labelled so, whose genotypes carry real linkage (founders
// from the coalescent, individuals as mosaics of them), analysed as a study would be:
//
//   the scan       every variant tested with principal components as covariates; genomic control;
//   the locus      the strongest signal fine-mapped by SuSiE (individual data and summary statistics)
//                  and by Wakefield's ABF, credible sets against the truth;
//   expression     cis-eQTLs with FastQTL's beta approximation and q-values; colocalisation; expression
//                  predicted from genotypes (PrediXcan) and its association with the trait (TWAS);
//   the whole      heritability by REML and Haseman–Elston, and polygenic scores on held-out people.
//
// Every estimate stands beside the truth that made it.

import { el, clear, esc } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, CATEGORIES } from '../palette.js';
import { btn, seg, select, fld, card, metrics, callout, empty, table, working, busy, onEnter } from '../ui/kit.js';

const kept = { p: { individuals: 1500, heritability: 0.4, causal: 20, genes: 12, seed: 1 }, r: null, error: null, running: false };
let view = null;

export function mount(section, ctx) {
  view = { ctx, section };
  if (/^\d{1,6}$/.test(ctx.params.seed || '')) kept.p.seed = Number(ctx.params.seed);
  view.controls = el('div', { class: 'controls' });
  view.stage = el('div', { class: 'stage' });
  section.append(
    el('p', { class: 'lede' }, 'Which variants make a trait, and through which genes? A scan tests every variant, but linkage makes its neighbours rise with it, and the ancestry of a cohort makes whole regions rise if it is not removed. Fine-mapping asks which variants carry the signal, as sets that cover them with a stated probability. Expression is itself a trait of genotype: its cis-eQTLs are mapped, the genes whose expression and the trait share a variant are found by colocalisation, and expression is predicted from genotype alone, so that a gene can be tested in people whose expression was never measured. Every number here stands beside the truth of the simulation that made it.'),
    el('div', { class: 'bench' }, view.controls, view.stage));
  buildControls();
  draw();
}

export function show(ctx) { view.ctx = ctx; }

function buildControls() {
  const p = kept.p;
  view.runBtn = btn('Simulate and analyse', { icon: 'play', kind: 'primary', onClick: study });
  if (kept.running) busy(view.runBtn, true, 'Analysing…');
  const again = btn('Another', { icon: 'refresh', title: 'another cohort: the next seed', onClick: () => { p.seed++; study(); } });
  const c = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'The cohort'),
      fld('people', seg([[1000, '1,000'], [1500, '1,500'], [3000, '3,000']], p.individuals, (v) => { p.individuals = Number(v); }, { block: true }), { hint: '500 of them with measured expression; the rest for the scan' }),
      fld('heritability', seg([[0.2, '0.2'], [0.4, '0.4'], [0.6, '0.6']], p.heritability, (v) => { p.heritability = Number(v); }, { block: true })),
      el('div', { class: 'ctl-row' },
        fld('causal variants', select([[5, '5'], [20, '20'], [80, '80']], p.causal, (v) => { p.causal = Number(v); })),
        fld('genes', select([[8, '8'], [12, '12'], [30, '30']], p.genes, (v) => { p.genes = Number(v); }), { help: 'with cis-eQTLs' }))),
    el('div', { class: 'run-bar' }, view.runBtn, again),
    el('p', { class: 'hint', style: { margin: '0 14px 12px' } }, 'A minute or two: predicting expression trains an elastic net for every gene, with nested cross-validation.'));
  onEnter(c, study);
  clear(view.controls).append(c);
}

async function study() {
  if (kept.running) return;
  kept.running = true; kept.error = null;
  busy(view.runBtn, true, 'Analysing…');
  view.ctx.setParams({ seed: kept.p.seed });
  draw();
  const p = kept.p;
  try {
    kept.r = await api.run('association.study', { individuals: p.individuals, heritability: p.heritability, causal: p.causal, genes: p.genes, regions: 4, seed: p.seed, permutations: 300 });
  } catch (e) {
    kept.error = e.message;
  } finally {
    kept.running = false;
    busy(view.runBtn, false);
    draw();
  }
}

const W = () => Math.max(700, (view.stage?.clientWidth || 1000) - 34);

function draw() {
  const box = clear(view.stage);
  if (kept.running) box.append(el('div', { class: 'results-head' }, working('Simulating the cohort, scanning, fine-mapping, mapping eQTLs, predicting expression, estimating heritability…')));
  if (kept.error) box.append(callout(kept.error, 'bad'));
  const r = kept.r;
  if (!r) { if (!kept.running && !kept.error) box.append(card({ body: empty({ icon: 'association', title: 'A cohort and its trait', text: 'Choose its size and the trait\'s heritability, then simulate it and analyse it as a study would.' }) })); return; }
  const h = r.heritability || {}, g = r.gwas || {};
  box.append(el('div', { class: 'results-head' }, el('span', { class: 'badge warn' }, 'simulated'),
    el('span', { class: 'when' }, `${fmt.count(r.individuals)} people in ${r.populations} populations, ${fmt.count(r.sites)} variants in ${r.regions.length} regions · ${fmt.seconds(Object.values(r.seconds || {}).reduce((a, b) => a + b, 0))}`)));
  box.append(metrics([
    [fmt.number(g.lambda, 3), `λ, genomic control (${fmt.number(g.lambda_no_pcs, 3)} without PCs)`],
    [`${g.causal_found} of ${r.causal}`, 'causal variants significant'],
    [fmt.number(h.greml, 3), `heritability by REML (truth ${fmt.number(h.true, 3)})`],
    [`${r.eqtl?.egenes ?? '–'} of ${r.eqtl?.genes?.length ?? '–'}`, 'genes with an eQTL'],
    [fmt.number(Math.max(...(r.prs?.ct || []).map((x) => x.r2)), 3), 'R² of the best polygenic score'],
  ]));
  box.append(manhattanCard(r));
  box.append(el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', gap: '16px' } }, qqCard(r), pcsCard(r)));
  box.append(locusCard(r));
  box.append(el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', gap: '16px' } }, eqtlCard(r), predixcanCard(r)));
  box.append(el('div', { style: { display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', gap: '16px' } }, heritabilityCard(r), prsCard(r)));
}

function manhattanCard(r) {
  const rows = r.gwas.manhattan;  // [chrom, pos, mlog10p, causal, mlog10p_logistic]
  const L = r.regions.map((x) => x.length);
  const offs = L.map((_, i) => L.slice(0, i).reduce((a, b) => a + b, 0));
  const total = L.reduce((a, b) => a + b, 0);
  const top = Math.max(10, ...rows.map((x) => x[2])) * 1.05;
  const width = W(), height = 260, m = { l: 50, r: 12, t: 12, b: 34 };
  const X = (c, p) => m.l + ((offs[c] + p) / total) * (width - m.l - m.r);
  const Y = (v) => m.t + (1 - v / top) * (height - m.t - m.b);
  const out = [];
  const gw = -Math.log10(5e-8);
  out.push(`<line x1="${m.l}" y1="${Y(gw)}" x2="${width - m.r}" y2="${Y(gw)}" stroke="${token('--T')}" stroke-dasharray="4 4" />`);
  out.push(`<text class="tick-label" x="${width - m.r}" y="${Y(gw) - 4}" text-anchor="end">5 × 10⁻⁸</text>`);
  for (const v of [0, 10, 20, 30, 40].filter((v) => v <= top)) out.push(`<text class="tick-label" x="${m.l - 6}" y="${Y(v) + 3.5}" text-anchor="end">${v}</text>`);
  for (const x of rows) {
    const col = CATEGORIES[[1, 3, 5, 2][x[0] % 4]];
    out.push(`<circle cx="${X(x[0], x[1]).toFixed(1)}" cy="${Y(x[2]).toFixed(1)}" r="${x[3] ? 3.6 : 1.7}" fill="${x[3] ? 'none' : col}" stroke="${x[3] ? token('--ink') : 'none'}" stroke-width="1.2" fill-opacity="0.75" />`);
  }
  r.regions.forEach((reg, c) => out.push(`<text class="tick-label" x="${X(c, reg.length / 2)}" y="${height - 12}" text-anchor="middle">region ${reg.name}</text>`));
  out.push(`<text class="tick-label" transform="translate(13 ${(m.t + height - m.b) / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">−log₁₀ p</text>`);
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">${out.join('')}</svg>`;
  return card({ title: 'The scan', sub: 'every variant tested, with sex, age and four principal components as covariates; the truly causal variants ringed; the dashed line is genome-wide significance', pad: true, body: fig });
}

function qqCard(r) {
  const q = r.gwas.qq, q0 = r.gwas.qq_no_pcs || [];
  const maxE = Math.max(...q.map((x) => x[0]));
  const maxO = Math.min(40, Math.max(...q.map((x) => x[1]), ...q0.map((x) => x[1])));
  const band = q.filter((x) => x.length >= 4);
  const fig = el('div');
  fig.innerHTML = chart({ width: 420, height: 300, x: [0, maxE * 1.02], y: [0, maxO * 1.05], xlabel: 'expected −log₁₀ p', ylabel: 'observed −log₁₀ p', series: [
    { points: band.map((x) => [x[0], x[3]]), colour: token('--hair-3'), dash: '2 3', width: 1 },
    { points: band.map((x) => [x[0], x[2]]), colour: token('--hair-3'), dash: '2 3', width: 1 },
    { points: [[0, 0], [maxE, maxE]], colour: token('--ink-faint'), width: 1 },
    { points: q0.map((x) => [x[0], Math.min(maxO, x[1])]), colour: token('--C'), line: false, dots: true, label: `without PCs, λ ${fmt.number(r.gwas.lambda_no_pcs, 3)}` },
    { points: q.map((x) => [x[0], Math.min(maxO, x[1])]), colour: token('--accent'), line: false, dots: true, label: `with PCs, λ ${fmt.number(r.gwas.lambda, 3)}` }] }).replace(/r="3"/g, 'r="2"');
  return card({ title: 'Observed against expected', sub: 'the p-values against their law under no association, with the 95% band; ancestry not removed lifts the whole line', pad: true, body: fig });
}

function pcsCard(r) {
  const pts = r.pcs.points;
  const pops = [...new Set(pts.map((x) => x[2]))];
  const fig = el('div');
  fig.innerHTML = chart({ width: 420, height: 300, x: [Math.min(...pts.map((x) => x[0])), Math.max(...pts.map((x) => x[0]))], y: [Math.min(...pts.map((x) => x[1])), Math.max(...pts.map((x) => x[1]))], xlabel: 'PC 1', ylabel: 'PC 2',
    series: pops.map((p, k) => ({ points: pts.filter((x) => x[2] === p).map((x) => [x[0], x[1]]), colour: CATEGORIES[k % CATEGORIES.length], line: false, dots: true, label: `population ${p + 1}` })) }).replace(/r="3"/g, 'r="1.8" fill-opacity="0.6"');
  return card({ title: 'The cohort\'s ancestry', sub: `the first two principal components of the genotypes (eigenvalues ${r.pcs.values.slice(0, 2).map((v) => fmt.number(v, 3)).join(', ')})`, pad: true, body: fig });
}

function locusCard(r) {
  const L = r.locus, col = (n) => L.columns.indexOf(n);
  const pos = col('position'), lp = col('mlog10p'), r2 = col('r2_lead'), ps = col('pip_susie'), pa = col('pip_abf'), pr = col('pip_rss'), cs = col('causal');
  const width = W(), m = { l: 56, r: 12 };
  const X = (p) => m.l + ((p - L.from) / (L.to - L.from)) * (width - m.l - m.r);
  const out = [];
  const top = Math.max(...L.sites.map((s) => s[lp])) * 1.08;
  const H1 = 190, Y1 = (v) => 10 + (1 - v / top) * (H1 - 20);
  // r² to the lead, in five steps from sage to rose, as locus zooms colour it
  const r2col = (v) => ['#A8B3C2', '#8FB09B', '#C9A043', '#D2786B', '#B5463E'][Math.min(4, Math.floor(v * 5))];
  for (const s of L.sites) out.push(`<circle cx="${X(s[pos]).toFixed(1)}" cy="${Y1(s[lp]).toFixed(1)}" r="${s[cs] ? 4 : 2.6}" fill="${r2col(s[r2])}" stroke="${s[cs] ? token('--ink') : 'none'}" stroke-width="1.2" fill-opacity="0.85" />`);
  out.push(`<text class="tick-label" x="${m.l - 8}" y="${Y1(top / 2)}" text-anchor="end" style="font-style:italic">−log₁₀ p</text>`);
  // PIPs, three tracks
  const tracks = [['SuSiE', ps, token('--accent')], ['SuSiE-RSS', pr, token('--C')], ['ABF', pa, token('--G')]];
  let y = H1 + 8;
  for (const [name, k, c] of tracks) {
    const h = 44;
    out.push(`<line x1="${m.l}" y1="${y + h}" x2="${width - m.r}" y2="${y + h}" stroke="${token('--hair-2')}" />`);
    out.push(`<text class="tick-label" x="${m.l - 8}" y="${y + h / 2 + 4}" text-anchor="end">${name}</text>`);
    for (const s of L.sites) if (s[k] > 0.01) out.push(`<line x1="${X(s[pos]).toFixed(1)}" y1="${y + h}" x2="${X(s[pos]).toFixed(1)}" y2="${(y + h - s[k] * h).toFixed(1)}" stroke="${c}" stroke-width="2" />`);
    y += h + 8;
  }
  // credible sets as bars
  (L.susie?.sets || []).forEach((set, k) => {
    const xs = set.positions.map(X);
    out.push(`<rect x="${Math.min(...xs) - 3}" y="${y}" width="${Math.max(6, Math.max(...xs) - Math.min(...xs) + 6)}" height="8" rx="2" fill="${CATEGORIES[k % CATEGORIES.length]}" fill-opacity="${set.covers ? 0.85 : 0.35}"><title>${esc(`credible set ${k + 1}: ${set.positions.length} variants, purity ${fmt.number(set.purity, 3)}, coverage ${fmt.number(set.coverage, 3)}${set.covers ? ', holds a causal variant' : ''}`)}</title></rect>`);
  });
  y += 20;
  out.push(`<text class="tick-label" x="${m.l}" y="${y + 6}">${esc(fmt.length(L.from))}</text><text class="tick-label" x="${width - m.r}" y="${y + 6}" text-anchor="end">${esc(fmt.length(L.to))}</text>`);
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = `<svg width="${width}" height="${y + 14}" viewBox="0 0 ${width} ${y + 14}">${out.join('')}</svg>`;
  const s = L.susie || {};
  return card({ title: 'The strongest locus, fine-mapped', sub: `coloured by r² to the lead variant, causal variants ringed; beneath, each variant's posterior probability of being causal by three methods, and SuSiE's 95% credible sets: ${s.covering} of ${(s.sets || []).length} hold a causal variant, ${s.causal_found} of the ${s.causal_in_region} causal variants of the region found`, pad: true, body: fig });
}

function eqtlCard(r) {
  const e = r.eqtl, co = new Map((r.coloc || []).map((c) => [c.gene, c]));
  const rows = e.genes.map((g) => ({ gene: g.name, causal: g.causal, h2: fmt.number(g.cis_h2, 2), p: g.p_beta < 1e-4 ? g.p_beta.toExponential(1) : fmt.number(g.p_beta, 3), q: fmt.number(g.q, 3),
    egene: g.egene ? 'yes' : '', lead: g.egene ? (g.lead_causal ? 'causal' : 'linked') : '', h4: co.get(g.name) ? fmt.percent(co.get(g.name).pp[4], 0) : '' }));
  return card({ title: 'Genes and their eQTLs', sub: `the lead variant of each gene's cis-window, its p by FastQTL's beta approximation of ${e.permutations} permutations, its q-value (π₀ ${fmt.number(e.pi0, 3)}); H4, colocalisation with the trait`, pad: true,
    body: table([{ key: 'gene', label: 'gene' }, { key: 'causal', label: 'causal', num: true }, { key: 'h2', label: 'cis h²', num: true }, { key: 'p', label: 'p', num: true }, { key: 'q', label: 'q', num: true }, { key: 'egene', label: 'eGene' }, { key: 'lead', label: 'lead' }, { key: 'h4', label: 'H4', num: true }], rows) });
}

function predixcanCard(r) {
  const pg = r.predixcan.genes.filter((g) => g.kept);
  const sc = r.predixcan.scatter;
  const fig = el('div');
  const a = chart({ width: 420, height: 260, x: [0, Math.max(0.1, ...r.predixcan.genes.map((g) => g.cis_h2)) * 1.1], y: [Math.min(0, ...r.predixcan.genes.map((g) => g.r2_nested)), Math.max(0.1, ...r.predixcan.genes.map((g) => g.r2_nested)) * 1.1],
    xlabel: 'true cis-heritability', ylabel: 'held-out R² of the prediction', series: [
      { points: [[0, 0], [1, 1]], colour: token('--hair-3'), dash: '3 3', width: 1 },
      { points: r.predixcan.genes.map((g) => [g.cis_h2, g.r2_nested]), colour: token('--accent'), line: false, dots: true }] });
  const b = sc?.points ? chart({ width: 420, height: 260, x: [Math.min(...sc.points.map((p) => p[0])), Math.max(...sc.points.map((p) => p[0]))], y: [Math.min(...sc.points.map((p) => p[1])), Math.max(...sc.points.map((p) => p[1]))],
    xlabel: `${sc.gene}, predicted from genotype`, ylabel: 'measured expression', series: [{ points: sc.points.map((p) => [p[0], p[1]]), colour: token('--C'), line: false, dots: true }] }).replace(/r="3"/g, 'r="1.6" fill-opacity="0.5"') : '';
  fig.innerHTML = a + b;
  const twas = pg.filter((g) => Number.isFinite(g.twas_z)).map((g) => ({ gene: g.name, z: fmt.number(g.twas_z, 3), sz: fmt.number(g.spredixcan_z, 3), r2: fmt.number(g.r2_nested, 3), mediating: g.mediating ? 'yes' : '' }));
  return card({ title: 'Expression predicted from genotype', sub: `an elastic net per gene (α ${r.predixcan.alpha}) from the variants within ${fmt.length(r.predixcan.window)}; its held-out R² against the truth it can reach; then the trait tested on the predicted expression (TWAS), from individuals and from summary statistics (S-PrediXcan)`, pad: true,
    body: [fig, twas.length ? table([{ key: 'gene', label: 'gene' }, { key: 'r2', label: 'R²', num: true }, { key: 'z', label: 'TWAS z', num: true }, { key: 'sz', label: 'S-PrediXcan z', num: true }, { key: 'mediating', label: 'mediates' }], twas) : null] });
}

function heritabilityCard(r) {
  const h = r.heritability;
  const pr = h.profile;
  const maxL = Math.max(...pr.map((x) => x[1]));
  const series = [{ points: pr.map((x) => [x[0], x[1] - maxL]), colour: token('--accent'), label: 'REML\'s profile' },
    { points: [[h.true, -12], [h.true, 0.5]], colour: token('--C'), dash: '3 3', width: 1.4, label: `truth ${fmt.number(h.true, 3)}` },
    { points: [[h.he, -12], [h.he, 0.5]], colour: token('--G'), dash: '2 3', width: 1.2, label: `Haseman–Elston ${fmt.number(h.he, 3)}` }];
  const fig = el('div');
  fig.innerHTML = chart({ width: 420, height: 260, x: [0, 1], y: [-12, 1], xlabel: 'h²', ylabel: 'log likelihood from its best', series });
  return card({ title: 'Heritability', sub: `GREML ${fmt.number(h.greml, 3)} ± ${fmt.number(h.greml_se, 2)}, Haseman–Elston ${fmt.number(h.he, 3)} ± ${fmt.number(h.he_se, 2)} on ${fmt.count(h.samples)} people; of the case–control trait, ${fmt.number(h.case_control?.liability, 3)} on the liability scale`, pad: true, body: fig });
}

function prsCard(r) {
  const p = r.prs;
  const ct = p.ct.filter((x) => x.threshold > 0);
  const fig = el('div');
  fig.innerHTML = chart({ width: 420, height: 260, x: [Math.min(...ct.map((x) => x.threshold)), 1], y: [0, Math.max(0.05, ...ct.map((x) => x.r2), p.ldpred_inf) * 1.15], log: 'x', xlabel: 'p-value threshold', ylabel: 'R² on held-out people',
    series: [{ points: ct.map((x) => [x.threshold, x.r2]), colour: token('--accent'), dots: true, label: 'clumping and thresholding' },
      { points: [[ct[0].threshold, p.ldpred_inf], [1, p.ldpred_inf]], colour: token('--C'), dash: '4 3', width: 1.4, label: `LDpred-inf ${fmt.number(p.ldpred_inf, 3)}` },
      { points: [[ct[0].threshold, p.guide], [1, p.guide]], colour: token('--hair-3'), dash: '2 3', width: 1, label: `expected (Daetwyler) ${fmt.number(p.guide, 3)}` }] });
  return card({ title: 'Polygenic scores', sub: `trained on ${fmt.count(p.train)} people, tested on ${fmt.count(p.test)}; with few causal variants, thresholding does better than LDpred's infinitesimal prior, as Vilhjálmsson et al. (2015) found`, pad: true, body: fig });
}
