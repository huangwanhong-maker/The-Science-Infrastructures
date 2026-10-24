// The Trees view: genomes and genes placed on the tree of their descent. A gene taken from each
// chosen genome is aligned (progressive multiple alignment), its distances corrected for
// multiple substitutions, and trees built from them by UPGMA and neighbour joining, or from the
// alignment itself by maximum likelihood under a substitution model with Γ rates; the bootstrap
// says how firmly each split is supported. Between every pair, the gene's codons say how selection
// has acted on it: dN/dS by Nei and Gojobori. Whole genomes are compared by their sketches (Mash
// distances), which need no alignment at all.

import { el, clear, esc, notice, toggles, field, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, warm, CATEGORIES } from '../palette.js';

const PRIMATES = ['ncbi:NC_012920.1', 'ncbi:NC_011137.1', 'ncbi:NC_001643.1', 'ncbi:NC_001644.1', 'ncbi:NC_011120.1', 'ncbi:NC_002083.1', 'ncbi:NC_001646.1', 'ncbi:NC_002082.1', 'ncbi:NC_005943.1', 'ncbi:NC_020006.2'];
const GENES = [['CYTB', 'cytochrome b'], ['COX1', 'cytochrome c oxidase I'], ['ND2', 'NADH dehydrogenase 2'], ['ND5', 'NADH dehydrogenase 5'], ['ATP6', 'ATP synthase 6']];

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, controls: el('div', { class: 'align-inputs' }), main: el('div'), selection: el('div'), sketches: el('div'),
           s: { genomes: new Set(PRIMATES), gene: 'CYTB', protein: false, method: 'ml', model: 'HKY85', gamma: true, bootstrap: 100 }, result: null };
  section.append(
    el('h2', { class: 'view-title' }, 'Trees'),
    el('p', { class: 'lede' }, 'A tree with lengths on its branches is a metric on its leaves, and from the distances alone neighbour joining finds it when the distances are those of a tree (Saitou and Nei 1987). Distances seen in sequences are not: substitutions hide one another, which a model corrects. A Markov chain of substitutions along the branches makes the tree a probability model of the alignment, whose likelihood Felsenstein\'s pruning computes as a sum over every ancestral state; maximum likelihood chooses the tree and its lengths, and the bootstrap, resampling columns, says which splits the data hold.'),
    el('div', { class: 'align-grid' }, view.controls, view.main),
    el('h3', { class: 'section-title' }, 'Selection on the gene'), view.selection,
    el('h3', { class: 'section-title' }, 'Whole genomes by their sketches'), view.sketches);
  buildControls();
  buildSketches();
}

export function show(ctx) {
  view.ctx = ctx;
  buildControls();
  if (!view.result) build();
}

function buildControls() {
  const s = view.s, box = clear(view.controls);
  const list = el('div', { class: 'stack small' });
  for (const g of view.ctx.catalogue.filter((x) => x.group === 'mitochondria')) {
    list.append(el('label', { class: 'field' }, el('input', { type: 'checkbox', checked: s.genomes.has(g.id), disabled: !g.fetched, onchange: (e) => { if (e.target.checked) s.genomes.add(g.id); else s.genomes.delete(g.id); } }),
      el('span', { style: { fontStyle: 'normal', color: 'var(--ink)' } }, g.title.replace(/ mitochondrion.*$/, '').replace(/^Homo sapiens$/, 'Homo sapiens (human)'))));
  }
  const gene = el('select', { onchange: () => { s.gene = gene.value; } }, ...GENES.map(([g, name]) => el('option', { value: g, selected: g === s.gene }, `${g}, ${name}`)));
  const model = el('select', { onchange: () => { s.model = model.value; } }, ...['JC69', 'K80', 'F81', 'HKY85', 'TN93', 'GTR'].map((m) => el('option', { value: m, selected: m === s.model }, m)));
  const boots = el('select', { onchange: () => { s.bootstrap = Number(boots.value); } }, ...[0, 50, 100, 200].map((b) => el('option', { value: b, selected: b === s.bootstrap }, b ? `${b} replicates` : 'none')));
  box.append(
    el('p', { class: 'eyebrow' }, 'the genomes (mitochondria)'), list,
    el('div', { class: 'row' }, field('gene', gene)),
    el('p', { class: 'faint small' }, 'Its bases: the models here are of nucleotides (JC69 to GTR); models of amino acids and codons are not yet written.'),
    el('p', { class: 'eyebrow' }, 'the tree'),
    el('div', { class: 'row' }, toggles([['upgma', 'UPGMA'], ['nj', 'neighbour joining'], ['ml', 'maximum likelihood']], s.method, (v) => { s.method = v; })),
    el('div', { class: 'row' }, field('model', model), el('label', { class: 'field' }, el('input', { type: 'checkbox', checked: s.gamma, onchange: (e) => { s.gamma = e.target.checked; } }), '+Γ, four rates')),
    el('div', { class: 'row' }, field('bootstrap', boots)),
    el('div', { class: 'row' }, el('button', { class: 'button primary', onclick: build }, 'Build the tree')));
}

async function build() {
  const s = view.s;
  if (s.genomes.size < 3) { notice('Choose at least three genomes.', true); return; }
  clear(view.main).append(el('p', { class: 'wait' }, `Aligning ${s.gene} from ${s.genomes.size} genomes, then building the tree${s.bootstrap ? ` and ${s.bootstrap} bootstrap replicates across the workers` : ''}…`));
  try {
    selection();
    view.result = await api.run('trees.gene', { genomes: [...s.genomes], gene: s.gene, protein: false, method: s.method, model: s.model, gamma: s.gamma, bootstrap: s.bootstrap, seed: 1 });
    draw();
  } catch (e) { clear(view.main).append(el('p', { class: 'error' }, e.message)); }
}

// ---- selection: dN/dS between every pair ---------------------------------------------------------------
// Each pair's codons, aligned through their proteins; the synonymous and non-synonymous sites and
// differences, Jukes–Cantor corrected (Nei and Gojobori 1986). Neutral change gives ω = dN/dS = 1;
// purifying selection, which removes most changes of the protein, ω well below it.
async function selection() {
  const s = view.s, box = clear(view.selection);
  box.append(el('p', { class: 'wait' }, `Aligning the codons of ${s.gene} pair by pair…`));
  let r;
  try { r = await api.run('trees.dnds', { genomes: [...s.genomes], gene: s.gene }); } catch (e) { clear(box).append(el('p', { class: 'error' }, e.message)); return; }
  clear(box);
  const pairs = r.pairs.filter((q) => q.dS > 0 && q.dN > 0 && Number.isFinite(q.dS) && Number.isFinite(q.dN));
  const short = (n) => n.replace(/^(\w)\w+ /, '$1. ');
  if (pairs.length) {
    const xs = pairs.map((q) => q.dS), ys = pairs.map((q) => q.dN);
    const x = [Math.min(...xs) / 1.6, Math.max(...xs) * 1.6], y = [Math.min(...ys, ...xs.map((v) => v * 0.01)) / 1.6, Math.max(...ys, ...xs) * 1.6];
    const line = (w, label, dash) => ({ points: [[x[0], x[0] * w], [x[1], x[1] * w]], colour: token('--ink-faint'), dash, width: 1, label });
    const omegas = pairs.map((q) => q.omega).filter(Number.isFinite).sort((a, b) => a - b);
    const median = omegas[Math.floor(omegas.length / 2)];
    const c = el('div');
    c.innerHTML = chart({ width: 460, height: 340, x, y, log: 'xy', xlabel: 'dS, synonymous substitutions a synonymous site', ylabel: 'dN, a non-synonymous site',
      series: [line(1, 'ω = 1, neutral', '5 3'), line(0.1, 'ω = 0.1', '2 3'), line(0.01, 'ω = 0.01', '1 4'),
        { points: pairs.map((q) => [q.dS, q.dN]), colour: token('--accent'), line: false, dots: true }] });
    const table = el('table', { class: 'data' }, el('tr', {}, ...['pair', 'codons', 'dS', 'dN', 'ω'].map((h) => el('th', { class: h === 'pair' ? '' : 'num' }, h))),
      ...[...r.pairs].sort((a, b) => a.dS - b.dS).map((q) => el('tr', {}, el('td', { class: 'small' }, `${short(r.names[q.a])} · ${short(r.names[q.b])}`), el('td', { class: 'num' }, fmt.count(q.codons)),
        el('td', { class: 'num' }, q.dS === null ? 'saturated' : fmt.number(q.dS, 3)), el('td', { class: 'num' }, q.dN === null ? 'saturated' : fmt.number(q.dN, 3)),
        el('td', { class: 'num' }, q.omega === null ? '–' : fmt.number(q.omega, 2)))));
    box.append(el('p', { class: 'note' }, `${r.gene} under the ${r.code} code (table ${r.table}), its codons aligned through the proteins. Of the possible changes at a codon, those that keep the amino acid are synonymous sites, the rest non-synonymous; the differences of each kind over its sites, corrected for multiple hits, are dS and dN (Nei and Gojobori 1986). The median ω here is ${fmt.number(median, 2)}: most changes of the protein were removed by selection before they could be seen, and the closest pairs show the highest ω, the trace of slightly harmful changes that selection has not yet had time to remove.`),
      el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '28px' } }, c, el('div', { style: { maxHeight: '340px', overflowY: 'auto' } }, table)));
  } else box.append(el('p', { class: 'note' }, 'No pair has both kinds of difference.'));
}

// ---- drawing a tree ----------------------------------------------------------------------------------------
// nodes: [{parent, length, label, support}], the root's parent −1. Rectangular, root at the left.
export function treeSvg(nodes, { width = 620, rowHeight = 22, labelWidth = 230, colour = () => token('--ink-2'), scaleLabel = 'substitutions a site' } = {}) {
  const n = nodes.length;
  const kids = Array.from({ length: n }, () => []);
  let root = 0;
  nodes.forEach((d, v) => { if (d.parent < 0) root = v; else kids[d.parent].push(v); });
  const depth = new Array(n).fill(0), y = new Array(n).fill(0);
  const order = [];
  const walk = (v) => {
    for (const c of kids[v]) { depth[c] = depth[v] + Math.max(0, nodes[c].length || 0); walk(c); }
    if (!kids[v].length) { y[v] = order.length; order.push(v); }
    else y[v] = (y[kids[v][0]] + y[kids[v][kids[v].length - 1]]) / 2;
  };
  walk(root);
  const maxDepth = Math.max(1e-9, ...depth);
  const W = width, H = order.length * rowHeight + 46, left = 14, right = labelWidth;
  const X = (d) => left + (d / maxDepth) * (W - left - right), Y = (k) => 16 + k * rowHeight;
  let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`;
  nodes.forEach((d, v) => {
    if (d.parent < 0) return;
    s += `<path d="M${X(depth[d.parent]).toFixed(1)},${Y(y[v]).toFixed(1)}H${X(depth[v]).toFixed(1)}" stroke="${colour(v)}" stroke-width="1.6" fill="none" />`;
  });
  nodes.forEach((d, v) => {
    if (!kids[v].length) return;
    const ys = kids[v].map((c) => Y(y[c]));
    s += `<path d="M${X(depth[v]).toFixed(1)},${Math.min(...ys).toFixed(1)}V${Math.max(...ys).toFixed(1)}" stroke="${token('--ink-2')}" stroke-width="1.6" />`;
    if (d.support !== undefined && d.support !== null && d.parent >= 0) {
      const sup = d.support;
      s += `<text class="tick-label" x="${(X(depth[v]) - 4).toFixed(1)}" y="${(Y(y[v]) - 5).toFixed(1)}" text-anchor="end" style="fill:${sup >= 0.95 ? token('--good') : sup >= 0.7 ? token('--ink-2') : token('--bad')}">${Math.round(sup * 100)}</text>`;
    }
  });
  order.forEach((v) => {
    s += `<circle cx="${X(depth[v]).toFixed(1)}" cy="${Y(y[v]).toFixed(1)}" r="3" fill="${colour(v)}" />`;
    s += `<text x="${(X(depth[v]) + 8).toFixed(1)}" y="${(Y(y[v]) + 4).toFixed(1)}" style="font-size:13px; font-style:italic; fill:${token('--ink')}">${esc(nodes[v].label || '')}</text>`;
  });
  // A scale bar.
  const bar = Math.pow(10, Math.floor(Math.log10(maxDepth / 3)));
  const by = H - 16;
  s += `<line x1="${left}" y1="${by}" x2="${X(bar)}" y2="${by}" stroke="${token('--ink')}" stroke-width="1.4" />`;
  s += `<text class="tick-label" x="${X(bar) + 6}" y="${by + 4}">${fmt.number(bar, 2)} ${esc(scaleLabel)}</text></svg>`;
  return { svg: s, order };
}

function draw() {
  const r = view.result, box = clear(view.main);
  const { svg, order } = treeSvg(r.tree.nodes, { width: Math.min(680, view.section.clientWidth - 460) });
  const fig = el('div');
  fig.innerHTML = svg;
  const facts = el('dl', { class: 'kv' });
  const row = (k, v) => facts.append(el('dt', {}, k), el('dd', {}, v));
  row('the alignment', `${r.alignment.rows.length} sequences × ${fmt.count(r.alignment.columns)} columns, ${fmt.count(r.alignment.variable)} variable, ${fmt.count(r.alignment.informative)} parsimony-informative`);
  row('aligned by', r.alignment.method);
  if (r.log_likelihood !== undefined) row('log-likelihood', fmt.number(r.log_likelihood, 7));
  if (r.parameters) row('model', Object.entries(r.parameters).map(([k, v]) => `${k} ${fmt.number(v, 3)}`).join(' · '));
  row('tree length', `${fmt.number(r.tree.length, 4)} substitutions a site`);
  if (r.compare) for (const [k, v] of Object.entries(r.compare)) row(`Robinson–Foulds to ${k}`, String(v));
  if (r.bootstrap) row('bootstrap', `${r.bootstrap.replicates} replicates in ${fmt.seconds(r.bootstrap.seconds)}${r.bootstrap.shards ? `, ${r.bootstrap.shards} shards` : ''}`);
  box.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '24px' } }, fig, el('div', { style: { maxWidth: '420px' } }, facts,
    el('p', { class: 'note' }, 'Support at a node is the fraction of bootstrap trees that contain its split: green at 95% or more, red below 70%.'))));
  // The distances, in the order of the tree's leaves.
  const D = r.distances, names = r.alignment.names;
  const leafIndex = order.map((v) => names.indexOf(r.tree.nodes[v].label));
  if (D && leafIndex.every((k) => k >= 0)) {
    const m = leafIndex.length, cell = Math.max(10, Math.min(30, Math.floor(420 / m)));
    const max = Math.max(...D.flat());
    const canvas = el('canvas', { width: m * cell, height: m * cell, class: 'dp' });
    const g = canvas.getContext('2d');
    leafIndex.forEach((a, i) => leafIndex.forEach((b, j) => { g.fillStyle = warm(1 - D[a][b] / Math.max(1e-9, max)); g.fillRect(j * cell, i * cell, cell, cell); }));
    canvas.addEventListener('pointermove', (e) => {
      const rect = canvas.getBoundingClientRect();
      const i = Math.floor((e.clientY - rect.top) / cell), j = Math.floor((e.clientX - rect.left) / cell);
      if (i < 0 || j < 0 || i >= m || j >= m) return;
      hover(`${esc(names[leafIndex[i]])} — ${esc(names[leafIndex[j]])}<br>distance <b>${fmt.number(D[leafIndex[i]][leafIndex[j]], 4)}</b>`, e);
    });
    canvas.addEventListener('pointerleave', () => hover(null));
    box.append(el('p', { class: 'eyebrow' }, `distances (${r.distance_model}), in the order of the tree; darker is closer`), canvas);
  }
  // The alignment, coloured by base or residue, its variable columns marked.
  if (r.alignment.rows.length && r.alignment.columns <= 6000) {
    const rows = r.alignment.rows, cols = r.alignment.columns;
    const cw = Math.max(1, Math.min(4, Math.floor((view.section.clientWidth - 100) / cols))), ch = 12;
    const canvas = el('canvas', { width: cols * cw, height: rows.length * ch, class: 'dp' });
    const g = canvas.getContext('2d');
    const colours = { A: token('--A'), C: token('--C'), G: token('--G'), T: token('--T') };
    order.forEach((v, k) => {
      const row = rows[names.indexOf(r.tree.nodes[v].label)];
      if (!row) return;
      for (let j = 0; j < cols; j++) {
        const c = row[j];
        const variable = r.alignment.variable_columns ? r.alignment.variable_columns[j] : true;
        g.fillStyle = c === '-' ? token('--paper-2') : r.protein ? (variable ? CATEGORIES['ARNDCQEGHILKMFPSTWYV'.indexOf(c) % CATEGORIES.length] || token('--N') : token('--paper-4')) : (variable ? colours[c] || token('--N') : token('--paper-4'));
        g.fillRect(j * cw, k * ch, cw, ch - 1);
      }
    });
    box.append(el('p', { class: 'eyebrow' }, 'the alignment, rows in the order of the tree: variable columns coloured, constant ones pale'), el('div', { style: { overflowX: 'auto' } }, canvas));
  }
}

// ---- sketches ----------------------------------------------------------------------------------------
function buildSketches() {
  const box = clear(view.sketches);
  const out = el('div');
  box.append(el('p', { class: 'note' }, 'A genome\'s k-mers, hashed, and the smallest thousand of them kept: two such sketches estimate the Jaccard index of the k-mer sets, and the Mash distance turns it into a rate of mutation (Ondov et al. 2016). Every fetched genome is sketched; the distances, computed as a task graph of blocks, give a tree by neighbour joining.'),
    el('button', { class: 'button', onclick: async () => {
      clear(out).append(el('p', { class: 'wait' }, 'Sketching every fetched genome…'));
      try {
        const ids = view.ctx.catalogue.filter((g) => g.fetched && g.group !== 'simulated').map((g) => g.id);
        const r = await api.run('trees.sketches', { genomes: ids, k: 21, size: 1000 });
        clear(out);
        const { svg } = treeSvg(r.tree.nodes, { width: Math.min(760, view.section.clientWidth - 60), scaleLabel: 'Mash distance' });
        const fig = el('div');
        fig.innerHTML = svg;
        out.append(el('p', { class: 'note' }, `${r.genomes} genomes sketched (the ${fmt.count(r.size)} smallest hashes of their canonical ${r.k}-mers) in ${fmt.seconds(r.sketch_seconds)}; ${fmt.count(r.pairs)} distances in ${fmt.seconds(r.distance_seconds)}, as ${r.tasks} tasks of blocks. Beyond a distance of about 0.3 the sketches share almost nothing and the distances saturate at 1: sketches tell close relatives apart, and say little about how far apart the distant ones are.`), fig);
      } catch (e) { clear(out).append(el('p', { class: 'error' }, e.message)); }
    } }, 'Sketch the genomes'), out);
}
