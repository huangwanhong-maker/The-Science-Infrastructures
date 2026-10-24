// The Trees view: genomes and genes placed on the tree of their descent. A gene taken from each
// chosen genome is aligned (progressive multiple alignment), its distances corrected for
// multiple substitutions, and trees built from them by UPGMA and neighbour joining, or from the
// alignment itself by maximum likelihood under a substitution model with Γ rates; the bootstrap
// says how firmly each split is supported. Between every pair, the gene's codons say how selection
// has acted on it: dN/dS by Nei and Gojobori. Whole genomes are compared by their sketches (Mash
// distances), which need no alignment at all.
//
// A bench: the genomes, the gene and the method at the left; the tree with its distances beside
// its leaves, the alignment and the selection on the gene at the right; below, the sketches. The
// address keeps the gene and the method (#trees&gene=COX1&method=nj).

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { token, warm, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, select, numberField, toggle, fld, card, metrics, callout, empty, kv, download, toast } from '../ui/kit.js';

const PRIMATES = ['ncbi:NC_012920.1', 'ncbi:NC_011137.1', 'ncbi:NC_001643.1', 'ncbi:NC_001644.1', 'ncbi:NC_011120.1', 'ncbi:NC_002083.1', 'ncbi:NC_001646.1', 'ncbi:NC_002082.1', 'ncbi:NC_005943.1', 'ncbi:NC_020006.2'];
const GENES = [['CYTB', 'cytochrome b'], ['COX1', 'cytochrome c oxidase I'], ['ND2', 'NADH dehydrogenase 2'], ['ND5', 'NADH dehydrogenase 5'], ['ATP6', 'ATP synthase 6']];
const METHODS = [
  ['upgma', 'UPGMA', 'Average linkage of the corrected distances: a molecular clock is assumed, and roots the tree.'],
  ['nj', 'Neighbour joining', 'The tree of the distances, exactly when they are a tree\'s (Saitou and Nei 1987); no clock, rooted at its midpoint.'],
  ['ml', 'Maximum likelihood', 'The tree and lengths under which the alignment is likeliest, by Felsenstein\'s pruning, searched from the joining tree.'],
];
const MODELS = [
  ['JC69', 'every substitution at one rate, the bases equally frequent (Jukes and Cantor 1969)'],
  ['K80', 'transitions and transversions each at a rate of their own (Kimura 1980)'],
  ['F81', 'one rate, the bases at their frequencies in the alignment (Felsenstein 1981)'],
  ['HKY85', 'K80\'s two rates with the bases at their frequencies (Hasegawa, Kishino and Yano 1985)'],
  ['TN93', 'the two kinds of transition and the transversions, each at its rate (Tamura and Nei 1993)'],
  ['GTR', 'each of the six exchanges at its own rate: the general time-reversible model (Tavaré 1986)'],
];
const BY = { upgma: 'by UPGMA', nj: 'by neighbour joining', ml: 'by maximum likelihood' };
const ROW = 22;   // the tree's row height, which the distances beside it share

let view = null;
// What was chosen and computed, kept when the view is mounted again (a change of theme redraws it).
const kept = { s: null, tree: null, dnds: null, sketch: null, sk: { k: 21, size: 1000 } };

export function mount(section, ctx) {
  kept.s ||= { genomes: new Set(PRIMATES), gene: 'CYTB', method: 'ml', model: 'HKY85', gamma: true, bootstrap: 100 };
  view = { ctx, section, s: kept.s, running: false, ticket: 0, dndsTicket: 0,
    picks: el('div', { class: 'picks' }), count: el('span', { class: 'count' }),
    head: el('div', { class: 'results-head' }), error: el('div'), main: el('div', { class: 'stack', style: { gap: '16px' } }), selection: el('div'),
    sketchBody: el('div') };
  view.stage = el('div', { class: 'stage' }, view.head, view.error, view.main, view.selection);
  section.append(
    el('p', { class: 'lede' }, 'A tree with lengths on its branches is a metric on its leaves, and from the distances alone neighbour joining finds it when the distances are those of a tree (Saitou and Nei 1987). Distances seen in sequences are not: substitutions hide one another, which a model corrects. A Markov chain of substitutions along the branches makes the tree a probability model of the alignment, whose likelihood Felsenstein\'s pruning computes as a sum over every ancestral state; maximum likelihood chooses the tree and its lengths, and the bootstrap, resampling columns, says which splits the data hold.'),
    el('div', { class: 'bench' }, controls(), view.stage),
    el('div', { class: 'section-head' }, el('h3', {}, 'Whole genomes, by their sketches'), el('span', { class: 'sub' }, 'Mash distances, which need no alignment')),
    sketchCard());
  drawPicks();
  if (kept.tree) drawTree();
  if (kept.dnds) drawSelection();
  if (kept.sketch) drawSketch();
}

export function show(ctx) {
  view.ctx = ctx;
  drawPicks();
  const changed = fromAddress(ctx.params);
  if (!kept.tree || changed) build();
}

// ---- the controls ----------------------------------------------------------------------------------
function controls() {
  const s = view.s;
  const gene = select(GENES.map(([g, name]) => [g, `${g}, ${name}`]), s.gene, (v) => { s.gene = v; });
  const model = select(MODELS.map(([m]) => [m, m]), s.model, (v) => { s.model = v; describe(); });
  const modelField = fld('substitution model', model, { hint: ' ' });
  const describe = () => { modelField.querySelector('.hint').textContent = MODELS.find(([m]) => m === s.model)?.[1] || ''; };
  describe();
  const gamma = toggle('+Γ, four categories of rate', s.gamma, (v) => { s.gamma = v; });
  gamma.hidden = s.method !== 'ml';
  const methods = choices(METHODS, s.method, (v) => { s.method = v; gamma.hidden = v !== 'ml'; });
  const boot = numberField({ value: s.bootstrap, min: 0, max: 1000, step: 50, unit: 'replicates', onChange: (v) => { s.bootstrap = v; } });
  view.inputs = { gene, model, methods, gamma, boot, describe };
  view.runBtn = btn('Build the tree', { icon: 'play', kind: 'primary', onClick: build });
  const box = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Genomes', view.count), view.picks),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Gene'),
      fld('the gene of each genome', gene, { hint: 'Its bases: the models here are of nucleotides (JC69 to GTR); models of amino acids and codons are not yet written.' })),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'Tree'), methods, modelField, gamma, fld('bootstrap', boot)),
    el('div', { class: 'run-bar' }, view.runBtn));
  box.addEventListener('keydown', (e) => {
    if (e.key !== 'Enter' || !(e.target.tagName === 'INPUT' || e.ctrlKey || e.metaKey)) return;
    e.preventDefault();
    if (e.target.tagName === 'INPUT') e.target.dispatchEvent(new Event('change'));
    build();
  });
  return el('div', { class: 'controls' }, box,
    quiet('These are RefSeq\'s mitochondria. With your own sequences, ', link('Build a tree', () => view.ctx.openTool('tree-build', toolValues())), ' in the Toolbox.'));
}

// The Toolbox's tree-build, opened with this view's choices and, once built, the gene's sequences.
function toolValues() {
  const s = view.s, v = { method: s.method, model: s.model, gamma: s.gamma, bootstrap: s.bootstrap, aligned: false };
  const a = kept.tree?.r.alignment;
  if (a) v.sequences = { text: a.names.map((n, k) => `>${n.replace(/\s+/g, '_')}\n${a.rows[k].replace(/-/g, '')}`).join('\n') + '\n' };
  return v;
}

// Exclusive choices with what each means: a column of cards, one of them chosen.
function choices(options, value, onPick) {
  const box = el('div', { class: 'choices', role: 'radiogroup' });
  const set = (v) => { for (const b of box.children) { const on = b.dataset.v === v; b.classList.toggle('on', on); b.setAttribute('aria-checked', String(on)); } };
  for (const [v, title, text] of options) {
    box.append(el('button', { type: 'button', class: 'choice', role: 'radio', 'data-v': v, onclick: () => { set(v); onPick(v); } },
      el('span', { class: 'ring' }), el('span', { class: 't' }, title), el('span', { class: 'd' }, text)));
  }
  set(value);
  box.set = set;
  return box;
}

// The genomes, as chips grouped as the catalogue groups them: the mitochondria, and your GenBank
// files (whose genes are named, as a gene tree needs).
const GROUPS = [['Mitochondria', (g) => g.group === 'mitochondria'], ['Your files', (g) => g.group === 'files' && g.kind === 'genbank']];
const speciesName = (g) => g.title.replace(/ mitochondrion.*$/, '').replace(/, complete genome.*$/, '').replace(/\.$/, '');

function drawPicks() {
  const s = view.s, box = clear(view.picks);
  for (const [label, test] of GROUPS) {
    const items = view.ctx.catalogue.filter(test);
    if (!items.length) continue;
    const ids = items.filter((g) => g.fetched).map((g) => g.id);
    const on = ids.filter((id) => s.genomes.has(id)).length;
    box.append(el('div', { class: 'pick-head' }, el('span', {}, label), el('span', { class: 'n' }, `${on} of ${items.length}`), el('span', { style: { flex: 1 } }),
      el('button', { class: 'link', disabled: on === ids.length || null, onclick: () => { ids.forEach((id) => s.genomes.add(id)); drawPicks(); } }, 'all'),
      el('button', { class: 'link', disabled: !on || null, onclick: () => { ids.forEach((id) => s.genomes.delete(id)); drawPicks(); } }, 'none')));
    const chips = el('div', { class: 'chips' });
    for (const g of items) {
      const chosen = s.genomes.has(g.id);
      chips.append(el('button', { type: 'button', class: `chip pick${chosen ? ' on' : ''}`, 'aria-pressed': String(chosen), disabled: !g.fetched || null, title: g.fetched ? g.title : `${g.title}: not fetched`,
        onclick: () => { if (s.genomes.has(g.id)) s.genomes.delete(g.id); else s.genomes.add(g.id); drawPicks(); } },
        el('span', { class: 'tick' }, icon('check', 'sm')), g.group === 'files' ? el('span', {}, speciesName(g)) : el('i', {}, speciesName(g))));
    }
    box.append(chips);
  }
  const n = chosenIds().length;
  view.count.textContent = `${n} chosen`;
  view.count.classList.toggle('short', n < 3);
  view.count.title = n < 3 ? 'a tree needs at least three' : '';
}

const chosenIds = () => [...view.s.genomes].filter((id) => view.ctx.catalogue.some((g) => g.id === id && g.fetched));

// The gene and the method from the address; true when they differ from what was built.
function fromAddress(p) {
  const s = view.s, before = JSON.stringify([s.gene, s.method, s.model, s.gamma, s.bootstrap]);
  if (GENES.some(([g]) => g === p.gene)) s.gene = p.gene;
  if (BY[p.method]) s.method = p.method;
  if (MODELS.some(([m]) => m === p.model)) s.model = p.model;
  if (p.gamma !== undefined) s.gamma = p.gamma !== '0';
  if (p.bootstrap !== undefined && Number.isFinite(Number(p.bootstrap))) s.bootstrap = Math.min(1000, Math.max(0, Math.round(Number(p.bootstrap))));
  const i = view.inputs;
  i.gene.value = s.gene; i.model.value = s.model; i.describe(); i.methods.set(s.method);
  i.gamma.querySelector('input').checked = s.gamma; i.gamma.hidden = s.method !== 'ml'; i.boot.input.value = s.bootstrap;
  return before !== JSON.stringify([s.gene, s.method, s.model, s.gamma, s.bootstrap]);
}

// ---- running ----------------------------------------------------------------------------------------
// A run in progress: its button disabled and turning, and a line saying what is being done with the
// seconds since it began. Returns the function that ends it.
function working(button, line, text) {
  const keep = [...button.childNodes];
  button.disabled = true;
  button.replaceChildren(el('span', { class: 'spin' }), el('span', {}, 'Working…'));
  const elapsed = el('span', { class: 'elapsed' });
  line.replaceChildren(el('span', { class: 'spin' }), el('span', {}, text), elapsed);
  const t0 = performance.now();
  const tick = () => { elapsed.textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`; };
  tick();
  const timer = setInterval(tick, 100);
  return () => { clearInterval(timer); button.replaceChildren(...keep); button.disabled = false; line.remove(); };
}

async function build() {
  if (view.running) return;
  const s = view.s, ids = chosenIds();
  if (ids.length < 3) { toast('Choose at least three genomes: a tree needs three leaves.', 'error'); return; }
  view.running = true;
  view.ctx.setParams({ gene: s.gene, method: s.method, model: s.model, gamma: s.gamma ? null : 0, bootstrap: s.bootstrap === 100 ? null : s.bootstrap });
  const text = `Aligning ${s.gene} from ${ids.length} genomes, then building the tree ${BY[s.method]}${s.bootstrap ? ` and ${s.bootstrap} bootstrap replicates` : ''}…`;
  const line = el('div', { class: 'wait' });
  clear(view.error);
  if (kept.tree) { view.head.append(line); view.main.classList.add('is-busy'); }
  else clear(view.main).append(card({ body: el('div', { style: { padding: '18px 18px 8px' } }, line, el('div', { class: 'skeleton', style: { height: '240px', margin: '6px 0 12px' } })) }));
  const stop = working(view.runBtn, line, text);
  const ticket = ++view.ticket;
  selection(ids);
  const t0 = performance.now();
  try {
    const r = await api.run('trees.gene', { genomes: ids, gene: s.gene, protein: false, method: s.method, model: s.model, gamma: s.gamma, bootstrap: s.bootstrap, seed: 1 });
    if (ticket !== view.ticket) return;
    kept.tree = { r, s: { ...s, genomes: ids }, at: new Date(), seconds: (performance.now() - t0) / 1000 };
    drawTree();
  } catch (e) {
    if (!kept.tree) clear(view.main);
    view.error.append(callout(e.message, 'bad'));
  } finally {
    stop();
    view.main.classList.remove('is-busy');
    view.running = false;
  }
}

// ---- the tree ----------------------------------------------------------------------------------------
function drawTree() {
  const { r, s, at, seconds } = kept.tree;
  const W = view.stage.clientWidth || 900;
  clear(view.head).append(el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Result'),
    el('span', { class: 'when' }, `${s.gene} of ${r.alignment.rows.length} genomes · ${at.toLocaleTimeString()} · ${fmt.seconds(seconds)}`));
  const box = clear(view.main);
  box.append(metrics([[fmt.count(r.alignment.rows.length), 'sequences'], [fmt.count(r.alignment.columns), 'columns aligned'], [fmt.count(r.alignment.variable), 'variable'],
    [fmt.count(r.alignment.informative), 'parsimony-informative'], r.log_likelihood !== undefined ? [fmt.number(r.log_likelihood, 7), 'log-likelihood'] : null,
    [fmt.number(r.tree.length, 4), 'tree length, a site']]));

  // The tree, and beside its leaves the distances between them, row by row.
  const D = r.distances, names = r.alignment.names;
  const leaves = r.tree.nodes.filter((d, v) => !r.tree.nodes.some((x) => x.parent === v)).length;
  const matrixW = D ? leaves * ROW + 2 : 0;
  const { svg, order } = treeSvg(r.tree.nodes, { width: Math.max(380, Math.min(720, W - 64 - matrixW)), rowHeight: ROW });
  const fig = el('div', { class: 'fig' });
  fig.innerHTML = svg;
  const leafIndex = order.map((v) => names.indexOf(r.tree.nodes[v].label));
  const parts = [fig];
  if (D && leafIndex.every((k) => k >= 0)) {
    const m = leafIndex.length, cell = ROW;
    const max = Math.max(...D.flat().filter(Number.isFinite));
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
    parts.push(el('div', { class: 'matrix' }, canvas, el('div', { class: 'cap' }, `${r.distance_model} distances; darker is closer`)));
  }
  const facts = [['aligned by', r.alignment.method]];
  if (r.parameters && Object.keys(r.parameters).length) facts.push(['the model fitted', Object.entries(r.parameters).map(([k, v]) => `${k} ${fmt.number(v, 3)}`).join(' · ')]);
  for (const [k, v] of Object.entries(r.compare || {})) facts.push([`Robinson–Foulds to ${k}`, String(v)]);
  if (r.bootstrap) facts.push(['bootstrap', `${r.bootstrap.replicates} replicates in ${fmt.seconds(r.bootstrap.seconds)}${r.bootstrap.threads ? ` on ${r.bootstrap.threads} threads` : ''}`]);
  if (r.seconds) facts.push(['computed', `the alignment in ${fmt.seconds(r.seconds.align)}, the tree in ${fmt.seconds(r.seconds.tree)}`]);
  const key = el('div', { class: 'support-key' }, el('span', {}, 'support of a split, the share of bootstrap trees holding it:'),
    ...[['--good', '95% or more'], ['--ink-2', '70 to 95%'], ['--bad', 'below 70%']].map(([c, t]) => el('span', { class: 'key-item' }, el('span', { class: 'swatch', style: { background: `var(${c})` } }), t)));
  const nwk = newick(r.tree.nodes);
  box.append(card({ title: `The tree of ${s.gene}`, sub: `${BY[s.method]}${s.method === 'ml' ? ` under ${s.model}${s.gamma ? '+Γ' : ''}` : ` on ${r.distance_model} distances`}${r.bootstrap ? `, ${r.bootstrap.replicates} bootstrap replicates` : ''}`,
    actions: [btn('Newick', { icon: 'download', size: 'sm', title: 'the tree in Newick', onClick: () => download(`${s.gene}-${s.method}.nwk`, nwk) })],
    body: [el('div', { class: 'tree-body' }, ...parts), el('div', { class: 'tree-facts' }, kv(facts), r.bootstrap ? key : null)] }));

  // The alignment, coloured by base or residue, its variable columns marked.
  if (r.alignment.rows.length && r.alignment.columns <= 6000) {
    const rows = r.alignment.rows, cols = r.alignment.columns;
    const cw = Math.max(1, Math.min(4, Math.floor((W - 260) / cols))), ch = 12;
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
    canvas.addEventListener('pointermove', (e) => {
      const rect = canvas.getBoundingClientRect();
      const k = Math.floor((e.clientY - rect.top) / ch), j = Math.floor((e.clientX - rect.left) / cw);
      if (k < 0 || j < 0 || k >= order.length || j >= cols) return;
      const name = r.tree.nodes[order[k]].label, row = rows[names.indexOf(name)];
      hover(`<i>${esc(name)}</i><br>column ${fmt.count(j + 1)}: <b>${esc(row?.[j] || '')}</b>${r.alignment.variable_columns && !r.alignment.variable_columns[j] ? ' · constant' : ''}`, e);
    });
    canvas.addEventListener('pointerleave', () => hover(null));
    const labels = el('div', { class: 'names' }, ...order.map((v) => el('div', { style: { height: `${ch}px`, lineHeight: `${ch}px` } }, r.tree.nodes[v].label)));
    const keyBases = el('span', { class: 'actions-key' }, ...['A', 'C', 'G', 'T'].map((b) => el('span', { class: 'key-item' }, el('span', { class: 'swatch', style: { background: `var(--${b})` } }), b)),
      el('span', { class: 'key-item' }, el('span', { class: 'swatch', style: { background: 'var(--paper-4)' } }), 'constant'));
    box.append(card({ title: 'The alignment', sub: 'rows in the order of the tree; variable columns coloured, constant ones pale', actions: [keyBases],
      body: el('div', { class: 'aln' }, labels, el('div', { class: 'scroll' }, canvas)) }));
  }
}

// A Newick string of the tree: lengths on the branches, supports (as percentages) on the inner nodes.
function newick(nodes) {
  const kids = nodes.map(() => []);
  let root = 0;
  nodes.forEach((d, v) => { if (d.parent < 0) root = v; else kids[d.parent].push(v); });
  const name = (s) => { const t = String(s || '').replace(/ /g, '_'); return /[(),:;'[\]\s]/.test(t) ? `'${t.replace(/'/g, "''")}'` : t; };
  const walk = (v) => {
    const d = nodes[v];
    let s = kids[v].length ? `(${kids[v].map(walk).join(',')})` : name(d.label);
    if (kids[v].length && d.parent >= 0 && d.support !== undefined && d.support !== null) s += Math.round(d.support * 100);
    if (d.parent >= 0 && Number.isFinite(d.length)) s += `:${+d.length.toPrecision(6)}`;
    return s;
  };
  return `${walk(root)};\n`;
}

// ---- selection: dN/dS between every pair ---------------------------------------------------------------
// Each pair's codons, aligned through their proteins; the synonymous and non-synonymous sites and
// differences, Jukes–Cantor corrected (Nei and Gojobori 1986). Neutral change gives ω = dN/dS = 1;
// purifying selection, which removes most changes of the protein, ω well below it.
async function selection(ids) {
  const s = view.s, ticket = ++view.dndsTicket;
  const line = el('div', { class: 'wait', style: { padding: '18px' } }, el('span', { class: 'spin' }), `Aligning the codons of ${s.gene} pair by pair…`);
  if (kept.dnds) view.selection.classList.add('is-busy');
  else clear(view.selection).append(card({ title: 'Selection on the gene', body: line }));
  try {
    const r = await api.run('trees.dnds', { genomes: ids, gene: s.gene });
    if (ticket !== view.dndsTicket) return;
    kept.dnds = r;
    drawSelection();
  } catch (e) {
    if (ticket !== view.dndsTicket) return;
    clear(view.selection).append(card({ title: 'Selection on the gene', body: el('div', { style: { padding: '14px 16px' } }, callout(e.message, 'bad')) }));
  } finally {
    if (ticket === view.dndsTicket) view.selection.classList.remove('is-busy');
  }
}

function drawSelection() {
  const r = kept.dnds, box = clear(view.selection);
  const pairs = r.pairs.filter((q) => q.dS > 0 && q.dN > 0 && Number.isFinite(q.dS) && Number.isFinite(q.dN));
  const short = (n) => n.replace(/^(\w)\w+ /, '$1. ');
  const sub = `${r.gene} under the ${r.code} code (table ${r.table}), its codons aligned through the proteins`;
  if (!pairs.length) { box.append(card({ title: 'Selection on the gene', sub, body: el('p', { class: 'note', style: { padding: '0 16px 14px' } }, 'No pair has both kinds of difference.') })); return; }
  const xs = pairs.map((q) => q.dS), ys = pairs.map((q) => q.dN);
  const x = [Math.min(...xs) / 1.6, Math.max(...xs) * 1.6], y = [Math.min(...ys, ...xs.map((v) => v * 0.01)) / 1.6, Math.max(...ys, ...xs) * 1.6];
  const line = (w, label, dash) => ({ points: [[x[0], x[0] * w], [x[1], x[1] * w]], colour: token('--ink-faint'), dash, width: 1, label });
  const omegas = pairs.map((q) => q.omega).filter(Number.isFinite).sort((a, b) => a - b);
  const median = omegas[Math.floor(omegas.length / 2)];
  const c = el('div', { class: 'fig' });
  c.innerHTML = chart({ width: 460, height: 340, x, y, log: 'xy', xlabel: 'dS, synonymous substitutions a synonymous site', ylabel: 'dN, a non-synonymous site',
    series: [line(1, 'ω = 1, neutral', '5 3'), line(0.1, 'ω = 0.1', '2 3'), line(0.01, 'ω = 0.01', '1 4'),
      { points: pairs.map((q) => [q.dS, q.dN]), colour: token('--accent'), line: false, dots: true }] });
  const table = el('table', { class: 'data wide' }, el('tr', {}, ...['pair', 'codons', 'dS', 'dN', 'ω'].map((h) => el('th', { class: h === 'pair' ? '' : 'num' }, h))),
    ...[...r.pairs].sort((a, b) => a.dS - b.dS).map((q) => el('tr', {}, el('td', { class: 'small' }, el('i', {}, `${short(r.names[q.a])} · ${short(r.names[q.b])}`)), el('td', { class: 'num' }, fmt.count(q.codons)),
      el('td', { class: 'num' }, q.dS === null ? 'saturated' : fmt.number(q.dS, 3)), el('td', { class: 'num' }, q.dN === null ? 'saturated' : fmt.number(q.dN, 3)),
      el('td', { class: 'num' }, q.omega === null ? '–' : fmt.number(q.omega, 2)))));
  box.append(card({ title: 'Selection on the gene', sub, actions: [el('span', { class: 'badge accent', title: 'the median of ω over the pairs' }, `median ω ${fmt.number(median, 2)}`)],
    body: el('div', { class: 'sel-body' },
      el('p', { class: 'note' }, `Of the possible changes at a codon, those that keep the amino acid are synonymous sites, the rest non-synonymous; the differences of each kind over its sites, corrected for multiple hits, are dS and dN (Nei and Gojobori 1986). The median ω here is ${fmt.number(median, 2)}: most changes of the protein were removed by selection before they could be seen, and the closest pairs show the highest ω, the trace of slightly harmful changes that selection has not yet had time to remove.`),
      el('div', { class: 'sel-row' }, c, el('div', { class: 'sel-table' }, table))) }));
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

// ---- sketches ----------------------------------------------------------------------------------------
// A genome's k-mers, hashed, and the smallest of them kept: two sketches estimate the Jaccard index
// of the k-mer sets, and the Mash distance turns it into a rate of mutation (Ondov et al. 2016).
const sketchIds = () => view.ctx.catalogue.filter((g) => g.fetched && g.group !== 'simulated').map((g) => g.id);

function sketchCard() {
  const sk = kept.sk;
  const k = numberField({ value: sk.k, min: 9, max: 31, step: 2, onChange: (v) => { sk.k = v; } });
  const size = numberField({ value: sk.size, min: 100, max: 100000, step: 500, unit: 'hashes', onChange: (v) => { sk.size = v; } });
  view.sketchBtn = btn('Sketch the genomes', { icon: 'play', kind: 'primary', onClick: sketch });
  const strip = el('div', { class: 'ctl-strip' }, fld('k, the length of a word', k), fld('sketch size', size), el('span', { class: 'grow' }),
    el('span', { class: 'toolbox-hint' }, icon('toolbox', 'sm'), el('span', {}, link('Your own genomes', () => view.ctx.openTool('sketch-compare', { genomes: sketchIds().join(', '), k: sk.k, size: sk.size })), ', in the Toolbox')),
    view.sketchBtn);
  strip.addEventListener('keydown', (e) => { if (e.key === 'Enter' && e.target.tagName === 'INPUT') { e.preventDefault(); e.target.dispatchEvent(new Event('change')); sketch(); } });
  view.sketchBody.append(empty({ icon: 'sketch', title: 'Not sketched yet', text: 'Every fetched genome is sketched, the simulated ones aside; the distances between the sketches, computed as a task graph of blocks, give a tree by neighbour joining.' }));
  return card({ cls: 'sketches', title: 'Mash distances and their tree', sub: 'every genome of the catalogue that is fetched: bacteria, viruses and mitochondria together', body: [strip, view.sketchBody] });
}

async function sketch() {
  if (view.sketching) return;
  view.sketching = true;
  const line = el('div', { class: 'wait', style: { padding: '14px 16px 0' } });
  if (kept.sketch) { view.sketchBody.prepend(line); [...view.sketchBody.children].slice(1).forEach((c) => c.classList.add('is-busy')); }
  else clear(view.sketchBody).append(line);
  const stop = working(view.sketchBtn, line, 'Sketching every fetched genome, then the distances between the sketches…');
  try {
    kept.sketch = await api.run('trees.sketches', { genomes: sketchIds(), k: kept.sk.k, size: kept.sk.size });
    drawSketch();
  } catch (e) {
    [...view.sketchBody.children].forEach((c) => c.classList.remove('is-busy'));
    view.sketchBody.prepend(el('div', { style: { padding: '14px 16px 0' } }, callout(e.message, 'bad')));
  } finally { stop(); view.sketching = false; }
}

function drawSketch() {
  const r = kept.sketch, W = view.section.clientWidth || 1200;
  const { svg } = treeSvg(r.tree.nodes, { width: Math.min(980, W - 80), labelWidth: 400, scaleLabel: 'Mash distance' });
  const fig = el('div', { class: 'fig' });
  fig.innerHTML = svg;
  clear(view.sketchBody).append(el('div', { class: 'sketch-body' },
    metrics([[fmt.count(r.genomes), 'genomes'], [`${r.k}`, 'k'], [fmt.count(r.size), 'hashes a sketch'], [fmt.count(r.pairs), 'distances'], [fmt.seconds(r.sketch_seconds), 'sketching'], [fmt.seconds(r.distance_seconds), `the distances, as ${r.tasks} task${r.tasks === 1 ? '' : 's'}`]]),
    fig,
    el('p', { class: 'note' }, `Each genome is the ${fmt.count(r.size)} smallest hashes of its canonical ${r.k}-mers. Beyond a distance of about 0.3 two sketches share almost nothing and the distances saturate at 1: sketches tell close relatives apart, and say little about how far apart the distant ones are.`)));
}

// ---- small pieces ----------------------------------------------------------------------------------------
const link = (text, onClick) => el('button', { type: 'button', class: 'link', onclick: onClick }, text);
const quiet = (...parts) => el('p', { class: 'toolbox-hint' }, icon('toolbox', 'sm'), el('span', {}, ...parts));
