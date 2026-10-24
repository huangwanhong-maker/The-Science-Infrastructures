// The Align view: two sequences aligned, and the alignment shown as what it is, a path through
// the edit graph. The score of every path is a product of weights in a semiring and the
// optimum a sum; every algorithm is a schedule of that sum, and the view draws each schedule's
// footprint: the full matrix of Gotoh with its traceback, the band, the wavefronts of WFA that
// rise by score and touch few cells, the posterior of a pair HMM, and the parametric polygon
// whose normal fan says which alignment wins for every pair of penalties.

import { el, clear, esc, notice, toggles, field, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import * as fmt from '../format.js';
import { BASE_COLOURS, token, warm, CATEGORIES } from '../palette.js';

const PRESETS = [
  { name: 'Durbin’s example (BLOSUM50, gap 8)', a: 'HEAGAWGHEE', b: 'PAWHEAE', alphabet: 'protein', matrix: 'BLOSUM50', open: 0, extend: 8, mode: 'global',
    note: 'Durbin, Eddy, Krogh and Mitchison (1998), Figs 2.5 and 2.6: the global alignment scores 1, HEAGAWGHE-E against --P-AW-HEAE, and the local 28, AWGHE against AW-HE.' },
  { name: 'Two reads of one place', a: 'ACGTTAGCATTCGGATCGAGGCTAACGTTAGCGATCGTAACGGATTAC', b: 'ACGTTAGCGATTCGGATCGAGCTAACGTTAGCGTCGTAACGGTATTAC', alphabet: 'dna', matrix: 'dna', open: 4, extend: 2, mode: 'global' },
  { name: 'A read placed on a reference', a: 'GATCGAGGCTAACGTTAGCG', b: 'TTGACACGTTAGCATTCGGATCGAGGCTAACGTTAGCGATCGTAACGGATTACCATG', alphabet: 'dna', matrix: 'dna', open: 4, extend: 2, mode: 'semiglobal' },
  { name: 'Cytochrome b of human and chimpanzee (proteins, from RefSeq)', genes: ['ncbi:NC_012920.1', 'ncbi:NC_001643.1'], gene: 'CYTB', protein: true, alphabet: 'protein', matrix: 'BLOSUM62', open: 11, extend: 1, mode: 'global',
    note: 'MT-CYB as RefSeq translates it in each mitochondrial genome, read by the engine from the fetched GenBank files.' },
  { name: 'Cytochrome b of human and gorilla (genes, from RefSeq)', genes: ['ncbi:NC_012920.1', 'ncbi:NC_011120.1'], gene: 'CYTB', protein: false, alphabet: 'dna', matrix: 'dna', open: 4, extend: 2, mode: 'global',
    note: 'The coding sequences of MT-CYB, read by the engine from the fetched GenBank files.' },
];

const ALGORITHMS = [
  ['gotoh', 'Gotoh', 'the full matrix, three states, and its traceback (Gotoh 1982)'],
  ['hirschberg', 'Hirschberg', 'linear space by divide and conquer (Hirschberg 1975; Myers and Miller 1988)'],
  ['banded', 'banded', 'only the diagonals of a band (Chao, Pearson and Miller 1992)'],
  ['wfa', 'WFA', 'wavefronts by penalty (Marco-Sola et al. 2021), global or semi-global'],
  ['biwfa', 'BiWFA', 'wavefronts from both ends to a breakpoint, in memory O(s) (Marco-Sola et al. 2023)'],
  ['myers', 'Myers', 'edit distance by bit vectors (Myers 1999)'],
  ['striped', 'striped', 'Farrar’s stripes on vector lanes (Farrar 2007), local only'],
];
// The modes each schedule is written for here; choosing a schedule moves the mode into them.
const MODES_OF = { wfa: ['global', 'semiglobal'], biwfa: ['global'], striped: ['local'] };

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, inputs: el('div', { class: 'align-inputs' }), result: el('div', { class: 'align-result' }), figure: el('div', { class: 'align-figure' }),
           tab: 'matrix', last: null, state: { ...PRESETS[0], algorithm: 'gotoh', band: 8 } };
  section.append(
    el('h2', { class: 'view-title' }, 'Align'),
    el('p', { class: 'lede' }, 'An alignment of two sequences is a path through their edit graph, from one corner to the other (or from anywhere to anywhere, for a local alignment); its score is the sum of its steps\' weights. Every algorithm here computes the best of all paths, which is a sum over paths in the tropical semiring, and each does it in its own order: the whole matrix, a band, the wavefronts of equal score, the columns as bit vectors, the lanes of a vector register. The same edit graph in other semirings counts the alignments, gives each pair its probability, or, in the polytope algebra, every optimal alignment for every scoring at once.'),
    el('div', { class: 'align-grid' }, view.inputs, el('div', {}, view.result, view.figure)));
  const p = ctx.params;
  if (p.a) Object.assign(view.state, { a: p.a, b: p.b || '', alphabet: p.alphabet || 'dna', matrix: p.matrix || 'dna', open: Number(p.open ?? 4), extend: Number(p.extend ?? 2), mode: p.mode || 'global', algorithm: p.algorithm || 'gotoh' });
  buildInputs();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.last) align();
}

function buildInputs() {
  const s = view.state, box = clear(view.inputs);
  const preset = el('select', { onchange: async () => {
    const q = PRESETS[Number(preset.value)];
    if (q.genes) {
      try {
        const [x, y] = await Promise.all(q.genes.map((g) => api.run('genome.gene', { genome: g, gene: q.gene })));
        Object.assign(s, q, { a: q.protein ? x.translation : x.sequence, b: q.protein ? y.translation : y.sequence });
      } catch (e) { notice(`${e.message} (fetch the NCBI genomes with Setup)`, true); return; }
    } else Object.assign(s, q);
    buildInputs();
    align();
  } },
    el('option', { value: '', disabled: true, selected: true }, 'an example…'), ...PRESETS.map((q, k) => el('option', { value: k }, q.name)));
  const a = el('textarea', { rows: 4, spellcheck: false, oninput: () => { s.a = a.value.replace(/\s+/g, '').toUpperCase(); } }, s.a);
  const b = el('textarea', { rows: 4, spellcheck: false, oninput: () => { s.b = b.value.replace(/\s+/g, '').toUpperCase(); } }, s.b);
  const alphabet = toggles([['dna', 'DNA'], ['protein', 'protein']], s.alphabet, (v) => { s.alphabet = v; s.matrix = v === 'dna' ? 'dna' : 'BLOSUM62'; buildInputs(); });
  const matrix = el('select', { onchange: () => { s.matrix = matrix.value; } });
  for (const m of s.alphabet === 'dna' ? ['dna'] : ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']) matrix.append(el('option', { value: m, selected: m === s.matrix }, m === 'dna' ? 'match +2, mismatch −4' : m));
  const open = el('input', { type: 'number', value: s.open, min: 0, max: 50, onchange: () => { s.open = Number(open.value); } });
  const extend = el('input', { type: 'number', value: s.extend, min: 1, max: 50, onchange: () => { s.extend = Number(extend.value); } });
  const mode = toggles([['global', 'global'], ['local', 'local'], ['semiglobal', 'semi-global', 'the first end to end, free ends on the second'], ['overlap', 'overlap', 'free end gaps on both']], s.mode, (v) => { s.mode = v; align(); });
  const algos = el('div', { class: 'stack' });
  for (const [k, label, title] of ALGORITHMS) {
    algos.append(el('label', { class: 'field', title }, el('input', { type: 'radio', name: 'algorithm', checked: s.algorithm === k, onchange: () => {
      s.algorithm = k;
      const modes = MODES_OF[k];
      if (modes && !modes.includes(s.mode)) { s.mode = modes[0]; buildInputs(); }
      align();
    } }), el('span', { style: { fontStyle: 'normal', color: 'var(--ink)' } }, label), el('span', { class: 'faint small' }, title)));
  }
  const random = el('button', { class: 'button small', onclick: async () => {
    try {
      const r = await api.run('align.random', { length: 300, divergence: 0.12, seed: Math.floor(Math.random() * 1e9) });
      Object.assign(s, { a: r.a, b: r.b, alphabet: 'dna', matrix: 'dna', open: 4, extend: 2 });
      buildInputs();
      align();
    } catch (e) { notice(e.message, true); }
  } }, 'two random relatives (300 bp, 12%)');
  box.append(
    el('div', { class: 'row' }, preset, random),
    el('label', { class: 'stack small' }, el('span', { class: 'faint' }, 'the first sequence (the query)'), a),
    el('label', { class: 'stack small' }, el('span', { class: 'faint' }, 'the second (the target)'), b),
    el('div', { class: 'row' }, alphabet, field('scores', matrix)),
    el('div', { class: 'row' }, field('gap open', open), field('extend', extend), el('span', { class: 'faint small' }, 'a gap of L costs open + L·extend')),
    el('div', { class: 'row' }, mode),
    el('p', { class: 'eyebrow', style: { marginTop: '10px' } }, 'the schedule'), algos,
    el('div', { class: 'row', style: { marginTop: '10px' } }, el('button', { class: 'button primary', onclick: align }, 'Align')),
    s.note ? el('p', { class: 'note' }, s.note) : null);
}

function request() {
  const s = view.state;
  return { a: s.a, b: s.b, alphabet: s.alphabet, matrix: s.matrix, gap_open: s.open, gap_extend: s.extend, mode: s.mode, algorithm: s.algorithm, band: s.band };
}

async function align() {
  const s = view.state;
  if (!s.a || !s.b) return;
  view.ctx.setParams({ a: s.a.length <= 400 ? s.a : null, b: s.b.length <= 400 ? s.b : null, alphabet: s.alphabet, matrix: s.matrix, open: s.open, extend: s.extend, mode: s.mode, algorithm: s.algorithm });
  clear(view.result).append(el('p', { class: 'wait' }, 'Aligning…'));
  try {
    const r = await api.run('align.pair', { ...request(), with_matrix: s.a.length * s.b.length <= 160000 });
    view.last = r;
    drawResult(r);
    drawFigure();
  } catch (e) { clear(view.result).append(el('p', { class: 'error' }, e.message)); }
}

function drawResult(r) {
  const box = clear(view.result);
  const kv = el('div', { class: 'row', style: { gap: '22px', alignItems: 'baseline' } },
    el('span', {}, el('span', { class: 'faint small' }, 'score '), el('b', { style: { fontSize: '22px', fontWeight: 500 } }, String(r.score))),
    r.penalty !== undefined ? el('span', { class: 'small' }, el('span', { class: 'faint' }, 'penalty '), String(r.penalty)) : null,
    r.distance !== undefined ? el('span', { class: 'small' }, el('span', { class: 'faint' }, 'edit distance '), String(r.distance)) : null,
    el('span', { class: 'small' }, el('span', { class: 'faint' }, 'identity '), fmt.percent(r.identity)),
    el('span', { class: 'small mono' }, r.cigar.length > 60 ? `${r.cigar.slice(0, 58)}…` : r.cigar),
    el('span', { class: 'small faint' }, `${fmt.count(r.cells)} ${r.penalty !== undefined ? 'offsets' : 'cells'} in ${fmt.seconds(r.seconds)}${r.gcups ? ` · ${fmt.number(r.gcups, 3)} GCUPS` : ''}`),
    r.memory !== undefined ? el('span', { class: 'small faint', title: 'the most bytes of wavefronts held at once' }, `${fmt.bytes(r.memory)} of wavefronts`) : null,
    r.lanes ? el('span', { class: 'small faint', title: `${r.lanes.passes} passes, ${fmt.count(r.lanes.lazy_f)} corrections of the lazy F loop` }, `${r.lanes.width} lanes of ${r.lanes.bits} bits (${r.lanes.isa})`) : null,
    r.unit_costs ? el('span', { class: 'small faint' }, `unit edit distance ${r.distance_unit}; the score shown rescores its path`) : null,
    r.reference ? el('span', { class: `badge ${r.reference.same ? 'good' : 'bad'}`, title: `Gotoh's full matrix gives ${r.reference.score}` }, r.reference.same ? 'the reference agrees' : `the reference gives ${r.reference.score}`) : null);
  box.append(kv);
  box.append(rows(r.rows, view.state.alphabet, r.query_begin, r.target_begin));
}

// The alignment's two rows and the line between, in blocks of 60 columns.
function rows(rw, alphabet, qb, tb) {
  const colours = BASE_COLOURS();
  const box = el('div', { class: 'alignment-rows' });
  const width = 60;
  let qi = qb, ti = tb;
  for (let k = 0; k < rw.query.length; k += width) {
    const q = rw.query.slice(k, k + width), m = rw.middle.slice(k, k + width), t = rw.target.slice(k, k + width);
    const line = (text, start, isMiddle) => {
      const span = el('div', { class: 'aline' });
      span.append(el('span', { class: 'apos' }, isMiddle ? '' : String(start + 1)));
      for (const c of text) span.append(el('span', { class: 'achar', style: isMiddle || alphabet !== 'dna' || c === '-' ? {} : { color: colours[c] || '' } }, c));
      return span;
    };
    box.append(el('div', { class: 'ablock' }, line(q, qi, false), line(m, 0, true), line(t, ti, false)));
    qi += [...q].filter((c) => c !== '-').length;
    ti += [...t].filter((c) => c !== '-').length;
  }
  return box;
}

// ---- the figure ----------------------------------------------------------------------------------------
const TABS = [['matrix', 'the edit graph'], ['algorithms', 'every schedule'], ['semirings', 'other semirings'], ['parametric', 'the parametric polygon'], ['poa', 'many sequences: a partial order']];

function drawFigure() {
  const box = clear(view.figure);
  box.append(el('div', { class: 'row', style: { margin: '18px 0 10px' } }, toggles(TABS.map(([k, l]) => [k, l]), view.tab, (v) => { view.tab = v; drawFigure(); })));
  const body = el('div');
  box.append(body);
  if (view.tab === 'matrix') drawMatrix(body);
  if (view.tab === 'algorithms') drawAlgorithms(body);
  if (view.tab === 'semirings') drawSemirings(body);
  if (view.tab === 'parametric') drawParametric(body);
  if (view.tab === 'poa') drawPoa(body);
}

// The DP matrix as a field, the path on it; WFA's wavefronts or the band over it.
function drawMatrix(body) {
  const r = view.last, M = r.matrix;
  if (!M) { body.append(el('p', { class: 'note' }, 'The matrix is drawn for sequences of up to 160,000 cells; these are longer.')); return; }
  const n = M.rows, m = M.cols;   // (n + 1) × (m + 1) values: query down, target across
  const size = Math.max(2, Math.min(18, Math.floor(Math.min(720 / (m + 1), 560 / (n + 1)))));
  const canvas = el('canvas', { width: (m + 1) * size + 1, height: (n + 1) * size + 1, class: 'dp' });
  const g = canvas.getContext('2d');
  const values = M.values;
  let lo = Infinity, hi = -Infinity;
  for (const v of values) if (Number.isFinite(v) && v > -1e8) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  // Under WFA's wavefronts, Gotoh's field is drawn faint: the cells the wavefronts never reach.
  if (M.fronts) g.globalAlpha = 0.32;
  for (let i = 0; i <= n; i++)
    for (let j = 0; j <= m; j++) {
      const v = values[i * (m + 1) + j];
      g.fillStyle = Number.isFinite(v) && v > -1e8 ? warm(0.08 + 0.85 * ((v - lo) / Math.max(1, hi - lo))) : token('--paper-3');
      g.fillRect(j * size, i * size, size, size);
    }
  g.globalAlpha = 1;
  // Each wavefront joins the furthest cells reached on its diagonals by paths of one penalty:
  // pale at penalty 0, deep at the optimum.
  if (M.fronts) {
    const F = M.fronts.length;
    // About forty fronts, evenly spaced in penalty, and always the last: more would fill the band.
    const stride = Math.max(1, Math.ceil(F / 40));
    M.fronts.forEach((w, k) => {
      if (k % stride && k !== F - 1) return;
      const t = F > 1 ? k / (F - 1) : 1;
      // A cool slate against the warm field, so that the path (ink) stays the darkest line.
      g.strokeStyle = `rgba(${Math.round(120 - 60 * t)}, ${Math.round(146 - 70 * t)}, ${Math.round(190 - 72 * t)}, ${0.3 + 0.5 * t})`;
      g.fillStyle = g.strokeStyle;
      g.lineWidth = Math.max(0.9, size / 7);
      // A diagonal no path of exactly this penalty reaches is null; with most penalties every
      // other one is, so the line bridges a single empty diagonal and breaks at two.
      g.beginPath();
      let gap = 2;
      for (const pt of w.points) {
        if (!pt) { gap++; continue; }
        const x = pt[1] * size + size / 2, y = pt[0] * size + size / 2;
        if (gap <= 1) g.lineTo(x, y); else g.moveTo(x, y);
        gap = 0;
      }
      g.stroke();
      for (const pt of w.points) if (pt) { g.beginPath(); g.arc(pt[1] * size + size / 2, pt[0] * size + size / 2, Math.max(1, size / 5), 0, 2 * Math.PI); g.fill(); }
    });
  }
  // The cells a schedule touched, when it says (the band, the wavefronts).
  if (M.touched) {
    g.fillStyle = 'rgba(58, 40, 48, 0.55)';
    for (const [i, j] of M.touched) g.fillRect(j * size + size * 0.35, i * size + size * 0.35, size * 0.3, size * 0.3);
  }
  if (M.band) {
    g.strokeStyle = token('--C');
    g.lineWidth = 1.5;
    for (const d of [M.band[0], M.band[1] + 1]) { g.beginPath(); g.moveTo(Math.max(0, d) * size, Math.max(0, -d) * size); g.lineTo((Math.min(m, n + d) + 0.0) * size, (Math.min(n, m - d)) * size); g.stroke(); }
  }
  // The path.
  g.strokeStyle = token('--ink');
  g.lineWidth = Math.max(M.fronts ? 2.4 : 1.5, size / 4);
  g.beginPath();
  M.path.forEach(([i, j], k) => { const x = j * size + size / 2, y = i * size + size / 2; if (k) g.lineTo(x, y); else g.moveTo(x, y); });
  g.stroke();
  // Letters along the edges when there is room.
  const wrap = el('div', { class: 'dp-wrap' });
  if (size >= 12) {
    const top = el('div', { class: 'dp-top', style: { paddingLeft: `${size}px` } }, ...[...view.state.b].map((c) => el('span', { style: { width: `${size}px` } }, c)));
    const left = el('div', { class: 'dp-left', style: { paddingTop: `${size}px` } }, ...[...view.state.a].map((c) => el('span', { style: { height: `${size}px`, lineHeight: `${size}px` } }, c)));
    wrap.append(el('div'), top, left, canvas);
  } else wrap.append(canvas);
  canvas.addEventListener('pointermove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const j = Math.floor((e.clientX - rect.left) / size), i = Math.floor((e.clientY - rect.top) / size);
    if (i < 0 || j < 0 || i > n || j > m) { hover(null); return; }
    const v = values[i * (m + 1) + j];
    hover(`cell (${i}, ${j}): ${esc(view.state.a[i - 1] || '·')} against ${esc(view.state.b[j - 1] || '·')}<br>best score of a path to here: <b>${Number.isFinite(v) && v > -1e8 ? v : '−∞'}</b>`, e);
  });
  canvas.addEventListener('pointerleave', () => hover(null));
  const caption = M.fronts
    ? `The wavefronts of the WFA: each line joins the furthest cells that paths of one penalty reach on each diagonal, from penalty 0 (pale) to the optimum, ${r.penalty} (deep); between two penalties the fronts slide along the diagonals while the letters match. ${M.fronts.length} wavefronts (${M.fronts.length > 40 ? `every ${Math.ceil(M.fronts.length / 40)}th drawn` : 'all drawn'}), ${fmt.count(r.record.offsets)} offsets computed and ${fmt.count(r.record.extended)} matches walked, against the ${fmt.count((n + 1) * (m + 1))} cells of the full matrix, which lies faint beneath as Gotoh's field of best scores. The penalties are the scoring's own, by Eizenga and Paten's exchange: mismatch ${r.penalties.mismatch}, gap ${r.penalties.gap_open} + ${r.penalties.gap_extend}·L.`
    : 'Each cell holds the best score of a path from the start to it; the dark line is the path traced back from the best end.';
  body.append(el('p', { class: 'note' }, caption), wrap);
}

// Every schedule on this pair: its value, its time, whether it agrees with the reference.
async function drawAlgorithms(body) {
  body.append(el('p', { class: 'wait' }, 'Running every schedule on this pair…'));
  let r;
  try { r = await api.run('align.compare', request()); } catch (e) { clear(body).append(el('p', { class: 'error' }, e.message)); return; }
  clear(body);
  // The rate at which each disposes of the whole edit graph: n·m over its time. A schedule that
  // computes fewer cells (wavefronts) or many at once (lanes, bit vectors) shows as faster.
  const rate = (x) => (x.seconds > 0 ? r.cells / x.seconds / 1e9 : 0);
  const gotoh = r.runs.find((x) => x.name === 'Gotoh');
  const table = el('table', { class: 'data' }, el('tr', {}, ...['schedule', 'computes', 'value', 'cells computed', 'time', 'n·m / time', 'agrees'].map((h) => el('th', { class: ['cells computed', 'time', 'n·m / time'].includes(h) ? 'num' : '' }, h))));
  for (const x of r.runs) {
    table.append(el('tr', {}, el('td', {}, x.name), el('td', { class: 'faint' }, x.computes), el('td', {}, String(x.value)), el('td', { class: 'num' }, fmt.count(x.cells)),
      el('td', { class: 'num' }, fmt.seconds(x.seconds)), el('td', { class: 'num', title: 'billions of cells of the full matrix a second' }, `${fmt.number(rate(x), 3)} G/s`),
      el('td', {}, x.agrees === null || x.agrees === undefined ? el('span', { class: 'faint small' }, 'another order') : el('span', { class: `badge ${x.agrees ? 'good' : 'bad'}` }, x.agrees ? 'yes' : 'no'))));
  }
  // Bars on a logarithmic axis, each against Gotoh's full matrix.
  const runs = r.runs.filter((x) => x.seconds > 0);
  const W = 640, row = 24, left = 190, right = 96, H = runs.length * row + 30;
  const lo = Math.log10(Math.min(...runs.map(rate)) / 1.6), hi = Math.log10(Math.max(...runs.map(rate)) * 1.3);
  const sx = (v) => left + ((Math.log10(v) - lo) / (hi - lo)) * (W - left - right);
  let svg = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" class="bars">`;
  for (let d = Math.ceil(lo); d <= Math.floor(hi); d++) {
    const x = sx(10 ** d);
    svg += `<line x1="${x}" y1="4" x2="${x}" y2="${H - 22}" stroke="${token('--hair')}" /><text class="tick-label" x="${x}" y="${H - 8}" text-anchor="middle">${10 ** d >= 1 ? 10 ** d : (10 ** d).toPrecision(1)}</text>`;
  }
  runs.forEach((x, k) => {
    const y = 6 + k * row, v = rate(x), speed = gotoh ? gotoh.seconds / x.seconds : null;
    svg += `<text class="tick-label" x="${left - 8}" y="${y + 13}" text-anchor="end" style="font-size:12px;fill:var(--ink-2)">${esc(x.name)}</text>`;
    svg += `<rect x="${left}" y="${y + 3}" width="${Math.max(1, sx(v) - left)}" height="${row - 9}" rx="2" fill="${CATEGORIES[k % CATEGORIES.length]}" fill-opacity="0.8"><title>${esc(x.name)}: ${fmt.number(v, 3)} billion cells of the full matrix a second</title></rect>`;
    if (speed) svg += `<text class="tick-label" x="${sx(v) + 6}" y="${y + 13}">${speed >= 1 ? `×${fmt.number(speed, 3)}` : `÷${fmt.number(1 / speed, 3)}`}</text>`;
  });
  svg += `<text class="tick-label" x="${(left + W - right) / 2}" y="${H + 0}" text-anchor="middle" style="font-style:italic"></text></svg>`;
  const c = el('div', { style: { marginTop: '14px' } }, el('p', { class: 'eyebrow' }, 'billions of cells of the full matrix a second, against Gotoh (×)'));
  c.insertAdjacentHTML('beforeend', svg);
  body.append(el('p', { class: 'note' }, 'Each schedule computes the optimum of the same sum, in its own order: scores maximised by Gotoh, Hirschberg and the stripes; penalties minimised by the wavefronts, whose global optimum is the same alignment (Eizenga and Paten 2022); unit costs by Myers\' bit vectors, a different scoring. Where the order is another (Myers; the ends-free WFA, which minimises the penalty rather than maximising the score), the value is not compared. Cells computed are the work done: the wavefronts\' offsets, Myers\' 64-row blocks, the stripes\' padded segments; n·m / time is the rate at which a schedule disposes of the whole edit graph.'), table, c);
}

async function drawSemirings(body) {
  body.append(el('p', { class: 'wait' }, 'Summing over the paths in other semirings…'));
  let r;
  try { r = await api.run('align.semiring', request()); } catch (e) { clear(body).append(el('p', { class: 'error' }, e.message)); return; }
  clear(body);
  const kv = el('dl', { class: 'kv' });
  const row = (k, v) => kv.append(el('dt', {}, k), el('dd', {}, v));
  if (r.alignments) row('all alignments (counting semiring)', el('span', { class: 'mono' }, r.alignments));
  if (r.optimal) row('optimal alignments (tropical with multiplicities)', el('span', { class: 'mono' }, r.optimal));
  if (r.delannoy) row('Delannoy number D(m, n)', el('span', { class: 'mono' }, r.delannoy));
  if (r.forward !== undefined) row('log P(x, y), the pair HMM’s forward sum', fmt.number(r.forward, 6));
  body.append(kv);
  if (r.posterior) {
    const { rows: n, cols: m, values } = r.posterior;
    const size = Math.max(2, Math.min(14, Math.floor(Math.min(640 / m, 420 / n))));
    const canvas = el('canvas', { width: m * size, height: n * size, class: 'dp' });
    const g = canvas.getContext('2d');
    for (let i = 0; i < n; i++)
      for (let j = 0; j < m; j++) { g.fillStyle = warm(Math.sqrt(values[i * m + j])); g.fillRect(j * size, i * size, size, size); }
    body.append(el('p', { class: 'note' }, 'The posterior probability that the i-th letter of the first sequence is aligned with the j-th of the second, under the pair HMM (Durbin et al. 1998, ch. 4): the probability semiring\'s forward and backward sums, divided. Where the alignment is certain the diagonal is dark; where it is not, the probability spreads.'), canvas);
    if (r.mea) body.append(el('p', { class: 'eyebrow', style: { marginTop: '12px' } }, `the alignment of maximum expected accuracy: ${fmt.number(r.mea.expected_pairs, 4)} pairs expected correct, accuracy ${fmt.percent(r.mea.accuracy)}`), rows(r.mea.rows, view.state.alphabet, 0, 0));
  }
}

async function drawParametric(body) {
  body.append(el('p', { class: 'wait' }, 'Propagating polygons through the edit graph…'));
  let r;
  try { r = await api.run('align.parametric', { a: view.state.a.slice(0, 120), b: view.state.b.slice(0, 120), alphabet: view.state.alphabet }); } catch (e) { clear(body).append(el('p', { class: 'error' }, e.message)); return; }
  clear(body);
  body.append(el('p', { class: 'note' }, `Every global alignment counts its mismatches x and its gap columns y; the scores match +1, mismatch −α, gap −β rank the alignments by αx + βy alone (the matches are fixed by x and y). The convex hull of all the (x, y) is the alignment polygon, computed in the polytope algebra, where ⊕ is the hull of a union and ⊗ a Minkowski sum (Pachter and Sturmfels 2004). Each vertex is the optimum for a cone of (α, β): the normal fan. Only the vertices whose cones meet the quadrant of positive penalties can win for a biologist's scoring; the others win where a mismatch or a gap is rewarded. ${r.vertices.length} vertices${view.state.a.length > 120 || view.state.b.length > 120 ? ', for the first 120 letters of each sequence' : ''}.`));
  const verts = r.vertices;
  const xs = verts.map((v) => v.mismatches), ys = verts.map((v) => v.gaps);
  const wrap = el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '28px' } });
  const poly = el('div');
  const W = 420, H = 360, pad = 40;
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  const sx = (x) => pad + ((x - x0) / Math.max(1, x1 - x0)) * (W - 2 * pad), sy = (y) => H - pad - ((y - y0) / Math.max(1, y1 - y0)) * (H - 2 * pad);
  let s = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`;
  s += `<polygon points="${verts.map((v) => `${sx(v.mismatches)},${sy(v.gaps)}`).join(' ')}" fill="${token('--accent')}" fill-opacity="0.1" stroke="${token('--accent')}" stroke-width="1.5" />`;
  verts.forEach((v, k) => { s += `<circle class="pvertex" data-k="${k}" cx="${sx(v.mismatches)}" cy="${sy(v.gaps)}" r="6" fill="${CATEGORIES[k % CATEGORIES.length]}" stroke="${token('--paper-2')}" stroke-width="1.5" style="cursor:pointer"><title>${v.mismatches} mismatches, ${v.gaps} gap columns</title></circle>`; });
  s += `<text class="tick-label" x="${W / 2}" y="${H - 8}" text-anchor="middle" style="font-style:italic">mismatches x</text><text class="tick-label" transform="translate(12 ${H / 2}) rotate(-90)" text-anchor="middle" style="font-style:italic">gap columns y</text></svg>`;
  poly.innerHTML = s;
  // The fan over the quadrant α, β ≥ 0: each point coloured by the vertex optimal there.
  const fan = el('div', { class: 'fan' });
  const G = r.grid, F = 300;
  const canvas = el('canvas', { width: G, height: G, style: { width: `${F}px`, height: `${F}px` } });
  const g = canvas.getContext('2d');
  const hex = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16));
  const img = g.createImageData(G, G);
  r.optimal.forEach((k, idx) => {
    const [cr, cg, cb] = hex(CATEGORIES[k % CATEGORIES.length]);
    img.data.set([Math.round(cr * 0.55 + 255 * 0.45), Math.round(cg * 0.55 + 255 * 0.45), Math.round(cb * 0.55 + 255 * 0.45), 255], idx * 4);
  });
  g.putImageData(img, 0, 0);
  canvas.style.imageRendering = 'auto';
  canvas.addEventListener('pointermove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const c = Math.floor(((e.clientX - rect.left) / rect.width) * G), row = Math.floor(((e.clientY - rect.top) / rect.height) * G);
    const k = r.optimal[row * G + c];
    if (k === undefined) return;
    const at = (k) => Math.exp(Math.log(r.low) + (Math.log(r.high) - Math.log(r.low)) * (k + 0.5) / G);
    const alpha = at(c), beta = at(G - 1 - row);
    hover(`α ${alpha.toFixed(2)}, β ${beta.toFixed(2)}<br>the optimum: vertex ${k + 1} (${verts[k].mismatches} mismatches, ${verts[k].gaps} gap columns)`, e);
  });
  canvas.addEventListener('pointerleave', () => hover(null));
  canvas.addEventListener('click', (e) => {
    const rect = canvas.getBoundingClientRect();
    const c = Math.floor(((e.clientX - rect.left) / rect.width) * G), row = Math.floor(((e.clientY - rect.top) / rect.height) * G);
    showVertex(r.optimal[row * G + c]);
  });
  fan.append(el('div', { class: 'fan-plot' }, el('span', { class: 'fan-y' }, `gap penalty β, ${fmt.number(r.low, 2)} to ${fmt.number(r.high, 3)}, logarithmic`), canvas, el('span', { class: 'fan-x' }, `mismatch penalty α, ${fmt.number(r.low, 2)} to ${fmt.number(r.high, 3)}, logarithmic`)));
  const detail = el('div', { class: 'small' }, el('p', { class: 'faint' }, 'Choose a vertex, or a point of the fan, to see its alignment.'));
  const showVertex = (k) => {
    const v = verts[k];
    if (!v) return;
    clear(detail).append(el('p', { class: 'eyebrow' }, el('span', { class: 'swatch', style: { background: CATEGORIES[k % CATEGORIES.length] } }), `vertex ${k + 1}: ${v.mismatches} mismatches, ${v.gaps} gap columns, ${v.matches} matches`), rows(v.rows, view.state.alphabet, 0, 0));
  };
  poly.addEventListener('click', (e) => { const c = e.target.closest('.pvertex'); if (c) showVertex(Number(c.dataset.k)); });
  wrap.append(el('div', {}, el('p', { class: 'eyebrow' }, 'the alignment polygon'), poly), el('div', {}, el('p', { class: 'eyebrow' }, 'its normal fan: who wins where'), fan));
  body.append(wrap, detail);
}

// ---- many sequences: partial-order alignment -------------------------------------------------------------
// Noisy reads of one sequence fused into a directed acyclic graph, read by read; the consensus is
// the heaviest bundle through it (Lee, Grasso and Sharlow 2002; Lee 2003). The graph is drawn by
// its columns (the nodes aligned to one another share one), each column's nodes stacked by weight.
const POA = { length: 150, reads: 10, error: 0.08, seed: 1 };

async function drawPoa(body) {
  const controls = el('div', { class: 'row' });
  const length = el('input', { type: 'number', value: POA.length, min: 20, max: 600, style: { width: '6em' }, onchange: () => { POA.length = Number(length.value); drawPoa(clear(body)); } });
  const reads = el('input', { type: 'number', value: POA.reads, min: 2, max: 40, style: { width: '5em' }, onchange: () => { POA.reads = Number(reads.value); drawPoa(clear(body)); } });
  const error = el('select', { onchange: () => { POA.error = Number(error.value); drawPoa(clear(body)); } },
    ...[0.02, 0.05, 0.08, 0.12, 0.15].map((e) => el('option', { value: e, selected: e === POA.error }, fmt.percent(e, 0))));
  controls.append(field('length', length), field('reads', reads), field('errors a base', error),
    el('button', { class: 'button small', onclick: () => { POA.seed++; drawPoa(clear(body)); } }, 'other reads'));
  body.append(controls, el('p', { class: 'wait' }, 'Fusing the reads into a graph…'));
  let r;
  try { r = await api.run('align.poa', POA); } catch (e) { body.lastChild.replaceWith(el('p', { class: 'error' }, e.message)); return; }
  body.lastChild.remove();
  const d = r.read_distances, mean = d.reduce((a, b) => a + b, 0) / d.length;
  body.append(el('p', { class: 'note' }, `${d.length} reads of a random sequence of ${POA.length} bp, each with ${fmt.percent(POA.error, 0)} errors a base (substitutions, insertions and deletions, as a long-read instrument makes them): between ${Math.min(...d)} and ${Math.max(...d)} edits from the truth, ${fmt.number(mean, 3)} on average. Each read is aligned to the graph of those before it, by dynamic programming over the graph in a topological order, and fused into it: where it agrees it adds weight to a node, where it differs it adds a branch, aligned to the nodes it differs from. The graph has ${fmt.count(r.nodes.length)} nodes and ${fmt.count(r.edges.length)} edges in ${r.columns} columns, made in ${fmt.seconds(r.seconds)}. Its heaviest path, the consensus, is `,
    el('b', {}, `${r.distance} ${r.distance === 1 ? 'edit' : 'edits'} from the truth`), '.'));

  // The graph.
  const byColumn = new Map();
  r.nodes.forEach(([, col, w], u) => { if (!byColumn.has(col)) byColumn.set(col, []); byColumn.get(col).push(u); });
  const onPath = new Set(r.path);
  const rank = new Array(r.nodes.length).fill(0);
  let depth = 0;
  for (const us of byColumn.values()) {
    us.sort((a, b) => (onPath.has(b) - onPath.has(a)) || r.nodes[b][2] - r.nodes[a][2]);
    us.forEach((u, k) => { rank[u] = k; });
    depth = Math.max(depth, us.length);
  }
  const dx = 24, dy = 30, pad = 18, maxW = Math.max(...r.nodes.map((x) => x[2]));
  const Hgraph = pad * 2 + (depth - 1) * dy + 20;
  const mid = pad + 10 + Math.floor((depth - 1) / 2) * dy;
  const X = (u) => pad + r.nodes[u][1] * dx, Y = (u) => mid + (rank[u] % 2 ? -1 : 1) * Math.ceil(rank[u] / 2) * dy;
  const W = pad * 2 + (r.columns - 1) * dx;
  const colours = BASE_COLOURS();
  const consensusEdge = new Set();
  for (let k = 1; k < r.path.length; k++) consensusEdge.add(`${r.path[k - 1]}:${r.path[k]}`);
  const maxE = Math.max(...r.edges.map((e) => e[2]));
  let g = `<svg width="${W}" height="${Hgraph}" viewBox="0 0 ${W} ${Hgraph}" class="poa">`;
  for (const [a, b, w] of r.edges) {
    const x1 = X(a), y1 = Y(a), x2 = X(b), y2 = Y(b), c = (x2 - x1) / 2;
    const main = consensusEdge.has(`${a}:${b}`);
    g += `<path d="M${x1},${y1} C${x1 + c},${y1} ${x2 - c},${y2} ${x2},${y2}" fill="none" stroke="${main ? token('--accent') : token('--ink-faint')}" stroke-opacity="${main ? 0.9 : 0.45}" stroke-width="${(0.6 + 3.4 * w / maxE).toFixed(2)}"><title>${w} ${w === 1 ? 'read' : 'reads'} pass this way</title></path>`;
  }
  r.nodes.forEach(([base, col, w], u) => {
    const main = onPath.has(u);
    g += `<g><title>${base} in column ${col + 1}: ${w} of ${d.length} reads${main ? '; on the consensus' : ''}</title><circle cx="${X(u)}" cy="${Y(u)}" r="8.5" fill="${colours[base] || colours.N}" fill-opacity="${(0.18 + 0.72 * w / maxW).toFixed(2)}" stroke="${main ? token('--accent') : 'none'}" stroke-width="1.6" />`;
    g += `<text x="${X(u)}" y="${Y(u) + 3.5}" text-anchor="middle" style="font:600 10px var(--mono);fill:var(--ink)">${base}</text></g>`;
  });
  g += '</svg>';
  const scroll = el('div', { class: 'poa-scroll' });
  scroll.innerHTML = g;

  // The multiple alignment the graph induces: one row a read, the consensus last.
  const msa = el('div', { class: 'alignment-rows poa-msa' });
  const cols = r.msa[0].length;
  const disagree = new Array(cols).fill(false);
  for (let j = 0; j < cols; j++) {
    const seen = new Set(r.msa.slice(0, -1).map((row) => row[j]));
    disagree[j] = seen.size > 1;
  }
  r.msa.forEach((row, k) => {
    const last = k === r.msa.length - 1;
    const line = el('div', { class: 'aline' }, el('span', { class: 'apos', style: { width: '8.5em' } }, last ? 'consensus' : `read ${k + 1} · ${d[k]}`));
    for (let j = 0; j < row.length; j++) {
      const c = row[j];
      line.append(el('span', { class: 'achar', style: { color: c === '-' ? 'var(--ink-faint)' : colours[c] || '', fontWeight: last ? 700 : 400, background: disagree[j] ? 'var(--paper-3)' : '' } }, c));
    }
    msa.append(line);
  });
  body.append(el('p', { class: 'eyebrow', style: { marginTop: '12px' } }, 'the graph: nodes by column, each read’s weight in the colour, the consensus ringed'), scroll,
    el('p', { class: 'eyebrow', style: { marginTop: '14px' } }, 'the multiple alignment it induces: each read with its edits from the truth, the consensus last; shaded, the columns where the reads disagree'), msa);
}
