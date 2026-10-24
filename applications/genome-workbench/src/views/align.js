// The Align view: two sequences aligned, and the alignment shown as what it is, a path through
// the edit graph. The score of every path is a product of weights in a semiring and the
// optimum a sum; every algorithm is a schedule of that sum, and the view draws each schedule's
// footprint: the full matrix of Gotoh with its traceback, the band, the wavefronts of WFA that
// rise by score and touch few cells, the posterior of a pair HMM, and the parametric polygon
// whose normal fan says which alignment wins for every pair of penalties.
//
// A bench: the instruments at the left (the two sequences, the scoring and the mode, the schedule
// chosen among cards), the alignment and its figures on the stage. The examples here carry what
// is known of them; the Toolbox's pairwise and multiple alignments take the user's own files.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import * as fmt from '../format.js';
import * as library from '../core/library.js';
import { BASE_COLOURS, token, warm, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, iconBtn, seg, select, numberField, fld, seqEditor, parseSequences, guessAlphabet, card, tabs, metrics, callout, empty, kv, download, copy, toast } from '../ui/kit.js';

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
const RANDOM_NOTE = 'Two descendants of one random sequence of 300 bp, about 12% apart, simulated by the engine from a new seed each time.';

// The schedules of the one sum: what each computes, its cost (n and m the lengths), its source, and
// the modes and letters it is written for here. Choosing a schedule moves the mode into its own.
const SCHEDULES = [
  { key: 'gotoh', name: 'Gotoh', short: 'the full matrix', cost: 'Time and space O(nm).', about: 'the full matrix, three states, and its traceback (Gotoh 1982)' },
  { key: 'hirschberg', name: 'Hirschberg', short: 'linear space', cost: 'Time O(nm), space O(n + m).', about: 'linear space by divide and conquer (Hirschberg 1975; Myers and Miller 1988)' },
  { key: 'banded', name: 'banded', short: 'a band of diagonals', cost: 'Time and space O(kn), for k diagonals.', about: 'only the diagonals of a band (Chao, Pearson and Miller 1992)' },
  { key: 'wfa', name: 'WFA', short: 'wavefronts by penalty', cost: 'Time O(ns), space O(s²), for a penalty s.', about: 'wavefronts by penalty (Marco-Sola et al. 2021), global or semi-global', modes: ['global', 'semiglobal'], tag: 'global, semi-global' },
  { key: 'biwfa', name: 'BiWFA', short: 'from both ends', cost: 'Time O(ns), space O(s).', about: 'wavefronts from both ends to a breakpoint, in memory O(s) (Marco-Sola et al. 2023)', modes: ['global'], tag: 'global' },
  { key: 'myers', name: 'Myers', short: 'bit vectors', cost: 'Time O(⌈m/w⌉ n), for words of w bits.', about: 'edit distance by bit vectors (Myers 1999)', dna: true, tag: 'DNA' },
  { key: 'striped', name: 'striped', short: 'vector lanes', cost: 'Time about nm / p, for p lanes.', about: 'Farrar’s stripes on vector lanes (Farrar 2007), local only', modes: ['local'], tag: 'local' },
];
const schedule = (k) => SCHEDULES.find((x) => x.key === k) || SCHEDULES[0];

const MODES = [
  ['global', 'global', 'end to end, both sequences'],
  ['local', 'local', 'the best-scoring pair of stretches, anywhere in each'],
  ['semiglobal', 'semi-global', 'the first end to end, free ends on the second'],
  ['overlap', 'overlap', 'free end gaps on both'],
];
const modeName = (k) => MODES.find(([v]) => v === k)?.[1] || k;
const MATRICES = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250'];

const TABS = [['matrix', 'The edit graph'], ['algorithms', 'Every schedule'], ['semirings', 'Other semirings'], ['parametric', 'The parametric polygon'], ['poa', 'Many sequences: a partial order']];

let view = null;
// Every request is numbered across mounts (the view is mounted again on a change of theme), so that
// only the latest one's answer is drawn.
let serial = 0;

export function mount(section, ctx) {
  view = { ctx, section, tab: 'matrix', last: null, shown: null, ticket: 0, figTicket: 0, figShown: null, timer: null,
           state: { ...PRESETS[0], example: 0, algorithm: 'gotoh', band: 8 } };
  const p = ctx.params;
  if (p.a) {
    const alphabet = p.alphabet === 'protein' ? 'protein' : 'dna';
    Object.assign(view.state, { a: p.a, b: p.b || '', alphabet, matrix: alphabet === 'dna' ? 'dna' : MATRICES.includes(p.matrix) ? p.matrix : 'BLOSUM62', open: Number(p.open ?? 4), extend: Number(p.extend ?? 2),
      mode: MODES.some(([k]) => k === p.mode) ? p.mode : 'global', algorithm: SCHEDULES.some((x) => x.key === p.algorithm) ? p.algorithm : 'gotoh', band: Number(p.band ?? 8) });
    const k = PRESETS.findIndex((q) => q.a === view.state.a && q.b === view.state.b);
    view.state.example = k >= 0 ? k : null;
    view.state.note = k >= 0 ? PRESETS[k].note || null : null;
  }
  if (TABS.some(([k]) => k === p.tab)) view.tab = p.tab;

  view.controls = el('div', { class: 'controls' });
  view.result = el('div');
  view.figure = figureCard();
  const bench = el('div', { class: 'bench wide-controls' }, view.controls, el('div', { class: 'stage' }, view.result, view.figure));
  section.append(
    el('h2', { class: 'view-title' }, 'Align'),
    el('p', { class: 'lede' }, 'An alignment of two sequences is a path through their edit graph, from one corner to the other (or from anywhere to anywhere, for a local alignment); its score is the sum of its steps\' weights. Every algorithm here computes the best of all paths, which is a sum over paths in the tropical semiring, and each does it in its own order: the whole matrix, a band, the wavefronts of equal score, the columns as bit vectors, the lanes of a vector register. The same edit graph in other semirings counts the alignments, gives each pair its probability, or, in the polytope algebra, every optimal alignment for every scoring at once.'),
    bench);
  // Ctrl Enter anywhere in the view aligns; Enter in a number of the controls does too. (On the
  // bench, which each mount makes anew, rather than on the section, which stays.)
  bench.addEventListener('keydown', (e) => {
    if (e.key !== 'Enter') return;
    const number = e.target.matches?.('.controls input[type="number"]');
    if (!(e.ctrlKey || e.metaKey) && !number) return;
    e.preventDefault();
    if (number) e.target.dispatchEvent(new Event('change'));
    align();
  });
  buildControls();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.last) align();
}

// ---- the instruments ---------------------------------------------------------------------------------
function buildControls() {
  const s = view.state;
  const fromLibrary = (anchor, use) => library.pick(anchor, ['fasta'], async (r) => use(await library.text(r.id)));
  view.edA = seqEditor({ value: s.a, rows: 3, library: fromLibrary, placeholder: 'the query: paste it, FASTA, or drop a file', onChange: (t) => takeText('a', t) });
  view.edB = seqEditor({ value: s.b, rows: 3, library: fromLibrary, placeholder: 'the target: paste it, FASTA, or drop a file', onChange: (t) => takeText('b', t) });
  view.exampleBox = el('div');
  view.sizeNote = el('div', { class: 'ctl-note' });
  view.scoringBox = el('div', { class: 'ctl-group' });
  view.schedBox = el('div', { class: 'ctl-group' });
  view.runBtn = btn('Align', { icon: 'play', kind: 'primary', onClick: align, title: 'align the two (Ctrl Enter)' });
  clear(view.controls).append(
    el('div', { class: 'card' },
      el('div', { class: 'ctl-group' }, el('p', { class: 'eyebrow' }, 'Two sequences'), view.exampleBox,
        fld('The first, the query', view.edA), fld('The second, the target', view.edB), view.sizeNote),
      view.schedBox,
      view.scoringBox,
      el('div', { class: 'run-bar' }, view.runBtn)),
    el('div', { class: 'toolbox-hint' }, icon('toolbox', 'sm'), el('span', {},
      'Ctrl Enter aligns. The examples here carry what is known of them; with your own files, ',
      el('button', { class: 'link', onclick: toToolbox, title: 'the Toolbox\'s pairwise alignment, with these two sequences and this scoring' }, 'take these two to the Toolbox'), '.')));
  drawExample();
  drawScoring();
  drawSchedules();
  sizeNote();
}

function drawExample() {
  const s = view.state;
  const options = [{ value: '', label: 'your own sequences', disabled: true },
    ...PRESETS.map((q, k) => ({ value: String(k), label: q.name, group: 'Examples' })),
    { value: 'random', label: 'Two random relatives, 300 bp, 12% apart', group: 'Simulated' }];
  const pick = select(options, s.example ?? '', (v) => takeExample(v, pick));
  const again = s.example === 'random' ? iconBtn('refresh', { title: 'two other relatives', onClick: () => takeExample('random', pick) }) : null;
  clear(view.exampleBox).append(fld('Example', el('div', { class: 'al-pick' }, pick, again)));
}

// The size of the edit graph, and a word when the full matrix would be very large.
function sizeNote() {
  const s = view.state, n = s.a.length, m = s.b.length;
  const big = n * m > 2e8 && s.algorithm === 'gotoh';
  view.sizeNote.textContent = n && m ? `${fmt.count(n)} × ${fmt.count(m)} letters: ${fmt.count(n * m)} cells in the full matrix.${big ? ' Gotoh holds every one; for a pair this large, Hirschberg, a band or the wavefronts.' : ''}` : '';
  view.sizeNote.classList.toggle('warn', big);
}

function drawScoring() {
  const s = view.state;
  const letters = seg([['dna', 'DNA'], ['protein', 'protein']], s.alphabet, (v) => setAlphabet(v), { block: true });
  const matrix = select(s.alphabet === 'dna' ? [['dna', '+2 / −4']] : MATRICES.map((m) => [m, m]), s.matrix, (v) => { s.matrix = v; stale(); }, { title: s.alphabet === 'dna' ? 'match +2, mismatch −4' : 'the substitution matrix' });
  const open = numberField({ value: s.open, min: 0, max: 50, onChange: (v) => { s.open = v; stale(); } });
  const extend = numberField({ value: s.extend, min: 1, max: 50, onChange: (v) => { s.extend = v; stale(); } });
  clear(view.scoringBox).append(el('p', { class: 'eyebrow' }, 'Scoring'),
    el('div', { class: 'ctl-row' }, fld('Letters', letters), fld('Scores', matrix)),
    el('div', { class: 'ctl-row' }, fld('Gap open', open), fld('Gap extend', extend)),
    el('div', { class: 'ctl-note' }, `${s.alphabet === 'dna' ? 'A match scores +2 and a mismatch −4; a' : 'A'} gap of L letters costs open + L · extend.`));
}

// The mode, and the schedules as cards: the name, what it does in a few words, its cost; the
// chosen one's description and source beneath, and its band when it has one. A card written for
// other modes is quieter, and choosing it moves the mode.
function drawSchedules() {
  const s = view.state;
  const mode = seg(MODES, s.mode, (v) => setMode(v), { block: true });
  const grid = el('div', { class: 'sched-grid', role: 'radiogroup', 'aria-label': 'the schedule' });
  for (const x of SCHEDULES) {
    const blocked = x.dna && s.alphabet !== 'dna';
    const elsewhere = x.modes && !x.modes.includes(s.mode);
    const on = s.algorithm === x.key;
    grid.append(el('button', { type: 'button', role: 'radio', 'aria-checked': String(on), class: `sched${on ? ' on' : ''}${elsewhere ? ' elsewhere' : ''}`, disabled: blocked || null,
      title: blocked ? 'Myers\' bit vectors here count unit edits on DNA' : elsewhere ? `written here for ${x.modes.map(modeName).join(' and ')} alignment: choosing it moves the mode` : x.about,
      onclick: () => pickSchedule(x.key) },
      el('span', { class: 'n' }, x.name, x.tag ? el('span', { class: 'tag' }, x.tag) : null),
      el('span', { class: 's' }, x.short)));
  }
  const cur = schedule(s.algorithm);
  const band = cur.key === 'banded'
    ? [fld('Band', numberField({ value: s.band, min: 0, max: 100000, unit: 'diagonals', onChange: (v) => { s.band = v; stale(); } })),
      el('div', { class: 'ctl-note' }, 'On each side of the diagonals between the corners: beyond the difference of the two lengths.')]
    : [];
  clear(view.schedBox).append(el('p', { class: 'eyebrow' }, 'Mode and schedule'),
    fld('Mode', mode), el('div', { class: 'ctl-note' }, `${modeName(s.mode)}: ${MODES.find(([v]) => v === s.mode)?.[2] || ''}.`),
    fld(el('span', { class: 'al-lab' }, 'Schedule', el('span', { class: 'aside' }, 'the order of one sum')), grid),
    el('div', { class: 'sched-about' }, el('b', {}, cur.name), `: ${cur.about}. `, el('span', { class: 'cost' }, cur.cost)), ...band);
}

function pickSchedule(k) {
  const s = view.state;
  s.algorithm = k;
  const modes = schedule(k).modes;
  if (modes && !modes.includes(s.mode)) s.mode = modes[0];
  drawSchedules();
  sizeNote();
  align();
}

function setMode(v) {
  const s = view.state;
  s.mode = v;
  const cur = schedule(s.algorithm);
  // A mode the schedule is not written for here moves the schedule to Gotoh's, which takes all four.
  if (cur.modes && !cur.modes.includes(v)) {
    s.algorithm = 'gotoh';
    toast(`${cur.name} is written here for ${cur.modes.map(modeName).join(' and ')} alignment; Gotoh's full matrix takes the ${modeName(v)} one.`);
  }
  drawSchedules();
  align();
}

function setAlphabet(v) {
  const s = view.state;
  if (s.alphabet === v) return;
  s.alphabet = v;
  s.matrix = v === 'dna' ? 'dna' : 'BLOSUM62';
  if (v !== 'dna' && s.algorithm === 'myers') s.algorithm = 'gotoh';
  drawScoring();
  drawSchedules();
  stale();
}

// The letters from an editor: a FASTA of two or more records, given to either, fills both; the
// letters' kind follows the sequences when they plainly say so.
const letters = (seq) => seq.replace(/-/g, '');
function takeText(which, text) {
  const s = view.state;
  const records = parseSequences(text);
  if (records.length >= 2) {
    s.a = letters(records[0].seq);
    s.b = letters(records[1].seq);
    view.edA.set(`>${records[0].name}\n${s.a}`);
    view.edB.set(`>${records[1].name}\n${s.b}`);
    toast(`${records[0].name} is the query and ${records[1].name} the target${records.length > 2 ? `; ${records.length - 2} more ${records.length === 3 ? 'record was' : 'records were'} left aside` : ''}.`);
  } else s[which] = letters(records[0]?.seq || '');
  if (s.example !== null) { s.example = null; s.note = null; drawExample(); }
  if (s.a.length + s.b.length >= 24 && s.a && s.b) {
    const kind = guessAlphabet(s.a + s.b);
    if (kind !== s.alphabet) { setAlphabet(kind); toast(`Read as ${kind === 'dna' ? 'DNA' : 'protein'}: the scores are ${s.alphabet === 'dna' ? 'match +2, mismatch −4' : s.matrix} now.`); }
  }
  sizeNote();
  stale();
}

async function takeExample(v, pick) {
  const s = view.state;
  pick.disabled = true;
  try {
    if (v === 'random') {
      const r = await api.run('align.random', { length: 300, divergence: 0.12, seed: Math.floor(Math.random() * 1e9) });
      Object.assign(s, { a: r.a, b: r.b, alphabet: 'dna', matrix: 'dna', open: 4, extend: 2, example: 'random', note: RANDOM_NOTE });
    } else {
      const q = PRESETS[Number(v)];
      if (q.genes) {
        const [x, y] = await Promise.all(q.genes.map((g) => api.run('genome.gene', { genome: g, gene: q.gene })));
        Object.assign(s, q, { a: q.protein ? x.translation : x.sequence, b: q.protein ? y.translation : y.sequence });
      } else Object.assign(s, q);
      s.example = Number(v);
      s.note = q.note || null;
    }
  } catch (e) {
    toast(`${e.message}${v === 'random' ? '' : ' (fetch the NCBI genomes with Setup)'}`, 'error');
    drawExample();
    return;
  }
  // An example keeps its own mode; a schedule not written for it, or for its letters, gives way to Gotoh's.
  const cur = schedule(s.algorithm);
  if ((cur.modes && !cur.modes.includes(s.mode)) || (cur.dna && s.alphabet !== 'dna')) s.algorithm = 'gotoh';
  view.edA.set(s.a);
  view.edB.set(s.b);
  drawExample();
  drawScoring();
  drawSchedules();
  sizeNote();
  align();
}

function toToolbox() {
  const s = view.state;
  view.ctx.openTool('align-pair', { a: s.a, b: s.b, mode: s.mode, algorithm: s.algorithm, alphabet: s.alphabet, gap_open: s.open, gap_extend: s.extend,
    ...(s.alphabet === 'dna' ? { match: 2, mismatch: -4 } : { matrix: s.matrix }) });
}

// ---- aligning ----------------------------------------------------------------------------------------
function request(st = view.state) {
  return { a: st.a, b: st.b, alphabet: st.alphabet, matrix: st.matrix, gap_open: st.open, gap_extend: st.extend, mode: st.mode, algorithm: st.algorithm, band: st.band };
}

// The inputs have changed since the alignment shown: it says so until the next.
function stale() {
  if (view.staleBadge) view.staleBadge.hidden = false;
}

// The run button while the engine works: disabled, turning, and the time so far once it is felt.
function busy(on) {
  const b = view.runBtn;
  clearInterval(view.timer);
  if (!on) { b.disabled = false; b.replaceChildren(icon('play', 'sm'), el('span', {}, 'Align')); return; }
  b.disabled = true;
  const label = el('span', {}, 'Aligning…');
  b.replaceChildren(el('span', { class: 'spin' }), label);
  const t0 = performance.now();
  const timer = view.timer = setInterval(() => {
    if (!b.isConnected) { clearInterval(timer); return; }
    const t = (performance.now() - t0) / 1000;
    if (t > 0.4) label.textContent = `Aligning… ${t.toFixed(1)} s`;
  }, 100);
}

async function align() {
  const s = view.state;
  if (!s.a || !s.b) { if (view.last) toast('Both sequences are needed to align.', 'error'); else drawEmpty(); return; }
  const ticket = view.ticket = ++serial;
  view.ctx.setParams({ a: s.a.length <= 400 ? s.a : null, b: s.b.length <= 400 ? s.b : null, alphabet: s.alphabet, matrix: s.matrix, open: s.open, extend: s.extend, mode: s.mode, algorithm: s.algorithm, band: s.algorithm === 'banded' ? s.band : null });
  busy(true);
  // The earlier alignment stays in view, dimmed, until the new one arrives.
  if (view.resultCard?.isConnected) { view.resultCard.classList.add('loading'); view.resultCard.body.classList.add('is-busy'); } else drawSkeleton();
  const asked = { ...s };
  try {
    const r = await api.run('align.pair', { ...request(asked), with_matrix: asked.a.length * asked.b.length <= 160000 });
    if (ticket !== view.ticket) return;
    view.last = r;
    view.shown = asked;
    drawResult(r);
    view.figure.hidden = false;
    drawFigure(true);
  } catch (e) {
    if (ticket !== view.ticket) return;
    drawError(e);
  } finally {
    if (ticket === view.ticket) busy(false);
  }
}

function drawEmpty() {
  view.figure.hidden = true;
  view.resultCard = null;
  clear(view.result).append(card({ body: empty({ icon: 'align', title: 'Two sequences make an alignment', text: 'Paste them, read them from files or from your data, or choose an example; then align them.' }) }));
}

function drawSkeleton() {
  const tiles = el('div', { class: 'al-skel' }, ...[0, 1, 2, 3, 4].map(() => el('div', { class: 'skeleton', style: { height: '62px' } })));
  clear(view.result).append(card({ title: 'The alignment', sub: 'aligning…', pad: true, body: [tiles, el('div', { class: 'skeleton', style: { height: '92px', marginTop: '14px' } })] }));
}

function drawError(e) {
  if (view.resultCard?.isConnected) {
    view.resultCard.classList.remove('loading');
    view.resultCard.body.classList.remove('is-busy');
    stale();
    view.result.querySelector('.al-error')?.remove();
    view.result.prepend(el('div', { class: 'al-error' }, callout(e.message, 'bad')));
  } else clear(view.result).append(callout(e.message, 'bad'));
}

function drawResult(r) {
  const sh = view.shown;
  const ref = r.reference;
  const refBadge = !ref ? null
    : ref.same === null || ref.same === undefined
      ? el('span', { class: 'badge', title: 'this schedule orders the alignments otherwise: its value is not compared' }, `another order; Gotoh gives ${ref.score}`)
      : el('span', { class: `badge ${ref.same ? 'good' : 'bad'}`, title: `Gotoh's full matrix gives ${ref.score}` }, ref.same ? 'the reference agrees' : `the reference gives ${ref.score}`);
  view.staleBadge = el('span', { class: 'badge warn', hidden: true, title: 'the inputs have changed since this alignment; Align to update it' }, 'of earlier inputs');
  const fasta = () => `>query\n${r.rows.query}\n>target\n${r.rows.target}\n`;
  const tiles = metrics([
    [String(r.score), r.unit_costs ? 'score of its path, rescored' : 'score'],
    r.penalty !== undefined ? [String(r.penalty), 'penalty, minimised'] : null,
    r.unit_costs ? [String(r.distance_unit), 'unit edit distance'] : null,
    r.distance !== undefined ? [String(r.distance), 'edit distance'] : null,
    [fmt.percent(r.identity), 'identity'],
    [fmt.count(r.cells), `${r.penalty !== undefined ? 'offsets' : 'cells'} computed`],
    [r.seconds > 0 && r.seconds < 1e-6 ? '< 1 µs' : fmt.seconds(r.seconds), r.gcups ? `${fmt.number(r.gcups, 3)} GCUPS` : 'time'],
    r.memory !== undefined ? [fmt.bytes(r.memory), 'of wavefronts held at most'] : null,
    r.lanes ? [`${r.lanes.width} × ${r.lanes.bits} bits`, `lanes (${r.lanes.isa}), ${r.lanes.passes} passes, ${fmt.count(r.lanes.lazy_f)} lazy-F corrections`] : null,
  ]);
  const cigar = el('div', { class: 'al-cigar' }, el('span', { class: 'eyebrow' }, 'CIGAR'),
    el('code', { title: r.cigar }, r.cigar.length > 72 ? `${r.cigar.slice(0, 70)}…` : r.cigar),
    iconBtn('copy', { title: 'copy the CIGAR', size: 'sm', onClick: () => copy(r.cigar) }),
    el('span', { class: 'small faint' }, `query ${fmt.count(r.query_begin + 1)}..${fmt.count(r.query_end)}, target ${fmt.count(r.target_begin + 1)}..${fmt.count(r.target_end)}`));
  const c = card({ title: 'The alignment', cls: 'al-card',
    sub: `${schedule(r.algorithm).name} · ${modeName(r.mode)} · ${sh.alphabet === 'dna' ? 'match +2, mismatch −4' : sh.matrix}, gap ${sh.open} + ${sh.extend}·L`,
    actions: [view.staleBadge, refBadge, iconBtn('download', { title: 'save the alignment as FASTA', size: 'sm', onClick: () => download('alignment.fasta', fasta()) })],
    pad: true,
    body: [tiles, cigar,
      r.unit_costs ? el('p', { class: 'note' }, 'Myers\' bit vectors count unit edits; the score shown rescores the path they found under the scoring chosen.') : null,
      rows(r.rows, sh.alphabet, r.query_begin, r.target_begin),
      sh.note ? el('p', { class: 'note al-source' }, icon('info', 'sm'), el('span', {}, sh.note)) : null] });
  c.prepend(el('div', { class: 'busy-line' }));
  view.resultCard = c;
  clear(view.result).append(c);
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
function figureCard() {
  view.tabBar = tabs(TABS, view.tab, (v) => { view.tab = v; view.ctx.setParams({ tab: v === 'matrix' ? null : v }); drawFigure(); });
  view.figHost = el('div', { class: 'fig-body' });
  return el('div', { class: 'card fig-card', hidden: true }, el('div', { class: 'busy-line' }), view.tabBar, view.figHost);
}

// The tab in hand, drawn afresh: a new tab replaces the old at once; the same tab, after a new
// alignment, keeps its drawing in view, dimmed, until the new one is ready. The partial order
// does not depend on the pair, and is kept across alignments.
async function drawFigure(afterAlign = false) {
  const tab = view.tab;
  if (afterAlign && tab === 'poa' && view.figShown === 'poa') return;
  const ticket = view.figTicket = ++serial;
  const host = view.figHost;
  const keep = view.figShown === tab && host.firstChild;
  const body = el('div');
  if (keep) host.firstChild.classList.add('is-busy'); else clear(host).append(body);
  view.figShown = tab;
  view.figure.classList.add('loading');
  const draw = { matrix: drawMatrix, algorithms: drawAlgorithms, semirings: drawSemirings, parametric: drawParametric, poa: drawPoa }[tab];
  try { await draw(body); } catch (e) { clear(body).append(callout(e.message, 'bad')); }
  if (ticket !== view.figTicket) return;
  if (keep) clear(host).append(body);
  view.figure.classList.remove('loading');
}

// A line that turns while the engine works, with the time so far once it is felt.
function waiting(text) {
  const t0 = performance.now(), elapsed = el('span', { class: 'elapsed' });
  const box = el('div', { class: 'wait' }, el('span', { class: 'spin' }), el('span', {}, text), elapsed);
  let seen = false;
  const tick = setInterval(() => {
    if (box.isConnected) seen = true; else if (seen || performance.now() - t0 > 6e5) { clearInterval(tick); return; }
    const t = (performance.now() - t0) / 1000;
    if (t > 0.5) elapsed.textContent = `${t.toFixed(1)} s`;
  }, 100);
  return box;
}

// A drawing with its caption and key beside it, or beneath when there is no room.
const split = (main, ...side) => el('div', { class: 'fig-split' }, el('div', { class: 'fig-main' }, main), el('div', { class: 'fig-side' }, ...side));
const keyItem = (swatch, label) => el('div', { class: 'key-item' }, swatch, el('span', {}, label));

// The DP matrix as a field, the path on it; WFA's wavefronts or the band over it.
function drawMatrix(body) {
  const r = view.last, M = r.matrix, sh = view.shown;
  if (!M) { body.append(callout('The matrix is drawn for sequences of up to 160,000 cells; these are longer. The other tabs still compute on them.')); return; }
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
    const top = el('div', { class: 'dp-top', style: { paddingLeft: `${size}px` } }, ...[...sh.b].map((c) => el('span', { style: { width: `${size}px` } }, c)));
    const left = el('div', { class: 'dp-left', style: { paddingTop: `${size}px` } }, ...[...sh.a].map((c) => el('span', { style: { height: `${size}px`, lineHeight: `${size}px` } }, c)));
    wrap.append(el('div'), top, left, canvas);
  } else wrap.append(canvas);
  canvas.addEventListener('pointermove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const j = Math.floor((e.clientX - rect.left) / size), i = Math.floor((e.clientY - rect.top) / size);
    if (i < 0 || j < 0 || i > n || j > m) { hover(null); return; }
    const v = values[i * (m + 1) + j];
    hover(`cell (${i}, ${j}): ${esc(sh.a[i - 1] || '·')} against ${esc(sh.b[j - 1] || '·')}<br>best score of a path to here: <b>${Number.isFinite(v) && v > -1e8 ? v : '−∞'}</b>`, e);
  });
  canvas.addEventListener('pointerleave', () => hover(null));
  const caption = M.fronts
    ? `The wavefronts of the WFA: each line joins the furthest cells that paths of one penalty reach on each diagonal, from penalty 0 (pale) to the optimum, ${r.penalty} (deep); between two penalties the fronts slide along the diagonals while the letters match. ${M.fronts.length} wavefronts (${M.fronts.length > 40 ? `every ${Math.ceil(M.fronts.length / 40)}th drawn` : 'all drawn'}), ${fmt.count(r.record.offsets)} offsets computed and ${fmt.count(r.record.extended)} matches walked, against the ${fmt.count((n + 1) * (m + 1))} cells of the full matrix, which lies faint beneath as Gotoh's field of best scores. The penalties are the scoring's own, by Eizenga and Paten's exchange: mismatch ${r.penalties.mismatch}, gap ${r.penalties.gap_open} + ${r.penalties.gap_extend}·L.`
    : 'Each cell holds the best score of a path from the start to it; the dark line is the path traced back from the best end.';
  // The key: the field's ramp, the path, and what the schedule adds over them.
  const line = (background) => el('span', { class: 'kline', style: { background } });
  const legend = el('div', { class: 'fig-legend' },
    el('div', {}, 'the best score of a path to each cell', el('div', { class: 'key-item' }, el('span', { class: 'kend' }, Number.isFinite(lo) ? String(lo) : '–'),
      el('span', { class: 'ramp', style: { background: `linear-gradient(90deg, ${[0.08, 0.3, 0.5, 0.7, 0.93].map(warm).join(', ')})`, opacity: M.fronts ? 0.45 : 1 } }),
      el('span', { class: 'kend' }, Number.isFinite(hi) ? String(hi) : '–'))),
    keyItem(line(token('--ink')), 'the optimal path, traced back'),
    M.band ? keyItem(line(token('--C')), `the band's edges: diagonals j − i from ${M.band[0]} to ${M.band[1]}`) : null,
    M.fronts ? keyItem(line('linear-gradient(90deg, rgba(120, 146, 190, 0.3), rgba(60, 76, 118, 0.8))'), `wavefronts, penalty 0 to ${r.penalty}`) : null,
    M.touched ? keyItem(el('span', { class: 'kdot' }), 'the cells the schedule touched') : null,
    el('div', { class: 'small faint' }, `${fmt.count(n + 1)} × ${fmt.count(m + 1)} cells, the query down, the target across; point at one to read it.`));
  body.append(split(wrap, el('p', { class: 'note', style: { marginTop: 0 } }, caption), legend));
}

// Every schedule on this pair: its value, its time, whether it agrees with the reference.
async function drawAlgorithms(body) {
  body.append(waiting('Running every schedule on this pair…'));
  let r;
  try { r = await api.run('align.compare', request(view.shown)); } catch (e) { clear(body).append(callout(e.message, 'bad')); return; }
  clear(body);
  // The rate at which each disposes of the whole edit graph: n·m over its time. A schedule that
  // computes fewer cells (wavefronts) or many at once (lanes, bit vectors) shows as faster.
  const rate = (x) => (x.seconds > 0 ? r.cells / x.seconds / 1e9 : 0);
  const gotoh = r.runs.find((x) => x.name === 'Gotoh');
  const table = el('table', { class: 'data wide' }, el('tr', {}, ...['schedule', 'computes', 'value', 'cells computed', 'time', 'n·m / time', 'agrees'].map((h) => el('th', { class: ['cells computed', 'time', 'n·m / time'].includes(h) ? 'num' : '' }, h))));
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
  const c = el('div', { class: 'al-bars' }, el('p', { class: 'eyebrow' }, 'billions of cells of the full matrix a second, against Gotoh (×)'));
  c.insertAdjacentHTML('beforeend', svg);
  body.append(el('p', { class: 'note', style: { margin: '0 0 12px' } }, 'Each schedule computes the optimum of the same sum, in its own order: scores maximised by Gotoh, Hirschberg and the stripes; penalties minimised by the wavefronts, whose global optimum is the same alignment (Eizenga and Paten 2022); unit costs by Myers\' bit vectors, a different scoring. Where the order is another (Myers; the ends-free WFA, which minimises the penalty rather than maximising the score), the value is not compared. Cells computed are the work done: the wavefronts\' offsets, Myers\' 64-row blocks, the stripes\' padded segments; n·m / time is the rate at which a schedule disposes of the whole edit graph.'),
    el('div', { style: { overflowX: 'auto' } }, table), c);
}

async function drawSemirings(body) {
  body.append(waiting('Summing over the paths in other semirings…'));
  let r;
  try { r = await api.run('align.semiring', request(view.shown)); } catch (e) { clear(body).append(callout(e.message, 'bad')); return; }
  clear(body);
  const mono = (x) => el('span', { class: 'mono' }, x);
  body.append(kv([
    ['all alignments (counting semiring)', r.alignments ? mono(r.alignments) : null],
    ['optimal alignments (tropical with multiplicities)', r.optimal ? mono(r.optimal) : null],
    ['Delannoy number D(m, n)', r.delannoy ? mono(r.delannoy) : null],
    ['log P(x, y), the pair HMM’s forward sum', r.forward !== undefined ? fmt.number(r.forward, 6) : null],
  ]));
  if (r.posterior) {
    const { rows: n, cols: m, values } = r.posterior;
    const size = Math.max(2, Math.min(14, Math.floor(Math.min(640 / m, 420 / n))));
    const canvas = el('canvas', { width: m * size, height: n * size, class: 'dp' });
    const g = canvas.getContext('2d');
    for (let i = 0; i < n; i++)
      for (let j = 0; j < m; j++) { g.fillStyle = warm(Math.sqrt(values[i * m + j])); g.fillRect(j * size, i * size, size, size); }
    body.append(el('p', { class: 'eyebrow', style: { marginTop: '16px' } }, 'the posterior of each pair of letters'),
      split(canvas, el('p', { class: 'note', style: { marginTop: 0 } }, 'The posterior probability that the i-th letter of the first sequence is aligned with the j-th of the second, under the pair HMM (Durbin et al. 1998, ch. 4): the probability semiring\'s forward and backward sums, divided. Where the alignment is certain the diagonal is dark; where it is not, the probability spreads.'),
        el('div', { class: 'fig-legend' }, el('div', {}, 'the probability, on a square-root scale', el('div', { class: 'key-item' }, el('span', { class: 'kend' }, '0'),
          el('span', { class: 'ramp', style: { background: `linear-gradient(90deg, ${[0, 0.25, 0.5, 0.75, 1].map(warm).join(', ')})` } }), el('span', { class: 'kend' }, '1'))))));
    if (r.mea) body.append(el('p', { class: 'eyebrow', style: { marginTop: '16px' } }, `the alignment of maximum expected accuracy: ${fmt.number(r.mea.expected_pairs, 4)} pairs expected correct, accuracy ${fmt.percent(r.mea.accuracy)}`), rows(r.mea.rows, view.shown.alphabet, 0, 0));
  }
}

async function drawParametric(body) {
  body.append(waiting('Propagating polygons through the edit graph…'));
  const sh = view.shown;
  let r;
  try { r = await api.run('align.parametric', { a: sh.a.slice(0, 120), b: sh.b.slice(0, 120), alphabet: sh.alphabet }); } catch (e) { clear(body).append(callout(e.message, 'bad')); return; }
  clear(body);
  body.append(el('p', { class: 'note', style: { margin: '0 0 12px' } }, `Every global alignment counts its mismatches x and its gap columns y; the scores match +1, mismatch −α, gap −β rank the alignments by αx + βy alone (the matches are fixed by x and y). The convex hull of all the (x, y) is the alignment polygon, computed in the polytope algebra, where ⊕ is the hull of a union and ⊗ a Minkowski sum (Pachter and Sturmfels 2004). Each vertex is the optimum for a cone of (α, β): the normal fan. Only the vertices whose cones meet the quadrant of positive penalties can win for a biologist's scoring; the others win where a mismatch or a gap is rewarded. ${r.vertices.length} vertices${sh.a.length > 120 || sh.b.length > 120 ? ', for the first 120 letters of each sequence' : ''}.`));
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
  const detail = el('div', { class: 'small al-vertex' }, el('p', { class: 'faint', style: { margin: 0 } }, 'Choose a vertex, or a point of the fan, to see its alignment.'));
  const showVertex = (k) => {
    const v = verts[k];
    if (!v) return;
    clear(detail).append(el('p', { class: 'eyebrow' }, el('span', { class: 'swatch', style: { background: CATEGORIES[k % CATEGORIES.length] } }), `vertex ${k + 1}: ${v.mismatches} mismatches, ${v.gaps} gap columns, ${v.matches} matches`), rows(v.rows, sh.alphabet, 0, 0));
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
  const out = el('div');
  let ticket = 0;
  // The reads drawn again: the earlier graph stays, dimmed, until the new one is made.
  const fill = async () => {
    const mine = ++ticket;
    if (out.firstChild) out.classList.add('is-busy'); else out.append(waiting('Fusing the reads into a graph…'));
    view.figure.classList.add('loading');
    let r;
    try { r = await api.run('align.poa', POA); } catch (e) { if (mine === ticket) { if (out.isConnected) view.figure.classList.remove('loading'); out.classList.remove('is-busy'); clear(out).append(callout(e.message, 'bad')); } return; }
    if (mine !== ticket) return;
    if (out.isConnected) view.figure.classList.remove('loading');
    out.classList.remove('is-busy');
    poaFigure(clear(out), r);
  };
  const reads = numberField({ value: POA.reads, min: 2, max: 40, onChange: (v) => { POA.reads = v; fill(); } });
  const length = numberField({ value: POA.length, min: 20, max: 600, step: 10, unit: 'bp', onChange: (v) => { POA.length = v; fill(); } });
  const error = select([0.02, 0.05, 0.08, 0.12, 0.15].map((e) => [e, fmt.percent(e, 0)]), POA.error, (v) => { POA.error = Number(v); fill(); });
  body.append(el('div', { class: 'fig-toolbar' },
    fld('Length', length, { inline: true }), fld('Reads', reads, { inline: true }), fld('Errors a base', error, { inline: true }),
    btn('Other reads', { icon: 'refresh', size: 'sm', onClick: () => { POA.seed++; fill(); } }),
    el('span', { class: 'grow' }),
    btn('Your own FASTA', { icon: 'toolbox', size: 'sm', kind: 'ghost', title: 'the Toolbox\'s multiple alignment, by partial orders, on sequences of your own', onClick: () => view.ctx.openTool('align-msa', { method: 'poa' }) })),
    out);
  await fill();
}

function poaFigure(body, r) {
  const d = r.read_distances, mean = d.reduce((a, b) => a + b, 0) / d.length;
  body.append(el('p', { class: 'note', style: { marginTop: 0 } }, `${d.length} reads of a random sequence of ${POA.length} bp, each with ${fmt.percent(POA.error, 0)} errors a base (substitutions, insertions and deletions, as a long-read instrument makes them): between ${Math.min(...d)} and ${Math.max(...d)} edits from the truth, ${fmt.number(mean, 3)} on average. Each read is aligned to the graph of those before it, by dynamic programming over the graph in a topological order, and fused into it: where it agrees it adds weight to a node, where it differs it adds a branch, aligned to the nodes it differs from. The graph has ${fmt.count(r.nodes.length)} nodes and ${fmt.count(r.edges.length)} edges in ${r.columns} columns, made in ${fmt.seconds(r.seconds)}. Its heaviest path, the consensus, is `,
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
  body.append(el('p', { class: 'eyebrow', style: { marginTop: '14px' } }, 'the graph: nodes by column, each read’s weight in the colour, the consensus ringed'), scroll,
    el('p', { class: 'eyebrow', style: { marginTop: '16px' } }, 'the multiple alignment it induces: each read with its edits from the truth, the consensus last; shaded, the columns where the reads disagree'), msa);
}
