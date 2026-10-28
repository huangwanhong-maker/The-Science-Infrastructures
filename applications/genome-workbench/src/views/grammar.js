// The Grammar view: the user's draft, "Towards a group context-sensitive grammar for genome
// analysis" (GE-CSG), made exact (domains/genomics/theory/group-grammars.md). Five sections, each
// computed by the engine on morphe and the genomics domain's grammar module:
//
//   1  the bases as a torsor of the Klein group: its square, its subgroups' cosets (IUPAC's codes),
//      its characters (the classifications) at a base point, and the chain of breakings;
//   2  parsing up to symmetry: a grammar and a group, the rule orbits by Burnside's lemma, Earley's
//      parse and tree, the equivariant parser beside the draft's orbit-only scan, spans shared, and
//      the canonical search of a monotone grammar;
//   3  the strand symmetry of a genome: the eight symmetries of words, and the Wedderburn channels;
//   4  the genetic code on (C₂)⁶: invariance profile, stabiliser and power spectrum;
//   5  the characters of a finite group, exactly, with its Wedderburn blocks.
//
// Every figure is drawn by draw/grammar.js. The state of each section (and its last answer) is kept
// across a remount, so that a change of theme redraws without asking again.

import { el, clear, hover, esc } from '../dom.js';
import * as api from '../api.js';
import { icon } from '../ui/icons.js';
import { btn, iconBtn, seg, select, toggle, fld, card, metrics, callout, empty, working, onEnter, dialog, seqEditor } from '../ui/kit.js';
import { CODES } from '../tools/common.js';
import * as G from '../draw/grammar.js';

let view = null;
// The latest request of each section, kept outside the view so that an answer arriving after a remount
// (a change of theme) is still drawn.
const tickets = {};

const NOTE = '/domains/genomics/theory/group-grammars.md';
const fmt = (x) => (Number.isFinite(x) ? x.toLocaleString('en-US') : '–');
const pct = (x, d = 1) => (Number.isFinite(x) ? `${(100 * x).toFixed(d)}%` : '–');
const sci = (x) => (x === 0 ? '0' : x < 1e-3 ? x.toExponential(1).replace('e', ' × 10^').replace(/\^(-?\d+)/, (m, e) => `${e.replace('-', '⁻').replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d])}`) : x.toPrecision(3));
const yes = (b, good = 'yes', bad = 'no') => el('span', { class: `badge ${b ? 'good' : 'bad'}` }, b ? good : bad);

// ---- the presets of §2 --------------------------------------------------------------------------------
// A word made to fold: a stem, a loop, and the stem's partner read back.
const revcomp = (s) => s.split('').reverse().map((c) => ({ A: 'T', C: 'G', G: 'C', T: 'A' }[c])).join('');
const PAIRS = ['A Q T', 'T Q A', 'C Q G', 'G Q C'];
const hairpinRules = (wobble) => [
  `S -> ${['A P T', 'T P A', 'C P G', 'G P C', ...(wobble ? ['G P T', 'T P G'] : [])].join(' | ')}`,
  `P -> ${[...PAIRS, ...(wobble ? ['G Q T', 'T Q G'] : [])].join(' | ')}`,
  `Q -> ${[...PAIRS, 'A L T', 'T L A', 'C L G', 'G L C', ...(wobble ? ['G Q T', 'T Q G', 'G L T', 'T L G'] : [])].join(' | ')}`,
  'L -> N N N | N N N N | N N N N N',
  'N -> A | C | G | T',
].join('\n');

const PRESETS = [
  { id: 'scan', title: 'S → aa | bb under a ↔ b', alphabet: 'ab', rules: 'S -> a a | b b', group: [['ba', false]], word: 'ab',
    about: 'The draft\'s example (§2.1): invariant, its two rules one orbit. The orbit-only scan matches a terminal by its orbit alone and accepts ab, which the language does not contain; the equivariant parser binds the element and does not.' },
  { id: 'brackets', title: 'Balanced brackets under the mirror', alphabet: '()', rules: 'S -> ( ) | ( S ) | S S', group: [[')(', true]], word: '(()())()',
    about: 'The mirror reverses a word and exchanges the brackets. It keeps every rule, so the three rules are three orbits: |R/G| = 3, above |R|/|G| = 1.5 (Burnside, Proposition 1.4). It reverses, so the equivariant parser has nothing to bind; the spans still share by the mirror.' },
  { id: 'hairpin', title: 'A DNA hairpin under the strand group', alphabet: 'ACGT', rules: hairpinRules(false), group: [['TGCA', false], ['ACGT', true]], word: `GCA${'TTCG'}${revcomp('GCA')}`,
    about: 'A stem of at least three Watson–Crick pairs closing a loop of three to five (§5). The grammar is invariant under complement κ, reversal ρ and their product κρ. The equivariant parser works with {1, κ}, the elements that act by letters alone. The word is made to fold: the stem GCA, a loop, the stem read back.' },
  { id: 'wobble', title: 'The hairpin with wobble pairs', alphabet: 'ACGT', rules: hairpinRules(true), group: [['TGCA', false], ['ACGT', true]], word: 'GCATTCGTGT',
    about: 'G·T pairs added (G·U in RNA): reversal still keeps the rules, but the complement does not, since κ(G·T) = C·A pairs in neither (Proposition 5.1). The chemistry breaks the symmetry, and the stabiliser records the part that survives. The word closes its stem on a G·T pair.' },
  { id: 'anbncn', title: 'aⁿbⁿcⁿ, monotone', alphabet: 'abc', rules: 'S -> a S B C | a B C\nC B -> B C\na B -> a b\nb B -> b b\nb C -> b c\nc C -> c c', group: [], word: 'aabbcc',
    about: 'The classic context-sensitive language (Hopcroft and Ullman 1979, §9). Earley\'s parser does not apply: membership is decided by searching the forms no longer than the word (§2.4), exponentially in the worst case.' },
  { id: 'wcd', title: 'w cⁿ dⁿ under a ↔ b, monotone', alphabet: 'abcd', rules: 'S -> a S C D | b S C D | a C D | b C D\nD C -> C D\na C -> a c\nb C -> b c\nc C -> c c\nc D -> c d\nd D -> d d', group: [['bacd', false]], word: 'abccdd',
    about: 'A monotone grammar with a symmetry: w over {a, b} followed by as many c\'s and d\'s as w has letters. It is invariant under a ↔ b, so the search may visit the forms\' classes instead of the forms: between |forms|/|G| and |forms| of them (Proposition 2.4). The complexity class does not change.' },
];

const GROUPS = [
  ['cyclic 4', 'C₄, cyclic of order 4'], ['klein', 'V₄ = C₂ × C₂, the Klein group'], ['cyclic 6', 'C₆'], ['symmetric 3', 'S₃'], ['dihedral 4', 'D₄, of order 8'], ['quaternion', 'Q₈, the quaternions'],
  ['alternating 4', 'A₄'], ['symmetric 4', 'S₄'], ['dihedral 5', 'D₅, of order 10'], ['cyclic 3 × cyclic 3', 'C₃ × C₃'], ['dihedral 4 × cyclic 2', 'D₄ × C₂'], ['alternating 5', 'A₅, the icosahedral group'], ['symmetric 5', 'S₅'],
];

// What the view holds between remounts.
const kept = {
  bp: 'A', subgroup: '⟨τ⟩', chain: 0, nuc: null,
  preset: 'scan', alphabet: 'ab', rules: 'S -> a a | b b', gens: [['ba', false]], start: '', word: 'ab', use: 'stabiliser', parse: null, parseError: null,
  source: 'genome', pasted: '', k: 4, sym: null, symError: null, symFor: null,
  code: 1, element: null, codeRes: null,
  group: 'dihedral 4', compare: 'quaternion', chars: {},
  started: false,
};

export function mount(section, ctx) {
  view = { ctx, section };
  const p = ctx.params;
  if (p.bp && 'ACGT'.includes(p.bp)) kept.bp = p.bp;
  if (p.preset && PRESETS.some((x) => x.id === p.preset) && p.preset !== kept.preset) usePreset(p.preset, false);
  if (p.k && +p.k >= 1 && +p.k <= 6) kept.k = +p.k;
  if (p.code && CODES.some(([v]) => String(v) === p.code)) kept.code = +p.code;
  if (p.group) kept.group = p.group;
  if (p.vs !== undefined) kept.compare = p.vs;
  section.append(
    el('p', { class: 'lede' }, 'W. Huang\'s draft, ', el('i', {}, 'Towards a Group Context-Sensitive Grammar for Genome Analysis'), ' (2026), takes a finite group as the primitive of a grammar of genomes. Here it is made exact, as the theory note states it. The bases are a torsor of the Klein group: its subgroups\' cosets are IUPAC\'s codes, its characters the classical classifications. A grammar invariant under a group can be parsed up to that symmetry if each rule\'s group element is bound as the parse goes; the draft\'s orbit-only scan accepts too much. Of the eight symmetries of words, a genome\'s k-mer statistics keep only the reverse complement: Chargaff\'s second rule, found as a group. The genetic code is a function on (C₂)⁶, and the vertebrate mitochondrial code keeps one symmetry more than the standard one. The Wedderburn blocks of a group algebra, computed exactly, show that the algebra alone does not determine the group.'),
    contents(),
    basesSection(), parseSection(), strandSection(), codeSection(), charactersSection());
  drawBases();
  drawParse();
  drawStrands();
  drawCode();
  drawCharacters();
}

export function show(ctx) {
  const genomeChanged = view.ctx.genome !== ctx.genome;
  view.ctx = ctx;
  if (!kept.started) {
    kept.started = true;
    runBases();
    runParse();
    runStrands();
    runCode();
    runCharacters();
  } else if (genomeChanged && kept.source === 'genome') runStrands();
  view.genomeName.textContent = genomeTitle();
}

const genomeTitle = () => (view.ctx.genomeInfo()?.title || 'the genome in hand').replace(/, complete genome.*$/, '').replace(/ \(.*\)$/, '');

// ---- the head of the page and of each section --------------------------------------------------------
function contents() {
  const go = (id) => view.section.querySelector(`#gr-${id}`)?.scrollIntoView({ behavior: 'smooth', block: 'start' });
  const chips = [['bases', '1 · the bases'], ['parse', '2 · parsing'], ['strands', '3 · strands'], ['code', '4 · the code'], ['characters', '5 · characters']]
    .map(([id, label]) => el('button', { class: 'chip', type: 'button', onclick: () => go(id) }, label));
  return el('div', { class: 'gr-contents' }, icon('grammar', 'sm'), el('span', { class: 'small faint' }, 'Group grammars for genomes, made exact:'), ...chips,
    el('span', { style: { flex: 1 } }),
    btn('Read the theory note', { icon: 'file', size: 'sm', kind: 'ghost', onClick: readNote }),
    el('a', { class: 'small link', href: NOTE, target: '_blank', rel: 'noopener', title: 'the note as Markdown, from the repository' }, 'as a file'));
}

function sectionHead(n, title, sub, ...actions) {
  return el('div', { class: 'section-head' }, el('h3', {}, `${n}. ${title}`), el('span', { class: 'sub' }, sub), actions.length ? el('div', { class: 'actions' }, ...actions) : null);
}

// A stage that keeps its last figures, dimmed, while a new answer is computed, with a line of what
// is being done and for how long above them.
function busy(box, text) {
  const kids = [...box.children];
  kids.forEach((k) => k.classList.add('stale'));
  const line = el('div', { class: 'gr-wait' }, working(text));
  box.prepend(line);
  return () => { line.remove(); kids.forEach((k) => k.classList.remove('stale')); };
}

// A quiet line that sends the user to the Toolbox with their own data.
function toolboxHint(text, id, label, ic = 'toolbox') {
  return el('p', { class: 'toolbox-hint' }, icon(ic, 'sm'), el('span', {}, text, el('button', { type: 'button', class: 'link', onclick: () => view.ctx.openTool(id) }, label), ', in the Toolbox.'));
}
// Headline numbers in one row that fills the stage.
const tiles = (items) => { const m = metrics(items); m.classList.add('gr-metrics'); return m; };
const supNum = (n) => String(n).replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);

// ---- 1. the bases -------------------------------------------------------------------------------------
function basesSection() {
  view.basesBody = el('div', { class: 'gr-three' });
  return el('section', { id: 'gr-bases', class: 'gr-section' },
    sectionHead(1, 'The bases as a group', 'a torsor of the Klein group: its subgroups are IUPAC\'s codes, its characters the classifications (§3)'),
    card({ body: view.basesBody }));
}

async function runBases() {
  const t = (tickets.bases = (tickets.bases || 0) + 1);
  try {
    const r = await api.run('grammar.nucleotides', { base_point: kept.bp });
    if (t !== tickets.bases) return;
    kept.nuc = r;
  } catch (e) {
    kept.nuc = { error: e.message };
  }
  drawBases();
}

function drawBases() {
  const box = clear(view.basesBody);
  const r = kept.nuc;
  if (!r) { box.append(working('the Klein group of the bases')); return; }
  if (r.error) { box.append(callout(r.error, 'bad')); return; }
  const sg = r.subgroups.find((s) => s.name === kept.subgroup) || r.subgroups[1];
  // the square
  const square = el('div', { class: 'gr-figure' });
  square.innerHTML = G.kleinSquareSvg({ basePoint: r.base_point, torsor: r.torsor, cosets: sg.cosets });
  const subChips = el('div', { class: 'chips' }, ...r.subgroups.map((s) => el('button', { class: `chip${s.name === sg.name ? ' on' : ''}`, type: 'button', title: `{${s.elements.join(', ')}}`, onclick: () => { kept.subgroup = s.name; drawBases(); } }, s.name)));
  const cosetText = sg.name === '1' ? 'The trivial subgroup: each base is a coset of its own.' : sg.name === 'V₄' ? 'The whole group: one coset, N, any base.'
    : `${sg.name} = {${sg.elements.join(', ')}}: its cosets are ${sg.cosets.map((c) => `${c.letter} = {${c.bases.split('').join(', ')}}`).join(' and ')}${(() => { const c = r.classifications.find((x) => `⟨${x.kernel}⟩` === sg.name); return c ? `, the ${c.name} classes` : ''; })()}.`;
  const legend = el('div', { class: 'gr-legend' },
    el('span', {}, el('i', { class: 'k-kappa' }), 'κ complement, A↔T and C↔G'), el('span', {}, el('i', { class: 'k-tau' }), 'τ transition, A↔G and C↔T'), el('span', {}, el('i', { class: 'k-sigma' }), 'σ = κτ, A↔C and G↔T'));
  const col1 = el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, 'The group on the bases'), square, legend,
    fld('A subgroup, and its cosets', subChips), el('p', { class: 'note' }, cosetText),
    el('p', { class: 'note' }, 'The action is free and transitive: the bases are a V₄-torsor (Proposition 3.1). By each base stands the element that carries the base point to it; the three-letter codes B, D, H, V are cosets of no subgroup (Proposition 3.2).'));
  // the characters
  const sign = (v) => el('td', { class: `gr-sign ${v > 0 ? 'plus' : 'minus'}` }, v > 0 ? '+1' : '−1');
  const rows = [{ name: 'trivial', kernel: 'V₄', on_group: [1, 1, 1, 1], on_bases: [1, 1, 1, 1], letters: ['N', ''] }, ...r.classifications];
  const groupTable = el('table', { class: 'data gr-chars' },
    el('tr', {}, el('th', {}, 'character'), ...r.elements.map((e) => el('th', { class: 'num', title: e.about }, e.name))),
    ...rows.map((c) => el('tr', { class: `${`⟨${c.kernel}⟩` === sg.name ? 'on' : ''}`, onclick: () => { if (c.kernel !== 'V₄') { kept.subgroup = `⟨${c.kernel}⟩`; drawBases(); } } },
      el('td', {}, c.name, c.kernel !== 'V₄' ? el('span', { class: 'faint small' }, `  kernel ⟨${c.kernel}⟩`) : null), ...c.on_group.map(sign))));
  const bpSeg = seg('ACGT'.split('').map((b) => [b, b]), kept.bp, (b) => { kept.bp = b; view.ctx.setParams({ bp: b }); runBases(); });
  const baseTable = el('table', { class: 'data gr-chars' },
    el('tr', {}, el('th', {}, `at the base point ${r.base_point}`), ...'ACGT'.split('').map((b) => el('th', { class: 'num' }, el('span', { class: `gr-base b-${b}` }, b)))),
    ...rows.map((c) => el('tr', { class: `${`⟨${c.kernel}⟩` === sg.name ? 'on' : ''}` }, el('td', {}, c.name, c.letters?.[1] ? el('span', { class: 'faint small' }, `  +1 on ${c.letters[0]}`) : null), ...c.on_bases.map(sign))));
  const col2 = el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, 'Its characters, the classifications'), groupTable,
    el('p', { class: 'note' }, 'Each nontrivial character is ±1 with a subgroup of order 2 for its kernel; together they form a group, χ', el('sub', {}, 'RY'), 'χ', el('sub', {}, 'WS'), ' = χ', el('sub', {}, 'MK'), ' ', yes(r.product_rule, 'holds', 'fails'), ', so two classifications name a base (Proposition 3.3).'),
    fld('The base point of the breaking', bpSeg), baseTable,
    el('p', { class: 'note' }, r.base_point === 'A' ? 'Through the torsor each character becomes a function on the bases, χ(g·σ₀) = χ(g). Choose another base point: the signs move, the partition does not.'
      : `Moving the base point from A to ${r.base_point} multiplies each character by its value at the element that carries A to ${r.base_point}. The signs flip where that value is −1; the partition into classes stays. Which class a character calls +1 is the one choice a breaking makes.`));
  // the chain of breakings
  const chain = r.chains[kept.chain] || r.chains[0];
  const chainSeg = seg(r.chains.map((c, k) => [k, c.middle.join('/'), c.name]), kept.chain, (k) => { kept.chain = k; drawBases(); });
  const tree = el('div', { class: 'gr-figure' });
  tree.innerHTML = G.breakingTreeSvg(chain);
  const col3 = el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, 'The chain of breakings'), fld('Through', chainSeg), tree,
    el('p', { class: 'note' }, `V₄ ⊇ ⟨${chain.kernel}⟩ ⊇ 1 gives the alphabets N, then ${chain.middle.join(' and ')} (the ${chain.name} classes), then the bases: each step the orbit map of a smaller subgroup. The residual symmetry on the middle alphabet is ${chain.residual}; on the bases it is V₄ again (§3.4).`),
    el('p', { class: 'note' }, 'The minimal idempotents of ℂ[V₄] are the characters, not the letters: the classifications arise as idempotents, the letters from the torsor\'s trivialisation, and the two are Fourier dual (Proposition 3.4).'));
  box.append(col1, col2, col3);
}

// ---- 2. parsing up to symmetry ---------------------------------------------------------------------------
function usePreset(id, run = true) {
  const p = PRESETS.find((x) => x.id === id);
  if (!p) return;
  Object.assign(kept, { preset: id, alphabet: p.alphabet, rules: p.rules, gens: p.group.map((g) => [...g]), start: '', word: p.word, use: 'stabiliser' });
  if (run && view) { view.ctx.setParams({ preset: id }); drawParseControls(); runParse(); }
}

function parseSection() {
  view.parseControls = el('div', { class: 'controls' });
  view.parseStage = el('div', { class: 'stage' });
  drawParseControls();
  return el('section', { id: 'gr-parse', class: 'gr-section' },
    sectionHead(2, 'Parsing up to symmetry', 'a grammar invariant under a group, parsed with one rule per orbit (§1, §2)'),
    el('div', { class: 'bench wide-controls' }, view.parseControls, view.parseStage));
}

function drawParseControls() {
  const box = clear(view.parseControls);
  const preset = PRESETS.find((p) => p.id === kept.preset);
  const presetSel = select([...PRESETS.map((p) => [p.id, p.title]), ...(preset ? [] : [['', 'your own grammar']])], preset ? kept.preset : '', (id) => usePreset(id));
  const alphabet = el('input', { type: 'text', value: kept.alphabet, spellcheck: false, class: 'mono', oninput: () => { kept.alphabet = alphabet.value.trim(); mine(); drawGens(); } });
  const rules = el('textarea', { rows: Math.min(9, Math.max(3, kept.rules.split('\n').length + 1)), spellcheck: false, class: 'gr-rules', oninput: () => { kept.rules = rules.value; mine(); } });
  rules.value = kept.rules;
  const start = el('input', { type: 'text', value: kept.start, placeholder: 'the first rule\'s left side', spellcheck: false, oninput: () => { kept.start = start.value.trim(); } });
  const word = el('input', { type: 'text', value: kept.word, spellcheck: false, class: 'mono', oninput: () => { kept.word = word.value; } });
  const gens = el('div', { class: 'stack gr-gens' });
  const drawGens = () => {
    clear(gens);
    kept.gens.forEach((g, k) => {
      const t = el('input', { type: 'text', value: g[0], spellcheck: false, class: 'mono', placeholder: kept.alphabet, oninput: () => { g[0] = t.value.trim(); mine(); } });
      gens.append(el('div', { class: 'gr-gen' }, el('span', { class: 'gr-gen-from mono', title: 'the letters of the alphabet, and their images' }, `${kept.alphabet} ↦`), t,
        toggle('reverses', g[1], (v) => { g[1] = v; mine(); }, { title: 'the generator also reverses words (an anti-automorphism, sign −1)' }),
        iconBtn('x', { title: 'remove this generator', size: 'sm', onClick: () => { kept.gens.splice(k, 1); mine(); drawGens(); } })));
    });
    if (!kept.gens.length) gens.append(el('p', { class: 'small faint', style: { margin: 0 } }, 'No generators: the trivial group.'));
  };
  drawGens();
  const mine = () => { if (kept.preset) { kept.preset = ''; presetSel.value = ''; aboutLine.textContent = 'Your own grammar.'; } };
  const aboutLine = el('div', { class: 'hint' }, preset ? preset.about : 'Your own grammar.');
  const runBtn = btn('Parse', { icon: 'play', kind: 'primary', onClick: () => runParse() });
  view.parseRunBtn = runBtn;
  const card_ = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'A grammar'), fld('Example', presetSel), aboutLine,
      el('div', { class: 'ctl-row' }, fld('Alphabet', alphabet, { help: 'its letters, as one word' }), fld('Start', start)),
      fld('Rules', rules, { help: 'one to a line: S -> a S b | ε', hint: 'Symbols are separated by spaces; a letter of the alphabet is a terminal, anything else a nonterminal.' })),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'A group on its letters'), gens,
      el('div', { class: 'row' }, btn('add a generator', { icon: 'plus', size: 'sm', kind: 'ghost', onClick: () => { kept.gens.push([kept.alphabet, false]); mine(); drawGens(); } })),
      el('div', { class: 'hint' }, 'Each generator is the images of the letters in order; a reversing one acts by anti-automorphism (Definition 1.1). They are closed to a group.')),
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'A word, and what to do without invariance'),
      fld('Word', word),
      fld('When the rules are not invariant, use', seg([['stabiliser', 'their stabiliser'], ['closure', 'their closure']], kept.use, (v) => { kept.use = v; }, { block: true }),
        { hint: 'The stabiliser keeps the language and loses symmetry; the closure keeps the group and changes the language.' })),
    el('div', { class: 'run-bar' }, runBtn));
  onEnter(card_, () => runParse());
  box.append(card_, toolboxHint('A grammar on words of your own: ', 'grammar-parse', 'Parse with a grammar'),
    toolboxHint('Every hairpin of a stretch of your own sequence, by the inside algorithm over all its spans: ', 'grammar-hairpins', 'Hairpins by grammar', 'annotation'));
}

async function runParse() {
  const t = (tickets.parse = (tickets.parse || 0) + 1);
  let request;
  try {
    request = { alphabet: kept.alphabet, rules: G.rulesFromText(kept.rules), group: kept.gens.filter((g) => g[0]), word: kept.word, use: kept.use };
    if (kept.start) request.start = kept.start;
  } catch (e) {
    kept.parse = null; kept.parseError = e.message; drawParse(); return;
  }
  const done = view.parseStage.children.length ? busy(view.parseStage, 'parsing') : null;
  if (!done) view.parseStage.replaceChildren(working('parsing'));
  const b = view.parseRunBtn;
  if (b) b.disabled = true;
  try {
    const r = await api.run('grammar.parse', request);
    if (t !== tickets.parse) return;
    kept.parse = r; kept.parseError = null;
  } catch (e) {
    if (t !== tickets.parse) return;
    kept.parseError = e.message;
  } finally {
    if (t === tickets.parse) { done?.(); if (b) b.disabled = false; }
  }
  drawParse();
}

function drawParse() {
  const box = clear(view.parseStage);
  if (kept.parseError) box.append(callout(kept.parseError, 'bad'));
  const r = kept.parse;
  if (!r) { if (!kept.parseError) box.append(card({ body: empty({ icon: 'grammar', title: 'A grammar and a group', text: 'Choose an example or write your own, then parse.' }) })); return; }
  box.append(...parseCards(r));
}

// The answer of grammar.parse as cards (the view's and the toolbox's).
export function parseCards(r) {
  const items = [];
  const box = { append: (...xs) => items.push(...xs.filter(Boolean)) };
  const o = r.orbits;
  box.append(tiles([
    [fmt(o.rules), `rules |R| of the grammar ${r.used === 'closure' ? 'closed under G' : 'used'}`],
    [fmt(o.group_order), `|G|, the elements that keep the rules${r.invariant ? '' : r.used === 'closure' ? ' (all of them, after the closure)' : ' (the stabiliser)'}`],
    [fmt(o.count), '|R/G|, the rule orbits: one representative of each is kept'],
    [fmt(o.burnside), '(1/|G|) Σ |Fix(g)| over the elements, Burnside\'s count of the orbits', o.burnside === o.count ? 'good' : ''],
    [o.ratio.toFixed(o.ratio % 1 ? 2 : 0), '|R|/|G|: the orbits number at least this (Proposition 1.4)'],
  ]));
  // invariance
  const chips = el('div', { class: 'chips' }, ...r.group.elements.map((e, k) => el('span', { class: `chip${r.keeps[k] ? ' on' : ''}`, title: `${e.letters}${e.reverses ? ', reversing' : ''}; order ${e.order}` },
    r.keeps[k] ? icon('check', 'sm') : icon('x', 'sm'), e.name)));
  const inv = [el('div', { class: 'row' }, yes(r.invariant, 'invariant', 'not invariant'), el('span', { class: 'small faint' }, `under the ${r.group.order} elements the generators make`)), chips];
  if (!r.invariant) {
    inv.push(el('p', { class: 'note' }, `The elements that keep the rules form the stabiliser {${r.stabiliser.join(', ')}}. Their closure under the whole group would add ${fmt(r.closure.added_count)} rule${r.closure.added_count === 1 ? '' : 's'}${r.closure.added.length ? `: ${r.closure.added.slice(0, 8).join(';  ')}${r.closure.added.length > 8 ? ' …' : ''}` : ''}.`));
    inv.push(el('p', { class: 'note' }, r.used === 'closure' ? 'What follows is computed on the closure, under the whole group: its language is the smallest invariant one containing the grammar\'s.' : 'What follows is computed on the grammar as written, under its stabiliser: the symmetry that survives.'));
  } else inv.push(el('p', { class: 'note' }, 'g·R = R for every element, so L(Γ) is invariant, and a nonterminal fixed by the group derives u exactly when it derives g·u (Proposition 1.3).'));
  box.append(card({ title: 'Invariance', sub: `the rules and their images under ${r.acting.order === r.group.order ? 'the group' : 'each element'}`, body: el('div', { class: 'card-body stack' }, ...inv) }));

  if (r.earley) drawParsers(box, r);
  else if (r.monotone_search) drawMonotone(box, r);
  else if (r.note) box.append(callout(r.note, 'warn'));

  // the orbits
  const list = el('div', { class: 'gr-orbits' }, ...o.list.map((x) => el('div', { class: 'gr-orbit' }, el('span', { class: 'badge accent' }, `${x.size}`), el('span', { class: 'mono' }, x.representative),
    x.size > 1 ? el('span', { class: 'faint small mono' }, x.members.filter((m) => m !== x.representative).join('   ')) : el('span', { class: 'faint small' }, 'fixed by every element'))));
  box.append(card({ title: 'The rule orbits', sub: `${fmt(o.count)} orbits of ${fmt(o.rules)} rules under ${fmt(o.group_order)} elements; the size of each beside it`, body: el('div', { class: 'card-body' }, list) }));
  return items;
}

function drawParsers(box, r) {
  const e = r.earley, q = r.equivariant;
  const wrong = q.scan.accepts && !e.accepts;
  const acc = (b) => el('span', { class: `badge ${b ? 'good' : ''}` }, b ? 'accepts' : 'rejects');
  const rows = [
    ['Earley\'s, on every rule', `${fmt(r.orbits.rules)} rules`, acc(e.accepts), fmt(e.items), 'the chart, Σ |Sⱼ|'],
    ['equivariant, one rule per orbit, its element bound', `${fmt(q.orbits)} orbit${q.orbits === 1 ? '' : 's'} under {${q.group.join(', ')}}`, acc(q.consistent.accepts), fmt(q.consistent.items), `standing for ${fmt(q.consistent.plain_equivalent)} of Earley's`],
    ['the draft\'s scan, a terminal matched by its orbit', `${fmt(q.orbits)} orbit${q.orbits === 1 ? '' : 's'}`, el('span', {}, acc(q.scan.accepts), wrong ? el('span', { class: 'badge bad', style: { marginLeft: '6px' } }, 'outside the language') : null), fmt(q.scan.items), ''],
  ];
  const t = el('table', { class: 'data wide gr-parsers' }, el('tr', {}, ...['parser', 'holds', `the word ${r.word || 'ε'}`, 'items', ''].map((h, k) => el('th', { class: k === 3 ? 'num' : '' }, h))),
    ...rows.map((row, k) => el('tr', { class: k === 2 && wrong ? 'gr-wrong' : '' }, ...row.map((c, j) => el('td', { class: j === 3 ? 'num' : j === 4 ? 'faint small' : '' }, c)))));
  const lines = [];
  if (q.order === 1) lines.push(el('p', { class: 'note' }, 'No element but 1 acts by letters alone here, so the equivariant parser has nothing to bind and the scan matches letters exactly: both are Earley\'s parser. Reversing elements still serve the rule orbits and the shared spans below.'));
  else {
    lines.push(el('p', { class: 'note' }, `The equivariant parser predicts one item per orbit and narrows each item's set of elements C at every scan to {g ∈ C : g·a = x}. Its chart stands for Earley's exactly: ${fmt(q.consistent.plain_equivalent)} plain items against ${fmt(e.items)} `, yes(q.chart_agrees, 'equal', 'unequal'), ' (Proposition 2.3). The worst case stays O(n³ |R|); what shrinks is the rule store and the predictions.'));
    lines.push(el('p', { class: 'note' }, 'The scan drops the binding: it accepts π⁻¹(L(Γ/G)), the words whose letters\' orbits spell a word of the orbit grammar, which contains L(Γ) and in general strictly (Proposition 2.1).'));
    if (q.survey) {
      const s = q.survey;
      lines.push(el('p', { class: 'note gr-survey' }, `Of all ${fmt(s.words)} words up to length ${s.up_to}: ${fmt(s.in_language)} in the language; the scan accepts ${fmt(s.by_scan)}, ${fmt(s.outside)} of them outside it`,
        s.examples.length ? el('span', {}, ', such as ', ...s.examples.flatMap((x, k) => [k ? ', ' : '', el('span', { class: 'mono' }, x)])) : '', '. The equivariant parser agrees with Earley\'s on every one ', yes(s.equivariant_agrees, 'yes', 'no'), '.'));
    }
  }
  box.append(card({ title: 'Three parsers on one word', sub: `the subgroup that acts by letters alone, {${q.group.join(', ')}}, binds the elements`, body: el('div', { class: 'card-body stack' }, t, ...lines) }));
  // the tree
  if (e.tree) {
    const fig = el('div', { class: 'gr-figure gr-scroll' });
    const tr = G.treeSvg(e.tree, { maxWidth: 1050 });
    fig.innerHTML = tr.svg;
    box.append(card({ title: 'The parse tree', sub: `one parse of ${r.word || 'ε'}, read back from Earley's chart; hover a node for its rule`, body: el('div', { class: 'card-body' }, fig) }));
  } else if (!e.accepts) {
    box.append(callout(`Earley's parser rejects ${r.word || 'the empty word'}: it is not in the language${r.equivariant?.scan?.accepts ? ', although the orbit-only scan accepts it' : ''}.`, wrong ? 'warn' : ''));
  }
  // spans shared
  if (r.spans && !r.spans.skipped) {
    const g = r.spans.by_group, c = r.spans.by_content;
    const fig = el('div', { class: 'gr-figure' });
    fig.innerHTML = G.hbarsSvg([
      { label: 'by the group', sub: `${fmt(g.classes)} classes of spans`, value: g.reused, note: `${fmt(g.reused)} of ${fmt(g.spans)} spans reused`, colour: 'var(--accent)' },
      { label: 'by content alone', sub: `${fmt(c.classes)} classes`, value: c.reused, note: `${fmt(c.reused)} of ${fmt(c.spans)}`, colour: 'var(--ink-faint)' },
    ], { width: 620, row: 38, label: 150, max: Math.max(1, g.spans) });
    box.append(card({ title: 'Spans decided once per class', sub: 'the inside algorithm, each span\'s original nonterminals looked up by its class (§2.3)', body: el('div', { class: 'card-body' }, fig,
      el('p', { class: 'note' }, 'A nonterminal fixed by the group derives a span exactly when it derives its images, so a span whose class was decided before is not decided again. With the trivial group the class is the content, and repeated substrings share; the group widens the sharing to their images. This is exact, and it is where the draft\'s shared parse caches hold as stated.')) }));
  } else if (r.spans?.skipped) {
    box.append(callout(`Spans are not shared here: ${r.spans.skipped} (the inside algorithm takes ε-free grammars).`));
  }
}

function drawMonotone(box, r) {
  const m = r.monotone_search;
  const fig = el('div', { class: 'gr-figure' });
  const rows = [];
  if (!m.classes.stopped) rows.push({ label: 'classes', sub: 'forms up to the group', value: m.classes.visited, note: `${fmt(m.classes.visited)} visited`, colour: 'var(--accent)' });
  if (!m.forms.stopped) rows.push({ label: 'forms', sub: 'every sentential form', value: m.forms.visited, note: `${fmt(m.forms.visited)} visited`, colour: 'var(--ink-faint)' });
  fig.innerHTML = rows.length ? G.hbarsSvg(rows, { width: 620, row: 38, label: 150 }) : '';
  const derives = m.classes.stopped ? m.forms.derives : m.classes.derives;
  const stopped = m.classes.stopped || m.forms.stopped;
  box.append(card({ title: 'The search of a monotone grammar', sub: 'every sentential form no longer than the word, breadth first from S (§2.4)', body: el('div', { class: 'card-body stack' },
    el('div', { class: 'row' }, el('span', { class: 'small' }, `${r.word}:`), derives === undefined ? el('span', { class: 'badge warn' }, 'undecided') : el('span', { class: `badge ${derives ? 'good' : ''}` }, derives ? 'derived' : 'not derived'),
      el('span', { class: 'small faint' }, `the grammar is monotone and not context-free, so Earley's parser does not apply`)),
    fig,
    stopped ? callout(`The search stopped after ${fmt(m.classes.budget || m.forms.budget)} forms; try a shorter word.`, 'warn') : null,
    el('p', { class: 'note' }, `A monotone grammar derives w only through forms no longer than w, so membership is decided by searching them, in exponential time in the worst case (PSPACE-complete in general). When the grammar is invariant and the group fixes S, the canonical search visits the classes of the derivable forms instead: between |forms|/|G| and |forms| of them, here ${m.classes.visited && m.forms.visited ? (m.forms.visited / m.classes.visited).toFixed(2) : '–'} times fewer under |G| = ${r.orbits.group_order} (Proposition 2.4). The wall is lowered, not restructured.`)) }));
}

// ---- 3. the strands ----------------------------------------------------------------------------------------
function strandSection() {
  view.genomeName = el('span');
  view.strandStage = el('div', { class: 'gr-two' });
  view.pasteBox = seqEditor({ value: kept.pasted, rows: 4, placeholder: 'paste a sequence or FASTA, or drop a file here', onChange: (v) => { kept.pasted = v; } });
  view.pasteBox.hidden = kept.source !== 'paste';
  const source = seg([['genome', el('span', {}, 'the genome in hand: ', view.genomeName)], ['paste', 'a sequence of your own']], kept.source, (v) => { kept.source = v; view.pasteBox.hidden = v !== 'paste'; if (v === 'genome') runStrands(); });
  const k = seg([1, 2, 3, 4, 5, 6].map((x) => [x, `k = ${x}`]), kept.k, (v) => { kept.k = v; view.ctx.setParams({ k: v }); runStrands(); });
  view.strandBtn = btn('Measure', { icon: 'play', kind: 'primary', size: 'sm', onClick: () => runStrands() });
  return el('section', { id: 'gr-strands', class: 'gr-section' },
    sectionHead(3, 'The strand symmetry of a genome', 'which of the eight symmetries of words its k-mer statistics keep (§4, §7)'),
    card({ body: [el('div', { class: 'ctl-strip' }, fld('Of', source), fld('Words of length', k), el('span', { class: 'grow' }), view.strandBtn), view.pasteBox, el('div', { class: 'card-body' }, view.strandStage)] }));
}

async function runStrands() {
  const t = (tickets.strands = (tickets.strands || 0) + 1);
  const params = { k: kept.k };
  if (kept.source === 'paste') {
    if (!kept.pasted.trim()) { kept.sym = null; kept.symError = 'Paste a sequence first, or choose the genome in hand.'; drawStrands(); return; }
    params.text = kept.pasted;
  } else {
    if (!view.ctx.genome) { kept.symError = 'No genome is in hand.'; drawStrands(); return; }
    params.genome = view.ctx.genome;
  }
  const done = kept.sym ? busy(view.strandStage, 'counting k-mers and their images') : null;
  if (!done) view.strandStage.replaceChildren(working('counting k-mers and their images'));
  view.strandBtn.disabled = true;
  try {
    const r = await api.run('grammar.symmetry', params);
    if (t !== tickets.strands) return;
    kept.sym = r; kept.symError = null;
  } catch (e) {
    if (t !== tickets.strands) return;
    kept.symError = e.message;
  } finally {
    if (t === tickets.strands) { done?.(); view.strandBtn.disabled = false; }
  }
  drawStrands();
}

function drawStrands() {
  const box = clear(view.strandStage);
  if (kept.symError) box.append(callout(kept.symError, 'bad'));
  const r = kept.sym;
  if (!r) { if (!kept.symError) box.append(working('counting k-mers and their images')); return; }
  box.append(...strandFigures(r, { hint: true }));
}

// The answer of grammar.symmetry as two figures and their reading, for a two-column grid.
export function strandFigures(r, { hint = false } = {}) {
  const items = [];
  const box = { append: (...xs) => items.push(...xs.filter(Boolean)) };
  const rc = r.elements.find((e) => e.short === 'κρ');
  const title = String(r.name).replace(/, complete genome.*$/, '').replace(/\.$/, '');
  const ABOUT8 = { '1': 'the identity', 'ρ': 'reversal', 'τ': 'transition', 'κ': 'complement', 'σ': 'A↔C, G↔T', 'τρ': 'transition, reversed', 'κρ': 'reverse complement', 'σρ': 'σ, reversed' };
  const order = [...r.elements].sort((a, b) => a.total_variation - b.total_variation);
  const tv = el('div', { class: 'gr-figure' });
  tv.innerHTML = G.hbarsSvg(order.map((e) => ({ label: e.short, sub: ABOUT8[e.short] || '', value: e.total_variation,
    note: e.total_variation.toFixed(4), colour: e.short === 'κρ' ? 'var(--accent)' : e.short === '1' ? 'var(--good)' : 'var(--ink-faint)', title: `${e.name}: TV(f, g·f) = ${e.total_variation.toFixed(5)}` })),
  { width: 560, row: 34, label: 130, max: Math.max(...r.elements.map((e) => e.total_variation)) * 1.05 });
  const name = (c) => `κ${c.kappa > 0 ? '+' : '−'} τ${c.tau > 0 ? '+' : '−'} ρ${c.rho > 0 ? '+' : '−'}`;
  const chs = [...r.channels].sort((a, b) => b.fraction - a.fraction);
  const ch = el('div', { class: 'gr-figure' });
  ch.innerHTML = G.hbarsSvg(chs.map((c) => ({ label: name(c), sub: c.rc > 0 ? 'even under κρ' : 'odd under κρ', value: c.fraction, note: sci(c.fraction),
    colour: c.rc > 0 ? 'var(--accent)' : 'var(--C)', title: `the channel with χ(κ) = ${c.kappa}, χ(τ) = ${c.tau}, χ(ρ) = ${c.rho}: ${c.fraction.toExponential(3)} of ‖f‖²` })),
  { width: 560, row: 34, label: 130, log: true, axis: 'the share of ‖f‖² in the channel, on a log scale', format: sci });
  box.append(
    el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, `The eight symmetries, by how far they move the ${r.k}-mer statistics`), tv,
      el('p', { class: 'note' }, 'Each bar is the total variation TV(f, g·f) = ½ Σ |f(w) − f(g⁻¹·w)| over the words w between the frequencies of one strand\'s k-mers and their image: 0 exactly when g keeps the statistics.')),
    el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, 'The eight channels, the frequency vector\'s isotypic parts'), ch,
      el('p', { class: 'note' }, `The central idempotents of the group\'s algebra split f into eight orthogonal parts, one per character; their energies sum to ‖f‖² (${(Math.abs(r.channel_sum - r.norm) / r.norm).toExponential(0)} apart, in double precision). Odd under κρ: ${sci(r.odd_rc)}; odd under κ: ${sci(r.odd_complement)}, ${r.odd_rc > 0 ? Math.round(r.odd_complement / r.odd_rc).toLocaleString('en-US') : '∞'} times more.`)));
  // Chargaff's second rule is an empirical fact of long genomes: say it only when this sequence shows it.
  const others = Math.min(...r.elements.filter((e) => e.short !== '1' && e.short !== 'κρ').map((e) => e.total_variation));
  const chargaff = rc && rc.total_variation < 0.05 && others > 4 * rc.total_variation;
  const reading = chargaff
    ? `Of the eight elements, only the identity and the reverse complement κρ nearly keep the statistics of ${title}: TV ${rc.total_variation.toFixed(4)} for κρ, ${others.toFixed(3)} at the least for any other. This is Chargaff's second parity rule, found as a group: the empirical symmetry group of the statistics is {1, κρ} ≅ C₂. Complement alone and reversal alone, which the draft lists as separate factors, are not symmetries of the statistics; only their product is.`
    : `Here the reverse complement is not singled out (TV ${rc ? rc.total_variation.toFixed(4) : '–'} for κρ, ${others.toFixed(3)} at the least for any other). Chargaff's second parity rule, that one strand's k-mer counts nearly equal their reverse complements', is a fact of long genomes; a short sequence, or long words over a small one, need not show it.`;
  box.append(el('div', { class: 'gr-span' }, callout(el('span', {}, reading), '', 'info'),
    el('p', { class: 'note' }, `${fmt(r.bases)} bases${r.records > 1 ? ` in ${r.records} records` : ''}; of the ${fmt(r.kmers)} words of length ${r.k}, ${r.canonical.reverse_complement} are canonical under κρ, ${r.canonical.strand} under the strand group {1, κ, ρ, κρ} and ${r.canonical.all_eight} under all eight, by Burnside's lemma with exact fixed-word counts.`),
    hint ? toolboxHint('The same for a file of your own: ', 'grammar-symmetry', 'Strand symmetry of a sequence') : null));
  return items;
}

// ---- 4. the genetic code ---------------------------------------------------------------------------------------
function codeSection() {
  view.codeStage = el('div', { class: 'gr-code' });
  const codes = select(CODES, kept.code, (v) => { kept.code = +v; kept.element = null; view.ctx.setParams({ code: v }); runCode(); });
  view.codeSelect = codes;
  const quick = el('div', { class: 'chips' }, ...[[1, 'the standard code'], [2, 'vertebrate mitochondria']].map(([v, l]) => el('button', { class: 'chip', type: 'button', onclick: () => { kept.code = v; kept.element = null; codes.value = v; view.ctx.setParams({ code: v }); runCode(); } }, l)));
  return el('section', { id: 'gr-code', class: 'gr-section' },
    sectionHead(4, 'The genetic code on (C₂)⁶', 'codons as a torsor of V₄³, a code as a function on it, and its symmetries (§6)'),
    card({ body: [el('div', { class: 'ctl-strip' }, fld('Genetic code', codes), fld('Compare', quick), el('span', { class: 'grow' }),
      el('span', { class: 'small faint', style: { alignSelf: 'center' } }, 'Click an element of the left grid to see where it fails.')), el('div', { class: 'card-body' }, view.codeStage)] }));
}

async function runCode() {
  const t = (tickets.code = (tickets.code || 0) + 1);
  const done = kept.codeRes ? busy(view.codeStage, 'the code on (C₂)⁶') : null;
  if (!done) view.codeStage.replaceChildren(working('the code on (C₂)⁶'));
  try {
    const params = { table: kept.code };
    if (kept.element !== null) params.element = kept.element;
    const r = await api.run('grammar.code', params);
    if (t !== tickets.code) return;
    kept.codeRes = r;
  } catch (e) {
    if (t !== tickets.code) return;
    kept.codeRes = { error: e.message };
  } finally {
    if (t === tickets.code) done?.();
  }
  drawCode();
}

function drawCode() {
  const box = clear(view.codeStage);
  const r = kept.codeRes;
  if (!r) { box.append(working('the code on (C₂)⁶')); return; }
  if (r.error) { box.append(callout(r.error, 'bad')); return; }
  box.append(...codeFigures(r, (bits) => { kept.element = bits; runCode(); }));
}

// The answer of grammar.code as its two grids, their reading, and the facts beside them; a click on
// an element of the profile calls onPick(its bits).
export function codeFigures(r, onPick = null) {
  const items = [];
  const box = { append: (...xs) => items.push(...xs.filter(Boolean)) };
  const stab = new Set(r.stabiliser.map((s) => s.bits));
  const ramp = (t) => `color-mix(in srgb, var(--accent) ${Math.round(8 + 80 * t)}%, var(--surface))`;
  const profile = el('div', { class: 'gr-figure gr-grid' });
  profile.innerHTML = G.codeGridSvg(r.elements, { names: ['1', 'σ', 'τ', 'κ'], selected: r.element.bits,
    colour: (c) => (c.keeps === 64 ? 'var(--accent)' : ramp(Math.max(0, (c.keeps - 16) / 48) ** 2)), text: (c) => c.keeps, mark: (c) => stab.has(c.bits),
    title: (c) => `${c.name}: keeps the amino acid of ${c.keeps} codons of 64${c.keeps === 64 ? ', a symmetry of the code' : ''}` });
  profile.addEventListener('click', (e) => {
    const g = e.target.closest('[data-bits]');
    if (g && onPick) onPick(+g.dataset.bits);
  });
  const maxP = Math.max(...r.characters.slice(1).map((c) => c.power));
  const spectrum = el('div', { class: 'gr-figure gr-grid' });
  spectrum.innerHTML = G.codeGridSvg(r.characters, { names: ['·', 'RY', 'MK', 'WS'],
    colour: (c) => (c.bits === 0 ? 'var(--paper-4)' : c.power === 0 ? 'var(--surface)' : ramp(Math.log1p(c.power) / Math.log1p(maxP))),
    text: () => '', mark: (c) => c.bits !== 0 && c.power === 0 && (r.element.keeps === 64) && (popcount(c.bits & r.element.bits) % 2 === 1),
    title: (c) => `${c.name}: P(u) = ${fmt(c.power)}${c.bits === 0 ? ' (the constant part)' : ''}` });
  const els = r.stabiliser.map((s) => el('span', { class: 'chip on' }, s.name));
  const side = el('div', { class: 'gr-col' },
    el('div', { class: 'eyebrow' }, `Code ${r.table}: ${r.name}`),
    fld('Its stabiliser, the elements that keep every amino acid', el('div', { class: 'chips' }, ...els)),
    el('p', { class: 'note' }, r.stabiliser.length === 1 ? 'Only the identity keeps the code: its stabiliser is trivial.' : `{${r.stabiliser.map((s) => s.name).join(', ')}}: the wobble of the third position is invariance under a subgroup.`),
    fld(`The element ${r.element.name}`, el('div', { class: 'stack' }, el('span', { class: 'small' }, `keeps ${r.element.keeps} codons of 64`),
      r.element.failures.length ? el('div', { class: 'chips' }, ...r.element.failures.slice(0, 24).map((f) => el('span', { class: 'chip mono' }, f)), r.element.failures.length > 24 ? el('span', { class: 'small faint' }, `and ${r.element.failures.length - 24} more pairs`) : null)
        : el('span', { class: 'badge good' }, 'it fails nowhere: a symmetry'))),
    metrics([[pct(r.energy.first_two), 'of the non-constant spectrum in characters of the first two bases'], [pct(r.energy.first), 'in the first base alone'], [pct(r.energy.third), 'in the third alone']]),
    el('p', { class: 'note' }, 'P = Ŵ I, the Walsh–Hadamard transform of the profile, exactly in integers ', yes(r.wiener_khinchin, 'holds', 'fails'), '; m keeps the code exactly when P vanishes on every character u with u·m odd ', yes(r.duality, 'holds', 'fails'), ' (Proposition 6.2).'));
  box.append(
    el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, 'The invariance profile I(m)'), profile,
      el('p', { class: 'note' }, 'For each of the 64 elements m, the codons c with aa(c ⊕ m) = aa(c). The outlined cells are the stabiliser.')),
    el('div', { class: 'gr-col' }, el('div', { class: 'eyebrow' }, 'The power spectrum P(u)'), spectrum,
      el('p', { class: 'note' }, 'P(u) = Σ f̂ₐ(u)² over the amino acids a, the amino acids\' indicators transformed; each position\'s character is one of the classifications RY, MK, WS, or trivial (·). Darker holds more.')),
    side);
  if (r.table === 1 || r.table === 2) {
    box.append(el('div', { class: 'gr-span' }, callout(r.table === 2 ? 'The vertebrate mitochondrial code reads ATA as Met and TGA as Trp. Its stabiliser is {1, τ₃}, and its spectrum vanishes on every character odd under τ₃ (Propositions 6.1 and 6.2). The two codes are two functions on one group, and this one keeps one symmetry more.'
      : 'The standard code\'s stabiliser is trivial. The third-position transition τ₃ keeps 60 codons of 64 and fails at exactly two pairs, ATA/ATG (Ile/Met) and TGA/TGG (stop/Trp) (Proposition 6.1). Most of the code\'s structure lies in characters that ignore the third base: the wobble, in Fourier terms.', '', 'info')));
  }
  return items;
}

const popcount = (x) => { let n = 0; while (x) { n += x & 1; x >>= 1; } return n; };

// ---- 5. characters -------------------------------------------------------------------------------------------
function charactersSection() {
  view.charStage = el('div', { class: 'gr-two' });
  const opts = GROUPS.map(([v, l]) => [v, l]);
  const a = select(opts.some(([v]) => v === kept.group) ? opts : [[kept.group, kept.group], ...opts], kept.group, (v) => { kept.group = v; view.ctx.setParams({ group: v }); runCharacters(); });
  const b = select([['', 'nothing'], ...opts], kept.compare, (v) => { kept.compare = v; view.ctx.setParams({ vs: v }); runCharacters(); });
  const own = el('input', { type: 'text', placeholder: 'or write one: dihedral 6, cyclic 2 × symmetric 3 …', style: { width: '300px' } });
  own.addEventListener('keydown', (e) => { if (e.key === 'Enter' && own.value.trim()) { kept.group = own.value.trim(); view.ctx.setParams({ group: kept.group }); runCharacters(); } });
  const pair = (x, y, label) => el('button', { class: 'chip', type: 'button', onclick: () => { kept.group = x; kept.compare = y; a.value = x; b.value = y; view.ctx.setParams({ group: x, vs: y }); runCharacters(); } }, label);
  return el('section', { id: 'gr-characters', class: 'gr-section' },
    sectionHead(5, 'The characters of a finite group', 'exactly, by Dixon and Schneider, with the blocks of ℂ[G] ≅ ⊕ Mχ₍₁₎(ℂ) (§7)'),
    card({ body: [el('div', { class: 'ctl-strip' }, fld('A group', a), fld('Beside', b), fld('Your own', own), fld('Two that share their degrees', el('div', { class: 'chips' }, pair('cyclic 4', 'klein', 'C₄ and C₂ × C₂'), pair('dihedral 4', 'quaternion', 'D₄ and Q₈')))),
      el('div', { class: 'card-body' }, view.charStage)] }));
}

async function runCharacters() {
  const t = (tickets.chars = (tickets.chars || 0) + 1);
  const want = [kept.group, kept.compare].filter(Boolean);
  const done = Object.keys(kept.chars).length ? busy(view.charStage, 'character tables') : null;
  if (!done) view.charStage.replaceChildren(working('character tables'));
  const out = {};
  await Promise.all(want.map(async (g) => {
    try { out[g] = kept.chars[g] && !kept.chars[g].error ? kept.chars[g] : await api.run('grammar.characters', { group: g }); } catch (e) { out[g] = { error: e.message, name: g }; }
  }));
  if (t !== tickets.chars) return;
  done?.();
  kept.chars = out;
  drawCharacters();
}

// The real blocks from the Frobenius–Schur indicators: +1 gives M_d(ℝ), −1 gives M_{d/2}(ℍ), and a
// pair of conjugate characters (indicator 0) one M_d(ℂ).
export function realBlocks(r) {
  const parts = [];
  let complexPairs = 0;
  for (const c of r.characters) {
    const ind = Math.round(c.indicator);
    if (ind === 1) parts.push(c.degree === 1 ? 'ℝ' : `M${G.sub(c.degree)}(ℝ)`);
    else if (ind === -1) parts.push(c.degree === 2 ? 'ℍ' : `M${G.sub(c.degree / 2)}(ℍ)`);
    else complexPairs++;
  }
  const cx = r.characters.filter((c) => Math.round(c.indicator) === 0);
  for (let k = 0; k < complexPairs / 2; k++) parts.push(cx[2 * k].degree === 1 ? 'ℂ' : `M${G.sub(cx[2 * k].degree)}(ℂ)`);
  return collapse(parts);
}
const collapse = (parts) => {
  const counts = new Map();
  for (const p of parts) counts.set(p, (counts.get(p) || 0) + 1);
  return [...counts].map(([p, n]) => (n > 1 && !p.startsWith('M') ? `${p}${supNum(n)}` : Array(n).fill(p).join(' ⊕ '))).join(' ⊕ ');
};
// ℚ[G] for an abelian group: each Galois orbit of characters of order d is the field ℚ(ζ_d).
function rationalBlocks(r) {
  const w = r.wedderburn;
  if (!w?.rational) return null;
  if (!r.abelian) return null;
  const parts = w.rational.map((b) => {
    const c = r.characters[b.characters[0]];
    let d = 1;
    for (const [re, im] of c.complex) {
      const turns = Math.atan2(im, re) / (2 * Math.PI);
      for (let q = 1; q <= r.exponent; q++) if (Math.abs(turns * q - Math.round(turns * q)) < 1e-6) { d = lcm(d, q); break; }
    }
    return d <= 2 ? 'ℚ' : d === 4 ? 'ℚ(i)' : `ℚ(ζ${G.sub(d)})`;
  });
  return collapse(parts);
}
const gcd = (a, b) => (b ? gcd(b, a % b) : a);
const lcm = (a, b) => (a * b) / gcd(a, b);

// One group's character table, its degrees and its blocks (the view's and the toolbox's).
export function characterCard(r) {
  if (r.error) return callout(`${r.name}: ${r.error}`, 'bad');
  const degrees = r.characters.map((c) => c.degree);
  const t = el('table', { class: 'data gr-table' },
    el('tr', {}, el('th', {}, ''), ...r.classes.map((c) => el('th', { class: 'num', title: `${c.size} element${c.size > 1 ? 's' : ''} of order ${c.order}` }, c.label, el('div', { class: 'faint small' }, `${c.size}`))), el('th', { class: 'num', title: 'the Frobenius–Schur indicator' }, 'ι')),
    ...r.characters.map((c, k) => el('tr', {}, el('td', {}, 'χ', el('sub', {}, String(k + 1))),
      ...c.values.map((v, j) => el('td', { class: `num mono ${v === '0' ? 'faint' : ''}`, title: `${v} ≈ ${c.complex[j][0].toFixed(4)}${c.complex[j][1] ? ` ${c.complex[j][1] < 0 ? '−' : '+'} ${Math.abs(c.complex[j][1]).toFixed(4)}i` : ''}` }, G.cyclotomic(v))),
      el('td', { class: 'num' }, el('span', { class: `badge ${c.indicator > 0.5 ? 'good' : c.indicator < -0.5 ? 'accent' : ''}` }, c.indicator > 0.5 ? '+1' : c.indicator < -0.5 ? '−1' : '0')))));
  const blocks = el('div', { class: 'gr-figure' });
  blocks.innerHTML = G.wedderburnSvg(degrees, { unit: Math.max(6, Math.min(22, Math.floor(220 / degrees.reduce((s, d) => s + d, 0)))) });
  const w = r.wedderburn;
  const complexBlocks = collapse(degrees.map((d) => (d === 1 ? 'ℂ' : `M${G.sub(d)}(ℂ)`)));
  const rat = rationalBlocks(r);
  return el('div', { class: 'gr-col' },
    el('div', { class: 'row', style: { alignItems: 'baseline' } }, el('h3', { class: 'gr-group-name' }, r.name), el('span', { class: 'small faint' }, `order ${r.order}, ${r.classes.length} classes, exponent ${r.exponent}${r.abelian ? ', abelian' : ''}`)),
    el('div', { class: 'gr-scroll' }, t),
    el('p', { class: 'note' }, `Degrees ${degrees.join(', ')}; Σ χ(1)² = ${fmt(r.sum_of_squares)} = |G| `, yes(r.sum_of_squares === r.order, 'holds', 'fails'), '; the rows orthonormal ', yes(r.orthonormal, 'yes', 'no'),
      r.symmetric ? el('span', {}, '; the degrees those of the hook length formula on the partitions of ', String(r.symmetric.n), ' ', yes(r.symmetric.degrees_agree, 'agree', 'differ')) : '', '.'),
    el('div', { class: 'gr-blocks' }, blocks, el('div', { class: 'stack', style: { gap: '4px' } },
      el('div', {}, `ℂ[${r.name}] ≅ ${complexBlocks}`), el('div', {}, `ℝ[${r.name}] ≅ ${realBlocks(r)}`),
      rat ? el('div', {}, `ℚ[${r.name}] ≅ ${rat}`) : w.rational ? el('div', {}, `ℚ[${r.name}]: blocks of dimension ${w.rational.map((x) => x.dimension).join(' + ')}, one for each Galois orbit of characters`) : null,
      el('div', { class: 'small faint' }, w.measured ? `each block's dimension measured as the rank of x ↦ x·eᵪ, exactly: ${w.blocks.map((b) => b.dimension).join(' + ')} = ${w.total}` : `|G| > 64: the dimensions d² are the theorem's, not measured`))));
}

function drawCharacters() {
  const box = clear(view.charStage);
  const a = kept.chars[kept.group], b = kept.compare ? kept.chars[kept.compare] : null;
  if (!a) { box.append(working('character tables')); return; }
  box.append(characterCard(a));
  if (b) box.append(characterCard(b));
  if (a && b && !a.error && !b.error) {
    const da = a.characters.map((c) => c.degree).sort((x, y) => x - y).join(','), db = b.characters.map((c) => c.degree).sort((x, y) => x - y).join(',');
    const ia = a.characters.map((c) => Math.round(c.indicator)).sort().join(','), ib = b.characters.map((c) => Math.round(c.indicator)).sort().join(',');
    const same = da === db;
    const text = same
      ? `${a.name} and ${b.name} have the same degrees {${da.replace(/,/g, ', ')}}, so ℂ[${a.name}] ≅ ℂ[${b.name}]: as an algebra, ℂ[G] determines only the multiset of degrees${a.order === b.order && (a.abelian !== b.abelian || ia !== ib || a.exponent !== b.exponent) ? ', and these groups differ' : ''}. ${ia !== ib ? `The indicators differ, so over ℝ the algebras do: ${realBlocks(a)} against ${realBlocks(b)}.` : ''} ${rationalBlocks(a) && rationalBlocks(b) && rationalBlocks(a) !== rationalBlocks(b) && a.abelian && b.abelian ? `Over ℚ, which exact arithmetic sees, they differ too: ${rationalBlocks(a)} against ${rationalBlocks(b)}.` : ''} The group is recovered as the group-like elements of the Hopf algebra ℂ[G], which needs the coproduct Δ(g) = g ⊗ g (§7; the draft's §4.2).`
      : `The degrees differ ({${da.replace(/,/g, ', ')}} against {${db.replace(/,/g, ', ')}}), so the two group algebras are not isomorphic.`;
    box.append(el('div', { class: 'gr-span' }, callout(text, '', 'info')));
  }
}

// ---- the theory note ------------------------------------------------------------------------------------------
async function readNote() {
  const body = el('div', { class: 'gr-note' }, working('the note'));
  dialog({ title: 'Group grammars for genomes', body, width: 900, cls: 'gr-note-dialog' });
  try {
    const r = await fetch(NOTE);
    if (!r.ok) throw new Error(`${r.status} ${r.statusText}`);
    body.innerHTML = markdown(await r.text());
  } catch (e) {
    body.replaceChildren(callout(`The note could not be read: ${e.message}`, 'bad'));
  }
}

// Markdown as the note is written: headings, paragraphs, nested lists, tables, emphasis, code and
// links; everything escaped first.
function markdown(text) {
  const inline = (s) => esc(s)
    .replace(/`([^`]+)`/g, '<code>$1</code>')
    .replace(/\*\*([^*]+)\*\*/g, '<b>$1</b>')
    .replace(/(^|[\s(])\*([^*\s][^*]*)\*/g, '$1<i>$2</i>')
    .replace(/\[([^\]]+)\]\(([^)]+)\)/g, (m, t, u) => `<a class="link" href="${/^https?:/.test(u) ? u : '#'}" target="_blank" rel="noopener">${t}</a>`);
  const lines = text.replace(/\r/g, '').split('\n');
  const out = [];
  let para = [], list = [];
  const flushPara = () => { if (para.length) out.push(`<p>${inline(para.join(' '))}</p>`); para = []; };
  const flushList = () => {
    if (!list.length) return;
    const html = [], stack = [];
    for (const { depth, ordered, text: t } of list) {
      while (stack.length > depth + 1) html.push(`</li></${stack.pop()}>`);
      if (stack.length === depth + 1) html.push('</li>');
      else while (stack.length < depth + 1) { const tag = ordered ? 'ol' : 'ul'; stack.push(tag); html.push(`<${tag}>`); }
      html.push(`<li>${inline(t)}`);
    }
    while (stack.length) html.push(`</li></${stack.pop()}>`);
    out.push(html.join(''));
    list = [];
  };
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    const h = line.match(/^(#{1,4})\s+(.*)$/);
    const li = line.match(/^(\s*)([-*]|\d+\.)\s+(.*)$/);
    if (h) { flushPara(); flushList(); out.push(`<h${h[1].length + 1}>${inline(h[2])}</h${h[1].length + 1}>`); }
    else if (/^\s*\|/.test(line)) {
      flushPara(); flushList();
      const rows = [];
      while (i < lines.length && /^\s*\|/.test(lines[i])) rows.push(lines[i++]);
      i--;
      const cells = (r) => r.trim().replace(/^\||\|$/g, '').split('|').map((c) => c.trim());
      const body = rows.filter((r) => !/^\s*\|[\s|:-]+\|\s*$/.test(r));
      out.push(`<table class="data">${body.map((r, k) => `<tr>${cells(r).map((c) => (k ? `<td>${inline(c)}</td>` : `<th>${inline(c)}</th>`)).join('')}</tr>`).join('')}</table>`);
    } else if (li) { flushPara(); list.push({ depth: Math.floor(li[1].length / 2), ordered: /\d/.test(li[2]), text: li[3] }); }
    else if (!line.trim()) { flushPara(); flushList(); }
    else if (list.length && /^\s+/.test(line)) list[list.length - 1].text += ` ${line.trim()}`;
    else { flushList(); para.push(line.trim()); }
  }
  flushPara(); flushList();
  return out.join('\n');
}
