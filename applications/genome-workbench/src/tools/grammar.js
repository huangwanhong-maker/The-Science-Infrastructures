// Tools of group grammars (domains/genomics/theory/group-grammars.md): the hairpins of a sequence by
// the inside algorithm over all its spans, the strand symmetry of a sequence's k-mers, the genetic
// code on (C₂)⁶, a grammar parsed up to a group, and the characters of a finite group. Their figures
// are the Grammar view's (views/grammar.js, draw/grammar.js).

import { el } from '../dom.js';
import { tool } from '../core/registry.js';
import * as api from '../api.js';
import { card, callout, metrics } from '../ui/kit.js';
import { tableCard } from '../ui/results.js';
import * as G from '../draw/grammar.js';
import { parseCards, strandFigures, codeFigures, characterCard } from '../views/grammar.js';
import { source, CODES, SAMPLES, n, pct, note } from './common.js';

// What the engine read of the input, when it stopped early or found lines it could not read.
const reading = (r, out) => {
  const g = r.reading;
  if (!g) return;
  if (g.truncated && g.stopped) out.append(note(`Of ${g.source}: ${g.stopped}.`));
  if (g.problem_count) out.append(note(`${n(g.problem_count)} ${g.problem_count === 1 ? 'line' : 'lines'} could not be read, the first at line ${g.problems[0].line}: ${g.problems[0].message}.`));
};

// Samples: transfer RNA genes of the human mitochondrion (RefSeq NC_012920.1), read from the record
// when asked, by their gene names; never typed in.
const MITO = 'ncbi:NC_012920.1';
async function trna(gene) {
  const g = await api.run('genome.gene', { genome: MITO, gene });
  return `>human_mt_${gene}\n${g.sequence}\n`;
}
async function stretch(genes) {
  const found = await Promise.all(genes.map((gene) => api.run('genome.gene', { genome: MITO, gene })));
  const begin = Math.min(...found.flatMap((g) => g.parts.map((p) => p[0]))), end = Math.max(...found.flatMap((g) => g.parts.map((p) => p[1])));
  const t = await api.run('genome.tracks', { genome: MITO, begin, end, bins: 4, letters_max: end - begin + 10 });
  return `>human_mt_${begin + 1}_${end}\n${t.letters}\n`;
}
const TRNAS = () => stretch(['TRNI', 'TRNQ', 'TRNM']);

// The 64 elements of V₄³ by name, as the engine names them: per position σ, τ or κ with the position.
const BY_BITS = ['', 'σ', 'τ', 'κ'];
const elementName = (m) => [0, 1, 2].map((p) => { const b = (m >> (2 * (2 - p))) & 3; return b ? `${BY_BITS[b]}${G.sub(p + 1)}` : ''; }).join('') || '1';
const ELEMENTS = Array.from({ length: 64 }, (_, m) => [m, elementName(m)]).sort((a, b) => a[1].length - b[1].length || a[0] - b[0]);

const GROUP_OPTIONS = [
  ['cyclic 4', 'C₄'], ['klein', 'C₂ × C₂, the Klein group'], ['cyclic 6', 'C₆'], ['symmetric 3', 'S₃'], ['dihedral 4', 'D₄, of order 8'], ['quaternion', 'Q₈'],
  ['alternating 4', 'A₄'], ['symmetric 4', 'S₄'], ['dihedral 5', 'D₅'], ['cyclic 3 × cyclic 3', 'C₃ × C₃'], ['alternating 5', 'A₅'], ['symmetric 5', 'S₅'],
];

tool({
  id: 'grammar-hairpins', title: 'Hairpins by grammar', group: 'grammar', icon: 'annotation', keywords: ['hairpin', 'stem', 'loop', 'rna', 'trna', 'secondary structure', 'inside', 'cyk', 'palindrome', 'wobble'],
  summary: 'Every hairpin of a stretch of up to 400 bases: the spans a hairpin grammar derives, from one inside pass over all of them, with the grammar\'s symmetry and graduated breaking over R/Y.',
  refs: 'Lari and Young (1990), Computer Speech and Language 4, 35–56 (inside–outside); Charniak et al. (2006), NAACL (coarse-to-fine parsing); the theory note, §2.5 and §5.',
  about: 'A hairpin is a stem of paired bases closing a loop: u v u′ with |u| = |u′|, the outermost bases of u and u′ paired, then the next, and so on inwards, every pair allowed. The grammar emits stem and loop from the outside in, and its inside values over every span say which spans are hairpins; the definition, checked by brute force on the same spans, must agree. With wobble pairs the grammar keeps reversal and loses the complement.',
  inputs: [
    { key: 'sequences', type: 'sequences', label: 'Sequence', required: true, kinds: ['fasta', 'genbank'], samples: [['three tRNA genes of the human mitochondrion', TRNAS], ['its tRNA-Phe', () => trna('TRNF')]] },
    { key: 'min_stem', type: 'number', label: 'A stem of at least', min: 2, max: 12, step: 1, unit: 'pairs', default: 4 },
    { key: 'max_span', type: 'number', label: 'Hairpins up to', min: 10, max: 200, step: 10, unit: 'bases', default: 80 },
    { key: 'min_loop', type: 'number', label: 'A loop from', min: 1, max: 12, step: 1, unit: 'bases', default: 3 },
    { key: 'max_loop', type: 'number', label: 'to', min: 1, max: 20, step: 1, unit: 'bases', default: 8 },
    { key: 'wobble', type: 'switch', label: 'Pairs', text: 'G·T wobble pairs too (G·U in RNA)', default: false },
    { key: 'from', type: 'number', label: 'From base', min: 1, max: 1e9, step: 1, default: 1 },
    { key: 'to', type: 'number', label: 'to base', min: 0, max: 1e9, step: 1, default: 0, hint: '0: to the end; a stretch of at most 400 bases' },
  ],
  example: async () => ({ sequences: { text: await TRNAS() }, min_stem: 4, min_loop: 3, max_loop: 8, wobble: true }),
  run: (v) => api.run('grammar.hairpins', { ...source(v.sequences), min_stem: v.min_stem, min_loop: v.min_loop, max_loop: Math.max(v.min_loop, v.max_loop), wobble: v.wobble, max_span: v.max_span,
    begin: Math.max(0, v.from - 1), ...(v.to > 0 ? { end: v.to } : {}) }),
  render(r, v, out) {
    const longest = r.maximal[r.longest ?? 0];
    out.append(metrics([[n(r.maximal_count), 'maximal hairpins, inside no other'], [n(r.spans), 'spans the grammar derives'], [r.agrees_with_definition ? 'agree' : 'differ', 'with the definition, by brute force on every span', r.agrees_with_definition ? 'good' : ''],
      [longest ? `${longest.stem} + ${longest.loop}` : '–', 'pairs and loop of the longest'], [`${n(r.inside_seconds * 1000, 1)} ms`, 'the inside pass']]));
    reading(r, out);
    if (!r.maximal.length) { out.append(callout(`No hairpin of ${r.spec.min_stem} pairs or more closing a loop of ${r.spec.min_loop} to ${r.spec.max_loop} in ${r.name} ${r.begin + 1}..${r.end}.`)); }
    else {
      let chosen = r.longest ?? 0;
      const arcs = el('div', { class: 'gr-figure gr-scroll' });
      const loop = el('div', { class: 'gr-figure' });
      const draw = () => {
        arcs.innerHTML = G.arcSvg(r.letters, r.maximal, { begin: r.begin, selected: chosen, width: 980 });
        const h = r.maximal[chosen];
        loop.innerHTML = G.stemLoopSvg(h);
        caption.textContent = `${h.begin + 1}..${h.end}: ${h.stem} pairs closing a loop of ${h.loop}; ${h.structure}`;
      };
      const caption = el('p', { class: 'note mono small' });
      arcs.addEventListener('click', (e) => { const g = e.target.closest('[data-hairpin]'); if (g) { chosen = +g.dataset.hairpin; draw(); } });
      draw();
      out.append(card({ title: 'The hairpins over their sequence', sub: `each pair an arc; ${r.name}, bases ${r.begin + 1}..${r.end}; click one to draw it`, body: el('div', { class: 'card-body stack' }, arcs) }));
      out.append(el('div', { class: 'gr-two' },
        card({ title: 'Stem and loop', sub: 'Watson–Crick pairs solid, wobble dotted', body: el('div', { class: 'card-body stack' }, loop, caption) }),
        r.tree ? card({ title: 'The Viterbi parse of the longest', sub: `stem rules weighing 1: the parse with the most pairs, ${r.tree_pairs}`, body: el('div', { class: 'card-body gr-scroll' }, (() => { const f = el('div', { class: 'gr-figure' }); f.innerHTML = G.treeSvg(r.tree, { slot: 18, level: 30, maxWidth: 640 }).svg; return f; })()) }) : null));
      out.append(tableCard('The maximal hairpins', [{ key: 'begin', label: 'from', num: true, format: (x) => n(x + 1) }, { key: 'end', label: 'to', num: true, format: (x) => n(x) }, { key: 'stem', label: 'pairs', num: true }, { key: 'loop', label: 'loop', num: true },
        { key: 'structure', label: 'structure', format: (x) => el('span', { class: 'mono small' }, x) }, { key: 'letters', label: 'bases', format: (x) => el('span', { class: 'mono small' }, x) }], r.maximal, { file: 'hairpins.csv' }));
    }
    const gm = r.grammar;
    out.append(card({ title: 'The grammar\'s symmetry', sub: `${gm.rules} rules, ${gm.nonterminals} nonterminals${r.spec.wobble ? ', with wobble pairs' : ''}`, body: el('div', { class: 'card-body stack' },
      el('div', { class: 'chips' }, ...gm.by_element.map((e) => el('span', { class: `chip${e.keeps ? ' on' : ''}` }, e.keeps ? '✓ ' : '✗ ', e.name))),
      note(`Of the eight symmetries of words, ${gm.stabiliser.length === 8 ? 'all keep the rules: a Watson–Crick pair is (x, κx), and every element of V₄ commutes with κ' : `{${gm.stabiliser.join(', ')}} keep the rules`}; ${gm.strand_invariant ? 'the strand group {1, κ, ρ, κρ} among them' : 'not the complement κ, since κ(G·T) = C·A pairs in neither, while reversal survives'} (Proposition 5.1). Under them the rules fall into ${n(gm.orbits.count)} orbits, as Burnside's count says (${n(gm.orbits.burnside)}), against |R|/|G| = ${(gm.rules / gm.orbits.group_order).toFixed(1)}.`)) }));
    const c = r.coarse;
    if (c && !c.skipped) {
      const f = el('div', { class: 'gr-figure' });
      f.innerHTML = G.hbarsSvg(c.by.map((b) => ({ label: `over ${b.name}`, value: b.fraction, note: `${pct(b.fraction)} of ${n(b.spans)} spans kept${b.exact ? '' : ', not exact'}`, colour: 'var(--accent)' })), { width: 600, row: 34, label: 110, max: 1 });
      out.append(card({ title: 'Graduated breaking', sub: `the grammar of a hairpin anywhere on bases ${c.begin + 1}..${c.end}, parsed first over a coarser alphabet (§2.5)`, body: el('div', { class: 'card-body stack' }, f,
        note(`Each coarse pass parses the image grammar on the word's image, and keeps the spans that some nonterminal both derives and is derived around. The fine pass over those spans alone finds every one of the ${n(c.hairpins)} hairpin spans, exactly as the full pass does: ${c.exact ? 'checked' : 'not so here'}. The prune loses nothing because the coarse grammar is an image, not an approximation.`)) }));
    } else if (c?.skipped) out.append(note(`Graduated breaking is not shown: ${c.skipped}.`));
  },
});

tool({
  id: 'grammar-symmetry', title: 'Strand symmetry of a sequence', group: 'grammar', icon: 'statistics', keywords: ['chargaff', 'reverse complement', 'symmetry', 'kmer', 'strand', 'parity', 'group', 'channels', 'wedderburn'],
  summary: 'Which of the eight symmetries of words (complement, transition, reversal and their products) a sequence\'s k-mer statistics keep, and the energy in each of its eight Wedderburn channels.',
  refs: 'Rudner, Karkas and Chargaff (1968), PNAS 60, 921–922 (the second parity rule); Serre (1977), Linear Representations of Finite Groups, §2.6.',
  inputs: [
    { key: 'sequences', type: 'sequences', label: 'Sequences', required: true, kinds: ['fasta', 'genbank'], samples: [['phage lambda, whole', SAMPLES.lambdaWhole], ['the start of phiX174', SAMPLES.phix]] },
    { key: 'k', type: 'choice', label: 'Words of length', options: [['1', '1'], ['2', '2'], ['3', '3'], ['4', '4'], ['5', '5'], ['6', '6']], default: '4', block: false },
  ],
  example: async () => ({ sequences: { text: await SAMPLES.lambdaWhole() }, k: '4' }),
  run: (v) => api.run('grammar.symmetry', { ...source(v.sequences), k: Number(v.k) }),
  render(r, v, out) {
    const rc = r.elements.find((e) => e.short === 'κρ');
    out.append(metrics([[n(r.bases), 'bases'], [n(r.kmers), `words of length ${r.k}`], [rc ? rc.total_variation.toFixed(4) : '–', 'TV under the reverse complement'], [r.odd_rc.toExponential(1), 'of ‖f‖² odd under κρ']]));
    reading(r, out);
    out.append(el('div', { class: 'gr-two' }, ...strandFigures(r)));
  },
});

tool({
  id: 'grammar-code', title: 'Code symmetry', group: 'grammar', icon: 'sequence', keywords: ['genetic code', 'codon', 'wobble', 'walsh', 'hadamard', 'spectrum', 'stabiliser', 'mitochondrial', 'symmetry'],
  summary: 'A genetic code of NCBI\'s as a function on (C₂)⁶: which elements keep the amino acid of every codon, where any one of them fails, and the code\'s power spectrum.',
  refs: 'Terras (1999), Fourier Analysis on Finite Groups and Applications, ch. 2; the theory note, §6.',
  inputs: [
    { key: 'table', type: 'select', label: 'Genetic code', options: CODES, default: 1, numeric: true },
    { key: 'element', type: 'select', label: 'The element whose failures to list', options: ELEMENTS, default: 2, numeric: true, help: 'τ₃: the transition at the third position' },
  ],
  example: { table: 2, element: 2 },
  run: (v) => api.run('grammar.code', { table: v.table, element: v.element }),
  render(r, v, out) {
    out.append(metrics([[r.stabiliser.map((s) => s.name).join(', '), 'its stabiliser'], [`${r.element.keeps} of 64`, `codons ${r.element.name} keeps`], [pct(r.energy.first_two), 'of the spectrum in the first two bases']]));
    out.append(card({ title: `Code ${r.table}, ${r.name}`, sub: 'the profile, the spectrum and what they say', body: el('div', { class: 'card-body' }, el('div', { class: 'gr-code' }, ...codeFigures(r))) }));
    const top = [...r.elements].filter((e) => e.bits).sort((a, b) => b.keeps - a.keeps).slice(0, 16);
    out.append(tableCard('The elements that keep the most codons', [{ key: 'name', label: 'element' }, { key: 'keeps', label: 'codons kept, of 64', num: true }, { key: 'bits', label: 'as six bits', format: (x) => el('span', { class: 'mono' }, x.toString(2).padStart(6, '0')) }], top, { file: 'code-profile.csv' }));
  },
});

tool({
  id: 'grammar-parse', title: 'Parse with a grammar', group: 'grammar', icon: 'grammar', keywords: ['grammar', 'earley', 'parse', 'context-free', 'context-sensitive', 'monotone', 'burnside', 'orbit', 'equivariant', 'symmetry'],
  summary: 'A grammar and a group acting on its letters: invariance, rule orbits by Burnside\'s lemma, Earley\'s parse, the equivariant parser beside the orbit-only scan, and spans shared across orbits.',
  refs: 'Earley (1970), Comm. ACM 13, 94–102; Hopcroft and Ullman (1979); the theory note, §1 and §2.',
  about: 'Rules are separated by semicolons (S -> a S b | ε; …), symbols by spaces; a letter of the alphabet is a terminal, anything else a nonterminal. Each generator of the group is the images of the alphabet\'s letters in order, followed by ρ when it also reverses words; generators are separated by semicolons. The Grammar view has a fuller editor and examples.',
  inputs: [
    { key: 'alphabet', type: 'text', label: 'Alphabet', required: true, default: 'ab', half: true },
    { key: 'word', type: 'text', label: 'Word', default: 'ab', half: true },
    { key: 'rules', type: 'text', label: 'Rules', required: true, default: 'S -> a a | b b', help: 'separated by ;' },
    { key: 'group', type: 'text', label: 'Group', default: 'ba', help: 'generators separated by ;, ρ to reverse' },
    { key: 'use', type: 'choice', label: 'When the rules are not invariant, use', options: [['stabiliser', 'their stabiliser'], ['closure', 'their closure']], default: 'stabiliser' },
  ],
  example: { alphabet: 'ab', word: 'ab', rules: 'S -> a a | b b', group: 'ba', use: 'stabiliser' },
  run: (v) => api.run('grammar.parse', { alphabet: v.alphabet.trim(), word: v.word, rules: G.rulesFromText(v.rules), group: G.generatorsFromText(v.group), use: v.use }),
  render(r, v, out) {
    out.append(...parseCards(r));
  },
});

tool({
  id: 'grammar-characters', title: 'Characters of a group', group: 'grammar', icon: 'grammar', keywords: ['character', 'character table', 'representation', 'wedderburn', 'group algebra', 'frobenius', 'schur', 'dixon', 'cyclotomic'],
  summary: 'A finite group\'s exact character table (Dixon–Schneider, in cyclotomic numbers), its Frobenius–Schur indicators, and the Wedderburn blocks of its group algebra.',
  refs: 'Dixon (1967), Numer. Math. 10, 446–450; Schneider (1990), J. Symbolic Comput. 9, 601–606; Serre (1977), §2.6 and §6.2.',
  inputs: [
    { key: 'group', type: 'select', label: 'A group', options: GROUP_OPTIONS, default: 'quaternion' },
    { key: 'own', type: 'text', label: 'Or any other', default: '', placeholder: 'dihedral 6, cyclic 2 × symmetric 3 …', help: 'cyclic, klein, dihedral, quaternion, symmetric, alternating, and products with ×' },
  ],
  example: { group: 'dihedral 4', own: '' },
  run: (v) => api.run('grammar.characters', { group: v.own?.trim() || v.group }),
  render(r, v, out) {
    out.append(metrics([[n(r.order), '|G|'], [n(r.classes.length), 'classes, and characters'], [r.characters.map((c) => c.degree).join(', '), 'their degrees'], [n(r.sum_of_squares), 'Σ χ(1)², which is |G|', r.sum_of_squares === r.order ? 'good' : '']]));
    out.append(card({ title: `The characters of ${r.name}`, sub: 'exact, in cyclotomic numbers; hover a value for its complex number', body: el('div', { class: 'card-body' }, characterCard(r)) }));
  },
});

