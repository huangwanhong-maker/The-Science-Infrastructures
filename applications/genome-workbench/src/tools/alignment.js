// Tools of alignment: two sequences by any schedule, many sequences at once, a query searched
// against many with Karlin–Altschul E-values, and the selection on a gene between two species
// from its codon alignment.

import { el } from '../dom.js';
import { tool } from '../core/registry.js';
import * as api from '../api.js';
import { alignmentView, tableCard, textCard } from '../ui/results.js';
import { card, kv } from '../ui/kit.js';
import { source, CODES, MATRICES, SAMPLES, n, pct, tiles, note } from './common.js';
import { parseSequences, guessAlphabet } from '../ui/kit.js';

const scoringInputs = [
  { key: 'alphabet', type: 'choice', label: 'Letters', options: [['dna', 'DNA'], ['protein', 'protein']], default: 'dna' },
  { key: 'matrix', type: 'select', label: 'Scores', options: MATRICES.map((m) => [m, m]), default: 'BLOSUM62', when: (v) => v.alphabet === 'protein' },
  { key: 'match', type: 'number', label: 'Match', min: 1, max: 10, default: 2, when: (v) => v.alphabet === 'dna' },
  { key: 'mismatch', type: 'number', label: 'Mismatch', min: -20, max: -1, default: -4, when: (v) => v.alphabet === 'dna' },
  { key: 'gap_open', type: 'number', label: 'Gap open', min: 0, max: 50, default: 4 },
  { key: 'gap_extend', type: 'number', label: 'Gap extend', min: 1, max: 20, default: 2 },
];
const scoring = (v) => ({ alphabet: v.alphabet, matrix: v.alphabet === 'dna' ? 'dna' : v.matrix, match: v.match, mismatch: v.mismatch, gap_open: v.gap_open, gap_extend: v.gap_extend });
// What the engine read of the input, when it stopped early or found lines it could not read.
const reading = (r, out) => {
  const g = r.reading;
  if (!g) return;
  if (g.truncated && g.stopped) out.append(note(`Of ${g.source}: ${g.stopped}.`));
  if (g.problem_count) out.append(note(`${n(g.problem_count)} ${g.problem_count === 1 ? 'line' : 'lines'} could not be read, the first at line ${g.problems[0].line}: ${g.problems[0].message}.`));
};

tool({
  id: 'align-pair', title: 'Pairwise alignment', group: 'alignment', icon: 'align', keywords: ['needleman', 'wunsch', 'smith', 'waterman', 'gotoh', 'wfa', 'global', 'local'],
  summary: 'Two sequences aligned by the schedule you choose, global, local, semi-global or overlap, with its score checked against Gotoh\'s.',
  refs: 'Gotoh (1982); Myers and Miller (1988); Farrar (2007); Marco-Sola et al. (2021, 2023).',
  inputs: [
    { key: 'a', type: 'sequence', label: 'First sequence (the query)', required: true, rows: 3 },
    { key: 'b', type: 'sequence', label: 'Second (the target)', required: true, rows: 3 },
    { group: 'How' },
    { key: 'mode', type: 'choice', label: 'Mode', options: [['global', 'global'], ['local', 'local'], ['semiglobal', 'semi-global'], ['overlap', 'overlap']], default: 'global' },
    { key: 'algorithm', type: 'select', label: 'Schedule', options: [['gotoh', 'Gotoh: the full matrix'], ['hirschberg', 'Myers–Miller: linear space'], ['banded', 'a band of diagonals'], ['wfa', 'WFA: wavefronts'], ['biwfa', 'BiWFA: wavefronts, O(s) memory'], ['striped', 'Farrar\'s stripes (local)'], ['myers', 'Myers\' bit vectors (edit distance)']], default: 'gotoh' },
    ...scoringInputs,
  ],
  example: async () => { const s = parseSequences(await SAMPLES.cox1()); return { a: s[0].seq, b: s[1].seq, mode: 'global', algorithm: 'wfa', alphabet: 'dna' }; },
  async run(v) {
    const mode = v.algorithm === 'striped' ? 'local' : ['wfa', 'biwfa'].includes(v.algorithm) && !['global', 'semiglobal'].includes(v.mode) ? 'global' : v.mode;
    const clean = (s) => (parseSequences(s)[0]?.seq || '');
    return api.run('align.pair', { a: clean(v.a), b: clean(v.b), mode, algorithm: v.algorithm, ...scoring(v) });
  },
  render(r, v, out) {
    out.append(tiles([[n(r.score), 'score'], [pct(r.identity), 'identity'], [n(r.distance), 'edit distance'], [r.reference?.same === false ? 'differs' : r.reference?.same ? 'agrees' : 'another order', 'with Gotoh\'s optimum'], [n(r.cells), `cells, in ${n(r.seconds * 1000, 2)} ms`]]));
    const rows = [{ name: 'query', row: r.rows.query }, { name: 'target', row: r.rows.target }];
    out.append(alignmentView(rows, { alphabet: v.alphabet, extra: [{ name: 'match', row: r.rows.middle, plain: true }], title: `The alignment, ${r.algorithm} (${r.mode})`, file: 'pair' }));
    out.append(textCard('Its CIGAR', 'alignment.cigar.txt', r.cigar, { sub: `query ${n(r.query_begin + 1)}..${n(r.query_end)}, target ${n(r.target_begin + 1)}..${n(r.target_end)}` }));
  },
});

tool({
  id: 'align-msa', title: 'Multiple alignment', group: 'alignment', icon: 'align', keywords: ['msa', 'muscle', 'clustal', 'progressive', 'consensus'],
  summary: 'Many sequences aligned together: progressively along a guide tree (MUSCLE\'s three stages), by the centre star, or by partial orders.',
  refs: 'Edgar (2004), Nucleic Acids Research 32, 1792–1797; Gusfield (1993); Lee, Grasso and Sharlow (2002).',
  inputs: [
    { key: 'sequences', type: 'sequences', label: 'Sequences', required: true, min: 2, kinds: ['fasta'], samples: [['cytochrome b of five primates', SAMPLES.cytb5], ['their proteins', SAMPLES.cytbProteins]] },
    { key: 'method', type: 'choice', label: 'Method', options: [['progressive', 'progressive'], ['centre_star', 'centre star'], ['poa', 'partial order']], default: 'progressive' },
    { key: 'alphabet', type: 'choice', label: 'Letters', options: [['auto', 'as they are'], ['dna', 'DNA'], ['protein', 'protein']], default: 'auto' },
  ],
  example: async () => ({ sequences: { text: await SAMPLES.cytbProteins() }, method: 'progressive', alphabet: 'protein' }),
  run: (v) => api.run('tool.align.msa', { ...source(v.sequences), method: v.method, alphabet: v.alphabet === 'auto' ? (v.sequences?.text ? guessAlphabet(parseSequences(v.sequences.text).map((r) => r.seq).join('')) : 'dna') : v.alphabet }),
  render(r, v, out) {
    out.append(tiles([[n(r.rows.length), 'sequences'], [n(r.columns), 'columns'], [n(r.sp), 'sum-of-pairs score'], [pct(r.identical_columns / Math.max(1, r.columns)), 'columns all alike']]));
    reading(r, out);
    out.append(alignmentView(r.rows, { alphabet: r.alphabet, extra: [{ name: 'consensus', row: r.consensus }, { name: 'conservation', row: r.conservation, plain: true }], title: `Aligned ${r.method === 'progressive' ? 'progressively' : r.method === 'centre_star' ? 'by the centre star' : 'by partial orders'}`, file: 'msa' }));
    if (r.tree) out.append(textCard('The guide tree', 'guide-tree.nwk', r.tree));
    out.append(note('Conservation in Clustal\'s marks: * the column is one residue, : of one strong group, . of one weak group.'));
  },
});

tool({
  id: 'align-search', title: 'Search a query against many', group: 'alignment', icon: 'search', keywords: ['blast', 'search', 'database', 'evalue', 'smith waterman', 'swipe'],
  summary: 'Smith–Waterman of one query against every sequence of a file, one target to each vector lane, ranked by Karlin–Altschul E-values.',
  refs: 'Rognes (2011), BMC Bioinformatics 12, 221; Karlin and Altschul (1990), PNAS 87, 2264–2268.',
  inputs: [
    { key: 'query', type: 'sequence', label: 'Query', required: true, rows: 3 },
    { key: 'targets', type: 'sequences', label: 'Targets', required: true, kinds: ['fasta'] },
    ...scoringInputs.map((i) => (i.key === 'gap_open' ? { ...i, default: 11 } : i.key === 'gap_extend' ? { ...i, default: 1 } : i)),
    { key: 'top', type: 'number', label: 'Show the best', min: 5, max: 500, step: 5, default: 50 },
  ],
  example: async () => { const p = parseSequences(await SAMPLES.cytbProteins()); return { query: p[0].seq.slice(0, 120), targets: { text: await SAMPLES.cytbProteins() }, alphabet: 'protein', matrix: 'BLOSUM62', gap_open: 11, gap_extend: 1 }; },
  run: (v) => api.run('tool.align.search', { query: parseSequences(v.query)[0]?.seq || '', ...source(v.targets), ...scoring(v), top: v.top }),
  render(r, v, out) {
    out.append(tiles([[n(r.searched), 'targets'], [n(r.cells), 'cells'], [`${n(r.gcups, 2)} GCUPS`, r.kernel], [r.statistics ? `λ ${n(r.statistics.lambda, 3)}, K ${n(r.statistics.K, 3)}` : '–', r.statistics?.gapped ? 'gapped statistics (BLAST\'s table)' : 'ungapped statistics']]));
    reading(r, out);
    if (r.note) out.append(note(`${r.note[0].toUpperCase()}${r.note.slice(1)}.`));
    out.append(tableCard('The best hits', [{ key: 'name', label: 'target' }, { key: 'length', label: 'length', num: true }, { key: 'score', label: 'score', num: true }, { key: 'bits', label: 'bits', num: true, format: (x) => n(x, 1) },
      { key: 'evalue', label: 'E-value', num: true, format: (x) => (x < 1e-3 ? x.toExponential(1) : n(x, 3)) }], r.hits, { file: 'hits.csv' }));
    out.append(note('E is the number of alignments at least this good expected by chance in a search of this size: K m n e^(−λS).'));
  },
});

tool({
  id: 'align-dnds', title: 'dN/dS of a gene', group: 'phylogeny', icon: 'trees', keywords: ['selection', 'omega', 'synonymous', 'nei', 'gojobori', 'codon'],
  summary: 'Two coding sequences aligned by their codons, and their synonymous and non-synonymous differences (Nei and Gojobori).',
  refs: 'Nei and Gojobori (1986), Molecular Biology and Evolution 3, 418–426.',
  inputs: [
    { key: 'a', type: 'sequence', label: 'First coding sequence', required: true, rows: 3, alphabet: 'dna' },
    { key: 'b', type: 'sequence', label: 'Second', required: true, rows: 3, alphabet: 'dna' },
    { key: 'table', type: 'select', label: 'Genetic code', options: CODES, default: 1, numeric: true },
  ],
  example: async () => { const s = parseSequences(await SAMPLES.cytb2()); return { a: s[0].seq, b: s[1].seq, table: 2 }; },
  run: (v) => api.run('tool.align.dnds', { a: parseSequences(v.a)[0]?.seq || '', b: parseSequences(v.b)[0]?.seq || '', table: v.table }),
  render(r, v, out) {
    out.append(tiles([[n(r.dN, 4), 'dN'], [n(r.dS, 4), 'dS'], [n(r.omega, 3), 'ω = dN/dS', Number.isFinite(r.omega) && r.omega < 0.5 ? 'good' : ''], [n(r.codons), 'codons compared']]));
    out.append(card({ title: 'Sites and differences', pad: true, body: kv([['synonymous sites S', n(r.S, 2)], ['non-synonymous sites N', n(r.N, 2)], ['synonymous differences', n(r.Sd, 2)], ['non-synonymous differences', n(r.Nd, 2)], ['pS, pN', `${n(r.pS, 4)}, ${n(r.pN, 4)}`], ['codon columns with a gap', n(r.gap_codons)]]) }));
    out.append(alignmentView([{ name: 'first', row: r.protein_a }, { name: 'second', row: r.protein_b }], { alphabet: 'protein', title: 'The proteins, as their codons were aligned', file: 'codon-protein' }));
    out.append(note(!Number.isFinite(r.omega) ? 'ω is not defined here: no synonymous difference (dS = 0), or the differences are saturated.' : r.omega < 1 ? 'ω below 1: most changes of the protein were removed by selection; ω near 1 would be neutral, above 1 favoured.' : 'ω above 1: changes of the protein outnumber what neutrality would give.'));
  },
});
