// Tools of annotation and indexes: genes, CpG islands and repeats found from the sequence alone;
// the suffix array and Burrows–Wheeler transform of a text; genomes compared by their sketches.

import { el } from '../dom.js';
import { tool } from '../core/registry.js';
import * as api from '../api.js';
import { tableCard, treeCard, figure } from '../ui/results.js';
import { token, warm } from '../palette.js';
import { source, CODES, SAMPLES, n, pct, tiles, note } from './common.js';

const oneSequence = (extra = {}) => ({ key: 'sequences', type: 'sequences', label: 'Sequence', required: true, kinds: ['fasta', 'genbank'], samples: [['phage lambda, whole', SAMPLES.lambdaWhole], ['its first 4 kb', SAMPLES.lambda]], ...extra });
// What the engine read of the input, when it stopped early or found lines it could not read.
const reading = (r, out) => {
  const g = r.reading;
  if (!g) return;
  if (g.truncated && g.stopped) out.append(note(`Of ${g.source}: ${g.stopped}.`));
  if (g.problem_count) out.append(note(`${n(g.problem_count)} ${g.problem_count === 1 ? 'line' : 'lines'} could not be read, the first at line ${g.problems[0].line}: ${g.problems[0].message}.`));
};

function track(L, items, colour, title) {
  const W = 820;
  let s = `<svg width="${W}" height="54" viewBox="0 0 ${W} 54"><line x1="10" y1="27" x2="${W - 10}" y2="27" stroke="${token('--hair-2')}"/>`;
  for (const [b, e, up, label] of items.slice(0, 6000)) {
    const x0 = 10 + (b / L) * (W - 20), x1 = 10 + (e / L) * (W - 20);
    s += `<rect x="${x0.toFixed(1)}" y="${up === false ? 30 : up === true ? 16 : 20}" width="${Math.max(1, x1 - x0).toFixed(1)}" height="${up === undefined ? 14 : 8}" rx="1.5" fill="${colour}" fill-opacity="0.8"><title>${label}</title></rect>`;
  }
  return figure(title, `${s}<text class="tick-label" x="10" y="52">1</text><text class="tick-label" x="${W - 10}" y="52" text-anchor="end">${n(L)}</text></svg>`);
}

tool({
  id: 'annot-genes', title: 'Find genes', group: 'annotation', icon: 'annotation', keywords: ['gene', 'prediction', 'genemark', 'prokaryote', 'cds'],
  summary: 'Genes of a prokaryote found from the sequence alone: periodic Markov chains of coding trained on the genome itself.',
  refs: 'After GeneMarkS: Besemer, Lomsadze and Borodovsky (2001), Nucleic Acids Research 29, 2607–2618.',
  inputs: [oneSequence(), { key: 'table', type: 'select', label: 'Genetic code', options: CODES, default: 11, numeric: true }],
  example: async () => ({ sequences: { text: await SAMPLES.lambdaWhole() }, table: 11 }),
  run: (v) => api.run('tool.annot.genes', { ...source(v.sequences), table: v.table }),
  render(r, v, out) {
    out.append(tiles([[n(r.count ?? r.genes.length), 'genes'], [n(r.length), 'bases'], [pct(r.coding), 'of the sequence coding'], [n(r.rounds), 'rounds of training']]));
    reading(r, out);
    out.append(track(r.length, r.genes.map((g) => [g.begin, g.end, g.strand > 0, `${g.begin + 1}..${g.end} ${g.strand > 0 ? '+' : '−'}, ${g.codons} codons, ${g.score} bits`]), token('--cds'), 'The genes, forward above and reverse below'));
    out.append(tableCard('The genes', [{ key: 'begin', label: 'from', num: true, format: (x) => n(x + 1) }, { key: 'end', label: 'to', num: true, format: (x) => n(x) }, { key: 'strand', label: 'strand', format: (x) => (x > 0 ? '+' : '−') },
      { key: 'codons', label: 'codons', num: true }, { key: 'score', label: 'bits', num: true, format: (x) => n(x, 1) }, { key: 'start', label: 'start', format: (x) => el('span', { class: 'mono' }, x) }, { key: 'protein', label: 'protein', format: (x) => el('span', { class: 'mono', title: x }, x.length > 48 ? `${x.slice(0, 46)}…` : x) }], r.genes, { file: 'genes.csv' }));
    if (r.length < 200000) out.append(note('The coding chains are trained on the sequence itself; a few hundred kilobases give them enough to learn from, less makes them guess.'));
  },
});

tool({
  id: 'annot-cpg', title: 'CpG islands', group: 'annotation', icon: 'annotation', keywords: ['cpg', 'island', 'hmm', 'methylation'],
  summary: 'Stretches rich in CpG, found by the two-state hidden Markov model of Durbin et al.',
  refs: 'Durbin, Eddy, Krogh and Mitchison (1998), Biological Sequence Analysis, ch. 3.',
  inputs: [oneSequence()],
  example: async () => ({ sequences: { text: await SAMPLES.lambdaWhole() } }),
  run: (v) => api.run('tool.annot.cpg', source(v.sequences)),
  render(r, v, out) {
    out.append(tiles([[n(r.count ?? r.islands.length), 'islands'], [n(r.length), 'bases']]));
    reading(r, out);
    out.append(track(r.length, r.islands.map((x) => [x.begin, x.end, undefined, `${x.begin + 1}..${x.end}, GC ${(x.gc * 100).toFixed(1)}%, CpG o/e ${x.oe}`]), token('--G'), 'The islands'));
    out.append(tableCard('The islands', [{ key: 'begin', label: 'from', num: true, format: (x) => n(x + 1) }, { key: 'end', label: 'to', num: true, format: (x) => n(x) }, { key: 'gc', label: 'GC', num: true, format: (x) => pct(x) }, { key: 'oe', label: 'CpG o/e', num: true, format: (x) => n(x, 2) }, { key: 'log_odds', label: 'bits', num: true, format: (x) => n(x, 1) }], r.islands, { file: 'cpg-islands.csv' }));
    out.append(note('The model\'s background is vertebrate DNA, depleted of CpG; a bacterium is not, and the model flags much of it.'));
  },
});

tool({
  id: 'annot-repeats', title: 'Repeats', group: 'annotation', icon: 'annotation', keywords: ['repeat', 'maximal', 'duplicate', 'inverted'],
  summary: 'Maximal repeats on both strands, at least so long, found by the suffix array and its LCP.',
  refs: 'Gusfield (1997), Algorithms on Strings, Trees, and Sequences, §7.12.',
  inputs: [oneSequence(), { key: 'min_length', type: 'number', label: 'At least', min: 10, max: 100000, step: 10, unit: 'bases', default: 30 }],
  example: async () => ({ sequences: { text: await SAMPLES.lambdaWhole() }, min_length: 20 }),
  run: (v) => api.run('tool.annot.repeats', { ...source(v.sequences), min_length: v.min_length }),
  render(r, v, out) {
    out.append(tiles([[n(r.count), 'repeats'], [n(r.longest), 'longest'], [n(r.length), 'bases']]));
    reading(r, out);
    out.append(tableCard('The repeats', [{ key: 'length', label: 'length', num: true }, { key: 'copies', label: 'copies', num: true }, { key: 'where', label: 'where (a reversed copy marked −)', format: (x) => el('span', { class: 'mono small' }, x) }, { key: 'word', label: 'the word', format: (x) => el('span', { class: 'mono', title: x }, x.length > 40 ? `${x.slice(0, 38)}…` : x) }], r.repeats, { file: 'repeats.csv' }));
  },
});

tool({
  id: 'index-sa', title: 'Suffix array and BWT', group: 'index', icon: 'index', keywords: ['suffix array', 'bwt', 'burrows wheeler', 'lcp', 'lyndon', 'fm index'],
  summary: 'A text\'s suffixes sorted (SA-IS), their longest common prefixes, its Burrows–Wheeler transform and Lyndon factors.',
  refs: 'Nong, Zhang and Chan (2011); Burrows and Wheeler (1994); Duval (1983).',
  inputs: [{ key: 'text', type: 'text', label: 'Text', required: true, default: 'mississippi', placeholder: 'any word, or a short sequence' }],
  example: { text: 'abracadabra' },
  run: (v) => api.run('index.word', { text: v.text }),
  render(r, v, out) {
    const rows = r.sa.map((p, k) => ({ rank: k, position: p, lcp: r.lcp?.[k] ?? '', suffix: v.text.slice(p) + '$' }));
    out.append(tiles([[el('span', { class: 'mono' }, r.bwt), 'its BWT'], [n(r.runs ?? 0), 'runs in the BWT'], [el('span', { class: 'mono' }, (r.lyndon || []).map(([b, e]) => v.text.slice(b, e)).join(' · ')), 'Lyndon factors']]));
    out.append(tableCard('The suffixes, sorted', [{ key: 'rank', label: 'rank', num: true }, { key: 'position', label: 'starts at', num: true, format: (x) => n(x + 1) }, { key: 'lcp', label: 'LCP', num: true }, { key: 'suffix', label: 'suffix', format: (x) => el('span', { class: 'mono' }, x) }], rows, { file: 'suffix-array.csv' }));
  },
});

tool({
  id: 'sketch-compare', title: 'Compare genomes by sketches', group: 'index', icon: 'sketch', keywords: ['mash', 'minhash', 'distance', 'ani', 'jaccard'],
  summary: 'Each genome\'s k-mers hashed and the smallest kept; any two sketches estimate their Jaccard index and Mash distance, and the distances a tree.',
  refs: 'Ondov et al. (2016), Genome Biology 17, 132.',
  inputs: [{ key: 'genomes', type: 'text', label: 'Genomes', default: (ctx) => ctx.catalogue.filter((g) => g.fetched && g.group === 'mitochondria').map((g) => g.id).join(', '), help: 'ids of the catalogue or of your files, separated by commas' },
    { key: 'k', type: 'number', label: 'k', min: 9, max: 31, step: 2, default: 21 }, { key: 'size', type: 'number', label: 'Sketch size', min: 100, max: 100000, step: 100, default: 1000 }],
  example: (ctx) => ({ genomes: ctx.catalogue.filter((g) => g.fetched && g.group === 'mitochondria').map((g) => g.id).join(', ') }),
  run: (v) => api.run('trees.sketches', { genomes: v.genomes.split(/[\s,]+/).filter(Boolean), k: v.k, size: v.size }),
  render(r, v, out) {
    out.append(tiles([[n(r.genomes), 'genomes'], [n(r.pairs), 'distances'], [`${n(r.sketch_seconds * 1000, 1)} ms`, 'sketching'], [`${n(r.distance_seconds * 1000, 1)} ms`, 'distances']]));
    out.append(treeCard('Their tree, by neighbour joining on Mash distances', r.tree.nodes, { scaleLabel: 'Mash distance' }));
    const m = r.names.length, cell = Math.max(14, Math.min(34, Math.floor(460 / m)));
    let s = `<svg width="${m * cell + 220}" height="${m * cell + 10}" viewBox="0 0 ${m * cell + 220} ${m * cell + 10}">`;
    const max = Math.max(1e-9, ...r.distances.flat().filter(Number.isFinite));
    r.distances.forEach((row, i) => row.forEach((d, j) => { s += `<rect x="${j * cell}" y="${i * cell}" width="${cell - 1}" height="${cell - 1}" fill="${warm(1 - (d ?? 1) / max)}"><title>${r.names[i]} — ${r.names[j]}: ${d === null ? '–' : d.toFixed(4)}</title></rect>`; }));
    r.names.forEach((nm, i) => { s += `<text x="${m * cell + 8}" y="${i * cell + cell * 0.65}" style="font-size:12px;font-style:italic;fill:${token('--ink')}">${nm}</text>`; });
    out.append(figure('The distances, darker closer', `${s}</svg>`));
  },
});
