// Tools of the user's files: reads checked, mapped to a genome and their variants called; a genome
// assembled from its reads; transcripts quantified from RNA-seq; trees built or read; a VCF's
// population; the statistics of an alignment score.

import { el } from '../dom.js';
import { tool } from '../core/registry.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import { token, CATEGORIES } from '../palette.js';
import { tableCard, textCard, treeCard, figure, sequences } from '../ui/results.js';
import { card, kv } from '../ui/kit.js';
import { source, MATRICES, SAMPLES, n, pct, tiles, note } from './common.js';

const readsFile = { key: 'reads_file', type: 'file', label: 'Reads', required: true, kinds: ['fastq', 'fasta'] };
const genome = { key: 'genome', type: 'genome', label: 'Genome' };
// What the engine read of the input, when it stopped early or found lines it could not read.
const reading = (r, out) => {
  const g = r.reading;
  if (!g) return;
  if (g.truncated && g.stopped) out.append(note(`Of ${g.source}: ${g.stopped}.`));
  if (g.problem_count) out.append(note(`${n(g.problem_count)} ${g.problem_count === 1 ? 'line' : 'lines'} could not be read, the first at line ${g.problems[0].line}: ${g.problems[0].message}.`));
};

tool({
  id: 'reads-qc', title: 'Check reads', group: 'reads', icon: 'reads', keywords: ['fastq', 'quality', 'qc', 'fastqc'],
  summary: 'A FASTQ file read through: its reads, their lengths, the quality along the read, GC and N.',
  inputs: [readsFile],
  hint: 'Drop a FASTQ file on the page to check it here.',
  run: (v) => api.run('tool.reads.qc', { reads_file: v.reads_file }),
  render(r, v, out) {
    out.append(tiles([[n(r.reads), 'reads'], [n(r.bases), 'bases'], [n(r.mean_length, 1), 'mean length'], [pct(r.gc), 'GC'], [pct(r.n_rate, 3), 'N']]));
    reading(r, out);
    // Long reads come in bins: of positions along the read (quality_bin) and of lengths (length_bin).
    const qb = r.quality_bin || 1, lb = r.length_bin || 1;
    const q = el('div');
    q.innerHTML = chart({ width: 480, height: 220, x: [1, Math.max(2, r.quality.length * qb)], y: [0, 42], xlabel: 'position in the read', ylabel: 'mean Phred quality',
      series: [{ points: r.quality.map((x, k) => [k * qb + (qb + 1) / 2, x]), colour: token('--accent'), area: true }] });
    const l = el('div');
    const lens = r.lengths.filter(([, c]) => c > 0);
    l.innerHTML = chart({ width: 480, height: 220, x: [Math.max(0, (lens[0]?.[0] || 0) - lb), (lens[lens.length - 1]?.[0] || 1) + lb], y: [0, Math.max(1, ...lens.map(([, c]) => c)) * 1.1], xlabel: 'length', ylabel: 'reads',
      bars: lens.map(([len, c]) => ({ x0: len - 0.45 * lb, x1: len + 0.45 * lb, y: c, colour: token('--C'), title: lb > 1 ? `${n(c)} reads of ${n(len - lb / 2)} to ${n(len + lb / 2)}` : `${n(c)} reads of ${len}` })) });
    out.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '14px' } }, figure('Quality along the read', q), figure('Lengths', l)));
  },
});

tool({
  id: 'reads-map', title: 'Map reads', group: 'reads', icon: 'reads', keywords: ['map', 'align', 'minimap', 'bwa', 'sam', 'coverage'],
  summary: 'Reads placed on a genome (minimizer anchors, chains, base-level alignment), with their mapping qualities, coverage and SAM.',
  refs: 'After Li (2018), Bioinformatics 34, 3094–3100.',
  inputs: [genome, readsFile, { key: 'kind', type: 'choice', label: 'Reads', options: [['short', 'short (Illumina)'], ['long', 'long (8% errors)']], default: 'short' }],
  run: (v) => api.run('tool.reads.map', { genome: v.genome, reads_file: v.reads_file, kind: v.kind }),
  render(r, v, out) {
    out.append(tiles([[n(r.reads), 'reads'], [pct(r.mapped / Math.max(1, r.reads)), 'mapped'], [pct(r.mapq20 / Math.max(1, r.reads)), 'at MAPQ ≥ 20'], [n(r.mean_depth, 1), 'mean depth'], [`${n(r.reads / Math.max(1e-9, r.seconds))}/s`, 'reads mapped']]));
    reading(r, out);
    const c = el('div');
    c.innerHTML = chart({ width: 820, height: 200, x: [0, r.length], y: [0, Math.max(1, ...r.coverage) * 1.15], xlabel: 'position', ylabel: 'depth',
      series: [{ points: r.coverage.map((d, k) => [((k + 0.5) * r.length) / r.coverage.length, d]), colour: token('--accent'), area: true }] });
    out.append(figure('Coverage along the genome', c));
    out.append(textCard('SAM', 'reads.sam', r.sam, { sub: r.sam_truncated ? `the first ${n(r.sam_records)} records` : `${n(r.sam_records)} records` }));
  },
});

tool({
  id: 'reads-call', title: 'Call variants', group: 'reads', icon: 'variant', keywords: ['variant', 'snp', 'snv', 'indel', 'vcf', 'genotype'],
  summary: 'Short reads mapped to a genome and the variants of a diploid called from them: genotype likelihoods, filters, VCF.',
  refs: 'Li (2011), Bioinformatics 27, 2987–2993; Tan, Abecasis and Kang (2015).',
  inputs: [genome, readsFile],
  run: (v) => api.run('tool.reads.call', { genome: v.genome, reads_file: v.reads_file }),
  render(r, v, out) {
    out.append(tiles([[n(r.called ?? r.calls.length), 'variants called'], [n(r.passed), 'passing the filters'], [n(r.snvs), 'SNVs'], [n(r.indels), 'indels'], [n(r.mean_depth, 1), 'mean depth']]));
    reading(r, out);
    out.append(tableCard('The calls', [{ key: 'contig', label: 'contig' }, { key: 'position', label: 'position', num: true, format: (x) => n(x + 1) }, { key: 'ref', label: 'ref', format: (x) => el('span', { class: 'mono' }, x) },
      { key: 'alt', label: 'alt', format: (x) => el('span', { class: 'mono' }, x) }, { key: 'genotype', label: 'genotype' }, { key: 'qual', label: 'QUAL', num: true, format: (x) => n(x, 1) }, { key: 'depth', label: 'depth', num: true }, { key: 'filter', label: 'filter' }], r.calls, { file: 'calls.csv' }));
    out.append(textCard('VCF', 'calls.vcf', r.vcf));
  },
});

tool({
  id: 'assembly-reads', title: 'Assemble reads', group: 'assembly', icon: 'assembly', keywords: ['assemble', 'de bruijn', 'contig', 'n50', 'velvet', 'spades'],
  summary: 'Reads counted into k-mers, the solid ones joined in a de Bruijn graph, its tips, bubbles and weak paths cleaned away; the contigs that remain.',
  refs: 'Pevzner, Tang and Waterman (2001); Zerbino and Birney (2008).',
  inputs: [readsFile, { key: 'k', type: 'number', label: 'k', min: 11, max: 31, step: 2, default: 31 }],
  run: (v) => api.run('tool.assembly.reads', { reads_file: v.reads_file, k: v.k }),
  render(r, v, out) {
    const s = r.stats;
    out.append(tiles([[n(s.contigs), 'contigs'], [n(s.total), 'bases in them'], [n(s.n50), 'N50'], [n(s.largest), 'largest'], [n(r.size_estimate), 'genome size, from the spectrum']]));
    reading(r, out);
    const c = el('div');
    const pts = r.spectrum.filter(([cc]) => cc <= Math.max(30, r.peak * 3));
    c.innerHTML = chart({ width: 520, height: 220, x: [0.5, Math.max(...pts.map((p) => p[0])) + 0.5], y: [0.8, Math.max(...pts.map((p) => p[1])) * 1.6], log: 'y', xlabel: 'times a k-mer is seen', ylabel: 'k-mers',
      bars: pts.map(([cc, k]) => ({ x0: cc - 0.45, x1: cc + 0.45, y: k, colour: cc < r.threshold ? token('--ink-faint') : token('--accent') })) });
    out.append(figure('The k-mer spectrum', c, { sub: `errors below ${r.threshold}, the peak at ${r.peak}` }));
    out.append(textCard('The contigs (FASTA)', 'contigs.fasta', r.contigs_fasta, { sub: `${n(s.contigs)} contigs` }));
  },
});

tool({
  id: 'expr-quant', title: 'Quantify transcripts', group: 'expression', icon: 'expression', keywords: ['rna-seq', 'kallisto', 'tpm', 'pseudoalignment', 'abundance'],
  summary: 'RNA-seq reads pseudoaligned to a transcriptome by k-mer equivalence classes, and each transcript\'s abundance by EM.',
  refs: 'Bray, Pimentel, Melsted and Pachter (2016), Nature Biotechnology 34, 525–527.',
  inputs: [{ key: 'transcripts_file', type: 'file', label: 'Transcripts (FASTA)', required: true, kinds: ['fasta'] }, readsFile,
    { key: 'fragment', type: 'number', label: 'Fragment length', min: 50, max: 1000, step: 10, unit: 'bases', default: 200 }],
  run: (v) => api.run('tool.expr.quant', { transcripts_file: v.transcripts_file, reads_file: v.reads_file, fragment: v.fragment }),
  render(r, v, out) {
    out.append(tiles([[n(r.transcripts.length), 'transcripts'], [n(r.fragments), 'reads'], [pct(r.aligned / Math.max(1, r.fragments)), 'pseudoaligned'], [n(r.classes), 'equivalence classes'], [n(r.rounds), 'EM rounds']]));
    reading(r, out);
    out.append(tableCard('Abundances', [{ key: 'name', label: 'transcript' }, { key: 'length', label: 'length', num: true }, { key: 'effective', label: 'effective length', num: true, format: (x) => n(x, 1) }, { key: 'count', label: 'estimated reads', num: true, format: (x) => n(x, 1) }, { key: 'tpm', label: 'TPM', num: true, format: (x) => n(x, 2) }], r.transcripts, { file: 'abundance.csv' }));
  },
});

tool({
  id: 'tree-build', title: 'Build a tree', group: 'phylogeny', icon: 'trees', keywords: ['phylogeny', 'neighbour joining', 'nj', 'upgma', 'maximum likelihood', 'bootstrap', 'newick'],
  summary: 'Sequences aligned (or already aligned) and their tree: by UPGMA or neighbour joining on corrected distances, or by maximum likelihood, with the bootstrap.',
  refs: 'Saitou and Nei (1987); Felsenstein (1981); Felsenstein (1985).',
  inputs: [
    { key: 'sequences', type: 'sequences', label: 'Sequences (DNA)', required: true, min: 3, kinds: ['fasta'], samples: [['cytochrome b of five primates', SAMPLES.cytb5]] },
    { key: 'aligned', type: 'switch', label: 'Aligned', text: 'they are aligned already', default: false },
    { key: 'method', type: 'choice', label: 'Tree', options: [['upgma', 'UPGMA'], ['nj', 'neighbour joining'], ['ml', 'likelihood']], default: 'nj' },
    { key: 'model', type: 'select', label: 'Model', options: ['JC69', 'K80', 'F81', 'HKY85', 'TN93', 'GTR'].map((m) => [m, m]), default: 'HKY85' },
    { key: 'gamma', type: 'switch', label: 'Rates', text: 'four Γ categories', default: true, when: (v) => v.method === 'ml' },
    { key: 'bootstrap', type: 'number', label: 'Bootstrap', min: 0, max: 1000, step: 50, unit: 'replicates', default: 100 },
  ],
  example: async () => ({ sequences: { text: await SAMPLES.cytb5() }, method: 'ml', model: 'HKY85', gamma: true, bootstrap: 100 }),
  run: (v) => api.run('tool.tree.build', { ...source(v.sequences), aligned: v.aligned, method: v.method, model: v.model, gamma: v.gamma, bootstrap: v.bootstrap }),
  render(r, v, out) {
    out.append(tiles([[n(r.taxa), 'sequences'], [n(r.columns), 'columns'], [n(r.variable), 'variable'], r.log_likelihood !== undefined ? [n(r.log_likelihood, 2), 'log-likelihood'] : null, ...Object.entries(r.parameters || {}).map(([k, x]) => [n(x, 3), k])]));
    out.append(treeCard(`The tree, ${{ upgma: 'by UPGMA', nj: 'by neighbour joining', ml: 'by maximum likelihood' }[v.method]}${v.bootstrap ? `, ${v.bootstrap} bootstrap replicates` : ''}`, r.tree.nodes, { newick: r.newick }));
    out.append(note('Support at a node is the fraction of bootstrap trees holding its split: green from 95%, red below 70%.'));
  },
});

tool({
  id: 'tree-read', title: 'Read a tree', group: 'phylogeny', icon: 'trees', keywords: ['newick', 'draw', 'tree', 'reroot', 'midpoint'],
  summary: 'A tree in Newick drawn, rooted at its midpoint or as written, with its leaves and length.',
  inputs: [{ key: 'newick', type: 'sequence', label: 'Newick', rows: 4, required: true, placeholder: '((A:0.1,B:0.2):0.05,C:0.3,D:0.4);  or drop a .nwk file' },
    { key: 'root', type: 'choice', label: 'Root', options: [['as', 'as written'], ['midpoint', 'at the midpoint']], default: 'midpoint' }],
  // the example's tree is computed, not typed: the five primates' cytochrome b by neighbour joining
  example: async () => ({ newick: (await api.run('tool.tree.build', { text: await SAMPLES.cytb5(), aligned: false, method: 'nj', model: 'HKY85', bootstrap: 0 })).newick, root: 'midpoint' }),
  run: (v) => api.run('tool.tree.read', { text: v.newick, root: v.root }),
  render(r, v, out) {
    out.append(tiles([[n(r.leaves), 'leaves'], [n(r.length, 4), 'tree length'], [r.binary ? 'yes' : 'no', 'binary']]));
    out.append(treeCard('The tree', r.tree.nodes, { newick: r.newick }));
  },
});

tool({
  id: 'pop-vcf', title: 'A population from its VCF', group: 'population', icon: 'populations', keywords: ['vcf', 'pca', 'diversity', 'tajima', 'structure', 'genotype'],
  summary: 'The genotypes of a VCF: principal components of the samples, and diversity (θ_W, π, Tajima\'s D).',
  refs: 'Patterson, Price and Reich (2006), PLoS Genetics 2, e190; Tajima (1989), Genetics 123, 585–595.',
  inputs: [{ key: 'vcf_file', type: 'file', label: 'Variants (VCF)', required: true, kinds: ['vcf'] }],
  run: (v) => api.run('tool.pop.vcf', { vcf_file: v.vcf_file }),
  render(r, v, out) {
    out.append(tiles([[n(r.sample_count ?? r.samples.length), 'samples'], [n(r.sites), 'biallelic sites'], [n(r.diversity.segregating), 'segregating'], [n(r.diversity.pi, 2), 'π'], [n(r.diversity.theta, 2), 'θ_W'], [n(r.diversity.tajima_d, 3), 'Tajima\'s D']]));
    reading(r, out);
    if (r.pca) {
      const xs = r.pca.pc1, ys = r.pca.pc2;
      const c = el('div');
      c.innerHTML = chart({ width: 520, height: 380, x: [Math.min(...xs) * 1.15, Math.max(...xs) * 1.15], y: [Math.min(...ys) * 1.15, Math.max(...ys) * 1.15], xlabel: `PC 1 (${pct(r.pca.explained[0])})`, ylabel: `PC 2 (${pct(r.pca.explained[1])})`,
        series: [{ points: xs.map((x, k) => [x, ys[k]]), colour: token('--accent'), dots: true, line: false }] });
      out.append(figure('The samples on their first two components', c));
    }
    out.append(card({ title: 'What was read', pad: true, body: kv(Object.entries(r.conversion || {}).map(([k, x]) => [k.replace(/_/g, ' '), String(x)])) }));
  },
});

tool({
  id: 'stats-evalue', title: 'Score statistics', group: 'statistics', icon: 'statistics', keywords: ['evalue', 'karlin', 'altschul', 'bit score', 'significance', 'blast'],
  summary: 'What a local alignment score means: λ, K and H of a scoring, and the score\'s bits, E-value and P-value in a search of a given size.',
  refs: 'Karlin and Altschul (1990), PNAS 87, 2264–2268; Altschul and Gish (1996).',
  inputs: [{ key: 'score', type: 'number', label: 'Raw score', min: 1, max: 100000, default: 100 },
    { key: 'matrix', type: 'select', label: 'Scores', options: MATRICES.map((m) => [m, m]), default: 'BLOSUM62' },
    { key: 'gap_open', type: 'number', label: 'Gap open', min: 0, max: 50, default: 11 }, { key: 'gap_extend', type: 'number', label: 'Gap extend', min: 1, max: 20, default: 1 },
    { key: 'm', type: 'number', label: 'Query length', min: 1, max: 1e9, default: 300 }, { key: 'n', type: 'number', label: 'Database length', min: 1, max: 1e12, step: 1000, default: 1e8, unit: 'residues' }],
  example: { score: 100, matrix: 'BLOSUM62', gap_open: 11, gap_extend: 1, m: 300, n: 1e8 },
  run: (v) => api.run('tool.stats.evalue', v),
  render(r, v, out) {
    out.append(tiles([[n(r.bits, 1), 'bits'], [r.evalue < 1e-3 ? r.evalue.toExponential(2) : n(r.evalue, 4), 'E-value'], [r.pvalue < 1e-3 ? r.pvalue.toExponential(2) : n(r.pvalue, 4), 'P-value']]));
    out.append(card({ title: r.gapped ? 'Gapped statistics, from BLAST\'s table' : 'Ungapped statistics, computed from the matrix', pad: true, body: kv([['λ', n(r.lambda, 4)], ['K', n(r.K, 4)], ['H (nats a pair)', n(r.H, 4)], ['scoring', `${v.matrix}, gaps ${v.gap_open} + L·${v.gap_extend}`]]) }));
    out.append(note('A score S in bits is (λS − ln K)/ln 2; E = m n 2^(−bits). Ungapped λ is the root of Σ p_i p_j e^(λ s_ij) = 1; gapped values come from simulation, as BLAST tabulates them.'));
  },
});

export { sequences, CATEGORIES };
