// Tools of sequences: what is in them, the other strand, their translations and open reading
// frames, motifs found on both strands, k-mers, codon usage.

import { el } from '../dom.js';
import { tool } from '../core/registry.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import { token } from '../palette.js';
import { sequences, tableCard, figure } from '../ui/results.js';
import { source, CODES, SAMPLES, n, pct, tiles, note } from './common.js';

const seqs = (extra = {}) => ({ key: 'sequences', type: 'sequences', label: 'Sequences', required: true, kinds: ['fasta'], samples: [['cytochrome b, human and chimpanzee', SAMPLES.cytb2], ['the start of phiX174', SAMPLES.phix]], ...extra });
// What the engine read of the input, when it stopped early or found lines it could not read.
const reading = (r, out) => {
  const g = r.reading;
  if (!g) return;
  if (g.truncated && g.stopped) out.append(note(`Of ${g.source}: ${g.stopped}.`));
  if (g.problem_count) out.append(note(`${n(g.problem_count)} ${g.problem_count === 1 ? 'line' : 'lines'} could not be read, the first at line ${g.problems[0].line}: ${g.problems[0].message}.`));
};

tool({
  id: 'seq-stats', title: 'Composition', group: 'sequence', icon: 'sequence', keywords: ['gc', 'content', 'length', 'n50', 'skew', 'statistics'],
  summary: 'Lengths, GC and its skew, CpG observed over expected, ambiguity and entropy, sequence by sequence.',
  inputs: [seqs()],
  example: async () => ({ sequences: { text: await SAMPLES.cytb2() } }),
  run: (v) => api.run('tool.seq.stats', source(v.sequences)),
  render(r, v, out) {
    const t = r.total;
    out.append(tiles([[n(t.records ?? r.records.length), 'sequences'], [n(t.length), 'bases'], [pct(t.gc), 'GC'], [n(t.n50), 'N50'], [n(t.longest), 'longest'], [n(t.ambiguous), 'not A, C, G or T']]));
    reading(r, out);
    out.append(tableCard('Each sequence', [
      { key: 'name', label: 'name' }, { key: 'length', label: 'length', num: true, format: (x) => n(x) }, { key: 'gc', label: 'GC', num: true, format: (x) => pct(x) },
      { key: 'gc_skew', label: 'GC skew', num: true, format: (x) => n(x, 3) }, { key: 'cpg_oe', label: 'CpG o/e', num: true, format: (x) => n(x, 2) },
      { key: 'entropy', label: 'entropy of 2-mers', num: true, format: (x) => n(x, 3) }, { key: 'ambiguous', label: 'N', num: true, format: (x) => n(x) }], r.records, { file: 'composition.csv' }));
    out.append(note('GC skew (G − C)/(G + C) changes sign at a bacterial origin of replication; CpG observed over expected is low in vertebrate DNA, where methylated CpG decays to TpG, and near 1 in CpG islands. Entropy of 2-mers is 4 bits for random sequence, less where it repeats.'));
  },
});

tool({
  id: 'seq-revcomp', title: 'Reverse complement', group: 'sequence', icon: 'sequence', keywords: ['complement', 'strand', 'reverse'],
  summary: 'The other strand, read 5′ to 3′; IUPAC ambiguity letters complemented too.',
  inputs: [seqs(), { key: 'what', type: 'choice', label: 'Give', options: [['revcomp', 'reverse complement'], ['complement', 'complement'], ['reverse', 'reverse']], default: 'revcomp' }],
  example: async () => ({ sequences: { text: await SAMPLES.phix() } }),
  run: (v) => api.run('tool.seq.revcomp', { ...source(v.sequences), what: v.what }),
  render: (r, v, out) => out.append(sequences(r.records, { title: { revcomp: 'Reverse complements', complement: 'Complements', reverse: 'Reversed' }[v.what], file: `${v.what}.fasta` })),
});

tool({
  id: 'seq-translate', title: 'Translate', group: 'sequence', icon: 'sequence', keywords: ['protein', 'frame', 'genetic code', 'codon'],
  summary: 'Proteins in one, three or six frames, under any of NCBI\'s 25 genetic codes.',
  inputs: [seqs(), { key: 'table', type: 'select', label: 'Genetic code', options: CODES, default: 1, numeric: true },
    { key: 'frames', type: 'choice', label: 'Frames', options: [['1', 'frame 1'], ['3', 'three forward'], ['6', 'all six']], default: '1' }],
  example: async () => ({ sequences: { text: await SAMPLES.cytb2() }, table: 2, frames: '1' }),
  run: (v) => api.run('tool.seq.translate', { ...source(v.sequences), table: v.table, frames: Number(v.frames) }),
  render: (r, v, out) => out.append(sequences(r.records.map((x) => ({ name: `${x.name} frame ${x.frame}`, seq: x.protein, note: `${n((x.protein.match(/\*/g) || []).length)} stops` })), { alphabet: 'protein', title: 'Translations', file: 'translations.fasta' })),
});

tool({
  id: 'seq-orfs', title: 'Open reading frames', group: 'sequence', icon: 'annotation', keywords: ['orf', 'gene', 'protein', 'start', 'stop'],
  summary: 'Every ORF on both strands from a start codon to a stop, at least so many codons, with its protein.',
  inputs: [seqs({ label: 'Sequence', samples: [['the start of lambda', SAMPLES.lambda]] }), { key: 'table', type: 'select', label: 'Genetic code', options: CODES, default: 11, numeric: true },
    { key: 'min_codons', type: 'number', label: 'At least', min: 10, max: 5000, step: 10, unit: 'codons', default: 60 },
    { key: 'any_start', type: 'switch', label: 'Starts', text: 'every start codon of the code, not only ATG', default: false }],
  example: async () => ({ sequences: { text: await SAMPLES.lambda() }, table: 11, min_codons: 50 }),
  run: (v) => api.run('tool.seq.orfs', { ...source(v.sequences), table: v.table, min_codons: v.min_codons, any_start: v.any_start }),
  render(r, v, out) {
    out.append(tiles([[n(r.orfs.length), 'ORFs'], [n(r.length), 'bases searched'], [n(r.orfs.reduce((a, o) => Math.max(a, o.codons), 0)), 'longest, codons']]));
    reading(r, out);
    if (r.orfs.length && r.orfs.length < 4000) {
      const L = r.length, W = 820;
      let s = `<svg width="${W}" height="70" viewBox="0 0 ${W} 70"><line x1="10" y1="35" x2="${W - 10}" y2="35" stroke="${token('--hair-2')}"/>`;
      for (const o of r.orfs) {
        const x0 = 10 + (o.begin / L) * (W - 20), x1 = 10 + (o.end / L) * (W - 20), y = o.strand > 0 ? 20 - (o.frame % 3) * 5 : 42 + (o.frame % 3) * 5;
        s += `<rect x="${x0.toFixed(1)}" y="${y}" width="${Math.max(1, x1 - x0).toFixed(1)}" height="4" rx="1" fill="${token('--accent')}" fill-opacity="0.75"><title>${o.begin + 1}..${o.end} ${o.strand > 0 ? '+' : '−'} frame ${o.frame}, ${o.codons} codons</title></rect>`;
      }
      out.append(figure('Where they are', `${s}<text class="tick-label" x="10" y="66">1</text><text class="tick-label" x="${W - 10}" y="66" text-anchor="end">${n(L)}</text></svg>`, { sub: 'forward strand above, reverse below, one line for each frame' }));
    }
    out.append(tableCard('The ORFs', [{ key: 'begin', label: 'from', num: true, format: (x) => n(x + 1) }, { key: 'end', label: 'to', num: true, format: (x) => n(x) }, { key: 'strand', label: 'strand', format: (x) => (x > 0 ? '+' : '−') },
      { key: 'frame', label: 'frame', num: true }, { key: 'codons', label: 'codons', num: true }, { key: 'protein', label: 'protein', format: (x) => el('span', { class: 'mono', title: x }, x.length > 60 ? `${x.slice(0, 58)}…` : x) }], r.orfs, { file: 'orfs.csv' }));
  },
});

tool({
  id: 'seq-find', title: 'Find motifs', group: 'index', icon: 'search', keywords: ['motif', 'pattern', 'site', 'iupac', 'restriction', 'search'],
  summary: 'Every occurrence of motifs (IUPAC letters allowed: N, R, Y, W…) on both strands, by one Aho–Corasick automaton.',
  refs: 'Aho and Corasick (1975), Communications of the ACM 18, 333–340.',
  inputs: [seqs({ label: 'Search in' }), { key: 'patterns', type: 'text', label: 'Motifs', placeholder: 'GAATTC, GGATCC, TATAWAW …', default: 'GAATTC, GGATCC, AAGCTT, CTGCAG, GCGGCCGC', help: 'separated by commas or spaces' },
    { key: 'both', type: 'switch', label: 'Strands', text: 'both strands', default: true }],
  example: async () => ({ sequences: { text: await SAMPLES.lambda() }, patterns: 'GAATTC, GGATCC, AAGCTT, CTGCAG, GGWCC' }),
  run: (v) => api.run('tool.seq.find', { ...source(v.sequences), patterns: v.patterns.split(/[\s,;]+/).filter(Boolean), both: v.both }),
  render(r, v, out) {
    out.append(tiles([[n(r.count), 'occurrences'], [n(r.patterns), 'words searched (IUPAC expanded)'], [n(r.searched), 'bases searched']]));
    reading(r, out);
    out.append(tableCard('By motif', [{ key: 'motif', label: 'motif' }, { key: 'count', label: 'occurrences', num: true }, { key: 'forward', label: 'forward', num: true }, { key: 'reverse', label: 'reverse', num: true }], r.by_motif, { file: 'motifs.csv' }));
    out.append(tableCard('Where', [{ key: 'record', label: 'sequence' }, { key: 'motif', label: 'motif' }, { key: 'begin', label: 'from', num: true, format: (x) => n(x + 1) }, { key: 'end', label: 'to', num: true, format: (x) => n(x) },
      { key: 'strand', label: 'strand', format: (x) => (x > 0 ? '+' : '−') }, { key: 'match', label: 'as found', format: (x) => el('span', { class: 'mono' }, x) }], r.hits, { file: 'hits.csv' }));
  },
});

tool({
  id: 'seq-kmers', title: 'k-mer spectrum', group: 'sequence', icon: 'statistics', keywords: ['kmer', 'spectrum', 'count', 'histogram'],
  summary: 'Canonical k-mers counted: how many are seen once, twice, …, and the commonest.',
  inputs: [seqs({ kinds: ['fasta', 'fastq'] }), { key: 'k', type: 'number', label: 'k', min: 3, max: 31, step: 2, default: 21 }, { key: 'top', type: 'number', label: 'The commonest', min: 5, max: 500, step: 5, default: 25 }],
  example: async () => ({ sequences: { text: await SAMPLES.lambda() }, k: 7 }),
  run: (v) => api.run('tool.seq.kmers', { ...source(v.sequences, v.sequences?.file && /fastq/.test(v.sequences.file) ? 'reads' : 'fasta'), k: v.k, top: v.top }),
  render(r, v, out) {
    out.append(tiles([[n(r.distinct), 'distinct k-mers'], [n(r.total), 'k-mers in all'], [pct(r.distinct / Math.max(1, Math.min(r.total, 4 ** v.k / 2))), 'of the possible, or of all']]));
    reading(r, out);
    const pts = r.spectrum.filter(([c]) => c <= 80);
    const c = el('div');
    c.innerHTML = chart({ width: 560, height: 240, x: [0.5, Math.max(5, ...pts.map((p) => p[0])) + 0.5], y: [0.8, Math.max(...pts.map((p) => p[1])) * 1.6], log: 'y', xlabel: 'times a k-mer is seen', ylabel: 'k-mers',
      bars: pts.map(([cc, k]) => ({ x0: cc - 0.45, x1: cc + 0.45, y: k, colour: token('--accent'), title: `${n(k)} k-mers seen ${cc} times` })) });
    out.append(figure('The spectrum', c));
    out.append(tableCard('The commonest', [{ key: 'kmer', label: 'k-mer', format: (x) => el('span', { class: 'mono' }, x) }, { key: 'count', label: 'seen', num: true }], r.top.map(([kmer, count]) => ({ kmer, count })), { file: 'kmers.csv' }));
  },
});

tool({
  id: 'seq-codons', title: 'Codon usage', group: 'sequence', icon: 'sequence', keywords: ['codon', 'rscu', 'enc', 'bias'],
  summary: 'Codons counted over coding sequences, their relative synonymous usage, and the effective number of codons.',
  refs: 'Sharp and Li (1987), Nucleic Acids Research 15, 1281–1295; Wright (1990), Gene 87, 23–29.',
  inputs: [seqs({ label: 'Coding sequences' }), { key: 'table', type: 'select', label: 'Genetic code', options: CODES, default: 1, numeric: true }],
  example: async () => ({ sequences: { text: await SAMPLES.cytb5() }, table: 2 }),
  run: (v) => api.run('tool.seq.codons', { ...source(v.sequences), table: v.table }),
  render(r, v, out) {
    out.append(tiles([[n(r.codons_counted), 'codons'], [n(r.enc, 1), 'effective number of codons (20 to 61)']]));
    reading(r, out);
    out.append(tableCard('Each codon', [{ key: 'codon', label: 'codon', format: (x) => el('span', { class: 'mono' }, x) }, { key: 'aa', label: 'amino acid' }, { key: 'count', label: 'count', num: true }, { key: 'rscu', label: 'RSCU', num: true, format: (x) => n(x, 2) }], r.codons, { file: 'codons.csv' }));
  },
});
