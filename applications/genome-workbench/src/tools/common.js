// What the tools share: sequences given as text or as a file of your data, turned into the
// engine's parameters; the genetic codes; small builders of results.

import { el } from '../dom.js';
import { metrics } from '../ui/kit.js';
import * as api from '../api.js';

// {text} or {file} (ui/form.js, the 'sequences' input), as the engine's parameters: the text, or a
// file of the library under `${name}_file` (the server turns its id into its path).
export function source(v, name = 'fasta') {
  if (v && typeof v === 'object' && v.file) return { [`${name}_file`]: v.file };
  return { text: typeof v === 'string' ? v : v?.text || '' };
}

export const CODES = [
  [1, 'Standard'], [2, 'Vertebrate mitochondrial'], [3, 'Yeast mitochondrial'], [4, 'Mold, protozoan, Mycoplasma'], [5, 'Invertebrate mitochondrial'],
  [6, 'Ciliate nuclear'], [9, 'Echinoderm mitochondrial'], [10, 'Euplotid nuclear'], [11, 'Bacterial, archaeal and plant plastid'], [12, 'Alternative yeast nuclear'],
  [13, 'Ascidian mitochondrial'], [14, 'Alternative flatworm mitochondrial'], [16, 'Chlorophycean mitochondrial'], [21, 'Trematode mitochondrial'], [22, 'Scenedesmus mitochondrial'],
  [23, 'Thraustochytrium mitochondrial'], [24, 'Rhabdopleuridae mitochondrial'], [25, 'Candidate division SR1 and Gracilibacteria'], [26, 'Pachysolen nuclear'], [27, 'Karyorelict nuclear'],
  [28, 'Condylostoma nuclear'], [29, 'Mesodinium nuclear'], [30, 'Peritrich nuclear'], [31, 'Blastocrithidia nuclear'], [33, 'Cephalodiscidae mitochondrial'],
].map(([v, l]) => [v, `${v} · ${l}`]);

export const MATRICES = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250'];

export const n = (x, d = 0) => (Number.isFinite(x) ? x.toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d }) : '–');
export const pct = (x, d = 1) => (Number.isFinite(x) ? `${(100 * x).toFixed(d)}%` : '–');

export const tiles = (items) => metrics(items);
export const note = (text) => el('p', { class: 'note', style: { margin: '0 2px' } }, text);

// Samples to try: real sequences, read when asked from the RefSeq genomes the recipe fetched
// (data/registry/ncbi-genomes.json), never typed in. Each is an async function giving FASTA.
const PRIMATES = [['ncbi:NC_012920.1', 'human'], ['ncbi:NC_001643.1', 'chimpanzee'], ['ncbi:NC_011120.1', 'gorilla'], ['ncbi:NC_002083.1', 'orangutan'], ['ncbi:NC_002082.1', 'gibbon']];
async function genes(gene, which, protein) {
  const out = [];
  for (const [id, name] of which) {
    const g = await api.run('genome.gene', { genome: id, gene });
    out.push(`>${name}_${gene}\n${protein ? g.translation : g.sequence}`);
  }
  return out.join('\n') + '\n';
}
async function stretch(id, name, begin, end) {
  const t = await api.run('genome.tracks', { genome: id, begin, end, bins: 4, letters_max: end - begin + 10 });
  return `>${name}_${begin + 1}_${end}\n${t.letters}\n`;
}
export const SAMPLES = {
  cytb2: () => genes('CYTB', PRIMATES.slice(0, 2), false),
  cytb5: () => genes('CYTB', PRIMATES, false),
  cytbProteins: () => genes('CYTB', PRIMATES, true),
  cox1: () => genes('COX1', PRIMATES.slice(0, 2), false),
  phix: () => stretch('ncbi:NC_001422.1', 'phiX174', 0, 1600),
  lambda: () => stretch('ncbi:NC_001416.1', 'lambda', 0, 4000),
  lambdaWhole: () => stretch('ncbi:NC_001416.1', 'lambda', 0, 48502),
};
