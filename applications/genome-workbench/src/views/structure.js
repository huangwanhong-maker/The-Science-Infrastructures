// The Rearrangements view: structural variants and copy number found from reads (genomics::structural,
// ADR 0064). A rearranged genome is a set of new adjacencies between oriented breakends, and each kind
// of rearrangement leaves its signature in three kinds of evidence:
//
//   pairs     mates farther apart than the library allows, or facing the wrong way;
//   splits    reads whose two parts map on either side of a junction, to the base;
//   depth     the number of fragments per bin, against the GC and mappability that shape it,
//             segmented into integer copy numbers; and at heterozygous sites, the B-allele fraction,
//             which tells a copy-neutral loss of heterozygosity from a deletion that depth alone
//             cannot.
//
// Every call is genotyped from the reads for and against it and said in the three notations: VCF,
// HGVS and ISCN. The individuals are simulated on real references (GRCh38's chromosomes 21 and 22;
// Escherichia coli K-12), labelled so, and their truth is drawn beneath every call.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { chart } from '../draw/chart.js';
import { genomeTracksSvg } from '../draw/tracks.js';
import { paint } from '../draw/ideogram.js';
import * as fmt from '../format.js';
import { token, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, seg, select, fld, card, metrics, callout, empty, table, working, busy, onEnter, download } from '../ui/kit.js';

const HUMAN = {
  genomes: ['data:ncbi-grch38-chromosomes/NC_000022.11.fa', 'data:ncbi-grch38-chromosomes/NC_000021.9.fa'],
  regions: [{ contig: 0, begin: 30000000, end: 35000000 }, { contig: 1, begin: 30000000, end: 31800000 }],
  ideogram: { file: 'data:ncbi-ideogram/ideogram_9606_GCF_000001305.16_850_V1.tsv' },
};
// The individuals, their variants placed by hand on the chromosomes' own numbering (0-based); none is
// observed in a person: they are planted, and the reads drawn from them.
const SAMPLES = [
  { id: 'human', title: 'GRCh38 22q12 and 21q21: a translocation, a large deletion, a copy-neutral LOH', request: { ...HUMAN, variants: [
    { type: 'TRA', contig: 0, begin: 32200000, end: 32200000, zygosity: 'het', haplotype: 1, partner: 1, partner_position: 30800000 },
    { type: 'DEL', contig: 0, begin: 30500000, end: 31700000, zygosity: 'het', haplotype: 0 },
    { type: 'LOH', contig: 0, begin: 32400000, end: 32900000, zygosity: 'het', haplotype: 0 },
    { type: 'DEL', contig: 1, begin: 30300000, end: 30305000, zygosity: 'hom' },
    { type: 'DUP', contig: 1, begin: 31200000, end: 31230000, zygosity: 'het', haplotype: 1 },
    { type: 'INV', contig: 1, begin: 31500000, end: 31510000, zygosity: 'het', haplotype: 0 }] } },
  { id: 'trisomy', title: 'GRCh38 22q12 and 21q21: a trisomy of the chromosome-21 region', request: { ...HUMAN, variants: [{ type: 'ANEUPLOIDY', contig: 1, begin: 30000000, end: 31800000, copies: 3 }] } },
  { id: 'ecoli', title: 'Escherichia coli K-12: thirty variants drawn at random', request: { genome: 'ncbi:NC_000913.3', random: { deletions: 10, tandem_duplications: 4, dispersed_duplications: 2, inversions: 4, insertions: 6, gains: 1, losses: 1, losses_of_heterozygosity: 1 } } },
];
const COVERAGES = [15, 30, 60];
const KIND_COLOUR = { deletion: '--C', copy_loss: '--C', tandem_duplication: '--T', copy_gain: '--T', dispersed_duplication: '--accent-2', inversion: '--A', insertion: '--G', translocation: '--accent', breakend: '--ink-faint' };
const kindColour = (k) => token(KIND_COLOUR[k] || '--ink-faint');
// The engine names kinds as VCF does; the page by their words.
const KIND_NAMES = { DEL: 'deletion', DUP: 'tandem_duplication', 'DUP:TANDEM': 'tandem_duplication', 'DUP:DISPERSED': 'dispersed_duplication', INV: 'inversion', INS: 'insertion', TRA: 'translocation', BND: 'breakend', CNV: 'copy_number' };
const kindName = (k) => KIND_NAMES[k] || String(k).toLowerCase();

const kept = { sample: 'human', coverage: 30, seed: 51, r: null, error: null, running: false, call: null, window: null };
let view = null;

export function mount(section, ctx) {
  view = { ctx, section };
  const p = ctx.params;
  if (SAMPLES.some((s) => s.id === p.sample)) kept.sample = p.sample;
  if (/^\d{1,6}$/.test(p.seed || '')) kept.seed = Number(p.seed);
  view.controls = el('div', { class: 'controls' });
  view.stage = el('div', { class: 'stage' });
  section.append(
    el('p', { class: 'lede' }, 'A rearranged genome differs from its reference by new adjacencies: two ends of the reference, each with its orientation, joined where they were not. Each kind of rearrangement leaves a signature in the reads. Mates of a pair lie farther apart than the library allows, or face the wrong way. A read is split, its two parts mapped either side of the junction. The depth of fragments steps up or down. Where depth cannot decide, as between a deletion of one copy and a copy-neutral loss of heterozygosity, the balance of alleles at heterozygous sites can. Every call is then said in the three international notations.'),
    el('div', { class: 'bench' }, view.controls, view.stage));
  buildControls();
  draw();
  if (!kept.r && !kept.running) runSample();
}

export function show(ctx) { view.ctx = ctx; }

function buildControls() {
  view.runBtn = btn('Read, map and call', { icon: 'play', kind: 'primary', onClick: runSample });
  if (kept.running) busy(view.runBtn, true, 'Reading…');
  const again = btn('Another', { icon: 'refresh', title: 'the same variants, other reads: the next seed', onClick: () => { kept.seed++; runSample(); } });
  const c = el('div', { class: 'card' },
    el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'The individual'),
      fld('sample', select(SAMPLES.map((s) => [s.id, s.title]), kept.sample, (v) => { kept.sample = v; }),
        { hint: 'A real reference; the individual on it is simulated, its variants planted, so that every call has its truth.' }),
      fld('coverage', seg(COVERAGES.map((v) => [v, `${v}×`]), kept.coverage, (v) => { kept.coverage = Number(v); }, { block: true }), { hint: 'paired reads of 2 × 150 bases from fragments of 450 ± 45, with a bias for GC' })),
    el('div', { class: 'run-bar' }, view.runBtn, again));
  onEnter(c, runSample);
  clear(view.controls).append(c);
}

async function runSample() {
  if (kept.running) return;
  kept.running = true; kept.error = null;
  busy(view.runBtn, true, 'Reading…');
  view.ctx.setParams({ sample: kept.sample, seed: kept.seed });
  draw();
  const s = SAMPLES.find((x) => x.id === kept.sample);
  try {
    kept.r = await api.run('structural.sample', { ...s.request, coverage: kept.coverage, seed: kept.seed });
    for (const c of kept.r.calls || []) c.kind = kindName(c.kind);
    kept.call = (kept.r.calls || []).find((c) => c.passed)?.id ?? null;
  } catch (e) {
    kept.error = e.message;
  } finally {
    kept.running = false;
    busy(view.runBtn, false);
    draw();
  }
}

// ---- drawing ------------------------------------------------------------------------------------------
const W = () => Math.max(700, (view.stage?.clientWidth || 1000) - 34);
const contigLabel = (c) => {
  const n = /NC_0000(\d\d)/.exec(c.accession || '')?.[1];
  const chrom = n ? String(Number(n)).replace(/^23$/, 'X').replace(/^24$/, 'Y') : null;
  return chrom ? `chr${chrom} ${fmt.number(c.offset / 1e6, 4)}–${fmt.number((c.offset + c.length) / 1e6, 4)} Mb` : c.name;
};
const chromOf = (c) => { const n = /NC_0000(\d\d)/.exec(c.accession || '')?.[1]; return n ? String(Number(n)) : null; };

function draw() {
  const box = clear(view.stage);
  if (kept.running) box.append(el('div', { class: 'results-head' }, working('Planting the variants, drawing paired reads, mapping both mates, and calling from pairs, splits, depth and alleles…')));
  if (kept.error) box.append(callout(kept.error, 'bad'));
  const r = kept.r;
  if (!r) { if (!kept.running && !kept.error) box.append(card({ body: empty({ icon: 'structure', title: 'A rearranged genome', text: 'Choose an individual, then read it.' }) })); return; }
  const ev = r.evaluation || {}, o = ev.overall || {};
  const sum = Object.values(r.seconds || {}).reduce((a, b) => a + (Number(b) || 0), 0);
  box.append(el('div', { class: 'results-head' }, el('span', { class: 'badge warn' }, 'simulated'),
    el('span', { class: 'when' }, `${fmt.count(r.pairs_simulated)} pairs on ${r.contigs.map(contigLabel).join(' and ')} · ${fmt.seconds(sum)} on ${r.threads} threads`)));
  box.append(metrics([
    [fmt.percent(o.precision, 1), 'precision of the calls'], [fmt.percent(o.recall, 1), 'recall of the variants'],
    [fmt.percent(o.genotype, 1), 'genotypes right'],
    [fmt.percent(ev.copy_number?.fraction, 1), 'bins at the right copy number'],
    [String((r.calls || []).filter((c) => c.passed).length), 'calls that pass'],
    [`${fmt.number(r.library?.median, 4)} ± ${fmt.number(r.library?.sd, 3)}`, 'fragment length, estimated'],
  ]));
  box.append(genomeCard(r));
  box.append(el('div', { style: { display: 'grid', gridTemplateColumns: 'minmax(0, 1fr) minmax(0, 1fr)', gap: '16px' } }, callCard(r), chromosomesCard(r)));
  box.append(descriptionCard(r));
  box.append(evaluationCard(r));
}

function genomeCard(r) {
  const contigs = r.contigs.map((c, i) => ({ name: String(i), length: c.length, label: contigLabel(c) }));
  const local = (contig, pos) => ({ contig: String(contig), pos: pos - r.contigs[contig].offset });
  const d = r.depth, bin = d.bin || 1000;
  const copyPts = d.begin.map((b, i) => ({ ...local(d.contig[i], b + bin / 2), value: d.copies[i], colour: token('--ink-faint') }));
  const segs = (r.segments || []).map((s) => ({ contig: String(s.contig), start: s.begin - r.contigs[s.contig].offset, end: s.end - r.contigs[s.contig].offset, value: s.copies,
    colour: s.copies < 2 ? token('--C') : s.copies > 2 ? token('--T') : token('--accent'), width: 2.6, title: `${s.copies} copies (${fmt.number(s.estimate, 3)}) over ${s.bins} bins` }));
  const bafPts = r.baf.position.map((p, i) => ({ ...local(r.baf.contig[i], p), value: r.baf.baf[i], colour: r.baf.genotype[i] === 2 || r.baf.genotype[i] === '1/1' ? token('--hair-3') : token('--C') }));
  const allelic = (r.allelic || []).filter((a) => !/^heterozygous|^normal|^diploid/.test(a.state)).map((a) => ({ contig: String(a.contig), start: a.begin - r.contigs[a.contig].offset, end: a.end - r.contigs[a.contig].offset,
    colour: /neutral/.test(a.state) ? token('--high') : /hemi|loss|del/.test(a.state) ? token('--C') : token('--T'), opacity: 0.55, title: `${a.state.replace(/_/g, ' ')}: ${a.copies} copies, ${a.sites} sites` }));
  const callSpans = (r.calls || []).filter((c) => c.passed).map((c, k) => ({ contig: String(c.contig), start: c.begin - r.contigs[c.contig].offset, end: Math.max(c.end, c.begin + 1) - r.contigs[c.contig].offset,
    colour: kindColour(c.kind), title: `${c.kind.replace(/_/g, ' ')} · ${fmt.length(c.length || 0)} · ${c.genotype}`, key: c.id, row: c.kind === 'translocation' ? 0 : 1 + (k % 2), cls: c.id === kept.call ? 'lit' : '' }));
  const truthSpans = (r.truth?.variants || []).map((v) => ({ contig: String(v.contig), start: v.begin - r.contigs[v.contig].offset, end: Math.max(v.end, v.begin + 1) - r.contigs[v.contig].offset, colour: token('--hair-3'), title: `truth: ${v.type} ${v.zygosity}`, opacity: 0.9 }));
  const arcs = (r.calls || []).filter((c) => c.passed).map((c) => (c.kind === 'translocation' || c.kind === 'dispersed_duplication')
    ? { a: local(c.contig, c.begin), b: local(c.partner, c.partner_position), colour: kindColour(c.kind), width: c.id === kept.call ? 2.6 : 1.4, title: `${c.kind}: ${c.hgvs || ''}` }
    : { a: local(c.contig, c.begin), b: local(c.contig, c.end), colour: kindColour(c.kind), width: c.id === kept.call ? 2.4 : 1.2, opacity: 0.7 });
  const ymax = Math.max(4, ...segs.map((s) => s.value + 1));
  const svg = genomeTracksSvg({ contigs, width: W(), arcs, arcHeight: 64, tracks: [
    { kind: 'points', label: 'copies', height: 120, y: [0, ymax], lines: [2], ticks: [0, 1, 2, 3, 4].filter((v) => v <= ymax), opacity: 0.4, points: copyPts, segments: segs },
    { kind: 'points', label: 'B-allele fraction', height: 90, y: [0, 1], lines: [0.5], ticks: [0, 0.5, 1], opacity: 0.45, points: bafPts },
    { kind: 'spans', label: 'allelic state', height: 18, spans: allelic },
    { kind: 'spans', label: 'calls', height: 46, spans: callSpans },
    { kind: 'spans', label: 'truth', height: 16, spans: truthSpans },
  ] });
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = svg;
  fig.addEventListener('click', (e) => {
    const t = e.target.closest('[data-key]');
    if (!t || t.dataset.key === '') return;
    kept.call = Number.isNaN(Number(t.dataset.key)) ? t.dataset.key : Number(t.dataset.key);
    draw();
  });
  const key = el('div', { class: 'paints' }, ...[['deletion', 'deletion, loss'], ['tandem_duplication', 'duplication, gain'], ['inversion', 'inversion'], ['insertion', 'insertion'], ['translocation', 'translocation'], ['dispersed_duplication', 'dispersed duplication']]
    .map(([k, l]) => el('span', {}, el('span', { class: 'sw', style: { background: kindColour(k) } }), l)),
    el('span', {}, el('span', { class: 'sw', style: { background: token('--high') } }), 'copy-neutral LOH'));
  return card({ title: 'The genome read three ways', sub: 'copies per bin of fragments, corrected for GC and mappability and segmented to whole numbers; the B-allele fraction at the heterozygous sites; and the calls, joined by arcs to their other ends. Click a call.', pad: true, body: [fig, key] });
}

const callById = (r, id) => (r.calls || []).find((c) => c.id === id);

function callCard(r) {
  const c = callById(r, kept.call);
  if (!c) return card({ title: 'A call', pad: true, body: empty({ icon: 'structure', title: 'Choose a call', text: 'Click one on the calls\' track.' }) });
  const t = Number.isInteger(c.truth) ? r.truth?.variants?.[c.truth] : c.truth;
  const facts = el('div', { class: 'forms' },
    form('kind', `${c.kind.replace(/_/g, ' ')}${c.inverted ? ', inverted' : ''}`, c.precise ? 'to the base, by split reads' : `imprecise: within ${fmt.count((c.ci_end?.[1] ?? 0) - (c.ci_begin?.[0] ?? 0))} bases`),
    form('genotype', `${c.genotype} · GQ ${c.gq}`, `${Math.round(c.alternate_support)} reads for it, ${Math.round(c.reference_support)} against; ${c.support?.pairs ?? 0} pairs, ${c.support?.splits ?? 0} splits, ${c.support?.clips ?? 0} clips`),
    form('size', c.length ? fmt.length(c.length) : '–', Number.isFinite(c.copies) ? `${fmt.number(c.copies, 3)} copies by depth` : ''),
    form('truth', t ? `${kindName(t.type).replace(/_/g, ' ')}, ${t.zygosity === 'het' ? 'heterozygous' : 'homozygous'}` : 'none: a false call', t && Number.isFinite(c.breakpoint_error) ? `breakpoints ${c.breakpoint_error === 0 ? 'exact' : `${c.breakpoint_error} bases off`}` : ''));
  const names = el('div', { class: 'forms' },
    form('HGVS', c.hgvs || '–', 'the sequence\'s description', true),
    form('ISCN', c.iscn || '–', 'the chromosomes\' description'),
    form('VCF', (c.vcf || []).map((l) => l.split('\t').slice(0, 5).join('  ')).join('\n') || '–', 'the record, first five columns'));
  if (c.filters?.length) names.append(form('filters', c.filters.join(', '), 'why it does not pass'));
  return card({ title: `The evidence for call ${c.id}`, sub: 'what was found, what it is, and what the reads say for and against it', pad: true, body: [facts, evidenceFigure(r, c), names] });
}

function form(k, v, note, lead = false) {
  return el('div', { class: `form${lead ? ' lead' : ''}` }, el('div', { class: 'k' }, k), el('div', { class: 'v', style: { whiteSpace: 'pre-wrap' } }, v || '–'), el('div', { class: 'note' }, note));
}

// Two windows about the call's two breakends on a broken axis: discordant pairs as mates joined by an
// arc across the break, split reads as their parts joined by a dashed line.
function evidenceFigure(r, c) {
  const ends = [{ contig: c.contig, at: c.begin }, (c.kind === 'translocation' || c.kind === 'dispersed_duplication') ? { contig: c.partner, at: c.partner_position } : { contig: c.contig, at: c.end }];
  const half = 900, w = 260, gap = 26, H = 210, top = 18;
  const X = (k, contig, pos) => (contig === ends[k].contig && Math.abs(pos - ends[k].at) <= half ? 8 + k * (w + gap) + ((pos - ends[k].at + half) / (2 * half)) * w : null);
  const place = (contig, pos) => { for (const k of [0, 1]) { const x = X(k, contig, pos); if (x !== null) return x; } return null; };
  const out = [];
  for (const k of [0, 1]) {
    const x0 = 8 + k * (w + gap);
    out.push(`<rect x="${x0}" y="${top}" width="${w}" height="${H - top - 24}" fill="${token('--band')}" />`);
    out.push(`<line x1="${x0 + w / 2}" y1="${top}" x2="${x0 + w / 2}" y2="${H - 24}" stroke="${kindColour(c.kind)}" stroke-dasharray="3 3" />`);
    out.push(`<text class="tick-label" x="${x0 + w / 2}" y="${H - 8}" text-anchor="middle">${esc(`${contigLabel(r.contigs[ends[k].contig]).split(' ')[0]}:${fmt.count(ends[k].at + 1)}`)}</text>`);
  }
  out.push(`<path d="M${8 + w + gap / 2 - 4},${top + 6} l8,-6 M${8 + w + gap / 2 - 4},${top + 14} l8,-6" stroke="${token('--ink-faint')}" />`);
  let row = 0;
  const rowY = (i) => top + 8 + (i % 22) * 7.4;
  for (const pr of (c.evidence?.pairs || []).slice(0, 22)) {
    const [m1, m2] = pr;
    const a = place(m1[0], (m1[1] + m1[2]) / 2), b = place(m2[0], (m2[1] + m2[2]) / 2);
    const y = rowY(row++);
    for (const [m, x] of [[m1, a], [m2, b]]) if (x !== null) out.push(`<path d="M${x - 6},${y} h12 l${m[3] ? -16 : 4},0" stroke="${token('--C')}" stroke-width="2.2" stroke-linecap="round" />`);
    if (a !== null && b !== null) out.push(`<path d="M${a},${y} Q${(a + b) / 2},${y - 10} ${b},${y}" fill="none" stroke="${token('--C')}" stroke-opacity="0.35" />`);
  }
  for (const sp of (c.evidence?.splits || []).slice(0, 22 - Math.min(22, row))) {
    const y = rowY(row++);
    const xs = sp.parts.map((p) => [place(p[0], p[1]), place(p[0], p[2])]).filter(([u, v]) => u !== null && v !== null);
    xs.forEach(([u, v]) => out.push(`<line x1="${u}" y1="${y}" x2="${Math.max(u + 2, v)}" y2="${y}" stroke="${token('--accent')}" stroke-width="2.6" stroke-linecap="round" />`));
    if (xs.length > 1) out.push(`<line x1="${xs[0][1]}" y1="${y}" x2="${xs[1][0]}" y2="${y}" stroke="${token('--accent')}" stroke-dasharray="2 3" stroke-opacity="0.6" />`);
  }
  const width = 16 + 2 * w + gap;
  const fig = el('div', { style: { overflowX: 'auto', margin: '10px 0' } });
  fig.innerHTML = `<svg width="${width}" height="${H}" viewBox="0 0 ${width} ${H}">${out.join('')}</svg>`;
  fig.append(el('div', { class: 'hint' }, 'Discordant pairs in blue, their mates joined across the break; split reads in mulberry, their parts joined by a dashed line; each window ', fmt.length(2 * half), ' about a breakend.'));
  return fig;
}

// The individual's chromosomes as they are, pieces of the reference in order and orientation.
function chromosomesCard(r) {
  const chs = r.truth?.chromosomes || [];
  if (!chs.length) return card({ title: 'Its chromosomes', pad: true, body: empty({ title: 'No truth to draw' }) });
  const longest = Math.max(...chs.map((c) => c.length));
  const w = 560, rowH = 30, top = 10;
  const out = [];
  chs.forEach((ch, i) => {
    const y = top + i * rowH;
    const nm = ch.name.replace(/NC_0000(\d\d)\.\d+/, (_, n) => `chr${Number(n)}`).replace(/#(\d+)/, ' · haplotype $1');
    out.push(`<text class="tick-label" x="0" y="${y + 13}">${esc(nm)}</text>`);
    let x = 150;
    for (const p of ch.pieces) {
      const len = p[0] === 'novel' ? p[1] : p[2] - p[1];
      const pw = Math.max(1.5, (len / longest) * (w - 158));
      if (p[0] === 'novel') {
        out.push(`<rect x="${x}" y="${y + 3}" width="${pw}" height="14" fill="${token('--G')}"><title>${esc(`inserted: ${fmt.length(len)}`)}</title></rect>`);
      } else {
        const chrom = chromOf(r.contigs[p[0]]);
        const fill = chrom ? paint(chrom) : CATEGORIES[p[0] % CATEGORIES.length];
        out.push(`<rect x="${x}" y="${y + 3}" width="${pw}" height="14" fill="${fill}" stroke="${token('--surface')}" stroke-width="0.6"><title>${esc(`${contigLabel(r.contigs[p[0]]).split(' ')[0]}:${fmt.count(p[1] + 1)}–${fmt.count(p[2])}${p[3] ? ', reversed' : ''}`)}</title></rect>`);
        if (p[3] && pw > 10) out.push(`<path d="M${x + pw - 3},${y + 10} h${-(pw - 6)} m4,-3 l-4,3 l4,3" fill="none" stroke="${token('--surface')}" stroke-width="1.2" />`);
      }
      x += pw;
    }
  });
  const fig = el('div', { style: { overflowX: 'auto' } });
  fig.innerHTML = `<svg width="${w}" height="${top + chs.length * rowH}" viewBox="0 0 ${w} ${top + chs.length * rowH}">${out.join('')}</svg>`;
  return card({ title: 'Its chromosomes, as planted', sub: 'each haplotype\'s chromosomes as pieces of the reference, in their order and orientation, painted by the chromosome each came from; an arrow marks a piece reversed', pad: true, body: fig });
}

function descriptionCard(r) {
  const d = r.descriptions || {};
  const body = [el('div', { class: 'forms' },
    form('ISCN · arr', d.arr || '–', 'the copy number of the sample, as an array would report it', true),
    form('ISCN · seq', d.seq || '–', 'the junctions, as sequencing reports them'))];
  if (d.seq_detailed) body.push(el('div', { class: 'long-form', style: { marginTop: '10px' } }, d.seq_detailed));
  const vcfBtn = btn('VCF', { icon: 'download', size: 'sm', onClick: () => download('structural-calls.vcf', r.vcf || '', 'text/vcf') });
  return card({ title: 'The sample, said in the notations', sub: d.ideogram ? 'with GRCh38\'s bands from NCBI\'s ideogram' : 'no ideogram: bands are named only on human chromosomes', actions: [vcfBtn], pad: true, body });
}

function evaluationCard(r) {
  const ev = r.evaluation || {};
  const rows = Object.entries(ev.by_kind || {}).map(([k, t]) => ({ kind: kindName(k).replace(/_/g, ' '), tp: t.tp, fp: t.fp, fn: t.fn, precision: fmt.percent(t.precision, 1), recall: fmt.percent(t.recall, 1) }));
  const errs = ev.breakpoint_errors || [];
  const hist = new Map();
  for (const e of errs) { const b = Math.min(20, Math.round(e)); hist.set(b, (hist.get(b) || 0) + 1); }
  const bars = [...hist].map(([b, n]) => ({ x0: b - 0.4, x1: b + 0.4, y: n, colour: token('--accent-2') }));
  const fig = el('div');
  fig.innerHTML = chart({ width: 360, height: 180, x: [-1, 21], y: [0, Math.max(1, ...hist.values()) * 1.15], xlabel: 'breakpoint error (bases)', ylabel: 'calls', bars });
  const lib = r.library?.histogram || [];
  const libFig = el('div');
  if (lib.length) libFig.innerHTML = chart({ width: 360, height: 180, x: [lib[0][0], lib[lib.length - 1][0]], y: [0, Math.max(...lib.map((x) => x[1])) * 1.1], xlabel: 'fragment length', ylabel: 'pairs', series: [{ points: lib, colour: token('--C'), area: true }] });
  return card({ title: 'Against the truth', sub: `matched as Truvari matches: the same kind, overlapping by half, ends within 500 bases; copy number right in ${fmt.percent(ev.copy_number?.fraction, 1)} of bins; LOH recall ${fmt.percent(ev.loh?.recall, 1)}`, pad: true,
    body: el('div', { style: { display: 'grid', gridTemplateColumns: 'minmax(0, 1fr) auto auto', gap: '16px', alignItems: 'start' } },
      table([{ key: 'kind', label: 'kind' }, { key: 'tp', label: 'right', num: true }, { key: 'fp', label: 'false', num: true }, { key: 'fn', label: 'missed', num: true }, { key: 'precision', label: 'precision' }, { key: 'recall', label: 'recall' }], rows), fig, libFig) });
}
