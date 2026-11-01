// The Describe view: the international notations of genomes that depart from their reference, each
// read by its grammar (morphe's BNF and Earley's parser) and given its meaning by the domain (ADR 0064).
//
//   Sequence variants · HGVS   a description parsed, its derivation drawn over its text and read in
//                              words; on a real reference, a change written in every form (HGVS by
//                              the 3′ rule, VCF's left form, SPDI's hull) over the interval of its
//                              placements, and projected through each transcript to its protein.
//   Chromosomes · ISCN         a karyotype, a microarray or a sequencing result parsed and
//                              interpreted: every chromosome painted as segments of the normal ones,
//                              the copy number at every band, what is lost and gained from where to
//                              where, which pieces are joined, and the detailed system.
//
// A rejected text keeps its viable prefix underlined and the place of the fault marked with what
// could have come: the chart of the parse says exactly how far the text is right.

import { el, clear, esc, hover } from '../dom.js';
import * as api from '../api.js';
import { icon } from '../ui/icons.js';
import { btn, seg, select, tabs, card, metrics, callout, empty, working, onEnter } from '../ui/kit.js';
import { derivationSvg } from '../draw/derivation.js';
import { karyogramSvg, imbalanceSvg, paint, CHROMOSOMES } from '../draw/ideogram.js';
import * as fmt from '../format.js';

// Real descriptions, each checked against Mutalyzer or NCBI's Variation Services by the domain's
// fixtures (domains/genomics/tests/fixtures/); "genome" is the data library's file of its reference.
const HBB = 'data:ncbi-human-genes/NC_000011.10_5222001-5232000.gb', CFTR = 'data:ncbi-human-genes/NC_000007.14_117479001-117670000.gb';
const BRCA1 = 'data:ncbi-human-genes/NC_000017.11_43043001-43172000.gb', MITO = 'ncbi:NC_012920.1';
export const HGVS_EXAMPLES = [
  { text: 'NC_000011.10:g.5227002T>A', name: 'HbS', genome: HBB, first: 5222001, about: 'the sickle-cell variant of HBB: c.20A>T, p.(Glu7Val)' },
  { text: 'NC_000007.14:g.117559592_117559594del', name: 'F508del', genome: CFTR, first: 117479001, about: "CFTR's commonest variant, a deletion in a repeat: three placements, three normal forms" },
  { text: 'NC_000017.11:g.43124030_43124031del', name: '185delAG', genome: BRCA1, first: 43043001, about: 'a frameshift of BRCA1: c.68_69del, p.(Glu23ValfsTer17)' },
  { text: 'NC_000017.11:g.43057065dup', name: '5382insC', genome: BRCA1, first: 43043001, about: 'a duplication of BRCA1 that the old name called an insertion: c.5266dup' },
  { text: 'NC_000011.10:g.5226929C>T', name: 'IVS-I-1', genome: HBB, first: 5222001, about: "the first base of HBB's first intron: c.92+1G>A, a splice donor" },
  { text: 'NC_000011.10:g.5225598T>A', name: 'stop lost', genome: HBB, first: 5222001, about: "HBB's stop codon read through: p.(Ter148TyrextTer21)" },
  { text: 'NC_000007.14:g.117548634_117548635del', name: 'CFTR 5T', genome: CFTR, first: 117479001, about: 'a deletion in an intron near a splice site: c.1210-7_1210-6del' },
  { text: 'NC_012920.1:m.8993T>G', name: 'NARP', genome: MITO, first: 1, about: 'a mitochondrial variant of ATP6, read with the vertebrate mitochondrial code: p.(Leu156Arg)' },
  { text: 'NC_012920.1:m.14484T>C', name: 'LHON', genome: MITO, first: 1, about: 'on the minus strand: ND6 p.(Met64Val)' },
  { text: 'NC_012920.1:m.3243A>G', name: 'MELAS', genome: MITO, first: 1, about: 'in a tRNA gene, MT-TL1: a variant of a non-coding transcript' },
];
// Descriptions to read by the grammar alone (HGVS 21.1's own examples, and a fault).
const HGVS_GRAMMAR = [
  ['NC_000023.11(NM_004006.2):c.(4071+1_4072-1)_(5154+1_5155-1)del', 'uncertain breakpoints, in introns'],
  ['NP_003997.2:p.[(Asn158Asp)(;)(Asn158Ile)]^[(Asn158Val)]', 'alleles and alternatives, predicted'],
  ['NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]', "a translocation's derivative chromosome"],
  ['NM_004006.2:c.4375dex', 'a fault: where the text stops being HGVS'],
];

const kept = {
  tab: 'iscn',
  hgvs: { text: '', parse: null, describe: null, error: null, busy: false, genome: null },
  iscn: { text: '47,XX,+21', edition: '2020', result: null, error: null, busy: false, line: 0, focus: null, painted: true },
  examples: null,
  ideogram: null,
};
let view = null;
let ticket = 0;

export function mount(section, ctx) {
  view = { ctx, section };
  const p = ctx.params;
  if (p.tab === 'hgvs' || p.tab === 'iscn') kept.tab = p.tab;
  if (typeof p.iscn === 'string' && p.iscn) kept.iscn.text = p.iscn;
  if (typeof p.hgvs === 'string' && p.hgvs) kept.hgvs.text = p.hgvs;
  view.tabs = tabs([['hgvs', 'Sequence variants · HGVS'], ['iscn', 'Chromosomes · ISCN']], kept.tab, (v) => { kept.tab = v; remember(); draw(); });
  view.body = el('div');
  section.append(
    el('p', { class: 'lede' }, 'A genome that departs from its reference is described in two international notations, both languages with grammars. HGVS writes a change of the sequence against a named reference: where it lies, what was deleted, inserted, duplicated or inverted, and what it does to the protein. ISCN writes a change of the chromosomes as a cytogeneticist sees them: how many there are, which bands are lost or gained, and which pieces of which chromosomes are joined, from where to where. Each text here is read by its grammar; its derivation is drawn over it, and its meaning is computed: the forms of a variant, and every chromosome of a karyotype rebuilt from the normal ones and painted.'),
    view.tabs, view.body);
  draw();
}

export function show(ctx) {
  view.ctx = ctx;
  const p = ctx.params;
  let changed = false;
  if ((p.tab === 'hgvs' || p.tab === 'iscn') && p.tab !== kept.tab) { kept.tab = p.tab; view.tabs.set(p.tab); changed = true; }
  if (typeof p.iscn === 'string' && p.iscn && p.iscn !== kept.iscn.text) { kept.iscn.text = p.iscn; kept.iscn.result = null; changed = true; }
  if (typeof p.hgvs === 'string' && p.hgvs && p.hgvs !== kept.hgvs.text) { kept.hgvs.text = p.hgvs; kept.hgvs.parse = null; changed = true; }
  if (changed) draw();
}

function remember() {
  if (view.ctx.view !== 'describe') return;
  view.ctx.setParams({ tab: kept.tab, iscn: kept.tab === 'iscn' ? kept.iscn.text : null, hgvs: kept.tab === 'hgvs' ? kept.hgvs.text || null : null });
}

function draw() {
  clear(view.body);
  if (kept.tab === 'iscn') drawIscn();
  else drawHgvs();
}

// ---- shared: a notation's line and its derivation ---------------------------------------------------------
function notationLine(value, placeholder, onRun, label = 'Read') {
  const input = el('input', { type: 'text', value, placeholder, spellcheck: false, autocomplete: 'off' });
  const run = btn(label, { icon: 'play', kind: 'primary', onClick: () => onRun(input.value.trim()) });
  input.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); onRun(input.value.trim()); } });
  const box = el('div', { class: 'notation' }, input, run);
  box.input = input;
  box.run = run;
  return box;
}

function derivationFigure(text, r) {
  const fig = el('div', { class: 'figure' });
  if (r.accepted && r.tree) {
    fig.append(el('div', { class: 'eyebrow' }, icon('grammar', 'sm'), 'Its derivation, read by the grammar'));
    const box = el('div', { style: { overflowX: 'auto' } });
    box.innerHTML = derivationSvg({ text, tree: r.tree });
    // A rule under the pointer: its name, and the text it derives.
    box.addEventListener('mousemove', (e) => {
      const g = e.target.closest('.dnode');
      if (!g) { hover(null); return; }
      const [name] = g.dataset.node.split(':');
      const b = Number(g.dataset.begin), en = Number(g.dataset.end);
      hover(`<b>${esc(name)}</b><br><code>${esc(text.slice(b, en))}</code>`, e);
    });
    box.addEventListener('mouseleave', () => hover(null));
    fig.append(box);
  } else {
    fig.append(el('div', { class: 'eyebrow' }, icon('warn', 'sm'), 'Not in the language'));
    const box = el('div', { style: { overflowX: 'auto' } });
    box.innerHTML = derivationSvg({ text, furthest: r.furthest ?? 0, expected: r.expected || [] });
    fig.append(box, el('p', { class: 'reading', style: { fontSize: '14.5px' } }, r.error || 'The text is not in the language.'));
  }
  return fig;
}

function diagnosticList(ds) {
  if (!ds || !ds.length) return null;
  return el('div', { style: { marginTop: '8px' } }, ...ds.map((d) => el('div', { class: 'diag' },
    el('span', { class: `badge ${/error/.test(d.level || d.severity) ? 'bad' : /warn/.test(d.level || d.severity) ? 'warn' : ''}` }, d.level || d.severity || 'note'),
    el('span', {}, d.message))));
}

// ---- ISCN ---------------------------------------------------------------------------------------------------
async function loadIscnBasics() {
  const jobs = [];
  if (!kept.examples) jobs.push(api.run('karyotype.examples', {}).then((r) => { kept.examples = r.examples; }).catch(() => { kept.examples = []; }));
  if (!kept.ideogram) jobs.push(api.run('karyotype.ideogram', { assembly: 'GRCh38', level: 850 }).then((r) => {
    kept.ideogram = { assembly: r.assembly, level: r.level, chromosomes: r.chromosomes.map((c) => ({ name: c.name, length: c.length, centromere: c.centromere,
      bands: c.bands.map(([name, start, end, stain]) => ({ name, start, end, stain })) })) };
  }).catch((e) => { kept.ideogramError = e.message; }));
  await Promise.all(jobs);
}

function drawIscn() {
  const s = kept.iscn;
  const line = notationLine(s.text, '46,XY,t(9;22)(q34;q11.2)', (t) => { if (t) { s.text = t; interpret(); } }, 'Interpret');
  view.iscnLine = line;
  const edition = seg([['2020', 'ISCN 2020'], ['2016', '2016']], s.edition, (v) => { s.edition = v; if (s.result) interpret(); });
  const ex = el('div', { class: 'examples' });
  const fillExamples = () => {
    clear(ex);
    for (const e of kept.examples || []) {
      ex.append(el('button', { class: 'chip', title: `${e.kind} — ${e.source}`, onclick: () => { s.text = e.text; line.input.value = e.text; interpret(); } },
        e.text.length > 44 ? `${e.text.slice(0, 42)}…` : e.text));
    }
  };
  view.iscnStage = el('div', { class: 'stack', style: { gap: '16px', marginTop: '16px' } });
  view.body.append(card({ pad: true, body: [
    el('div', { class: 'ctl-row', style: { justifyContent: 'space-between', alignItems: 'center', marginBottom: '8px' } },
      el('div', { class: 'eyebrow', style: { margin: 0 } }, 'A karyotype, an array or a sequencing result'), edition),
    line, ex,
  ] }), view.iscnStage);
  if (kept.examples) fillExamples();
  loadIscnBasics().then(() => { fillExamples(); if (!s.result && !s.busy) interpret(); else drawIscnResult(); });
}

async function interpret() {
  const s = kept.iscn;
  const t = ++ticket;
  s.busy = true; s.error = null;
  remember();
  if (view.iscnLine) { view.iscnLine.input.value = s.text; view.iscnLine.input.classList.remove('bad'); }
  clear(view.iscnStage).append(working('Reading the description by ISCN\'s grammar and rebuilding its chromosomes…'));
  try {
    const r = await api.run('karyotype.interpret', { text: s.text, edition: s.edition, assembly: 'GRCh38', level: 850 });
    if (t !== ticket) return;
    s.result = r; s.line = 0; s.focus = null;
  } catch (e) {
    if (t !== ticket) return;
    s.error = e.message; s.result = null;
  } finally {
    if (t === ticket) s.busy = false;
  }
  drawIscnResult();
}

// A chromosome's copies, grouped under the chromosome whose centromere they carry; markers and
// pieces of unknown origin apart.
function copiesOf(lineR) {
  const copies = {};
  for (const n of CHROMOSOMES) copies[n] = [];
  const odd = [];
  for (const ch of lineR.chromosomes || []) {
    const segments = (ch.segments || []).filter((x) => x.kind === 'reference').map((x) => ({ chrom: x.chrom, start: x.start, end: x.end, reverse: !!x.reverse }));
    const entry = { key: ch.key, label: ch.label, abnormal: !!ch.abnormal, segments, long: ch.long_form, marker: ch.marker };
    if (ch.centromere_of && copies[ch.centromere_of]) copies[ch.centromere_of].push(entry);
    else odd.push(entry);
  }
  return { copies, odd };
}

function drawIscnResult() {
  const s = kept.iscn;
  const box = clear(view.iscnStage);
  if (s.error) { box.append(callout(s.error, 'bad')); return; }
  const r = s.result;
  if (!r) return;
  if (view.iscnLine) view.iscnLine.input.classList.toggle('bad', !r.accepted);
  box.append(derivationFigure(s.text, r));
  if (!r.accepted) return;
  const head = el('div', { class: 'figure' },
    el('div', { class: 'eyebrow' }, 'Written canonically'),
    el('div', { class: 'reading' }, el('code', {}, r.canonical || s.text)),
    diagnosticList(r.diagnostics));
  box.append(head);
  if (!kept.ideogram) { box.append(callout(kept.ideogramError || 'The ideogram of GRCh38 is not in the cache: python data/acquisition/fetch.py ncbi-ideogram', 'warn')); return; }
  const lines = r.lines || [];
  const techniques = (r.techniques || []).filter((x) => (x.changes || []).length || x.line);
  if (!lines.length && !techniques.length) return;
  if (lines.length > 1) {
    const total = lines.reduce((a, L) => a + (L.cells || 0), 0);
    box.append(el('div', { class: 'ctl-row', style: { alignItems: 'center', gap: '10px' } },
      el('span', { class: 'eyebrow', style: { margin: 0 } }, `${lines.length} cell lines`),
      seg(lines.map((L, i) => [i, `${i + 1}${L.cells ? ` · ${L.cells} of ${total} cells` : ''}`]), s.line, (v) => { s.line = Number(v); s.focus = null; drawIscnResult(); })));
  }
  if (lines.length) box.append(...lineCards(lines[Math.min(s.line, lines.length - 1)], r));
  for (const tq of techniques) box.append(...techniqueCards(tq));
}

function lineCards(L, r) {
  const s = kept.iscn;
  const out = [];
  const { copies, odd } = copiesOf(L);
  const changes = L.changes || [];
  const losses = changes.filter((c) => c.copies < c.normal), gains = changes.filter((c) => c.copies > c.normal);
  const fusions = L.lgf?.fusion_labels || [];
  const abnormal = (L.chromosomes || []).filter((c) => c.abnormal);
  out.push(metrics([
    [String(L.counted ?? (L.chromosomes || []).length), 'chromosomes counted'],
    [L.balanced ? 'balanced' : 'unbalanced', L.balanced ? 'no band lost or gained' : 'bands lost or gained', L.balanced ? 'good' : ''],
    [String(abnormal.length), abnormal.length === 1 ? 'abnormal chromosome' : 'abnormal chromosomes'],
    [String(fusions.length), fusions.length === 1 ? 'new junction' : 'new junctions'],
    [fmt.length(losses.reduce((a, c) => a + (c.end - c.start), 0)), 'lost'],
    [fmt.length(gains.reduce((a, c) => a + (c.end - c.start), 0)), 'gained'],
  ]));
  // The karyogram, painted.
  const kary = el('div', { class: `karyo${s.focus ? ' focus' : ''}`, style: { overflowX: 'auto' } });
  const width = Math.max(760, Math.min(1080, (view.section.clientWidth || 1000) - 60));
  kary.innerHTML = karyogramSvg({ ideogram: kept.ideogram, copies, width, painted: s.painted, highlight: s.focus });
  kary.addEventListener('click', (e) => {
    const g = e.target.closest('.rod');
    const k = g?.dataset.key || null;
    s.focus = s.focus === k ? null : k;
    drawIscnResult();
  });
  const involved = [...new Set(abnormal.flatMap((c) => (c.segments || []).filter((x) => x.kind === 'reference').map((x) => x.chrom)))];
  const legend = el('div', { class: 'paints' }, ...involved.map((n) => el('span', {}, el('span', { class: 'sw', style: { background: paint(n) } }), `chromosome ${n}`)));
  out.push(card({ title: 'The karyotype, rebuilt and painted', sub: 'every chromosome as segments of the normal ones, coloured by the chromosome each piece came from; a rose tick marks a junction',
    actions: [seg([[true, 'painted'], [false, 'banded']], s.painted, (v) => { s.painted = v === true || v === 'true'; drawIscnResult(); })],
    pad: true, body: [kary, involved.length ? legend : null, odd.length ? el('p', { class: 'hint' }, `Not drawn: ${odd.map((o) => o.label).join(', ')} (no centromere of known origin).`) : null] }));
  // Gains and losses on the reference.
  if (changes.length) {
    const imb = el('div', { style: { overflowX: 'auto' } });
    imb.innerHTML = imbalanceSvg({ ideogram: kept.ideogram, segments: changes.map((c) => ({ chrom: c.chrom, start: c.start, end: c.end, copies: c.copies, normal: c.normal, label: `${c.from_band}–${c.to_band}` })), width, sex: L.male ? 'XY' : 'XX' });
    out.push(card({ title: 'What is lost and gained', sub: 'losses to the left of each chromosome in blue, gains to the right in rose, as many bars as copies away from the normal', pad: true, body: imb }));
  }
  // In words: from where to where.
  const list = el('div', { class: 'change-list' });
  for (const c of losses) list.append(changeRow('loss', `${c.normal - c.copies === 1 ? 'one copy lost' : `${c.normal - c.copies} copies lost`}`, c));
  for (const c of gains) list.append(changeRow('gain', `${c.copies - c.normal === 1 ? 'one copy gained' : `${c.copies - c.normal} copies gained`}`, c));
  for (const f of fusions) list.append(el('div', { class: 'change' }, el('span', { class: 'what fusion' }, 'joined'), el('span', {}, el('code', {}, f))));
  for (const ch of abnormal) if (ch.long_form) list.append(el('div', { class: 'change' }, el('span', { class: 'what' }, esc(ch.label)), el('span', {}, el('code', {}, ch.long_form))));
  if (list.children.length) out.push(card({ title: 'From where to where', sub: 'the bands and GRCh38 bases of every loss and gain, every new junction, and each abnormal chromosome in ISCN\'s detailed form', pad: true, body: list }));
  if (L.detailed) out.push(card({ title: 'The detailed system', sub: 'the whole karyotype with every abnormal chromosome written pter to qter', pad: true, body: el('div', { class: 'long-form' }, L.detailed) }));
  const lineDiag = diagnosticList(L.diagnostics);
  if (lineDiag) out.push(card({ title: 'Remarks on its sense', pad: true, body: lineDiag }));
  return out;
}

function changeRow(kind, what, c) {
  const span = c.end - c.start;
  const uncertain = c.outer && c.inner && (c.outer[0] !== c.inner[0] || c.outer[1] !== c.inner[1]);
  return el('div', { class: 'change' }, el('span', { class: `what ${kind}` }, what),
    el('span', {}, el('b', {}, c.from_band === c.to_band ? c.from_band : `${c.from_band} → ${c.to_band}`), ' · ',
      el('code', {}, `chr${c.chrom}:${fmt.count(c.start + 1)}–${fmt.count(c.end)}`), ` (${fmt.length(span)})`,
      uncertain ? el('span', { class: 'hint' }, ` · breakpoints known to within ${fmt.length(Math.max(c.inner[0] - c.outer[0], c.outer[1] - c.inner[1]))}`) : null));
}

function techniqueCards(tq) {
  const out = [];
  const changes = tq.changes || [];
  const width = Math.max(760, Math.min(1080, (view.section.clientWidth || 1000) - 60));
  if (changes.length) {
    const imb = el('div', { style: { overflowX: 'auto' } });
    imb.innerHTML = imbalanceSvg({ ideogram: kept.ideogram, segments: changes.map((c) => ({ chrom: c.chrom, start: c.start, end: c.end, copies: c.copies, normal: c.normal })), width });
    const list = el('div', { class: 'change-list' });
    for (const c of changes) list.append(changeRow(c.copies < c.normal ? 'loss' : 'gain', c.copies < c.normal ? `${c.copies} ${c.copies === 1 ? 'copy' : 'copies'} left` : `${c.copies} copies`, c));
    out.push(card({ title: `${tq.technique}[${tq.build || tq.assembly}]`, sub: 'the result\'s copy number on the reference', pad: true, body: [imb, list, diagnosticList(tq.diagnostics)] }));
  }
  if (tq.line) out.push(...lineCards(tq.line, null));
  return out;
}

// ---- HGVS ---------------------------------------------------------------------------------------------------
function drawHgvs() {
  const s = kept.hgvs;
  const line = notationLine(s.text, 'NC_012920.1:m.3243A>G', (t) => { if (t) { s.text = t; readHgvs(); } }, 'Read');
  view.hgvsLine = line;
  const pick = (text, genome = null) => { s.text = text; s.genome = genome; line.input.value = text; readHgvs(); };
  const ex = el('div', { class: 'examples' }, ...HGVS_EXAMPLES.map((e) => el('button', { class: 'chip', title: e.about, onclick: () => pick(e.text, e) }, e.name, el('i', {}, e.text.replace(/^.*?:/, '')))));
  const gx = el('div', { class: 'examples' }, el('span', { class: 'hint', style: { alignSelf: 'center' } }, 'the grammar alone:'),
    ...HGVS_GRAMMAR.map(([t, about]) => el('button', { class: 'chip', title: about, onclick: () => pick(t) }, t.length > 46 ? `${t.slice(0, 44)}…` : t)));
  view.hgvsStage = el('div', { class: 'stack', style: { gap: '16px', marginTop: '16px' } });
  view.body.append(card({ pad: true, body: [el('div', { class: 'eyebrow' }, 'A sequence variant, as HGVS writes it'), line, ex, gx] }), view.hgvsStage);
  if (s.parse) drawHgvsResult();
  else if (s.text) readHgvs();
  else view.hgvsStage.append(card({ body: empty({ icon: 'variant', title: 'A description of a variant', text: 'Type one, or choose an example: a substitution in the mitochondrion, a deletion in a repeat, a frameshift.' }) }));
}

async function readHgvs() {
  const s = kept.hgvs;
  const t = ++ticket;
  s.busy = true; s.error = null; s.describe = null;
  remember();
  if (view.hgvsLine) { view.hgvsLine.input.value = s.text; view.hgvsLine.input.classList.remove('bad'); }
  clear(view.hgvsStage).append(working('Reading the description by HGVS\'s grammar…'));
  try {
    s.parse = await api.run('nomenclature.parse', { text: s.text });
    if (t !== ticket) return;
    drawHgvsResult();
    const ref = s.genome || HGVS_EXAMPLES.find((e) => e.text === s.text) || genomeFor(s.parse);
    if (s.parse.accepted && ref && /:[gmo]\./.test(s.text)) {
      view.hgvsMore.replaceChildren(working('Placing the change on its reference, in every form, through every transcript…'));
      s.describe = await api.run('nomenclature.describe', { genome: ref.genome, first: ref.first || 1, contig: 0, description: s.text });
      if (t !== ticket) return;
      drawDescribe();
    }
  } catch (e) {
    if (t !== ticket) return;
    s.error = e.message;
    if (view.hgvsMore) view.hgvsMore.replaceChildren(callout(e.message, 'warn'));
    else clear(view.hgvsStage).append(callout(e.message, 'bad'));
  } finally {
    if (t === ticket) s.busy = false;
  }
}

// The reference a description names, when the library holds it.
function genomeFor(parse) {
  const acc = parse?.value?.descriptions?.[0]?.reference?.accession || /^([A-Z]{2}_\d+\.\d+)/.exec(kept.hgvs.text)?.[1];
  if (!acc) return null;
  const known = HGVS_EXAMPLES.find((e) => e.text.startsWith(`${acc}:`));
  if (known) return known;
  if (view.ctx.catalogue?.some((g) => g.id === `ncbi:${acc}` && g.fetched)) return { genome: `ncbi:${acc}`, first: 1 };
  return null;
}

function drawHgvsResult() {
  const s = kept.hgvs;
  const box = clear(view.hgvsStage);
  const r = s.parse;
  if (!r) return;
  if (view.hgvsLine) view.hgvsLine.input.classList.toggle('bad', !r.accepted);
  box.append(derivationFigure(s.text, r));
  if (!r.accepted) return;
  box.append(el('div', { class: 'figure' },
    el('div', { class: 'eyebrow' }, 'What it says'),
    el('p', { class: 'reading' }, r.reading || ''),
    el('div', { class: 'forms' }, el('div', { class: 'form lead' }, el('div', { class: 'k' }, 'written canonically'), el('div', { class: 'v' }, r.canonical || s.text))),
    diagnosticList(r.diagnostics)));
  view.hgvsMore = el('div', { class: 'stack', style: { gap: '16px' } });
  box.append(view.hgvsMore);
  if (s.describe) drawDescribe();
}

// The change on its reference: the forms, the interval of placements on the letters, the transcripts.
function drawDescribe() {
  const d = kept.hgvs.describe;
  const box = clear(view.hgvsMore);
  if (!d) return;
  const forms = el('div', { class: 'forms' },
    form('HGVS', d.hgvs, 'the 3′ rule: as far towards the 3′ end as the change can go', true),
    form('VCF', d.vcf ? `${d.vcf.pos} ${d.vcf.ref} ${d.vcf.alt}` : '–', 'left-aligned and padded (Tan, Abecasis and Kang 2015)'),
    form('SPDI', d.spdi, 'fully justified: the whole interval it could occupy (Holmes et al. 2020)'));
  const cards = [forms];
  if (d.placements && d.context) cards.push(placementsFigure(d));
  box.append(card({ title: 'One change, three normal forms', sub: 'a change in a repeat can be written at every place of an interval; each standard picks a canonical one', pad: true, body: cards }));
  const tx = (d.transcripts || []);
  if (tx.length) {
    box.append(card({ title: 'Through its transcripts to the protein', sub: `${tx.length} ${tx.length === 1 ? 'transcript' : 'transcripts'} near the change`, pad: true,
      body: el('div', { class: 'change-list' }, ...tx.map((t) => el('div', { class: 'change', style: { gridTemplateColumns: '170px 1fr' } },
        el('span', {}, el('b', {}, t.gene || t.id), el('div', { class: 'hint' }, t.id)),
        el('span', {},
          el('div', {}, el('code', {}, t.c || '–'), t.p ? ['  ', el('code', {}, t.p)] : null),
          el('div', { style: { marginTop: '4px', display: 'flex', gap: '6px', flexWrap: 'wrap' } },
            ...(t.consequences || []).map((c) => el('span', { class: `badge ${t.impact === 'HIGH' ? 'bad' : t.impact === 'MODERATE' ? 'warn' : ''}` }, c.replace(/_/g, ' '))),
            t.impact ? el('span', { class: 'badge accent' }, t.impact.toLowerCase()) : null),
          t.note ? el('div', { class: 'hint' }, t.note) : null,
          geneModel(t, d))))) }));
  }
}

// A transcript drawn along the reference: exons as boxes, the coding part full height, introns as a
// line, the strand by chevrons, and the change as a rose mark on it.
function geneModel(t, d) {
  const exons = (t.exons || []).map(([b, e]) => [Math.min(b, e), Math.max(b, e)]).sort((a, b) => a[0] - b[0]);
  if (!exons.length) return null;
  const lo = exons[0][0], hi = exons[exons.length - 1][1], at = d.change?.position ?? d.placements?.first;
  const span = Math.max(1, hi - lo), W = 560, H = 34, pad = 6;
  const X = (p) => pad + ((p - lo) / span) * (W - 2 * pad);
  const cds = t.cds ? [Math.min(...t.cds), Math.max(...t.cds)] : null;
  const out = [`<line x1="${X(lo)}" y1="17" x2="${X(hi)}" y2="17" stroke="var(--hair-3)" />`];
  for (let p = lo + span / 14; p < hi; p += span / 7) {
    const x = X(p), dir = t.strand === '-' || t.strand === -1 ? -1 : 1;
    out.push(`<path d="M${x - 2 * dir},14 l${3 * dir},3 l${-3 * dir},3" fill="none" stroke="var(--ink-faint)" />`);
  }
  for (const [b, e] of exons) {
    out.push(`<rect x="${X(b)}" y="12" width="${Math.max(1.5, X(e) - X(b))}" height="10" fill="var(--accent-2)" fill-opacity="0.55" />`);
    if (cds) {
      const cb = Math.max(b, cds[0]), ce = Math.min(e, cds[1]);
      if (ce > cb) out.push(`<rect x="${X(cb)}" y="8" width="${Math.max(1.5, X(ce) - X(cb))}" height="18" fill="var(--accent)" />`);
    }
  }
  if (Number.isFinite(at) && at >= lo - span * 0.05 && at <= hi + span * 0.05) out.push(`<path d="M${X(at)},3 v28" stroke="var(--T)" stroke-width="2" /><circle cx="${X(at)}" cy="3" r="2.6" fill="var(--T)" />`);
  const fig = el('div', { style: { marginTop: '8px', overflowX: 'auto' } });
  fig.innerHTML = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${out.join('')}</svg>`;
  return fig;
}

function form(k, v, note, lead = false) {
  return el('div', { class: `form${lead ? ' lead' : ''}` }, el('div', { class: 'k' }, k), el('div', { class: 'v' }, v || '–'), el('div', { class: 'note' }, note));
}

// The letters about the change, with the interval of its placements: VCF's place at its left end,
// HGVS's at its right, SPDI's the whole hull.
function placementsFigure(d) {
  const { begin, letters } = d.context;
  const pl = d.placements;
  // Positions are numbered as the reference's bases (letters[i] is base begin + i); the hull is closed.
  const CW = 11, pad = 12, n = letters.length;
  const W = pad * 2 + n * CW, H = 104;
  const x = (pos) => pad + (pos - begin) * CW;
  const out = [];
  const [h0, h1] = pl.hull || [pl.first, pl.last];
  out.push(`<rect x="${x(h0)}" y="24" width="${Math.max(2, (h1 - h0 + 1) * CW)}" height="26" rx="4" fill="var(--accent-softer)" stroke="var(--accent-ring)" />`);
  for (let i = 0; i < n; i++) {
    const b = letters[i];
    out.push(`<text x="${pad + i * CW + CW / 2}" y="42" text-anchor="middle" style="font: 600 14px var(--mono); fill: var(--${'ACGT'.includes(b) ? b : 'ink-faint'})">${esc(b)}</text>`);
  }
  // Each form's placement: the bases it says are changed (an insertion's place between two bases).
  const len = Math.max(0, pl.length || 0);
  const span = (from, label, y, colour, below) => {
    const x0 = x(from), x1 = len ? x(from + len) : x0;
    const tick = below ? 6 : -6;
    return `<path d="M${x0 + 1},${y - tick} v${tick} H${Math.max(x0 + 3, x1 - 1)} v${-tick}" fill="none" stroke="${colour}" stroke-width="1.6" stroke-linejoin="round" />` +
      `<text x="${(x0 + Math.max(x0 + 3, x1)) / 2}" y="${below ? y + 15 : y - 9}" text-anchor="middle" style="font: italic 11.5px var(--serif); fill: ${colour}">${esc(label)}</text>`;
  };
  if (Number.isFinite(pl.first)) out.push(span(pl.first, 'VCF, left', 60, 'var(--C)', true));
  if (Number.isFinite(pl.last)) out.push(span(pl.last, 'HGVS, 3′', 86, 'var(--accent)', true));
  out.push(`<text x="${x(h0)}" y="14" style="font: italic 11.5px var(--serif); fill: var(--ink-faint)">SPDI, the hull of every placement</text>`);
  const fig = el('div', { style: { overflowX: 'auto', marginTop: '12px' } });
  fig.innerHTML = `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${out.join('')}</svg>`;
  return fig;
}
