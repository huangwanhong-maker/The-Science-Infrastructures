// The Index view: how a text becomes searchable. For a short word: its suffixes sorted (the
// suffix array) with their longest common prefixes, its rotations sorted, whose last column is
// the Burrows–Wheeler transform; the LF map that walks the text backwards through the BWT, and
// backward search, which narrows an interval of rows one letter of the pattern at a time; and
// the word's Lyndon factorisation, through which the bijective transform needs no sentinel. For
// a genome: its FMD-index (the genome and its reverse complement, so that one search finds both
// strands), patterns counted and located in it, and a read's super-maximal exact matches, the
// seeds from which a mapper starts.

import { el, clear, esc, notice, toggles, field } from '../dom.js';
import * as api from '../api.js';
import * as fmt from '../format.js';
import { BASE_COLOURS, token, CATEGORIES } from '../palette.js';

const WORDS = ['mississippi', 'banana', 'abracadabra', 'GATTACAGATTACA', 'ACGTACGTTTGCA'];
let view = null;

export function mount(section, ctx) {
  view = { ctx, section, word: el('div'), genome: el('div'), text: ctx.params.text || 'mississippi', pattern: ctx.params.pattern || 'ssi', step: null, data: null };
  section.append(
    el('h2', { class: 'view-title' }, 'Index'),
    el('p', { class: 'lede' }, 'A genome is searched through an index of its suffixes. Sort every suffix of a text and you have its suffix array; sort its rotations instead and their last column is the Burrows–Wheeler transform, a permutation of the text that groups equal letters, which compresses, and from which the text and the array can be recovered. Searching it backwards, one letter of the pattern at a time, narrows a range of rows whose size is the number of occurrences: the FM-index. For DNA the index holds the genome and its reverse complement together, so a search finds both strands, and extends a match in either direction.'),
    el('h3', { class: 'section-title' }, 'A word'), view.word,
    el('h3', { class: 'section-title' }, 'A genome'), view.genome);
  buildWord();
  buildGenome();
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.data) computeWord();
  drawGenomeHeader();
}

// ---- a word ----------------------------------------------------------------------------------------
function buildWord() {
  const box = clear(view.word);
  const text = el('input', { value: view.text, style: { width: '260px', font: '14px var(--mono)' }, onkeydown: (e) => { if (e.key === 'Enter') { view.text = text.value.trim(); computeWord(); } } });
  const pattern = el('input', { value: view.pattern, style: { width: '120px', font: '14px var(--mono)' }, onkeydown: (e) => { if (e.key === 'Enter') { view.pattern = pattern.value.trim(); computeWord(); } } });
  const examples = el('span', { class: 'row' }, ...WORDS.map((w) => el('button', { class: 'link small', onclick: () => { view.text = w; text.value = w; computeWord(); } }, w)));
  box.append(el('div', { class: 'row' }, field('text', text), field('pattern', pattern), el('button', { class: 'button primary', onclick: () => { view.text = text.value.trim(); view.pattern = pattern.value.trim(); computeWord(); } }, 'Index it'), examples));
  view.wordBody = el('div');
  box.append(view.wordBody);
}

async function computeWord() {
  if (!view.text) return;
  view.ctx.setParams({ text: view.text, pattern: view.pattern });
  try {
    view.data = await api.run('index.word', { text: view.text, pattern: view.pattern });
    view.step = view.data.search ? view.data.search.length : null;
    drawWord();
  } catch (e) { clear(view.wordBody).append(el('p', { class: 'error' }, e.message)); }
}

function drawWord() {
  const d = view.data, box = clear(view.wordBody);
  const n = d.sa.length;
  // The interval of the backward search after `step` letters (all rows before the first).
  const interval = view.step === null || !d.search ? null : view.step === 0 ? [0, n] : [d.search[view.step - 1].lo, d.search[view.step - 1].hi];
  const table = el('table', { class: 'data bwt' });
  table.append(el('tr', {}, el('th', { class: 'num' }, 'row'), el('th', { class: 'num' }, 'SA'), el('th', { class: 'num' }, 'LCP'), el('th', {}, 'sorted rotations: F … L'), el('th', {}, 'the suffix')));
  for (let r = 0; r < n; r++) {
    const rot = d.rotations ? d.rotations[r] : null;
    const inside = interval && r >= interval[0] && r < interval[1];
    const rotCell = el('td', { class: 'mono rot' });
    if (rot) {
      [...rot].forEach((c, k) => rotCell.append(el('span', { class: k === 0 ? 'first' : k === rot.length - 1 ? 'last' : '' }, c)));
    }
    table.append(el('tr', { class: inside ? 'inside' : '' },
      el('td', { class: 'num faint' }, String(r)), el('td', { class: 'num' }, String(d.sa[r])), el('td', { class: 'num faint' }, String(d.lcp[r] ?? '')),
      rotCell, el('td', { class: 'mono faint' }, d.text.slice(d.sa[r]))));
  }
  const facts = el('dl', { class: 'kv' });
  const row = (k, v) => facts.append(el('dt', {}, k), el('dd', {}, v));
  row('BWT (the last column)', el('span', { class: 'mono' }, d.bwt));
  row('runs of the BWT', `${d.runs}, against ${n} letters`);
  row('distinct substrings', `${fmt.count(d.distinct)} = n(n + 1)/2 − Σ LCP`);
  if (d.lyndon) {
    const lw = el('span', { class: 'mono' });
    d.lyndon.forEach(([b, e], k) => lw.append(el('span', { style: { background: `${CATEGORIES[k % CATEGORIES.length]}33`, borderBottom: `2px solid ${CATEGORIES[k % CATEGORIES.length]}`, padding: '0 1px', marginRight: '2px' } }, d.word.slice(b, e))));
    row('Lyndon factorisation', lw);
  }
  if (d.bijective) row('bijective BWT (no sentinel)', el('span', { class: 'mono' }, d.bijective));
  if (d.inverse !== undefined) row('the text from the BWT by LF', el('span', { class: 'mono' }, d.inverse));
  // The search, step by step.
  const steps = el('div');
  if (d.search) {
    const ctrl = el('div', { class: 'row' },
      el('span', { class: 'faint small' }, `backward search for “${d.pattern}”:`),
      el('button', { class: 'button small', disabled: view.step <= 0, onclick: () => { view.step = Math.max(0, view.step - 1); drawWord(); } }, '‹'),
      el('button', { class: 'button small', disabled: view.step >= d.search.length, onclick: () => { view.step = Math.min(d.search.length, view.step + 1); drawWord(); } }, '›'),
      el('span', { class: 'small' }, view.step === 0 ? 'every row: the empty pattern' : `after “${d.pattern.slice(d.pattern.length - view.step)}”: rows [${interval[0]}, ${interval[1]}), ${interval[1] - interval[0]} occurrence${interval[1] - interval[0] === 1 ? '' : 's'}`));
    steps.append(ctrl, el('p', { class: 'note' }, 'Each step prepends a letter c to the matched suffix of the pattern: the new interval is C(c) + occ(c, lo) to C(c) + occ(c, hi), where C(c) counts the letters smaller than c and occ(c, i) the c’s among the first i letters of the BWT (Ferragina and Manzini 2000). The rows inside are the rotations that begin with what has been matched.'));
  }
  box.append(el('div', { class: 'row', style: { alignItems: 'flex-start', gap: '30px' } }, el('div', { class: 'bwt-wrap' }, table), el('div', { style: { maxWidth: '520px' } }, facts, steps)));
}

// ---- a genome ----------------------------------------------------------------------------------------
function buildGenome() {
  const box = clear(view.genome);
  view.genomeHead = el('div');
  const pattern = el('input', { placeholder: 'a pattern, ACGT…', style: { width: '240px', font: '13px var(--mono)' } });
  const searchOut = el('div');
  const read = el('textarea', { rows: 3, placeholder: 'a read; or simulate one', style: { width: '100%', maxWidth: '760px' } });
  const smemOut = el('div');
  box.append(view.genomeHead,
    el('div', { class: 'row' }, field('count and locate', pattern), el('button', { class: 'button', onclick: async () => {
      try {
        const r = await api.run('index.count', { genome: view.ctx.genome, pattern: pattern.value.trim() });
        clear(searchOut).append(el('p', { class: 'note' }, `${fmt.count(r.count)} occurrence${r.count === 1 ? '' : 's'} on both strands${r.palindrome ? ` (the pattern is its own reverse complement, so each of its ${fmt.count(r.count / 2)} sites is found once on each strand)` : ''}, found by ${r.steps} backward steps in ${fmt.seconds(r.seconds)}${r.located?.length ? `; the first ${r.located.length} located:` : ''}`),
          r.located ? el('p', { class: 'mono small' }, r.located.map(([c, p, rev]) => `${fmt.position(p)}${rev ? '−' : '+'}`).join('  ')) : null);
        drawGenomeHeader(r.index);
      } catch (e) { notice(e.message, true); }
    } }, 'Search')), searchOut,
    el('p', { class: 'eyebrow', style: { marginTop: '18px' } }, 'the super-maximal exact matches of a read'),
    read,
    el('div', { class: 'row' },
      el('button', { class: 'button small', onclick: async () => {
        try { const r = await api.run('index.sample_read', { genome: view.ctx.genome, length: 250, error: 0.04, seed: Math.floor(Math.random() * 1e9) }); read.value = r.read; } catch (e) { notice(e.message, true); }
      } }, 'simulate a read with 4% errors'),
      el('button', { class: 'button', onclick: async () => {
        try {
          const r = await api.run('index.smems', { genome: view.ctx.genome, read: read.value.replace(/\s+/g, ''), min_length: 17 });
          drawSmems(smemOut, r, read.value.replace(/\s+/g, ''));
          drawGenomeHeader(r.index);
        } catch (e) { notice(e.message, true); }
      } }, 'Find its SMEMs')),
    smemOut);
}

function drawGenomeHeader(index) {
  if (!view.genomeHead) return;
  const g = view.ctx.genomeInfo();
  clear(view.genomeHead).append(el('p', { class: 'note' }, index
    ? `The FMD-index of ${g?.title || view.ctx.genome}: ${fmt.length(index.length)} and its reverse complement, ${fmt.count(index.bwt_length)} letters of BWT in ${fmt.length(index.bytes).replace('bp', 'bytes')} (sampled suffix array every ${index.sa_sample}), built in ${fmt.seconds(index.build_seconds)}.`
    : `The genome chosen above (${g?.title || view.ctx.genome}) is indexed when first searched.`));
}

function drawSmems(box, r, readText) {
  clear(box);
  const n = readText.length;
  const W = Math.min(900, view.section.clientWidth - 80), rowH = 15;
  const x = (i) => 10 + (i / n) * (W - 20);
  const colours = BASE_COLOURS();
  let s = `<svg width="${W}" height="${40 + r.smems.length * rowH}" viewBox="0 0 ${W} ${40 + r.smems.length * rowH}">`;
  for (let i = 0; i < n; i++) s += `<rect x="${x(i).toFixed(2)}" y="6" width="${Math.max(0.6, (W - 20) / n).toFixed(2)}" height="12" fill="${colours[readText[i]] || colours.N}" />`;
  r.smems.forEach((m, k) => {
    const y = 28 + k * rowH;
    s += `<rect x="${x(m.begin).toFixed(1)}" y="${y}" width="${(x(m.end) - x(m.begin)).toFixed(1)}" height="9" rx="2" fill="${m.hits === 1 ? token('--accent') : token('--C')}" fill-opacity="${m.hits === 1 ? 0.85 : 0.5}"><title>${m.end - m.begin} bases at ${m.begin + 1}, ${m.hits} hit${m.hits === 1 ? '' : 's'}</title></rect>`;
    s += `<text class="tick-label" x="${(x(m.end) + 4).toFixed(1)}" y="${y + 8}">${m.end - m.begin} bp × ${m.hits}</text>`;
  });
  s += '</svg>';
  const fig = el('div');
  fig.innerHTML = s;
  box.append(el('p', { class: 'note' }, `${r.smems.length} SMEMs of at least ${r.min_length} bases, found in ${fmt.seconds(r.seconds)}: each a match that cannot be extended in either direction and is contained in no longer one (Li 2012). Mulberry: unique in the genome; blue: repeated. A mapper chains them and aligns the gaps between.`), fig);
}
