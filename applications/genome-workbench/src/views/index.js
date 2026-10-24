// The Index view: how a text becomes searchable. For a short word: its suffixes sorted (the
// suffix array) with their longest common prefixes, its rotations sorted, whose last column is
// the Burrows–Wheeler transform; the LF map that walks the text backwards through the BWT, and
// backward search, which narrows an interval of rows one letter of the pattern at a time; and
// the word's Lyndon factorisation, through which the bijective transform needs no sentinel. For
// a genome: its FMD-index (the genome and its reverse complement, so that one search finds both
// strands), patterns counted and located in it, and a read's super-maximal exact matches, the
// seeds from which a mapper starts.
//
// Two benches, one under the other: the word's controls beside its table, its search and its
// transform; the genome's beside its index, its counts and the seeds of a read. Enter runs either;
// ← and → step the backward search. The word and its pattern live in the address (text, pattern).

import { el, clear } from '../dom.js';
import * as api from '../api.js';
import * as fmt from '../format.js';
import { BASE_COLOURS, token, CATEGORIES } from '../palette.js';
import { icon } from '../ui/icons.js';
import { btn, iconBtn, fld, card, metrics, kv, callout, empty, seqEditor, numberField, toast } from '../ui/kit.js';

const WORDS = ['mississippi', 'banana', 'abracadabra', 'GATTACAGATTACA', 'ACGTACGTTTGCA'];
let view = null;
let keysBound = false;

export function mount(section, ctx) {
  view = { ctx, section, text: ctx.params.text || 'mississippi', pattern: ctx.params.pattern ?? 'ssi', step: null, data: null,
           genomeOf: null, index: null, minLength: 17 };
  section.append(
    el('h2', { class: 'view-title' }, 'Index'),
    el('p', { class: 'lede' }, 'A genome is searched through an index of its suffixes. Sort every suffix of a text and you have its suffix array; sort its rotations instead and their last column is the Burrows–Wheeler transform, a permutation of the text that groups equal letters, which compresses, and from which the text and the array can be recovered. Searching it backwards, one letter of the pattern at a time, narrows a range of rows whose size is the number of occurrences: the FM-index. For DNA the index holds the genome and its reverse complement together, so a search finds both strands, and extends a match in either direction.'),
    el('div', { class: 'section-head ix-first' }, el('h3', {}, 'A word'), el('span', { class: 'sub' }, 'its suffixes and rotations sorted, and a pattern searched backwards')),
    wordBench(),
    el('div', { class: 'section-head' }, el('h3', {}, 'A genome'), el('span', { class: 'sub' }, 'its FMD-index: patterns counted and located, and the seeds of a read')),
    genomeBench());
  if (!keysBound) { window.addEventListener('keydown', onKey); keysBound = true; }
}

export function show(ctx) {
  view.ctx = ctx;
  if (!view.data) computeWord();
  // Results belong to the genome they were found in: another genome starts afresh.
  if (view.genomeOf !== ctx.genome) resetGenome();
  drawGenomeHeader();
}

// A button that shows it is working while it runs, and stays where it is.
function busy(b, on, label) {
  if (on) {
    b.keep = [...b.childNodes];
    b.disabled = true;
    b.replaceChildren(el('span', { class: 'spin' }), el('span', {}, label));
  } else {
    b.disabled = false;
    if (b.keep) b.replaceChildren(...b.keep);
  }
}

const mono = (attrs) => el('input', { type: 'text', class: 'ix-mono', spellcheck: false, autocomplete: 'off', ...attrs });
const toolboxLink = (text, label, onClick) => el('p', { class: 'toolbox-hint' }, icon('toolbox', 'sm'),
  el('span', {}, text, el('button', { class: 'link', onclick: onClick }, label), '.'));

// ---- a word ----------------------------------------------------------------------------------------
function wordBench() {
  const text = mono({ value: view.text, placeholder: 'any word', 'aria-label': 'the text' });
  const pattern = mono({ value: view.pattern, placeholder: 'a pattern', 'aria-label': 'the pattern' });
  const go = () => { view.text = text.value.trim(); view.pattern = pattern.value.trim(); computeWord(); };
  for (const i of [text, pattern]) i.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); go(); } });
  view.textInput = text;
  view.patternInput = pattern;
  view.examples = el('div', { class: 'chips' }, ...WORDS.map((w) => el('button', { type: 'button', class: `chip ix-word-chip${w === view.text ? ' on' : ''}`, 'data-w': w, onclick: () => { text.value = w; go(); } }, w)));
  view.wordRun = btn('Index it', { icon: 'play', kind: 'primary', onClick: go });
  const controls = el('div', { class: 'controls' },
    el('div', { class: 'card' },
      el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'the word'),
        fld('text', text, { hint: 'any letters; the sentinel $ is added at its end' }),
        fld('pattern', pattern, { help: 'searched backwards' })),
      el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'or one of these'), view.examples),
      el('div', { class: 'run-bar' }, view.wordRun)),
    el('p', { class: 'small faint ix-keys' }, el('kbd', {}, 'Enter'), ' indexes it; ', el('kbd', {}, '←'), ' ', el('kbd', {}, '→'), ' step the search.'),
    toolboxLink('A longer text, or a file of yours: ', 'the suffix array and BWT in the Toolbox', () => view.ctx.openTool('index-sa', { text: view.textInput.value.trim() || view.text })));
  view.wordStage = el('div', { class: 'stage' }, card({ body: el('div', { class: 'wait ix-wait' }, el('span', { class: 'spin' }), 'Sorting the suffixes…') }));
  return el('div', { class: 'bench' }, controls, view.wordStage);
}

async function computeWord() {
  if (!view.text) { toast('A text is needed', 'error'); return; }
  view.ctx.setParams({ text: view.text, pattern: view.pattern });
  for (const c of view.examples.children) c.classList.toggle('on', c.dataset.w === view.text);
  busy(view.wordRun, true, 'Indexing…');
  view.wordStage.classList.add('ix-stale');
  try {
    view.data = await api.run('index.word', { text: view.text, pattern: view.pattern });
    view.step = view.data.search ? view.data.search.length : null;
    drawWord();
  } catch (e) {
    view.wordStage.querySelector(':scope > .callout')?.remove();
    view.wordStage.prepend(callout(e.message, 'bad'));
  } finally {
    busy(view.wordRun, false);
    view.wordStage.classList.remove('ix-stale');
  }
}

function drawWord() {
  const d = view.data, box = clear(view.wordStage);
  const n = d.sa.length;
  const steps = d.search?.length || 0;
  // The interval of the backward search after `step` letters (all rows before the first).
  const interval = view.step === null || !d.search ? null : view.step === 0 ? [0, n] : [d.search[view.step - 1].lo, d.search[view.step - 1].hi];
  const matched = interval ? d.pattern.slice(d.pattern.length - view.step) : '';
  const found = steps ? d.search[steps - 1].hi - d.search[steps - 1].lo : null;
  box.append(metrics([[fmt.count(n), 'rows: its letters and $'], [fmt.count(d.runs), 'runs in the BWT'], [fmt.count(d.distinct), 'distinct substrings'],
    found !== null ? [fmt.count(found), `occurrence${found === 1 ? '' : 's'} of “${d.pattern}”`, found ? 'good' : ''] : null]));

  // The table: rows of the suffix array, their LCPs, the sorted rotations from F to L, the suffixes.
  const table = el('table', { class: 'data bwt' });
  table.append(el('tr', {}, el('th', { class: 'num' }, 'row'), el('th', { class: 'num' }, 'SA'), el('th', { class: 'num' }, 'LCP'), el('th', {}, 'sorted rotations: F … L'), el('th', {}, 'the suffix')));
  for (let r = 0; r < n; r++) {
    const rot = d.rotations ? d.rotations[r] : null;
    const inside = interval && r >= interval[0] && r < interval[1];
    const rotCell = el('td', { class: 'mono rot' });
    if (rot) {
      // Inside the interval, the rotation begins with what has been matched: underlined.
      [...rot].forEach((c, k) => rotCell.append(el('span', { class: [k === 0 ? 'first' : k === rot.length - 1 ? 'last' : '', inside && k < matched.length ? 'match' : ''].join(' ').trim() || null }, c)));
    }
    table.append(el('tr', { class: inside ? 'inside' : '' },
      el('td', { class: 'num faint' }, String(r)), el('td', { class: 'num' }, String(d.sa[r])), el('td', { class: 'num faint' }, String(d.lcp[r] ?? '')),
      rotCell, el('td', { class: 'mono faint' }, d.text.slice(d.sa[r]))));
  }
  const tableCard = card({ title: 'The sorted rotations', sub: interval ? `the rows of the interval shaded; F in blue, L (the BWT) in mulberry` : 'F in blue, L (the BWT) in mulberry', body: el('div', { class: 'bwt-wrap' }, table) });

  // The search, step by step: the pattern's letters, those matched so far lit, and the interval.
  const side = el('div', { class: 'ix-side' });
  if (d.search && d.pattern) {
    const m = d.pattern.length;
    const setStep = (k) => { view.step = Math.max(0, Math.min(steps, k)); drawWord(); };
    const letters = el('div', { class: 'ix-letters', role: 'group', 'aria-label': 'the pattern, matched from its end' },
      ...[...d.pattern].map((c, j) => {
        const k = m - j;   // the step at which this letter has been matched
        return el('button', { type: 'button', class: `ix-letter${k <= view.step ? ' on' : ''}${k === view.step + 1 ? ' next' : ''}`, disabled: k > steps || null, title: k > steps ? 'the search stopped before this letter' : `after ${k} letter${k === 1 ? '' : 's'}`, onclick: () => setStep(k === view.step ? k - 1 : k) }, c);
      }));
    const prev = iconBtn('chevron', { title: 'a letter fewer (←)', size: 'sm', onClick: () => setStep(view.step - 1) });
    prev.classList.add('flip');
    prev.disabled = view.step <= 0;
    const next = iconBtn('chevron', { title: 'a letter more (→)', size: 'sm', onClick: () => setStep(view.step + 1) });
    next.disabled = view.step >= steps;
    const size = interval[1] - interval[0];
    const now = view.step === 0 ? 'every row: the empty pattern' : `after “${matched}”: rows [${interval[0]}, ${interval[1]}), ${size} occurrence${size === 1 ? '' : 's'}`;
    const list = el('table', { class: 'data ix-steps' }, el('tr', {}, el('th', { class: 'num' }, 'step'), el('th', {}, 'letter'), el('th', {}, 'interval'), el('th', { class: 'num' }, 'rows')),
      ...[{ letter: '', lo: 0, hi: n }, ...d.search].map((s, k) => el('tr', { class: k === view.step ? 'on' : '', onclick: () => setStep(k) },
        el('td', { class: 'num faint' }, String(k)), el('td', { class: 'mono' }, s.letter || '—'), el('td', { class: 'mono' }, `[${s.lo}, ${s.hi})`), el('td', { class: 'num' }, String(s.hi - s.lo)))));
    side.append(card({ title: `Backward search for “${d.pattern}”`, sub: `step ${view.step} of ${steps}`, actions: [prev, next], body: el('div', { class: 'card-body stack ix-search' },
      letters, el('div', { class: 'ix-now' }, now), list,
      el('p', { class: 'note' }, 'Each step prepends a letter c to the matched suffix of the pattern: the new interval is C(c) + occ(c, lo) to C(c) + occ(c, hi), where C(c) counts the letters smaller than c and occ(c, i) the c’s among the first i letters of the BWT (Ferragina and Manzini 2000). The rows inside are the rotations that begin with what has been matched.')) }));
  }

  // The transform and what follows from it.
  const pairs = [['BWT (the last column)', el('span', { class: 'mono' }, d.bwt)], ['runs of the BWT', `${d.runs}, against ${n} letters`], ['distinct substrings', `${fmt.count(d.distinct)} = n(n + 1)/2 − Σ LCP`]];
  if (d.lyndon) {
    const lw = el('span', { class: 'mono' });
    d.lyndon.forEach(([b, e], k) => lw.append(el('span', { style: { background: `${CATEGORIES[k % CATEGORIES.length]}33`, borderBottom: `2px solid ${CATEGORIES[k % CATEGORIES.length]}`, padding: '0 1px', marginRight: '2px' } }, d.word.slice(b, e))));
    pairs.push(['Lyndon factorisation', lw]);
  }
  if (d.bijective) pairs.push(['bijective BWT (no sentinel)', el('span', { class: 'mono' }, d.bijective)]);
  if (d.inverse !== undefined) pairs.push(['the text from the BWT by LF', el('span', { class: 'mono' }, d.inverse)]);
  side.append(card({ title: 'Its transform', sub: 'and what can be read back from it', pad: true, body: kv(pairs) }));
  box.append(el('div', { class: 'ix-word' }, tableCard, side));
}

// ---- a genome ----------------------------------------------------------------------------------------
function genomeBench() {
  view.genomeHead = el('div', { class: 'ix-genome' });
  const pattern = mono({ placeholder: 'ACGT…', 'aria-label': 'a pattern to count and locate' });
  view.genomePattern = pattern;
  view.countBtn = btn('Count and locate', { icon: 'search', kind: 'primary', onClick: () => count() });
  pattern.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); count(); } });
  view.read = seqEditor({ rows: 4, alphabet: 'dna', placeholder: 'a read, pasted; or simulate one from the genome' });
  view.read.textarea.addEventListener('keydown', (e) => { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); smems(); } });
  const minLength = numberField({ value: view.minLength, min: 8, max: 100, step: 1, unit: 'bases', onChange: (v) => { view.minLength = v; } });
  view.smemBtn = btn('Find its SMEMs', { icon: 'play', kind: 'primary', onClick: () => smems() });
  view.simBtn = btn('Simulate', { icon: 'spark', title: '250 bases drawn from the genome, with 4% errors', onClick: () => simulate() });
  const controls = el('div', { class: 'controls' },
    el('div', { class: 'card' }, el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'the genome in hand'), view.genomeHead)),
    el('div', { class: 'card' },
      el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'count and locate'),
        fld('pattern', pattern, { hint: 'found on both strands by one backward search' })),
      el('div', { class: 'run-bar' }, view.countBtn)),
    el('div', { class: 'card' },
      el('div', { class: 'ctl-group' }, el('div', { class: 'eyebrow' }, 'the seeds of a read'),
        fld('read', view.read), fld('matches of at least', minLength)),
      el('div', { class: 'run-bar' }, view.smemBtn, view.simBtn)),
    toolboxLink('Your own sequences: ', 'find motifs in them', () => view.ctx.openTool('seq-find', view.genomePattern.value.trim() ? { patterns: view.genomePattern.value.trim() } : null)),
    toolboxLink('Your own reads, mapped from their seeds: ', 'map reads', () => view.ctx.openTool('reads-map')));
  view.indexCard = el('div');
  view.countCard = el('div');
  view.smemCard = el('div');
  view.genomeStage = el('div', { class: 'stage' }, view.indexCard, view.countCard, view.smemCard);
  return el('div', { class: 'bench' }, controls, view.genomeStage);
}

function resetGenome() {
  view.genomeOf = view.ctx.genome;
  view.index = null;
  drawIndexCard();
  clear(view.countCard).append(card({ body: empty({ icon: 'search', title: 'Nothing counted yet',
    text: 'A pattern is counted in as many backward steps as it has letters, whatever the length of the genome; its places are then read from the sampled suffix array.',
    action: btn('Count GAATTC', { icon: 'play', onClick: () => { view.genomePattern.value = 'GAATTC'; count(); } }) }) }));
  clear(view.smemCard).append(card({ body: empty({ icon: 'index', title: 'No read yet',
    text: 'Paste a read, or simulate one from the genome, and find its super-maximal exact matches: the seeds from which a mapper starts.',
    action: btn('Simulate one and find them', { icon: 'spark', onClick: async () => { if (await simulate()) smems(); } }) }) }));
}

function drawGenomeHeader() {
  if (!view.genomeHead) return;
  const g = view.ctx.genomeInfo();
  clear(view.genomeHead).append(
    el('div', { class: 'ix-genome-name' }, icon(g?.group === 'files' ? 'file' : 'genome', 'sm'), el('span', {}, g ? g.title.replace(/, complete genome.*$/, '') : view.ctx.genome || 'no genome')),
    el('div', { class: 'small faint' }, view.index ? `indexed: ${fmt.length(view.index.length)} and its reverse complement` : 'indexed when first searched; another genome is chosen at the top right'));
}

function drawIndexCard() {
  const g = view.ctx.genomeInfo();
  const name = g?.title?.replace(/, complete genome.*$/, '') || view.ctx.genome;
  const i = view.index;
  clear(view.indexCard).append(card({ title: 'The FMD-index', sub: `of ${name}, with its reverse complement`, pad: true,
    body: i ? metrics([[fmt.length(i.length), 'the genome'], [fmt.count(i.bwt_length), 'letters of BWT, both strands'], [fmt.bytes(i.bytes), 'in memory'], [`every ${i.sa_sample}`, 'suffix array, sampled'], [fmt.seconds(i.build_seconds), 'to build']])
      : el('p', { class: 'note ix-flush' }, 'It is built when the genome is first searched, and kept by the engine: a count or a read below builds it.') }));
  drawGenomeHeader();
}

function noteIndex(index) {
  if (!index) return;
  view.index = index;
  drawIndexCard();
}

async function count() {
  const q = view.genomePattern.value.trim();
  if (!q) { toast('A pattern is needed', 'error'); view.genomePattern.focus(); return; }
  if (!view.ctx.genome) { toast('No genome is in hand', 'error'); return; }
  const genome = view.ctx.genome;
  busy(view.countBtn, true, 'Counting…');
  view.countCard.classList.add('ix-stale');
  try {
    const r = await api.run('index.count', { genome, pattern: q });
    if (view.ctx.genome !== genome) return;
    const head = el('p', { class: 'note ix-flush' }, `${fmt.count(r.count)} occurrence${r.count === 1 ? '' : 's'} on both strands${r.palindrome ? ` (the pattern is its own reverse complement, so each of its ${fmt.count(r.count / 2)} sites is found once on each strand)` : ''}, found by ${r.steps} backward steps in ${fmt.seconds(r.seconds)}${r.located?.length ? `; the first ${r.located.length} located, each a way into the Genome view:` : ''}`);
    // Each place opens the Genome view there.
    const places = r.located?.length ? el('div', { class: 'chips ix-places' }, ...r.located.map(([c, p, rev]) => el('button', { type: 'button', class: 'chip', title: 'show it in the Genome view',
      onclick: () => view.ctx.go('genome', { ...(c ? { contig: c } : {}), begin: Math.max(0, p - 120), end: p + 120 + q.length }) }, el('span', { class: 'mono' }, `${fmt.position(p)}${rev ? '−' : '+'}`)))) : null;
    clear(view.countCard).append(card({ title: `“${q}”`, sub: 'counted and located', pad: true, body: el('div', { class: 'stack' },
      metrics([[fmt.count(r.count), 'occurrences', r.count ? 'good' : ''], [fmt.count(r.steps), 'backward steps'], [fmt.seconds(r.seconds), 'to count']]), head, places) }));
    noteIndex(r.index);
  } catch (e) {
    toast(e.message, 'error');
  } finally {
    busy(view.countBtn, false);
    view.countCard.classList.remove('ix-stale');
  }
}

async function simulate() {
  if (!view.ctx.genome) { toast('No genome is in hand', 'error'); return false; }
  busy(view.simBtn, true, 'Drawing…');
  try {
    const r = await api.run('index.sample_read', { genome: view.ctx.genome, length: 250, error: 0.04, seed: Math.floor(Math.random() * 1e9) });
    view.read.set(r.read);
    return true;
  } catch (e) { toast(e.message, 'error'); return false; } finally { busy(view.simBtn, false); }
}

async function smems() {
  const readText = view.read.value().replace(/\s+/g, '');
  if (!readText) { toast('A read is needed: paste one, or simulate one', 'error'); return; }
  const genome = view.ctx.genome;
  busy(view.smemBtn, true, 'Finding…');
  view.smemCard.classList.add('ix-stale');
  try {
    const r = await api.run('index.smems', { genome, read: readText, min_length: view.minLength });
    if (view.ctx.genome !== genome) return;
    drawSmems(view.smemCard, r, readText);
    noteIndex(r.index);
  } catch (e) {
    toast(e.message, 'error');
  } finally {
    busy(view.smemBtn, false);
    view.smemCard.classList.remove('ix-stale');
  }
}

function drawSmems(box, r, readText) {
  clear(box);
  const n = readText.length;
  // The read spans the figure less a margin at the right for the last match's label.
  const W = Math.max(320, Math.min(1100, view.genomeStage.clientWidth - 36)), rowH = 15, span = W - 96;
  const x = (i) => 10 + (i / n) * span;
  const colours = BASE_COLOURS();
  let s = `<svg width="${W}" height="${40 + r.smems.length * rowH}" viewBox="0 0 ${W} ${40 + r.smems.length * rowH}">`;
  for (let i = 0; i < n; i++) s += `<rect x="${x(i).toFixed(2)}" y="6" width="${Math.max(0.6, span / n).toFixed(2)}" height="12" fill="${colours[readText[i]] || colours.N}" />`;
  r.smems.forEach((m, k) => {
    const y = 28 + k * rowH;
    s += `<rect x="${x(m.begin).toFixed(1)}" y="${y}" width="${(x(m.end) - x(m.begin)).toFixed(1)}" height="9" rx="2" fill="${m.hits === 1 ? token('--accent') : token('--C')}" fill-opacity="${m.hits === 1 ? 0.85 : 0.5}"><title>${m.end - m.begin} bases at ${m.begin + 1}, ${m.hits} hit${m.hits === 1 ? '' : 's'}</title></rect>`;
    s += `<text class="tick-label" x="${(x(m.end) + 4).toFixed(1)}" y="${y + 8}">${m.end - m.begin} bp × ${m.hits}</text>`;
  });
  s += '</svg>';
  const fig = el('div', { class: 'ix-fig' });
  fig.innerHTML = s;
  const unique = r.smems.filter((m) => m.hits === 1).length;
  const key = (colour, text) => el('span', { class: 'k' }, el('span', { class: 'swatch', style: { background: colour } }), text);
  box.append(card({ title: 'The super-maximal exact matches', sub: `of a read of ${fmt.count(n)} bases`, body: el('div', {},
    el('div', { class: 'card-body stack' },
      metrics([[fmt.count(r.smems.length), `SMEMs of at least ${r.min_length} bases`], [fmt.count(unique), 'unique in the genome'], [fmt.seconds(r.seconds), 'to find them']]),
      el('p', { class: 'note ix-flush' }, `Each a match that cannot be extended in either direction and is contained in no longer one (Li 2012). A mapper chains them and aligns the gaps between.`),
      fig),
    el('div', { class: 'key-row' }, key(token('--accent'), 'unique in the genome'), key(token('--C'), 'repeated'), el('span', { class: 'k faint' }, 'the read above them, base by base'))) }));
}

// ← and → step the backward search while the view is shown and nothing is being typed.
function onKey(e) {
  if (!view || view.section.hidden || !view.data?.search || e.ctrlKey || e.metaKey || e.altKey) return;
  const a = document.activeElement;
  if (/^(INPUT|TEXTAREA|SELECT)$/.test(a?.tagName || '') || a?.isContentEditable) return;
  if (document.querySelector('.dialog-veil, .popover')) return;
  const steps = view.data.search.length;
  if (e.key === 'ArrowLeft' && view.step > 0) view.step--;
  else if (e.key === 'ArrowRight' && view.step < steps) view.step++;
  else return;
  e.preventDefault();
  drawWord();
}
