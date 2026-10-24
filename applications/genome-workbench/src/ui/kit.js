// The workbench's kit of components: each a function returning a DOM element, styled by
// styles.css §4. Views and tools are built from these, so that a control looks and behaves the
// same everywhere and a new tool needs no styling of its own.
//
//   btn, iconBtn, seg, select, numberField, slider, toggle, fld      controls
//   seqEditor, dropzone                                                inputs of sequences and files
//   card, tabs, metrics, callout, empty, table, kv                     surfaces and results
//   toast, popover, menu, dialog, vhead                                overlays and heads
//   working, busy, onEnter                                             work in progress
//   download, copy, readFile, fmtBytes                                 helpers

import { el, clear, esc, toast } from '../dom.js';
import { icon } from './icons.js';

// ---- controls ----------------------------------------------------------------------------------------
export function btn(label, { icon: ic = null, kind = '', size = '', onClick = null, title = null, disabled = false, type = 'button' } = {}) {
  const b = el('button', { class: `btn ${kind} ${size}`.trim(), type, title, disabled, onclick: onClick });
  if (ic) b.append(icon(ic, 'sm'));
  if (label) b.append(el('span', {}, label));
  return b;
}

// A button that shows it is working while its action's promise is pending.
export function actionBtn(label, action, opts = {}) {
  const b = btn(label, { ...opts, onClick: async () => {
    if (b.disabled) return;
    const keep = [...b.childNodes];
    b.disabled = true;
    clear(b).append(el('span', { class: 'spin' }), el('span', {}, opts.busy || label));
    try { await action(); } finally { clear(b).append(...keep); b.disabled = false; }
  } });
  return b;
}

export function iconBtn(name, { title = '', onClick = null, pressed = null, size = '' } = {}) {
  const b = el('button', { class: `icon-btn ${size}`.trim(), title, 'aria-label': title, onclick: onClick });
  if (pressed !== null) b.setAttribute('aria-pressed', String(pressed));
  b.append(icon(name, size === 'sm' ? 'sm' : ''));
  return b;
}

// A segmented control: options are [value, label, title?]; onPick(value).
export function seg(options, value, onPick, { block = false } = {}) {
  const box = el('div', { class: `seg${block ? ' block' : ''}`, role: 'radiogroup' });
  const set = (v) => { for (const b of box.children) b.classList.toggle('on', b.dataset.v === String(v)); };
  for (const [v, label, title] of options) {
    box.append(el('button', { type: 'button', role: 'radio', 'data-v': String(v), title: title || null, onclick: () => { set(v); onPick(v); } }, label));
  }
  set(value);
  box.set = set;
  return box;
}

// A select: options are [value, label] or {value, label, group, disabled}.
export function select(options, value, onChange, { title = null, width = null } = {}) {
  const s = el('select', { title, style: width ? { width } : null, onchange: () => onChange(s.value) });
  let group = null, groupEl = null;
  for (const o of options) {
    const opt = Array.isArray(o) ? { value: o[0], label: o[1] } : o;
    if (opt.group && opt.group !== group) { group = opt.group; groupEl = el('optgroup', { label: group }); s.append(groupEl); }
    (opt.group ? groupEl : s).append(el('option', { value: opt.value, disabled: opt.disabled || null, selected: String(opt.value) === String(value) || null }, opt.label));
  }
  return s;
}

// A number with steppers and a unit.
export function numberField({ value, min = null, max = null, step = 1, unit = '', onChange = () => {}, width = null }) {
  const input = el('input', { type: 'number', value, min, max, step, style: width ? { width } : null });
  const clamp = (v) => Math.min(max ?? Infinity, Math.max(min ?? -Infinity, v));
  const commit = () => { const v = clamp(Number(input.value)); if (Number.isFinite(v)) { input.value = v; onChange(v); } };
  input.addEventListener('change', commit);
  // rounded to the step's decimals, so that 0.1 + 0.05 is written 0.15
  const decimals = (String(step).split('.')[1] || '').length;
  const stepBy = (k) => { input.value = +clamp(Number(input.value || 0) + k * step).toFixed(decimals); commit(); };
  const box = el('div', { class: 'unit-input' }, input,
    el('button', { type: 'button', class: 'step', title: 'less', onclick: () => stepBy(-1) }, '−'),
    el('button', { type: 'button', class: 'step', title: 'more', onclick: () => stepBy(1) }, '+'));
  if (unit) box.append(el('span', { class: 'unit' }, unit));
  box.input = input;
  return box;
}

// A slider with its value written beside it.
export function slider({ value, min, max, step = 1, format = (v) => String(v), onInput = null, onChange = () => {} }) {
  const out = el('output', {}, format(value));
  const input = el('input', { type: 'range', min, max, step, value });
  const paint = () => { input.style.setProperty('--fill', `${(100 * (Number(input.value) - min)) / (max - min)}%`); out.textContent = format(Number(input.value)); };
  input.addEventListener('input', () => { paint(); onInput?.(Number(input.value)); });
  input.addEventListener('change', () => onChange(Number(input.value)));
  paint();
  const box = el('div', { class: 'slider' }, input, out);
  box.input = input;
  return box;
}

export function toggle(label, checked, onChange, { title = null } = {}) {
  const input = el('input', { type: 'checkbox', checked: checked || null, onchange: () => onChange(input.checked) });
  const box = el('label', { class: 'switch', title }, input, el('span', { class: 'track' }), el('span', {}, label));
  box.setValue = (v) => { input.checked = !!v; };
  box.input = input;
  return box;
}

// A labelled field: fld('gap open', control, {help, hint, inline}).
export function fld(label, control, { help = null, hint = null, inline = false } = {}) {
  return el('div', { class: `fld${inline ? ' inline' : ''}` },
    el('label', { class: 'lab' }, label, help ? el('span', { class: 'help' }, help) : null),
    control, hint ? el('div', { class: 'hint' }, hint) : null);
}

// ---- work in progress ---------------------------------------------------------------------------------
// A line saying what is being computed, with the seconds it has taken so far; its timer stops once
// the line has left the page.
export function working(text) {
  const t0 = performance.now(), secs = el('span', { class: 'elapsed' }, '0.0 s');
  const box = el('div', { class: 'wait' }, el('span', { class: 'spin' }), el('span', {}, text), secs);
  let seen = false;
  const timer = setInterval(() => {
    if (box.isConnected) seen = true;
    else if (seen) { clearInterval(timer); return; }
    secs.textContent = `${((performance.now() - t0) / 1000).toFixed(1)} s`;
  }, 100);
  return box;
}

// A run button while its job runs: disabled, with a spinner and what it is doing.
export function busy(b, on, label) {
  if (!b) return;
  if (on) { b.keep = b.keep || [...b.childNodes]; b.disabled = true; b.replaceChildren(el('span', { class: 'spin' }), el('span', {}, label)); }
  else if (b.keep) { b.replaceChildren(...b.keep); b.keep = null; b.disabled = false; }
}

// Enter anywhere in an instrument's card runs it (not in a text area, nor on a button).
export function onEnter(box, run) {
  box.addEventListener('keydown', (e) => {
    if (e.key !== 'Enter' || e.target.matches('textarea, button')) return;
    if (e.target.tagName === 'INPUT') e.target.dispatchEvent(new Event('change'));
    e.preventDefault();
    run();
  });
}

// ---- sequences and files -----------------------------------------------------------------------------
// A text of letters with what it holds written under it: its length, its letters' kind, its GC;
// FASTA is read (the first record, or all with `multi`). Buttons: a file, the library, clear.
export function seqEditor({ value = '', placeholder = 'paste a sequence, FASTA, or drop a file here', rows = 4, onChange = () => {}, multi = false, library = null, alphabet = 'auto' } = {}) {
  const ta = el('textarea', { rows, spellcheck: false, placeholder });
  ta.value = value;
  const stats = el('span', { class: 'stats' });
  const describe = () => {
    const parsed = parseSequences(ta.value);
    const n = parsed.reduce((a, r) => a + r.seq.length, 0);
    if (!n) { stats.textContent = 'empty'; return parsed; }
    const all = parsed.map((r) => r.seq).join('');
    const kind = alphabet !== 'auto' ? alphabet : guessAlphabet(all);
    const gc = kind === 'dna' ? ` · GC ${((100 * (all.match(/[GCgc]/g) || []).length) / n).toFixed(1)}%` : '';
    stats.textContent = `${parsed.length > 1 ? `${parsed.length} sequences · ` : ''}${n.toLocaleString('en-US')} ${kind === 'protein' ? 'residues' : 'bases'}${gc}`;
    return parsed;
  };
  ta.addEventListener('input', () => { describe(); onChange(ta.value); });
  const pick = el('input', { type: 'file', hidden: true, onchange: async () => { const f = pick.files[0]; if (f) { ta.value = await readFile(f, 8 * 1024 * 1024); describe(); onChange(ta.value); } pick.value = ''; } });
  const acts = el('span', { class: 'acts' },
    iconBtn('upload', { title: 'read a file', size: 'sm', onClick: () => pick.click() }),
    library ? iconBtn('data', { title: 'from the library', size: 'sm', onClick: (e) => library(e.currentTarget, (text) => { ta.value = text; describe(); onChange(ta.value); }) }) : null,
    iconBtn('x', { title: 'clear', size: 'sm', onClick: () => { ta.value = ''; describe(); onChange(''); } }));
  const box = el('div', { class: 'seq-editor' }, ta, el('div', { class: 'bar' }, stats, acts), pick);
  box.addEventListener('dragover', (e) => { e.preventDefault(); e.stopPropagation(); });
  box.addEventListener('drop', async (e) => {
    const f = e.dataTransfer?.files?.[0];
    if (!f) return;
    e.preventDefault(); e.stopPropagation();
    ta.value = await readFile(f, 8 * 1024 * 1024);
    describe(); onChange(ta.value);
  });
  describe();
  box.textarea = ta;
  box.value = () => ta.value;
  box.set = (v) => { ta.value = v; describe(); };
  box.sequences = () => parseSequences(ta.value).filter((r) => r.seq.length);
  return box;
}

// FASTA records, or one bare sequence: [{name, seq}] with the letters cleaned of spaces and digits.
export function parseSequences(text) {
  const t = String(text || '');
  if (!t.trim()) return [];
  if (!t.trimStart().startsWith('>')) return [{ name: 'sequence', seq: t.replace(/[^A-Za-z*-]/g, '').toUpperCase() }];
  const out = [];
  let cur = null;
  for (const line of t.split(/\r?\n/)) {
    if (line.startsWith('>')) { cur = { name: line.slice(1).trim().split(/\s+/)[0] || `seq${out.length + 1}`, seq: '' }; out.push(cur); }
    else if (cur) cur.seq += line.replace(/[^A-Za-z*-]/g, '').toUpperCase();
  }
  return out;
}

export function guessAlphabet(s) {
  const sample = s.slice(0, 5000).toUpperCase();
  if (!sample) return 'dna';
  const nuc = (sample.match(/[ACGTUN]/g) || []).length;
  return nuc / sample.length > 0.9 ? 'dna' : 'protein';
}

export function dropzone({ label = 'Drop files here, or choose them', hint = 'FASTA · FASTQ · GenBank · VCF · Newick · SAM, plain or gzipped', multiple = true, onFiles }) {
  const pick = el('input', { type: 'file', multiple: multiple || null, hidden: true, onchange: () => { if (pick.files.length) onFiles([...pick.files]); pick.value = ''; } });
  const z = el('div', { class: 'dropzone', role: 'button', tabindex: 0, onclick: () => pick.click(), onkeydown: (e) => { if (e.key === 'Enter' || e.key === ' ') pick.click(); } },
    icon('upload', 'lg'), el('div', {}, el('b', {}, label), el('div', { class: 'small' }, hint)), pick);
  z.addEventListener('dragover', (e) => { e.preventDefault(); e.stopPropagation(); z.classList.add('over'); });
  z.addEventListener('dragleave', () => z.classList.remove('over'));
  z.addEventListener('drop', (e) => { e.preventDefault(); e.stopPropagation(); z.classList.remove('over'); if (e.dataTransfer.files.length) onFiles([...e.dataTransfer.files]); });
  return z;
}

// The text of a file (gzip unpacked by the browser where it can), at most `limit` bytes.
export async function readFile(file, limit = Infinity) {
  let blob = file.size > limit ? file.slice(0, limit) : file;
  if (/\.gz$/i.test(file.name) && 'DecompressionStream' in window) {
    blob = await new Response(file.stream().pipeThrough(new DecompressionStream('gzip'))).blob();
    if (blob.size > limit) blob = blob.slice(0, limit);
  }
  return blob.text();
}

// ---- surfaces -----------------------------------------------------------------------------------------
export function card({ title = null, sub = null, actions = [], body = [], pad = false, cls = '' } = {}) {
  const c = el('div', { class: `card ${cls}`.trim() });
  if (title || actions.length) c.append(el('div', { class: 'card-head' }, el('div', {}, title ? el('h3', {}, title) : null, sub ? el('div', { class: 'sub' }, sub) : null), el('div', { class: 'actions' }, ...actions)));
  const b = el('div', { class: pad ? 'card-body' : '' }, ...[].concat(body));
  c.append(b);
  c.body = b;
  return c;
}

export function tabs(options, value, onPick) {
  const box = el('div', { class: 'tabs', role: 'tablist' });
  const set = (v) => { for (const b of box.children) b.classList.toggle('on', b.dataset.v === String(v)); };
  for (const [v, label] of options) box.append(el('button', { role: 'tab', 'data-v': String(v), onclick: () => { set(v); onPick(v); } }, label));
  set(value);
  box.set = set;
  return box;
}

export function metrics(items) {
  return el('div', { class: 'metrics' }, ...items.filter(Boolean).map(([v, k, cls]) => el('div', { class: `metric ${cls || ''}` }, el('div', { class: 'v' }, v), el('div', { class: 'k' }, k))));
}

export function callout(text, kind = '', ic = null) {
  return el('div', { class: `callout ${kind}`.trim() }, icon(ic || (kind === 'bad' ? 'warn' : kind === 'warn' ? 'warn' : 'info'), 'sm'), el('div', {}, text));
}

export function empty({ icon: ic = 'spark', title, text = '', action = null }) {
  return el('div', { class: 'empty' }, el('div', {}, icon(ic), el('h3', {}, title), text ? el('p', {}, text) : null, action));
}

// A table: columns [{key, label, num, format, sort}], rows of objects; sortable by a click.
export function table(columns, rows, { wide = true, max = 2000 } = {}) {
  let sortKey = null, dir = 1;
  const t = el('table', { class: `data${wide ? ' wide' : ''}` });
  const draw = () => {
    clear(t);
    t.append(el('tr', {}, ...columns.map((c) => el('th', { class: `${c.num ? 'num' : ''} sortable`, onclick: () => { dir = sortKey === c.key ? -dir : 1; sortKey = c.key; draw(); } }, c.label, sortKey === c.key ? (dir > 0 ? ' ↑' : ' ↓') : ''))));
    let rs = rows;
    if (sortKey) {
      const col = columns.find((c) => c.key === sortKey);
      const val = col.sort || ((r) => r[sortKey]);
      rs = [...rows].sort((a, b) => { const x = val(a), y = val(b); return (x > y ? 1 : x < y ? -1 : 0) * dir; });
    }
    for (const r of rs.slice(0, max)) t.append(el('tr', {}, ...columns.map((c) => { const v = c.format ? c.format(r[c.key], r) : r[c.key]; return el('td', { class: c.num ? 'num' : '' }, v instanceof Node ? v : String(v ?? '')); })));
  };
  draw();
  return rows.length > max ? el('div', {}, t, el('p', { class: 'note' }, `the first ${max.toLocaleString('en-US')} of ${rows.length.toLocaleString('en-US')} rows`)) : t;
}

export function kv(pairs) {
  const d = el('dl', { class: 'kv' });
  for (const [k, v] of pairs) if (v !== null && v !== undefined) d.append(el('dt', {}, k), el('dd', {}, v instanceof Node ? v : String(v)));
  return d;
}

// ---- overlays -----------------------------------------------------------------------------------------
// A popover beside an anchor; closes on a click elsewhere or Escape. Returns {box, close}.
export function popover(anchor, content, { width = null, place = 'below' } = {}) {
  const box = el('div', { class: 'popover', style: width ? { width: `${width}px` } : null }, content);
  document.body.append(box);
  const r = anchor.getBoundingClientRect();
  const w = box.offsetWidth, h = box.offsetHeight;
  let x = Math.min(r.left, window.innerWidth - w - 10), y = place === 'above' ? r.top - h - 6 : r.bottom + 6;
  if (y + h > window.innerHeight - 10) y = Math.max(10, r.top - h - 6);
  box.style.left = `${Math.max(10, x)}px`;
  box.style.top = `${y}px`;
  const close = () => { box.remove(); document.removeEventListener('pointerdown', away, true); document.removeEventListener('keydown', esc_); };
  const away = (e) => { if (!box.contains(e.target) && !anchor.contains(e.target)) close(); };
  const esc_ = (e) => { if (e.key === 'Escape') close(); };
  setTimeout(() => { document.addEventListener('pointerdown', away, true); document.addEventListener('keydown', esc_); });
  return { box, close };
}

// A menu of items [{label, icon, meta, onClick} | 'sep' | {group}].
export function menu(anchor, items, opts = {}) {
  const list = el('div');
  const p = popover(anchor, list, opts);
  for (const it of items) {
    if (it === 'sep') { list.append(el('div', { class: 'menu-sep' })); continue; }
    if (it.group) { list.append(el('div', { class: 'menu-group' }, it.group)); continue; }
    list.append(el('button', { class: 'menu-item', disabled: it.disabled || null, onclick: () => { p.close(); it.onClick?.(); } }, it.icon ? icon(it.icon, 'sm') : null, el('span', {}, it.label), it.meta ? el('span', { class: 'meta' }, it.meta) : null));
  }
  return p;
}

export function dialog({ title, body, foot = [], width = null, onClose = null, cls = '' }) {
  const veil = el('div', { class: 'dialog-veil' });
  const box = el('div', { class: `dialog ${cls}`.trim(), role: 'dialog', 'aria-modal': 'true', style: width ? { width: `min(${width}px, 94vw)` } : null });
  const close = () => { veil.remove(); document.removeEventListener('keydown', key); onClose?.(); };
  const key = (e) => { if (e.key === 'Escape') close(); };
  if (title) box.append(el('div', { class: 'dialog-head' }, el('h2', {}, title), el('span', { style: { flex: 1 } }), iconBtn('x', { title: 'close', onClick: close })));
  box.append(el('div', { class: 'dialog-body' }, body));
  if (foot.length) box.append(el('div', { class: 'dialog-foot' }, ...foot));
  veil.append(box);
  veil.addEventListener('pointerdown', (e) => { if (e.target === veil) close(); });
  document.addEventListener('keydown', key);
  document.body.append(veil);
  return { box, close };
}

// A view's head: its title and tagline, an "about" toggle for the long text, and its actions.
export function vhead({ title, tagline = '', about = '', actions = [] }) {
  const aboutBox = el('div', { class: 'about', hidden: true }, about);
  const h = el('div', { class: 'vhead' }, el('div', { class: 'titles' }, el('h2', {}, title), tagline ? el('p', { class: 'tagline' }, tagline) : null),
    el('div', { class: 'actions' }, ...actions, about ? iconBtn('info', { title: 'about this view', onClick: (e) => { aboutBox.hidden = !aboutBox.hidden; e.currentTarget.setAttribute('aria-pressed', String(!aboutBox.hidden)); } }) : null));
  return [h, aboutBox];
}

// ---- helpers -------------------------------------------------------------------------------------------
export function download(name, text, type = 'text/plain') {
  const a = el('a', { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
}

export async function copy(text) {
  try { await navigator.clipboard.writeText(text); toast('copied', 'good', 1800); } catch { toast('the browser would not copy', 'error'); }
}

export function fmtBytes(b) {
  if (!Number.isFinite(b)) return '–';
  const u = ['B', 'KB', 'MB', 'GB'];
  let k = 0;
  while (b >= 1024 && k < u.length - 1) { b /= 1024; k++; }
  return `${b.toFixed(k === 0 || b >= 100 ? 0 : 1)} ${u[k]}`;
}

export { esc, toast };
