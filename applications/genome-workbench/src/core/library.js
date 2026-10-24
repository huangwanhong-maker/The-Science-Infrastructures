// The user's files, as the page knows them: the server's library (files.py) read once and kept
// current, files added (dropped anywhere, chosen, or pasted), and a picker that offers the files
// of a kind wherever an input takes one.

import * as api from '../api.js';
import { el, toast } from '../dom.js';
import { icon } from '../ui/icons.js';
import { popover, fmtBytes } from '../ui/kit.js';

let records = [];
let loaded = null;
const listeners = new Set();

export function subscribe(f) { listeners.add(f); return () => listeners.delete(f); }
const emit = () => { for (const f of listeners) f(records); };

export async function refresh() {
  try { records = (await api.files()).files || []; } catch { records = []; }
  emit();
  return records;
}
export const ready = () => (loaded ||= refresh());
export const all = () => records;
export const ofKinds = (kinds) => records.filter((r) => !kinds || kinds.includes(r.kind));
export const byId = (id) => records.find((r) => r.id === id);

// Files added, one after another, with a toast for each that follows its progress.
export async function add(files) {
  const added = [];
  for (const f of files) {
    const t = toast(`Adding ${f.name}…`, '', 600000);
    const line = t.querySelector('div');
    try {
      const r = await api.upload(f, (p) => { line.textContent = `Adding ${f.name}: ${Math.round(p * 100)}%`; });
      line.textContent = `${r.name} added: ${r.about}${describe(r) ? `, ${describe(r)}` : ''}.`;
      t.classList.add('good');
      added.push(r);
    } catch (e) {
      line.textContent = `${f.name}: ${e.message}`;
      t.classList.add('error');
    }
    setTimeout(() => t.remove(), 5000);
  }
  await refresh();
  return added;
}

// The whole text of a small file (a FASTA of a few sequences, a tree).
export const text = async (id) => (await api.fileText(id)).text;

export async function remove(id) {
  await api.removeFile(id);
  await refresh();
}

// What is in a file, in a few words.
export function describe(r) {
  const s = r.summary || {};
  const n = (x) => Number(x || 0).toLocaleString('en-US');
  switch (r.kind) {
    case 'fasta': return `${n(s.sequences)} sequence${s.sequences === 1 ? '' : 's'}, ${n(s.bases)} bases${s.gc ? `, GC ${(s.gc * 100).toFixed(1)}%` : ''}`;
    case 'fastq': return `${n(s.reads)} reads, ${n(s.bases)} bases (mean ${s.mean_length})`;
    case 'genbank': return `${n(s.records)} record${s.records === 1 ? '' : 's'}${s.organism ? `, ${s.organism}` : ''}`;
    case 'vcf': return `${n(s.samples)} samples, ${n(s.sites)} sites`;
    case 'newick': return `${n(s.trees)} tree${s.trees === 1 ? '' : 's'}, ${n(s.leaves)} leaves`;
    case 'sam': return `${n(s.records)} records`;
    default: return fmtBytes(r.bytes);
  }
}

// A picker beside an anchor: the library's files of some kinds, and a way to add one.
export function pick(anchor, kinds, onPick) {
  const list = el('div', { style: { maxHeight: '360px', overflowY: 'auto' } });
  const p = popover(anchor, list, { width: 380 });
  const fill = () => {
    list.replaceChildren();
    const rs = ofKinds(kinds);
    list.append(el('div', { class: 'menu-group' }, rs.length ? 'your files' : 'no files of this kind yet'));
    for (const r of rs) {
      list.append(el('button', { class: 'menu-item', onclick: () => { p.close(); onPick(r); } }, icon('file', 'sm'),
        el('span', { style: { minWidth: 0, overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' } }, r.name), el('span', { class: 'meta' }, `${r.kind} · ${fmtBytes(r.bytes)}`)));
    }
    const chooser = el('input', { type: 'file', multiple: true, hidden: true, onchange: async () => {
      const added = await add([...chooser.files]);
      const fit = added.filter((r) => !kinds || kinds.includes(r.kind));
      if (fit.length) { p.close(); onPick(fit[0]); } else fill();
    } });
    list.append(el('div', { class: 'menu-sep' }), el('button', { class: 'menu-item', onclick: () => chooser.click() }, icon('upload', 'sm'), el('span', {}, 'add a file from this computer…')), chooser);
  };
  ready().then(fill);
  fill();
  return p;
}
