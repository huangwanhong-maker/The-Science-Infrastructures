// A tool's inputs, declared as a schema, made into a form. Each input is
//
//   {key, type, label, help, hint, default, when(values), ...}
//
// with the types:
//   sequence     one sequence, pasted or read from a file          (alphabet: 'dna' | 'protein' | 'auto')
//   sequences    several (FASTA), pasted, or a file of your data   (kinds: the file kinds it takes)
//   file         a file of your data                               (kinds)
//   genome       a genome of the catalogue or of your files
//   number       {min, max, step, unit}
//   range        {min, max, step, format}
//   choice       a few options shown at once                       (options: [value, label][])
//   select       many options                                       (options)
//   switch       on or off
//   text         a line of text
// A schema may group inputs ({group: 'title'} entries) and hide one until another has a value
// (when). form() returns {element, values(), set(values), validate()}.

import { el } from '../dom.js';
import { icon } from './icons.js';
import { seqEditor, numberField, slider, seg, select, toggle, fld, btn, parseSequences } from './kit.js';
import * as library from '../core/library.js';

export function form(schema, ctx, initial = {}, onChange = () => {}) {
  const values = {};
  const rows = [];
  const element = el('div', { class: 'tool-form' });
  const notify = () => { for (const r of rows) if (r.input.when) r.node.hidden = !r.input.when(values); onChange(values); };
  let group = null;
  for (const input of schema) {
    if (input.group) { group = el('div', { class: 'eyebrow', style: { marginTop: element.children.length ? '6px' : '0' } }, input.group); element.append(group); continue; }
    const v0 = initial[input.key] ?? (typeof input.default === 'function' ? input.default(ctx) : input.default);
    values[input.key] = v0;
    const set = (v) => { values[input.key] = v; notify(); };
    const control = make(input, v0, set, ctx);
    const node = fld(input.label, control, { help: input.help || null, hint: input.hint || null, inline: input.type === 'switch' });
    // Numbers stand two to a row; anything may ask for half a row, or a whole one, by `half`.
    if (input.half ?? input.type === 'number') node.classList.add('half');
    rows.push({ input, node, control });
    element.append(node);
  }
  notify();
  return {
    element,
    values: () => ({ ...values }),
    set(vs) { for (const r of rows) if (r.input.key in vs && r.control.setValue) { r.control.setValue(vs[r.input.key]); values[r.input.key] = vs[r.input.key]; } notify(); },
    validate() {
      for (const r of rows) {
        if (r.node.hidden || !r.input.required) continue;
        const v = values[r.input.key];
        const emptyText = typeof v === 'string' ? !v.trim() : v && typeof v === 'object' && 'text' in v ? !String(v.text).trim() : v === null || v === undefined || v === '';
        if (emptyText) return `${r.input.label} is needed`;
        if (r.input.type === 'sequences' && v?.text !== undefined && r.input.min && parseSequences(v.text).length < r.input.min) return `${r.input.label}: at least ${r.input.min} sequences`;
      }
      return null;
    },
  };
}

function make(input, v0, set, ctx) {
  switch (input.type) {
    case 'sequence': {
      const ed = seqEditor({ value: v0 || '', rows: input.rows || 4, alphabet: input.alphabet || 'auto', placeholder: input.placeholder, onChange: set,
        library: (anchor, use) => library.pick(anchor, ['fasta'], async (r) => use(await library.text(r.id))) });
      ed.setValue = (v) => ed.set(v || '');
      if (input.samples) ed.append(samples(input.samples, (s) => { ed.set(s); set(s); }));
      return ed;
    }
    case 'sequences': return sourceInput(input, v0, set, true);
    case 'file': return filePicker(input, v0, set);
    case 'genome': {
      const opts = ctx.catalogue.filter((g) => g.fetched).map((g) => ({ value: g.id, label: g.title.replace(/, complete genome.*$/, '').replace(/ \(.*\)$/, ''), group: { bacteria: 'Bacteria', viruses: 'Viruses', mitochondria: 'Mitochondria', simulated: 'Simulated', files: 'Your files' }[g.group] || g.group }));
      const v = v0 || ctx.genome;
      set(v);
      const s = select(opts, v, set);
      s.setValue = (x) => { s.value = x; };
      return s;
    }
    case 'number': {
      const n = numberField({ value: v0, min: input.min, max: input.max, step: input.step || 1, unit: input.unit || '', onChange: set });
      n.setValue = (x) => { n.input.value = x; };
      return n;
    }
    case 'range': {
      const r = slider({ value: v0, min: input.min, max: input.max, step: input.step || 1, format: input.format, onChange: set, onInput: set });
      r.setValue = (x) => { r.input.value = x; r.input.dispatchEvent(new Event('input')); };
      return r;
    }
    case 'choice': {
      const c = seg(input.options, v0, set, { block: input.block !== false && input.options.length <= 4 });
      c.setValue = (x) => c.set(x);
      return c;
    }
    case 'select': {
      const s = select(input.options, v0, (x) => set(input.numeric ? Number(x) : x));
      s.setValue = (x) => { s.value = x; };
      return s;
    }
    case 'switch': {
      const t = toggle(input.text || '', !!v0, set);
      t.setValue = (x) => { t.querySelector('input').checked = !!x; };
      return t;
    }
    default: {
      const t = el('input', { type: 'text', value: v0 ?? '', placeholder: input.placeholder || '', oninput: () => set(t.value) });
      t.setValue = (x) => { t.value = x ?? ''; };
      return t;
    }
  }
}

function samples(list, use) {
  return el('div', { class: 'chips', style: { padding: '6px 8px', borderTop: '1px solid var(--hair)', background: 'var(--surface-2)' } },
    el('span', { class: 'small faint', style: { alignSelf: 'center' } }, 'try:'), ...list.map(([label, text]) => {
      const b = el('button', { class: 'chip', type: 'button', onclick: async () => {
        b.disabled = true;
        try { use(typeof text === 'function' ? await text() : text); } catch (e) { b.title = e.message; b.classList.add('bad'); } finally { b.disabled = false; }
      } }, label);
      return b;
    }));
}

// Sequences pasted, or a file of your data: the value is {text} or {file: id}.
function sourceInput(input, v0, set) {
  let mode = v0?.file ? 'file' : 'paste';
  let text = v0?.text ?? (typeof v0 === 'string' ? v0 : '');
  let file = v0?.file || null;
  const ed = seqEditor({ value: text, rows: input.rows || 6, multi: true, alphabet: input.alphabet || 'auto', placeholder: input.placeholder || 'paste FASTA (>name, then the sequence), or drop a file', onChange: (t) => { text = t; push(); } });
  const pickBox = filePicker({ ...input, kinds: input.kinds || ['fasta'] }, file, (f) => { file = f; push(); });
  const push = () => set(mode === 'file' ? { file } : { text });
  const switcher = seg([['paste', 'paste'], ['file', 'from your data']], mode, (m) => { mode = m; ed.hidden = m !== 'paste'; pickBox.hidden = m !== 'file'; push(); });
  ed.hidden = mode !== 'paste';
  pickBox.hidden = mode !== 'file';
  const box = el('div', { class: 'stack' }, el('div', { class: 'row' }, switcher), ed, pickBox);
  if (input.samples) ed.append(samples(input.samples, (s) => { ed.set(s); text = s; push(); }));
  box.setValue = (v) => { if (v?.file) { mode = 'file'; file = v.file; switcher.set('file'); } else { mode = 'paste'; text = v?.text ?? v ?? ''; ed.set(text); switcher.set('paste'); } ed.hidden = mode !== 'paste'; pickBox.hidden = mode !== 'file'; };
  push();
  return box;
}

// A file of your data, of some kinds: a chip that shows the file chosen and opens the picker.
function filePicker(input, v0, set) {
  let current = v0;
  const chip = el('button', { type: 'button', class: 'btn block', style: { justifyContent: 'flex-start' } });
  const draw = () => {
    const r = current ? library.byId(current) : null;
    chip.replaceChildren(icon('file', 'sm'), el('span', { style: { overflow: 'hidden', textOverflow: 'ellipsis', minWidth: 0 } }, r ? r.name : current ? current : `choose ${input.kinds ? input.kinds.join(' or ').toUpperCase() : 'a file'}…`),
      r ? el('span', { class: 'kind', style: { marginLeft: 'auto' } }, r.kind) : null);
  };
  chip.addEventListener('click', () => library.pick(chip, input.kinds || null, (r) => { current = r.id; draw(); set(r.id); }));
  library.ready().then(() => {
    if (!current && input.autopick !== false) { const r = library.ofKinds(input.kinds)[0]; if (r) { current = r.id; set(r.id); } }
    draw();
  });
  draw();
  chip.setValue = (v) => { current = v; draw(); };
  return chip;
}

export { btn };
