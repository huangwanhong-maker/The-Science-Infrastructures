// What a tool's result is drawn with, shared by every tool: sequences, multiple alignments, tables,
// trees, files to save, figures. A tool's render(result, values, out) picks from these; out is the
// element the result goes in.

import { el, esc, hover } from '../dom.js';
import { icon } from './icons.js';
import { btn, iconBtn, table, download, copy, card } from './kit.js';
import { token, BASE_COLOURS } from '../palette.js';
import { treeSvg } from '../views/trees.js';
import * as library from '../core/library.js';

// Amino acids coloured by their chemistry, in the workbench's hand.
const PROTEIN_GROUPS = [
  ['AVLIMC', '#C9A043'],   // hydrophobic
  ['FWY', '#B07D62'],      // aromatic
  ['KRH', '#5D7FB2'],      // positive
  ['DE', '#C55A72'],       // negative
  ['STNQ', '#6E9C7A'],     // polar
  ['G', '#8C7BB8'], ['P', '#A8577E'],
];
export function residueColour(c, alphabet) {
  if (c === '-' || c === '.') return null;
  if (alphabet === 'dna') return BASE_COLOURS()[c] || token('--N');
  for (const [set, colour] of PROTEIN_GROUPS) if (set.includes(c)) return colour;
  return token('--N');
}

// Files to save; and the first of a kind another tool takes can be kept in Your data instead,
// so that one tool's result becomes the next one's input.
const KEEPABLE = /\.(fasta|fa|fna|faa|fastq|fq|vcf|nwk|newick|sam|gb|gbk)$/i;
export function downloads(items) {
  const list = items.filter(Boolean);
  const keep = list.find(([, name]) => KEEPABLE.test(name));
  const text = (t) => (typeof t === 'function' ? t() : t);
  return el('div', { class: 'row', style: { gap: '6px' } },
    ...list.map(([label, name, t, type]) => btn(label, { icon: 'download', size: 'sm', onClick: () => download(name, text(t), type || 'text/plain') })),
    keep ? btn('keep', { icon: 'plus', size: 'sm', kind: 'ghost', title: `keep ${keep[1]} in Your data, for other tools`, onClick: () => library.add([new File([text(keep[2])], keep[1], { type: 'text/plain' })]) }) : null);
}

export const fasta = (records, width = 70) => records.map((r) => `>${r.name}\n${(r.seq.match(new RegExp(`.{1,${width}}`, 'g')) || ['']).join('\n')}`).join('\n') + '\n';

// Named sequences, each in a block with its positions, to copy or save.
export function sequences(records, { alphabet = 'dna', title = 'Sequences', file = 'sequences.fasta' } = {}) {
  const body = el('div', { class: 'stack' });
  for (const r of records.slice(0, 200)) {
    const block = el('div', { class: 'seq-block' });
    const s = r.seq || '';
    const width = 60;
    const lines = [];
    for (let i = 0; i < s.length && i < 600000; i += width) lines.push(`<span class="pos">${String(i + 1).padStart(8, ' ')}  </span>${esc(s.slice(i, i + width))}`);
    block.innerHTML = lines.join('\n') || '<span class="faint">empty</span>';
    body.append(el('div', {}, el('div', { class: 'row', style: { marginBottom: '4px' } }, el('b', { style: { fontWeight: 500 } }, r.name), el('span', { class: 'small faint' }, `${s.length.toLocaleString('en-US')} ${alphabet === 'protein' ? 'residues' : 'bases'}${r.note ? ` · ${r.note}` : ''}`),
      el('span', { style: { flex: 1 } }), iconBtn('copy', { title: 'copy', size: 'sm', onClick: () => copy(s) })), block));
  }
  return card({ title, sub: `${records.length} sequence${records.length === 1 ? '' : 's'}`, actions: [downloads([['FASTA', file, () => fasta(records)]])], body, pad: true });
}

// A multiple alignment: names fixed at the left, columns coloured, a ruler and a line of
// conservation or consensus underneath, scrolled sideways.
export function alignmentView(rows, { alphabet = 'dna', extra = [], title = 'The alignment', file = 'alignment' } = {}) {
  const n = rows[0]?.row.length || 0;
  const cw = 9, ch = 15, nameW = Math.min(180, 8 + 7 * Math.max(...rows.map((r) => r.name.length), 8));
  const all = [...rows, ...extra];
  const canvas = el('canvas', { width: n * cw * 2, height: (all.length + 1) * ch * 2, style: { width: `${n * cw}px`, height: `${(all.length + 1) * ch}px`, display: 'block' } });
  const g = canvas.getContext('2d');
  g.scale(2, 2);
  g.font = `11.5px ${getComputedStyle(document.body).getPropertyValue('--mono') || 'monospace'}`;
  g.textAlign = 'center';
  g.textBaseline = 'middle';
  // the ruler
  g.fillStyle = token('--ink-faint');
  for (let j = 9; j < n; j += 10) g.fillText(String(j + 1), j * cw + cw / 2, ch / 2);
  all.forEach((r, i) => {
    const y = (i + 1) * ch;
    for (let j = 0; j < n; j++) {
      const c = r.row[j];
      const col = r.plain ? null : residueColour(c, alphabet);
      if (col) { g.fillStyle = col; g.globalAlpha = 0.22; g.fillRect(j * cw, y + 1, cw, ch - 2); g.globalAlpha = 1; }
      g.fillStyle = col && !r.plain ? token('--ink') : token('--ink-faint');
      g.fillText(c, j * cw + cw / 2, y + ch / 2 + 0.5);
    }
  });
  canvas.addEventListener('pointermove', (e) => {
    const b = canvas.getBoundingClientRect();
    const j = Math.floor((e.clientX - b.left) / cw), i = Math.floor((e.clientY - b.top) / ch) - 1;
    if (i < 0 || i >= all.length || j < 0 || j >= n) { hover(null); return; }
    const r = all[i];
    const pos = r.row.slice(0, j + 1).replace(/[-.]/g, '').length;
    hover(`<b>${esc(r.name)}</b> · column ${j + 1}${r.plain ? '' : ` · residue ${pos}`}<br>${esc(r.row[j])}`, e);
  });
  canvas.addEventListener('pointerleave', () => hover(null));
  const names = el('div', { style: { width: `${nameW}px`, flex: 'none', font: '12px var(--mono)', paddingTop: `${ch}px` } },
    ...all.map((r) => el('div', { style: { height: `${ch}px`, lineHeight: `${ch}px`, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis', color: r.plain ? 'var(--ink-faint)' : 'var(--ink)', fontStyle: r.plain ? 'italic' : 'normal' }, title: r.name }, r.name)));
  const scroller = el('div', { style: { overflowX: 'auto', flex: 1, minWidth: 0 } }, canvas);
  const text = () => fasta(rows.map((r) => ({ name: r.name, seq: r.row })), 60);
  return card({ title, sub: `${rows.length} sequences × ${n.toLocaleString('en-US')} columns`, actions: [downloads([['aligned FASTA', `${file}.fasta`, text]])],
    body: el('div', { style: { display: 'flex', gap: '8px', padding: '12px 14px' } }, names, scroller) });
}

// Rows of a table, sortable, with a CSV to save.
export function tableCard(title, columns, rows, { sub = null, file = 'table.csv' } = {}) {
  const csv = () => [columns.map((c) => c.label).join(','), ...rows.map((r) => columns.map((c) => { const v = r[c.key]; const s = String(v ?? ''); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; }).join(','))].join('\n') + '\n';
  return card({ title, sub: sub ?? `${rows.length.toLocaleString('en-US')} rows`, actions: [downloads([['CSV', file, csv]])], body: el('div', { style: { overflow: 'auto', maxHeight: '520px', padding: '4px 14px 10px' } }, table(columns, rows)) });
}

// A text file: its first lines shown, and all of it to save.
export function textCard(title, name, text, { sub = null, lines = 40 } = {}) {
  const shown = text.split('\n').slice(0, lines).join('\n');
  return card({ title, sub: sub ?? `${text.length.toLocaleString('en-US')} characters`, actions: [downloads([['save', name, text]]), iconBtn('copy', { title: 'copy', size: 'sm', onClick: () => copy(text) })],
    body: el('pre', { class: 'file', style: { border: 0, borderRadius: 0, maxHeight: '340px' } }, shown + (text.split('\n').length > lines ? '\n…' : '')) });
}

// A tree with lengths and supports, rooted at its left.
export function treeCard(title, nodes, { newick = null, scaleLabel = 'substitutions a site', width = 760 } = {}) {
  const fig = el('div', { style: { padding: '10px 14px', overflowX: 'auto' } });
  fig.innerHTML = treeSvg(nodes, { width, scaleLabel }).svg;
  // drawn again at the card's own width once it is on the page, so that a narrow card is not cut
  requestAnimationFrame(() => {
    const w = fig.clientWidth - 28;
    if (fig.isConnected && w > 0 && Math.abs(w - width) > 20) fig.innerHTML = treeSvg(nodes, { width: Math.max(420, w), scaleLabel }).svg;
  });
  return card({ title, sub: `${nodes.filter((n) => !nodes.some((m) => m.parent === nodes.indexOf(n))).length} leaves`, actions: newick ? [downloads([['Newick', 'tree.nwk', newick]])] : [], body: fig });
}

// A figure of the workbench's own drawing (an SVG string, or an element).
export function figure(title, content, { sub = null, actions = [] } = {}) {
  const body = el('div', { style: { padding: '10px 14px', overflowX: 'auto' } });
  if (typeof content === 'string') body.innerHTML = content; else body.append(content);
  return card({ title, sub, actions, body });
}

export { icon };
