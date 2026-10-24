// Your data: the files you have given the workbench and the genomes it knows. Files are added by
// dropping them anywhere, or here; each is recognised by its content and summarised, and can be
// opened as a genome (FASTA, GenBank), handed to any tool that takes its kind, previewed or
// removed. Below them, the catalogue's genomes: NCBI's, fetched by the recipe, and the simulated
// ones, which carry their truth.

import { el, clear } from '../dom.js';
import { icon } from '../ui/icons.js';
import { btn, dropzone, menu, dialog, empty, fmtBytes, callout } from '../ui/kit.js';
import * as library from '../core/library.js';
import * as registry from '../core/registry.js';
import * as api from '../api.js';

let view = null;

export function mount(section, ctx) {
  view = { ctx, section, files: el('div'), genomes: el('div'), head: el('div') };
  section.append(view.head, view.files, view.genomes);
  library.subscribe(() => { if (!section.hidden) draw(); });
  draw();
}

export function show(ctx) {
  view.ctx = ctx;
  library.refresh();
  if (ctx.params.add) { ctx.setParams({ add: null }); setTimeout(() => view.section.querySelector('.dropzone')?.click(), 50); }
}

function draw() {
  const rs = library.all();
  clear(view.head).append(
    dropzone({ label: 'Add files: drop them here, or anywhere on the page', onFiles: async (fs) => { const added = await library.add(fs); if (added.some((r) => r.genome)) await view.ctx.refreshCatalogue(); } }));
  // your files
  const box = clear(view.files);
  const total = rs.reduce((a, r) => a + (r.bytes || 0), 0);
  box.append(el('div', { class: 'section-head' }, el('h3', {}, 'Your files'), el('span', { class: 'sub' }, 'kept in the data cache on this computer, under their content\'s hash'),
    rs.length ? el('span', { class: 'actions small faint' }, `${rs.length} ${rs.length === 1 ? 'file' : 'files'} · ${fmtBytes(total)}`) : null));
  if (!rs.length) {
    box.append(el('div', { class: 'card' }, empty({ icon: 'upload', title: 'No files yet', text: 'Drop a FASTA, FASTQ, GenBank, VCF, Newick or SAM file anywhere on the page. A genome becomes one of the genomes you can open; reads, variants and trees go to the tools that take them.' })));
  } else {
    const grid = el('div', { class: 'lib-grid' });
    for (const r of rs) grid.append(fileCard(r));
    box.append(grid);
  }
  // the genomes
  const g = clear(view.genomes);
  const cat = view.ctx.catalogue;
  const groups = [['bacteria', 'Bacteria'], ['viruses', 'Viruses'], ['mitochondria', 'Mitochondria'], ['simulated', 'Simulated, with their truth']];
  g.append(el('div', { class: 'section-head' }, el('h3', {}, 'Genomes'), el('span', { class: 'sub' }, 'RefSeq\'s, fetched by the recipe, and genomes simulated from a seed')));
  if (cat.some((x) => !x.fetched)) g.append(callout(el('span', {}, 'Some of NCBI\'s genomes are not fetched yet. ', el('button', { class: 'link', onclick: () => view.ctx.go('run') }, 'Run the recipe'), ' to fetch them (a few megabytes).'), 'warn'));
  for (const [key, label] of groups) {
    const items = cat.filter((x) => x.group === key);
    if (!items.length) continue;
    g.append(el('p', { class: 'eyebrow', style: { marginTop: '16px' } }, label));
    const grid = el('div', { class: 'lib-grid' });
    for (const x of items) {
      grid.append(el('div', { class: `lib-card${x.id === view.ctx.genome ? ' current' : ''}` },
        el('div', { class: 'row', style: { gap: '8px', flexWrap: 'nowrap' } }, icon(key === 'simulated' ? 'spark' : 'genome', 'sm'), el('span', { class: 'title' }, x.title.replace(/, complete genome.*$/, ''))),
        el('div', { class: 'facts' }, [x.accession, x.size ? fmtBytes(x.size) : null, x.fetched ? null : 'not fetched'].filter(Boolean).join(' · ')),
        el('div', { class: 'acts' }, btn(x.id === view.ctx.genome ? 'In hand' : 'Open', { size: 'sm', kind: x.id === view.ctx.genome ? '' : 'primary', disabled: !x.fetched, onClick: () => view.ctx.setGenome(x.id, true) }))));
    }
    g.append(grid);
  }
}

function fileCard(r) {
  const tools = registry.toolsFor(r.kind);
  return el('div', { class: `lib-card${r.id === view.ctx.genome ? ' current' : ''}` },
    el('div', { class: 'row', style: { gap: '8px', flexWrap: 'nowrap' } }, el('span', { class: 'kind' }, r.kind), el('span', { class: 'title', title: r.name }, r.name)),
    el('div', { class: 'facts' }, library.describe(r)),
    el('div', { class: 'facts' }, `${fmtBytes(r.bytes)} · added ${new Date(r.added).toLocaleString()}`),
    el('div', { class: 'acts' },
      r.genome ? btn('Open as genome', { size: 'sm', kind: 'primary', icon: 'genome', onClick: async () => { await view.ctx.refreshCatalogue(); view.ctx.setGenome(r.id, true); } }) : null,
      tools.length ? btn('Use in…', { size: 'sm', icon: 'toolbox', onClick: (e) => menu(e.currentTarget, tools.map((t) => ({ label: t.title, icon: t.icon || 'toolbox', meta: t.summary?.slice(0, 40), onClick: () => view.ctx.openTool(t.id, presetFor(t, r)) }))) }) : null,
      btn('Preview', { size: 'sm', kind: 'ghost', onClick: () => preview(r) }),
      btn('', { size: 'sm', kind: 'ghost', icon: 'trash', title: 'remove', onClick: () => remove(r) })));
}

// The values that open a tool on this file: the first input that takes its kind.
function presetFor(t, r) {
  for (const i of t.inputs || []) {
    if (i.type === 'genome' && r.genome) return { [i.key]: r.id };
    if ((i.kinds || []).includes(r.kind)) return { [i.key]: i.type === 'sequences' ? { file: r.id } : r.id };
  }
  return {};
}

async function preview(r) {
  const pre = el('pre', { class: 'file', style: { maxHeight: '52vh' } }, 'Reading…');
  dialog({ title: r.name, body: el('div', { class: 'stack' }, el('div', { class: 'small faint' }, `${r.about} · ${library.describe(r)}`), pre), width: 860 });
  try { pre.textContent = (await api.fileHead(r.id)).head || '(empty)'; } catch (e) { pre.textContent = e.message; }
}

function remove(r) {
  const d = dialog({ title: 'Remove this file?', body: el('p', {}, `${r.name} will be removed from the data cache. Results already made from it stay where they are.`),
    foot: [btn('Keep it', { onClick: () => d.close() }), btn('Remove', { kind: 'primary', icon: 'trash', onClick: async () => { d.close(); await library.remove(r.id); await view.ctx.refreshCatalogue(); } })] });
}
