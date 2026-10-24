// The command palette (Ctrl+K, or /): one box that finds anything the workbench offers, its views,
// its tools, the genomes and the user's files, and its commands; and that reads what is typed as
// an action where it can (a range of the genome, a sequence to search, a gene's name).

import { el, esc } from '../dom.js';
import { icon } from '../ui/icons.js';
import { dialog } from '../ui/kit.js';
import * as registry from './registry.js';
import * as library from './library.js';

// A fuzzy match of q in s: its score (higher is better) and the positions matched, or null.
export function match(q, s) {
  q = q.toLowerCase(); const t = s.toLowerCase();
  if (!q) return { score: 0, at: [] };
  const direct = t.indexOf(q);
  if (direct >= 0) return { score: 100 - direct + (direct === 0 || /\W/.test(t[direct - 1]) ? 40 : 0) + q.length * 2, at: [...Array(q.length).keys()].map((k) => direct + k) };
  let i = 0, score = 0, last = -2;
  const at = [];
  for (let j = 0; j < t.length && i < q.length; j++) {
    if (t[j] !== q[i]) continue;
    score += j === last + 1 ? 6 : 1;
    if (j === 0 || /\W/.test(t[j - 1])) score += 5;
    at.push(j); last = j; i++;
  }
  // a scattered match of a longer query is noise, not a find
  return i === q.length && (q.length < 3 || score >= 4 * q.length) ? { score, at } : null;
}

function marked(s, at) {
  const set = new Set(at);
  return [...s].map((c, k) => (set.has(k) ? `<mark>${esc(c)}</mark>` : esc(c))).join('');
}

// The palette over the page. ctx: the shell's context (go, openTool, setGenome, catalogue …).
export function open(ctx, initial = '') {
  const input = el('input', { placeholder: 'Go to a view, a tool, a genome or a file; a range like 1200-4800; a sequence; a gene…', autocomplete: 'off', spellcheck: false });
  const list = el('div', { class: 'palette-list' });
  const body = el('div', {}, el('div', { class: 'palette-input' }, icon('search'), input, el('kbd', {}, 'Esc')), list,
    el('div', { class: 'palette-foot' }, el('span', {}, el('kbd', {}, '↑'), ' ', el('kbd', {}, '↓'), ' to move'), el('span', {}, el('kbd', {}, 'Enter'), ' to choose'), el('span', { style: { marginLeft: 'auto' } }, 'Ctrl K anywhere')));
  const d = dialog({ title: null, body, cls: 'palette' });
  d.box.querySelector('.dialog-body').style.padding = '0';
  let items = [], active = 0;

  const sources = () => {
    const out = [];
    for (const v of registry.views) out.push({ group: 'Views', icon: v.icon, title: v.title, meta: v.tagline || '', run: () => ctx.go(v.id) });
    for (const t of registry.tools) out.push({ group: 'Tools', icon: t.icon || 'toolbox', title: t.title, meta: t.summary || '', keys: (t.keywords || []).join(' '), run: () => ctx.openTool(t.id) });
    for (const g of ctx.catalogue.filter((x) => x.fetched)) out.push({ group: 'Genomes', icon: g.group === 'files' ? 'file' : 'genome', title: g.title.replace(/, complete genome.*$/, ''), meta: g.id.replace(/^ncbi:/, ''), run: () => ctx.setGenome(g.id, true) });
    for (const r of library.all()) out.push({ group: 'Your files', icon: 'file', title: r.name, meta: `${r.kind} · ${library.describe(r)}`, run: () => ctx.go('data', { file: r.id }) });
    for (const c of registry.commands) out.push({ group: 'Commands', icon: c.icon || 'command', title: c.title, meta: c.meta || '', keys: (c.keywords || []).join(' '), run: () => c.run(ctx) });
    return out;
  };

  // What the query itself can mean.
  const readings = (q) => {
    const out = [];
    const range = q.replace(/,/g, '').match(/^(\d+)\s*(?:-|\.\.|–)\s*(\d+)$/);
    if (range && ctx.genome) out.push({ group: 'In this genome', icon: 'genome', title: `Show ${Number(range[1]).toLocaleString('en-US')}..${Number(range[2]).toLocaleString('en-US')}`, meta: 'the genome view', run: () => ctx.go('genome', { begin: Number(range[1]) - 1, end: Number(range[2]) }) });
    if (/^[ACGTNacgtn]{6,}$/.test(q) && ctx.genome) out.push({ group: 'In this genome', icon: 'search', title: `Find ${q.toUpperCase()} on both strands`, meta: 'the genome view', run: () => ctx.go('genome', { find: q.toUpperCase() }) });
    else if (/^[A-Za-z][A-Za-z0-9_.-]{1,24}$/.test(q) && ctx.genome) out.push({ group: 'In this genome', icon: 'search', title: `Find the gene ${q}`, meta: 'by name or locus tag', run: () => ctx.go('genome', { find: q }), gene: true });
    return out;
  };

  const draw = () => {
    const q = input.value.trim();
    const scored = [];
    for (const it of sources()) {
      const m = match(q, `${it.title} ${it.keys || ''}`) || (q.length > 2 ? match(q, it.meta) : null);
      if (m) scored.push({ ...it, score: m.score + (it.group === 'Views' ? 4 : 0), at: m.at.filter((k) => k < it.title.length) });
    }
    scored.sort((a, b) => b.score - a.score);
    // Each group once, the groups in the order of their best find.
    const best = q ? scored.slice(0, 40) : scored.filter((s) => s.group !== 'Genomes' && s.group !== 'Your files').slice(0, 30);
    const order = [...new Set(best.map((s) => s.group))];
    const found = order.flatMap((g) => best.filter((s) => s.group === g));
    // A word read as a gene's name leads only when it looks like one (lacZ, CYTB, ND5, b0001)
    // or nothing else was found; otherwise it waits at the end.
    const rs = readings(q);
    const lead = rs.filter((r) => !r.gene || /[A-Z0-9]/.test(q.slice(1)) || !found.length);
    items = [...lead, ...found, ...rs.filter((r) => !lead.includes(r))];
    active = Math.min(active, Math.max(0, items.length - 1));
    list.replaceChildren();
    let group = null;
    items.forEach((it, k) => {
      if (it.group !== group) { group = it.group; list.append(el('div', { class: 'menu-group' }, group)); }
      const b = el('button', { class: `menu-item${k === active ? ' active' : ''}`, onclick: () => choose(k), onmousemove: () => { if (active !== k) { active = k; paint(); } } }, icon(it.icon || 'spark', 'sm'));
      const t = el('span', {});
      t.innerHTML = it.at?.length ? marked(it.title, it.at) : esc(it.title);
      b.append(t, it.meta ? el('span', { class: 'meta', style: { maxWidth: '45%', overflow: 'hidden', textOverflow: 'ellipsis' } }, it.meta) : null);
      list.append(b);
    });
    if (!items.length) list.append(el('div', { class: 'empty' }, el('p', {}, 'Nothing by that name.')));
  };
  const paint = () => {
    [...list.querySelectorAll('.menu-item')].forEach((b, k) => b.classList.toggle('active', k === active));
    list.querySelectorAll('.menu-item')[active]?.scrollIntoView({ block: 'nearest' });
  };
  const choose = (k) => { const it = items[k]; if (!it) return; d.close(); it.run(); };
  input.addEventListener('input', () => { active = 0; draw(); });
  input.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown') { active = Math.min(items.length - 1, active + 1); paint(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { active = Math.max(0, active - 1); paint(); e.preventDefault(); }
    else if (e.key === 'Enter') { choose(active); e.preventDefault(); }
  });
  input.value = initial;
  library.ready().then(draw);
  draw();
  setTimeout(() => input.focus());
}
