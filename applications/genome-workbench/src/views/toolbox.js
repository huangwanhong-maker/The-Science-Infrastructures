// The Toolbox: every instrument of the workbench, on the user's own sequences and files. The
// catalogue at the left (searchable, by group) comes from the registry; a tool opens with its form
// (built from its inputs, ui/form.js) beside its results, and every run is kept, so that an earlier
// result can be brought back. A new tool joins by being declared (src/tools/), with no change here.

import { el, clear } from '../dom.js';
import { icon } from '../ui/icons.js';
import { btn, callout, empty, card, toast } from '../ui/kit.js';
import { form } from '../ui/form.js';
import * as registry from '../core/registry.js';
import { match } from '../core/palette.js';
import * as fmt from '../format.js';

let view = null;

// The tools used last, remembered in this browser (a convenience: nothing depends on it).
const RECENT = 'genome-workbench:recent-tools';
const recent = () => { try { return JSON.parse(localStorage.getItem(RECENT) || '[]').filter((id) => registry.toolById(id)); } catch { return []; } };
const remember = (id) => { try { localStorage.setItem(RECENT, JSON.stringify([id, ...recent().filter((x) => x !== id)].slice(0, 5))); } catch { /* private window */ } };

export function mount(section, ctx) {
  view = { ctx, section, list: el('div', { class: 'groups' }), pane: el('div', { class: 'tool-pane' }), query: '', current: null, runs: new Map(), forms: new Map() };
  const search = el('input', { type: 'search', placeholder: `Search ${registry.tools.length} tools…`, oninput: () => { view.query = search.value; drawList(); } });
  section.append(el('div', { class: 'toolbox' }, el('div', { class: 'toolbox-list' }, el('div', { class: 'search-row' }, search), view.list), view.pane));
  drawList();
}

export function show(ctx) {
  view.ctx = ctx;
  const pending = ctx.pendingTool;
  ctx.pendingTool = null;
  const id = pending?.id || ctx.params.tool || view.current || registry.tools[0]?.id;
  if (id) open(id, pending?.values || null);
}

function drawList() {
  const box = clear(view.list);
  const q = view.query.trim();
  let shown = 0;
  const item = (t, ic) => el('button', { class: `tool-item${t.id === view.current ? ' on' : ''}`, 'data-id': t.id, onclick: () => { view.ctx.setParams({ tool: t.id }, false); open(t.id); } },
    icon(t.icon || ic, 'sm'), el('div', {}, el('div', { class: 't' }, t.title), t.summary ? el('div', { class: 'd' }, t.summary) : null));
  const used = q ? [] : recent().slice(0, 4).map(registry.toolById);
  if (used.length > 1) {
    box.append(el('div', { class: 'rail-group', style: { margin: '14px 0 4px', paddingLeft: '10px' } }, 'Used lately'));
    for (const t of used) box.append(item(t, 'history'));
  }
  for (const [g, label, ic] of registry.TOOL_GROUPS) {
    let ts = registry.tools.filter((t) => t.group === g);
    if (q) ts = ts.filter((t) => match(q, `${t.title} ${t.summary || ''} ${(t.keywords || []).join(' ')}`));
    if (!ts.length) continue;
    box.append(el('div', { class: 'rail-group', style: { margin: '14px 0 4px', paddingLeft: '10px' } }, label));
    for (const t of ts) {
      shown++;
      box.append(item(t, ic));
    }
  }
  if (!shown) box.append(el('p', { class: 'note', style: { padding: '12px' } }, 'No tool by that name.'));
}

function open(id, values = null) {
  const t = registry.toolById(id);
  if (!t) return;
  view.current = id;
  for (const b of view.list.querySelectorAll('.tool-item')) b.classList.toggle('on', b.dataset.id === id);
  const pane = clear(view.pane);
  const group = registry.TOOL_GROUPS.find(([g]) => g === t.group);
  pane.append(el('div', { class: 'tool-head' }, el('div', { class: 'glyph' }, icon(t.icon || group?.[2] || 'toolbox', 'lg')),
    el('div', {}, el('div', { class: 'eyebrow', style: { marginBottom: '2px' } }, group?.[1] || ''), el('h2', {}, t.title), t.summary ? el('p', {}, t.summary) : null, t.refs ? el('div', { class: 'refs' }, t.refs) : null)));
  if (t.about) pane.append(el('p', { class: 'note', style: { maxWidth: '900px', margin: '-4px 0 16px' } }, t.about));
  // the form, kept per tool so that what was typed survives a visit elsewhere
  let f = view.forms.get(id);
  if (!f) { f = form(t.inputs || [], view.ctx, values || {}); view.forms.set(id, f); }
  else if (values) f.set(values);
  const results = el('div', { class: 'stage' });
  const runBtn = btn(t.action || 'Run', { icon: 'play', kind: 'primary' });
  const go = async () => {
    const problem = f.validate();
    if (problem) { toast(problem, 'error'); return; }
    runBtn.disabled = true;
    const keep = [...runBtn.childNodes];
    runBtn.replaceChildren(el('span', { class: 'spin' }), el('span', {}, 'Running…'));
    const values = f.values();
    const t0 = performance.now();
    remember(id);
    try {
      const result = await t.run(values, view.ctx);
      const run = { values, result, at: new Date(), seconds: (performance.now() - t0) / 1000 };
      const runs = view.runs.get(id) || [];
      runs.unshift(run);
      if (runs.length > 12) runs.length = 12;
      view.runs.set(id, runs);
      drawResults(t, results, run);
    } catch (e) {
      clear(results).append(callout(e.message, 'bad'));
    } finally {
      runBtn.disabled = false;
      runBtn.replaceChildren(...keep);
    }
  };
  runBtn.onclick = go;
  f.element.onkeydown = (e) => { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) go(); };
  const controls = el('div', { class: 'controls' }, el('div', { class: 'card' }, f.element, el('div', { class: 'run-bar' }, runBtn,
    t.example ? btn('Example', { icon: 'spark', onClick: async () => { f.set(typeof t.example === 'function' ? await t.example(view.ctx) : t.example); go(); } }) : null)),
    el('p', { class: 'small faint', style: { margin: '0 4px' } }, 'Ctrl Enter runs it.'));
  pane.append(el('div', { class: 'bench wide-controls' }, controls, results));
  const last = view.runs.get(id)?.[0];
  if (last) drawResults(t, results, last);
  else results.append(card({ body: empty({ icon: t.icon || group?.[2] || 'toolbox', title: 'Nothing run yet', text: t.hint || 'Fill in the form and run it; your files are in Your data, and anything you drop on the page joins them.', action: t.example ? btn('Run the example', { icon: 'spark', onClick: async () => { f.set(typeof t.example === 'function' ? await t.example(view.ctx) : t.example); go(); } }) : null }) }));
}

function drawResults(t, box, run) {
  clear(box);
  const runs = view.runs.get(t.id) || [];
  box.append(el('div', { class: 'results-head' }, el('span', { class: 'eyebrow', style: { margin: 0 } }, 'Result'), el('span', { class: 'when' }, `${run.at.toLocaleTimeString()} · ${fmt.seconds(run.seconds)}`),
    runs.length > 1 ? el('div', { class: 'history', style: { marginLeft: 'auto' } }, el('span', { class: 'small faint' }, 'earlier:'), ...runs.map((r, k) => el('button', { class: `chip${r === run ? ' on' : ''}`, onclick: () => drawResults(t, box, r) }, `${k + 1} · ${r.at.toLocaleTimeString()}`))) : null));
  const out = el('div', { class: 'stack', style: { gap: '14px' } });
  box.append(out);
  try { t.render(run.result, run.values, out, view.ctx); } catch (e) { out.append(callout(`The result could not be drawn: ${e.message}`, 'bad')); }
}
