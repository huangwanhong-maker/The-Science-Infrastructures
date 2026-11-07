// The toolbox (ADR 0056): the quantitative economist's instruments, each a form built from its
// declaration, run by the engine on a table (the library's, one's own, or one pasted here) and answered
// in blocks: key figures, tables, charts, matrices, notes, and values that can be carried to the map as
// a layer. The catalogue, the form and the blocks are the visualization layer's (ADR 0080): what stands
// here is the atlas's own rail, tables and map.

import * as api from './api.js';
import { el, clear, put, fill } from './dom.js';
import { GROUPS, TOOLS, tool as toolOf } from '/visualization/js/toolbox/catalogue.js';
import { field, freshForm, paramControl, parsePaste, roleControl as roleOf, rolesOf } from '/visualization/js/toolbox/form.js';
import { blockView } from '/visualization/js/toolbox/answer.js';
import { formatNumber } from '/visualization/js/toolbox/charts.js';
import { shared, emit, on } from './state.js';
import { openDataDialog } from './data.js';
import { onTheme } from './palette.js';

const state = {
  wired: false, tables: null, current: null, form: null, answer: null, filter: '',
};

// ---- the form of a tool ----

function columnsOf(form) {
  if (form.table === '__paste') return parsePaste(form.paste).columns;
  return state.tables?.find((t) => t.id === form.table)?.columns || [];
}

const roleControl = (role, form, columns) => roleOf(role, form, columns, () => renderForm());

function tableChooser(form) {
  const options = (state.tables || []).map((t) => el('option', { value: t.id, text: `${t.title}${t.own ? ' (yours)' : ''} · ${t.rows.toLocaleString('en-GB')} rows`, selected: t.id === form.table }));
  const select = el('select', { onchange: (e) => {
    if (e.target.value === '__add') { openDataDialog('table'); e.target.value = form.table; return; }
    form.table = e.target.value; form.where = {}; renderForm();
  } }, options, el('option', { value: '__paste', text: 'a table pasted here…', selected: form.table === '__paste' }), el('option', { value: '__add', text: '+ a table of your own…' }));
  const cols = columnsOf(form);
  const whereCol = Object.keys(form.where)[0] || '';
  const whereRow = el('div', { class: 'where' },
    el('select', { onchange: (e) => { const val = form.where[whereCol]; form.where = e.target.value ? { [e.target.value]: val ?? '' } : {}; renderForm(); } },
      el('option', { value: '', text: 'every row' }), cols.map((c) => el('option', { value: c.name, text: `where ${c.name} =`, selected: c.name === whereCol }))),
    whereCol ? el('input', { type: 'text', value: form.where[whereCol] ?? '', placeholder: 'value', oninput: (e) => {
      const kind = cols.find((c) => c.name === whereCol)?.kind;
      form.where = { [whereCol]: kind === 'number' && e.target.value !== '' && Number.isFinite(Number(e.target.value)) ? Number(e.target.value) : e.target.value };
    } }) : null);
  return el('div', { class: 'tb-card' }, el('div', { class: 'heading' }, el('span', { text: 'the data' })),
    field('Table', select),
    form.table === '__paste' ? field('Paste a table (a header row; comma, semicolon or tab separated)', el('textarea', { rows: 6, class: 'mono', oninput: (e) => { form.paste = e.target.value; },
      onchange: () => renderForm() }, form.paste)) : null,
    form.table !== '__paste' ? field('Rows', whereRow) : null,
    form.table !== '__paste' ? el('p', { class: 'note', text: state.tables?.find((t) => t.id === form.table)?.source || '' }) : null);
}

// ---- running and the answer -------------------------------------------------------------------------------

async function run() {
  const t = state.current, form = state.form;
  const status = document.getElementById('tb-status');
  const out = document.getElementById('tb-answer');
  const asked = rolesOf(t, form);
  if (asked.missing) { status.textContent = `choose ${asked.missing}`; return; }
  const body = { tool: t.id, roles: asked.roles, params: form.params };
  if (!t.noTable) {
    if (form.table === '__paste') body.columns = parsePaste(form.paste).data;
    else { body.table = form.table; if (Object.keys(form.where).length) body.where = form.where; }
  }
  status.textContent = 'computing…';
  status.classList.remove('error');
  try {
    state.answer = await api.tool(body);
    status.textContent = '';
    renderAnswer();
  } catch (error) {
    status.textContent = error.message;
    status.classList.add('error');
    clear(out);
  }
}

function placesOn(b) {
  if (b.kind !== 'map') return null;
  const entries = Object.entries(b.values || {}).filter(([, v]) => Number.isFinite(v)).sort((a, c) => c[1] - a[1]);
  return el('div', { class: 'tb-block' }, el('div', { class: 'heading' }, el('span', { text: b.title })),
    el('p', {}, `${entries.length} places. Highest: ${entries.slice(0, 5).map(([k, v]) => `${k} ${formatNumber(v)}`).join(', ')}.`),
    el('button', { class: 'primary', text: 'Show on the map', onclick: () => toMap(b) }));
}
const block = (b) => blockView(b, { name: state.current.id, places: placesOn });

// Values for places become a layer of the map, as a plug-in's do.
function toMap(b) {
  const key = `toolbox:${state.current.id}#${Date.now()}`;
  shared.layers = shared.layers.filter((l) => l.file !== `toolbox:${state.current.id}`);
  shared.layers.push({ key, file: `toolbox:${state.current.id}`, title: b.title, plugin: state.current.title, language: 'tool', year: state.form.params.year,
    dataset: 'library:trade', result: { kind: 'values', title: b.title, values: b.values, scale: b.scale || 'sequential', unit: b.unit || '' } });
  shared.colour = key;
  shared.flows = 'library:trade';
  emit('layers');
  location.hash = `#map${state.form.params.year ? `&year=${state.form.params.year}` : ''}`;
}

function renderAnswer() {
  const out = document.getElementById('tb-answer');
  const A = state.answer;
  if (!out || !A) return;
  fill(out, el('h3', { class: 'tb-answer-title', text: A.title || state.current.title }), (A.blocks || []).map(block));
}

// ---- the view ------------------------------------------------------------------------------------------------

function renderRail() {
  const rail = document.getElementById('tb-rail');
  const q = state.filter.toLowerCase();
  fill(rail,
    el('input', { type: 'search', class: 'tb-search', placeholder: 'find a tool…', value: state.filter, oninput: (e) => { state.filter = e.target.value; renderRail(); } }),
    GROUPS.map(([g, title]) => {
      const tools = TOOLS.filter((t) => t.group === g && (!q || `${t.title} ${t.summary}`.toLowerCase().includes(q)));
      if (!tools.length) return null;
      return el('div', { class: 'tb-group' }, el('div', { class: 'heading' }, el('span', { text: title })),
        tools.map((t) => el('button', { class: `tb-tool${state.current?.id === t.id ? ' on' : ''}`, onclick: () => choose(t.id) },
          el('span', { text: t.title }), el('span', { class: `badge ${t.foundation}`, text: t.foundation }))));
    }));
}

function renderForm() {
  const t = state.current, form = state.form;
  const main = document.getElementById('tb-main');
  const columns = columnsOf(form);
  fill(main,
    el('div', { class: 'tb-head' }, el('h2', {}, t.title, ' ', el('span', { class: `badge ${t.foundation}`, text: t.foundation })),
      el('p', { class: 'lede', text: t.summary }), t.about ? el('p', { class: 'note', text: t.about }) : null),
    el('div', { class: 'tb-form' },
      t.noTable ? null : tableChooser(form),
      (t.roles || []).length ? el('div', { class: 'tb-card' }, el('div', { class: 'heading' }, el('span', { text: 'the variables' })),
        t.roles.map((r) => roleControl(r, form, columns))) : null,
      (t.params || []).length ? el('div', { class: 'tb-card' }, el('div', { class: 'heading' }, el('span', { text: t.noTable ? 'the model' : 'the choices' })),
        el('div', { class: 'params' }, t.params.map((q) => paramControl(q, form)))) : null),
    el('div', { class: 'tb-run' }, el('button', { class: 'primary', text: 'Run', onclick: run }),
      el('button', { class: 'quiet', text: 'the example again', onclick: () => { state.form = freshForm(t); renderForm(); } }),
      el('span', { class: 'note', id: 'tb-status', role: 'status' })),
    el('div', { id: 'tb-answer', class: 'tb-answer' }));
  renderAnswer();
}

function choose(id) {
  const t = toolOf(id) || TOOLS[0];
  state.current = t;
  state.form = freshForm(t);
  state.answer = null;
  if (location.hash.slice(1).split('&')[0] === 'toolbox') history.replaceState(null, '', `#toolbox&tool=${t.id}`);
  renderRail();
  renderForm();
  // A tool with an example runs it at once.
  if (t.preset || t.noTable) run();
}

// The address may name a tool: #toolbox&tool=ols.
export async function showToolbox(data, params = {}) {
  if (!state.wired) {
    state.wired = true;
    onTheme(() => renderAnswer());
    // A table of one's own just added: the tables again, and the form on it.
    on('dataset-chosen', async (kind, id) => {
      if (kind !== 'table') return;
      state.tables = (await api.tables()).tables;
      if (state.current && !state.current.noTable && id) { state.form.table = id; state.form.where = {}; state.form.roles = {}; }
      renderForm();
    });
  }
  try {
    state.tables = (await api.tables()).tables;
  } catch (error) {
    state.tables = [];
    document.getElementById('status').textContent = error.message;
  }
  if (params.tool && toolOf(params.tool)) choose(params.tool);
  else if (!state.current) choose('summary');
  else { renderRail(); renderForm(); }
}
