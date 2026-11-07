// The toolbox of quantitative economics on the place open (ADR 0080). Questions in plain words first, each a tool's form filled in; or
// any tool, on any table here. What a tool answers for places goes on the map through the line of time, and stands beside its
// numbers. One thing at a time: choosing a question or a tool folds the list away.

import { el, fill } from '/visualization/js/toolbox/dom.js';
import { GROUPS, TOOLS, tool as toolOf } from '/visualization/js/toolbox/catalogue.js';
import { field, freshForm, paramControl, roleControl, rolesOf } from '/visualization/js/toolbox/form.js';
import { blockView } from '/visualization/js/toolbox/answer.js';
import { formatNumber, lineChart } from '/visualization/js/toolbox/charts.js';
import { tabs } from '/src/quiet.js';
import { findTable, pasted, tablesHere, worldMosaic } from '../tables.js';
import { questionsHere } from '../questions.js';
import { drawables, ranking, recordOf, styleOf, whereIs } from '../answers.js';

const DEG = Math.PI / 180;
const memory = { tool: null, form: null, question: null, table: null, answer: null, drawn: null, chosen: null };
const REGRESSIONS = new Set(['ols', 'glm', 'panel']);

export async function render(root, ctx) {
  const open = ctx.spec?.open || 'questions';
  root.append(el('div', { class: 'view-head' }, el('h1', { text: 'Quantitative economics' }),
    el('span', { class: 'lede', text: 'Regressions, series, inequality and games, on this place’s data. What is found of places stands on the map.' })));
  const wrap = el('div', { class: 'tb' });
  root.append(wrap);

  // ---- running ----------------------------------------------------------------------------------------------------------------------------------------

  const status = el('p', { class: 'tb-status', role: 'status' });
  const answerEl = el('div', { class: 'tb-answer' });

  // Only the columns a tool's roles name, and those it is split by, go to the engine, in the rows the filter keeps.
  function requestOf(t, form, table, each, byRow) {
    const asked = rolesOf(t, form);
    if (asked.missing) throw new Error(`Choose ${asked.missing}.`);
    const body = { tool: t.id, roles: asked.roles, params: form.params };
    if (t.noTable) return { body, kept: null };
    const wanted = new Set([...Object.values(asked.roles).flat(), ...each]);
    for (const c of wanted) if (!table.data[c]) throw new Error(`The table has no column ${c}.`);
    const conditions = Object.entries(form.where || {}).filter(([, v]) => v !== '' && v !== null && v !== undefined);
    const n = table.rows;
    const kept = [];
    for (let i = 0; i < n; i++) {
      if (conditions.every(([c, v]) => table.data[c] && (table.data[c][i] === v || String(table.data[c][i]) === String(v)))) kept.push(i);
    }
    if (!kept.length) throw new Error('No row of the table is kept by that filter.');
    body.columns = Object.fromEntries([...wanted].map((c) => [c, kept.map((i) => table.data[c][i])]));
    if (each.length) body.each = each;
    if (byRow) {
      body.by_row = true;
      if (table.data.place_lon) body.places = { lon: kept.map((i) => table.data.place_lon[i]), lat: kept.map((i) => table.data.place_lat[i]) };
      else if (table.places?.kind === 'row') body.places = { lon: kept.map((i) => table.data[table.places.lon][i]), lat: kept.map((i) => table.data[table.places.lat][i]) };
    }
    return { body, kept };
  }

  async function run(t, form, table, { each = [], byRow = false, draw = null, read = null, title = null } = {}) {
    status.textContent = 'computing…';
    status.classList.remove('error');
    answerEl.replaceChildren();
    try {
      const { body, kept } = requestOf(t, form, table, each, byRow);
      const A = await ctx.ask('tool', body);
      // The engine's rows are of the rows kept: back to the table's.
      for (const b of A.blocks || []) if (b.kind === 'rows' && kept) b.rows = b.rows.map((i) => kept[i]);
      memory.answer = { A, t, table, draw, read, title, source: table?.source || '' };
      status.textContent = '';
      await showAnswer();
    } catch (error) {
      status.textContent = String(error.message || error);
      status.classList.add('error');
    }
  }

  // ---- the answer ----------------------------------------------------------------------------------------------------------------------------------------

  async function drawOnMap(d, table, source) {
    const record = await recordOf(d, table, source, worldMosaic);
    memory.drawn?.remove();
    memory.drawn = ctx.st.show(record, styleOf(d, record, table));
    memory.chosen = d.title;
    // Go there if it is not in view.
    const [s, w, n, e] = memory.drawn.record.box, v = ctx.earth.view;
    const across = Math.hypot(n - s, e - w) * 111e3, away = Math.hypot((s + n) / 2 - v.lat / DEG, ((((w + e) / 2 - v.lon / DEG) + 540) % 360) - 180) * 111e3;
    if (away > Math.max(v.range * 1.4, across) || across > v.range * 5) ctx.st.fit(memory.drawn.id);
    else ctx.st.sceneChanged?.();
  }

  function mapStrip(list, table, source) {
    const strip = el('div', { class: 'tb-map' });
    const paint = (d) => {
      const R = ranking(d, table);
      const row = (x) => el('tr', { class: whereIs(x.id, table) ? 'go' : '', onclick: () => {
        const c = whereIs(x.id, table);
        if (c) ctx.earth.flyTo({ lon: c[0] * DEG, lat: c[1] * DEG, height: ctx.earth.tree.heightAt(c[0] * DEG, c[1] * DEG) || 0, range: Math.min(ctx.earth.view.range, 2.6e3 * 6) }, 1.4);
      } }, el('td', { text: x.name }), el('td', { class: 'num', text: formatNumber(x.value) }));
      const words = el('div', { class: 'tb-rank' },
        el('table', { class: 'results' }, el('thead', {}, el('tr', {}, el('th', { text: 'Highest' }), el('th', {}))), el('tbody', {}, R.high.map(row))),
        el('table', { class: 'results' }, el('thead', {}, el('tr', {}, el('th', { text: 'Lowest' }), el('th', {}))), el('tbody', {}, R.low.map(row))));
      return el('div', {}, R.last !== undefined ? el('p', { class: 'note', text: `The list is of the last moment, ${R.last}; ${R.count.toLocaleString('en')} places. The line of time above the clock moves the map.` }) : el('p', { class: 'note', text: `${R.count.toLocaleString('en')} places.` }), words);
    };
    const body = el('div', {});
    const choose = async (d) => {
      memory.chosen = d.title;
      body.replaceChildren(el('p', { class: 'note', text: 'drawing…' }));
      try {
        await drawOnMap(d, table, source);
        body.replaceChildren(paint(d));
      } catch (error) {
        body.replaceChildren(el('p', { class: 'note error', text: String(error.message || error) }));
      }
    };
    const start = list.find((d) => memory.answer?.title && d.title === memory.answer.title) || list.find((d) => memory.answer?.draw && d.title.includes(memory.answer.draw)) || list[0];
    const picker = list.length > 6
      ? el('select', { class: 'tb-pick', onchange: (e) => choose(list[+e.target.value]) }, list.map((d, k) => el('option', { value: k, selected: d === start, text: d.title })))
      : list.length > 1 ? el('div', { class: 'tb-keys' }, list.map((d) => el('button', { class: d === start ? 'on' : '', title: d.title, onclick: (e) => { for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget); choose(d); } }, d.title.split(':')[0].slice(0, 42)))) : null;
    strip.append(el('div', { class: 'heading' }, el('span', { text: 'On the map' })), list.length === 1 ? el('p', { class: 'tb-drawn', text: start.title }) : picker, body);
    choose(start);
    return strip;
  }

  // Each part's numbers through time only: a line for a chosen number, there being no place to draw.
  function eachOverTime(b) {
    const body = el('div', {});
    const xs = b.keys.map((k) => Number(k[0]));
    const order = xs.map((_, i) => i).sort((p, q) => xs[p] - xs[q]);
    const paint = (m) => body.replaceChildren(lineChart({ title: m.name, x: order.map((i) => xs[i]), series: [{ name: m.name, y: order.map((i) => m.values[i]) }], xlabel: b.columns[0], ylabel: '' }));
    const start = b.measures.find((m) => memory.answer?.draw && m.name.includes(memory.answer.draw)) || b.measures[0];
    paint(start);
    return el('div', { class: 'tb-block' }, b.measures.length > 1 ? el('select', { class: 'tb-pick', onchange: (e) => paint(b.measures[+e.target.value]) }, b.measures.map((m, k) => el('option', { value: k, selected: m === start, text: m.name }))) : null, body);
  }

  function leftOut(b) {
    if (!b.left?.length) return null;
    const count = b.left.reduce((a, x) => a + x.count, 0);
    return el('p', { class: 'note', text: `${count.toLocaleString('en')} left out: ${b.left.slice(0, 3).map((x) => `${x.count.toLocaleString('en')}, ${x.reason}`).join('; ')}.` });
  }

  async function showAnswer() {
    const held = memory.answer;
    if (!held) return;
    const { A, t, table, source, read } = held;
    const blocks = A.blocks || [];
    const list = drawables(blocks, table);
    // A question names what it draws first in plain words.
    const named = held.title && list.find((d) => held.draw && d.title.includes(held.draw));
    if (named) named.title = held.title;
    if (!list.length) { memory.drawn?.remove(); memory.drawn = null; }
    if (blocks.some((b) => b.kind === 'map') && !table) {
      // A tool with no table: its places are economies. Their names come with the mosaic's table.
      const economies = await findTable(ctx, 'economies')?.load().catch(() => null);
      held.table = economies || null;
    }
    const each = blocks.find((b) => b.kind === 'each');
    fill(answerEl,
      el('h3', { class: 'tb-answer-title', text: A.title || t.title }),
      read ? el('p', { class: 'tb-read', text: read }) : null,
      list.length ? mapStrip(list, held.table, source) : null,
      each ? leftOut(each) : null,
      each && !list.length ? eachOverTime(each) : null,
      blocks.map((b) => blockView(b, { name: t.id })));
  }

  // ---- questions -------------------------------------------------------------------------------------------------------------------------------------------

  const questionsPanel = () => {
    const box = el('div', {});
    const list = () => {
      const groups = questionsHere(ctx.area?.id);
      fill(box, groups.map((g) => el('div', { class: 'tb-group' }, el('div', { class: 'heading' }, el('span', { text: g.title })),
        g.items.map((q) => el('button', { class: 'tb-q', onclick: () => ask(q) }, el('span', { class: 'tb-ask', text: q.ask }), el('span', { class: 'tb-by', text: toolOf(q.tool)?.title || q.tool }))))));
    };
    const ask = async (q) => {
      memory.question = q;
      const t = toolOf(q.tool);
      fill(box, el('button', { class: 'tb-back', onclick: list }, '← all questions'), el('p', { class: 'tb-asked', text: q.ask }), status, answerEl);
      let table = null;
      if (q.table) {
        const entry = findTable(ctx, q.table);
        if (!entry) { status.textContent = 'That table is not of this place.'; return; }
        status.textContent = 'reading the table…';
        try { table = await entry.load(); } catch (error) { status.textContent = String(error.message || error); return; }
      }
      const form = { table: q.table || null, paste: '', where: { ...(q.filter || {}) }, roles: JSON.parse(JSON.stringify(q.roles || {})), params: { ...Object.fromEntries((t.params || []).map((p) => [p.name, p.default])), ...q.params } };
      memory.tool = t;
      memory.form = form;
      memory.each = q.each || [];
      memory.table = table;
      await run(t, form, table, { each: q.each || [], byRow: !!q.byRow, draw: q.draw, read: q.read, title: q.title || null });
    };
    list();
    return box;
  };

  // ---- tools -----------------------------------------------------------------------------------------------------------------------------------------------

  const toolsPanel = () => {
    const box = el('div', {});
    const list = () => {
      fill(box, el('div', { class: 'tb-tools' }, GROUPS.map(([g, title]) => {
        const tools = TOOLS.filter((t) => t.group === g);
        return el('details', { class: 'tb-fold', open: g === open || open === 'questions' ? true : null }, el('summary', { text: title }),
          tools.map((t) => el('button', { class: 'tb-tool', onclick: () => choose(t) }, el('span', { class: 'tb-ask', text: t.title }), el('span', { class: 'tb-by', text: t.summary }))));
      })));
    };

    const choose = async (t) => {
      memory.tool = t;
      memory.form = freshForm(t);
      memory.each = [];
      memory.table = null;
      memory.answer = null;
      answerEl.replaceChildren();
      const here = tablesHere(ctx);
      // A preset's table if this place has it, else what stands first.
      if (!t.noTable && memory.form.table !== '__paste' && !here.some((e) => e.id === memory.form.table)) memory.form.table = here[0]?.id;
      await loadTable();
      drawForm();
      if (t.preset || t.noTable) go();
    };

    const loadTable = async () => {
      const form = memory.form;
      if (memory.tool.noTable) { memory.table = null; return; }
      if (form.table === '__paste') { memory.table = pasted(form.paste || ''); return; }
      const entry = findTable(ctx, form.table);
      memory.table = entry ? await entry.load() : null;
    };

    const go = () => {
      const t = memory.tool;
      if (memory.form.table === '__paste') memory.table = pasted(memory.form.paste || '');
      const each = memory.each || [];
      run(t, memory.form, memory.table, { each, byRow: REGRESSIONS.has(t.id) && !each.length && !!memory.table?.places });
    };

    const dataCard = () => {
      const form = memory.form, T = memory.table;
      const here = tablesHere(ctx);
      const groups = [...new Set(here.map((e) => e.group))];
      const select = el('select', { onchange: async (e) => {
        form.table = e.target.value;
        form.where = {};
        form.roles = {};
        memory.each = [];
        status.textContent = 'reading the table…';
        try { await loadTable(); status.textContent = ''; } catch (error) { status.textContent = String(error.message || error); memory.table = null; }
        if (memory.tool.roles?.some((r) => r.role === 'lon') && memory.table?.data.place_lon) Object.assign(form.roles, { lon: 'place_lon', lat: 'place_lat' });
        drawForm();
      } }, groups.map((g) => el('optgroup', { label: g }, here.filter((e) => e.group === g).map((e) => el('option', { value: e.id, selected: e.id === form.table, text: e.title })))),
      el('option', { value: '__paste', selected: form.table === '__paste', text: 'A table pasted here…' }));
      const cols = T?.columns || [];
      const whereCol = Object.keys(form.where)[0] || '';
      const whereRow = el('div', { class: 'where' },
        el('select', { onchange: (e) => { const val = form.where[whereCol]; form.where = e.target.value ? { [e.target.value]: val ?? '' } : {}; drawForm(); } },
          el('option', { value: '', text: 'every row' }), cols.map((c) => el('option', { value: c.name, selected: c.name === whereCol, text: `where ${c.name} =` }))),
        whereCol ? el('input', { type: 'text', value: form.where[whereCol] ?? '', placeholder: 'value', oninput: (e) => { form.where = { [whereCol]: e.target.value }; } }) : null);
      // For each: only for what a table's rows can be split by.
      const P = T?.places;
      const splits = [];
      if (P?.column && !memory.tool.noTable) {
        splits.push(['', 'the whole table at once']);
        splits.push([P.column, 'each place']);
        if (P.time) { splits.push([`${P.column},${P.time}`, 'each place, at each moment']); splits.push([P.time, 'each moment']); }
      }
      const eachRow = splits.length && memory.tool.id !== 'moran' ? field('Run it for', el('select', { onchange: (e) => { memory.each = e.target.value ? e.target.value.split(',') : []; } },
        splits.map(([v, label]) => el('option', { value: v, selected: (memory.each || []).join(',') === v, text: label })))) : null;
      return el('div', { class: 'tb-card' }, el('div', { class: 'heading' }, el('span', { text: 'The data' })),
        field('Table', select),
        form.table === '__paste' ? field('Paste a table (a header row; comma, semicolon or tab separated)', el('textarea', { rows: 6, class: 'mono', oninput: (e) => { form.paste = e.target.value; }, onchange: () => { memory.table = pasted(form.paste); drawForm(); } }, form.paste)) : null,
        form.table !== '__paste' ? field('Rows', whereRow) : null,
        eachRow,
        T?.source ? el('p', { class: 'note', text: T.source }) : null);
    };

    const drawForm = () => {
      const t = memory.tool, form = memory.form;
      const columns = memory.table?.columns || [];
      fill(box, el('button', { class: 'tb-back', onclick: list }, '← all tools'),
        el('div', { class: 'tb-head' }, el('h2', { text: t.title }), el('p', { class: 'lede', text: t.summary }), t.about ? el('p', { class: 'note', text: t.about }) : null),
        t.noTable ? null : dataCard(),
        (t.roles || []).length ? el('div', { class: 'tb-card' }, el('div', { class: 'heading' }, el('span', { text: 'The variables' })), t.roles.map((r) => roleControl(r, form, columns, drawForm))) : null,
        (t.params || []).length ? el('div', { class: 'tb-card' }, el('div', { class: 'heading' }, el('span', { text: t.noTable ? 'The model' : 'The choices' })), el('div', { class: 'params' }, t.params.map((q) => paramControl(q, form)))) : null,
        el('div', { class: 'tb-run' }, el('button', { class: 'key on', onclick: go }, 'Run'), el('button', { class: 'key', onclick: () => choose(t) }, 'the example again')),
        status, answerEl);
    };

    if (memory.form && memory.tool && open !== 'questions' && memory.tool.group === open) { drawForm(); if (memory.answer) showAnswer(); } else list();
    return box;
  };

  const panels = tabs(ctx.dom.h, [{ name: 'Questions', panel: questionsPanel }, { name: 'Tools', panel: toolsPanel }], { lazy: true });
  wrap.append(panels.el);
  panels.show(open === 'questions' ? 'Questions' : 'Tools');
}
