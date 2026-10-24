// Plug-ins (ADR 0046): Python scripts and C++ programs in the atlas's plugins/ directory, run by the
// server through the API's science.networks on the network the map shows (its places, its flows in
// the year, and the engine's decomposition of them). Each describes itself and its parameters, from
// which its form is built; what it returns becomes layers: values colour the places, arcs are drawn
// between them, sketches over them, and tables and notes go to the panel.

import * as api from './api.js';
import { el, clear, put, fill } from './dom.js';
import { shared, emit } from './state.js';

const state = { list: null, described: new Map(), running: new Map(), runs: 0 };

function paramField(p, values) {
  const id = `param-${p.name}`;
  const set = (v) => { values[p.name] = v; };
  if (!(p.name in values)) values[p.name] = p.default;
  let input;
  if (p.kind === 'choice') {
    input = el('select', { id, onchange: (e) => set(e.target.value) },
      (p.options || []).map((o) => (o && typeof o === 'object' ? [o.value, o.label ?? String(o.value)] : [o, String(o)]))
        .map(([v, label]) => el('option', { value: v, text: label, selected: v === values[p.name] })));
  } else if (p.kind === 'boolean') {
    input = el('input', { id, type: 'checkbox', checked: Boolean(values[p.name]), onchange: (e) => set(e.target.checked) });
    return el('label', { class: 'toggles', title: p.help || '' }, el('label', {}, input, p.label || p.name));
  } else {
    input = el('input', { id, type: 'number', value: values[p.name], step: p.step ?? (p.kind === 'integer' ? 1 : 'any'),
      min: p.min ?? undefined, max: p.max ?? undefined,
      onchange: (e) => set(p.kind === 'integer' ? parseInt(e.target.value, 10) : Number(e.target.value)) });
  }
  return el('label', { class: 'field', title: p.help || '' }, el('span', { text: p.label || p.name }), input);
}

function card(file, language, context) {
  const values = {};
  const body = el('div', { class: 'body' });
  const titleEl = el('span', { class: 'title', text: state.described.get(file)?.title || file.replace(/^plugins\//, '') });
  const details = el('details', { class: 'plugin' },
    el('summary', {}, el('span', { class: `badge ${language === 'cpp' ? 'cpp' : 'py'}`, text: language === 'cpp' ? 'C++' : 'Py' }), titleEl),
    body);
  const open = async () => {
    fill(body, el('p', { class: 'note', text: language === 'cpp' ? 'compiling and describing…' : 'describing…' }));
    let d = state.described.get(file);
    try {
      if (!d) {
        d = await api.describePlugin(file);
        state.described.set(file, d);
        titleEl.textContent = d.title || titleEl.textContent;
      }
    } catch (error) {
      fill(body, el('p', { class: 'error', text: error.message }));
      return;
    }
    const status = el('div', { class: 'note' });
    const bar = el('div', { class: 'progress', hidden: true }, el('span', { style: { width: '0%' } }));
    const run = el('button', { class: 'primary', text: 'Run on this map' });
    const stop = el('button', { class: 'quiet', text: 'stop', hidden: true });
    run.onclick = async () => {
      const { dataset, year } = context();
      run.disabled = true;
      bar.hidden = false;
      status.classList.remove('error');
      status.textContent = 'starting…';
      try {
        const { job } = await api.runPlugin(file, values, dataset, year);
        stop.hidden = false;
        stop.onclick = () => api.cancelJob(job).catch(() => {});
        await follow(job, d, bar, status);
      } catch (error) {
        status.textContent = error.message;
        status.classList.add('error');
      }
      run.disabled = false;
      stop.hidden = true;
    };
    fill(body, el('p', { text: d.description || '' }),
      (d.params || []).map((p) => paramField(p, values)),
      el('div', { style: { display: 'flex', gap: '8px', alignItems: 'center' } }, run, stop), bar, status);
  };
  details.addEventListener('toggle', () => { if (details.open && !body.firstChild) open(); });
  return details;
}

async function follow(id, description, bar, status) {
  for (;;) {
    const j = await api.job(id);
    bar.firstChild.style.width = `${Math.round(100 * (j.progress || 0))}%`;
    status.textContent = j.text || '';
    if (j.state === 'done') {
      adopt(j, description);
      status.textContent = `done in ${j.seconds} s: ${(j.results || []).length} result${(j.results || []).length === 1 ? '' : 's'}`;
      return;
    }
    if (j.state === 'failed' || j.state === 'cancelled') {
      status.textContent = j.state === 'cancelled' ? 'stopped' : j.message || 'failed';
      status.classList.toggle('error', j.state === 'failed');
      return;
    }
    await new Promise((r) => setTimeout(r, 600));
  }
}

// A run's results become layers; the first values colour the places and the first arcs are drawn.
function adopt(j, d) {
  const run = ++state.runs;
  shared.layers = shared.layers.filter((l) => l.file !== j.file);   // a plug-in's newer run replaces its older
  let colour = null, arcs = null;
  (j.results || []).forEach((result, i) => {
    const key = `${j.file}#${run}.${i}`;
    shared.layers.push({ key, file: j.file, title: result.title || d.title, plugin: d.title, language: d.language, year: j.year, dataset: shared.flows, result });
    if (result.kind === 'values' && !colour) colour = key;
    if (result.kind === 'arcs' && !arcs) arcs = key;
  });
  if (colour) shared.colour = colour;
  if (arcs) shared.arcs = arcs;
  emit('layers');
}

export async function pluginsSection(context) {
  const box = el('div', {});
  const list = el('div', {});
  const reload = async () => {
    fill(list, el('p', { class: 'note', text: 'looking in plugins/…' }));
    try {
      state.list = (await api.plugins()).plugins;
      for (const p of state.list) if (p.description) state.described.set(p.file, p.description);
      fill(list, state.list.length ? state.list.map((p) => card(p.file, p.language, context))
        : el('p', { class: 'note', text: 'No plug-ins yet: put a Python script or a C++ program in applications/economy-atlas/plugins/ (its README shows how).' }));
    } catch (error) {
      fill(list, el('p', { class: 'error', text: error.message }));
    }
  };
  put(box, el('div', { class: 'heading' }, el('span', { text: 'plug-ins' }), el('button', { class: 'note', text: 'look again', onclick: reload })),
    el('p', { class: 'note', text: 'Python and C++, run on the flows shown; their results become layers.' }), list);
  await reload();
  return box;
}
