// The editor of the network: settings, populations, projections, devices and recorders, each
// a card of fields written as the specification writes them ("20 ms", or a distribution as
// JSON). A field changes the specification at once; a change of structure (an item added or
// removed, a model changed) draws the editor again.

import { el, clear } from './dom.js';
import { colourOf } from './palette.js';
import { RULES, SYNAPSES, DEVICES, DEVICE_COMMON, RECORDERS, INTEGRATORS, ruleOf, setRule, addPopulation, renamePopulation, removePopulation, removeProjection, unique } from './spec.js';
import { format, hint, parseField, fieldText, wellFormed } from './units.js';

export function renderEditor(root, spec, catalogue, { onChange, onStructure, onFocusModel }) {
  clear(root);
  const changed = () => onChange?.();
  const restructure = () => { onStructure?.(); renderEditor(root, spec, catalogue, { onChange, onStructure, onFocusModel }); };
  const names = () => spec.populations.map((p) => p.name);
  const modelOf = (name) => catalogue.get(spec.populations.find((p) => p.name === name)?.model);

  // ---- fields ----
  const quantity = (label, get, set, unit, doc, placeholder) => {
    const input = el('input', { type: 'text', value: fieldText(get()), placeholder: placeholder ?? hint(unit), title: doc ?? null, spellcheck: 'false' });
    const mark = () => input.classList.toggle('bad', !wellFormed(input.value, unit));
    input.addEventListener('input', () => {
      mark();
      if (!wellFormed(input.value, unit)) return;
      try { set(parseField(input.value)); changed(); } catch { input.classList.add('bad'); }
    });
    mark();
    return field(label, input, doc);
  };
  const integer = (label, get, set, { min = 0, doc } = {}) => {
    const input = el('input', { type: 'number', min, step: 1, value: get() ?? '', title: doc ?? null });
    input.addEventListener('input', () => { const v = input.value === '' ? null : Math.round(Number(input.value)); set(v); changed(); });
    return field(label, input, doc);
  };
  const list = (label, get, set, kind) => {
    const input = el('input', { type: 'text', value: (get() ?? []).join(', '), placeholder: kind === 'index' ? '0, 1, 2' : `values in ${hint(kind)}, by commas`, spellcheck: 'false' });
    input.addEventListener('input', () => {
      const parts = input.value.split(',').map((s) => s.trim()).filter(Boolean);
      set(kind === 'index' ? parts.map(Number).filter(Number.isFinite) : parts.map((x) => parseField(x)));
      changed();
    });
    return field(label, input);
  };
  const choice = (label, options, get, set, { structural = false } = {}) => {
    const select = el('select', {}, options.map(([value, text]) => el('option', { value, text, selected: value === get() ? true : null })));
    select.addEventListener('change', () => { set(select.value); structural ? restructure() : changed(); });
    return field(label, select);
  };

  const remover = (what, act) => el('button', { class: 'remove', title: `Remove ${what}`, 'aria-label': `Remove ${what}`, onclick: () => { act(); restructure(); } }, '×');

  // ---- the run ----
  const run = section('Run');
  run.body.append(el('div', { class: 'grid' },
    quantity('duration', () => spec.duration, (v) => { spec.duration = v; }, 's', 'simulated time'),
    quantity('resolution', () => spec.resolution, (v) => { spec.resolution = v; }, 's', 'the step h of the grid of time'),
    integer('seed', () => spec.seed, (v) => { spec.seed = v; }, { doc: 'the seed of every random stream' }),
    choice('integrator', INTEGRATORS.map((m) => [m, m === 'default' ? "each model's own" : m]), () => spec.integrator ?? 'default', (v) => { if (v === 'default') delete spec.integrator; else spec.integrator = v; })));
  root.append(run.node);

  // ---- populations ----
  const pops = section('Populations', spec.populations.length);
  spec.populations.forEach((p, i) => {
    const model = catalogue.get(p.model);
    const card = el('div', { class: 'card' });
    const name = el('input', { type: 'text', value: p.name, class: 'title-input', 'aria-label': 'name', spellcheck: 'false' });
    name.addEventListener('change', () => {
      const to = name.value.trim() || p.name;
      if (to !== p.name) {
        renamePopulation(spec, p.name, unique(names().filter((n) => n !== p.name), to));
        restructure();
      }
    });
    card.append(el('div', { class: 'card-head' },
      el('i', { class: 'swatch', style: `background:${colourOf(i)}` }), name,
      remover(`population ${p.name}`, () => removePopulation(spec, p.name))));
    const body = el('div', { class: 'grid' },
      choice('model', catalogue.models.map((m) => [m.name, m.name]), () => p.model, (v) => { p.model = v; delete p.parameters; delete p.initial; delete p.bias; onFocusModel?.(v); }, { structural: true }),
      integer('neurons', () => p.size, (v) => { p.size = v; }, { min: 1 }));
    card.append(body);
    if (model) {
      const details = el('details', { class: 'more' }, el('summary', {}, 'parameters', p.parameters && Object.keys(p.parameters).length ? el('span', { class: 'count', text: ` ${Object.keys(p.parameters).length} set` }) : null));
      const grid = el('div', { class: 'grid' });
      for (const q of model.parameters)
        grid.append(quantity(q.name, () => p.parameters?.[q.name], (v) => {
          p.parameters ??= {};
          if (v === null) delete p.parameters[q.name]; else p.parameters[q.name] = v;
        }, q.unit, q.doc, format(q.default, q.unit, { ascii: true })));
      details.append(grid);
      const initial = el('details', { class: 'more' }, el('summary', {}, 'initial state and bias'));
      const ig = el('div', { class: 'grid' });
      for (const c of model.components)
        ig.append(quantity(c.name, () => p.initial?.[c.name], (v) => {
          p.initial ??= {};
          if (v === null) delete p.initial[c.name]; else p.initial[c.name] = v;
        }, c.unit, `${c.doc}: a value, or a distribution such as {"uniform": ["0 mV", "20 mV"]}`, 'at rest'));
      if (Array.isArray(p.bias)) {
        ig.append(field('bias', el('div', { class: 'fixed' }, `${p.bias.length} values, ${p.bias[0]} to ${p.bias[p.bias.length - 1]} `,
          el('button', { class: 'link', onclick: () => { delete p.bias; restructure(); } }, 'clear'))));
      } else {
        ig.append(quantity('bias', () => p.bias, (v) => { if (v === null) delete p.bias; else p.bias = v; }, model.input,
          'a constant input to each neuron: a value, a distribution, or a list with one value per neuron'));
      }
      initial.append(ig);
      card.append(details, initial);
    }
    pops.body.append(card);
  });
  pops.body.append(el('button', { class: 'add-item', onclick: () => { addPopulation(spec, 'lif_psc_exp'); restructure(); } }, '+ population'));
  root.append(pops.node);

  // ---- projections ----
  const projs = section('Projections', spec.projections.length);
  spec.projections.forEach((p, i) => {
    const card = el('div', { class: 'card' });
    const target = modelOf(p.target);
    const receptor = target?.receptors.find((r) => r.name === p.receptor);
    const r = ruleOf(p);
    card.append(el('div', { class: 'card-head' }, el('span', { class: 'title', text: `${i + 1}. ${p.source} → ${p.target}` }), remover(`projection ${i + 1}`, () => removeProjection(spec, i))));
    const grid = el('div', { class: 'grid' },
      choice('from', names().map((n) => [n, n]), () => p.source, (v) => { p.source = v; }, { structural: true }),
      choice('to', names().map((n) => [n, n]), () => p.target, (v) => {
        p.target = v;
        const m = modelOf(v);
        if (m && !m.receptors.some((x) => x.name === p.receptor)) p.receptor = m.receptors[0]?.name;
      }, { structural: true }),
      choice('receptor', (target?.receptors ?? []).map((x) => [x.name, x.name]), () => p.receptor, (v) => { p.receptor = v; }, { structural: true }),
      choice('rule', [...Object.entries(RULES).map(([k, v]) => [k, v.label]), ...(r.kind === 'pairs' ? [['pairs', 'listed pairs']] : [])], () => r.kind, (v) => setRule(p, v), { structural: true }));
    if (RULES[r.kind]?.value) grid.append(quantity(RULES[r.kind].value, () => r.value, (v) => setRule(p, r.kind, v), '1'));
    grid.append(
      quantity('weight', () => p.weight, (v) => { p.weight = v; }, receptor?.unit, receptor ? `jumps ${receptor.component}: ${receptor.doc}` : null),
      quantity('delay', () => p.delay, (v) => { p.delay = v; }, 's', 'at least one step'));
    const syn = p.synapse?.model ?? 'static';
    grid.append(choice('synapse', Object.keys(SYNAPSES).map((k) => [k, k]), () => syn, (v) => { p.synapse = v === 'static' ? undefined : { model: v }; if (!p.synapse) delete p.synapse; }, { structural: true }));
    card.append(grid);
    if (SYNAPSES[syn]?.length) {
      const sg = el('div', { class: 'grid' });
      for (const [k, u] of SYNAPSES[syn])
        sg.append(quantity(k, () => p.synapse?.[k], (v) => { p.synapse ??= { model: syn }; if (v === null) delete p.synapse[k]; else p.synapse[k] = v; }, u === 'w' ? receptor?.unit : u));
      card.append(el('details', { class: 'more', open: true }, el('summary', {}, `${syn} synapse`), sg));
    }
    projs.body.append(card);
  });
  if (spec.populations.length)
    projs.body.append(el('button', {
      class: 'add-item', onclick: () => {
        const s = spec.populations[0], t = spec.populations[spec.populations.length - 1];
        spec.projections.push({ source: s.name, target: t.name, receptor: catalogue.get(t.model)?.receptors[0]?.name, rule: { fixed_probability: 0.1 }, delay: '1.5 ms' });
        restructure();
      },
    }, '+ projection'));
  root.append(projs.node);

  // ---- devices ----
  const devs = section('Devices', spec.devices.length);
  spec.devices.forEach((d, i) => {
    const card = el('div', { class: 'card' });
    const target = modelOf(d.target);
    const receptor = target?.receptors.find((r) => r.name === d.receptor);
    card.append(el('div', { class: 'card-head' }, el('span', { class: 'title', text: `${DEVICES[d.model]?.label ?? d.model} → ${d.target}` }), remover('this device', () => spec.devices.splice(i, 1))));
    const grid = el('div', { class: 'grid' },
      choice('kind', Object.entries(DEVICES).map(([k, v]) => [k, v.label]), () => d.model, (v) => {
        for (const [k] of DEVICES[d.model]?.fields ?? []) delete d[k];
        d.model = v;
        if (DEVICES[v].fields.some(([, u]) => u === 'receptor')) d.receptor = modelOf(d.target)?.receptors[0]?.name;
      }, { structural: true }),
      choice('target', names().map((n) => [n, n]), () => d.target, (v) => { d.target = v; }, { structural: true }));
    for (const [k, u] of DEVICES[d.model]?.fields ?? []) {
      if (u === 'receptor') grid.append(choice('receptor', (target?.receptors ?? []).map((x) => [x.name, x.name]), () => d.receptor, (v) => { d.receptor = v; }, { structural: true }));
      else if (u.startsWith('list')) grid.append(list(k, () => d[k], (v) => { d[k] = v; }, u === 'list in' ? target?.input : u.slice(5)));
      else grid.append(quantity(k, () => d[k], (v) => { d[k] = v; }, u === 'w' ? receptor?.unit : u === 'in' ? target?.input : u));
    }
    card.append(grid);
    const more = el('div', { class: 'grid' });
    for (const [k, u] of DEVICE_COMMON)
      more.append(u === 'index' ? integer(k, () => d[k], (v) => { if (v === null) delete d[k]; else d[k] = v; }) : quantity(k, () => d[k], (v) => { if (v === null) delete d[k]; else d[k] = v; }, u));
    card.append(el('details', { class: 'more' }, el('summary', {}, 'when and to whom'), more));
    devs.body.append(card);
  });
  if (spec.populations.length)
    devs.body.append(el('button', { class: 'add-item', onclick: () => { spec.devices.push({ model: 'dc', target: spec.populations[0].name }); restructure(); } }, '+ device'));
  root.append(devs.node);

  // ---- recorders ----
  const recs = section('Recorders', spec.recorders.length);
  spec.recorders.forEach((r, i) => {
    const card = el('div', { class: 'card' });
    const what = r.model === 'weights' ? `projection ${Number(r.projection) + 1}` : r.population;
    card.append(el('div', { class: 'card-head' }, el('span', { class: 'title', text: `${RECORDERS[r.model]?.label ?? r.model} of ${what}` }), remover('this recorder', () => spec.recorders.splice(i, 1))));
    const grid = el('div', { class: 'grid' },
      choice('kind', Object.entries(RECORDERS).map(([k, v]) => [k, v.label]), () => r.model, (v) => {
        for (const [k] of RECORDERS[r.model]?.fields ?? []) delete r[k];
        r.model = v;
        if (v === 'weights') r.projection = 0; else r.population = spec.populations[0]?.name;
        if (v === 'state') r.neurons = [0];
      }, { structural: true }));
    const model = modelOf(r.population);
    for (const [k, u] of RECORDERS[r.model]?.fields ?? []) {
      if (u === 'population') grid.append(choice('population', names().map((n) => [n, n]), () => r.population, (v) => { r.population = v; delete r.components; }, { structural: true }));
      else if (u === 'projection') grid.append(choice('projection', spec.projections.map((p, j) => [String(j), `${j + 1}. ${p.source} → ${p.target}`]), () => String(r.projection), (v) => { r.projection = Number(v); }, { structural: true }));
      else if (u === 'index') grid.append(integer(k, () => r[k], (v) => { if (v === null) delete r[k]; else r[k] = v; }));
      else if (u === 'list index') grid.append(list(k, () => r[k], (v) => { r[k] = v; }, 'index'));
      else if (u === 'components') {
        const box = el('div', { class: 'checks' });
        for (const c of model?.components ?? []) {
          const on = !r.components || r.components.includes(c.name);
          const input = el('input', { type: 'checkbox', checked: on ? true : null });
          input.addEventListener('change', () => {
            const all = model.components.map((x) => x.name);
            const set = new Set(r.components ?? all);
            if (input.checked) set.add(c.name); else set.delete(c.name);
            r.components = all.filter((x) => set.has(x));
            changed();
          });
          box.append(el('label', { title: c.doc }, input, c.name));
        }
        grid.append(field('components', box));
      } else grid.append(quantity(k, () => r[k], (v) => { if (v === null) delete r[k]; else r[k] = v; }, u));
    }
    card.append(grid);
    recs.body.append(card);
  });
  if (spec.populations.length)
    recs.body.append(el('button', { class: 'add-item', onclick: () => { spec.recorders.push({ model: 'spikes', population: spec.populations[0].name }); restructure(); } }, '+ recorder'));
  root.append(recs.node);
}

function section(title, count) {
  const body = el('div', { class: 'section-body' });
  const node = el('section', { class: 'editor-section' }, el('h3', {}, title, count !== undefined ? el('span', { class: 'count', text: ` ${count}` }) : null), body);
  return { node, body };
}

function field(label, control, doc) {
  return el('label', { class: 'field', title: doc ?? null }, el('span', { class: 'label', text: label }), control);
}
