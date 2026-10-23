// The Fields panel: a lattice gauge field on the complex. A group element sits on each edge (a link
// variable); the holonomy of a face is the product along its attaching word; the Wilson action
// weighs each face by β (1 − Re tr U / N). The engine samples the field by Monte Carlo and returns
// its last configuration, the history of the mean plaquette, the exact value where one is known
// (two-dimensional Yang-Mills on a closed surface), Wilson loops along the generators of the
// fundamental group, and, for U(1) on a closed orientable surface, the topological charge.

import { el, clear } from './dom.js';
import { lineChart } from './chart.js';
import { number, withError, wordText } from './format.js';

// The form's values, kept between complexes.
export const form = { kind: 'U1', n: 3, beta: 2, sweeps: 400, thermalise: 100, update: 'heatbath', seed: 1 };

export function groupOf(f) {
  if (f.kind === 'Z') return { kind: 'Z', n: Math.max(2, f.n | 0) };
  if (f.kind === 'SU') return { kind: 'SU', n: Math.max(2, f.n | 0) };
  return { kind: f.kind };
}

export function renderFields(host, { complex, analysis, names, result, running, error, showing, onRun, onSelect, onShow, selection }) {
  clear(host);
  if (!complex) return;
  if (!(complex.faces || []).length) {
    host.append(el('p', { class: 'note', text: 'A gauge field needs 2-cells: its curvature is read on them. This complex has none.' }));
    return;
  }
  host.append(el('p', { class: 'note', text: 'A group element on each edge; the holonomy of each face, the product along its attaching word, enters the Wilson action β Σ (1 − Re tr U / N). Gauge transformations act at the vertices.' }));
  const field = (label, input, wide) => el('label', { class: `field${wide ? ' wide' : ''}` }, el('span', { class: 'label', text: label }), input);
  const kind = el('select', { onchange: () => { form.kind = kind.value; n.disabled = !['Z', 'SU'].includes(form.kind); } },
    ...[['Z', 'ℤ/N'], ['U1', 'U(1)'], ['SU2', 'SU(2)'], ['SU', 'SU(N)']].map(([v, t]) => el('option', { value: v, selected: form.kind === v ? true : null }, t)));
  const n = el('input', { type: 'number', min: 2, max: 8, value: form.n, disabled: !['Z', 'SU'].includes(form.kind) ? true : null, oninput: () => { form.n = Number(n.value); } });
  const num = (key, attrs) => { const i = el('input', { type: 'number', value: form[key], ...attrs, oninput: () => { form[key] = Number(i.value); } }); return i; };
  const update = el('select', { onchange: () => { form.update = update.value; } },
    ...[['heatbath', 'heat bath'], ['metropolis', 'Metropolis']].map(([v, t]) => el('option', { value: v, selected: form.update === v ? true : null }, t)));
  const run = el('button', { class: 'button primary', disabled: running ? true : null, onclick: () => onRun({ group: groupOf(form), beta: form.beta, sweeps: form.sweeps, thermalise: form.thermalise, update: form.update, seed: form.seed }) }, running ? 'Running…' : 'Run');
  host.append(el('div', { class: 'grid' },
    field('group', kind), field('N', n),
    field('β, the coupling', num('beta', { step: 0.1, min: 0 })), field('update', update),
    field('sweeps measured', num('sweeps', { min: 1, step: 50 })), field('sweeps to thermalise', num('thermalise', { min: 0, step: 50 })),
    field('seed', num('seed', { min: 0 }))),
  el('div', { class: 'row' }, run, error ? el('span', { class: 'error', text: error }) : null));

  if (!result) return;
  const r = result;
  host.append(el('label', { class: 'check' }, el('input', { type: 'checkbox', checked: showing ? true : null, onchange: (e) => onShow(e.target.checked) }), 'colour the complex by this field'));
  const tile = (label, value) => el('div', { class: 'tile' }, el('div', { class: 'tile-label', text: label }), el('div', { class: 'tile-value', text: value }));
  host.append(el('div', { class: 'tiles' },
    tile(`⟨Re tr U / N⟩, ${r.group}`, withError(r.mean?.plaquette, r.mean?.error)),
    tile('exact', r.exact ? number(r.exact.plaquette, 5) : '—'),
    tile('autocorrelation time', r.mean?.tau !== undefined ? `${number(r.mean.tau, 3)} sweeps` : '—'),
    tile(r.charge !== null && r.charge !== undefined ? 'topological charge' : 'acceptance', r.charge !== null && r.charge !== undefined ? number(r.charge) : number(r.acceptance, 3))));
  if (r.history?.plaquette?.length > 1) {
    host.append(lineChart(r.history.plaquette, { exact: r.exact?.plaquette ?? null, mean: r.mean?.plaquette ?? null, label: 'sweeps' }));
    host.append(el('p', { class: 'note', text: `The mean plaquette after each sweep${r.exact ? '; dashed, the exact value' : ''}; the lilac line, the mean.` }));
  }
  if (r.loops?.length) {
    host.append(el('h4', { text: 'Wilson loops' }));
    const gens = analysis?.fundamental_group?.simplified?.loops ? analysis.fundamental_group.simplified : analysis?.fundamental_group;
    const table = el('table', { class: 'data' }, el('tr', {}, el('th', { text: 'loop' }), el('th', { text: '⟨Re tr W / N⟩' })));
    r.loops.forEach((l, i) => {
      const name = gens?.loops?.findIndex((w) => JSON.stringify(w) === JSON.stringify(l.loop));
      const label = name >= 0 ? gens.generators[name] : wordText(l.loop, names, { powers: false });
      const key = ['wilson', i];
      const on = selection && JSON.stringify(selection.key) === JSON.stringify(key);
      table.append(el('tr', {}, el('td', {}, el('button', { class: `chip${on ? ' on' : ''}`, onclick: () => onSelect(on ? null : { kind: 'loop', word: l.loop, key }) }, label)), el('td', { text: withError(l.mean, l.error) })));
    });
    host.append(table);
  }
}
