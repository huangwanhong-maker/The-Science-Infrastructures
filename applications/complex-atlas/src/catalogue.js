// The catalogue: the families of complexes the engine builds (spheres and surfaces, presentation
// complexes, lattices, projective spaces, triangulations, coloured graphs), grouped as it groups
// them; a family chosen shows its account and its parameters, and builds a complex.

import { el, clear } from './dom.js';

export function renderCatalogue(host, { families, selected, onSelect, onBuild, busy }) {
  clear(host);
  if (!families) {
    host.append(el('p', { class: 'note pad', text: 'The catalogue comes from the engine. Until it is built, draw the examples below, or your own complex.' }));
    return;
  }
  const groups = new Map();
  for (const f of families) (groups.get(f.group) || groups.set(f.group, []).get(f.group)).push(f);
  for (const [group, list] of groups) {
    host.append(el('div', { class: 'family', text: group }));
    for (const f of list) {
      host.append(el('button', { class: `item${selected?.id === f.id ? ' selected' : ''}`, onclick: () => onSelect(f) }, f.title));
      if (selected?.id === f.id) host.append(card(f, selected.params, onBuild, busy));
    }
  }
}

function card(f, values, onBuild, busy) {
  const inputs = {};
  const fields = (f.params || []).map((p) => {
    const value = values?.[p.name] ?? p.default;
    const input = p.type === 'text'
      ? el('input', { type: 'text', value: value ?? '', spellcheck: 'false' })
      : el('input', { type: 'number', value, min: p.min, max: p.max, step: p.type === 'int' ? 1 : 'any' });
    inputs[p.name] = [p, input];
    return el('label', { class: `field${p.type === 'text' ? ' wide' : ''}` }, el('span', { class: 'label', text: p.name }), input);
  });
  const params = () => Object.fromEntries(Object.entries(inputs).map(([k, [p, i]]) => [k, p.type === 'text' ? i.value : Number(i.value)]));
  const build = el('button', { class: 'button primary', disabled: busy ? true : null, onclick: () => onBuild(f, params()) }, 'Build');
  const form = el('div', { class: 'family-card' }, el('p', { text: f.description || '' }), fields.length ? el('div', { class: 'grid' }, fields) : null, el('div', { class: 'row' }, build));
  form.addEventListener('keydown', (e) => { if (e.key === 'Enter') onBuild(f, params()); });
  return form;
}
