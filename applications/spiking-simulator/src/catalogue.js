// The catalogue of models, as the engine describes them: what each is, its source, its
// parameters with their defaults and units, its state and its receptors.

import { el, clear } from './dom.js';
import { format, hint } from './units.js';

// The models by family, in the order the catalogue shows them.
export const FAMILIES = [
  ['Integrate-and-fire', ['lif_psc_delta', 'lif_psc_exp', 'lif_psc_alpha', 'lif_cond_exp', 'alif', 'adlif', 'glif', 'mat2']],
  ['Escape noise', ['lif_escape', 'gif']],
  ['Spike initiation', ['adex', 'izhikevich', 'izhikevich_2007', 'qif']],
  ['Conductances and excitability', ['hodgkin_huxley', 'wang_buzsaki', 'morris_lecar', 'fitzhugh_nagumo']],
  ['Resonators', ['resonate_and_fire', 'brf']],
  ['Quantum', ['quantum_emitter']],
  ['Relays', ['parrot']],
];

export class Catalogue {
  constructor(models) {
    this.models = models;
    this.byName = new Map(models.map((m) => [m.name, m]));
  }

  get(name) { return this.byName.get(name); }

  // The families with the models the engine has; any not listed go last.
  families() {
    const listed = new Set(FAMILIES.flatMap(([, names]) => names));
    const out = FAMILIES.map(([title, names]) => [title, names.map((n) => this.get(n)).filter(Boolean)]).filter(([, ms]) => ms.length);
    const rest = this.models.filter((m) => !listed.has(m.name));
    if (rest.length) out.push(['Other', rest]);
    return out;
  }

  parameter(model, name) { return this.get(model)?.parameters.find((p) => p.name === name); }
  receptor(model, name) { return this.get(model)?.receptors.find((r) => r.name === name); }
}

// The panel: a filter, the families, and the chosen model's card.
export function renderCatalogue(root, catalogue, { selected, onSelect, onAdd }) {
  clear(root);
  const filter = el('input', { type: 'search', placeholder: 'Find a model', 'aria-label': 'Find a model', class: 'find' });
  const list = el('div', { class: 'models' });
  const card = el('div', { class: 'model-card' });
  root.append(el('h3', { text: 'Models' }), el('div', { class: 'pad' }, filter), list, card);

  const draw = () => {
    clear(list);
    const q = filter.value.trim().toLowerCase();
    for (const [title, models] of catalogue.families()) {
      const shown = models.filter((m) => !q || m.name.includes(q) || m.summary.toLowerCase().includes(q) || m.citation.toLowerCase().includes(q));
      if (!shown.length) continue;
      list.append(el('div', { class: 'family', text: title }));
      for (const m of shown) {
        list.append(el('button', {
          class: 'model' + (m.name === selected ? ' selected' : ''), 'aria-pressed': m.name === selected ? 'true' : 'false',
          onclick: () => { selected = m.name; onSelect?.(m.name); draw(); show(); },
        }, el('span', { class: 'name', text: m.name }), el('span', { class: 'method', text: m.method })));
      }
    }
  };

  const show = () => {
    clear(card);
    const m = catalogue.get(selected);
    if (!m) { card.append(el('p', { class: 'note', text: 'Choose a model to read what it is and what it takes.' })); return; }
    card.append(
      el('div', { class: 'card-head' },
        el('h4', { text: m.name }),
        el('button', { class: 'add', onclick: () => onAdd?.(m.name), title: 'Add a population of this model' }, '+ population')),
      el('p', { class: 'summary', text: m.summary }),
      el('p', { class: 'citation', text: m.citation }),
      el('p', { class: 'facts' },
        el('span', {}, 'integrated ', el('b', { text: m.method === 'exact' ? 'exactly' : m.method })),
        el('span', {}, 'input in ', el('b', { text: hint(m.input) }))),
      table('Parameters', ['name', 'default', ''], m.parameters.map((p) => [p.name, format(p.default, p.unit), p.doc])),
      table('State', ['component', 'unit', ''], m.components.map((c) => [c.name, hint(c.unit), c.doc])),
      table('Receptors', ['receptor', 'jumps', ''], m.receptors.map((r) => [r.name, `${r.component} (${hint(r.unit)})`, r.doc])),
    );
  };

  filter.addEventListener('input', draw);
  draw();
  show();
  return { select(name) { selected = name; draw(); show(); } };
}

function table(title, head, rows) {
  if (!rows.length) return el('div');
  return el('div', { class: 'table' },
    el('div', { class: 'table-title', text: title }),
    el('table', {},
      el('thead', {}, el('tr', {}, head.map((h) => el('th', { text: h })))),
      el('tbody', {}, rows.map((r) => el('tr', {}, r.map((c, i) => el('td', { class: i === 2 ? 'doc' : null, text: c })))))));
}
