// The complex: its counts and Euler characteristic, its homology and fundamental group (morphe's,
// computed by the engine), the loops of H₁ told as the beings and occasions they pass, and the
// closed shells of H₂ with the Bianchi identity on each; and the complex saved for the Complex Atlas.

import { el, clear } from './dom.js';
import { shortName, yearText, loopWords, signText, plural, listText } from './format.js';
import { renderModels } from './models.js';

const SUP = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
const groupText = (text) => String(text).replace(/\^(\d+)/g, (_, n) => [...n].map((d) => SUP[d]).join(''));

export function renderComplex(host, data, act) {
  clear(host);
  const inv = data.invariants, [V, E, F] = inv.counts;
  const H = (k) => inv.homology.find((h) => h.dimension === k);
  const save = () => {
    const blob = new Blob([JSON.stringify(data.complex.atlas, null, 1)], { type: 'application/json' });
    const a = el('a', { href: URL.createObjectURL(blob), download: `${data.canon.id}.complex.json` });
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  };
  const occ = data.occasions;
  host.append(
    el('h2', { class: 'view-title', text: 'The complex' }),
    el('p', { class: 'lede narrow', text: 'The occasions are the vertices; the life of each being between two of its occasions is an edge; each tie between two beings is a chain of lenses, one 2-cell between each two of their meetings, its boundary running forward along one life and back along the other. The invariants are exact, from morphe.' }),
    el('div', { class: 'complex-grid' },
      el('section', { class: 'invariants' },
        el('div', { class: 'figures big' },
          figure(V, 'vertices', 'occasions'), figure(E, 'edges', 'lives between occasions'), figure(F, 'faces', 'lenses of ties'),
          figure(String(inv.euler).replace('-', '−'), 'χ', 'V − E + F')),
        el('table', { class: 'groups' }, el('tbody', {},
          row('H₀', groupText(H(0)?.text ?? '—'), plural(inv.components, 'connected history', 'connected histories')),
          row('H₁', groupText(H(1)?.text ?? '—'), 'loops of encounters no recorded tie fills'),
          row('H₂', groupText(H(2)?.text ?? '—'), 'closed shells of ties'),
          row('π₁', groupText(inv.pi1.text), `from ${plural(inv.pi1.generators, 'generator')} and ${plural(inv.pi1.relators, 'relator')}, simplified to ${inv.pi1.simplified_generators} and ${inv.pi1.simplified_relators}`))),
        el('p', { class: 'note', text: 'H₁ is where flat connections wind: a circle of lives meeting at occasions, with no tie recorded across it. Each is a silence of the record, and a place to look in the sources next.' }),
        el('div', { class: 'row' }, el('button', { class: 'button', onclick: save }, 'Save the complex for the Complex Atlas'),
          el('span', { class: 'faint', text: 'JSON with positions, to open in the Complex Atlas as your own complex' }))),
      el('section', { class: 'loops' },
        el('h3', { text: `The loops of H₁ (${inv.loops.length})` }),
        inv.loops.length ? el('ol', { class: 'loop-list' }, ...inv.loops.map((loop, i) => {
          const w = loopWords(loop, data);
          return el('li', { class: data.view?.loop === i ? 'on' : '' },
            el('p', { class: 'loop-path', text: w.path }),
            w.turns.length ? el('p', { class: 'loop-turns' }, el('span', { class: 'where', text: 'where ' }), listText(w.turns)) : null,
            el('span', { class: 'row-links' }, el('button', { class: 'link', onclick: () => act.showLoop(i, 'foam') }, 'Show it in the foam'), ' · ',
              el('button', { class: 'link', onclick: () => act.showLoop(i, 'loom') }, 'on the loom')));
        })) : el('p', { class: 'note', text: 'None: every circle of encounters is filled by ties.' }),
        el('h3', { text: `Closed shells (${inv.shells.length} minimal, of H₂ of rank ${inv.shell_rank})` }),
        el('p', { class: 'note', text: 'Three beings pairwise tied, all present at two consecutive meetings of their three ties: three lenses closing into a sphere. On each, the Bianchi identity of ℤ₂: the sign of the triangle at the second meeting is its sign at the first, times −1 for each lens that flips.' }),
        inv.shells.length ? el('table', { class: 'shells' },
          el('thead', {}, el('tr', {}, el('th', { text: 'Beings' }), el('th', { text: 'Between' }), el('th', { class: 'num', text: 'Before' }), el('th', { class: 'num', text: 'Flips' }), el('th', { class: 'num', text: 'After' }), el('th', { text: '' }))),
          el('tbody', {}, ...inv.shells.map((s) => {
            const ok = s.product_to === s.product_from * (s.flips % 2 ? -1 : 1);
            return el('tr', {},
              el('td', { text: listText(s.beings.map((b) => data.beings[b].short)) }),
              el('td', {}, linkTo(act, occ[s.from]), ' → ', linkTo(act, occ[s.to])),
              el('td', { class: `num sign ${s.product_from < 0 ? 'neg' : 'pos'}`, text: signText(s.product_from) }),
              el('td', { class: 'num', text: String(s.flips) }),
              el('td', { class: `num sign ${s.product_to < 0 ? 'neg' : 'pos'}`, text: signText(s.product_to) }),
              el('td', { class: ok ? 'faint' : 'error', text: ok ? '✓' : 'fails' }));
          }))) : el('p', { class: 'note', text: 'None.' }))),
    renderModels(data, act));
}

function figure(value, name, about) {
  return el('div', { class: 'figure' }, el('span', { class: 'value', text: String(value) }), el('span', { class: 'name', text: name }), el('span', { class: 'about', text: about }));
}
function row(name, value, about) {
  return el('tr', {}, el('th', { text: name }), el('td', { class: 'group', text: value }), el('td', { class: 'faint', text: about }));
}
function linkTo(act, o) {
  return el('button', { class: 'link', title: o.label, onclick: () => act.pick({ kind: 'occasion', id: o.id }) }, `${o.label.length > 34 ? `${o.label.slice(0, 32)}…` : o.label} (${yearText(o.date)})`);
}
