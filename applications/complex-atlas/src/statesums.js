// The State sums panel: the Turaev-Viro invariant of a triangulated closed 3-manifold. Each edge is
// coloured by a spin j in {0, 1/2, ..., (r − 2)/2}; a colouring weighs each edge by its quantum
// dimension and each tetrahedron by its quantum 6j symbol at q = e^{iπ/r}; the sum, normalised by
// the total quantum dimension at each vertex, does not depend on the triangulation (Turaev and Viro
// 1992). It is the Ponzano-Regge model of three-dimensional quantum gravity, made finite.

import { el, clear } from './dom.js';
import { number } from './format.js';

export const form = { r: 4 };

// Whether the complex is a triangulation of dimension 3: its facets are tetrahedra, or it is the
// 2-complex of a 4-coloured graph, whose bicoloured cycles are the edges of the triangulation it encodes.
const isColoured = (c) => c?.coloured_graph?.colours === 4;
export const isThreeDimensional = (c) => isColoured(c) || (Array.isArray(c?.simplices) && c.simplices.length > 0 && c.simplices.every((s) => s.length === 4));

export function renderSums(host, { complex, result, running, error, onRun }) {
  clear(host);
  if (!complex) return;
  host.append(el('p', { class: 'note', text: 'The Turaev–Viro state sum: every edge coloured by a spin up to (r − 2)/2, every tetrahedron weighed by its quantum 6j symbol at q = e^{iπ/r}. It is a topological invariant of the 3-manifold, the Ponzano–Regge model of three-dimensional gravity made finite at a root of unity.' }));
  if (!isThreeDimensional(complex)) {
    host.append(el('p', { class: 'note', text: 'It needs a triangulated closed 3-manifold, given by its tetrahedra or by a 4-coloured graph: try ∂Δ⁴ or a melonic graph with D = 3, or facets of four vertices each.' }));
    return;
  }
  const edges = isColoured(complex) ? complex.faces.length : complex.edges.length;
  const r = el('input', { type: 'number', min: 3, max: 9, value: form.r, oninput: () => { form.r = Number(r.value); size.textContent = estimate(); } });
  const estimate = () => `up to ${number((form.r - 1) ** edges, 3)} colourings of ${edges} edges`;
  const size = el('span', { class: 'note', text: estimate() });
  host.append(el('div', { class: 'row' }, el('label', { class: 'field' }, el('span', { class: 'label', text: 'r, the root of unity' }), r), size),
    el('div', { class: 'row' }, el('button', { class: 'button primary', disabled: running ? true : null, onclick: () => onRun(form.r) }, running ? 'Summing…' : 'Sum'), error ? el('span', { class: 'error', text: error }) : null));
  if (result) {
    const tile = (label, value) => el('div', { class: 'tile' }, el('div', { class: 'tile-label', text: label }), el('div', { class: 'tile-value', text: value }));
    host.append(el('div', { class: 'tiles' }, tile(`TV at r = ${result.r}`, number(result.value, 8)), tile('admissible colourings', number(result.colourings)),
      ...(result.sphere !== undefined ? [tile('the 3-sphere, w⁻²', number(result.sphere, 8))] : [])));
  }
}
