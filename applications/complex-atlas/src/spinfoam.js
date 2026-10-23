// The Spin foam panel. A closed triangulated 3-manifold (or a 4-coloured graph, the Feynman graph
// of a group field theory in three dimensions) carries a spin foam on its dual 2-complex: a spin j
// on each 2-cell (dual to an edge of the triangulation, j + ½ its length), and at each vertex
// (dual to a tetrahedron) the amplitude of its boundary spin network, the 6j symbol of the six
// spins, quantum at q = exp(iπ/r) (Ponzano and Regge 1968; Turaev and Viro 1992). Summed over all
// spins the amplitudes give the Turaev–Viro invariant; at large spins the 6j symbol oscillates
// as the cosine of the Regge action of the tetrahedron (the Ponzano–Regge formula), the bridge
// from the spin foam to discrete general relativity.
//
// The engine samples spins (Metropolis with weight |A|, the sign kept) or sets them uniform; the
// page draws the foam with its faces coloured by spin and its vertices by their amplitudes, and a
// vertex opens its spin network and its quantum tetrahedron.

import { el, clear } from './dom.js';
import { histogram } from './plot.js';
import { spinNetworkSvg, TETRAHEDRON_EDGES } from './graphs.js';
import { number, spinText } from './format.js';
import { sequentialColour, cssColour } from '../../../visualization/js/index.js';
import { TetrahedronView } from './views/tetrahedron.js';

export const foamForm = { r: 5, mode: 'sample', two_j: 2, sweeps: 2000, seed: 1 };

// Whether a complex carries a spin foam here: a closed triangulated 3-manifold given by its
// tetrahedra, or the 2-complex of a 4-coloured graph.
export const foamable = (c) => Boolean(c && (c.coloured_graph?.colours === 4 || (Array.isArray(c.simplices) && c.simplices.length > 0 && c.simplices.every((s) => s.length === 4))));

// The colour of a spin among spins up to max.
export const spinColour = (twoJ, max) => sequentialColour(0.2 + 0.8 * (max ? twoJ / max : 0));

// What the realisation shows of a spin foam.
export function foamOverlay(result) {
  const max = Math.max(1, ...result.face_spins);
  const amps = result.vertices.map((v) => Math.abs(v.amplitude ?? 0)), top = Math.max(1e-300, ...amps);
  return {
    faceColours: result.face_spins.map((j) => spinColour(j, max)),
    faceLabels: result.face_spins.map((j) => `j = ${spinText(j)}`),
    vertexScales: amps.map((a) => 0.9 + 1.5 * Math.sqrt(a / top)),
    vertexColours: result.vertices.map((v) => ((v.amplitude ?? 0) < 0 ? [0.725, 0.635, 0.878] : [0.949, 0.761, 0.486])),
    maxSpin: max,
  };
}

export function renderFoam(host, { complex, foam, active, running, error, onRun, onBack, onVertex, selected }) {
  clear(host);
  if (!complex) return;
  host.append(el('p', { class: 'note', text: 'A spin foam: a spin on each 2-cell of the dual 2-complex (an edge of the triangulation, of length j + ½), and at each vertex (a tetrahedron) the 6j symbol of its six spins, quantum at q = exp(iπ/r). Summed over spins: the Turaev–Viro invariant; at large spins: the Regge action of three-dimensional gravity.' }));
  if (!foam && !foamable(complex)) {
    host.append(el('p', { class: 'note', text: 'It needs a closed triangulated 3-manifold (its tetrahedra) or a 4-coloured graph: try ∂Δ⁴ among the examples, the suspended 3-sphere or a melonic graph in the catalogue, or a graph of Boulatov\'s GFT in the Generation mode.' }));
    return;
  }
  if (!active) {
    const field = (label, input) => el('label', { class: 'field' }, el('span', { class: 'label', text: label }), input);
    const num = (key, attrs) => { const i = el('input', { type: 'number', value: foamForm[key], ...attrs, oninput: () => { foamForm[key] = Number(i.value); } }); return i; };
    const mode = el('select', { onchange: () => { foamForm.mode = mode.value; } },
      ...[['sample', 'sampled (Metropolis on |A|)'], ['uniform', 'uniform spin']].map(([v, t]) => el('option', { value: v, selected: foamForm.mode === v ? true : null }, t)));
    host.append(el('div', { class: 'grid' },
      field('r, the root of unity', num('r', { min: 3, max: 9 })), field('spins', mode),
      field('uniform 2j', num('two_j', { min: 0, max: 7 })), field('sweeps', num('sweeps', { min: 1, step: 500 })),
      field('seed', num('seed', { min: 0 }))),
    el('div', { class: 'row' }, el('button', { class: 'button primary', disabled: running ? true : null, onclick: () => onRun({ ...foamForm }) }, running ? 'Sampling…' : 'Lay the spin foam'),
      error ? el('span', { class: 'error', text: error }) : null));
  }
  if (!foam) return;
  const r = foam;
  if (active) host.append(el('div', { class: 'row' }, el('button', { class: 'button', onclick: onBack }, '← The triangulation'),
    el('span', { class: 'note', text: `r = ${r.r ?? foamForm.r}; click a vertex for its tetrahedron` })));
  const tile = (label, value) => el('div', { class: 'tile' }, el('div', { class: 'tile-label', text: label }), el('div', { class: 'tile-value', text: value }));
  host.append(el('div', { class: 'tiles' },
    tile('amplitude of these spins', number(r.amplitude, 5)),
    tile('Turaev–Viro, all spins', r.tv !== null && r.tv !== undefined ? number(r.tv, 6) : '—'),
    tile('vertices · faces', `${r.vertices.length} · ${r.face_spins.length}`),
    tile('sign of A, mean', r.statistics?.sign_mean !== undefined ? number(r.statistics.sign_mean, 3) : r.sign > 0 ? '+' : '−')));
  if (r.statistics) host.append(el('p', { class: 'fact', text: `${number(r.statistics.samples)} samples, acceptance ${number(r.statistics.acceptance, 3)}.` }));
  host.append(el('h4', { text: 'The spins' }), histogram(r.face_spins, { label: (j) => `j = ${spinText(j)}` }));
  // A legend of the faces' colours.
  const max = Math.max(1, ...r.face_spins);
  host.append(el('div', { class: 'spin-ramp' }, ...Array.from({ length: max + 1 }, (_, j) => el('span', { style: `background:${cssColour(spinColour(j, max))}`, title: `j = ${spinText(j)}` }, spinText(j)))));
  if (r.regge) {
    const d = r.regge.deficits || [];
    host.append(el('h4', { text: 'Its geometry (Regge)' }),
      el('p', { class: 'fact', text: `With lengths j + ½: Regge action Σ ℓ ε = ${number(r.regge.action, 5)}; deficit angles from ${number(Math.min(...d), 4)} to ${number(Math.max(...d), 4)} rad.` }),
      el('p', { class: 'note', text: 'A positive deficit angle is positive curvature, concentrated on the edge; flat space has none.' }));
  }
  host.append(el('h4', { text: 'Its vertices' }));
  const gens = el('div', { class: 'gens' });
  r.vertices.forEach((v, i) => gens.append(el('button', { class: `chip${selected === i ? ' on' : ''}`, title: `{6j} = ${number(v.amplitude, 5)}`, onclick: () => onVertex(selected === i ? null : i) }, `v${i} · ${number(v.amplitude, 3)}`)));
  host.append(gens);
}

// The inset of a vertex: its spin network, its quantum tetrahedron, its amplitude.
let tetraView = null;
export function renderInset(host, { vertex, index, max, onClose }) {
  if (tetraView) { tetraView.dispose(); tetraView = null; }
  clear(host);
  if (!vertex) { host.hidden = true; return; }
  host.hidden = false;
  const geo = vertex.tetrahedron || {};
  const tetra = el('div', { class: 'tetra' });
  host.append(
    el('div', { class: 'inset-head' }, el('b', { text: `Vertex v${index}: a quantum tetrahedron` }), el('button', { class: 'close', title: 'Close', onclick: onClose }, '×')),
    el('div', { class: 'inset-body' },
      el('div', {}, el('div', { class: 'inset-label', text: 'its boundary spin network' }),
        spinNetworkSvg(vertex.spins, { size: 190, spinText: (j) => spinText(j), colourOf: (j) => cssColour(spinColour(j, max)) })),
      el('div', {}, el('div', { class: 'inset-label', text: geo.euclidean === false ? 'its tetrahedron (not Euclidean: flattened)' : 'its tetrahedron, edges j + ½' }), tetra)),
    el('dl', { class: 'facts' },
      el('dt', { text: 'spins' }), el('dd', { text: TETRAHEDRON_EDGES.map(([a, b], k) => `${a}${b}: ${spinText(vertex.spins[k])}`).join(' · ') }),
      el('dt', { text: '{6j}' }), el('dd', { text: number(vertex.amplitude, 6) }),
      el('dt', { text: 'Ponzano–Regge' }), el('dd', { text: vertex.asymptotic !== null && vertex.asymptotic !== undefined ? `${number(vertex.asymptotic, 6)}  (cos(Σ (j + ½) θ + π/4) / √(12πV), exact as j → ∞)` : 'no Euclidean tetrahedron: the 6j decays' }),
      el('dt', { text: 'volume' }), el('dd', { text: number(geo.volume, 5) }),
      el('dt', { text: 'dihedral angles' }), el('dd', { text: (geo.dihedral || []).map((t) => number(t, 4)).join(', ') })));
  if (geo.lengths) {
    tetraView = new TetrahedronView(tetra);
    tetraView.show(geo.lengths, { colours: vertex.spins.map((j) => spinColour(j, max)), labels: vertex.spins.map((j) => spinText(j + 1)) });  // the length ℓ = j + ½, from twice the spin
  }
}
