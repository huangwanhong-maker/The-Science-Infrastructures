// The record and its complex (ADR 0067): what the district's record holds and where it comes from; its event complex
// K and refinement K̂, their cells and homology; the network of the moment at the clock's date; how the spacetime is
// drawn; and what was left out of the record, with its reason. One thing at a time, numbers first.

import { cutAt, ensureRecord, setView, state } from '../shared.js';
import { figures, fold, seg, tabs, toggle, hint } from '/src/quiet.js';

export async function render(root, ctx) {
  const { h, busy, failure } = ctx.dom;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The record and its complex'),
    h('span', { class: 'lede' }, 'A district\'s events as a 2-complex: volumes as its edges, the interfaces between them as its faces, events as its vertices.')));
  const body = h('div', {});
  root.append(body);
  busy(body, 'building the record of the district…');
  let R;
  try { R = await ensureRecord(ctx); } catch (e) { return failure(body, e); }
  const n = (x) => (x === undefined || x === null ? '—' : Number(x).toLocaleString('en-GB'));

  // ---- what the record holds ----
  const holds = () => {
    const kinds = Object.entries(R.kinds).filter(([k]) => !k.startsWith('observation')).sort((a, b) => b[1] - a[1]);
    const dropped = Object.entries(R.dropped).sort((a, b) => b[1] - a[1]);
    return h('div', {},
      figures(h, [
        { value: `${R.window.fromDate.slice(0, 4)}–${R.window.toDate.slice(0, 4)}`, label: 'the window' },
        { value: n(R.network.buildings + R.network.streets + R.network.junctions), label: 'volumes (couplings)' },
        { value: n(R.primitives.length), label: 'interfaces (primitives)' },
        { value: n(R.events.length), label: 'events' },
        { value: n(R.censuses.reduce((a, c) => a + c.buildings, 0)), label: 'census observations' },
      ]),
      h('table', { class: 'data' }, h('tbody', {},
        h('tr', {}, h('td', {}, 'buildings'), h('td', {}, n(R.network.buildings)), h('td', { class: 'faint' }, 'the massing\'s 3-cells')),
        h('tr', {}, h('td', {}, 'street spaces'), h('td', {}, n(R.network.streets)), h('td', { class: 'faint' }, 'one for each segment of the walking network')),
        h('tr', {}, h('td', {}, 'junction spaces'), h('td', {}, n(R.network.junctions)), h('td', {})),
        h('tr', {}, h('td', {}, 'party walls'), h('td', {}, n(R.network.partyWalls)), h('td', { class: 'faint' }, 'pairs of buildings sharing a wall')),
        h('tr', {}, h('td', {}, 'frontages'), h('td', {}, n(R.network.frontages)), h('td', { class: 'faint' }, `${n(R.network.facadesFronting)} of ${n(R.network.facades)} outer wall pieces face a street`)),
        h('tr', {}, h('td', {}, 'cross-sections'), h('td', {}, n(R.network.crossSections)), h('td', { class: 'faint' }, 'a street space meeting a junction')))),
      h('h3', {}, 'Events by kind'),
      h('table', { class: 'data' }, h('tbody', {}, ...kinds.map(([k, v]) => h('tr', {}, h('td', {}, k.replace('/', ': ')), h('td', {}, n(v)))))),
      dropped.length ? fold(h, `${n(dropped.reduce((a, [, v]) => a + v, 0))} facts left out of the record, with their reasons`,
        h('table', { class: 'data' }, h('tbody', {}, ...dropped.map(([k, v]) => h('tr', {}, h('td', {}, k), h('td', {}, n(v))))))) : null,
      fold(h, 'where the record comes from', h('ul', {}, ...(R.sources || []).map((s) => h('li', {}, `${s.title || s.id}${s.credit ? ` — ${s.credit}` : ''}${s.retrieved ? `, retrieved ${s.retrieved}` : ''}`)))));
  };

  // ---- the complex ----
  const complex = () => {
    const row = (name, T) => h('tr', {}, h('td', {}, name), h('td', {}, n(T.vertices)), h('td', {}, n(T.edges)), h('td', {}, n(T.faces)), h('td', {}, n(T.euler)),
      h('td', {}, `${n(T.b0)} ${T.b1 === undefined ? '?' : n(T.b1)} ${T.b2 === undefined ? '?' : n(T.b2)}`));
    return h('div', {},
      h('table', { class: 'data' }, h('thead', {}, h('tr', {}, ...['', 'vertices', 'edges', 'faces', 'χ', 'b₀ b₁ b₂'].map((x) => h('th', {}, x)))),
        h('tbody', {}, row('K, the event complex', R.topology.K), row('its refinement (K-hat)', R.topology.refined))),
      hint(h, `K: the events are its vertices (with the window's two ends), each volume's life between two events it attends an edge, each interface's life between two events that meet it a face, attached up one side and down the other. K̂ splits each event into its volumes there, joined by rungs along the interfaces it meets, with ${n(R.topology.triangles)} triangles where three interfaces meet pairwise at one event (two buildings sharing a wall and fronting one street: a hinge).`),
      hint(h, 'A closed shell (b₂) is a set of interfaces whose lenses close up round a hinge between two events that meet them all; a loop (b₁) a circuit of volumes met pairwise at events with no interface across some pair. The Betti numbers are given where the complex is small enough for an exact elimination over ℤ/2.'));
  };

  // ---- the moment ----
  const moment = h('div', {});
  const paintMoment = async () => {
    const t = ctx.clock.instant / 86400000 + 40587;
    moment.replaceChildren(h('p', { class: 'busy' }, 'the network of the moment…'));
    const S = await cutAt(ctx, t);
    if (!S) return;
    const types = {};
    for (const k of S.primitives) types[R.primitives[k].type] = (types[R.primitives[k].type] || 0) + 1;
    const vols = {};
    for (const c of S.couplings) vols[R.couplings[c].type] = (vols[R.couplings[c].type] || 0) + 1;
    moment.replaceChildren(
      figures(h, [{ value: S.date, label: 'the clock\'s date' }, { value: n(S.couplings.length), label: 'volumes standing' }, { value: n(S.primitives.length), label: 'interfaces' }]),
      h('div', { class: 'ld-line' }, Object.entries(vols).map(([k, v]) => `${n(v)} ${k}`).join(' · ')),
      h('div', { class: 'ld-line' }, Object.entries(types).map(([k, v]) => `${n(v)} ${k.replace('-', ' ')}${v === 1 ? '' : 's'}`).join(' · ')),
      hint(h, 'The section of the complex at the date: the volumes as nodes, the interfaces as links between them. It is the boundary between the past of the record and its future, the spin network of the moment of the design (§2.3).'));
  };

  // ---- the drawing ----
  const drawing = () => {
    const V = state.view;
    return h('div', {},
      h('div', { class: 'ld-line' }, 'the complex drawn'),
      seg(h, [['refined', 'refined: each volume at each event'], ['events', 'events: each event one vertex']], V.mode, (m) => setView({ mode: m })),
      h('div', { class: 'ld-line' }, 'time rising'),
      seg(h, [[4, '4 m a year'], [10, '10 m'], [20, '20 m']], V.scale, (s) => setView({ scale: s })),
      seg(h, [['ground', 'from the ground'], ['above', 'above the tallest building']], V.lift, (l) => setView({ lift: l })),
      h('div', { class: 'ld-line' }, 'which cells'),
      ...[['building', 'buildings'], ['street', 'street spaces'], ['junction', 'junction spaces'], ['party-wall', 'party walls'], ['frontage', 'frontages'], ['cross-section', 'cross-sections']]
        .map(([k, name]) => toggle(h, name, V.show[k] !== false, (on) => setView({ show: { ...state.view.show, [k]: on } }))),
      toggle(h, 'the censuses on the buildings\' worldlines', V.censuses, (on) => setView({ censuses: on })),
      h('div', { class: 'ld-line' }, 'the city seen through'),
      seg(h, [[1, 'solid'], [0.35, 'half'], [0.12, 'a ghost']], ctx.cityOpacity ? ctx.cityOpacity() : 1, (a) => ctx.setCityOpacity && ctx.setCityOpacity(a)),
      hint(h, 'The layer "the district\'s events in spacetime" draws it; switch it on in the rail. Pick a bead, a worldline or a sheet on the map for what it is.'));
  };

  const T = tabs(h, [
    { name: 'The record', panel: holds },
    { name: 'The complex', panel: complex },
    { name: 'The moment', panel: () => { paintMoment(); return moment; } },
    { name: 'The drawing', panel: drawing },
  ], { remember: 'event-spacetime' });
  body.replaceChildren(T.el);
}
