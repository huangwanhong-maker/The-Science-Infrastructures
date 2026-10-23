// The invariants of the complex, as the engine computes them: its cells and Euler characteristic,
// homology with generators (each can be drawn on the complex), cohomology, Betti numbers over
// finite fields, the fundamental group as a presentation with its generators as loops, the
// surface it is, whether it is a manifold, whether it collapses, and, for the dual of a coloured
// graph, the data of the tensor model.

import { el, clear } from './dom.js';
import { groupText, presentationText, wordText, surfaceText, normalFormNames, countsText, chainText, superscript } from './format.js';

const H = (k) => `H${'₀₁₂₃₄₅₆₇₈₉'[k] ?? k}`;

export function renderInvariants(host, { complex, analysis, names, selection, onSelect, state }) {
  clear(host);
  if (!complex) return;
  const chip = (label, sel, title) => {
    const on = selection && sel && JSON.stringify(selection.key) === JSON.stringify(sel.key);
    return el('button', { class: `chip${on ? ' on' : ''}`, title, onclick: () => onSelect(on ? null : sel) }, label);
  };
  if (state === 'loading') host.append(el('p', { class: 'note', text: 'Computing the invariants…' }));
  if (state && state.error) host.append(el('p', { class: 'error', text: state.error }));
  const a = analysis;
  const counts = a?.counts || null;
  host.append(el('h4', { text: 'The cells' }));
  host.append(el('div', { class: 'counts', text: countsText(counts || [complex.vertices, complex.edges.length, (complex.faces || []).length, ...(complex.cells || []).map((c) => c.length)].filter((n, k, all) => k < 2 || n || all.slice(k).some(Boolean))) }));
  if (a) host.append(el('div', { class: 'euler' }, 'Euler characteristic χ = ', el('b', { text: String(a.euler).replace('-', '−') }), a.components > 1 ? `, in ${a.components} components` : ''));
  if (!a) {
    if (state !== 'loading' && !(state && state.error)) host.append(el('p', { class: 'note', text: 'The invariants are computed exactly by the engine, the foundations\' C++ (morphe\'s topology). Build it to see them; the views on the left need no engine.' }));
    return;
  }

  // Homology.
  host.append(el('h4', { text: 'Homology' }));
  const table = el('table', { class: 'data' }, el('tr', {}, el('th', { text: '' }), el('th', { text: 'over ℤ' }), el('th', { text: 'generators' })));
  for (const h of a.homology || []) {
    const gens = el('div', { class: 'gens' });
    (h.generators || []).forEach((g, i) => {
      const label = g.order ? `ℤ/${g.order}` : `ℤ ${i + 1}`;
      gens.append(chip(label, { kind: 'chain', dim: h.dimension, chain: g.chain, key: ['h', h.dimension, i] }, chainText(g.chain, h.dimension, names)));
    });
    table.append(el('tr', {}, el('td', { text: H(h.dimension) }), el('td', { class: 'group', text: groupText(h.rank, h.torsion) }), el('td', {}, gens)));
  }
  host.append(table);
  if (a.cohomology?.length) host.append(el('p', { class: 'fact', text: `Cohomology: ${a.cohomology.map((h) => `H${superscript(h.dimension)} = ${groupText(h.rank, h.torsion)}`).join(', ')}` }));
  if (a.betti) {
    const fields = Object.entries(a.betti).map(([p, b]) => `over 𝔽${'₀₁₂₃₄₅₆₇₈₉'[p] ?? p}: ${b.join(', ')}`);
    if (fields.length) host.append(el('p', { class: 'fact', text: `Betti numbers ${fields.join('; ')}` }));
  }

  // The fundamental group.
  const pg = a.fundamental_group;
  if (pg) {
    host.append(el('h4', { text: 'The fundamental group' }));
    const big = pg.generators.length > 24 || pg.relators.length > 24;
    host.append(el('div', { class: 'presentation' }, big
      ? el('span', { class: 'faint', text: `From a spanning tree: ${pg.generators.length} generators, ${pg.relators.length} relators.` })
      : presentationText(pg.generators, pg.relators)));
    if (pg.simplified) host.append(el('div', { class: 'presentation' }, el('span', { class: 'faint', text: 'simplified ' }), presentationText(pg.simplified.generators, pg.simplified.relators)));
    if (pg.abelianisation) host.append(el('p', { class: 'fact', text: `Abelianised: ${groupText(pg.abelianisation.rank, pg.abelianisation.torsion)}` }));
    // The loops of the simplified generators when there are any (the generators that remain
    // after Tietze's moves are among the original ones), else of the raw ones.
    const shown = pg.simplified?.loops?.length ? pg.simplified : pg;
    if (shown.loops?.length) {
      const gens = el('div', { class: 'gens' });
      shown.loops.slice(0, 48).forEach((loop, i) => gens.append(chip(shown.generators[i] ?? `x${i + 1}`, { kind: 'loop', word: loop, key: ['pi', i] }, `the loop ${wordText(loop, names, { powers: false })}`)));
      host.append(el('p', { class: 'note', text: `Each ${shown === pg ? '' : 'remaining '}generator is a loop at the base vertex: through the spanning tree, along its edge, and back.` }), gens);
    }
  }

  // What kind of space it is.
  const facts = el('dl', { class: 'facts' });
  const fact = (k, v) => facts.append(el('dt', { text: k }), el('dd', {}, v));
  if (a.surface) {
    const s = a.surface;
    fact('surface', surfaceText(s));
    fact('normal form', `${wordText(s.normal_form || [], normalFormNames(s))}${s.boundary_components ? `, and ${s.boundary_components} boundary circles` : ''}`);
  }
  if (a.manifold) fact('manifold', a.manifold.manifold ? `a combinatorial manifold of dimension ${a.manifold.dimension}` : a.manifold.pseudomanifold ? `a pseudomanifold of dimension ${a.manifold.dimension}, not a manifold` : 'not a pseudomanifold');
  if (a.collapse) {
    const c = a.collapse;
    fact('collapse', c.collapsible ? `collapsible, to a point in ${c.steps} elementary collapses` : `not collapsible: ${c.steps} elementary collapses leave ${countsText(c.remaining)}`);
  }
  const g = a.coloured_graph;
  if (g) {
    fact('coloured graph', `${g.colours} colours, ${g.vertices} vertices, ${g.faces} faces, ${g.bubbles} bubbles`);
    fact('jackets', g.jackets?.length ? `genera ${g.jackets.join(', ')}` : '—');
    fact('degree', `ω = ${g.degree}${g.melonic ? ', melonic' : ''}`);
  }
  if (facts.childNodes.length) host.append(el('h4', { text: 'The space' }), facts);
}
