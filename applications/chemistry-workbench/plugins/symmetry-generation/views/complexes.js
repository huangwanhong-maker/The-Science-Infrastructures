// Case Study II's complexes: the six-carbon cycle, the same cycle with one 2-cell attached, and the
// two-vertex complex of the field theory, with their boundary matrices and homology from morphe.

export async function render(root, ctx) {
  const { h, busy, failure } = ctx.dom;
  const { draw, matrix } = ctx.complex;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Three complexes of Case Study II'),
    h('span', { class: 'lede' }, 'A carbon skeleton, a representational choice that fills it, and a Feynman graph of the tetrahedral field theory: three different objects, whose homology morphe computes from their integer boundary matrices.')));
  const body = h('div', {});
  root.append(body);
  busy(body);
  let c;
  try {
    c = await ctx.ask('complexes');
  } catch (error) {
    return failure(body, error);
  }
  const card = (X, title, text, names) => {
    const homology = X.homology.map((g, k) => `H${'₀₁₂₃'[k]} = ${g.group}`).join(', ');
    const B = X.boundaries;
    const edges = X.edges.map((_, i) => `e${i + 1}`), verts = Array.from({ length: X.vertices }, (_, i) => (names ? names[i] : `v${i + 1}`)),
      faces = X.faces.map((_, i) => `f${i + 1}`);
    return h('div', { class: 'panel' }, h('h3', {}, title), h('p', { class: 'caption' }, text),
      draw(X, { width: 300, height: 210, vertexNames: names }),
      h('dl', { class: 'facts' }, h('dt', {}, 'cells'), h('dd', {}, `(${X.counts.join(', ')}), χ = ${X.euler}`), h('dt', {}, 'homology'), h('dd', {}, homology),
        B.length > 1 ? [h('dt', {}, '∂₁∂₂ = 0'), h('dd', {}, X.boundaryOfBoundary ? 'yes' : 'no')] : null),
      h('p', { class: 'caption' }, '∂₁ (vertices × edges)'), matrix(B[0], verts, edges),
      B.length > 1 ? [h('p', { class: 'caption' }, '∂₂ (edges × faces)'), matrix(B[1], edges, faces)] : null);
  };
  body.replaceChildren(h('div', { class: 'grid cols-3' },
    card(c.ring, 'the six-cycle', 'Six carbons, six edges oriented round the ring: one independent cycle, H₁ = ℤ.'),
    card(c.filled, 'with one 2-cell', 'Attaching f along e₁ + … + e₆ kills H₁. It changes the model, not the molecule: aromaticity is not this face.'),
    card(c.dipole, 'the two-vertex complex', 'Two tetrahedral vertices joined factor to factor: each pair of propagators bounds a face, ∂f = e_i − e_j. It is the dual 2-skeleton of S³ cut into two tetrahedra, so H₂ = ℤ³ = ℤ^(4 − 1).', ['u', 'v'])));
  body.append(h('p', { class: 'note' }, 'The two-vertex complex is one of the Feynman graphs of the order-λ² term: of the 24 ways to join two vertices factor to factor, only the 4 automorphisms of the vertex give it; the others close their strands into four or three faces. The Field view lists them all.'));
}
