// Case Study I in three dimensions: the two configurations X₊ and X₋ of CHBrClF, mirror images,
// with the manuscript's sign τ, the CIP descriptor it corresponds to, and the congruences between
// them; and the two 1,2-dichloroethenes, one constitutional graph, one dipole allowed and one not.

export async function render(root, ctx) {
  const { h, busy, failure } = ctx.dom;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'The stereocentre of Case Study I'),
    h('span', { class: 'lede' }, 'One constitutional graph, the star K₁,₄ with H, F, Cl and Br; two embeddings, mirror images of each other, which no rotation carries one onto the other.')));
  const body = h('div', {});
  root.append(body);
  busy(body);
  let g, cat;
  try {
    [g, cat] = await Promise.all([ctx.ask('grammar', { example: 'stereocentre' }), ctx.ask('catalogue')]);
  } catch (error) {
    return failure(body, error);
  }
  const { view } = await import('/src/draw/molecule.js');
  const geo = g.geometry;
  const canvas = h('div', { class: 'canvas3d' }, h('span', { class: 'tag' }, 'X₊ (left) and X₋ (right), mirror images in the plane between them; drag to turn'));
  body.replaceChildren(h('div', { class: 'grid', style: 'grid-template-columns: minmax(0, 1.6fr) minmax(0, 1fr)' }, canvas,
    h('dl', { class: 'facts panel' },
      h('dt', {}, 'τ(X₊)'), h('dd', {}, `${geo.tauPlus > 0 ? '+1' : '−1'}  →  ${geo.cipPlus}`),
      h('dt', {}, 'τ(X₋)'), h('dd', {}, `${geo.tauMinus > 0 ? '+1' : '−1'}  →  ${geo.cipMinus}`),
      h('dt', {}, 'congruences X₊ → X₋'), h('dd', {}, `${geo.congruences}, with det = ${geo.congruenceDet} (a reflection, no rotation)`),
      h('dt', {}, 'point group'), h('dd', {}, `${geo.pointGroup}: chiral`),
      h('dt', {}, 'τ after F ↔ Cl'), h('dd', {}, `${geo.tauAfterSwap > 0 ? '+1' : '−1'}: relabelling flips τ as a reflection does`),
      h('dt', {}, 'grammar'), h('dd', {}, `${g.terminalClasses} terminal classes, ${g.orbits} orbit of p, ${g.constitutionalClasses} constitutional class`))));
  body.append(h('p', { class: 'note' }, 'τ = sgn det[r_F − r_H, r_Cl − r_H, r_Br − r_H] is defined relative to the ligands\' order. The manuscript\'s order H, F, Cl, Br is increasing CIP priority, and then τ = +1 is exactly the descriptor R: X₊ is (R)-bromochlorofluoromethane.'));
  const dispose = view(canvas, [{ molecule: geo.plus, offset: [-2.4, 0, 0], title: `X₊ (${geo.cipPlus})` }, { molecule: geo.minus, offset: [2.4, 0, 0], title: `X₋ (${geo.cipMinus})` }], { mirror: true, view: [0, 0.12, 1] });
  ctx.onDispose(dispose);

  // The dichloroethenes.
  const z = cat.find((m) => m.id === 'dichloroethene-z'), e = cat.find((m) => m.id === 'dichloroethene-e');
  const pair = h('div', { class: 'canvas3d', style: 'height: 320px' }, h('span', { class: 'tag' }, '(Z) left, (E) right; the arrow is the direction a dipole may take'));
  body.append(h('div', { class: 'section' }, h('h2', {}, 'Two configurations, one graph: 1,2-dichloroethene (App. I.1)'),
    h('div', { class: 'grid', style: 'grid-template-columns: minmax(0, 1.6fr) minmax(0, 1fr)' }, pair,
      h('dl', { class: 'facts panel' },
        h('dt', {}, '(Z)'), h('dd', {}, `${z.pointGroup.schoenflies}: ${z.pointGroup.dipole} fixed direction, a dipole allowed`),
        h('dt', {}, '(E)'), h('dd', {}, `${e.pointGroup.schoenflies}: ${e.pointGroup.dipole} fixed directions, no dipole (centrosymmetric)`),
        h('dt', {}, 'constitutional graph'), h('dd', {}, 'the same for both')))));
  const dz = view(pair, [{ molecule: z, offset: [-2.2, 0, 0] }, { molecule: e, offset: [2.2, 0, 0] }], { dipole: { origin: [-2.2, 0, 0], direction: [0, 1, 0], length: 1.5 }, view: [0, 0.12, 1] });
  ctx.onDispose(dz);
}
