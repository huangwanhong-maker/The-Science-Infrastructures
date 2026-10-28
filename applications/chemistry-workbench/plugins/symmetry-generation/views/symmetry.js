// Point groups: a molecule in three dimensions with its symmetry elements, its group found by physis
// from the positions of its nuclei, the character table computed exactly by morphe, and the
// responses the group allows (Neumann's principle): a polar vector (a dipole), an axial vector, and
// chirality.

export async function render(root, ctx) {
  const { h, busy, failure, fmt } = ctx.dom;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Point groups and what they allow'),
    h('span', { class: 'lede' }, 'A response that transforms by a representation ρ must lie in the subspace fixed by ρ(H). Its dimension is the mean of the character over the group: that number decides whether a molecule may carry a dipole, and the absence of improper elements whether it is chiral.')));
  const body = h('div', {});
  root.append(body);
  busy(body);
  let cat;
  try {
    cat = await ctx.ask('catalogue');
  } catch (error) {
    return failure(body, error);
  }
  const { view } = await import('/src/draw/molecule.js');
  const stored = () => { try { return sessionStorage.getItem('symmetry-molecule'); } catch { return null; } };
  let current = (ctx.args && ctx.args[0]) || stored() || 'dichloroethene-e';
  const keys = h('div', { class: 'keys' });
  const show = h('div', {});
  body.replaceChildren(keys, show);
  let dispose = null;
  const draw = () => {
    if (dispose) dispose();
    for (const k of keys.children) k.classList.toggle('on', k.dataset.id === current);
    const m = cat.find((x) => x.id === current) || cat[0];
    const pg = m.pointGroup;
    const canvas = h('div', { class: 'canvas3d' }, h('span', { class: 'tag' }, `${m.name} · ${m.formula} · axes violet, mirrors apricot`));
    const facts = h('dl', { class: 'facts panel' },
      h('dt', {}, 'point group'), h('dd', {}, h('b', {}, pg.schoenflies), pg.linear ? ' (on a line: infinite)' : `, order ${pg.order}`),
      h('dt', {}, 'chiral'), h('dd', {}, pg.chiral ? 'yes: no improper element' : 'no'),
      h('dt', {}, 'polar vectors fixed'), h('dd', {}, `${pg.dipole}: ${pg.dipole > 0 ? 'a permanent dipole is allowed' : 'no permanent dipole'}`),
      h('dt', {}, 'axial vectors fixed'), h('dd', {}, String(pg.axial)));
    const parts = [h('div', { class: 'grid', style: 'grid-template-columns: minmax(0, 1.6fr) minmax(0, 1fr)' }, canvas, h('div', {}, facts, m.note ? h('p', { class: 'caption' }, m.note) : null))];
    if (pg.characters) {
      const C = pg.characters;
      const head = h('tr', {}, h('th', {}, ''), ...C.classes.map((c) => h('th', {}, `${c.size > 1 ? c.size : ''}${c.operation}`)));
      const rows = C.rows.map((row, i) => h('tr', {}, h('td', {}, `χ${i + 1}`), ...row.map((v) => h('td', { class: 're', title: v.exact }, Math.abs(v.im) > 1e-9 ? `${fmt(v.re, 3)}${v.im >= 0 ? '+' : '−'}${fmt(Math.abs(v.im), 3)}i` : fmt(Math.round(v.re * 1e9) / 1e9, 4)))));
      const vec = h('tr', { class: 'hl' }, h('td', {}, 'vector'), ...C.vector.map((x) => h('td', {}, fmt(Math.round(x * 1e9) / 1e9, 4))));
      const ax = h('tr', {}, h('td', {}, 'axial'), ...C.axial.map((x) => h('td', {}, fmt(Math.round(x * 1e9) / 1e9, 4))));
      parts.push(h('div', { class: 'section' }, h('h2', {}, 'Character table'),
        h('p', { class: 'caption' }, 'The irreducible characters of the group found, computed exactly by Dixon and Schneider\'s method in morphe (hover a value for its exact form), with the characters of the vector and axial representations, tr Q and det Q · tr Q, whose means over the group give the fixed dimensions above.'),
        h('div', { style: 'overflow-x:auto' }, h('table', { class: 'data chartable' }, head, ...rows, vec, ax))));
    }
    show.replaceChildren(...parts);
    dispose = view(canvas, [{ molecule: m }], { elements: pg });
  };
  for (const m of cat) {
    const b = h('button', { class: 'key', 'data-id': m.id }, m.name.replace(/, X.*$/, ' X' + (m.id.endsWith('plus') ? '₊' : m.id.endsWith('minus') ? '₋' : '')));
    b.addEventListener('click', () => { current = m.id; try { sessionStorage.setItem('symmetry-molecule', current); } catch { /* private window */ } draw(); });
    keys.append(b);
  }
  draw();
  ctx.onDispose(() => dispose && dispose());
}
