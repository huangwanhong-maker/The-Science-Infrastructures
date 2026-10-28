// Formalism II: the tetrahedral interaction on a real field over a finite group. Its value on a field
// the reader sets (over ℤ₂, eight numbers), its translation symmetries (G × G, not only the diagonal),
// its values along the eight characters of ℤ₂³ and the odd order parameters, with and without the
// gauge projection; the action along the characters (−1)^a and (−1)^b; and the Wick expansion: the
// Feynman graphs at order n, each with its 2-complex, its weight and its homology.

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Formalism II: the tetrahedral field'),
    h('span', { class: 'lede' }, 'V[ψ] = Σ ψ(a,b,c) ψ(c,d,e) ψ(e,b,f) ψ(f,d,a), Boulatov\'s vertex: four triangles glued into a tetrahedron. Over ℤ₂ a field is eight numbers.')));
  const psi = [1, 1, 1, 1, 1, 1, 1, 1];
  const fieldBox = h('div', {});
  const landauBox = h('div', {});
  const wickBox = h('div', {});
  root.append(fieldBox, landauBox, wickBox);

  const presets = [['constant 1', [1, 1, 1, 1, 1, 1, 1, 1]], ['(−1)^a, the manuscript\'s ψ₁', [1, 1, 1, 1, -1, -1, -1, -1]], ['(−1)^b', [1, 1, -1, -1, 1, 1, -1, -1]], ['a negative value', [-2, -2, -2, -2, 2, -2, 1, -2]]];
  const showField = async () => {
    busy(fieldBox);
    let f;
    try {
      f = await ctx.ask('field', { group: 'Z2', psi });
    } catch (error) {
      return failure(fieldBox, error);
    }
    const inputs = f.points.map((p, i) => {
      const inp = h('input', { type: 'number', value: psi[i], step: 1, style: 'width: 4.2em' });
      inp.addEventListener('change', () => { psi[i] = Math.round(+inp.value || 0); showField(); });
      return h('label', { style: 'display:inline-flex; gap:5px; align-items:center; margin: 2px 10px 2px 0' }, h('span', { class: 'mono faint' }, `ψ${p}`), inp);
    });
    const keys = h('div', { class: 'keys' }, ...presets.map(([name, v]) => {
      const b = h('button', { class: 'key' }, name);
      b.addEventListener('click', () => { v.forEach((x, i) => { psi[i] = x; }); showField(); });
      return b;
    }));
    const translations = f.translations.map((t) => h('span', { class: `chip ${t.invariant ? 'pass' : ''}`, style: 'margin: 2px' }, `(${t.g.join('')})`));
    const modeRows = f.modes.map((m) => h('tr', { class: m.s.join('') === '100' ? 'hl' : '' },
      h('td', {}, `(−1)^{${['a', 'b', 'c'].filter((_, i) => m.s[i]).join('+') || '0'}}`), h('td', {}, m.odd ? 'odd' : 'even'),
      h('td', {}, fmt(m.V)), h('td', {}, fmt(m.eta)), h('td', {}, fmt(m.etaProjected))));
    fieldBox.replaceChildren(h('div', { class: 'grid cols-2' },
      h('div', { class: 'panel' }, h('h3', {}, 'the field ψ on ℤ₂³'), keys, h('div', {}, ...inputs),
        h('p', { style: 'font-size:1.15rem; margin-top:8px' }, 'V[ψ] = ', h('b', {}, fmt(f.value))),
        h('p', { class: 'caption' }, `Slots joined by strands: {${f.slotClasses.map((c) => c.join(', ')).join('}, {')}}. V is invariant under translating the classes independently: of the eight slot translations (g₁g₂g₃), those that keep V on a generic field are marked, (g, h, g): the group G × G, of which the manuscript's simultaneous translation (111) is the diagonal.`),
        h('div', {}, ...translations)),
      h('div', { class: 'panel' }, h('h3', {}, 'characters and order parameters'),
        h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'character χ'), h('th', {}, 'under U_p'), h('th', {}, 'V[χ]'), h('th', {}, 'η_χ[ψ]'), h('th', {}, 'η_χ[Pψ]')), ...modeRows),
        h('p', { class: 'caption' }, 'η_χ[ψ] = 2⁻³ Σ χ(x) ψ(x). The highlighted row is the manuscript\'s η: V vanishes identically along it. The last column applies the gauge projection P: over the abelian ℤ₂ it annihilates every odd η.'))));
  };

  const showLandau = async () => {
    let L;
    try {
      L = await ctx.ask('landau', { mu: -1, lambda: 0.05 });
    } catch (error) {
      return failure(landauBox, error);
    }
    // (−1)^{a+b+c} has the same V as (−1)^b, so its curve would lie on it: two are drawn.
    const curves = L.curves.slice(0, 2);
    const colors = [css('--series-2'), css('--series-1')];
    landauBox.replaceChildren(h('div', { class: 'section' }, h('h2', {}, 'The action along a character: S(t) = (μ/2)·8t² + λ V[χ] t⁴ (μ = −1, λ = 0.05)'),
      lines(curves.map((c, i) => ({ name: c.mode, x: c.t, y: c.S, color: colors[i], dash: c.V === 0 ? '6 4' : null })), { xlabel: 't', ylabel: 'S', height: 240 }),
      legend(curves.map((c, i) => ({ name: `${c.mode}: V = ${c.V}${c.minimum ? `, minima ±${fmt(c.minimum, 3)}` : ', no quartic term'}`, color: colors[i] }))),
      h('p', { class: 'note' }, 'Along (−1)^a the interaction contributes nothing, so §5.3\'s β cannot come from it; along the odd character (−1)^b it gives 64λt⁴, and S = 4μη² + 64λη⁴ is the Landau potential with α = 4μ, β = 64λ.')));
  };

  let n = 2, group = 'Z2';
  const showWick = async () => {
    busy(wickBox, 'enumerating the Wick pairings…');
    let w;
    try {
      w = await ctx.ask('wick', { n, group });
    } catch (error) {
      return failure(wickBox, error);
    }
    const nk = h('div', { class: 'keys' }, ...[1, 2, 3].map((k) => {
      const b = h('button', { class: `key ${k === n ? 'on' : ''}` }, `order ${k}`);
      b.addEventListener('click', () => { n = k; showWick(); });
      return b;
    }), ...['1', 'Z2', 'Z3', 'S3'].map((g) => {
      const b = h('button', { class: `key ${g === group ? 'on' : ''}` }, g === '1' ? 'G = 1' : `G = ${g.replace('Z', 'ℤ').replace('S3', 'S₃')}`);
      b.addEventListener('click', () => { group = g; showWick(); });
      return b;
    }));
    const classes = [...w.classes].sort((a, b) => b.faces - a.faces || b.pairings - a.pairings);
    const rows = classes.map((c, i) => {
      const tr = h('tr', { style: 'cursor:pointer' }, h('td', {}, String(i + 1)), h('td', {}, c.connected ? 'connected' : 'disconnected'), h('td', {}, String(c.faces)),
        h('td', {}, String(c.pairings)), h('td', {}, String(c.automorphisms)), h('td', {}, c.complex.homology.map((g) => g.group).join(', ')), h('td', {}, c.amplitude), h('td', {}, c.bf));
      tr.addEventListener('click', () => {
        for (const r of tbody.children) r.classList.remove('hl');
        tr.classList.add('hl');
        drawing.replaceChildren(h('h3', {}, `class ${i + 1}: ${c.faces} faces`), ctx.complex.draw(c.complex, { width: 320, height: 230, vertexNames: Array.from({ length: c.complex.vertices }, (_, k) => `u${k + 1}`) }));
      });
      return tr;
    });
    const tbody = h('tbody', {}, ...rows);
    const drawing = h('div', { class: 'panel' }, h('p', { class: 'caption' }, 'Click a class to draw its 2-complex.'));
    const totals = w.moment ? `Σ over pairings of |G|^F = ${w.sum}; the Gaussian moment ⟨Vⁿ⟩ by brute force = ${w.moment}: ${w.sum === w.moment ? 'equal' : 'DIFFERENT'}.` : `Σ over pairings of |G|^F = ${w.sum}.`;
    wickBox.replaceChildren(h('div', { class: 'section' }, h('h2', {}, 'The Wick expansion: Feynman graphs and their 2-complexes'),
      h('p', { class: 'muted' }, 'Each pairing of the 4n field factors of n vertices closes the strands into faces; each face carries a free sum over G, so a pairing weighs |G|^F. Classes are taken up to permuting the vertices and the 4 automorphisms of each.'), nk,
      h('div', { class: 'grid', style: 'grid-template-columns: minmax(0, 1.4fr) minmax(0, 1fr)' },
        h('div', { style: 'overflow-x:auto' }, h('table', { class: 'data' },
          h('thead', {}, h('tr', {}, h('th', {}, '#'), h('th', {}, ''), h('th', {}, 'faces F'), h('th', {}, 'pairings'), h('th', {}, '|Aut Γ|'), h('th', {}, 'H₀, H₁, H₂'), h('th', {}, '|G|^F'), h('th', {}, 'BF (gauge)'))), tbody)),
        drawing),
      h('p', { class: 'caption' }, totals),
      w.gluingFaces ? h('p', { class: 'note' }, `Of the 24 gluings of two vertices factor to factor: ${Object.entries(w.gluingFaces).map(([f, k]) => `${k} with ${f} faces`).join(', ')}. The manuscript's dipole is the six-face one.`) : null));
    if (rows.length) rows[0].click();
  };
  await Promise.all([showField(), showLandau(), showWick()]);
}
