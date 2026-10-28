// Symmetry selection (§7, App. I.4): Frank's mirror-symmetric rule family, deterministic and
// stochastic. Mass action from a slight excess: the excess m and its equation; from an exact racemate
// it stays racemic (Lemma 4.1). The Markov chain from a racemic seed: a law of outcomes exactly
// symmetric, every outcome asymmetric (Prop. 4.2 and its remark). And the cubic normal form.

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Selection: a symmetric law, asymmetric outcomes'),
    h('span', { class: 'lede' }, 'A + L → 2L, A + R → 2R, L + R → D: the mirror L ↔ R carries the rule family onto itself. Physis\'s reaction kinetics integrate it as mass action and solve its Markov chain exactly.')));
  const P = { k: 1, ki: 1, k0: 0, l0: 0.011, r0: 0.01, t: 30, precursor: 30 };
  const sliders = h('div', { class: 'panel', style: 'align-self: start' });
  const out = h('div', {});
  root.append(h('div', { class: 'grid', style: 'grid-template-columns: 300px minmax(0, 1fr)' }, sliders, out));
  const control = (key, name, min, max, step) => {
    const input = h('input', { type: 'range', min, max, step, value: P[key] });
    const o = h('output', {}, fmt(P[key]));
    input.addEventListener('input', () => { o.textContent = fmt(+input.value); });
    input.addEventListener('change', () => { P[key] = +input.value; run(); });
    return h('label', { class: 'slider' }, name, input, o);
  };
  sliders.append(h('h3', {}, 'rates and start'),
    control('k', 'k (autocat.)', 0.1, 3, 0.1), control('ki', 'kᵢ (inhibition)', 0, 3, 0.1), control('k0', 'k₀ (uncatalysed)', 0, 0.2, 0.01),
    control('l0', 'ℓ₀', 0, 0.05, 0.001), control('r0', 'r₀', 0, 0.05, 0.001), control('t', 'time', 5, 200, 5),
    h('h3', { style: 'margin-top:10px' }, 'Markov chain'), control('precursor', 'A molecules', 4, 60, 1),
    h('p', { class: 'caption' }, 'The chain starts from one L, one R and the A molecules; the ODE from concentrations a₀ = 1, ℓ₀, r₀.'));
  const run = async () => {
    busy(out);
    let f, pf;
    try {
      [f, pf] = await Promise.all([ctx.ask('frank', P), ctx.ask('pitchfork', { b: 1 })]);
    } catch (error) {
      return failure(out, error);
    }
    const c1 = css('--series-1'), c2 = css('--series-2'), c3 = css('--series-3'), c4 = css('--series-4');
    const o = f.ode;
    const traj = lines([{ name: 'a (precursor)', x: o.t, y: o.a, color: c3 }, { name: 'ℓ', x: o.t, y: o.l, color: c1 }, { name: 'r', x: o.t, y: o.r, color: c2 }, { name: 'd', x: o.t, y: o.d, color: c4, dash: '5 4' }],
      { xlabel: 't', ylabel: 'concentration', height: 230 });
    const mchart = lines([{ name: 'm = (ℓ − r)/(ℓ + r)', x: o.t, y: o.m, color: c1 }], { xlabel: 't', ylabel: 'm', ylim: [-1, 1], height: 200 });
    const ch = f.chain;
    const hist = bars(ch.excess.map(([k, p]) => ({ x: k, y: p, color: k > 0 ? c1 : c2, label: `n_L − n_R = ${k}: ${fmt(p, 4)}` })),
      { xlabel: 'n_L − n_R when the chain stops', ylabel: 'probability', height: 230 });
    const paths = lines(ch.paths.map((p, i) => ({ x: p.t, y: p.m, color: p.m[p.m.length - 1] > 0 ? c1 : c2, width: 1.2, opacity: 0.75, name: `path ${i + 1}` })),
      { xlabel: 't', ylabel: 'm on sample paths', ylim: [-1, 1], height: 220 });
    const pfData = [];
    for (const b of pf.branches) for (const e of b.equilibria) pfData.push({ a: b.a, m: e.m, stable: e.stable });
    const stableUp = pfData.filter((d) => d.stable && d.m > 0), stableDown = pfData.filter((d) => d.stable && d.m < 0), stableZero = pfData.filter((d) => d.stable && d.m === 0), unstable = pfData.filter((d) => !d.stable && d.m === 0);
    const pitch = lines([
      { x: stableZero.map((d) => d.a), y: stableZero.map((d) => d.m), color: c1, name: 'stable' },
      { x: unstable.map((d) => d.a), y: unstable.map((d) => d.m), color: c3, dash: '4 4', name: 'unstable' },
      { x: stableUp.map((d) => d.a), y: stableUp.map((d) => d.m), color: c1, name: 'stable' },
      { x: stableDown.map((d) => d.a), y: stableDown.map((d) => d.m), color: c1, name: 'stable' }], { xlabel: 'a', ylabel: 'm*', height: 200 });
    out.replaceChildren(
      h('div', { class: 'grid cols-2' },
        h('div', {}, h('h3', {}, 'Mass action'), traj, legend([{ name: 'a', color: c3 }, { name: 'ℓ', color: c1 }, { name: 'r', color: c2 }, { name: 'd', color: c4 }])),
        h('div', {}, h('h3', {}, 'The excess m'), mchart,
          h('p', { class: 'caption' }, `ṁ = (kᵢ s/2) m (1 − m²) − (2k₀ a/s) m holds along the solution to ${fmt(o.identityError, 2)} (relative); a + ℓ + r + 2d is kept to ${fmt(o.conservationError, 2)}. m(t_end) = ${fmt(o.mFinal, 5)}.`))),
      h('p', { class: 'note' }, P.l0 === P.r0 ? 'With ℓ₀ = r₀ the solution stays exactly racemic: an equivariant flow keeps the stabilizer of its start (Lemma 4.1, which holds with equality for a flow).' :
        'δ = ℓ − r obeys δ̇ = k a δ: its sign never changes, so a deterministic run selects nothing; it amplifies the sign it was given. Without k₀, m tends to ±1 as the minor enantiomer is consumed.'),
      h('div', { class: 'grid cols-2 section' },
        h('div', {}, h('h3', {}, 'The Markov chain: law of outcomes'), hist,
          h('p', { class: 'caption' }, `${ch.states} reachable states; generator commutes with L ↔ R: ${ch.commutes ? 'exactly' : 'no'}. The law is symmetric to ${fmt(ch.asymmetry, 2)}; P(both extinct) = ${fmt(ch.extinct, 4)}. Every outcome but extinction has one enantiomer only.`)),
        h('div', {}, h('h3', {}, 'Sample paths (Gillespie)'), paths, h('p', { class: 'caption' }, 'Each path picks a sign; the ensemble does not.'))),
      h('div', { class: 'section' }, h('h3', {}, 'The cubic normal form ṁ = am − bm³ (§7.3)'), pitch,
        h('p', { class: 'caption' }, 'For a > 0 the origin loses stability to ±√(a/b), where f′ = −2a. The vector field stays odd for every a.')));
  };
  run();
}
