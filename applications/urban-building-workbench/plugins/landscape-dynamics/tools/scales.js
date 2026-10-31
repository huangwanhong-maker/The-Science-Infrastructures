// Scales and memory: the district's own frequencies and flows, computed from its cells. How fast differences of
// crowding between its parts even out, and along which division they last longest; the flows that go round its
// holes, which nothing at a junction can start or stop; what a description of the streets that leaves the shops
// out must remember — those who went in come out later, where they went in; and that the walkers counted one by
// one and the flow drawn on the map are one description. The map shows each in turn.

import { dir, ensureDistrict, mid, request, state } from '../shared.js';
import { brief, figures, fold, hint, swatches, tabs } from '../ui.js';

// How long a passer-by stays in a holding of each use, in hours (declared); uses no passer-by walks into are left out.
const STAY = {
  'food-and-drink': 0.6, 'other-retail': 0.25, 'clothing-and-fashion': 0.33, 'personal-services': 0.75, finance: 0.25, health: 0.5, jewellery: 0.25,
  'entertainment-and-culture': 1.5, 'community-and-public': 0.5, other: 0.33,
};

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars } = ctx.chart;
  const { gis, G } = ctx;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Scales and memory'),
    h('span', { class: 'lede' }, 'How fast the district evens out, what goes round its holes, what its streets remember.')));
  const head = h('div', {});
  root.append(head);
  busy(head, 'the district…');
  try { await ensureDistrict(ctx); } catch (e) { return failure(head, e); }
  const D = state.district;
  const shops = (state.candidates || []).filter((c) => STAY[c.use]).map((c) => ({ at: c.xy, stay: STAY[c.use], use: c.use, name: c.name }));
  let from = null, r = null, view = 'slowest';
  const night = ctx.theme === 'night';
  const ids = { edges: 'landscape-dynamics/tool-scales:edges', junctions: 'landscape-dynamics/tool-scales:junctions', arrows: 'landscape-dynamics/tool-scales:arrows', shops: 'landscape-dynamics/tool-scales:shops' };
  const clear = () => { for (const id of Object.values(ids)) gis.remove(id); };
  let alive = true;
  ctx.onDispose(() => { alive = false; clear(); });

  const out = h('div', {});
  head.replaceChildren(out);
  const setView = (v) => { view = v; paint(); };

  // ---- the map ----
  const street = night ? [0.42, 0.5, 0.58] : [0.7, 0.76, 0.8];
  function paint() {
    if (!r || !alive) return;
    clear();
    const cells = D.cells;
    gis.addRibbons(ids.edges, cells.edges, { line: (f) => f.line, width: () => 2.4, color: () => street, z: 0.35, order: 3, pick: (f) => `${f.name || 'a street'}: ${Math.round(f.length)} m` });
    if (view === 'slowest') {
      const mode = r.slowest.mode, span = Math.max(...Object.values(mode).map(Math.abs)) || 1;
      const div = G.diverging(span, 'sea-brick');
      gis.addPoints(ids.junctions, cells.vertices, { position: (v) => v.at, radius: () => 4.2, height: (v) => 4 + 30 * Math.abs(mode[v.id] || 0) / span, color: (v) => div(mode[v.id] || 0), shadow: true,
        pick: (v) => `${mode[v.id] >= 0 ? 'one side' : 'the other side'} of the slowest division` });
    } else if (view === 'harmonic') {
      const k = Math.min(harmonicPick, r.harmonic.length - 1);
      if (k >= 0) {
        const flow = r.harmonic[k].flow;
        const edges = cells.edges.filter((e) => Math.abs(flow[e.id] || 0) > 0.02);
        gis.addArrows(ids.arrows, edges, { position: (e) => mid(e.line), direction: (e) => dir(e.line, flow[e.id]), length: (e) => 6 + 16 * Math.abs(flow[e.id]), width: () => 4,
          color: () => (night ? [0.89, 0.67, 0.54] : [0.79, 0.55, 0.42]), z: 0.7, opacity: 0.95 });
        gis.restyle(ids.edges, { color: (e) => (Math.abs(flow[e.id] || 0) > 0.02 ? (night ? [0.62, 0.74, 0.84] : [0.36, 0.56, 0.73]) : street) });
      }
    } else if (view === 'memory') {
      const m = r.memory, J = new Map(m.junctions.map((j) => [j.id, j]));
      const top = Math.max(...m.junctions.map((j) => j.pull), 1e-9);
      const seq = G.sequential(0, Math.max(...m.junctions.map((j) => j.delay), 0.1), 'gold', { from: 0.1, to: 1 });
      gis.addPoints(ids.junctions, cells.vertices, { position: (v) => v.at, radius: (v) => (v.id === m.from ? 5.5 : 3.6), height: (v) => 3 + 40 * ((J.get(v.id) || { pull: 0 }).pull / top),
        color: (v) => (v.id === m.from ? (night ? [0.55, 0.72, 0.86] : [0.25, 0.48, 0.66]) : (J.get(v.id) || {}).pull ? seq(J.get(v.id).delay) : street), shadow: true,
        pick: (v) => { const j = J.get(v.id); return j && j.shops ? `${j.shops} shop${j.shops === 1 ? '' : 's'}, a stay of ${fmt(60 * j.delay, 3)} minutes on average${v.id === m.from ? ' — the walkers set down here' : ''}` : 'no shop here'; } });
      gis.addPoints(ids.shops, shops, { position: (s) => s.at, radius: () => 1.3, height: () => 2.5, color: () => (night ? [0.86, 0.72, 0.5] : [0.66, 0.52, 0.3]), pick: (s) => `${s.name} — ${s.use.replace(/-/g, ' ')}` });
    }
  }
  let harmonicPick = 0;
  // A junction picked on the map, in the memory's view, is where the walkers are set down.
  gis.onPick((hit) => {
    if (!alive || view !== 'memory' || !hit || hit.layer !== ids.junctions) return;
    from = hit.feature.id;
    run();
  });

  // ---- the tabs, each drawing its view on the map ----
  const modesPanel = () => {
    const s0 = r.spectra[0], rates = s0.values.slice(s0.harmonic, s0.harmonic + 12);
    return h('div', {},
      figures(h, [{ value: `${fmt(60 / r.slowest.rate, 2)} min`, label: 'the slowest to even out' }, { value: rates[1] ? `${fmt(60 / rates[1], 2)} min` : '—', label: 'the next' }]),
      hint(h, 'On the map: the district\'s slowest division, one side sea, the other brick.'),
      bars(rates.map((v, k) => ({ x: k, y: 60 / v, color: k ? css('--series-3') : css('--accent-deep'), label: `mode ${k + 1}: ${fmt(60 / v, 3)} minutes` })), { labels: rates.map((_, k) => `${k + 1}`), ylabel: 'minutes, the slowest modes first', height: 150 }),
      brief(ctx, { title: 'How fast unevenness of crowding fades', kind: 'proved',
        statement: 'A walker hopping junction to junction: any unevenness is a sum of the district\'s modes, each fading at its own rate. The slowest splits the district along its weakest joins.',
        method: `The Hodge Laplacian on the junctions, each street weighted by its walking rate (1.3 m/s over its length); its eigenvalues and slowest eigenvector (physis).${r.slowest.pieces > 1 ? ` The largest of the ${r.slowest.pieces} pieces the box cuts the network into.` : ''}` }));
  };
  const holesPanel = () => {
    const pickHole = h('select', {}, r.harmonic.map((_, k) => h('option', { value: k, selected: k === harmonicPick || null }, `round hole ${k + 1} of ${r.harmonic.length}`)));
    pickHole.addEventListener('change', () => { harmonicPick = +pickHole.value; paint(); });
    return h('div', {},
      figures(h, r.spectra.map((s) => ({ value: s.harmonic, label: `${['junctions', 'streets', 'blocks'][s.degree] || 'cells'}: harmonic = b${['₀', '₁', '₂'][s.degree] || ''}` }))),
      r.harmonic.length ? h('div', { class: 'ld-row' }, h('span', { class: 'ld-hint' }, 'On the map:'), pickHole) : null,
      brief(ctx, { title: r.harmonic.length ? `${r.harmonic.length} flows only a loop can drive` : 'No flow goes round without a cause', kind: 'proved',
        statement: r.harmonic.length ? 'Round each hole runs a flow that gathers walkers at no junction and circulates round no block: no attraction can make it, nothing at one junction undo it.' : 'Every block is a cell: there is no hole.',
        method: 'The kernel of the Hodge Laplacian on the streets, against the Betti numbers of the complex (Hodge\'s theorem; physis, morphe).' }));
  };
  const memoryPanel = () => {
    const m = r.memory, tauMin = m.tau.map((t) => 60 * t);
    const fromJ = m.junctions.find((j) => j.id === m.from) || { shops: 0, pull: 0, delay: 0 };
    return h('div', {},
      figures(h, [{ value: shops.length, label: 'shops' }, { value: fmt(m.inside, 3), label: `of 100 inside at most, after ${fmt(60 * m.insideAt, 2)} min` }, { value: `${fmt(60 * fromJ.delay, 3)} min`, label: 'a stay, where they start' }]),
      hint(h, 'A hundred walkers set down at the blue junction; pick another on the map.'),
      lines([{ name: 'on the streets, exactly', x: tauMin, y: m.onStreet, color: css('--accent-deep'), width: 2 }, { name: 'with the shops\' stays forgotten', x: tauMin, y: m.onStreetMarkov, color: css('--minor'), dash: '5 4' }],
        { xlabel: 'minutes', ylabel: 'walkers on the streets, of 100', height: 150 }),
      swatches(h, [{ color: css('--accent-deep'), label: 'exact', kind: 'line' }, { color: css('--minor'), label: 'stays forgotten', kind: 'line' }]),
      brief(ctx, { title: m.compartments ? 'The streets remember the shops' : 'No shop to remember', kind: 'proved',
        statement: m.compartments ? 'Those who go in come out later, where they went in: the streets settle in minutes, then go on changing for an hour. Forget the stays and only the minutes are left.' : 'Street Oceanography\'s holdings are not at hand here.',
        method: 'Junctions and shops as one master equation in detailed balance, by Crank–Nicolson; the memory kernel B exp(Dτ) C from the exact elimination of the shops (physis).',
        limits: 'Stays and the rate of entering declared by use; walkers stay in the district.' }),
      fromJ.shops ? fold(h, 'the memory at that junction', lines([{ name: 'returning', x: tauMin, y: m.kernel, color: css('--minor-deep'), width: 2 }], { xlabel: 'minutes since they went in', ylabel: 'returning, per hour', height: 120 })) : null);
  };
  const particlesPanel = () => h('div', {},
    figures(h, [{ value: fmt(r.operators.deviation, 2), label: 'operator against master equation' }, { value: r.operators.states, label: 'counts reachable' }]),
    brief(ctx, { title: 'One by one, or as a flow: one description', kind: 'proved',
      statement: `${r.operators.walkers} walkers among three junctions as a reaction network: its creation–annihilation operator is the master equation exactly, and the mean count follows the flow with no correction.`,
      method: 'The Doi–Peliti operator of the network on its Fock space, against the generator of the master equation (physis chemistry).' }));

  const T = tabs(h, [{ name: 'Modes', panel: () => modesPanel() }, { name: 'Holes', panel: () => holesPanel() }, { name: 'Memory', panel: () => memoryPanel() }, { name: 'Particles', panel: () => particlesPanel() }],
    { remember: 'landscape-scales', lazy: true, onShow: (n) => setView({ Modes: 'slowest', Holes: 'harmonic', Memory: 'memory', Particles: 'none' }[n]) });
  let current = null;
  async function run() {
    const on = out.querySelector('.ld-tab.on');
    current = on ? on.dataset.name : current;
    if (!r) out.replaceChildren(hint(h, 'the district\'s Laplacians, its harmonic flows, the shops eliminated exactly…'));
    try { r = await ctx.ask('scales', request({ shops: shops.map((s) => ({ at: s.at, stay: s.stay })), from })); } catch (e) { return failure(out, e); }
    out.replaceChildren(T.el);
    T.show(current);
  }
  run();
}
