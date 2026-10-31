// Scales and memory: the district's own frequencies and flows, computed from its cells. How fast differences of
// crowding between its parts even out, and along which division they last longest; the flows that go round its
// holes, which nothing at a junction can start or stop; what a description of the streets that leaves the shops
// out must remember — those who went in come out later, where they went in; and that the walkers counted one by
// one and the flow drawn on the map are one description. The map shows each in turn.

import { dir, ensureDistrict, mid, request, state } from '../shared.js';

// How long a passer-by stays in a holding of each use, in hours (declared); uses no passer-by walks into are left out.
const STAY = {
  'food-and-drink': 0.6, 'other-retail': 0.25, 'clothing-and-fashion': 0.33, 'personal-services': 0.75, finance: 0.25, health: 0.5, jewellery: 0.25,
  'entertainment-and-culture': 1.5, 'community-and-public': 0.5, other: 0.33,
};

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, legend } = ctx.chart;
  const { gis, G } = ctx;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Scales and memory'),
    h('span', { class: 'lede' }, 'The district\'s own frequencies, its flows round the holes, and what the streets remember of the shops along them.')));
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

  const mapKeys = h('div', { class: 'keys' });
  const setView = (v) => { view = v; for (const b of mapKeys.children) b.classList.toggle('on', b.dataset.v === v); paint(); };
  for (const [v, name] of [['slowest', 'on the map: the slowest division'], ['harmonic', 'the flows round the holes'], ['memory', 'the shops\' memory']])
    mapKeys.append(h('button', { class: `key${v === view ? ' on' : ''}`, 'data-v': v, onclick: () => setView(v) }, name));
  const out = h('div', {});
  head.replaceChildren(mapKeys, out);

  // ---- the map ----
  const street = night ? [0.42, 0.5, 0.58] : [0.7, 0.76, 0.8];
  function paint() {
    if (!r || !alive) return;
    clear();
    const cells = D.cells;
    gis.addRibbons(ids.edges, cells.edges, { line: (f) => f.line, width: () => 2.4, color: () => street, z: 0.35, order: 3, pick: (f) => `${f.name || 'a segment'}: ${Math.round(f.length)} m` });
    if (view === 'slowest') {
      const mode = r.slowest.mode, span = Math.max(...Object.values(mode).map(Math.abs)) || 1;
      const div = G.diverging(span, 'sea-brick');
      gis.addPoints(ids.junctions, cells.vertices, { position: (v) => v.at, radius: () => 4.2, height: (v) => 4 + 30 * Math.abs(mode[v.id] || 0) / span, color: (v) => div(mode[v.id] || 0), shadow: true,
        pick: (v) => `junction #${v.id}: ${mode[v.id] >= 0 ? 'one side' : 'the other side'} of the district's slowest division (${fmt(mode[v.id], 3)})` });
    } else if (view === 'harmonic') {
      const k = Math.min(harmonicPick, r.harmonic.length - 1);
      if (k >= 0) {
        const flow = r.harmonic[k].flow;
        const edges = cells.edges.filter((e) => Math.abs(flow[e.id] || 0) > 0.02);
        gis.addArrows(ids.arrows, edges, { position: (e) => mid(e.line), direction: (e) => dir(e.line, flow[e.id]), length: (e) => 6 + 16 * Math.abs(flow[e.id]), width: () => 4,
          color: () => (night ? [0.89, 0.67, 0.54] : [0.79, 0.55, 0.42]), z: 0.7, opacity: 0.95 });
        gis.restyle(ids.edges, { color: (e) => (Math.abs(flow[e.id] || 0) > 0.02 ? (night ? [0.62, 0.74, 0.84] : [0.36, 0.56, 0.73]) : street) });
      }
    } else {
      const m = r.memory, J = new Map(m.junctions.map((j) => [j.id, j]));
      const top = Math.max(...m.junctions.map((j) => j.pull), 1e-9);
      const seq = G.sequential(0, Math.max(...m.junctions.map((j) => j.delay), 0.1), 'gold', { from: 0.1, to: 1 });
      gis.addPoints(ids.junctions, cells.vertices, { position: (v) => v.at, radius: (v) => (v.id === m.from ? 5.5 : 3.6), height: (v) => 3 + 40 * ((J.get(v.id) || { pull: 0 }).pull / top),
        color: (v) => (v.id === m.from ? (night ? [0.55, 0.72, 0.86] : [0.25, 0.48, 0.66]) : (J.get(v.id) || {}).pull ? seq(J.get(v.id).delay) : street), shadow: true,
        pick: (v) => { const j = J.get(v.id); return j && j.shops ? `junction #${v.id}: ${j.shops} shop${j.shops === 1 ? '' : 's'} open on it, drawing ${fmt(j.pull, 3)} passers-by an hour each, keeping them ${fmt(60 * j.delay, 3)} minutes on average${v.id === m.from ? ' — the walkers set down here' : ''}` : `junction #${v.id}: no shop opens on it`; } });
      gis.addPoints(ids.shops, shops, { position: (s) => s.at, radius: () => 1.3, height: () => 2.5, color: () => (night ? [0.86, 0.72, 0.5] : [0.66, 0.52, 0.3]), pick: (s) => `${s.name} — ${s.use.replace(/-/g, ' ')}, a stay of about ${Math.round(60 * s.stay)} minutes (declared)` });
    }
  }
  let harmonicPick = 0;
  // A junction picked on the map, in the memory's view, is where the walkers are set down.
  gis.onPick((hit) => {
    if (!alive || view !== 'memory' || !hit || hit.layer !== ids.junctions) return;
    from = hit.feature.id;
    run();
  });

  async function run() {
    busy(out, 'the district\'s Laplacians, its harmonic flows, and the shops eliminated exactly…');
    try { r = await ctx.ask('scales', request({ shops: shops.map((s) => ({ at: s.at, stay: s.stay })), from })); } catch (e) { return failure(out, e); }
    const s0 = r.spectra[0];
    const rates = s0.values.slice(s0.harmonic, s0.harmonic + 12);
    const m = r.memory;
    const pickHole = h('select', {}, r.harmonic.map((_, k) => h('option', { value: k, selected: k === harmonicPick || null }, `the flow round hole ${k + 1} of ${r.harmonic.length}`)));
    pickHole.addEventListener('change', () => { harmonicPick = +pickHole.value; setView('harmonic'); });
    const harmonicKeys = h('div', { class: 'keys' }, h('label', { class: 'pick' }, 'on the map:', pickHole));
    const tauMin = m.tau.map((t) => 60 * t);
    const fromJ = m.junctions.find((j) => j.id === m.from) || { shops: 0, pull: 0, delay: 0 };
    out.replaceChildren(
      h('div', { class: 'section' }, h('h2', {}, 'How fast the district evens out'),
        bars(rates.map((v, k) => ({ x: k, y: 60 / v, color: k ? css('--series-3') : css('--accent-deep'), label: `mode ${k + 1}: ${fmt(60 / v, 3)} minutes` })), { labels: rates.map((_, k) => `${k + 1}`), ylabel: 'minutes to fall to 1/e, the slowest modes first', height: 170 }),
        ctx.finding({ title: `A lone walker forgets which side of the district it started on in about ${fmt(60 / r.slowest.rate, 2)} minutes`, kind: 'proved',
          statement: `A walker hopping from junction to junction along the streets at walking pace: an unevenness of crowding between the junctions is a sum of the district's modes, each fading at its own rate. The slowest divides the district in two along its weakest joins — drawn on the map, one side sea, the other brick — and fades in ${fmt(60 / r.slowest.rate, 3)} minutes.${r.slowest.pieces > 1 ? ` (The box cuts the network into ${r.slowest.pieces} pieces; this is the largest, with ${r.slowest.largest} of the ${r.slowest.junctions} junctions. The bars count every piece's modes.)` : ''} A district whose two sides are joined by few or long streets keeps its sides apart for longer.`,
          method: 'The Hodge Laplacian on the junctions, each segment weighted by its walking rate (1.3 m/s over its length): the generator of the hop; its eigenvalues and slowest eigenvector (physis).' })),
      h('div', { class: 'section' }, h('h2', {}, 'Flows round the holes'),
        r.harmonic.length ? harmonicKeys : null,
        h('table', { class: 'data' }, h('thead', {}, h('tr', {}, h('th', { class: 'left' }, 'cells'), h('th', {}, 'how many'), h('th', {}, 'Betti number'), h('th', {}, 'harmonic'))),
          h('tbody', {}, r.spectra.map((s) => h('tr', {}, h('td', { class: 'left' }, ['junctions', 'segments', 'blocks'][s.degree] || `${s.degree}-cells`), h('td', {}, String(s.cells)), h('td', {}, String(r.betti[s.degree] ?? 0)), h('td', {}, String(s.harmonic)))))),
        ctx.finding({ title: r.harmonic.length ? `${r.harmonic.length} flow${r.harmonic.length === 1 ? '' : 's'} that only a loop can drive` : 'No flow goes round without a cause at a junction', kind: 'proved',
          statement: r.harmonic.length
            ? `The district has ${r.betti[1]} hole${r.betti[1] === 1 ? '' : 's'} — block${r.betti[1] === 1 ? '' : 's'} with a dead end or a loop of lanes inside, left open. Round each runs a flow that neither gathers nor leaves walkers at any junction and circulates round no block: nothing that pushes walkers toward a place (a store, a station) can make it, nor anything at a single junction undo it; only something that itself goes round, a habit of walking the loop. Drawn on the map as arrows. In every degree the flows of this kind are exactly as many as the holes of that dimension.`
            : 'Every block is a cell: there is no hole, and any flow that conserves walkers at every junction circulates round some block. In every degree the harmonic flows are exactly as many as the holes of that dimension: none here.',
          method: 'The kernel of the Hodge Laplacian on the segments, weighted as above; the flux G φ of each harmonic cochain; against the Betti numbers of the complex (Hodge\'s theorem; physis, morphe).' })),
      h('div', { class: 'section' }, h('h2', {}, 'What the streets remember of the shops'),
        h('p', { class: 'caption' }, `${shops.length} holdings of the district a passer-by walks into (Street Oceanography's holdings), each opening on its nearest junction, drawing ${fmt(m.visit)} passers-by an hour and keeping each for a declared stay by its use. ${m.compartments} groups of shops alike at one junction. A hundred walkers are set down at the junction drawn blue; pick another on the map.`),
        lines([{ name: 'the exact streets', x: tauMin, y: m.onStreet, color: css('--accent-deep'), width: 2 }, { name: 'shops returning walkers at once', x: tauMin, y: m.onStreetMarkov, color: css('--minor'), dash: '5 4' },
          { name: 'at the junction where they were set down', x: tauMin, y: m.here, color: css('--series-2'), width: 1.4 }],
        { xlabel: 'minutes after a hundred walkers are set down', ylabel: 'walkers on the streets', height: 190, ylim: [0, 100] }),
        legend([{ name: 'on the streets, exactly', color: css('--accent-deep') }, { name: 'with the shops\' stays forgotten', color: css('--minor') }, { name: 'at the junction where they were set down', color: css('--series-2') }]),
        ctx.finding({ title: m.compartments ? `Leaving the shops out, the streets must remember: they settle in minutes, then go on for an hour more` : 'No shop in the district to remember', kind: 'proved',
          statement: m.compartments
            ? `Described by its streets alone, the district obeys an exact equation with a memory: walkers who went into a shop come out, where they went in, after its stay, so the streets' future depends on their past. Set down at one junction, a hundred walkers spread over the district within minutes, and the shops take in up to ${fmt(m.inside, 3)} of them at once, ${fmt(60 * m.insideAt, 2)} minutes after; they come back out over the shops' stays, and the streets go on changing for an hour or more as visitors go in and out. Forget the stays, returning everyone at once, and the streets keep all hundred and settle in minutes: one time scale where the district has two, the walk's and the shops'. At the junction where they were set down, ${fromJ.shops} shop${fromJ.shops === 1 ? '' : 's'} keep their visitors ${fmt(60 * fromJ.delay, 3)} minutes on average.`
            : 'Street Oceanography\'s holdings are not at hand, or none lies in the district: the streets have nothing to remember.',
          method: 'Walkers at the junctions and visitors in the shops as one master equation in detailed balance, propagated by Crank–Nicolson (calx\'s θ-method), against the streets alone with the shops returning their visitors at once; the memory at the junction, K(τ) = B exp(Dτ) C, from the exact elimination of its shops (physis).',
          against: `the walkers kept, streets and shops together, ${m.kept === 0 ? 'exactly' : `to ${fmt(Math.abs(m.kept), 2)} of the hundred`}`,
          limits: 'The stays and the rate of entering are declared, not measured; every shop is open; walkers stay in the district.' }),
        fromJ.shops ? lines([{ name: 'the kernel', x: tauMin, y: m.kernel, color: css('--minor-deep'), width: 2 }], { xlabel: 'minutes since they went in', ylabel: 'returning, per hour', height: 130 }) : null,
        fromJ.shops ? h('p', { class: 'caption' }, 'The memory at that junction: the rate at which those who went in some minutes ago come back out. A single kind of shop makes an exponential; shops of different stays, a sum of them.') : null),
      h('div', { class: 'section' }, h('h2', {}, 'Walkers one by one, and as a flow'),
        ctx.finding({ title: `Counted one by one or followed as a flow, the walkers are one description (deviation ${fmt(r.operators.deviation, 2)})`, kind: 'proved',
          statement: `${r.operators.walkers} walkers hopping among three junctions of the district (#${r.operators.junctions.join(', #')}) as a reaction network, each hop a reaction: written with creation and annihilation operators on the space of all counts (${r.operators.fock} states), its operator is the master equation on the ${r.operators.states} counts they can reach, exactly; and because a hop takes one walker at a time, the mean count obeys the flow's own equation with no correction. The flows drawn on the map are the means of walkers counted one by one.`,
          method: 'The Doi–Peliti operator of the network, against the generator of its master equation (physis chemistry).' })));
    paint();
  }
  run();
}
