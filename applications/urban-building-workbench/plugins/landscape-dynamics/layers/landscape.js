// The landscape running (a model): a day of the district as one process, drawn at the clock's hour. Walkers stand
// at the junctions as columns and walk along the streets each way, as many as cross it, an arrow where more walk
// one way than the other; the blocks are tinted by their warmth; the two stores stand tall while open; the lane
// appears when it is built — the cells of the stratum change — and turns deep when it is recognised. The legend keeps the hour's state and the
// events so far. The run is the engine's (urban-science, landscape/process.hpp), forced by the clock's date.

import { dir, ensureDistrict, frameAt, mid, on, runDay, state } from '../shared.js';

const SPEED = 1.3;          // m/s, walking
const PER_DOT = 25;         // persons an hour crossing a street one way for each walker drawn

export async function layer(ctx) {
  const { dom, gis, G } = ctx;
  const { h, fmt } = dom;
  const ids = { blocks: ctx.layerId(':blocks'), edges: ctx.layerId(':edges'), crowd: ctx.layerId(':crowd'), arrows: ctx.layerId(':arrows'), stores: ctx.layerId(':stores'), walkers: ctx.layerId(':walkers') };
  const status = h('div', { class: 'faint' }, 'running the day…');
  const now = h('div', {});
  const recent = h('div', { class: 'faint' });
  const legend = h('div', { class: 'card' }, h('b', {}, 'The landscape running (a model)'), ctx.chart.ramp('gold', 14, 34, 'block warmth, °C'), now, recent, status);
  state.active.add('landscape');
  const warmth = G.sequential(14, 34, 'gold', { from: 0.05, to: 1 });
  let drawnStratum = -1, frame = null, stop = null, clock = ctx.clock;
  const night = ctx.theme === 'night';
  const lineColour = night ? [0.6, 0.7, 0.8] : [0.36, 0.56, 0.73], laneColour = [0.79, 0.55, 0.42], deep = [0.6, 0.33, 0.24];

  const stratumCells = () => state.run.strata[frame.stratum].cells;
  function build() {
    const cells = stratumCells();
    gis.addFlats(ids.blocks, cells.blocks, { ring: (f) => f.ring, color: (f) => warmth(frame.temperature[f.id]) || [0.8, 0.8, 0.8], z: 0.12,
      pick: (f) => `a block at ${fmt(frame.temperature[f.id], 3)} °C` });
    gis.addRibbons(ids.edges, cells.edges, { line: (f) => f.line, width: (f) => (f.type === 'lane' ? 5 : 2.4), z: 0.3, order: 3,
      color: (f) => (f.id === frame.lane ? (frame.recognised ? deep : laneColour) : lineColour), pick: (f) => `${f.name || 'a segment'}: ${fmt((frame.forth[f.id] || 0) + (frame.back[f.id] || 0), 3)} walkers an hour, ${fmt(Math.abs(frame.flow[f.id] || 0), 2)} more one way than the other` });
    gis.addPoints(ids.crowd, cells.vertices, { position: (f) => f.at, radius: () => 2.2, height: (f) => Math.max(0.3, (frame.walkers[f.id] || 0) / 4), color: () => (night ? [0.62, 0.78, 0.86] : [0.25, 0.43, 0.6]), opacity: 0.9,
      pick: (f) => `${fmt(frame.walkers[f.id] || 0, 3)} walkers at this junction` });
    gis.addArrows(ids.arrows, cells.edges.filter((e) => e.length > 25), { position: (f) => mid(f.line), direction: (f) => dir(f.line, frame.flow[f.id] || 0), length: (f) => (Math.abs(frame.flow[f.id] || 0) < 2 ? 0.01 : Math.min(18, 5 + Math.abs(frame.flow[f.id]) / 3)),
      width: (f) => (Math.abs(frame.flow[f.id] || 0) < 2 ? 0.01 : 3.4), color: () => laneColour, z: 0.6, opacity: 0.85 });
    gis.addPoints(ids.stores, state.run.strata[frame.stratum].stores, { position: (s) => s.at, radius: () => 3.4, height: (s, k) => ((k === 0 ? frame.a : frame.b) ? 16 : 3), color: (s, k) => ((k === 0 ? frame.a : frame.b) ? laneColour : [0.7, 0.7, 0.7]), shadow: true,
      pick: (s, k) => `${s.name}: ${(k === 0 ? frame.a : frame.b) ? 'open' : 'shut'}` });
    drawnStratum = frame.stratum;
  }
  function restyle() {
    gis.restyle(ids.blocks, {});
    gis.restyle(ids.edges, {});
    gis.restyle(ids.crowd, {});
    gis.restyle(ids.arrows, {});
    gis.restyle(ids.stores, {});
  }
  // Walkers along the streets, each way as many as cross it that way, keeping to the left of their way, placed
  // afresh every frame of the animation.
  const xyz = new Float32Array(3 * 6000), rgb = new Float32Array(3 * 6000);
  function walk(ms) {
    if (!frame) return;
    const cells = stratumCells();
    let n = 0;
    const t = ms / 1000;
    for (const e of cells.edges) {
      const L = e.length || 1, a = e.line[0], b = e.line[e.line.length - 1], d = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
      const ux = (b[0] - a[0]) / d, uy = (b[1] - a[1]) / d;
      const c = e.id === frame.lane ? laneColour : night ? [0.85, 0.9, 0.95] : [0.18, 0.3, 0.42];
      for (const [way, count] of [[1, frame.forth[e.id] || 0], [-1, frame.back[e.id] || 0]]) {
        const k = Math.min(24, Math.floor(count / PER_DOT));
        for (let i = 0; i < k && n < 6000; i++) {
          let s = ((t * SPEED + (i / k) * L) % L + L) % L;
          if (way < 0) s = L - s;
          const p = G.pointAlong(e.line, s);
          // The left of the direction of travel: (−u_y, u_x) going forth, the opposite coming back.
          xyz[3 * n] = p[0] - way * uy * 1.1, xyz[3 * n + 1] = p[1] + way * ux * 1.1, xyz[3 * n + 2] = 0.3;
          rgb[3 * n] = c[0], rgb[3 * n + 1] = c[1], rgb[3 * n + 2] = c[2];
          n++;
        }
      }
    }
    gis.setSwarm(ids.walkers, n, xyz, rgb);
  }
  function show() {
    if (!state.run) return;
    frame = frameAt(state.run, clock.hour);
    if (!frame) return;
    if (frame.stratum !== drawnStratum || !gis.layers.get(ids.blocks)) build();
    else restyle();
    const S = state.run.strata[frame.stratum];
    const names = S.stores.map((s) => s.name);
    now.replaceChildren(
      h('div', {}, `${String(Math.floor(clock.hour)).padStart(2, '0')}:${String(Math.round((clock.hour % 1) * 60)).padStart(2, '0')} — ${names[0]} ${frame.a ? 'open' : 'shut'}, ${names[1]} ${frame.b ? 'open' : 'shut'}`),
      h('div', {}, frame.lane ? `the lane is there${frame.recognised ? ', recognised (it cannot be removed)' : ''}` : frame.permit ? 'a permit for the lane, not yet built' : 'no permit for a lane'),
      h('div', { class: 'faint' }, `${S.counts[0]} junctions, ${S.counts[1]} segments, ${S.counts[2]} blocks; b₁ = ${S.betti[1]}`));
    const past = state.run.events.filter((e) => e.time <= clock.hour).slice(-4).reverse();
    recent.replaceChildren(...past.map((e) => h('div', {}, `${hhmm(e.time)} ${e.mechanism}${e.match && !e.mechanism.includes('store') ? '' : `: ${e.match}`}`)));
    status.textContent = `${state.run.events.length} events through ${state.run.date}; walkers kept to ${fmt(Math.abs(state.run.balance.residual), 2)}; ${state.run.entriesSource}`;
  }
  const off = on((what) => { if (what === 'run') { drawnStratum = -1; show(); } });
  try {
    await ensureDistrict(ctx);
    if (!state.run) await runDay(ctx);
    gis.addSwarm(ids.walkers, { capacity: 6000, radius: 0.5, height: 1.7 });
    stop = gis.onFrame(walk);
    show();
  } catch (e) {
    status.textContent = String(e.message || e);
  }
  return {
    legend,
    setTime(c) { clock = c; show(); },
    dispose() { off(); if (stop) stop(); state.active.delete('landscape'); },
  };
}

const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.floor((x % 1) * 60)).padStart(2, '0')}`;
