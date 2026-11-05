// Walkers on the district's streets (a model): every walker of the urban-science domain's agent-based model, at the
// clock's hour, moving. The engine runs the day (or a tool's scenario) and keeps, at every quarter hour, a minute of
// the walkers' places every ten seconds; here each walker is moved between them in real time and the minute plays
// round while the clock stands. A walker keeps to the left of its footpath, as Melbourne's walkers do (29 of 41
// footpaths keep left, Street Oceanography's sign test), at a place across it of its own; one standing at an event
// stands anywhere across the street. Each is drawn at about one and a half times a person's size.
//
// Beneath them each segment is a ribbon of its walkable width, coloured by its density (walkers per m², sea foam to
// deep sea at 0.6 m⁻² and above, on a square-root scale); where a tool compares two runs, by the difference instead (sea where fewer
// walk, brick rose where more). At the mouth of a segment that is full, a column of the walkers waiting to enter.
// Events of a simulated record are drawn as pale violet glass, never as observed ones.

import { colourMap } from '/visualization/styles/index.js';
import { ask, DAY_NAMES, ensurePrepared, hhmm, MODE_COLOURS, MODES, on, rgb, specAt, state, STATUS, STATUS_COLOURS } from '../shared.js';

const SEA = colourMap('sea').colour, DIFF = colourMap('sea-brick').at;
const RECORD = [0.6, 0.56, 0.74];   // a simulated record's events: pale violet glass, never an observed event's colour
const LOS = [[0.31, 'A'], [0.43, 'B'], [0.72, 'C'], [1.08, 'D'], [2.0, 'E'], [Infinity, 'F']];   // walkers per m² (Fruin 1971; F above 2, as the model counts it)
const losOf = (d) => LOS.find(([x]) => d <= x)[1];
const hash = (i) => (((i * 2654435761) >>> 0) % 10007) / 10007;

export async function layer(ctx) {
  const { gis, dom } = ctx;
  const { h, fmt } = dom;
  const ids = { ribbons: ctx.layerId(':density'), swarm: ctx.layerId(':walkers'), queues: ctx.layerId(':queues'), record: ctx.layerId(':record') };
  let clock = ctx.clock, alive = true, quarter = null, loading = null, wanted = null, loadedAt = 0, stop = null, poll = null;

  const runLine = h('div', { class: 'ld-quiet' });
  const counts = h('div', { class: 'ld-line', style: 'font-size:0.88rem;margin-top:3px' });
  const status = h('div', { class: 'ld-quiet' }, 'preparing the district\'s walkers: its network, counters and streams…');
  const keyBox = h('div', {});
  const seg = h('div', { class: 'ld-seg', style: 'margin:5px 0 3px' });
  const legend = h('div', { class: 'card hf-legend' }, h('b', {}, 'Walkers (a model)'), runLine, seg, keyBox, counts, status);

  const drawKey = () => {
    seg.replaceChildren(...[['pattern', 'by pattern'], ['status', 'by what they do']].map(([k, name]) =>
      h('button', { class: state.colourBy === k ? 'on' : '', onclick: () => { state.colourBy = k; drawKey(); } }, name)));
    const items = state.colourBy === 'pattern' ? MODES.map((m, i) => [MODE_COLOURS[i], m]) : [[STATUS_COLOURS[1], 'walking'], [STATUS_COLOURS[3], 'waiting at a red signal'],
      [STATUS_COLOURS[2], 'waiting to enter a full street'], [STATUS_COLOURS[4], 'standing at the event']];
    const ramp = state.compare && quarter && quarter.scenario
      ? h('div', { class: 'ld-quiet', style: 'margin-top:4px' }, h('div', {}, `walkers on each street, ${state.compare.label}`),
        h('div', { style: `height:7px;border-radius:2px;margin:3px 0;background:linear-gradient(to right, ${[-1, -0.5, 0, 0.5, 1].map((u, i) => `${hex(DIFF(u))} ${25 * i}%`).join(', ')})` }),
        h('div', { style: 'display:flex;justify-content:space-between' }, h('span', {}, 'fewer'), h('span', {}, 'as many'), h('span', {}, 'more')))
      : h('div', { class: 'ld-quiet', style: 'margin-top:4px' }, h('div', {}, 'the street beneath: walkers per m²'),
        h('div', { style: `height:7px;border-radius:2px;margin:3px 0;background:linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map((u) => `${hex(SEA(0.14 + 0.86 * u))} ${100 * u}%`).join(', ')})` }),
        h('div', { style: 'display:flex;justify-content:space-between' }, h('span', {}, '0'), h('span', {}, '0.15'), h('span', {}, '0.6 and more')));
    keyBox.replaceChildren(...[h('div', { class: 'ld-swatches' }, items.map(([c, n]) => h('span', {}, h('i', { class: 'ld-sw dot', style: `background:${c}` }), n))), ramp,
      state.record ? h('div', { class: 'ld-swatches', style: 'margin-top:4px' }, h('span', {}, h('i', { class: 'ld-sw block', style: `background:${hex(RECORD)};opacity:0.6` }), 'a simulated record’s crowds and queues (never observed)')) : null].filter(Boolean));
    if (quarter) paintStatic();
  };
  const hex = (c) => `rgb(${c.map((v) => Math.round(255 * v)).join(',')})`;

  // ---- the quarter's frames ----
  // Walkers by id across the frames: [mode, status, x, y, heading, edge] per frame (NaN where absent).
  function unpack(F) {
    const K = F.frames.length, byId = new Map();
    F.frames.forEach((fr, k) => {
      for (let i = 0; i < fr.length; i += 7) {
        const id = fr[i];
        let w = byId.get(id);
        if (!w) { w = { id, mode: fr[i + 2], u: hash(id), s: new Int8Array(K).fill(-1), x: new Float32Array(K).fill(NaN), y: new Float32Array(K).fill(NaN), hd: new Float32Array(K), e: new Int32Array(K) }; byId.set(id, w); }
        w.s[k] = fr[i + 3], w.x[k] = fr[i + 4] / 10, w.y[k] = fr[i + 5] / 10, w.hd[k] = (fr[i + 6] * Math.PI) / 180, w.e[k] = fr[i + 1];
      }
    });
    return { ...F, K, walkers: [...byId.values()] };
  }

  async function load() {
    if (!alive) return;
    const want = specAt(clock), hour = clock.hour, sig = `${JSON.stringify(want.spec)}@${Math.floor(hour * 4)}`;
    if (quarter && quarter.sig === sig) return;
    wanted = { sig, want, hour };
    if (loading) return;
    loading = (async () => {
      while (wanted && alive) {
        const job = wanted;
        wanted = null;
        try {
          if (!quarter || quarter.sig !== job.sig) status.textContent = `fetching the walkers of ${hhmm(Math.floor(job.hour * 4) / 4)}…`;
          const F = await ask(ctx, 'frame', { ...job.want.spec, hour: job.hour });
          if (F.pending) {
            status.textContent = `simulating ${job.want.label}: the walkers of ${hhmm(job.hour)} come once the run reaches it (now ${hhmm(Math.max(job.want.spec.from, F.reached))})…`;
            clearTimeout(poll);
            poll = setTimeout(() => { if (!wanted) { wanted = job; load(); } }, 1500);
            continue;
          }
          if (F.outside) continue;
          const was = quarter && quarter.scenario;
          quarter = { ...unpack(F), sig: job.sig, label: job.want.label, scenario: job.want.scenario };
          loadedAt = performance.now();
          if (was !== quarter.scenario) drawKey();
          runLine.textContent = `${job.want.scenario ? job.want.label : `the ${DAY_NAMES[clock.weekday]} as built`}, ${hhmm(F.hour)}–${hhmm(F.hour + 0.25)}`;
          status.textContent = F.done ? '' : `the run has reached ${hhmm(F.reached)}`;
          paintStatic();
        } catch (e) {
          status.textContent = String(e.message || e);
        }
      }
      loading = null;
    })();
  }

  // The ribbons' colours, the queues' columns and the counts: once a quarter.
  function paintStatic() {
    const P = state.prep;
    if (!P || !quarter) return;
    const D = quarter.density;
    const cmp = quarter.scenario ? state.compare : null;   // a comparison is drawn with the run it compares, inside its hours
    gis.restyle(ids.ribbons, {
      color: (s, k) => {
        if (cmp) { const d = cmp.b[k] - cmp.a[k], m = Math.max(cmp.a[k], cmp.b[k]); return m < 0.5 ? [0.86, 0.85, 0.83] : DIFF(Math.max(-1, Math.min(1, d / Math.max(4, 0.6 * m)))); }
        return SEA(0.14 + 0.86 * Math.min(1, Math.sqrt((D[k] || 0) / 0.6)));
      },
    });
    const qs = quarter.queues.filter(([, n]) => n >= 3).map(([state2, n]) => {
      const g = P.segments[state2 >> 1], line = state2 % 2 ? [...g.line].reverse() : g.line;
      const a = line[0], b = line[1], L = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
      return { at: [a[0] + (3 * (b[0] - a[0])) / L, a[1] + (3 * (b[1] - a[1])) / L], n, name: g.name };
    });
    gis.addPoints(ids.queues, qs, { position: (f) => f.at, radius: () => 1.6, height: (f) => 2 + 0.6 * f.n, color: () => rgb(STATUS_COLOURS[2]), opacity: 0.9,
      pick: (f) => `${f.n} waiting to enter ${f.name || 'a street'}: it is full` });
    const k0 = quarter.walkers.filter((w) => w.s[0] >= 0);
    const by = (st) => k0.filter((w) => w.s[0] === st).length;
    const crowd = k0.filter((w) => w.mode === 5).length;
    counts.textContent = `${k0.length.toLocaleString('en')} walkers in the district · ${by(3).toLocaleString('en')} at red · ${by(2).toLocaleString('en')} queued${crowd ? ` · ${crowd.toLocaleString('en')} in the event's crowd` : ''}`;
    paintRecord();
  }

  function paintRecord() {
    const R = state.record;
    if (!R) { gis.remove(ids.record); return; }
    const evs = R.events.filter((e) => e.kind === 'crowd' || e.kind === 'queue');
    gis.addPoints(ids.record, evs, { position: (e) => e.at, radius: (e) => (e.kind === 'crowd' ? 4.5 : 2.6), height: (e) => (e.kind === 'crowd' ? 9 : 7), base: () => 0.2, color: () => RECORD,
      opacity: 0.32, pick: (e) => `simulated, not observed: a ${e.kind} ${e.kind === 'crowd' ? 'on' : 'at'} ${e.where}, ${hhmm(e.from)}–${hhmm(e.to)}\n${e.source}` });
  }

  // ---- every frame: the walkers between their places ----
  const cap = 16000, xyz = new Float32Array(3 * cap), col = new Float32Array(3 * cap);
  // By night the walkers are lifted a fifth of the way toward white, so that they stand out on the dark streets.
  const lift = ctx.theme === 'night' ? (c) => c.map((v) => v + (1 - v) * 0.2) : (c) => c;
  const modeRgb = MODE_COLOURS.map((x) => lift(rgb(x))), statusRgb = STATUS_COLOURS.map((x) => lift(rgb(x)));
  function frame(now) {
    if (!quarter || !state.prep) return;
    const segs = state.prep.segments, K = quarter.K, span = quarter.step * (K - 1);
    const t = (((now - loadedAt) / 1000) % span + span) % span, k = Math.min(K - 2, Math.floor(t / quarter.step)), u = t / quarter.step - k;
    let n = 0;
    for (const w of quarter.walkers) {
      if (n >= cap) break;
      const a = w.s[k] >= 0, b = w.s[k + 1] >= 0;
      if (!a && !b) continue;
      const j = a ? k : k + 1, st = w.s[j];
      let x = w.x[j], y = w.y[j];
      if (a && b) x += u * (w.x[k + 1] - x), y += u * (w.y[k + 1] - y);
      else if (a ? u > 0.85 : u < 0.15) continue;   // born or arrived between the frames: shown near its own frame
      const width = (segs[w.e[j]] || {}).width || 4, hd = w.hd[j];
      const off = st === 4 ? (w.u - 0.5) * 0.9 * width : (0.12 + 0.75 * w.u) * 0.5 * width;
      xyz[3 * n] = x - Math.sin(hd) * off, xyz[3 * n + 1] = y + Math.cos(hd) * off, xyz[3 * n + 2] = 0.14;
      const c = state.colourBy === 'pattern' ? modeRgb[w.mode] : statusRgb[st];
      col[3 * n] = c[0], col[3 * n + 1] = c[1], col[3 * n + 2] = c[2];
      n++;
    }
    gis.setSwarm(ids.swarm, n, xyz, col);
  }

  const off = on((what) => {
    if (what === 'scenario') { quarter = null; load(); }
    if (what === 'compare' || what === 'colour') { drawKey(); paintStatic(); }
    if (what === 'record') { paintRecord(); drawKey(); }
  });

  (async () => {
    try {
      const P = await ensurePrepared(ctx);
      if (!alive) return;
      gis.addRibbons(ids.ribbons, P.segments, {
        line: (s) => s.line, width: (s) => Math.max(2, s.width), color: () => SEA(0.1), z: 0.1, opacity: ctx.theme === 'night' ? 0.5 : 0.9, order: 2,
        pick: (s) => {
          const k = P.segments.indexOf(s), d = quarter ? quarter.density[k] || 0 : 0;
          return `${s.name || 'a way'}: ${fmt(d * s.width * s.length, 2)} walkers on ${Math.round(s.length)} m, ${fmt(d, 2)} per m² (level ${losOf(d)})`;
        },
      });
      gis.addSwarm(ids.swarm, { capacity: cap, radius: 0.7, height: 2.6 });
      // By night the walkers are drawn unlit, small lanterns on the dark streets.
      if (ctx.theme === 'night' && ctx.THREE) { const m = gis.layers.get(ids.swarm).object; m.material.dispose(); m.material = new ctx.THREE.MeshBasicMaterial(); }
      stop = gis.onFrame(frame);
      status.textContent = '';
      drawKey();
      load();
    } catch (e) {
      status.textContent = String(e.message || e);
    }
  })();

  return {
    legend,
    setTime(c) { clock = c; load(); },
    describe(hit) {
      if (hit.layer !== ids.ribbons || !state.prep) return null;
      const P = state.prep, k = hit.index, s = P.segments[k], d = quarter ? quarter.density[k] || 0 : 0;
      const cmp = state.compare;
      const sel = {
        kind: 'segment', edge: k, name: s.name, title: s.name || 'A way',
        lines: [`${Math.round(s.length)} m, ${fmt(s.width, 2)} m walkable (${s.kind})`, `at ${hhmm(clock.hour)}: ${fmt(d * s.width * s.length, 2)} walkers, ${fmt(d, 2)} per m², Fruin's level ${losOf(d)}`,
          cmp ? `${cmp.label}: ${fmt(cmp.a[k], 3)} → ${fmt(cmp.b[k], 3)} walkers on it, on average` : null].filter(Boolean),
      };
      sel.actions = [{ title: 'close it or make it one-way', run: () => { state.picked = { edge: k, name: s.name }; ctx.openTool('regime', sel); } }];
      return sel;
    },
    dispose() {
      alive = false;
      off();
      clearTimeout(poll);
      if (stop) stop();
      for (const id of Object.values(ids)) gis.remove(id);
    },
  };
}
