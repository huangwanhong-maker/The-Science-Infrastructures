// Waves on Swanston Street, a model: the crowd walking north from Flinders Street, second by second, as a
// kinematic wave (Lighthill–Whitham–Richards by Godunov's scheme on physis's network; the domain's
// flow/kinematic.hpp). The street's blocks and cross streets come from the city model; a signal stands at
// each street and none at a lane; walkers enter at the mean rate the sensors at its start count northward,
// at that hour of the clock's weekday. The engine runs the day from 06:00 and keeps four minutes from the
// clock's hour, two seconds apart; the page plays them in a loop in real time.
//   - the water: the density along the street (people per m²), on its footways or, with the carriageway
//     opened, across it;
//   - the walkers: as many as the density holds, moving at the model's speed (the flow over the density),
//     stopping at red, queues spilling back from the crossings as waves;
//   - the signals: a lamp each side of every signalled crossing.
// Nobody turns into a cross street in the model, so its flow further up overstates what the sensors there
// count: the difference is what the cross streets take.

import { sensorBlocks } from './blocks.js';
import { entering, streetModel } from './street-model.js';

const H = 40;            // metres of water per person per m²
const STREET = 'melbourne-swanston';
const FOOTWAY = 10, CARRIAGEWAY = 12;   // m, walkable on an ordinary day, and the carriageway between (declared)

export async function layer(ctx) {
  const { gis, ask, dom, G } = ctx;
  const { h, fmt } = dom;
  const slabId = ctx.layerId(':water'), swarmId = ctx.layerId(':walkers'), lampId = ctx.layerId(':lamps');
  let clock = ctx.clock, model = null, sensors = [], run = null, key = '', opened = false, scale = 1, stop = null, pending = null;
  const status = h('div', { class: 'faint' }, 'building the street from the city model…');
  const openBox = h('input', { type: 'checkbox' });
  const scaleKeys = h('div', { class: 'keys' }, [1, 2, 4].map((k) => h('button', { class: `key${k === 1 ? ' on' : ''}`, onclick: () => { scale = k; for (const b of scaleKeys.children) b.classList.toggle('on', b.textContent === `×${k}`); refresh(true); } }, `×${k}`)));
  const legend = h('div', { class: 'card' }, h('b', {}, 'Waves on Swanston Street (a model)'),
    ctx.chart.ramp('sea', 0, 1, 'people per m²'),
    h('div', { class: 'faint' }, 'Walkers going north from Flinders Street at the rate counted there, each a person, stopping at red; the water is their density.'),
    h('label', { class: 'pick' }, openBox, ' open the carriageway to walkers'), h('div', { class: 'faint' }, 'the crowd, against the counts:'), scaleKeys, status);
  openBox.addEventListener('change', () => { opened = openBox.checked; refresh(true); });

  // ---- the street and its walkers ----
  const walkers = { x: new Float32Array(8000), lat: new Float32Array(8000), n: 0 };
  const xyz = new Float32Array(3 * 8000), rgb = new Float32Array(3 * 8000);
  let tau = 0, last = 0, inflow = 0, stations = [], toAxis = (x) => x;
  const colourOf = G.sequential(0, 1, 'sea', { from: 0.12, to: 1 });
  const lateral = () => {
    const w = run.width[0], half = (CARRIAGEWAY + FOOTWAY) / 2;
    if (w > FOOTWAY + 0.1) return (Math.random() * 2 - 1) * half;
    const side = Math.random() < 0.5 ? -1 : 1;
    return side * (CARRIAGEWAY / 2 + Math.random() * (FOOTWAY / 2));
  };
  function seed() {
    walkers.n = 0;
    const rho = run.density[0], I = run.interfaces;
    for (let c = 0; c < rho.length; c++) {
      const dx = I[c + 1] - I[c], people = rho[c] * run.width[0] * dx * scaleWalk;
      let k = Math.floor(people + Math.random());
      while (k-- > 0 && walkers.n < 8000) { walkers.x[walkers.n] = I[c] + Math.random() * dx; walkers.lat[walkers.n] = lateral(); walkers.n++; }
    }
  }
  let scaleWalk = 1;   // walkers drawn per person (below one when the street holds more than the swarm can)
  const rowAt = (t) => { const k = t / (run.times[1] - run.times[0]); const i = Math.min(run.times.length - 2, Math.floor(k)); return [i, k - i]; };
  const speedAt = (i, u, x) => {
    const I = run.interfaces;
    let c = 0, lo = 0, hi = I.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (I[m] <= x) lo = m; else hi = m; }
    c = lo;
    const a = run.speed[i], b = run.speed[i + 1];
    const v0 = a[c] + u * (b[c] - a[c]), v1 = a[c + 1] + u * (b[c + 1] - a[c + 1]);
    const f = (x - I[c]) / (I[c + 1] - I[c] || 1);
    return v0 + f * (v1 - v0);
  };
  function frame(now) {
    if (!run) return;
    const dt = last ? Math.min(0.1, (now - last) / 1000) : 0;
    last = now;
    const span = run.times[run.times.length - 1] - run.times[0];
    tau += dt;
    if (tau >= span) { tau = 0; seed(); }
    const [i, u] = rowAt(tau);
    // Walkers move at the model's speed; those past the end leave; new ones enter at the entry's flow.
    let n = 0;
    for (let k = 0; k < walkers.n; k++) {
      const x = walkers.x[k] + speedAt(i, u, walkers.x[k]) * dt;
      if (x >= run.length) continue;
      walkers.x[n] = x, walkers.lat[n] = walkers.lat[k], n++;
    }
    walkers.n = n;
    inflow += (run.flow[i][0] + u * (run.flow[i + 1][0] - run.flow[i][0])) * dt * scaleWalk;
    while (inflow >= 1 && walkers.n < 8000) { walkers.x[walkers.n] = Math.random() * 2; walkers.lat[walkers.n] = lateral(); walkers.n++; inflow -= 1; }
    for (let k = 0; k < walkers.n; k++) {
      const p = G.pointAlong(model.axis, toAxis(walkers.x[k]));
      xyz[3 * k] = p.x - p.dy * walkers.lat[k], xyz[3 * k + 1] = p.y + p.dx * walkers.lat[k], xyz[3 * k + 2] = 0.12;
      rgb[3 * k] = 0.2, rgb[3 * k + 1] = 0.33, rgb[3 * k + 2] = 0.47;
    }
    gis.setSwarm(swarmId, walkers.n, xyz, rgb);
    // The water and the lamps follow the row.
    const rho = run.density[i].map((v, c) => v + u * (run.density[i + 1][c] - v));
    current = rho;
    gis.restyle(slabId, {});
    const g = run.green[u < 0.5 ? i : i + 1];
    lamps.forEach((l) => { l.on = g[l.c]; });
    gis.restyle(lampId, {});
  }
  let current = null, lamps = [];

  function buildDrawing() {
    gis.remove(slabId), gis.remove(lampId);
    const Lm = run.length;
    toAxis = (x) => (x * model.length) / Lm;
    stations = run.position.map((x) => { const p = G.pointAlong(model.axis, toAxis(x)); return [p.x, p.y]; });
    const wide = run.width[0] > FOOTWAY + 0.1;
    const sides = wide ? [{ line: stations, width: CARRIAGEWAY + FOOTWAY }]
      : [-1, 1].map((sd) => { const e = G.edges(stations, CARRIAGEWAY + FOOTWAY / 2); return { line: sd > 0 ? e.left : e.right, width: FOOTWAY / 2 }; });
    current = run.density[0];
    gis.addSlabs(slabId, sides, {
      line: (f) => f.line, width: (f) => f.width, z: 0.25, opacity: 0.55, order: 5,
      heights: () => current.map((v) => Math.max(0.05, v * H)),
      colors: () => current.map((v) => colourOf(v) || [0.85, 0.82, 0.8]),
      pick: () => 'Swanston Street (a model)\nthe density of walkers going north',
    });
    lamps = [];
    run.crossings.forEach((x, c) => {
      if (!run.signalled[c]) return;
      for (const sd of [-1, 1]) {
        const p = G.pointAlong(model.axis, toAxis(x));
        const off = sd * (CARRIAGEWAY / 2 + FOOTWAY / 2 + 1.5);
        lamps.push({ c, xy: [p.x - p.dy * off, p.y + p.dx * off], on: 1, name: model.crossings[c] ? model.crossings[c].name : '' });
      }
    });
    gis.addPoints(lampId, lamps, { shape: 'box', position: (l) => l.xy, radius: () => 1.6, height: () => 4, base: () => 0, color: (l) => (l.on ? [0.42, 0.7, 0.45] : [0.85, 0.3, 0.27]),
      pick: (l) => `the crossing at ${l.name}: ${l.on ? 'green' : 'red'} (a 120 s cycle, green half of it; declared)` });
  }

  async function refresh(force = false) {
    if (!model) return;
    const hour = Math.round(clock.hour * 4) / 4;
    const k = `${clock.weekday}/${hour}/${opened}/${scale}`;
    if (!force && k === key) return;
    key = k;
    clearTimeout(pending);
    pending = setTimeout(async () => {
      const { demand, used } = entering(ctx, model, sensors, clock.weekday);
      if (!used.length) { status.textContent = 'no sensor counts walkers at the street\'s start'; return; }
      status.textContent = 'running the street from 06:00…';
      try {
        const r = await ask('streetwave', { hour: Math.max(6, Math.min(23.9, hour)), blockLengths: model.blockLengths, signalled: model.signalled, demand: demand.map((d) => d * scale),
          footway: FOOTWAY, carriageway: CARRIAGEWAY, rule: opened, ruleFrom: 0, ruleTo: 24, window: 240, step: 2 });
        if (key !== k) return;
        const people = r.density[0].reduce((a, v, c) => a + v * r.width[0] * (r.interfaces[c + 1] - r.interfaces[c]), 0);
        scaleWalk = people > 6000 ? 6000 / people : 1;
        const fresh = !run || run.width[0] !== r.width[0] || run.position.length !== r.position.length;
        run = r;
        if (fresh) buildDrawing();
        tau = 0, inflow = 0;
        seed();
        const ph = (demand[Math.floor(hour)] || 0) * 3600 * scale;
        status.textContent = `${hour < 6 ? 'before 06:00: ' : ''}${Math.round(ph).toLocaleString('en')} an hour enter northward (${used.map((u) => u.name).join(', ')}${scale > 1 ? `, ×${scale}` : ''}); ${Math.round(people).toLocaleString('en')} on the street; ${model.crossings.filter((c) => c.signalled).length} signalled crossings of ${model.crossings.length}${scaleWalk < 1 ? `; one walker drawn for ${fmt(1 / scaleWalk, 2)} people` : ''}`;
      } catch (e) {
        status.textContent = String(e.message || e);
      }
    }, 250);
  }

  (async () => {
    try {
      const [street, cur] = await Promise.all([ask('street', { street: STREET }), ask('currents')]);
      model = streetModel(ctx, street);
      sensors = sensorBlocks(ctx, cur.sensors);
      const mid = G.pointAlong(model.axis, model.length * 0.35);
      gis.flyTo([mid.x, mid.y], 420);
      gis.addSwarm(swarmId, { capacity: 8000, radius: 0.6, height: 2.5 });
      stop = gis.onFrame(frame);
      refresh(true);
    } catch (e) {
      status.textContent = String(e.message || e);
    }
  })();

  return {
    legend,
    setTime(c) { clock = c; refresh(); },
    describe(hit) {
      const lines = [];
      if (run && model) {
        const mean = (arr) => arr.reduce((a, v) => a + v, 0) / arr.length;
        for (const f of sensors.filter((x) => x.street === model.name)) {
          const at = G.nearestOnLine(model.axis, f.foot);
          if (at.distance > 40) continue;
          const xi = (at.s * run.length) / model.length;
          let c = 0;
          while (c + 1 < run.interfaces.length && run.interfaces[c + 1] <= xi) c++;
          const modelled = 3600 * mean(run.flow.map((row) => row[c]));
          const dir = G.pointAlong(model.axis, at.s), w = f.toward[0][0] * dir.dx + f.toward[0][1] * dir.dy > 0 ? 0 : 1;
          const counted = (w ? f.way2 : f.way1)[24 * clock.weekday + Math.floor(clock.hour)] * scale;
          lines.push(`${f.name}, ${Math.round(at.s)} m up: the model carries ${Math.round(modelled).toLocaleString('en')} an hour north, the sensor counts ${Math.round(counted).toLocaleString('en')}`);
        }
      }
      return { kind: 'model', title: 'Swanston Street, walking north (a model)', lines: lines.length ? [...lines, 'The model lets nobody turn off; the difference is what the cross streets take.'] : ['the density of walkers going north'],
        actions: [{ title: 'crowds and circulation', run: () => ctx.openTool('flow') }] };
    },
    dispose() { if (stop) stop(); clearTimeout(pending); },
  };
}
