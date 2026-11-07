// How far in so many minutes, from wherever is picked: the ways in view, at the speeds of a mode, to the engine; the times come back
// per point, and the ways are drawn to fill as a wave from the place picked. Speeds are declared here, and they are assumptions: nothing of
// signals, gradients or traffic is in them, so that a time is a lower bound on the journey.

const KMH = (v) => v / 3.6;
const MODES = {
  foot: { name: 'on foot', walk: 1.34, speed: (cls, rail) => (rail || ['motorway', 'motorway_link', 'trunk', 'trunk_link'].includes(cls) ? 0 : cls === 'steps' ? 0.6 : 1.34), upto: 60,
    note: 'At 1.34 m/s (Weidmann 1993), slower on steps; motorways closed. No signals, no gradients.' },
  road: { name: 'by road', walk: 1.34, upto: 30, note: 'At a free-flow speed by class of road, 25 km/h on a residential street to 90 on a motorway; no traffic, no signals, no turns or one-way streets: a lower bound on a journey.',
    speed: (cls, rail) => (rail ? 0 : { motorway: KMH(90), motorway_link: KMH(50), trunk: KMH(70), trunk_link: KMH(40), primary: KMH(50), primary_link: KMH(35), secondary: KMH(40), secondary_link: KMH(30), tertiary: KMH(35), tertiary_link: KMH(25), residential: KMH(25), unclassified: KMH(25), living_street: KMH(10), service: KMH(15) }[cls] ?? 0) },
};

export async function layer(ctx) {
  const { gis, frame, earth, world, dom: { h } } = ctx;
  const id = ctx.layerId(':wave');
  const mark = ctx.layerId(':here');
  let mode = 'foot', from = null, dead = false, token = 0, minutes = MODES.foot.upto, shown = 0, playing = null;
  const status = h('div', { class: 'ld-quiet' }, 'Pick a place on the map.');
  const readout = h('div', { class: 'reach-readout' });
  const slider = h('input', { type: 'range', min: 1, max: 60, step: 1, value: 60, 'aria-label': 'minutes', oninput: () => { stop(); set(+slider.value); } });
  const label = h('span', { class: 'mono' }, '');
  const modeKeys = h('div', { class: 'keys-tight' }, Object.entries(MODES).map(([k, M]) => h('button', { class: k === mode ? 'on' : '', onclick: (e) => {
    mode = k;
    for (const b of e.currentTarget.parentElement.children) b.classList.toggle('on', b === e.currentTarget);
    note.textContent = M.note;
    slider.max = M.upto;
    if (from) compute();
  } }, M.name)));
  const note = h('div', { class: 'ld-quiet' }, MODES.foot.note);
  const set = (m) => { shown = m; label.textContent = `${Math.round(m)} min`; slider.value = m; gis.setWave(id, m * 60); };
  const stop = () => { if (playing) { cancelAnimationFrame(playing); playing = null; } };
  function play(until) {
    stop();
    const t0 = performance.now(), seconds = 5;
    const step = (now) => {
      const u = Math.min(1, (now - t0) / (seconds * 1000));
      set(until * u);
      playing = u < 1 && !dead ? requestAnimationFrame(step) : null;
    };
    playing = requestAnimationFrame(step);
  }

  async function compute() {
    const me = ++token;
    const M = MODES[mode];
    const all = world.lines();
    if (!all.length) { status.textContent = 'The ways here are still coming: try again in a moment.'; return; }
    status.textContent = `the minutes over ${all.length.toLocaleString('en')} ways…`;
    const ways = [];
    const kept = [];
    for (const l of all) {
      const speed = M.speed(l.cls, l.rail);
      const row = [speed];
      for (let i = 0; i < l.xyz.length; i += 3) row.push(Math.round(l.xyz[i] * 10) / 10, Math.round(l.xyz[i + 1] * 10) / 10);
      ways.push(row);
      kept.push(l);
    }
    let a;
    try {
      a = await ctx.ask('reach', { ways, from: [from], walk: M.walk, limits: Array.from({ length: 6 }, (_, k) => ((k + 1) * M.upto * 60) / 6) });
    } catch (e) { status.textContent = String(e.message || e); return; }
    if (me !== token || dead) return;
    // The wave: every way a ribbon, each of its points carrying its time.
    const far = Math.max(2.6, earth.view.range * 0.0018);
    gis.addWave(id, kept, a.times, { width: (l) => Math.max(far, l.width), colours: ctx.theme === 'night' ? { near: [1.0, 0.72, 0.84], far: [0.42, 0.26, 0.62], crest: [1.0, 0.97, 0.9] } : { near: [0.62, 0.12, 0.38], far: [0.56, 0.42, 0.74], crest: [0.98, 0.3, 0.5] }, upto: M.upto * 60 });
    readout.replaceChildren(...a.limits.map((lim, k) => h('div', {}, h('span', {}, `${Math.round(lim / 60)}′`), h('b', {}, `${(a.metres[k] / 1000).toFixed(a.metres[k] < 1e5 ? 1 : 0)} km`))));
    status.textContent = `${Math.round(a.gap)} m to the nearest way; ways in view only${a.stitched ? `; ${a.stitched.toLocaleString('en')} ends joined across tiles` : ''}.`;
    slider.max = M.upto;
    play(M.upto);
  }

  const legend = h('div', { class: 'card' }, h('b', {}, 'How far in so many minutes'), modeKeys,
    h('label', { class: 'view-slider reach-slider' }, slider, label), readout, note, status,
    h('div', { class: 'ld-quiet' }, 'Ways: Overture Maps, OpenStreetMap'));
  ctx.onDispose(() => { dead = true; stop(); });
  return {
    legend,
    // A click anywhere sets the place.
    pickAt(lat, lon) {
      const p = frame.toLocal(lat, lon);
      from = [Math.round(p[0] * 10) / 10, Math.round(p[1] * 10) / 10];
      const z = earth.city ? earth.view.height - earth.city.height : 0;
      gis.addPoints(mark, [{ at: [p[0], p[1]] }], { position: (f) => f.at, radius: () => Math.max(2, earth.view.range * 0.0016), base: () => z, height: () => Math.max(30, earth.view.range * 0.05), color: () => (ctx.theme === 'night' ? [1, 0.75, 0.86] : [0.76, 0.28, 0.5]), pick: false });
      compute();
      return { kind: 'place', title: 'From here', lines: [`${MODES[mode].name}`], at: { lat, lon } };
    },
    dispose() { dead = true; stop(); gis.remove(id); gis.remove(mark); },
  };
}
