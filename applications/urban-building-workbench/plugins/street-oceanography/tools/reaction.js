// Stores that work together: which uses keep each other's company beyond the blocks they share; two stores as
// a pathway in space and time — one raises a need (a bridal store, rings) that the other meets (a jeweller) —
// solved exactly for a visitor who wanders or searches, at the separation and hours of two stores picked on the
// map, and drawn on the map: the first store's visitors spreading along the street minute by minute, the
// water rose where the need has been raised, an arc joining the two stores as thick as their coupling, a
// column at the second rising with those who have bought; and the spacing of clusters of a use along a street.

// A store's day: its first opening and last closing among the intervals that start on the weekday.
const dayHours = (hours, weekday) => {
  const ivs = (hours || []).filter(([a]) => a >= 24 * weekday && a < 24 * (weekday + 1));
  if (!ivs.length) return null;
  return [Math.min(...ivs.map(([a]) => a)) - 24 * weekday, Math.max(...ivs.map(([, b]) => b)) - 24 * weekday];
};

export async function render(root, ctx) {
  const { h, svg, busy, failure, fmt, css, hover } = ctx.dom;
  const { lines, legend } = ctx.chart;
  const { colourMap } = await import('/visualization/styles/index.js');
  const diverging = colourMap('lavender-rust').at;
  const rgb = (c) => `rgb(${c.map((v) => Math.round(v * 255)).join(',')})`;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Stores that work together')));

  // ---- two stores as a pathway ----
  const P = { diffusivity: 1, search: 0.15, leaveMinutes: 30, storeMinutes: 15, catalysisMinutes: 20, purchaseMinutes: 30, bridalOpen: 11, bridalClose: 19, ringOpen: 11, ringClose: 19, windowSeparation: 30 };
  const picked = [];
  const pickNote = h('p', { class: 'note' }, 'Pick on the map the store that raises a need (a bridal store), then the one that meets it (a jeweller), on the same street — or keep the declared example.');
  const panel = h('div', {});
  const res = h('div', {});
  const slider = (key, name, min, max, step, unit = '') => {
    const input = h('input', { type: 'range', min, max, step, value: P[key] });
    const o = h('output', {}, `${fmt(P[key])}${unit}`);
    input.addEventListener('input', () => { o.textContent = `${fmt(+input.value)}${unit}`; });
    input.addEventListener('change', () => { P[key] = +input.value; runPath(); if (plumeArgs) drawPlume(...plumeArgs); });
    return h('label', { class: 'slider' }, name, input, o);
  };
  panel.append(slider('diffusivity', 'wandering', 0.25, 4, 0.25, ' m²/s'), slider('search', 'walks toward it at', 0, 0.6, 0.05, ' m/s'),
    slider('leaveMinutes', 'stays on the street', 10, 90, 5, ' min'), slider('storeMinutes', 'stays in a store', 5, 40, 5, ' min'),
    slider('catalysisMinutes', 'the need arises in', 5, 60, 5, ' min'), slider('purchaseMinutes', 'buys in', 5, 90, 5, ' min'));
  const mapNote = h('p', { class: 'caption' });
  root.append(h('h2', {}, 'Two stores as a pathway'), pickNote, h('div', { class: 'panel' }, panel), mapNote, res);
  const takePick = (s) => {
    if (!s || s.kind !== 'holding') return false;
    if (picked.length === 2 || (picked.length === 1 && picked[0].street !== s.street)) picked.length = 0;
    picked.push(s);
    if (picked.length === 1) pickNote.textContent = `First: ${s.name} (${s.use}). Now the store that meets the need.`;
    if (picked.length === 2) {
      const [a, b] = picked, wd = ctx.clock.weekday;
      const sep = Math.max(10, Math.abs(a.s - b.s));
      P.windowSeparation = sep;
      const ha = dayHours(a.hours, wd), hb = dayHours(b.hours, wd);
      if (ha) [P.bridalOpen, P.bridalClose] = [ha[0], Math.min(24, ha[1])];
      if (hb) [P.ringOpen, P.ringClose] = [hb[0], Math.min(24, hb[1])];
      pickNote.textContent = `${a.name} (${a.use}) raises the need; ${b.name} (${b.use}) meets it, ${Math.round(sep)} m along the street. On ${['Mondays', 'Tuesdays', 'Wednesdays', 'Thursdays', 'Fridays', 'Saturdays', 'Sundays'][wd]}: ${ha ? `${fmt(ha[0], 3)}–${fmt(ha[1], 3)} h` : 'no hours declared'} and ${hb ? `${fmt(hb[0], 3)}–${fmt(hb[1], 3)} h` : 'no hours declared'}.`;
      runPath();
      drawPlume(a, b);
    }
    return true;
  };
  takePick(ctx.selection);
  ctx.onSelect(takePick);
  const runPath = async () => {
    busy(res, 'solving the visitor\'s walk exactly…');
    let d;
    try { d = await ctx.ask('pathway', P); } catch (e) { return failure(res, e); }
    const c = d.coupling, w = d.walk, win = d.window;
    const interp = (xs, ys, x) => { for (let i = 1; i < xs.length; i++) if (xs[i] >= x) { const t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]); return ys[i - 1] + t * (ys[i] - ys[i - 1]); } return ys[ys.length - 1]; };
    const atSep = interp(c.separation, c.yield, P.windowSeparation);
    const atHours = win.yield[win.offset.indexOf(0)];
    res.replaceChildren(
      lines([{ name: 'searching', x: c.separation, y: c.yield, color: css('--accent') }, { name: 'wandering', x: c.separation, y: w.yield, color: css('--minor'), dash: '5 4' }],
        { xlabel: 'distance between the stores (m)', ylabel: 'chance the need is met', height: 210, marks: [{ x: P.windowSeparation, label: `${Math.round(P.windowSeparation)} m` }] }),
      legend([{ name: `walks toward it at ${fmt(P.search)} m/s`, color: css('--accent') }, { name: 'wanders', color: css('--minor') }]),
      lines([{ name: 'chance', x: win.offset, y: win.yield, color: css('--series-2') }],
        { xlabel: 'hours the second store\'s day is shifted', ylabel: `chance at ${Math.round(win.separation)} m`, height: 170 }),
      ctx.finding({ title: `${fmt(100 * atSep, 3)}% of the first store's visitors find what they now need (${fmt(100 * atHours, 3)}% with the declared hours)`, kind: 'model',
        statement: `Together the two stores sell something neither sells alone, but only near each other and while both are open. The chance halves every ${fmt(w.length * Math.LN2, 3)} m for a visitor who wanders (e-folding ${fmt(w.length, 3)} m, as a killed random walk gives exactly: ${fmt(d.scales.walkClosedForm, 3)} m) and every ${fmt(c.length * Math.LN2, 3)} m for one who knows where the second store is — legibility matters. Shift the second store's hours away and the pathway closes at any distance. Per visitor, the first store creates ${fmt(d.capture.created, 4)} of the second's takings against ${fmt(d.capture.own, 4)} of its own, which rents set by its own takings never count.`,
        method: 'The visitor\'s walk among frontages 10 m apart, its time in each store, the need and the purchase as a Markov jump process, solved exactly (absorption probabilities; the dwell is phase-type); hours by stepping the process through the day.',
        limits: 'Every rate is declared; values of sales are illustrative.' }));
  };
  runPath();

  // ---- on the map: the first store's visitors spreading along the street ----
  const IDS = ['reaction:plume', 'reaction:bond', 'reaction:bought'];
  let stopPlume = null, plumeArgs = null;
  const clearMap = () => { if (stopPlume) stopPlume(); stopPlume = null; for (const id of IDS) ctx.gis.remove(id); };
  ctx.onDispose(clearMap);
  const SEA = [0.31, 0.54, 0.55], ROSE = [0.8, 0.5, 0.36];
  async function drawPlume(a, b) {
    plumeArgs = [a, b];
    clearMap();
    let st, d;
    try {
      st = await ctx.ask('street', { street: a.street });
      d = await ctx.ask('plume', { ...P, length: st.length, bridal: a.s, ring: b.s, horizon: 3600, steps: 60 });
    } catch (e) { mapNote.textContent = String(e.message || e); return; }
    if (plumeArgs[0] !== a || plumeArgs[1] !== b) return;
    const { G } = ctx;
    const axis = st.axis.map(([lat, lon]) => ctx.frame.toLocal(lat, lon).slice(0, 2));
    const toAxis = (x) => (x * G.lineLength(axis)) / st.length;
    const at = (x) => { const q = G.pointAlong(axis, toAxis(x)); return [q.x, q.y]; };
    const stations = d.position.map(at);
    let k = 0;
    const heights = () => d.position.map((_, i) => Math.min(45, 2 + 400 * (d.errand[k][i] + d.need[k][i])));
    const colours = () => d.position.map((_, i) => {
      const e = d.errand[k][i], m = d.need[k][i], u = e + m > 1e-12 ? m / (e + m) : 0;
      return SEA.map((c, j) => c + u * (ROSE[j] - c));
    });
    ctx.gis.addSlabs(IDS[0], [{ line: stations }], { line: (f) => f.line, width: () => 14, z: 0.3, opacity: 0.72, order: 6, heights, colors: colours,
      pick: () => 'the first store\'s visitors along the street (a model): sea glass with the errand they came for, rose with the need the store raised' });
    const A = at(d.bridal), B = at(d.ring);
    ctx.gis.addArcs(IDS[1], [{ A, B }], { from: (f) => f.A, to: (f) => f.B, height: () => 18 + Math.abs(d.ring - d.bridal) / 3, radius: () => 0.5 + 8 * d.yield, base: 3,
      color: () => ROSE, pick: () => `the two stores' coupling: ${fmt(100 * d.yield, 3)}% of the first store's visitors buy at the second` });
    ctx.gis.addPoints(IDS[2], [{ xy: B }], { shape: 'cylinder', position: (f) => f.xy, radius: () => 5, height: () => 1 + 300 * d.bought[k], color: () => ROSE,
      pick: () => `bought so far: ${fmt(100 * d.bought[k], 3)}%` });
    ctx.gis.flyTo([(A[0] + B[0]) / 2, (A[1] + B[1]) / 2], Math.max(260, 2.5 * Math.hypot(B[0] - A[0], B[1] - A[1])));
    let last = 0;
    const show = () => {
      const on = d.errand[k].reduce((x, v) => x + v, 0) + d.need[k].reduce((x, v) => x + v, 0);
      mapNote.textContent = `On the map, ${Math.round(d.times[k] / 60)} minutes after a visitor leaves ${a.name}: ${fmt(100 * on, 3)}% still on the street (rose where the need is raised), ${fmt(100 * d.bought[k], 3)}% have bought at ${b.name}, ${fmt(100 * d.left[k], 3)}% have left. It loops every 15 seconds.`;
    };
    show();
    stopPlume = ctx.gis.onFrame((now) => {
      if (now - last < 250) return;
      last = now;
      k = (k + 1) % d.times.length;
      ctx.gis.restyle(IDS[0], {});
      ctx.gis.restyle(IDS[2], {});
      show();
    });
  }

  // ---- which uses keep company ----
  const coloc = h('div', { class: 'section' });
  root.append(coloc);
  const where = 'melbourne';
  const box = h('div', {});
  coloc.append(h('h2', {}, 'Which uses keep company'), box);
  busy(box, 'comparing with uses shuffled within each block…');
  try {
    const d = await ctx.ask('affinity', { where, cell: 100, permutations: 499 });
    const n = d.uses.length, idx = d.uses.map((_, i) => i).filter((i) => d.count[i] >= 3);
    const cell = 26, W = 170 + cell * idx.length, H = 26 + cell * idx.length + 120;
    const m = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
    idx.forEach((a, r) => {
      m.append(svg('text', { x: 162, y: 26 + r * cell + cell / 2 + 4, 'text-anchor': 'end', class: 'mark-label' }, d.uses[a]));
      idx.forEach((b, col) => {
        const q = d.clq[a * n + b], la = Number.isFinite(q) && q > 0 ? Math.log(q) : NaN;
        const rect = svg('rect', { x: 170 + col * cell, y: 26 + r * cell, width: cell - 2, height: cell - 2, rx: 3, fill: Number.isFinite(la) ? rgb(diverging(Math.max(-1, Math.min(1, la / 1.5)))) : css('--sunken') });
        hover(rect, `${d.uses[a]} → ${d.uses[b]}: ${fmt(q, 3)} times as often as chance\np = ${fmt(d.pUpper[a * n + b], 2)}; within blocks p = ${fmt(d.stratifiedPUpper[a * n + b], 2)}`);
        m.append(rect);
        if (d.pUpper[a * n + b] < 0.01 && d.stratifiedPUpper[a * n + b] < 0.05) m.append(svg('text', { x: 170 + col * cell + cell / 2 - 1, y: 26 + r * cell + cell / 2 + 4, 'text-anchor': 'middle', class: 'mark-label' }, '◆'));
      });
      m.append(svg('text', { x: 170 + r * cell + cell / 2, y: 26 + idx.length * cell + 8, transform: `rotate(55 ${170 + r * cell + cell / 2} ${26 + idx.length * cell + 8})`, class: 'mark-label' }, d.uses[a]));
    });
    const bi = d.uses.indexOf('bridal-and-wedding'), ji = d.uses.indexOf('jewellery');
    box.replaceChildren(m, ctx.finding({ title: 'Neighbours, against chance and against shared blocks', kind: 'measured',
      statement: `Each row: how often a holding's nearest neighbour has the column's use, against chance (rose: more often). ◆ marks pairs that keep company even when uses are shuffled only within each 100-metre block — closer than sharing good blocks explains. ${bi >= 0 && ji >= 0 && d.count[bi] > 0 && d.count[ji] > 0 ? `Bridal → jewellery: ${fmt(d.clq[bi * n + ji], 3)}× (p ${fmt(d.pUpper[bi * n + ji], 2)}; within blocks ${fmt(d.stratifiedPUpper[bi * n + ji], 2)}).` : ''}`,
      method: `Co-location quotients of ${d.points} holdings (${d.where}), ${d.permutations} shuffles.`, limits: 'Nearest neighbours only; the bridal class is read from trading names.' }));
  } catch (e) { failure(box, e); }
}
