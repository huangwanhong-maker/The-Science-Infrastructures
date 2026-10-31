// Crowds and circulation: the crowd on Swanston Street, its blocks and crossings from the city model and its
// walkers entering at the counted rate, as a kinematic wave on an ordinary day and with the carriageway opened
// to walkers, its numbers kept exactly; and movement among the intersections of a district of blocks, which
// can circulate round the blocks only when something other than an attraction drives it, with the time a
// walker takes to leave.

import { sensorBlocks } from '../layers/blocks.js';
import { entering, streetModel } from '../layers/street-model.js';

export async function render(root, ctx) {
  const { h, svg, busy, failure, fmt, css, hover } = ctx.dom;
  const { lines, raster, ramp } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Crowds and circulation'),
    h('span', { class: 'lede' }, 'Models with declared parameters: what a rule or a street pattern would do to a crowd, to set beside what the records show. The waves layer draws the first on the map.')));

  const slider = (P, key, name, min, max, step, unit, run) => {
    const input = h('input', { type: 'range', min, max, step, value: P[key] });
    const o = h('output', {}, `${fmt(P[key])}${unit}`);
    input.addEventListener('input', () => { o.textContent = `${fmt(+input.value)}${unit}`; });
    input.addEventListener('change', () => { P[key] = +input.value; run(); });
    return h('label', { class: 'slider' }, name, input, o);
  };

  // ---- Swanston Street ----
  const [street, cur] = await Promise.all([ctx.ask('street', { street: 'melbourne-swanston' }), ctx.ask('currents')]);
  const model = streetModel(ctx, street), sensors = sensorBlocks(ctx, cur.sensors);
  const P = { scale: 1, footway: 10, carriageway: 12, green: 0.5 };
  const res = h('div', {});
  const runDay = async () => {
    busy(res, 'advancing the waves…');
    const { demand, used } = entering(ctx, model, sensors, ctx.clock.weekday);
    let d;
    try {
      d = await ctx.ask('streetday', { blockLengths: model.blockLengths, signalled: model.signalled, demand: demand.map((x) => x * P.scale), footway: P.footway,
        carriageway: P.carriageway, green: P.green, from: 7, to: 21, ruleTo: 21, snapshot: 300 });
    } catch (e) { return failure(res, e); }
    const hi = Math.max(d.ordinary.maxDensity, d.rule.maxDensity, 0.05);
    const pic = (o, name) => h('div', {}, h('h3', {}, name), raster(o.density, { map: 'orchid', lo: 0, hi, height: 200 }),
      h('p', { class: 'caption' }, `Densest ${fmt(o.maxDensityRuleHours, 3)} per m² from 12:00; longest queue ${fmt(o.queue, 3)} m; ${fmt(o.meanTime / 60, 3)} min on the street on average.`));
    const day = ['Mondays', 'Tuesdays', 'Wednesdays', 'Thursdays', 'Fridays', 'Saturdays', 'Sundays'][ctx.clock.weekday];
    res.replaceChildren(h('div', { class: 'grid cols-2' }, pic(d.ordinary, 'An ordinary day'), pic(d.rule, 'Carriageway opened from noon')), ramp('orchid', 0, hi, 'people per m²'),
      h('p', { class: 'caption' }, `Rows: time from 07:00 (top) to 21:00; columns: along Swanston Street from Flinders Street (left) to La Trobe Street, ${Math.round(model.length)} m, crossing ${model.crossings.map((c) => c.name).join(', ')}.`),
      ctx.finding({ title: 'The same arrivals on two streets', kind: 'model',
        statement: `Walkers enter northward at the mean rate counted at ${used.map((u) => u.name).join(' and ')} on ${day}${P.scale > 1 ? `, ${P.scale} times over` : ''}. Opening the carriageway turns ${fmt(P.footway)} m of footway into ${fmt(P.footway + P.carriageway)} m and stops the crossings. From noon the densest point goes from ${fmt(d.ordinary.maxDensityRuleHours, 3)} to ${fmt(d.rule.maxDensityRuleHours, 3)} people per m², and the mean time on the street from ${fmt(d.ordinary.meanTime / 60, 3)} to ${fmt(d.rule.meanTime / 60, 3)} min. Red crossings make queues whose tails run back down the street as waves. Nobody is lost or created: the count closes to ${fmt(Math.abs(d.ordinary.budgetResidual), 2)} people.`,
        method: 'Godunov\'s scheme for the Lighthill–Whitham–Richards law on a network of the street\'s blocks and crossings, from the city model (physis continuum); walkers at 1.34 m/s when free, jammed at 5.4 per m² (Weidmann 1993); a signal at each street, none at a lane.',
        limits: 'Widths and signals are declared; nobody turns into a cross street, so the flow further up overstates what the sensors there count. The waves layer draws this street on the map, second by second.' }));
  };
  root.append(h('div', { class: 'section' }, h('h2', {}, 'Swanston Street, ordinary and opened'),
    h('div', { class: 'panel' }, slider(P, 'scale', 'the crowd, against the counts', 1, 6, 1, '×', runDay), slider(P, 'footway', 'footways', 4, 20, 1, ' m', runDay),
      slider(P, 'carriageway', 'carriageway opened', 0, 20, 1, ' m', runDay), slider(P, 'green', 'crossings green', 0.3, 0.8, 0.05, '', runDay)),
    res));
  runDay();

  // ---- the district ----
  const Q = { nx: 5, ny: 4, attraction: 2, circulation: 0.6, holeCirculation: 1.2, holes: [7], station: 0 };
  const dres = h('div', {});
  const holeKey = h('input', { type: 'checkbox', checked: true });
  const runDist = async () => {
    busy(dres);
    let d;
    try { d = await ctx.ask('district', Q); } catch (e) { return failure(dres, e); }
    const W = 460, H = 360, pad = 34, nx = d.nx, ny = d.ny;
    const X = (v) => pad + (v % (nx + 1)) * ((W - 2 * pad) / nx), Y = (v) => H - pad - Math.floor(v / (nx + 1)) * ((H - 2 * pad) / ny);
    const pic = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: '100%', style: `max-width:${W}px` });
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const hole = d.holes.includes(j * nx + i);
      pic.append(svg('rect', { x: X(i) + 6, y: Y((j + 1) * (nx + 1)) + 6, width: (W - 2 * pad) / nx - 12, height: (H - 2 * pad) / ny - 12, rx: 4, fill: hole ? css('--minor-soft') : css('--sunken'), stroke: hole ? css('--minor') : 'none' }));
    }
    const maxJ = Math.max(1e-12, ...d.edges.map((e) => Math.abs(e.current)));
    for (const e of d.edges) {
      const [a, b] = e.current >= 0 ? [e.tail, e.head] : [e.head, e.tail];
      const line = svg('line', { x1: X(a), y1: Y(a), x2: X(b), y2: Y(b), stroke: css('--accent'), 'stroke-width': 0.8 + (5 * Math.abs(e.current)) / maxJ, 'stroke-linecap': 'round', opacity: 0.85 });
      hover(line, `net flow ${fmt(e.current, 3)} /s`);
      pic.append(line);
      const mx = (X(a) + X(b)) / 2, my = (Y(a) + Y(b)) / 2, ang = Math.atan2(Y(b) - Y(a), X(b) - X(a));
      if (Math.abs(e.current) / maxJ > 0.15) pic.append(svg('path', { d: `M${mx + 6 * Math.cos(ang)},${my + 6 * Math.sin(ang)} L${mx - 5 * Math.cos(ang) + 4 * Math.sin(ang)},${my - 5 * Math.sin(ang) - 4 * Math.cos(ang)} L${mx - 5 * Math.cos(ang) - 4 * Math.sin(ang)},${my - 5 * Math.sin(ang) + 4 * Math.cos(ang)} Z`, fill: css('--accent-deep') }));
    }
    const maxPi = Math.max(...d.pi);
    d.pi.forEach((p, v) => { const c = svg('circle', { cx: X(v), cy: Y(v), r: 3 + (8 * p) / maxPi, fill: css('--minor'), opacity: 0.8 }); hover(c, `share of walkers here: ${fmt(p, 3)}`); pic.append(c); });
    const tr = d.transit;
    dres.replaceChildren(pic, h('p', { class: 'caption' }, 'Lines: the net flow between intersections (width), arrows its direction; dots: where walkers are.'),
      ctx.finding({ title: d.detailedBalance ? 'No net circulation' : 'An attraction alone never makes a crowd go round', kind: 'proved',
        statement: `However strongly the station pulls, with no drive round the blocks every flow balances its reverse and no net current circles a block (Kolmogorov's criterion). Only the part of the drive that is not a pull toward places — round the blocks, round the park — makes walkers circulate, and the station's pull changes none of the loops' drives (to ${fmt(d.affinityGapWithoutAttraction, 2)}). Here the irreversibility (entropy production) is ${fmt(d.entropyProduction, 4)} /s.`,
        method: 'A Markov jump process on the district\'s streets and blocks; its drive split exactly into pull, circulation round blocks and circulation round holes (a Hodge decomposition; physis).',
        against: 'A net current round a block in observed flows: a pull alone cannot make one.', limits: 'Rates are declared; real movement needs trajectories.' }),
      lines([{ name: 'time to leave', x: tr.time.map((t) => t / 60), y: tr.density.map((v) => v * 60), color: css('--series-2') }],
        { xlabel: 'minutes in the district', ylabel: 'probability per minute', height: 180, marks: [{ x: tr.median / 60, label: `median ${fmt(tr.median / 60, 3)} min` }] }),
      h('p', { class: 'caption' }, `From the station to either far corner: mean ${fmt(tr.mean / 60, 3)} min, 90% within ${fmt(tr.p90 / 60, 3)} min (an exact phase-type law).`));
  };
  holeKey.addEventListener('change', () => { Q.holes = holeKey.checked ? [7] : []; runDist(); });
  root.append(h('div', { class: 'section' }, h('h2', {}, 'Circulation in a district of blocks'),
    h('div', { class: 'panel' }, slider(Q, 'attraction', 'pull of the station', 0, 4, 0.25, '', runDist), slider(Q, 'circulation', 'drive round blocks', 0, 2, 0.1, '', runDist),
      slider(Q, 'holeCirculation', 'drive round the park', 0, 3, 0.1, '', runDist), h('label', { class: 'pick' }, holeKey, ' a park in the second row')),
    dres));
  runDist();
}
