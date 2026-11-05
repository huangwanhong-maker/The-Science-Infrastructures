// The district's walkers: who is on its streets at the clock's hour and through the day, the counters' readings
// beside what they count, and how the walkers' departures were set from the counters (the district-economy flow
// study's patterns, calibrated each way by NNLS). One thing at a time, in three tabs.

import { DAY_NAMES, ensurePrepared, hhmm, MODE_COLOURS, MODES, baseSpec, runAndWait, state } from '../shared.js';
import { brief, figures, fold, hint, tabs } from '/src/quiet.js';

export async function render(root, ctx) {
  const { h, fmt: fmt0 } = ctx.dom;
  // Whole numbers from 10 up (the shell's fmt strips the zeros of a whole number toPrecision gives: 90 → "9").
  const fmt = (x, d = 3) => (typeof x === 'number' && Math.abs(x) >= 10 ? Math.round(x).toLocaleString('en') : fmt0(x, d));
  const { lines, bars } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Walkers of the district'), h('span', { class: 'lede' }, 'Every walker an agent, set by the City\'s counters.')));
  const body = h('div', {});
  root.append(body);
  body.append(hint(h, 'preparing the district: its network, counters and streams…'));
  let alive = true;
  ctx.onDispose(() => { alive = false; });
  let P;
  try { P = await ensurePrepared(ctx); } catch (e) { return ctx.dom.failure(body, e); }
  const clock = ctx.clock, spec = baseSpec(clock);
  const progress = hint(h, `simulating the ${DAY_NAMES[clock.weekday]}…`);
  body.replaceChildren(progress);
  let R;
  try {
    R = await runAndWait(ctx, spec, (st) => { progress.textContent = `simulating the ${DAY_NAMES[clock.weekday]}: ${hhmm(Math.max(0, st.reached))} reached, ${st.peak.toLocaleString('en')} walkers at most so far`; }, () => alive);
  } catch (e) { return ctx.dom.failure(body, e); }
  if (!alive) return;
  const hour = Math.min(23, Math.floor(clock.hour)), H = R.hours[hour];
  const xs = R.hours.map((x) => x.hour + 0.5);

  const now = () => {
    const trips = H.trips.map((t, m) => ({ x: m, y: t, color: MODE_COLOURS[m], label: `${MODES[m]}: ${Math.round(t).toLocaleString('en')} trips begun` })).filter((d) => d.y > 0);
    return h('div', {},
      figures(h, [{ value: Math.round(H.walkers).toLocaleString('en'), label: `walkers on the streets, ${hhmm(hour)}–${hhmm(hour + 1)}, on average` },
        { value: H.trips.reduce((a, b) => a + b, 0).toLocaleString('en'), label: 'trips begun in the hour' },
        { value: `${fmt(R.totals.travelMinutes, 3)} min`, label: 'a trip, on average over the day' },
        { value: `${fmt(60 * R.totals.delayMinutes, 3)} s`, label: 'of it waiting at signals and in crowds' }]),
      lines([{ name: 'walkers on the streets', x: xs, y: R.hours.map((x) => x.walkers), color: ctx.dom.css('--accent'), width: 2 }],
        { height: 170, xlabel: `hour of the ${DAY_NAMES[clock.weekday]}`, ylabel: 'walkers at once', xlim: [0, 24], marks: [{ x: clock.hour, label: hhmm(clock.hour) }] }),
      h('div', { class: 'ld-label' }, `the trips begun, ${hhmm(hour)}–${hhmm(hour + 1)}, by pattern`),
      bars(trips.map((d, i) => ({ ...d, x: i })), { height: 170, labels: trips.map((d) => MODES[d.x]) }),
      brief(ctx, { kind: 'model', title: `${R.totals.departed.toLocaleString('en')} trips in the day, ${R.totals.failed ? `${R.totals.failed.toLocaleString('en')} with no way` : 'every one arrived'}`,
        statement: `The walkers spent ${Math.round(R.totals.personHours).toLocaleString('en')} person-hours on the district's streets, ${Math.round(100 * R.totals.waitingHours / R.totals.personHours)}% of them waiting at red or at a full street's mouth; ${fmt(R.totals.crowdedHours, 3)} person-hours were in crowds at Fruin's level E (above 1.08 walkers per m²).`,
        method: 'Each walker departs by its stream\'s hourly Poisson rate, chooses each next segment by the logit over every walk to its target (θ = 0.02 per metre, widths by Hoddle\'s classes of street), walks at its own speed (Weidmann\'s 1.34 ± 0.26 m/s) slowed by the density ahead (Greenshields), waits for its phase where two main streets cross (two phases of a 120 s cycle), and enters a segment within its capacity: the urban-science domain\'s movement module, run synchronously, the same on any number of threads.',
        limits: 'Signals stand where two of Hoddle\'s main streets cross, with offsets fixed and arbitrary: the City\'s timings are not in the data. Walkers inside buildings, in trams and in cars are not here. Visits dwell nowhere: a visit and its walk back are two trips, as in the flow study.' }));
  };

  const counters = () => {
    const list = h('div', {});
    for (const c of R.counters) {
      const info = P.counters.find((x) => x.id === c.id);
      const sum = (a) => a[0].map((y, i) => y + a[1][i]);
      list.append(h('div', { class: 'ld-card', style: 'margin:8px 0' },
        h('div', { class: 'ld-row', style: 'justify-content:space-between' }, h('span', { class: 'ld-title' }, c.name),
          h('button', { class: 'ld-link', onclick: () => ctx.gis.flyTo(info.at, 260) }, 'show it')),
        h('div', { class: 'ld-sub' }, `${info.street}; the calibration's R² over the week ${fmt(info.r2, 2)}`),
        lines([{ name: 'counted', x: xs, y: sum(c.observed), color: '#C98B6B', width: 2 },
          { name: 'expected (the calibrated flow)', x: xs, y: sum(c.expected), color: '#93A3AE', width: 1.4, dash: '4 3' },
          { name: 'the walkers\' reading', x: xs, y: sum(c.simulated), color: '#5B8FB9', width: 2 }],
        { height: 130, xlim: [0, 24], ylabel: 'an hour', marks: [{ x: clock.hour }] })));
    }
    return h('div', {}, ctx.chart.legend([{ name: 'counted (mean of the hour of the week, 2025–26)', color: '#C98B6B' }, { name: 'expected', color: '#93A3AE' }, { name: 'the walkers\' reading', color: '#5B8FB9' }]), list);
  };

  const how = () => {
    const C = P.calibration;
    const week = Array.from({ length: 168 }, (_, k) => k + 0.5);
    return h('div', {},
      figures(h, [{ value: fmt(C.r2Way, 3), label: 'R², each way at every counter and hour' }, { value: fmt(C.r2Gross, 3), label: 'R², both ways together' },
        { value: String(P.counters.length), label: 'counters that count each way' }, { value: String(P.entrances.length), label: 'entrances where streets leave the box' }]),
      h('div', { class: 'ld-label' }, 'trips an hour by pattern through the week, as calibrated'),
      lines(P.patterns.map((p) => ({ name: p.name, x: week, y: p.trips, color: MODE_COLOURS[p.mode], width: 1.5 })), { height: 190, xlabel: 'hour of the week (Monday 00:00 = 0)', ylabel: 'trips an hour', xlim: [0, 168] }),
      ctx.chart.legend(P.patterns.map((p) => ({ name: p.name, color: MODE_COLOURS[p.mode] }))),
      brief(ctx, { kind: 'measured', title: 'The walkers\' departures, set by the counters each way',
        statement: `Through trips, arrivals at work, departures and trips between buildings are fitted hour by hour of the week to what ${P.counters.length} counters count each way; visits to the district's ${Object.values(P.establishments).reduce((a, b) => a + b, 0).toLocaleString('en')} shops, cafés, services and places of culture are set at a tenth of the trips (${fmt(100 * C.visitShare, 3)}% as fitted), since counts cannot tell them apart.`,
        method: 'The district-economy flow study\'s patterns and calibration (research/experiments/2026-11-05-a-district-economy-as-a-flow, §§2–3), on the movement module\'s calibrate: tyche\'s NNLS each hour, each footpath\'s handedness fitted over the week, alternating, three rounds. Entrances take the gross flow of the nearest counter outside the box on their street; buildings carry their floor; CLUE\'s 2024 establishments draw visits by their kind\'s typical hours.',
        limits: 'A stream\'s law over its pairs is fixed: the through trips from an entrance go to the others in proportion to their mean flows (the study weights them by the hour\'s). Homes are not told from offices. The fit is in sample.' }),
      fold(h, 'the counters left out', ...(P.leftOut.length ? P.leftOut.map((s) => h('p', {}, s)) : [h('p', {}, 'none')])),
      fold(h, 'the entrances', ...P.entrances.map((e) => h('p', {}, `${e.street || 'a way'}: ${Math.round(e.mean).toLocaleString('en')} an hour on average (${e.source}${e.counter >= 0 ? `, ${e.distance} m away` : ''})`))));
  };

  const T = tabs(h, [{ name: 'Now', panel: now }, { name: 'Counters', panel: counters }, { name: 'How they are set', panel: how }], { remember: 'hf-walkers' });
  body.replaceChildren(T.el);
}
