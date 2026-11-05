// The City's counters in the district, the walkers' reading beside the count (a model beside a record): at each
// counter two columns, sand for what the counter counts on average at the clock's hour of the week (both ways, over the
// reference year 2025–26, public holidays left out) and sea blue for what the run's walkers passing its line would make
// it read (its footpath's share of the street, and its side's share of each way, as calibrated), a metre for every
// hundred an hour. The run is the day as built until it is done, then the run the map draws.

import { ask, ensurePrepared, hhmm, on, rgb, runAndWait, specAt, state, baseSpec } from '../shared.js';

const COUNTED = '#C98B6B', SIMULATED = '#5B8FB9', EXPECTED = '#93A3AE';

export async function layer(ctx) {
  const { gis, dom } = ctx;
  const { h, fmt } = dom;
  const id = ctx.layerId(':columns');
  let clock = ctx.clock, alive = true, results = new Map();
  const status = h('div', { class: 'ld-quiet' }, 'waiting for the day\'s walkers…');
  const legend = h('div', { class: 'card' }, h('b', {}, 'Counters: counted and simulated'),
    h('div', { class: 'ld-swatches' }, h('span', {}, h('i', { class: 'ld-sw block', style: `background:${COUNTED}` }), 'counted (the mean of this hour of the week)'),
      h('span', {}, h('i', { class: 'ld-sw block', style: `background:${SIMULATED}` }), 'the walkers\' reading')),
    h('div', { class: 'ld-quiet' }, 'a metre for every hundred an hour, both ways'), status);

  const resultFor = (spec) => results.get(JSON.stringify(spec));
  async function fetchResult(spec) {
    const k = JSON.stringify(spec);
    if (results.has(k)) return results.get(k);
    results.set(k, null);
    const R = await runAndWait(ctx, spec, (st) => { status.textContent = `the walkers reach ${hhmm(Math.max(spec.from, st.reached))}; the readings come when the run is done`; }, () => alive);
    results.set(k, R);
    return R;
  }

  const features = () => {
    const P = state.prep;
    const out = [];
    for (const c of P.counters) {
      const [ex, ey] = [Math.cos(0.6), Math.sin(0.6)];
      out.push({ c, which: 'counted', at: [c.at[0] - 2.2 * ex, c.at[1] - 2.2 * ey] }, { c, which: 'simulated', at: [c.at[0] + 2.2 * ex, c.at[1] + 2.2 * ey] });
    }
    return out;
  };
  // A counter's numbers at the clock's hour from a run's result (its window holding the hour).
  const at = (R, c) => {
    if (!R) return null;
    const k = R.counters.findIndex((x) => x.id === c.id), i = Math.floor(clock.hour) - R.status.from;
    if (k < 0 || i < 0 || i >= R.counters[k].observed[0].length) return null;
    const x = R.counters[k];
    return { counted: x.observed[0][i] + x.observed[1][i], expected: x.expected[0][i] + x.expected[1][i], simulated: x.simulated[0][i] + x.simulated[1][i], x };
  };

  function paint() {
    if (!state.prep) return;
    const want = specAt(clock);
    const R = resultFor(want.spec) || resultFor(baseSpec(clock));
    gis.addPoints(id, features(), {
      position: (f) => f.at, radius: () => 2.4, shape: 'box',
      height: (f) => { const v = at(R, f.c); return v ? Math.max(0.2, (f.which === 'counted' ? v.counted : v.simulated) / 100) : 0.2; },
      color: (f) => rgb(f.which === 'counted' ? COUNTED : SIMULATED), opacity: 0.95,
      pick: (f) => { const v = at(R, f.c); return v ? `${f.c.name}, ${hhmm(Math.floor(clock.hour))}: counted ${Math.round(v.counted)} an hour, the walkers' reading ${Math.round(v.simulated)} (the calibration expects ${Math.round(v.expected)})` : f.c.name; },
    });
    if (R) status.textContent = `${R === resultFor(want.spec) ? want.label : 'the day as built'}, ${hhmm(Math.floor(clock.hour))}–${hhmm(Math.floor(clock.hour) + 1)}`;
  }

  async function refresh() {
    try {
      await ensurePrepared(ctx);
      paint();
      const want = specAt(clock);
      const tasks = [fetchResult(baseSpec(clock))];
      if (want.scenario) tasks.push(fetchResult(want.spec));
      await Promise.all(tasks);
      if (alive) paint();
    } catch (e) {
      status.textContent = String(e.message || e);
    }
  }
  const off = on((w) => { if (w === 'scenario') refresh(); });
  refresh();

  return {
    legend,
    setTime(c) { const day = c.weekday !== clock.weekday; clock = c; if (day) refresh(); else paint(); },
    describe(hit) {
      if (hit.layer !== id) return null;
      const f = hit.feature, want = specAt(clock), R = resultFor(want.spec) || resultFor(baseSpec(clock)), v = at(R, f.c);
      if (!v) return { title: f.c.name, lines: ['the walkers\' readings come when the run is done'] };
      const xs = v.x.observed[0].map((_, i) => R.status.from + i + 0.5);
      const sum = (a) => a[0].map((y, i) => y + a[1][i]);
      const chart = ctx.chart.lines([
        { name: 'counted', x: xs, y: sum(v.x.observed), color: COUNTED, width: 2 },
        { name: 'the calibration\'s expectation', x: xs, y: sum(v.x.expected), color: EXPECTED, width: 1.4, dash: '4 3' },
        { name: 'the walkers\' reading', x: xs, y: sum(v.x.simulated), color: SIMULATED, width: 2 },
      ], { height: 150, width: 360, xlabel: 'hour', ylabel: 'people an hour, both ways', marks: [{ x: clock.hour, label: hhmm(clock.hour) }] });
      return { kind: 'counter', title: f.c.name, lines: [`${f.c.street}; counts toward ${f.c.toward.join(' and ')}; it reads ${Math.round(100 * f.c.share * (1 + f.c.hand))}% of the street's walkers toward ${f.c.toward[0]} and ${Math.round(100 * f.c.share * (1 - f.c.hand))}% toward ${f.c.toward[1]}`,
        `at ${hhmm(Math.floor(clock.hour))}: counted ${Math.round(v.counted)}, the walkers' reading ${Math.round(v.simulated)}, expected ${Math.round(v.expected)}`], extra: dom.h('div', {}, chart) };
    },
    dispose() { alive = false; off(); gis.remove(id); },
  };
}
