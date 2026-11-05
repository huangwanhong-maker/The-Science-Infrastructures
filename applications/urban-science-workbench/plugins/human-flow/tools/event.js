// An event's crowd: one of the City's event permits on the district's streets (2014–2018, the data library's record
// of the Hoddle Grid's events), with what the counters on its street counted that day against the same weekday in
// the weeks about it, and a crowd brought to its street: walkers coming from the district's entrances and buildings,
// standing there a while and going back. The walkers run again with the crowd, and its street's regime if one is
// chosen, beside the same hours without it; what they did is written as a record of events, every one simulated. And
// the street's regime decided online as the crowd comes, by POMCP on the walkers themselves, where the crowd's size is
// not known in advance.

import { ask, DAY_NAMES, emit, ensurePrepared, hhmm, mjdOf, runAndWait, setScenario, state, streetOf } from '../shared.js';
import { brief, figures, fold, hint, seg, slider, tabs } from '/src/quiet.js';

const weekdayOf = (date) => { const [y, m, d] = date.split('-').map(Number); return (new Date(Date.UTC(y, m - 1, d)).getUTCDay() + 6) % 7; };

export async function render(root, ctx) {
  const { h, fmt: fmt0 } = ctx.dom;
  // Whole numbers from 10 up (the shell's fmt strips the zeros of a whole number toPrecision gives: 90 → "9").
  const fmt = (x, d = 3) => (typeof x === 'number' && Math.abs(x) >= 10 ? Math.round(x).toLocaleString('en') : fmt0(x, d));
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'An event\'s crowd'), h('span', { class: 'lede' }, 'A permit on a street, its crowd, and the walkers around it.')));
  const body = h('div', {});
  root.append(body);
  body.append(hint(h, 'preparing the district\'s walkers…'));
  let alive = true;
  ctx.onDispose(() => { alive = false; });
  let P, permits;
  try {
    P = await ensurePrepared(ctx);
    permits = (await ask(ctx, 'permits')).permits;
  } catch (e) { return ctx.dom.failure(body, e); }
  if (!alive) return;
  const opt = { permit: null, street: 'Bourke Street Mall', crowd: 2000, from: 12, hours: 3, dwell: 40, regime: 'open', date: null, all: false };
  const out = h('div', {});

  // ---- the permits ----
  const listBox = h('div', {});
  const drawList = () => {
    const shown = permits.filter((p) => opt.all || /^Public Event/.test(p.kind)).slice().sort((a, b) => b.counts.length - a.counts.length || (a.from < b.from ? 1 : -1));
    listBox.replaceChildren(
      h('div', { class: 'ld-row' }, seg(h, [['public', 'public events'], ['all', 'every permit']], opt.all ? 'all' : 'public', (k) => { opt.all = k === 'all'; drawList(); }),
        h('span', { class: 'ld-quiet' }, `${shown.length} on the district's streets`)),
      h('ul', { class: 'ld-list stack', style: 'max-height:300px;overflow:auto' }, shown.slice(0, 80).map((p) => {
        const c = p.counts[0];
        const change = c && c.usual > 0 ? (100 * (c.day - c.usual)) / c.usual : NaN;
        return h('li', { style: `cursor:pointer;${opt.permit === p ? 'background:var(--accent-soft)' : ''}`, onclick: () => choose(p) },
          h('div', { class: 'ld-row' }, h('span', {}, `${p.from} · ${p.kind.replace(/^Public Event - /, '')}`), h('span', { class: 'faint' }, p.streets.join(', '))),
          h('div', { class: 'ld-quiet' }, p.location.length > 110 ? `${p.location.slice(0, 110)}…` : p.location),
          c ? h('div', { class: 'ld-quiet' }, `${c.counter} counted ${Math.round(c.day).toLocaleString('en')} that day, ${Number.isFinite(change) ? `${change >= 0 ? '+' : '−'}${fmt(Math.abs(change), 2)}%` : ''} on its usual ${DAY_NAMES[weekdayOf(p.from)]} (${Math.round(c.usual).toLocaleString('en')})`) : null);
      })));
  };
  const choose = (p) => {
    opt.permit = p, opt.date = p.from;
    const st = p.streets.find((s) => streetOf(s));
    if (st) opt.street = st;
    drawList(), drawSetup();
  };

  // ---- the crowd ----
  const setupBox = h('div', {});
  const drawSetup = () => {
    const weekday = opt.date ? weekdayOf(opt.date) : ctx.clock.weekday;
    const streets = P.streets.filter((s) => s.length >= 60);
    setupBox.replaceChildren(
      h('div', { class: 'ld-label' }, 'the crowd'),
      h('div', { class: 'ld-row' }, h('select', { onchange: (e) => { opt.street = e.target.value; } }, streets.map((s) => h('option', { value: s.name, selected: s.name === opt.street ? true : null }, s.name))),
        h('span', { class: 'ld-quiet' }, `on a ${DAY_NAMES[weekday]}${opt.date ? ` (as ${opt.date})` : ''}`),
        opt.date ? h('button', { class: 'ld-link', onclick: () => ctx.clock.set({ date: opt.date, hour: opt.from + 0.5 }) }, 'set the clock to that day') : null),
      slider(h, fmt, { label: 'walkers an hour', min: 250, max: 8000, step: 250, value: opt.crowd, onChange: (v) => { opt.crowd = v; } }),
      slider(h, fmt, { label: 'staying', min: 5, max: 120, step: 5, value: opt.dwell, unit: ' min', onChange: (v) => { opt.dwell = v; } }),
      slider(h, fmt, { label: 'from', min: 7, max: 21, step: 1, value: opt.from, unit: ':00', onChange: (v) => { opt.from = v; } }),
      slider(h, fmt, { label: 'for', min: 1, max: 4, step: 1, value: opt.hours, unit: ' h', onChange: (v) => { opt.hours = v; } }),
      h('div', { class: 'ld-label' }, 'the street meanwhile'),
      seg(h, [['open', 'open'], ['one-way', 'one-way (phased)'], ['closed', 'closed to walkers passing through']], opt.regime, (k) => { opt.regime = k; }),
      hint(h, 'The permits give dates, not hours or crowds: the crowd\'s size, its hours and how long it stays are yours to set.'),
      h('div', { class: 'ld-row', style: 'margin-top:6px' }, h('button', { class: 'primary', onclick: () => runEvent() }, 'Bring the crowd')));
  };

  const specs = () => {
    const weekday = opt.date ? weekdayOf(opt.date) : ctx.clock.weekday, from = opt.from - 1, to = Math.min(24, opt.from + opt.hours);
    // With the street's regime changed, the run without the crowd keeps the street's own walkers apart alike (open).
    const base = { weekday, from, to, ...(opt.regime !== 'open' ? { regime: { street: opt.street, kind: 'open', from: opt.from, to } } : {}) };
    const event = { street: opt.street, crowd: opt.crowd, from: opt.from, to, dwell: opt.dwell };
    const withCrowd = { ...base, event, ...(opt.regime !== 'open' ? { regime: { street: opt.street, kind: opt.regime, from: opt.from, to } } : {}) };
    return { base, withCrowd, to };
  };

  async function runEvent() {
    const S = specs();
    out.replaceChildren(hint(h, 'running the walkers without the crowd and with it…'));
    try {
      const prog = (what) => (st) => { if (alive && out.firstChild) out.firstChild.textContent = `${what}: ${hhmm(Math.max(S.base.from, st.reached))} reached`; };
      const RA = await runAndWait(ctx, S.base, prog('without the crowd'), () => alive);
      const RB = await runAndWait(ctx, S.withCrowd, prog('with the crowd'), () => alive);
      const C = await ask(ctx, 'compare', { a: S.base, b: S.withCrowd, over: opt.from, until: S.to });
      const date = opt.date || ctx.clock.date;
      const rec = await ask(ctx, 'record', { ...S.withCrowd, mjd0: mjdOf(date, S.base.from, ctx.area.timezone) });
      if (!alive) return;
      const text = `${opt.crowd.toLocaleString('en')} an hour to ${opt.street}, ${hhmm(opt.from)}–${hhmm(S.to)}${opt.regime !== 'open' ? `, the street ${opt.regime === 'closed' ? 'closed to walkers passing through' : 'one-way'}` : ''}`;
      state.compare = { a: C.a, b: C.b, label: `the crowd against the same hours without it` };
      state.record = rec;
      setScenario({ spec: S.withCrowd, label: text });
      emit('compare'), emit('record');
      if (ctx.clock.weekday !== S.base.weekday) ctx.clock.set({ date: date, hour: opt.from + 1 });
      else if (ctx.clock.hour < opt.from || ctx.clock.hour >= S.to) ctx.clock.set({ hour: opt.from + 1 });
      out.replaceChildren(result(C, RA, RB, rec, text, S));
    } catch (e) { ctx.dom.failure(out, e); }
  }

  function result(C, RA, RB, rec, text, S) {
    // The other streams draw alike with the crowd and without it (draws keyed by stream and minute), so the crowd's
    // departures are the difference in departures; with the street's regime changed as well, nearly so.
    const came = Math.max(0, RB.totals.departed - RA.totals.departed);
    // A junction named by the streets that meet there.
    const placeOf = (e) => {
      if (e.where !== 'a junction') return e.where;
      let j = 0, best = Infinity;
      P.junctions.forEach((q, k) => { const d = Math.hypot(q[0] - e.at[0], q[1] - e.at[1]); if (d < best) best = d, j = k; });
      const names = [...new Set(P.segments.filter((g) => (g.u === j || g.v === j) && g.name).map((g) => g.name))];
      return names.length ? names.slice(0, 2).join(' at ') : 'a junction';
    };
    const waiting = (e) => +((e.source.match(/at most (\d+) waiting/) || [])[1] || 0);
    const crowds = rec.events.filter((e) => e.kind === 'crowd');
    const queues = rec.events.filter((e) => e.kind === 'queue' && e.from >= opt.from - 1e-6).sort((a, b) => waiting(b) - waiting(a));
    const longest = crowds.slice().sort((a, b) => b.to - b.from - (a.to - a.from))[0];
    const T = tabs(h, [
      { name: 'Around it', panel: () => h('div', {},
        figures(h, [{ value: came.toLocaleString('en'), label: `came for the event${opt.regime !== 'open' ? ' (about)' : ''}` },
          { value: `+${fmt(C.windowB.personHours - C.windowA.personHours, 3)}`, label: 'person-hours on the streets' },
          { value: `${fmt(C.windowA.crowdedHours, 3)} → ${fmt(C.windowB.crowdedHours, 3)}`, label: 'person-hours in crowds (Fruin\'s E)' },
          { value: `${fmt(60 * (RB.totals.delayMinutes - RA.totals.delayMinutes), 2)} s`, label: 'more waiting a trip' }]),
        brief(ctx, { kind: 'model', title: text,
          statement: 'The crowd\'s walkers are drawn in rose, standing on the street while they stay; beneath them each street by how many more walk it than without the event. They come from the district\'s entrances and buildings as the visits do, walk the logit\'s ways there and back, and may enter the street even when it is closed to walkers passing through.',
          method: 'Two runs with the same walkers\' draws; the first hour of each is its warm-up. The crowd is a stream of its own: Poisson departures at the rate set over the hours set, an exponential stay of the mean set.',
          limits: 'The permit says when and where, not how many: the crowd is declared. A crowd standing counts in the street\'s density as a crowd does; nobody avoids the district because of the event.' })) },
      { name: 'The record', panel: () => h('div', {},
        figures(h, [{ value: String(crowds.length), label: 'crowds (ten walkers or more at level E for a minute)' }, { value: String(queues.length), label: 'queues of twenty or more at a street\'s mouth, in the event\'s hours' }]),
        h('p', { class: 'note' }, `Every event of this record is simulated (${rec.simulated ? 'checked' : 'not all'}), and is drawn on the map in pale violet glass, never as an observed event.${rec.problems.length ? ` Its problems: ${rec.problems.slice(0, 2).join('; ')}` : ' The record is sound.'}`),
        longest ? hint(h, `the longest crowd: on ${longest.where}, ${hhmm(longest.from)}–${hhmm(longest.to)}`) : null,
        h('ul', { class: 'ld-list', style: 'max-height:260px;overflow:auto' }, [...crowds, ...queues].slice(0, 60).map((e) =>
          h('li', { style: 'cursor:pointer', onclick: () => ctx.gis.flyTo(e.at, 220) },
            h('span', {}, e.kind === 'crowd' ? `a crowd on ${placeOf(e)}, ${hhmm(e.from)}–${hhmm(e.to)}` : `${waiting(e)} waiting at ${placeOf(e)}, ${hhmm(e.from)}–${hhmm(e.to)}`),
            h('span', { class: 'faint' }, 'simulated'))))) },
      { name: 'Its regime, online', panel: () => advicePanel(S) },
    ], { remember: 'hf-event' });
    return T.el;
  }

  // POMCP's advice: the street's regime every ten minutes as the crowd comes, its size hidden.
  function advicePanel(S) {
    const box = h('div', {});
    const params = { epochs: 3, simulations: 12, truth: 1.5 };
    const go = async () => {
      box.replaceChildren(hint(h, 'planning on the walkers…'));
      const spec = { ...S.base, event: S.withCrowd.event, regime: { street: opt.street, kind: 'open', from: opt.from, to: S.to } };
      try {
        let A = await ask(ctx, 'advice', { ...spec, ...params, epochMinutes: 10 });
        while (A.pending && alive) {
          box.firstChild.textContent = `planning on the walkers: ${A.epoch} of ${params.epochs} decisions made…`;
          await new Promise((r) => setTimeout(r, 2000));
          A = await ask(ctx, 'advice', { ...spec, ...params, epochMinutes: 10 });
        }
        if (!alive) return;
        const slow = Math.max(...A.steps.map((s) => s.planningSeconds));
        const best = A.fixed.slice().sort((a, b) => b.reward - a.reward)[0];
        box.replaceChildren(
          figures(h, [{ value: A.steps.map((s) => s.action).join(' → '), label: 'the regimes chosen, ten minutes each' },
            { value: `${fmt(slow, 2)} s`, label: 'to plan a decision, at most' },
            { value: `${fmt(600 / slow, 2)}×`, label: 'faster than the street\'s time' }]),
          h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'from'), h('th', {}, 'chosen'), h('th', {}, 'the count'), h('th', {}, 'crowd ½ · 1 · 2, believed after')),
            ...A.steps.map((s) => h('tr', {}, h('td', {}, hhmm(s.hour)), h('td', {}, s.action), h('td', {}, ['low', 'as expected', 'high'][s.observation] || String(s.observation)),
              h('td', {}, s.belief.map((b) => fmt(b, 2)).join(' · '))))),
          brief(ctx, { kind: 'model', title: `POMCP chose ${A.steps[0].action} first; in hindsight ${best.regime} held throughout did best`,
            statement: `Planning on ${A.simulations} simulations of the walkers a decision, it is fast enough to advise as the crowd comes (${fmt(slow, 2)} s for ten minutes of street). Held throughout, the regimes would have cost ${A.fixed.map((f) => `${f.regime} ${Math.round(-f.reward).toLocaleString('en')}`).join(', ')} (person-hours, crowds weighted ten times); POMCP's sequence ${Math.round(-A.reward).toLocaleString('en')}.`,
            method: `tyche's POMCP on the movement module's RegimeModel: the hidden state is every walker and the crowd's size (½, 1 or 2 times what was declared, equally likely); it observes ${A.sensed.join(' and ')} each epoch (low, as expected or high); the reward is minus the person-hours present and ten times those above Fruin's level F. The world's crowd is ${A.truth} times what was declared. Each particle is a copy of the walkers, whose streams keep the pairs carrying 90% of their weight.`,
            limits: `${A.simulations} simulations a decision is a small budget: the choice is an estimate, not an optimum, and the regimes' costs differ by little here. Planning and acting took ${fmt(A.seconds, 3)} s in all.` }));
      } catch (e) { ctx.dom.failure(box, e); }
    };
    return h('div', {}, hint(h, 'The street\'s regime decided every ten minutes as the crowd arrives, when nobody knows yet how large it will be.'),
      slider(h, fmt, { label: 'decisions', min: 1, max: 6, step: 1, value: params.epochs, onChange: (v) => { params.epochs = v; } }),
      slider(h, fmt, { label: 'simulations each', min: 4, max: 48, step: 4, value: params.simulations, onChange: (v) => { params.simulations = v; } }),
      slider(h, fmt, { label: 'the crowd, truly', min: 0.5, max: 2, step: 0.25, value: params.truth, unit: '×', onChange: (v) => { params.truth = v; } }),
      h('div', { class: 'ld-row', style: 'margin:6px 0' }, h('button', { class: 'primary', onclick: go }, 'Ask for the regime')), box);
  }

  drawList(), drawSetup();
  body.replaceChildren(h('div', { class: 'ld-label' }, 'the City\'s event permits on the district\'s streets, 2014–2018'), listBox, setupBox, out,
    fold(h, 'where the permits come from', h('p', {}, 'City of Melbourne, event permits 2014–2018 (CC BY 4.0), matched to the district\'s streets in the data library\'s record of the Hoddle Grid\'s events; the counts are the City\'s Pedestrian Counting System\'s, that day against the median of the same weekday in the four weeks either side.')));
}
