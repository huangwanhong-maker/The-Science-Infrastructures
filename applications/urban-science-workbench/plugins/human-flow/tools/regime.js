// A street closed or made one-way: the walkers run again over a few hours with the street's regime changed, and
// beside the same hours as built (the same walkers' draws, so the two runs differ only by the regime), how they
// redistribute over the district: drawn on the map street by street, and summed here. Closed means closed to walkers
// passing through (those bound for the street still reach it); one-way means walked one way only, from the end named,
// as crowds are phased at events. A street's stubs where it leaves the district stay open: walkers enter there.

import { ask, DAY_NAMES, emit, ensurePrepared, hhmm, runAndWait, setScenario, state, streetOf } from '../shared.js';
import { brief, figures, hint, seg, slider } from '/src/quiet.js';

const KINDS = [['closed', 'closed to walkers passing through'], ['one-way', 'one-way (phased)'], ['open', 'open']];

// A compass name of a junction's place on a street: its end farthest that way.
function endName(P, street, j) {
  const xs = street.ends.map((e) => P.junctions[e]);
  const p = P.junctions[j], cx = xs.reduce((a, q) => a + q[0], 0) / xs.length, cy = xs.reduce((a, q) => a + q[1], 0) / xs.length;
  const dx = p[0] - cx, dy = p[1] - cy;
  return Math.abs(dx) > Math.abs(dy) ? (dx > 0 ? 'the east end' : 'the west end') : (dy > 0 ? 'the north end' : 'the south end');
}

export async function render(root, ctx) {
  const { h, fmt: fmt0 } = ctx.dom;
  // Whole numbers from 10 up (the shell's fmt strips the zeros of a whole number toPrecision gives: 90 → "9").
  const fmt = (x, d = 3) => (typeof x === 'number' && Math.abs(x) >= 10 ? Math.round(x).toLocaleString('en') : fmt0(x, d));
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'A street closed or one-way'), h('span', { class: 'lede' }, 'Change a street; the walkers find other ways.')));
  const body = h('div', {});
  root.append(body);
  body.append(hint(h, 'preparing the district\'s walkers…'));
  let alive = true;
  ctx.onDispose(() => { alive = false; });
  let P;
  try { P = await ensurePrepared(ctx); } catch (e) { return ctx.dom.failure(body, e); }
  const streets = P.streets.filter((s) => s.length >= 80 && s.edges.length >= 2);
  const sel = ctx.selection && ctx.selection.kind === 'segment' ? ctx.selection : null;
  const opt = {
    street: sel && sel.name && streetOf(sel.name) ? sel.name : 'Swanston Street', block: sel ? sel.edge : null, useBlock: !!sel, kind: 'closed', entry: null,
    from: Math.min(21, Math.max(7, Math.floor(ctx.clock.hour))), hours: 3,
  };
  const out = h('div', {});
  const controls = h('div', {});
  const draw = () => {
    const st = streetOf(opt.street);
    const select = h('select', { onchange: (e) => { opt.street = e.target.value; opt.useBlock = false; opt.entry = null; draw(); } },
      streets.map((s) => h('option', { value: s.name, selected: s.name === opt.street && !opt.useBlock ? true : null }, `${s.name} (${Math.round(s.length)} m)`)));
    const where = opt.block !== null
      ? seg(h, [['street', 'the whole street'], ['block', 'the block picked']], opt.useBlock ? 'block' : 'street', (k) => { opt.useBlock = k === 'block'; draw(); })
      : hint(h, 'or pick a block of street on the map (the walkers layer), and choose it here');
    const ends = st ? st.ends : [];
    controls.replaceChildren(...[
      h('div', { class: 'ld-label' }, 'the street'), h('div', { class: 'ld-row' }, select), where,
      h('div', { class: 'ld-label' }, 'its regime'), seg(h, KINDS, opt.kind, (k) => { opt.kind = k; draw(); }),
      opt.kind === 'one-way' && !opt.useBlock && ends.length ? h('div', { class: 'ld-row' }, h('span', { class: 'ld-quiet' }, 'walked from'),
        seg(h, ends.slice(0, 4).map((j) => [String(j), endName(P, st, j)]), String(opt.entry ?? ends[0]), (k) => { opt.entry = +k; })) : null,
      slider(h, fmt, { label: 'from', min: 6, max: 22, step: 1, value: opt.from, unit: ':00', onChange: (v) => { opt.from = v; } }),
      slider(h, fmt, { label: 'for', min: 1, max: 4, step: 1, value: opt.hours, unit: ' h', onChange: (v) => { opt.hours = v; } }),
      h('div', { class: 'ld-row', style: 'margin-top:8px' }, h('button', { class: 'primary', onclick: () => compareOne() }, 'Run the walkers'),
        h('button', { class: 'plain', onclick: () => compareAll() }, 'all three regimes'))].filter(Boolean));
  };

  const specs = () => {
    const weekday = ctx.clock.weekday, from = opt.from - 1, to = Math.min(24, opt.from + opt.hours);
    const regime = (kind) => ({ ...(opt.useBlock && opt.block !== null ? { edges: [opt.block] } : { street: opt.street }), kind, from: opt.from, to, ...(kind === 'one-way' && opt.entry !== null && !opt.useBlock ? { entry: opt.entry } : {}) });
    // As built, the street open, its own walkers kept apart as under the change, so that both runs draw alike.
    const base = { weekday, from, to, regime: regime('open') };
    return { base, with: (kind) => (kind === 'open' ? base : { ...base, regime: regime(kind) }), where: opt.useBlock ? `the block of ${opt.street} picked` : opt.street, to };
  };
  const label = (S, kind) => `${S.where} ${kind === 'closed' ? 'closed to walkers passing through' : kind === 'one-way' ? 'one-way' : 'open'}, ${hhmm(opt.from)}–${hhmm(S.to)}`;

  // Both runs, then the comparison drawn on the map and summed here.
  async function compareOne() {
    const S = specs(), b = S.with(opt.kind);
    out.replaceChildren(hint(h, 'running the walkers as built and with the change…'));
    try {
      const prog = (what) => (st) => { if (alive) out.firstChild.textContent = `${what}: ${hhmm(Math.max(S.base.from, st.reached))} reached`; };
      const RA = await runAndWait(ctx, S.base, prog('as built'), () => alive);
      const RB = await runAndWait(ctx, b, prog('with the change'), () => alive);
      const C = await ask(ctx, 'compare', { a: S.base, b, over: opt.from, until: S.to });
      if (!alive) return;
      const text = label(S, opt.kind);
      state.compare = { a: C.a, b: C.b, label: `${text} against as built` };
      setScenario({ spec: b, label: text });
      emit('compare');
      if (ctx.clock.hour < opt.from || ctx.clock.hour >= S.to || ctx.clock.weekday !== S.base.weekday) ctx.clock.set({ hour: opt.from + 0.5 });
      out.replaceChildren(summary(C, RA, RB, text));
    } catch (e) { ctx.dom.failure(out, e); }
  }

  function byStreet(C) {
    const m = new Map();
    P.segments.forEach((s, k) => {
      const n = s.name || 'unnamed ways';
      const e = m.get(n) || { name: n, a: 0, b: 0 };
      e.a += C.a[k], e.b += C.b[k];
      m.set(n, e);
    });
    return [...m.values()].map((e) => ({ ...e, d: e.b - e.a })).sort((x, y) => x.d - y.d);
  }

  function summary(C, RA, RB, text) {
    const WA = C.windowA, WB = C.windowB, d = WB.personHours - WA.personHours;
    const streetsD = byStreet(C), lost = streetsD.slice(0, 4), gained = streetsD.slice(-4).reverse();
    const row = (e) => h('li', {}, h('span', {}, e.name), h('span', { class: 'faint' }, `${e.d > 0 ? '+' : '−'}${fmt(Math.abs(e.d), 3)} walkers at once (${fmt(e.a, 3)} → ${fmt(e.b, 3)})`));
    const failed = RB.totals.failed - RA.totals.failed;
    return h('div', {},
      figures(h, [{ value: `${d >= 0 ? '+' : '−'}${Math.round(Math.abs(d)).toLocaleString('en')}`, label: 'person-hours on the streets over those hours' },
        { value: `${d >= 0 ? '+' : '−'}${fmt(Math.abs(100 * d / WA.personHours), 2)}%`, label: 'of them' },
        { value: `${fmt(60 * (RB.totals.delayMinutes - RA.totals.delayMinutes), 2)} s`, label: 'more waiting a trip, over the run' },
        { value: failed.toLocaleString('en'), label: 'walkers left with no way' }]),
      h('div', { class: 'ld-label' }, 'fewer walkers'), h('ul', { class: 'ld-list' }, lost.map(row)),
      h('div', { class: 'ld-label' }, 'more walkers'), h('ul', { class: 'ld-list' }, gained.map(row)),
      brief(ctx, { kind: 'model', title: text,
        statement: `The map draws the walkers with the change, and beneath them each street by how many more (brick) or fewer (sea) walk it than as built over those hours. Walkers in crowds at Fruin's level E: ${fmt(WA.crowdedHours, 3)} → ${fmt(WB.crowdedHours, 3)} person-hours.`,
        method: 'The same walkers\' draws in both runs (common random numbers): a walker departs, chooses and walks as in the other run until the regime changes the ways open to it, then chooses again from its next junction. The first hour of each run is its warm-up and is not compared.',
        limits: 'One realisation each: differences of a few walkers on a street are the runs\' noise. The regime changes routes, not demand: nobody stays away because a street is closed.' }));
  }

  async function compareAll() {
    const S = specs();
    out.replaceChildren(hint(h, 'running the walkers under each regime…'));
    try {
      const rows = [];
      for (const [kind] of KINDS) {
        const spec = kind === 'open' ? S.base : S.with(kind);
        const R = await runAndWait(ctx, spec, (st) => { if (alive) out.firstChild.textContent = `${kind}: ${hhmm(Math.max(S.base.from, st.reached))} reached`; }, () => alive);
        const C = await ask(ctx, 'compare', { a: S.base, b: spec, over: opt.from, until: S.to });
        rows.push({ kind, R, W: C.windowB });
      }
      if (!alive) return;
      const base = rows.find((r) => r.kind === 'open');
      out.replaceChildren(h('table', { class: 'data' },
        h('tr', {}, h('th', {}, ''), h('th', {}, 'person-hours'), h('th', {}, 'in crowds (E)'), h('th', {}, 'waiting a trip'), h('th', {}, 'no way')),
        ...rows.map((r) => h('tr', {}, h('td', {}, KINDS.find((k) => k[0] === r.kind)[1]), h('td', {}, Math.round(r.W.personHours).toLocaleString('en')),
          h('td', {}, fmt(r.W.crowdedHours, 3)), h('td', {}, `${fmt(60 * r.R.totals.delayMinutes, 3)} s`), h('td', {}, String(r.R.totals.failed - base.R.totals.failed))))),
      hint(h, `${S.where}, ${hhmm(opt.from)}–${hhmm(S.to)} on a ${DAY_NAMES[S.base.weekday]}; person-hours and crowds over those hours, waiting over the run.`));
    } catch (e) { ctx.dom.failure(out, e); }
  }

  ctx.onSelect((s) => { if (s && s.kind === 'segment') { opt.block = s.edge; opt.useBlock = true; if (s.name && streetOf(s.name)) opt.street = s.name; draw(); } });
  draw();
  body.replaceChildren(controls, out);
}
