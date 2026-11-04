// Opening hours: where and when a street is open, from the hours its holdings declare; the rhythm of those
// hours over the week (a shop closing at five and a bar at two are different frequencies); whether the street
// opens and closes in one order; and which kinds of people its open shops draw at each hour (a model).

const DAYS = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'];
function profile(f) { const p = []; for (let k = 0; k < 168; k++) p.push(f(Math.floor(k / 24), (k % 24) + 0.5)); return p; }
const bump = (x, mu, s) => Math.exp(-0.5 * ((x - mu) / s) ** 2);
const TYPES = [
  { name: 'office workers', presence: profile((d, t) => (d < 5 ? 600 * bump(t, 12.7, 0.9) + 250 * bump(t, 18, 1.2) : 40 * bump(t, 13, 3))),
    affinity: { 'food-and-drink': 2.0, 'personal-services': 0.5, 'other-retail': 0.3 }, distanceCost: 0.008 },
  { name: 'couples planning a wedding', presence: profile((d, t) => (d >= 5 ? 60 * bump(t, 14.5, 2.2) : 12 * bump(t, 17.5, 1.6))),
    affinity: { 'bridal-and-wedding': 3.0, jewellery: 2.6, 'clothing-and-fashion': 0.8, 'food-and-drink': 0.6 }, distanceCost: 0.002 },
  { name: 'evening visitors', presence: profile((d, t) => (d >= 4 ? 300 : 140) * bump(t, 20.5, 1.8)),
    affinity: { 'food-and-drink': 2.4, 'entertainment-and-culture': 2.0 }, distanceCost: 0.004 },
];

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, raster, ramp, legend, category } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Opening hours'),
    h('span', { class: 'lede' }, 'A shop stays open while staying open pays: its declared hours are its own reading of the street. Pick a shop on the map, or choose a street.')));
  let streets;
  try { streets = (await ctx.ask('streets')).filter((s) => s.declared > 0 && s.area === ctx.area.id); } catch (e) { return failure(root, e); }
  if (!streets.length) return failure(root, new Error('no street with declared hours here'));
  const sel = h('select', {}, streets.map((s) => h('option', { value: s.id }, `${s.name} (${s.declared} of ${s.holdings} declare hours)`)));
  const day = h('select', {}, DAYS.map((d, i) => h('option', { value: i }, d)));
  day.value = String(ctx.clock.weekday);
  root.append(h('div', { class: 'keys' }, h('label', { class: 'pick' }, sel), h('label', { class: 'pick' }, 'order of opening on ', day)));
  const out = h('div', {});
  root.append(out);
  const choose = (s) => { if (s && s.kind === 'holding' && streets.some((x) => x.id === s.street)) { sel.value = s.street; return true; } return false; };
  choose(ctx.selection);
  ctx.onSelect((s) => { if (choose(s)) run(); });

  const run = async () => {
    busy(out);
    let d, v;
    try { [d, v] = await Promise.all([ctx.ask('street', { street: sel.value }), ctx.ask('visits', { street: sel.value, types: TYPES })]); } catch (e) { return failure(out, e); }
    const of = d.openFraction;
    const hx = d.forcing.map((_, k) => k + 0.5);
    const top = Object.entries(d.forcingByUse).sort((a, b) => b[1].reduce((s, x) => s + x, 0) - a[1].reduce((s, x) => s + x, 0)).slice(0, 5);
    const marks = [24, 48, 72, 96, 120, 144].map((x, i) => ({ x, label: DAYS[i + 1] }));
    const forcing = lines([{ name: 'all', x: hx, y: d.forcing, color: css('--ink-2'), width: 2 }, ...top.map(([u, f]) => ({ name: u, x: hx, y: f, color: category(d.uses.indexOf(u)), width: 1.3 }))],
      { xlabel: 'hour of the week', ylabel: 'open', height: 210, xlim: [0, 168], marks });
    const power = d.lines.filter((l) => l.m > 0).reduce((a, l) => { a[l.daily ? 0 : 1] += l.amplitude ** 2; return a; }, [0, 0]);
    const lineBars = bars(d.lines.filter((l) => l.m > 0).map((l) => ({ x: l.m, y: l.amplitude, color: l.daily ? css('--accent') : css('--minor'),
      label: `${l.m} a week (${fmt(l.frequency, 3)} a day): ${fmt(l.amplitude, 3)}${l.daily ? ', a harmonic of the day' : ', how days differ'}` })),
    { xlabel: 'cycles a week (multiples of 7: the day)', ylabel: 'amplitude', height: 180 });
    const a = d.ascentDescent[+day.value];
    const order = a.holdings >= 3 ? lines([], { xlabel: 'opens (h)', ylabel: 'closes (h)', height: 220, width: 380,
      xlim: [Math.min(...a.open) - 0.5, Math.max(...a.open) + 0.5], ylim: [Math.min(...a.close) - 0.5, Math.max(...a.close) + 0.5],
      points: a.open.map((o, i) => ({ x: o, y: a.close[i], label: `opens ${fmt(o, 3)}, closes ${fmt(a.close[i], 3)}`, color: css('--accent') })) }) : h('p', { class: 'faint' }, 'Too few declare that day.');
    const tc = [css('--series-1'), css('--series-2'), css('--series-4')];
    const ranked = d.holdings.map((x, i) => ({ x, i, total: v.holdings[i].totals.reduce((s, t) => s + t, 0) })).filter((r) => r.total > 0).sort((p, q) => q.total - p.total).slice(0, 10);
    out.replaceChildren(
      h('h3', {}, `${d.name}: open by place and hour`), raster(of.values, { map: 'orchid', lo: 0, hi: 1, height: Math.min(300, 40 + 2 * of.values.length) }),
      h('div', { style: 'display:grid;grid-template-columns:repeat(7,1fr);font-size:0.8rem;color:var(--ink-faint)' }, DAYS.map((x) => h('span', {}, x))),
      ramp('orchid', 0, 1, 'share open'), h('p', { class: 'caption' }, `Rows: along the street, 0 to ${fmt(d.length, 4)} m (top to bottom); columns: half-hours of the week. ${d.declared} of ${d.holdings.length} holdings declare hours.`),
      h('div', { class: 'section' }, h('h3', {}, 'Holdings open through the week'), forcing, legend([{ name: 'all', color: css('--ink-2') }, ...top.map(([u]) => ({ name: u, color: category(d.uses.indexOf(u)) }))]),
        lineBars, h('p', { class: 'caption' }, `${fmt(100 * power[1] / (power[0] + power[1]), 3)}% of the rhythm is how days differ (weekends, late nights); the rest repeats every day.`)),
      h('div', { class: 'section' }, h('h3', {}, `${DAYS[+day.value]}: who opens first, and who closes first`), order,
        ctx.finding({ title: a.holdings >= 3 ? (a.p < 0.05 ? (a.tau > 0 ? 'First to open, first to close' : 'First to open, last to close') : 'Opening and closing in unrelated orders') : '—', kind: 'measured',
          statement: a.holdings >= 3 ? `Kendall's τ = ${fmt(a.tau, 2)} (p = ${fmt(a.p, 2)}); openings spread over ${fmt(a.sdOpen, 2)} h and closings over ${fmt(a.sdClose, 2)} h. A street whose morning and evening are not mirror images moves between states rather than following the light.` : '',
          method: 'The declared hours of that day; p by shuffling closings among holdings.', limits: 'Declared, not observed.' })),
      h('div', { class: 'section' }, h('h3', {}, 'Which kinds of people its open shops draw'),
        lines(v.weekly.map((w, k) => ({ name: v.types[k], x: hx, y: w, color: tc[k] })), { xlabel: 'hour of the week', ylabel: 'visits an hour', height: 200, xlim: [0, 168], marks }),
        legend(v.types.map((t, k) => ({ name: t, color: tc[k] }))),
        h('table', { class: 'data' }, h('tr', {}, h('th', {}, 'holding'), h('th', {}, 'use'), ...v.types.map((t) => h('th', {}, t))),
          ranked.map((r) => h('tr', {}, h('td', {}, r.x.name || '(unnamed)'), h('td', { class: 'left' }, d.uses[r.x.use]), ...v.holdings[r.i].share.map((s) => h('td', {}, `${fmt(100 * s, 3)}%`))))),
        ctx.finding({ title: 'One street, a different street to each kind of visitor', kind: 'model',
          statement: `Office workers at lunch and after work, couples planning a wedding at weekends, evening visitors: each chooses among the shops open when they are there, by what they want and how far they must walk. What is open when a kind of person is present is what that person can reach. ${v.unserved.map((u, k) => `${fmt(u, 4)} ${v.types[k]}`).join(', ')} a week find nothing open.`,
          method: 'A choice model (logit) over open holdings; who is present when, and what they want, are declared.', limits: 'Nobody\'s actual visits are measured here.' })));
  };
  sel.addEventListener('change', run);
  day.addEventListener('change', run);
  run();
}
