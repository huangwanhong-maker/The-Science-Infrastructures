// The rhythm of a street: a sensor's hourly counts read on the calendars — the week (with the day as its
// harmonics), the synodic month, the year. Its spectrum against red noise (what slow, irregular variation alone
// would give), its lines classified (a calendar's own, a calendar modulated by the season, or the street's own),
// the harmonic constants of a standard list of constituents, how selectively it answers the calendars, whose
// clock it keeps when daylight saving begins and ends, and the coupling left in what the calendars do not explain.

const DAYS = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'];
const VERDICT = { forced: ['calendar', 'forced'], combination: ['calendar × season', 'combination'], undecided: ['calendar × season', 'undecided'], own: ['the street\'s own', 'own'] };

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, raster, ramp, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Rhythm of a street'),
    h('span', { class: 'lede' }, 'A street varies with the hour, the day, the week and the year. Pick a sensor on the map, or choose one here.')));
  let sensors = [];
  try { sensors = (await ctx.ask('sensors')).filter((s) => s.locationType === 'Outdoor').sort((a, b) => a.name.localeCompare(b.name)); }
  catch (e) { return failure(root, e); }
  const sel = h('select', {}, sensors.map((s) => h('option', { value: s.id }, `${s.name} (${s.id})`)));
  const era = h('select', {}, h('option', { value: 'recent' }, 'Oct 2024 – Oct 2026'), h('option', { value: 'history' }, '2009 – 2022'));
  const order = h('select', {}, [0, 1, 2, 3].map((d) => h('option', { value: d, selected: d === 2 }, `${d}`)));
  const range = h('select', {}, h('option', { value: 7.2 }, '0 – 7 a day'), h('option', { value: 1.25 }, '0 – 1.25 a day'), h('option', { value: 0.3 }, '0 – 0.3 a day'));
  root.append(h('div', { class: 'keys' }, h('label', { class: 'pick' }, 'sensor ', sel), h('label', { class: 'pick' }, 'years ', era),
    h('label', { class: 'pick', title: 'how many steps of month and year a line may combine and still count as the calendars\'' }, 'seasonal steps ≤ ', order), h('label', { class: 'pick' }, 'show ', range)));
  const out = h('div', {});
  root.append(out);
  const fromSel = ctx.selection && ctx.selection.kind === 'sensor' ? String(ctx.selection.id) : null;
  sel.value = fromSel || String((sensors.find((s) => /Bourke Street Mall \(North\)/.test(s.name)) || sensors[0] || {}).id);
  ctx.onSelect((s) => { if (s && s.kind === 'sensor') { sel.value = String(s.id); run(); } });

  const bicBox = h('div', {});
  const run = async () => {
    busy(out, 'computing the spectrum and the harmonic constants…');
    let d, clock;
    try {
      [d, clock] = await Promise.all([ctx.ask('rhythm', { sensor: +sel.value, era: era.value, order: +order.value }),
        ctx.ask('clock', { sensor: +sel.value }).catch((e) => ({ error: e.message }))]);
    } catch (e) { return failure(out, e); }
    const fmax = +range.value;
    const sp = d.spectrum;
    const keep = sp.frequency.map((_, i) => i).filter((i) => sp.frequency[i] <= fmax && sp.frequency[i] > 0);
    const at = (a) => keep.map((i) => a[i]);
    const f = at(sp.frequency);
    const marks = [];
    for (let k = 1; k <= 7; k++) if (k <= fmax) marks.push({ x: k, label: `${k}/day`, color: css('--accent') });
    if (fmax < 2) for (let m = 1; m <= 8; m++) if (m / 7 <= fmax && m % 7) marks.push({ x: m / 7, label: m <= 3 ? `${m}/week` : '', dash: '2 4' });
    if (fmax < 0.5) marks.push({ x: 1 / 29.530589, label: 'lunar month', color: css('--minor') });
    const spec = lines([
      { x: f, y: at(sp.lower), y2: at(sp.upper), color: css('--series-1'), opacity: 0.14 },
      { name: 'spectrum', x: f, y: at(sp.density), color: css('--series-1'), width: 1.2 },
      { name: 'red noise', x: f, y: at(sp.background), color: css('--series-3'), dash: '5 4' },
      { name: '95% level of red noise', x: f, y: at(sp.level), color: css('--minor') }],
    { logy: true, xlabel: 'cycles per day', ylabel: 'spectral density of ln(1 + count)', height: 300, marks, xlim: [0, fmax] });
    const peaks = h('table', { class: 'data' }, h('tr', {}, ['frequency', 'period', '× red noise', 'line', 'kind'].map((x) => h('th', {}, x))),
      d.peaks.slice(0, 20).map((p) => h('tr', {}, h('td', {}, `${fmt(p.cpd, 4)}/day`), h('td', {}, p.cpd >= 1 ? `${fmt(24 / p.cpd, 3)} h` : `${fmt(1 / p.cpd, 4)} d`),
        h('td', {}, fmt(p.ratio, 3)), h('td', { class: 'left' }, p.name || '—'), h('td', { class: 'left' }, h('span', { class: `chip ${VERDICT[p.verdict][1]}` }, VERDICT[p.verdict][0])))));
    const own = d.peaks.filter((p) => p.verdict === 'own');
    const hm = d.harmonics;
    const harm = h('table', { class: 'data' }, h('tr', {}, ['constituent', 'calendar', 'period', 'amplitude', '±', 'phase °'].map((x) => h('th', {}, x))),
      hm.terms.slice(0, 14).map((t) => h('tr', {}, h('td', {}, t.name), h('td', { class: 'left' }, t.driving), h('td', {}, t.cpd >= 1 ? `${fmt(24 / t.cpd, 4)} h` : `${fmt(1 / t.cpd, 4)} d`),
        h('td', {}, fmt(t.amplitude, 3)), h('td', {}, fmt(t.amplitudeSe, 2)), h('td', {}, fmt((t.phase * 180) / Math.PI, 3)))));
    const S = d.selectivity;
    const shares = bars(S.drivings.map((x, i) => ({ x: i, y: x.share, color: i === 0 ? css('--accent') : css('--minor'), label: `${x.driving}: ${fmt(100 * x.share, 3)}%` })),
      { labels: S.drivings.map((x) => x.driving), ylabel: 'share of the calendars\' variance', height: 190, width: 300 });
    const prof = lines([{ name: 'mean ln(1 + count)', x: d.profile.map((_, k) => k + 0.5), y: d.profile, color: css('--series-1') }],
      { xlabel: 'hour of the week (Monday 00:00 = 0)', ylabel: 'mean ln(1 + count)', height: 180, xlim: [0, 168], marks: [24, 48, 72, 96, 120, 144].map((x, i) => ({ x, label: DAYS[i + 1] })) });
    const record = lines([{ name: 'daily total', x: d.daily.day.map((x) => 1970 + x / 365.2425), y: d.daily.total, color: css('--series-2'), width: 0.9 }],
      { xlabel: 'year', ylabel: 'pedestrians a day', height: 170 });
    let clockBlock = h('p', { class: 'faint' }, clock.error || '');
    if (!clock.error) {
      const cs = clock.shifts, yr = (c) => +c.date.slice(0, 4) + (+c.date.slice(5, 7) - 0.5) / 12;
      clockBlock = h('div', {}, lines([{ name: 'civil share', x: cs.map(yr), y: cs.map((c) => c.civilShare), color: css('--accent') }],
        { xlabel: 'change of daylight saving', ylabel: 'share following the clock', height: 180, ylim: [-0.5, 1.5],
          points: cs.map((c) => ({ x: yr(c), y: c.civilShare, color: c.change > 0 ? css('--accent') : css('--minor'),
            label: `${c.date}: clocks ${c.change > 0 ? 'forward' : 'back'}; the daily peak moved ${fmt(c.shift, 3)} ± ${fmt(c.se, 2)} h in UTC` })) }),
        ctx.finding({ title: `It keeps the clock: ${fmt(clock.civilShare, 3)} ± ${fmt(clock.civilShareSe, 2)}`, kind: 'measured',
          statement: 'When daylight saving begins or ends, the civil clock moves an hour against the sun. If the street\'s day moves with the clock (1), work and timetables set it; if it stays with the sun (0), daylight does.',
          method: 'The UTC time of the daily peak, fitted on weekdays three weeks before and after each change since 2009 (holidays left out), weighted over the changes.',
          limits: 'Seasonal drift within the six weeks is not removed.' }));
    }
    out.replaceChildren(
      h('h3', {}, d.name), record, prof,
      h('div', { class: 'section' }, h('h2', {}, 'Spectrum'), spec,
        legend([{ name: 'spectrum and its 95% interval', color: css('--series-1') }, { name: 'red noise', color: css('--series-3') }, { name: '95% level', color: css('--minor') }]),
        h('p', { class: 'caption' }, `${d.span.from} to ${d.span.to} (${fmt(d.span.hours / 24, 4)} days); multitaper, ${d.tapers} tapers, resolution ${fmt(d.bandwidthCpd, 3)}/day.`)),
      h('div', { class: 'section' }, h('h3', {}, 'Its lines'), peaks,
        ctx.finding({ title: own.length ? `${own.length} line${own.length > 1 ? 's' : ''} the calendars do not explain` : 'Every strong line belongs to the calendars', kind: 'measured',
          statement: `The week's lines (the day and its harmonics, and how weekends and weeknights differ) are the calendar's. A line combining the week with the year is the week changing with the season — which the street\'s response or the hours themselves (longer in summer) can produce. ${own.length ? `${own.map((p) => `${fmt(p.cpd, 4)}/day`).join(', ')} lie off every combination with at most ${d.order} seasonal steps` : ''}. With that many steps allowed, ${fmt(100 * d.chance, 3)}% of the band would be "explained" by chance.`,
          method: 'Multitaper spectrum of ln(1 + count); lines are local maxima above the 95% level of red noise fitted to the record.',
          against: 'Red noise: slow irregular variation with no rhythm of its own.', limits: 'Counts only; nothing faster than the hour.' })),
      h('div', { class: 'section' }, h('h3', {}, 'Harmonic constants'), harm,
        h('p', { class: 'caption' }, `${hm.constituents} constituents the record can tell apart, explaining ${fmt(100 * hm.varianceExplained, 3)}% of the variance; phases from the record's start (UTC).`)),
      h('div', { class: 'section' }, h('h3', {}, 'How selectively it answers the calendars'), shares,
        ctx.finding({ title: `Selectivity ${fmt(S.S, 3)}`, kind: 'measured',
          statement: `0 when the day, the week and the year share the street's rhythm equally, 1 when one carries it all. The calendars explain ${fmt(100 * S.forced, 3)}% of the variance.`,
          method: 'The variance of the fitted constituents, summed by calendar.' })),
      h('div', { class: 'section' }, h('h2', {}, 'Whose clock does it keep?'), clockBlock),
      h('div', { class: 'section' }, h('h2', {}, 'Coupling in what is left'), h('button', { class: 'plain', onclick: bic }, 'compute the bicoherence (a minute)'), bicBox));
  };

  const bic = async () => {
    busy(bicBox, 'computing…');
    let b;
    try { b = await ctx.ask('bicoherence', { sensor: +sel.value, era: era.value }); } catch (e) { return failure(bicBox, e); }
    const top = Math.max(b.level * 3, 0.05);
    bicBox.replaceChildren(raster(b.b2.slice().reverse(), { map: 'plum', lo: 0, hi: top, height: 280 }), ramp('plum', 0, top, 'b²'),
      ctx.finding({ title: `${fmt(100 * b.shareAbove, 3)}% of frequency pairs coupled (5% expected by chance)`, kind: 'measured',
        statement: `After the calendars' lines are removed, pairs of frequencies whose phases move together show a street that answers irregular disturbance (weather, events) nonlinearly. ${b.above.length ? `Strongest: ${b.above.slice(0, 4).map(([x, y, v]) => `(${fmt(x, 3)}, ${fmt(y, 3)})/day b² ${fmt(v, 2)}`).join('; ')}.` : ''}`,
        method: `Bicoherence over ${b.segments} two-week segments of the residual of the harmonic fit.`, limits: 'Quadratic coupling only.' }));
  };
  for (const el of [sel, era, order, range]) el.addEventListener('change', run);
  run();
}
