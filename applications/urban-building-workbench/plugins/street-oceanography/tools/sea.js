// Tides and currents of the city: every Melbourne sensor's tide and current read together, as an
// oceanographer reads a coast's tide gauges and current meters. When the city's high water comes and which
// streets have two crests a day; the commuters' flood toward the centre of the counts and its ebb, and the
// people the measured cross-sections take in; the arrow of time on the streets (the walkers' entropy
// production, bounded below by the measured edges); which way walkers keep; and the surge on the clock's date.
// The tide layer and the currents layer draw the same on the map.

import { sensorBlocks } from '../layers/blocks.js';

const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.round((x % 1) * 60) % 60).padStart(2, '0')}`;
const DAYS = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun'];

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { lines, bars, legend } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Tides and currents of the city'),
    h('span', { class: 'lede' }, 'Melbourne\'s pedestrian sensors read as tide gauges and current meters: the counts since October 2024, the totals at every sensor, and each way at those that count each way.')));
  const box = h('div', {});
  root.append(box);
  busy(box, 'reading every street\'s tide and current…');
  let tides, cur;
  try { [tides, cur] = await Promise.all([ctx.ask('tides'), ctx.ask('currents')]); } catch (e) { return failure(box, e); }
  box.replaceChildren();

  // ---- high water ----
  const wd = tides.map((t) => t.weekday.highAt), we = tides.map((t) => t.weekend.highAt);
  const bins = (xs) => Array.from({ length: 48 }, (_, k) => xs.filter((x) => x >= k / 2 && x < (k + 1) / 2).length);
  const bw = bins(wd), be = bins(we);
  const evening = wd.filter((x) => x >= 16.5 && x < 18.5).length, noon = wd.filter((x) => x >= 11.5 && x < 14.5).length;
  box.append(h('div', { class: 'section' }, h('h2', {}, 'When the water is high'),
    bars(bw.map((y, k) => ({ x: k / 2 + 0.25, y, color: css('--minor'), label: `${hhmm(k / 2)}–${hhmm(k / 2 + 0.5)}: ${y} streets on weekdays, ${be[k]} at weekends` })),
      { xlabel: 'hour of high water on the mean weekday', ylabel: 'streets', height: 180 }),
    ctx.finding({ title: `High water at the evening commute on ${evening} streets, at midday on ${noon}`, kind: 'measured',
      statement: `Of ${tides.length} streets, ${evening} reach their weekday high water between 16:30 and 18:30, the hour the offices empty toward the stations, and ${noon} between 11:30 and 14:30, when the centre's own walkers are out. At weekends high water comes later on ${tides.filter((t) => t.weekend.highAt > t.weekday.highAt + 0.5).length} streets and earlier on ${tides.filter((t) => t.weekend.highAt < t.weekday.highAt - 0.5).length}.`,
      method: 'Each street\'s tide: least squares at the calendars\' constituents (the week with its days, and the year) on its hourly counts in civil time; its mean weekday, minute by minute, and its highest minute.',
      limits: 'A mean over two years; a street with two crests has its higher one counted. The tide layer\'s "time of high water" draws it on the map.' })));

  // ---- the kind of tide ----
  const FORMS = ['semidiurnal', 'mixed, mainly semidiurnal', 'mixed, mainly diurnal', 'diurnal'];
  const counts = FORMS.map((f) => tides.filter((t) => t.formClass === f).length);
  const two = tides.filter((t) => t.form < 3).sort((a, b) => a.form - b.form).slice(0, 6);
  box.append(h('div', { class: 'section' }, h('h2', {}, 'Two crests a day, or one'),
    bars(counts.map((y, k) => ({ x: k, y, color: css('--accent'), label: `${FORMS[k]}: ${y}` })), { labels: ['semidiurnal', 'mixed, semi', 'mixed, diurnal', 'diurnal'], ylabel: 'streets', height: 170 }),
    ctx.finding({ title: `${counts[3]} diurnal streets, ${counts[2]} mixed, ${counts[1] + counts[0]} mainly semidiurnal`, kind: 'measured',
      statement: `Tides are classed by the ratio F of their diurnal species to their semidiurnal. Above 3, one crest a day dominates: most of the centre's streets, whose day is one broad swell from morning to evening. Below 3 the second harmonic, a morning and an evening crest, takes a larger share: ${two.map((t) => `${t.name} (F = ${fmt(t.form, 3)})`).join(', ')}.`,
      method: 'F = ΣA of the constituents of one day (D1 and its weekly sidebands) over ΣA of two (D2 and its sidebands), Courtier\'s (1938) classes.',
      limits: 'The sea\'s thresholds, borrowed: a convention, not a law of streets.' })));

  // ---- the commuters' currents ----
  const meters = sensorBlocks(ctx, cur.sensors).filter((f) => f.toward.length === 2 && f.toward.every(Boolean));
  let wsum = 0;
  const centre = [0, 0];
  for (const f of meters) centre[0] += f.transport * f.foot[0], centre[1] += f.transport * f.foot[1], wsum += f.transport;
  centre[0] /= wsum, centre[1] /= wsum;
  const sign = (f) => { const d = f.toward[0]; return d[0] * (centre[0] - f.foot[0]) + d[1] * (centre[1] - f.foot[1]) > 0 ? 1 : -1; };
  const inwardDay = Array.from({ length: 24 }, (_, hr) => {
    let s = 0;
    for (const f of meters) for (let d = 0; d < 5; d++) s += (sign(f) * f.net[24 * d + hr]) / 5;
    return s;
  });
  const taken = [];
  inwardDay.reduce((acc, v, k) => { const a = k >= 5 ? acc + v : 0; taken.push(a); return a; }, 0);
  const prism = Math.max(...taken), prismAt = taken.indexOf(prism) + 1;
  const xs = Array.from({ length: 24 }, (_, k) => k + 0.5);
  box.append(h('div', { class: 'section' }, h('h2', {}, 'The commuters\' flood and ebb'),
    lines([{ name: 'net toward the centre', x: xs, y: inwardDay, color: css('--accent'), width: 2 }, { name: 'taken in since 05:00', x: xs.map((x) => x + 0.5), y: taken, color: css('--minor'), dash: '5 3' }],
      { xlabel: 'hour of the mean weekday', ylabel: 'people an hour (dashed: people)', height: 220, xlim: [0, 24], marks: [{ x: 12, label: 'noon' }] }),
    legend([{ name: 'net current toward the centre of the counts, summed over the sensors', color: css('--accent') }, { name: 'what the measured cross-sections have taken in since 05:00', color: css('--minor') }]),
    ctx.finding({ title: `A flood of ${Math.round(prism).toLocaleString('en')} people in by ${hhmm(prismAt)}, and out again by night`, kind: 'measured',
      statement: `On the mean weekday the ${meters.length} counting lines carry a net current toward the centre of the city's counts from early morning, strongest at ${hhmm(xs[inwardDay.indexOf(Math.max(...inwardDay))])}, and away from it after ${hhmm(xs[inwardDay.findIndex((v, k) => k > 12 && v < 0)] || 17)}. By ${hhmm(prismAt)} they have taken in ${Math.round(prism).toLocaleString('en')} more people than they let out: the prism of the city's tide through these cross-sections, as the water an estuary takes in between low and high tide. By midnight the same lines have let out ${Math.round(Math.abs(taken[23] - prism)).toLocaleString('en')}, ${taken[23] < 0 ? `${Math.round(-taken[23]).toLocaleString('en')} more than they took in: the day's arrivals come partly by ways no sensor counts (from the trains beneath, say) and leave past the sensors` : 'not all of it'}.`,
      method: 'Each sensor\'s mean count each way at every hour of the week (its current meter); the sign of its net current toward the count-weighted centre of the sensors; summed, and accumulated from 05:00.',
      limits: 'The sensors are a sample of the cross-sections, not a cordon: people also arrive by train beneath the streets, so this is what the measured lines carry, not the city\'s whole prism. A walker may cross two lines.' })));

  // ---- the arrow of time ----
  const sigma = Array.from({ length: 168 }, (_, k) => meters.reduce((s, f) => s + (Number.isFinite(f.irreversibility[k]) ? f.irreversibility[k] : 0), 0));
  const peak = Math.max(...sigma), peakAt = sigma.indexOf(peak), least = Math.min(...sigma), leastAt = sigma.indexOf(least);
  const marks = [24, 48, 72, 96, 120, 144].map((x, i) => ({ x, label: DAYS[i + 1] }));
  box.append(h('div', { class: 'section' }, h('h2', {}, 'The arrow of time on the streets'),
    lines([{ name: 'entropy production, measured edges', x: sigma.map((_, k) => k + 0.5), y: sigma, color: css('--series-2'), width: 1.8 }],
      { xlabel: 'hour of the week', ylabel: 'k_B an hour', height: 200, xlim: [0, 168], marks }),
    ctx.finding({ title: `Most irreversible at ${DAYS[Math.floor(peakAt / 24)]} ${hhmm((peakAt % 24) + 0.5)}, nearly reversible at ${hhmm((leastAt % 24) + 0.5)}`, kind: 'proved',
      statement: `Take the walkers as independent, each a Markov jump process among places. Their entropy production is a sum over the edges of the network of streets, each term (n₊ − n₋) ln(n₊/n₋) ≥ 0 for the counts n₊ and n₋ each way (Schnakenberg). The measured edges' terms therefore bound the whole from below: ${Math.round(peak).toLocaleString('en')} k_B an hour at the peak, when the morning commute runs one way, against ${Math.round(least).toLocaleString('en')} at the quietest hour. A street where as many walk each way, whatever their number, adds nothing: the arrow of time on a street is its one-wayness, not its crowd.`,
      method: 'Each counting line\'s mean count each way at every hour of the week; the edge term of Schnakenberg\'s entropy production (the domain\'s flow/currents.hpp, held equal to physis\'s on a driven ring), summed over the lines.',
      limits: 'A lower bound, from the hour\'s means; the independence of walkers is the model\'s. Hours with fewer than one walker either way are left out.' })));

  // ---- handedness ----
  const hd = cur.handedness;
  box.append(h('div', { class: 'section' }, h('h2', {}, 'Which way walkers keep'),
    ctx.finding({ title: `${hd.left} of ${hd.footpaths} footpaths keep left`, kind: 'measured',
      statement: `A sensor on one footpath counts along it, and over two years more walkers pass it one way than the other: its residual current. If walkers chose their side of the street as left-hand traffic chooses its lane, the residual would run east on north footpaths, west on south ones, north on west ones and south on east ones. ${hd.left} of the ${hd.footpaths} footpaths whose sensor names its side do so. With no handedness the count would be binomial(${hd.footpaths}, ½): P(≥ ${hd.left}) = ${fmt(hd.pLeft, 2)} (two-sided ${fmt(hd.pTwoSided, 2)}). The mirror symmetry that exchanges a street's two sides is broken, and broken the way the road's traffic breaks it.`,
      method: 'Each sensor\'s side of the street from its name ("Town Hall (West)"), the directions it counts, and the sign of its mean net current; a sign test on tyche\'s binomial law.',
      limits: 'Footpaths of one street share their walkers and are not independent; a station entrance or a crossing on one side can make a residual of its own.' })));

  // ---- the surge on the date ----
  const surgeBox = h('div', { class: 'section' });
  box.append(surgeBox);
  const showSurge = async () => {
    const date = ctx.clock.date;
    busy(surgeBox, `the counts of ${date}…`);
    let d;
    try { d = await ctx.ask('sea', { date }); } catch (e) { return failure(surgeBox, e); }
    const pred = Array(24).fill(0), obs = Array(24).fill(0), per = [];
    for (const s of d.sensors) {
      if (!s.observed.some(Number.isFinite)) continue;
      let po = 0, oo = 0;
      for (let hr = 0; hr < 24; hr++) {
        const p = (s.level[4 * hr] + 2 * s.level[4 * hr + 2] + s.level[4 * hr + 4]) / 4;
        if (!Number.isFinite(s.observed[hr])) continue;
        pred[hr] += p, obs[hr] += s.observed[hr], po += p, oo += s.observed[hr];
      }
      const t = tides.find((x) => x.id === s.id);
      if (po > 500) per.push({ name: t ? t.name : String(s.id), ratio: oo / po });
    }
    if (!d.observed) {
      surgeBox.replaceChildren(h('h2', {}, `The surge on ${date}`), h('p', { class: 'note' }, `${date} lies outside the record (October 2024 to October 2026): set the clock's date within it to see a day's surge.`));
      return;
    }
    per.sort((a, b) => a.ratio - b.ratio);
    const total = obs.reduce((a, v) => a + v, 0) / pred.reduce((a, v) => a + v, 0) - 1;
    surgeBox.replaceChildren(h('h2', {}, `The surge on ${date}`),
      lines([{ name: 'the tide (the calendars)', x: xs, y: pred, color: css('--minor'), width: 2 }, { name: 'observed', x: xs, y: obs, color: css('--accent'), width: 2 }],
        { xlabel: 'hour of the day', ylabel: 'people an hour, all sensors', height: 200, xlim: [0, 24] }),
      legend([{ name: 'the tide', color: css('--minor') }, { name: 'observed', color: css('--accent') }]),
      ctx.finding({ title: `${total >= 0 ? '+' : '−'}${fmt(Math.abs(100 * total), 2)}% against the tide`, kind: 'measured',
        statement: `${d.observed} streets recorded${d.holiday ? `, on ${d.holiday}` : ''}${d.weather ? `, ${d.weather.rain} mm of rain and ${d.weather.tmin}–${d.weather.tmax} °C` : ''}. Furthest above their tide: ${per.slice(-3).reverse().map((x) => `${x.name} (${fmt(100 * (x.ratio - 1), 2)}%)`).join(', ')}; furthest below: ${per.slice(0, 3).map((x) => `${x.name} (${fmt(100 * (x.ratio - 1), 2)}%)`).join(', ')}.`,
        method: 'Each street\'s count that day against its tide, the prediction from the calendars alone; the surge is what they leave, as a storm surge is the sea\'s departure from its astronomical tide.',
        limits: 'The tide carries no holidays, weather or events: those are the surge. Change the clock\'s date to read another day.' }));
  };
  showSurge();
  let lastDate = ctx.clock.date;
  const onClock = () => { if (ctx.clock.date !== lastDate) { lastDate = ctx.clock.date; showSurge(); } };
  ctx.clock.on(onClock);
  ctx.onDispose(() => { const i = ctx.clock.listeners.indexOf(onClock); if (i >= 0) ctx.clock.listeners.splice(i, 1); });
}
