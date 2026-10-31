// Weather and holidays: how much rain, warmth and cold, and public holidays change a street's counts, with the
// hour of the week and the season taken out first; holidays compared with the same weekday a week away; and
// every outdoor sensor side by side, since every street here receives the same weather and the same holidays.

export async function render(root, ctx) {
  const { h, busy, failure, fmt, css } = ctx.dom;
  const { bars } = ctx.chart;
  root.append(h('div', { class: 'view-head' }, h('h1', {}, 'Weather and holidays'),
    h('span', { class: 'lede' }, 'Every street receives the same rain and the same holidays; what it does with them is its own. Pick a sensor on the map, or choose one here.')));
  let sensors;
  try { sensors = (await ctx.ask('sensors')).filter((s) => s.locationType === 'Outdoor').sort((a, b) => a.name.localeCompare(b.name)); }
  catch (e) { return failure(root, e); }
  const sel = h('select', {}, sensors.map((s) => h('option', { value: s.id }, `${s.name} (${s.id})`)));
  const era = h('select', {}, h('option', { value: 'recent' }, 'Oct 2024 – Oct 2026'), h('option', { value: 'history' }, '2009 – 2022'));
  root.append(h('div', { class: 'keys' }, h('label', { class: 'pick' }, 'sensor ', sel), h('label', { class: 'pick' }, 'years ', era)));
  const out = h('div', {});
  const all = h('div', { class: 'section' });
  root.append(out, all);
  sel.value = ctx.selection && ctx.selection.kind === 'sensor' ? String(ctx.selection.id) : String((sensors.find((s) => /Southbank/.test(s.name)) || sensors[0]).id);
  ctx.onSelect((s) => { if (s && s.kind === 'sensor') { sel.value = String(s.id); run(); } });

  const run = async () => {
    busy(out, 'fitting…');
    let d;
    try { d = await ctx.ask('response', { sensor: +sel.value, era: era.value }); } catch (e) { return failure(out, e); }
    const rain = d.forcings[0];
    const lagChart = bars(rain.lags.map((l, i) => ({ x: l, y: 100 * Math.expm1(rain.beta[i]), color: css('--series-2'), label: `${l} h after: ${fmt(100 * Math.expm1(rain.beta[i]), 3)}% per mm/h (± ${fmt(100 * rain.se[i], 2)})` })),
      { xlabel: 'hours since the rain', ylabel: '% change per mm/h', height: 190, width: 360 });
    const table = h('table', { class: 'data' }, h('tr', {}, ['', 'per', 'change', '±', 'z'].map((x) => h('th', {}, x))),
      d.forcings.map((f) => h('tr', {}, h('td', {}, f.name), h('td', { class: 'left' }, f.unit), h('td', {}, `${fmt(f.percent, 3)}%`), h('td', {}, `${fmt(100 * f.totalSe, 2)}%`), h('td', {}, fmt(f.z, 3)))));
    const ho = d.holidays;
    out.replaceChildren(h('h3', {}, d.name), table, lagChart,
      ctx.finding({ title: `Rain: ${fmt(rain.percent, 3)}% per mm/h`, kind: 'measured',
        statement: `A covered arcade and an open promenade receive the same rain; this is how much of it ${d.name} turns into fewer (or more) people, over the hours after it falls.`,
        method: `Regression of ln(1 + count) on rain (this hour and the three before), the temperature's departure from its usual value for the month and hour, and public holidays, with the hour of the week, the season and a trend taken out; errors clustered by day. ${fmt(d.observations)} hours.`,
        limits: 'Hourly reanalysis weather at the city centre (Open-Meteo), not a gauge on the street; what is fitted describes the street as it was.' }),
      ctx.finding({ title: `A weekday holiday: ${fmt(ho.percent, 3)}% (p = ${fmt(ho.p, 2)})`, kind: 'measured',
        statement: `${ho.pairs} weekday holidays, each paired with the same weekday a week before or after.`,
        method: 'Mean of the paired differences of ln(1 + daily total); p by random sign flips within pairs.', limits: 'Weather on the paired days is not removed.' }));
  };

  const everyone = async () => {
    busy(all, 'fitting every outdoor sensor (once, about a minute)…');
    let rows;
    try { rows = (await ctx.ask('summary')).filter((r) => Number.isFinite(r.rain)); } catch (e) { return failure(all, e); }
    rows.sort((a, b) => a.rain - b.rain);
    const rainBars = bars(rows.map((r, i) => ({ x: i, y: 100 * Math.expm1(r.rain), color: Math.abs(r.rain) > 2 * r.rainSe ? css('--series-2') : css('--series-3'),
      label: `${r.name}: ${fmt(100 * Math.expm1(r.rain), 3)}% per mm/h (± ${fmt(200 * r.rainSe, 2)}%)` })), { ylabel: '% per mm/h of rain', height: 240, labels: rows.map((r) => r.name.slice(0, 20)) });
    const hol = rows.filter((r) => Number.isFinite(r.holiday)).sort((a, b) => a.holiday - b.holiday);
    const holBars = bars(hol.map((r, i) => ({ x: i, y: 100 * Math.expm1(r.holiday), color: r.holidayP < 0.05 ? css('--accent') : css('--series-3'), label: `${r.name}: ${fmt(100 * Math.expm1(r.holiday), 3)}% (p = ${fmt(r.holidayP, 2)})` })),
      { ylabel: '% on a weekday holiday', height: 240, labels: hol.map((r) => r.name.slice(0, 20)) });
    all.replaceChildren(h('h2', {}, 'Every street under the same weather and holidays'), h('h3', {}, 'Rain'), rainBars, h('h3', {}, 'Weekday holidays'), holBars,
      h('p', { class: 'caption' }, 'Coloured where the effect is clear of its uncertainty. Some streets fill on holidays and others empty: one holiday, read by two streets.'));
  };
  sel.addEventListener('change', run);
  era.addEventListener('change', run);
  run();
  all.append(h('button', { class: 'plain', onclick: everyone }, 'Compare every outdoor sensor'));
}
