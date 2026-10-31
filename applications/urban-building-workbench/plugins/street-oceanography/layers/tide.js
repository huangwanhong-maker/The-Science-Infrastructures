// The tide on the streets. Each Melbourne sensor's block of street is drawn as a channel of water standing
// to the level of its pedestrians, one metre for every hundred people an hour, as a tide gauge reads a port.
// The level comes from the engine's tide: the calendars' constituents fitted to the sensor's counts since
// October 2024 (the domain's rhythm/tide.hpp), evaluated through the clock's civil day. It rises and falls
// as the clock runs. Four readings:
//   - the tide: the level the calendars predict at the clock's hour;
//   - observed and surge: on a date the record holds, the count that hour, coloured by its departure from
//     the tide, with a tide mark floating at the predicted level;
//   - time of high water: the hour at which the mean weekday's tide is highest, as a port's establishment;
//   - kind of tide: two crests a day or one, by the ratio of the diurnal to the semidiurnal species
//     (Courtier's classes), the water standing at the daily range.

import { sensorBlocks } from './blocks.js';

const PER = 0.01;   // metres of water per person an hour
const MODES = [
  { id: 'tide', title: 'the tide', about: 'The level the calendars alone predict at this hour: the week with its days and the year.' },
  { id: 'surge', title: 'observed and surge', about: 'The count that hour on this date, coloured by how far it stands above or below the tide; a dark tide mark floats at the prediction.' },
  { id: 'high', title: 'time of high water', about: 'The hour of high water on the mean weekday, the water standing at its height.' },
  { id: 'form', title: 'kind of tide', about: 'Two crests a day (commuters\' streets) or one: the ratio of the day\'s first harmonic species to its second, in Courtier\'s classes; the water stands at the daily range.' },
];
const FORMS = ['semidiurnal', 'mixed, mainly semidiurnal', 'mixed, mainly diurnal', 'diurnal'];
const hhmm = (x) => `${String(Math.floor(x)).padStart(2, '0')}:${String(Math.round((x % 1) * 60) % 60).padStart(2, '0')}`;
const lerp = (arr, x) => {
  if (!arr) return NaN;
  const i = Math.max(0, Math.min(arr.length - 1, Math.floor(x))), u = x - i, a = arr[i], b = arr[Math.min(arr.length - 1, i + 1)];
  if (Number.isFinite(a) && Number.isFinite(b)) return a + Math.max(0, Math.min(1, u)) * (b - a);
  return Number.isFinite(a) ? a : NaN;
};

export async function layer(ctx) {
  const { gis, ask, dom, G } = ctx;
  const { h, fmt } = dom;
  const id = ctx.layerId(), markId = ctx.layerId(':marks');
  let mode = 'tide', clock = ctx.clock, sea = null, seaDate = null, features = [], drawn = false;
  const byId = new Map();
  const status = h('div', { class: 'faint' }, 'fitting every street\'s tide (the first time, about ten seconds)…');
  const ramp = h('div', {});
  const about = h('div', { class: 'faint' });
  const keys = h('div', { class: 'keys' }, MODES.map((m) => h('button', { class: `key${m.id === mode ? ' on' : ''}`, 'data-mode': m.id, onclick: () => setMode(m.id) }, m.title)));
  const legend = h('div', { class: 'card' }, h('b', {}, 'The tide on the streets'), keys, ramp, about, status);

  const predicted = (f) => lerp(sea && sea.get(f.id)?.level, clock.hour * 4);
  const observed = (f) => lerp(sea && sea.get(f.id)?.observed, Math.max(0, clock.hour - 0.5));
  const surge = (f) => { const o = observed(f), p = predicted(f); return Number.isFinite(o) && Number.isFinite(p) ? Math.log2((o + 30) / (p + 30)) : NaN; };
  let seaColour, maxHigh = 1;
  const formColour = G.sequential(0, 3, 'orchid', { from: 0.2, to: 0.95 });
  const hourColour = (x) => G.sequential(0, 1, 'warm', { from: 0, to: 1 })(x);
  const surgeColour = G.diverging(1, 'sea-brick');
  const grey = [0.8, 0.78, 0.76];
  const height = (f) => {
    if (mode === 'high') return f.weekday.high * PER;
    if (mode === 'form') return Math.max(0.3, f.weekday.range * PER);
    if (mode === 'surge' && Number.isFinite(observed(f))) return observed(f) * PER;
    return (predicted(f) || 0) * PER;
  };
  const colour = (f) => {
    if (mode === 'high') return hourColour(f.weekday.highAt / 24) || grey;
    if (mode === 'form') return formColour(FORMS.indexOf(f.formClass)) || grey;
    if (mode === 'surge') return Number.isFinite(surge(f)) ? surgeColour(surge(f)) : grey;
    return seaColour(predicted(f)) || grey;
  };
  const pick = (f) => {
    const p = predicted(f), o = observed(f);
    if (mode === 'high') return `${f.name}\nhigh water on weekdays at ${hhmm(f.weekday.highAt)}, ${fmt(f.weekday.high, 3)} an hour`;
    if (mode === 'form') return `${f.name}\n${f.formClass} (F = ${fmt(f.form, 3)}); daily range ${fmt(f.weekday.range, 3)} an hour`;
    return `${f.name}\ntide ${Number.isFinite(p) ? Math.round(p).toLocaleString('en') : '—'} an hour${Number.isFinite(o) ? `; observed ${Math.round(o).toLocaleString('en')} (${o >= p ? '+' : '−'}${Math.round(Math.abs(100 * (o - p) / Math.max(p, 1)))}%)` : ''}`;
  };

  function showLegend() {
    for (const b of keys.children) b.classList.toggle('on', b.dataset.mode === mode);
    const m = MODES.find((x) => x.id === mode);
    about.textContent = m.about;
    const { ramp: r, legend: lg } = ctx.chart;
    if (mode === 'tide') ramp.replaceChildren(r('sea', 0, Math.round(maxHigh), 'people an hour'));
    else if (mode === 'surge') ramp.replaceChildren(h('div', { class: 'legend' }, h('span', {}, 'half the tide'), h('span', { class: 'ramp', style: `background: linear-gradient(to right, ${[-1, -0.5, 0, 0.5, 1].map((t) => `rgb(${surgeColour(t).map((c) => Math.round(255 * c)).join(',')})`).join(',')})` }), h('span', {}, 'twice')));
    else if (mode === 'high') ramp.replaceChildren(r('warm', 0, 24, 'hour of the day'));
    else ramp.replaceChildren(lg(FORMS.map((name, i) => ({ name, color: `rgb(${formColour(i).map((c) => Math.round(255 * c)).join(',')})`, dot: true }))));
    const s = seaDate && sea && sea.meta;
    status.textContent = s ? (s.observed
      ? `${s.date}: ${s.observed} streets recorded${s.holiday ? `; ${s.holiday}` : ''}${s.weather ? `; ${s.weather.rain} mm of rain, ${s.weather.tmin}–${s.weather.tmax} °C` : ''}. Water one metre deep for every hundred people an hour.`
      : `${s.date} lies outside the record (October 2024 to October 2026): the tide alone. Water one metre deep for every hundred people an hour.`) : status.textContent;
  }

  function draw() {
    if (!features.length) return;
    if (!drawn) {
      gis.addSlabs(id, features, { line: (f) => f.line, width: (f) => f.width, height, color: colour, opacity: 0.5, z: 0.2, pick });
      gis.addSlabs(markId, features, { line: (f) => f.line, width: (f) => f.width + 0.6, base: (f) => 0.2 + (predicted(f) || 0) * PER, height: () => 0.5, color: () => [0.24, 0.19, 0.22], opacity: 0.9, pick: false });
      drawn = true;
    } else {
      gis.restyle(id, {});
      gis.restyle(markId, {});
    }
    gis.setVisible(markId, mode === 'surge');
  }

  function setMode(m) { mode = m; showLegend(); draw(); }

  async function loadDay(date) {
    if (seaDate === date) return;
    seaDate = date;
    const d = await ask('sea', { date });
    if (seaDate !== date) return;
    sea = new Map(d.sensors.map((x) => [x.id, x]));
    sea.meta = d;
    showLegend();
    draw();
  }

  // The tides first, then the day; the legend is shown at once.
  (async () => {
    try {
      const tides = await ask('tides');
      features = sensorBlocks(ctx, tides);
      for (const f of features) byId.set(f.id, f);
      maxHigh = Math.max(1, ...features.map((f) => f.weekday.high));
      seaColour = G.sequential(0, maxHigh, 'sea', { from: 0.12, to: 1 });
      status.textContent = '';
      await loadDay(clock.date);
    } catch (e) {
      status.textContent = String(e.message || e);
    }
  })();
  showLegend();

  return {
    legend,
    setTime(c) {
      clock = c;
      if (!features.length) return;
      if (c.date !== seaDate) loadDay(c.date);
      else draw();
    },
    describe(hit) {
      const f = hit.feature;
      const p = predicted(f), o = observed(f);
      const day = sea && sea.get(f.id);
      const xs = day ? day.level.map((_, k) => k / 4) : [];
      const chart = day ? ctx.chart.lines([{ name: 'the tide', x: xs, y: day.level, color: ctx.dom.css('--accent'), width: 2 }],
        { height: 150, width: 360, xlim: [0, 24], xlabel: 'hour of the day', ylabel: 'people an hour', marks: [{ x: clock.hour, label: hhmm(clock.hour) }],
          points: day.observed.map((v, k) => (Number.isFinite(v) ? { x: k + 0.5, y: v, r: 2.6, color: ctx.dom.css('--minor-deep'), label: `${k}:00–${k + 1}:00 observed ${Math.round(v)}` } : null)).filter(Boolean) }) : null;
      const sel = {
        kind: 'sensor', id: f.id, name: f.name, title: f.name,
        lines: [
          `tide ${Number.isFinite(p) ? Math.round(p).toLocaleString('en') : '—'} an hour now${Number.isFinite(o) ? `; observed ${Math.round(o).toLocaleString('en')}, ${o >= p ? 'a surge of' : 'below the tide by'} ${Math.round(Math.abs(100 * (o - p) / Math.max(p, 1)))}%` : ''}`,
          `high water ${hhmm(f.weekday.highAt)} on weekdays (${Math.round(f.weekday.high).toLocaleString('en')}), ${hhmm(f.weekend.highAt)} at weekends (${Math.round(f.weekend.high).toLocaleString('en')}); low water ${hhmm(f.weekday.lowAt)}`,
          `${f.formClass} (F = ${fmt(f.form, 3)}); the year swings it by ±${Math.round(f.annual)} an hour; ${Math.round(100 * f.varianceExplained)}% of its hours' variance is the calendars'`,
        ],
        extra: chart ? dom.h('div', {}, chart, dom.h('div', { class: 'faint' }, 'line: the tide through the day; dots: the counts observed that day')) : null,
      };
      sel.actions = [
        { title: 'its rhythm', run: () => ctx.openTool('rhythm', sel) },
        { title: 'weather and holidays', run: () => ctx.openTool('response', sel) },
        { title: 'tides and currents of the city', run: () => ctx.openTool('sea', sel) },
      ];
      return sel;
    },
    dispose() { gis.remove(markId); },
  };
}
