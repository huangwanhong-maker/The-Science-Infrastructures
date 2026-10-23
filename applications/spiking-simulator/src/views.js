// The visualizer and profiler: what a run's results show.
//
// - activity: the recorded spikes as a raster, population by population, over their rates;
//   drag across either to look closer in time, double-click to see the whole run;
// - statistics: each population's rate, irregularity (CV, LV), count variability (Fano),
//   synchrony (chi, correlation) and rhythm, as the engine's analysis gives them, with the
//   spectra of the population rates, the distributions of intervals and rates, and, where each
//   neuron was given its own current, the rate against the current (an f-I curve);
// - states: the recorded components, one chart each;
// - weights: the mean and spread of plastic weights over time;
// - profile: where the engine's time went, phase by phase and thread by thread;
// - network: the populations, projections and devices drawn.

import { el, clear } from './dom.js';
import { Figure, line, band, bars, dots, extent, decades } from './charts.js';
import { colourOf, POPULATION_COLOURS } from './palette.js';
import { compact, duration } from './results.js';
import { short, format, parseField } from './units.js';

const ms = (t) => short(t * 1e3);
const WAIT = 'rgba(240, 226, 238, 0.16)';

// A unit with a prefix that suits the values: [factor, label] for SI values in `unit`.
function display(unit, values) {
  if (!['V', 'A', 'S', 'F', 's', 'Hz'].includes(unit)) return [1, unit === '1' ? '' : unit];
  let m = 0;
  for (const v of values) if (Number.isFinite(v)) m = Math.max(m, Math.abs(v));
  const prefixes = [[1e-12, 'p'], [1e-9, 'n'], [1e-6, 'µ'], [1e-3, 'm'], [1, ''], [1e3, 'k']];
  const [f, p] = [...prefixes].reverse().find(([s]) => m >= s) ?? prefixes[0];
  return [1 / f, p + unit];
}

const colourFor = (results, name) => colourOf(Math.max(0, results.populations.findIndex((p) => p.name === name)));

// ---- activity ----------------------------------------------------------------------------------------
export function activity(root, results, spec, figures, { domain } = {}) {
  clear(root);
  const T = results.time;
  // Time in ms for short runs, in s for long ones.
  const long = T >= 5;
  const tl = long ? 'time (s)' : 'time (ms)', tf = long ? (t) => short(t) : ms, tt = (t) => (long ? `${short(Number(t.toPrecision(5)))} s` : `${ms(t)} ms`);
  const view = { domain: domain && domain[1] > domain[0] ? [Math.max(0, domain[0]), Math.min(T, domain[1])] : [0, T] };
  const recs = results.spikes;
  if (!recs.length && !results.rates.length) { root.append(el('p', { class: 'note', text: 'Nothing was recorded as spikes or rates: add a spike or rate recorder.' })); return; }
  root.append(el('p', { class: 'note', text: 'Drag across a chart to look closer in time; double-click to see the whole run.' }));
  const all = [];
  const select = (d) => { view.domain = d; for (const f of all) f.render(); };
  const reset = () => select([0, T]);

  if (recs.length) {
    let offset = 0;
    const bands = recs.map((r) => { const b = { r, from: offset, to: offset + r.count, colour: colourFor(results, r.population) }; offset += r.count; return b; });
    const rows = offset;
    const fig = new Figure(root, {
      title: 'Spikes', subtitle: `${compact(recs.reduce((n, r) => n + r.times.length, 0))} spikes of ${compact(rows)} recorded neurons`,
      height: Math.min(520, Math.max(240, rows * 0.7 + 60)),
      legend: bands.map((b) => ({ label: `${b.r.population} (${compact(b.r.count)})`, colour: b.colour, kind: 'dot' })),
      x: { label: tl, format: tf }, y: { label: 'neuron', ticks: () => bands.map((b) => (b.from + b.to) / 2), format: (v) => bands.find((b) => (b.from + b.to) / 2 === v)?.r.population ?? '', grid: false },
      margin: { left: 96 },
      onSelect: select, onReset: reset,
      hover: (t, row) => {
        const b = bands.find((q) => row >= q.from && row < q.to);
        if (!b) return null;
        const i = Math.floor(row - b.from) + b.r.first;
        const n = b.r.senders.reduce((c, s) => c + (s === i ? 1 : 0), 0);
        return { title: tt(t), rows: [{ colour: b.colour, value: `neuron ${i}`, label: `of ${b.r.population}` }, { value: `${n} spikes`, label: `${short(n / T)} Hz over the run` }] };
      },
    });
    fig.set(view.domain, [rows, 0], (ctx, x, y, a) => {
      const [t0, t1] = view.domain;
      fig.dx = view.domain;
      const rh = Math.max(1, (a.bottom - a.top) / rows);
      const w = Math.max(1, ((a.right - a.left) * results.resolution) / (t1 - t0));
      for (const b of bands) {
        ctx.fillStyle = b.colour;
        const { times, senders } = b.r;
        for (let k = 0; k < times.length; k++) {
          const t = times[k];
          if (t < t0 || t > t1) continue;
          ctx.fillRect(x(t) - w / 2, y(b.from + senders[k] - b.r.first + 1), w, Math.max(1, rh * 0.85));
        }
        ctx.fillStyle = 'rgba(240, 226, 238, 0.12)';
        ctx.fillRect(a.left, Math.round(y(b.to)) , a.right - a.left, 1);
      }
    });
    const render = fig.render.bind(fig);
    fig.render = () => { fig.dx = view.domain; render(); };
    all.push(fig);
    figures.push(fig);
  }

  // The rates: the populations' own rate recorders, or else the recorded neurons' rate.
  const series = results.rates.length
    ? results.rates.map((r) => ({ name: r.population, times: r.times, rate: r.rate, bin: r.bin, whole: true }))
    : (results.analysis ?? []).map((a) => ({ name: a.population, times: a.population_rate.times, rate: a.population_rate.rate, bin: a.population_rate.bin, whole: false }));
  if (series.length) {
    const fig = new Figure(root, {
      title: 'Population rate', subtitle: series[0].whole ? `each population in bins of ${ms(series[0].bin)} ms` : `the recorded neurons, in bins of ${ms(series[0].bin)} ms`,
      height: 200, legend: series.length > 1 ? series.map((s) => ({ label: s.name, colour: colourFor(results, s.name) })) : [],
      x: { label: tl, format: tf }, y: { label: 'Hz' }, onSelect: select, onReset: reset,
      margin: { left: recs.length ? 96 : 58 },
      hover: (t) => ({ title: tt(t), rows: series.map((s) => { const k = Math.min(s.rate.length - 1, Math.max(0, Math.floor(t / s.bin))); return { colour: colourFor(results, s.name), value: `${short(s.rate[k])} Hz`, label: s.name }; }) }),
    });
    const visibleMax = () => {
      let m = 0;
      for (const s of series) for (let k = 0; k < s.times.length; k++) if (s.times[k] >= view.domain[0] && s.times[k] <= view.domain[1]) m = Math.max(m, s.rate[k]);
      return m || 1;
    };
    const render = fig.render.bind(fig);
    fig.render = () => { fig.dx = view.domain; fig.dy = [0, visibleMax() * 1.08]; render(); };
    fig.set(view.domain, [0, visibleMax() * 1.08], (ctx, x, y) => {
      for (const s of series) line(ctx, x, y, s.times, s.rate, colourFor(results, s.name), series.length > 2 ? 1.5 : 2);
    });
    all.push(fig);
    figures.push(fig);
  }
}

// ---- statistics -----------------------------------------------------------------------------------------
const NOTES = {
  rate: 'spikes per neuron per second over the window analysed: mean and spread over the neurons',
  cv: 'coefficient of variation of the intervals between spikes: 1 for a Poisson process, 0 for a clock',
  lv: 'local variation of successive intervals (Shinomoto et al. 2003): 1 for Poisson, 0 for a clock, blind to slow changes of rate',
  fano: 'variance over mean of the counts in windows: 1 for Poisson',
  chi: 'Golomb\'s synchrony: 1 when all fire together, about 1/√N when independent',
  corr: 'mean correlation of counts in windows over a sample of pairs',
  peak: 'the strongest rhythm in the population\'s rate (the peak of its spectrum)',
};

export function statistics(root, results, spec, figures, catalogue) {
  clear(root);
  const A = results.analysis ?? [];
  if (!A.length) { root.append(el('p', { class: 'note', text: 'Statistics are computed from recorded spikes: add a spike recorder.' })); return; }
  const window = A[0];
  root.append(el('p', { class: 'note', text: `Over ${ms(window.from)}–${ms(window.to)} ms; counts in windows of ${ms(window.options.window)} ms for the Fano factor and correlations, in bins of ${ms(window.options.sync_bin)} ms for synchrony, the rate in bins of ${ms(window.options.bin)} ms. A specification sets them under "analysis".` }));
  const cell = (v, digits = 3) => (Number.isFinite(v) ? short(Number(v.toPrecision(digits))) : '—');
  root.append(el('table', { class: 'data stats' },
    el('thead', {}, el('tr', {},
      el('th', { text: 'population' }), el('th', { text: 'neurons' }),
      ...[['rate (Hz)', NOTES.rate], ['CV', NOTES.cv], ['LV', NOTES.lv], ['Fano', NOTES.fano], ['χ', NOTES.chi], ['correlation', NOTES.corr], ['rhythm (Hz)', NOTES.peak]].map(([t, n]) => el('th', { title: n }, t, el('sup', { text: ' ?' }))))),
    el('tbody', {}, A.map((a) => el('tr', {},
      el('td', {}, el('i', { class: 'swatch', style: `background:${colourFor(results, a.population)}` }), a.population),
      el('td', { text: `${compact(a.active)} of ${compact(a.count)} fired` }),
      el('td', { text: `${cell(a.rate_mean)} ± ${cell(a.rate_sd, 2)}` }),
      el('td', { text: cell(a.cv) }), el('td', { text: cell(a.lv) }), el('td', { text: cell(a.fano) }),
      el('td', { text: cell(a.chi) }), el('td', { text: cell(a.correlation, 2) }),
      el('td', { text: Number.isFinite(a.peak_frequency) && a.spectrum.segments ? short(a.peak_frequency) : '—' }))))));
  root.append(el('p', { class: 'note', text: 'For a Poisson process CV, LV and the Fano factor are 1, and for a clock 0; χ is near 1/√N for independent neurons and 1 for neurons that fire together. Hover a heading for its definition.' }));

  const grid = el('div', { class: 'figs' });
  root.append(grid);
  // The spectra of the populations' rates.
  const withSpectrum = A.filter((a) => a.spectrum.frequency.length > 2);
  if (withSpectrum.length) {
    const lo = withSpectrum[0].spectrum.frequency[1];
    const hi = Math.max(...withSpectrum.map((a) => a.spectrum.frequency[a.spectrum.frequency.length - 1]));
    const ys = withSpectrum.flatMap((a) => a.spectrum.density.slice(1)).filter((v) => v > 0);
    const yd = extent(ys.map(Math.log10));
    const fig = new Figure(grid, {
      title: 'Spectrum of the population rate', subtitle: 'power spectral density by Welch\'s method',
      height: 230, legend: withSpectrum.length > 1 ? withSpectrum.map((a) => ({ label: a.population, colour: colourFor(results, a.population) })) : [],
      x: { label: 'Hz', log: true }, y: { label: 'Hz²/Hz', log: true, format: (v) => v.toExponential(0) },
      hover: (f) => ({ title: `${short(f)} Hz`, rows: withSpectrum.map((a) => { const k = nearest(a.spectrum.frequency, f); return { colour: colourFor(results, a.population), value: short(a.spectrum.density[k]), label: a.population }; }) }),
    });
    fig.set([lo, hi], [10 ** yd[0], 10 ** yd[1]], (ctx, x, y) => {
      for (const a of withSpectrum) line(ctx, x, y, a.spectrum.frequency.slice(1), a.spectrum.density.slice(1).map((v) => (v > 0 ? v : NaN)), colourFor(results, a.population));
    });
    figures.push(fig);
  }
  // The intervals between spikes, each population in its own chart.
  for (const a of A) {
    const counts = a.isi.counts;
    if (!counts.some((c) => c > 0)) continue;
    const edges = a.isi.edges;
    const first = counts.findIndex((c) => c > 0), last = counts.length - 1 - [...counts].reverse().findIndex((c) => c > 0);
    const fig = new Figure(grid, {
      title: `Intervals between spikes: ${a.population}`, subtitle: `${compact(counts.reduce((s, c) => s + c, 0))} intervals, logarithmic bins`, height: 200,
      x: { label: 'interval (ms)', log: true, format: (v) => short(v * 1e3) }, y: { label: 'count' },
      hover: (t) => { const k = edges.findIndex((e, i) => i < counts.length && t >= e && t < edges[i + 1]); return k < 0 ? null : { title: `${ms(edges[k])}–${ms(edges[k + 1])} ms`, rows: [{ colour: colourFor(results, a.population), value: String(counts[k]), label: 'intervals' }] }; },
    });
    fig.set([edges[Math.max(0, first - 1)], edges[Math.min(edges.length - 1, last + 2)]], [0, Math.max(...counts) * 1.08], (ctx, x, y) => bars(ctx, x, y, edges, counts, colourFor(results, a.population)));
    figures.push(fig);
  }
  // The rates against each neuron's own current: an f-I curve.
  for (const a of A) {
    const pop = spec.populations.find((p) => p.name === a.population);
    if (!pop || !Array.isArray(pop.bias)) continue;
    const bias = pop.bias.slice(a.first, a.first + a.count).map((q) => toSI(q));
    if (bias.some((v) => !Number.isFinite(v))) continue;
    const unit = catalogue?.get(pop.model)?.input ?? 'A';
    const [f, label] = display(unit, bias);
    const xs = bias.map((v) => v * f);
    const fig = new Figure(grid, {
      title: `Rate against current: ${a.population}`, subtitle: 'each neuron given its own constant current (the f–I curve)', height: 220,
      x: { label }, y: { label: 'Hz' },
      hover: (I) => { const k = nearest(xs, I); return { title: `${short(xs[k])} ${label}`, rows: [{ colour: colourFor(results, a.population), value: `${short(a.rates[k])} Hz`, label: `neuron ${a.first + k}` }] }; },
    });
    fig.set(extent(xs), [0, Math.max(1, ...a.rates) * 1.08], (ctx, x, y) => {
      line(ctx, x, y, xs, a.rates, colourFor(results, a.population), 1.5, 0.55);
      dots(ctx, x, y, xs, a.rates, colourFor(results, a.population), 2.5);
    });
    figures.push(fig);
  }
  // The distribution of rates over the neurons.
  for (const a of A) {
    const { edges, counts } = a.rate_histogram;
    if (a.count < 10 || !counts.some((c) => c > 0)) continue;
    const fig = new Figure(grid, {
      title: `Rates over the neurons: ${a.population}`, subtitle: `${compact(a.count)} neurons`, height: 190,
      x: { label: 'Hz' }, y: { label: 'neurons' },
      hover: (r) => { const k = edges.findIndex((e, i) => i < counts.length && r >= e && r < edges[i + 1]); return k < 0 ? null : { title: `${short(edges[k])}–${short(edges[k + 1])} Hz`, rows: [{ colour: colourFor(results, a.population), value: String(counts[k]), label: 'neurons' }] }; },
    });
    fig.set([edges[0], edges[edges.length - 1]], [0, Math.max(...counts) * 1.08], (ctx, x, y) => bars(ctx, x, y, edges, counts, colourFor(results, a.population)));
    figures.push(fig);
  }
}

function nearest(xs, v) {
  let best = 0;
  for (let i = 1; i < xs.length; i++) if (Math.abs(xs[i] - v) < Math.abs(xs[best] - v)) best = i;
  return best;
}

// A quantity as written in a specification, in SI units (the units the editor's examples use).
const SI = { p: 1e-12, n: 1e-9, u: 1e-6, 'µ': 1e-6, m: 1e-3, '': 1, k: 1e3, M: 1e6 };
export function toSI(q) {
  if (typeof q === 'number') return q;
  const v = parseField(q);
  if (typeof v === 'number') return v;
  const m = String(v ?? '').match(/^([-+]?[\d.]+(?:e[-+]?\d+)?)\s*([pnuµmkM]?)(A|V|S|F|s|Hz)$/);
  return m ? Number(m[1]) * SI[m[2]] : NaN;
}

// ---- states ------------------------------------------------------------------------------------------------
export function states(root, results, spec, figures) {
  clear(root);
  if (!results.states.length) { root.append(el('p', { class: 'note', text: 'No state was recorded: add a state recorder.' })); return; }
  for (const s of results.states) {
    const nc = s.components.length;
    s.components.forEach((c, j) => {
      const series = s.neurons.map((n, k) => s.times.map((_, i) => s.values[k][i * nc + j]));
      const [f, label] = display(s.units[j], series.flat());
      const ys = series.map((v) => v.map((x) => x * f));
      const colours = s.neurons.map((_, k) => (s.neurons.length === 1 ? colourFor(results, s.population) : POPULATION_COLOURS[k % POPULATION_COLOURS.length]));
      const fig = new Figure(root, {
        title: `${c} of ${s.population}`, subtitle: label ? `in ${label}` : 'dimensionless', height: 190,
        legend: s.neurons.length > 1 ? s.neurons.slice(0, 8).map((n, k) => ({ label: `neuron ${n}`, colour: colours[k] })) : [],
        x: { label: 'time (ms)', format: ms }, y: { label },
        hover: (t) => { const i = nearest(s.times, t); return { title: `${ms(s.times[i])} ms`, rows: ys.slice(0, 8).map((v, k) => ({ colour: colours[k], value: `${short(v[i])}${label ? ` ${label}` : ''}`, label: `neuron ${s.neurons[k]}` })) }; },
      });
      fig.set([s.times[0] ?? 0, s.times[s.times.length - 1] ?? 1], extent(ys.slice(0, 8).flat()), (ctx, x, y) => {
        ys.slice(0, 8).forEach((v, k) => line(ctx, x, y, s.times, v, colours[k], k === 0 ? 2 : 1.25, k === 0 ? 1 : 0.8));
      });
      figures.push(fig);
    });
  }
}

// ---- weights -------------------------------------------------------------------------------------------------
export function weights(root, results, spec, figures, catalogue) {
  clear(root);
  if (!results.weights.length) { root.append(el('p', { class: 'note', text: 'No weights were recorded: add a weight recorder on a plastic projection.' })); return; }
  for (const w of results.weights) {
    const p = spec.projections[w.projection];
    const model = catalogue?.get(spec.populations.find((q) => q.name === p?.target)?.model);
    const unit = model?.receptors.find((r) => r.name === p?.receptor)?.unit ?? '1';
    const [f, label] = display(unit, w.mean);
    const mean = w.mean.map((v) => v * f), lo = w.mean.map((v, i) => (v - w.sd[i]) * f), hi = w.mean.map((v, i) => (v + w.sd[i]) * f);
    const colour = colourFor(results, p?.source);
    const fig = new Figure(root, {
      title: `Weights of ${p ? `${p.source} → ${p.target}` : `projection ${w.projection + 1}`}`, subtitle: `mean and one standard deviation${p?.synapse?.model ? `, ${p.synapse.model} synapses` : ''}`, height: 230,
      x: { label: 'time (s)' }, y: { label },
      hover: (t) => { const i = nearest(w.times, t); return { title: `${short(w.times[i])} s`, rows: [{ colour, value: `${short(mean[i])} ${label}`, label: 'mean' }, { value: `± ${short(hi[i] - mean[i])} ${label}`, label: 'standard deviation' }] }; },
    });
    fig.set([0, results.time], extent([...lo, ...hi], { zero: true }), (ctx, x, y) => {
      band(ctx, x, y, w.times, lo, hi, colour);
      line(ctx, x, y, w.times, mean, colour);
    });
    figures.push(fig);
  }
}

// ---- profile ----------------------------------------------------------------------------------------------------
const PHASES = [
  ['update', 'update of the neurons', 'devices, flows and the recording of states, in parallel'],
  ['exchange', 'exchange of spikes', 'the window\'s spikes gathered and sorted, on one thread'],
  ['deliver', 'delivery to synapses', 'every spike to its targets\' rings, with plasticity, in parallel'],
  ['record', 'recording', 'spikes and weights kept, on one thread'],
];
const PHASE_COLOURS = { update: POPULATION_COLOURS[0], exchange: POPULATION_COLOURS[3], deliver: POPULATION_COLOURS[2], record: POPULATION_COLOURS[1] };

export function profile(root, results) {
  clear(root);
  const st = results.statistics, p = st.profile;
  if (!p) { root.append(el('p', { class: 'note', text: 'This engine gives no profile: build it again.' })); return; }
  const steps = Math.round(results.time / results.resolution);
  const tiles = [
    ['computed in', duration(st.run_seconds)], ['real-time factor', short(st.realtime_factor)],
    ['synaptic events / s', compact(st.events / st.run_seconds)], ['spikes / s', compact(st.spikes / st.run_seconds)],
    ['windows', `${compact(p.windows)} of ${short(steps / Math.max(1, p.windows))} steps`], ['threads', String(st.threads)],
  ];
  root.append(el('div', { class: 'tiles' }, tiles.map(([l, v]) => el('div', { class: 'tile' }, el('div', { class: 'tile-label', text: l }), el('div', { class: 'tile-value', text: v })))));
  root.append(el('p', { class: 'note', text: `The engine advances in windows of the shortest delay (${st.min_delay_steps} steps): within one, no neuron can be reached by a spike emitted in it, so the threads update their own neurons without waiting; between windows the spikes are exchanged and delivered. Building the network took ${duration(st.build_seconds)} before the run.` }));

  // The run's time, phase by phase.
  const total = PHASES.reduce((s, [k]) => s + p[k], 0);
  const strip = el('div', { class: 'strip' });
  for (const [k, name] of PHASES) {
    const share = p[k] / total;
    strip.append(el('div', { class: 'seg', style: `flex:${Math.max(share, 0.0001)};background:${PHASE_COLOURS[k]}`, title: `${name}: ${duration(p[k])} (${short(share * 100)} %)` }));
  }
  root.append(el('h3', { text: 'Where the time went' }), strip,
    el('table', { class: 'data' }, el('tbody', {}, PHASES.map(([k, name, about]) => el('tr', {},
      el('td', {}, el('i', { class: 'swatch', style: `background:${PHASE_COLOURS[k]}` }), name),
      el('td', { text: duration(p[k]) }), el('td', { text: `${short((p[k] / total) * 100)} %` }), el('td', { class: 'faint', text: about }))))));

  // Thread by thread: busy, and waiting for the slowest.
  root.append(el('h3', { text: 'Thread by thread' }),
    el('div', { class: 'legend' }, ...[['update', PHASE_COLOURS.update], ['delivery', PHASE_COLOURS.deliver], ['waiting for the slowest thread', WAIT]].map(([l, c]) => el('span', { class: 'key' }, el('i', { class: 'k-box', style: `background:${c}` }), l))));
  const longest = p.update + p.deliver;
  const rows = el('div', { class: 'threads' });
  p.threads.forEach((t, i) => {
    const parts = [[t.update, PHASE_COLOURS.update, 'update'], [p.update - t.update, WAIT, 'waiting in the update'], [t.deliver, PHASE_COLOURS.deliver, 'delivery'], [p.deliver - t.deliver, WAIT, 'waiting in the delivery']];
    rows.append(el('div', { class: 'thread' }, el('span', { class: 'thread-name', text: `thread ${i + 1}` }),
      el('div', { class: 'bar' }, parts.map(([v, c, n]) => el('div', { class: 'seg', style: `flex:${Math.max(v, 0) / longest};background:${c}`, title: `${n}: ${duration(Math.max(v, 0))}` }))),
      el('span', { class: 'faint', text: `${short((1 - (t.update + t.deliver) / longest) * 100)} % waiting` })));
  });
  root.append(rows);
  const busy = p.threads.map((t) => t.update + t.deliver);
  const balance = Math.min(...busy) / Math.max(...busy);
  const waiting = 1 - busy.reduce((a, b) => a + b, 0) / busy.length / longest;
  root.append(el('p', { class: 'note', text: `The least busy thread worked ${short(balance * 100)} % as long as the busiest. The neurons are divided among the threads in equal numbers, and each thread delivers the spikes that reach its own neurons.` }));
  // Waiting has two causes: windows too short for their meeting's cost, or threads busy for
  // unequal times (uneven work, or other programs sharing the processor).
  const perWindow = (longest / Math.max(1, p.windows)) * 1e6;
  if (st.threads > 1 && waiting > 0.4 && perWindow < 1000) {
    root.append(el('p', { class: 'advice', text: `The threads waited ${short(waiting * 100)} % of the time: each window took ${short(perWindow)} µs, and at that size the threads' meeting at its end costs as much as its work. Fewer threads, or a longer shortest delay (fewer windows), would serve this network better.` }));
  } else if (st.threads > 1 && waiting > 0.25) {
    root.append(el('p', { class: 'advice', text: `The threads waited ${short(waiting * 100)} % of the time, in windows of ${short(perWindow / 1000)} ms, long enough that their meeting costs little: they were busy for unequal times. Other programs sharing the processor do this, as do hyperthreads sharing a core; the same run on an idle machine, or on as many threads as there are cores, shows how much of it is the network's own.` }));
  }
}

// ---- network --------------------------------------------------------------------------------------------------------
export function network(root, spec, catalogue) {
  clear(root);
  const pops = spec.populations;
  if (!pops.length) return;
  const W = 760, H = Math.max(420, 220 + pops.length * 40), cx = W / 2, cy = H / 2;
  const R = pops.length === 1 ? 0 : Math.min(W * 0.62, H) / 2 - 100;
  // On a ring, the first population on the left when there are two, at the top otherwise.
  const start = pops.length === 2 ? Math.PI : -Math.PI / 2;
  const pos = new Map(pops.map((p, i) => {
    const a = start + (2 * Math.PI * i) / pops.length;
    return [p.name, { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a), r: 16 + 9 * Math.log10(Math.max(1, p.size)), colour: colourOf(i), p }];
  }));
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.setAttribute('class', 'diagram');
  const node = (tag, attrs, text) => {
    const n = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    if (text !== undefined) n.textContent = text;
    svg.append(n);
    return n;
  };
  const inhibitory = (pr) => /inhib/.test(pr.receptor ?? '') || (typeof pr.weight === 'string' && pr.weight.trim().startsWith('-')) || (typeof pr.weight === 'number' && pr.weight < 0);
  spec.projections.forEach((pr) => {
    const a = pos.get(pr.source), b = pos.get(pr.target);
    if (!a || !b) return;
    const inh = inhibitory(pr);
    const label = `${pr.source} → ${pr.target}: ${Object.keys(pr.rule ?? {}).filter((k) => !['autapses', 'multapses'].includes(k)).map((k) => `${k.replace(/_/g, ' ')} ${pr.rule[k] === true ? '' : pr.rule[k]}`).join(' ')}; weight ${typeof pr.weight === 'object' ? JSON.stringify(pr.weight) : pr.weight}; delay ${pr.delay ?? 'one step'}${pr.synapse?.model ? `; ${pr.synapse.model}` : ''}`;
    let path, end;
    if (a === b) {
      const out = pops.length === 1 ? Math.PI / 2 : Math.atan2(a.y - cy, a.x - cx);
      const ux = Math.cos(out), uy = Math.sin(out), s = a.r + 34;
      const p1 = [a.x + ux * a.r + -uy * 10, a.y + uy * a.r + ux * 10], p2 = [a.x + ux * a.r - -uy * 10, a.y + uy * a.r - ux * 10];
      path = `M${p1[0]},${p1[1]} C${a.x + ux * s - uy * 40},${a.y + uy * s + ux * 40} ${a.x + ux * s + uy * 40},${a.y + uy * s - ux * 40} ${p2[0]},${p2[1]}`;
      end = p2;
    } else {
      const dx = b.x - a.x, dy = b.y - a.y, d = Math.hypot(dx, dy), ux = dx / d, uy = dy / d;
      const bend = 18, mx = (a.x + b.x) / 2 - uy * bend, my = (a.y + b.y) / 2 + ux * bend;
      const s = [a.x + ux * a.r, a.y + uy * a.r], e = [b.x - ux * (b.r + 6), b.y - uy * (b.r + 6)];
      path = `M${s[0]},${s[1]} Q${mx},${my} ${e[0]},${e[1]}`;
      end = e;
    }
    const g = node('path', { d: path, fill: 'none', stroke: a.colour, 'stroke-width': 1.6, opacity: 0.8 });
    node('title', {}, label);
    g.append(svg.lastChild);
    if (inh) node('circle', { cx: end[0], cy: end[1], r: 4.5, fill: a.colour });
    else {
      const t = document.createElementNS(NS, 'title');
      t.textContent = label;
      node('circle', { cx: end[0], cy: end[1], r: 2.5, fill: a.colour }).append(t);
    }
  });
  for (const [name, q] of pos) {
    node('circle', { cx: q.x, cy: q.y, r: q.r, fill: '#281d2c', stroke: q.colour, 'stroke-width': 2 });
    // The labels inside the circle when they fit; beneath it when they do not.
    const sub = `${compact(q.p.size)} · ${q.p.model}`;
    const inside = 2 * q.r - 8 >= 6.2 * Math.max(sub.length, name.length * 1.4);
    node('text', { x: q.x, y: inside ? q.y - 2 : q.y + q.r + 18, 'text-anchor': 'middle', class: 'd-name' }, name);
    node('text', { x: q.x, y: inside ? q.y + 14 : q.y + q.r + 33, 'text-anchor': 'middle', class: 'd-sub' }, sub);
    const devs = spec.devices.filter((d) => d.target === name);
    devs.forEach((d, k) => {
      // Devices come in across the population's outward direction, clear of its own loop.
      // From above, unless the population's own loop points up; then across its outward direction.
      const out = pops.length === 1 ? Math.PI / 2 : Math.atan2(q.y - cy, q.x - cx);
      const up = Math.abs(Math.sin(out) + 1) < 0.3;
      const ang = (up ? out - Math.PI / 2 : -Math.PI / 2) + (k - (devs.length - 1) / 2) * 0.5;
      const dx = Math.cos(ang), dy = Math.sin(ang), d0 = q.r + 46;
      const shown = (q, u) => (Number.isFinite(toSI(q)) ? format(toSI(q), u) : String(q ?? ''));
      const what = d.model === 'poisson' ? `Poisson ${shown(d.rate, 'Hz')}` : d.model === 'dc' ? `current ${shown(d.amplitude, 'A')}` : d.model === 'sine' ? `sine ${shown(d.frequency, 'Hz')}` : d.model === 'spike_times' ? `${(d.times ?? []).length} spikes at times` : d.model.replace(/_/g, ' ');
      node('line', { x1: q.x + dx * d0, y1: q.y + dy * d0, x2: q.x + dx * (q.r + 5), y2: q.y + dy * (q.r + 5), stroke: 'rgba(240,226,238,0.35)', 'stroke-width': 1, 'stroke-dasharray': '3 3' });
      node('text', { x: q.x + dx * (d0 + 6), y: q.y + dy * (d0 + 6) + 4, 'text-anchor': dx > 0.3 ? 'start' : dx < -0.3 ? 'end' : 'middle', class: 'd-device' }, what);
    });
  }
  root.append(el('p', { class: 'note', text: 'Each population is a circle, its area growing with its size; a projection is drawn in its source\'s colour and ends in a point (excitatory) or a dot (inhibitory). Hover a projection for its rule, weight and delay; dashed lines are devices.' }), svg);
}
