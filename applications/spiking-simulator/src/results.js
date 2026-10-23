// What a run gave: the engine's figures (simulated time, wall time, spikes, synaptic events),
// each population's recorded activity, and the results as a file.

import { el, clear } from './dom.js';
import { colourOf } from './palette.js';
import { short } from './units.js';

// A count in words people read: 1,284 / 12.9 thousand / 3.1 million.
export function count(n) {
  if (!Number.isFinite(n)) return '—';
  if (Math.abs(n) < 10000) return Math.round(n).toLocaleString('en');
  if (Math.abs(n) < 1e6) return `${short(n / 1e3)} thousand`;
  if (Math.abs(n) < 1e9) return `${short(n / 1e6)} million`;
  return `${short(n / 1e9)} billion`;
}

// A count in a tile: 1,284 / 12.9K / 3.1M.
export function compact(n) {
  if (!Number.isFinite(n)) return '—';
  if (Math.abs(n) < 10000) return Math.round(n).toLocaleString('en');
  if (Math.abs(n) < 1e6) return `${short(n / 1e3)}K`;
  if (Math.abs(n) < 1e9) return `${short(n / 1e6)}M`;
  return `${short(n / 1e9)}B`;
}

export function duration(s) {
  if (!Number.isFinite(s)) return '—';
  if (s < 1e-3) return `${short(s * 1e6)} µs`;
  if (s < 1) return `${short(s * 1e3)} ms`;
  if (s < 120) return `${short(s)} s`;
  return `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
}

// Each population's recorded spikes: the neurons recorded, the spikes, the mean rate.
export function activity(results) {
  const T = results.time;
  return results.populations.map((p, i) => {
    const recs = results.spikes.filter((s) => s.population === p.name);
    const neurons = recs.reduce((n, r) => n + r.count, 0);
    const spikes = recs.reduce((n, r) => n + r.times.length, 0);
    return { name: p.name, model: p.model, size: p.size, colour: colourOf(i), neurons, spikes, rate: neurons > 0 && T > 0 ? spikes / (neurons * T) : NaN };
  });
}

export function renderResults(root, results, { spec, onDownload } = {}) {
  clear(root);
  if (!results) {
    root.append(el('div', { class: 'empty' },
      el('p', { class: 'lead', text: 'Compose a network on the left, or begin from an example, and run it.' }),
      el('p', { class: 'note', text: 'The engine is the neuroscience domain\'s C++ simulator: exact where the dynamics are linear, the same spikes on any number of threads.' })));
    return;
  }
  const st = results.statistics;
  const tiles = [
    ['simulated', duration(results.time)],
    ['computed in', duration(st.run_seconds)],
    ['real-time factor', short(st.realtime_factor)],
    ['threads', String(st.threads)],
    ['neurons', compact(st.neurons)],
    ['synapses', compact(st.synapses)],
    ['spikes', compact(st.spikes)],
    ['synaptic events', compact(st.events)],
  ];
  root.append(el('div', { class: 'tiles' }, tiles.map(([label, value]) => el('div', { class: 'tile' }, el('div', { class: 'tile-label', text: label }), el('div', { class: 'tile-value', text: value })))));
  const eventsPerSecond = st.run_seconds > 0 ? st.events / st.run_seconds : NaN;
  root.append(el('p', { class: 'note' },
    `${count(st.spikes)} spikes and ${count(st.events)} synaptic events in ${duration(st.run_seconds)} of computing (${count(eventsPerSecond)} events a second), after ${duration(st.build_seconds)} building the network; the shortest delay is ${st.min_delay_steps} steps of ${short(results.resolution * 1e3)} ms.`));

  const rows = activity(results);
  root.append(el('h3', { text: 'Recorded spikes' }),
    el('table', { class: 'data' },
      el('thead', {}, el('tr', {}, ['population', 'model', 'neurons', 'recorded', 'spikes', 'mean rate'].map((h) => el('th', { text: h })))),
      el('tbody', {}, rows.map((r) => el('tr', {},
        el('td', {}, el('i', { class: 'swatch', style: `background:${r.colour}` }), r.name),
        el('td', { text: r.model }), el('td', { text: count(r.size) }), el('td', { text: r.neurons ? count(r.neurons) : '—' }),
        el('td', { text: r.neurons ? count(r.spikes) : '—' }), el('td', { text: r.neurons ? `${short(r.rate)} Hz` : '—' }))))));

  const other = [];
  for (const s of results.states) other.push(`the state of ${s.neurons.length} ${s.neurons.length === 1 ? 'neuron' : 'neurons'} of ${s.population} (${s.components.join(', ')}), ${count(s.times.length)} samples`);
  for (const r of results.rates) other.push(`the rate of ${r.population} in bins of ${short(r.bin * 1e3)} ms`);
  for (const w of results.weights) other.push(`the weights of projection ${w.projection + 1}, ${count(w.times.length)} samples`);
  if (other.length) root.append(el('h3', { text: 'Also recorded' }), el('ul', { class: 'plain' }, other.map((t) => el('li', { text: t }))));

  root.append(el('div', { class: 'actions' },
    el('button', { class: 'button', onclick: () => onDownload?.('results') }, 'Download the results'),
    el('button', { class: 'button', onclick: () => onDownload?.('spec') }, 'Download the specification')));
  if (spec) root.append(el('p', { class: 'note', text: 'Times are in seconds and every quantity in SI units; each recorded component names its unit.' }));
}
