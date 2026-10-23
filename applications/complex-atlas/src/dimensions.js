// The emergent dimensions of a discrete space: the spectral dimension d_s(t) = −2 d ln P / d ln t
// from the return probability P(t) of a lazy random walk, and the Hausdorff dimension
// d_H(r) = d ln |B(r)| / d ln r from the growth of balls, on the complex's 1-skeleton or on its
// dual graph (physis spacetime/dimensions.hpp). On a lattice both are the lattice's dimension at
// intermediate scales; on a random planar quadrangulation (pure 2D quantum gravity) d_H = 4 and
// d_s = 2; on branched polymers (melonic graphs) d_H = 2 and d_s = 4/3 (Ambjørn, Durhuus and
// Jonsson 1997; Le Gall 2013; Gurau and Ryan 2014). Finite graphs show these only between their
// smallest and largest scales.

import { el, clear } from './dom.js';
import { xyChart, chartLegend } from './plot.js';
import { number } from './format.js';

export const dimensionsForm = { graph: 'skeleton' };

// expected: {ds, dh, label} marked on the charts, or null.
export function renderDimensions(host, { complex, result, running, error, onRun, expected = null, compact = false }) {
  clear(host);
  if (!complex) return;
  if (!compact) host.append(el('h4', { text: 'Emergent dimensions' }));
  host.append(el('p', { class: 'note', text: 'How many dimensions the complex has as a space: by diffusion (the spectral dimension, from the return probability of a random walk) and by counting (the Hausdorff dimension, from the growth of balls).' }));
  const graph = el('select', { onchange: () => { dimensionsForm.graph = graph.value; } },
    ...[['skeleton', 'on its 1-skeleton'], ['dual', 'on its dual graph']].map(([v, t]) => el('option', { value: v, selected: dimensionsForm.graph === v ? true : null }, t)));
  host.append(el('div', { class: 'row' }, graph, el('button', { class: 'button', disabled: running ? true : null, onclick: () => onRun(dimensionsForm.graph) }, running ? 'Measuring…' : 'Measure'),
    error ? el('span', { class: 'error', text: error }) : null));
  if (!result) return;
  const s = result.spectral, h = result.hausdorff;
  if (s?.t?.length) {
    host.append(xyChart([{ x: s.t, y: s.ds, colour: 'rose', label: 'd_s' }], {
      xlabel: 'diffusion time t', ylabel: 'd_s(t)', logx: true, ylim: [0, Math.max(4.5, ...s.ds.filter(Number.isFinite)) ],
      hlines: expected?.ds ? [{ y: expected.ds, label: `d_s = ${number(expected.ds, 3)} ${expected.label || ''}` }] : [],
    }));
  }
  if (h?.r?.length) {
    host.append(xyChart([{ x: h.r, y: h.dh, colour: 'lilac', label: 'd_H' }], {
      xlabel: 'radius r', ylabel: 'd_H(r)', logx: true, ylim: [0, Math.max(4.5, ...h.dh.filter(Number.isFinite))],
      hlines: expected?.dh ? [{ y: expected.dh, label: `d_H = ${number(expected.dh, 3)} ${expected.label || ''}` }] : [],
    }));
  }
  host.append(chartLegend([{ label: 'spectral, from diffusion', colour: 'rose' }, { label: 'Hausdorff, from balls', colour: 'lilac' }, ...(expected ? [{ label: 'the continuum value', colour: 'gold', dash: '4 3' }] : [])]));
  const peak = (arr) => Math.max(...arr.filter(Number.isFinite));
  host.append(el('p', { class: 'fact', text: `Plateau: d_s ≈ ${number(peak(s?.ds || [NaN]), 3)}, d_H ≈ ${number(peak(h?.dh || [NaN]), 3)} (the largest values between the smallest and largest scales).` }));
}
