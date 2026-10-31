// The colours of a GIS's features, from the visualization layer's validated maps (styles/palettes.js): a
// category by its index in the categorical map whose neighbours stay apart in every vision; a magnitude on a
// sequential map over a declared domain; a signed value on a diverging map about zero. Colours are
// [r, g, b] in [0, 1], as the maps give them, so that a renderer can write them into its buffers.

import { colourMap } from '../../styles/index.js';

export function categorical(name = 'categorical-apart') {
  const m = colourMap(name).colour;
  return (i) => m(i);
}

export function sequential(lo, hi, name = 'plum', { from = 0.12, to = 1 } = {}) {
  const m = colourMap(name).colour;
  return (v) => (Number.isFinite(v) ? m(from + (to - from) * Math.min(1, Math.max(0, (v - lo) / ((hi - lo) || 1)))) : null);
}

export function diverging(span, name = 'lavender-rust') {
  const m = colourMap(name).at;
  return (v) => (Number.isFinite(v) ? m(Math.max(-1, Math.min(1, v / (span || 1)))) : null);
}

export const hexOf = (rgb) => `#${rgb.map((c) => Math.round(Math.min(1, Math.max(0, c)) * 255).toString(16).padStart(2, '0')).join('')}`;
export const rgbOfHex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);
