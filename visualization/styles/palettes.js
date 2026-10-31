// The colour maps of the warm family the infrastructure's pictures share (ADR 0005; the infrastructure's
// ADR 0053; CLAUDE.md §7, visual style), each made in OKLab (colour.js) and validated for colour-vision
// deficiencies by the checks (cvd.js).
//
// The kind of a map is decided by the reduction it draws (formalism Definition 3.10, its table of
// fibres): sequential for a magnitude, one hue from light to dark, its lightness monotone; diverging
// for a signed value about a natural zero, two hues apart through a light, neutral middle; cyclic for
// an angle; categorical for a label, never ordered.
//
//   plum           sequential   plum night, Old Rose, cream: magnitudes on a night ground
//   gold           sequential   pale gold to umber: the Economy Atlas's gold; the Moon's crust
//   orchid         sequential   pale orchid to deep violet: Cosmopoiesis's orchid
//   lavender-rust  diverging    lavender through a pale rose to rust
//   orchid-gold    diverging    violet inward, gold outward, through a pale middle
//   sea            sequential   sea foam to deep sea glass: a level of water, a crowd read as a sea
//   sea-brick      diverging    sea glass below, brick rose above, through a pale stone: a level against its tide
//   warm           cyclic       plum, Old Rose, gold, blush lilac, and round again
//   categorical    categorical  ten warm hues, in the order the Complex Atlas has always drawn its indices
//   categorical-apart            the same ten, ordered so that neighbours stay apart in every vision
//
// The ten hues are pastel and cannot all be told apart by a reader with a deficiency: under
// deuteranopia two of them come within 0.003 of each other in OKLab. In the order `categorical-apart`
// each hue differs from the next by at least 0.158 in normal vision and 0.071 under each deficiency
// (the checks measure it, against floors of 0.15 and 0.06): of all 10! orders, the one whose least
// difference, as a fraction of its floor, is greatest. Identity must still be carried by a label or a
// legend as well as by hue.

import { gradientCss, hexRgb, ramp, table } from './colour.js';

export const MAPS = {
  plum: { kind: 'sequential', stops: ['#2a1f3a', '#6b4a86', '#c08081', '#f0d9b5'], about: 'plum night to Old Rose to cream' },
  gold: { kind: 'sequential', stops: ['#f6e3b4', '#e9aa40', '#965208'], about: 'pale gold to umber' },
  orchid: { kind: 'sequential', stops: ['#ecdef7', '#b892dd', '#5e3a8b'], about: 'pale orchid to deep violet' },
  'lavender-rust': { kind: 'diverging', stops: ['#5d5aa0', '#a99ad0', '#f2e6ee', '#d38e8a', '#9c4a45'], about: 'lavender through a pale rose to rust' },
  'orchid-gold': { kind: 'diverging', stops: ['#5e3a8b', '#b892dd', '#f4ede4', '#e9aa40', '#965208'], about: 'violet below the middle, gold above it' },
  sea: { kind: 'sequential', stops: ['#e4f0ea', '#9cc7bf', '#4f8a8b', '#22505a'], about: 'sea foam to deep sea glass' },
  'sea-brick': { kind: 'diverging', stops: ['#2e6b70', '#8fbdb6', '#f1eae3', '#d69a8e', '#9a4a3f'], about: 'sea glass below the middle, brick rose above it' },
  warm: { kind: 'cyclic', stops: ['#6b4a86', '#c08081', '#e9b872', '#d9b8d6'], about: 'plum, Old Rose, gold, blush lilac, and round again' },
  categorical: {
    kind: 'categorical',
    stops: ['#e0a3b6', '#e9b872', '#9fc2b8', '#b9a2e0', '#e48f78', '#c9d38a', '#8fb3dd', '#d987c6', '#f0cf9a', '#a6d6c4'],
    about: 'ten warm hues',
  },
  'categorical-apart': {
    kind: 'categorical',
    stops: ['#e0a3b6', '#c9d38a', '#e48f78', '#a6d6c4', '#b9a2e0', '#e9b872', '#8fb3dd', '#f0cf9a', '#d987c6', '#9fc2b8'],
    about: 'ten warm hues, neighbours apart in every vision',
  },
};

// The default map of each kind.
export const DEFAULT_MAP = { sequential: 'plum', diverging: 'lavender-rust', cyclic: 'warm', categorical: 'categorical-apart' };

// A map by name: colour(u) for u in [0, 1] (for a diverging map u = (t + 1)/2 of t in [−1, 1], its
// middle at u = ½; for a categorical map the hue of index ⌊u⌋, cycled), and its table and gradient.
export function colourMap(name) {
  const m = MAPS[name];
  if (!m) throw new Error(`no colour map "${name}": one of ${Object.keys(MAPS).join(', ')}`);
  const colour = m.kind === 'categorical'
    ? (i) => hexRgb(m.stops[((Math.floor(i) % m.stops.length) + m.stops.length) % m.stops.length])
    : ramp(m.stops, { cyclic: m.kind === 'cyclic' });
  return {
    name, kind: m.kind, about: m.about, stops: m.stops, colour,
    at: m.kind === 'diverging' ? (t) => colour((t + 1) / 2) : colour,
    table: (steps = 256, alpha) => table(colour, steps, alpha),
    gradient: (options) => gradientCss(colour, options),
  };
}

// A part of a sequential map, u in [0, 1] taking it from `from` to `to`: on a night ground the darkest
// part of the plum map is left out, on a day ground its palest, so that every value stands out.
export const within = (colour, from, to) => (u) => colour(from + (to - from) * Math.min(1, Math.max(0, u)));

// The maps the provisional JavaScript has always offered, kept under their names (js/complex.js
// re-exports them).
const plum = colourMap('plum').colour, lavender = colourMap('lavender-rust').colour, warm = colourMap('warm').colour;
const hues = colourMap('categorical').colour;
export const sequentialColour = plum;                                     // magnitudes in [0, 1]
export const divergingColour = (t) => lavender((t + 1) / 2);               // signed values in [−1, 1]
export const cyclicColour = (theta) => warm(theta / (2 * Math.PI));        // angles in radians, any period of 2π
export const categoricalColour = (i) => hues(i);                           // labels, by index
