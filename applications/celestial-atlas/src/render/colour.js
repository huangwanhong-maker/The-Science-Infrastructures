// Colours: the atlas palette, and the colour of a star from its B-V index.
//
// The colour of a star is drawn by the shared visualization layer (visualization/js/stars.js).

// Dusty Orchid (#CA8DBE) as the main colour: plum nights with orchid accents;
// orchid-tinted paper with plum ink for the charts.
export const PALETTE = {
  night: '#18121a',
  text: '#f6eef5',
  textDim: '#cdb6c8',
  accent: '#CA8DBE',
  blush: '#e8c3dc',
  accentDeep: '#8e4f82',
  roseDeep: '#a7709b',
  gold: '#f6c47e',
  goldDeep: '#c98a2e',
  paper: '#fbf6fa',
  ink: '#3a2a3a',
  inkSoft: '#86707f',
};

// rgba() of a palette colour.
export function rgba(hex, alpha) {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${n >> 16 & 255},${n >> 8 & 255},${n & 255},${alpha})`;
}

// The colour of a star from its colour index is the shared visualization layer's.
export { bvToTemperature, temperatureToRgb, starRgb, bpRpToTemperature, starColour } from '../../../../visualization/js/index.js';

// Linear interpolation through colour stops [[t, '#rrggbb'], ...] sorted by t.
export function ramp(stops, t) {
  if (!Number.isFinite(t)) return hexToRgb(stops[0][1]);
  if (t <= stops[0][0]) return hexToRgb(stops[0][1]);
  for (let k = 1; k < stops.length; k++) {
    if (t <= stops[k][0]) {
      const [t0, c0] = stops[k - 1], [t1, c1] = stops[k];
      const a = hexToRgb(c0), b = hexToRgb(c1), f = (t - t0) / (t1 - t0);
      return a.map((x, i) => x + (b[i] - x) * f);
    }
  }
  return hexToRgb(stops[stops.length - 1][1]);
}

// Galaxy colours. Tidal index: field galaxies lavender-blue, bound galaxies
// orchid, deepening with the strength of the tie. Group size: single galaxies
// lavender-blue, groups orchid, rich clusters gold. Morphology (de Vaucouleurs
// T): early types gold, spirals orchid, irregulars blue.
export const TIDAL_STOPS = [[-2, '#9fb0ff'], [0, '#eadcf0'], [1.5, '#dca4cf'], [3.5, '#c06aae']];
export const TYPE_STOPS = [[-5, '#f6b26b'], [0, '#f6d29a'], [3, '#f0d6e8'], [7, '#dca4cf'], [9, '#a9c4ff'], [11, '#8fb3ff']];
export const GROUP_STOPS = [[0, '#9fb0ff'], [0.5, '#eadcf0'], [1.2, '#dca4cf'], [2, '#f6c47e']];

export function hexToRgb(hex) {
  const n = parseInt(hex.slice(1), 16);
  return [(n >> 16 & 255) / 255, (n >> 8 & 255) / 255, (n & 255) / 255];
}
