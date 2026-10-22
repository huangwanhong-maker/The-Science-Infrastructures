// Colours: the atlas palette, and the colour of a star from its B-V index.
//
// Sources:
//   Ballesteros (2012), EPL 97, 34008: effective temperature from B-V.
//   Helland (2012), "How to convert temperature (K) to RGB": blackbody fit.

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

export function bvToTemperature(bv) {
  return 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62));
}

const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));

// Approximate sRGB of a blackbody of temperature T (K), components in [0, 1].
export function temperatureToRgb(T) {
  const t = clamp(T, 1000, 40000) / 100;
  const r = t <= 66 ? 255 : 329.698727446 * (t - 60) ** -0.1332047592;
  const g = t <= 66 ? 99.4708025861 * Math.log(t) - 161.1195681661
    : 288.1221695283 * (t - 60) ** -0.0755148492;
  const b = t >= 66 ? 255 : t <= 19 ? 0 : 138.5177312231 * Math.log(t - 10) - 305.0447927307;
  return [clamp(r, 0, 255) / 255, clamp(g, 0, 255) / 255, clamp(b, 0, 255) / 255];
}

// Star colour from B-V; unknown indices are drawn as a Sun-like star.
export function starRgb(bv) {
  const index = Number.isFinite(bv) ? clamp(bv, -0.4, 2.0) : 0.65;
  return temperatureToRgb(bvToTemperature(index));
}

// Effective temperature from Gaia BP-RP, for dwarfs, without reddening
// correction (Mucciarelli & Bellazzini 2020, RNAAS 4, 52; solar metallicity).
// The quadratic is calibrated for 0.39 < BP-RP < 1.5 and extrapolated outside.
export function bpRpToTemperature(bpRp) {
  const c = clamp(bpRp, -0.5, 5);
  return 5040 / (0.4929 + 0.5092 * c - 0.0353 * c * c);
}

// Star colour for a catalogue's photometric system: 'johnson' (B-V) or 'gaia' (BP-RP).
export function starColour(photometry, index) {
  if (photometry === 'gaia') {
    return temperatureToRgb(Number.isFinite(index) ? bpRpToTemperature(index) : 5800);
  }
  return starRgb(index);
}

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
