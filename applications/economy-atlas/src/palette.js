// The atlas's colours: a diverging ramp from slate (value sent out, a negative curvature) through
// paper to terracotta (value drawn in, a positive curvature), and the colours of the parts of a flow.

const hex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
const SLATE = hex('#3f6b7c'), PALE = hex('#f3ecdf'), TERRA = hex('#b9612e');

function mix(a, b, t) {
  return a.map((x, i) => Math.round(x + (b[i] - x) * t));
}

// t in [−1, 1]: −1 slate, 0 paper, +1 terracotta; the middle eased so that small values stay pale.
export function diverging(t) {
  const u = Math.max(-1, Math.min(1, t));
  const s = Math.sign(u) * Math.pow(Math.abs(u), 0.7);
  const c = s < 0 ? mix(PALE, SLATE, -s) : mix(PALE, TERRA, s);
  return `rgb(${c[0]}, ${c[1]}, ${c[2]})`;
}

export const rampCss = () => `linear-gradient(90deg, ${[-1, -0.5, -0.2, 0, 0.2, 0.5, 1].map(diverging).join(', ')})`;

export const PART = {
  reciprocal: '#d8c3a0',   // exchanged both ways
  driven: '#5a4b3d',       // carried by the potential, from where value is drawn to where it is drawn to
  local: '#b98a2e',        // circulating round the map's triangles
  global: '#c9a255',       // circulating round its holes, through long-distance links
};

// A categorical sequence for currencies and sectors, warm and gentle.
export const SERIES = ['#b9612e', '#3f6b7c', '#8a6a2f', '#7a4f6d', '#5d7f4f', '#a8743a', '#4c5f8c', '#9b4f3f', '#6d7a3a', '#866447',
  '#3e7a72', '#a05a7a', '#7b6a9b', '#b08a4a', '#5a6e7f', '#8f4e2c', '#4f8a6a', '#9a7a3a', '#6a5a8a', '#7f3f4f'];
export const seriesColour = (i) => SERIES[((i % SERIES.length) + SERIES.length) % SERIES.length];
