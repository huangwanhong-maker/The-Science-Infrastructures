// The workbench's colours, read from the style sheet's tokens so that the page and its drawings
// agree: the four bases, the kinds of features, and the ramps of the data.

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

export function baseColour(b) { return css(`--${'ACGT'.includes(b) ? b : 'N'}`); }
export const BASE_COLOURS = () => ({ A: css('--A'), C: css('--C'), G: css('--G'), T: css('--T'), N: css('--N') });

// A feature's colour by its key.
export function featureColour(key) {
  if (key === 'CDS') return css('--cds');
  if (key === 'rRNA') return css('--rrna');
  if (key === 'tRNA') return css('--trna');
  if (/RNA/.test(key)) return css('--ncrna');
  if (/repeat|mobile|insertion/.test(key)) return css('--mobile');
  return css('--other');
}
export const token = css;

// Interpolation between colours given as hex.
const hex = (c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16));
function mix(a, b, t) {
  const x = hex(a), y = hex(b);
  return `rgb(${x.map((v, i) => Math.round(v + (y[i] - v) * t)).join(',')})`;
}
function through(stops, t) {
  t = Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0));
  const k = Math.min(stops.length - 2, Math.floor(t * (stops.length - 1)));
  return mix(stops[k], stops[k + 1], t * (stops.length - 1) - k);
}
// A warm sequential ramp: paper, rose, mulberry, plum (heat maps of scores and counts).
export const warm = (t) => through(['#FBF4F4', '#EBC7CF', '#C9839F', '#A8577E', '#5E2E4A'], t);
// A diverging ramp: sage below, paper at the middle, mulberry above; t in [−1, 1].
export const diverging = (t) => through(['#4F7D5D', '#9CC0A6', '#FBF4F4', '#D29BB6', '#8C3F66'], (t + 1) / 2);
// Categories (populations, clusters, workers): distinct and gentle.
export const CATEGORIES = ['#A8577E', '#5D7FB2', '#C9A043', '#6E9C7A', '#8C7BB8', '#C55A72', '#B07D62', '#4E8C8F', '#9C8A4A', '#7A6B9E'];
