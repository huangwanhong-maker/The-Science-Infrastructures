// The atlas's colours, read from the theme's tokens (styles.css) when drawing, so that the map follows
// day and night; outside a page (in the checks) the day's values. Every data colour was validated
// for colour vision deficiencies on its surface (the dataviz method: OKLab ΔE ≥ 8 between neighbours
// under protanopia and deuteranopia):
//   diverging, for a polarity: steel blue (value sent out) through a neutral to rose red (drawn in);
//   sequential, for a magnitude: the neutral to deep gold;
//   the parts of a flow: indigo (driven), rose (local, round the triangles), gold (global, round the holes).

const DAY = {
  sends: '#3A6EA5', draws: '#C4506E', neutral: '#EFE9DD', driven: '#4B4F8F', local: '#C4506E', global: '#A98A2C',
  reciprocal: '#E4D4A8', gold: '#C9A85C', 'gold-ink': '#7E621F', rose: '#C4506E', ink: '#2B2620', paper: '#FBF8F1',
};

let cache = null;
const listeners = new Set();

export function token(name) {
  if (!cache) {
    cache = { ...DAY };
    if (typeof document !== 'undefined') {
      const style = getComputedStyle(document.documentElement);
      for (const k of Object.keys(DAY)) {
        const v = style.getPropertyValue(`--${k}`).trim();
        if (v) cache[k] = v;
      }
    }
  }
  return cache[name];
}

// The theme changed: forget the colours read, and tell the views to draw again.
export function themeChanged() {
  cache = null;
  for (const f of listeners) f();
}
export const onTheme = (f) => listeners.add(f);

const rgb = (h) => {
  const s = h.replace('#', '');
  const x = s.length === 3 ? s.split('').map((c) => c + c).join('') : s;
  return [0, 2, 4].map((i) => parseInt(x.slice(i, i + 2), 16));
};
const mix = (a, b, t) => a.map((x, i) => Math.round(x + (b[i] - x) * t));
const css = (c) => `rgb(${c[0]}, ${c[1]}, ${c[2]})`;

// t in [−1, 1]: −1 value sent out, 0 neutral, +1 value drawn in; eased so that small values stay quiet.
export function diverging(t) {
  const u = Math.max(-1, Math.min(1, Number.isFinite(t) ? t : 0));
  const s = Math.sign(u) * Math.pow(Math.abs(u), 0.7);
  const mid = rgb(token('neutral'));
  return css(s < 0 ? mix(mid, rgb(token('sends')), -s) : mix(mid, rgb(token('draws')), s));
}

// t in [0, 1]: a magnitude, from the neutral to deep gold.
export function sequential(t) {
  const u = Math.max(0, Math.min(1, Number.isFinite(t) ? t : 0));
  return css(mix(rgb(token('neutral')), rgb(token('gold-ink')), Math.pow(u, 0.8)));
}

export const rampCss = (kind = 'diverging') => kind === 'sequential'
  ? `linear-gradient(90deg, ${[0, 0.25, 0.5, 0.75, 1].map(sequential).join(', ')})`
  : `linear-gradient(90deg, ${[-1, -0.5, -0.2, 0, 0.2, 0.5, 1].map(diverging).join(', ')})`;

// The parts of a flow, and the colour of what a plug-in adds.
export const PART = {
  get reciprocal() { return token('reciprocal'); },
  get driven() { return token('driven'); },
  get local() { return token('local'); },
  get global() { return token('global'); },
  get plugin() { return token('rose'); },
};

// A sequence for currencies, sectors and positions: identity only, muted, the chosen one drawn heavier.
export const SERIES = ['#A98A2C', '#3A6EA5', '#C4506E', '#4B4F8F', '#6E8B4E', '#9B6A3A', '#3E7A72', '#8C5A86', '#B07A3A', '#5A6E8F',
  '#A0504A', '#4F8A6A', '#7B6A9B', '#B09040', '#5A7A9A', '#8F4E3C', '#3F6E5A', '#9A7A3A', '#6A5A8A', '#7F4F5F'];
export const seriesColour = (i) => SERIES[((i % SERIES.length) + SERIES.length) % SERIES.length];
