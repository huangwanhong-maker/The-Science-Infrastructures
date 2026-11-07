// The colours of the scene, by theme and by the Sun: space and the graded Earth, the sky's zenith and horizon at the
// Sun's altitude, the haze, the ground beneath a city. The Earth's imagery is graded into the workbench's palette, its
// land mauve and its seas violet, so that what is laid on it (always in the visualization layer's validated maps) reads
// as data against a quiet ground. Colours are written in sRGB here and given to the shaders as linear light.

const lin = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
const rgb = (hex) => [1, 3, 5].map((i) => lin(parseInt(hex.slice(i, i + 2), 16) / 255));
const mix = (a, b, t) => a.map((v, i) => v + (b[i] - v) * t);
const smooth = (a, b, x) => { const t = Math.min(1, Math.max(0, (x - a) / (b - a))); return t * t * (3 - 2 * t); };

// The sky at the sine of the Sun's altitude s: stops of night, the Sun at the horizon, and full day.
function sky(stops) {
  const S = Object.fromEntries(Object.entries(stops).map(([k, v]) => [k, Object.fromEntries(Object.entries(v).map(([q, h]) => [q, rgb(h)]))]));
  return (s) => {
    const t = smooth(-0.22, 0.02, s), d = smooth(0.02, 0.42, s);
    const at = (q) => mix(mix(S.night[q], S.dusk[q], t), S.day[q], d);
    return { zenith: at('zenith'), horizon: at('horizon'), below: at('below'), sun: at('sun') };
  };
}

const THEMES = {
  night: {
    space: rgb('#0B0611'), grain: 1,
    grade: { amount: 0.94, shadow: rgb('#170E22'), middle: rgb('#5B3C6E'), light: rgb('#E6B9D2'), seaDeep: rgb('#0D0A22'), seaLight: rgb('#2B2355'), nightTint: rgb('#FFC2DA'), glint: [0.075, 0.04, 0.085] },
    air: { day: [0.46, 0.30, 0.82], dusk: [0.95, 0.42, 0.62] },
    sky: sky({
      night: { zenith: '#0B0611', horizon: '#1E1230', below: '#120A1C', sun: '#000000' },
      dusk: { zenith: '#261A49', horizon: '#E0799F', below: '#2C1A3B', sun: '#FFC9B8' },
      day: { zenith: '#47397E', horizon: '#D3AFDC', below: '#49365A', sun: '#FFF0E0' },
    }),
    tint: rgb('#1C1129'), tintAmount: 0.42,
    shadow: '#080410', shadowOpacity: 0.55, lightScale: 0.78,
  },
  day: {
    space: rgb('#F2E8F3'), grain: 0,
    grade: { amount: 0.94, shadow: rgb('#6E5480'), middle: rgb('#B59CBF'), light: rgb('#F8F0F6'), seaDeep: rgb('#B9B2E2'), seaLight: rgb('#DFDAF3'), nightTint: rgb('#C2477F'), glint: [0.30, 0.26, 0.34] },
    air: { day: [0.62, 0.44, 0.88], dusk: [0.92, 0.48, 0.66] },
    sky: sky({
      night: { zenith: '#2C1E46', horizon: '#57406F', below: '#3B2B4F', sun: '#000000' },
      dusk: { zenith: '#9C86CE', horizon: '#F4B9CF', below: '#C5AECB', sun: '#FFD8C8' },
      day: { zenith: '#CDBDEB', horizon: '#F7EBF6', below: '#E5D8E8', sun: '#FFFFFF' },
    }),
    tint: rgb('#EFE4EF'), tintAmount: 0.5,
    shadow: '#3C2650', shadowOpacity: 0.28, lightScale: 1,
  },
};

export const look = (theme) => THEMES[theme === 'day' ? 'day' : 'night'];

// A colour of the day ribbon at an altitude of the Sun (degrees), for the theme: the night's plum, the rose of the
// horizon, and day's pale lilac.
const RIBBON = {
  night: [[-90, '#140B1F'], [-18, '#190F27'], [-9, '#2C1B47'], [-3, '#7E4674'], [1, '#D66F96'], [6, '#EBA6BE'], [20, '#DCC3E2'], [65, '#F4E9F3']],
  day: [[-90, '#3B2A55'], [-18, '#443263'], [-9, '#5F4A83'], [-3, '#A36A98'], [1, '#E58DAE'], [6, '#F3BBCF'], [20, '#EBDDEF'], [65, '#FBF5FA']],
};
const hexTo = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
export function ribbonColour(theme, alt) {
  const stops = RIBBON[theme === 'day' ? 'day' : 'night'];
  for (let i = 1; i < stops.length; i++) {
    if (alt <= stops[i][0]) {
      const t = (alt - stops[i - 1][0]) / (stops[i][0] - stops[i - 1][0]);
      const a = hexTo(stops[i - 1][1]), b = hexTo(stops[i][1]);
      return `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * t)).join(',')})`;
    }
  }
  return stops[stops.length - 1][1];
}
