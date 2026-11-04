// The colours of the body's systems, softened from anatomy's conventions (arteries red, veins blue, nerves
// yellow, bone ivory, muscle red-brown, lymphoid green), and what each shows at first: the body's interior
// within a glass skin, the muscles and fascia off so that the organs are seen.
//
// Identity is carried by the systems' names in the rail and on the card as well as by colour; the arteries'
// red and the veins' blue stay apart in normal vision and under each deficiency, by the visualization layer's
// simulation (styles/cvd.js; the atlas's checks), and the veins are the darker.

export const SYSTEM_STYLE = {
  skin: { colour: '#e4c5b0', show: 'glass' },
  skeleton: { colour: '#e8e0cc', show: 'solid' },
  muscles: { colour: '#b9665a', show: 'off' },
  fascia: { colour: '#d9d4c6', show: 'off' },
  heart: { colour: '#b2505b', show: 'solid' },
  arteries: { colour: '#cc4a52', show: 'solid' },
  veins: { colour: '#3d5a9e', show: 'solid' },
  brain: { colour: '#e6c4ba', show: 'solid' },
  nerves: { colour: '#e2bd52', show: 'solid' },
  eye: { colour: '#d6d3e6', show: 'solid' },
  ear: { colour: '#e0b6a2', show: 'solid' },
  respiratory: { colour: '#dc9ea4', show: 'solid' },
  digestive: { colour: '#d8a57a', show: 'solid' },
  urinary: { colour: '#a95f57', show: 'solid' },
  genital: { colour: '#c991a8', show: 'solid' },
  endocrine: { colour: '#a68cc8', show: 'solid' },
  lymphoid: { colour: '#8db07c', show: 'solid' },
};

export const UNKNOWN = { colour: '#b9bec5', show: 'solid' };

export function systemColour(key) {
  return (SYSTEM_STYLE[key] || UNKNOWN).colour;
}

// A colour a little darker, for outlines on the plate.
export function shade(hex, f = 0.72) {
  const n = parseInt(hex.slice(1), 16);
  const c = [(n >> 16) & 255, (n >> 8) & 255, n & 255].map((x) => Math.round(x * f));
  return `#${c.map((x) => x.toString(16).padStart(2, '0')).join('')}`;
}
