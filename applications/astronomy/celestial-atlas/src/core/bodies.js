// Solar-system bodies shown by the atlas: NAIF identifiers, sizes, masses,
// orientation of the rotation pole, and display colours.
//
// Sources:
//   Archinal et al. (2018), IAU WGCCRE report: mean radii and pole directions (J2000).
//   Park et al. (2021), DE440: GM values (km^3/s^2), planetary systems.
//
// In DE440 the outer planets appear as barycentres of their systems
// (Mars = 4, ..., Pluto = 9). The barycentre lies within a few hundred km of
// the planet's centre, which is negligible at the scales shown.

export const GM_SUN = 132712440041.279419;

export const BODIES = [
  { key: 'sun', name: 'Sun', id: 10, radius: 695700, gm: GM_SUN,
    pole: [286.13, 63.87], colour: '#ffcf7a' },
  { key: 'mercury', name: 'Mercury', id: 199, radius: 2439.4, gm: 22031.868551,
    pole: [281.0103, 61.4155], colour: '#b9a89c' },
  { key: 'venus', name: 'Venus', id: 299, radius: 6051.8, gm: 324858.592,
    pole: [272.76, 67.16], colour: '#ecd2a2' },
  { key: 'earth', name: 'Earth', id: 399, radius: 6371.0, gm: 403503.235502,
    pole: [0, 90], colour: '#7fa7c9' },
  { key: 'moon', name: 'Moon', id: 301, radius: 1737.4, gm: 4902.800118,
    pole: [269.9949, 66.5392], colour: '#cfc6bd', parent: 'earth' },
  { key: 'mars', name: 'Mars', id: 4, radius: 3389.5, gm: 42828.375816,
    pole: [317.269, 54.432], colour: '#d8775a' },
  { key: 'jupiter', name: 'Jupiter', id: 5, radius: 69911, gm: 126712764.1,
    pole: [268.057, 64.495], colour: '#d9b48f' },
  { key: 'saturn', name: 'Saturn', id: 6, radius: 58232, gm: 37940584.8418,
    pole: [40.589, 83.537], colour: '#e6cc93',
    rings: { inner: 74500, outer: 136780 } },
  { key: 'uranus', name: 'Uranus', id: 7, radius: 25362, gm: 5794556.4,
    pole: [257.311, -15.175], colour: '#a9d3d4' },
  { key: 'neptune', name: 'Neptune', id: 8, radius: 24622, gm: 6836527.10058,
    pole: [299.36, 43.46], colour: '#7f98d6' },
  { key: 'pluto', name: 'Pluto', id: 9, radius: 1188.3, gm: 975.5,
    pole: [132.993, -6.163], colour: '#c9ab8f' },
];

export const BY_KEY = Object.fromEntries(BODIES.map((b) => [b.key, b]));

// Bodies that orbit the Sun (everything except the Sun and the Moon).
export const PLANETS = BODIES.filter((b) => b.key !== 'sun' && !b.parent);

// Gravitational parameter of the two-body problem for a body's orbit.
export function orbitalMu(body) {
  if (body.parent) return BY_KEY[body.parent].gm + body.gm;
  return GM_SUN + body.gm;
}
