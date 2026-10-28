// The geographic tile pyramid of the globes: which part of a body each tile covers.
//
// Level z has 2^(z+1) columns and 2^z rows of tiles, from 180 W and 90 N, in east longitude and
// latitude (geodetic on the Earth, planetocentric on the Moon): NASA Trek's scheme, and the data
// library's (data/tiles/pyramid.py, whose functions these are, held equal by the checks). A tile's
// four children are the quarters of it; its parent the tile of the level above that contains it.
// A terrain tile carries 65 x 65 samples at the corners of its 64 x 64 cells, so that neighbours
// sample the same points along their shared edge.

export const GRID = 65;

export const key = (z, x, y) => `${z}/${x}/${y}`;

// [west, south, east, north] in degrees.
export function tileBounds(z, x, y) {
  const span = 180 / 2 ** z;
  const west = -180 + x * span, north = 90 - y * span;
  return [west, north - span, west + span, north];
}

export function children(z, x, y) {
  return [[z + 1, 2 * x, 2 * y], [z + 1, 2 * x + 1, 2 * y], [z + 1, 2 * x, 2 * y + 1], [z + 1, 2 * x + 1, 2 * y + 1]];
}
export function parent(z, x, y) {
  return z === 0 ? null : [z - 1, x >> 1, y >> 1];
}
export const roots = () => [[0, 0, 0], [0, 1, 0]];

// The tile of level z containing a point (degrees).
export function tileAt(z, lonDeg, latDeg) {
  const span = 180 / 2 ** z;
  const x = Math.min(2 ** (z + 1) - 1, Math.max(0, Math.floor((((lonDeg + 180) % 360) + 360) % 360 / span)));
  const y = Math.min(2 ** z - 1, Math.max(0, Math.floor((90 - latDeg) / span)));
  return [z, x, y];
}

// The angular size of a tile's cell, degrees.
export const cellDegrees = (z, n = GRID) => 180 / 2 ** z / (n - 1);

// Where an ancestor's image lies inside it: the part [u0, v0, size] of the ancestor's texture
// that covers the tile (u to the east, v to the south), for drawing a tile with the image of an
// ancestor when its own is not loaded or does not exist.
export function subregion(z, x, y, za) {
  const d = z - za;
  if (d < 0) throw new Error('subregion: the ancestor is below the tile');
  const size = 1 / 2 ** d;
  return [(x - ((x >> d) << d)) * size, (y - ((y >> d) << d)) * size, size];
}
