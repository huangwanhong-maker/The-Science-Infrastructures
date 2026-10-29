// Colouring and drapes: the part of a renderer that is arithmetic, without a device (ADR 0005; the
// infrastructure's ADR 0053). A value becomes a colour through its expression's map and scale; a field
// sampled on a sphere or a plane becomes an image whose pixels are its cells.
//
// Images are RGBA bytes, row by row from the top. A sphere's image is its grid as it lies: row i the
// latitude 90° − (i + ½)·180°/nlat from the north, column j the longitude −180° + (j + ½)·360°/nlon from
// the west, so that it is an equirectangular map, north up. A plane's image has its v axis up: row r is
// the cells j = nv − 1 − r. A cell without a value is transparent, never a colour (formalism §3.2).

import { colourMap, within } from '../styles/palettes.js';

// The colour of a value under an expression and a scale: [r, g, b, a] in [0, 1], or null where the value
// has no place. With opacity 'rising', a sequential map rises from clear at its low end, and a diverging
// map from clear at its middle, so that a weak field lets the ground show through.
export function colouring(expression, scale, { alpha = 1, range = null } = {}) {
  const map = colourMap(expression.colour.map);
  const rising = expression.colour.opacity === 'rising';
  const diverging = map.kind === 'diverging';
  const colour = range ? within(map.colour, range[0], range[1]) : map.colour;
  return (v) => {
    const u = scale.u(v);
    if (!Number.isFinite(u)) return null;
    const c = colour(u);
    const a = rising ? alpha * Math.sqrt(diverging ? Math.abs(2 * u - 1) : u) : alpha;
    return [c[0], c[1], c[2], a];
  };
}

// The image of scalar values on a sphere or a plane sampling: { width, height, data }.
export function drapeImage(values, sampling, colourOf) {
  let width, height, pixel;
  if (sampling.kind === 'sphere') {
    width = sampling.nlon;
    height = sampling.nlat;
    pixel = (index) => index;
  } else if (sampling.kind === 'plane') {
    width = sampling.nu;
    height = sampling.nv;
    pixel = (index) => {
      const i = index % width, j = Math.floor(index / width);
      return (height - 1 - j) * width + i;
    };
  } else {
    throw new Error(`a drape lies on a sphere or a plane, not a ${sampling.kind}`);
  }
  if (values.length !== width * height) throw new Error(`a drape of ${width} × ${height} cells has ${values.length} values`);
  const data = new Uint8ClampedArray(4 * width * height);
  for (let index = 0; index < values.length; index++) {
    const c = colourOf(values[index], index);
    if (!c) continue;
    const p = 4 * pixel(index);
    data[p] = Math.round(255 * c[0]);
    data[p + 1] = Math.round(255 * c[1]);
    data[p + 2] = Math.round(255 * c[2]);
    data[p + 3] = Math.round(255 * c[3]);
  }
  return { width, height, data };
}

// Curves as line segments with a colour at each end, for a renderer's lines: the positions (x, y, z
// interleaved, two points a segment) and the colours (r, g, b, a). A segment with an end that has no
// colour (no value there) is left out.
export function curveSegments(curves, colourOf) {
  let count = 0;
  for (const c of curves) if (c.seeded) count += Math.max(0, c.points.length / 3 - 1);
  const positions = new Float32Array(6 * count), colours = new Float32Array(8 * count);
  let s = 0;
  for (const c of curves) {
    if (!c.seeded) continue;
    const n = c.points.length / 3;
    let previous = n ? colourOf(c.magnitude[0]) : null;
    for (let i = 1; i < n; i++) {
      const here = colourOf(c.magnitude[i]);
      if (previous && here) {
        positions.set(c.points.subarray(3 * (i - 1), 3 * i + 3), 6 * s);
        colours.set(previous, 8 * s);
        colours.set(here, 8 * s + 4);
        s++;
      }
      previous = here;
    }
  }
  return { positions: positions.subarray(0, 6 * s), colours: colours.subarray(0, 8 * s), segments: s };
}
