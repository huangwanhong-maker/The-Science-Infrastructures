// The base map: an equirectangular image of land and sea, tinted by elevation
// and depth and shaded by the relief, drawn under every data layer.
//
// Inputs are kept apart, as products: elevation from a field product (ETOPO
// 2022) and the outlines of land, lakes, glaciers and ice shelves from the
// geography product, rasterised at the image's resolution so that coasts are
// as sharp as the vector data allow.
//
// The relief is shaded by Lambert's law under a light from the north-west, 45
// degrees high, the cartographic convention. Slopes are taken from the elevation
// grid by central differences on the sphere and exaggerated, since a coarse grid
// smooths real slopes: the exaggeration is stated for a 10' grid and scaled with
// the grid's spacing (a 5' grid resolves steeper slopes and needs half of it).
// The darkening is limited, so that the steepest shaded slopes stay coloured.
// The shading is kept in the alpha channel, so that data layers drawn over the
// map can be shaded by the same relief.
//
// Sources: Horn (1981), "Hill shading and the reflectance map", Proc. IEEE 69(1);
// Imhof (1982), Cartographic Relief Presentation, ch. 8 and 13 (hypsometric tints).

import { maskRaster } from '../core/raster.js';

const hex = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];

// Tints: [elevation or depth in metres, colour], warm and quiet, so that data
// laid over them keep their own colours.
export const STYLES = {
  relief: {
    ocean: [[0, '#5d7ba3'], [50, '#4f6c97'], [200, '#3f5888'], [1000, '#30467a'], [3000, '#24366a'], [5000, '#1b2958'], [7000, '#141f45'], [11000, '#0e1532']],
    land: [[0, '#68745a'], [200, '#76805f'], [600, '#8e8c69'], [1200, '#a39072'], [2000, '#a8877d'], [3000, '#b0919a'], [4500, '#c9b4c3'], [6500, '#ece2ee']],
    ice: [[0, '#d9d4e6'], [2000, '#e6e2ef'], [4000, '#f3f0f7']],
    lake: '#56739d',
    shelf: '#c6c7dc',
    shade: { land: 0.85, ocean: 0.45, ice: 0.6, exaggeration: 14 },     // exaggeration for a 10' grid
  },
  plain: {
    ocean: [[0, '#241d33'], [11000, '#1b1527']],
    land: [[0, '#4a3c4f'], [6500, '#5a4b5f']],
    ice: [[0, '#5d5468'], [4000, '#685f72']],
    lake: '#241d33',
    shelf: '#4d4659',
    shade: { land: 0.55, ocean: 0.12, ice: 0.4, exaggeration: 14 },
  },
};

function ramp(stops, max, step) {
  const n = Math.ceil(max / step) + 1;
  const out = new Uint8Array(3 * n);
  const s = stops.map(([z, c]) => [z, hex(c)]);
  for (let i = 0; i < n; i++) {
    const z = i * step;
    let k = 0;
    while (k < s.length - 2 && z > s[k + 1][0]) k++;
    const [z0, c0] = s[k], [z1, c1] = s[Math.min(k + 1, s.length - 1)];
    const t = z1 === z0 ? 0 : Math.max(0, Math.min(1, (z - z0) / (z1 - z0)));
    for (let j = 0; j < 3; j++) out[3 * i + j] = Math.round(c0[j] + (c1[j] - c0[j]) * t);
  }
  return out;
}

// Lambert shading of the elevation grid, 1 on flat ground, about 0.3 to 1.6 on slopes.
export function hillshade(field, values, exaggeration) {
  const { nlon, nlat, dlon, dlat, lat0, periodic } = field.grid;
  const R = 6371008.8, rad = Math.PI / 180;
  const az = 315 * rad, alt = 45 * rad;
  const lx = Math.sin(az) * Math.cos(alt), ly = Math.cos(az) * Math.cos(alt), lz = Math.sin(alt);
  const out = new Float32Array(nlon * nlat);
  const dy = R * dlat * rad;
  for (let r = 0; r < nlat; r++) {
    const lat = lat0 - r * dlat;
    const dx = Math.max(R * Math.cos(lat * rad) * dlon * rad, dy * 0.05);
    const rn = Math.max(0, r - 1), rs = Math.min(nlat - 1, r + 1);
    for (let c = 0; c < nlon; c++) {
      const cw = c === 0 ? (periodic ? nlon - 1 : 0) : c - 1;
      const ce = c === nlon - 1 ? (periodic ? 0 : nlon - 1) : c + 1;
      const zx = (values[r * nlon + ce] - values[r * nlon + cw]) / (2 * dx);
      const zy = (values[rn * nlon + c] - values[rs * nlon + c]) / ((rs - rn) * dy);
      const nx = -exaggeration * zx, ny = -exaggeration * zy;
      const shade = (nx * lx + ny * ly + lz) / Math.sqrt(nx * nx + ny * ny + 1);
      out[r * nlon + c] = shade / lz;
    }
  }
  return out;
}

// The base image, width x height RGBA, row 0 at the north. Alpha holds the
// shading as 127.5 x shade (so 128 is flat ground).
export function baseImage({ elevation, geography, coords, width, height, style = 'relief' }) {
  const st = STYLES[style];
  const land = maskRaster(coords, geography.layers.land, width, height);
  maskRaster(coords, geography.layers.lakes, width, height, land, 2);
  maskRaster(coords, geography.layers.glaciers, width, height, land, 3);
  const shelves = maskRaster(coords, geography.layers.ice_shelves, width, height);
  const out = new Uint8ClampedArray(width * height * 4);

  const step = 10;
  const oceanRamp = ramp(st.ocean, 11000, step), landRamp = ramp(st.land, 9000, step), iceRamp = ramp(st.ice, 5000, step);
  const lake = hex(st.lake), shelf = hex(st.shelf);
  let z = null, shade = null;
  // Bilinear interpolation of the elevation grid at the image's pixel centres,
  // with the columns and rows of each pixel computed once (the grid has no gaps).
  const eg = elevation && elevation.grid;
  const c0 = new Int32Array(width), c1 = new Int32Array(width), fx = new Float32Array(width);
  const r0 = new Int32Array(height), r1 = new Int32Array(height), fy = new Float32Array(height);
  if (elevation) {
    z = elevation.slice(elevation.meta.variables[0].name, 0);
    shade = hillshade(elevation, z, st.shade.exaggeration * (eg.dlat / (10 / 60)));
    for (let i = 0; i < width; i++) {
      const x = ((((-180 + ((i + 0.5) * 360) / width - eg.lon0) / eg.dlon) % eg.nlon) + eg.nlon) % eg.nlon;
      c0[i] = Math.floor(x); c1[i] = (c0[i] + 1) % eg.nlon; fx[i] = x - c0[i];
    }
    for (let j = 0; j < height; j++) {
      const y = Math.max(0, Math.min(eg.nlat - 1, (eg.lat0 - (90 - ((j + 0.5) * 180) / height)) / eg.dlat));
      r0[j] = Math.min(eg.nlat - 2, Math.floor(y)); r1[j] = r0[j] + 1; fy[j] = y - r0[j];
    }
  }
  const at = (v, i, j) => {
    const a = r0[j] * eg.nlon, b = r1[j] * eg.nlon, x = fx[i], y = fy[j];
    return (v[a + c0[i]] * (1 - x) + v[a + c1[i]] * x) * (1 - y) + (v[b + c0[i]] * (1 - x) + v[b + c1[i]] * x) * y;
  };
  for (let j = 0; j < height; j++) {
    for (let i = 0; i < width; i++) {
      const k = j * width + i;
      const kind = land[k];                     // 0 sea, 1 land, 2 lake, 3 glacier
      const h = z ? at(z, i, j) : 0;
      const s = shade ? at(shade, i, j) : 1;
      let r, gr, b, strength;
      if (kind === 2) {
        [r, gr, b] = lake; strength = st.shade.ocean * 0.3;
      } else if (kind === 3 || (kind === 0 && shelves[k])) {
        const q = 3 * Math.min(iceRamp.length / 3 - 1, Math.max(0, Math.round(h / step)));
        [r, gr, b] = kind === 3 ? [iceRamp[q], iceRamp[q + 1], iceRamp[q + 2]] : shelf;
        strength = st.shade.ice;
      } else if (kind === 1) {
        const q = 3 * Math.min(landRamp.length / 3 - 1, Math.max(0, Math.round(h / step)));
        r = landRamp[q]; gr = landRamp[q + 1]; b = landRamp[q + 2];
        strength = st.shade.land;
      } else {
        const q = 3 * Math.min(oceanRamp.length / 3 - 1, Math.max(0, Math.round(-h / step)));
        r = oceanRamp[q]; gr = oceanRamp[q + 1]; b = oceanRamp[q + 2];
        strength = st.shade.ocean;
      }
      const f = Math.max(0.32, 1 + ((Number.isFinite(s) ? s : 1) - 1) * strength);
      out[4 * k] = r * f; out[4 * k + 1] = gr * f; out[4 * k + 2] = b * f;
      out[4 * k + 3] = Math.max(0, Math.min(255, Math.round(127.5 * (1 + ((Number.isFinite(s) ? s : 1) - 1) * (kind === 1 || kind === 3 ? 1 : 0.35)))));
    }
  }
  return { data: out, land };
}
