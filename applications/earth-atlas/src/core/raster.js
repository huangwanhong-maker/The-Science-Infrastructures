// Polygons drawn onto an equirectangular raster by scanline filling.
//
// Provisional stand-in for core/domains/geography (GIS; laboratory ADR 0014).
//
// The raster covers longitude -180 to 180 and latitude 90 to -90 in width x
// height pixels; the centre of pixel (i, j) is at longitude -180 + (i + 0.5) 360 / width
// and latitude 90 - (j + 0.5) 180 / height. A pixel belongs to a polygon when its
// centre is inside by the even-odd rule, so holes need no special treatment and
// neighbouring polygons that share an edge never both claim a pixel. Rings are
// runs of (longitude, latitude) pairs in a Float32Array, given as [first pair,
// number of pairs], as in the geography product.
//
// Source: Foley, van Dam, Feiner and Hughes (1990), Computer Graphics: Principles
// and Practice, 2nd ed., sec. 3.6, "Filling polygons".

// Fills the pixels inside the rings (one feature, even-odd) with `value`.
export function fillRings(coords, runs, width, height, value, target) {
  const sx = width / 360, sy = height / 180;
  let ymin = Infinity, ymax = -Infinity;
  for (const [first, count] of runs) {
    for (let k = first; k < first + count; k++) {
      const y = (90 - coords[2 * k + 1]) * sy;
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
    }
  }
  const j0 = Math.max(0, Math.ceil(ymin - 0.5)), j1 = Math.min(height - 1, Math.floor(ymax - 0.5));
  if (j1 < j0) return 0;
  const rows = Array.from({ length: j1 - j0 + 1 }, () => []);
  for (const [first, count] of runs) {
    for (let k = 0; k < count; k++) {
      const a = first + k, b = first + ((k + 1) % count);
      const xa = (coords[2 * a] + 180) * sx, ya = (90 - coords[2 * a + 1]) * sy;
      const xb = (coords[2 * b] + 180) * sx, yb = (90 - coords[2 * b + 1]) * sy;
      if (ya === yb) continue;
      const [ylo, yhi, xlo, xhi] = ya < yb ? [ya, yb, xa, xb] : [yb, ya, xb, xa];
      // Rows whose centres y = j + 0.5 satisfy ylo <= y < yhi (half-open, so a
      // vertex shared by two edges is counted once).
      const ja = Math.max(j0, Math.ceil(ylo - 0.5)), jb = Math.min(j1, Math.ceil(yhi - 0.5) - 1);
      const slope = (xhi - xlo) / (yhi - ylo);
      for (let j = ja; j <= jb; j++) rows[j - j0].push(xlo + (j + 0.5 - ylo) * slope);
    }
  }
  let filled = 0;
  for (let j = j0; j <= j1; j++) {
    const xs = rows[j - j0].sort((p, q) => p - q);
    for (let m = 0; m + 1 < xs.length; m += 2) {
      const i0 = Math.max(0, Math.ceil(xs[m] - 0.5)), i1 = Math.min(width - 1, Math.ceil(xs[m + 1] - 0.5) - 1);
      for (let i = i0; i <= i1; i++) target[j * width + i] = value;
      filled += Math.max(0, i1 - i0 + 1);
    }
  }
  return filled;
}

// An index raster of features: pixel = 1 + index of the feature containing it, 0 elsewhere.
export function indexRaster(coords, features, width, height) {
  const target = new Uint16Array(width * height);
  features.forEach((f, i) => fillRings(coords, f.rings, width, height, i + 1, target));
  return target;
}

// A mask raster: 1 inside any of the features, 0 elsewhere.
export function maskRaster(coords, features, width, height, target = new Uint8Array(width * height), value = 1) {
  for (const f of features) fillRings(coords, f.rings, width, height, value, target);
  return target;
}

// Fraction of the sphere's area covered by the pixels with a given test.
export function areaFraction(raster, width, height, test) {
  let covered = 0, total = 0;
  for (let j = 0; j < height; j++) {
    const w = Math.cos(((90 - ((j + 0.5) * 180) / height) * Math.PI) / 180);
    for (let i = 0; i < width; i++) {
      total += w;
      if (test(raster[j * width + i])) covered += w;
    }
  }
  return covered / total;
}

// Area of a feature (rings, even-odd) as a fraction of the sphere, by the shoelace
// formula in the plane of the cylindrical equal-area projection (x = longitude
// in radians, y = sin latitude), where areas are proportional to those on the
// sphere. Edges are straight in that plane, a negligible difference from the
// raster's straight edges in longitude and latitude for edges of a fraction of
// a degree. Rings may run either way; the even-odd area is taken as the absolute
// sum of signed ring areas, which holds when holes run against their outer ring.
export function featureAreaFraction(coords, runs) {
  let sum = 0;
  for (const [first, count] of runs) {
    let ring = 0;
    for (let k = 0; k < count; k++) {
      const a = first + k, b = first + ((k + 1) % count);
      const xa = (coords[2 * a] * Math.PI) / 180, ya = Math.sin((coords[2 * a + 1] * Math.PI) / 180);
      const xb = (coords[2 * b] * Math.PI) / 180, yb = Math.sin((coords[2 * b + 1] * Math.PI) / 180);
      ring += xa * yb - xb * ya;
    }
    sum += ring / 2;
  }
  return Math.abs(sum) / (4 * Math.PI);
}

// Water in a grey shaded relief: pixels of the water's grey (within 1) whose
// eight neighbours have it too, grown back by one pixel, so that a lone pixel of
// that grey in a shadow is not taken for water. grey: Uint8Array (w x h).
export function waterMask(grey, w, h, level) {
  const at = new Uint8Array(w * h);
  for (let i = 0; i < w * h; i++) at[i] = Math.abs(grey[i] - level) <= 1 ? 1 : 0;
  const core = new Uint8Array(w * h);
  for (let y = 1; y < h - 1; y++) {
    for (let x = 1; x < w - 1; x++) {
      const i = y * w + x;
      if (at[i] && at[i - 1] && at[i + 1] && at[i - w] && at[i + w] && at[i - w - 1] && at[i - w + 1] && at[i + w - 1] && at[i + w + 1]) core[i] = 1;
    }
  }
  const out = new Uint8Array(w * h);
  for (let y = 0; y < h; y++) {
    for (let x = 0; x < w; x++) {
      if (!at[y * w + x]) continue;
      let wet = 0;
      for (let dy = -1; dy <= 1 && !wet; dy++) {
        for (let dx = -1; dx <= 1 && !wet; dx++) {
          const yy = y + dy, xx = x + dx;
          if (yy >= 0 && yy < h && xx >= 0 && xx < w && core[yy * w + xx]) wet = 1;
        }
      }
      out[y * w + x] = wet;
    }
  }
  return out;
}
