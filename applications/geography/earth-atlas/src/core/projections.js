// Map projections of the sphere onto the plane, forward and inverse.
//
// Provisional stand-in for core/domains/geography (cartography; laboratory ADR 0014).
//
// Coordinates on the plane are in units of the sphere's radius. Each projection
// has a central meridian lambda0 (degrees); a longitude is first reduced to
// [-180, 180) about it. The fragment shader in src/render/surface.js transcribes
// the inverses; keep the two in step.
//
// Sources:
//   Equal Earth: Savric, B., Patterson, T., and Jenny, B. (2018), "The Equal Earth
//   map projection", Int. J. Geogr. Inf. Sci. 33(3), 454-465, eq. 1-3; the inverse
//   by Newton's method on the parametric latitude, as in the paper, sec. 3.
//   Equirectangular (plate carree): Snyder (1987), USGS Prof. Paper 1395, eq. 12-1.

import { DEG, wrapLongitude } from './sphere.js';

const A1 = 1.340264, A2 = -0.081106, A3 = 0.000893, A4 = 0.003796;
const M = Math.sqrt(3) / 2;

export const equalEarth = {
  id: 'equal-earth',
  name: 'Equal Earth',
  equalArea: true,
  // Plane coordinates (x, y) for a place, radians of longitude and latitude.
  raw(lambda, phi) {
    const t = Math.asin(M * Math.sin(phi)), t2 = t * t, t6 = t2 * t2 * t2;
    return [
      (lambda * Math.cos(t)) / (M * (A1 + 3 * A2 * t2 + t6 * (7 * A3 + 9 * A4 * t2))),
      t * (A1 + A2 * t2 + t6 * (A3 + A4 * t2)),
    ];
  },
  rawInverse(x, y) {
    let t = y, t2 = t * t, t6 = t2 * t2 * t2;
    for (let i = 0; i < 12; i++) {
      const fy = t * (A1 + A2 * t2 + t6 * (A3 + A4 * t2)) - y;
      const fpy = A1 + 3 * A2 * t2 + t6 * (7 * A3 + 9 * A4 * t2);
      const delta = fy / fpy;
      t -= delta; t2 = t * t; t6 = t2 * t2 * t2;
      if (Math.abs(delta) < 1e-12) break;
    }
    const lambda = (M * x * (A1 + 3 * A2 * t2 + t6 * (7 * A3 + 9 * A4 * t2))) / Math.cos(t);
    return [lambda, Math.asin(Math.max(-1, Math.min(1, Math.sin(t) / M)))];
  },
};

export const equirectangular = {
  id: 'equirectangular',
  name: 'Equirectangular',
  equalArea: false,
  raw: (lambda, phi) => [lambda, phi],
  rawInverse: (x, y) => [x, y],
};

export const PROJECTIONS = { [equalEarth.id]: equalEarth, [equirectangular.id]: equirectangular };

// Half-width and half-height of the projected sphere.
export function extent(p) {
  return [p.raw(Math.PI, 0)[0], p.raw(0, Math.PI / 2)[1]];
}

export function forward(p, lon, lat, lambda0 = 0) {
  return p.raw(wrapLongitude(lon - lambda0) * DEG, lat * DEG);
}

// The place at plane coordinates (x, y), or null outside the map.
export function inverse(p, x, y, lambda0 = 0) {
  const [, ymax] = extent(p);
  if (Math.abs(y) > ymax) return null;
  const [lambda, phi] = p.rawInverse(x, y);
  if (!(Math.abs(lambda) <= Math.PI + 1e-9)) return null;
  return [wrapLongitude(lambda / DEG + lambda0), phi / DEG];
}

// The map's outline: the projected meridians at lambda0 +/- 180, as a closed polyline.
export function outline(p, n = 90) {
  const pts = [];
  for (let i = 0; i <= n; i++) pts.push(p.raw(Math.PI, Math.PI / 2 - (Math.PI * i) / n));
  for (let i = 0; i <= n; i++) pts.push(p.raw(-Math.PI, -Math.PI / 2 + (Math.PI * i) / n));
  return pts;
}

// Ratio of an area on the map to the same area on the sphere, at a place
// (the Jacobian determinant over cos(phi)): constant for an equal-area projection.
export function areaScale(p, lon, lat, h = 1e-5) {
  const l = lon * DEG, f = lat * DEG;
  const [x1, y1] = p.raw(l + h, f), [x0, y0] = p.raw(l - h, f);
  const [x3, y3] = p.raw(l, f + h), [x2, y2] = p.raw(l, f - h);
  const det = ((x1 - x0) * (y3 - y2) - (x3 - x2) * (y1 - y0)) / (4 * h * h);
  return Math.abs(det) / Math.cos(f);
}
