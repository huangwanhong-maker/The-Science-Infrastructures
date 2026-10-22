// Osculating two-body orbits from a state vector.
//
// Provisional stand-in for core/domains/astronomy and physis/mechanics
// (laboratory ADR 0011).
//
// The osculating orbit is the conic the body would follow if, from this instant,
// only the central body attracted it. It is drawn from the eccentricity vector
// e = (v x h)/mu - r/|r| and the angular momentum h = r x v, without passing
// through angular elements, so it has no singularity at small e or i.
// Source: Vallado, Fundamentals of Astrodynamics and Applications (2013), sec. 2.5.

import { cross, dot, norm, scale, sub, unit } from './frames.js';

export function osculating(r, v, mu) {
  const h = cross(r, v);
  const hNorm = norm(h);
  const rNorm = norm(r);
  const eVec = sub(scale(cross(v, h), 1 / mu), scale(r, 1 / rNorm));
  const e = norm(eVec);
  const p = (hNorm * hNorm) / mu;
  const energy = dot(v, v) / 2 - mu / rNorm;
  const a = -mu / (2 * energy);
  const P = e > 1e-10 ? scale(eVec, 1 / e) : unit(r);
  const W = scale(h, 1 / hNorm);
  const Q = cross(W, P);
  const period = a > 0 ? 2 * Math.PI * Math.sqrt(a ** 3 / mu) : Infinity;
  return { a, e, p, P, Q, W, period };
}

// Points of the osculating ellipse relative to the focus, in the frame of r and v.
// Returns a Float64Array of 3 * count coordinates, or null for an open orbit.
export function orbitPoints(r, v, mu, count = 256) {
  const o = osculating(r, v, mu);
  if (!(o.e < 1)) return null;
  const out = new Float64Array(3 * count);
  for (let k = 0; k < count; k++) {
    const nu = (2 * Math.PI * k) / count;
    const radius = o.p / (1 + o.e * Math.cos(nu));
    const c = radius * Math.cos(nu), s = radius * Math.sin(nu);
    out[3 * k] = c * o.P[0] + s * o.Q[0];
    out[3 * k + 1] = c * o.P[1] + s * o.Q[1];
    out[3 * k + 2] = c * o.P[2] + s * o.Q[2];
  }
  return out;
}
