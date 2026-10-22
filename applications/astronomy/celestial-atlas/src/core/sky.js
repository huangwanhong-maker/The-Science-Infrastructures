// Apparent places of the Sun, the Moon and the planets at one instant: as seen
// from the Earth's centre and from an observer on its surface.
//
// Provisional stand-in for core/domains/astronomy (laboratory ADR 0011).
//
// Applied: light-time; precession to the mean equator of date; topocentric
// parallax (about 1° for the Moon); refraction for apparent altitudes.
// Not applied: annual aberration (≤ 20.5"), nutation (≤ 17"), gravitational
// light deflection, polar motion. The resulting directions are good to well
// under an arcminute, which is ample for charts and views.

import { BODIES, BY_KEY, orbitalMu } from './bodies.js';
import { gmst } from './time.js';
import {
  LIGHT_KM_S, apply, multiply, rotY, rotZ, precessionMatrix, geodeticToTerrestrial,
  terrestrialToMeanOfDate, add, sub, scale, norm, dot, cross, toSpherical, refraction,
} from './frames.js';
import { orbitPoints } from './kepler.js';

export function angleBetween(a, b) {
  return Math.atan2(norm(cross(a, b)), dot(a, b));
}

// Horizontal coordinates from a vector already rotated into the horizon frame.
export function altAzOf(h) {
  const r = norm(h);
  let az = Math.atan2(h[1], -h[0]);
  if (az < 0) az += 2 * Math.PI;
  return { alt: Math.asin(h[2] / r), az };
}

export function snapshot(spk, scales, observer) {
  const et = scales.et;
  const P = precessionMatrix(scales.jdTT);
  const g = gmst(scales.jdUT1, scales.jdTT);
  const lst = g + observer.lon;
  const dateToHorizon = multiply(rotY(Math.PI / 2 - observer.lat), rotZ(lst));
  const icrsToHorizon = multiply(dateToHorizon, P);
  const observerDate = terrestrialToMeanOfDate(
    geodeticToTerrestrial(observer.lat, observer.lon, observer.height || 0), g);

  const ssb = {};
  for (const b of BODIES) ssb[b.key] = spk.stateSSB(b.id, et);
  const earth = ssb.earth.p;
  const sun = ssb.sun;
  const embSsb = spk.stateSSB(3, et);
  const emb = { p: sub(embSsb.p, sun.p), v: sub(embSsb.v, sun.v) };

  const bodies = {};
  for (const b of BODIES) {
    const s = ssb[b.key];
    const entry = { body: b, ssb: s, helio: { p: sub(s.p, sun.p), v: sub(s.v, sun.v) } };
    if (b.key !== 'earth') {
      let tau = 0, geo = sub(s.p, earth);
      for (let i = 0; i < 2; i++) {
        tau = norm(geo) / LIGHT_KM_S;
        geo = sub(spk.stateSSB(b.id, et - tau).p, earth);
      }
      const radec = toSpherical(geo);
      const topo = sub(apply(P, geo), observerDate);
      const hz = altAzOf(apply(dateToHorizon, topo));
      Object.assign(entry, {
        geo,
        lightTime: tau,
        distKm: norm(geo),
        ra: radec.lon,
        dec: radec.lat,
        alt: hz.alt,
        az: hz.az,
        altApparent: hz.alt + refraction(hz.alt),
      });
    }
    bodies[b.key] = entry;
  }

  const sunGeo = bodies.sun.geo;
  for (const b of BODIES) {
    const e = bodies[b.key];
    if (!e.geo || b.key === 'sun') continue;
    e.elongation = angleBetween(e.geo, sunGeo);
    const toSun = sub(sun.p, add(earth, e.geo));
    e.phaseAngle = angleBetween(toSun, scale(e.geo, -1));
    e.illuminated = (1 + Math.cos(e.phaseAngle)) / 2;
  }

  return { scales, et, precession: P, gmst: g, lst, icrsToHorizon, dateToHorizon, observer, bodies, emb };
}

// Osculating orbit of a body, as points in km (ICRF) relative to its primary.
// The Earth's orbit is that of the Earth-Moon barycentre, which is free of the
// monthly wobble the Moon imposes on the Earth itself.
export function osculatingOrbit(snap, key, count = 256) {
  const body = BY_KEY[key];
  if (key === 'sun') return null;
  if (body.parent) {
    const e = snap.bodies[key].ssb;
    const parent = snap.bodies[body.parent].ssb;
    return orbitPoints(sub(e.p, parent.p), sub(e.v, parent.v), orbitalMu(body), count);
  }
  const state = key === 'earth' ? snap.emb : snap.bodies[key].helio;
  return orbitPoints(state.p, state.v, orbitalMu(body), count);
}
