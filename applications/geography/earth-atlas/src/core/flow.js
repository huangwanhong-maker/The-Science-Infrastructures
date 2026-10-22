// Particles carried by a vector field on the sphere, to make a flow visible.
//
// Provisional stand-in for core/domains (kinematics of a flow; laboratory ADR 0014).
//
// A particle at the unit vector p moves with dp/dt = w(p) / R, where w = u e + v n
// is the velocity (u east, v north, m/s) as a vector tangent to the sphere. Each
// step of length dt follows the midpoint (second-order Runge-Kutta) rule: the
// velocity is sampled half a step ahead at m, carried back to p by parallel
// transport, and p is turned along it by the angle |w| dt / R, which keeps it
// exactly on the sphere. The scheme is second-order accurate, so particles
// follow curved streamlines without drifting across them; it is exact for a
// steady flow along a great circle; and it holds at the poles, where east and
// north are undefined but w is not.
//
// A particle lives for a random number of steps, and dies early where the field
// has no value; it is then reborn at a random place chosen by the caller
// (usually within the visible part of the map). Each step records the segment
// the particle travelled, for drawing as a fading trail.
//
// The time step is a visual choice: a particle moves dt seconds of real flow per
// frame, so its speed on the screen is proportional to the flow's.
//
// Source: Hairer, Lubich and Wanner (2006), Geometric Numerical Integration,
// 2nd ed., sec. IV.5 (integration on manifolds by projection); Hairer, Norsett
// and Wanner (1993), Solving Ordinary Differential Equations I, sec. II.1.

import { DEG, EARTH_RADIUS_M } from './sphere.js';

// The velocity w = u e + v n at (lon, lat) as an Earth-fixed vector, into out.
function tangent(lon, lat, u, v, out) {
  const l = lon * DEG, f = lat * DEG;
  const sl = Math.sin(l), cl = Math.cos(l), sf = Math.sin(f), cf = Math.cos(f);
  out[0] = -sl * u - sf * cl * v;
  out[1] = cl * u - sf * sl * v;
  out[2] = cf * v;
  return out;
}

export class ParticleSystem {
  constructor(count, { minLife = 40, maxLife = 110, random = Math.random } = {}) {
    this.count = count;
    this.random = random;
    this.minLife = minLife;
    this.maxLife = maxLife;
    this.lon = new Float64Array(count);
    this.lat = new Float64Array(count);
    this.age = new Int32Array(count);
    this.life = new Int32Array(count);
    // Per step: the segment travelled (lon0, lat0, lon1, lat1), its speed (m/s),
    // and whether the particle was alive along it.
    this.segments = new Float32Array(4 * count);
    this.speed = new Float32Array(count);
    this.alive = new Uint8Array(count);
    this._w = [0, 0, 0];
  }

  reset(spawn) {
    for (let i = 0; i < this.count; i++) {
      this.reborn(i, spawn);
      this.age[i] = Math.floor(this.random() * this.life[i]);      // stagger the first deaths
    }
  }

  reborn(i, spawn) {
    const [lon, lat] = spawn(this.random);
    this.lon[i] = lon; this.lat[i] = lat;
    this.age[i] = 0;
    this.life[i] = this.minLife + Math.floor(this.random() * (this.maxLife - this.minLife));
  }

  // sample(lon, lat) -> [u, v] in m/s, or null where the field has no value.
  step(sample, dt, spawn) {
    const w = this._w, s = this.segments, k = dt / EARTH_RADIUS_M;
    for (let i = 0; i < this.count; i++) {
      const lon = this.lon[i], lat = this.lat[i];
      this.age[i]++;
      const v1 = this.age[i] > this.life[i] ? null : sample(lon, lat);
      if (!v1) { this.dead(i, spawn); continue; }
      const l = lon * DEG, f = lat * DEG, cf = Math.cos(f);
      const px = cf * Math.cos(l), py = cf * Math.sin(l), pz = Math.sin(f);
      // Half a step ahead (to first order), where the velocity is sampled again.
      tangent(lon, lat, v1[0], v1[1], w);
      let mx = px + w[0] * k / 2, my = py + w[1] * k / 2, mz = pz + w[2] * k / 2;
      const mr = Math.hypot(mx, my, mz);
      mx /= mr; my /= mr; mz /= mr;
      const mlon = Math.atan2(my, mx) / DEG, mlat = Math.asin(Math.max(-1, Math.min(1, mz))) / DEG;
      const v2 = sample(mlon, mlat);
      if (!v2) { this.dead(i, spawn); continue; }
      // The midpoint velocity, carried back to p by parallel transport along the
      // great circle from m to p: t = w - (w.p) / (1 + m.p) (m + p).
      tangent(mlon, mlat, v2[0], v2[1], w);
      const f2 = (w[0] * px + w[1] * py + w[2] * pz) / (1 + mx * px + my * py + mz * pz);
      const tx = w[0] - f2 * (mx + px), ty = w[1] - f2 * (my + py), tz = w[2] - f2 * (mz + pz);
      const speed = Math.hypot(tx, ty, tz);
      let qx = px, qy = py, qz = pz;
      if (speed > 0) {
        const angle = speed * k, c = Math.cos(angle), sn = Math.sin(angle) / speed;
        qx = c * px + sn * tx; qy = c * py + sn * ty; qz = c * pz + sn * tz;
      }
      const qlon = Math.atan2(qy, qx) / DEG, qlat = Math.asin(Math.max(-1, Math.min(1, qz))) / DEG;
      s[4 * i] = lon; s[4 * i + 1] = lat; s[4 * i + 2] = qlon; s[4 * i + 3] = qlat;
      this.speed[i] = Math.hypot(v2[0], v2[1]);
      this.alive[i] = 1;
      this.lon[i] = qlon; this.lat[i] = qlat;
    }
  }

  dead(i, spawn) {
    this.alive[i] = 0;
    this.speed[i] = 0;
    this.reborn(i, spawn);
  }
}

// A random place, uniform in area, within a band of latitude and an interval of
// longitude (which may wrap round the antimeridian).
export function uniformSpawn(latMin = -90, latMax = 90, lonMin = -180, lonMax = 180) {
  const s0 = Math.sin((latMin * Math.PI) / 180), s1 = Math.sin((latMax * Math.PI) / 180);
  const width = ((lonMax - lonMin) % 360 + 360) % 360 || 360;
  return (random) => {
    const lat = (Math.asin(s0 + (s1 - s0) * random()) * 180) / Math.PI;
    let lon = lonMin + width * random();
    if (lon >= 180) lon -= 360;
    return [lon, lat];
  };
}
