// Rigid motions sampled by an engine, interpolated for drawing at the frame rate.
//
// A page cannot ask the engine for the orientation of the Earth sixty times a second; it asks for
// a window of samples (a body's rotation to the world's axes and its position, every few minutes)
// and interpolates between them here. Rotations are unit quaternions (w, x, y, z) of the active
// rotation from the body's frame to the world's, interpolated along the great arc between samples
// (slerp; Shoemake 1985, "Animating rotation with quaternion curves", SIGGRAPH); positions by the
// cubic Hermite polynomial through the samples' positions and velocities, which is exact for any
// motion of constant acceleration over a step.
//
// A body that spins fast (the Earth: a turn a day) is interpolated with its spin factored out: the
// samples are R(tᵢ)·R_z(−ω tᵢ), slowly varying (precession, nutation, polar motion and the
// difference between ω and the true rate), and R(t) = S(t)·R_z(ω t). The interpolation is then
// good between samples hours apart, where plain slerp between rotations half a turn apart would
// take the shorter way round. With samples ten minutes apart the Earth's rotation is interpolated
// to well under a milliarcsecond.

// ---- quaternions --------------------------------------------------------------------------------

export const quatMul = (a, b) => [
  a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
  a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
  a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
  a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
];
export const quatConj = (q) => [q[0], -q[1], -q[2], -q[3]];
export const quatAboutZ = (angle) => [Math.cos(angle / 2), 0, 0, Math.sin(angle / 2)];

export function quatNormalise(q) {
  const n = Math.hypot(q[0], q[1], q[2], q[3]);
  return [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
}

// The great arc from a to b, at fraction t; the shorter way (q and −q are one rotation).
export function quatSlerp(a, b, t) {
  let d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
  let bb = b;
  if (d < 0) { d = -d; bb = [-b[0], -b[1], -b[2], -b[3]]; }
  if (d > 0.9999995) return quatNormalise(a.map((x, k) => x + (bb[k] - x) * t));
  const theta = Math.acos(d), s = Math.sin(theta);
  const wa = Math.sin((1 - t) * theta) / s, wb = Math.sin(t * theta) / s;
  return [wa * a[0] + wb * bb[0], wa * a[1] + wb * bb[1], wa * a[2] + wb * bb[2], wa * a[3] + wb * bb[3]];
}

// The 3 x 3 matrix (row-major arrays) of the rotation a unit quaternion performs.
export function quatToMatrix(q) {
  const [w, x, y, z] = q;
  return [
    [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
    [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
    [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
  ];
}

export const apply = (R, v) => [
  R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
  R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
  R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
];
export const applyTransposed = (R, v) => [
  R[0][0] * v[0] + R[1][0] * v[1] + R[2][0] * v[2],
  R[0][1] * v[0] + R[1][1] * v[1] + R[2][1] * v[2],
  R[0][2] * v[0] + R[1][2] * v[1] + R[2][2] * v[2],
];

// The angle (rad) of the rotation between two unit quaternions: 2 atan2(|v|, |w|) of a⁻¹b, exact to
// rounding for small angles too (an arccosine of the dot product loses half the digits there).
export function quatAngle(a, b) {
  const q = quatMul(quatConj(a), b);
  return 2 * Math.atan2(Math.hypot(q[1], q[2], q[3]), Math.abs(q[0]));
}

// ---- sampled motions ------------------------------------------------------------------------------

// The index i of the step [t_i, t_{i+1}] holding t (clamped to the window), and the fraction.
function locate(times, t) {
  const n = times.length;
  if (n < 2 || t <= times[0]) return [0, 0];
  if (t >= times[n - 1]) return [n - 2, 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (times[mid] <= t) lo = mid; else hi = mid;
  }
  return [lo, (t - times[lo]) / (times[lo + 1] - times[lo])];
}

// A rotation sampled at times (s), with an optional spin ω (rad/s) about the body's z axis factored out.
export class SampledRotation {
  constructor(times, quaternions, { spin = 0 } = {}) {
    this.times = times;
    this.spin = spin;
    this.slow = quaternions.map((q, i) => quatMul(q, quatAboutZ(-spin * times[i])));
    // Keep consecutive samples in one hemisphere, so that each step is the shorter arc.
    for (let i = 1; i < this.slow.length; i++) {
      const a = this.slow[i - 1], b = this.slow[i];
      if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] < 0) this.slow[i] = b.map((x) => -x);
    }
  }
  covers(t) { return this.times.length > 1 && t >= this.times[0] && t <= this.times[this.times.length - 1]; }
  // The quaternion at time t.
  at(t) {
    const [i, f] = locate(this.times, t);
    const s = this.slow.length > 1 ? quatSlerp(this.slow[i], this.slow[i + 1], f) : this.slow[0];
    return quatMul(s, quatAboutZ(this.spin * t));
  }
  matrix(t) { return quatToMatrix(this.at(t)); }
}

// A path sampled at times (s) with positions (m) and, if known, velocities (m/s).
export class SampledPath {
  constructor(times, positions, velocities = null) {
    Object.assign(this, { times, positions, velocities });
  }
  covers(t) { return this.times.length > 1 && t >= this.times[0] && t <= this.times[this.times.length - 1]; }
  at(t) {
    const [i, f] = locate(this.times, t);
    const p0 = this.positions[i], p1 = this.positions[Math.min(i + 1, this.positions.length - 1)];
    if (!this.velocities || this.positions.length < 2) return p0.map((x, k) => x + (p1[k] - x) * f);
    const h = this.times[i + 1] - this.times[i];
    const v0 = this.velocities[i], v1 = this.velocities[i + 1];
    // Cubic Hermite basis.
    const f2 = f * f, f3 = f2 * f;
    const h00 = 2 * f3 - 3 * f2 + 1, h10 = f3 - 2 * f2 + f, h01 = -2 * f3 + 3 * f2, h11 = f3 - f2;
    return [0, 1, 2].map((k) => h00 * p0[k] + h10 * h * v0[k] + h01 * p1[k] + h11 * h * v1[k]);
  }
}

// The Earth's nominal rate of rotation, rad per second of UT1 (the derivative of the Earth rotation
// angle, IERS Conventions 2010, eq. 5.15): the spin factored out of its sampled orientation.
export const EARTH_SPIN = 2 * Math.PI * 1.00273781191135448 / 86400;
