// A flight of the camera from one body to another: from a camera given in the world's axes to the
// view of a globe on another body, the bodies moving as it goes.
//
// The camera is carried relative to the target's centre, which moves with the instant: its distance
// from the centre interpolated in its logarithm, so that a flight from the Earth to the Moon spends
// its time evenly over the decades of distance, and raised in mid-flight to see both bodies (as a
// virtual globe's flight rises to see the ground it covers: van Wijk and Nuij 2003, in spirit); its
// direction from the centre turned along the great circle between the start and the end; its view
// turned from the start's to the end's, looking at the target in mid-flight.
//
//   from: { position, forward, up } in the world's axes (m), relative to the target's centre at the start
//   to:   { position, forward, up } likewise, at the end (the target's globe camera, turned to the world)
//   peak: the distance (m) to rise to, at least both ends' distances
//   s:    the flight's progress, 0 to 1 (eased by the caller or by ease())

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a) => { const n = norm(a) || 1; return [a[0] / n, a[1] / n, a[2] / n]; };
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

// The spherical interpolation of two unit vectors (any axis when they are opposite).
export function slerpVec(a, b, s) {
  const c = Math.min(1, Math.max(-1, dot(a, b)));
  const angle = Math.acos(c);
  if (angle < 1e-9) return a.slice();
  if (Math.PI - angle < 1e-6) {
    // Opposite: turn about any axis orthogonal to a.
    const axis = unit(Math.abs(a[0]) < 0.9 ? cross(a, [1, 0, 0]) : cross(a, [0, 1, 0]));
    const t = angle * s;
    const ax = cross(axis, a);
    return [0, 1, 2].map((k) => a[k] * Math.cos(t) + ax[k] * Math.sin(t));
  }
  const sa = Math.sin(angle);
  const wa = Math.sin((1 - s) * angle) / sa, wb = Math.sin(s * angle) / sa;
  return [0, 1, 2].map((k) => wa * a[k] + wb * b[k]);
}

export const ease = (t) => (t <= 0 ? 0 : t >= 1 ? 1 : t * t * (3 - 2 * t));

export function transferCamera(from, to, peak, s) {
  const r0 = norm(from.position), r1 = norm(to.position);
  const logR = (1 - s) * Math.log(Math.max(r0, 1)) + s * Math.log(Math.max(r1, 1));
  const rise = Math.sin(Math.PI * s);
  const r = Math.exp(logR) + rise * Math.max(0, peak - Math.exp(logR));
  const dir = slerpVec(unit(from.position), unit(to.position), s);
  const position = dir.map((x) => x * r);
  // The view: from the start's to the end's, looking at the target in mid-flight.
  const toward = dir.map((x) => -x);
  const base = slerpVec(unit(from.forward), unit(to.forward), s);
  const w = rise;
  const forward = unit([0, 1, 2].map((k) => (1 - w) * base[k] + w * toward[k]));
  let up = slerpVec(unit(from.up), unit(to.up), s);
  const along = dot(up, forward);
  up = unit([0, 1, 2].map((k) => up[k] - along * forward[k]));
  return { position, forward, up };
}
