// The camera of a globe, in the manner of a virtual globe: a point looked at, a range, a heading
// and a tilt, in the body-fixed frame of the body in focus; and the gestures that change them.
//
// The view is held as a target (geodetic longitude, latitude and height), the range from it to the
// camera, the heading of the view (from north through east) and its tilt from the vertical. In the
// target's local frame (east, north, up) the camera stands at
//   range · (−sin t sin h, −sin t cos h, cos t),
// looking at the target, its up vector cos t (sin h, cos h, 0) + sin t (0, 0, 1). Holding the view in
// the body-fixed frame keeps the ground still under the camera while the body turns in space.
//
// Gestures: dragging grabs the ground (the point under the pointer stays under it: the target
// is turned about the body's centre by the rotation that carries the point now under the pointer
// to the point first grabbed); the wheel zooms toward the pointer; dragging with the right button
// (or with Ctrl) turns the heading and the tilt; a double click flies toward the point; the arrow
// keys pan, + and − zoom, N turns north up and U looks straight down. A flight follows the great
// circle between targets, rising and falling in range (van Wijk and Nuij 2003, "Smooth and
// efficient zooming and panning", in spirit).

import { cartesian, geodetic, rayHit } from './figure.js';

const DEG = Math.PI / 180;
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

// The local frame at a geodetic point: rows east, north, up (body-fixed components).
export function localFrame(lon, lat) {
  const sl = Math.sin(lon), cl = Math.cos(lon), sp = Math.sin(lat), cp = Math.cos(lat);
  return [[-sl, cl, 0], [-sp * cl, -sp * sl, cp], [cp * cl, cp * sl, sp]];
}

// The camera's body-fixed position, view direction and up vector for a view.
export function cameraOf(fig, view) {
  const { lon, lat, height, range, heading, tilt } = view;
  const target = cartesian(fig, lon, lat, height);
  const [E, N, U] = localFrame(lon, lat);
  const st = Math.sin(tilt), ct = Math.cos(tilt), sh = Math.sin(heading), ch = Math.cos(heading);
  const off = [-st * sh, -st * ch, ct]; // in east, north, up
  const position = [0, 1, 2].map((k) => target[k] + range * (off[0] * E[k] + off[1] * N[k] + off[2] * U[k]));
  const forward = [0, 1, 2].map((k) => -(off[0] * E[k] + off[1] * N[k] + off[2] * U[k]));
  const up = [0, 1, 2].map((k) => ct * (sh * E[k] + ch * N[k]) + st * U[k]);
  return { position, forward, up, target };
}

// Rotate vector v about the unit axis by angle (Rodrigues).
function rotate(v, axis, angle) {
  const c = Math.cos(angle), s = Math.sin(angle), d = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2];
  const x = [axis[1] * v[2] - axis[2] * v[1], axis[2] * v[0] - axis[0] * v[2], axis[0] * v[1] - axis[1] * v[0]];
  return [0, 1, 2].map((k) => v[k] * c + x[k] * s + axis[k] * d * (1 - c));
}
function unit(v) {
  const n = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / n, v[1] / n, v[2] / n];
}

export class GlobeControls {
  // element: the canvas; rayAt(clientX, clientY) → { origin, direction } in the body-fixed frame of
  // the body in focus; heightAt(lon, lat) → the terrain's height there (m) or 0.
  constructor(element, { fig, rayAt, heightAt = () => 0, onChange = () => {} }) {
    Object.assign(this, { element, fig, rayAt, heightAt, onChange });
    this.view = { lon: 0, lat: 0, height: 0, range: 2.6 * fig.a, heading: 0, tilt: 0 };
    this.flight = null;
    this.minRange = 2;
    this.maxRange = 60 * fig.a;
    this.attach();
  }

  // A new body in focus; maxRange how far the camera may stand from it (by default 60 radii; a view
  // of a planetary system or of the Galaxy asks for more).
  setFigure(fig, { maxRange = 60 * fig.a } = {}) {
    this.fig = fig;
    this.maxRange = maxRange;
  }

  // The point of the ground under a pointer: the ray's hit on the figure raised to the height of the
  // terrain there (twice refined), or null when the ray misses the body.
  pick(clientX, clientY) {
    const { origin, direction } = this.rayAt(clientX, clientY);
    let hit = rayHit(this.fig, origin, direction, 0);
    if (!hit) return null;
    for (let i = 0; i < 2; i++) {
      const g = geodetic(this.fig, ...hit);
      const h = this.heightAt(g.lon, g.lat);
      const again = rayHit(this.fig, origin, direction, h);
      if (!again) break;
      hit = again;
    }
    return hit;
  }

  camera() {
    return cameraOf(this.fig, this.view);
  }

  setView(view) {
    Object.assign(this.view, view);
    this.limit();
    this.onChange(this.view);
  }

  limit() {
    const v = this.view;
    v.lat = clamp(v.lat, -89.999 * DEG, 89.999 * DEG);
    v.lon = Math.atan2(Math.sin(v.lon), Math.cos(v.lon));
    v.tilt = clamp(v.tilt, 0, 88 * DEG);
    v.range = clamp(v.range, this.minRange, this.maxRange);
    // Keep the camera above the ground beneath it.
    const cam = cameraOf(this.fig, v).position;
    const g = geodetic(this.fig, ...cam);
    const ground = this.heightAt(g.lon, g.lat);
    if (g.h < ground + 2) {
      // Lower the tilt, then lengthen the range, until it clears.
      for (let i = 0; i < 40 && geodetic(this.fig, ...cameraOf(this.fig, v).position).h < ground + 2; i++) {
        if (v.tilt > 0.02) v.tilt *= 0.9;
        else v.range *= 1.05;
      }
    }
  }

  // Fly to a view over `seconds`.
  flyTo(target, seconds = 2.5) {
    const from = { ...this.view };
    const to = { ...from, ...target };
    const a = cartesian(this.fig, from.lon, from.lat, 0), b = cartesian(this.fig, to.lon, to.lat, 0);
    const ua = unit(a), ub = unit(b);
    // The arc between them as atan2(|a × b|, a · b), exact for small arcs; a flight to where one is
    // already looking has no axis (acos of a rounded dot product gave an arc of 1.5e-8 rad about the
    // normalised zero cross product: a NaN camera).
    const c = [ua[1] * ub[2] - ua[2] * ub[1], ua[2] * ub[0] - ua[0] * ub[2], ua[0] * ub[1] - ua[1] * ub[0]];
    const sinAngle = Math.hypot(c[0], c[1], c[2]), cosAngle = ua[0] * ub[0] + ua[1] * ub[1] + ua[2] * ub[2];
    const angle = Math.atan2(sinAngle, cosAngle);
    const axis = sinAngle > 1e-15 ? unit(c)
      : cosAngle < 0 ? unit(Math.abs(ua[2]) < 0.9 ? [-ua[1], ua[0], 0] : [0, -ua[2], ua[1]])   // antipodes: any axis across
        : [0, 0, 1];
    // Rise by the ground covered, so that the flight sees where it goes.
    const peak = Math.max(from.range, to.range, angle * this.fig.a * 1.1);
    let dh = to.heading - from.heading;
    dh = Math.atan2(Math.sin(dh), Math.cos(dh));
    this.flight = { start: performance.now(), seconds, from, to, ua, axis, angle, peak, dh };
  }

  // Advance a flight; called each frame. Returns whether one is under way.
  tick(now = performance.now()) {
    const f = this.flight;
    if (!f) return false;
    const t = clamp((now - f.start) / (f.seconds * 1000), 0, 1);
    const s = t * t * (3 - 2 * t); // ease in and out
    const p = rotate(f.ua, f.axis, f.angle * s);
    const g = geodetic(this.fig, ...p.map((x) => x * this.fig.a));
    const rise = Math.sin(Math.PI * s);
    const logR = (1 - s) * Math.log(f.from.range) + s * Math.log(f.to.range);
    const range = Math.exp(logR) + rise * Math.max(0, f.peak - Math.exp(logR)) * 0.6;
    this.view.lon = g.lon;
    this.view.lat = g.lat;
    this.view.height = (1 - s) * f.from.height + s * f.to.height;
    this.view.range = range;
    this.view.heading = f.from.heading + f.dh * s;
    this.view.tilt = (1 - s) * f.from.tilt + s * f.to.tilt;
    this.limit();
    this.onChange(this.view);
    if (t >= 1) this.flight = null;
    return true;
  }

  // ---- gestures ---------------------------------------------------------------------------------
  attach() {
    const el = this.element;
    let drag = null;
    el.addEventListener('contextmenu', (e) => e.preventDefault());
    el.addEventListener('pointerdown', (e) => {
      el.setPointerCapture(e.pointerId);
      this.flight = null;
      const rotateMode = e.button === 2 || e.ctrlKey;
      const grabbed = rotateMode ? null : this.pick(e.clientX, e.clientY);
      drag = { x: e.clientX, y: e.clientY, rotateMode, grabbed: grabbed ? unit(grabbed) : null, view: { ...this.view } };
    });
    el.addEventListener('pointermove', (e) => {
      if (!drag) return;
      if (drag.rotateMode) {
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
        this.view.heading = drag.view.heading - dx * 0.006;
        this.view.tilt = drag.view.tilt - dy * 0.005;
        this.limit();
        this.onChange(this.view);
        return;
      }
      if (!drag.grabbed) return;
      // Where the grabbed point's sphere meets the pointer's ray now.
      const { origin, direction } = this.rayAt(e.clientX, e.clientY);
      const r = Math.hypot(...cartesian(this.fig, this.view.lon, this.view.lat, this.view.height));
      const hit = rayHitSphere(origin, direction, r);
      if (!hit) return;
      const now = unit(hit);
      // Turn the target by the rotation taking `now` to `grabbed`.
      const axisRaw = [now[1] * drag.grabbed[2] - now[2] * drag.grabbed[1], now[2] * drag.grabbed[0] - now[0] * drag.grabbed[2], now[0] * drag.grabbed[1] - now[1] * drag.grabbed[0]];
      const sinA = Math.hypot(...axisRaw), cosA = now[0] * drag.grabbed[0] + now[1] * drag.grabbed[1] + now[2] * drag.grabbed[2];
      if (sinA < 1e-12) return;
      const axis = axisRaw.map((x) => x / sinA), angle = Math.atan2(sinA, cosA);
      const t = cartesian(this.fig, this.view.lon, this.view.lat, this.view.height);
      const moved = rotate(t, axis, angle);
      const g = geodetic(this.fig, ...moved);
      // Turning the target turns the local north: carry the heading with it.
      const northBefore = localFrame(this.view.lon, this.view.lat)[1];
      const northMoved = rotate(northBefore, axis, angle);
      const [E2, N2] = localFrame(g.lon, g.lat);
      const twist = Math.atan2(northMoved[0] * E2[0] + northMoved[1] * E2[1] + northMoved[2] * E2[2], northMoved[0] * N2[0] + northMoved[1] * N2[1] + northMoved[2] * N2[2]);
      this.view.lon = g.lon;
      this.view.lat = g.lat;
      this.view.heading -= twist;
      this.limit();
      this.onChange(this.view);
    });
    const end = (e) => {
      drag = null;
      try { el.releasePointerCapture(e.pointerId); } catch { /* already released */ }
    };
    el.addEventListener('pointerup', end);
    el.addEventListener('pointercancel', end);
    el.addEventListener('wheel', (e) => {
      e.preventDefault();
      this.flight = null;
      const factor = Math.exp(Math.sign(e.deltaY) * Math.min(Math.abs(e.deltaY), 200) * 0.0018);
      const towards = this.pick(e.clientX, e.clientY);
      const before = this.view.range;
      this.view.range = clamp(before * factor, this.minRange, this.maxRange);
      if (towards && factor < 1) {
        // Move the target toward the point under the pointer by the share of the range removed.
        const g = geodetic(this.fig, ...towards);
        const share = 1 - this.view.range / before;
        const a = unit(cartesian(this.fig, this.view.lon, this.view.lat, 0)), b = unit(cartesian(this.fig, g.lon, g.lat, 0));
        const mid = unit(a.map((x, k) => x + (b[k] - x) * share));
        const m = geodetic(this.fig, ...mid.map((x) => x * this.fig.a));
        this.view.lon = m.lon;
        this.view.lat = m.lat;
        this.view.height = this.heightAt(m.lon, m.lat);
      }
      this.limit();
      this.onChange(this.view);
    }, { passive: false });
    el.addEventListener('dblclick', (e) => {
      const p = this.pick(e.clientX, e.clientY);
      if (!p) return;
      const g = geodetic(this.fig, ...p);
      this.flyTo({ lon: g.lon, lat: g.lat, height: g.h, range: this.view.range * 0.35 }, 1.2);
    });
    el.tabIndex = 0;
    el.addEventListener('keydown', (e) => {
      const step = this.view.range / this.fig.a * 0.08;
      const v = this.view;
      const k = e.key;
      // Forward is along the heading, right is a quarter turn clockwise from it.
      const move = (forward, right) => {
        const c = Math.max(Math.cos(v.lat), 0.05), sh = Math.sin(v.heading), ch = Math.cos(v.heading);
        v.lat += step * (forward * ch - right * sh);
        v.lon += (step * (forward * sh + right * ch)) / c;
      };
      if (k === 'ArrowLeft') move(0, -1);
      else if (k === 'ArrowRight') move(0, 1);
      else if (k === 'ArrowUp') move(1, 0);
      else if (k === 'ArrowDown') move(-1, 0);
      else if (k === '+' || k === '=') v.range *= 0.8;
      else if (k === '-' || k === '_') v.range *= 1.25;
      else if (k === 'n' || k === 'N') v.heading = 0;
      else if (k === 'u' || k === 'U') v.tilt = 0;
      else return;
      e.preventDefault();
      this.limit();
      this.onChange(v);
    });
  }
}

// The farther intersection of a ray with a sphere about the origin when the origin is outside it,
// the nearer when it is in front: the point the grab follows.
function rayHitSphere(o, d, r) {
  const a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2], b = 2 * (o[0] * d[0] + o[1] * d[1] + o[2] * d[2]), c = o[0] * o[0] + o[1] * o[1] + o[2] * o[2] - r * r;
  let disc = b * b - 4 * a * c;
  if (disc < 0) {
    // Off the limb: the point of the ray closest to the sphere, so that the drag continues.
    const t = -b / (2 * a);
    const p = [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]];
    const n = Math.hypot(...p);
    return p.map((x) => (x / n) * r);
  }
  disc = Math.sqrt(disc);
  const t = (-b - disc) / (2 * a);
  return t > 0 ? [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]] : null;
}
