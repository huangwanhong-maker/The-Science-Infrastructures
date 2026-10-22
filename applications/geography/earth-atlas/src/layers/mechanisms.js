// The mechanisms layer: how the ground broke in the world's earthquakes, from
// the Global CMT catalogue, drawn as beach balls coloured by the kind of
// faulting: thrust (one side pushed up over the other, at plates that
// converge), normal (one side dropped, where the crust is stretched), and
// strike-slip (the sides slid past each other). Either the past year up to the
// time shown, or every earthquake of magnitude 6.5 and more since 1976.

import { BeachBallLayer } from '../render/beachballs.js';
import { mechanism, tensorNED } from '../core/moment.js';

export const KIND_COLOURS = { thrust: '#e0697f', normal: '#7ea4dd', 'strike-slip': '#86c29b', oblique: '#c08ac4' };
const hex = (h) => [parseInt(h.slice(1, 3), 16) / 255, parseInt(h.slice(3, 5), 16) / 255, parseInt(h.slice(5, 7), 16) / 255];
const YEAR = 365.25 * 86400000;

export class MechanismsLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.shown = false;
    this.mode = 'year';
    this.balls = new BeachBallLayer();
    this.balls.mesh.visible = false;
    stage.overlay.add(this.balls.mesh);
    this.drawn = [];
  }

  // Reads the catalogue once; true if there is one.
  async load() {
    if (this.product) return true;
    const entry = this.library.list('events').find((e) => e.kind === 'mechanism');
    if (!entry) return false;
    if (!this.loading) this.loading = this.library.eventSet(entry.name).then((p) => this.prepare(p));
    await this.loading;
    return true;
  }

  async show(on) {
    this.shown = on;
    if (on && !(await this.load())) this.shown = false;
    this.balls.mesh.visible = this.shown && Boolean(this.product);
    this.stale = true;
  }

  prepare(product) {
    this.product = product;
    this.f = Object.fromEntries(product.fields.map((k, i) => [k, i]));
    const f = this.f;
    // Each event's tensor (north, east, down), computed once; its class when first needed.
    this.events = product.events.map((r) => {
      const el = { mrr: r[f.mrr], mtt: r[f.mtt], mpp: r[f.mpp], mrt: r[f.mrt], mrp: r[f.mrp], mtp: r[f.mtp] };
      return { row: r, t: r[f.time_ms], mw: r[f.mw], lon: r[f.lon], lat: r[f.lat], m: tensorNED(el), elements: el, kind: null };
    });
    this.times = this.events.map((e) => e.t);
  }

  // Every event as the section figure wants it.
  all() {
    if (!this.events) return [];
    const f = this.f;
    return this.events.map((e) => {
      if (!e.kind) e.kind = mechanism(e.elements).kind;
      return { lon: e.lon, lat: e.lat, depth: e.row[f.depth_km], mw: e.mw, kind: e.kind, t: e.t, name: e.row[f.region] };
    });
  }

  setMode(mode) { this.mode = mode; this.stale = true; }

  // The events shown at a time: in the year before it, or the large ones before it.
  selected(time) {
    if (!this.events) return [];
    // The first event after the time shown, by bisection.
    let lo = 0, hi = this.times.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (this.times[m] <= time) lo = m + 1; else hi = m; }
    if (this.mode === 'year') {
      let a = lo;
      while (a > 0 && this.times[a - 1] > time - YEAR) a--;
      return this.events.slice(a, lo);
    }
    return this.events.slice(0, lo).filter((e) => e.mw >= 6.5);
  }

  radius(mw) { return Math.max(2.5, Math.min(20, 2.6 + (mw - 5) * 3.3)); }

  // Instances for the current view and time: called when the view has settled or the time moved.
  layout(time) {
    this.stale = false;
    this.time = time;
    if (!this.shown || !this.product) { this.balls.clear(); this.drawn = []; return; }
    const s = this.stage, proj = s.projector();
    const list = this.selected(time).slice().sort((a, b) => a.mw - b.mw);
    const n = list.length;
    const A = { aPos: new Float32Array(3 * n), aSize: new Float32Array(n), aM1: new Float32Array(4 * n), aM2: new Float32Array(2 * n), aNorth: new Float32Array(2 * n), aColor: new Float32Array(4 * n) };
    const p = [0, 0, 0];
    const drawn = [];
    let k = 0;
    for (const e of list) {
      const a = s.screen(e.lon, e.lat), b = s.screen(e.lon, Math.min(89.9, e.lat + 0.05));
      if (!a || !b || a[2] <= 0.3 || !proj.project(e.lon, e.lat, p, 1.0025)) continue;
      const dx = b[0] - a[0], dy = a[1] - b[1], len = Math.hypot(dx, dy) || 1;
      if (!e.kind) e.kind = mechanism(e.elements).kind;
      const m = e.m, c = hex(KIND_COLOURS[e.kind]);
      const age = this.mode === 'year' ? (time - e.t) / YEAR : 0;
      const size = this.radius(e.mw);
      A.aPos.set(p, 3 * k); A.aSize[k] = size;
      A.aM1.set([m[0][0], m[0][1], m[0][2], m[1][1]], 4 * k); A.aM2.set([m[1][2], m[2][2]], 2 * k);
      A.aNorth.set([dx / len, dy / len], 2 * k);
      A.aColor.set([c[0], c[1], c[2], Math.max(0.55, 1 - 0.45 * age) * a[2]], 4 * k);
      drawn.push({ e, x: a[0], y: a[1], r: size });
      k++;
    }
    this.balls.set(A, k);
    this.drawn = drawn;
  }

  frame(resolution) { if (this.shown) this.balls.frame(resolution); }

  // The topmost drawn event under a point of the screen.
  pick(x, y) {
    for (let i = this.drawn.length - 1; i >= 0; i--) {
      const d = this.drawn[i];
      if (Math.hypot(d.x - x, d.y - y) <= d.r + 2) return d.e;
    }
    return null;
  }

  find(name) {
    if (!this.events) return null;
    const i = this.f.name;
    return this.events.find((e) => e.row[i] === name) || null;
  }

  // Earthquakes near a place large enough to move it measurably: magnitude 6 and
  // more within 10^(0.5 Mw - 1.4) km (40 km at magnitude 6, 400 km at 8, 1,260 at 9).
  near(lon, lat, arcKm) {
    if (!this.events) return [];
    const out = [];
    for (const e of this.events) {
      if (e.mw < 6) continue;
      const km = arcKm(lon, lat, e.lon, e.lat);
      if (km <= 10 ** (0.5 * e.mw - 1.4)) out.push({ t: e.t, mw: e.mw, km, name: e.row[this.f.region] });
    }
    return out;
  }

  legend() {
    return {
      title: 'Earthquake mechanisms', count: `${this.drawn.length.toLocaleString('en-GB')} shown`,
      keys: Object.entries(KIND_COLOURS).map(([k, c]) => ({ label: k, colour: c, ball: true })),
      note: `${this.mode === 'year' ? 'the year before the time shown' : 'magnitude 6.5 and more, 1976 to the time shown'} · coloured quadrants: compression · Global CMT`,
    };
  }
}
