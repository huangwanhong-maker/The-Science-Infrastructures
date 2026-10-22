// The events layer: things that happened at a place and an instant, up to the
// time shown. Earthquakes are drawn by magnitude (size) and depth (colour), and
// fade with age over the window of the product. Disasters (GDACS) are drawn by
// their alert level (colour) and kind (size), from their start until a week
// after their last update.

import { MarkerLayer } from '../render/markers.js';

export const DEPTH_CLASSES = [
  { max: 70, label: 'shallow, < 70 km', colour: [246, 196, 126] },
  { max: 300, label: 'intermediate, 70–300 km', colour: [236, 143, 167] },
  { max: Infinity, label: 'deep, > 300 km', colour: [167, 148, 226] },
];

export const radiusOf = (mag) => Math.max(2.5, 3.2 + 3.3 * (mag - 4.5));

export const ALERTS = { Green: [143, 209, 180], Orange: [246, 196, 126], Red: [232, 132, 143] };
export const DISASTER_SIZE = { TC: 9, DR: 9, VO: 8, FL: 7, EQ: 6, WF: 3.5 };
const WEEK = 7 * 86400000;

export class EventsLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.markers = new MarkerLayer({ fill: 0.36 });
    stage.overlay.add(this.markers.mesh);
    this.def = null;
    this.shown = [];
  }

  async select(def) {
    this.def = def;
    this.product = null;
    this.markers.clear();
    if (!def) return;
    const p = await this.library.eventSet(def.product);
    if (this.def !== def) return;
    this.product = p;
    this.index = Object.fromEntries(p.fields.map((f, i) => [f, i]));
    this.key = null;
  }

  window() {
    if (!this.product) return null;
    const w = this.product.window;
    return [Date.parse(w.start), Date.parse(w.end)];
  }

  update(time, force = false) {
    if (!this.product) return;
    const key = `${Math.floor(time / 60000)}/${this.stage.mode}/${this.stage.lambda0}`;
    if (key === this.key && !force) return;
    this.key = key;
    const [w0, w1] = this.window();
    const span = Math.max(1, w1 - w0);
    const f = this.index, proj = this.stage.projector();
    const pos = [], col = [], size = [];
    this.shown = [];
    const p = [0, 0, 0];
    if (this.product.kind === 'disaster') {
      for (const e of this.product.events) {
        const t = e[f.time_ms], end = e[f.end_ms];
        if (t > time || time > end + WEEK) continue;
        const age = Math.max(0, time - end) / WEEK;
        if (!proj.project(e[f.lon], e[f.lat], p, 1.002)) continue;
        const c = ALERTS[e[f.alert]] || ALERTS.Green;
        pos.push(p[0], p[1], p[2]);
        col.push(c[0] / 255, c[1] / 255, c[2] / 255, 1 - 0.6 * age);
        size.push(DISASTER_SIZE[e[f.type]] || 6);
        this.shown.push(e);
      }
      this.markers.set(new Float32Array(pos), new Float32Array(col), new Float32Array(size));
      return;
    }
    for (const e of this.product.events) {
      const t = e[f.time_ms];
      if (t > time) continue;
      const age = (time - t) / span;
      if (age > 1) continue;
      proj.project(e[f.lon], e[f.lat], p, 1.002);
      const cls = DEPTH_CLASSES.find((c) => e[f.depth_km] < c.max);
      pos.push(p[0], p[1], p[2]);
      col.push(cls.colour[0] / 255, cls.colour[1] / 255, cls.colour[2] / 255, 1 - 0.65 * age);
      size.push(radiusOf(e[f.mag]));
      this.shown.push(e);
    }
    this.markers.set(new Float32Array(pos), new Float32Array(col), new Float32Array(size));
  }

  // The event nearest a screen point, within its radius plus a few pixels.
  pick(x, y) {
    if (!this.product) return null;
    const f = this.index;
    let best = null, bestD = Infinity;
    for (const e of this.shown) {
      const s = this.stage.screen(e[f.lon], e[f.lat]);
      if (!s) continue;
      const d = Math.hypot(s[0] - x, s[1] - y);
      const r = this.product.kind === 'disaster' ? DISASTER_SIZE[e[f.type]] || 6 : radiusOf(e[f.mag]);
      if (d < r + 4 && d < bestD) { best = e; bestD = d; }
    }
    return best && Object.fromEntries(this.product.fields.map((k, i) => [k, best[i]]));
  }

  frame(resolution) { this.markers.frame(resolution); }

  set visible(v) { this.markers.visible = v; }
}
