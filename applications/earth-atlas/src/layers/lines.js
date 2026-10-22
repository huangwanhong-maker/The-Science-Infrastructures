// The lines layer: isolines of a field at an instant, such as isobars, with
// the field's highs and lows marked.

import { isolines, extrema } from '../core/contour.js';
import { SegmentLayer, pairSegments } from '../render/segments.js';
import { displayUnit, toShown } from '../core/units.js';

export class LinesLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.minor = new SegmentLayer({ color: '#fbefe9', width: 1.0, opacity: 0.42, order: 2 });
    this.major = new SegmentLayer({ color: '#fbefe9', width: 1.6, opacity: 0.62, order: 2 });
    stage.overlay.add(this.minor.mesh, this.major.mesh);
    this.def = null;
    this.marks = [];
    this.key = null;
  }

  async select(def) {
    this.def = def;
    this.product = null;
    this.key = null;
    this.minor.clear(); this.major.clear();
    this.marks = [];
    if (!def) return;
    const product = await this.library.field(def.product);
    if (this.def !== def) return;
    this.product = product;
  }

  // Recomputes the lines for an instant (and when the projection changes).
  update(time, force = false) {
    if (!this.product) return;
    const b = this.product.bracket(time);
    const key = `${b.i0}/${b.i1}/${b.w.toFixed(1)}/${this.stage.mode}/${this.stage.lambda0}`;
    if (key === this.key && !force) return;
    this.key = key;
    const d = this.def, g = this.product.grid;
    const a = this.product.slice(d.variable, Math.max(0, b.i0));
    let values = a;
    if (b.w > 0) {
      const c = this.product.slice(d.variable, b.i1);
      values = new Float32Array(a.length);
      for (let i = 0; i < a.length; i++) values[i] = a[i] + (c[i] - a[i]) * b.w;
    }
    let lo = Infinity, hi = -Infinity;
    for (const x of values) if (x === x) { if (x < lo) lo = x; if (x > hi) hi = x; }
    const minor = [], major = [];
    for (let k = Math.ceil((lo - d.base) / d.step); d.base + k * d.step <= hi; k++) {
      const seg = isolines(values, g, d.base + k * d.step);
      (k % 5 === 0 ? major : minor).push(seg);
    }
    const proj = this.stage.projector();
    const join = (list) => {
      const n = list.reduce((s, x) => s + x.length, 0);
      const all = new Float32Array(n);
      let o = 0;
      for (const x of list) { all.set(x, o); o += x.length; }
      return pairSegments(all, (lon, lat, out) => proj.project(lon, lat, out, 1.0012), proj.breaks);
    };
    const m1 = join(minor), m2 = join(major);
    this.minor.set(m1.a, m1.b, m1.info);
    this.major.set(m2.a, m2.b, m2.info);
    const unit = displayUnit(d.units, d.labelUnits);
    this.marks = d.extrema
      ? extrema(values, g, d.extrema).map((e) => ({ ...e, text: `${e.kind === 'high' ? 'H' : 'L'}`, value: Math.round(toShown(e.value, unit)) }))
      : [];
  }

  frame(resolution) { this.minor.frame(resolution, 0); this.major.frame(resolution, 0); }

  set visible(v) { this.minor.visible = this.major.visible = v; }
}
