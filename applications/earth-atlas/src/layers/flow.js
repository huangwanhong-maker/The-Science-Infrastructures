// The flow layer: particles carried by a vector field, drawn as fading trails.
//
// Kinds of flow:
//   vector  the two components of a velocity in the product (wind, currents);
//   waves   the velocity at which the energy of ocean waves travels, derived
//           from the peak direction and period (src/core/derived.js).

import { ParticleSystem } from '../core/flow.js';
import { bilinear } from '../core/grid.js';
import { deepWaterGroupSpeed } from '../core/derived.js';
import { SegmentLayer } from '../render/segments.js';

const TRAIL = 16;            // segments kept per particle
const PX_PER_FRAME = 0.16;   // screen speed, pixels per frame for 1 m/s

export class FlowLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.layer = new SegmentLayer({ capacity: 1, color: '#fff3e6', width: 1.35, opacity: 0.82, trail: TRAIL - 1, order: 3 });
    stage.overlay.add(this.layer.mesh);
    this.def = null;
    this.frameNo = 0;
    this.derived = new Map();
    this.time = null;
  }

  async select(def) {
    this.def = def;
    this.product = null;
    this.layer.clear();
    if (!def) return;
    const product = await this.library.field(def.product);
    if (this.def !== def) return;
    this.product = product;
    this.reset();
  }

  reset() {
    const s = this.stage;
    const n = Math.round(Math.max(1500, Math.min(9000, (s.coverage() * s.width * s.height) / 105)));
    this.particles = new ParticleSystem(n);
    this.particles.reset(s.spawner());
    this.layer.allocate(n * TRAIL);
    this.layer.geometry.instanceCount = n * TRAIL;
    this.layer.info.array.fill(0);
    this.layer.info.needsUpdate = true;
    this.a = new Float32Array(3 * n); this.b = new Float32Array(3 * n); this.info = new Float32Array(4 * n);
    this.frameNo = TRAIL;
  }

  // The velocity components at the time indices of a bracket, as grids.
  grids(i) {
    if (this.def.kind === 'vector') return [this.product.slice(this.def.variables[0], i), this.product.slice(this.def.variables[1], i)];
    const key = `${this.def.product}/${i}`;
    if (!this.derived.has(key)) {
      const dir = this.product.slice(this.def.variables[0], i), per = this.product.slice(this.def.variables[1], i);
      const u = new Float32Array(dir.length), v = new Float32Array(dir.length);
      for (let k = 0; k < dir.length; k++) {
        const c = deepWaterGroupSpeed(per[k]), a = (dir[k] * Math.PI) / 180;
        u[k] = -c * Math.sin(a); v[k] = -c * Math.cos(a);
      }
      this.derived.set(key, [u, v]);
      if (this.derived.size > 6) this.derived.delete(this.derived.keys().next().value);
    }
    return this.derived.get(key);
  }

  update(time) {
    if (!this.product) return;
    this.time = time;
    const b = this.product.bracket(time);
    const g = this.product.grid;
    const [u0, v0] = this.grids(Math.max(0, b.i0));
    const [u1, v1] = this.grids(Math.max(0, b.i1));
    const w = b.w;
    const out = [0, 0];
    this.sample = (lon, lat) => {
      let u = bilinear(u0, g, lon, lat), v = bilinear(v0, g, lon, lat);
      if (w > 0) {
        u += (bilinear(u1, g, lon, lat) - u) * w;
        v += (bilinear(v1, g, lon, lat) - v) * w;
      }
      if (!(u === u && v === v)) return null;
      out[0] = u; out[1] = v;
      return out;
    };
  }

  step() {
    if (!this.product || !this.sample || !this.particles) return;
    const s = this.stage, ps = this.particles;
    // Particles move PX_PER_FRAME pixels a frame for one unit of the field, times the quantity's pace.
    const dt = (PX_PER_FRAME * (this.def.pace || 1)) / s.pixelsPerMetre();
    ps.step(this.sample, dt, this.spawn || (this.spawn = s.spawner()));
    const proj = s.projector();
    const p = [0, 0, 0], q = [0, 0, 0];
    const n = ps.count;
    for (let i = 0; i < n; i++) {
      const seg = ps.segments;
      let alpha = 0;
      if (ps.alive[i]) {
        proj.project(seg[4 * i], seg[4 * i + 1], p);
        proj.project(seg[4 * i + 2], seg[4 * i + 3], q);
        alpha = proj.breaks(p, q) ? 0 : Math.min(1, 0.3 + ps.speed[i] / 14);
      }
      this.a.set(p, 3 * i); this.b.set(q, 3 * i);
      this.info[4 * i] = this.frameNo; this.info[4 * i + 1] = alpha;
    }
    const slot = this.frameNo % TRAIL;
    this.layer.write(slot * n, n, this.a, this.b, this.info);
    this.frameNo++;
  }

  // Called when the view changes: new places to be born in.
  viewChanged() { this.spawn = null; }

  frame(resolution) { this.layer.frame(resolution, this.frameNo - 1); }

  set visible(v) { this.layer.visible = v; }
}
