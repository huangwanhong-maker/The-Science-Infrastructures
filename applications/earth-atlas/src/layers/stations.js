// The GPS stations layer: how fast and which way the ground moves, measured at
// some 20,000 continuous GNSS stations (MIDAS velocities, Nevada Geodetic
// Laboratory). Each station's horizontal velocity is an arrow and its vertical
// rate colours a dot beneath it (red rising, blue sinking).
//
// Two frames: IGS20 (aligned with ITRF2020), in which whole plates move by
// centimetres a year; or relative to the plate the station stands on, in which
// a rigid plate stands still and what remains is the deformation near plate
// boundaries, the rebound of the land from the ice ages, and the stations' own
// errors. The plate's motion subtracted is that of ITRF2020's plate motion
// model, in the stations' own frame, for the 13 plates it covers, and that of
// NNR-MORVEL56 for the others (where a frame difference of a millimetre or two
// a year remains).
//
// Arrows keep their length on screen, so many millimetres a year are many
// pixels at any zoom; where stations crowd, one per cell of the screen is
// drawn, the one with the longest record.

import { SegmentLayer, pairSegments } from '../render/segments.js';
import { MarkerLayer } from '../render/markers.js';
import { colourAt } from '../render/colormaps.js';
import { velocityAt, itrfVelocityAt } from '../core/plates.js';
import { compassPoint } from '../core/derived.js';
import { flowToward } from '../core/ocean.js';

const CELL = 16;                                 // px: one station per cell where they crowd
const M_PER_DEG = 111195;
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

export class StationsLayer {
  constructor(stage, library, plates) {
    this.stage = stage;
    this.library = library;
    this.plates = plates;
    this.shown = false;
    this.mode = 'plate';
    this.arrows = new SegmentLayer({ color: '#fff3e6', width: 1.25, opacity: 0.92, order: 4 });
    this.dots = new MarkerLayer({ fill: 0.95, order: 5 });
    for (const l of [this.arrows, this.dots]) { l.mesh.visible = false; stage.overlay.add(l.mesh); }
    this.drawn = [];
  }

  // px per mm/yr: whole plates move far in IGS20; deformation is small.
  get scale() { return this.mode === 'plate' ? 3.2 : 0.75; }

  async show(on) {
    this.shown = on;
    if (on && !this.product) {
      const entry = this.library.list('stations')[0];
      if (!entry) { this.shown = false; return; }
      this.product = await this.library.stationSet(entry.name);
      const f = Object.fromEntries(this.product.fields.map((k, i) => [k, i]));
      this.f = f;
      this.list = this.product.stations;
      // Each station's plate and the plate's motion there, for the plate-fixed frame.
      this.residual = this.list.map((r) => {
        const p = this.plates.plates ? this.plates.plateAt(r[f.lon], r[f.lat]) : null;
        if (!p || !p.omega) return null;
        const itrf = itrfVelocityAt(p.code, r[f.lon], r[f.lat]);
        const v = itrf || velocityAt(p.omega, r[f.lon], r[f.lat]);
        return { plate: p, model: itrf ? 'ITRF2020-PMM' : 'NNR-MORVEL56', motion: v, east: r[f.ve] - v.east, north: r[f.vn] - v.north };
      });
      this.order = this.list.map((_, i) => i).sort((a, b) => (this.list[b][f.last] - this.list[b][f.first]) - (this.list[a][f.last] - this.list[a][f.first]));
    }
    for (const l of [this.arrows, this.dots]) l.mesh.visible = this.shown && Boolean(this.product);
    this.stale = true;
  }

  setMode(mode) { this.mode = mode; this.stale = true; }

  velocity(i) {
    const r = this.list[i], f = this.f;
    if (this.mode === 'plate') { const x = this.residual[i]; return x ? [x.east, x.north] : null; }
    return [r[f.ve], r[f.vn]];
  }

  // Arrows for the current view: called when the view has settled.
  layout() {
    this.stale = false;
    if (!this.shown || !this.product) return;
    const s = this.stage, f = this.f, rows = this.list;
    const taken = new Set(), drawn = [];
    for (const i of this.order) {
      const r = rows[i];
      const at = s.screen(r[f.lon], r[f.lat]);
      if (!at || at[2] <= 0.2 || at[0] < -20 || at[1] < -20 || at[0] > s.width + 20 || at[1] > s.height + 20) continue;
      const key = `${Math.floor(at[0] / CELL)},${Math.floor(at[1] / CELL)}`;
      if (taken.has(key)) continue;
      const v = this.velocity(i);
      if (!v) continue;
      taken.add(key);
      drawn.push({ i, x: at[0], y: at[1] });
    }
    this.drawn = drawn;
    const mPerPx = 1 / s.pixelsPerMetre();
    const pairs = [], pos = [], col = [], size = [];
    const proj = s.projector(), p = [0, 0, 0];
    for (const d of drawn) {
      const r = rows[d.i], [ve, vn] = this.velocity(d.i);
      const lon = r[f.lon], lat = r[f.lat];
      const speed = Math.hypot(ve, vn);
      const length = speed * this.scale * mPerPx;                   // metres on the ground for the arrow
      const cos = Math.max(0.05, Math.cos((lat * Math.PI) / 180));
      const tip = [lon + (ve / Math.max(speed, 1e-9)) * length / (M_PER_DEG * cos), lat + (vn / Math.max(speed, 1e-9)) * length / M_PER_DEG];
      if (speed * this.scale >= 2) {
        pairs.push(lon, lat, tip[0], tip[1]);
        // The head: two strokes back from the tip at 26 degrees.
        const head = Math.min(9, Math.max(3.5, speed * this.scale * 0.3)) * mPerPx;
        const back = Math.atan2(-vn, -ve);
        for (const turn of [0.45, -0.45]) {
          const a = back + turn;
          pairs.push(tip[0], tip[1], tip[0] + (Math.cos(a) * head) / (M_PER_DEG * cos), tip[1] + (Math.sin(a) * head) / M_PER_DEG);
        }
      }
      if (proj.project(lon, lat, p, 1.002)) {
        const t = 0.5 + Math.max(-1, Math.min(1, r[f.vu] / 6)) * 0.5;
        const [cr, cg, cb] = colourAt('balance', t);
        pos.push(p[0], p[1], p[2]); col.push(cr / 255, cg / 255, cb / 255, 1); size.push(3.2);
      }
    }
    const lift = (lon, lat, out) => proj.project(lon, lat, out, 1.0016);
    const seg = pairSegments(new Float32Array(pairs), lift, proj.breaks);
    this.arrows.set(seg.a, seg.b, seg.info);
    this.dots.set(new Float32Array(pos), new Float32Array(col), new Float32Array(size));
  }

  frame(resolution) {
    if (!this.shown) return;
    this.arrows.frame(resolution, 0);
    this.dots.frame(resolution);
  }

  // The drawn station nearest a point on the screen, within r pixels.
  pick(x, y, r = 9) {
    let best = null, bestD = r;
    for (const d of this.drawn) {
      const dist = Math.hypot(d.x - x, d.y - y);
      if (dist < bestD) { bestD = dist; best = d.i; }
    }
    return best === null ? null : this.station(best);
  }

  station(i) {
    const r = this.list[i], f = this.f;
    return { index: i, id: r[f.id], lon: r[f.lon], lat: r[f.lat], height: r[f.height_m], ve: r[f.ve], vn: r[f.vn], vu: r[f.vu],
      se: r[f.se], sn: r[f.sn], su: r[f.su], first: r[f.first], last: r[f.last], residual: this.residual[i] };
  }

  find(id) {
    if (!this.list) return null;
    const i = this.list.findIndex((r) => r[this.f.id] === id);
    return i < 0 ? null : this.station(i);
  }

  // Rows for the place card: [term, html, kind].
  rows(st) {
    const out = [['GPS station', `${esc(st.id)} <span class="dir">${st.first.toFixed(1)}–${st.last.toFixed(1)}</span>`]];
    const h = Math.hypot(st.ve, st.vn);
    out.push(['Moves', `${h.toFixed(1)} mm/yr <span class="dir">toward ${compassPoint(flowToward(st.ve, st.vn))}, ±${Math.max(st.se, st.sn).toFixed(1)}</span>`]);
    out.push(['', 'in IGS20, the frame of ITRF2020 (MIDAS; Blewitt et al. 2016)', 'gloss']);
    if (st.residual) {
      const x = st.residual, r = Math.hypot(x.east, x.north);
      out.push([`Relative to the ${esc(st.residual.plate.name)} plate`, `${r.toFixed(1)} mm/yr <span class="dir">toward ${compassPoint(flowToward(x.east, x.north))}</span>`]);
      const model = x.model === 'ITRF2020-PMM' ? 'ITRF2020’s plate motion model' : 'NNR-MORVEL56, as ITRF2020’s model does not cover this plate';
      out.push(['', `${r < (x.model === 'ITRF2020-PMM' ? 1.5 : 3) ? 'within the plate’s rigid motion: the ground here is hardly deforming' : 'what the plate’s rigid motion leaves: the ground deforms here'} (${model})`, 'gloss']);
    }
    out.push(['Vertical', `${st.vu >= 0 ? 'rising' : 'sinking'} ${Math.abs(st.vu).toFixed(1)} mm/yr <span class="dir">±${st.su.toFixed(1)}</span>`]);
    out.push(['Daily positions', `<button type="button" class="action" data-action="station-series">Show the record</button>`]);
    return out;
  }

  legend() {
    return {
      title: 'GPS stations', count: `${this.drawn.length.toLocaleString('en-GB')} shown`,
      keys: [{ label: `${this.mode === 'plate' ? 5 : 20} mm/yr`, colour: '#fff3e6', length: (this.mode === 'plate' ? 5 : 20) * this.scale }],
      dots: [['rising', colourAt('balance', 1)], ['steady', colourAt('balance', 0.5)], ['sinking', colourAt('balance', 0)]],
      note: `${this.mode === 'plate' ? 'relative to each station’s plate (ITRF2020’s model, or NNR-MORVEL56)' : 'in IGS20 (ITRF2020)'} · MIDAS velocities, Nevada Geodetic Laboratory`,
    };
  }
}
