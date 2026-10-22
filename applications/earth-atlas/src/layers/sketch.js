// Sketches: polylines, markers, labels and moving points given in longitude
// and latitude, drawn over the globe as they are. They come from plug-ins,
// whose results are data; this layer is the one place that turns them into
// geometry. Polylines are grouped by style into segment layers; markers and the
// moving points at the time shown share one marker layer.

import { SegmentLayer, polylineSegments } from '../render/segments.js';
import { MarkerLayer } from '../render/markers.js';

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// '#rrggbb' -> [r, g, b] in 0..1.
function rgb(colour) {
  const m = /^#([0-9a-f]{6})$/i.exec(colour || '');
  const n = m ? parseInt(m[1], 16) : 0xe8c3dc;
  return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
}

// Runs of a polyline that do not cross the antimeridian: { coords, runs }.
function splitAtDateline(coords) {
  const flat = Float32Array.from(coords), runs = [];
  let first = 0;
  const n = flat.length / 2;
  for (let i = 1; i < n; i++) {
    if (Math.abs(flat[2 * i] - flat[2 * i - 2]) > 180) { if (i - first > 1) runs.push([first, i - first]); first = i; }
  }
  if (n - first > 1) runs.push([first, n - first]);
  return { coords: flat, runs };
}

// Where a moving point is at `time`: linear between its times, nowhere outside them.
function positionAt(m, time) {
  const t = m.times, n = t.length;
  if (!n || time < t[0] || time > t[n - 1]) return null;
  let i = 0;
  while (i + 1 < n - 1 && t[i + 1] <= time) i++;
  const w = t[i + 1] > t[i] ? (time - t[i]) / (t[i + 1] - t[i]) : 0;
  let lon0 = m.coords[2 * i], lon1 = m.coords[2 * i + 2];
  if (lon1 - lon0 > 180) lon0 += 360; else if (lon0 - lon1 > 180) lon1 += 360;
  return [lon0 + (lon1 - lon0) * w, m.coords[2 * i + 1] + (m.coords[2 * i + 3] - m.coords[2 * i + 1]) * w];
}

const concat = (parts) => {
  const out = new Float32Array(parts.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of parts) { out.set(p, at); at += p.length; }
  return out;
};

export class SketchLayer {
  constructor(stage) {
    this.stage = stage;
    this.sketches = [];
    this.lines = new Map();          // style key -> SegmentLayer
    this.markers = new MarkerLayer({ capacity: 512, fill: 0.55, order: 5 });
    stage.overlay.add(this.markers.mesh);
    this.time = 0;
    this.shown = true;
  }

  // Replaces the sketches drawn.
  set(sketches) {
    this.sketches = sketches || [];
    this.rebuild();
    this.update(this.time, true);
  }

  clear() { this.set([]); }

  // Geometry for the projection in force.
  rebuild() {
    const proj = this.stage.projector();
    const lift = (lon, lat, out) => proj.project(lon, lat, out, 1.0024);
    const groups = new Map();
    for (const s of this.sketches) {
      for (const p of s.polylines || []) {
        const key = `${p.colour || '#e8c3dc'}|${p.width || 1.2}|${p.dashed ? 1 : 0}`;
        if (!groups.has(key)) groups.set(key, []);
        groups.get(key).push(p);
      }
    }
    for (const [key, layer] of this.lines) if (!groups.has(key)) layer.set(new Float32Array(0), new Float32Array(0), new Float32Array(0));
    for (const [key, list] of groups) {
      let layer = this.lines.get(key);
      if (!layer) {
        const [colour, width, dashed] = key.split('|');
        layer = new SegmentLayer({ color: colour, width: Number(width), opacity: 0.95, dash: dashed === '1' ? 8 : 0, order: 3 });
        layer.visible = this.shown;
        this.stage.overlay.add(layer.mesh);
        this.lines.set(key, layer);
      }
      const a = [], b = [], info = [];
      for (const p of list) {
        const { coords, runs } = splitAtDateline(p.coords);
        const seg = polylineSegments(coords, runs, lift, proj.breaks, { alpha: p.alpha ?? 1, subdivide: 0.5 });
        a.push(seg.a); b.push(seg.b); info.push(seg.info);
      }
      layer.set(concat(a), concat(b), concat(info));
    }
  }

  // Markers at an instant: the fixed ones, and each moving point where it is then.
  update(time, force = false) {
    if (!force && time === this.time) return;
    this.time = time;
    const proj = this.stage.projector(), p = [0, 0, 0], pos = [], col = [], size = [];
    const put = (lon, lat, colour, s) => {
      if (!proj.project(lon, lat, p, 1.0026)) return;
      const [r, g, b] = rgb(colour);
      pos.push(p[0], p[1], p[2]); col.push(r, g, b, 1); size.push(s);
    };
    for (const s of this.sketches) {
      for (const m of s.markers || []) put(m.lon, m.lat, m.colour, m.size || 3);
      for (const m of s.moving || []) { const at = positionAt(m, time); if (at) put(at[0], at[1], m.colour, m.size || 3); }
    }
    this.markers.set(new Float32Array(pos), new Float32Array(col), new Float32Array(size));
  }

  frame(resolution) {
    for (const layer of this.lines.values()) layer.frame(resolution, 0);
    this.markers.frame(resolution);
  }

  // Items for the label layer.
  labelItems() {
    return this.sketches.flatMap((s, k) => (s.labels || []).map((l, i) => ({
      key: `sketch${k}-${i}`, text: l.text, html: esc(l.text), lon: l.lon, lat: l.lat, className: 'sketch', font: '12px serif', size: 18, always: true,
    })));
  }

  set visible(v) {
    this.shown = v;
    for (const layer of this.lines.values()) layer.visible = v;
    this.markers.visible = v;
  }

  get visible() { return this.shown; }
}
