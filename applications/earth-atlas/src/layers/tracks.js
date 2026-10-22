// The tracks layer: tropical cyclones as their centres move, from several
// sources at once, so that their forecasts can be compared:
// - tracks found here in a model's sea-level pressure (src/core/cyclones.js),
//   one colour per model;
// - a tracks product, such as GDACS's: observed positions (solid) and the
//   official forecast (dashed).
// At the time shown, each cyclone's centre is marked, by interpolation between
// its positions, with its central pressure where the source gives one, and the
// name of the official cyclone it follows.

import { SegmentLayer, polylineSegments } from '../render/segments.js';
import { MarkerLayer } from '../render/markers.js';
import { trackProduct, arcDeg, positionAt } from '../core/cyclones.js';

const hex = (h) => [parseInt(h.slice(1, 3), 16) / 255, parseInt(h.slice(3, 5), 16) / 255, parseInt(h.slice(5, 7), 16) / 255];
const OFFICIAL = '#f6ece8';

// The name of the official cyclone a track follows: the one whose position at a
// common time is nearest, within 5 degrees.
function nameOf(track, official) {
  let best = null, bestD = 5;
  for (const o of official) {
    for (const p of track.points) {
      const q = positionAt(o.points, p.t);
      if (!q) continue;
      const d = arcDeg(p.lon, p.lat, q.lon, q.lat);
      if (d < bestD) { bestD = d; best = o.name; }
    }
  }
  return best;
}

export class TracksLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.sources = [];                // [{ def, tracks, lines: SegmentLayer[], dots }]
    this.cache = new Map();           // def id -> tracks
    this.marks = [];
  }

  // defs: the sources to show; isOcean(lon, lat): the land-sea test for genesis.
  async select(defs, isOcean) {
    const ticket = (this.ticket = {});
    const sources = [];
    for (const def of defs) {
      let tracks = this.cache.get(def.id);
      if (!tracks) {
        if (def.kind === 'official') {
          const p = await this.library.trackSet(def.product);
          tracks = p.tracks.map((tr) => ({ name: tr.name, alert: tr.alert, official: true,
            points: tr.points.map(([t, lon, lat, observed]) => ({ t, lon, lat, observed: Boolean(observed) })) }));
        } else {
          const product = await this.library.field(def.product);
          tracks = trackProduct(product, def.variable, { isOcean });
        }
        this.cache.set(def.id, tracks);
      }
      sources.push({ def, tracks });
    }
    if (this.ticket !== ticket) return;
    const official = sources.filter((s) => s.def.kind === 'official').flatMap((s) => s.tracks);
    for (const s of sources) if (s.def.kind !== 'official') for (const tr of s.tracks) tr.name = nameOf(tr, official);
    this.clear();
    this.sources = sources.map((s) => {
      const colour = s.def.kind === 'official' ? OFFICIAL : s.def.colour;
      const solid = new SegmentLayer({ color: colour, width: 1.6, opacity: 0.85, order: 3 });
      const dashed = new SegmentLayer({ color: colour, width: 1.3, opacity: 0.75, dash: 0.01, order: 3 });
      const dots = new MarkerLayer({ fill: 0.9 });
      for (const l of [solid.mesh, dashed.mesh, dots.mesh]) this.stage.overlay.add(l);
      return { ...s, colour, solid, dashed, dots };
    });
    this.rebuild();
  }

  clear() {
    for (const s of this.sources) {
      for (const l of [s.solid, s.dashed, s.dots]) { this.stage.overlay.remove(l.mesh); l.mesh.geometry.dispose(); }
    }
    this.sources = [];
    this.marks = [];
  }

  // Lines in world coordinates for the current projection.
  rebuild() {
    const proj = this.stage.projector();
    const lift = (lon, lat, out) => proj.project(lon, lat, out, 1.0016);
    for (const s of this.sources) {
      const solid = [], dashed = [], solidRuns = [], dashedRuns = [];
      const add = (coords, runs, pts) => { if (pts.length > 1) { runs.push([coords.length / 2, pts.length]); for (const p of pts) coords.push(p.lon, p.lat); } };
      for (const tr of s.tracks) {
        if (s.def.kind === 'official') {
          const obs = tr.points.filter((p) => p.observed), fc = tr.points.filter((p) => !p.observed);
          add(solid, solidRuns, obs);
          add(dashed, dashedRuns, obs.length ? [obs[obs.length - 1], ...fc] : fc);
        } else {
          add(solid, solidRuns, tr.points);
        }
      }
      const seg = (coords, runs) => polylineSegments(new Float32Array(coords), runs, lift, proj.breaks, { subdivide: 1 });
      const a = seg(solid, solidRuns), b = seg(dashed, dashedRuns);
      s.solid.set(a.a, a.b, a.info);
      s.dashed.set(b.a, b.b, b.info);
    }
    this.lastTime = null;
  }

  // Positions as dots, and each cyclone's centre at the time shown.
  update(time, force = false) {
    if (!force && this.lastTime === time) return;
    this.lastTime = time;
    const proj = this.stage.projector();
    const p = [0, 0, 0];
    this.marks = [];
    for (const s of this.sources) {
      const [r, g, b] = hex(s.colour);
      const pos = [], col = [], size = [];
      for (const tr of s.tracks) {
        for (const q of tr.points) {
          if (!proj.project(q.lon, q.lat, p, 1.0018)) continue;
          pos.push(p[0], p[1], p[2]); col.push(r, g, b, q.t <= time ? 0.85 : 0.4); size.push(1.6);
        }
        const now = positionAt(tr.points, time);
        if (now && proj.project(now.lon, now.lat, p, 1.002)) {
          // The centre: larger for a deeper cyclone (6 px at 1000 hPa, 14 px at 920 hPa).
          const hPa = now.pa !== undefined ? now.pa / 100 : undefined;
          pos.push(p[0], p[1], p[2]); col.push(r, g, b, 1);
          size.push(hPa !== undefined ? Math.max(5, Math.min(15, 6 + (1000 - hPa) / 10)) : 7);
          this.marks.push({ lon: now.lon, lat: now.lat, hPa, name: tr.name, source: s.def, colour: s.colour });
        }
      }
      s.dots.set(new Float32Array(pos), new Float32Array(col), new Float32Array(size));
    }
  }

  frame(resolution) {
    for (const s of this.sources) { s.solid.frame(resolution, 0); s.dashed.frame(resolution, 0); s.dots.frame(resolution); }
  }

  // The span of times covered, for the timeline.
  window() {
    const ts = this.sources.flatMap((s) => s.tracks.flatMap((tr) => [tr.points[0].t, tr.points[tr.points.length - 1].t]));
    return ts.length ? [Math.min(...ts), Math.max(...ts)] : null;
  }
}
