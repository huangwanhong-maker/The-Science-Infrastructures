// The plates layer: the tectonic plates of PB2002 (Bird 2003), their boundaries
// drawn by class (spreading ridges and rifts, transform faults, convergent
// boundaries, and subduction zones with teeth on the overriding plate), their
// names, and their motion, from the rotations of NNR-MORVEL56 (Argus et al.
// 2011; src/core/plates.js).
//
// The motion is offered as a field product computed here (east and north
// velocity, mm/yr, every half degree), so that the flow and surface layers draw
// it like any other vector field: plates carried as particles, or coloured by
// speed.

import { SegmentLayer, polylineSegments } from '../render/segments.js';
import { ComputedField } from '../core/grid.js';
import { fillRings, featureAreaFraction } from '../core/raster.js';
import { poleOf, omega, velocityAt, subtract, acrossStep, azimuth, inRings, morvelCode } from '../core/plates.js';
import { arcDeg } from '../core/cyclones.js';
import { compassPoint } from '../core/derived.js';
import { flowToward } from '../core/ocean.js';
import { travel as offsetAlong, relativeTo } from '../core/section.js';

export const MOTION_PRODUCT = 'plate-motion-nnr-morvel56';
export const GROUPS = [
  { id: 'divergent', label: 'spreading ridge or rift', classes: ['OSR', 'CRB'], colour: '#f6c47e', width: 1.7 },
  { id: 'transform', label: 'transform fault', classes: ['OTF', 'CTF'], colour: '#e9d6e4', width: 1.3 },
  { id: 'convergent', label: 'convergent boundary', classes: ['OCB', 'CCB'], colour: '#e8a0a0', width: 1.7 },
  { id: 'subduction', label: 'subduction zone', classes: ['SUB'], colour: '#ef8a98', width: 1.9 },
];
const W = 1440, H = 720;                         // the plates' raster: a quarter of a degree
const KM = 111.195;                              // km per degree of great circle (mean Earth radius)
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// A point at distances east and north (km) of a place, on a sphere; enough for teeth tens of km long.
const offset = (lon, lat, east, north) => [lon + east / (KM * Math.cos((lat * Math.PI) / 180)), lat + north / KM];

export class PlatesLayer {
  constructor(stage, library) {
    this.stage = stage;
    this.library = library;
    this.shown = false;
    this.lines = GROUPS.map((g) => {
      const l = new SegmentLayer({ color: g.colour, width: g.width, opacity: 0.9, order: 3 });
      l.mesh.visible = false;
      stage.overlay.add(l.mesh);
      return l;
    });
    this.teeth = new SegmentLayer({ color: GROUPS[3].colour, width: 1.4, opacity: 0.9, order: 3 });
    this.teeth.mesh.visible = false;
    stage.overlay.add(this.teeth.mesh);
  }

  // Reads the plates, rasterises them, and computes their motion.
  async load() {
    if (this.product !== undefined) return this.product;
    const set = await this.library.plates();
    if (!set) { this.product = null; return null; }
    const { meta, coords } = set;
    this.meta = meta;
    this.coords = coords;
    this.plates = meta.plates.map((p) => {
      let w = 180, e = -180, s = 90, n = -90;
      for (const [first, count] of p.rings) {
        for (let k = first; k < first + count; k++) {
          const x = coords[2 * k], y = coords[2 * k + 1];
          w = Math.min(w, x); e = Math.max(e, x); s = Math.min(s, y); n = Math.max(n, y);
        }
      }
      const pole = poleOf(p.code);
      return { ...p, box: [w, s, e, n], area: featureAreaFraction(coords, p.rings) * 4 * Math.PI, pole, omega: pole ? omega(pole) : null };
    });
    this.classes = meta.classes.map((c) => c.code);
    this.raster = new Uint8Array(W * H);
    this.plates.forEach((p, i) => fillRings(coords, p.rings, W, H, i + 1, this.raster));
    this.labels = this.labelPoints();
    this.product = this.motionField();
    this.library.addField(this.product);
    return this.product;
  }

  // The plate at a place: exact, by its outline.
  plateAt(lon, lat) {
    for (const p of this.plates || []) {
      const [w, s, e, n] = p.box;
      if (lon < w || lon > e || lat < s || lat > n) continue;
      if (inRings(this.coords, p.rings, lon, lat)) return p;
    }
    // On a shared edge, the raster decides.
    return this.plateFromRaster(lon, lat);
  }

  plateFromRaster(lon, lat) {
    const i = Math.min(W - 1, Math.max(0, Math.floor(((lon + 180) / 360) * W)));
    const j = Math.min(H - 1, Math.max(0, Math.floor(((90 - lat) / 180) * H)));
    const k = this.raster[j * W + i];
    return k ? this.plates[k - 1] : null;
  }

  // East and north velocity (mm/yr, no-net-rotation frame) every half degree.
  motionField() {
    const grid = { lon0: -180, dlon: 0.5, nlon: 720, lat0: 90, dlat: 0.5, nlat: 361, periodic: true };
    const u = new Float32Array(grid.nlon * grid.nlat), v = new Float32Array(grid.nlon * grid.nlat);
    for (let r = 0; r < grid.nlat; r++) {
      const lat = Math.max(-89.99, Math.min(89.99, grid.lat0 - r * grid.dlat));
      for (let c = 0; c < grid.nlon; c++) {
        const lon = grid.lon0 + c * grid.dlon;
        const p = this.plateFromRaster(lon, lat);
        const vel = p && p.omega ? velocityAt(p.omega, lon, lat) : { east: NaN, north: NaN };
        u[r * grid.nlon + c] = vel.east; v[r * grid.nlon + c] = vel.north;
      }
    }
    const stats = (a) => { const s = [...a].filter(Number.isFinite).sort((x, y) => x - y); return { p02: s[Math.floor(0.02 * s.length)], p98: s[Math.floor(0.98 * s.length)] }; };
    const meta = {
      schema: 'field-grid/1', name: MOTION_PRODUCT, computed: true,
      title: 'Motion of the tectonic plates, NNR-MORVEL56 on PB2002', model: 'NNR-MORVEL56', domain: 'geophysics',
      grid, times: [], time_scale: 'UTC',
      variables: [
        { name: 'plate_u', long_name: 'eastward velocity of the plate', standard_name: 'eastward_plate_velocity', units: 'mm yr-1', ...stats(u) },
        { name: 'plate_v', long_name: 'northward velocity of the plate', standard_name: 'northward_plate_velocity', units: 'mm yr-1', ...stats(v) },
      ],
      vectors: [{ name: 'plate_velocity', u: 'plate_u', v: 'plate_v', units: 'mm yr-1', standard_name: 'plate_velocity' }],
      source: {
        dataset: 'pb2002', provider: 'computed in the page: the rotations of NNR-MORVEL56 (Argus, Gordon and DeMets 2011) applied to the plates of PB2002 (Bird 2003)',
        credit: 'NNR-MORVEL56 (Argus et al. 2011) on PB2002 (Bird 2003)',
        licence: this.meta.source.licence, citation: `Argus, D. F., Gordon, R. G. and DeMets, C. (2011), Geochem. Geophys. Geosyst. 12, Q11001; ${this.meta.source.citation}`,
      },
    };
    return new ComputedField(meta, { plate_u: [u], plate_v: [v] });
  }

  // For each plate, the place farthest from its edges (a chamfer distance
  // transform of the raster, east-west steps shortened by the cosine of
  // latitude), where its name is written.
  labelPoints() {
    const d = new Float32Array(W * H).fill(1e9);
    const id = this.raster;
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const k = id[j * W + i], e = id[j * W + ((i + 1) % W)], s = j + 1 < H ? id[(j + 1) * W + i] : k;
        if (k !== e || k !== s) d[j * W + i] = 0;
      }
    }
    const cosAt = (j) => Math.max(0.05, Math.cos(((90 - ((j + 0.5) * 180) / H) * Math.PI) / 180));
    for (let pass = 0; pass < 2; pass++) {
      const order = pass === 0 ? [...Array(H).keys()] : [...Array(H).keys()].reverse();
      for (const j of order) {
        const cx = cosAt(j), dir = pass === 0 ? 1 : -1, back = j - dir;
        for (let n = 0; n < W; n++) {
          const i = pass === 0 ? n : W - 1 - n, at = j * W + i;
          const im = (i - dir + W) % W;
          let best = Math.min(d[at], d[j * W + im] + cx);
          if (back >= 0 && back < H) best = Math.min(best, d[back * W + i] + 1, d[back * W + im] + Math.hypot(1, cx));
          d[at] = best;
        }
      }
    }
    const far = this.plates.map(() => ({ d: -1, i: 0, j: 0 }));
    for (let j = 0; j < H; j++) {
      for (let i = 0; i < W; i++) {
        const k = id[j * W + i];
        if (k && d[j * W + i] > far[k - 1].d) far[k - 1] = { d: d[j * W + i], i, j };
      }
    }
    return this.plates.map((p, k) => ({ plate: p, lon: -180 + ((far[k].i + 0.5) * 360) / W, lat: 90 - ((far[k].j + 0.5) * 180) / H }));
  }

  show(on) {
    this.shown = on;
    if (on && !this.built && this.plates) this.rebuild();
    for (const l of [...this.lines, this.teeth]) l.mesh.visible = on && Boolean(this.plates);
  }

  // Lines in world coordinates for the current projection.
  rebuild() {
    if (!this.plates) return;
    this.built = true;
    const proj = this.stage.projector();
    const lift = (lon, lat, out) => proj.project(lon, lat, out, 1.0012);
    const fields = this.meta.steps;
    GROUPS.forEach((g, gi) => {
      const coords = [], runs = [];
      for (const s of fields) {
        if (!g.classes.includes(this.classes[s[4]])) continue;
        let [lon1, lat1, lon2, lat2] = s;
        // A step across the antimeridian is drawn the short way.
        if (lon2 - lon1 > 180) lon2 -= 360; else if (lon1 - lon2 > 180) lon2 += 360;
        runs.push([coords.length / 2, 2]);
        coords.push(lon1, lat1, lon2, lat2);
      }
      const seg = polylineSegments(new Float32Array(coords), runs, lift, proj.breaks, { subdivide: 0.5 });
      this.lines[gi].set(seg.a, seg.b, seg.info);
    });
    // Teeth: on every other subduction step, a triangle pointing toward the
    // overriding plate, set on its side of the trench.
    const coords = [], runs = [];
    let n = 0;
    for (const s of fields) {
      if (this.classes[s[4]] !== 'SUB' || !s[7] || n++ % 2) continue;
      let [lon1, lat1, lon2, lat2] = s;
      if (lon2 - lon1 > 180) lon2 -= 360; else if (lon1 - lon2 > 180) lon2 += 360;
      const lon = (lon1 + lon2) / 2, lat = (lat1 + lat2) / 2;
      const az = (azimuth(lon1, lat1, lon2, lat2) * Math.PI) / 180;
      const length = arcDeg(lon1, lat1, lon2, lat2) * KM;
      const half = Math.min(9, length * 0.22), tall = Math.min(14, length * 0.32);
      const side = s[7] === 2 ? 1 : -1;                     // right of travel (2) or left (1)
      const t = [Math.sin(az), Math.cos(az)], nrm = [side * Math.cos(az), -side * Math.sin(az)];
      const a = offset(lon, lat, -t[0] * half, -t[1] * half);
      const b = offset(lon, lat, t[0] * half, t[1] * half);
      const c = offset(lon, lat, nrm[0] * tall, nrm[1] * tall);
      runs.push([coords.length / 2, 4]);
      coords.push(...a, ...c, ...b, ...a);
    }
    const seg = polylineSegments(new Float32Array(coords), runs, lift, proj.breaks);
    this.teeth.set(seg.a, seg.b, seg.info);
  }

  frame(resolution) {
    if (!this.shown) return;
    for (const l of [...this.lines, this.teeth]) l.frame(resolution, 0);
  }

  // The boundary step nearest a place: { step, km }, by distance to the step's arc.
  nearestStep(lon, lat) {
    let best = null, bestD = Infinity;
    for (const s of this.meta.steps) {
      const d = Math.min(arcDeg(lon, lat, s[0], s[1]), arcDeg(lon, lat, s[2], s[3]),
        arcDeg(lon, lat, (s[0] + s[2]) / 2 + (Math.abs(s[2] - s[0]) > 180 ? 180 : 0), (s[1] + s[3]) / 2));
      if (d < bestD) { bestD = d; best = s; }
    }
    return best ? { step: best, km: bestD * KM } : null;
  }

  // The relative motion across a step: the right-hand plate relative to the left.
  across(step) {
    const [lon1, lat1, lon2, lat2, , left, right] = step;
    const L = this.plates[left], R = this.plates[right];
    if (!L.omega || !R.omega) return null;
    let mlon = (lon1 + lon2) / 2;
    if (Math.abs(lon2 - lon1) > 180) mlon += 180;
    const mlat = (lat1 + lat2) / 2;
    const v = velocityAt(subtract(R.omega, L.omega), mlon, mlat);
    return { left: L, right: R, v, ...acrossStep(v, azimuth(lon1, lat1, lon2, lat2)) };
  }

  // Where a line of samples [{ s, lon, lat }] crosses plate boundaries: [{ s, label, colour }].
  crossings(points) {
    const out = [];
    for (let i = 1; i < points.length; i++) {
      const a = this.plateFromRaster(points[i - 1].lon, points[i - 1].lat), b = this.plateFromRaster(points[i].lon, points[i].lat);
      if (!a || !b || a === b) continue;
      const lon = (points[i - 1].lon + points[i].lon) / 2, lat = (points[i - 1].lat + points[i].lat) / 2;
      const near = this.nearestStep(lon, lat);
      const code = near ? this.classes[near.step[4]] : null;
      const group = GROUPS.find((g) => g.classes.includes(code)) || GROUPS[1];
      out.push({ s: (points[i - 1].s + points[i].s) / 2, label: code ? this.meta.classes[near.step[4]].name : `${a.name} | ${b.name}`, colour: group.colour });
    }
    return out;
  }

  // A section across the nearest boundary: perpendicular to it, toward the
  // overriding plate at a subduction zone (so that the sinking plate dips to
  // the right), from 400 km before the boundary to 1,200 km beyond; else east
  // through the place. { lon, lat, az, half, why, mark }.
  sectionAt(lon, lat) {
    const near = this.plates ? this.nearestStep(lon, lat) : null;
    if (!near || near.km > 900) return { lon, lat, az: 90, half: 700, why: 'east–west through the place chosen', mark: 0 };
    const [lon1, lat1, lon2, lat2, k, , , over] = near.step;
    let l2 = lon2;
    if (l2 - lon1 > 180) l2 -= 360; else if (lon1 - l2 > 180) l2 += 360;
    const az = azimuth(lon1, lat1, lon2, lat2);
    const code = this.classes[k];
    const across = code === 'SUB' && over === 1 ? az - 90 : az + 90;      // toward the overriding plate
    const mlon = (lon1 + l2) / 2, mlat = (lat1 + lat2) / 2;
    // The centre 400 km beyond the boundary, so that the section runs from 400 km before it to 1,200 km after.
    const centre = offsetAlong(mlon, mlat, across, 400);
    const name = this.meta.classes[k].name;
    const rel = relativeTo(centre[0], centre[1], (across + 360) % 360, lon, lat);
    return { lon: centre[0], lat: centre[1], az: (across + 360) % 360, half: 800, mark: rel.along,
      why: code === 'SUB' ? `across the ${name}, toward the overriding plate` : `across the nearest boundary, a ${name}` };
  }

  // Rows for the place card: [term, html, kind].
  rows(lon, lat) {
    if (!this.plates) return [];
    const p = this.plateAt(lon, lat);
    if (!p) return [];
    const rows = [['Tectonic plate', `${esc(p.name)} <span class="dir">${esc(morvelCode(p.code))}</span>`]];
    if (p.omega) {
      const v = velocityAt(p.omega, lon, lat), speed = Math.hypot(v.east, v.north);
      rows.push(['Plate moves', `${speed.toFixed(1)} mm/yr <span class="dir">toward ${compassPoint(flowToward(v.east, v.north))}</span>`]);
      const pole = p.pole;
      rows.push(['', `relative to the mean of all plates (no net rotation; NNR-MORVEL56), turning ${pole.rate.toFixed(3)}°/Myr anticlockwise about ${Math.abs(pole.lat).toFixed(1)}° ${pole.lat >= 0 ? 'N' : 'S'}, ${Math.abs(pole.lon).toFixed(1)}° ${pole.lon >= 0 ? 'E' : 'W'}`, 'gloss']);
    }
    const near = this.nearestStep(lon, lat);
    if (near && near.km < 1500) {
      const kind = this.meta.classes[near.step[4]];
      const a = this.across(near.step);
      rows.push(['Nearest plate boundary', near.km < 1 ? 'here' : `${near.km.toFixed(0)} km away`]);
      if (a) {
        const other = a.left === p ? a.right : a.right === p ? a.left : null;
        const pair = other ? `${esc(other.name)} and ${esc(p.name)}` : `${esc(a.left.name)} and ${esc(a.right.name)}`;
        const opening = a.divergence >= 0 ? `open at ${a.divergence.toFixed(1)} mm/yr` : `close at ${(-a.divergence).toFixed(1)} mm/yr`;
        const slide = Math.abs(a.rightLateral) >= 0.5 ? ` and slide ${Math.abs(a.rightLateral).toFixed(1)} mm/yr past each other, ${a.rightLateral > 0 ? 'right' : 'left'}-laterally` : '';
        let text = `${esc(kind.name.charAt(0).toUpperCase() + kind.name.slice(1))}: ${pair} ${opening}${slide}.`;
        if (kind.code === 'SUB') {
          const over = near.step[7] === 1 ? a.left : a.right, under = near.step[7] === 1 ? a.right : a.left;
          text += ` ${esc(under.name)} sinks beneath ${esc(over.name)}.`;
        }
        rows.push(['', text, 'gloss']);
      }
    }
    return rows;
  }

  // Names of the plates, larger ones first and smaller ones closer in.
  labelItems(zoom) {
    if (!this.shown || !this.labels) return [];
    const minArea = zoom < 2.6 ? 0.25 : zoom < 3.4 ? 0.05 : zoom < 4.4 ? 0.01 : 0;
    return this.labels.filter((l) => l.plate.area >= minArea).sort((a, b) => b.plate.area - a.plate.area)
      .map((l) => ({ key: `plate${l.plate.code}`, text: `${l.plate.name} plate`, lon: l.lon, lat: l.lat, className: 'plate',
        font: 'italic 600 12px "Source Serif 4", serif', tracking: 1.6, size: 14 }));
  }

  legend() {
    return { title: 'Plate boundaries', keys: GROUPS.map((g) => ({ label: g.label, colour: g.colour })),
      note: 'Bird (2003); teeth on the overriding plate; motion: NNR-MORVEL56 (Argus et al. 2011)' };
  }
}
