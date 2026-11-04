// A map of the sphere in a projection, in the plane of an SVG drawing: the geography domain's
// projections (Equal Earth, equirectangular), x east and y down, `scale` drawing units for the
// projection's unit. Pure arithmetic, checked in Node.
//
// The raw projection is used, not the domain's `forward`, which wraps +180° to −180°: Natural Earth
// splits its polygons at the antimeridian, and a ring's vertices on +180° must stay on the right edge.
// A great-circle arc is interpolated on the sphere (slerp), lifted a little so that the arcs between
// the same places part, and split where it crosses the map's edge.

import { PROJECTIONS, extent, outline as outlineOf, inverse } from '../../../domains/geography/js/index.js';

const DEG = Math.PI / 180;
const f1 = (v) => v.toFixed(1);

export function planar(id = 'equal-earth', scale = 300) {
  const p = PROJECTIONS[id];
  if (!p) throw new Error(`no projection ${id}`);
  const [xm, ym] = extent(p);
  const W = xm * scale, H = ym * scale;
  const project = (lon, lat) => {
    const [x, y] = p.raw(Math.max(-180, Math.min(180, lon)) * DEG, Math.max(-90, Math.min(90, lat)) * DEG);
    return [x * scale, -y * scale];
  };
  const ring = (r) => r.map(([lon, lat], i) => { const [x, y] = project(lon, lat); return `${i ? 'L' : 'M'}${f1(x)},${f1(y)}`; }).join('') + 'Z';
  return {
    id, projection: p, scale, W, H,
    name: p.name || id,
    project,
    // The place at a point of the drawing, or null off the map.
    unproject(x, y) { return inverse(p, x / scale, -y / scale); },
    viewBox(pad = 8) { return `${f1(-W - pad)} ${f1(-H - pad)} ${f1(2 * W + 2 * pad)} ${f1(2 * H + 2 * pad)}`; },
    // The outline of the sphere on the map, closed.
    outline(n = 160) { return outlineOf(p, n).map(([x, y], i) => `${i ? 'L' : 'M'}${f1(x * scale)},${f1(-y * scale)}`).join('') + 'Z'; },
    // Meridians and parallels every `step` degrees (parallels short of the poles).
    graticule(step = 30, every = 3) {
      let d = '';
      for (let lon = -180; lon <= 180; lon += step)
        for (let lat = -90, i = 0; lat <= 90; lat += every, i++) { const [x, y] = project(lon, lat); d += `${i ? 'L' : 'M'}${f1(x)},${f1(y)}`; }
      for (let lat = -90 + step; lat < 90; lat += step)
        for (let lon = -180, i = 0; lon <= 180; lon += every, i++) { const [x, y] = project(lon, lat); d += `${i ? 'L' : 'M'}${f1(x)},${f1(y)}`; }
      return d;
    },
    ring,
    // A GeoJSON Polygon or MultiPolygon as a path (holes kept: draw with the even-odd or non-zero rule).
    feature(geometry) {
      const polygons = geometry?.type === 'Polygon' ? [geometry.coordinates] : geometry?.type === 'MultiPolygon' ? geometry.coordinates : [];
      return polygons.map((poly) => poly.map(ring).join('')).join('');
    },
    // A polyline of [lon, lat] (or a flat [lon, lat, lon, lat, …]) as a path.
    line(coords) {
      const pts = Array.isArray(coords[0]) ? coords : Array.from({ length: Math.floor(coords.length / 2) }, (_, i) => [coords[2 * i], coords[2 * i + 1]]);
      return pts.map(([lon, lat], i) => { const [x, y] = project(lon, lat); return `${i ? 'L' : 'M'}${f1(x)},${f1(y)}`; }).join('');
    },
    // The great-circle arc from a to b ({latitude, longitude} in degrees), lifted by `bulge` of its
    // angle, split where it crosses the map's edge; '' for two places at one point.
    greatCircle(a, b, bulge = 0.12) {
      const vec = (lat, lon) => [Math.cos(lat * DEG) * Math.cos(lon * DEG), Math.cos(lat * DEG) * Math.sin(lon * DEG), Math.sin(lat * DEG)];
      const u = vec(a.latitude, a.longitude), v = vec(b.latitude, b.longitude);
      const omega = Math.acos(Math.max(-1, Math.min(1, u[0] * v[0] + u[1] * v[1] + u[2] * v[2])));
      if (omega < 1e-6) return '';
      const n = Math.max(8, Math.ceil(omega * 40));
      let d = '', prev = null;
      for (let i = 0; i <= n; i++) {
        const t = i / n;
        const s0 = Math.sin((1 - t) * omega) / Math.sin(omega), s1 = Math.sin(t * omega) / Math.sin(omega);
        const q = [s0 * u[0] + s1 * v[0], s0 * u[1] + s1 * v[1], s0 * u[2] + s1 * v[2]];
        const lat = Math.asin(Math.max(-1, Math.min(1, q[2]))) / DEG, lon = Math.atan2(q[1], q[0]) / DEG;
        let [x, y] = project(lon, lat);
        y -= bulge * scale * Math.sin(Math.PI * t) * omega * 0.35;
        const jump = prev && Math.abs(x - prev[0]) > W;
        d += `${!prev || jump ? 'M' : 'L'}${f1(x)},${f1(y)}`;
        prev = [x, y];
      }
      return d;
    },
  };
}
