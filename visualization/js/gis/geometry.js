// The planar geometry of a city GIS, on the local frame's east and north (m): rings and their orientation,
// containment, polylines and the ribbons that draw them with a width, and extents. Pure arithmetic, checked in
// Node; the three.js scene builds its meshes from it.

// The signed area of a ring (shoelace formula): positive when counter-clockwise seen from above.
export function ringArea(ring) {
  let a = 0;
  for (let i = 0, n = ring.length; i < n; i++) {
    const [x0, y0] = ring[i], [x1, y1] = ring[(i + 1) % n];
    a += x0 * y1 - x1 * y0;
  }
  return a / 2;
}

// The ring with the orientation asked for (counter-clockwise outer rings, clockwise holes), and without a
// closing point repeating the first.
export function oriented(ring, ccw = true) {
  let r = ring.length > 1 && ring[0][0] === ring[ring.length - 1][0] && ring[0][1] === ring[ring.length - 1][1] ? ring.slice(0, -1) : ring.slice();
  if ((ringArea(r) > 0) !== ccw) r = r.reverse();
  return r;
}

// Whether a point lies inside a ring (even–odd rule, by crossings of a ray toward +x).
export function inRing([x, y], ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i], [xj, yj] = ring[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

export function inPolygon(p, outer, holes = []) {
  return inRing(p, outer) && !holes.some((h) => inRing(p, h));
}

// The centroid of a ring's area.
export function centroid(ring) {
  let a = 0, cx = 0, cy = 0;
  for (let i = 0, n = ring.length; i < n; i++) {
    const [x0, y0] = ring[i], [x1, y1] = ring[(i + 1) % n];
    const c = x0 * y1 - x1 * y0;
    a += c, cx += (x0 + x1) * c, cy += (y0 + y1) * c;
  }
  if (Math.abs(a) < 1e-12) {
    const n = ring.length;
    return [ring.reduce((s, p) => s + p[0], 0) / n, ring.reduce((s, p) => s + p[1], 0) / n];
  }
  return [cx / (3 * a), cy / (3 * a)];
}

export function extent(points) {
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const [x, y] of points) x0 = Math.min(x0, x), y0 = Math.min(y0, y), x1 = Math.max(x1, x), y1 = Math.max(y1, y);
  return { x0, y0, x1, y1 };
}

// The length of a polyline, and the point and direction at an arc length along it.
export function lineLength(line) {
  let s = 0;
  for (let i = 1; i < line.length; i++) s += Math.hypot(line[i][0] - line[i - 1][0], line[i][1] - line[i - 1][1]);
  return s;
}

export function pointAlong(line, s) {
  let left = s;
  for (let i = 1; i < line.length; i++) {
    const dx = line[i][0] - line[i - 1][0], dy = line[i][1] - line[i - 1][1], L = Math.hypot(dx, dy);
    if (left <= L || i === line.length - 1) {
      const t = L > 0 ? Math.min(1, Math.max(0, left / L)) : 0;
      return { x: line[i - 1][0] + t * dx, y: line[i - 1][1] + t * dy, dx: L > 0 ? dx / L : 1, dy: L > 0 ? dy / L : 0 };
    }
    left -= L;
  }
  return { x: line[0][0], y: line[0][1], dx: 1, dy: 0 };
}

// The point of a polyline nearest to p: its arc length s, distance, position, unit direction there, and side
// (+1 when p lies to the left of the line's direction, −1 to the right).
export function nearestOnLine(line, [px, py]) {
  let best = { s: 0, distance: Infinity, x: line[0][0], y: line[0][1], dx: 1, dy: 0, side: 1 };
  let s0 = 0;
  for (let i = 1; i < line.length; i++) {
    const [ax, ay] = line[i - 1], [bx, by] = line[i];
    const dx = bx - ax, dy = by - ay, L2 = dx * dx + dy * dy, L = Math.sqrt(L2);
    const t = L2 > 0 ? Math.min(1, Math.max(0, ((px - ax) * dx + (py - ay) * dy) / L2)) : 0;
    const x = ax + t * dx, y = ay + t * dy, d = Math.hypot(px - x, py - y);
    if (d < best.distance) {
      const cross = dx * (py - ay) - dy * (px - ax);
      best = { s: s0 + t * L, distance: d, x, y, dx: L > 0 ? dx / L : 1, dy: L > 0 ? dy / L : 0, side: cross >= 0 ? 1 : -1 };
    }
    s0 += L;
  }
  return best;
}

// The part of a polyline between arc lengths s0 and s1 (clamped to the line), as a polyline.
export function subLine(line, s0, s1) {
  const L = lineLength(line);
  const a = Math.max(0, Math.min(L, Math.min(s0, s1))), b = Math.max(0, Math.min(L, Math.max(s0, s1)));
  const out = [];
  const first = pointAlong(line, a);
  out.push([first.x, first.y]);
  let s = 0;
  for (let i = 1; i < line.length; i++) {
    s += Math.hypot(line[i][0] - line[i - 1][0], line[i][1] - line[i - 1][1]);
    if (s > a && s < b) out.push([line[i][0], line[i][1]]);
  }
  const last = pointAlong(line, b);
  out.push([last.x, last.y]);
  return out;
}

// A polyline resampled at n + 1 points equally spaced in arc length (its ends kept).
export function resample(line, n) {
  const L = lineLength(line), out = [];
  for (let k = 0; k <= n; k++) {
    const p = pointAlong(line, (k / n) * L);
    out.push([p.x, p.y]);
  }
  return out;
}

// The left and right edges of a ribbon of a width along a polyline, at its vertices: offsets along the
// bisector of the neighbouring segments' normals, mitred, the miter limited to 3 half-widths (as ribbon()).
export function edges(line, width) {
  const n = line.length, half = width / 2;
  const normals = [];
  for (let i = 1; i < n; i++) {
    const dx = line[i][0] - line[i - 1][0], dy = line[i][1] - line[i - 1][1], L = Math.hypot(dx, dy) || 1;
    normals.push([-dy / L, dx / L]);
  }
  const left = [], right = [];
  for (let i = 0; i < n; i++) {
    let o;
    if (i === 0) o = normals[0].map((v) => v * half);
    else if (i === n - 1) o = normals[n - 2].map((v) => v * half);
    else {
      const a = normals[i - 1], b = normals[i];
      const m = [a[0] + b[0], a[1] + b[1]], ml = Math.hypot(m[0], m[1]);
      if (ml < 1e-9) o = a.map((v) => v * half);
      else {
        const cos = (m[0] * a[0] + m[1] * a[1]) / ml;
        const len = Math.min(3, 1 / Math.max(cos, 1e-6)) * half;
        o = [(m[0] / ml) * len, (m[1] / ml) * len];
      }
    }
    left.push([line[i][0] + o[0], line[i][1] + o[1]]);
    right.push([line[i][0] - o[0], line[i][1] - o[1]]);
  }
  return { left, right };
}

// The triangles of a ribbon of a width along a polyline: two per segment, the corners mitred (the miter
// limited to 3 half-widths, beyond which the joint is left as the segments' own ends). Returns a flat list of
// (x, y) for the triangles' vertices, six numbers a triangle.
export function ribbon(line, width) {
  const out = [];
  const n = line.length;
  if (n < 2 || !(width > 0)) return out;
  const half = width / 2;
  const normals = [];
  for (let i = 1; i < n; i++) {
    const dx = line[i][0] - line[i - 1][0], dy = line[i][1] - line[i - 1][1], L = Math.hypot(dx, dy) || 1;
    normals.push([-dy / L, dx / L]);
  }
  const offsetAt = (i) => {
    if (i === 0) return normals[0].map((v) => v * half);
    if (i === n - 1) return normals[n - 2].map((v) => v * half);
    const a = normals[i - 1], b = normals[i];
    const m = [a[0] + b[0], a[1] + b[1]], ml = Math.hypot(m[0], m[1]);
    if (ml < 1e-9) return a.map((v) => v * half);
    const cos = (m[0] * a[0] + m[1] * a[1]) / ml;
    const len = Math.min(3, 1 / Math.max(cos, 1e-6)) * half;
    return [(m[0] / ml) * len, (m[1] / ml) * len];
  };
  const offs = line.map((_, i) => offsetAt(i));
  for (let i = 1; i < n; i++) {
    const [ax, ay] = line[i - 1], [bx, by] = line[i];
    const [oax, oay] = offs[i - 1], [obx, oby] = offs[i];
    const l0 = [ax + oax, ay + oay], r0 = [ax - oax, ay - oay], l1 = [bx + obx, by + oby], r1 = [bx - obx, by - oby];
    out.push(...r0, ...r1, ...l1, ...r0, ...l1, ...l0);
  }
  return out;
}
