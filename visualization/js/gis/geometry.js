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

// The triangles of a polygon with holes, by ear clipping, each hole first bridged to the outer ring at a point it can see (Eberly,
// "Triangulation by ear clipping", 2002). Returns indices into [...outer, ...holes.flat()], three a triangle; their areas sum to
// the polygon's. Rings may wind either way and need not be closed; degenerate input gives fewer triangles, never a failure.
export function triangulate(outer, holes = []) {
  const pts = [...outer];
  const indices = (ring, start, ccw) => {
    const ids = ring.map((_, i) => start + i);
    return (ringArea(ring) > 0) === ccw ? ids : ids.reverse();
  };
  let loop = indices(outer, 0, true);
  const inner = holes.map((h) => { const start = pts.length; pts.push(...h); return indices(h, start, false); })
    .map((ids) => { let m = 0; ids.forEach((id, k) => { if (pts[id][0] > pts[ids[m]][0]) m = k; }); return { ids, m }; })
    .sort((a, b) => pts[b.ids[b.m]][0] - pts[a.ids[a.m]][0]);
  const cross = (a, b, c) => (pts[b][0] - pts[a][0]) * (pts[c][1] - pts[a][1]) - (pts[b][1] - pts[a][1]) * (pts[c][0] - pts[a][0]);
  const inside = (p, a, b, c) => {
    const [x, y] = pts[p], [ax, ay] = pts[a], [bx, by] = pts[b], [cx, cy] = pts[c];
    if ((x === ax && y === ay) || (x === bx && y === by) || (x === cx && y === cy)) return false;
    const d1 = (bx - ax) * (y - ay) - (by - ay) * (x - ax), d2 = (cx - bx) * (y - by) - (cy - by) * (x - bx), d3 = (ax - cx) * (y - cy) - (ay - cy) * (x - cx);
    return d1 >= 0 && d2 >= 0 && d3 >= 0;
  };
  for (const { ids, m } of inner) {
    const M = pts[ids[m]];
    let best = Infinity, at = -1;
    for (let i = 0; i < loop.length; i++) {
      const a = pts[loop[i]], b = pts[loop[(i + 1) % loop.length]];
      if ((a[1] > M[1]) === (b[1] > M[1])) continue;
      const x = a[0] + ((M[1] - a[1]) / (b[1] - a[1])) * (b[0] - a[0]);
      if (x >= M[0] && x < best) { best = x; at = a[0] > b[0] ? i : (i + 1) % loop.length; }
    }
    if (at < 0) {   // the hole lies outside: join it at the nearest point, as bad input allows
      let near = Infinity;
      loop.forEach((id, i) => { const d = Math.hypot(pts[id][0] - M[0], pts[id][1] - M[1]); if (d < near) { near = d; at = i; } });
    } else {
      // A reflex point of the ring may lie between: the one nearest in angle from the hole's point is visible.
      const I = pts.length;
      pts.push([best, M[1]]);
      const P = loop[at];
      let angle = Infinity, length = Infinity;
      loop.forEach((id, i) => {
        if (id === P) return;
        const prev = loop[(i + loop.length - 1) % loop.length], next = loop[(i + 1) % loop.length];
        if (cross(prev, id, next) >= 0) return;
        const up = pts[P][1] > M[1];
        const [a, b, c] = up ? [ids[m], I, P] : [ids[m], P, I];
        if (!inside(id, a, b, c)) return;
        const dx = pts[id][0] - M[0], dy = pts[id][1] - M[1], tilt = Math.abs(dy) / (Math.hypot(dx, dy) || 1), far = Math.hypot(dx, dy);
        if (tilt < angle || (tilt === angle && far < length)) { angle = tilt; length = far; at = i; }
      });
      pts.pop();
    }
    const cycle = [...ids.slice(m), ...ids.slice(0, m), ids[m]];
    loop = [...loop.slice(0, at + 1), ...cycle, ...loop.slice(at)];
  }
  // Clip ears from a linked loop.
  const out = [];
  const n = loop.length;
  if (n < 3) return out;
  const prev = loop.map((_, i) => (i + n - 1) % n), next = loop.map((_, i) => (i + 1) % n);
  let left = n, i = 0, tried = 0;
  while (left > 3 && tried <= left) {
    const a = prev[i], c = next[i];
    const turn = cross(loop[a], loop[i], loop[c]);
    let ear = turn > 0;
    if (ear) for (let k = next[c]; k !== a; k = next[k]) {
      if (inside(loop[k], loop[a], loop[i], loop[c])) { ear = false; break; }
    }
    if (ear || turn === 0) {
      if (ear) out.push(loop[a], loop[i], loop[c]);
      next[a] = c, prev[c] = a;
      left--;
      i = c;
      tried = 0;
    } else { i = c; tried++; }
  }
  if (left === 3) { if (cross(loop[prev[i]], loop[i], loop[next[i]]) > 0) out.push(loop[prev[i]], loop[i], loop[next[i]]); return out; }
  // What self-crossing input leaves: a fan of what turns the right way.
  for (let k = next[i]; next[k] !== i; k = next[k]) if (cross(loop[i], loop[k], loop[next[k]]) > 0) out.push(loop[i], loop[k], loop[next[k]]);
  return out;
}
