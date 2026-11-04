// The view of a map: a scale k and a translation (x, y) of the drawing, applied as
// translate(x, y) scale(k). Zooming about a point keeps that point where it was; fitting frames a set of
// points. Pure arithmetic, checked in Node.

export const identity = Object.freeze({ k: 1, x: 0, y: 0 });

export const transform = (v) => `translate(${v.x.toFixed(2)},${v.y.toFixed(2)}) scale(${v.k.toFixed(4)})`;

// The view at scale k, the drawing's point under (px, py) kept fixed; at k = 1, the whole map.
export function zoomAt(v, k, px, py, { min = 1, max = 18 } = {}) {
  const k2 = Math.max(min, Math.min(max, k));
  if (k2 === 1 && min === 1) return { ...identity };
  return { k: k2, x: px - (px - v.x) * (k2 / v.k), y: py - (py - v.y) * (k2 / v.k) };
}

// A point of the screen's drawing coordinates in the map's own plane.
export const toPlane = (v, x, y) => [(x - v.x) / v.k, (y - v.y) / v.k];

// The view framing points [[x, y], …] of the plane on a map of half-size (W, H), with a margin; the
// whole map when they span most of it.
export function fit(points, W, H, { margin = 0.8, max = 18, least = [40, 30] } = {}) {
  if (!points.length) return { ...identity };
  const xs = points.map((p) => p[0]), ys = points.map((p) => p[1]);
  const bw = Math.max(...xs) - Math.min(...xs), bh = Math.max(...ys) - Math.min(...ys);
  if (bw > 1.2 * W || bh > 1.2 * H) return { ...identity };
  const k = Math.max(1, Math.min(max, (margin * 2 * W) / Math.max(bw, least[0]), (margin * 2 * H) / Math.max(bh, least[1])));
  const cx = (Math.max(...xs) + Math.min(...xs)) / 2, cy = (Math.max(...ys) + Math.min(...ys)) / 2;
  return { k, x: -cx * k, y: -cy * k };
}

// The view centred on a point of the plane at the present scale.
export const centreOn = (v, x, y) => ({ k: v.k, x: -x * v.k, y: -y * v.k });
