// The layout of the foam, apart from its drawing so that it can be checked without a browser: the
// place of each being on the ground plane, and the height of each occasion.

import { colourIndex, COLOURS } from './palette.js';

export const HEIGHT = 40;  // the span of time, in scene units
export const RING = 13;    // the radius of the ring of beings

// A small deterministic generator, so that the layout of a canon is always the same.
function mulberry(seed) {
  return () => {
    seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// The place of each being on the ground plane: a force layout in which beings that meet or are tied
// attract, all repel, and each is drawn toward the sector of its discipline on the ring.
export function placeBeings(data) {
  const n = data.beings.length;
  const rand = mulberry(20261024);
  const weight = new Map();
  const add = (a, b, w) => { if (a === b) return; const k = a < b ? `${a},${b}` : `${b},${a}`; weight.set(k, (weight.get(k) || 0) + w); };
  for (const o of data.occasions) {
    const here = o.present.filter((p) => !p.derived).map((p) => p.being);
    for (let i = 0; i < here.length; i++) for (let j = i + 1; j < here.length; j++) add(here[i], here[j], 1 / Math.max(1, here.length - 1));
  }
  for (const t of data.ties) add(t.a, t.b, 1.5 + t.attested.length * 0.3);
  const sectors = COLOURS.length + 1;
  const home = data.beings.map((b) => {
    const angle = (2 * Math.PI * (colourIndex(b.colour) + 0.5)) / sectors;
    return { angle, x: Math.cos(angle) * RING * (0.6 + 0.4 * rand()) + rand() - 0.5, z: Math.sin(angle) * RING * (0.6 + 0.4 * rand()) + rand() - 0.5 };
  });
  const pairs = [...weight.entries()].map(([k, w]) => { const [a, b] = k.split(',').map(Number); return { a, b, w }; });
  for (let step = 0; step < 500; step++) {
    const cool = 1 - step / 500;
    const fx = new Float64Array(n), fz = new Float64Array(n);
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) {
        const dx = home[i].x - home[j].x, dz = home[i].z - home[j].z;
        const d2 = Math.max(0.25, dx * dx + dz * dz);
        const f = 28 / d2;
        fx[i] += f * dx; fz[i] += f * dz; fx[j] -= f * dx; fz[j] -= f * dz;
      }
    for (const { a, b, w } of pairs) {
      const dx = home[b].x - home[a].x, dz = home[b].z - home[a].z;
      const d = Math.hypot(dx, dz) || 1e-6;
      const f = 0.06 * Math.min(4, w) * (d - 3.2);
      fx[a] += (f * dx) / d; fz[a] += (f * dz) / d; fx[b] -= (f * dx) / d; fz[b] -= (f * dz) / d;
    }
    for (let i = 0; i < n; i++) {
      const tx = Math.cos(home[i].angle) * RING, tz = Math.sin(home[i].angle) * RING;
      fx[i] += 0.05 * (tx - home[i].x);
      fz[i] += 0.05 * (tz - home[i].z);
      const lim = 1.2 * cool + 0.05;
      home[i].x += Math.max(-lim, Math.min(lim, fx[i] * 0.4));
      home[i].z += Math.max(-lim, Math.min(lim, fz[i] * 0.4));
    }
  }
  // Scaled so that the farthest being stands on the ring.
  const far = Math.max(1e-6, ...home.map((h) => Math.hypot(h.x, h.z)));
  for (const h of home) { h.x *= RING / far; h.z *= RING / far; }
  return home;
}

// The height of each occasion: by calendar time, or by its rank among the occasions ("events").
export function heights(data, time = 'calendar') {
  const order = data.complex.vertices;
  const y = new Array(data.occasions.length).fill(0);
  if (time === 'events') {
    order.forEach((o, k) => { y[o] = (HEIGHT * k) / Math.max(1, order.length - 1); });
  } else {
    const y0 = data.occasions[order[0]].year, y1 = data.occasions[order[order.length - 1]].year;
    for (const o of order) y[o] = (HEIGHT * (data.occasions[o].year - y0)) / Math.max(1e-9, y1 - y0);
  }
  return y;
}

