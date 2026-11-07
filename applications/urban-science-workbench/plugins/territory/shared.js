// What the territory's layers share: the world's cells, Space2Stats' tiles, the colours of values on the scene, a compact
// legend. Colours come from the visualization layer's validated maps; a map that runs from pale to deep is turned over by
// night, so that much of something stands out against the plum in light and, by day, against the lilac in depth.

import { colourMap } from '/visualization/styles/index.js';

export function tone(name, theme, { from = 0.1, to = 1, keep = false } = {}) {
  const m = colourMap(name).colour;
  const darkLow = name === 'plum';          // plum runs from its night upward already
  const flip = !keep && ((theme === 'night') !== darkLow);
  return (u) => {
    const t = Math.min(1, Math.max(0, u));
    return m(flip ? 1 - (from + (to - from) * t) : from + (to - from) * t);
  };
}
const hex = (c) => `rgb(${c.map((v) => Math.round(v * 255)).join(',')})`;
export const gradient = (colour, n = 12) => `linear-gradient(90deg, ${Array.from({ length: n }, (_, i) => hex(colour(i / (n - 1)))).join(',')})`;

export function card(h, title, { colour, low, high, note, extra = [], source }) {
  return h('div', { class: 'card' }, h('b', {}, title), ...extra,
    colour ? h('span', { class: 'ramp', style: `background:${gradient(colour)}` }) : null,
    low !== undefined ? h('div', { class: 'ramp-ends' }, h('span', {}, low), h('span', {}, high)) : null,
    note ? h('div', { class: 'ld-quiet' }, note) : null,
    source ? h('div', { class: 'ld-quiet' }, source) : null);
}

export const compact = (v) => {
  if (!Number.isFinite(v)) return '—';
  const a = Math.abs(v);
  if (a >= 1e9) return `${(v / 1e9).toFixed(a >= 1e10 ? 0 : 1)} bn`;
  if (a >= 1e6) return `${(v / 1e6).toFixed(a >= 1e7 ? 0 : 1)} m`;
  if (a >= 1e4) return `${Math.round(v / 1e3)}k`;
  if (a >= 100) return Math.round(v).toLocaleString('en');
  if (a >= 1) return v.toFixed(1);
  return v.toFixed(2);
};

// The world's cells at a resolution: kept for the page's life.
const worlds = new Map();
export function world(ctx, res) {
  if (!worlds.has(res)) {
    worlds.set(res, ctx.ask('earth_people', { res }).then((w) => {
      w.centres = w.rings.map((r) => r.reduce((a, p) => [a[0] + p[0] / r.length, a[1] + p[1] / r.length], [0, 0]));
      return w;
    }).catch((e) => { worlds.delete(res); throw e; }));
  }
  return worlds.get(res);
}

// The cell nearest a point, among many, through a coarse index on the sphere.
export function nearest(centres, step = 4) {
  const buckets = new Map();
  const key = (lon, lat) => `${Math.floor((((lon + 540) % 360)) / step)},${Math.floor((lat + 90) / step)}`;
  centres.forEach((c, k) => { const q = key(c[0], c[1]); if (!buckets.has(q)) buckets.set(q, []); buckets.get(q).push(k); });
  return (lon, lat, limit = Infinity) => {
    let best = -1, d = Infinity;
    const x = Math.floor(((lon + 540) % 360) / step), y = Math.floor((lat + 90) / step), span = Math.ceil(360 / step);
    for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) {
      for (const k of buckets.get(`${(x + i + span) % span},${y + j}`) || []) {
        const c = centres[k], dx = (((c[0] - lon + 540) % 360) - 180) * Math.cos((lat * Math.PI) / 180), dy = c[1] - lat, q = dx * dx + dy * dy;
        if (q < d) { d = q; best = k; }
      }
    }
    return d <= limit * limit ? best : -1;
  };
}

// The tiles of Space2Stats about the camera: 4 x 4 degrees each, asked once and kept.
const tiles = new Map();
export function tileAt(lat, lon) {
  const s = Math.floor(lat / 4) * 4, w = Math.floor(lon / 4) * 4, key = `${s},${w}`;
  if (!tiles.has(key)) {
    tiles.set(key, fetch(`/world/space2stats?lat=${s + 2}&lon=${w + 2}`).then(async (r) => {
      const t = await r.json();
      if (!r.ok) throw new Error(t.message || `the service answered ${r.status}`);
      t.centres = t.rings.map((q) => q.reduce((a, p) => [a[0] + p[0] / q.length, a[1] + p[1] / q.length], [0, 0]));
      return t;
    }).catch((e) => { setTimeout(() => tiles.delete(key), 90000); throw e; }));   // a failing service is asked again after a while, not at every frame
  }
  return tiles.get(key);
}
export function tilesInView(view, limit = 12) {
  const DEG = 180 / Math.PI, lat = view.lat * DEG, lon = view.lon * DEG;
  const span = Math.max(1.5, (view.range * 1.15) / 111000);
  const slat = Math.min(span, 18), slon = Math.min(span / Math.max(0.25, Math.cos(view.lat)), 26);
  const out = [];
  for (let la = Math.floor((lat - slat) / 4) * 4; la <= lat + slat; la += 4) for (let lo = Math.floor((lon - slon) / 4) * 4; lo <= lon + slon; lo += 4) {
    if (la < -60 || la > 80) continue;
    out.push([la, lo, Math.hypot(la + 2 - lat, ((lo + 2 - lon + 540) % 360) - 180)]);
  }
  return out.sort((a, b) => a[2] - b[2]).slice(0, limit);
}
