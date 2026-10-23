// The loom's lanes and its scale of time. Pure: the application's checks run these in Node.

import { COLOURS, colourIndex } from './palette.js';
import { roundYears } from './format.js';

// One lane per being, grouped by the being's discipline in the palette's fixed order, then by its
// first occasion, then by name; the beings of hidden disciplines have no lane. Returns the lanes in
// order with their heights, the groups with their extents, and the whole height.
export function laneLayout(data, hidden = new Set(), { lane = 19, gap = 22, top = 0 } = {}) {
  const firstOf = (b) => {
    const os = data.beings[b].occasions || [];
    return os.length ? Math.min(...os.map((o) => data.occasions[o].first)) : Infinity;
  };
  const beings = data.beings.map((_, b) => b).filter((b) => !hidden.has(data.beings[b].colour || 'unknown'));
  beings.sort((a, b) => colourIndex(data.beings[a].colour) - colourIndex(data.beings[b].colour)
    || firstOf(a) - firstOf(b) || String(data.beings[a].name).localeCompare(String(data.beings[b].name)));
  // Within each discipline, beings move toward the mean place of those they are tied to
  // (barycentric ordering, a few sweeps), so that most lenses span few lanes; the groups keep
  // the palette's order, so that every legend reads the same.
  const partners = new Map(beings.map((b) => [b, []]));
  for (const t of data.ties || []) {
    const w = 1 + (t.attested?.length || 0);
    if (partners.has(t.a) && partners.has(t.b)) { partners.get(t.a).push([t.b, w]); partners.get(t.b).push([t.a, w]); }
  }
  for (let sweep = 0; sweep < 6; sweep++) {
    const place = new Map(beings.map((b, i) => [b, i]));
    const centre = (b) => {
      const ps = partners.get(b);
      if (!ps.length) return place.get(b);
      const w = ps.reduce((s, [, x]) => s + x, 0);
      return (ps.reduce((s, [q, x]) => s + place.get(q) * x, 0) + place.get(b)) / (w + 1);
    };
    const keys = new Map(beings.map((b) => [b, centre(b)]));
    beings.sort((a, b) => colourIndex(data.beings[a].colour) - colourIndex(data.beings[b].colour) || keys.get(a) - keys.get(b) || a - b);
  }
  const y = new Map(), groups = [];
  let at = top, last = null;
  for (const b of beings) {
    const c = data.beings[b].colour || 'unknown';
    if (c !== last) {
      if (last !== null) at += gap - lane;
      at += gap * 0.9;
      groups.push({ colour: c, top: at - lane / 2, bottom: at - lane / 2, first: b });
      last = c;
    }
    y.set(b, at);
    groups[groups.length - 1].bottom = at + lane / 2;
    at += lane;
  }
  return { order: beings, y, groups, height: at };
}

// Ticks of years for a span of time drawn over a width: a step of 1, 2, 5, 10, 20, 25, 50 or 100
// years, the smallest that leaves at least `spacing` pixels between labels.
export function yearTicks(t0, t1, width, spacing = 70) {
  const span = Math.max(t1 - t0, 1e-9);
  const step = [1, 2, 5, 10, 20, 25, 50, 100, 200, 500].find((s) => (s / span) * width >= spacing) || 1000;
  return { step, ticks: roundYears(t0, t1, step) };
}

// The span of a canon on the line of time: its occasions' first and last days, with a margin.
export function timeSpan(data, margin = 0.6) {
  const ys = data.occasions.flatMap((o) => [o.year, o.yearLast ?? o.year]);
  if (!ys.length) return [0, 1];
  return [Math.min(...ys) - margin, Math.max(...ys) + margin];
}

export { COLOURS };
