// A street as the crowd model takes it: its axis from the street's record of holdings, its cross streets where
// the city model's named ways cross the axis (crossings within 25 m of one another are one), a signal at
// each street and none at a lane or a mall, and the blocks between; and the walkers entering at its start at the
// rate the sensors there count toward its far end.

const NOT_CROSSING = new Set(['footway', 'path', 'steps', 'corridor', 'cycleway', 'busway']);
const UNSIGNALLED = /\b(Lane|Place|Alley|Arcade|Walk|Mall|Way)\b/;

// The arc lengths along `axis` at which polyline `line` crosses it.
function crossings(axis, line) {
  const out = [];
  let s0 = 0;
  for (let i = 1; i < axis.length; i++) {
    const [ax, ay] = axis[i - 1], [bx, by] = axis[i];
    const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy);
    for (let j = 1; j < line.length; j++) {
      const [cx, cy] = line[j - 1], [ex, ey] = line[j];
      const fx = ex - cx, fy = ey - cy, den = dx * fy - dy * fx;
      if (Math.abs(den) < 1e-9) continue;
      const t = ((cx - ax) * fy - (cy - ay) * fx) / den, u = ((cx - ax) * dy - (cy - ay) * dx) / den;
      if (t >= 0 && t <= 1 && u >= 0 && u <= 1) out.push(s0 + t * L);
    }
    s0 += L;
  }
  return out;
}

export function streetModel(ctx, street, { crossing = 20 } = {}) {
  const { G, frame, model } = ctx;
  const axis = street.axis.map(([lat, lon]) => frame.toLocal(lat, lon).slice(0, 2));
  const length = G.lineLength(axis);
  const hits = [];
  for (const r of (model && model.roads) || []) {
    if (!r.name || r.name === street.name || NOT_CROSSING.has(r.kind) || !r.line) continue;
    for (const s of crossings(axis, r.line)) hits.push({ s, name: r.name });
  }
  hits.sort((a, b) => a.s - b.s);
  const groups = [];
  for (const x of hits) {
    const g = groups[groups.length - 1];
    if (g && x.s - g.last < 25) { g.items.push(x); g.last = x.s; } else groups.push({ items: [x], last: x.s });
  }
  let cross = groups.map((g) => {
    const names = {};
    for (const x of g.items) names[x.name] = (names[x.name] || 0) + 1;
    const name = Object.entries(names).sort((a, b) => b[1] - a[1])[0][0];
    return { s: g.items.reduce((a, x) => a + x.s, 0) / g.items.length, name, signalled: !UNSIGNALLED.test(name) };
  }).filter((c) => c.s > 30 && c.s < length + 10);
  // The street ends at a crossing: its last cross street, or its end, unsignalled.
  if (!cross.length || cross[cross.length - 1].s < length - 30) cross.push({ s: length, name: 'the street\'s end', signalled: false });
  const blockLengths = [], signalled = [];
  let from = 0;
  for (const c of cross) {
    blockLengths.push(Math.max(10, c.s - crossing / 2 - from));
    signalled.push(c.signalled ? 1 : 0);
    from = c.s + crossing / 2;
  }
  return { axis, length, crossings: cross, blockLengths, signalled, crossing, name: street.name };
}

// The walkers entering at the street's start each second, hour by hour of a weekday (0 = Monday): the mean
// counts toward the axis's direction at the sensors on the street within `within` metres of its start.
export function entering(ctx, model, sensors, weekday, { within = 160 } = {}) {
  const { G } = ctx;
  const start = model.axis[0], dir = G.pointAlong(model.axis, 0);
  const used = [];
  const demand = new Array(24).fill(0);
  for (const f of sensors) {
    if (f.street !== model.name || !f.toward || f.toward.length !== 2) continue;
    const at = G.nearestOnLine(model.axis, f.foot);
    if (at.s > within || at.distance > 40 || Math.hypot(f.foot[0] - start[0], f.foot[1] - start[1]) > within + 40) continue;
    const w = f.toward[0][0] * dir.dx + f.toward[0][1] * dir.dy > 0 ? 0 : 1;
    const week = w ? f.way2 : f.way1;
    for (let hr = 0; hr < 24; hr++) demand[hr] += week[24 * weekday + hr] / 3600;
    used.push({ id: f.id, name: f.name, toward: f.directions[w] });
  }
  return { demand, used };
}
