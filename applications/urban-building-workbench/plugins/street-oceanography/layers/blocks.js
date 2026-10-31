// The street a sensor stands on, for drawing what it measures there: the nearest mapped way (a named one if
// one lies within 35 m), the part of it within `half` metres of the sensor's foot, and the half of the
// street on the sensor's side where its name gives one ("Town Hall (West)"), with the unit direction of the
// way at the foot and, for each of its counts, the direction along the way it runs (`toward`, unit vectors)
// and its compass name (`directions`).

const SIDES = ['North', 'South', 'East', 'West'];
const COMPASS = { North: [0, 1], South: [0, -1], East: [1, 0], West: [-1, 0] };

export function sideOf(name) {
  const m = /\((North|South|East|West)\)/.exec(name || '');
  return m ? m[1] : null;
}

export function sensorBlocks(ctx, sensors, { half = 60 } = {}) {
  const { G, frame, model } = ctx;
  const roads = (model && model.roads || []).filter((r) => r.line && r.line.length >= 2 && !['steps', 'corridor', 'service'].includes(r.kind));
  const out = [];
  for (const s of sensors) {
    const p = frame.toLocal(s.lat, s.lon).slice(0, 2);
    let best = null, named = null;
    for (const r of roads) {
      const n = G.nearestOnLine(r.line, p);
      if (!best || n.distance < best.n.distance) best = { r, n };
      if (r.name && n.distance <= 35 && (!named || n.distance < named.n.distance)) named = { r, n };
    }
    const pick = named || best;
    if (!pick || pick.n.distance > 60) continue;
    const { r, n } = pick;
    const width = Math.max(4, r.width_m || 8);
    let line = G.subLine(r.line, n.s - half, n.s + half);
    const side = sideOf(s.name);
    // A sensor on one footpath: the half of the street on its side.
    // Drawn at least 8 m wide (5 m for one footpath), so that the water reads at the city's scale.
    let drawn = Math.max(8, width);
    if (side) {
      const e = G.edges(line, width / 2);
      line = n.side > 0 ? e.left : e.right;
      drawn = Math.max(5, width / 2);
    }
    // The way's direction at the foot, and toward which compass point each of the sensor's counts runs.
    const axis = [n.dx, n.dy];
    const names = s.directions || s.toward || [];
    const toward = names.map((name) => {
      const c = COMPASS[name];
      if (!c) return null;
      const dot = axis[0] * c[0] + axis[1] * c[1];
      return dot >= 0 ? axis : [-axis[0], -axis[1]];
    });
    out.push({ ...s, xy: p, line, width: drawn, street: r.name || r.kind, side, axis, directions: names, toward, foot: [n.x, n.y] });
  }
  return out;
}

export { SIDES };
