// The district of an area, shared by the plug-ins that work on one: its box (the box its entry in areas.json names,
// else the middle of its first view), the walking network cut from the city model in that box, and the footprints
// of the buildings standing in it.
//
// The network: the city model's ways where people walk (streets, lanes, malls; not motorways, not the separate
// footway lines drawn beside streets, which would double every street), cut to the box, split at every point two
// ways share and at their ends: the junctions; each piece between two junctions a segment. A way that closes on
// itself is split once more so that no segment is a loop. Ways that cross without a shared point are left to the
// engines, which drop the shorter of each crossing pair and name it.
//
// A building stands in the box when its footprint's centre does; it is sent as its key (the city model's id), name,
// ring (local metres), base and top.

const WALK = new Set(['primary', 'secondary', 'tertiary', 'residential', 'unclassified', 'living_street', 'pedestrian', 'service', 'primary_link', 'secondary_link', 'tertiary_link']);

// The district of an area: the box its entry in areas.json names (`district`), else the middle of its first view.
export function districtOf(area) {
  if (area.district) return area.district;
  const [s, w, n, e] = area.view.box, cy = (s + n) / 2, cx = (w + e) / 2, dy = 0.0028, dx = dy / Math.cos((cy * Math.PI) / 180);
  return { south: cy - dy, west: cx - dx, north: cy + dy, east: cx + dx };
}

export function boxOf(frame, b) {
  const [x0, y0] = frame.toLocal(b.south, b.west), [x1, y1] = frame.toLocal(b.north, b.east);
  return { x0: Math.min(x0, x1), y0: Math.min(y0, y1), x1: Math.max(x0, x1), y1: Math.max(y0, y1), geo: b };
}
export const inside = (p, B) => p[0] >= B.x0 && p[0] <= B.x1 && p[1] >= B.y0 && p[1] <= B.y1;

export function networkOf(model, B) {
  const runs = [];
  for (const r of model.roads || []) {
    if (!WALK.has(r.kind) || !r.line || r.line.length < 2) continue;
    let run = [];
    for (const p of r.line) {
      if (inside(p, B)) run.push(p);
      else { if (run.length >= 2) runs.push({ line: run, kind: r.kind, name: r.name || '' }); run = []; }
    }
    if (run.length >= 2) runs.push({ line: run, kind: r.kind, name: r.name || '' });
  }
  const key = (p) => `${Math.round(p[0] * 100)},${Math.round(p[1] * 100)}`;
  const uses = new Map();
  for (const r of runs) for (const p of new Set(r.line.map(key))) uses.set(p, (uses.get(p) || 0) + 1);
  const junctions = [], index = new Map();
  const junction = (p) => { const k = key(p); if (!index.has(k)) { index.set(k, junctions.length); junctions.push([+p[0].toFixed(2), +p[1].toFixed(2)]); } return index.get(k); };
  const segments = [], seen = new Set();
  const push = (line, kind, name) => {
    if (line.length < 2) return;
    let u = junction(line[0]), v = junction(line[line.length - 1]);
    if (u === v) {
      // A closed piece: split it at its middle point so that no segment is a loop.
      if (line.length < 3) return;
      const m = Math.floor(line.length / 2);
      push(line.slice(0, m + 1), kind, name);
      push(line.slice(m), kind, name);
      return;
    }
    const sig = `${Math.min(u, v)}-${Math.max(u, v)}-${line.length}-${key(line[Math.floor(line.length / 2)])}`;
    if (seen.has(sig)) return;
    seen.add(sig);
    const kindOf = kind === 'service' ? 'laneway' : kind === 'pedestrian' || kind === 'living_street' ? 'mall' : 'street';
    segments.push({ u, v, line: line.map((p) => [+p[0].toFixed(2), +p[1].toFixed(2)]), kind: kindOf, name });
  };
  for (const r of runs) {
    let start = 0;
    for (let k = 1; k < r.line.length; k++) {
      const last = k === r.line.length - 1;
      if (last || uses.get(key(r.line[k])) > 1) { push(r.line.slice(start, k + 1), r.kind, r.name); start = k; }
    }
  }
  // Drop junctions no segment uses (renumbering).
  const used = new Set(segments.flatMap((s) => [s.u, s.v]));
  const renumber = new Map();
  const J = [];
  junctions.forEach((p, i) => { if (used.has(i)) { renumber.set(i, J.length); J.push(p); } });
  for (const s of segments) s.u = renumber.get(s.u), s.v = renumber.get(s.v);
  return { junctions: J, segments };
}

export function footprintsIn(model, B) {
  const out = [];
  for (const b of (model && model.buildings) || []) {
    const f = b.footprint;
    if (!f || f.length < 3) continue;
    const cx = f.reduce((a, p) => a + p[0], 0) / f.length, cy = f.reduce((a, p) => a + p[1], 0) / f.length;
    if (cx >= B.x0 && cx <= B.x1 && cy >= B.y0 && cy <= B.y1) out.push({ key: b.id, name: b.name || '', ring: f, base: b.min_height_m || 0, top: b.height_m || 10 });
  }
  return out;
}
