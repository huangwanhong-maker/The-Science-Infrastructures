// An answer carried to the map: values by place, and by time where a table has it, made a record of the workbench's line of time
// (ADR 0079), so that the card, the hover and the moment work as for any other record. Pure: the world's mosaic is handed in.

// The mean of a ring's corners: near enough for finding neighbours and for flying there.
export const centreOf = (ring) => ring.reduce((a, p) => [a[0] + p[0] / ring.length, a[1] + p[1] / ring.length], [0, 0]);

const COUNTER_RADIUS = 46;   // metres: a point of a table, as a small hexagon

export function hexAbout(lon, lat, radius = COUNTER_RADIUS) {
  const k = 111320, c = Math.cos((lat * Math.PI) / 180);
  return Array.from({ length: 6 }, (_, i) => {
    const a = (i * Math.PI) / 3;
    return [lon + (radius * Math.cos(a)) / (k * c), lat + (radius * Math.sin(a)) / k];
  });
}

// A label of a table's time column, as milliseconds since 1970 in UTC.
export function moment(label, kind) {
  if (kind === 'ms') return Number(label);
  if (kind === 'year') return Date.UTC(Number(label), 0, 1);
  if (kind === 'month') { const [y, m] = String(label).split('-').map(Number); return Date.UTC(y, m - 1, 15); }
  if (kind === 'date') return Date.parse(`${label}T00:00:00Z`);
  const n = Number(label);
  return Number.isFinite(n) ? (n > 1e9 ? n : Date.UTC(n, 0, 1)) : Date.parse(label);
}

// What can be drawn of an answer's values: a list of { title, scale, unit, values: [{ id, time?, value }] }.
export function drawables(answerBlocks, table) {
  const out = [];
  const P = table?.places;
  for (const b of answerBlocks) {
    if (b.kind === 'map') {
      out.push({ title: b.title, scale: b.scale, unit: b.unit, economies: true, values: Object.entries(b.values || {}).map(([id, value]) => ({ id, value })) });
    } else if (b.kind === 'rows' && P) {
      const ids = P.kind === 'row' ? null : table.data[P.column];
      const times = P.time ? table.data[P.time] : null;
      out.push({
        title: b.title, scale: b.scale, unit: b.unit,
        values: b.rows.map((row, k) => ({ row, id: ids ? String(ids[row]) : String(row), time: times ? times[row] : undefined, value: b.values[k] })),
      });
    } else if (b.kind === 'each' && P && b.columns.includes(P.column ?? P.lon)) {
      const at = b.columns.indexOf(P.column), when = P.time ? b.columns.indexOf(P.time) : -1;
      for (const m of b.measures) {
        out.push({ title: m.name, scale: 'auto', unit: '', values: b.keys.map((key, k) => ({ id: key[at], time: when >= 0 ? key[when] : undefined, value: m.values[k] })) });
      }
    }
  }
  return out;
}

// A value that does not change in time for any place (a unit's effect, given on each of its rows) stands at one moment: the
// values once for each place, without a time; otherwise as they are.
export function stillOrNot(finite) {
  const first = new Map();
  const still = finite.every((v) => { const a = first.get(v.id); if (a === undefined) { first.set(v.id, v); return true; } return Math.abs(a.value - v.value) <= 1e-12 * Math.max(1, Math.abs(a.value)); });
  return still && finite.some((v) => v.time !== undefined) ? [...first.values()].map((v) => ({ id: v.id, row: v.row, value: v.value })) : finite;
}

// Rings for a place: its polygon, its small hexagon, or an economy's cells of the world's mosaic.
async function ringsFor(P, economies, mosaic) {
  if (economies || P?.kind === 'economies') { const by = await mosaic(); return (id) => by.get(id) || []; }
  if (P.kind === 'cells') return (id) => (P.rings[id] ? [P.rings[id]] : []);
  if (P.kind === 'points') return (id) => (P.at[id] ? [hexAbout(P.at[id][0], P.at[id][1])] : []);
  return () => [];
}

// A record of the line of time from what is to be drawn. `mosaic()` answers the world's cells by economy, for values of economies.
export async function recordOf(d, table, source, mosaic = async () => new Map()) {
  const P = table?.places;
  const named = P?.names || {};
  const kind = P?.time_kind;
  const finite = d.values.filter((v) => Number.isFinite(v.value));
  if (!finite.length) throw new Error('nothing here has a value to draw');
  const name = d.title.length > 70 ? `${d.title.slice(0, 68)}…` : d.title;

  // Each row its own point.
  if (P?.kind === 'row') {
    const lon = table.data[P.lon], lat = table.data[P.lat], time = P.time ? table.data[P.time] : null;
    return {
      kind: 'events', title: name, credit: source,
      time: finite.map((v) => (time ? moment(time[v.row], kind) : Date.UTC(2000, 0, 1))),
      lat: finite.map((v) => lat[v.row]), lon: finite.map((v) => lon[v.row]), values: { [name]: finite.map((v) => v.value) },
    };
  }

  const values = stillOrNot(finite);

  const rings = await ringsFor(P, d.economies, mosaic);
  const times = [...new Set(values.map((v) => v.time))];
  const stamped = times.map((label) => ({ label, t: label === undefined ? Date.UTC(2000, 0, 1) : moment(label, kind) })).sort((a, b) => a.t - b.t);
  const step = new Map(stamped.map((s, i) => [s.label, i]));
  const ids = [...new Set(values.map((v) => v.id))];
  const cells = [], ofId = new Map();
  for (const id of ids) {
    const at = [];
    rings(id).forEach((ring, k) => {
      at.push(cells.length);
      cells.push({ id: `${id}:${k}`, name: named[id] || id, ring });
    });
    ofId.set(id, at);
  }
  if (!cells.length) throw new Error('none of these places can be drawn');
  const frames = stamped.map(() => new Array(cells.length).fill(null));
  for (const v of values) for (const c of ofId.get(v.id)) frames[step.get(v.time)][c] = v.value;
  return { kind: 'cells', title: name, credit: source, time: stamped.map((s) => s.t), cells, values: { [name]: frames }, units: d.unit ? { [name]: d.unit } : {} };
}

// Which colours and height suit it: a value about zero, or a magnitude; and how tall the table's places stand.
export function styleOf(d, record, table) {
  const flat = record.kind === 'cells' ? record.values[record.title].flat().filter(Number.isFinite) : record.values[record.title];
  const low = Math.min(...flat), high = Math.max(...flat);
  const signed = d.scale === 'diverging' || (low < 0 && high > 0);
  return { toolbox: true, map: signed ? 'lavender-rust' : 'plum', scale: signed ? 'diverging' : high / Math.max(low, 1e-9) > 40 ? 'sqrt' : 'linear', credit: record.credit, ...(table?.places?.top ? { top: table.places.top } : {}), ...(table?.places?.base ? { base: table.places.base } : {}) };
}

// The place an entry names, for reading: the highest and lowest at the last moment drawn.
export function ranking(d, table, count = 5) {
  const names = table?.places?.names || {};
  const finite = stillOrNot(d.values.filter((v) => Number.isFinite(v.value)));
  const last = finite.reduce((a, v) => (v.time !== undefined && (a === undefined || String(v.time) > String(a)) ? v.time : a), undefined);
  const now = finite.filter((v) => v.time === last).sort((a, b) => b.value - a.value);
  const row = (v) => ({ id: v.id, row: v.row, name: names[v.id] || v.id, value: v.value });
  return { last, count: now.length, high: now.slice(0, count).map(row), low: now.slice(-count).reverse().map(row) };
}

// Where a place is, to go there.
export function whereIs(id, table) {
  const P = table?.places;
  if (!P) return null;
  if (P.kind === 'cells') return P.rings[id] ? centreOf(P.rings[id]) : null;
  return P.at?.[id] ?? null;
}
