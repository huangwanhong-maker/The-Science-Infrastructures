// Data in space and time: what happened where and when, in the forms a plug-in may give or be given (ADR 0079, and api/python/science/spacetime.py, which sets out
// the same contract for Python and C++). Space is degrees on WGS 84, time milliseconds since 1970-01-01 in UTC. A record is a field on a part of the Earth
// times a line of time, in one of five forms, by what it is sampled on: events (points), tracks (curves), cells (a partition, at steps), grid (a raster, at
// steps) and flows (between places).
//
// This module holds the arithmetic, checkable in Node and with nothing of a page in it, so that a script run in a worker can import it:
//   events, tracks, cells, grid, flows, table, series, note, param, SpaceTimeAnalysis    the contract, as in Python
//   readTable, guessColumns, recordFromTable, recordFromGeoJSON, readRecords             a file of one's own made into records
//   parseTime                                                                            a date as a file writes it, or a refusal with its reason
//   prepared                                                                             a record checked and put in columns, with what was left out and why
//   histogram, timeTicks, formatTime, spanOf                                             a line of time and how to draw it

export const SECOND = 1000, MINUTE = 60000, HOUR = 3600000, DAY = 86400000, YEAR = 31556952000;
export const KINDS = ['events', 'tracks', 'cells', 'grid', 'flows'];

// ---- the contract ---------------------------------------------------------------------------------------------------------------------------------

const common = (kind, title, { units, show, note: text } = {}) => ({ kind, title, units: units || {}, ...(show ? { show } : {}), ...(text ? { note: text } : {}) });
export const events = (title, { time, lat, lon, values = {}, labels, ...rest }) => ({ ...common('events', title, rest), time, lat, lon, values, ...(labels ? { labels } : {}) });
export const tracks = (title, items, rest = {}) => ({ ...common('tracks', title, rest), tracks: items });
export const cells = (title, { time, cells: items, values, ...rest }) => ({ ...common('cells', title, rest), time, cells: items, values });
export const grid = (title, { time, box, width, height, values, ...rest }) => ({ ...common('grid', title, rest), time, box, width, height, values });
export const flows = (title, { places, flows: items, ...rest }) => ({ ...common('flows', title, rest), places, flows: items });
export const table = (title, columns, rows, text = '') => ({ kind: 'table', title, columns, rows, note: text });
export const series = (title, lines, { xlabel = '', ylabel = '', time_axis = false, note: text = '' } = {}) => ({ kind: 'series', title, lines, xlabel, ylabel, time_axis, note: text });
export const note = (text) => ({ kind: 'note', text });
export const param = (name, label, value, { kind = typeof value === 'boolean' ? 'boolean' : 'number', min, max, step, options, help } = {}) => ({ name, label, kind, default: value, min, max, step, options, help });
export class SpaceTimeAnalysis {
  id = ''; title = ''; description = ''; params = [];
  describe() { return { id: this.id, title: this.title, description: this.description, kind: 'spacetime-analysis', language: 'javascript', params: this.params }; }
}

// ---- reading a date ---------------------------------------------------------------------------------------------------------------------------------

const ISO = /^(\d{4})(?:-(\d{1,2})(?:-(\d{1,2})(?:[T ](\d{1,2}):(\d{2})(?::(\d{2})(?:\.(\d+))?)?\s*(Z|[+-]\d{2}(?::?\d{2})?)?)?)?)?$/i;

// A date as a file writes it, to ms UTC: { t, assumed } or { why }. ISO 8601 with its zone is read as it says. Without one it is taken as UTC, and said to be
// (`assumed`). A number is taken by its size: a year, seconds, or milliseconds since 1970. A form that could be read two ways (03/04/2020) is refused.
export function parseTime(value) {
  if (value === null || value === undefined || value === '') return { why: 'no time' };
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) return { why: 'no time' };
    if (Number.isInteger(value) && value >= 1000 && value <= 3000) return { t: Date.UTC(value, 0, 1), assumed: true };
    if (Math.abs(value) > 1e11) return { t: value };
    if (Math.abs(value) > 1e8) return { t: value * 1000 };
    return { why: `the number ${value} is no year, seconds or milliseconds` };
  }
  const s = String(value).trim();
  if (/^-?\d+(\.\d+)?$/.test(s) && !/^\d{4}$/.test(s)) return parseTime(Number(s));
  const m = ISO.exec(s);
  if (m) {
    const [, y, mo = 1, d = 1, h = 0, mi = 0, se = 0, frac = '', zone] = m;
    const year = +y, month = +mo, day = +d;
    if (month < 1 || month > 12 || day < 1 || day > 31 || +h > 24 || +mi > 59 || +se > 60) return { why: `"${s}" is no date` };
    let t = Date.UTC(year, month - 1, day, +h, +mi, +se, Math.round(+`0.${frac || 0}` * 1000));
    if (year < 100) { const dt = new Date(t); dt.setUTCFullYear(year); t = dt.getTime(); }
    if (zone && zone.toUpperCase() !== 'Z') {
      const sign = zone[0] === '-' ? -1 : 1, digits = zone.slice(1).replace(':', '');
      t -= sign * (+digits.slice(0, 2) * HOUR + +(digits.slice(2) || 0) * MINUTE);
    }
    return { t, assumed: !zone && m[4] !== undefined };
  }
  if (/^\d{1,2}[/.-]\d{1,2}[/.-]\d{2,4}/.test(s)) return { why: `"${s}" can be read two ways: write it as ISO 8601 (2020-03-04)` };
  const parsed = Date.parse(s);
  if (Number.isFinite(parsed) && /[a-z]{3}/i.test(s)) return { t: parsed, assumed: !/(gmt|utc|z|[+-]\d{2}:?\d{2})\s*$/i.test(s) };
  return { why: `"${s}" is not a date this reads` };
}

// ---- a file of one's own ------------------------------------------------------------------------------------------------------------------------------

// A delimited text as columns and rows of strings: commas, semicolons or tabs, whichever divides its first line most; quotes as RFC 4180 has them.
export function readTable(text) {
  const clean = text.replace(/^﻿/, '');
  const first = clean.slice(0, clean.indexOf('\n') < 0 ? clean.length : clean.indexOf('\n'));
  const delimiter = [',', ';', '\t'].map((d) => [d, first.split(d).length]).sort((a, b) => b[1] - a[1])[0][0];
  const rows = [];
  let row = [], cell = '', quoted = false;
  for (let i = 0; i < clean.length; i++) {
    const c = clean[i];
    if (quoted) {
      if (c === '"' && clean[i + 1] === '"') { cell += '"'; i++; } else if (c === '"') quoted = false; else cell += c;
    } else if (c === '"' && !cell) quoted = true;
    else if (c === delimiter) { row.push(cell); cell = ''; } else if (c === '\n' || c === '\r') {
      if (c === '\r' && clean[i + 1] === '\n') i++;
      row.push(cell); cell = '';
      if (row.some((x) => x !== '')) rows.push(row);
      row = [];
    } else cell += c;
  }
  if (cell || row.length) { row.push(cell); if (row.some((x) => x !== '')) rows.push(row); }
  const columns = (rows.shift() || []).map((c) => c.trim());
  return { columns, rows };
}

const NAMES = {
  lat: /^(lat|latitude|y|lat_?deg|latitude_?deg)$/i,
  lon: /^(lon|lng|long|longitude|x|lon_?deg|longitude_?deg)$/i,
  time: /^(time|date|datetime|date_?time|timestamp|t|when|utc|time_?utc|time_?ms|occurred|year|started?|observed)(_at)?$/i,
  until: /^(until|end|end_?time|finish(ed)?|stop(ped)?)$/i,
  track: /^(track|track_?id|trajectory|traj_?id|trip|trip_?id|vehicle|vessel|mmsi|device|device_?id|animal|individual|individual_?id|storm|sid)$/i,
  label: /^(name|label|title|place|description|desc|location)$/i,
};
// Which columns are which, by their names and then their contents. `numbers` are the columns that read as numbers in most rows.
export function guessColumns({ columns, rows }) {
  const out = { lat: null, lon: null, time: null, until: null, track: null, label: null, numbers: [] };
  for (const key of ['lat', 'lon', 'time', 'until', 'track', 'label']) out[key] = columns.find((c) => NAMES[key].test(c)) ?? null;
  const sample = rows.slice(0, 400);
  // A column of ids or names whose values repeat says which thing each row is: those rows are a track.
  if (!out.track) {
    const k = columns.findIndex((c) => /^(id|uid|name|object|object_?id|agent|agent_?id)$/i.test(c));
    if (k >= 0 && rows.length && new Set(rows.map((r) => r[k])).size <= rows.length / 2) out.track = columns[k];
  }
  columns.forEach((c, k) => {
    if ([out.lat, out.lon, out.time, out.until].includes(c)) return;
    const known = sample.filter((r) => r[k] !== '' && r[k] !== undefined);
    if (known.length && known.filter((r) => Number.isFinite(Number(r[k]))).length >= 0.9 * known.length) out.numbers.push(c);
  });
  if (!out.time) {
    const dated = columns.find((c, k) => c !== out.lat && c !== out.lon && !out.numbers.includes(c) && sample.length && sample.filter((r) => Number.isFinite(parseTime(r[k]).t)).length >= 0.8 * sample.length);
    if (dated) out.time = dated;
  }
  return out;
}

const coordinate = (cell) => (cell === undefined || cell === null || String(cell).trim() === '' ? NaN : Number(cell));
// A table as a record: its rows events, or tracks where a column says which thing each row is. Rows that cannot be placed are left out, and counted by reason in `left`.
export function recordFromTable(table_, mapping, title = 'A table') {
  const at = Object.fromEntries(table_.columns.map((c, k) => [c, k]));
  const need = ['lat', 'lon', 'time'].filter((k) => !mapping[k] || !(mapping[k] in at));
  if (need.length) throw new Error(`Say which column is the ${need.join(', the ')}.`);
  const left = {};
  const drop = (why) => { left[why] = (left[why] || 0) + 1; };
  let assumed = 0;
  const numbers = (mapping.numbers || []).filter((c) => c in at);
  const put = [];
  for (const r of table_.rows) {
    // An empty cell is no coordinate: Number('') is 0, which would put the row on the equator and stretch the whole record to reach it.
    const lat = coordinate(r[at[mapping.lat]]), lon = coordinate(r[at[mapping.lon]]);
    if (!Number.isFinite(lat) || !Number.isFinite(lon) || Math.abs(lat) > 90 || Math.abs(lon) > 360) { drop('no place on the Earth'); continue; }
    const time = parseTime(r[at[mapping.time]]);
    if (time.why) { drop(time.why.includes('two ways') ? time.why : 'no time that reads'); continue; }
    if (time.assumed) assumed++;
    put.push({ lat, lon, t: time.t, r });
  }
  const valueOf = (x, c) => { const v = Number(x.r[at[c]]); return x.r[at[c]] === '' || !Number.isFinite(v) ? null : v; };
  const report = { left, assumedUTC: assumed, rows: table_.rows.length };
  if (mapping.track && mapping.track in at) {
    const by = new Map();
    for (const x of put) { const id = x.r[at[mapping.track]]; if (!by.has(id)) by.set(id, []); by.get(id).push(x); }
    const items = [...by].map(([id, xs]) => {
      xs.sort((a, b) => a.t - b.t);
      return { id, name: id, time: xs.map((x) => x.t), lat: xs.map((x) => x.lat), lon: xs.map((x) => x.lon), values: Object.fromEntries(numbers.map((c) => [c, xs.map((x) => valueOf(x, c))])) };
    });
    return { ...tracks(title, items), report };
  }
  return {
    ...events(title, {
      time: put.map((x) => x.t), lat: put.map((x) => x.lat), lon: put.map((x) => x.lon),
      values: Object.fromEntries(numbers.map((c) => [c, put.map((x) => valueOf(x, c))])),
      labels: mapping.label && mapping.label in at ? put.map((x) => x.r[at[mapping.label]]) : undefined,
    }),
    report,
  };
}

// GeoJSON as a record (RFC 7946): Points with a time are events; LineStrings with `times` are tracks, or with a fourth coordinate; Polygons with a time and an id are cells.
const TIME_KEYS = ['time', 'date', 'datetime', 'timestamp', 'when', 'start', 'mag_time'];
export function recordFromGeoJSON(data, title = 'A GeoJSON') {
  const features = data.type === 'FeatureCollection' ? data.features : data.type === 'Feature' ? [data] : [];
  const timeOf = (p) => { for (const k of TIME_KEYS) if (p && p[k] !== undefined && p[k] !== null) { const t = parseTime(p[k]); if (!t.why) return t.t; } return null; };
  const points = features.filter((f) => f.geometry?.type === 'Point'), lines = features.filter((f) => /LineString$/.test(f.geometry?.type || '')), areas = features.filter((f) => /Polygon$/.test(f.geometry?.type || ''));
  const numeric = (fs) => { const names = new Set(); for (const f of fs.slice(0, 200)) for (const [k, v] of Object.entries(f.properties || {})) if (typeof v === 'number' && !TIME_KEYS.includes(k)) names.add(k); return [...names]; };
  if (points.length >= lines.length && points.length >= areas.length) {
    const kept = points.filter((f) => timeOf(f.properties) !== null);
    if (!kept.length) throw new Error('These points have no time: a property named time, date or timestamp is read.');
    const names = numeric(kept);
    return { ...events(title, { time: kept.map((f) => timeOf(f.properties)), lon: kept.map((f) => f.geometry.coordinates[0]), lat: kept.map((f) => f.geometry.coordinates[1]),
      values: Object.fromEntries(names.map((n) => [n, kept.map((f) => (typeof f.properties[n] === 'number' ? f.properties[n] : null))])),
      labels: kept.map((f) => f.properties?.name ?? f.properties?.title ?? f.properties?.place ?? '') }), report: { left: points.length > kept.length ? { 'no time': points.length - kept.length } : {}, rows: points.length } };
  }
  if (lines.length >= areas.length) {
    const items = [];
    lines.forEach((f, k) => {
      const parts = f.geometry.type === 'LineString' ? [f.geometry.coordinates] : f.geometry.coordinates;
      const times = f.properties?.times || f.properties?.coordTimes;
      parts.forEach((c, j) => {
        const t = times ? (Array.isArray(times[0]) ? times[j] : times).map((x) => parseTime(x).t) : c.map((p) => parseTime(p[3]).t);
        if (t.some((x) => !Number.isFinite(x))) return;
        items.push({ id: f.id ?? f.properties?.id ?? `${k}${parts.length > 1 ? `.${j}` : ''}`, name: f.properties?.name, time: t, lon: c.map((p) => p[0]), lat: c.map((p) => p[1]) });
      });
    });
    if (!items.length) throw new Error('These lines have no times: a property `times` beside their coordinates, or a fourth coordinate, is read.');
    return tracks(title, items);
  }
  // Polygons: a cell for each id, a step for each time.
  const stamps = [...new Set(areas.map((f) => timeOf(f.properties)).filter((t) => t !== null))].sort((a, b) => a - b);
  if (!stamps.length) throw new Error('These polygons have no time: a property named time, date or timestamp is read.');
  const ids = [...new Set(areas.map((f, k) => f.id ?? f.properties?.id ?? f.properties?.name ?? k))];
  const names = numeric(areas);
  const place = (f) => ids.indexOf(f.id ?? f.properties?.id ?? f.properties?.name ?? areas.indexOf(f));
  const values = Object.fromEntries(names.map((n) => [n, stamps.map(() => new Array(ids.length).fill(null))]));
  const rings = new Array(ids.length).fill(null);
  for (const f of areas) {
    const k = place(f), step = stamps.indexOf(timeOf(f.properties));
    if (!rings[k]) rings[k] = f.geometry.type === 'Polygon' ? f.geometry.coordinates[0] : f.geometry.coordinates[0][0];
    if (step >= 0) for (const n of names) if (typeof f.properties[n] === 'number') values[n][step][k] = f.properties[n];
  }
  return cells(title, { time: stamps, cells: ids.map((id, k) => ({ id: String(id), ring: rings[k].slice(0, -1) })), values });
}

// A file's text as records, by what it is: our own JSON ({kind} or {records} or a list), GeoJSON, JSON lines, or a delimited table (given back to be mapped).
export function readRecords(text, name = 'A file') {
  const title = name.replace(/\.[a-z0-9]+$/i, '');
  const head = text.trimStart()[0];
  if (head === '{' || head === '[') {
    let data = null;
    try { data = JSON.parse(text); } catch {
      const lines = text.split(/\r?\n/).filter((l) => l.trim());
      try { const rows = lines.map((l) => JSON.parse(l)); const columns = [...new Set(rows.flatMap((r) => Object.keys(r)))]; return { table: { columns, rows: rows.map((r) => columns.map((c) => (r[c] === undefined || r[c] === null ? '' : String(r[c])))) }, title }; } catch { throw new Error(`${name} is not JSON that this reads.`); }
    }
    if (data.type === 'FeatureCollection' || data.type === 'Feature') return { records: [recordFromGeoJSON(data, title)], title };
    if (KINDS.includes(data.kind)) return { records: [{ title, ...data }], title };
    if (Array.isArray(data.records)) return { records: data.records, title };
    if (Array.isArray(data) && data.every((r) => KINDS.includes(r?.kind))) return { records: data, title };
    if (Array.isArray(data) && data.length && typeof data[0] === 'object') { const columns = [...new Set(data.flatMap((r) => Object.keys(r)))]; return { table: { columns, rows: data.map((r) => columns.map((c) => (r[c] === undefined || r[c] === null ? '' : String(r[c])))) }, title }; }
    throw new Error(`${name} is JSON, but of no form this reads: a record of ${KINDS.join(', ')}, GeoJSON, or a list of rows.`);
  }
  return { table: readTable(text), title };
}

// ---- a record in columns ---------------------------------------------------------------------------------------------------------------------------------

const finiteRange = (values) => { let a = Infinity, b = -Infinity; for (const v of values) if (Number.isFinite(v)) { if (v < a) a = v; if (v > b) b = v; } return a <= b ? [a, b] : null; };

// A record checked and put in columns for drawing: time as Float64Array, its span and box, each variable's range; what could not be drawn is left out, and counted in `left`. A form that cannot be read at all is refused with its reason.
export function prepared(record) {
  if (!record || !KINDS.includes(record.kind)) throw new Error(`A record is one of ${KINDS.join(', ')}.`);
  // What the file's reader left out counts with what is left out here.
  const left = { ...(record.report?.left || {}) };
  const drop = (why, n = 1) => { left[why] = (left[why] || 0) + n; };
  const out = { kind: record.kind, title: record.title || 'A record', units: record.units || {}, show: record.show || {}, note: record.note || '', left, report: record.report || null };
  const inPlace = (la, lo) => Number.isFinite(la) && Number.isFinite(lo) && Math.abs(la) <= 90 && Math.abs(lo) <= 360;
  let south = 90, north = -90, west = 360, east = -360;
  const see = (la, lo) => { if (la < south) south = la; if (la > north) north = la; if (lo < west) west = lo; if (lo > east) east = lo; };
  const stats = (cols) => Object.fromEntries(Object.entries(cols).map(([k, v]) => [k, finiteRange(v.flat ? v.flat(2) : v)]));
  if (record.kind === 'events') {
    const n = record.time?.length ?? 0;
    if (record.lat?.length !== n || record.lon?.length !== n) throw new Error('An event has a time, a latitude and a longitude: these columns are not of one length.');
    const order = [];
    for (let k = 0; k < n; k++) {
      if (!Number.isFinite(record.time[k])) { drop('no time'); continue; }
      if (!inPlace(record.lat[k], record.lon[k])) { drop('no place on the Earth'); continue; }
      order.push(k);
    }
    order.sort((a, b) => record.time[a] - record.time[b]);
    out.time = Float64Array.from(order, (k) => record.time[k]);
    out.lat = Float64Array.from(order, (k) => record.lat[k]);
    out.lon = Float64Array.from(order, (k) => record.lon[k]);
    out.values = Object.fromEntries(Object.entries(record.values || {}).map(([name, c]) => [name, Float32Array.from(order, (k) => (c[k] === null || c[k] === undefined ? NaN : c[k]))]));
    out.labels = record.labels ? order.map((k) => record.labels[k]) : null;
    for (let i = 0; i < out.time.length; i++) see(out.lat[i], out.lon[i]);
    out.count = out.time.length;
    out.times = out.time;
  } else if (record.kind === 'tracks') {
    out.tracks = [];
    for (const t of record.tracks || []) {
      const keep = [];
      for (let k = 0; k < (t.time?.length ?? 0); k++) { if (Number.isFinite(t.time[k]) && inPlace(t.lat[k], t.lon[k])) keep.push(k); else drop('a point with no time or place'); }
      if (keep.length < 2) { drop('a track of fewer than two points'); continue; }
      keep.sort((a, b) => t.time[a] - t.time[b]);
      out.tracks.push({ id: t.id, name: t.name, time: Float64Array.from(keep, (k) => t.time[k]), lat: Float64Array.from(keep, (k) => t.lat[k]), lon: Float64Array.from(keep, (k) => t.lon[k]),
        values: Object.fromEntries(Object.entries(t.values || {}).map(([name, c]) => [name, Float32Array.from(keep, (k) => (c[k] === null || c[k] === undefined ? NaN : c[k]))])) });
      for (const k of keep) see(t.lat[k], t.lon[k]);
    }
    const names = [...new Set(out.tracks.flatMap((t) => Object.keys(t.values)))];
    out.values = Object.fromEntries(names.map((n) => [n, out.tracks.flatMap((t) => (t.values[n] ? [...t.values[n]] : []))]));
    out.count = out.tracks.length;
    out.times = Float64Array.from(out.tracks.flatMap((t) => [...t.time]));
  } else if (record.kind === 'cells') {
    const steps = record.time?.length ?? 0, n = record.cells?.length ?? 0;
    if (!steps || !n) throw new Error('Cells have a line of time and places on it.');
    out.time = Float64Array.from(record.time);
    out.cells = record.cells;
    for (const c of record.cells) { if (!Array.isArray(c.ring) || c.ring.length < 3) throw new Error(`The cell ${c.id} has no ring of [lon, lat].`); for (const [lo, la] of c.ring) see(la, lo); }
    out.values = {};
    for (const [name, frames] of Object.entries(record.values || {})) {
      if (frames.length !== steps) throw new Error(`The values of ${name} have ${frames.length} steps, and its line of time ${steps}.`);
      const a = new Float32Array(steps * n);
      frames.forEach((row, s) => { if (row.length !== n) throw new Error(`A step of ${name} has ${row.length} values for ${n} cells.`); row.forEach((v, k) => { a[s * n + k] = v === null || v === undefined ? NaN : v; }); });
      out.values[name] = a;
    }
    out.count = n;
    out.times = out.time;
  } else if (record.kind === 'grid') {
    const steps = record.time?.length ?? 0, w = record.width, h = record.height;
    if (!steps || !(w > 0) || !(h > 0) || record.box?.length !== 4) throw new Error('A grid has a line of time, a box [south, west, north, east], a width and a height.');
    out.time = Float64Array.from(record.time);
    out.width = w, out.height = h, out.box = record.box;
    [south, west, north, east] = record.box;
    out.values = {};
    for (const [name, frames] of Object.entries(record.values || {})) {
      if (frames.length !== steps) throw new Error(`The values of ${name} have ${frames.length} steps, and its line of time ${steps}.`);
      out.values[name] = frames.map((f) => { if (f.length !== w * h) throw new Error(`A step of ${name} has ${f.length} values for a grid of ${w} × ${h}.`); return Float32Array.from(f, (v) => (v === null || v === undefined ? NaN : v)); });
    }
    out.count = w * h;
    out.times = out.time;
  } else {
    const at = new Map();
    for (const p of record.places || []) { if (inPlace(p.lat, p.lon)) { at.set(p.id, p); see(p.lat, p.lon); } else drop('a place with no position'); }
    out.places = at;
    out.flows = [];
    for (const f of record.flows || []) {
      if (!at.has(f.from) || !at.has(f.to)) { drop('a flow between unknown places'); continue; }
      if (!Number.isFinite(f.time)) { drop('a flow with no time'); continue; }
      out.flows.push(f);
    }
    out.flows.sort((a, b) => a.time - b.time);
    out.values = { value: out.flows.map((f) => f.value ?? null) };
    out.count = out.flows.length;
    out.times = Float64Array.from(out.flows, (f) => f.time);
  }
  if (!out.count) throw new Error(`Nothing in ${out.title} could be placed in space and time${Object.keys(left).length ? `: ${Object.entries(left).map(([k, v]) => `${v} with ${k}`).join(', ')}` : ''}.`);
  out.ranges = record.kind === 'cells' ? Object.fromEntries(Object.entries(out.values).map(([k, v]) => [k, finiteRange(v)])) : record.kind === 'grid' ? Object.fromEntries(Object.entries(out.values).map(([k, v]) => [k, finiteRange(v.flatMap((f) => [...f]))])) : stats(out.values);
  out.span = spanOf(out.times);
  out.box = out.box || [south, west, north, east];
  out.variables = Object.keys(out.values).filter((k) => out.ranges[k]);
  return out;
}

export const spanOf = (times) => { let a = Infinity, b = -Infinity; for (const t of times) { if (t < a) a = t; if (t > b) b = t; } return [a, b]; };

// The first index with times[k] >= t, in a sorted column.
export function lowerBound(times, t) {
  let a = 0, b = times.length;
  while (a < b) { const m = (a + b) >> 1; if (times[m] < t) a = m + 1; else b = m; }
  return a;
}

// ---- a line of time ------------------------------------------------------------------------------------------------------------------------------------

// How many times fall in each of n equal bins over a span.
export function histogram(times, [t0, t1], n) {
  const out = new Float32Array(n);
  const width = t1 > t0 ? (t1 - t0) / n : 1;
  for (const t of times) { const k = Math.min(n - 1, Math.floor((t - t0) / width)); if (k >= 0) out[k]++; }
  return out;
}

const STEPS = [
  [SECOND, 'second'], [5 * SECOND, 'second'], [15 * SECOND, 'second'], [MINUTE, 'minute'], [5 * MINUTE, 'minute'], [15 * MINUTE, 'minute'], [HOUR, 'hour'], [3 * HOUR, 'hour'], [6 * HOUR, 'hour'], [12 * HOUR, 'hour'],
  [DAY, 'day'], [2 * DAY, 'day'], [7 * DAY, 'week'], [14 * DAY, 'week'],
];
// Round moments within a span, about `most` of them, in UTC: [{ t, unit }]. Months and years are counted on the calendar.
export function timeTicks(t0, t1, most = 8) {
  const span = t1 - t0, want = span / Math.max(1, most);
  if (!(span > 0)) return [{ t: t0, unit: 'second' }];
  if (want >= 25 * DAY) {
    const months = want < 40 * DAY ? 1 : want < 80 * DAY ? 2 : want < 100 * DAY ? 3 : want < 200 * DAY ? 6 : 0;
    const out = [];
    const d = new Date(t0);
    if (months) {
      let y = d.getUTCFullYear(), m = d.getUTCMonth() - (d.getUTCMonth() % months);
      for (;;) {
        const t = Date.UTC(y, m);
        if (t > t1) break;
        if (t >= t0) out.push({ t, unit: m === 0 ? 'year' : 'month' });
        m += months;
        if (m >= 12) { m -= 12; y++; }
      }
      return out;
    }
    const years = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000].find((s) => s * YEAR >= want) || 1000;
    for (let y = Math.ceil(d.getUTCFullYear() / years) * years; ; y += years) {
      const t = Date.UTC(y, 0);
      if (t > t1) break;
      if (t >= t0) out.push({ t, unit: 'year' });
    }
    return out;
  }
  const [step, unit] = STEPS.find(([s]) => s >= want) || [28 * DAY, 'week'];
  const out = [];
  for (let t = Math.ceil(t0 / step) * step; t <= t1; t += step) out.push({ t, unit });
  return out;
}

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
const pad = (v) => String(v).padStart(2, '0');
// A moment as fits a scale, in UTC: a tick's label, or (with unit 'full', fitting a span) the date of the cursor.
export function formatTime(t, unit = 'day', span = 0) {
  const d = new Date(t), y = d.getUTCFullYear(), mo = MONTHS[d.getUTCMonth()], day = d.getUTCDate(), hm = `${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}`;
  if (unit === 'full') return span > 40 * YEAR ? String(y) : span > 2 * YEAR ? `${mo} ${y}` : span > 20 * DAY ? `${day} ${mo} ${y}` : span > DAY ? `${day} ${mo} ${hm}` : `${hm}:${pad(d.getUTCSeconds())}`;
  if (unit === 'year') return String(y);
  if (unit === 'month') return mo;
  if (unit === 'week' || unit === 'day') return `${day} ${mo}`;
  if (unit === 'hour' || unit === 'minute') return hm;
  return `${hm}:${pad(d.getUTCSeconds())}`;
}

// How long a window reads: "3 days".
export function formatSpan(ms) {
  const say = (v, unit) => `${v >= 10 ? Math.round(v) : Math.round(v * 10) / 10} ${unit}${v === 1 ? '' : 's'}`;
  if (ms >= 2 * YEAR) return say(ms / YEAR, 'year');
  if (ms >= 60 * DAY) return say(ms / (YEAR / 12), 'month');
  if (ms >= 2 * DAY) return say(ms / DAY, 'day');
  if (ms >= 2 * HOUR) return say(ms / HOUR, 'hour');
  if (ms >= 2 * MINUTE) return say(ms / MINUTE, 'minute');
  return say(ms / SECOND, 'second');
}
