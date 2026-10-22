// Data mapped onto the brain: a table of values by region, coordinates with
// values, and channel values at electrodes, parsed from text someone pastes.
//
// A table is lines of fields separated by commas, tabs or semicolons; a
// header line is skipped when its fields are not numbers. By region: the
// first field names the region (by name, acronym, label or a FreeSurfer key,
// with the hemisphere written any way) and the second is the value. By
// coordinate: x, y, z, then optionally a value and a name. By electrode: the
// electrode's name (Fz, C3, ...) and the value.

import { colourAt } from '../render/colormaps.js';

const splitLine = (line) => line.split(/\t|;|,(?=(?:[^"]*"[^"]*")*[^"]*$)/).map((f) => f.trim().replace(/^"|"$/g, ''));
const isNumber = (s) => s !== '' && Number.isFinite(Number(s));

export function parseTable(text) {
  const rows = [];
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith('#')) continue;
    rows.push(splitLine(line));
  }
  if (rows.length && !rows[0].some(isNumber)) rows.shift();     // a header
  return rows;
}

// Values by region: {values: Map<Region, number>, unmatched: string[], count}.
export function valuesByRegion(parcellation, text) {
  const values = new Map(), unmatched = [];
  for (const row of parseTable(text)) {
    if (row.length < 2) continue;
    const [name, value] = [row[0], Number(row[row.length - 1])];
    if (!Number.isFinite(value)) { unmatched.push(name); continue; }
    const region = parcellation.match(isNumber(name) ? Number(name) : name);
    if (!region) { unmatched.push(name); continue; }
    values.set(region, value);
  }
  return { values, unmatched, count: values.size };
}

// Points: [{position: [x, y, z], value, name}].
export function parsePoints(text) {
  const points = [], bad = [];
  for (const row of parseTable(text)) {
    const nums = row.slice(0, 3).map(Number);
    if (row.length < 3 || nums.some((x) => !Number.isFinite(x))) { bad.push(row.join(', ')); continue; }
    const rest = row.slice(3);
    const value = rest.length && isNumber(rest[0]) ? Number(rest[0]) : null;
    const name = rest.find((f) => !isNumber(f)) || null;
    points.push({ position: nums, value, name });
  }
  return { points, bad };
}

// Channel values: {values: Map<name, number>, unknown: string[]}.
export function channelValues(text, names) {
  const known = new Map(names.map((n) => [n.toLowerCase(), n]));
  const values = new Map(), unknown = [];
  for (const row of parseTable(text)) {
    if (row.length < 2) continue;
    const name = known.get(row[0].toLowerCase().replace(/^t3$/, 't7').replace(/^t4$/, 't8').replace(/^t5$/, 'p7').replace(/^t6$/, 'p8'));
    const value = Number(row[row.length - 1]);
    if (!name || !Number.isFinite(value)) { unknown.push(row[0]); continue; }
    values.set(name, value);
  }
  return { values, unknown };
}

// A scale for values: [low, high] from the data, symmetric about zero for a diverging map.
export function valueRange(values, { diverging = false } = {}) {
  const v = [...values].filter(Number.isFinite);
  if (!v.length) return [0, 1];
  let lo = Math.min(...v), hi = Math.max(...v);
  if (diverging) { const m = Math.max(Math.abs(lo), Math.abs(hi)) || 1; lo = -m; hi = m; }
  if (hi === lo) { hi = lo + 1; }
  return [lo, hi];
}

// The colour of a value on a map, as [r, g, b] bytes.
export const colourOfValue = (value, [lo, hi], map) => colourAt(map, (value - lo) / (hi - lo)).map(Math.round);

// A few significant figures for a legend.
export function tidy(x) {
  if (!Number.isFinite(x)) return '—';
  const a = Math.abs(x);
  if (a === 0) return '0';
  if (a >= 1000) return x.toFixed(0);
  if (a >= 10) return x.toFixed(1);
  if (a >= 1) return x.toFixed(2);
  return x.toPrecision(2);
}
