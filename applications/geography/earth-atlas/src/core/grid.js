// Values on a regular longitude-latitude grid: where a place falls on the grid,
// interpolation in space and time, and statistics.
//
// Provisional stand-in for core/inference (interpolation; laboratory ADR 0014).
//
// A grid is { lon0, dlon, nlon, lat0, dlat, nlat, periodic }: the centre of
// column c is at longitude lon0 + c dlon, the centre of row r at latitude
// lat0 - r dlat (rows run north to south), as written by
// data/products/field_grid.py. Values are row-major Float32Arrays with NaN for
// no value.
//
// Interpolation is bilinear between the four surrounding grid points. Where some
// of the four have no value (along a coast, for an ocean field), the weights of
// those that have one are renormalised, provided they carry at least half of the
// total weight; otherwise the place has no value. The fragment shader in
// src/render/surface.js applies the same rule.

export const MISSING = -32768;

// Fractional column and row of a place, or null outside the grid. A place within
// half a cell of the first or last row is clamped onto it.
export function gridPosition(grid, lon, lat) {
  let x = (lon - grid.lon0) / grid.dlon;
  if (grid.periodic) x = ((x % grid.nlon) + grid.nlon) % grid.nlon;
  else if (x < -0.5 || x > grid.nlon - 0.5) return null;
  let y = (grid.lat0 - lat) / grid.dlat;
  if (y < -0.5 || y > grid.nlat - 0.5) return null;
  return [x, Math.max(0, Math.min(grid.nlat - 1, y))];
}

export function bilinear(values, grid, lon, lat) {
  const pos = gridPosition(grid, lon, lat);
  if (!pos) return NaN;
  const [x, y] = pos;
  const nlon = grid.nlon;
  let c0 = Math.floor(x), r0 = Math.floor(y);
  const fx = x - c0, fy = y - r0;
  let c1 = c0 + 1;
  const r1 = Math.min(r0 + 1, grid.nlat - 1);
  if (c1 >= nlon) c1 = grid.periodic ? c1 - nlon : nlon - 1;
  if (c0 < 0) c0 = 0;
  const v00 = values[r0 * nlon + c0], v01 = values[r0 * nlon + c1];
  const v10 = values[r1 * nlon + c0], v11 = values[r1 * nlon + c1];
  const w00 = (1 - fx) * (1 - fy), w01 = fx * (1 - fy), w10 = (1 - fx) * fy, w11 = fx * fy;
  let sum = 0, weight = 0;
  if (v00 === v00) { sum += w00 * v00; weight += w00; }
  if (v01 === v01) { sum += w01 * v01; weight += w01; }
  if (v10 === v10) { sum += w10 * v10; weight += w10; }
  if (v11 === v11) { sum += w11 * v11; weight += w11; }
  return weight >= 0.5 ? sum / weight : NaN;
}

// The bracketing times of an instant in a sorted list of times (ms), and the
// weight of the later one. Outside the list the nearest end is taken, and
// `offset` says by how much (ms; negative before the first time).
export function timeBracket(times, ms, nearest = false) {
  const n = times.length;
  if (n === 0) return { i0: -1, i1: -1, w: 0, shown: null, offset: 0 };
  if (ms <= times[0]) return { i0: 0, i1: 0, w: 0, shown: times[0], offset: ms - times[0] };
  if (ms >= times[n - 1]) return { i0: n - 1, i1: n - 1, w: 0, shown: times[n - 1], offset: ms - times[n - 1] };
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (times[mid] <= ms) lo = mid; else hi = mid;
  }
  const w = (ms - times[lo]) / (times[hi] - times[lo]);
  if (nearest) {
    const i = w < 0.5 ? lo : hi;
    return { i0: i, i1: i, w: 0, shown: times[i], offset: ms - times[i] };
  }
  return { i0: lo, i1: hi, w, shown: ms, offset: 0 };
}

// A field product: a grid, times, and int16-coded variables (see field_grid.py).
export class FieldProduct {
  constructor(meta, buffer) {
    this.meta = meta;
    this.grid = meta.grid;
    this.times = (meta.times || []).map((t) => Date.parse(t));
    this.codes = new Int16Array(buffer);
    this.variables = new Map(meta.variables.map((v) => [v.name, v]));
    this.cache = new Map();
    this.cacheLimit = 48;
    // A climatology describes a typical year: its times are the months of one
    // year (meta.climatology.year), and any instant is shown by its month.
    this.climatology = meta.climatology || null;
  }

  has(name) { return this.variables.has(name); }

  // Decoded values of a variable at a time index (0 for a field without times).
  slice(name, ti = 0) {
    const key = `${name}@${ti}`;
    const hit = this.cache.get(key);
    if (hit) { this.cache.delete(key); this.cache.set(key, hit); return hit; }
    const v = this.variables.get(name);
    if (!v) throw new Error(`no variable ${name} in ${this.meta.name}`);
    const n = this.grid.nlon * this.grid.nlat;
    const start = v.offsets[Math.max(0, ti)] / 2;
    const out = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const q = this.codes[start + i];
      out[i] = q === MISSING ? NaN : q * v.scale + v.offset;
    }
    this.cache.set(key, out);
    if (this.cache.size > this.cacheLimit) this.cache.delete(this.cache.keys().next().value);
    return out;
  }

  bracket(ms, name) {
    if (this.climatology && this.times.length === 12) {
      const i = new Date(ms).getUTCMonth();
      return { i0: i, i1: i, w: 0, shown: this.times[i], offset: 0, month: i };
    }
    const v = name && this.variables.get(name);
    const nearest = Boolean(v && v.cell_methods && /time:\s*(mean|sum|maximum|minimum)/.test(v.cell_methods));
    return timeBracket(this.times, ms, nearest);
  }

  // A variable's value at a place and instant: bilinear in space, linear in time.
  sample(name, ms, lon, lat) {
    const b = this.bracket(ms, name);
    if (b.i0 < 0) return bilinear(this.slice(name, 0), this.grid, lon, lat);
    const a = bilinear(this.slice(name, b.i0), this.grid, lon, lat);
    if (b.w === 0) return a;
    const c = bilinear(this.slice(name, b.i1), this.grid, lon, lat);
    return a + (c - a) * b.w;
  }

  // A direction (degrees) at a place and instant, interpolated through its unit
  // vector, since the mean of 350 and 10 degrees is 0 and not 180.
  sampleDirection(name, ms, lon, lat) {
    const b = this.bracket(ms, name);
    const at = (i) => {
      const [sn, cs] = this.unitVectors(name, i);
      return [bilinear(sn, this.grid, lon, lat), bilinear(cs, this.grid, lon, lat)];
    };
    let [x, y] = at(Math.max(0, b.i0));
    if (b.w > 0) {
      const [x1, y1] = at(b.i1);
      x += (x1 - x) * b.w; y += (y1 - y) * b.w;
    }
    if (!(x === x && y === y) || (x === 0 && y === 0)) return NaN;
    return ((Math.atan2(x, y) * 180) / Math.PI + 360) % 360;
  }

  unitVectors(name, ti) {
    const key = `${name}@${ti}:unit`;
    let hit = this.cache.get(key);
    if (!hit) {
      const d = this.slice(name, ti);
      const sn = new Float32Array(d.length), cs = new Float32Array(d.length);
      for (let i = 0; i < d.length; i++) {
        const a = (d[i] * Math.PI) / 180;
        sn[i] = Math.sin(a); cs[i] = Math.cos(a);
      }
      hit = [sn, cs];
      this.cache.set(key, hit);
    }
    return hit;
  }

  // A variable's value at a place for each time of the product.
  series(name, lon, lat) {
    return this.times.map((t, i) => bilinear(this.slice(name, i), this.grid, lon, lat));
  }
}

// A field computed in the application (by a model of our own) rather than read
// from the library: the same interface, with values held as Float32Arrays.
// slices: { variable name: [Float32Array per time] }.
export class ComputedField extends FieldProduct {
  constructor(meta, slices) {
    super(meta, new ArrayBuffer(0));
    this.slices = slices;
  }

  slice(name, ti = 0) {
    const list = this.slices[name];
    if (!list) throw new Error(`no variable ${name} in ${this.meta.name}`);
    return list[Math.max(0, Math.min(list.length - 1, ti))];
  }
}

// Summary statistics of a variable, as written by data/products/field_grid.py.
export function statistics(slices) {
  const all = [];
  for (const s of slices) for (let i = 0; i < s.length; i += 7) if (s[i] === s[i]) all.push(s[i]);
  all.sort((a, b) => a - b);
  const at = (p) => all[Math.min(all.length - 1, Math.floor(p * (all.length - 1)))];
  return { min: all[0], max: all[all.length - 1], p02: at(0.02), p98: at(0.98) };
}

// Area-weighted (cos latitude) mean of a gridded field, ignoring missing values.
export function areaMean(values, grid) {
  let sum = 0, weight = 0;
  for (let r = 0; r < grid.nlat; r++) {
    const w = Math.cos(((grid.lat0 - r * grid.dlat) * Math.PI) / 180);
    for (let c = 0; c < grid.nlon; c++) {
      const v = values[r * grid.nlon + c];
      if (v === v) { sum += w * v; weight += w; }
    }
  }
  return weight > 0 ? sum / weight : NaN;
}

// A class product: a category per cell (data/products/class_raster.py), on a
// grid like a field's. Classes are not interpolated: a place takes the class of
// the cell that contains it. values: Uint8Array, 0 for no value.
export class ClassRaster {
  constructor(meta, values) {
    this.meta = meta;
    this.grid = meta.grid;
    this.values = values instanceof Uint8Array ? values : new Uint8Array(values);
    this.classes = new Map(meta.classes.map((c) => [c.value, c]));
  }

  // The index of the cell containing a place, or -1 outside the grid.
  cell(lon, lat) {
    const g = this.grid;
    let c = Math.round((lon - g.lon0) / g.dlon);
    if (g.periodic) c = ((c % g.nlon) + g.nlon) % g.nlon;
    else if (c < 0 || c >= g.nlon) return -1;
    const r = Math.round((g.lat0 - lat) / g.dlat);
    if (r < 0 || r >= g.nlat) return -1;
    return r * g.nlon + c;
  }

  valueAt(lon, lat) {
    const k = this.cell(lon, lat);
    return k < 0 ? 0 : this.values[k];
  }

  // The class at a place ({ value, code, name, colour }), or null.
  classAt(lon, lat) {
    return this.classes.get(this.valueAt(lon, lat)) || null;
  }
}
