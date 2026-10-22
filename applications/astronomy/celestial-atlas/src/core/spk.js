// Reader and evaluator for JPL SPK kernels with type-2 (Chebyshev position) segments,
// such as the DE440 planetary ephemeris.
//
// Provisional stand-in for the ephemeris evaluator planned in
// core/domains/astronomy (laboratory ADR 0011).
//
// Sources: NAIF, "DAF Required Reading" and "SPK Required Reading" (type 2).
// Units: positions in km, velocities in km/s, time in seconds past J2000 TDB,
// frame ICRF (SPK frame 1, "J2000").

const RECORD_BYTES = 1024;

function ascii(bytes, start, length) {
  return String.fromCharCode(...bytes.subarray(start, start + length));
}

// Chebyshev series and its derivative with respect to s, for coefficients
// c[off .. off + n - 1] at s in [-1, 1].
function chebyshev(c, off, n, s) {
  let t0 = 1, t1 = s, d0 = 0, d1 = 1;
  let value = c[off] + (n > 1 ? c[off + 1] * s : 0);
  let deriv = n > 1 ? c[off + 1] : 0;
  for (let k = 2; k < n; k++) {
    const t2 = 2 * s * t1 - t0;
    const d2 = 2 * t1 + 2 * s * d1 - d0;
    value += c[off + k] * t2;
    deriv += c[off + k] * d2;
    t0 = t1; t1 = t2; d0 = d1; d1 = d2;
  }
  return [value, deriv];
}

export class SPK {
  constructor(buffer) {
    const bytes = new Uint8Array(buffer);
    const view = new DataView(buffer);
    const id = ascii(bytes, 0, 8);
    if (!id.startsWith('DAF/SPK')) throw new Error(`not an SPK kernel: "${id}"`);
    const format = ascii(bytes, 88, 8);
    if (format !== 'LTL-IEEE') throw new Error(`unsupported binary format ${format}`);
    if (new Uint8Array(new Uint16Array([1]).buffer)[0] !== 1) {
      throw new Error('a little-endian platform is required');
    }

    const nd = view.getInt32(8, true);
    const ni = view.getInt32(12, true);
    const summaryDoubles = nd + Math.ceil(ni / 2);
    this.data = new Float64Array(buffer, 0, Math.floor(buffer.byteLength / 8));
    this.segments = [];

    let record = view.getInt32(76, true);   // first summary record (FWARD)
    while (record > 0) {
      const base = (record - 1) * RECORD_BYTES;
      const next = view.getFloat64(base, true);
      const count = view.getFloat64(base + 16, true);
      for (let i = 0; i < count; i++) {
        const off = base + 24 + i * summaryDoubles * 8;
        const ints = off + nd * 8;
        const segment = {
          startEt: view.getFloat64(off, true),
          endEt: view.getFloat64(off + 8, true),
          target: view.getInt32(ints, true),
          center: view.getInt32(ints + 4, true),
          frame: view.getInt32(ints + 8, true),
          type: view.getInt32(ints + 12, true),
          begin: view.getInt32(ints + 16, true),   // 1-based double addresses
          end: view.getInt32(ints + 20, true),
        };
        if (segment.type === 2) {
          const d = this.data, e = segment.end;
          segment.init = d[e - 4];
          segment.intlen = d[e - 3];
          segment.rsize = d[e - 2];
          segment.n = d[e - 1];
          segment.ncoef = (segment.rsize - 2) / 3;
        }
        this.segments.push(segment);
      }
      record = next;
    }

    this.coverage = {
      startEt: Math.max(...this.segments.map((s) => s.startEt)),
      endEt: Math.min(...this.segments.map((s) => s.endEt)),
    };
  }

  static async fromUrl(url, onProgress) {
    const response = await fetch(url);
    if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
    const total = Number(response.headers.get('Content-Length')) || 0;
    if (!onProgress || !response.body || !total) {
      return new SPK(await response.arrayBuffer());
    }
    const buffer = new Uint8Array(total);
    const reader = response.body.getReader();
    let done = 0;
    for (;;) {
      const { value, done: finished } = await reader.read();
      if (finished) break;
      buffer.set(value, done);
      done += value.length;
      onProgress(done / total);
    }
    return new SPK(buffer.buffer);
  }

  segmentFor(target, et) {
    for (const s of this.segments) {
      if (s.target === target && et >= s.startEt && et <= s.endEt) return s;
    }
    return null;
  }

  covers(et) {
    return et >= this.coverage.startEt && et <= this.coverage.endEt;
  }

  // State of one segment's target relative to its centre.
  segmentState(segment, et) {
    if (segment.type !== 2) throw new Error(`SPK type ${segment.type} is not supported`);
    let index = Math.floor((et - segment.init) / segment.intlen);
    index = Math.min(Math.max(index, 0), segment.n - 1);
    const off = segment.begin - 1 + index * segment.rsize;
    const d = this.data;
    const mid = d[off], radius = d[off + 1];
    const s = (et - mid) / radius;
    const n = segment.ncoef;
    const x = chebyshev(d, off + 2, n, s);
    const y = chebyshev(d, off + 2 + n, n, s);
    const z = chebyshev(d, off + 2 + 2 * n, n, s);
    return {
      p: [x[0], y[0], z[0]],
      v: [x[1] / radius, y[1] / radius, z[1] / radius],
    };
  }

  // State of a body relative to the solar-system barycentre (NAIF id 0),
  // following the chain of segments, e.g. Moon = (3 wrt 0) + (301 wrt 3).
  stateSSB(target, et) {
    const p = [0, 0, 0], v = [0, 0, 0];
    let body = target;
    while (body !== 0) {
      const segment = this.segmentFor(body, et);
      if (!segment) throw new Error(`no ephemeris for body ${body} at ${et} s`);
      const s = this.segmentState(segment, et);
      for (let i = 0; i < 3; i++) { p[i] += s.p[i]; v[i] += s.v[i]; }
      body = segment.center;
    }
    return { p, v };
  }

  // State of target relative to observer, both by NAIF id.
  state(target, observer, et) {
    const a = this.stateSSB(target, et);
    const b = this.stateSSB(observer, et);
    return {
      p: [a.p[0] - b.p[0], a.p[1] - b.p[1], a.p[2] - b.p[2]],
      v: [a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]],
    };
  }
}
