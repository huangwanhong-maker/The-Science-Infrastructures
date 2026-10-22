// The ocean beneath a place: quantities of a column of currents and temperature
// (time x depth), such as a current profiler moored there would record.
//
// Provisional stand-in for core/domains/oceanography (laboratory ADR 0023).
//
// Conventions: currents are given by the direction they flow toward (the
// oceanographer's convention, opposite to the wind's); depth is positive down,
// in metres; times are milliseconds since 1970-01-01 UTC.
//
// Sources: de Boyer Montégut, C., Madec, G., Fischer, A. S., Lazar, A. and
// Iudicone, D. (2004), Mixed layer depth over the global ocean: an examination
// of profile data and a profile-based climatology, J. Geophys. Res. 109,
// C12003 (the 0.2 °C criterion from 10 m); Gill (1982), Atmosphere-Ocean
// Dynamics, §8.3 (inertial oscillations); IERS Conventions (2010), Table 1.1
// (the Earth's rate of rotation).

export const OMEGA = 7.292115e-5;              // rad s-1, the Earth's nominal mean angular velocity

// The Coriolis parameter f = 2 Ω sin φ (s-1).
export const coriolis = (latDeg) => 2 * OMEGA * Math.sin((latDeg * Math.PI) / 180);

// The inertial period 2π / |f| (hours): the period of a current turning freely
// on the rotating Earth, half a sidereal day divided by |sin φ|.
export const inertialPeriodHours = (latDeg) => (2 * Math.PI) / Math.abs(coriolis(latDeg)) / 3600;

// The direction a current flows toward, in degrees clockwise from north.
export function flowToward(u, v) {
  if (!Number.isFinite(u) || !Number.isFinite(v) || (u === 0 && v === 0)) return NaN;
  return ((Math.atan2(u, v) * 180) / Math.PI + 360) % 360;
}

// The value of a profile at depth z, by linear interpolation between the levels
// around it; NaN outside the levels with values.
export function atDepth(profile, depths, z) {
  for (let k = 0; k + 1 < depths.length; k++) {
    if (z < depths[k] || z > depths[k + 1]) continue;
    const a = profile[k], b = profile[k + 1];
    if (!Number.isFinite(a) || !Number.isFinite(b)) return z === depths[k] ? a : z === depths[k + 1] ? b : NaN;
    const w = (z - depths[k]) / (depths[k + 1] - depths[k]);
    return a + (b - a) * w;
  }
  return z === depths[depths.length - 1] ? profile[depths.length - 1] : NaN;
}

// The depth of the mixed layer: the shallowest depth below the reference depth
// at which temperature differs from its value there by `delta` (de Boyer
// Montégut et al. 2004: 0.2 °C from 10 m), interpolated linearly between
// levels. { depth, reachesBottom }: when the criterion is not met above the
// deepest level with values, the layer is mixed to that level at least.
export function mixedLayerDepth(temps, depths, { reference = 10, delta = 0.2 } = {}) {
  const t0 = atDepth(temps, depths, reference);
  if (!Number.isFinite(t0)) return { depth: NaN, reachesBottom: false };
  let last = reference;
  for (let k = 0; k < depths.length; k++) {
    if (depths[k] <= reference) continue;
    const t = temps[k];
    if (!Number.isFinite(t)) break;
    if (Math.abs(t - t0) >= delta) {
      // Between the level above (or the reference) and this one.
      const zA = Math.max(reference, depths[k - 1]);
      const tA = depths[k - 1] >= reference ? temps[k - 1] : t0;
      const dA = Math.abs(tA - t0), dB = Math.abs(t - t0);
      const w = dB === dA ? 1 : (delta - dA) / (dB - dA);
      return { depth: zA + (depths[k] - zA) * w, reachesBottom: false };
    }
    last = depths[k];
  }
  return { depth: last, reachesBottom: true };
}

// A function (t, z) -> value over a column, linear in time and in depth between
// the stored times and levels. values[i][k] at times[i], depths[k] (null or NaN
// where missing). Where one of the four neighbours is missing, the valid ones
// carry the weight if they hold at least half of it, as the map's fields do.
export function columnSampler(times, depths, values) {
  const nt = times.length, nz = depths.length;
  const bracket = (axis, n, x) => {
    if (x <= axis[0]) return [0, 0, 0];
    if (x >= axis[n - 1]) return [n - 1, n - 1, 0];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (axis[m] <= x) lo = m; else hi = m; }
    return [lo, hi, (x - axis[lo]) / (axis[hi] - axis[lo])];
  };
  const get = (i, k) => { const v = values[i][k]; return v === null || v === undefined ? NaN : v; };
  return (t, z) => {
    if (z < depths[0] || z > depths[nz - 1]) return NaN;
    const [i0, i1, a] = bracket(times, nt, t);
    const [k0, k1, b] = bracket(depths, nz, z);
    let sum = 0, weight = 0;
    for (const [i, wi] of [[i0, 1 - a], [i1, a]]) {
      for (const [k, wk] of [[k0, 1 - b], [k1, b]]) {
        const w = wi * wk, v = get(i, k);
        if (w > 0 && Number.isFinite(v)) { sum += w * v; weight += w; }
      }
    }
    return weight >= 0.5 ? sum / weight : NaN;
  };
}

// A symmetric range [-r, r] for a signed quantity: the q-th quantile of |x|,
// rounded up to a multiple of `step`, at least `floor`.
export function symmetricRange(values, { q = 0.98, step = 0.1, floor = 0.1 } = {}) {
  const a = values.filter(Number.isFinite).map(Math.abs).sort((x, y) => x - y);
  if (!a.length) return floor;
  const x = a[Math.min(a.length - 1, Math.floor(q * (a.length - 1)))];
  return Math.max(floor, Math.ceil(x / step - 1e-9) * step);
}

// The deepest level with a value at any time (index), or -1.
export function deepestLevel(values) {
  let deepest = -1;
  for (const row of values) for (let k = row.length - 1; k > deepest; k--) if (row[k] !== null && Number.isFinite(row[k])) { deepest = k; break; }
  return deepest;
}
