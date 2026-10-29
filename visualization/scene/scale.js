// Scales: where a value falls on a colour map (the infrastructure's ADR 0053; formalism Definition 3.10,
// the channels of a scalar).
//
// A scale maps a value v to u in [0, 1], the place on the map, and has its domain and its ticks:
//
//   linear     u = (v − lo)/(hi − lo);
//   log        u = log(v/lo)/log(hi/lo), for positive values spanning decades, lo raised to a floor
//              when one is given (the values below it drawn at the floor);
//   diverging  u = ½ + ½ (v − c)/m about a centre c with a natural meaning (zero for a signed
//              component, a reference for a departure), m = max |v − c| so that the two arms are
//              alike; with a floor f, logarithmic in |v − c| from f to m on each arm
//              (u = ½ ± ½ log(|v − c|/f)/log(m/f), values within f of the centre at it).
//
// "auto" chooses diverging when a centre is given, log when the values are positive and span two
// decades or more, linear otherwise. A value that is not finite has no place (NaN), and is not drawn.

const clamp01 = (u) => Math.min(1, Math.max(0, u));

function finiteRange(values) {
  let lo = Infinity, hi = -Infinity;
  for (const v of values) if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
  return hi >= lo ? [lo, hi] : null;
}

// "Nice" steps of 1, 2 and 5 times a power of ten, about `count` of them across [lo, hi].
export function niceTicks(lo, hi, count = 5) {
  if (!(hi > lo)) return Number.isFinite(lo) ? [lo] : [];
  const raw = (hi - lo) / Math.max(1, count), p = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * p).find((s) => s >= raw) || 10 * p;
  const out = [];
  for (let k = Math.ceil(lo / step - 1e-9); k * step <= hi + 1e-9 * step; k++) out.push(Math.abs(k * step) < 1e-12 * step ? 0 : k * step);
  return out;
}

export function scaleOf(values, { kind = 'auto', centre = null, floor = null, domain = null } = {}) {
  const range = domain || finiteRange(values);
  if (!range) return null;
  let [lo, hi] = range;
  if (kind === 'auto') kind = centre !== null && centre !== undefined ? 'diverging' : lo > 0 && hi / lo >= 100 ? 'log' : 'linear';
  if (kind === 'diverging') {
    const c = centre ?? 0;
    const m = Math.max(Math.abs(lo - c), Math.abs(hi - c)) || Math.abs(c) * 1e-9 || 1;
    const f = floor > 0 && floor < m ? floor : null;
    const u = f
      ? (v) => {
        if (!Number.isFinite(v)) return NaN;
        const d = v - c, a = Math.abs(d);
        return a <= f ? 0.5 : 0.5 + 0.5 * Math.sign(d) * clamp01(Math.log(a / f) / Math.log(m / f));
      }
      : (v) => (Number.isFinite(v) ? clamp01(0.5 + (0.5 * (v - c)) / m) : NaN);
    const ticks = (count = 5) => {
      if (!f) return niceTicks(c - m, c + m, count).map((value) => ({ value, u: u(value) }));
      const out = [{ value: c, u: 0.5 }];
      for (let e = Math.ceil(Math.log10(f) - 1e-9); 10 ** e <= m * (1 + 1e-9); e++) for (const s of [1, -1]) out.push({ value: c + s * 10 ** e, u: u(c + s * 10 ** e) });
      return out.sort((a, b) => a.value - b.value);
    };
    return { kind, domain: [c - m, c + m], centre: c, floor: f, u, t: (v) => 2 * u(v) - 1, ticks };
  }
  if (kind === 'log') {
    const base = floor > 0 ? Math.max(lo, floor) : lo;
    if (base > 0 && hi > base) {
      const a = Math.log10(base), s = Math.log10(hi) - a;
      const u = (v) => (Number.isFinite(v) ? clamp01((Math.log10(Math.max(v, base)) - a) / s) : NaN);
      const ticks = () => {
        const out = [];
        for (let e = Math.ceil(a - 1e-9); e <= Math.log10(hi) + 1e-9; e++) out.push({ value: 10 ** e, u: u(10 ** e) });
        if (out.length < 3) for (const m of [2, 5]) for (let e = Math.floor(a); e <= Math.log10(hi); e++) { const v = m * 10 ** e; if (v >= base && v <= hi) out.push({ value: v, u: u(v) }); }
        return out.sort((x, y) => x.value - y.value);
      };
      return { kind, domain: [base, hi], centre: null, floor: floor > 0 ? floor : null, u, t: u, ticks };
    }
    kind = 'linear';
  }
  const s = hi - lo || 1;
  const u = (v) => (Number.isFinite(v) ? clamp01((v - lo) / s) : NaN);
  return { kind: 'linear', domain: [lo, hi], centre: null, floor: null, u, t: u, ticks: (count = 5) => niceTicks(lo, hi, count).map((value) => ({ value, u: u(value) })) };
}
