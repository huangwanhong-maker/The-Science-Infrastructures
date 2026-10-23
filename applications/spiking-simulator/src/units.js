// Quantities as the engine writes and reads them: values in SI units, shown with the prefix
// that suits them ("250 pF", "20 ms"), and written back as the specification writes them.

const PREFIXES = [[1e9, 'G'], [1e6, 'M'], [1e3, 'k'], [1, ''], [1e-3, 'm'], [1e-6, 'µ'], [1e-9, 'n'], [1e-12, 'p'], [1e-15, 'f']];
// What a specification may write for "µ".
const ASCII = { 'µ': 'u' };

// The engine names dimensions by symbols; a few compound ones read better in the units a
// neuroscientist writes them in.
const READABLE = { 'kg^-2 m^-4 s^6 A^3': ['nS/mV', 1e-6], 'A/s': ['pA/ms', 1e-9], 'V/s': ['mV/ms', 1] };

// A short number: four significant digits, no trailing zeros.
export function short(x) {
  if (!Number.isFinite(x)) return String(x);
  if (x === 0) return '0';
  const a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) return x.toExponential(3).replace(/\.?0+e/, 'e');
  return String(Number(x.toPrecision(4)));
}

// A value in SI units, with the prefix that brings it between 1 and 1000.
export function format(value, unit, { ascii = false } = {}) {
  if (value === null || value === undefined || Number.isNaN(value)) return '';
  if (unit === '1' || unit === '' || unit === undefined) return short(value);
  if (READABLE[unit]) return `${short(value / READABLE[unit][1])} ${READABLE[unit][0]}`;
  if (!['V', 'A', 'F', 'S', 's', 'Hz', 'Ω'].includes(unit)) return `${short(value)} ${unit}`;
  if (value === 0) return `0 ${unit}`;
  const a = Math.abs(value);
  let [scale, prefix] = PREFIXES.find(([s]) => a >= s * 0.9995) ?? PREFIXES[PREFIXES.length - 1];
  if (unit === 's' && scale > 1) [scale, prefix] = [1, ''];          // seconds, not kiloseconds
  if (unit === 'Hz' && scale < 1) [scale, prefix] = [1, ''];         // hertz, not millihertz
  const p = ascii ? (ASCII[prefix] ?? prefix) : prefix;
  return `${short(value / scale)} ${p}${unit}`;
}

// The unit a field is written in, for its placeholder.
export function hint(unit) {
  if (!unit || unit === '1') return 'a number';
  return READABLE[unit]?.[0] ?? unit;
}

const NUMBER = /^[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?$/i;
const QUANTITY = /^[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?\s+\S+$/i;

// What a field's text is in a specification: a distribution (JSON), a bare number, or a
// quantity with its unit; null when empty. The engine checks the dimension.
export function parseField(text) {
  const t = String(text ?? '').trim();
  if (!t) return null;
  if (t.startsWith('{') || t.startsWith('[')) return JSON.parse(t);
  if (NUMBER.test(t)) return Number(t);
  return t.replace(/µ/g, 'u');
}

// A field's text from its value in a specification.
export function fieldText(value) {
  if (value === null || value === undefined) return '';
  if (typeof value === 'object') return JSON.stringify(value);
  return String(value);
}

// Whether a field's text is well formed (the engine will judge its dimension).
export function wellFormed(text, unit) {
  const t = String(text ?? '').trim();
  if (!t) return true;
  if (t.startsWith('{') || t.startsWith('[')) {
    try { JSON.parse(t); return true; } catch { return false; }
  }
  if (NUMBER.test(t)) return !unit || unit === '1';
  return QUANTITY.test(t.replace(/µ/g, 'u'));
}

// Seconds from a duration as a specification writes it ("2 s", "500 ms", 1.5).
export function seconds(q) {
  if (typeof q === 'number') return q;
  const m = String(q ?? '').trim().match(/^([-+]?[\d.]+(?:e[-+]?\d+)?)\s*(s|ms|us|µs|min)?$/i);
  if (!m) return NaN;
  return Number(m[1]) * ({ s: 1, ms: 1e-3, us: 1e-6, 'µs': 1e-6, min: 60 }[m[2] ?? 's']);
}
