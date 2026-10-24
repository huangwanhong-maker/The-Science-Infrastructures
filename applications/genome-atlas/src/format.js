// How numbers and places are written: lengths in bp, kb and Mb; positions 1-based, as biologists
// read them (the engine's are 0-based and half-open); fractions as percentages.

const group = (n) => Math.round(n).toLocaleString('en-US');

// A length: 5,386 bp; 48.5 kb; 4.64 Mb.
export function length(n) {
  if (!Number.isFinite(n)) return '–';
  if (n < 10000) return `${group(n)} bp`;
  if (n < 1e6) return `${(n / 1000).toFixed(n < 1e5 ? 1 : 0)} kb`;
  return `${(n / 1e6).toFixed(n < 1e7 ? 2 : 1)} Mb`;
}

// A 0-based position written 1-based: 3,925,598.
export const position = (p) => group(p + 1);
// A half-open interval [b, e) written as biologists do: 101..250.
export const interval = (b, e) => `${group(b + 1)}..${group(e)}`;
export const count = (n) => group(n);
export const percent = (x, digits = 1) => (Number.isFinite(x) ? `${(100 * x).toFixed(digits)}%` : '–');
export function seconds(s) {
  if (!Number.isFinite(s)) return '–';
  if (s < 1e-3) return `${(s * 1e6).toFixed(0)} µs`;
  if (s < 1) return `${(s * 1e3).toFixed(s < 0.01 ? 2 : 0)} ms`;
  return `${s.toFixed(s < 10 ? 2 : 1)} s`;
}
export function bytes(b) {
  if (!Number.isFinite(b)) return '–';
  const units = ['B', 'KB', 'MB', 'GB', 'TB'];
  let k = 0;
  while (b >= 1024 && k < units.length - 1) { b /= 1024; k++; }
  return `${b.toFixed(k === 0 || b >= 100 ? 0 : 1)} ${units[k]}`;
}
export function number(x, digits = 3) {
  if (!Number.isFinite(x)) return '–';
  const a = Math.abs(x);
  if (a !== 0 && (a < 1e-3 || a >= 1e6)) return x.toExponential(2).replace('e+', 'e');
  return Number(x.toPrecision(digits)).toLocaleString('en-US', { maximumFractionDigits: 6 });
}

// A tick spacing of 1, 2 or 5 times a power of ten for about n ticks over a span.
export function tickStep(span, n = 8) {
  const raw = Math.max(span, 1e-12) / n, p = 10 ** Math.floor(Math.log10(raw)), f = raw / p;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
}
// A tick's label on a genome: 0, 500, 1 kb, 1.5 kb, 2 Mb.
export function tickLabel(p, step, total = 0) {
  if (step >= 1e6 || p >= 1e7 || (total >= 2e6 && step >= 1e5)) return `${+(p / 1e6).toFixed(step >= 1e6 ? 0 : step >= 1e5 ? 1 : 2)} Mb`;
  if (step >= 1000 || p >= 1e4) return `${+(p / 1000).toFixed(step >= 1000 ? 0 : step >= 100 ? 1 : 2)} kb`;
  return group(p);
}
