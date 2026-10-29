// The legend of an expression: what is shown, in what unit and frame, on what scale, from what source and
// of what provenance, and where it could not be read and why (the infrastructure's ADR 0053). Data, not
// HTML: an application writes it in its own style.
//
// Units are the engine's SI symbols (physis units: "T", "K", "Pa", "m s^-2"); a value is written with the
// SI prefix that brings it between 1 and 1000 (nT for the Moon's crust, μT for the Earth's field), or as
// given for units without prefixes (K).

import { colourMap } from '../styles/palettes.js';

const PRETTY = { 'm s^-2': 'm/s²', 'm s^-1': 'm/s', 'kg m^-3': 'kg/m³', 'm^2': 'm²', 'm^3': 'm³', 'mol/m^3': 'mol/m³', 1: '' };
export const unitText = (unit) => PRETTY[unit] ?? unit;

const PREFIXES = [[1e-12, 'p'], [1e-9, 'n'], [1e-6, 'μ'], [1e-3, 'm'], [1, ''], [1e3, 'k'], [1e6, 'M'], [1e9, 'G']];
const PREFIXED = new Set(['T', 'Pa', 'm', 's', 'V', 'A', 'W', 'J', 'N', 'Hz']);

// The prefix for values of a size: the largest whose multiple does not exceed it.
export function prefixFor(size, unit) {
  if (!PREFIXED.has(unit) || !(size > 0)) return { factor: 1, symbol: unitText(unit) };
  let best = PREFIXES[0];
  for (const p of PREFIXES) if (size >= p[0] * 0.9995) best = p;
  return { factor: best[0], symbol: `${best[1]}${unit}` };
}

// A value written for reading, in a unit with its prefix: `figures` significant figures (three by
// default), the value first rounded to them, so that 99.999… nT is written 100 nT.
export function quantityText(v, unit, prefix = null, figures = 3) {
  if (!Number.isFinite(v)) return '—';
  const p = prefix || prefixFor(Math.abs(v), unit);
  const x = Number((v / p.factor).toPrecision(figures)), a = Math.abs(x);
  const decimals = a === 0 ? 0 : Math.max(0, figures - 1 - Math.floor(Math.log10(a)));
  const digits = a === 0 ? '0' : x.toFixed(Math.min(decimals, 12));
  return `${digits.replace('-', '−')}${p.symbol ? ` ${p.symbol}` : ''}`;
}

// The labels of ticks with the fewest significant figures, from three, that tell them apart.
function tickLabels(ticks, unit, prefix) {
  for (let figures = 3; figures <= 8; figures++) {
    const labels = ticks.map((t) => quantityText(t.value, unit, prefix, figures));
    if (new Set(labels).size === labels.length) return labels;
  }
  return ticks.map((t) => quantityText(t.value, unit, prefix, 8));
}

// The legend of an expression drawn from samples on a scale.
export function legendOf(expression, scale, samples = null) {
  const { field, colour } = expression;
  const map = colourMap(colour.map);
  const top = scale ? Math.max(Math.abs(scale.domain[0]), Math.abs(scale.domain[1])) : 0;
  const prefix = prefixFor(top, field.fibre.unit);
  const raw = scale ? scale.ticks(4) : [];
  const labels = tickLabels(raw, field.fibre.unit, prefix);
  const ticks = raw.map((t, k) => ({ ...t, label: labels[k] }));
  const counts = samples?.counts || null;
  const refused = samples ? samples.reasons.map((r) => ({ text: r.text, count: r.count, share: r.count / samples.size })) : [];
  return {
    title: expression.title,
    shows: expression.shows,
    unit: prefix.symbol,
    fibreUnit: field.fibre.unit,
    frame: expression.frame,
    scale: scale ? { kind: scale.kind, floor: scale.floor, centre: scale.centre, domain: scale.domain } : null,
    ticks,
    map: { name: map.name, kind: map.kind, gradient: map.gradient() },
    provenance: field.provenance,
    parts: field.parts || [],
    source: field.source,
    counts,
    refused,
    notes: samples ? samples.notes : [],
  };
}
