// Units: from the units in which data are stored (SI or as the source gives
// them, in UDUNITS notation) to the units in which they are shown.
//
// Provisional stand-in for physis/units (laboratory ADR 0014).
//
// Data stay in their stored units throughout; conversion happens only when a
// number is shown or a colour scale is set (CLAUDE.md, section 4.1). Every
// conversion is affine, x_shown = a x + b, so a colour scale in shown units is
// the same scale in stored units with its ends converted back.

const AFFINE = {
  'K': { unit: '°C', a: 1, b: -273.15 },
  'degree_C': { unit: '°C', a: 1, b: 0 },
  'Pa': { unit: 'hPa', a: 0.01, b: 0 },
  'kg m-2 s-1': { unit: 'mm/h', a: 3600, b: 0 },       // 1 kg of water per m^2 is a depth of 1 mm
  'm s-1': { unit: 'm/s', a: 1, b: 0 },
  'm': { unit: 'm', a: 1, b: 0 },
  's': { unit: 's', a: 1, b: 0 },
  'degree': { unit: '°', a: 1, b: 0 },
  '%': { unit: '%', a: 1, b: 0 },
  '1': { unit: '', a: 1, b: 0 },
  'W m-2': { unit: 'W/m²', a: 1, b: 0 },
  'kg m-2': { unit: 'mm', a: 1, b: 0 },                // of liquid water: 1 kg per m^2 is 1 mm
  'J kg-1': { unit: 'J/kg', a: 1, b: 0 },
  's-1': { unit: '10⁻⁵ s⁻¹', a: 1e5, b: 0 },           // vorticity, in the unit of synoptic charts
  'mm yr-1': { unit: 'mm/yr', a: 1, b: 0 },            // plate motion
  'km': { unit: 'km', a: 1, b: 0 },
};

// Alternative units the reader may choose for a quantity stored in given units.
export const ALTERNATIVES = {
  'm s-1': [{ unit: 'km/h', a: 3.6, b: 0 }, { unit: 'kn', a: 3600 / 1852, b: 0 }],
  'K': [{ unit: '°F', a: 1.8, b: -459.67 }],
  'degree_C': [{ unit: '°F', a: 1.8, b: 32 }],
};

export function displayUnit(stored, chosen) {
  const base = AFFINE[stored] || { unit: stored, a: 1, b: 0 };
  if (!chosen || chosen === base.unit) return base;
  return (ALTERNATIVES[stored] || []).find((u) => u.unit === chosen) || base;
}

export const toShown = (x, u) => u.a * x + u.b;
export const toStored = (y, u) => (y - u.b) / u.a;

// A number written for reading: significant figures suited to its size.
const significant = (v, digits) => Number(v.toPrecision(digits)).toLocaleString('en-GB', { maximumFractionDigits: 6 });

export function format(x, digits = 3) {
  if (!Number.isFinite(x)) return '—';
  const ax = Math.abs(x);
  if (ax >= 1e12) return `${significant(x / 1e12, digits)} T`;
  if (ax >= 1e9) return `${significant(x / 1e9, digits)} bn`;
  if (ax >= 1e6) return `${significant(x / 1e6, digits)} M`;
  if (ax >= 1e4) return Math.round(x).toLocaleString('en-GB');
  if (ax >= 100) return x.toFixed(0);
  if (ax >= 10) return x.toFixed(1);
  if (ax >= 1) return x.toFixed(digits - 1 > 2 ? 2 : digits - 1);
  if (ax === 0) return '0';
  return x.toPrecision(2);
}
