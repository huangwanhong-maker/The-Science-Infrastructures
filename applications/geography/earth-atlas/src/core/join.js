// Joining data about places to the geometry of places, by key.
//
// Provisional stand-in for core/domains/geography (GIS; laboratory ADR 0014).
//
// Base geography and thematic data are separate products (laboratory ADR 0014):
// the countries of a geography product carry several systems of identifiers
// (Natural Earth's adm0_a3, ISO 3166-1 alpha-3 and alpha-2, the World Bank's
// codes), and an indicator product is keyed by one system. A join assigns each
// country the record whose key matches one of the country's identifiers, tried
// in a stated order, and reports what remained unmatched on either side, so that
// gaps in a map are explained rather than silent.

// Identifier fields of a country, tried in this order against three-letter keys,
// then against two-letter keys.
export const THREE_LETTER = ['wb_a3', 'iso_a3_eh', 'iso_a3', 'adm0_a3'];
export const TWO_LETTER = ['iso_a2_eh', 'iso_a2', 'wb_a2'];

// countries: [{ keys }]; records: { code: { iso2 } }.
// Returns { codes: [code or null per country], unmatched: [codes], via: {field: count} }.
export function joinCountries(countries, records) {
  const byIso2 = new Map();
  for (const [code, r] of Object.entries(records)) if (r.iso2) byIso2.set(r.iso2, code);
  const used = new Set();
  const via = {};
  const codes = countries.map((c) => {
    for (const field of THREE_LETTER) {
      const k = c.keys[field];
      if (k && records[k]) { used.add(k); via[field] = (via[field] || 0) + 1; return k; }
    }
    for (const field of TWO_LETTER) {
      const k = c.keys[field];
      if (k && byIso2.has(k)) { const code = byIso2.get(k); used.add(code); via[field] = (via[field] || 0) + 1; return code; }
    }
    return null;
  });
  const unmatched = Object.keys(records).filter((k) => !used.has(k)).sort();
  return { codes, unmatched, via };
}

// The value of a record for a year, or for the latest earlier year that has one
// (within `lookback` years). Returns { value, year } or null.
export function valueForYear(record, years, year, lookback = 0) {
  const first = years[0];
  for (let y = Math.min(year, years[years.length - 1]); y >= Math.max(first, year - lookback); y--) {
    const v = record.values[y - first];
    if (v !== null && v !== undefined) return { value: v, year: y };
  }
  return null;
}

// The latest year for which at least `share` of the records have a value.
export function latestWellCoveredYear(records, years, share = 0.8) {
  const list = Object.values(records);
  for (let i = years.length - 1; i >= 0; i--) {
    const n = list.reduce((s, r) => s + (r.values[i] !== null && r.values[i] !== undefined ? 1 : 0), 0);
    if (n >= share * list.length) return years[i];
  }
  return years[years.length - 1];
}
