// Tropical cyclones found in the sea-level pressure of field products, by the
// meteorology domain's tracker (domains/meteorology/js/cyclones.js, after
// TempestExtremes), for data/products/verification.py, which matches them to
// the observed storms of IBTrACS.
//
// The tracker is the one the Earth Atlas runs in its page
// (applications/earth-atlas/src/layers/tracks.js), with the same test of
// genesis over the sea, read from a relief product (sea where the height is
// below 0 m). One option differs: a run may begin in the middle of a storm's
// life, poleward of where storms form, so the latitude of a track's first
// position may be up to `genesisLat` (default 40 degrees here, 25 in the page).
//
// A product is read as data/products/field_grid.py writes it: int16 codes,
// value = q * scale + offset, q = -32768 for no value.
//
// Usage, from the repository root:
//   node data/products/track_cyclones.mjs <fields-*.json> <variable> [--relief <fields-*.json>]
//   node data/products/track_cyclones.mjs < job.json > tracks.json
// where job.json is { "relief": path | null, "options": {...}, "products": [{ "path", "variable" }] }.
// The second form prints { "<path>": { "variable", "tracks": [{ "minPa", "points": [{ "t", "lon", "lat", "pa" }] }] } },
// with t in ms since 1970-01-01T00:00:00Z (UTC) and pa in Pa.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';

import { trackProduct } from '../../domains/meteorology/js/index.js';
import { bilinear } from '../../foundations/physis/js/index.js';

const MISSING = -32768;

// A field product with the interface trackProduct reads: times (ms), grid, slice(name, ti).
export function readProduct(path) {
  const meta = JSON.parse(readFileSync(path, 'utf8'));
  const raw = readFileSync(join(dirname(path), meta.file));
  const codes = new Int16Array(raw.buffer, raw.byteOffset, raw.byteLength / 2);
  const variables = new Map(meta.variables.map((v) => [v.name, v]));
  const n = meta.grid.nlon * meta.grid.nlat;
  return {
    meta, grid: meta.grid, times: (meta.times || []).map((t) => Date.parse(t)),
    slice(name, ti = 0) {
      const v = variables.get(name);
      if (!v) throw new Error(`no variable ${name} in ${meta.name}`);
      const start = v.offsets[Math.max(0, ti)] / 2;
      const out = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        const q = codes[start + i];
        out[i] = q === MISSING ? NaN : q * v.scale + v.offset;
      }
      return out;
    },
  };
}

// The test of genesis over the sea, from a relief product's first variable.
export function oceanTest(path) {
  if (!path) return null;
  const relief = readProduct(path);
  const height = relief.slice(relief.meta.variables[0].name, 0);
  return (lon, lat) => bilinear(height, relief.grid, lon, lat) < 0;
}

export function tracksOf(path, variable, options) {
  const tracks = trackProduct(readProduct(path), variable, options);
  return tracks.map((tr) => ({ minPa: Math.round(tr.minPa),
    points: tr.points.map((p) => ({ t: p.t, lon: +p.lon.toFixed(3), lat: +p.lat.toFixed(3), pa: Math.round(p.pa) })) }));
}

function main() {
  const args = process.argv.slice(2);
  if (args.length) {
    const r = args.indexOf('--relief');
    const relief = r >= 0 ? args.splice(r, 2)[1] : null;
    const [path, variable] = args;
    const tracks = tracksOf(path, variable, { genesisLat: 40, isOcean: oceanTest(relief) });
    for (const tr of tracks) {
      const a = tr.points[0], b = tr.points[tr.points.length - 1];
      console.log(`${new Date(a.t).toISOString().slice(0, 13)}Z to ${new Date(b.t).toISOString().slice(0, 13)}Z: `
        + `${tr.points.length} positions, from ${a.lon} ${a.lat}, lowest ${(tr.minPa / 100).toFixed(1)} hPa`);
    }
    return;
  }
  const job = JSON.parse(readFileSync(0, 'utf8'));
  const options = { genesisLat: 40, ...(job.options || {}), isOcean: oceanTest(job.relief) };
  const out = {};
  for (const { path, variable } of job.products) {
    try {
      out[path] = { variable, tracks: tracksOf(path, variable, options) };
    } catch (error) {
      out[path] = { variable, error: String(error.message || error) };
    }
  }
  process.stdout.write(JSON.stringify(out));
}

main();
