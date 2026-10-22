// Tropical cyclones found and followed in a forecast's fields of sea-level
// pressure, after TempestExtremes (Ullrich and Zarzycki 2017; Zarzycki and
// Ullrich 2017).
//
// Provisional stand-in for core/domains/meteorology (laboratory ADR 0022).
//
// 1. Detection, at each time ("DetectNodes"). A candidate is a grid point
//    whose pressure is the lowest within `mergeDeg` of it (so that two
//    candidates are never closer). It is kept if a closed contour surrounds
//    it: walking outward along any of `rays` great circles, the pressure rises
//    by at least `closedPa` before `closedDeg` of arc. Zarzycki and Ullrich
//    use 200 Pa within 5.5 degrees; their warm-core test (the thickness from
//    300 to 500 hPa) needs fields that not every product carries, and is
//    replaced here by a requirement on the pressure itself (`maxPa`).
// 2. Linking, from time to time ("StitchNodes"). A track takes, at each next
//    time, the nearest candidate within `stepDeg` per 6 hours (8 degrees in
//    TempestExtremes), each candidate going to one track only, nearest pairs
//    first. A track is kept if it lasts at least `minHours`, deepens to
//    `peakPa` or below, and begins equatorward of `genesisLat` and, when a
//    test `isOcean(lon, lat)` is given, over the sea, as a tropical cyclone
//    does; the summer heat lows over land, and the lows of the westerlies,
//    are left out.
//
// Positions are refined by fitting a paraboloid to the 3 x 3 neighbourhood of
// the minimum (a tenth of a cell is typical on a 1 degree grid).
//
// References:
//   Ullrich, P. A. and Zarzycki, C. M. (2017), TempestExtremes: a framework for
//   scale-insensitive pointwise feature tracking on unstructured grids,
//   Geosci. Model Dev. 10, 1069-1090.
//   Zarzycki, C. M. and Ullrich, P. A. (2017), Assessing sensitivities in
//   algorithmic detection of tropical cyclones in climate data, Geophys. Res.
//   Lett. 44, 1141-1149.

const DEG = Math.PI / 180;

export const DEFAULTS = { mergeDeg: 6, closedPa: 200, closedDeg: 5.5, rays: 16, maxPa: 100500,
  stepDeg: 8, minHours: 24, genesisLat: 25, maxLat: 50, peakPa: 100000, isOcean: null };

// Great-circle angle (degrees) between two places in degrees.
export function arcDeg(lon1, lat1, lon2, lat2) {
  const s = Math.sin(((lat2 - lat1) * DEG) / 2) ** 2
    + Math.cos(lat1 * DEG) * Math.cos(lat2 * DEG) * Math.sin(((lon2 - lon1) * DEG) / 2) ** 2;
  return (2 * Math.asin(Math.min(1, Math.sqrt(s)))) / DEG;
}

// The place reached along a great circle from (lon, lat) at a bearing (radians) after d degrees.
function travel(lon, lat, bearing, d) {
  const p = lat * DEG, l = lon * DEG, a = d * DEG;
  const p2 = Math.asin(Math.sin(p) * Math.cos(a) + Math.cos(p) * Math.sin(a) * Math.cos(bearing));
  const l2 = l + Math.atan2(Math.sin(bearing) * Math.sin(a) * Math.cos(p), Math.cos(a) - Math.sin(p) * Math.sin(p2));
  return [((l2 / DEG + 540) % 360) - 180, p2 / DEG];
}

// Candidates at one time: [{ lon, lat, pa }]. values: Float32Array on a regular
// grid { lon0, dlon, nlon, lat0, dlat, nlat, periodic }, rows north to south.
export function detect(values, grid, options = {}) {
  const o = { ...DEFAULTS, ...options };
  const { nlon, nlat, dlon, dlat, lon0, lat0 } = grid;
  const at = (r, c) => values[r * nlon + ((c % nlon) + nlon) % nlon];
  // The pressure at a place, by its nearest grid point (enough for the contour test).
  const sample = (lon, lat) => {
    const r = Math.round((lat0 - lat) / dlat);
    if (r < 0 || r >= nlat) return NaN;
    return at(r, Math.round((lon - lon0) / dlon));
  };
  const rMerge = Math.ceil(o.mergeDeg / dlat);
  const out = [];
  for (let r = 1; r < nlat - 1; r++) {
    const lat = lat0 - r * dlat;
    if (Math.abs(lat) > o.maxLat) continue;
    const cMerge = Math.ceil(o.mergeDeg / (dlon * Math.max(0.2, Math.cos(lat * DEG))));
    for (let c = 0; c < nlon; c++) {
      const v = at(r, c);
      if (!(v <= o.maxPa)) continue;
      // The lowest within mergeDeg (ties go to the first in grid order).
      let lowest = true;
      for (let dr = -rMerge; dr <= rMerge && lowest; dr++) {
        const rr = r + dr;
        if (rr < 0 || rr >= nlat) continue;
        for (let dc = -cMerge; dc <= cMerge; dc++) {
          if (!dr && !dc) continue;
          const w = at(rr, c + dc);
          if (w < v || (w === v && (dr < 0 || (dr === 0 && dc < 0)))) {
            if (arcDeg(lon0 + c * dlon, lat, lon0 + (c + dc) * dlon, lat0 - rr * dlat) <= o.mergeDeg) { lowest = false; break; }
          }
        }
      }
      if (!lowest) continue;
      const lon = lon0 + c * dlon;
      // A closed contour: along every ray, the pressure rises by closedPa within closedDeg.
      let closed = true;
      for (let k = 0; k < o.rays && closed; k++) {
        const bearing = (2 * Math.PI * k) / o.rays;
        let rises = false;
        for (let d = Math.min(dlat, dlon) / 2; d <= o.closedDeg + 1e-9; d += Math.min(dlat, dlon) / 2) {
          const [ln, lt] = travel(lon, lat, bearing, d);
          if (sample(ln, lt) - v >= o.closedPa) { rises = true; break; }
        }
        closed = rises;
      }
      if (!closed) continue;
      // A paraboloid through the 3 x 3 neighbourhood places the minimum between grid points.
      const fx = (at(r, c + 1) - at(r, c - 1)) / 2, fxx = at(r, c + 1) - 2 * v + at(r, c - 1);
      const fy = (at(r + 1, c) - at(r - 1, c)) / 2, fyy = at(r + 1, c) - 2 * v + at(r - 1, c);
      const ox = fxx > 0 ? Math.max(-0.5, Math.min(0.5, -fx / fxx)) : 0;
      const oy = fyy > 0 ? Math.max(-0.5, Math.min(0.5, -fy / fyy)) : 0;
      const pa = v - (fxx > 0 ? (fx * fx) / (2 * fxx) : 0) - (fyy > 0 ? (fy * fy) / (2 * fyy) : 0);
      out.push({ lon: ((lon + ox * dlon + 540) % 360) - 180, lat: lat - oy * dlat, pa });
    }
  }
  return out;
}

// Tracks from candidates at successive times: [{ points: [{ t, lon, lat, pa }] }].
// times: ms; nodes: an array (per time) of candidate arrays.
export function stitch(times, nodes, options = {}) {
  const o = { ...DEFAULTS, ...options };
  const open = [], done = [];
  for (let i = 0; i < times.length; i++) {
    const here = nodes[i].map((n) => ({ ...n, t: times[i], taken: false }));
    const hours = i ? (times[i] - times[i - 1]) / 3600000 : 6;
    const reach = (o.stepDeg * hours) / 6;
    // Nearest pairs first: each candidate to one track, each track to one candidate.
    const pairs = [];
    open.forEach((track, k) => {
      const last = track.points[track.points.length - 1];
      here.forEach((n, j) => {
        const d = arcDeg(last.lon, last.lat, n.lon, n.lat);
        if (d <= reach) pairs.push([d, k, j]);
      });
    });
    pairs.sort((a, b) => a[0] - b[0]);
    const continued = new Set();
    for (const [, k, j] of pairs) {
      if (continued.has(k) || here[j].taken) continue;
      continued.add(k);
      here[j].taken = true;
      open[k].points.push({ t: here[j].t, lon: here[j].lon, lat: here[j].lat, pa: here[j].pa });
    }
    // Tracks not continued end; untaken candidates begin tracks.
    for (let k = open.length - 1; k >= 0; k--) if (!continued.has(k) && i > 0) done.push(...open.splice(k, 1));
    for (const n of here) if (!n.taken) open.push({ points: [{ t: n.t, lon: n.lon, lat: n.lat, pa: n.pa }] });
  }
  done.push(...open);
  return done.map((track) => ({ ...track, minPa: Math.min(...track.points.map((p) => p.pa)) })).filter((track) => {
    const p = track.points;
    return (p[p.length - 1].t - p[0].t) / 3600000 >= o.minHours && track.minPa <= o.peakPa
      && Math.abs(p[0].lat) <= o.genesisLat && (!o.isOcean || o.isOcean(p[0].lon, p[0].lat));
  })
    .sort((a, b) => a.minPa - b.minPa);
}

// A position at an instant, between the two positions that bracket it (longitude across the date line).
export function positionAt(points, t) {
  if (!points.length || t < points[0].t || t > points[points.length - 1].t) return null;
  let k = 0;
  while (k < points.length - 2 && points[k + 1].t < t) k++;
  const a = points[k], b = points[Math.min(k + 1, points.length - 1)];
  const w = b.t === a.t ? 0 : (t - a.t) / (b.t - a.t);
  let dl = b.lon - a.lon;
  if (dl > 180) dl -= 360; else if (dl < -180) dl += 360;
  const pa = a.pa !== undefined && b.pa !== undefined ? a.pa + (b.pa - a.pa) * w : undefined;
  return { lon: ((a.lon + dl * w + 540) % 360) - 180, lat: a.lat + (b.lat - a.lat) * w, pa };
}

// Tracks in a field product: the variable of sea-level pressure, all its times.
export function trackProduct(product, name, options = {}) {
  const nodes = product.times.map((_, i) => detect(product.slice(name, i), product.grid, options));
  return stitch(product.times, nodes, options);
}
