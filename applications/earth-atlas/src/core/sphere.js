// The Earth as a sphere: positions, local directions, distances and motion
// along great circles.
//
// Provisional stand-in for core/domains/geodesy (laboratory ADR 0014). The
// sphere of mean radius serves drawing, sampling and particle motion; it is not
// a geodetic model (WGS 84 differs by up to 0.19 degrees of latitude and 21 km of
// radius, which no layer of this atlas resolves).
//
// Conventions: longitude east-positive and latitude north-positive, in degrees
// at the interface. Earth-centred, Earth-fixed unit vectors: x toward (0 E, 0 N),
// y toward (90 E, 0 N), z toward the north pole.
//
// Sources:
//   IUGG mean radius R1 = (2a + b) / 3 for WGS 84: Moritz (2000), Geodetic Reference System 1980.
//   Spherical trigonometry of the direct problem: Snyder (1987), Map Projections:
//   A Working Manual, USGS Prof. Paper 1395, eq. 5-5 and 5-6.

export const DEG = Math.PI / 180;
export const EARTH_RADIUS_M = 6371008.8;

export const wrapLongitude = (lon) => ((((lon + 180) % 360) + 360) % 360) - 180;

export function toVector(lon, lat, out = [0, 0, 0]) {
  const l = lon * DEG, p = lat * DEG, c = Math.cos(p);
  out[0] = c * Math.cos(l);
  out[1] = c * Math.sin(l);
  out[2] = Math.sin(p);
  return out;
}

export function fromVector(v) {
  const r = Math.hypot(v[0], v[1], v[2]);
  return [Math.atan2(v[1], v[0]) / DEG, Math.asin(Math.max(-1, Math.min(1, v[2] / r))) / DEG];
}

// Unit vectors pointing east and north at a place (undefined at the poles, where
// east is taken along the 90 E meridian's tangent).
export function eastNorth(lon, lat) {
  const l = lon * DEG, p = lat * DEG;
  const sl = Math.sin(l), cl = Math.cos(l), sp = Math.sin(p), cp = Math.cos(p);
  return [[-sl, cl, 0], [-sp * cl, -sp * sl, cp]];
}

// Great-circle (central) angle between two places, radians (haversine form).
export function centralAngle(lon1, lat1, lon2, lat2) {
  const p1 = lat1 * DEG, p2 = lat2 * DEG;
  const a = Math.sin((p2 - p1) / 2) ** 2 + Math.cos(p1) * Math.cos(p2) * Math.sin((lon2 - lon1) * DEG / 2) ** 2;
  return 2 * Math.asin(Math.min(1, Math.sqrt(a)));
}

export const distanceM = (lon1, lat1, lon2, lat2) => EARTH_RADIUS_M * centralAngle(lon1, lat1, lon2, lat2);

// The place reached from (lon, lat) by moving east_m metres east and north_m
// metres north along the great circle of that initial bearing. Valid through the
// poles. Returns [lon, lat].
export function move(lon, lat, eastM, northM, out = [0, 0]) {
  const d = Math.hypot(eastM, northM) / EARTH_RADIUS_M;
  if (d === 0) { out[0] = lon; out[1] = lat; return out; }
  const bearing = Math.atan2(eastM, northM);
  const p1 = lat * DEG;
  const sp1 = Math.sin(p1), cp1 = Math.cos(p1), sd = Math.sin(d), cd = Math.cos(d);
  const sp2 = Math.max(-1, Math.min(1, sp1 * cd + cp1 * sd * Math.cos(bearing)));
  const dl = Math.atan2(Math.sin(bearing) * sd * cp1, cd - sp1 * sp2);
  out[0] = wrapLongitude(lon + dl / DEG);
  out[1] = Math.asin(sp2) / DEG;
  return out;
}

// Area of the sphere between two latitudes, as a fraction of the whole.
export const zoneFraction = (lat1, lat2) => Math.abs(Math.sin(lat2 * DEG) - Math.sin(lat1 * DEG)) / 2;
