// A vertical section through the Earth along a great circle: the places along
// it, and where other places (earthquakes) lie relative to it.
//
// Provisional stand-in for core/domains/geodesy (laboratory ADR 0024).
//
// On a sphere of radius R, for a great circle leaving point 1 at azimuth θ12
// and a point 3 at angular distance δ13 and azimuth θ13 from point 1, the
// cross-track distance is R asin(sin δ13 sin(θ13 - θ12)) (positive to the
// right of the direction of travel) and the along-track distance R atan2(
// sin δ13 cos(θ13 - θ12), cos δ13) — the second form is exact on both sides
// of point 1, where the more usual acos(cos δ13 / cos δxt) loses the sign.

const DEG = Math.PI / 180;
export const R_KM = 6371.0088;          // mean Earth radius (IUGG), km

// The point reached from (lon, lat) after s km at initial azimuth az (degrees).
export function travel(lon, lat, az, s) {
  const p = lat * DEG, l = lon * DEG, a = az * DEG, d = s / R_KM;
  const p2 = Math.asin(Math.sin(p) * Math.cos(d) + Math.cos(p) * Math.sin(d) * Math.cos(a));
  const l2 = l + Math.atan2(Math.sin(a) * Math.sin(d) * Math.cos(p), Math.cos(d) - Math.sin(p) * Math.sin(p2));
  return [((l2 / DEG + 540) % 360) - 180, p2 / DEG];
}

// Angular distance (radians) and initial azimuth (radians) from point 1 to point 2.
function inverse(lon1, lat1, lon2, lat2) {
  const p1 = lat1 * DEG, p2 = lat2 * DEG, dl = (lon2 - lon1) * DEG;
  const h = Math.sin((p2 - p1) / 2) ** 2 + Math.cos(p1) * Math.cos(p2) * Math.sin(dl / 2) ** 2;
  const delta = 2 * Math.asin(Math.min(1, Math.sqrt(h)));
  const theta = Math.atan2(Math.sin(dl) * Math.cos(p2), Math.cos(p1) * Math.sin(p2) - Math.sin(p1) * Math.cos(p2) * Math.cos(dl));
  return { delta, theta };
}

// Where a place lies relative to the section from (lon0, lat0) at azimuth az:
// { along, across } in km (across positive to the right).
export function relativeTo(lon0, lat0, az, lon, lat) {
  const { delta, theta } = inverse(lon0, lat0, lon, lat);
  const dTheta = theta - az * DEG;
  return {
    across: Math.asin(Math.sin(delta) * Math.sin(dTheta)) * R_KM,
    along: Math.atan2(Math.sin(delta) * Math.cos(dTheta), Math.cos(delta)) * R_KM,
  };
}

// A section: its centre, azimuth (degrees) and half-length (km); its ends and
// samples every `step` km: [{ s, lon, lat }] from -half to +half.
export function sectionLine(lon, lat, az, half, step = 5) {
  const out = [];
  for (let s = -half; s <= half + 1e-9; s += step) {
    const [x, y] = s >= 0 ? travel(lon, lat, az, s) : travel(lon, lat, az + 180, -s);
    out.push({ s, lon: x, lat: y });
  }
  return out;
}
