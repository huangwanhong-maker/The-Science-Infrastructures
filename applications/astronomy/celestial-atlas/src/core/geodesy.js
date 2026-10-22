// Places on the Earth: geodetic and geocentric latitude on the WGS 84
// ellipsoid, and distances along the sphere.
//
// Provisional stand-in for core/domains/geodesy (laboratory ADR 0021).
//
// For a point on the ellipsoid's surface, the geocentric latitude psi (the
// angle of the radius vector) and the geodetic latitude phi (the angle of the
// normal) satisfy tan(psi) = (1 - e^2) tan(phi) (Torge and Müller 2012, eq. 4.11);
// they differ by up to 0.19 degrees, at 45 degrees of latitude.
//
// Distances use the haversine formula on a sphere of the Earth's mean radius,
// R1 = 6371.0088 km (Moritz 1980), within 0.5 % of the ellipsoidal distance.

export const WGS84_E2 = 0.00669437999014;       // first eccentricity squared (NIMA TR8350.2)
export const MEAN_RADIUS_KM = 6371.0088;

export const geocentricLatitude = (phi) => Math.atan((1 - WGS84_E2) * Math.tan(phi));
export const geodeticLatitude = (psi) => Math.atan(Math.tan(psi) / (1 - WGS84_E2));

// Great-circle distance in km between two places given in radians.
export function distanceKm(lat1, lon1, lat2, lon2) {
  const s = Math.sin((lat2 - lat1) / 2) ** 2 + Math.cos(lat1) * Math.cos(lat2) * Math.sin((lon2 - lon1) / 2) ** 2;
  return 2 * MEAN_RADIUS_KM * Math.asin(Math.min(1, Math.sqrt(s)));
}

// The nearest of a list of places ({ lat, lon } in degrees) to a place in
// radians, with its distance in km.
export function nearestPlace(places, lat, lon) {
  let best = null, bestKm = Infinity;
  for (const p of places) {
    const d = distanceKm(lat, lon, p.lat * Math.PI / 180, p.lon * Math.PI / 180);
    if (d < bestKm) { best = p; bestKm = d; }
  }
  return best ? { place: best, km: bestKm } : null;
}
