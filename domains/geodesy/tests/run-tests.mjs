// Checks of the geodesy domain: the Earth as a sphere, latitudes on the WGS 84
// ellipsoid, distances, and directions in the local horizon.
// Run from the repository root:  node domains/geodesy/tests/run-tests.mjs

import { suite, seeded } from '../../../development/tools/checks.mjs';
import {
  DEG, EARTH_RADIUS_M, toVector, fromVector, eastNorth, centralAngle, distanceM, move, wrapLongitude,
  geocentricLatitude, geodeticLatitude, distanceKm, nearestPlace, compassPoint, meanDirection,
} from '../js/index.js';

const { check, assert, done } = suite();
const random = seeded(20261022);

// ---- the sphere ------------------------------------------------------------
{
  let worst = 0;
  for (let i = 0; i < 1000; i++) {
    const lon = random() * 360 - 180, lat = Math.asin(2 * random() - 1) / DEG;
    const [l2, p2] = fromVector(toVector(lon, lat));
    worst = Math.max(worst, Math.abs(wrapLongitude(l2 - lon)) * Math.cos(lat * DEG), Math.abs(p2 - lat));
  }
  check('sphere: vector round trip', worst, 1e-10, 'deg');
  const [e, n] = eastNorth(35, 40);
  const p = toVector(35, 40);
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  check('sphere: east, north and up are orthonormal', Math.max(Math.abs(dot(e, n)), Math.abs(dot(e, p)), Math.abs(dot(n, p)), Math.abs(dot(e, e) - 1), Math.abs(dot(n, n) - 1)), 1e-12);
  // One degree of a great circle is R pi / 180 = 111.195 km.
  check('sphere: one degree of latitude', distanceM(10, 0, 10, 1) - (EARTH_RADIUS_M * Math.PI) / 180, 1e-6, 'm');
  // Paris to Tokyo, 9 713 km on the sphere of mean radius (independently: 0.152 rad x R).
  const angle = centralAngle(2.3522, 48.8566, 139.6917, 35.6895);
  check('sphere: haversine equals the spherical law of cosines', angle - Math.acos(dot(toVector(2.3522, 48.8566), toVector(139.6917, 35.6895))), 1e-12, 'rad');
  // Moving along the equator by a quarter of the circumference.
  const q = move(0, 0, (Math.PI / 2) * EARTH_RADIUS_M, 0);
  check('sphere: a quarter turn east along the equator', Math.hypot(q[0] - 90, q[1]), 1e-9, 'deg');
  // Moving north through the pole comes down the opposite meridian.
  const r = move(30, 89, 0, (2 * Math.PI / 180) * EARTH_RADIUS_M);
  check('sphere: through the north pole', Math.hypot(wrapLongitude(r[0] - (-150)), r[1] - 89), 1e-9, 'deg');
}

// ---- the ellipsoid, and where an observer stands ------------------------------
{
  let worst = 0, maxDiff = 0, atLat = 0;
  for (let d = -89.5; d <= 89.5; d += 0.5) {
    const phi = d * DEG, psi = geocentricLatitude(phi);
    worst = Math.max(worst, Math.abs(geodeticLatitude(psi) - phi));
    if (Math.abs(phi - psi) > maxDiff) { maxDiff = Math.abs(phi - psi); atLat = d; }
  }
  check('geodetic latitude from geocentric and back', worst / DEG, 1e-12, 'deg');
  // The two latitudes differ most at 45 degrees, by 11.5 arcminutes (Torge and Müller 2012, sec. 4.1).
  check('geodetic minus geocentric latitude at most 0.1924 deg, at 45 deg', Math.abs(maxDiff / DEG - 0.1924) + Math.abs(Math.abs(atLat) - 45) / 100, 0.0006, 'deg');
  // London to New York along the sphere: 5570 km (a standard example of the haversine formula).
  check('distance London - New York', distanceKm(51.5074 * DEG, -0.1278 * DEG, 40.7128 * DEG, -74.006 * DEG) - 5570, 3, 'km');
  const near = nearestPlace([{ lat: 35.68, lon: 139.77 }, { lat: 48.86, lon: 2.35 }], 35.0 * DEG, 139.0 * DEG);
  // From 35 N, 139 E to Tokyo: 0.68 degrees of latitude and 0.77 of longitude, 103.0 km.
  check('the nearest place is Tokyo, 103.0 km away', Math.abs(near.place.lon - 139.77) + Math.abs(near.km - 103.0), 0.1, 'km');
}

// ---- directions in the local horizon ------------------------------------------
{
  assert('compass points', compassPoint(225) === 'SW' && compassPoint(359) === 'N' && compassPoint(100) === 'E');
  check('mean of 350 and 10 degrees is 0', Math.min(meanDirection([350, 10]), 360 - meanDirection([350, 10])), 1e-9, 'deg');
}

done();
