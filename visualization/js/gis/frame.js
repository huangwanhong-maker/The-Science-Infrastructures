// The local frame of a city-scale GIS: east, north and up in metres from an origin on the WGS 84 ellipsoid.
//
// A place's geodetic latitude, longitude and height go to body-fixed Cartesian coordinates (the globe's
// figure.js, the geodesy domain's formulas) and from there to the origin's local frame by the rotation whose
// rows are east, north and up at the origin (the globe's localFrame). The map is exact: no projection's
// distortion enters, and a building 3 km from the origin stands where it stands. Heights above the
// ellipsoid are taken as heights above the ground at the origin, which over a city of a few kilometres
// departs from the tangent plane by d²/2R: 0.7 m at 3 km — the reason the frame's z is "up from the origin's
// plane", not "height above the ground".

import { cartesianDeg, figure, geodetic } from '../globe/figure.js';
import { localFrame } from '../globe/controls.js';

const DEG = Math.PI / 180;

export function localFrameAt(latDeg, lonDeg, heightM = 0, body = 'earth') {
  const fig = figure(body);
  const o = cartesianDeg(fig, lonDeg, latDeg, heightM);
  const R = localFrame(lonDeg * DEG, latDeg * DEG);   // rows: east, north, up
  return {
    origin: { lat: latDeg, lon: lonDeg, height: heightM },
    figure: fig,
    // Geodetic degrees and metres to local east, north, up (m).
    toLocal(lat, lon, h = 0) {
      const p = cartesianDeg(fig, lon, lat, h);
      const d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
      return R.map((row) => row[0] * d[0] + row[1] * d[1] + row[2] * d[2]);
    },
    // Local east, north, up (m) to geodetic degrees and metres.
    toGeodetic(x, y, z = 0) {
      const p = [0, 1, 2].map((k) => o[k] + R[0][k] * x + R[1][k] * y + R[2][k] * z);
      const g = geodetic(fig, p[0], p[1], p[2]);
      return { lat: g.lat / DEG, lon: g.lon / DEG, height: g.h };
    },
  };
}
