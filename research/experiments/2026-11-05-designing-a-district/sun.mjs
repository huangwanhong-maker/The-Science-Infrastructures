// The Sun's altitude and azimuth at a place and at given instants, from the astronomy domain's subsolar point
// (domains/astronomy/js/sun.js), as the Urban Science Workbench's clock takes it (plugins/landscape-dynamics/shared.js,
// sunByHour). Low-precision theory: 0.01° in right ascension and declination (the module's own statement).
//
//   node sun.mjs '{"lat": -37.81, "lon": 144.96, "instants": ["2026-06-21T02:15:00Z", ...]}'
//
// prints [{"instant", "altitude_deg", "azimuth_deg"}, ...] (azimuth from north through east).

import { subsolarPoint, sunAltAz } from '../../../domains/astronomy/js/index.js';

const request = JSON.parse(process.argv[2]);
const out = request.instants.map((iso) => {
  const s = sunAltAz(request.lon, request.lat, subsolarPoint(Date.parse(iso)));
  return { instant: iso, altitude_deg: s.altitude, azimuth_deg: s.azimuth };
});
process.stdout.write(JSON.stringify(out));
