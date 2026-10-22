// Checks of what is the Celestial Atlas's own: the globe on which the observer's
// place is chosen. The science it shows is checked with its domains (metrology,
// astronomy, geodesy) and the shared visualization layer.
// Run from the repository root:  node applications/celestial-atlas/tests/run-tests.mjs

import { suite, seeded } from '../../../development/tools/checks.mjs';
import { orthographic, orthographicInverse } from '../src/render/globe-picker.js';

const { check, done } = suite();

// The globe: the orthographic projection and its inverse (Snyder 1987).
{
  const rnd = seeded(4);
  let err = 0;
  for (let k = 0; k < 2000; k++) {
    const lon0 = (rnd() * 2 - 1) * Math.PI, lat0 = (rnd() - 0.5) * 0.9 * Math.PI;
    const lon = (rnd() * 2 - 1) * Math.PI, lat = Math.asin(rnd() * 2 - 1);
    const p = orthographic(lon, lat, lon0, lat0);
    if (!p.visible) continue;
    const q = orthographicInverse(p.x, p.y, lon0, lat0);
    const dl = Math.abs(((q.lon - lon + 3 * Math.PI) % (2 * Math.PI)) - Math.PI);
    err = Math.max(err, Math.abs(q.lat - lat), Math.cos(lat) * dl);
  }
  check('the globe: orthographic projection and its inverse (Snyder 1987)', err, 1e-9, 'rad');
}

done();
