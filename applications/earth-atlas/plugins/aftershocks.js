// Aftershocks: the probability of aftershocks of magnitude 5, 6 and 7 or more
// over the next day, week and month, by the model of Reasenberg and Jones
// (1989; geophysics's aftershocks.js), with their generic Californian
// parameters unless others are set; a Plugin. The mainshock is the largest
// earthquake of the library's month within a radius of the chosen place, or
// the largest of the month. The zone of aftershocks is drawn as a circle whose
// radius is the subsurface rupture length of Wells and Coppersmith (1994), as
// the rupture may run that far from the epicentre in either direction; the
// earthquakes of the catalogue within it since the mainshock are compared
// with the number the model expects.

import { GENERIC_CALIFORNIA, expectedAftershocks, aftershockProbability, subsurfaceRuptureLengthKm } from '../../../domains/geophysics/js/index.js';
import { distanceM, move } from '../../../domains/geodesy/js/index.js';

const DAY = 86400000;

export function create(api) {
  class Aftershocks extends api.Plugin {
    id = 'aftershocks';
    title = 'Aftershocks';
    domain = 'geophysics';
    description = 'The chance of aftershocks of magnitude 5, 6 and 7 or more over the next day, week and month after the largest recent earthquake near the chosen place, by the model of Reasenberg and Jones (1989), and the zone in which they fall. Not a warning service: aftershock forecasts are issued by the USGS and national geological surveys.';
    parameters = {
      radius: { label: 'Mainshock within', default: 300, min: 10, max: 3000, step: 10, unit: 'km' },
      a: { label: 'Productivity a', default: GENERIC_CALIFORNIA.a, min: -4, max: 0, step: 0.01 },
      b: { label: 'Gutenberg–Richter b', default: GENERIC_CALIFORNIA.b, min: 0.5, max: 1.5, step: 0.01 },
      p: { label: 'Omori decay p', default: GENERIC_CALIFORNIA.p, min: 0.5, max: 2, step: 0.01 },
      c: { label: 'Omori delay c', default: GENERIC_CALIFORNIA.c, min: 0.001, max: 1, step: 0.001, unit: 'day' },
    };

    needs() { return []; }

    async run(inputs, params, context) {
      const sets = context.list('events').filter((e) => e.kind === 'earthquake');
      const entry = sets.find((e) => /month/.test(e.name)) || sets[0];
      if (!entry) return [api.note('No earthquakes are in the library.')];
      const set = await context.events(entry.name), f = Object.fromEntries(set.fields.map((k, i) => [k, i]));
      const quakes = set.events.map((e) => ({ t: e[f.time_ms], lon: e[f.lon], lat: e[f.lat], depth: e[f.depth_km], mag: e[f.mag], type: e[f.mag_type], place: e[f.place] }));
      if (!quakes.length) return [api.note(`${entry.title} holds no earthquakes.`)];
      const largest = (list) => list.reduce((a, b) => (b.mag > a.mag ? b : a));
      const model = { a: Number(params.a), b: Number(params.b), p: Number(params.p), c: Number(params.c) };
      // The mainshock.
      let main = null, how = '';
      if (context.pinned) {
        const [lon, lat] = context.pinned, near = quakes.filter((q) => distanceM(lon, lat, q.lon, q.lat) <= Number(params.radius) * 1000);
        if (near.length) { main = largest(near); how = `the largest within ${params.radius} km of the chosen place`; }
      }
      if (!main) { main = largest(quakes); how = context.pinned ? `none within ${params.radius} km of the chosen place, so the largest of the month` : 'the largest of the month'; }
      // From the instant shown, or from the mainshock if the instant precedes it.
      const since = Math.max(0, (context.time - main.t) / DAY);
      const windows = [['Next day', 1], ['Next week', 7], ['Next month', 30]];
      const percent = (x) => (x >= 0.995 ? '> 99 %' : x >= 0.1 ? `${Math.round(100 * x)} %` : x >= 0.001 ? `${(100 * x).toPrecision(2)} %` : '< 0.1 %');
      const rows = [5, 6, 7].map((M) => [`M ≥ ${M}`, ...windows.map(([, days]) => percent(aftershockProbability(since, since + days, M, main.mag, model)))]);
      // The zone, and the catalogue's earthquakes within it since the mainshock.
      const L = subsurfaceRuptureLengthKm(main.mag), end = Date.parse(set.window && set.window.end) || Math.max(...quakes.map((q) => q.t));
      const within = quakes.filter((q) => q !== main && q.t > main.t && distanceM(main.lon, main.lat, q.lon, q.lat) <= L * 1000);
      const least = set.selection && set.selection.min_magnitude !== undefined ? set.selection.min_magnitude : Math.min(...quakes.map((q) => q.mag));
      const expected = expectedAftershocks(0, Math.max(0, (end - main.t) / DAY), least, main.mag, model);
      const ring = [];
      for (let k = 0; k <= 96; k++) { const az = (2 * Math.PI * k) / 96; ring.push(...move(main.lon, main.lat, L * 1000 * Math.sin(az), L * 1000 * Math.cos(az))); }
      const named = `M ${main.mag} (${main.type}), ${api.iso(main.t).slice(0, 16).replace('T', ' ')} UTC, ${main.depth} km deep, ${main.place}`;
      const moment = /^mw/i.test(main.type || '');
      return [
        api.table(`Chance of at least one aftershock, after ${named}`, ['Magnitude', ...windows.map(([w]) => w)], rows,
          `Mainshock: ${how}; from ${since > 0 ? `the instant shown, ${since.toFixed(1)} days after it` : 'the mainshock'}. Reasenberg and Jones (1989), a = ${model.a}, b = ${model.b}, p = ${model.p}, c = ${model.c} day${model.a === GENERIC_CALIFORNIA.a && model.b === GENERIC_CALIFORNIA.b && model.p === GENERIC_CALIFORNIA.p && model.c === GENERIC_CALIFORNIA.c ? ' (generic California)' : ''}.`),
        api.sketch({
          title: `Aftershock zone of ${named}: ${Math.round(L)} km (Wells and Coppersmith 1994), and the catalogue's earthquakes in it since`,
          polylines: [{ coords: ring, colour: '#f6c47e', width: 1.4, alpha: 0.9, dashed: true }],
          markers: [{ lon: main.lon, lat: main.lat, colour: '#f6a07e', size: 7 }, ...within.map((q) => ({ lon: q.lon, lat: q.lat, colour: '#e8c3dc', size: 2 + 0.6 * (q.mag - least) }))],
        }),
        api.note(`Observed: ${within.length} earthquakes of M ${least} or more within ${Math.round(L)} km since the mainshock, to ${api.iso(end).slice(0, 10)}; the model expects ${expected.toFixed(1)}. The catalogue holds every earthquake there, not only aftershocks, and misses small ones in the first hours.`),
        api.note(`The generic parameters were fitted to California's sequences; elsewhere, and from one sequence to the next, the rate differs by up to ten times (Page et al. 2016), and the parameters set here may be a region's.${moment ? '' : ` The mainshock's magnitude is ${main.type}, not moment magnitude, for which the model and the rupture length were fitted.`} For M at or above the mainshock's, the chance is that of a larger earthquake to come, which the model counts among the aftershocks.`),
      ];
    }
  }
  return new Aftershocks();
}
