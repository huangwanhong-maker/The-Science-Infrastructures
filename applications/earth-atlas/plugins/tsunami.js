// Tsunami travel times: the first arrival of the long wave, c = √(g h), from
// a source over the sea of the library's relief, by a shortest path on the
// grid (oceanography's tsunami.js); a Plugin. The source is the chosen place,
// or else the latest shallow earthquake (depth below 70 km) of the magnitude
// set or more in the library's earthquakes. The result is a field of hours with
// hourly contours, the source marked, and the arrivals near coastal places.

import { TravelTimes } from '../../../domains/oceanography/js/index.js';
import { distanceM, EARTH_RADIUS_M, DEG } from '../../../domains/geodesy/js/index.js';

const HOUR = 3600000;

// Depths (m, positive in the sea; NaN on land) on a grid coarser by an integer
// factor f, from the mean altitude over each block of f × f cells.
function depthsOn(z, grid, f) {
  if (f <= 1) return { depth: Float32Array.from(z, (x) => (x < 0 ? -x : NaN)), grid };
  const nlon = Math.floor(grid.nlon / f), nlat = Math.floor(grid.nlat / f), depth = new Float32Array(nlon * nlat);
  for (let r = 0; r < nlat; r++) {
    for (let c = 0; c < nlon; c++) {
      let sum = 0, count = 0;
      for (let i = 0; i < f; i++) for (let j = 0; j < f; j++) { const x = z[(r * f + i) * grid.nlon + c * f + j]; if (x === x) { sum += x; count++; } }
      depth[r * nlon + c] = count && sum / count < 0 ? -sum / count : NaN;
    }
  }
  return { depth, grid: { ...grid, lon0: grid.lon0 + ((f - 1) / 2) * grid.dlon, dlon: grid.dlon * f, nlon, lat0: grid.lat0 - ((f - 1) / 2) * grid.dlat, dlat: grid.dlat * f, nlat,
    periodic: Boolean(grid.periodic) && nlon * f === grid.nlon } };
}

// The latest earthquake of magnitude ≥ least and depth < 70 km in the library's
// earthquakes (the month's, if there is such a set), or null.
async function latestShallow(context, least) {
  const sets = context.list('events').filter((e) => e.kind === 'earthquake');
  const entry = sets.find((e) => /month/.test(e.name)) || sets[0];
  if (!entry) return null;
  const set = await context.events(entry.name), f = Object.fromEntries(set.fields.map((k, i) => [k, i]));
  const found = set.events.filter((e) => e[f.mag] >= least && e[f.depth_km] < 70);
  if (!found.length) return null;
  const e = found.reduce((a, b) => (b[f.time_ms] > a[f.time_ms] ? b : a));
  return { lon: e[f.lon], lat: e[f.lat], mag: e[f.mag], magType: e[f.mag_type], depthKm: e[f.depth_km], place: e[f.place], time: e[f.time_ms], credit: (set.source && set.source.credit) || entry.title };
}

export function create(api) {
  class Tsunami extends api.Plugin {
    id = 'tsunami';
    title = 'Tsunami travel times';
    domain = 'oceanography';
    description = 'The first arrival of a tsunami’s long wave from a source over the sea, in hours with hourly contours: from the chosen place, or else from the latest shallow earthquake of the magnitude set or more. It says when a wave could first arrive, not whether one was made, nor its height. Not a warning service: tsunami warnings are issued by the Pacific Tsunami Warning Center, the regional centres and national agencies.';
    parameters = {
      spacing: { label: 'Grid', default: 10, choices: [[10, '10′ (18 km)'], [20, '20′ (37 km)'], [30, '30′ (56 km)']] },
      stencil: { label: 'Stencil', default: 32, choices: [[32, '32 neighbours, within 2.1 %'], [16, '16 neighbours, within 4.2 %']] },
      magnitude: { label: 'Earthquakes of magnitude at least', default: 6.5, min: 5, max: 9.5, step: 0.1 },
    };

    needs() { return [{ name: 'relief', label: 'Bathymetry', standard_name: 'altitude' }]; }

    async run(inputs, params, context) {
      const started = performance.now();
      const { product, names: [name] } = inputs.relief;
      // The source.
      let source = null, what = '';
      if (context.pinned) {
        source = { lon: context.pinned[0], lat: context.pinned[1] };
        what = `the chosen place, ${source.lat.toFixed(2)}°, ${source.lon.toFixed(2)}°`;
      } else {
        const quake = await latestShallow(context, Number(params.magnitude));
        if (!quake) return [api.note(`No earthquake of magnitude ${params.magnitude} or more, shallower than 70 km, is among the library's earthquakes. Choose a place in the sea on the map and run again.`)];
        source = quake;
        what = `the earthquake of ${api.iso(quake.time)}, M ${quake.mag} (${quake.magType}), ${quake.depthKm} km deep, ${quake.place} (${quake.credit})`;
      }
      // The sea on the grid chosen.
      const factor = Math.max(1, Math.round(Number(params.spacing) / (product.grid.dlon * 60)));
      const { depth, grid } = depthsOn(product.slice(name, 0), product.grid, factor);
      let wet = 0;
      for (const x of depth) if (x > 0) wet++;
      const solver = new TravelTimes(depth, grid, [source.lon, source.lat], { neighbours: Number(params.stencil) });
      if (!solver.source) return [api.note(`The source, ${what}, is more than 100 km from the sea of ${inputs.relief.model} on this grid.`)];
      while (!solver.advance(150000)) await this.breathe(context, solver.settled / wet, `${Math.round((100 * solver.settled) / wet)} % of the sea`);
      const seconds = solver.result(), hours = Float32Array.from(seconds, (s) => s / 3600);
      const elapsed = (performance.now() - started) / 1000;
      const made = api.makeProduct({
        name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain, grid, times: [],
        variables: [{ name: 'tsunami_travel_time', long_name: 'Tsunami travel time, first arrival of the long wave', units: 'h', contours: { step: 1, base: 0 }, slices: [hours] }],
        source: { dataset: inputs.relief.meta.name, provider: `computed in the page from ${inputs.relief.model}`, credit: `${this.title} (plug-in), over ${inputs.relief.model}` },
      });
      const s = solver.source;
      const markers = [{ lon: source.lon, lat: source.lat, colour: '#f6a07e', size: 7 }];
      if (s.movedKm > 1) markers.push({ lon: s.lon, lat: s.lat, colour: '#f6c47e', size: 3.5 });
      const results = [
        api.field(made, { surface: 'tsunami_travel_time', lines: 'tsunami_travel_time' }),
        api.sketch({ title: `Source: ${what}`, markers, labels: [{ lon: source.lon, lat: source.lat, text: source.mag ? `M ${source.mag}` : 'source' }] }),
        api.note(`From ${what}${s.movedKm > 1 ? `, started at the nearest sea, ${s.movedKm.toFixed(0)} km away, ${Math.round(s.depth)} m deep` : ''}. Over ${inputs.relief.model} on a grid of ${Math.round(grid.dlon * 60)}′ (${grid.nlon} × ${grid.nlat}, ${wet.toLocaleString('en-GB')} points of sea), ${params.stencil} neighbours: ${elapsed.toFixed(1)} s.`),
        api.note('The long wave travels at √(g h) over the depth h, so the times are those of linear long waves from a point. A great earthquake breaks hundreds of kilometres of fault, and the coast facing the rupture is reached sooner; observed arrivals across an ocean come up to a few per cent later (the water is compressible and the Earth elastic); near the coast the times depend on how finely the depth is known.'),
      ];
      // First arrivals at the sea nearest to coastal places (within 40 km).
      const geography = await (context.geography ? context.geography() : null);
      const places = (geography && geography.meta && geography.meta.layers && geography.meta.layers.places) || [];
      if (places.length) {
        const rowsOf = Math.ceil(40000 / (EARTH_RADIUS_M * grid.dlat * DEG)), arrivals = [];
        for (const p of places) {
          if (p.class === 'Scientific station') continue;
          const r0 = Math.round((grid.lat0 - p.lat) / grid.dlat), c0 = Math.round((p.lon - grid.lon0) / grid.dlon);
          const cols = Math.min(grid.nlon >> 1, Math.ceil(rowsOf / Math.max(0.05, Math.cos(p.lat * DEG))));
          let best = Infinity;
          for (let r = Math.max(0, r0 - rowsOf); r <= Math.min(grid.nlat - 1, r0 + rowsOf); r++) {
            for (let c = c0 - cols; c <= c0 + cols; c++) {
              const cc = ((c % grid.nlon) + grid.nlon) % grid.nlon, h = hours[r * grid.nlon + cc];
              if (h === h && h < best && distanceM(p.lon, p.lat, grid.lon0 + cc * grid.dlon, grid.lat0 - r * grid.dlat) <= 40000) best = h;
            }
          }
          if (best < Infinity) arrivals.push({ p, h: best });
        }
        arrivals.sort((a, b) => a.h - b.h);
        const clock = (h) => `${Math.floor(h)} h ${String(Math.round((h % 1) * 60)).padStart(2, '0')} min`;
        results.push(api.table('First arrivals near coastal places', ['Place', 'After the source', 'At (UTC)', 'Distance'],
          arrivals.slice(0, 15).map(({ p, h }) => [`${p.names.en || p.names.default}${p.country ? `, ${p.country}` : ''}`, clock(h),
            source.time ? api.iso(source.time + h * HOUR).slice(0, 16).replace('T', ' ') : '—', `${Math.round(distanceM(source.lon, source.lat, p.lon, p.lat) / 1000)} km`]),
          `The sea within 40 km of each place, the earliest of its points; ${arrivals.length} places reached.`));
      }
      return results;
    }
  }
  return new Tsunami();
}
