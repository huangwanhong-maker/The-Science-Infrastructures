// Volcanic ash: parcels released over an eruption column, carried by the
// forecast wind and falling at the settling speed of their particles until
// they reach the ground (meteorology's settling.js); a Plugin. The wind is
// interpolated linearly in height between 10 m, 850, 500 and 250 hPa, the
// levels standing at their heights in the standard atmosphere, and held at
// its 250 hPa value above; the particles are spheres falling at Stokes' speed
// in the standard atmosphere's air. The vent is the chosen place, or else the
// latest volcano among GDACS's events; the ground is the library's relief.

import { pressureAltitude, stokesSettling, levelBracket, settlingPath } from '../../../domains/meteorology/js/index.js';
import { distanceM } from '../../../domains/geodesy/js/index.js';

const HOUR = 3600000;
// The levels of the wind: need, level, height (m) in the standard atmosphere.
const LEVELS = [['wind10', '10 m', 10], ['wind850', '850 hPa', pressureAltitude(85000)], ['wind500', '500 hPa', pressureAltitude(50000)], ['wind250', '250 hPa', pressureAltitude(25000)]];
// Release heights are coloured from coral (the lowest) to plum (the highest).
const RAMP = ['#f6a07e', '#f6c47e', '#f4dbe7', '#e3b5d3', '#ca8dbe', '#a2629f'];
function colourAt(f) {
  const x = Math.max(0, Math.min(1, f)) * (RAMP.length - 1), k = Math.min(RAMP.length - 2, Math.floor(x)), w = x - k;
  const rgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  const a = rgb(RAMP[k]), b = rgb(RAMP[k + 1]);
  return `#${a.map((v, i) => Math.round(v + (b[i] - v) * w).toString(16).padStart(2, '0')).join('')}`;
}

// The latest volcano among the library's disasters (GDACS), or null.
async function latestVolcano(context) {
  for (const entry of context.list('events').filter((e) => e.kind === 'disaster')) {
    const set = await context.events(entry.name), f = Object.fromEntries(set.fields.map((k, i) => [k, i]));
    const found = set.events.filter((e) => e[f.type] === 'VO');
    if (found.length) {
      const e = found.reduce((a, b) => (b[f.time_ms] > a[f.time_ms] ? b : a));
      return { lon: e[f.lon], lat: e[f.lat], name: e[f.name], time: e[f.time_ms], credit: (set.source && set.source.credit) || entry.title };
    }
  }
  return null;
}

export function create(api) {
  class VolcanicAsh extends api.Plugin {
    id = 'volcanic-ash';
    title = 'Volcanic ash';
    domain = 'meteorology';
    description = 'Parcels of ash released over an eruption column at the chosen place, or at the latest volcano in GDACS’s events, carried by the forecast wind at their height and falling at the settling speed of their particles until they reach the ground. Not a warning service: the Volcanic Ash Advisory Centres issue advisories for aviation, and national agencies warn of ashfall.';
    parameters = {
      top: { label: 'Top of the column (above sea level)', default: 12, min: 1, max: 20, step: 0.5, unit: 'km' },
      count: { label: 'Parcels up the column', default: 12, min: 1, max: 40, step: 1 },
      diameter: { label: 'Particle diameter', default: 25, min: 1, max: 200, step: 1, unit: 'µm' },
      density: { label: 'Particle density', default: 2500, min: 500, max: 3500, step: 50, unit: 'kg/m³' },
      duration: { label: 'Eruption lasting (a release each hour)', default: 0, min: 0, max: 48, step: 1, unit: 'h' },
      hours: { label: 'Hours ahead', default: 72, min: 6, max: 240, step: 6, unit: 'h' },
      step: { label: 'Time step', default: 0.5, min: 0.1, max: 3, step: 0.1, unit: 'h' },
    };

    needs() {
      return LEVELS.map(([name, level]) => ({ name, label: `Wind at ${level}`, standard_names: ['eastward_wind', 'northward_wind'], level, times: 2 }));
    }

    async run(inputs, params, context) {
      // The vent.
      let vent = null, what = '';
      if (context.pinned) {
        vent = { lon: context.pinned[0], lat: context.pinned[1] };
        what = `the chosen place, ${vent.lat.toFixed(2)}°, ${vent.lon.toFixed(2)}°`;
      } else {
        vent = await latestVolcano(context);
        if (!vent) return [api.note('No volcano is among the library’s disasters (GDACS’s events of the month). Choose a vent on the map and run again.')];
        what = `${vent.name} (${vent.credit}, ${api.iso(vent.time).slice(0, 10)})`;
      }
      // The ground: the coarsest relief in the library, or sea level.
      const reliefs = context.fields().filter((m) => !(m.times && m.times.length) && (m.variables || []).some((v) => v.standard_name === 'altitude'));
      let ground = () => 0, groundFrom = 'sea level';
      if (reliefs.length) {
        const meta = reliefs.sort((a, b) => a.grid.nlon * a.grid.nlat - b.grid.nlon * b.grid.nlat)[0];
        const relief = await context.field(meta.name), name = meta.variables.find((v) => v.standard_name === 'altitude').name;
        ground = (lon, lat) => { const z = relief.sample(name, 0, lon, lat); return z > 0 ? z : 0; };
        groundFrom = meta.title || meta.name;
      }
      const base = ground(vent.lon, vent.lat), top = Number(params.top) * 1000;
      if (top <= base) return [api.note(`The top of the column, ${params.top} km, is not above the vent, ${Math.round(base)} m.`)];
      // The wind at a height: linear between the two levels about it.
      const heights = LEVELS.map(([, , z]) => z), winds = LEVELS.map(([name]) => inputs[name]);
      const flow = (lon, lat, z, ms) => {
        const [i, j, w] = levelBracket(z, heights), a = winds[i], b = winds[j];
        const ua = a.product.sample(a.names[0], ms, lon, lat), va = a.product.sample(a.names[1], ms, lon, lat);
        if (!(ua === ua && va === va)) return null;
        if (w === 0) return [ua, va];
        const ub = b.product.sample(b.names[0], ms, lon, lat), vb = b.product.sample(b.names[1], ms, lon, lat);
        return ub === ub && vb === vb ? [ua + (ub - ua) * w, va + (vb - va) * w] : null;
      };
      const d = Number(params.diameter) * 1e-6, rho = Number(params.density);
      const settle = (z) => stokesSettling(d, rho, z).w;
      // The span all four winds cover.
      const first = Math.max(...winds.map((x) => x.product.times[0])), last = Math.min(...winds.map((x) => x.product.times[x.product.times.length - 1]));
      if (!(last > first)) return [api.note('The four winds chosen share no span of time.')];
      const t0 = Math.max(first, Math.min(context.time, last)), dt = Number(params.step) * 3600;
      const releases = Math.max(0, Math.floor(Number(params.duration))) + 1, n = Math.max(1, Math.round(Number(params.count)));
      const paths = [];
      for (let r = 0; r < releases; r++) {
        const t = t0 + r * HOUR;
        if (t >= last) break;
        const steps = Math.max(0, Math.min(Math.round((Number(params.hours) * 3600) / dt), Math.floor((last - t) / 1000 / dt)));
        for (let k = 1; k <= n; k++) {
          const z = base + ((top - base) * k) / n;
          paths.push({ z, t, k, path: settlingPath({ lon: vent.lon, lat: vent.lat, z, t }, flow, settle, { dt, steps, ground }) });
        }
        await this.breathe(context, (r + 1) / releases, `release ${r + 1} of ${releases}`);
      }
      const colour = (k) => colourAt(n === 1 ? 1 : (k - 1) / (n - 1));
      const when = (t) => `${api.iso(Math.round(t / 60000) * 60000).slice(0, 16).replace('T', ' ')} UTC`;
      const landed = paths.filter((p) => p.path.landed);
      const results = [api.sketch({
        title: `${paths.length} parcels of ${params.diameter} µm ash from ${what}, released at ${((base + (top - base) / n) / 1000).toFixed(1)}-${params.top} km above sea level from ${when(t0)}${releases > 1 ? ` for ${releases - 1} h` : ''}: coral at the lowest, plum at the highest`,
        polylines: paths.map((p) => ({ coords: p.path.coords, colour: colour(p.k), width: 1.3, alpha: 0.85 })),
        moving: paths.map((p) => ({ times: p.path.times, coords: p.path.coords, colour: colour(p.k), size: 3 })),
        markers: [{ lon: vent.lon, lat: vent.lat, colour: '#f6a07e', size: 7 }, ...landed.map((p) => ({ lon: p.path.landed.lon, lat: p.path.landed.lat, colour: colour(p.k), size: 2.4 }))],
      })];
      // The first release, height by height.
      const rows = paths.filter((p) => p.t === t0).map((p) => {
        const s = stokesSettling(d, rho, p.z), c = p.path.coords, m = c.length;
        const far = distanceM(vent.lon, vent.lat, c[m - 2], c[m - 1]) / 1000;
        return [`${(p.z / 1000).toFixed(1)} km`, `${(100 * s.w).toFixed(1)} cm/s`, s.reynolds.toFixed(2),
          p.path.landed ? `${((p.path.landed.t - p.t) / HOUR).toFixed(0)} h` : `aloft at ${(p.path.heights[p.path.heights.length - 1] / 1000).toFixed(1)} km`, `${Math.round(far)} km`];
      });
      results.push(api.table(`Parcels released at ${when(t0)}`, ['Released at', 'Settling', 'Reynolds', 'Landed after', 'From the vent'], rows,
        `Stokes' law holds for Reynolds numbers below about 0.1 (within 3 %); at 1 it overstates the speed by some 15 %. Ground: ${groundFrom}.`));
      results.push(api.note(`The wind of ${winds.map((w) => w.model).filter((x, i, all) => all.indexOf(x) === i).join(' and ')} at 10 m and at 850, 500 and 250 hPa, placed at ${heights.map((z) => Math.round(z)).join(', ')} m (the standard atmosphere), linear between, and the 250 hPa wind above. The parcels follow the wind without rising or sinking air, without turbulence to spread them, and without the aggregation that brings fine ash down sooner: the paths show where the wind carries ash of one size, not how much, nor its concentration.`));
      return results;
    }
  }
  return new VolcanicAsh();
}
