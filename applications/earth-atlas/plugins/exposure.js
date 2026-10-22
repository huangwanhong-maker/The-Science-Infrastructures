// People exposed to a hazard in the forecast chosen: a Plugin. The hazard is
// heat, the greatest 2 m temperature over the forecast above a threshold (35 °C
// by default), or wind, the greatest 10 m wind speed above one (17 m/s, a gale:
// force 8 on the Beaufort scale is 17.2-20.7 m/s at 10 m). The greatest value
// of each cell is taken over the forecast's times from the time shown to its
// end, sampled bilinearly at the centres of the cells of the population grid
// (GHS-POP at 0.25°; the forecasts are at 1°), and a cell's people are exposed
// when it exceeds the threshold there. The exposed are then summed over the
// countries and the first-order divisions by the geography domain's zonal
// statistics: a cell's people divided among the regions it meets in proportion
// to the areas they share, and those of a coast finer than the regions' cells
// given to the nearest region within 0.3°.
//
// Caveats. The forecasts are instants six hours apart, so that the day's
// greatest temperature, near mid-afternoon, falls between them for some
// longitudes and is underestimated by a degree or more, and gusts are not winds
// of the hour. The people are those of 2025 (GHS-POP R2023A), where they live,
// not where they are at the time. Exposure is not harm: shelter, health and
// acclimatisation decide that.

import { overlapTable, zonalStatistics } from '../../../domains/geography/js/index.js';

const REGIONS = new WeakMap();  // regions products read once, by library: name -> { meta, ids }

export function create(api) {
  class Exposure extends api.Plugin {
    id = 'exposure';
    title = 'People exposed';
    domain = 'geography';
    description = 'People where the forecast’s greatest 2 m temperature, or 10 m wind speed, exceeds a threshold; by country and by province.';
    parameters = {
      hazard: { label: 'Hazard', default: 'heat', choices: [['heat', 'Heat: greatest 2 m temperature'], ['wind', 'Wind: greatest 10 m wind speed']] },
      heat: { label: 'Heat above', default: 35, min: 20, max: 50, step: 0.5, unit: '°C' },
      wind: { label: 'Wind above', default: 17, min: 5, max: 60, step: 0.5, unit: 'm/s' },
    };

    needs() {
      return [
        { name: 'temperature', label: 'Temperature at 2 m', standard_name: 'air_temperature', level: '2 m', times: 2 },
        { name: 'wind', label: 'Wind at 10 m', standard_names: ['eastward_wind', 'northward_wind'], level: '10 m', times: 2 },
        { name: 'people', label: 'People', standard_name: 'population' },
      ];
    }

    // A regions product with its raster of ids, read once.
    async regions(context, name) {
      if (!REGIONS.has(context.library)) REGIONS.set(context.library, new Map());
      const read = REGIONS.get(context.library);
      if (!read.has(name)) {
        read.set(name, (async () => {
          const meta = await context.library.regions(name);
          const ids = new Uint16Array(await (await fetch(`data/products/${meta.file}`)).arrayBuffer());
          return { meta, ids };
        })());
      }
      return read.get(name);
    }

    // The greatest value of each cell over the times from `i0` on: of the temperature, or of the wind's speed.
    greatest(inputs, hazard, i0) {
      const source = hazard === 'heat' ? inputs.temperature : inputs.wind, { product } = source;
      const out = new Float32Array(product.grid.nlon * product.grid.nlat).fill(-Infinity);
      for (let k = i0; k < product.times.length; k++) {
        if (hazard === 'heat') {
          const t = product.slice(source.names[0], k);
          for (let i = 0; i < out.length; i++) if (t[i] > out[i]) out[i] = t[i];
        } else {
          const u = product.slice(source.names[0], k), v = product.slice(source.names[1], k);
          for (let i = 0; i < out.length; i++) { const s = Math.hypot(u[i], v[i]); if (s > out[i]) out[i] = s; }
        }
      }
      for (let i = 0; i < out.length; i++) if (out[i] === -Infinity) out[i] = NaN;
      return { values: out, source };
    }

    async run(inputs, params, context) {
      const hazard = params.hazard === 'wind' ? 'wind' : 'heat';
      const threshold = hazard === 'heat' ? Number(params.heat ?? 35) : Number(params.wind ?? 17);
      const limit = hazard === 'heat' ? threshold + 273.15 : threshold;          // in the stored units, K or m/s
      const started = Date.now();
      // The greatest value over the forecast, from the time shown to its end.
      const forecast = (hazard === 'heat' ? inputs.temperature : inputs.wind).product;
      const i0 = api.indexAtOrBefore(forecast.times, context.time), t0 = forecast.times[i0], t1 = forecast.times[forecast.times.length - 1];
      const { values, source } = this.greatest(inputs, hazard, i0);
      const quantity = hazard === 'heat'
        ? { name: 'tmax2m', long_name: 'Greatest air temperature at 2 m over the forecast', standard_name: 'air_temperature', units: 'K', level: '2 m', cell_methods: 'time: maximum' }
        : { name: 'wspd10m_max', long_name: 'Greatest wind speed at 10 m over the forecast', standard_name: 'wind_speed', units: 'm s-1', level: '10 m', cell_methods: 'time: maximum' };
      const onForecast = api.makeProduct({ name: `${this.productName}-hazard`, title: this.title, grid: forecast.grid, times: [t0], variables: [{ ...quantity, slices: [values] }] });
      await this.breathe(context, 0.2, 'the greatest values over the forecast');
      // On the population grid: the hazard at each cell's centre, and the people where it exceeds the threshold.
      const people = inputs.people.product, grid = people.grid;
      const count = people.slice(inputs.people.names[0], 0);
      const there = api.sampleOnGrid(onForecast, quantity.name, 0, grid);
      const exposed = new Float32Array(count.length);
      let total = 0, everyone = 0;
      for (let i = 0; i < count.length; i++) {
        const n = count[i] === count[i] ? count[i] : 0;
        everyone += n;
        exposed[i] = there[i] > limit ? n : 0;
        total += exposed[i];
      }
      await this.breathe(context, 0.4, 'by country');
      // By country and by province: the exposed, and everyone, for the share.
      const by = {};
      for (const name of ['countries', 'admin1']) {
        const { meta, ids } = await this.regions(context, name);
        const table = overlapTable(grid, meta.grid), options = { count: meta.regions.length, reach: 3, table };
        const hit = zonalStatistics(exposed, grid, ids, meta.grid, options), all = zonalStatistics(count, grid, ids, meta.grid, options);
        const values = {}, rows = [];
        for (const r of meta.regions) {
          if (hit.sum[r.id] >= 0.5) { values[r.code] = hit.sum[r.id]; rows.push({ region: r, exposed: hit.sum[r.id], people: all.sum[r.id] }); }
        }
        by[name] = { meta, values, rows: rows.sort((a, b) => b.exposed - a.exposed) };
        await this.breathe(context, name === 'countries' ? 0.7 : 0.95, 'by province');
      }
      const what = hazard === 'heat' ? `greatest 2 m temperature above ${threshold} °C` : `greatest 10 m wind above ${threshold} m/s`;
      const clause = hazard === 'heat' ? `the greatest 2 m temperature is above ${threshold} °C` : `the greatest 10 m wind is above ${threshold} m/s`;
      const when = (ms) => api.iso(ms).slice(0, 16).replace('T', ' ');
      const window = `${when(t0)} to ${when(t1)} UTC`;
      const made = api.makeProduct({
        name: this.productName, title: `${this.title} (${hazard})`, model: `${this.title} (plug-in)`, domain: this.domain,
        grid, times: [t0],
        variables: [
          { ...quantity, long_name: `${quantity.long_name}, ${window}`, slices: [there] },
          { name: 'exposed', long_name: `People where ${clause}`, standard_name: 'population', units: '1', cell_methods: 'area: sum', slices: [exposed] },
        ],
        source: { dataset: source.meta.name, provider: `computed in the page from ${source.model} and ${inputs.people.model}`, credit: `${this.title} (plug-in), from ${source.model} and ${inputs.people.model}` },
      });
      const fmt = (n) => Math.round(n).toLocaleString('en-GB');
      const top = by.countries.rows.slice(0, 15).map((x) => [x.region.name, fmt(x.exposed), `${((100 * x.exposed) / x.people).toFixed(1)} %`]);
      const seconds = ((Date.now() - started) / 1000).toFixed(1);
      return [
        api.field(made, { surface: 'exposed' }),
        { kind: 'regions', by: 'countries', title: `People exposed by country: ${what}`, units: 'people', values: by.countries.values, colormap: 'matter', scale: 'log' },
        { kind: 'regions', by: 'admin1', title: `People exposed by province: ${what}`, units: 'people', values: by.admin1.values, colormap: 'matter', scale: 'log' },
        api.table(`The ${Math.min(15, top.length)} countries with the most people exposed`, ['Country', 'People exposed', 'Of its people'], top,
          `${source.model}, ${window}; people of ${inputs.people.meta.epoch || 'the population product'} from ${inputs.people.model}.`),
        api.note(`${fmt(total)} people (${((100 * total) / everyone).toFixed(2)} % of ${fmt(everyone)}) live where ${clause}, in ${by.countries.rows.length} countries and ${by.admin1.rows.length} provinces, from ${window} (${seconds} s). Six-hourly instants miss part of the day’s peak; exposure is not harm.`),
      ];
    }
  }
  return new Exposure();
}
