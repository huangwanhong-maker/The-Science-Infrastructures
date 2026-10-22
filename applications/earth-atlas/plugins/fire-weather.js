// Fire weather over the days of a forecast: the Canadian Forest Fire Weather
// Index System (Van Wagner 1987; meteorology's fireweather.js), a Diagnostic
// whose times are days. At each place and day it takes the forecast time
// nearest local noon (mean solar time, within half the forecast's spacing),
// the temperature, humidity and wind there, and the rain of the 24 hours
// before, accumulated from the precipitation rate; and it carries the three
// moisture codes from day to day, from starting values on the first day.
// Over the sea, judged by the relief in the library, nothing is computed.

import { fireWeatherDay, STANDARD_START, localNoon, nearestIndex, accumulated, forecastNoonDays } from '../../../domains/meteorology/js/index.js';

const HOUR = 3600000, DAY = 24 * HOUR;
const sameGrid = (a, b) => a.nlon === b.nlon && a.nlat === b.nlat && a.lon0 === b.lon0 && a.lat0 === b.lat0 && a.dlon === b.dlon && a.dlat === b.dlat;

export function create(api) {
  // A need's variable at each of `times`, on `grid` (sampled when its own grid differs); null where it has no such time.
  const series = (input, name, times, grid) => times.map((t) => {
    const i = input.product.times.indexOf(t);
    if (i < 0) return null;
    return sameGrid(input.product.grid, grid) ? input.product.slice(name, i) : api.sampleOnGrid(input.product, name, i, grid);
  });
  // The sea on `grid`, from the coarsest relief in the library (altitude below 0 m), or null.
  async function seaMask(context, grid) {
    const reliefs = context.fields().filter((m) => !(m.times && m.times.length) && (m.variables || []).some((v) => v.standard_name === 'altitude'));
    if (!reliefs.length) return null;
    const meta = reliefs.sort((a, b) => a.grid.nlon * a.grid.nlat - b.grid.nlon * b.grid.nlat)[0];
    const relief = await context.field(meta.name), z = api.sampleOnGrid(relief, meta.variables.find((v) => v.standard_name === 'altitude').name, 0, grid);
    return Uint8Array.from(z, (x) => (x < 0 ? 1 : 0));
  }

  class FireWeather extends api.Diagnostic {
    id = 'fire-weather';
    title = 'Fire weather';
    domain = 'meteorology';
    description = 'The Canadian Fire Weather Index System over the forecast days: the fine fuel, duff and drought codes from the weather at local noon and the rain of the day before, and the spread, buildup and fire weather indices from them. The codes start from the values set here on the first day. Not a warning service: fire danger is rated officially by national and regional agencies.';
    parameters = {
      ffmc: { label: 'Fine Fuel Moisture Code on the first day', default: STANDARD_START.ffmc, min: 0, max: 101, step: 1 },
      dmc: { label: 'Duff Moisture Code on the first day', default: STANDARD_START.dmc, min: 0, max: 400, step: 1 },
      dc: { label: 'Drought Code on the first day', default: STANDARD_START.dc, min: 0, max: 1200, step: 5 },
    };

    needs() {
      return [
        { name: 'temperature', label: 'Temperature at 2 m', standard_name: 'air_temperature', level: '2 m', times: 8 },
        { name: 'humidity', label: 'Relative humidity at 2 m', standard_name: 'relative_humidity', level: '2 m', times: 8 },
        { name: 'wind', label: 'Wind at 10 m', standard_names: ['eastward_wind', 'northward_wind'], level: '10 m', times: 8 },
        { name: 'rain', label: 'Precipitation', standard_name: 'precipitation_flux', times: 8 },
      ];
    }

    outputs() {
      const at = { cell_methods: 'time: point (local noon)', units: '1' };
      return [
        { name: 'fire_weather_index', long_name: 'Fire Weather Index', ...at },
        { name: 'initial_spread_index', long_name: 'Initial Spread Index', ...at },
        { name: 'buildup_index', long_name: 'Buildup Index', ...at },
        { name: 'fine_fuel_moisture_code', long_name: 'Fine Fuel Moisture Code', ...at },
        { name: 'duff_moisture_code', long_name: 'Duff Moisture Code', ...at },
        { name: 'drought_code', long_name: 'Drought Code', ...at },
      ];
    }

    async run(inputs, params, context) {
      const climatology = Object.values(inputs).find((x) => x.meta.climatology);
      if (climatology) return [api.note(`${climatology.model} is a climatology of monthly normals, not a forecast: choose a forecast for every source.`)];
      const product = inputs.temperature.product, grid = product.grid, times = product.times, n = grid.nlon * grid.nlat;
      // The days every source covers.
      const days = Object.values(inputs).map((x) => forecastNoonDays(x.product.times)).reduce((a, b) => a.filter((d) => b.includes(d)));
      if (!days.length) return [api.note('The sources share no whole day of fire weather: it needs local noon at every longitude and the 24 hours of rain before it.')];
      const T = series(inputs.temperature, inputs.temperature.names[0], times, grid), H = series(inputs.humidity, inputs.humidity.names[0], times, grid);
      const U = series(inputs.wind, inputs.wind.names[0], times, grid), V = series(inputs.wind, inputs.wind.names[1], times, grid);
      const P = series(inputs.rain, inputs.rain.names[0], times, grid);
      const sea = await seaMask(context, grid);
      const spacing = times[1] - times[0];
      const codes = { ffmc: new Float64Array(n).fill(Number(params.ffmc)), dmc: new Float64Array(n).fill(Number(params.dmc)), dc: new Float64Array(n).fill(Number(params.dc)) };
      const outputs = this.outputs(), lists = Object.fromEntries(outputs.map((o) => [o.name, []]));
      const keys = { fire_weather_index: 'fwi', initial_spread_index: 'isi', buildup_index: 'bui', fine_fuel_moisture_code: 'ffmc', duff_moisture_code: 'dmc', drought_code: 'dc' };
      for (const [di, d] of days.entries()) {
        const month = new Date(d).getUTCMonth() + 1;
        const out = Object.fromEntries(outputs.map((o) => [o.name, new Float32Array(n).fill(NaN)]));
        for (let c = 0; c < grid.nlon; c++) {
          const noon = localNoon(d, grid.lon0 + c * grid.dlon), k = nearestIndex(times, noon, spacing / 2);
          if (k < 0 || !T[k] || !H[k] || !U[k] || !V[k]) continue;
          for (let r = 0; r < grid.nlat; r++) {
            const i = r * grid.nlon + c;
            if (sea && sea[i]) continue;
            const rain = accumulated(times, (j) => (P[j] ? P[j][i] : NaN), noon - DAY, noon).amount;
            const weather = { T: T[k][i] - 273.15, H: Math.min(100, Math.max(0, H[k][i])), W: 3.6 * Math.hypot(U[k][i], V[k][i]), rain };
            if (!(weather.T === weather.T && weather.H === weather.H && weather.W === weather.W && rain === rain)) continue;
            const day = fireWeatherDay({ ffmc: codes.ffmc[i], dmc: codes.dmc[i], dc: codes.dc[i] }, weather, month, grid.lat0 - r * grid.dlat);
            codes.ffmc[i] = day.ffmc; codes.dmc[i] = day.dmc; codes.dc[i] = day.dc;
            for (const o of outputs) out[o.name][i] = day[keys[o.name]];
          }
        }
        for (const o of outputs) lists[o.name].push(out[o.name]);
        await this.breathe(context, (di + 1) / days.length, `day ${di + 1} of ${days.length}`);
      }
      const dayTimes = days.map((d) => d + 12 * HOUR);
      const made = api.makeProduct({
        name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain, grid, times: dayTimes,
        variables: outputs.map((o) => ({ ...o, slices: lists[o.name] })),
        source: { dataset: inputs.temperature.meta.name, provider: `computed in the page from ${inputs.temperature.model}`, credit: `${this.title} (plug-in), from ${inputs.temperature.model}` },
      });
      const first = api.iso(days[0]).slice(0, 10);
      const results = [
        api.field(made, { surface: 'fire_weather_index' }),
        api.note(`${days.length} days, ${first} to ${api.iso(days[days.length - 1]).slice(0, 10)}, each at the forecast time nearest local noon; the rain is that of the 24 hours before, from the precipitation rate${sea ? '; the sea is left out' : ''}.`),
        api.note(`The codes start on ${first} from FFMC ${params.ffmc}, DMC ${params.dmc} and DC ${params.dc} (the standard spring start is 85, 6 and 15). The fine fuel forgets its start within two or three days; the duff and drought codes carry a difference in their start unchanged through dry days, and lose it only as rain wets them (20 mm halves it in the duff, and takes a tenth of it from the deep layers). So the DMC, the DC, and with them the BUI and the FWI, describe the forecast days and not the season's dryness, unless the starting codes are an agency's current values.`),
      ];
      if (context.pinned) {
        const [lon, lat] = context.pinned;
        const rows = dayTimes.map((t) => [api.iso(t).slice(0, 10), ...['fine_fuel_moisture_code', 'duff_moisture_code', 'drought_code', 'initial_spread_index', 'buildup_index', 'fire_weather_index']
          .map((name) => { const v = made.sample(name, t, lon, lat); return v === v ? v.toFixed(1) : '—'; })]);
        results.push(api.table(`At ${lat.toFixed(2)}°, ${lon.toFixed(2)}°`, ['Day', 'FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI'], rows,
          'Interpolated between grid points. Classes of danger depend on the region and its fuels.'));
      }
      return results;
    }
  }
  return new FireWeather();
}
