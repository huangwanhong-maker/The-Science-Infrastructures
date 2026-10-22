// Heat over the days of a forecast: the daily maxima of the wet-bulb
// temperature (Stull 2011) and of the heat index (the NWS's), and the excess
// of the day's greatest temperature over the month's normal daily maximum
// (TerraClimate, 1991-2020, when the library holds it). A Diagnostic whose
// times are days: each place's local solar day, from midnight to midnight of
// mean solar time, with the maxima taken over the forecast's samples within it.

import { wetBulbStull, withinStullRange, heatIndex, sampledMaximum, forecastLocalDays, localDay } from '../../../domains/meteorology/js/index.js';

const HOUR = 3600000;
const sameGrid = (a, b) => a.nlon === b.nlon && a.nlat === b.nlat && a.lon0 === b.lon0 && a.lat0 === b.lat0 && a.dlon === b.dlon && a.dlat === b.dlat;

export function create(api) {
  // A need's variable at each of `times`, on `grid`; null where it has no such time.
  const series = (input, name, times, grid) => times.map((t) => {
    const i = input.product.times.indexOf(t);
    if (i < 0) return null;
    return sameGrid(input.product.grid, grid) ? input.product.slice(name, i) : api.sampleOnGrid(input.product, name, i, grid);
  });
  // The library's monthly normals carrying a variable that passes `test`: { product, name, meta }, or null.
  async function normals(context, test) {
    const meta = context.fields().find((m) => m.climatology && (m.variables || []).some(test));
    return meta ? { product: await context.field(meta.name), name: meta.variables.find(test).name, meta } : null;
  }

  class Heat extends api.Diagnostic {
    id = 'heat';
    title = 'Heat';
    domain = 'meteorology';
    description = 'The daily maxima of the wet-bulb temperature and the heat index over the forecast days, and how far the day’s greatest temperature exceeds the month’s normal daily maximum. Not a warning service: heat warnings are issued by national weather services.';

    needs() {
      return [
        { name: 'temperature', label: 'Temperature at 2 m', standard_name: 'air_temperature', level: '2 m', times: 8 },
        { name: 'humidity', label: 'Relative humidity at 2 m', standard_name: 'relative_humidity', level: '2 m', times: 8 },
      ];
    }

    outputs() {
      return [
        { name: 'wet_bulb_max', long_name: 'Daily maximum wet-bulb temperature (Stull 2011)', standard_name: 'wet_bulb_temperature', units: 'K', level: '2 m', cell_methods: 'time: maximum (local day)' },
        { name: 'heat_index_max', long_name: 'Daily maximum heat index (NWS)', standard_name: 'heat_index_of_air_temperature', units: 'K', level: '2 m', cell_methods: 'time: maximum (local day)' },
        { name: 'tmax_excess', long_name: 'Daily maximum temperature above the normal', units: 'K', level: '2 m', difference: true, cell_methods: 'time: maximum (local day)' },
      ];
    }

    async run(inputs, params, context) {
      const climatology = Object.values(inputs).find((x) => x.meta.climatology);
      if (climatology) return [api.note(`${climatology.model} is a climatology of monthly normals, not a forecast: choose a forecast for every source.`)];
      const product = inputs.temperature.product, grid = product.grid, times = product.times, n = grid.nlon * grid.nlat;
      // The days both sources cover.
      const days = Object.values(inputs).map((x) => forecastLocalDays(x.product.times)).reduce((a, b) => a.filter((d) => b.includes(d)));
      if (!days.length) return [api.note('The sources share no whole local day at every longitude.')];
      const T = series(inputs.temperature, inputs.temperature.names[0], times, grid), H = series(inputs.humidity, inputs.humidity.names[0], times, grid);
      const normal = await normals(context, (v) => v.standard_name === 'air_temperature' && /maximum within days/.test(v.cell_methods || ''));
      const outputs = this.outputs(), lists = Object.fromEntries(outputs.map((o) => [o.name, []]));
      for (const [di, d] of days.entries()) {
        const wet = new Float32Array(n).fill(NaN), hi = new Float32Array(n).fill(NaN), excess = new Float32Array(n).fill(NaN);
        // The month's normal daily maximum (°C), on the forecast's grid; the month of the day at 0° E.
        const tmax = normal ? api.sampleOnGrid(normal.product, normal.name, normal.product.bracket(d + 12 * HOUR).i0, grid) : null;
        for (let c = 0; c < grid.nlon; c++) {
          const [t0, t1] = localDay(d, grid.lon0 + c * grid.dlon);
          for (let r = 0; r < grid.nlat; r++) {
            const i = r * grid.nlon + c;
            const at = (k) => (T[k] && H[k] ? [T[k][i] - 273.15, Math.min(H[k][i], 100)] : [NaN, NaN]);
            // Stull's fit where it is valid; at saturation (RH above 99 %) as at 99 %, within 0.2 °C.
            const w = sampledMaximum(times, (k) => {
              const [t, rh] = at(k);
              if (!(t === t && rh === rh)) return NaN;
              return withinStullRange(t, Math.min(rh, 99)) ? wetBulbStull(t, Math.min(rh, 99)) : -Infinity;
            }, t0, t1);
            wet[i] = w > -Infinity ? w + 273.15 : NaN;
            hi[i] = sampledMaximum(times, (k) => (T[k] && H[k] ? heatIndex(T[k][i], Math.min(H[k][i], 100)) : NaN), t0, t1);
            if (tmax) excess[i] = sampledMaximum(times, (k) => (T[k] ? T[k][i] : NaN), t0, t1) - 273.15 - tmax[i];
          }
        }
        lists.wet_bulb_max.push(wet); lists.heat_index_max.push(hi); lists.tmax_excess.push(excess);
        await this.breathe(context, (di + 1) / days.length, `day ${di + 1} of ${days.length}`);
      }
      const shown = normal ? outputs : outputs.filter((o) => o.name !== 'tmax_excess');
      const dayTimes = days.map((d) => d + 12 * HOUR);
      const made = api.makeProduct({
        name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain, grid, times: dayTimes,
        variables: shown.map((o) => ({ ...o, long_name: o.name === 'tmax_excess' ? `${o.long_name} (${normal.meta.title || normal.meta.model})` : o.long_name, slices: lists[o.name] })),
        source: { dataset: inputs.temperature.meta.name, provider: `computed in the page from ${inputs.temperature.model}`, credit: `${this.title} (plug-in), from ${inputs.temperature.model}` },
      });
      const spacing = (times[1] - times[0]) / HOUR;
      const results = [
        api.field(made, { surface: 'wet_bulb_max' }),
        api.note(`${days.length} local days, ${api.iso(days[0]).slice(0, 10)} to ${api.iso(days[days.length - 1]).slice(0, 10)}. The maxima are of the forecast's samples every ${spacing} h, which miss the afternoon peak between them: the greatest temperature by 1-2 °C${normal ? ', and so the excess over the normal, which is a true daily maximum, by as much' : ''}. The wet-bulb temperature is Stull's fit at sea-level pressure, within 1 °C between -20 and 50 °C and 5 and 99 % (left out beyond); a wet-bulb temperature of 35 °C is the bound of human heat loss in theory (Sherwood and Huber 2010), and healthy people reach their limits near 31 °C. The heat index is the NWS's, for a person in the shade; in the sun it is up to 8 °C higher.`),
      ];
      if (!normal) results.push(api.note('No monthly normals of the daily maximum temperature are in the library, so the excess over the normal is not computed.'));
      else results.push(api.note(`The normal, from ${normal.meta.title || normal.meta.model}, is read at the forecast's grid points from its own finer grid: in mountains the two grids stand at different heights, and the excess carries their difference (about 6.5 °C per km).`));
      if (context.pinned) {
        const [lon, lat] = context.pinned, c = (x) => (x === x ? (x - 273.15).toFixed(1) : '—');
        const rows = dayTimes.map((t) => [api.iso(t).slice(0, 10), c(made.sample('wet_bulb_max', t, lon, lat)), c(made.sample('heat_index_max', t, lon, lat)),
          normal ? (() => { const x = made.sample('tmax_excess', t, lon, lat); return x === x ? `${x > 0 ? '+' : ''}${x.toFixed(1)}` : '—'; })() : '—']);
        results.push(api.table(`At ${lat.toFixed(2)}°, ${lon.toFixed(2)}° (°C)`, ['Day', 'Wet bulb', 'Heat index', 'Above normal'], rows, 'Daily maxima of the local solar day, interpolated between grid points.'));
      }
      return results;
    }
  }
  return new Heat();
}
