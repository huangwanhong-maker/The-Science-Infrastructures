// Heavy rain in a forecast: the precipitation of the 24 and 72 hours before
// each forecast time, accumulated from the precipitation rate (meteorology's
// accumulated()), and each as a fraction of the month's normal precipitation
// (TerraClimate, 1991-2020, when the library holds it). A Diagnostic at the
// forecast's own times, from the first whose 24 hours the forecast covers.

import { accumulated } from '../../../domains/meteorology/js/index.js';

const HOUR = 3600000;

export function create(api) {
  // The library's monthly normals carrying a variable that passes `test`: { product, name, meta }, or null.
  async function normals(context, test) {
    const meta = context.fields().find((m) => m.climatology && (m.variables || []).some(test));
    return meta ? { product: await context.field(meta.name), name: meta.variables.find(test).name, meta } : null;
  }

  class HeavyRain extends api.Diagnostic {
    id = 'heavy-rain';
    title = 'Heavy rain';
    domain = 'meteorology';
    description = 'The rain of the last 24 and 72 hours at each forecast time, from the precipitation rate, and each as a fraction of the month’s normal precipitation. The fraction compares an amount with an ordinary month; it is not a return period. Not a warning service: flood and rain warnings are issued by national weather and hydrological services.';
    parameters = {
      least: { label: 'Fraction only where the month’s normal is at least', default: 1, min: 0, max: 50, step: 0.5, unit: 'mm' },
    };

    needs() { return [{ name: 'rain', label: 'Precipitation', standard_name: 'precipitation_flux', times: 5 }]; }

    outputs() {
      return [
        { name: 'rain_24h', long_name: 'Precipitation in the last 24 hours', standard_name: 'precipitation_amount', units: 'kg m-2', cell_methods: 'time: sum (24 h)' },
        { name: 'rain_72h', long_name: 'Precipitation in the last 72 hours', standard_name: 'precipitation_amount', units: 'kg m-2', cell_methods: 'time: sum (72 h)' },
        { name: 'rain_24h_ratio', long_name: 'Precipitation in 24 hours, as a fraction of the month’s normal', units: '1', cell_methods: 'time: sum (24 h)' },
        { name: 'rain_72h_ratio', long_name: 'Precipitation in 72 hours, as a fraction of the month’s normal', units: '1', cell_methods: 'time: sum (72 h)' },
      ];
    }

    async run(inputs, params, context) {
      const { product, names: [name] } = inputs.rain, grid = product.grid, all = product.times, n = grid.nlon * grid.nlat;
      const spacing = all[1] - all[0], begins = all[0] - spacing;            // the first rate is the mean over the interval before it
      const times = all.filter((t) => t - 24 * HOUR >= begins);
      if (!times.length) return [api.note('The forecast is too short for 24 hours of rain.')];
      const rates = all.map((_, k) => product.slice(name, k));
      const normal = await normals(context, (v) => v.standard_name === 'precipitation_amount');
      const least = Number(params.least ?? 1);
      const outputs = normal ? this.outputs() : this.outputs().filter((o) => !o.name.endsWith('ratio'));
      const lists = Object.fromEntries(outputs.map((o) => [o.name, []]));
      let month = -1, ppt = null;
      for (const [ti, t] of times.entries()) {
        const r24 = new Float32Array(n), r72 = new Float32Array(n).fill(NaN), whole72 = t - 72 * HOUR >= begins;
        for (let i = 0; i < n; i++) {
          r24[i] = accumulated(all, (k) => rates[k][i], t - 24 * HOUR, t).amount;
          if (whole72) r72[i] = accumulated(all, (k) => rates[k][i], t - 72 * HOUR, t).amount;
        }
        lists.rain_24h.push(r24); lists.rain_72h.push(r72);
        if (normal) {
          // The normal of the month in which the window ends.
          const m = normal.product.bracket(t).i0;
          if (m !== month) { month = m; ppt = api.sampleOnGrid(normal.product, normal.name, m, grid); }
          const f24 = new Float32Array(n), f72 = new Float32Array(n);
          for (let i = 0; i < n; i++) {
            const ok = ppt[i] >= least && ppt[i] > 0;
            f24[i] = ok ? r24[i] / ppt[i] : NaN;
            f72[i] = ok ? r72[i] / ppt[i] : NaN;
          }
          lists.rain_24h_ratio.push(f24); lists.rain_72h_ratio.push(f72);
        }
        await this.breathe(context, (ti + 1) / times.length, `time ${ti + 1} of ${times.length}`);
      }
      const made = api.makeProduct({
        name: this.productName, title: this.title, model: `${this.title} (plug-in)`, domain: this.domain, grid, times,
        variables: outputs.map((o) => ({ ...o, slices: lists[o.name] })),
        source: { dataset: inputs.rain.meta.name, provider: `computed in the page from ${inputs.rain.model}`, credit: `${this.title} (plug-in), from ${inputs.rain.model}` },
      });
      const results = [
        api.field(made, { surface: 'rain_24h' }),
        api.note(`From ${inputs.rain.model}'s precipitation rate every ${spacing / HOUR} h, each value taken as the mean over the ${spacing / HOUR} hours before it: a sum of samples, not the model's own accumulation, so a shower between samples may be missed or counted for the whole interval.`),
      ];
      if (normal) {
        results.push(api.note(`What the fraction says, and what it does not. It compares the forecast amount with the precipitation of an ordinary month at the place (${normal.meta.title || normal.meta.model}): 0.5 in 24 hours is half a normal month's rain in a day, which is unusual in most climates. It does not say how rare the amount is: rarity needs the distribution of the heaviest days at the place, a climatology of extremes from which return periods are drawn, and the library holds none. The same fraction is ordinary where rain comes in a few heavy storms (a dry season's edge, a desert) and exceptional where it falls steadily; and an amount of small fraction can still flood steep or saturated ground, or a river fed from far away. Where the normal is below ${least} mm, the fraction is left out.`));
      } else {
        results.push(api.note('No monthly normals of precipitation are in the library, so the fractions are not computed.'));
      }
      if (context.pinned) {
        const [lon, lat] = context.pinned, fmt = (x, d = 0) => (x === x ? x.toFixed(d) : '—');
        const rows = times.filter((t) => new Date(t).getUTCHours() % 12 === 0).map((t) => [api.iso(t).slice(0, 13).replace('T', ' ') + ' UTC',
          fmt(made.sample('rain_24h', t, lon, lat), 1), fmt(made.sample('rain_72h', t, lon, lat), 1),
          normal ? fmt(made.sample('rain_24h_ratio', t, lon, lat), 2) : '—', normal ? fmt(made.sample('rain_72h_ratio', t, lon, lat), 2) : '—']);
        results.push(api.table(`At ${lat.toFixed(2)}°, ${lon.toFixed(2)}°`, ['Ending', '24 h (mm)', '72 h (mm)', '24 h / month', '72 h / month'], rows, 'Interpolated between grid points.'));
      }
      return results;
    }
  }
  return new HeavyRain();
}
