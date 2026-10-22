// The barotropic vorticity equation, run through the ForecastModel class: the
// model the atlas also runs in a worker (ADR 0015), here as a plug-in, stepped
// in the page from the 500 hPa wind of the model chosen. The domain's model
// does the dynamics (domains/meteorology, spectral T42 or less); the plug-in
// puts the wind on the model's Gaussian grid, records the balanced height and
// the vorticity on a regular grid, and the base class verifies the height
// against the source's later analyses, beside persistence.

import { BarotropicModel } from '../../../domains/meteorology/js/index.js';
import { bilinear } from '../../../foundations/physis/js/index.js';

const OUTPUT = { lon0: -178.75, dlon: 2.5, nlon: 144, lat0: 88.75, dlat: 2.5, nlat: 72, periodic: true };

export function create(api) {
  class Barotropic extends api.ForecastModel {
    id = 'barotropic';
    title = 'Barotropic model';
    description = 'The barotropic vorticity equation from the 500 hPa wind, the height carried in balance, verified against the source.';
    dt = 1800;
    every = 6;
    parameters = {
      hours: { label: 'Forecast length', default: 48, min: 6, max: 168, step: 6, unit: 'h' },
      truncation: { label: 'Truncation', default: 42, choices: [[21, 'T21'], [31, 'T31'], [42, 'T42']] },
    };

    needs() {
      return [
        { name: 'wind', label: 'Wind at 500 hPa', standard_names: ['eastward_wind', 'northward_wind'], level: '500 hPa', times: 2 },
        { name: 'height', label: 'Height at 500 hPa', standard_name: 'geopotential_height', level: '500 hPa', times: 2 },
      ];
    }

    grid() { return OUTPUT; }

    outputs() {
      return [
        { name: 'z@500hPa', long_name: 'Geopotential height at 500 hPa', standard_name: 'geopotential_height', units: 'm', level: '500 hPa' },
        { name: 'zeta@500hPa', long_name: 'Relative vorticity at 500 hPa', standard_name: 'atmosphere_relative_vorticity', units: 's-1', level: '500 hPa' },
      ];
    }

    verify() { return { output: 'z@500hPa', need: 'height', label: 'height', unit: 'm' }; }
    show() { return { surface: 'z@500hPa', lines: 'z@500hPa' }; }

    // The wind on the model's Gaussian grid (rows south to north), the analysed
    // height on the output grid, and the balanced height at the start: the
    // forecast height is the analysis plus the change of the balanced height.
    async initial(inputs, index, params) {
      const model = new BarotropicModel({ truncation: Number(params.truncation) || 42 });
      const { nlat, nlon, mu } = model.sh;
      const { product, names: [un, vn] } = inputs.wind;
      const su = product.slice(un, index), sv = product.slice(vn, index);
      const u = new Float64Array(nlat * nlon), v = new Float64Array(nlat * nlon);
      for (let j = 0; j < nlat; j++) {
        const lat = (Math.asin(mu[j]) * 180) / Math.PI;
        for (let k = 0; k < nlon; k++) {
          const lon = (((360 * k) / nlon + 180) % 360) - 180;
          u[j * nlon + k] = bilinear(su, product.grid, lon, lat) || 0;
          v[j * nlon + k] = bilinear(sv, product.grid, lon, lat) || 0;
        }
      }
      model.setWind(u, v);
      const height = inputs.height, at = api.indexAtOrBefore(height.product.times, product.times[index]);
      const zStart = api.sampleOnGrid(height.product, height.names[0], at, OUTPUT);
      const z0 = Float64Array.from(model.sh.synthesiseOn(model.balancedHeight(model.streamFunction()), OUTPUT));
      return { model, zStart, z0 };
    }

    advance(state, seconds) { state.model.step(seconds); }

    record(state) {
      const { model, zStart, z0 } = state;
      const balanced = model.sh.synthesiseOn(model.balancedHeight(model.streamFunction()), OUTPUT);
      const z = new Float32Array(balanced.length);
      for (let i = 0; i < z.length; i++) z[i] = zStart[i] + balanced[i] - z0[i];
      return { 'z@500hPa': z, 'zeta@500hPa': Float32Array.from(model.sh.synthesiseOn(model.zeta, OUTPUT)) };
    }
  }
  return new Barotropic();
}
