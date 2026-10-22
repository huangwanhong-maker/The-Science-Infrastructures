// Relative vorticity and divergence of the wind, by centred differences on
// the sphere: a Diagnostic. The wind may be at any level the library holds;
// the user chooses the source and level in the panel.
//
//   ζ = (1 / (a cos φ)) [ ∂v/∂λ − ∂(u cos φ)/∂φ ],    δ = (1 / (a cos φ)) [ ∂u/∂λ + ∂(v cos φ)/∂φ ]
//
// (Holton and Hakim 2013, §4.1), with a the mean radius of the Earth. Rows of
// the grid run north to south; the poles and, on a grid that does not wrap, the
// edge columns are left without a value.

import { EARTH_RADIUS_M, DEG } from '../../../domains/geodesy/js/index.js';

export function create(api) {
  class Vorticity extends api.Diagnostic {
    id = 'vorticity';
    title = 'Vorticity and divergence';
    description = 'Relative vorticity and divergence of the wind, by centred differences on the sphere.';

    needs() { return [{ name: 'wind', label: 'Wind', standard_names: ['eastward_wind', 'northward_wind'] }]; }

    outputs() {
      return [
        { name: 'zeta', long_name: 'Relative vorticity', standard_name: 'atmosphere_relative_vorticity', units: 's-1' },
        { name: 'delta', long_name: 'Divergence of the wind', standard_name: 'divergence_of_wind', units: 's-1' },
      ];
    }

    compute(inputs, params, index) {
      const { product, names: [un, vn] } = inputs.wind;
      const u = product.slice(un, index), v = product.slice(vn, index);
      const { nlon, nlat, dlon, dlat, lat0, periodic } = product.grid;
      const zeta = new Float32Array(nlon * nlat).fill(NaN), delta = new Float32Array(nlon * nlat).fill(NaN);
      const dx = 2 * dlon * DEG * EARTH_RADIUS_M, dy = 2 * dlat * DEG * EARTH_RADIUS_M;
      for (let r = 1; r < nlat - 1; r++) {
        const lat = lat0 - r * dlat, cos = Math.cos(lat * DEG);
        if (cos < 1e-6) continue;
        const cosN = Math.cos((lat + dlat) * DEG), cosS = Math.cos((lat - dlat) * DEG);
        for (let c = 0; c < nlon; c++) {
          const e = c + 1 < nlon ? c + 1 : periodic ? 0 : -1, w = c > 0 ? c - 1 : periodic ? nlon - 1 : -1;
          if (e < 0 || w < 0) continue;
          const i = r * nlon + c, iN = (r - 1) * nlon + c, iS = (r + 1) * nlon + c;
          const dvdl = (v[r * nlon + e] - v[r * nlon + w]) / dx, dudl = (u[r * nlon + e] - u[r * nlon + w]) / dx;
          const ducos = (u[iN] * cosN - u[iS] * cosS) / dy, dvcos = (v[iN] * cosN - v[iS] * cosS) / dy;
          zeta[i] = (dvdl - ducos) / cos;
          delta[i] = (dudl + dvcos) / cos;
        }
      }
      return { zeta, delta };
    }
  }
  return new Vorticity();
}
