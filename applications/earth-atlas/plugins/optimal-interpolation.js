// Optimal interpolation of a background field toward a network of
// observations: a DataAssimilation plug-in. The observations are the values of
// a second model at the points of a lattice over the globe, so the analysis
// also shows where two models disagree. With innovations d = y − H x_b,
//
//   x_a = x_b + B Hᵀ (H B Hᵀ + R)⁻¹ d,
//
// B a Gaussian of great-circle distance with scale L and variance σ_b², and
// R = σ_o² I (Kalnay 2003, §5.4). The observation-space system is solved by
// Cholesky; the increment at a grid point sums the observations within 4 L.

import { toVector, EARTH_RADIUS_M } from '../../../domains/geodesy/js/index.js';
import { bilinear } from '../../../foundations/physis/js/index.js';

export function create(api) {
  class OptimalInterpolation extends api.DataAssimilation {
    id = 'optimal-interpolation';
    title = 'Optimal interpolation';
    description = 'A background analysed toward observations taken from a second model at a lattice of points; the increment is shown, and the analysis joins the models of the quantity.';
    parameters = {
      count: { label: 'Observations', default: 300, min: 10, max: 1500, step: 10 },
      length: { label: 'Correlation length L', default: 500, min: 50, max: 3000, step: 50, unit: 'km' },
      sigmaB: { label: 'Background error σb', default: 2, min: 0.1, max: 20, step: 0.1 },
      sigmaO: { label: 'Observation error σo', default: 1, min: 0.01, max: 20, step: 0.1 },
    };

    needs() {
      return [
        { name: 'background', label: 'Background', standard_name: 'air_temperature', level: '2 m' },
        { name: 'truth', label: 'Observed from', standard_name: 'air_temperature', level: '2 m', distinct: true },
      ];
    }

    // The second model's values at a lattice of points, at the background's time.
    async observations(inputs, index, params) {
      const { product, names: [name] } = inputs.truth, t = inputs.background.product.times[index];
      return api.lattice(Math.round(params.count))
        .map((p) => ({ ...p, value: product.sample(name, t, p.lon, p.lat) }))
        .filter((o) => o.value === o.value);
    }

    async analyse(background, grid, obs, params, context) {
      const n = obs.length, L = params.length * 1000, sb2 = params.sigmaB ** 2, so2 = params.sigmaO ** 2;
      const reach = Math.cos(Math.min(Math.PI, (4 * L) / EARTH_RADIUS_M));          // cos of the angle beyond which B is taken as 0
      const covariance = (cosAngle) => sb2 * Math.exp(-((EARTH_RADIUS_M * Math.acos(Math.min(1, Math.max(-1, cosAngle)))) ** 2) / (2 * L * L));
      const points = obs.map((o) => toVector(o.lon, o.lat));
      const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
      // The observation-space system (H B Hᵀ + R) w = d.
      const A = new Float64Array(n * n), d = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j <= i; j++) A[i * n + j] = A[j * n + i] = covariance(dot(points[i], points[j]));
        A[i * n + i] += so2;
        d[i] = obs[i].value - bilinear(background, grid, obs[i].lon, obs[i].lat);
      }
      const w = choleskySolve(A, d, n);
      // The analysis: the background plus B Hᵀ w, over the observations within reach.
      const out = new Float32Array(background.length), { nlon, nlat, lon0, dlon, lat0, dlat } = grid;
      for (let r = 0; r < nlat; r++) {
        const lat = lat0 - r * dlat;
        for (let c = 0; c < nlon; c++) {
          const i = r * nlon + c, b = background[i];
          if (b !== b) { out[i] = NaN; continue; }
          const p = toVector(lon0 + c * dlon, lat);
          let s = 0;
          for (let k = 0; k < n; k++) { const cosAngle = dot(p, points[k]); if (cosAngle > reach) s += covariance(cosAngle) * w[k]; }
          out[i] = b + s;
        }
        if (r % 20 === 19) await this.breathe(context, r / nlat, `row ${r + 1} of ${nlat}`);
      }
      return out;
    }
  }
  return new OptimalInterpolation();
}

// Solves A x = b for a symmetric positive-definite A (n × n, row-major) by
// Cholesky factorisation. (A linear solver belongs to calx; this is its stand-in.)
function choleskySolve(A, b, n) {
  const L = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i * n + j];
      for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = i === j ? Math.sqrt(Math.max(s, 1e-12)) : s / L[j * n + j];
    }
  }
  const y = new Float64Array(n), x = new Float64Array(n);
  for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k]; y[i] = s / L[i * n + i]; }
  for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k]; x[i] = s / L[i * n + i]; }
  return x;
}
