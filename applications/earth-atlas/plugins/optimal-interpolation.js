// Optimal interpolation of a background field toward observations: a
// DataAssimilation plug-in. The observations are real ones when an observation
// set is chosen (weather stations, buoys, or our own sensors: each station's
// report nearest the background's time); otherwise the values of a second model
// at a lattice of points, which shows where two models disagree. With
// innovations d = y − H x_b,
//
//   x_a = x_b + B Hᵀ (H B Hᵀ + R)⁻¹ d,
//
// B a Gaussian of great-circle distance with scale L and variance σ_b², and
// R = σ_o² I (Kalnay 2003, §5.4). The system is solved locally (Lorenc 1981):
// at each grid point, with the nearest observations within 4 L, at most a few
// dozen, by Cholesky. Observations whose innovation exceeds four times its
// expected spread, √(σ_b² + σ_o²), are rejected as gross errors.
//
// Caveat. No correction is made for the height of a station above the model's
// ground, and a sensor's errors are taken to be those of a weather station.

import { toVector, EARTH_RADIUS_M } from '../../../domains/geodesy/js/index.js';
import { bilinear } from '../../../foundations/physis/js/index.js';

const HOUR = 3600000;

export function create(api) {
  class OptimalInterpolation extends api.DataAssimilation {
    id = 'optimal-interpolation';
    title = 'Optimal interpolation';
    description = 'A background analysed toward observations: weather stations, buoys or your own sensors, or else a second model at a lattice of points. The increment is shown, and the analysis joins the models of the quantity.';
    parameters = {
      count: { label: 'Lattice points (second model)', default: 300, min: 10, max: 1500, step: 10 },
      length: { label: 'Correlation length L', default: 500, min: 50, max: 3000, step: 50, unit: 'km' },
      sigmaB: { label: 'Background error σb', default: 2, min: 0.1, max: 20, step: 0.1 },
      sigmaO: { label: 'Observation error σo', default: 1, min: 0.01, max: 20, step: 0.1 },
      window: { label: 'Reports within', default: 1.5, min: 0.25, max: 12, step: 0.25, unit: 'h' },
      nearest: { label: 'Observations per point', default: 24, min: 1, max: 64, step: 1 },
    };

    needs() {
      return [
        { name: 'background', label: 'Background', standard_name: 'air_temperature', level: '2 m' },
        { name: 'stations', label: 'Observed by', kind: 'observations', variable: 'air_temperature', optional: true },
        { name: 'truth', label: 'Or by a model', standard_name: 'air_temperature', level: '2 m', distinct: true, optional: true },
      ];
    }

    // Real observations near the analysis's time t, else the second model at a lattice of points.
    async observations(inputs, t, params) {
      if (inputs.stations) {
        const set = inputs.stations.observations, ix = Object.fromEntries(set.fields.map((f, i) => [f, i]));
        const col = ix[inputs.stations.variable], best = new Map();
        for (const row of set.rows) {
          const value = row[col];
          if (value === null || !Number.isFinite(value)) continue;
          const dt = Math.abs(Date.parse(row[ix.time]) - t);
          if (dt > (params.window ?? this.parameters.window.default) * HOUR) continue;
          const prior = best.get(row[ix.id]);
          if (!prior || dt < prior.dt) best.set(row[ix.id], { lon: row[ix.lon], lat: row[ix.lat], value, dt, id: row[ix.id] });
        }
        this.source = `${best.size} reports of ${inputs.stations.model} within ${params.window ?? this.parameters.window.default} h`;
        return [...best.values()];
      }
      if (!inputs.truth) throw new Error('choose observations, or a second model to observe from');
      const { product, names: [name] } = inputs.truth;
      this.source = `${inputs.truth.model} at ${Math.round(params.count)} lattice points`;
      return api.lattice(Math.round(params.count))
        .map((p) => ({ ...p, value: product.sample(name, t, p.lon, p.lat) }))
        .filter((o) => o.value === o.value);
    }

    async analyse(background, grid, all, params, context) {
      const value = (k) => params[k] ?? this.parameters[k].default;      // a parameter, or its default when a caller omits it
      const L = value('length') * 1000, sb2 = value('sigmaB') ** 2, so2 = value('sigmaO') ** 2, K = Math.max(1, Math.round(value('nearest')));
      const reach = Math.cos(Math.min(Math.PI, (4 * L) / EARTH_RADIUS_M));       // cos of the angle beyond which B is taken as 0
      const covariance = (c) => sb2 * Math.exp(-((EARTH_RADIUS_M * Math.acos(Math.min(1, Math.max(-1, c)))) ** 2) / (2 * L * L));
      // Innovations, and the gross-error check.
      const limit = 4 * Math.sqrt(sb2 + so2);
      const obs = [];
      let rejected = 0;
      for (const o of all) {
        const d = o.value - bilinear(background, grid, o.lon, o.lat);
        if (!(d === d)) continue;
        if (Math.abs(d) > limit) { rejected++; continue; }
        obs.push({ ...o, d, p: toVector(o.lon, o.lat) });
      }
      // Observations binned by 5° of latitude and longitude, to find the near ones fast.
      const bins = new Map(), key = (lon, lat) => `${Math.floor((lat + 90) / 5)},${Math.floor((((lon + 180) % 360) + 360) % 360 / 5)}`;
      for (const o of obs) { const k = key(o.lon, o.lat); if (!bins.has(k)) bins.set(k, []); bins.get(k).push(o); }
      const span = Math.ceil(((4 * L) / EARTH_RADIUS_M) * (180 / Math.PI) / 5) + 1;
      const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
      const out = new Float32Array(background.length), { nlon, nlat, lon0, dlon, lat0, dlat } = grid;
      const A = new Float64Array(K * K), d = new Float64Array(K);
      for (let r = 0; r < nlat; r++) {
        const lat = lat0 - r * dlat, bi = Math.floor((lat + 90) / 5);
        for (let c = 0; c < nlon; c++) {
          const i = r * nlon + c, b = background[i];
          if (b !== b) { out[i] = NaN; continue; }
          const lon = lon0 + c * dlon, p = toVector(lon, lat), bj = Math.floor((((lon + 180) % 360) + 360) % 360 / 5);
          const near = [];
          for (let di = -span; di <= span; di++) {
            const ii = bi + di;
            if (ii < 0 || ii >= 36) continue;
            // Near the poles every longitude is near.
            const width = Math.abs(lat) > 60 ? 36 : span;
            for (let dj = -width; dj <= width; dj++) {
              const list = bins.get(`${ii},${(((bj + dj) % 72) + 72) % 72}`);
              if (!list) continue;
              for (const o of list) { const cs = dot(p, o.p); if (cs > reach) near.push([cs, o]); }
            }
            if (Math.abs(lat) > 60) break;
          }
          if (!near.length) { out[i] = b; continue; }
          near.sort((x, y) => y[0] - x[0]);
          const n = Math.min(K, near.length), used = near.slice(0, n).map(([, o]) => o);
          for (let a = 0; a < n; a++) {
            for (let e = 0; e <= a; e++) A[a * n + e] = A[e * n + a] = covariance(dot(used[a].p, used[e].p));
            A[a * n + a] += so2;
            d[a] = used[a].d;
          }
          const w = choleskySolve(A, d, n);
          let s = 0;
          for (let a = 0; a < n; a++) s += covariance(near[a][0]) * w[a];
          out[i] = b + s;
        }
        if (r % 15 === 14) await this.breathe(context, r / nlat, `row ${r + 1} of ${nlat}`);
      }
      // The diagnostics of an analysis: observation minus background, and minus analysis.
      const rms = (xs) => Math.sqrt(xs.reduce((a, x) => a + x * x, 0) / Math.max(1, xs.length));
      const oa = obs.map((o) => o.value - bilinear(out, grid, o.lon, o.lat)).filter((x) => x === x);
      this.stats = { used: obs.length, rejected, omb: rms(obs.map((o) => o.d)), oma: rms(oa) };
      return out;
    }

    async run(inputs, params, context) {
      const results = await super.run(inputs, params, context);
      const s = this.stats, unit = inputs.background.variables[0].units === 'K' ? 'K' : inputs.background.variables[0].units;
      results.push(api.table(`Optimal interpolation of ${inputs.background.model}: ${this.source}`, ['', 'Value'], [
        ['Observations used', String(s.used)],
        ['Rejected (gross error)', String(s.rejected)],
        ['RMS observation − background', `${s.omb.toFixed(2)} ${unit}`],
        ['RMS observation − analysis', `${s.oma.toFixed(2)} ${unit}`],
      ], 'The analysis draws close to the observations where they are dense, and stays with the background elsewhere.'));
      return results;
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
