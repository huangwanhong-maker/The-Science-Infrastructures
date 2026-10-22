// The shallow-water equations on the sphere, by the spectral transform method:
// the laboratory's second model, one step beyond the barotropic vorticity
// equation (src/core/barotropic.js).
//
// Provisional stand-in for core/domains/meteorology (laboratory ADR 0019).
//
// Hypotheses: a single layer of fluid, hydrostatic and incompressible, of
// geopotential Phi = g h, on a rotating sphere, without friction or heating.
// Its flow may diverge, so it carries gravity waves as well as Rossby waves,
// and its Rossby waves feel a finite radius of deformation sqrt(Phi) / f,
// which slows the ultra-long waves that the barotropic model moves westward too
// fast. In vorticity-divergence form (Hack and Jakob 1992, eqs. 2.10-2.12):
//
//   d zeta / dt  = - div( eta V ),                              eta = zeta + f   (1)
//   d delta / dt =   k . curl( eta V ) - lap( Phi + |V|^2 / 2 )                  (2)
//   d Phi / dt   = - div( Phi V )                                                (3)
//
// with V = k x grad(psi) + grad(chi), zeta = lap(psi), delta = lap(chi).
// The model holds zeta, delta and Phi as spherical-harmonic coefficients,
// evaluates the products on the Gaussian grid, steps by the classical
// fourth-order Runge-Kutta method, and damps the smallest scales by lap^2, as
// the barotropic model does.
//
// Started from an analysis, the model is not in balance, and gravity waves of
// a few hours' period ring through the forecast. A digital filter
// (Lynch and Huang 1992) removes them: the model is run backward for half a
// span without damping, then forward for the whole span, and the states are
// averaged with the weights of a low-pass filter (a Lanczos-windowed sinc)
// whose cut-off period is a few hours. The average, valid at the start, keeps
// the slow, balanced part of the state.
//
// Sources:
//   Hack, J. J. and Jakob, R. (1992), Description of a global shallow water
//   model based on the spectral transform method, NCAR/TN-343+STR.
//   Williamson, D. L., et al. (1992), A standard test set for numerical
//   approximations to the shallow water equations in spherical geometry,
//   J. Comput. Phys. 102, 211-224: test case 2 (tests/run-tests.mjs).
//   Lynch, P. and Huang, X.-Y. (1992), Initialization of the HIRLAM model using
//   a digital filter, Mon. Wea. Rev. 120, 1019-1034.

import { SphericalHarmonics, EARTH_RADIUS, OMEGA } from './spectral.js';
import { GRAVITY } from './barotropic.js';

export class ShallowWaterModel {
  // tilt: the angle (radians) by which the axis of rotation leans from the
  // grid's pole toward longitude 180 degrees; 0 for the Earth, other values to
  // test flow across the grid's poles (Williamson et al. 1992, test case 2).
  constructor({ truncation = 42, radius = EARTH_RADIUS, omega = OMEGA, dampingHours = 8, tilt = 0 } = {}) {
    this.sh = new SphericalHarmonics(truncation);
    this.radius = radius;
    this.omega = omega;
    this.damping = dampingHours > 0 ? 1 / (dampingHours * 3600) : 0;
    const S = this.sh.size;
    this.state = { zeta: this.sh.zeros(), delta: this.sh.zeros(), phi: this.sh.zeros() };
    this.time = 0;
    const { nlat, nlon, mu } = this.sh;
    this.f = new Float64Array(nlat * nlon);
    this.cos2 = new Float64Array(nlat * nlon);
    for (let j = 0; j < nlat; j++) {
      const c = Math.sqrt(1 - mu[j] * mu[j]);
      for (let k = 0; k < nlon; k++) {
        const lambda = (2 * Math.PI * k) / nlon;
        this.f[j * nlon + k] = 2 * omega * (mu[j] * Math.cos(tilt) - Math.cos(lambda) * c * Math.sin(tilt));
        this.cos2[j * nlon + k] = 1 - mu[j] * mu[j];
      }
    }
    this.S = S;
  }

  // The state from winds u, v (m/s) and geopotential height z (m) on the
  // model's Gaussian grid (rows south to north). With `divergent` false, the
  // divergence of the wind is set to zero, as at the level of non-divergence.
  setState(u, v, z, { divergent = false } = {}) {
    const sh = this.sh, n = u.length;
    const A = new Float64Array(n), B = new Float64Array(n), P = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const c = Math.sqrt(this.cos2[i]);
      A[i] = u[i] * c; B[i] = v[i] * c; P[i] = GRAVITY * z[i];
    }
    const zeta = sh.curl(A, B, this.radius);
    zeta.re[0] = 0; zeta.im[0] = 0;
    const delta = divergent ? sh.divergence(A, B, this.radius) : sh.zeros();
    delta.re[0] = 0; delta.im[0] = 0;
    this.state = { zeta, delta, phi: sh.analyse(P) };
    this.time = 0;
  }

  // U = u cos(phi), V = v cos(phi) on the Gaussian grid, from the stream
  // function and velocity potential:
  //   U = [-(1 - mu^2) dpsi/dmu + dchi/dlambda] / a,   V = [dpsi/dlambda + (1 - mu^2) dchi/dmu] / a.
  wind(psi, chi) {
    const sh = this.sh, a = this.radius;
    const U = sh.fromFourier(sh.toFourier(psi, sh.H)), Uc = sh.fromFourier(sh.toFourier(chi, sh.P, true));
    const V = sh.fromFourier(sh.toFourier(psi, sh.P, true)), Vc = sh.fromFourier(sh.toFourier(chi, sh.H));
    for (let i = 0; i < U.length; i++) { U[i] = (-U[i] + Uc[i]) / a; V[i] = (V[i] + Vc[i]) / a; }
    return { U, V };
  }

  // The tendencies (1)-(3) of a state.
  tendency({ zeta, delta, phi }) {
    const sh = this.sh, a = this.radius;
    const psi = sh.invertLaplacian(zeta, a), chi = sh.invertLaplacian(delta, a);
    const { U, V } = this.wind(psi, chi);
    const q = sh.synthesise(zeta), P = sh.synthesise(phi);
    const n = U.length;
    const etaU = new Float64Array(n), etaV = new Float64Array(n), phiU = new Float64Array(n), phiV = new Float64Array(n), E = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const eta = q[i] + this.f[i];
      etaU[i] = eta * U[i]; etaV[i] = eta * V[i];
      phiU[i] = P[i] * U[i]; phiV[i] = P[i] * V[i];
      E[i] = P[i] + (U[i] * U[i] + V[i] * V[i]) / (2 * this.cos2[i]);
    }
    const dzeta = sh.divergence(etaU, etaV, a);
    const ddelta = sh.curl(etaU, etaV, a);
    const e = sh.analyse(E);
    const dphi = sh.divergence(phiU, phiV, a);
    for (let k = 0; k < this.S; k++) {
      const lap = sh.laplacianEigen(sh.n[k], a);
      dzeta.re[k] = -dzeta.re[k]; dzeta.im[k] = -dzeta.im[k];
      ddelta.re[k] -= lap * e.re[k]; ddelta.im[k] -= lap * e.im[k];
      dphi.re[k] = -dphi.re[k]; dphi.im[k] = -dphi.im[k];
    }
    dzeta.re[0] = dzeta.im[0] = ddelta.re[0] = ddelta.im[0] = dphi.re[0] = dphi.im[0] = 0;
    return { zeta: dzeta, delta: ddelta, phi: dphi };
  }

  // One step of length dt (s; negative to run backward): Runge-Kutta 4, then
  // the damping unless `damp` is false.
  step(dt, damp = true) {
    const S = this.S, x = this.state;
    const keys = ['zeta', 'delta', 'phi'];
    const add = (a, b, s) => {
      const out = {};
      for (const key of keys) {
        const re = new Float64Array(S), im = new Float64Array(S);
        for (let k = 0; k < S; k++) { re[k] = a[key].re[k] + s * b[key].re[k]; im[k] = a[key].im[k] + s * b[key].im[k]; }
        out[key] = { re, im };
      }
      return out;
    };
    const k1 = this.tendency(x);
    const k2 = this.tendency(add(x, k1, dt / 2));
    const k3 = this.tendency(add(x, k2, dt / 2));
    const k4 = this.tendency(add(x, k3, dt));
    for (const key of keys) {
      for (let k = 0; k < S; k++) {
        x[key].re[k] += (dt / 6) * (k1[key].re[k] + 2 * k2[key].re[k] + 2 * k3[key].re[k] + k4[key].re[k]);
        x[key].im[k] += (dt / 6) * (k1[key].im[k] + 2 * k2[key].im[k] + 2 * k3[key].im[k] + k4[key].im[k]);
      }
    }
    if (damp && this.damping) {
      const M = this.sh.M, top = M * (M + 1);
      const factor = (n) => (n === 0 ? 1 : Math.exp(-Math.abs(dt) * this.damping * ((n * (n + 1)) / top) ** 2));
      for (const key of keys) this.sh.scale(x[key], factor);
    }
    this.time += dt;
  }

  // Digital-filter initialisation over a span of `spanHours` (half backward,
  // half forward), with a cut-off period of `cutoffHours`.
  initialise({ dt = 900, spanHours = 6, cutoffHours = 6 } = {}) {
    const N = Math.round((spanHours * 3600) / (2 * dt));
    const weights = lanczosWeights(N, (2 * Math.PI * dt) / (cutoffHours * 3600));
    const start = this.copy(this.state), t0 = this.time;
    for (let s = 0; s < N; s++) this.step(-dt, false);
    const sum = { zeta: this.sh.zeros(), delta: this.sh.zeros(), phi: this.sh.zeros() };
    const accumulate = (w) => {
      for (const key of ['zeta', 'delta', 'phi']) {
        for (let k = 0; k < this.S; k++) { sum[key].re[k] += w * this.state[key].re[k]; sum[key].im[k] += w * this.state[key].im[k]; }
      }
    };
    accumulate(weights[0]);
    for (let s = 1; s <= 2 * N; s++) { this.step(dt, false); accumulate(weights[s]); }
    this.state = sum;
    this.time = t0;
    return { steps: 3 * N, unfiltered: start };
  }

  copy(state) {
    const out = {};
    for (const key of ['zeta', 'delta', 'phi']) out[key] = { re: Float64Array.from(state[key].re), im: Float64Array.from(state[key].im) };
    return out;
  }

  // Global means of mass (Phi) and of total energy (Phi |V|^2 / 2 + Phi^2 / 2),
  // which equations (1)-(3) conserve; computed on the Gaussian grid.
  invariants(state = this.state) {
    const sh = this.sh;
    const { U, V } = this.wind(sh.invertLaplacian(state.zeta, this.radius), sh.invertLaplacian(state.delta, this.radius));
    const P = sh.synthesise(state.phi);
    let mass = 0, energy = 0, w = 0;
    for (let j = 0; j < sh.nlat; j++) {
      for (let k = 0; k < sh.nlon; k++) {
        const i = j * sh.nlon + k, wj = sh.w[j];
        mass += wj * P[i];
        energy += wj * (P[i] * (U[i] * U[i] + V[i] * V[i]) / (2 * this.cos2[i]) + (P[i] * P[i]) / 2);
        w += wj;
      }
    }
    return { mass: mass / w, energy: energy / w };
  }

  // The global root-mean-square of d Phi / dt (m^2 s^-3): a measure of gravity-wave noise.
  noise(state = this.state) {
    const d = this.sh.synthesise(this.tendency(state).phi);
    let s = 0, w = 0;
    for (let j = 0; j < this.sh.nlat; j++) for (let k = 0; k < this.sh.nlon; k++) { s += this.sh.w[j] * d[j * this.sh.nlon + k] ** 2; w += this.sh.w[j]; }
    return Math.sqrt(s / w);
  }
}

// The weights h_k, k = -N..N (returned for k = -N first), of a low-pass filter
// with cut-off frequency theta (radians per step): sin(k theta) / (k pi),
// windowed by Lanczos's sigma factors, normalised to sum to one.
export function lanczosWeights(N, theta) {
  const h = [];
  for (let k = -N; k <= N; k++) {
    const ideal = k === 0 ? theta / Math.PI : Math.sin(k * theta) / (k * Math.PI);
    const x = (k * Math.PI) / (N + 1);
    const sigma = k === 0 ? 1 : Math.sin(x) / x;
    h.push(ideal * sigma);
  }
  const total = h.reduce((a, b) => a + b, 0);
  return h.map((x) => x / total);
}
