// The barotropic vorticity equation on the sphere: the model of the first
// numerical weather prediction (Charney, Fjortoft and von Neumann 1950), run
// by the spectral transform method.
//
// Provisional stand-in for core/domains/meteorology (laboratory ADR 0015).
//
// Hypotheses: the flow at the level of non-divergence (about 500 hPa) is
// horizontal, non-divergent and barotropic. Its relative vorticity zeta is then
// carried by the flow, together with the planetary vorticity f = 2 Omega sin(phi):
//
//   d zeta / dt = - V . grad(zeta + f) = - div( (zeta + f) V ),   V = k x grad(psi),   zeta = lap(psi).   (1)
//
// The model holds zeta as spherical-harmonic coefficients (src/core/spectral.js),
// evaluates the flux (zeta + f) V on the Gaussian grid, and steps (1) in time by
// the classical fourth-order Runge-Kutta method. A scale-selective damping of
// order four (lap^2) removes the enstrophy that cascades to the truncation.
//
// Heights are diagnosed from the stream function by the linear balance
// equation, g lap(z) = div(f grad psi) (Holton and Hakim 2013, eq. 6.62), and
// added, as changes, to the analysed height of the first time, so that a
// forecast starts exactly from the analysis.
//
// Sources:
//   Charney, J. G., Fjortoft, R., and von Neumann, J. (1950), "Numerical
//   integration of the barotropic vorticity equation", Tellus 2, 237-254.
//   Haurwitz, B. (1940), "The motion of atmospheric disturbances on the
//   spherical earth", J. Mar. Res. 3, 254-267: the Rossby-Haurwitz wave (tests).
//   Williamson, D. L., et al. (1992), J. Comput. Phys. 102, 211-224, test case 6.
//   Durran (2010), Numerical Methods for Fluid Dynamics, sec. 6.4.

import { SphericalHarmonics, EARTH_RADIUS, OMEGA } from './spectral.js';

export const GRAVITY = 9.80665;

export class BarotropicModel {
  constructor({ truncation = 42, radius = EARTH_RADIUS, omega = OMEGA, dampingHours = 8 } = {}) {
    this.sh = new SphericalHarmonics(truncation);
    this.radius = radius;
    this.omega = omega;
    this.damping = dampingHours > 0 ? 1 / (dampingHours * 3600) : 0;
    this.zeta = this.sh.zeros();
    this.time = 0;
    const { nlat, nlon, mu } = this.sh;
    this.f = new Float64Array(nlat * nlon);
    for (let j = 0; j < nlat; j++) for (let k = 0; k < nlon; k++) this.f[j * nlon + k] = 2 * omega * mu[j];
  }

  // Vorticity from winds u, v (m/s) on the model's Gaussian grid (rows south to north).
  setWind(u, v) {
    const { nlat, nlon, mu } = this.sh;
    const A = new Float64Array(nlat * nlon), B = new Float64Array(nlat * nlon);
    for (let j = 0; j < nlat; j++) {
      const c = Math.sqrt(1 - mu[j] * mu[j]);
      for (let k = 0; k < nlon; k++) { A[j * nlon + k] = u[j * nlon + k] * c; B[j * nlon + k] = v[j * nlon + k] * c; }
    }
    this.zeta = this.sh.curl(A, B, this.radius);
    this.zeta.re[0] = 0; this.zeta.im[0] = 0;          // the global mean of vorticity is zero
    this.time = 0;
  }

  setStreamFunction(psi) {
    this.zeta = this.sh.scale({ re: Float64Array.from(psi.re), im: Float64Array.from(psi.im) },
      (n) => this.sh.laplacianEigen(n, this.radius));
    this.time = 0;
  }

  streamFunction(zeta = this.zeta) { return this.sh.invertLaplacian(zeta, this.radius); }

  // d zeta / dt, equation (1).
  tendency(zeta) {
    const sh = this.sh;
    const psi = sh.invertLaplacian(zeta, this.radius);
    const { U, V } = sh.windFromStreamFunction(psi, this.radius);
    const q = sh.synthesise(zeta);
    for (let i = 0; i < q.length; i++) {
      const a = q[i] + this.f[i];
      U[i] *= a; V[i] *= a;
    }
    const d = sh.divergence(U, V, this.radius);
    for (let k = 0; k < sh.size; k++) { d.re[k] = -d.re[k]; d.im[k] = -d.im[k]; }
    d.re[0] = 0; d.im[0] = 0;
    return d;
  }

  // One step of length dt (s): fourth-order Runge-Kutta, then the damping.
  step(dt) {
    const z = this.zeta, S = this.sh.size;
    const add = (a, b, s) => {
      const out = { re: new Float64Array(S), im: new Float64Array(S) };
      for (let k = 0; k < S; k++) { out.re[k] = a.re[k] + s * b.re[k]; out.im[k] = a.im[k] + s * b.im[k]; }
      return out;
    };
    const k1 = this.tendency(z);
    const k2 = this.tendency(add(z, k1, dt / 2));
    const k3 = this.tendency(add(z, k2, dt / 2));
    const k4 = this.tendency(add(z, k3, dt));
    for (let k = 0; k < S; k++) {
      z.re[k] += (dt / 6) * (k1.re[k] + 2 * k2.re[k] + 2 * k3.re[k] + k4.re[k]);
      z.im[k] += (dt / 6) * (k1.im[k] + 2 * k2.im[k] + 2 * k3.im[k] + k4.im[k]);
    }
    if (this.damping) {
      const M = this.sh.M, top = M * (M + 1);
      this.sh.scale(z, (n) => Math.exp(-dt * this.damping * ((n * (n + 1)) / top) ** 2));
    }
    this.time += dt;
  }

  // Kinetic energy and enstrophy per unit area (global means), from the spectrum:
  // the mean of f^2 over the sphere is (1/2) sum_n (|a_n^0|^2 + 2 sum_{m>0} |a_n^m|^2).
  invariants(zeta = this.zeta) {
    const sh = this.sh;
    let e = 0, z = 0;
    for (let k = 0; k < sh.size; k++) {
      const n = sh.n[k];
      if (n === 0) continue;
      const w = sh.m[k] === 0 ? 0.5 : 1;
      const a2 = zeta.re[k] ** 2 + zeta.im[k] ** 2;
      z += w * a2 / 2;
      e += (w * a2) / (2 * ((n * (n + 1)) / this.radius ** 2));
    }
    return { energy: e, enstrophy: z };
  }

  // The height change implied by the stream function under linear balance.
  balancedHeight(psi) {
    const sh = this.sh;
    const { U, V } = sh.windFromStreamFunction(psi, this.radius);
    const A = new Float64Array(U.length), B = new Float64Array(U.length);
    for (let i = 0; i < U.length; i++) { A[i] = this.f[i] * V[i]; B[i] = -this.f[i] * U[i]; }
    const rhs = sh.divergence(A, B, this.radius);      // div(f grad psi)
    return sh.scale(rhs, (n) => (n === 0 ? 0 : 1 / (GRAVITY * sh.laplacianEigen(n, this.radius))));
  }
}

// Root-mean-square difference of two fields on a regular grid, weighted by the
// cosine of latitude, between two latitudes.
export function rmsDifference(a, b, grid, latMin = -90, latMax = 90) {
  let s = 0, w = 0;
  for (let r = 0; r < grid.nlat; r++) {
    const lat = grid.lat0 - r * grid.dlat;
    if (lat < latMin || lat > latMax) continue;
    const c = Math.cos((lat * Math.PI) / 180);
    for (let k = 0; k < grid.nlon; k++) {
      const d = a[r * grid.nlon + k] - b[r * grid.nlon + k];
      if (d === d) { s += c * d * d; w += c; }
    }
  }
  return w > 0 ? Math.sqrt(s / w) : NaN;
}
