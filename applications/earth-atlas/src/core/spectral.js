// The spectral transform method on the sphere: spherical harmonics, Gaussian
// grids, and the transforms between them, for triangular truncation T_M.
//
// Provisional stand-in for calx (numerical methods) and physis (discretisation)
// (laboratory ADR 0015).
//
// A real field on the sphere is written
//
//   f(lambda, mu) = sum_{m=-M..M} sum_{n=|m|..M} a_n^m Pbar_n^m(mu) e^{i m lambda},   mu = sin(phi),
//
// with a_n^{-m} = conj(a_n^m), so only m >= 0 is stored. Pbar_n^m are the
// associated Legendre functions normalised so that the integral of their square
// over mu in [-1, 1] is 1. Coefficients are stored as Float64Arrays of real and
// imaginary parts, indexed by index(m, n).
//
// Analysis integrates over longitude by a discrete Fourier transform (by the
// radix-2 fast Fourier transform when the number of longitudes is a power of
// two, as on the grids used here) and over latitude by Gaussian quadrature, which is exact for the products of a
// quadratic grid (nlat >= (3M + 1) / 2, nlon >= 3M + 1), so the transforms of
// products of two fields carry no aliasing.
//
// Sources:
//   Durran (2010), Numerical Methods for Fluid Dynamics, 2nd ed., sec. 6.4 (the
//   spectral transform method on the sphere).
//   Hack and Jakob (1992), Description of a global shallow water model based on
//   the spectral transform method, NCAR/TN-343+STR, sec. 3: the recurrences for
//   Pbar_n^m and H_n^m = (1 - mu^2) dPbar_n^m/dmu, and the transforms of
//   divergence and curl.
//   Press et al. (2007), Numerical Recipes, 3rd ed., sec. 4.6: Gauss-Legendre nodes.

export const EARTH_RADIUS = 6371220;      // m, as in Williamson et al. (1992)
export const OMEGA = 7.292e-5;            // s^-1, Earth's rotation rate (Williamson et al. 1992)

// Nodes (mu, ascending from the south) and weights of n-point Gauss-Legendre quadrature.
export function gaussLegendre(n) {
  const mu = new Float64Array(n), w = new Float64Array(n);
  for (let i = 0; i < Math.ceil(n / 2); i++) {
    let x = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let dp = 0;
    for (let iter = 0; iter < 100; iter++) {
      let p0 = 1, p1 = x;
      for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      dp = (n * (x * p1 - p0)) / (x * x - 1);
      const dx = p1 / dp;
      x -= dx;
      if (Math.abs(dx) < 1e-15) break;
    }
    mu[n - 1 - i] = x; mu[i] = -x;
    w[i] = w[n - 1 - i] = 2 / ((1 - x * x) * dp * dp);
  }
  return { mu, w };
}

// The discrete Fourier transform of n points (n a power of two), in place:
// X_m = sum_k x_k e^{sign i 2 pi m k / n}, by the iterative radix-2 algorithm
// (Cooley and Tukey 1965; Press et al. 2007, sec. 12.2).
export class FFT {
  constructor(n) {
    this.n = n;
    this.rev = new Uint32Array(n);
    const bits = Math.log2(n);
    for (let i = 0; i < n; i++) {
      let r = 0;
      for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b);
      this.rev[i] = r;
    }
    this.cos = new Float64Array(n / 2);
    this.sin = new Float64Array(n / 2);
    for (let k = 0; k < n / 2; k++) { this.cos[k] = Math.cos((2 * Math.PI * k) / n); this.sin[k] = Math.sin((2 * Math.PI * k) / n); }
  }

  transform(re, im, sign) {
    const n = this.n, rev = this.rev;
    for (let i = 0; i < n; i++) {
      const j = rev[i];
      if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let size = 2; size <= n; size *= 2) {
      const half = size / 2, step = n / size;
      for (let start = 0; start < n; start += size) {
        for (let k = 0; k < half; k++) {
          const c = this.cos[k * step], s = sign * this.sin[k * step];
          const a = start + k, b = a + half;
          const tr = re[b] * c - im[b] * s, ti = re[b] * s + im[b] * c;
          re[b] = re[a] - tr; im[b] = im[a] - ti;
          re[a] += tr; im[a] += ti;
        }
      }
    }
  }
}

export class SphericalHarmonics {
  // M: truncation; nlat, nlon: the Gaussian grid (defaults: the quadratic grid).
  constructor(M, nlat = 2 * Math.ceil((3 * M + 1) / 4), nlon = 2 * Math.ceil((3 * M + 1) / 2)) {
    this.M = M;
    this.nlat = nlat;
    this.nlon = nlon;
    this.size = ((M + 1) * (M + 2)) / 2;              // coefficients with 0 <= m <= n <= M
    const { mu, w } = gaussLegendre(nlat);
    this.mu = mu;
    this.w = w;
    this.lat = Array.from(mu, (x) => (Math.asin(x) * 180) / Math.PI);
    this.lon = Array.from({ length: nlon }, (_, k) => (360 * k) / nlon);
    this.n = new Int32Array(this.size);
    this.m = new Int32Array(this.size);
    for (let m = 0; m <= M; m++) for (let n = m; n <= M; n++) { const k = this.index(m, n); this.n[k] = n; this.m[k] = m; }
    [this.P, this.H] = this.legendre(mu);
    // Fourier tables.
    this.cos = new Float64Array((M + 1) * nlon);
    this.sin = new Float64Array((M + 1) * nlon);
    for (let m = 0; m <= M; m++) for (let k = 0; k < nlon; k++) {
      this.cos[m * nlon + k] = Math.cos((2 * Math.PI * m * k) / nlon);
      this.sin[m * nlon + k] = Math.sin((2 * Math.PI * m * k) / nlon);
    }
    this.fft = (nlon & (nlon - 1)) === 0 && nlon >= 2 * M + 2 ? new FFT(nlon) : null;
  }

  index(m, n) { return m * (this.M + 1) - (m * (m - 1)) / 2 + (n - m); }

  zeros() { return { re: new Float64Array(this.size), im: new Float64Array(this.size) }; }

  // Pbar_n^m and H_n^m at the given mu, for 0 <= m <= n <= M: arrays indexed [j * size + index(m, n)].
  legendre(mus) {
    const M = this.M, S = this.size, J = mus.length;
    const P = new Float64Array(J * S), H = new Float64Array(J * S);
    const eps = (n, m) => (n <= m ? 0 : Math.sqrt((n * n - m * m) / (4 * n * n - 1)));
    for (let j = 0; j < J; j++) {
      const x = mus[j], c = Math.sqrt(Math.max(0, 1 - x * x));
      let pmm = Math.sqrt(0.5);                        // Pbar_0^0
      for (let m = 0; m <= M; m++) {
        if (m > 0) pmm *= Math.sqrt((2 * m + 1) / (2 * m)) * c;
        // Pbar_n^m for n = m .. M + 1 by the three-term recurrence.
        const col = new Float64Array(M + 3);
        col[0] = pmm;
        col[1] = Math.sqrt(2 * m + 3) * x * pmm;
        for (let n = m + 2; n <= M + 1; n++) {
          col[n - m] = (x * col[n - 1 - m] - eps(n - 1, m) * col[n - 2 - m]) / eps(n, m);
        }
        for (let n = m; n <= M; n++) {
          const k = j * S + this.index(m, n);
          P[k] = col[n - m];
          // (1 - mu^2) dPbar_n^m/dmu = -n eps_{n+1}^m Pbar_{n+1}^m + (n + 1) eps_n^m Pbar_{n-1}^m
          H[k] = -n * eps(n + 1, m) * col[n + 1 - m] + (n + 1) * eps(n, m) * (n > m ? col[n - 1 - m] : 0);
        }
      }
    }
    return [P, H];
  }

  // Fourier coefficients F_m (m = 0..M) of each latitude row of a grid field
  // (rows south to north, as this.mu): F_m = (1/nlon) sum_k f_k e^{-i m lambda_k}.
  fourier(grid) {
    const { M, nlat, nlon } = this;
    const re = new Float64Array(nlat * (M + 1)), im = new Float64Array(nlat * (M + 1));
    if (this.fft) {
      const xr = new Float64Array(nlon), xi = new Float64Array(nlon);
      for (let j = 0; j < nlat; j++) {
        xr.set(grid.subarray ? grid.subarray(j * nlon, (j + 1) * nlon) : grid.slice(j * nlon, (j + 1) * nlon));
        xi.fill(0);
        this.fft.transform(xr, xi, -1);
        for (let m = 0; m <= M; m++) { re[j * (M + 1) + m] = xr[m] / nlon; im[j * (M + 1) + m] = xi[m] / nlon; }
      }
      return { re, im };
    }
    for (let j = 0; j < nlat; j++) {
      const row = j * nlon;
      for (let m = 0; m <= M; m++) {
        let a = 0, b = 0;
        const t = m * nlon;
        for (let k = 0; k < nlon; k++) { const f = grid[row + k]; a += f * this.cos[t + k]; b -= f * this.sin[t + k]; }
        re[j * (M + 1) + m] = a / nlon; im[j * (M + 1) + m] = b / nlon;
      }
    }
    return { re, im };
  }

  // A grid field from Fourier coefficients: f_k = F_0 + 2 Re sum_{m>0} F_m e^{i m lambda_k}.
  fromFourier(F, out = new Float64Array(this.nlat * this.nlon)) {
    const { M, nlat, nlon } = this;
    if (this.fft) {
      // f_k = sum over m = -M..M of F_m e^{i m lambda_k}, with F_{-m} = conj(F_m).
      const xr = new Float64Array(nlon), xi = new Float64Array(nlon);
      for (let j = 0; j < nlat; j++) {
        xr.fill(0); xi.fill(0);
        xr[0] = F.re[j * (M + 1)];
        for (let m = 1; m <= M; m++) {
          const a = F.re[j * (M + 1) + m], b = F.im[j * (M + 1) + m];
          xr[m] = a; xi[m] = b; xr[nlon - m] = a; xi[nlon - m] = -b;
        }
        this.fft.transform(xr, xi, 1);
        out.set(xr, j * nlon);
      }
      return out;
    }
    for (let j = 0; j < nlat; j++) {
      for (let k = 0; k < nlon; k++) {
        let s = F.re[j * (M + 1)];
        for (let m = 1; m <= M; m++) {
          const a = F.re[j * (M + 1) + m], b = F.im[j * (M + 1) + m];
          s += 2 * (a * this.cos[m * nlon + k] - b * this.sin[m * nlon + k]);
        }
        out[j * nlon + k] = s;
      }
    }
    return out;
  }

  // Spectral coefficients of a grid field.
  analyse(grid) {
    const F = this.fourier(grid), a = this.zeros();
    const { M, nlat, size } = this;
    for (let j = 0; j < nlat; j++) {
      const wj = this.w[j];
      for (let k = 0; k < size; k++) {
        const m = this.m[k], p = this.P[j * size + k] * wj;
        a.re[k] += p * F.re[j * (M + 1) + m];
        a.im[k] += p * F.im[j * (M + 1) + m];
      }
    }
    return a;
  }

  // Fourier coefficients of sum_n a_n^m T_n^m(mu_j), with T = P (or H, given as table).
  toFourier(a, table = this.P, derivLambda = false) {
    const { M, nlat, size } = this;
    const F = { re: new Float64Array(nlat * (M + 1)), im: new Float64Array(nlat * (M + 1)) };
    for (let j = 0; j < nlat; j++) {
      for (let k = 0; k < size; k++) {
        const m = this.m[k], t = table[j * size + k];
        let re = a.re[k] * t, im = a.im[k] * t;
        if (derivLambda) { const r = -m * im; im = m * re; re = r; }       // multiply by i m
        F.re[j * (M + 1) + m] += re;
        F.im[j * (M + 1) + m] += im;
      }
    }
    return F;
  }

  synthesise(a) { return this.fromFourier(this.toFourier(a)); }

  // Coefficients of div = [dA/dlambda + (1 - mu^2) dB/dmu] / (a (1 - mu^2)) (sign +1),
  // or of curl = [dB/dlambda - (1 - mu^2) dA/dmu] / (a (1 - mu^2)) (sign -1, with A, B swapped),
  // for grid fields A = X cos(phi), B = Y cos(phi) of a vector (X, Y).
  divergence(A, B, radius = EARTH_RADIUS) {
    const FA = this.fourier(A), FB = this.fourier(B), d = this.zeros();
    const { M, nlat, size } = this;
    for (let j = 0; j < nlat; j++) {
      const c = this.w[j] / ((1 - this.mu[j] * this.mu[j]) * radius);
      for (let k = 0; k < size; k++) {
        const m = this.m[k], i = j * (M + 1) + m, p = this.P[j * size + k], h = this.H[j * size + k];
        // i m A_m P - B_m H, integrated (by parts in mu for the second term)
        d.re[k] += c * (-m * FA.im[i] * p - FB.re[i] * h);
        d.im[k] += c * (m * FA.re[i] * p - FB.im[i] * h);
      }
    }
    return d;
  }

  // Coefficients of the vertical component of the curl of (X, Y), given A = X cos(phi), B = Y cos(phi).
  curl(A, B, radius = EARTH_RADIUS) {
    const FA = this.fourier(A), FB = this.fourier(B), z = this.zeros();
    const { M, nlat, size } = this;
    for (let j = 0; j < nlat; j++) {
      const c = this.w[j] / ((1 - this.mu[j] * this.mu[j]) * radius);
      for (let k = 0; k < size; k++) {
        const m = this.m[k], i = j * (M + 1) + m, p = this.P[j * size + k], h = this.H[j * size + k];
        // i m B_m P + A_m H
        z.re[k] += c * (-m * FB.im[i] * p + FA.re[i] * h);
        z.im[k] += c * (m * FB.re[i] * p + FA.im[i] * h);
      }
    }
    return z;
  }

  // U = u cos(phi) and V = v cos(phi) on the grid for the non-divergent flow of
  // stream function psi: u = -(1/a) dpsi/dphi, v = (1 / (a cos phi)) dpsi/dlambda.
  windFromStreamFunction(psi, radius = EARTH_RADIUS) {
    const U = this.fromFourier(this.toFourier(psi, this.H));
    const V = this.fromFourier(this.toFourier(psi, this.P, true));
    for (let i = 0; i < U.length; i++) { U[i] *= -1 / radius; V[i] /= radius; }
    return { U, V };
  }

  // In place: multiply each coefficient by f(n).
  scale(a, f) {
    for (let k = 0; k < this.size; k++) { const s = f(this.n[k]); a.re[k] *= s; a.im[k] *= s; }
    return a;
  }

  laplacianEigen(n, radius = EARTH_RADIUS) { return (-n * (n + 1)) / (radius * radius); }

  // The stream function of a vorticity (the inverse Laplacian; the n = 0 term is set to 0).
  invertLaplacian(zeta, radius = EARTH_RADIUS) {
    const psi = { re: Float64Array.from(zeta.re), im: Float64Array.from(zeta.im) };
    return this.scale(psi, (n) => (n === 0 ? 0 : 1 / this.laplacianEigen(n, radius)));
  }

  // Values of a spectral field on a regular latitude-longitude grid (rows north
  // to south from lat0, columns east from lon0), for display and verification.
  synthesiseOn(a, grid) {
    const mus = Array.from({ length: grid.nlat }, (_, r) => Math.sin(((grid.lat0 - r * grid.dlat) * Math.PI) / 180));
    const [P] = this.legendre(mus);
    const { M, size } = this;
    const out = new Float32Array(grid.nlat * grid.nlon);
    const cosT = [], sinT = [];
    for (let m = 0; m <= M; m++) {
      cosT.push(Float64Array.from({ length: grid.nlon }, (_, c) => Math.cos((m * (grid.lon0 + c * grid.dlon) * Math.PI) / 180)));
      sinT.push(Float64Array.from({ length: grid.nlon }, (_, c) => Math.sin((m * (grid.lon0 + c * grid.dlon) * Math.PI) / 180)));
    }
    for (let r = 0; r < grid.nlat; r++) {
      const Fr = new Float64Array(M + 1), Fi = new Float64Array(M + 1);
      for (let k = 0; k < size; k++) { const m = this.m[k], p = P[r * size + k]; Fr[m] += a.re[k] * p; Fi[m] += a.im[k] * p; }
      for (let c = 0; c < grid.nlon; c++) {
        let s = Fr[0];
        for (let m = 1; m <= M; m++) s += 2 * (Fr[m] * cosT[m][c] - Fi[m] * sinT[m][c]);
        out[r * grid.nlon + c] = s;
      }
    }
    return out;
  }
}
