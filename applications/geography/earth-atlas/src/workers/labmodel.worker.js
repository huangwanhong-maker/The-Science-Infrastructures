// Runs one of the laboratory's models away from the page's main thread, and
// returns its fields every few hours on a regular grid:
//   'barotropic'     the barotropic vorticity equation (src/core/barotropic.js);
//   'shallow-water'  the shallow-water equations (src/core/shallowwater.js),
//                    started by a digital filter.
//
// Message in:  { model, u, v, z: Float64Array on the model's Gaussian grid (rows
//                south to north; z, the height, for the shallow-water model only),
//                hours, dt (s), every (h), grid: the output grid, truncation }
// Messages out: { progress: 0..1 } ..., then { done: true, times: [s], fields }
//   fields: { dz (barotropic: the balanced height change since the start) or
//   z (shallow water: the height), u, v (m/s), zeta (1/s) }, each a list of
//   Float32Arrays, one per time.

import { BarotropicModel, GRAVITY } from '../core/barotropic.js';
import { ShallowWaterModel } from '../core/shallowwater.js';

self.onmessage = (event) => {
  const { model: kind = 'barotropic', u, v, z, hours, dt, every, grid, truncation } = event.data;
  const shallow = kind === 'shallow-water';
  const model = shallow ? new ShallowWaterModel({ truncation }) : new BarotropicModel({ truncation });
  const sh = model.sh;
  const steps = Math.round((hours * 3600) / dt);
  let done = 0;
  const total = steps + (shallow ? Math.round((6 * 3600) / dt) * 1.5 : 0);
  const progress = () => { done++; if (done % 4 === 0) self.postMessage({ progress: done / total }); };
  let z0 = null;
  if (shallow) {
    model.setState(u, v, z);
    // The digital filter: half a span backward, a span forward (Lynch and Huang 1992).
    const step = model.step.bind(model);
    model.step = (h, damp) => { step(h, damp); progress(); };
    model.initialise({ dt, spanHours: 6, cutoffHours: 6 });
    model.step = step;
  } else {
    model.setWind(u, v);
    z0 = sh.synthesiseOn(model.balancedHeight(model.streamFunction()), grid);
  }
  const out = shallow ? { z: [], u: [], v: [], zeta: [] } : { dz: [], u: [], v: [], zeta: [] };
  const times = [];
  const record = () => {
    const zeta = shallow ? model.state.zeta : model.zeta;
    const psi = sh.invertLaplacian(zeta, model.radius);
    const chi = shallow ? sh.invertLaplacian(model.state.delta, model.radius) : null;
    if (shallow) {
      const h = sh.synthesiseOn(model.state.phi, grid);
      for (let i = 0; i < h.length; i++) h[i] /= GRAVITY;
      out.z.push(h);
    } else {
      const h = sh.synthesiseOn(model.balancedHeight(psi), grid);
      for (let i = 0; i < h.length; i++) h[i] -= z0[i];
      out.dz.push(h);
    }
    const { uu, vv } = windOn(sh, psi, chi, grid, model.radius);
    out.u.push(uu); out.v.push(vv);
    out.zeta.push(sh.synthesiseOn(zeta, grid));
    times.push(model.time);
  };
  record();
  const recordEvery = Math.max(1, Math.round((every * 3600) / dt));
  for (let s = 1; s <= steps; s++) {
    model.step(dt);
    progress();
    if (s % recordEvery === 0) record();
  }
  const transfer = Object.values(out).flat().map((a) => a.buffer);
  self.postMessage({ done: true, times, fields: out }, transfer);
};

// The wind on a regular grid from the stream function psi and, if given, the
// velocity potential chi: u = [-(1 - mu^2) dpsi/dmu + dchi/dlambda] / (a cos phi),
// v = [dpsi/dlambda + (1 - mu^2) dchi/dmu] / (a cos phi).
function windOn(sh, psi, chi, grid, radius) {
  const mus = Array.from({ length: grid.nlat }, (_, r) => Math.sin(((grid.lat0 - r * grid.dlat) * Math.PI) / 180));
  const [P, H] = sh.legendre(mus);
  const { M, size } = sh;
  const uu = new Float32Array(grid.nlat * grid.nlon), vv = new Float32Array(grid.nlat * grid.nlon);
  for (let r = 0; r < grid.nlat; r++) {
    const c = Math.sqrt(Math.max(1e-12, 1 - mus[r] * mus[r]));
    const Ur = new Float64Array(M + 1), Ui = new Float64Array(M + 1), Vr = new Float64Array(M + 1), Vi = new Float64Array(M + 1);
    for (let k = 0; k < size; k++) {
      const m = sh.m[k], p = P[r * size + k], h = H[r * size + k];
      Ur[m] -= (psi.re[k] * h) / radius; Ui[m] -= (psi.im[k] * h) / radius;          // -(1/a) sum psi H
      Vr[m] += (-m * psi.im[k] * p) / radius; Vi[m] += (m * psi.re[k] * p) / radius; // (1/a) sum i m psi P
      if (chi) {
        Ur[m] += (-m * chi.im[k] * p) / radius; Ui[m] += (m * chi.re[k] * p) / radius; // (1/a) sum i m chi P
        Vr[m] += (chi.re[k] * h) / radius; Vi[m] += (chi.im[k] * h) / radius;          // (1/a) sum chi H
      }
    }
    for (let col = 0; col < grid.nlon; col++) {
      const lam = ((grid.lon0 + col * grid.dlon) * Math.PI) / 180;
      let U = Ur[0], V = Vr[0];
      for (let m = 1; m <= M; m++) {
        const cm = Math.cos(m * lam), sm = Math.sin(m * lam);
        U += 2 * (Ur[m] * cm - Ui[m] * sm);
        V += 2 * (Vr[m] * cm - Vi[m] * sm);
      }
      uu[r * grid.nlon + col] = U / c;
      vv[r * grid.nlon + col] = V / c;
    }
  }
  return { uu, vv };
}
