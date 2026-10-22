// Settling particles carried by the wind: the standard atmosphere, the
// terminal velocity of a small sphere (Stokes' law), the wind interpolated in
// height between levels, and the path of a falling parcel.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).
//
// The standard atmosphere is ICAO's (1993), the same as the U.S. Standard
// Atmosphere (1976) below 32 km: from T0 = 288.15 K and p0 = 101 325 Pa at sea
// level, temperature falls by 6.5 K per km to 216.65 K at 11 km, and is
// constant to 20 km; pressure follows hydrostatically, p = p0 (T/T0)^(g0/(R L))
// below 11 km and p = p11 exp(-g0 (z - 11 km)/(R T11)) above; the viscosity of
// air follows Sutherland's law, μ = β T^(3/2) / (T + S). Heights are
// geopotential; the geometric height differs by z²/r, 23 m at 12 km.
//
// A pressure level's height in the standard atmosphere (its pressure
// altitude): 850 hPa at 1457 m, 500 hPa at 5574 m, 250 hPa at 10 363 m. A
// forecast's level stands at its own geopotential height, within a few hundred
// metres of these in middle latitudes.
//
// Stokes' law. A sphere of diameter d and density ρp, falling through air of
// density ρa and viscosity μ at small Reynolds number Re = ρa w d / μ, reaches
// the terminal velocity w = (ρp - ρa) g d² / (18 μ). Valid for Re ≲ 0.1, where
// the drag is within 3 % of Stokes' (Schiller and Naumann's correction is
// 0.15 Re^0.687); at Re = 1 it overstates w by about 15 %. For volcanic glass
// (2500 kg m-3) Re = 0.1 at 27 µm near the ground and 34 µm at 10 km.
// The slip of air at the particle (Cunningham's correction, 1 + 2.5 λ/d with
// λ the mean free path) adds less than 3 % for d > 25 µm below 12 km and is
// neglected; so are the shapes of real ash, which fall more slowly than spheres
// of the same volume.
//
// References:
//   ICAO (1993), Manual of the ICAO Standard Atmosphere, Doc 7488/3, 3rd ed.
//   NOAA, NASA and USAF (1976), U.S. Standard Atmosphere, 1976, eq. 51 (viscosity).
//   Stokes, G. G. (1851), On the effect of the internal friction of fluids on
//   the motion of pendulums, Trans. Camb. Phil. Soc. 9, 8-106.
//   Schiller, L. and Naumann, A. (1933), Über die grundlegenden Berechnungen
//   bei der Schwerkraftaufbereitung, Z. Ver. Dtsch. Ing. 77, 318-320.
//   Seinfeld, J. H. and Pandis, S. N. (2016), Atmospheric Chemistry and
//   Physics, 3rd ed., §9.3-9.5 (settling, slip correction).

import { EARTH_RADIUS_M, DEG, wrapLongitude } from '../../geodesy/js/index.js';

export const G0 = 9.80665;               // standard gravity, m s-2
export const R_AIR = 287.05287;          // specific gas constant of dry air (ICAO), J kg-1 K-1
const T0 = 288.15, P0 = 101325, LAPSE = 0.0065, Z11 = 11000, T11 = 216.65;
const P11 = P0 * (T11 / T0) ** (G0 / (R_AIR * LAPSE));   // 22 632 Pa
const SUTHERLAND_BETA = 1.458e-6, SUTHERLAND_S = 110.4;

// The standard atmosphere at geopotential height z (m), from 0 to 20 km:
// { T (K), p (Pa), rho (kg m-3), mu (Pa s) }.
export function standardAtmosphere(z) {
  const h = Math.max(0, Math.min(z, 20000));
  const T = h <= Z11 ? T0 - LAPSE * h : T11;
  const p = h <= Z11 ? P0 * (T / T0) ** (G0 / (R_AIR * LAPSE)) : P11 * Math.exp((-G0 * (h - Z11)) / (R_AIR * T11));
  return { T, p, rho: p / (R_AIR * T), mu: (SUTHERLAND_BETA * T ** 1.5) / (T + SUTHERLAND_S) };
}

// The height (m) at which the standard atmosphere has pressure p (Pa), to 20 km.
export function pressureAltitude(p) {
  if (p >= P11) return (T0 / LAPSE) * (1 - (p / P0) ** ((R_AIR * LAPSE) / G0));
  return Z11 + ((R_AIR * T11) / G0) * Math.log(P11 / p);
}

// Stokes' terminal velocity of a sphere of diameter d (m) and density rhoP (kg m-3)
// at height z (m) in the standard atmosphere: { w (m/s, downward), reynolds }.
export function stokesSettling(d, rhoP, z) {
  const { rho, mu } = standardAtmosphere(z);
  const w = ((rhoP - rho) * G0 * d * d) / (18 * mu);
  return { w, reynolds: (rho * w * d) / mu };
}

// The wind at height z, linear in height between levels [{ z, u, v }] in
// increasing order of height; below the lowest and above the highest, theirs.
export function windAtHeight(z, levels) {
  const n = levels.length;
  if (z <= levels[0].z) return [levels[0].u, levels[0].v];
  if (z >= levels[n - 1].z) return [levels[n - 1].u, levels[n - 1].v];
  let k = 0;
  while (z > levels[k + 1].z) k++;
  const a = levels[k], b = levels[k + 1], w = (z - a.z) / (b.z - a.z);
  return [a.u + (b.u - a.u) * w, a.v + (b.v - a.v) * w];
}

// The two levels about a height, and the weight of the upper: [lower, upper, w],
// with the indices of `heights` (increasing); at or beyond the ends, the end twice.
export function levelBracket(z, heights) {
  const n = heights.length;
  if (z <= heights[0]) return [0, 0, 0];
  if (z >= heights[n - 1]) return [n - 1, n - 1, 0];
  let k = 0;
  while (z > heights[k + 1]) k++;
  return [k, k + 1, (z - heights[k]) / (heights[k + 1] - heights[k])];
}

// The path of a parcel carried by the wind and falling at its settling speed,
// from { lon, lat, z (m), t (ms UTC) }, by a fourth-order Runge–Kutta step on the
// sphere of mean radius, until it reaches the ground or `steps` steps of dt s:
//   flow(lon, lat, z, ms) -> [u, v] (m/s) or null outside the data;
//   settle(z) -> w (m/s, downward); ground(lon, lat) -> height of the ground (m).
// Returns { times, coords: [lon, lat, ...], heights, landed: { lon, lat, t } | null };
// the landing is placed between the last two steps, linearly in height.
export function settlingPath(start, flow, settle, { dt = 1800, steps = 96, ground = () => 0 } = {}) {
  const rate = (lon, lat, z, ms) => {
    const uv = flow(wrapLongitude(lon), Math.max(-89.99, Math.min(89.99, lat)), z, ms);
    if (!uv) return null;
    return [(uv[0] / (EARTH_RADIUS_M * Math.cos(lat * DEG))) / DEG, (uv[1] / EARTH_RADIUS_M) / DEG, -settle(z)];
  };
  let { lon, lat, z, t } = start;
  const times = [t], coords = [lon, lat], heights = [z];
  for (let s = 0; s < steps; s++) {
    const k1 = rate(lon, lat, z, t); if (!k1) break;
    const k2 = rate(lon + 0.5 * dt * k1[0], lat + 0.5 * dt * k1[1], z + 0.5 * dt * k1[2], t + 500 * dt); if (!k2) break;
    const k3 = rate(lon + 0.5 * dt * k2[0], lat + 0.5 * dt * k2[1], z + 0.5 * dt * k2[2], t + 500 * dt); if (!k3) break;
    const k4 = rate(lon + dt * k3[0], lat + dt * k3[1], z + dt * k3[2], t + 1000 * dt); if (!k4) break;
    let nlon = lon + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    let nlat = lat + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    const nz = z + (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    if (nlat > 90) { nlat = 180 - nlat; nlon += 180; } else if (nlat < -90) { nlat = -180 - nlat; nlon += 180; }
    nlon = wrapLongitude(nlon);
    const floor = ground(nlon, nlat), nt = t + 1000 * dt;
    if (nz <= floor) {
      // Landed within the step: where the height met the ground, linearly.
      const g0 = ground(lon, lat), f = Math.max(0, Math.min(1, (z - g0) / Math.max(1e-9, (z - g0) - (nz - floor))));
      let dl = nlon - lon;
      if (dl > 180) dl -= 360; else if (dl < -180) dl += 360;
      const landed = { lon: wrapLongitude(lon + dl * f), lat: lat + (nlat - lat) * f, t: t + f * 1000 * dt };
      times.push(landed.t); coords.push(landed.lon, landed.lat); heights.push(g0 + (floor - g0) * f);
      return { times, coords, heights, landed };
    }
    lon = nlon; lat = nlat; z = nz; t = nt;
    times.push(t); coords.push(lon, lat); heights.push(z);
  }
  return { times, coords, heights, landed: null };
}
