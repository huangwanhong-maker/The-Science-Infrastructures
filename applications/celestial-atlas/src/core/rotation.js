// The orientation of the Sun, the Moon and the planets: the direction of each
// body's north pole and the angle of its prime meridian, at an instant.
//
// Provisional stand-in for core/domains/astronomy (laboratory ADR 0016).
//
// The IAU model gives, in the ICRF, the right ascension alpha0 and declination
// delta0 of the north pole and the angle W of the prime meridian, measured
// eastward along the body's equator from the node of the equator on the ICRF
// equator. A vector in the ICRF becomes body-fixed by
//
//   B = Rz(W) Rx(90 deg - delta0) Rz(90 deg + alpha0),   v_body = B v_ICRF,   (1)
//
// with x toward the prime meridian, z toward the north pole, and east
// longitude atan2(y, x).
//
// Times: d is days and T Julian centuries from J2000.0, in TDB.
//
// Source: Archinal, B. A., et al. (2018), "Report of the IAU Working Group on
// Cartographic Coordinates and Rotational Elements: 2015", Celest. Mech. Dyn.
// Astron. 130:22, Tables 1-2 (planets and the Sun) and Table 3 (the Moon).

import { DEG, rotX, rotZ, multiply } from './frames.js';

const sin = (x) => Math.sin(x * DEG), cos = (x) => Math.cos(x * DEG);

const MODELS = {
  sun: () => [286.13, 63.87, (d) => 84.176 + 14.1844000 * d],
  mercury: (d, T) => {
    const M = [174.7910857 + 4.092335 * d, 349.5821714 + 8.184670 * d, 164.3732571 + 12.277005 * d,
      339.1643429 + 16.369340 * d, 153.9554286 + 20.461675 * d];
    return [281.0103 - 0.0328 * T, 61.4155 - 0.0049 * T,
      () => 329.5988 + 6.1385108 * d + 0.01067257 * sin(M[0]) - 0.00112309 * sin(M[1]) - 0.00011040 * sin(M[2])
        - 0.00002539 * sin(M[3]) - 0.00000571 * sin(M[4])];
  },
  venus: () => [272.76, 67.16, (d) => 160.20 - 1.4813688 * d],
  earth: (d, T) => [0 - 0.641 * T, 90 - 0.557 * T, () => 190.147 + 360.9856235 * d],
  moon: (d, T) => {
    const E = [125.045 - 0.0529921 * d, 250.089 - 0.1059842 * d, 260.008 + 13.0120009 * d, 176.625 + 13.3407154 * d,
      357.529 + 0.9856003 * d, 311.589 + 26.4057084 * d, 134.963 + 13.0649930 * d, 276.617 + 0.3287146 * d,
      34.226 + 1.7484877 * d, 15.134 - 0.1589763 * d, 119.743 + 0.0036096 * d, 239.961 + 0.1643573 * d,
      25.053 + 12.9590088 * d];
    const s = (i) => sin(E[i - 1]), c = (i) => cos(E[i - 1]);
    return [
      269.9949 + 0.0031 * T - 3.8787 * s(1) - 0.1204 * s(2) + 0.0700 * s(3) - 0.0172 * s(4) + 0.0072 * s(6)
        - 0.0052 * s(10) + 0.0043 * s(13),
      66.5392 + 0.0130 * T + 1.5419 * c(1) + 0.0239 * c(2) - 0.0278 * c(3) + 0.0068 * c(4) - 0.0029 * c(6)
        + 0.0009 * c(7) + 0.0008 * c(10) - 0.0009 * c(13),
      () => 38.3213 + 13.17635815 * d - 1.4e-12 * d * d + 3.5610 * s(1) + 0.1208 * s(2) - 0.0642 * s(3) + 0.0158 * s(4)
        + 0.0252 * s(5) - 0.0066 * s(6) - 0.0047 * s(7) - 0.0046 * s(8) + 0.0028 * s(9) + 0.0052 * s(10)
        + 0.0040 * s(11) + 0.0019 * s(12) - 0.0044 * s(13),
    ];
  },
  mars: (d, T) => [
    317.269202 - 0.10927547 * T + 0.000068 * sin(198.991226 + 19139.4819985 * T) + 0.000238 * sin(226.292679 + 38280.8511281 * T)
      + 0.000052 * sin(249.663391 + 57420.7251593 * T) + 0.000009 * sin(266.183510 + 76560.6367950 * T)
      + 0.419057 * sin(79.398797 + 0.5042615 * T),
    54.432516 - 0.05827105 * T + 0.000051 * cos(122.433576 + 19139.9407476 * T) + 0.000141 * cos(43.058401 + 38280.8753272 * T)
      + 0.000031 * cos(57.663379 + 57420.7517205 * T) + 0.000005 * cos(79.476401 + 76560.6495004 * T)
      + 1.591274 * cos(166.325722 + 0.5042615 * T),
    () => 176.049863 + 350.891982443297 * d + 0.000145 * sin(129.071773 + 19140.0328244 * T)
      + 0.000157 * sin(36.352167 + 38281.0473591 * T) + 0.000040 * sin(56.668646 + 57420.9295360 * T)
      + 0.000001 * sin(67.364003 + 76560.2552215 * T) + 0.000001 * sin(104.792680 + 95700.4387578 * T)
      + 0.584542 * sin(95.391654 + 0.5042615 * T),
  ],
  jupiter: (d, T) => {
    const J = [99.360714 + 4850.4046 * T, 175.895369 + 1191.9605 * T, 300.323162 + 262.5475 * T,
      114.012305 + 6070.2476 * T, 49.511251 + 64.3000 * T];
    return [
      268.056595 - 0.006499 * T + 0.000117 * sin(J[0]) + 0.000938 * sin(J[1]) + 0.001432 * sin(J[2]) + 0.000030 * sin(J[3]) + 0.002150 * sin(J[4]),
      64.495303 + 0.002413 * T + 0.000050 * cos(J[0]) + 0.000404 * cos(J[1]) + 0.000617 * cos(J[2]) - 0.000013 * cos(J[3]) + 0.000926 * cos(J[4]),
      () => 284.95 + 870.5360000 * d,
    ];
  },
  saturn: (d, T) => [40.589 - 0.036 * T, 83.537 - 0.004 * T, () => 38.90 + 810.7939024 * d],
  uranus: () => [257.311, -15.175, (d) => 203.81 - 501.1600928 * d],
  neptune: (d, T) => {
    const N = 357.85 + 52.316 * T;
    return [299.36 + 0.70 * sin(N), 43.46 - 0.51 * cos(N), () => 249.978 + 541.1397757 * d - 0.48 * sin(N)];
  },
  pluto: () => [132.993, -6.163, (d) => 302.695 + 56.3625225 * d],
};

// Pole and prime meridian, degrees, at a TDB Julian date.
export function orientation(key, jdTdb) {
  const model = MODELS[key];
  if (!model) return null;
  const d = jdTdb - 2451545.0, T = d / 36525;
  const [alpha0, delta0, W] = model(d, T);
  return { alpha0, delta0, W: ((W(d) % 360) + 360) % 360 };
}

// The matrix of (1): ICRF to body-fixed.
export function bodyFixedMatrix(key, jdTdb) {
  const o = orientation(key, jdTdb);
  if (!o) return null;
  return multiply(rotZ(o.W * DEG), multiply(rotX((90 - o.delta0) * DEG), rotZ((90 + o.alpha0) * DEG)));
}

// Longitudes are counted positive to the west on bodies that turn in the
// prograde sense, by the IAU convention (east on the Sun, the Earth, the Moon).
export const WEST_POSITIVE = new Set(['mercury', 'mars', 'jupiter', 'saturn', 'neptune']);
