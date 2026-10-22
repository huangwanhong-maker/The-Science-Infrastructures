// Relations between galaxies: supergalactic positions, the tidal index and main
// disturber of each galaxy, and the groups that tidal relations define.
//
// Provisional stand-in for core/domains/astronomy (laboratory ADR 0011, 0013).
//
// Sources:
//   Karachentsev and Makarov (1999), IAU Symp. 186, 109; Karachentsev, Makarov
//   and Kaisina (2013), AJ 145, 101: the tidal index
//       Θ1 = max_k log10(M_k / D_ik³) + C,     C = −10.96,
//   where M_k is the mass of neighbour k in solar masses, estimated from its
//   K-band luminosity with M/L_K = 1, and D_ik the distance in space (Mpc).
//   The neighbour that attains the maximum is the main disturber. Θ1 = 0 when
//   the Keplerian period of the galaxy around its main disturber equals the age
//   of the universe: Θ1 > 0 marks a galaxy bound in a group, Θ1 < 0 one in the
//   field. Θ5 sums the five strongest terms instead of taking the largest.

import { ICRS_TO_SUPERGALACTIC, apply, fromSpherical, DEG } from './frames.js';

export const TIDAL_C = -10.96;

// Supergalactic position (in the units of dist) from ICRS right ascension and
// declination in degrees.
export function supergalactic(raDeg, decDeg, dist) {
  const u = apply(ICRS_TO_SUPERGALACTIC, fromSpherical(raDeg * DEG, decDeg * DEG));
  return [u[0] * dist, u[1] * dist, u[2] * dist];
}

// Tidal indices for galaxies at positions (3n coordinates, Mpc) with K-band
// luminosities logL (log10 of solar units; NaN where unknown). A galaxy of
// unknown luminosity still receives an index, but disturbs no other galaxy.
export function tidalIndices(positions, logL, { massToLight = 1, constant = TIDAL_C } = {}) {
  const n = logL.length;
  const theta1 = new Float64Array(n).fill(NaN);
  const theta5 = new Float64Array(n).fill(NaN);
  const disturber = new Int32Array(n).fill(-1);
  const logMass = Float64Array.from(logL, (l) => l + Math.log10(massToLight));
  const strongest = [];
  for (let i = 0; i < n; i++) {
    const xi = positions[3 * i], yi = positions[3 * i + 1], zi = positions[3 * i + 2];
    let best = -Infinity, bestK = -1;
    strongest.length = 0;
    for (let k = 0; k < n; k++) {
      if (k === i || !Number.isFinite(logMass[k])) continue;
      const dx = positions[3 * k] - xi, dy = positions[3 * k + 1] - yi, dz = positions[3 * k + 2] - zi;
      const d2 = dx * dx + dy * dy + dz * dz;
      if (d2 === 0) continue;
      const f = logMass[k] - 1.5 * Math.log10(d2);
      if (f > best) { best = f; bestK = k; }
      if (strongest.length < 5) {
        strongest.push(f);
        strongest.sort((a, b) => b - a);
      } else if (f > strongest[4]) {
        strongest[4] = f;
        strongest.sort((a, b) => b - a);
      }
    }
    if (bestK < 0) continue;
    theta1[i] = best + constant;
    disturber[i] = bestK;
    let sum = 0;
    for (const f of strongest) sum += 10 ** (f - best);
    theta5[i] = best + Math.log10(sum) + constant;
  }
  return { theta1, theta5, disturber };
}

// Luminosities of main disturbers that the catalogue gives without one (in the
// UNGC, the Milky Way), inferred from the published tidal indices of the
// galaxies they disturb, by inverting the definition:
//     log10 M_k = Θ1(i) − C + 3 log10 D_ik,
// and taking the median over those galaxies. Returns the completed luminosities
// (as log10 L_K with M/L_K = 1) and a record of each inference.
export function inferDisturberLuminosities(positions, logL, publishedDisturber, publishedTheta1,
  { constant = TIDAL_C } = {}) {
  const estimates = new Map();
  for (let i = 0; i < logL.length; i++) {
    const k = publishedDisturber[i];
    if (k < 0 || Number.isFinite(logL[k]) || !Number.isFinite(publishedTheta1[i])) continue;
    const dx = positions[3 * k] - positions[3 * i];
    const dy = positions[3 * k + 1] - positions[3 * i + 1];
    const dz = positions[3 * k + 2] - positions[3 * i + 2];
    const d = Math.hypot(dx, dy, dz);
    if (d === 0) continue;
    if (!estimates.has(k)) estimates.set(k, []);
    estimates.get(k).push(publishedTheta1[i] - constant + 3 * Math.log10(d));
  }
  const filled = Float64Array.from(logL);
  const inferred = [];
  for (const [k, values] of estimates) {
    values.sort((a, b) => a - b);
    const mid = values.length >> 1;
    const median = values.length % 2 ? values[mid] : (values[mid - 1] + values[mid]) / 2;
    filled[k] = median;
    inferred.push({ index: k, logL: median, count: values.length, range: [values[0], values[values.length - 1]] });
  }
  return { logL: filled, inferred };
}

// Groups defined by tidal relations: a galaxy joins its main disturber when its
// tidal index exceeds the threshold. Returns, for every galaxy, the index of its
// group, and the groups sorted by size, each with its most luminous member.
export function tidalGroups(disturber, theta1, logL, threshold = 0) {
  const n = disturber.length;
  const parent = Int32Array.from({ length: n }, (_, i) => i);
  const find = (i) => {
    while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
  };
  for (let i = 0; i < n; i++) {
    const k = disturber[i];
    if (k >= 0 && theta1[i] > threshold) {
      const a = find(i), b = find(k);
      if (a !== b) parent[a] = b;
    }
  }
  const byRoot = new Map();
  for (let i = 0; i < n; i++) {
    const r = find(i);
    if (!byRoot.has(r)) byRoot.set(r, []);
    byRoot.get(r).push(i);
  }
  const groups = [...byRoot.values()].map((members) => {
    let leader = members[0];
    for (const m of members) {
      if ((Number.isFinite(logL[m]) ? logL[m] : -Infinity) > (Number.isFinite(logL[leader]) ? logL[leader] : -Infinity)) leader = m;
    }
    return { members, leader, size: members.length };
  }).sort((a, b) => b.size - a.size);
  const groupOf = new Int32Array(n);
  groups.forEach((g, index) => { for (const m of g.members) groupOf[m] = index; });
  return { groupOf, groups };
}

// Readable names from the catalogue's compact form: MESSIER081 -> M 81,
// NGC0253 -> NGC 253, ESO499-037 -> ESO 499-037; other names are kept.
const PREFIXES = ['UGCA', 'UGC', 'NGC', 'IC', 'PGC', 'DDO', 'KKH', 'KKR', 'KKs', 'KK', 'ESO', 'AGC', 'LEDA', 'CGCG', 'MCG'];
const PREFIX_PATTERN = new RegExp(`^(${PREFIXES.join('|')})0*(\\d.*)$`);

export function displayName(name) {
  const messier = /^MESSIER0*(\d+)$/.exec(name);
  if (messier) return `M ${messier[1]}`;
  const m = PREFIX_PATTERN.exec(name);
  return m ? `${m[1]} ${m[2]}` : name;
}
