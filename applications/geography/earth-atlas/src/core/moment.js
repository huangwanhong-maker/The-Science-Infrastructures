// Earthquake mechanisms from moment tensors.
//
// Provisional stand-in for core/domains/geophysics (seismology; laboratory ADR 0024).
//
// The moment tensor M of a source is symmetric; the catalogue gives it in the
// basis (r up, t south, p east). Its eigenvectors are the T (tension), N
// (null) and P (pressure) axes, in decreasing order of eigenvalue. For a
// double couple, M0 (n d^T + d n^T) with |n| = |d| = 1 and n . d = 0, the
// eigenvalues are M0, 0, -M0 and T = (n + d)/√2, P = (n - d)/√2; so the fault
// plane has normal (T + P)/√2 and slip (T - P)/√2, and the auxiliary plane the
// two exchanged. A P wave leaving the source in the direction x starts with a
// push (compression) where x^T M x > 0: the dark quadrants of a beach ball.
//
// Conventions: directions in (north, east, down); strike clockwise from north
// with the fault dipping to the right of it, dip from the horizontal, rake from
// the strike direction in the fault plane, positive for the hanging wall going
// up (Aki and Richards 2002, §4.2, Fig. 4.13). Projections are of the lower
// hemisphere, equal-area (Lambert, the Schmidt net), north up.
//
// Sources: Aki, K. and Richards, P. G. (2002), Quantitative Seismology, 2nd ed.,
// §3.3 and §4.2; Jost, M. L. and Herrmann, R. B. (1989), A student's guide to
// and review of moment tensors, Seismol. Res. Lett. 60, 37-57; Frohlich, C.
// (1992), Triangle diagrams: ternary graphs to display similarity and diversity
// of earthquake focal mechanisms, Phys. Earth Planet. Inter. 75, 193-198.

const DEG = Math.PI / 180;

// The tensor in (north, east, down) from the catalogue's (r up, t south, p east):
// N = -t, E = p, D = -r.
export function tensorNED({ mrr, mtt, mpp, mrt, mrp, mtp }) {
  return [
    [mtt, -mtp, mrt],
    [-mtp, mpp, -mrp],
    [mrt, -mrp, mrr],
  ];
}

// Eigenvalues and eigenvectors of a symmetric 3 x 3 matrix by Jacobi rotations:
// [{ value, vector }] in decreasing order of value.
export function eigenSymmetric(m) {
  const a = m.map((r) => r.slice()), v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
    if (off < 1e-15) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      if (Math.abs(a[p][q]) < 1e-18) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 3; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < 3; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
      }
      for (let k = 0; k < 3; k++) {
        const vkp = v[k][p], vkq = v[k][q];
        v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
      }
    }
  }
  return [0, 1, 2].map((i) => ({ value: a[i][i], vector: [v[0][i], v[1][i], v[2][i]] })).sort((x, y) => y.value - x.value);
}

// Plunge (degrees below the horizontal) and azimuth (degrees from north) of an
// axis in (north, east, down), taken pointing downward.
export function plungeAzimuth([n, e, d]) {
  if (d < 0) { n = -n; e = -e; d = -d; }
  const h = Math.hypot(n, e);
  return { plunge: Math.atan2(d, h) / DEG, azimuth: ((Math.atan2(e, n) / DEG) + 360) % 360 };
}

// Strike, dip and rake (degrees) of a plane with normal n and slip d (NED).
export function strikeDipRake(n, d) {
  if (n[2] > 0) { n = n.map((x) => -x); d = d.map((x) => -x); }       // the normal up, into the hanging wall
  const dip = Math.acos(Math.max(-1, Math.min(1, -n[2])));
  const strike = Math.atan2(-n[0], n[1]);
  const sd = Math.sin(dip);
  const rake = sd < 1e-9 ? Math.atan2(-d[0] * Math.sin(strike) + d[1] * Math.cos(strike), d[0] * Math.cos(strike) + d[1] * Math.sin(strike))
    : Math.atan2(-d[2] / sd, d[0] * Math.cos(strike) + d[1] * Math.sin(strike));
  return { strike: ((strike / DEG) + 360) % 360, dip: dip / DEG, rake: rake / DEG };
}

// The mechanism of a tensor given by its six elements (any scale): axes, planes,
// the fraction of it that is not a double couple, and its class.
export function mechanism(elements) {
  const m = tensorNED(elements);
  const [T, N, P] = eigenSymmetric(m);
  const t = T.vector, p = P.vector;
  const s = Math.SQRT1_2;
  const n1 = t.map((x, i) => s * (x + p[i])), d1 = t.map((x, i) => s * (x - p[i]));
  const planes = [strikeDipRake(n1, d1), strikeDipRake(d1, n1)];
  const axes = { T: plungeAzimuth(t), N: plungeAzimuth(N.vector), P: plungeAzimuth(p) };
  // epsilon = -lambda_min / |lambda_max| of absolute value, by deviatoric eigenvalue (Jost and Herrmann 1989, eq. 23).
  const values = [T.value, N.value, P.value];
  const big = Math.max(Math.abs(T.value), Math.abs(P.value));
  const epsilon = big ? N.value / big : 0;
  return { axes, planes, values, clvd: Math.abs(epsilon) * 2, kind: kindOf(axes) };
}

// Frohlich's (1992) classes by the plunges of the axes: thrust when T is
// steeper than 50°, normal when P is steeper than 60°, strike-slip when N is
// steeper than 60°; otherwise oblique.
export function kindOf(axes) {
  if (axes.N.plunge > 60) return 'strike-slip';
  if (axes.T.plunge > 50) return 'thrust';
  if (axes.P.plunge > 60) return 'normal';
  return 'oblique';
}

// x^T M x for a direction x (north, east, down) and a tensor M in NED: positive
// where the first motion of P waves is a push (compression).
export function radiation(m, x) {
  let s = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) s += x[i] * m[i][j] * x[j];
  return s;
}

// A direction (NED) in the lower hemisphere from a point (x east, y north) of
// the unit disc of the equal-area projection, or null outside it.
export function fromEqualArea(x, y) {
  const r = Math.hypot(x, y);
  if (r > 1) return null;
  const theta = 2 * Math.asin(Math.min(1, r * Math.SQRT1_2));       // from the downward vertical
  const az = Math.atan2(x, y);
  return [Math.sin(theta) * Math.cos(az), Math.sin(theta) * Math.sin(az), Math.cos(theta)];
}

// The point (x east, y north) of the unit disc for a direction (NED), lower hemisphere.
export function toEqualArea([n, e, d]) {
  if (d < 0) { n = -n; e = -e; d = -d; }
  const len = Math.hypot(n, e, d);
  const theta = Math.acos(Math.max(-1, Math.min(1, d / len)));
  const r = Math.SQRT2 * Math.sin(theta / 2);
  const h = Math.hypot(n, e) || 1;
  return [r * (e / h), r * (n / h)];
}

// Points of the great circle of a plane (strike, dip) on the projection, from
// one end of its strike to the other.
export function planeTrace({ strike, dip }, steps = 60) {
  const s = strike * DEG, dp = dip * DEG;
  const along = [Math.cos(s), Math.sin(s), 0];                                   // the strike direction
  const down = [-Math.sin(s) * Math.cos(dp), Math.cos(s) * Math.cos(dp), Math.sin(dp)];   // down the dip
  const out = [];
  for (let k = 0; k <= steps; k++) {
    const a = (k / steps) * Math.PI;
    out.push(toEqualArea(along.map((x, i) => Math.cos(a) * x + Math.sin(a) * down[i])));
  }
  return out;
}

// The angle (degrees) between two planes given by strike, dip and rake, counting
// a plane and its reverse (strike + 180, dip, -rake) as the same.
export function planeDifference(a, b) {
  const normal = ({ strike, dip }) => { const s = strike * DEG, d = dip * DEG; return [-Math.sin(d) * Math.sin(s), Math.sin(d) * Math.cos(s), -Math.cos(d)]; };
  const na = normal(a), nb = normal(b);
  const c = Math.abs(na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2]);
  return Math.acos(Math.min(1, c)) / DEG;
}
