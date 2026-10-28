// The figure of a body for drawing: an ellipsoid of revolution (the Earth, WGS 84) or a sphere
// (the Moon, 1737.4 km), and the conversions the renderer needs between geodetic coordinates and
// body-fixed Cartesian coordinates, in metres.
//
// These are the formulas of the geodesy domain's C++ (geodesy/figure/ellipsoid.hpp; Torge and
// Müller 2012, eq. 4.24), repeated here because the page draws at the frame rate and cannot ask
// the engine for every vertex; the checks hold the two equal. Heights are along the normal of the
// ellipsoid (geodetic); on a sphere geodetic and geocentric latitudes coincide.

export const FIGURES = {
  earth: { name: 'WGS 84', a: 6378137.0, f: 1 / 298.257223563 },
  moon: { name: 'Moon, IAU mean sphere', a: 1737400.0, f: 0 },
};

export function figure(nameOrFigure) {
  const f = typeof nameOrFigure === 'string' ? FIGURES[nameOrFigure] : nameOrFigure;
  if (!f) throw new Error(`no figure for ${nameOrFigure}`);
  const e2 = f.f * (2 - f.f);
  return { ...f, b: f.a * (1 - f.f), e2 };
}

const DEG = Math.PI / 180;

// Body-fixed x, y, z (m) of a geodetic longitude, latitude (radians) and height (m).
export function cartesian(fig, lon, lat, h, out = [0, 0, 0], offset = 0) {
  const sp = Math.sin(lat), cp = Math.cos(lat);
  const n = fig.a / Math.sqrt(1 - fig.e2 * sp * sp);
  out[offset] = (n + h) * cp * Math.cos(lon);
  out[offset + 1] = (n + h) * cp * Math.sin(lon);
  out[offset + 2] = (n * (1 - fig.e2) + h) * sp;
  return out;
}
export const cartesianDeg = (fig, lonDeg, latDeg, h, out, offset) => cartesian(fig, lonDeg * DEG, latDeg * DEG, h, out, offset);

// The outward normal at a geodetic longitude and latitude (radians).
export function normal(lon, lat, out = [0, 0, 0], offset = 0) {
  const cp = Math.cos(lat);
  out[offset] = cp * Math.cos(lon);
  out[offset + 1] = cp * Math.sin(lon);
  out[offset + 2] = Math.sin(lat);
  return out;
}

// Geodetic longitude, latitude (radians) and height of a body-fixed point: Vermeille (2002) with a
// Newton step, as geodesy's C++.
export function geodetic(fig, x, y, z) {
  const rho2 = x * x + y * y, rho = Math.sqrt(rho2);
  const lon = rho2 > 0 ? Math.atan2(y, x) : 0;
  if (fig.f === 0) return { lon, lat: Math.atan2(z, rho), h: Math.hypot(rho, z) - fig.a };
  const a2 = fig.a * fig.a, e2 = fig.e2, e4 = e2 * e2;
  const p = rho2 / a2, q = ((1 - e2) / a2) * z * z, r = (p + q - e4) / 6;
  let phi;
  if (r > 0 || q > 0) {
    const s = (e4 * p * q) / (4 * r * r * r);
    const t = Math.cbrt(1 + s + Math.sqrt(s * (2 + s)));
    const u = r * (1 + t + 1 / t);
    const v = Math.sqrt(u * u + e4 * q);
    const w = (e2 * (u + v - q)) / (2 * v);
    const k = Math.sqrt(u + v + w * w) - w;
    const d = (k * rho) / (k + e2);
    phi = 2 * Math.atan2(z, d + Math.hypot(d, z));
  } else {
    phi = Math.atan2(z, rho);
  }
  for (let i = 0; i < 2; i++) {
    const sp = Math.sin(phi), cp = Math.cos(phi), w2 = 1 - e2 * sp * sp, n = fig.a / Math.sqrt(w2);
    const f = rho * sp - z * cp - e2 * n * sp * cp;
    const df = rho * cp + z * sp - e2 * n * (cp * cp - sp * sp) - (e2 * e2 * n * sp * sp * cp * cp) / w2;
    if (df === 0) break;
    phi -= f / df;
  }
  const sp = Math.sin(phi), cp = Math.cos(phi), n = fig.a / Math.sqrt(1 - e2 * sp * sp);
  return { lon, lat: phi, h: rho * cp + z * sp - (fig.a * fig.a) / n };
}

// The first intersection of a ray (origin o, direction d, body-fixed) with the ellipsoid raised by
// a height h, or null: the ellipsoid scaled to the unit sphere, as Cesium's rayEllipsoid.
export function rayHit(fig, o, d, h = 0) {
  const ra = fig.a + h, rb = fig.b + h;
  const ox = o[0] / ra, oy = o[1] / ra, oz = o[2] / rb, dx = d[0] / ra, dy = d[1] / ra, dz = d[2] / rb;
  const A = dx * dx + dy * dy + dz * dz, B = 2 * (ox * dx + oy * dy + oz * dz), C = ox * ox + oy * oy + oz * oz - 1;
  const disc = B * B - 4 * A * C;
  if (disc < 0) return null;
  const sq = Math.sqrt(disc);
  let t = (-B - sq) / (2 * A);
  if (t < 0) t = (-B + sq) / (2 * A);
  if (t < 0) return null;
  return [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]];
}

// Whether a body-fixed point is hidden from a viewpoint by the ellipsoid shrunk by `inset` metres
// (the lowest terrain): Cesium's horizon test in the ellipsoid's scaled space (Ring 2013, "Horizon
// culling", cesium.com blog), conservative for any point above that shrunk ellipsoid.
export function occluded(fig, viewer, point, inset = 0) {
  const ra = fig.a - inset, rb = fig.b - inset;
  const vx = viewer[0] / ra, vy = viewer[1] / ra, vz = viewer[2] / rb;
  const px = point[0] / ra, py = point[1] / ra, pz = point[2] / rb;
  const vh2 = vx * vx + vy * vy + vz * vz - 1;
  if (vh2 <= 0) return false; // inside the occluder: nothing is hidden by it
  const tx = px - vx, ty = py - vy, tz = pz - vz;
  const tDotV = -(tx * vx + ty * vy + tz * vz);
  return tDotV > vh2 && (tDotV * tDotV) / (tx * tx + ty * ty + tz * tz) > vh2;
}

// A point standing for a set of points in the horizon test (Cesium's horizon culling point; Ring
// 2013, "Computing the horizon occlusion point", cesium.com blog): a point P along the set's mean
// direction such that, if P is hidden by the ellipsoid shrunk by `inset`, every point of the set is.
// In the shrunk ellipsoid's scaled space, a point at distance m from the centre and angle α from the
// direction d bounds P's distance by 1 / cos(α + β), cos β = 1/m; when α + β reaches a right angle
// no such P exists (the set wraps round the horizon) and null is returned: never cull it.
export function horizonCullingPoint(fig, points, inset = 0) {
  const ra = fig.a - inset, rb = fig.b - inset;
  let dx = 0, dy = 0, dz = 0;
  for (const p of points) { dx += p[0] / ra; dy += p[1] / ra; dz += p[2] / rb; }
  const dn = Math.hypot(dx, dy, dz);
  if (!(dn > 0)) return null;
  dx /= dn; dy /= dn; dz /= dn;
  let largest = 0;
  for (const p of points) {
    const x = p[0] / ra, y = p[1] / ra, z = p[2] / rb;
    let m2 = x * x + y * y + z * z;
    let m = Math.sqrt(m2);
    const ux = x / m, uy = y / m, uz = z / m;
    m2 = Math.max(1, m2);
    m = Math.max(1, m);
    const cosA = ux * dx + uy * dy + uz * dz;
    const sinA = Math.hypot(uy * dz - uz * dy, uz * dx - ux * dz, ux * dy - uy * dx);
    const cosB = 1 / m, sinB = Math.sqrt(m2 - 1) * cosB;
    const denominator = cosA * cosB - sinA * sinB;
    if (denominator <= 1e-12) return null;
    largest = Math.max(largest, 1 / denominator);
  }
  return [dx * largest * ra, dy * largest * ra, dz * largest * rb];
}
