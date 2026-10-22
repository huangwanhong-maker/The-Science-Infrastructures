// Rotations between the ICRS, the ecliptic, the mean equator of date, the
// galactic frame and the local horizon; topocentric observers; refraction.
//
// Provisional stand-in for core/metrology/frames (laboratory ADR 0011).
//
// Conventions: matrices are row-major arrays of nine numbers and describe
// passive rotations of the coordinate frame, as in SOFA: v' = M v.
//
// Sources:
//   IAU 2006 obliquity of the ecliptic at J2000: 84381.406 arcsec.
//   Lieske et al. (1977); Meeus, Astronomical Algorithms (1998), eq. 21.2–21.4: precession.
//   ESA (1997), The Hipparcos Catalogue, vol. 1, sec. 1.5.3: ICRS to galactic.
//   Meeus (1998), ch. 13: horizontal coordinates. Saemundsson (1986): refraction.
//   NIMA TR8350.2: WGS 84 ellipsoid.

export const DEG = Math.PI / 180;
export const ARCSEC = Math.PI / 648000;
export const AU_KM = 149597870.7;                 // IAU 2012, exact
export const LIGHT_KM_S = 299792.458;             // exact
export const OBLIQUITY_J2000 = 84381.406 * ARCSEC;

// ---- 3 x 3 rotations ----------------------------------------------------

export function rotX(a) {
  const c = Math.cos(a), s = Math.sin(a);
  return [1, 0, 0, 0, c, s, 0, -s, c];
}

export function rotY(a) {
  const c = Math.cos(a), s = Math.sin(a);
  return [c, 0, -s, 0, 1, 0, s, 0, c];
}

export function rotZ(a) {
  const c = Math.cos(a), s = Math.sin(a);
  return [c, s, 0, -s, c, 0, 0, 0, 1];
}

export function multiply(a, b) {
  const m = new Array(9);
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      m[3 * i + j] = a[3 * i] * b[j] + a[3 * i + 1] * b[3 + j] + a[3 * i + 2] * b[6 + j];
    }
  }
  return m;
}

export function transpose(m) {
  return [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
}

export function apply(m, v) {
  return [
    m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
    m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
    m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
  ];
}

// ---- vectors -------------------------------------------------------------

export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const scale = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);
export const unit = (a) => scale(a, 1 / norm(a));

// Unit vector from spherical angles (radians).
export function fromSpherical(lon, lat) {
  const c = Math.cos(lat);
  return [c * Math.cos(lon), c * Math.sin(lon), Math.sin(lat)];
}

// Spherical angles (radians) of a vector; longitude in [0, 2π).
export function toSpherical(v) {
  const r = norm(v);
  let lon = Math.atan2(v[1], v[0]);
  if (lon < 0) lon += 2 * Math.PI;
  return { lon, lat: Math.asin(v[2] / r), r };
}

// ---- fixed frames ---------------------------------------------------------

// ICRS (equator and equinox of J2000) to the ecliptic of J2000; frame bias ignored (< 0.03").
export const ICRS_TO_ECLIPTIC = rotX(OBLIQUITY_J2000);
export const ECLIPTIC_TO_ICRS = transpose(ICRS_TO_ECLIPTIC);

// ICRS to galactic coordinates (Hipparcos, ESA 1997).
export const ICRS_TO_GALACTIC = [
  -0.0548755604162154, -0.8734370902348850, -0.4838350155487132,
  +0.4941094278755837, -0.4448296299600112, +0.7469822444972189,
  -0.8676661490190047, -0.1980763734312015, +0.4559837761750669,
];

// Galactic to supergalactic coordinates (de Vaucouleurs, de Vaucouleurs and
// Corwin 1976): the north supergalactic pole lies at galactic (l, b) =
// (47.37°, +6.32°), and supergalactic longitude is zero at (137.37°, 0°).
// The two directions are orthogonal, so the rows below are orthonormal.
const SG_POLE = fromSpherical(47.37 * DEG, 6.32 * DEG);
const SG_ORIGIN = fromSpherical(137.37 * DEG, 0);
export const GALACTIC_TO_SUPERGALACTIC = [...SG_ORIGIN, ...cross(SG_POLE, SG_ORIGIN), ...SG_POLE];
export const ICRS_TO_SUPERGALACTIC = multiply(GALACTIC_TO_SUPERGALACTIC, ICRS_TO_GALACTIC);

// ---- precession -----------------------------------------------------------

// Precession from J2000 to the mean equator and equinox of date (IAU 1976).
// Accurate to well under 1" within a few centuries of J2000.
export function precessionMatrix(jdTT) {
  const t = (jdTT - 2451545.0) / 36525;
  const zeta = (2306.2181 * t + 0.30188 * t * t + 0.017998 * t ** 3) * ARCSEC;
  const z = (2306.2181 * t + 1.09468 * t * t + 0.018203 * t ** 3) * ARCSEC;
  const theta = (2004.3109 * t - 0.42665 * t * t - 0.041833 * t ** 3) * ARCSEC;
  return multiply(rotZ(-z), multiply(rotY(theta), rotZ(-zeta)));
}

// ---- observers and the horizon -------------------------------------------

const WGS84_A = 6378.137;               // km
const WGS84_F = 1 / 298.257223563;
const WGS84_E2 = WGS84_F * (2 - WGS84_F);

// Terrestrial position of a geodetic site, km (WGS 84).
export function geodeticToTerrestrial(latRad, lonRad, heightKm = 0) {
  const s = Math.sin(latRad), c = Math.cos(latRad);
  const n = WGS84_A / Math.sqrt(1 - WGS84_E2 * s * s);
  return [
    (n + heightKm) * c * Math.cos(lonRad),
    (n + heightKm) * c * Math.sin(lonRad),
    (n * (1 - WGS84_E2) + heightKm) * s,
  ];
}

// Terrestrial vector to the mean equator of date, rotating by GMST
// (nutation and polar motion ignored: < 20").
export function terrestrialToMeanOfDate(v, gmstRad) {
  return apply(rotZ(-gmstRad), v);
}

// Horizontal coordinates of a direction given in the mean equator of date.
// Azimuth from north through east, radians.
export function horizontal(vDate, lstRad, latRad) {
  const h = apply(rotY(Math.PI / 2 - latRad), apply(rotZ(lstRad), vDate));
  const r = norm(h);
  const north = -h[0], east = h[1], up = h[2];
  let az = Math.atan2(east, north);
  if (az < 0) az += 2 * Math.PI;
  return { alt: Math.asin(up / r), az, r };
}

// Atmospheric refraction for a true altitude, radians (Saemundsson 1986;
// standard pressure and temperature). Zero far below the horizon.
export function refraction(altRad) {
  const h = altRad / DEG;
  if (h < -1.5) return 0;
  const arcmin = 1.02 / Math.tan((h + 10.3 / (h + 5.11)) * DEG);
  return Math.max(0, arcmin) / 60 * DEG;
}
