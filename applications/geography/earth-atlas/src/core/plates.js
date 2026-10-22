// The motion of the tectonic plates. Each plate is taken to be rigid, so that its
// motion over the sphere is a rotation: an angular velocity ω, a vector through
// the Earth's centre whose direction is the plate's pole of rotation and whose
// length is its rate. A point r of the plate moves with velocity v = ω × r
// (Euler's theorem; Cox and Hart 1986, ch. 4). Angular velocities add as
// vectors, so the motion of plate B relative to plate A is ω_B − ω_A: the Lie
// algebra so(3) is a vector space even though the rotations it generates do not
// commute.
//
// Provisional stand-in for core/domains/geodesy (laboratory ADR 0024).
//
// Conventions: angular velocities in radians per million years, so that with r
// in kilometres v is in km/Myr, which is mm/yr; longitudes east and latitudes north,
// geodetic, on WGS 84 (NIMA 2000); velocities as east and north components.
//
// Sources: Argus, D. F., Gordon, R. G. and DeMets, C. (2011), Geologically
// current motion of 56 plates relative to the no-net-rotation reference frame,
// Geochem. Geophys. Geosyst. 12, Q11001, Table 1; Bird, P. (2003), An updated
// digital model of plate boundaries, G-cubed 4(3), 1027.

const DEG = Math.PI / 180;

// NNR-MORVEL56 (Argus et al. 2011, Table 1): [name, code, pole latitude °N,
// pole longitude °E, rate °/Myr, 95 % confidence of the rate °/Myr]. Each plate
// rotates anticlockwise about its pole, relative to the no-net-rotation frame.
// Lower-case codes are the 25 plates of MORVEL; upper-case codes are the 31
// plates of Bird (2003) added to them.
export const NNR_MORVEL56 = [
  ['Amur', 'am', 63.17, -122.82, 0.297, 0.020], ['Antarctica', 'an', 65.42, -118.11, 0.250, 0.008],
  ['Arabia', 'ar', 48.88, -8.49, 0.559, 0.016], ['Australia', 'au', 33.86, 37.94, 0.632, 0.017],
  ['Capricorn', 'cp', 44.44, 23.09, 0.608, 0.019], ['Caribbean', 'ca', 35.20, -92.62, 0.286, 0.023],
  ['Cocos', 'co', 26.93, -124.31, 1.198, 0.045], ['Eurasia', 'eu', 48.85, -106.50, 0.223, 0.009],
  ['India', 'in', 50.37, -3.29, 0.544, 0.010], ['Juan de Fuca', 'jf', -38.31, 60.04, 0.951, 0.256],
  ['Lwandle', 'lw', 51.89, -69.52, 0.286, 0.026], ['Macquarie', 'mq', 49.19, 11.05, 1.144, 0.274],
  ['Nazca', 'nz', 46.23, -101.06, 0.696, 0.029], ['North America', 'na', -4.85, -80.64, 0.209, 0.013],
  ['Nubia', 'nb', 47.68, -68.44, 0.292, 0.007], ['Pacific', 'pa', -63.58, 114.70, 0.651, 0.011],
  ['Philippine Sea', 'ps', -46.02, -31.36, 0.910, 0.050], ['Rivera', 'ri', 20.25, -107.29, 4.536, 0.630],
  ['Sandwich', 'sw', -29.94, -36.87, 1.362, 0.744], ['Scotia', 'sc', 22.52, -106.15, 0.146, 0.016],
  ['Somalia', 'sm', 49.95, -84.52, 0.339, 0.011], ['South America', 'sa', -22.62, -112.83, 0.109, 0.011],
  ['Sunda', 'su', 50.06, -95.02, 0.337, 0.020], ['Sur', 'sr', -32.50, -111.32, 0.107, 0.028],
  ['Yangtze', 'yz', 63.03, -116.62, 0.334, 0.013],
  ['Aegean Sea', 'AS', 19.43, 122.87, 0.124, 0.050], ['Altiplano', 'AP', -6.58, -83.98, 0.488, 0.053],
  ['Anatolia', 'AT', 40.11, 26.66, 1.210, 0.049], ['Balmoral Reef', 'BR', -63.74, 142.06, 0.490, 0.045],
  ['Banda Sea', 'BS', -1.49, 121.64, 2.475, 0.058], ['Birds Head', 'BH', -40.00, 100.50, 0.799, 0.089],
  ['Burma', 'BU', -6.13, -78.10, 2.229, 0.060], ['Caroline', 'CL', -72.78, 72.05, 0.607, 0.047],
  ['Conway Reef', 'CR', -20.40, 170.53, 3.923, 0.048], ['Easter', 'EA', 24.97, 67.53, 11.334, 0.081],
  ['Futuna', 'FT', -16.33, 178.07, 5.101, 0.047], ['Galapagos', 'GP', 2.53, 81.18, 5.487, 0.063],
  ['Juan Fernandez', 'JZ', 34.25, 70.74, 22.368, 0.087], ['Kermadec', 'KE', 39.99, 6.46, 2.347, 0.048],
  ['Manus', 'MN', -3.67, 150.27, 51.569, 0.061], ['Maoke', 'MO', 14.25, 92.67, 0.774, 0.049],
  ['Mariana', 'MA', 11.05, 137.84, 1.306, 0.070], ['Molucca Sea', 'MS', 2.15, -56.09, 3.566, 0.058],
  ['New Hebrides', 'NH', 0.57, -6.60, 2.469, 0.045], ["Niuafo'ou", 'NI', -3.29, -174.49, 3.314, 0.067],
  ['North Andes', 'ND', 17.73, -122.68, 0.116, 0.051], ['North Bismarck', 'NB', -45.04, 127.64, 0.856, 0.045],
  ['Okhotsk', 'OK', 30.30, -92.28, 0.229, 0.051], ['Okinawa', 'ON', 36.12, 137.92, 2.539, 0.050],
  ['Panama', 'PM', 31.35, -113.90, 0.317, 0.055], ['Shetland', 'SL', 50.71, -143.47, 0.268, 0.048],
  ['Solomon Sea', 'SS', -2.87, 130.62, 1.703, 0.049], ['South Bismarck', 'SB', 6.88, -31.89, 8.111, 0.091],
  ['Timor', 'TI', -4.44, 113.50, 1.864, 0.067], ['Tonga', 'TO', 25.87, 4.48, 8.942, 0.100],
  ['Woodlark', 'WL', 0.10, 128.52, 1.744, 0.049],
];

const POLES = new Map(NNR_MORVEL56.map(([name, code, lat, lon, rate, ci]) => [code, { name, code, lat, lon, rate, ci }]));

// PB2002's plates in NNR-MORVEL56. The 18 MORVEL plates that Bird also has are
// written in upper case by Bird and lower case by MORVEL; Bird's Africa is
// MORVEL's Nubia (MORVEL's Lwandle lies within it), his Somalia and Yangtze are
// MORVEL's sm and yz; the other 31 plates are Bird's own, with his codes. (A
// case-blind match would join MORVEL's Nubia, nb, to Bird's North Bismarck, NB.)
// MORVEL's Capricorn, Lwandle, Macquarie and Sur have no outline in PB2002;
// there the motion of the plate containing them is used.
export const PB2002_TO_MORVEL = {
  AF: 'nb', SO: 'sm', YA: 'yz', AM: 'am', AN: 'an', AR: 'ar', AU: 'au', CA: 'ca', CO: 'co', EU: 'eu', IN: 'in',
  JF: 'jf', NA: 'na', NZ: 'nz', PA: 'pa', PS: 'ps', RI: 'ri', SA: 'sa', SC: 'sc', SU: 'su', SW: 'sw',
};

export const morvelCode = (pb2002) => PB2002_TO_MORVEL[pb2002] || pb2002;
export const poleOf = (pb2002) => POLES.get(morvelCode(pb2002)) || null;

// The plate motion model of ITRF2020 (Altamimi, Métivier, Rebischung, Collilieux,
// Chanard and Altamimi 2023, Geophys. Res. Lett. 50, e2023GL106373), estimated
// from GNSS velocities in the frame of ITRF2020 itself: for the 13 plates it
// covers, the Cartesian angular velocity (ωx, ωy, ωz) in milliarcseconds a year,
// as transcribed in PROJ's data file ITRF2020 (MIT licence), and the model's
// origin rate bias (mm/yr), to be added to ω × r to predict ITRF2020 velocities.
// It describes the present motion in the stations' own frame, so it, and not the
// geological NNR-MORVEL56, is the plate to subtract from a GPS velocity: stations
// in stable interiors then move by a few tenths of a millimetre a year.
export const ITRF2020_PMM = {
  AM: [-0.131, -0.551, 0.837], AN: [-0.269, -0.312, 0.678], AR: [1.129, -0.146, 1.438], AU: [1.487, 1.175, 1.223],
  CA: [0.207, -1.422, 0.726], EU: [-0.085, -0.519, 0.753], IN: [1.137, 0.013, 1.444], NZ: [-0.327, -1.561, 1.605],
  NA: [0.045, -0.666, -0.098], AF: [0.090, -0.585, 0.717], PA: [-0.404, 1.021, -2.154], SA: [-0.261, -0.282, -0.157],
  SO: [-0.081, -0.719, 0.864],
};
export const ITRF2020_ORIGIN_RATE_BIAS = [0.37, 0.35, 0.74];          // mm/yr, X, Y, Z
const MAS_PER_YR = (Math.PI / 180 / 3600 / 1000) * 1e6;               // 1 mas/yr in rad/Myr

// The angular velocity (rad/Myr) of a PB2002 plate in ITRF2020's model, or null.
export function itrfOmega(pb2002) {
  const w = ITRF2020_PMM[pb2002];
  return w ? w.map((x) => x * MAS_PER_YR) : null;
}

// WGS 84: semi-major axis (m) and the square of the first eccentricity.
const A = 6378137, FLAT = 1 / 298.257223563, E2 = FLAT * (2 - FLAT);

// The Earth-centred, Earth-fixed position (m) of a point at geodetic longitude
// and latitude (degrees) and height h (m) on WGS 84.
export function ecef(lon, lat, h = 0) {
  const l = lon * DEG, p = lat * DEG;
  const N = A / Math.sqrt(1 - E2 * Math.sin(p) ** 2);
  return [(N + h) * Math.cos(p) * Math.cos(l), (N + h) * Math.cos(p) * Math.sin(l), (N * (1 - E2) + h) * Math.sin(p)];
}

// The angular velocity (rad/Myr) of a pole { lat, lon, rate (°/Myr) }.
export function omega(pole) {
  const r = pole.rate * DEG, p = pole.lat * DEG, l = pole.lon * DEG;
  return [r * Math.cos(p) * Math.cos(l), r * Math.cos(p) * Math.sin(l), r * Math.sin(p)];
}

// The pole { lat, lon, rate } of an angular velocity (rad/Myr).
export function poleOfOmega([x, y, z]) {
  const r = Math.hypot(x, y, z);
  if (r === 0) return { lat: 90, lon: 0, rate: 0 };
  return { lat: Math.asin(z / r) / DEG, lon: Math.atan2(y, x) / DEG, rate: r / DEG };
}

const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const subtract = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

// The velocity (mm/yr, east and north) of a point turning with angular velocity w (rad/Myr).
export function velocityAt(w, lon, lat) {
  const v = cross(w, ecef(lon, lat).map((x) => x / 1000));       // km x rad/Myr = mm/yr
  const l = lon * DEG, p = lat * DEG;
  return {
    east: -Math.sin(l) * v[0] + Math.cos(l) * v[1],
    north: -Math.sin(p) * Math.cos(l) * v[0] - Math.sin(p) * Math.sin(l) * v[1] + Math.cos(p) * v[2],
  };
}

// The velocity (mm/yr) that ITRF2020's model predicts for a point of a plate,
// origin rate bias included; null for a plate it does not cover.
export function itrfVelocityAt(pb2002, lon, lat) {
  const w = itrfOmega(pb2002);
  if (!w) return null;
  const v = velocityAt(w, lon, lat);
  const [tx, ty, tz] = ITRF2020_ORIGIN_RATE_BIAS, l = lon * DEG, p = lat * DEG;
  return {
    east: v.east - Math.sin(l) * tx + Math.cos(l) * ty,
    north: v.north - Math.sin(p) * Math.cos(l) * tx - Math.sin(p) * Math.sin(l) * ty + Math.cos(p) * tz,
  };
}

// The motion of the right-hand plate relative to the left-hand one across a
// boundary step of azimuth `azimuth` (degrees from north, the direction of
// travel along the step), as Bird (2003) gives it: the divergent component
// (positive opening, negative convergence) and the right-lateral component.
// v: the velocity of the right-hand plate relative to the left (east, north).
export function acrossStep(v, azimuth) {
  const a = azimuth * DEG;
  const t = [Math.sin(a), Math.cos(a)];          // along the step
  const n = [Math.cos(a), -Math.sin(a)];         // to the right of it, into the right-hand plate
  const along = v.east * t[0] + v.north * t[1];
  // Seen from either plate, the other moves to the observer's right when the
  // right-hand plate moves against the direction of travel.
  return { divergence: v.east * n[0] + v.north * n[1], rightLateral: -along };
}

// The initial azimuth (degrees from north) of the great circle from (lon1, lat1) to (lon2, lat2).
export function azimuth(lon1, lat1, lon2, lat2) {
  const p1 = lat1 * DEG, p2 = lat2 * DEG, dl = (lon2 - lon1) * DEG;
  const y = Math.sin(dl) * Math.cos(p2), x = Math.cos(p1) * Math.sin(p2) - Math.sin(p1) * Math.cos(p2) * Math.cos(dl);
  return ((Math.atan2(y, x) / DEG) + 360) % 360;
}

// Whether a point lies within a plate's rings (even-odd, in longitude and
// latitude: the outlines are cut at the antimeridian and run along the poles).
export function inRings(coords, rings, lon, lat) {
  let inside = false;
  for (const [first, count] of rings) {
    for (let k = 0, j = count - 1; k < count; j = k++) {
      const xa = coords[2 * (first + k)], ya = coords[2 * (first + k) + 1];
      const xb = coords[2 * (first + j)], yb = coords[2 * (first + j) + 1];
      if ((ya > lat) !== (yb > lat) && lon < ((xb - xa) * (lat - ya)) / (yb - ya) + xa) inside = !inside;
    }
  }
  return inside;
}
