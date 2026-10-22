// Quantities derived from the variables of a dataset: wind speed and direction,
// and the propagation of ocean waves.
//
// Provisional stand-in for core/domains/meteorology and core/domains/oceanography
// (laboratory ADR 0014).
//
// Conventions:
//   u is the eastward and v the northward component of a velocity.
//   A meteorological wind direction is the direction FROM which the wind blows,
//   in degrees clockwise from north: a wind from the west (u > 0, v = 0) is 270.
//   Wave directions in WaveWatch III follow the same "from" convention.
//
// Sources:
//   WMO (2018), Guide to Instruments and Methods of Observation (WMO-No. 8),
//   vol. I, ch. 5: wind direction.
//   Deep-water dispersion of surface gravity waves, omega^2 = g k, gives the phase
//   speed c = g T / (2 pi) and the group speed c_g = c / 2 = g T / (4 pi), at which
//   wave energy travels: Holthuijsen (2007), Waves in Oceanic and Coastal Waters,
//   eq. 5.4.33 and 5.4.36. Valid where the depth exceeds half a wavelength.
//   Standard gravity: CGPM (1901).

export const STANDARD_GRAVITY = 9.80665;

export const windSpeed = (u, v) => Math.hypot(u, v);

export function windFromDirection(u, v) {
  if (u === 0 && v === 0) return NaN;
  return ((Math.atan2(-u, -v) * 180) / Math.PI + 360) % 360;
}

const POINTS = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'];
export const compassPoint = (degrees) => POINTS[Math.round((((degrees % 360) + 360) % 360) / 22.5) % 16];

export const deepWaterGroupSpeed = (periodS) => (STANDARD_GRAVITY * periodS) / (4 * Math.PI);
export const deepWaterWavelength = (periodS) => (STANDARD_GRAVITY * periodS * periodS) / (2 * Math.PI);

// The velocity at which wave energy travels, from the peak direction (degrees,
// "from") and the peak period (s): [u, v] in m/s, pointing where the waves go.
export function waveEnergyVelocity(fromDirection, periodS) {
  const c = deepWaterGroupSpeed(periodS);
  const a = (fromDirection * Math.PI) / 180;
  return [-c * Math.sin(a), -c * Math.cos(a)];
}

// Mean of directions (degrees) with weights, through their unit vectors.
export function meanDirection(directions, weights) {
  let x = 0, y = 0;
  directions.forEach((d, i) => {
    const w = weights ? weights[i] : 1;
    x += w * Math.sin((d * Math.PI) / 180);
    y += w * Math.cos((d * Math.PI) / 180);
  });
  return x === 0 && y === 0 ? NaN : ((Math.atan2(x, y) * 180) / Math.PI + 360) % 360;
}
