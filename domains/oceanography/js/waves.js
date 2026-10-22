// The propagation of ocean waves: the speed and direction in which their
// energy travels.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).
//
// Conventions:
//   u is the eastward and v the northward component of a velocity.
//   Wave directions in WaveWatch III are the directions FROM which the waves come,
//   in degrees clockwise from north, as wind directions are.
//
// Sources:
//   Deep-water dispersion of surface gravity waves, omega^2 = g k, gives the phase
//   speed c = g T / (2 pi) and the group speed c_g = c / 2 = g T / (4 pi), at which
//   wave energy travels: Holthuijsen (2007), Waves in Oceanic and Coastal Waters,
//   eq. 5.4.33 and 5.4.36. Valid where the depth exceeds half a wavelength.
//   Standard gravity: CGPM (1901).

export const STANDARD_GRAVITY = 9.80665;

export const deepWaterGroupSpeed = (periodS) => (STANDARD_GRAVITY * periodS) / (4 * Math.PI);
export const deepWaterWavelength = (periodS) => (STANDARD_GRAVITY * periodS * periodS) / (2 * Math.PI);

// The velocity at which wave energy travels, from the peak direction (degrees,
// "from") and the peak period (s): [u, v] in m/s, pointing where the waves go.
export function waveEnergyVelocity(fromDirection, periodS) {
  const c = deepWaterGroupSpeed(periodS);
  const a = (fromDirection * Math.PI) / 180;
  return [-c * Math.sin(a), -c * Math.cos(a)];
}
