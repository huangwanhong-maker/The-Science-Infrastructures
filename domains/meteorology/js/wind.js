// The wind: its speed, and the direction from which it blows.
//
// Provisional JavaScript, standing in for the C++ planned for this domain (ADR 0014).
//
// Conventions:
//   u is the eastward and v the northward component of a velocity.
//   A meteorological wind direction is the direction FROM which the wind blows,
//   in degrees clockwise from north: a wind from the west (u > 0, v = 0) is 270.
//
// Source:
//   WMO (2018), Guide to Instruments and Methods of Observation (WMO-No. 8),
//   vol. I, ch. 5: wind direction.

export const windSpeed = (u, v) => Math.hypot(u, v);

export function windFromDirection(u, v) {
  if (u === 0 && v === 0) return NaN;
  return ((Math.atan2(-u, -v) * 180) / Math.PI + 360) % 360;
}
