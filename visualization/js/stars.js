// The colour of a star from its colour index, for drawing skies: from the Johnson B−V index or
// Gaia's BP−RP, through the effective temperature, to the sRGB of a blackbody at it.
//
// Sources:
//   Ballesteros (2012), EPL 97, 34008: effective temperature from B−V.
//   Mucciarelli and Bellazzini (2020), RNAAS 4, 52: effective temperature from BP−RP, for dwarfs,
//     without reddening correction, solar metallicity; calibrated for 0.39 < BP−RP < 1.5 and
//     extrapolated outside.
//   Helland (2012), "How to convert temperature (K) to RGB": a fit to the blackbody's colour.
//
// First written for the Celestial Atlas, and moved here when the globes' sky needed it.

const clamp = (x, lo, hi) => Math.min(hi, Math.max(lo, x));

export function bvToTemperature(bv) {
  return 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62));
}

// Approximate sRGB of a blackbody of temperature T (K), components in [0, 1].
export function temperatureToRgb(T) {
  const t = clamp(T, 1000, 40000) / 100;
  const r = t <= 66 ? 255 : 329.698727446 * (t - 60) ** -0.1332047592;
  const g = t <= 66 ? 99.4708025861 * Math.log(t) - 161.1195681661
    : 288.1221695283 * (t - 60) ** -0.0755148492;
  const b = t >= 66 ? 255 : t <= 19 ? 0 : 138.5177312231 * Math.log(t - 10) - 305.0447927307;
  return [clamp(r, 0, 255) / 255, clamp(g, 0, 255) / 255, clamp(b, 0, 255) / 255];
}

// Star colour from B−V; unknown indices are drawn as a Sun-like star.
export function starRgb(bv) {
  const index = Number.isFinite(bv) ? clamp(bv, -0.4, 2.0) : 0.65;
  return temperatureToRgb(bvToTemperature(index));
}

export function bpRpToTemperature(bpRp) {
  const c = clamp(bpRp, -0.5, 5);
  return 5040 / (0.4929 + 0.5092 * c - 0.0353 * c * c);
}

// Star colour for a catalogue's photometric system: 'johnson' (B−V) or 'gaia' (BP−RP).
export function starColour(photometry, index) {
  if (photometry === 'gaia') {
    return temperatureToRgb(Number.isFinite(index) ? bpRpToTemperature(index) : 5800);
  }
  return starRgb(index);
}
