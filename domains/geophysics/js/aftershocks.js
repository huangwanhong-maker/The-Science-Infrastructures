// Aftershocks: the rate and probability of aftershocks after a mainshock, by
// the model of Reasenberg and Jones (1989), and the size of the zone in which
// they fall, from the length of the rupture (Wells and Coppersmith 1994).
//
// Provisional JavaScript, standing in for the C++ planned for this domain (seismology; ADR 0024).
//
// The model joins two empirical laws: the Gutenberg-Richter law (the number
// of earthquakes of magnitude ≥ M falls by 10^b per unit of M) and the
// modified Omori law (the rate of aftershocks decays as (t + c)^-p). The rate
// of aftershocks of magnitude ≥ M, at time t (days) after a mainshock of
// magnitude Mm, is (Reasenberg and Jones 1989)
//   λ(t, M) = 10^(a + b (Mm - M)) (t + c)^-p,
// and the number expected in [t1, t2] is its integral,
//   N = 10^(a + b (Mm - M)) [(t1 + c)^(1-p) - (t2 + c)^(1-p)] / (p - 1)   (p ≠ 1),
//   N = 10^(a + b (Mm - M)) ln[(t2 + c) / (t1 + c)]                       (p = 1).
// Aftershocks being a nonhomogeneous Poisson process, the probability of at
// least one is P = 1 - exp(-N). Their generic parameters, fitted to 62
// sequences of mainshocks of magnitude ≥ 5 in California from 1933 to 1987,
// are a = -1.67, b = 0.91, p = 1.08, c = 0.05 day.
//
// Range of validity: California's sequences; elsewhere productivity and decay
// differ by region (subduction zones and mid-ocean ridges among them), and
// from sequence to sequence by a factor of ten in rate (Page et al. 2016).
// Early in a sequence the catalogue misses small aftershocks, which the model
// counts; for M near or above Mm, P is the chance of a later, larger
// earthquake, which the model includes as an aftershock.
//
// Wells and Coppersmith (1994), Table 2A: subsurface rupture length, all slip
// types, log10 RLD (km) = -2.44 + 0.59 M, with standard deviation 0.16 in
// log10 RLD, for 4.8 ≤ M ≤ 8.1 (moment magnitude).
//
// References:
//   Reasenberg, P. A. and Jones, L. M. (1989), Earthquake hazard after a
//   mainshock in California, Science 243, 1173-1176.
//   Wells, D. L. and Coppersmith, K. J. (1994), New empirical relationships
//   among magnitude, rupture length, rupture width, rupture area, and surface
//   displacement, Bull. Seismol. Soc. Am. 84, 974-1002.
//   Page, M. T., van der Elst, N., Hardebeck, J., Felzer, K. and Michael, A. J.
//   (2016), Three ingredients for improved global aftershock forecasts:
//   tectonic region, time-dependent catalog incompleteness, and intersequence
//   variability, Bull. Seismol. Soc. Am. 106, 2290-2301.

export const GENERIC_CALIFORNIA = { a: -1.67, b: 0.91, p: 1.08, c: 0.05 };

// The rate (per day) of aftershocks of magnitude ≥ M at t days after a mainshock of magnitude Mm.
export function aftershockRate(t, M, Mm, { a, b, p, c } = GENERIC_CALIFORNIA) {
  return 10 ** (a + b * (Mm - M)) * (t + c) ** -p;
}

// The number of aftershocks of magnitude ≥ M expected between t1 and t2 days after the mainshock.
export function expectedAftershocks(t1, t2, M, Mm, { a, b, p, c } = GENERIC_CALIFORNIA) {
  const k = 10 ** (a + b * (Mm - M));
  if (Math.abs(p - 1) < 1e-12) return k * Math.log((t2 + c) / (t1 + c));
  return (k * ((t1 + c) ** (1 - p) - (t2 + c) ** (1 - p))) / (p - 1);
}

// The probability of at least one aftershock of magnitude ≥ M between t1 and t2 days.
export const aftershockProbability = (t1, t2, M, Mm, params = GENERIC_CALIFORNIA) => 1 - Math.exp(-expectedAftershocks(t1, t2, M, Mm, params));

// The subsurface rupture length (km) of an earthquake of moment magnitude M, Wells and Coppersmith (1994).
export const subsurfaceRuptureLengthKm = (M) => 10 ** (-2.44 + 0.59 * M);
