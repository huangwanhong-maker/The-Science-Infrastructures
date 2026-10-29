// Colour-vision deficiencies simulated, so that a map's colours are validated for the readers who see
// them otherwise (the infrastructure's ADR 0053; CLAUDE.md §7, visual style).
//
// Machado, Oliveira and Fernandes (2009), "A physiologically-based model for simulation of color vision
// deficiency", IEEE TVCG 15(6), 1291–1298: a matrix on linear sRGB for each deficiency; at severity 1,
// protanopia, deuteranopia and tritanopia (their Table, as colorspacious 1.1.2 tabulates it; the checks
// hold these to colorspacious's own simulation, tests/fixtures/cvd.json, made by
// tests/make_cvd_fixtures.py). The simulated colour is not clamped, as colorspacious's is not; for
// display, clamp it.
//
// A map is validated by what its use needs: a sequential map's lightness must stay monotone under each
// deficiency (its order then reads without hue), the two arms of a diverging map must stay apart, and
// the adjacent hues of a categorical map must differ by a noticeable step.

import { deltaE, oklab, toGamma, toLinear } from './colour.js';

export const MACHADO_2009 = {
  protan: [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
  deutan: [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.01182, 0.04294, 0.968881]],
  tritan: [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.3039]],
};
export const DEFICIENCIES = Object.keys(MACHADO_2009);

// The colour as a reader with the deficiency `kind` sees it (unclamped sRGB). Linear values below the
// transfer function's knee, negative ones among them, take its linear part, as colorspacious's do.
export function simulate(rgb, kind) {
  const M = MACHADO_2009[kind];
  if (!M) throw new Error(`no deficiency "${kind}": one of ${DEFICIENCIES.join(', ')}`);
  const lin = rgb.map(toLinear);
  const sim = M.map((row) => row[0] * lin[0] + row[1] * lin[1] + row[2] * lin[2]);
  return sim.map(toGamma);
}

const clamp01 = (rgb) => rgb.map((x) => Math.min(1, Math.max(0, x)));

// The lightness of a colour function along u, under each deficiency and in normal vision: whether it
// is monotone (in one direction) and by how much its least step falls short.
export function lightnessProfile(colour, { steps = 64 } = {}) {
  const out = {};
  for (const kind of ['normal', ...DEFICIENCIES]) {
    const L = Array.from({ length: steps }, (_, k) => {
      const c = colour(k / (steps - 1));
      return oklab(kind === 'normal' ? c : clamp01(simulate(c, kind)))[0];
    });
    const rising = L[L.length - 1] >= L[0];
    let worst = Infinity;
    for (let k = 1; k < steps; k++) worst = Math.min(worst, rising ? L[k] - L[k - 1] : L[k - 1] - L[k]);
    out[kind] = { monotone: worst >= 0, least_step: worst, from: L[0], to: L[L.length - 1] };
  }
  return out;
}

// The least difference between adjacent colours of a list, in normal vision and under each deficiency.
export function separation(colours) {
  const out = {};
  for (const kind of ['normal', ...DEFICIENCIES]) {
    const seen = colours.map((c) => (kind === 'normal' ? c : clamp01(simulate(c, kind))));
    let least = Infinity;
    for (let i = 1; i < seen.length; i++) least = Math.min(least, deltaE(seen[i - 1], seen[i]));
    out[kind] = least;
  }
  return out;
}
